"""
Principled Small-N Treatment for Mahalanobis S_MD
=================================================

Problem:
  Pure Mahalanobis uses the empirical Lor4D centroid μ(N) for *all* features.
  At very small N (notably N=16), the Lor4D vs Lor5D separation in d_eff is
  weak; μ_d(N) can drift away from the physical value d*=4, allowing Lor5D
  to intrude.

"Principled" fix (no tunable hyperparameters):
  Anchor ONLY the dimension component of the reference manifold to the
  physical prior d*=4, while keeping (c,w) empirical:

    μ_anchor(N) = (4.0, μ_c(N), μ_w(N))

  Keep Σ(N) from the Lor4D ensemble, so the metric remains data-driven.

This yields an "anchored Mahalanobis" action:
  S_MD^anchor = (I - μ_anchor)^T Σ^{-1} (I - μ_anchor)

Goal:
  Test whether this removes the N=16 Lor5D intruder and restores Lor4D #1,
  while preserving the strong large-N behavior.
"""
from __future__ import annotations

import time
from collections import defaultdict
from pathlib import Path

import numpy as np

from expanded_family_robustness import ALL_FAMILIES, compute_features


def mahalanobis_score(feat: np.ndarray, mu: np.ndarray, cov_inv: np.ndarray) -> float:
    delta = feat - mu
    return float(delta @ cov_inv @ delta)


def oas_shrinkage_cov(X: np.ndarray) -> np.ndarray:
    """Oracle Approximating Shrinkage (OAS) covariance toward scaled identity.

    Parameter-free shrinkage estimator, useful when sample size is small.
    Reference: Chen et al. (2010). Here implemented for small p (p=3).
    """
    X = np.asarray(X, dtype=float)
    n, p = X.shape
    if n <= 1:
        return np.eye(p)
    Xc = X - np.mean(X, axis=0, keepdims=True)
    S = (Xc.T @ Xc) / float(n)  # biased covariance
    trS = float(np.trace(S))
    if trS <= 0:
        return np.eye(p)
    mu = trS / float(p)
    F = mu * np.eye(p)
    # Compute shrinkage intensity alpha
    # phi = (1/n) * sum ||x_i x_i^T - S||_F^2
    # For efficiency with small p, compute explicitly.
    phi = 0.0
    for i in range(n):
        xi = Xc[i : i + 1].T  # (p,1)
        outer = xi @ xi.T
        diff = outer - S
        phi += float(np.sum(diff * diff))
    phi /= float(n)
    gamma = float(np.sum((S - F) * (S - F)))
    if gamma <= 1e-18:
        alpha = 1.0
    else:
        alpha = min(1.0, phi / (gamma * float(n)))
    return (1.0 - alpha) * S + alpha * F


def main() -> None:
    # Keep N=16 in scope, add a couple of near-boundary points.
    N_VALUES = [16, 20, 28, 48, 64, 96]
    REPS = 20
    SEED_BASE = 42

    n_fam = len(ALL_FAMILIES)
    total_samples = n_fam * len(N_VALUES) * REPS
    print("=" * 80)
    print("Mahalanobis Small-N Prior Test (Anchored d*=4)")
    print(f"  Families: {n_fam}")
    print(f"  N values: {N_VALUES}")
    print(f"  Reps per (family, N): {REPS}")
    print(f"  Total samples: {total_samples}")
    print("=" * 80)

    t0 = time.time()

    data: dict[int, dict[str, list[np.ndarray]]] = defaultdict(lambda: defaultdict(list))
    for fi, (fam_name, gen_fn) in enumerate(ALL_FAMILIES.items()):
        print(f"  [{fi+1:2d}/{n_fam}] {fam_name:15s}", end="", flush=True)
        for N in N_VALUES:
            ok = 0
            for rep in range(REPS):
                seed = SEED_BASE + hash(fam_name) % 10000 + N * 100 + rep
                seed = seed % (2**31)
                try:
                    poset = gen_fn(N, seed=seed)
                    feat = compute_features(poset, N)
                    data[N][fam_name].append(feat)
                    ok += 1
                except Exception:
                    pass
            print(f"  N={N}:{ok}", end="", flush=True)
        print()

    print(f"\nGeneration done in {time.time() - t0:.1f}s")

    results_baseline: dict[int, list[tuple[str, float]]] = {}
    results_anchor: dict[int, list[tuple[str, float]]] = {}
    results_oas: dict[int, list[tuple[str, float]]] = {}

    for N in N_VALUES:
        lor4d_feats = np.array(data[N]["Lor4D"])
        if len(lor4d_feats) == 0:
            print(f"  WARNING: no Lor4D data at N={N}")
            continue

        mu = np.mean(lor4d_feats, axis=0)
        cov = np.cov(lor4d_feats.T)
        cov_inv = np.linalg.inv(cov + 1e-12 * np.eye(3))
        cov_oas = oas_shrinkage_cov(lor4d_feats)
        cov_oas_inv = np.linalg.inv(cov_oas + 1e-12 * np.eye(3))

        mu_anchor = mu.copy()
        mu_anchor[0] = 4.0  # physical prior d*=4

        fam_scores_base = []
        fam_scores_anchor = []
        fam_scores_oas = []
        for fam_name in ALL_FAMILIES:
            feats = data[N][fam_name]
            if not feats:
                continue
            feats_arr = np.array(feats)
            avg_base = np.mean([mahalanobis_score(f, mu, cov_inv) for f in feats_arr])
            avg_anchor = np.mean([mahalanobis_score(f, mu_anchor, cov_inv) for f in feats_arr])
            avg_oas = np.mean([mahalanobis_score(f, mu, cov_oas_inv) for f in feats_arr])
            fam_scores_base.append((fam_name, avg_base))
            fam_scores_anchor.append((fam_name, avg_anchor))
            fam_scores_oas.append((fam_name, avg_oas))

        results_baseline[N] = sorted(fam_scores_base, key=lambda x: x[1])
        results_anchor[N] = sorted(fam_scores_anchor, key=lambda x: x[1])
        results_oas[N] = sorted(fam_scores_oas, key=lambda x: x[1])

    # ── Report ───────────────────────────────────────────────────────────
    report: list[str] = []
    report.append("# Mahalanobis Small-N Prior Test — Anchored d*=4\n")
    report.append(f"Date: {time.strftime('%Y-%m-%d %H:%M')}")
    report.append(f"Families: {n_fam}")
    report.append(f"N values: {N_VALUES}")
    report.append(f"Reps: {REPS}\n")

    report.append("## 1. Definition\n")
    report.append("- Baseline (full): `S = (I-mu)^T Sigma^{-1} (I-mu)` with mu,Sigma from Lor4D ensemble at each N.")
    report.append("- Anchored (full): same Sigma, but `mu_anchor = (4.0, mu_c, mu_w)` (physical d*=4 prior).")
    report.append("- Anchored (diag): assume feature independence, `Sigma^{-1} -> diag(1/var_i)` from Lor4D.\n")

    report.append("## 2. Lor4D Rank Summary\n")
    report.append("| N | Baseline winner | Lor4D rank | OAS winner | Lor4D rank | Anchored(full) winner | Lor4D rank | Anchored(diag) winner | Lor4D rank |")
    report.append("|---|---|---:|---|---:|---|---:|---|---:|")
    for N in N_VALUES:
        if N not in results_baseline:
            continue
        ranked_b = results_baseline[N]
        ranked_o = results_oas[N]
        ranked_a = results_anchor[N]
        # diag variant computed on demand from cached data
        lor4d_feats = np.array(data[N]["Lor4D"])
        mu = np.mean(lor4d_feats, axis=0)
        mu_anchor = mu.copy()
        mu_anchor[0] = 4.0
        var = np.var(lor4d_feats, axis=0, ddof=1)
        inv_var = 1.0 / np.maximum(var, 1e-12)
        fam_scores_diag = []
        for fam_name in ALL_FAMILIES:
            feats = data[N][fam_name]
            if not feats:
                continue
            feats_arr = np.array(feats)
            deltas = feats_arr - mu_anchor
            scores = np.sum(inv_var * (deltas ** 2), axis=1)
            fam_scores_diag.append((fam_name, float(np.mean(scores))))
        ranked_d = sorted(fam_scores_diag, key=lambda x: x[1])

        b_rank = next(i + 1 for i, (f, _) in enumerate(ranked_b) if f == "Lor4D")
        o_rank = next(i + 1 for i, (f, _) in enumerate(ranked_o) if f == "Lor4D")
        a_rank = next(i + 1 for i, (f, _) in enumerate(ranked_a) if f == "Lor4D")
        d_rank = next(i + 1 for i, (f, _) in enumerate(ranked_d) if f == "Lor4D")
        report.append(
            f"| {N} | {ranked_b[0][0]} | {b_rank} | {ranked_o[0][0]} | {o_rank} | {ranked_a[0][0]} | {a_rank} | {ranked_d[0][0]} | {d_rank} |"
        )

    report.append("\n## 3. N=16 Top-5 (Baseline vs Anchored)\n")
    if 16 in results_baseline:
        report.append("### Baseline Top 5 (N=16)\n")
        report.append("| Rank | Family | Score |")
        report.append("|---:|---|---:|")
        for i, (fam, score) in enumerate(results_baseline[16][:5], start=1):
            report.append(f"| {i} | {fam} | {score:.4f} |")
        report.append("\n### Anchored Top 5 (N=16)\n")
        report.append("| Rank | Family | Score |")
        report.append("|---:|---|---:|")
        for i, (fam, score) in enumerate(results_anchor[16][:5], start=1):
            report.append(f"| {i} | {fam} | {score:.4f} |")
        # Anchored-diagonal Top 5
        lor4d_feats = np.array(data[16]["Lor4D"])
        mu = np.mean(lor4d_feats, axis=0)
        mu_anchor = mu.copy()
        mu_anchor[0] = 4.0
        var = np.var(lor4d_feats, axis=0, ddof=1)
        inv_var = 1.0 / np.maximum(var, 1e-12)
        fam_scores_diag = []
        for fam_name in ALL_FAMILIES:
            feats = data[16][fam_name]
            if not feats:
                continue
            feats_arr = np.array(feats)
            deltas = feats_arr - mu_anchor
            scores = np.sum(inv_var * (deltas ** 2), axis=1)
            fam_scores_diag.append((fam_name, float(np.mean(scores))))
        ranked_d = sorted(fam_scores_diag, key=lambda x: x[1])
        report.append("\n### Anchored (diag) Top 5 (N=16)\n")
        report.append("| Rank | Family | Score |")
        report.append("|---:|---|---:|")
        for i, (fam, score) in enumerate(ranked_d[:5], start=1):
            report.append(f"| {i} | {fam} | {score:.4f} |")

    # Quick diagnostics: how far μ_d is from 4 at each N
    report.append("\n## 4. mu_d(N) Drift (Lor4D empirical centroid)\n")
    report.append("| N | mu_d | |mu_d-4| | mu_c | mu_w |")
    report.append("|---|---:|---:|---:|---:|")
    for N in N_VALUES:
        lor4d_feats = np.array(data[N]["Lor4D"])
        if len(lor4d_feats) == 0:
            continue
        mu = np.mean(lor4d_feats, axis=0)
        report.append(f"| {N} | {mu[0]:.4f} | {abs(mu[0]-4.0):.4f} | {mu[1]:.4f} | {mu[2]:.4f} |")

    report.append(f"\nRuntime: {time.time() - t0:.1f}s\n")

    full_report = "\n".join(report)
    print("\n" + full_report)

    outdir = Path(__file__).parent / "outputs_carlip"
    outdir.mkdir(exist_ok=True)
    outpath = outdir / "mahalanobis_small_n_prior_results.md"
    outpath.write_text(full_report, encoding="utf-8")
    print(f"\nSaved to {outpath}")


if __name__ == "__main__":
    main()
