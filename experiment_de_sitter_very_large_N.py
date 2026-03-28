"""
Very-large-N de Sitter convergence experiment.
Extends to N=768 and N=1024 to verify that ranking stability and
scaling exponent acceleration persist deep into the large-N regime.

Key optimizations:
  - 15 reps (sufficient for mean estimation)
  - 3 H values: 0.0 (control), 0.1, 0.3
  - Standard family reference cached per N
  - Per-N timing for diagnostics
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import linregress

from generators import Poset, transitive_closure
from curvature_layer2_de_sitter_scan import (
    sprinkle_de_sitter_like_diamond,
    poset_from_de_sitter_points,
)
from expanded_family_robustness import (
    ALL_FAMILIES,
    compute_features,
    mahalanobis_score,
)

SEED_BASE = 42
REPS = 15
N_GRID = [256, 512, 768, 1024]
H_GRID = [0.0, 0.1, 0.3]

OUT_DIR = Path("outputs_carlip")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def generate_ds4d(n, hubble, seed):
    pts = sprinkle_de_sitter_like_diamond(n, d_spatial=3, hubble=hubble, seed=seed)
    return poset_from_de_sitter_points(pts, hubble)


def feature_ensemble(gen_fn, N, reps, seed_base):
    feats = []
    for r in range(reps):
        poset = gen_fn(N, seed=seed_base + r)
        f = compute_features(poset, N)
        feats.append(f)
    return np.array(feats)


def ds4d_feature_ensemble(N, hubble, reps, seed_base):
    feats = []
    for r in range(reps):
        poset = generate_ds4d(N, hubble, seed=seed_base + r)
        f = compute_features(poset, N)
        feats.append(f)
    return np.array(feats)


def main():
    if hasattr(sys.stdout, "reconfigure"):
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    t0 = time.time()
    report = ["# Very-Large-N de Sitter Convergence Experiment\n"]
    report.append(f"N grid: {N_GRID}")
    report.append(f"H grid: {H_GRID}")
    report.append(f"REPS per (family, N): {REPS}\n")

    # ================================================================
    # Part 1: Rankings at very large N
    # ================================================================
    report.append("=" * 72)
    report.append("# Part 1: S_MD Rankings at Very Large N")
    report.append("=" * 72)
    report.append("")
    report.append("| N | H | S_MD(dS4D) | S_MD(Lor4D) | Delta | Lor4D_rank | dS4D_rank | total | t_N(s) |")
    report.append("|---|---|------------|-------------|-------|------------|-----------|-------|--------|")

    cached = {}

    for N in N_GRID:
        tN = time.time()
        print(f"[Part 1] N={N} — computing {len(ALL_FAMILIES)} families x {REPS} reps ...", flush=True)

        ref_feats = feature_ensemble(ALL_FAMILIES["Lor4D"], N, REPS, SEED_BASE + 9000)
        mu = ref_feats.mean(axis=0)
        cov = np.cov(ref_feats.T) + 1e-8 * np.eye(3)
        cov_inv = np.linalg.inv(cov)

        std_scores = {}
        for name, gen_fn in ALL_FAMILIES.items():
            feats = feature_ensemble(gen_fn, N, REPS, SEED_BASE + 5000)
            std_scores[name] = mahalanobis_score(feats.mean(axis=0), mu, cov_inv)

        cached[N] = (mu, cov_inv, std_scores)
        lor4d_score = std_scores["Lor4D"]

        t_ref = time.time() - tN
        print(f"  Reference built in {t_ref:.0f}s", flush=True)

        for H in H_GRID:
            ds_feats = ds4d_feature_ensemble(N, H, REPS, SEED_BASE + 8000)
            ds_mu = ds_feats.mean(axis=0)
            ds_score = mahalanobis_score(ds_mu, mu, cov_inv)

            all_scores = {**std_scores, f"dS4D_H{H}": ds_score}
            sorted_all = sorted(all_scores.items(), key=lambda x: x[1])
            names = [n for n, _ in sorted_all]
            lor4d_rank = names.index("Lor4D") + 1
            ds_rank = names.index(f"dS4D_H{H}") + 1
            delta = ds_score - lor4d_score

            elapsed_N = time.time() - tN
            report.append(
                f"| {N} | {H} | {ds_score:.4f} | {lor4d_score:.4f} | "
                f"{delta:+.4f} | #{lor4d_rank} | #{ds_rank} | {len(all_scores)} | {elapsed_N:.0f} |"
            )

        print(f"  N={N} done in {time.time() - tN:.0f}s", flush=True)

    report.append("")

    # ================================================================
    # Part 2: H_c at very large N (reusing cache)
    # ================================================================
    report.append("=" * 72)
    report.append("# Part 2: Critical Hubble H_c")
    report.append("=" * 72)
    report.append("")
    report.append("| N | H_c (dS4D <= #2) | S_MD at H_c | Nearest non-dS competitor | Competitor S_MD |")
    report.append("|---|------------------|-------------|--------------------------|-----------------|")

    for N in N_GRID:
        mu, cov_inv, std_scores = cached[N]
        non_lor4d = {k: v for k, v in std_scores.items() if k != "Lor4D"}
        nearest_name = min(non_lor4d, key=non_lor4d.get)
        nearest_score = non_lor4d[nearest_name]

        h_c = None
        h_c_score = None
        for H in H_GRID:
            ds_feats = ds4d_feature_ensemble(N, H, REPS, SEED_BASE + 8000)
            ds_mu = ds_feats.mean(axis=0)
            ds_score = mahalanobis_score(ds_mu, mu, cov_inv)

            all_scores_h = {**std_scores, f"dS4D_H{H}": ds_score}
            sorted_h = sorted(all_scores_h.items(), key=lambda x: x[1])
            names_h = [n for n, _ in sorted_h]
            ds_rank = names_h.index(f"dS4D_H{H}") + 1
            if ds_rank <= 2:
                h_c = H
                h_c_score = ds_score

        if h_c is not None:
            report.append(f"| {N} | {h_c} | {h_c_score:.4f} | {nearest_name} | {nearest_score:.4f} |")
        else:
            report.append(f"| {N} | <{H_GRID[0]} | — | {nearest_name} | {nearest_score:.4f} |")

    report.append("")

    # ================================================================
    # Part 3: Scaling exponent beta(H) across full N range
    # ================================================================
    report.append("=" * 72)
    report.append("# Part 3: Scaling Exponent beta(H)")
    report.append("=" * 72)
    report.append("")
    report.append("| H | beta | SE(beta) | R^2 | N_points | interpretation |")
    report.append("|---|------|----------|-----|----------|----------------|")

    for H in H_GRID:
        if H == 0.0:
            report.append(f"| {H} | — | — | — | — | flat reference |")
            continue

        log_n_list = []
        log_s_list = []
        for N in N_GRID:
            mu, cov_inv, _ = cached[N]
            ds_feats = ds4d_feature_ensemble(N, H, REPS, SEED_BASE + 8000)
            ds_mu = ds_feats.mean(axis=0)
            ds_score = mahalanobis_score(ds_mu, mu, cov_inv)
            if ds_score > 0:
                log_n_list.append(np.log(N))
                log_s_list.append(np.log(ds_score))

        if len(log_n_list) >= 2:
            slope, intercept, r, p, se = linregress(log_n_list, log_s_list)
            if abs(slope) < 0.3:
                interp = "flat (indistinguishable)"
            elif abs(slope) < 0.8:
                interp = "sub-linear growth"
            elif abs(slope) < 1.5:
                interp = "~linear divergence"
            else:
                interp = "super-linear divergence"
            report.append(f"| {H} | {slope:.2f} | {se:.2f} | {r**2:.3f} | {len(log_n_list)} | {interp} |")
        else:
            report.append(f"| {H} | — | — | — | {len(log_n_list)} | insufficient data |")

    report.append("")

    # ================================================================
    # Part 4: Comparison with previous experiments
    # ================================================================
    report.append("=" * 72)
    report.append("# Part 4: Convergence Summary")
    report.append("=" * 72)
    report.append("")
    report.append("Previous beta values:")
    report.append("  Extended (N=16-256):  H=0.1: 0.91, H=0.3: 0.97")
    report.append("  Large-N (N=256-512):  H=0.1: 0.76, H=0.3: 2.53")
    report.append("")
    report.append("If beta here is similar to Large-N, the curvature convexity is confirmed.")
    report.append("If beta here is even larger, the acceleration effect strengthens.")
    report.append("")

    elapsed = time.time() - t0
    report.append(f"\nTotal elapsed: {elapsed:.0f}s")
    print(f"\nTotal elapsed: {elapsed:.0f}s", flush=True)

    out_path = OUT_DIR / "experiment_de_sitter_very_large_N.md"
    out_path.write_text("\n".join(report), encoding="utf-8")
    print(f"Report saved to {out_path}", flush=True)


if __name__ == "__main__":
    main()
