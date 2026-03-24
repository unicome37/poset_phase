"""Conjecture E — PC-Space Curvature Proxy Design (§4.1.26).

Following §4.1.25's discovery that curvature information lives in the PC2+
shape subspace (11.4% of variance), this experiment:

  1. Per-dimension PCA: fit PCA within each d separately (since loadings
     differ across dimensions)
  2. Density subtraction: project each realization onto PC2+ (remove PC1)
  3. Optimal proxy: find the linear combination in PC2+ that maximizes
     Spearman(proxy, H²) — this is the "PC-optimal curvature proxy"
  4. Back-project: convert the optimal PC2+ direction back to raw p_k
     coefficients → a closed-form curvature observable
  5. Compare: PC-optimal vs bdg_d2c vs KL vs raw PC2 score
  6. N-convergence: does the optimal proxy improve with larger N?

Data: loads §4.1.24 CSV (360 realizations with C0–C4).
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats
from scipy.optimize import minimize_scalar


def load_data(csv_path: str) -> list[dict]:
    rows = []
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            d = {}
            for k, v in row.items():
                try:
                    d[k] = float(v)
                except ValueError:
                    d[k] = v
            rows.append(d)
    return rows


def build_distribution_matrix(rows: list[dict], max_k: int = 4):
    n = len(rows)
    X = np.zeros((n, max_k + 1))
    for i, r in enumerate(rows):
        for k in range(max_k + 1):
            X[i, k] = float(r.get(f"C{k}", 0))
    row_sums = X.sum(axis=1, keepdims=True)
    row_sums = np.clip(row_sums, 1.0, None)
    P = X / row_sums
    return P, X


def pca_manual(X_centered):
    """PCA via SVD on already-centered data. Returns scores, components, explained_ratio."""
    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
    scores = U * S
    components = Vt
    ev = (S ** 2) / (X_centered.shape[0] - 1)
    total = ev.sum()
    er = ev / total if total > 0 else ev
    return scores, components, er


def find_optimal_angle_2d(scores_2d, h2):
    """Find angle θ in 2D PC subspace that maximizes |Spearman(cos θ · PC_a + sin θ · PC_b, H²)|."""
    def neg_abs_spearman(theta):
        proxy = np.cos(theta) * scores_2d[:, 0] + np.sin(theta) * scores_2d[:, 1]
        rho, _ = sp_stats.spearmanr(h2, proxy)
        return -abs(rho)

    # Grid search then refine
    best_theta = 0.0
    best_val = 0.0
    for t in np.linspace(0, np.pi, 180):
        val = -neg_abs_spearman(t)
        if val > best_val:
            best_val = val
            best_theta = t

    # Refine with Brent
    result = minimize_scalar(neg_abs_spearman,
                             bounds=(best_theta - 0.05, best_theta + 0.05),
                             method='bounded')
    return result.x, -result.fun


def find_optimal_direction_nd(scores_nd, h2, n_restarts=50):
    """Find unit vector in n-D PC subspace maximizing |Spearman(proxy, H²)|.
    Uses random restarts with coordinate descent."""
    n_dims = scores_nd.shape[1]
    if n_dims == 0:
        return np.array([]), 0.0

    if n_dims == 1:
        rho, _ = sp_stats.spearmanr(h2, scores_nd[:, 0])
        return np.array([1.0 if rho > 0 else -1.0]), abs(rho)

    best_rho = 0.0
    best_w = np.zeros(n_dims)

    rng = np.random.RandomState(42)
    for _ in range(n_restarts):
        w = rng.randn(n_dims)
        w /= np.linalg.norm(w)

        for _iter in range(20):
            proxy = scores_nd @ w
            rho, _ = sp_stats.spearmanr(h2, proxy)
            if abs(rho) > best_rho:
                best_rho = abs(rho)
                best_w = w.copy() * np.sign(rho)

            # Gradient-free: try perturbing each coordinate
            improved = False
            for dim in range(n_dims):
                for delta in [0.1, -0.1, 0.05, -0.05]:
                    w2 = w.copy()
                    w2[dim] += delta
                    w2 /= np.linalg.norm(w2)
                    proxy2 = scores_nd @ w2
                    rho2, _ = sp_stats.spearmanr(h2, proxy2)
                    if abs(rho2) > abs(rho):
                        w = w2
                        rho = rho2
                        improved = True
                        if abs(rho2) > best_rho:
                            best_rho = abs(rho2)
                            best_w = w2.copy() * np.sign(rho2)
            if not improved:
                break

    return best_w, best_rho


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="outputs_unified_functional/conjecture_e_shape_ratio_family.csv")
    ap.add_argument("--report", default="outputs_unified_functional/conjecture_e_pc_space_proxy.md")
    ap.add_argument("--max-k", type=int, default=4)
    args = ap.parse_args()

    rows = load_data(args.csv)
    print(f"Loaded {len(rows)} rows")

    P, X_raw = build_distribution_matrix(rows, args.max_k)
    n_samples, n_features = P.shape

    dims = np.array([r["d"] for r in rows])
    ns = np.array([r["N"] for r in rows])
    hubbles = np.array([r["hubble"] for r in rows])
    h2 = hubbles ** 2
    kl_vals = np.array([r.get("S10_kl_vs_flat", float("nan")) for r in rows])
    bdg_d2c_vals = np.array([r.get("S8_bdg_d2c_shape", 0.0) for r in rows])

    unique_dims = sorted(set(dims))
    unique_ns = sorted(set(ns))

    lines: list[str] = []
    lines.append("# Conjecture E — PC-Space Curvature Proxy (§4.1.26)\n")
    lines.append(f"Loaded {n_samples} realizations, p_0 through p_{args.max_k}\n")

    # ══════════════════════════════════════════════════════════════
    # PART 1: Per-dimension PCA + optimal proxy in PC2+ subspace
    # ══════════════════════════════════════════════════════════════
    lines.append("\n## Part 1: Per-dimension optimal PC2+ proxy\n")
    lines.append("Strategy: within each d, fit PCA on standardized p_k, then find the unit "
                 "vector in PC2+ subspace that maximizes |Spearman(proxy, H²)|.\n")

    all_results = []

    for d_val in unique_dims:
        d_mask = dims == d_val
        P_d = P[d_mask]
        h2_d = h2[d_mask]
        ns_d = ns[d_mask]
        kl_d = kl_vals[d_mask]
        bdg_d = bdg_d2c_vals[d_mask]

        # Standardize
        p_mean = P_d.mean(axis=0)
        p_std = P_d.std(axis=0, ddof=0)
        p_std = np.where(p_std < 1e-12, 1.0, p_std)
        P_std = (P_d - p_mean) / p_std

        # PCA
        scores, components, er = pca_manual(P_std)

        lines.append(f"\n### d = {int(d_val)}\n")
        lines.append(f"PC1 variance: {er[0]*100:.1f}%\n")

        # PC2+ scores (remove PC1)
        scores_shape = scores[:, 1:]  # all PCs except PC1
        components_shape = components[1:]

        # Find optimal direction in full PC2+ subspace
        # Do this pooled across all N first
        opt_w_pooled, opt_rho_pooled = find_optimal_direction_nd(scores_shape, h2_d)

        lines.append(f"**Pooled optimal direction** (across all N): |ρ| = {opt_rho_pooled:.3f}\n")
        lines.append(f"Weights on PC2..PC{n_features}: {', '.join(f'{w:+.3f}' for w in opt_w_pooled)}\n")

        # Convert optimal direction back to raw p_k coefficients
        # opt_w is in PC2+ space → multiply by components_shape to get standardized p_k
        if len(opt_w_pooled) > 0:
            coef_std = opt_w_pooled @ components_shape  # in standardized p_k space
            coef_raw = coef_std / p_std  # in raw p_k space
            coef_raw_norm = coef_raw / np.max(np.abs(coef_raw))  # normalize for readability

            lines.append(f"Back-projected to raw p_k coefficients (normalized):")
            lines.append(f"  [{', '.join(f'{c:+.3f}' for c in coef_raw_norm)}]\n")

            # Also express in terms of raw C_k (unnormalized)
            # p_k = C_k / total → proxy ∝ Σ coef_raw_k * C_k / total
            lines.append(f"Closed-form proxy ∝ Σ_k w_k · p_k where w = [{', '.join(f'{c:+.3f}' for c in coef_raw_norm)}]\n")

        # Per-N analysis
        lines.append("#### Per-N slice results\n")
        lines.append("| N | PC2-opt ρ | p | bdg_d2c ρ | p | KL ρ | p | PC2 alone ρ | p |")
        lines.append("|---|----------|---|----------|---|------|---|------------|---|")

        for N_val in unique_ns:
            n_mask = ns_d == N_val
            if n_mask.sum() < 5:
                continue

            h2_sub = h2_d[n_mask]
            scores_sub = scores_shape[n_mask]

            # Optimal proxy score (using pooled weights)
            if len(opt_w_pooled) > 0:
                proxy_sub = scores_sub @ opt_w_pooled
                rho_opt, p_opt = sp_stats.spearmanr(h2_sub, proxy_sub)
            else:
                rho_opt, p_opt = 0.0, 1.0

            # bdg_d2c
            bdg_sub = bdg_d[n_mask]
            rho_bdg, p_bdg = sp_stats.spearmanr(h2_sub, bdg_sub)

            # KL
            kl_sub = kl_d[n_mask]
            valid = ~np.isnan(kl_sub)
            if valid.sum() >= 5:
                rho_kl, p_kl = sp_stats.spearmanr(h2_sub[valid], kl_sub[valid])
            else:
                rho_kl, p_kl = float("nan"), 1.0

            # PC2 alone
            if scores_sub.shape[1] >= 1:
                rho_pc2, p_pc2 = sp_stats.spearmanr(h2_sub, scores_sub[:, 0])
            else:
                rho_pc2, p_pc2 = 0.0, 1.0

            lines.append(f"| {int(N_val)} | {rho_opt:+.3f} | {p_opt:.1e} | "
                         f"{rho_bdg:+.3f} | {p_bdg:.1e} | "
                         f"{rho_kl:+.3f} | {p_kl:.1e} | "
                         f"{rho_pc2:+.3f} | {p_pc2:.1e} |")

            # Also find per-N optimal (may differ from pooled)
            opt_w_n, opt_rho_n = find_optimal_direction_nd(scores_sub, h2_sub, n_restarts=30)

            all_results.append({
                "d": int(d_val), "N": int(N_val),
                "rho_opt_pooled": rho_opt, "p_opt_pooled": p_opt,
                "rho_opt_local": opt_rho_n,
                "rho_bdg": rho_bdg, "p_bdg": p_bdg,
                "rho_kl": rho_kl, "p_kl": p_kl,
                "rho_pc2": rho_pc2, "p_pc2": p_pc2,
            })

        # Per-N optimal (local fit)
        lines.append("\n#### Per-N locally-optimized proxy\n")
        lines.append("| N | local-opt ρ | pooled-opt ρ | bdg_d2c ρ | KL ρ |")
        lines.append("|---|-----------|-------------|----------|------|")
        for r in all_results:
            if r["d"] == int(d_val):
                lines.append(f"| {r['N']} | {r['rho_opt_local']:+.3f} | "
                             f"{r['rho_opt_pooled']:+.3f} | {r['rho_bdg']:+.3f} | "
                             f"{r['rho_kl']:+.3f} |")

    # ══════════════════════════════════════════════════════════════
    # PART 2: Summary comparison
    # ══════════════════════════════════════════════════════════════
    lines.append("\n## Part 2: Summary — PC-optimal vs bdg_d2c vs KL\n")

    # Count significant positive correlations
    lines.append("### Significant positive slices (ρ > 0, p < 0.05)\n")
    lines.append("| Method | Count | Slices |")
    lines.append("|--------|-------|--------|")

    methods = {
        "PC2-opt (pooled)": ("rho_opt_pooled", "p_opt_pooled"),
        "PC2-opt (local)": ("rho_opt_local", None),
        "bdg_d2c": ("rho_bdg", "p_bdg"),
        "KL(p||p_flat)": ("rho_kl", "p_kl"),
        "PC2 alone": ("rho_pc2", "p_pc2"),
    }

    for name, (rho_key, p_key) in methods.items():
        sig_slices = []
        for r in all_results:
            rho = r[rho_key]
            if p_key and not np.isnan(r.get(p_key, float("nan"))):
                p = r[p_key]
                if rho > 0 and p < 0.05:
                    sig_slices.append(f"d={r['d']}N={r['N']}")
            elif p_key is None:
                # For local-opt we don't have p-values directly, count if |rho| > 0.3
                if rho > 0.3:
                    sig_slices.append(f"d={r['d']}N={r['N']}")
        lines.append(f"| {name} | {len(sig_slices)}/9 | {', '.join(sig_slices)} |")

    # Mean rho comparison
    lines.append("\n### Mean |ρ| across all 9 slices\n")
    lines.append("| Method | Mean |ρ| | Min ρ | Max ρ |")
    lines.append("|--------|---------|-------|-------|")

    for name, (rho_key, _) in methods.items():
        rhos = [r[rho_key] for r in all_results if not np.isnan(r[rho_key])]
        if rhos:
            lines.append(f"| {name} | {np.mean(np.abs(rhos)):.3f} | "
                         f"{min(rhos):+.3f} | {max(rhos):+.3f} |")

    # ══════════════════════════════════════════════════════════════
    # PART 3: N-convergence
    # ══════════════════════════════════════════════════════════════
    lines.append("\n## Part 3: N-convergence of PC-optimal proxy\n")
    lines.append("Does the shape-space curvature signal improve with N?\n")
    lines.append("| d | N=128 ρ | N=256 ρ | N=512 ρ | Trend |")
    lines.append("|---|---------|---------|---------|-------|")

    for d_val in unique_dims:
        rhos_by_n = {}
        for r in all_results:
            if r["d"] == int(d_val):
                rhos_by_n[r["N"]] = r["rho_opt_pooled"]
        r128 = rhos_by_n.get(128, float("nan"))
        r256 = rhos_by_n.get(256, float("nan"))
        r512 = rhos_by_n.get(512, float("nan"))
        if r512 > r256 > r128:
            trend = "↑ improving"
        elif r512 < r128:
            trend = "↓ degrading"
        else:
            trend = "~ mixed"
        lines.append(f"| {int(d_val)} | {r128:+.3f} | {r256:+.3f} | {r512:+.3f} | {trend} |")

    # ══════════════════════════════════════════════════════════════
    # PART 4: Physical interpretation
    # ══════════════════════════════════════════════════════════════
    lines.append("\n## Part 4: Interpretation\n")
    lines.append("""### Two-line architecture for pure-causal curvature encoding

**Line A — Wall / Admissibility (PC1 = density mode)**
- Encodes causal connectivity density
- Anti-monotone with curvature (raw C_k) or positive (normalized p_k)
- Already implemented as sigmoid wall in F7/F8a
- 88.6% of total variance → dominant, robust, cross-dimensional

**Line B — Bulk Curvature Proxy (PC2+ = shape residual)**
- Encodes shape distortion of C_k distribution
- Positive correlation with H² (when density mode is removed)
- Only 11.4% of variance → subleading, weaker signal
- bdg_d2c is an impure projection (71.7% shape, 28.3% density)
- PC-optimal proxy is a pure shape projection

### Key insight
> Wall succeeded first because it rides the dominant mode.
> Bulk proxy is harder because it must extract a subleading residual.
> This is not an engineering accident — it is dictated by the information hierarchy
> of causal interval distributions under de Sitter expansion.
""")

    # Save report
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_text = "\n".join(lines) + "\n"
    report_path.write_text(report_text, encoding="utf-8")
    print(f"\nSaved: {report_path}")

    # Print summary
    print("\n=== SUMMARY ===")
    for line in lines[-50:]:
        print(line)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
