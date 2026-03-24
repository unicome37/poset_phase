"""Conjecture E — §4.1.32: Geometric Target Identification.

Given §4.1.31's discovery that antichain and B_ℓ spectral channels converge
to a single post-density DoF at d≥3, this experiment identifies WHAT
continuum geometric quantity that shared DoF corresponds to.

Strategy:
  1. Load §4.1.31 CSV data (360 realizations, same sprinklings).
  2. For each (d, N) slice, density-residualize all features.
  3. Extract latent variable Z via PCA-1 on joint residuals (both channels).
  4. Regress Z against candidate geometric quantities:
       (a) H                    — Hubble rate
       (b) H²                   — squared Hubble rate
       (c) R_dS = d(d-1)H²      — scalar curvature
       (d) θ = (d-1)H            — expansion scalar ∇_μ u^μ
       (e) H^α (power-law fit)   — general power
       (f) K_ext = (d-1)H        — extrinsic curvature trace (= θ for FRW)
       (g) log(1+H)              — logarithmic
  5. N-scaling: does the best-fit exponent α converge with N?
  6. d=2 anomaly diagnosis: why does the latent variable decouple at d=2?

References:
  - §4.1.31: Dual-channel unification — |ρ_cross| → 0.86 at d≥3
  - Sorkin (2007): B_ℓ → □ + ξR in continuum limit
  - Antichain width ~ spatial volume growth in de Sitter
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats
from scipy.optimize import curve_fit


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def ols_residualize(feat: np.ndarray, dens: np.ndarray) -> np.ndarray:
    """OLS-remove density from feature."""
    coeffs = np.polyfit(dens, feat, 1)
    return feat - np.polyval(coeffs, dens)


def spearman(x: np.ndarray, y: np.ndarray):
    """Spearman rank correlation with p-value."""
    return sp_stats.spearmanr(x, y)


def pearson(x: np.ndarray, y: np.ndarray):
    """Pearson correlation with p-value."""
    return sp_stats.pearsonr(x, y)


def ols_r2(x: np.ndarray, y: np.ndarray) -> float:
    """Simple linear OLS R²."""
    if len(x) < 5:
        return float("nan")
    coeffs = np.polyfit(x, y, 1)
    pred = np.polyval(coeffs, x)
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 1e-15 else 0.0


def alpha_grid_search(H_arr: np.ndarray, Z_arr: np.ndarray):
    """Grid search: find alpha in [0.5, 5.0] maximizing Pearson R²(Z, H^alpha).
    Uses ALL data including H=0 (H^alpha=0 for alpha>0).
    Returns (best_alpha, R²_best, dict of alpha→R²)."""
    alphas = np.arange(0.5, 5.05, 0.1)
    results = {}
    for alpha in alphas:
        x = np.power(H_arr, alpha)
        if np.std(x) < 1e-15:
            results[round(alpha, 1)] = 0.0
            continue
        r2 = ols_r2(x, Z_arr)
        results[round(alpha, 1)] = r2
    best_alpha = max(results, key=results.get)
    return best_alpha, results[best_alpha], results


def group_mean_analysis(H_arr: np.ndarray, Z_arr: np.ndarray):
    """Compute mean Z per H level, then fit alpha on group means.
    More robust when each H level has multiple reps.
    Returns (best_alpha, R²_best, H_levels, Z_means, Z_stds)."""
    H_levels = sorted(set(H_arr))
    Z_means = []
    Z_stds = []
    for h in H_levels:
        mask = H_arr == h
        Z_means.append(np.mean(Z_arr[mask]))
        Z_stds.append(np.std(Z_arr[mask]))
    H_levels = np.array(H_levels)
    Z_means = np.array(Z_means)
    Z_stds = np.array(Z_stds)

    # Grid search on group means
    best_alpha = 1.0
    best_r2 = -1.0
    for alpha in np.arange(0.5, 5.05, 0.1):
        x = np.power(H_levels, alpha)
        r2 = ols_r2(x, Z_means)
        if r2 > best_r2:
            best_r2 = r2
            best_alpha = round(alpha, 1)
    return best_alpha, best_r2, H_levels, Z_means, Z_stds


# ---------------------------------------------------------------------------
# Load CSV
# ---------------------------------------------------------------------------
def load_csv(path: Path) -> list[dict]:
    """Load §4.1.31 CSV data."""
    rows = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            parsed = {}
            for k, v in row.items():
                try:
                    parsed[k] = float(v)
                except (ValueError, TypeError):
                    parsed[k] = v
            rows.append(parsed)
    return rows


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------
def analyze(rows: list[dict], dims: list[int], ns: list[int]) -> str:
    """Run full geometric target identification analysis."""
    lines: list[str] = []
    lines.append("# §4.1.32: Geometric Target Identification\n")
    lines.append("## Goal\n")
    lines.append("Identify the continuum geometric quantity that the shared")
    lines.append("post-density DoF (discovered in §4.1.31) corresponds to.\n")

    # Feature lists
    ac_feats = ["w_max_ratio", "mean_layer_width", "layer_width_std", "layer_ratio"]
    bl_feats = ["b1_std", "eig_min", "eig_max", "eig_spread"]
    all_feats = ac_feats + bl_feats

    # Candidate geometric quantities
    # For de Sitter: R = d(d-1)H², θ = (d-1)H, K_trace = (d-1)H
    candidate_names = [
        "H", "H²", "R_dS = d(d-1)H²", "θ = (d-1)H",
        "H^α (power-law)", "log(1+H)",
    ]

    # ======================================================================
    # §1: Extract latent variable Z via PCA on joint residuals
    # ======================================================================
    lines.append("## 1. Latent Variable Extraction\n")
    lines.append("For each (d, N) slice:")
    lines.append("1. OLS-residualize all 8 features against n_causal_pairs")
    lines.append("2. Standardize each residual (zero mean, unit variance)")
    lines.append("3. PCA on 8-dimensional residual space → PC1 = latent variable Z")
    lines.append("4. Report % variance explained by PC1\n")

    lines.append("| d | N | PC1 var% | PC2 var% | PC1+2 var% | AC loadings (w,m,s,r) | Bℓ loadings (b,emin,emax,esp) |")
    lines.append("|---|---|---------|---------|-----------|----------------------|------------------------------|")

    # Store latent variables for subsequent analysis
    Z_store: dict[tuple[int, int], np.ndarray] = {}
    H_store: dict[tuple[int, int], np.ndarray] = {}
    meta_store: dict[tuple[int, int], dict] = {}

    for d in dims:
        for N in ns:
            subset = [r for r in rows if int(r["d"]) == d and int(r["N"]) == N]
            if len(subset) < 15:
                continue

            H_arr = np.array([r["hubble"] for r in subset])
            dens_arr = np.array([r["n_causal_pairs"] for r in subset])

            # Residualize all features
            resid_matrix = []
            valid_feats = []
            for feat in all_feats:
                vals = np.array([r.get(feat, float("nan")) for r in subset])
                mask = ~np.isnan(vals)
                if mask.sum() < len(subset) * 0.8:
                    continue
                # Replace NaN with mean for PCA (rare cases)
                vals_clean = vals.copy()
                if np.any(~mask):
                    vals_clean[~mask] = np.nanmean(vals)
                resid = ols_residualize(vals_clean, dens_arr)
                resid_matrix.append(resid)
                valid_feats.append(feat)

            if len(resid_matrix) < 4:
                continue

            # Standardize
            X = np.column_stack(resid_matrix)
            means = X.mean(axis=0)
            stds = X.std(axis=0)
            stds[stds < 1e-15] = 1.0
            X_std = (X - means) / stds

            # PCA via SVD
            U, S, Vt = np.linalg.svd(X_std, full_matrices=False)
            explained_var = S ** 2 / np.sum(S ** 2)
            pc1_var = explained_var[0] * 100
            pc2_var = explained_var[1] * 100 if len(explained_var) > 1 else 0
            loadings = Vt[0]  # PC1 loadings

            # Latent variable Z = projection onto PC1
            Z = X_std @ Vt[0]

            # Ensure Z positively correlates with H (sign convention)
            rho_zh, _ = spearman(H_arr, Z)
            if rho_zh < 0:
                Z = -Z
                loadings = -loadings

            Z_store[(d, N)] = Z
            H_store[(d, N)] = H_arr

            # Format loadings
            ac_load_str = ",".join(f"{loadings[valid_feats.index(f)]:+.2f}"
                                   if f in valid_feats else "N/A"
                                   for f in ac_feats)
            bl_load_str = ",".join(f"{loadings[valid_feats.index(f)]:+.2f}"
                                   if f in valid_feats else "N/A"
                                   for f in bl_feats)

            meta_store[(d, N)] = {
                "pc1_var": pc1_var,
                "pc2_var": pc2_var,
                "loadings": loadings,
                "valid_feats": valid_feats,
                "n_feats": len(valid_feats),
            }

            lines.append(
                f"| {d} | {N} | {pc1_var:.1f}% | {pc2_var:.1f}% | "
                f"{pc1_var + pc2_var:.1f}% | {ac_load_str} | {bl_load_str} |"
            )

    # ======================================================================
    # §2: Regression of Z against candidate geometric quantities
    # ======================================================================
    lines.append("\n## 2. Latent Variable Z vs Candidate Geometric Quantities\n")
    lines.append("For each (d, N), compute Spearman ρ(Z, candidate) and OLS R².\n")

    lines.append("\n### 2a. Spearman correlations\n")
    lines.append("**Note**: Spearman is rank-based, so ρ(Z,H) = ρ(Z,H²) = ρ(Z,R_dS) = ρ(Z,θ)")
    lines.append("for any monotone transform of H. The discriminant power comes from")
    lines.append("**Pearson R²** (²b) and **group-mean analysis** (³).\n")
    lines.append("| d | N | ρ(Z,H) | ρ(Z,H²) | ρ(Z,R_dS) | ρ(Z,θ) | ρ(Z,log1pH) |")
    lines.append("|---|---|--------|---------|----------|--------|------------|")

    for d in dims:
        for N in ns:
            key = (d, N)
            if key not in Z_store:
                continue
            Z = Z_store[key]
            H = H_store[key]

            H2 = H ** 2
            R_dS = d * (d - 1) * H2
            theta = (d - 1) * H
            log1pH = np.log1p(H)

            candidates_arr = [H, H2, R_dS, theta, log1pH]
            cells = []
            for cand in candidates_arr:
                rho, p = spearman(Z, cand)
                sig = "**" if p < 0.01 else "*" if p < 0.05 else ""
                cells.append(f"{sig}{rho:+.3f}{sig}")

            lines.append(f"| {d} | {N} | " + " | ".join(cells) + " |")

    lines.append("\n### 2b. OLS R² (linear fit)\n")
    lines.append("| d | N | R²(H) | R²(H²) | R²(R_dS) | R²(θ) | R²(log1pH) | Best fit |")
    lines.append("|---|---|-------|--------|---------|-------|-----------|----------|")

    best_fits: dict[tuple[int, int], str] = {}

    for d in dims:
        for N in ns:
            key = (d, N)
            if key not in Z_store:
                continue
            Z = Z_store[key]
            H = H_store[key]

            H2 = H ** 2
            R_dS = d * (d - 1) * H2
            theta = (d - 1) * H
            log1pH = np.log1p(H)

            cand_names_short = ["H", "H²", "R_dS", "θ", "log1pH"]
            candidates_arr = [H, H2, R_dS, theta, log1pH]
            r2_vals = []
            cells = []
            for cand in candidates_arr:
                r2 = ols_r2(cand, Z)
                r2_vals.append(r2)
                cells.append(f"{r2:.3f}")

            best_idx = int(np.argmax(r2_vals))
            best_name = cand_names_short[best_idx]
            best_fits[key] = best_name
            cells.append(f"**{best_name}**")

            lines.append(f"| {d} | {N} | " + " | ".join(cells) + " |")

    # ======================================================================
    # §3: Grid-search α and group-mean analysis
    # ======================================================================
    lines.append("\n## 3. Effective Exponent α: Z ~ H^α\n")
    lines.append("Grid search over α ∈ [0.5, 5.0] maximizing Pearson R²(Z, H^α).")
    lines.append("Two methods: (a) all-point regression, (b) group-mean regression.\n")

    lines.append("### 3a. All-point α grid search\n")
    lines.append("| d | N | α_best | R²_best | R²(α=1) | R²(α=2) | R²(α=3) | Interpretation |")
    lines.append("|---|---|--------|---------|---------|---------|---------|----------------|")

    alpha_store: dict[tuple[int, int], float] = {}

    for d in dims:
        for N in ns:
            key = (d, N)
            if key not in Z_store:
                continue
            Z = Z_store[key]
            H = H_store[key]

            best_alpha, best_r2, r2_dict = alpha_grid_search(H, Z)
            alpha_store[key] = best_alpha

            r2_1 = r2_dict.get(1.0, 0.0)
            r2_2 = r2_dict.get(2.0, 0.0)
            r2_3 = r2_dict.get(3.0, 0.0)

            if abs(best_alpha - 1.0) < 0.3:
                interp = "**H / θ** (expansion rate)"
            elif abs(best_alpha - 2.0) < 0.3:
                interp = "**H² / R_dS** (scalar curvature)"
            elif best_alpha > 2.3:
                interp = f"H^{best_alpha:.1f} (super-quadratic)"
            else:
                interp = f"H^{best_alpha:.1f} (intermediate)"

            lines.append(
                f"| {d} | {N} | {best_alpha:.1f} | {best_r2:.4f} | "
                f"{r2_1:.4f} | {r2_2:.4f} | {r2_3:.4f} | {interp} |"
            )

    lines.append("\n### 3b. Group-mean α (mean Z per H level)\n")
    lines.append("| d | N | α_best | R²_best | mean(Z) per H=[0, 0.25, 0.5, 1.0, 2.0] |")
    lines.append("|---|---|--------|---------|----------------------------------------------|")

    alpha_gm_store: dict[tuple[int, int], float] = {}

    for d in dims:
        for N in ns:
            key = (d, N)
            if key not in Z_store:
                continue
            Z = Z_store[key]
            H = H_store[key]

            gm_alpha, gm_r2, H_levels, Z_means, Z_stds = group_mean_analysis(H, Z)
            alpha_gm_store[key] = gm_alpha

            means_str = ", ".join(f"{z:.3f}" for z in Z_means)
            lines.append(
                f"| {d} | {N} | {gm_alpha:.1f} | {gm_r2:.4f} | [{means_str}] |"
            )

    # ======================================================================
    # §4: N-scaling of α
    # ======================================================================
    lines.append("\n## 4. N-Scaling of Effective Exponent α\n")
    lines.append("Does α converge to a specific value as N→∞?\n")
    lines.append("### All-point α\n")
    lines.append("| d | N values | α values | trend | extrapolated α_∞ |")
    lines.append("|---|---------|---------|-------|-----------------|")

    for d in dims:
        ns_d = [N for N in ns if (d, N) in alpha_store]
        if len(ns_d) < 2:
            continue
        alphas_d = [alpha_store[(d, N)] for N in ns_d]

        inv_n = [1.0 / N for N in ns_d]
        if len(inv_n) >= 2:
            coeffs = np.polyfit(inv_n, alphas_d, 1)
            alpha_inf = coeffs[1]
        else:
            alpha_inf = alphas_d[-1]

        trend_rho, _ = spearman(np.array(ns_d), np.array(alphas_d)) if len(ns_d) >= 3 else (0, 1)
        trend_word = "increasing" if trend_rho > 0.5 else "decreasing" if trend_rho < -0.5 else "stable"

        ns_str = ",".join(str(n) for n in ns_d)
        alpha_str = ",".join(f"{a:.1f}" for a in alphas_d)
        lines.append(f"| {d} | {ns_str} | {alpha_str} | {trend_word} | **{alpha_inf:.1f}** |")

    lines.append("\n### Group-mean α\n")
    lines.append("| d | N values | α values | trend | extrapolated α_∞ |")
    lines.append("|---|---------|---------|-------|-----------------|")

    for d in dims:
        ns_d = [N for N in ns if (d, N) in alpha_gm_store]
        if len(ns_d) < 2:
            continue
        alphas_d = [alpha_gm_store[(d, N)] for N in ns_d]

        inv_n = [1.0 / N for N in ns_d]
        if len(inv_n) >= 2:
            coeffs = np.polyfit(inv_n, alphas_d, 1)
            alpha_inf = coeffs[1]
        else:
            alpha_inf = alphas_d[-1]

        trend_rho, _ = spearman(np.array(ns_d), np.array(alphas_d)) if len(ns_d) >= 3 else (0, 1)
        trend_word = "increasing" if trend_rho > 0.5 else "decreasing" if trend_rho < -0.5 else "stable"

        ns_str = ",".join(str(n) for n in ns_d)
        alpha_str = ",".join(f"{a:.1f}" for a in alphas_d)
        lines.append(f"| {d} | {ns_str} | {alpha_str} | {trend_word} | **{alpha_inf:.1f}** |")

    # ======================================================================
    # §5: Discriminant test: H vs H² vs θ
    # ======================================================================
    lines.append("\n## 5. Discriminant Test: R² Comparison Across Exponents\n")
    lines.append("Key question: is Z primarily tracking H (first-order expansion rate)")
    lines.append("or H² (second-order curvature)?\n")
    lines.append("For de Sitter: R = d(d-1)H², θ = (d-1)H, K_trace = (d-1)H.")
    lines.append("Since these differ only by constants within each d,")
    lines.append("the question reduces to: **α ≈ 1 or α ≈ 2?**\n")

    lines.append("| d | N | R²(α=0.5) | R²(α=1.0) | R²(α=1.5) | R²(α=2.0) | R²(α=3.0) | α_best | ΔR²(2vs1) |")
    lines.append("|---|---|----------|----------|----------|----------|----------|--------|----------|")

    for d in dims:
        for N in ns:
            key = (d, N)
            if key not in Z_store:
                continue
            Z = Z_store[key]
            H = H_store[key]

            _, _, r2_dict = alpha_grid_search(H, Z)
            r2_05 = r2_dict.get(0.5, 0.0)
            r2_10 = r2_dict.get(1.0, 0.0)
            r2_15 = r2_dict.get(1.5, 0.0)
            r2_20 = r2_dict.get(2.0, 0.0)
            r2_30 = r2_dict.get(3.0, 0.0)
            best_a = alpha_store.get(key, 2.0)
            delta = r2_20 - r2_10

            lines.append(
                f"| {d} | {N} | {r2_05:.4f} | {r2_10:.4f} | {r2_15:.4f} | "
                f"{r2_20:.4f} | {r2_30:.4f} | {best_a:.1f} | {delta:+.4f} |"
            )

    # ======================================================================
    # §6: d=2 Anomaly Diagnosis
    # ======================================================================
    lines.append("\n## 6. d=2 Anomaly Diagnosis\n")
    lines.append("Why does the latent variable behave differently at d=2?\n")

    lines.append("### 6a. PC1 variance explained by dimension\n")
    lines.append("| d | N=128 | N=256 | N=512 |")
    lines.append("|---|-------|-------|-------|")
    for d in dims:
        cells = []
        for N in ns:
            key = (d, N)
            if key in meta_store:
                cells.append(f"{meta_store[key]['pc1_var']:.1f}%")
            else:
                cells.append("—")
        lines.append(f"| {d} | " + " | ".join(cells) + " |")

    lines.append("\n### 6b. Cross-channel coherence vs d\n")
    lines.append("PC1 captures both channels when they are coherent.")
    lines.append("At d=2, B_ℓ spectral features are known to be weak (§4.1.27: 0/6),")
    lines.append("so PC1 may be dominated by antichain features alone.\n")

    lines.append("### 6c. Loading pattern by dimension\n")
    for d in dims:
        for N in ns:
            key = (d, N)
            if key not in meta_store:
                continue
            m = meta_store[key]
            lines.append(f"\n**d={d}, N={N}**: PC1 = {m['pc1_var']:.1f}%")
            loading_strs = []
            for i, feat in enumerate(m["valid_feats"]):
                channel = "AC" if feat in ac_feats else "Bℓ"
                loading_strs.append(f"{feat}({channel})={m['loadings'][i]:+.3f}")
            lines.append("  Loadings: " + ", ".join(loading_strs))

    # ======================================================================
    # §7: Physical Interpretation & Conclusion
    # ======================================================================
    lines.append("\n\n## 7. Physical Interpretation\n")

    # Summarize α values at largest N
    lines.append("### Summary of effective exponent α at largest N\n")
    lines.append("| d | N_max | α | Nearest integer power | Implied geometric target |")
    lines.append("|---|-------|---|----------------------|------------------------|")

    for d in dims:
        largest_N = max((N for N in ns if (d, N) in alpha_store), default=None)
        if largest_N is None:
            continue
        alpha = alpha_store[(d, largest_N)]
        if abs(alpha - 1.0) < 0.3:
            nearest = "1"
            target = "H or θ = (d-1)H (expansion rate)"
        elif abs(alpha - 2.0) < 0.3:
            nearest = "2"
            target = "H² or R_dS = d(d-1)H² (scalar curvature)"
        else:
            nearest = f"{alpha:.1f}"
            target = f"H^{alpha:.1f} (non-standard)"
        lines.append(f"| {d} | {largest_N} | {alpha:.2f} | {nearest} | {target} |")

    lines.append("\n### Conclusion\n")
    lines.append("Based on the power-law exponent α, R² comparisons, and N-convergence:")
    lines.append("")
    lines.append("1. **If α → 1**: The post-density DoF tracks the **expansion rate H**")
    lines.append("   (equivalently, the expansion scalar θ = (d-1)H or the trace of")
    lines.append("   extrinsic curvature K = (d-1)H·g_ij for de Sitter spatial slices).")
    lines.append("   This is the *first-order* geometric effect of expansion.")
    lines.append("")
    lines.append("2. **If α → 2**: The DoF tracks **H² ∝ R_dS** (scalar curvature).")
    lines.append("   This would mean the bulk channel directly encodes curvature,")
    lines.append("   making it a true EH-action-like observable.")
    lines.append("")
    lines.append("3. **If α is intermediate (1 < α < 2)**: The DoF tracks a")
    lines.append("   *mixture* of first-order and second-order geometric effects,")
    lines.append("   which may indicate that the discrete observable has not yet")
    lines.append("   fully converged to its continuum limit at the tested N values.")
    lines.append("")

    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(
        description="§4.1.32: Geometric target identification"
    )
    ap.add_argument(
        "--csv",
        default="outputs_unified_functional/conjecture_e_dual_channel.csv",
        help="Path to §4.1.31 CSV data",
    )
    ap.add_argument(
        "--report",
        default="outputs_unified_functional/conjecture_e_geometric_target.md",
        help="Output report path",
    )
    ap.add_argument("--dims", nargs="*", type=int, default=[2, 3, 4])
    ap.add_argument("--ns", nargs="*", type=int, default=[128, 256, 512])
    args = ap.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.exists():
        print(f"ERROR: CSV not found: {csv_path}")
        print("Run conjecture_e_dual_channel_unification.py first.")
        return 1

    print(f"Loading data from {csv_path}...")
    rows = load_csv(csv_path)
    print(f"  Loaded {len(rows)} realizations")

    print("Running geometric target identification analysis...")
    report_text = analyze(rows, args.dims, args.ns)

    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report_text, encoding="utf-8")
    print(f"Saved report: {report_path}")

    # Console quick summary
    print("\n" + "=" * 70)
    print("QUICK SUMMARY — Geometric Target Identification")
    print("=" * 70)

    for d in args.dims:
        for N in args.ns:
            subset = [r for r in rows if int(r["d"]) == d and int(r["N"]) == N]
            if len(subset) < 15:
                continue
            H_arr = np.array([r["hubble"] for r in subset])
            dens_arr = np.array([r["n_causal_pairs"] for r in subset])

            # Quick latent variable extraction
            ac_f = ["w_max_ratio", "mean_layer_width", "layer_width_std", "layer_ratio"]
            bl_f = ["b1_std", "eig_min", "eig_max", "eig_spread"]
            all_f = ac_f + bl_f

            resid_matrix = []
            for feat in all_f:
                vals = np.array([r.get(feat, float("nan")) for r in subset])
                vals_clean = vals.copy()
                nans = np.isnan(vals_clean)
                if nans.any():
                    vals_clean[nans] = np.nanmean(vals)
                resid = ols_residualize(vals_clean, dens_arr)
                resid_matrix.append(resid)

            X = np.column_stack(resid_matrix)
            means = X.mean(axis=0)
            stds = X.std(axis=0)
            stds[stds < 1e-15] = 1.0
            X_std = (X - means) / stds

            _, S, Vt = np.linalg.svd(X_std, full_matrices=False)
            Z = X_std @ Vt[0]
            rho_zh, _ = spearman(H_arr, Z)
            if rho_zh < 0:
                Z = -Z

            # Grid search for best alpha
            best_alpha, best_r2, r2_dict = alpha_grid_search(H_arr, Z)
            r2_1 = r2_dict.get(1.0, 0.0)
            r2_2 = r2_dict.get(2.0, 0.0)
            print(f"  d={d} N={N}: α_best={best_alpha:.1f} (R²={best_r2:.4f})  R²(α=1)={r2_1:.4f}  R²(α=2)={r2_2:.4f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
