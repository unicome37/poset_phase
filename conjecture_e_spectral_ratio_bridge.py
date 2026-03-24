"""Conjecture E — T4: Spectral Ratio Bridge (EH Bridge Reanalysis).

Tests whether eigenvalue RATIOS or GAPS of the B_ℓ d'Alembertian matrix
track H² (scalar curvature R) rather than H (expansion rate).

Hypothesis: Individual eigenvalues scale as density × f(H), so after
density removal they track H. But eigenvalue RATIOS cancel the density
scale factor, potentially exposing the H² dependence directly.

Method:
  1. Load existing §4.1.31 dual-channel CSV (360 realizations).
  2. Compute new ratio/gap features from raw eigenvalues.
  3. OLS-residualize against density (n_causal_pairs).
  4. Power-law α grid scan: R²(residual, H^α) for α ∈ [0.25, 8.0].
  5. Compare α values with §4.1.32 baselines (α≈1 at d=4 for raw features).

New features tested:
  - eig_ratio        = eig_min / eig_max
  - eig_abs_ratio    = |eig_min| / |eig_max|
  - eig_gap_norm     = eig_gap / eig_spread
  - eig_gap_over_max = eig_gap / |eig_max|
  - eig_neg_pos_gap  = |eig_min| - eig_max  (asymmetry)
  - eig_neg_pos_prod = |eig_min| * eig_max   (product → H²?)
  - b1_std_over_mean = b1_std / |b1_mean|    (coefficient of variation)
  - b1_std_sq        = b1_std²                (variance → H²?)
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats


# ---------------------------------------------------------------------------
# Helpers (same as §4.1.32)
# ---------------------------------------------------------------------------
def ols_residualize(feat: np.ndarray, dens: np.ndarray) -> np.ndarray:
    coeffs = np.polyfit(dens, feat, 1)
    return feat - np.polyval(coeffs, dens)


def ols_r2(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.std(x) < 1e-15 or np.std(y) < 1e-15:
        return 0.0
    coeffs = np.polyfit(x, y, 1)
    pred = np.polyval(coeffs, x)
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return max(0.0, 1.0 - ss_res / ss_tot) if ss_tot > 1e-15 else 0.0


def alpha_grid(H_arr: np.ndarray, Y_arr: np.ndarray) -> dict:
    alpha_range = np.arange(0.25, 8.05, 0.25)
    results = {}
    for alpha in alpha_range:
        x = np.power(np.abs(H_arr), alpha)
        results[round(alpha, 2)] = ols_r2(x, Y_arr)
    return results


def group_means(H_arr: np.ndarray, Y_arr: np.ndarray):
    H_levels = sorted(set(np.round(H_arr, 6)))
    means, counts = [], []
    H_out = []
    for h in H_levels:
        mask = np.abs(H_arr - h) < 1e-6
        if mask.sum() > 0:
            means.append(np.mean(Y_arr[mask]))
            counts.append(int(mask.sum()))
            H_out.append(h)
    return np.array(H_out), np.array(means), np.array(counts)


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
def load_csv(path: Path) -> list[dict]:
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = []
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
# Compute new spectral ratio features
# ---------------------------------------------------------------------------
NEW_FEATURES = [
    "eig_ratio",
    "eig_abs_ratio",
    "eig_gap_norm",
    "eig_gap_over_max",
    "eig_neg_pos_gap",
    "eig_neg_pos_prod",
    "b1_std_over_mean",
    "b1_std_sq",
]

# Baseline features from §4.1.32 for comparison
BASELINE_FEATURES = ["w_max_ratio", "b1_std", "eig_spread", "mean_layer_width"]


def compute_new_features(row: dict) -> dict:
    """Compute new ratio/gap features from existing eigenvalue data."""
    eig_min = row.get("eig_min", float("nan"))
    eig_max = row.get("eig_max", float("nan"))
    eig_gap = row.get("eig_gap", float("nan"))
    eig_spread = row.get("eig_spread", float("nan"))
    b1_std = row.get("b1_std", float("nan"))
    b1_mean = row.get("b1_mean", float("nan"))

    out = {}

    # Eigenvalue ratios
    if abs(eig_max) > 1e-15:
        out["eig_ratio"] = eig_min / eig_max
        out["eig_abs_ratio"] = abs(eig_min) / abs(eig_max)
        out["eig_gap_over_max"] = eig_gap / abs(eig_max)
    else:
        out["eig_ratio"] = float("nan")
        out["eig_abs_ratio"] = float("nan")
        out["eig_gap_over_max"] = float("nan")

    # Normalized gap
    if abs(eig_spread) > 1e-15:
        out["eig_gap_norm"] = eig_gap / eig_spread
    else:
        out["eig_gap_norm"] = float("nan")

    # Asymmetry and product
    out["eig_neg_pos_gap"] = abs(eig_min) - eig_max
    out["eig_neg_pos_prod"] = abs(eig_min) * eig_max

    # B1 coefficient of variation
    if abs(b1_mean) > 1e-15:
        out["b1_std_over_mean"] = b1_std / abs(b1_mean)
    else:
        out["b1_std_over_mean"] = float("nan")

    # B1 variance (squared std)
    out["b1_std_sq"] = b1_std ** 2

    return out


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------
def main():
    data_path = Path(__file__).parent / "outputs_unified_functional" / "conjecture_e_dual_channel.csv"
    out_dir = Path(__file__).parent / "outputs_unified_functional"

    print(f"Loading data from {data_path} ...")
    raw_rows = load_csv(data_path)
    print(f"  Loaded {len(raw_rows)} rows.")

    # Augment with new features
    for row in raw_rows:
        row.update(compute_new_features(row))

    dims = sorted(set(int(r["d"]) for r in raw_rows))
    ns = sorted(set(int(r["N"]) for r in raw_rows))
    hubbles = sorted(set(r["hubble"] for r in raw_rows))

    all_features = NEW_FEATURES + BASELINE_FEATURES
    lines: list[str] = []

    lines.append("# T4: Spectral Ratio Bridge — EH Reanalysis\n")
    lines.append("## Goal\n")
    lines.append("Test whether eigenvalue RATIOS/GAPS of B_ℓ achieve α ≈ 2")
    lines.append("(tracking R = d(d-1)H²) instead of α ≈ 1 (tracking H).\n")
    lines.append(f"- Data source: §4.1.31 dual-channel CSV ({len(raw_rows)} realizations)")
    lines.append(f"- Dimensions: {dims}")
    lines.append(f"- Sizes: {ns}")
    lines.append(f"- Hubble values: {hubbles}")
    lines.append(f"- New features: {NEW_FEATURES}")
    lines.append(f"- Baseline features: {BASELINE_FEATURES}\n")

    # ==================================================================
    # Section 1: Raw correlation with H² (before density removal)
    # ==================================================================
    lines.append("## 1. Raw Spearman ρ(feature, H²) — pooled by d\n")
    lines.append("| d | " + " | ".join(all_features) + " |")
    lines.append("|---|" + "|".join(["------"] * len(all_features)) + "|")

    for d in dims:
        subset = [r for r in raw_rows if int(r["d"]) == d]
        h2_arr = np.array([r["H2"] for r in subset])
        cells = []
        for feat in all_features:
            fa = np.array([r.get(feat, float("nan")) for r in subset])
            mask = ~np.isnan(fa)
            if mask.sum() < 10:
                cells.append("N/A")
                continue
            rho, p = sp_stats.spearmanr(h2_arr[mask], fa[mask])
            sig = "**" if p < 0.01 else "*" if p < 0.05 else ""
            cells.append(f"{sig}{rho:+.3f}{sig}")
        lines.append(f"| {d} | " + " | ".join(cells) + " |")

    # ==================================================================
    # Section 2: Density-residual Spearman ρ per (d, N)
    # ==================================================================
    lines.append("\n## 2. Density-Residual ρ(residual, H²) per (d, N)\n")
    lines.append("After OLS-removing n_causal_pairs.\n")

    # Store residuals for α scan
    residual_store: dict[tuple[int, int], dict[str, tuple[np.ndarray, np.ndarray]]] = {}

    for d in dims:
        for N in ns:
            subset = [r for r in raw_rows if int(r["d"]) == d and int(r["N"]) == N]
            if len(subset) < 15:
                continue

            h2_arr = np.array([r["H2"] for r in subset])
            H_arr = np.array([r["hubble"] for r in subset])
            dens_arr = np.array([r["n_causal_pairs"] for r in subset])

            residual_store[(d, N)] = {}

            for feat in all_features:
                fa = np.array([r.get(feat, float("nan")) for r in subset])
                mask = ~np.isnan(fa)
                if mask.sum() < 15:
                    continue
                resid = ols_residualize(fa[mask], dens_arr[mask])
                residual_store[(d, N)][feat] = (H_arr[mask], resid)

    # Print table
    lines.append("| d | N | " + " | ".join(all_features) + " |")
    lines.append("|---|---|" + "|".join(["------"] * len(all_features)) + "|")

    for d in dims:
        for N in ns:
            key = (d, N)
            if key not in residual_store:
                continue
            cells = []
            for feat in all_features:
                if feat not in residual_store[key]:
                    cells.append("N/A")
                    continue
                H_arr, resid = residual_store[key][feat]
                h2 = H_arr ** 2
                rho, p = sp_stats.spearmanr(h2, resid)
                sig = "**" if abs(rho) > 0.3 and p < 0.01 else ""
                cells.append(f"{sig}{rho:+.3f}{sig}")
            lines.append(f"| {d} | {N} | " + " | ".join(cells) + " |")

    # ==================================================================
    # Section 3: α GRID SCAN (the key test)
    # ==================================================================
    lines.append("\n## 3. Power-Law α Grid Scan (KEY TEST)\n")
    lines.append("For each feature, find the α that maximizes R²(group-mean residual, H^α).")
    lines.append("Target: α ≈ 2 at d=4 → feature tracks R = d(d-1)H².")
    lines.append("Baseline from §4.1.32: w_max_ratio α≈1.25, b1_std α≈1.00 at d=4.\n")

    # Pooled analysis: for each d, pool all N values and compute group means per H
    lines.append("### 3a. Pooled group-mean analysis (all N combined)\n")
    lines.append("| d | Feature | α_best | R²_best | α_2nd | R²_2nd | α=1.0 R² | α=2.0 R² |")
    lines.append("|---|---------|--------|---------|-------|--------|----------|----------|")

    pooled_results: dict[int, dict[str, dict]] = {}

    for d in dims:
        pooled_results[d] = {}
        for feat in all_features:
            # Collect all residuals across N for this (d, feat)
            all_H = []
            all_resid = []
            for N in ns:
                key = (d, N)
                if key in residual_store and feat in residual_store[key]:
                    H_arr, resid = residual_store[key][feat]
                    all_H.extend(H_arr.tolist())
                    all_resid.extend(resid.tolist())

            if len(all_H) < 15:
                lines.append(f"| {d} | {feat} | N/A | N/A | N/A | N/A | N/A | N/A |")
                continue

            H_arr = np.array(all_H)
            resid_arr = np.array(all_resid)

            # Remove H=0 (undefined for H^α)
            nonzero = H_arr > 1e-10
            if nonzero.sum() < 10:
                lines.append(f"| {d} | {feat} | N/A | N/A | N/A | N/A | N/A | N/A |")
                continue

            H_nz = H_arr[nonzero]
            R_nz = resid_arr[nonzero]

            # Group means
            H_gm, R_gm, _ = group_means(H_nz, R_nz)

            if len(H_gm) < 3:
                lines.append(f"| {d} | {feat} | N/A | N/A | N/A | N/A | N/A | N/A |")
                continue

            # α grid scan on group means
            ag = alpha_grid(H_gm, R_gm)
            sorted_ag = sorted(ag.items(), key=lambda x: -x[1])

            alpha_best, r2_best = sorted_ag[0]
            alpha_2nd, r2_2nd = sorted_ag[1] if len(sorted_ag) > 1 else (0, 0)
            r2_at_1 = ag.get(1.0, 0.0)
            r2_at_2 = ag.get(2.0, 0.0)

            pooled_results[d][feat] = {
                "alpha_best": alpha_best,
                "r2_best": r2_best,
                "r2_at_1": r2_at_1,
                "r2_at_2": r2_at_2,
            }

            lines.append(
                f"| {d} | {feat} | **{alpha_best:.2f}** | {r2_best:.4f} "
                f"| {alpha_2nd:.2f} | {r2_2nd:.4f} "
                f"| {r2_at_1:.4f} | {r2_at_2:.4f} |"
            )

    # ==================================================================
    # Section 3b: Per-(d, N) α scan for N-scaling
    # ==================================================================
    lines.append("\n### 3b. Per-(d, N) α scan — N-scaling of α_eff\n")
    lines.append("Does α_eff increase with N at d=4? (Tests T5 hypothesis)\n")

    focus_feats = ["eig_ratio", "eig_abs_ratio", "eig_neg_pos_prod",
                   "b1_std_sq", "w_max_ratio", "b1_std"]

    lines.append("| d | N | " + " | ".join(f"α({f})" for f in focus_feats) + " |")
    lines.append("|---|---|" + "|".join(["------"] * len(focus_feats)) + "|")

    for d in dims:
        for N in ns:
            key = (d, N)
            if key not in residual_store:
                continue
            cells = []
            for feat in focus_feats:
                if feat not in residual_store[key]:
                    cells.append("N/A")
                    continue
                H_arr, resid = residual_store[key][feat]
                nonzero = H_arr > 1e-10
                if nonzero.sum() < 8:
                    cells.append("N/A")
                    continue
                H_nz = H_arr[nonzero]
                R_nz = resid[nonzero]
                H_gm, R_gm, _ = group_means(H_nz, R_nz)
                if len(H_gm) < 3:
                    cells.append("N/A")
                    continue
                ag = alpha_grid(H_gm, R_gm)
                alpha_best = max(ag, key=ag.get)
                cells.append(f"{alpha_best:.2f}")
            lines.append(f"| {d} | {N} | " + " | ".join(cells) + " |")

    # ==================================================================
    # Section 4: Direct correlation with R_dS
    # ==================================================================
    lines.append("\n## 4. Direct Correlation: Residual vs R_dS = d(d-1)H²\n")
    lines.append("| d | Feature | Spearman(resid, R_dS) | Pearson R²(resid, R_dS) | Pearson R²(resid, H) | R²(R)/R²(H) |")
    lines.append("|---|---------|----------------------|------------------------|---------------------|-------------|")

    for d in dims:
        for feat in all_features:
            all_H = []
            all_resid = []
            for N in ns:
                key = (d, N)
                if key in residual_store and feat in residual_store[key]:
                    H_arr, resid = residual_store[key][feat]
                    all_H.extend(H_arr.tolist())
                    all_resid.extend(resid.tolist())

            if len(all_H) < 15:
                continue

            H_arr = np.array(all_H)
            resid_arr = np.array(all_resid)
            R_dS = d * (d - 1) * H_arr ** 2

            rho_R, _ = sp_stats.spearmanr(R_dS, resid_arr)
            r2_R = ols_r2(R_dS, resid_arr)
            r2_H = ols_r2(H_arr, resid_arr)
            ratio = r2_R / r2_H if r2_H > 1e-10 else float("inf")

            lines.append(
                f"| {d} | {feat} | {rho_R:+.3f} | {r2_R:.4f} | {r2_H:.4f} | {ratio:.2f}× |"
            )

    # ==================================================================
    # Section 5: Verdict
    # ==================================================================
    lines.append("\n## 5. Verdict\n")

    # Check if any new feature achieves α ≈ 2 at d=4
    d4_hits = []
    if 4 in pooled_results:
        for feat, res in pooled_results[4].items():
            if feat in NEW_FEATURES and 1.75 <= res["alpha_best"] <= 2.25:
                d4_hits.append((feat, res["alpha_best"], res["r2_best"]))

    if d4_hits:
        lines.append("### ✅ α ≈ 2 CANDIDATES FOUND at d=4\n")
        lines.append("| Feature | α_best | R²_best |")
        lines.append("|---------|--------|---------|")
        for feat, alpha, r2 in d4_hits:
            lines.append(f"| **{feat}** | **{alpha:.2f}** | **{r2:.4f}** |")
        lines.append("")
        lines.append("These features are candidates for a finite-N R-tracking observable.")
        lines.append("If confirmed, this would substantially advance E-bulk-second-order.")
    else:
        lines.append("### ❌ No new feature achieves α ≈ 2 at d=4\n")
        lines.append("All spectral ratio/gap features still track H (α ≈ 1), not R (α ≈ 2).")
        lines.append("This **confirms** the §4.1.33 conclusion: the H→R bridge is a")
        lines.append("continuum-limit phenomenon, not achievable by algebraic ratios")
        lines.append("of finite-N eigenvalue features.\n")

    # Summary comparison table
    lines.append("### Summary: α comparison (d=4, pooled)\n")
    lines.append("| Feature | Type | α_best | R²_best | vs baseline |")
    lines.append("|---------|------|--------|---------|-------------|")

    if 4 in pooled_results:
        for feat in all_features:
            if feat in pooled_results[4]:
                res = pooled_results[4][feat]
                ftype = "NEW" if feat in NEW_FEATURES else "baseline"
                vs = ""
                if ftype == "NEW":
                    if res["alpha_best"] > 1.5:
                        vs = "↑ closer to R"
                    elif res["alpha_best"] < 0.75:
                        vs = "↓ below H"
                    else:
                        vs = "≈ baseline"
                lines.append(
                    f"| {feat} | {ftype} | {res['alpha_best']:.2f} "
                    f"| {res['r2_best']:.4f} | {vs} |"
                )

    # N-scaling verdict
    lines.append("\n### N-scaling at d=4\n")
    d4_scaling = {}
    for feat in focus_feats:
        alphas = []
        for N in ns:
            key = (4, N)
            if key not in residual_store or feat not in residual_store[key]:
                continue
            H_arr, resid = residual_store[key][feat]
            nonzero = H_arr > 1e-10
            if nonzero.sum() < 8:
                continue
            H_nz = H_arr[nonzero]
            R_nz = resid[nonzero]
            H_gm, R_gm, _ = group_means(H_nz, R_nz)
            if len(H_gm) < 3:
                continue
            ag = alpha_grid(H_gm, R_gm)
            alpha_best = max(ag, key=ag.get)
            alphas.append((N, alpha_best))

        if len(alphas) >= 2:
            ns_arr = np.array([a[0] for a in alphas], dtype=float)
            al_arr = np.array([a[1] for a in alphas])
            rho_trend, p_trend = sp_stats.spearmanr(ns_arr, al_arr)
            d4_scaling[feat] = (alphas, rho_trend, p_trend)

    if d4_scaling:
        lines.append("| Feature | N values | α values | ρ(N, α) | Trend |")
        lines.append("|---------|----------|----------|---------|-------|")
        for feat, (alphas, rho_trend, p_trend) in d4_scaling.items():
            ns_str = "/".join(str(a[0]) for a in alphas)
            al_str = "/".join(f"{a[1]:.2f}" for a in alphas)
            trend = "↑ converging to 2" if rho_trend > 0.5 else "→ flat" if abs(rho_trend) < 0.3 else "↓ diverging"
            lines.append(f"| {feat} | {ns_str} | {al_str} | {rho_trend:+.2f} | {trend} |")
    else:
        lines.append("Insufficient data for N-scaling analysis at d=4.")

    lines.append("\n---\n")
    lines.append("*Generated by conjecture_e_spectral_ratio_bridge.py*")
    lines.append(f"*Data: {len(raw_rows)} realizations from §4.1.31*")

    # Write report
    report_path = out_dir / "conjecture_e_spectral_ratio_bridge.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nReport written to {report_path}")


if __name__ == "__main__":
    main()
