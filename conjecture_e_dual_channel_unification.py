"""Conjecture E — §4.1.31: Dual-Channel Unification Experiment.

Tests whether the antichain (transverse) and B_ℓ (spectral) channels —
the two supported DDT escape paths — converge to the same geometric content.

Design:
  - Extract BOTH antichain and B_ℓ features from the SAME de Sitter sprinklings.
  - Extended to larger N (up to 1024) to test N-scaling convergence.
  - Cross-channel analysis: correlate density-residuals of antichain features
    with density-residuals of B_ℓ features.
  - If the two channels' residuals are highly correlated, they measure the
    same post-density geometric content → unified E-bulk proxy.
  - If uncorrelated, they are independent signals → E-bulk is multi-dimensional.

Key Questions:
  Q1: Do both channels independently correlate with H² at each (d, N)?
  Q2: After density removal, do both channels' residuals correlate with each other?
  Q3: How do cross-channel correlations scale with N?
  Q4: Can a combined proxy (antichain + B_ℓ) explain more curvature variance than either alone?

References:
  - §4.1.27: B_ℓ spectral — 6/18 beyond density
  - §4.1.28: Antichain — 21/21 beyond density
  - §4.1.30: DDT verdict — vertical closed, horizontal+spectral open
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

from curvature_layer2_de_sitter_scan import sprinkle_de_sitter_like_diamond
from curvature_layer2_recovery_fast import build_causal_matrix_de_sitter

# Re-use infrastructure from §4.1.27 and §4.1.28
from conjecture_e_sorkin_dalembertian import (
    build_dalembertian_matrix,
    compute_b1_features,
)
from conjecture_e_antichain_structure import (
    compute_antichain_features,
)


# ---------------------------------------------------------------------------
# Joint data row: both channels on the same realization
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ExpRow:
    d: int
    N: int
    hubble: float
    rep: int
    R_dS: float
    H2: float
    n_causal_pairs: int
    # --- Antichain features (transverse channel) ---
    w_max_ratio: float
    n_layers: float
    layer_ratio: float
    mean_layer_width: float
    layer_width_std: float
    layer_width_cv: float
    max_layer_width_ratio: float
    layer_entropy: float
    # --- B_ℓ spectral features ---
    b1_mean: float
    b1_std: float
    b1_neg_frac: float
    eig_min: float
    eig_max: float
    eig_gap: float
    eig_spread: float


def run_single(
    d: int, N: int, hubble: float, rep: int, seed: int
) -> ExpRow:
    """Run one realization: sprinkle once, extract both channel features."""
    points = sprinkle_de_sitter_like_diamond(
        N, d - 1, hubble=hubble,
        seed=seed + d * 100000 + N * 100 + rep + int(hubble * 1000),
    )
    causal = build_causal_matrix_de_sitter(points, hubble)

    n_causal_pairs = int(np.sum(causal))
    R_dS = d * (d - 1) * hubble ** 2
    H2 = hubble ** 2

    # --- Antichain features ---
    ac_feats = compute_antichain_features(causal, compute_dilworth=(N <= 600))

    # --- B_ℓ spectral features ---
    # Build d'Alembertian (without density prefactor — preserves rank correlations)
    B = build_dalembertian_matrix(causal, d, rho=1.0)
    b1_feats = compute_b1_features(B)

    return ExpRow(
        d=d, N=N, hubble=hubble, rep=rep,
        R_dS=R_dS, H2=H2,
        n_causal_pairs=n_causal_pairs,
        # Antichain
        w_max_ratio=ac_feats.get("w_max_ratio", float("nan")),
        n_layers=ac_feats.get("n_layers", float("nan")),
        layer_ratio=ac_feats.get("layer_ratio", float("nan")),
        mean_layer_width=ac_feats.get("mean_layer_width", float("nan")),
        layer_width_std=ac_feats.get("layer_width_std", float("nan")),
        layer_width_cv=ac_feats.get("layer_width_cv", float("nan")),
        max_layer_width_ratio=ac_feats.get("max_layer_width_ratio", float("nan")),
        layer_entropy=ac_feats.get("layer_entropy", float("nan")),
        # B_ℓ spectral
        b1_mean=b1_feats.get("b1_mean", float("nan")),
        b1_std=b1_feats.get("b1_std", float("nan")),
        b1_neg_frac=b1_feats.get("b1_neg_frac", float("nan")),
        eig_min=b1_feats.get("eig_min", float("nan")),
        eig_max=b1_feats.get("eig_max", float("nan")),
        eig_gap=b1_feats.get("eig_gap", float("nan")),
        eig_spread=b1_feats.get("eig_spread", float("nan")),
    )


# ---------------------------------------------------------------------------
# Density residualization helper
# ---------------------------------------------------------------------------
def ols_residualize(feat_a: np.ndarray, dens_a: np.ndarray) -> np.ndarray:
    """OLS-remove density from feature: feat = a*dens + b → residual."""
    coeffs = np.polyfit(dens_a, feat_a, 1)
    return feat_a - np.polyval(coeffs, dens_a)


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------
def generate_report(
    rows: list[ExpRow],
    dims: list[int],
    ns: list[int],
    hubbles: list[float],
) -> str:
    lines: list[str] = []
    lines.append("# §4.1.31: Dual-Channel Unification Experiment\n")
    lines.append("## Design\n")
    lines.append(f"- Dimensions: {dims}")
    lines.append(f"- Sizes: {ns}")
    lines.append(f"- Hubble values: {hubbles}")
    lines.append(f"- Total realizations: {len(rows)}")
    lines.append("- Both antichain and B_ℓ features extracted from the **same** sprinklings")
    lines.append("- Core question: do the two supported DDT escape channels measure the same geometric content?\n")

    ac_feats = [
        "w_max_ratio", "layer_ratio", "mean_layer_width",
        "layer_width_std", "layer_width_cv", "max_layer_width_ratio",
        "layer_entropy",
    ]
    bl_feats = ["b1_mean", "b1_std", "b1_neg_frac", "eig_min", "eig_max", "eig_gap", "eig_spread"]

    # ======================================================================
    # Q1: Independent channel performance vs H² (replication check)
    # ======================================================================
    lines.append("## Q1: Independent Channel Performance vs H²\n")
    lines.append("Replication of §4.1.27/28 on the same sprinklings.\n")

    lines.append("### Antichain features vs H² (pooled by d)\n")
    lines.append("| d | " + " | ".join(ac_feats) + " |")
    lines.append("|---|" + "|".join(["------"] * len(ac_feats)) + "|")
    for d in dims:
        h2_arr = []
        vals = {f: [] for f in ac_feats}
        for r in rows:
            if r.d != d:
                continue
            h2_arr.append(r.H2)
            for f in ac_feats:
                vals[f].append(getattr(r, f))
        if len(h2_arr) < 10:
            continue
        h2_a = np.array(h2_arr)
        cells = []
        for f in ac_feats:
            fa = np.array(vals[f])
            mask = ~np.isnan(fa)
            if mask.sum() < 10:
                cells.append("N/A")
                continue
            rho, p = sp_stats.spearmanr(h2_a[mask], fa[mask])
            sig = "**" if p < 0.01 else "*" if p < 0.05 else ""
            cells.append(f"{sig}{rho:+.3f}{sig}")
        lines.append(f"| {d} | " + " | ".join(cells) + " |")

    lines.append("\n### B_ℓ spectral features vs H² (pooled by d)\n")
    lines.append("| d | " + " | ".join(bl_feats) + " |")
    lines.append("|---|" + "|".join(["------"] * len(bl_feats)) + "|")
    for d in dims:
        h2_arr = []
        vals = {f: [] for f in bl_feats}
        for r in rows:
            if r.d != d:
                continue
            h2_arr.append(r.H2)
            for f in bl_feats:
                vals[f].append(getattr(r, f))
        if len(h2_arr) < 10:
            continue
        h2_a = np.array(h2_arr)
        cells = []
        for f in bl_feats:
            fa = np.array(vals[f])
            mask = ~np.isnan(fa)
            if mask.sum() < 10:
                cells.append("N/A")
                continue
            rho, p = sp_stats.spearmanr(h2_a[mask], fa[mask])
            sig = "**" if p < 0.01 else "*" if p < 0.05 else ""
            cells.append(f"{sig}{rho:+.3f}{sig}")
        lines.append(f"| {d} | " + " | ".join(cells) + " |")

    # ======================================================================
    # Q2: Density-Residual Analysis per channel (per d,N slice)
    # ======================================================================
    lines.append("\n## Q2: Density-Residual Analysis per Channel per (d, N)\n")
    lines.append("After OLS-removing n_causal_pairs, Spearman ρ(residual, H²).\n")

    # Best representatives: antichain → w_max_ratio, mean_layer_width, layer_width_std
    # B_ℓ → b1_std, eig_min, eig_max, eig_spread
    ac_key = ["w_max_ratio", "mean_layer_width", "layer_width_std", "layer_ratio"]
    bl_key = ["b1_std", "eig_min", "eig_max", "eig_spread"]
    all_key = ac_key + bl_key

    lines.append("| d | N | " + " | ".join(all_key) + " |")
    lines.append("|---|---|" + "|".join(["------"] * len(all_key)) + "|")

    # Store residual arrays for cross-channel analysis
    residuals_store: dict[tuple[int, int], dict[str, np.ndarray]] = {}
    h2_store: dict[tuple[int, int], np.ndarray] = {}

    for d in dims:
        for N in ns:
            subset = [r for r in rows if r.d == d and r.N == N]
            if len(subset) < 15:
                continue

            h2_a = np.array([r.H2 for r in subset])
            dens_a = np.array([float(r.n_causal_pairs) for r in subset])
            h2_store[(d, N)] = h2_a

            cells = []
            resid_dict: dict[str, np.ndarray] = {}
            for feat in all_key:
                fa = np.array([getattr(r, feat) for r in subset])
                mask = ~np.isnan(fa)
                if mask.sum() < 15:
                    cells.append("N/A")
                    continue
                # Use only non-NaN entries
                resid = ols_residualize(fa[mask], dens_a[mask])
                rho_r, p_r = sp_stats.spearmanr(h2_a[mask], resid)
                resid_dict[feat] = resid
                sig = "**" if abs(rho_r) > 0.3 and p_r < 0.01 else ""
                cells.append(f"{sig}{rho_r:+.3f}{sig}")

            residuals_store[(d, N)] = resid_dict
            lines.append(f"| {d} | {N} | " + " | ".join(cells) + " |")

    # Count beyond-density per channel per (d, N)
    lines.append("\n### Summary: Beyond-density count per (d, N)\n")
    lines.append("| d | N | AC beyond | B_ℓ beyond | Total |")
    lines.append("|---|---|-----------|-----------|-------|")

    for d in dims:
        for N in ns:
            key = (d, N)
            if key not in residuals_store:
                continue
            rd = residuals_store[key]
            h2_a = h2_store[key]
            ac_beyond = 0
            bl_beyond = 0
            for feat in ac_key:
                if feat in rd:
                    rho_r, p_r = sp_stats.spearmanr(h2_a[:len(rd[feat])], rd[feat])
                    if abs(rho_r) > 0.3 and p_r < 0.01:
                        ac_beyond += 1
            for feat in bl_key:
                if feat in rd:
                    rho_r, p_r = sp_stats.spearmanr(h2_a[:len(rd[feat])], rd[feat])
                    if abs(rho_r) > 0.3 and p_r < 0.01:
                        bl_beyond += 1
            lines.append(f"| {d} | {N} | {ac_beyond}/{len(ac_key)} | {bl_beyond}/{len(bl_key)} | {ac_beyond + bl_beyond}/{len(all_key)} |")

    # ======================================================================
    # Q3: CROSS-CHANNEL CORRELATION (the key test)
    # ======================================================================
    lines.append("\n## Q3: Cross-Channel Residual Correlation (KEY TEST)\n")
    lines.append("After removing density from both channels independently,")
    lines.append("do their residuals correlate with each other?\n")
    lines.append("- High |ρ| → same underlying geometric content (unified E-bulk)")
    lines.append("- Low |ρ| → independent signals (multi-dimensional E-bulk)\n")

    # For each (d, N), correlate best antichain residual with best B_ℓ residual
    lines.append("### Best-representative cross-correlation\n")
    lines.append("AC representative: w_max_ratio (strongest antichain)")
    lines.append("B_ℓ representative: b1_std (strongest B_ℓ spectral)\n")
    lines.append("| d | N | ρ(AC_resid, Bℓ_resid) | p-value | interpretation |")
    lines.append("|---|---|----------------------|---------|----------------|")

    cross_results: list[tuple[int, int, float, float]] = []

    for d in dims:
        for N in ns:
            key = (d, N)
            if key not in residuals_store:
                continue
            rd = residuals_store[key]
            # Pick best available from each channel
            ac_best = None
            for ac_cand in ["w_max_ratio", "mean_layer_width", "layer_width_std"]:
                if ac_cand in rd and len(rd[ac_cand]) >= 15:
                    ac_best = rd[ac_cand]
                    break
            bl_best = None
            for bl_cand in ["b1_std", "eig_min", "eig_max", "eig_spread"]:
                if bl_cand in rd and len(rd[bl_cand]) >= 15:
                    bl_best = rd[bl_cand]
                    break
            if ac_best is None or bl_best is None:
                continue
            # Ensure same length
            min_len = min(len(ac_best), len(bl_best))
            rho_cross, p_cross = sp_stats.spearmanr(ac_best[:min_len], bl_best[:min_len])

            if abs(rho_cross) > 0.5:
                interp = "**UNIFIED** (same content)"
            elif abs(rho_cross) > 0.3:
                interp = "partial overlap"
            else:
                interp = "independent signals"

            cross_results.append((d, N, rho_cross, p_cross))
            lines.append(f"| {d} | {N} | {rho_cross:+.3f} | {p_cross:.2e} | {interp} |")

    # Full cross-correlation matrix for each (d, N)
    lines.append("\n### Full cross-correlation matrix (d=4 slices)\n")
    for d in dims:
        for N in ns:
            key = (d, N)
            if key not in residuals_store:
                continue
            rd = residuals_store[key]
            available = [f for f in all_key if f in rd and len(rd[f]) >= 15]
            if len(available) < 4:
                continue

            lines.append(f"\n#### d={d}, N={N}\n")
            lines.append("| | " + " | ".join(available) + " |")
            lines.append("|---|" + "|".join(["------"] * len(available)) + "|")
            for fi in available:
                cells = []
                for fj in available:
                    min_len = min(len(rd[fi]), len(rd[fj]))
                    if min_len < 10:
                        cells.append("N/A")
                        continue
                    rho_ij, _ = sp_stats.spearmanr(rd[fi][:min_len], rd[fj][:min_len])
                    cells.append(f"{rho_ij:+.3f}")
                lines.append(f"| {fi} | " + " | ".join(cells) + " |")

    # ======================================================================
    # Q4: N-Scaling of Cross-Channel Correlation
    # ======================================================================
    lines.append("\n## Q4: N-Scaling of Cross-Channel Correlation\n")
    lines.append("Does the cross-channel correlation strengthen or weaken with N?\n")
    lines.append("| d | N | ρ_cross(AC, Bℓ) |")
    lines.append("|---|---|-----------------|")

    for d in dims:
        for cr in cross_results:
            if cr[0] == d:
                lines.append(f"| {d} | {cr[1]} | {cr[2]:+.3f} |")

    # Check trend
    for d in dims:
        d_results = [(cr[1], cr[2]) for cr in cross_results if cr[0] == d]
        if len(d_results) >= 3:
            ns_arr = [x[0] for x in d_results]
            rhos_arr = [x[1] for x in d_results]
            trend_rho, trend_p = sp_stats.spearmanr(ns_arr, [abs(r) for r in rhos_arr])
            trend_word = "strengthening" if trend_rho > 0 else "weakening" if trend_rho < 0 else "flat"
            lines.append(f"\n**d={d} trend**: ρ(N, |cross_corr|) = {trend_rho:+.3f} (p={trend_p:.2e}) → **{trend_word}**")

    # ======================================================================
    # Q5: Combined Proxy — Does adding both channels explain more?
    # ======================================================================
    lines.append("\n\n## Q5: Combined Proxy Analysis\n")
    lines.append("Can a linear combination of AC + B_ℓ residuals explain more H² variance than either alone?\n")
    lines.append("| d | N | R²(AC only) | R²(Bℓ only) | R²(combined) | ΔR² |")
    lines.append("|---|---|------------|-------------|--------------|------|")

    for d in dims:
        for N in ns:
            key = (d, N)
            if key not in residuals_store:
                continue
            rd = residuals_store[key]
            h2_a = h2_store[key]

            # Find best available features
            ac_best_name = None
            for ac_cand in ["w_max_ratio", "mean_layer_width", "layer_width_std"]:
                if ac_cand in rd and len(rd[ac_cand]) >= 15:
                    ac_best_name = ac_cand
                    break
            bl_best_name = None
            for bl_cand in ["b1_std", "eig_min", "eig_max", "eig_spread"]:
                if bl_cand in rd and len(rd[bl_cand]) >= 15:
                    bl_best_name = bl_cand
                    break

            if ac_best_name is None or bl_best_name is None:
                continue

            ac_r = rd[ac_best_name]
            bl_r = rd[bl_best_name]
            min_len = min(len(ac_r), len(bl_r), len(h2_a))
            ac_r = ac_r[:min_len]
            bl_r = bl_r[:min_len]
            h2_s = h2_a[:min_len]

            if min_len < 15:
                continue

            # Single-feature OLS R²
            def ols_r2(x: np.ndarray, y: np.ndarray) -> float:
                coeffs = np.polyfit(x, y, 1)
                pred = np.polyval(coeffs, x)
                ss_res = np.sum((y - pred) ** 2)
                ss_tot = np.sum((y - np.mean(y)) ** 2)
                return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

            # Multivariate OLS R²
            def multi_ols_r2(X: np.ndarray, y: np.ndarray) -> float:
                X_aug = np.column_stack([X, np.ones(len(y))])
                try:
                    beta, _, _, _ = np.linalg.lstsq(X_aug, y, rcond=None)
                    pred = X_aug @ beta
                    ss_res = np.sum((y - pred) ** 2)
                    ss_tot = np.sum((y - np.mean(y)) ** 2)
                    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
                except np.linalg.LinAlgError:
                    return float("nan")

            r2_ac = ols_r2(ac_r, h2_s)
            r2_bl = ols_r2(bl_r, h2_s)
            X_combined = np.column_stack([ac_r, bl_r])
            r2_both = multi_ols_r2(X_combined, h2_s)
            delta_r2 = r2_both - max(r2_ac, r2_bl)

            lines.append(
                f"| {d} | {N} | {r2_ac:.3f} | {r2_bl:.3f} | {r2_both:.3f} | {delta_r2:+.3f} |"
            )

    # ======================================================================
    # Conclusion
    # ======================================================================
    lines.append("\n## Conclusion\n")

    # Classify cross-channel results
    unified_count = sum(1 for cr in cross_results if abs(cr[2]) > 0.5)
    partial_count = sum(1 for cr in cross_results if 0.3 < abs(cr[2]) <= 0.5)
    independent_count = sum(1 for cr in cross_results if abs(cr[2]) <= 0.3)
    total_cross = len(cross_results)

    lines.append(f"**Cross-channel correlation summary**: {total_cross} (d,N) slices tested")
    lines.append(f"- Unified (|ρ| > 0.5): {unified_count}/{total_cross}")
    lines.append(f"- Partial overlap (0.3 < |ρ| ≤ 0.5): {partial_count}/{total_cross}")
    lines.append(f"- Independent (|ρ| ≤ 0.3): {independent_count}/{total_cross}\n")

    if unified_count > total_cross * 0.6:
        lines.append("**VERDICT: CHANNELS UNIFIED** — antichain and B_ℓ spectral residuals")
        lines.append("are highly correlated after density removal, suggesting they measure")
        lines.append("the same underlying geometric content. E-bulk has a single dominant")
        lines.append("post-density degree of freedom.\n")
    elif unified_count + partial_count > total_cross * 0.5:
        lines.append("**VERDICT: PARTIAL UNIFICATION** — the two channels share significant")
        lines.append("common content but also carry independent information.")
        lines.append("E-bulk may have 1–2 post-density degrees of freedom.\n")
    else:
        lines.append("**VERDICT: INDEPENDENT CHANNELS** — antichain and B_ℓ spectral")
        lines.append("residuals are largely uncorrelated. The two DDT escape paths access")
        lines.append("different geometric information. E-bulk is multi-dimensional.\n")

    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(
        description="§4.1.31: Dual-channel unification experiment"
    )
    ap.add_argument("--dims", nargs="*", type=int, default=[2, 3, 4])
    ap.add_argument("--ns", nargs="*", type=int, default=[128, 256, 512])
    ap.add_argument("--hubbles", nargs="*", type=float,
                    default=[0.0, 0.25, 0.5, 1.0, 2.0])
    ap.add_argument("--reps", type=int, default=8)
    ap.add_argument("--seed", type=int, default=2031)
    ap.add_argument("--out",
                    default="outputs_unified_functional/conjecture_e_dual_channel.csv")
    ap.add_argument("--report",
                    default="outputs_unified_functional/conjecture_e_dual_channel.md")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[ExpRow] = []
    total = len(args.dims) * len(args.ns) * len(args.hubbles) * args.reps
    done = 0

    print(f"Dual-channel unification experiment: {total} realizations")
    print(f"  dims={args.dims}, ns={args.ns}, hubbles={args.hubbles}, reps={args.reps}")
    print()

    for d in args.dims:
        for N in args.ns:
            for hubble in args.hubbles:
                for rep in range(args.reps):
                    row = run_single(d, N, hubble, rep, args.seed)
                    rows.append(row)
                    done += 1
                    if done % 10 == 0 or done == total:
                        print(
                            f"  [{done:4d}/{total}] d={d} N={N:4d} H={hubble:.2f} "
                            f"w_max_r={row.w_max_ratio:.3f} "
                            f"b1_std={row.b1_std:.4f} "
                            f"pairs={row.n_causal_pairs}"
                        )

    # Save CSV
    fieldnames = [f.name for f in fields(ExpRow)]
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow({fn: getattr(row, fn) for fn in fieldnames})
    print(f"\nSaved CSV: {out_path}")

    # Generate report
    report_text = generate_report(rows, args.dims, args.ns, args.hubbles)
    report_path = Path(args.report)
    report_path.write_text(report_text, encoding="utf-8")
    print(f"Saved report: {report_path}")

    # Console summary
    print("\n" + "=" * 70)
    print("QUICK SUMMARY — Cross-Channel Correlations")
    print("=" * 70)

    for d in args.dims:
        for N in args.ns:
            subset = [r for r in rows if r.d == d and r.N == N]
            if len(subset) < 15:
                continue
            h2_a = np.array([r.H2 for r in subset])
            dens_a = np.array([float(r.n_causal_pairs) for r in subset])

            # Residualize w_max_ratio and b1_std
            ac_a = np.array([r.w_max_ratio for r in subset])
            bl_a = np.array([r.b1_std for r in subset])
            ac_mask = ~np.isnan(ac_a)
            bl_mask = ~np.isnan(bl_a)
            both_mask = ac_mask & bl_mask

            if both_mask.sum() < 15:
                continue

            ac_resid = ols_residualize(ac_a[both_mask], dens_a[both_mask])
            bl_resid = ols_residualize(bl_a[both_mask], dens_a[both_mask])

            rho_cross, p_cross = sp_stats.spearmanr(ac_resid, bl_resid)
            print(f"  d={d} N={N}: ρ(AC_resid, Bℓ_resid) = {rho_cross:+.3f}  (p={p_cross:.2e})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
