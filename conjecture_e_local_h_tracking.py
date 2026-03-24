"""Conjecture E — §4.1.35: Density-Matched Local H(t) Tracking.

Redesigned Phase B from §4.1.34 to avoid the density/finite-size confound.

Three independent methodological fixes:

Method 1 — Independent Patch Sprinklings (cleanest):
  Sprinkle N points independently in a narrow time window [t_c - δ, t_c + δ]
  centered at different epochs t_c (early vs late).
  Each patch is a fresh, independent sprinkling with its OWN N elements.
  No sub-poset extraction → no density imbalance.
  Compare antichain features across patches with different local H(t_c) = p/t_c.

Method 2 — Density-Matched Sub-sampling:
  From the full [t_min, t_max] sprinkling, split into early/late by median time,
  then DOWNSAMPLE the larger bin to match the smaller bin's element count.
  Random sub-sampling repeated K times → robust mean.

Method 3 — Density-Residualized Local Features:
  Pool all local bins across reps, regress out n_elements from each feature,
  then test residual vs H_local.

Physical prediction (if local H-tracking holds):
  Early epoch (large H = p/t_early) → wider antichains (w_max_ratio ↑)
  Late epoch (small H = p/t_late) → narrower antichains (w_max_ratio ↓)

References:
  - §4.1.34: Beyond de Sitter (Phase B confound identified)
  - §4.1.28: Antichain structure (de Sitter, 21/21 beyond density)
  - §4.1.31: Dual-channel unification
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

from conjecture_e_antichain_structure import compute_antichain_features
from conjecture_e_beyond_de_sitter import (
    power_law_scale_factor,
    build_causal_matrix_power_law,
    ols_residualize,
)


# ---------------------------------------------------------------------------
# Method 1: Independent Patch Sprinkling
# ---------------------------------------------------------------------------

def sprinkle_patch(
    n: int,
    d_spatial: int,
    p: float,
    t_center: float,
    t_half_width: float,
    t0: float = 1.0,
    seed: int | None = None,
) -> np.ndarray:
    """Sprinkle N points in a narrow time window [t_c - δ, t_c + δ].

    Volume weighting still applies: accept prob ∝ a(t)^{d_spatial}.
    But because the window is narrow, the acceptance rate variation is small,
    so element counts are approximately equal across patches.
    """
    t_min = max(t_center - t_half_width, 0.01)
    t_max = t_center + t_half_width

    rng = np.random.default_rng(seed)
    points = np.empty((0, 1 + d_spatial), dtype=float)

    a_tmin = power_law_scale_factor(np.array([t_min]), p, t0)[0]
    a_tmax = power_law_scale_factor(np.array([t_max]), p, t0)[0]
    a_max = max(abs(a_tmin), abs(a_tmax))

    while points.shape[0] < n:
        batch_size = max(500, 20 * n)
        t = rng.uniform(t_min, t_max, batch_size)
        x = rng.random((batch_size, d_spatial)) - 0.5
        a = power_law_scale_factor(t, p, t0)
        accept_prob = np.clip((a / a_max) ** d_spatial, 0.0, 1.0)
        accept = rng.random(batch_size) < accept_prob
        accepted = np.column_stack([t[accept], x[accept]])
        points = np.vstack([points, accepted])

    return points[:n]


@dataclass(frozen=True)
class PatchRow:
    d: int
    N: int
    p_frw: float
    rep: int
    t_center: float
    H_local: float         # p / t_center
    n_causal_pairs: int
    w_max_ratio: float
    n_layers: float
    layer_ratio: float
    mean_layer_width: float
    layer_width_std: float
    layer_width_cv: float
    max_layer_width_ratio: float
    layer_entropy: float


def run_patch(
    d: int, N: int, p_frw: float, rep: int,
    t_center: float, t_half_width: float,
    base_seed: int,
) -> PatchRow:
    """Run one independent patch sprinkling."""
    seed = (base_seed
            + d * 1000000
            + N * 1000
            + rep * 10
            + int(t_center * 100))

    points = sprinkle_patch(
        N, d - 1, p_frw,
        t_center=t_center,
        t_half_width=t_half_width,
        t0=1.0,
        seed=seed,
    )
    causal = build_causal_matrix_power_law(points, p_frw, t0=1.0)
    n_causal = int(np.sum(causal))

    ac = compute_antichain_features(causal, compute_dilworth=(N <= 300))

    return PatchRow(
        d=d, N=N, p_frw=p_frw, rep=rep,
        t_center=t_center,
        H_local=p_frw / t_center,
        n_causal_pairs=n_causal,
        w_max_ratio=ac.get("w_max_ratio", float("nan")),
        n_layers=ac.get("n_layers", float("nan")),
        layer_ratio=ac.get("layer_ratio", float("nan")),
        mean_layer_width=ac.get("mean_layer_width", float("nan")),
        layer_width_std=ac.get("layer_width_std", float("nan")),
        layer_width_cv=ac.get("layer_width_cv", float("nan")),
        max_layer_width_ratio=ac.get("max_layer_width_ratio", float("nan")),
        layer_entropy=ac.get("layer_entropy", float("nan")),
    )


# ---------------------------------------------------------------------------
# Method 2: Density-Matched Sub-sampling
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MatchedRow:
    d: int
    N: int
    p_frw: float
    rep: int
    bin_label: str          # "early" or "late"
    n_matched: int          # number of elements after matching
    H_local: float
    n_causal_pairs: int
    w_max_ratio: float
    n_layers: float
    layer_ratio: float
    mean_layer_width: float


def run_matched_bins(
    d: int, N: int, p_frw: float, rep: int,
    base_seed: int, n_subsample: int = 5,
) -> list[MatchedRow]:
    """Sprinkle full [0.1, 1.0], split early/late, downsample to match."""
    from conjecture_e_beyond_de_sitter import sprinkle_power_law_frw

    seed = (base_seed
            + d * 100000
            + N * 100
            + rep
            + int(p_frw * 1000))
    points = sprinkle_power_law_frw(
        N, d - 1, p_frw,
        t_min=0.1, t_max=1.0, t0=1.0, seed=seed,
    )
    causal = build_causal_matrix_power_law(points, p_frw, t0=1.0)
    t_vals = points[:, 0]

    t_median = np.median(t_vals)
    early_idx = np.where(t_vals <= t_median)[0]
    late_idx = np.where(t_vals > t_median)[0]

    n_early = len(early_idx)
    n_late = len(late_idx)
    n_match = min(n_early, n_late)

    if n_match < 10:
        return []

    rng = np.random.default_rng(seed + 999)
    results = []

    for label, full_idx in [("early", early_idx), ("late", late_idx)]:
        # Accumulate over n_subsample random sub-samples
        w_max_list = []
        n_layers_list = []
        lr_list = []
        mlw_list = []
        ncp_list = []

        for k in range(n_subsample):
            if len(full_idx) > n_match:
                sub_idx = rng.choice(full_idx, size=n_match, replace=False)
            else:
                sub_idx = full_idx.copy()
            sub_idx = np.sort(sub_idx)

            sub_causal = causal[np.ix_(sub_idx, sub_idx)]
            n_cp = int(np.sum(sub_causal))

            ac = compute_antichain_features(sub_causal, compute_dilworth=(n_match <= 300))
            w_max_list.append(ac.get("w_max_ratio", float("nan")))
            n_layers_list.append(ac.get("n_layers", float("nan")))
            lr_list.append(ac.get("layer_ratio", float("nan")))
            mlw_list.append(ac.get("mean_layer_width", float("nan")))
            ncp_list.append(n_cp)

        t_sub = t_vals[full_idx]
        h_local = float(np.mean(p_frw / t_sub))

        results.append(MatchedRow(
            d=d, N=N, p_frw=p_frw, rep=rep,
            bin_label=label,
            n_matched=n_match,
            H_local=h_local,
            n_causal_pairs=int(np.mean(ncp_list)),
            w_max_ratio=float(np.nanmean(w_max_list)),
            n_layers=float(np.nanmean(n_layers_list)),
            layer_ratio=float(np.nanmean(lr_list)),
            mean_layer_width=float(np.nanmean(mlw_list)),
        ))

    return results


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def generate_report(
    patch_rows: list[PatchRow],
    matched_rows: list[MatchedRow],
    dims: list[int],
    ns: list[int],
    p_values: list[float],
    t_centers: list[float],
) -> str:
    lines: list[str] = []
    lines.append("# §4.1.35: Density-Matched Local H(t) Tracking\n")
    lines.append("## Motivation\n")
    lines.append("§4.1.34 Phase B found 0/54 cells in the expected direction, but this was")
    lines.append("a **methodological confound**: early bins have fewer elements than late bins")
    lines.append("due to volume weighting ∝ t^{p(d-1)}, so finite-size effects dominated.\n")
    lines.append("This experiment redesigns the local test with three fixes:\n")
    lines.append("1. **Independent Patch Sprinklings**: Equal-N sprinklings at different epochs")
    lines.append("2. **Density-Matched Sub-sampling**: Downsample larger bin to match smaller")
    lines.append("3. **Density-Residualized Comparison**: Regress out n_elements from features\n")

    lines.append("## Design\n")
    lines.append(f"- Dimensions: {dims}")
    lines.append(f"- Sizes: {ns}")
    lines.append(f"- Power-law exponents p: {p_values}")
    lines.append(f"- Patch time centers: {t_centers}")
    lines.append(f"- Total patch rows (Method 1): {len(patch_rows)}")
    lines.append(f"- Total matched rows (Method 2): {len(matched_rows)}\n")

    ac_feats = [
        "w_max_ratio", "n_layers", "layer_ratio",
        "mean_layer_width", "layer_width_std",
    ]

    # ==================================================================
    # METHOD 1: Independent Patch Sprinklings
    # ==================================================================
    lines.append("---\n")
    lines.append("## Method 1: Independent Patch Sprinklings\n")
    lines.append("Each epoch t_c gets its own fresh N-element sprinkling in [t_c±δ].")
    lines.append("No sub-poset extraction → no density imbalance.\n")
    lines.append("**Physical prediction**: Higher H(t_c) = p/t_c → higher w_max_ratio\n")

    # M1-A: Per-p correlation of features with H_local across epochs
    lines.append("### M1-A: Feature vs H_local (Spearman, pooled across epochs)\n")
    lines.append("| d | N | p | n_rows | w_max_ratio ρ | n_layers ρ | layer_ratio ρ | mean_layer_width ρ |")
    lines.append("|---|---|---|--------|--------------|-----------|--------------|-------------------|")

    m1_summary = {d: {"positive": 0, "significant": 0, "total": 0} for d in dims}

    for d in dims:
        for N in ns:
            for p_frw in p_values:
                sub = [r for r in patch_rows
                       if r.d == d and r.N == N and abs(r.p_frw - p_frw) < 0.01]
                if len(sub) < 6:
                    continue

                h_arr = np.array([r.H_local for r in sub])
                cells = []
                for feat in ["w_max_ratio", "n_layers", "layer_ratio", "mean_layer_width"]:
                    fa = np.array([getattr(r, feat) for r in sub])
                    mask = ~np.isnan(fa) & ~np.isnan(h_arr)
                    if mask.sum() < 6:
                        cells.append("N/A")
                        continue
                    rho, pv = sp_stats.spearmanr(h_arr[mask], fa[mask])
                    sig = "**" if pv < 0.01 else "*" if pv < 0.05 else ""
                    cells.append(f"{sig}{rho:+.3f}{sig}")

                    if feat == "w_max_ratio":
                        m1_summary[d]["total"] += 1
                        if rho > 0:
                            m1_summary[d]["positive"] += 1
                        if pv < 0.05 and rho > 0:
                            m1_summary[d]["significant"] += 1

                lines.append(f"| {d} | {N} | {p_frw:.2f} | {len(sub)} | " + " | ".join(cells) + " |")

    lines.append("\n**Method 1 — w_max_ratio direction summary:**\n")
    for d in dims:
        s = m1_summary[d]
        lines.append(f"- d={d}: {s['positive']}/{s['total']} positive direction, "
                     f"{s['significant']}/{s['total']} significantly positive (p<0.05)")

    # M1-B: Paired early vs late comparison
    lines.append("\n### M1-B: Paired Early vs Late (same d, N, p, rep)\n")
    lines.append("| d | N | p | n_pairs | mean Δw (early−late) | fraction early>late | p-value (sign test) |")
    lines.append("|---|---|---|---------|---------------------|--------------------|--------------------|")

    m1b_total_pos = 0
    m1b_total_n = 0

    for d in dims:
        for N in ns:
            for p_frw in p_values:
                # Find paired early (smallest t_center) vs late (largest t_center)
                sub = [r for r in patch_rows
                       if r.d == d and r.N == N and abs(r.p_frw - p_frw) < 0.01]
                if len(sub) < 4:
                    continue

                t_min_c = min(r.t_center for r in sub)
                t_max_c = max(r.t_center for r in sub)

                early_reps = {r.rep: r for r in sub if abs(r.t_center - t_min_c) < 0.01}
                late_reps = {r.rep: r for r in sub if abs(r.t_center - t_max_c) < 0.01}
                common = sorted(set(early_reps) & set(late_reps))

                if len(common) < 3:
                    continue

                dw_list = []
                for rep in common:
                    ew = early_reps[rep].w_max_ratio
                    lw = late_reps[rep].w_max_ratio
                    if not (np.isnan(ew) or np.isnan(lw)):
                        dw_list.append(ew - lw)

                if len(dw_list) < 3:
                    continue

                n_pos = sum(1 for dw in dw_list if dw > 0)
                n_pairs = len(dw_list)
                dw_mean = np.mean(dw_list)

                # Sign test (binomial)
                pv_sign = sp_stats.binomtest(n_pos, n_pairs, 0.5).pvalue

                m1b_total_pos += n_pos
                m1b_total_n += n_pairs

                sig = "**" if pv_sign < 0.01 else "*" if pv_sign < 0.05 else ""
                lines.append(
                    f"| {d} | {N} | {p_frw:.2f} | {n_pairs} | {dw_mean:+.4f} | "
                    f"{n_pos}/{n_pairs} ({100*n_pos/n_pairs:.0f}%) | {sig}{pv_sign:.3f}{sig} |"
                )

    if m1b_total_n > 0:
        lines.append(f"\n**Overall paired direction: {m1b_total_pos}/{m1b_total_n} "
                     f"({100*m1b_total_pos/m1b_total_n:.1f}%) early > late**\n")

    # M1-C: Density-residualized patch comparison
    lines.append("### M1-C: Density-Residualized Patch Features vs H_local\n")
    lines.append("OLS-remove n_causal_pairs from each feature, then Spearman vs H_local.\n")
    lines.append("| d | N | Feature | ρ_raw | ρ_resid | Beyond density? |")
    lines.append("|---|---|---------|-------|---------|-----------------|")

    m1c_beyond = {d: [0, 0] for d in dims}

    for d in dims:
        for N in ns:
            sub = [r for r in patch_rows if r.d == d and r.N == N]
            if len(sub) < 10:
                continue
            h_arr = np.array([r.H_local for r in sub])
            dens_arr = np.array([r.n_causal_pairs for r in sub], dtype=float)

            for feat in ac_feats:
                fa = np.array([getattr(r, feat) for r in sub])
                mask = ~np.isnan(fa) & ~np.isnan(h_arr)
                if mask.sum() < 10:
                    continue

                rho_raw, _ = sp_stats.spearmanr(h_arr[mask], fa[mask])
                resid = ols_residualize(fa[mask], dens_arr[mask])
                rho_res, pv_res = sp_stats.spearmanr(h_arr[mask], resid)
                beyond = abs(rho_res) > 0.3 and pv_res < 0.05
                m1c_beyond[d][1] += 1
                if beyond:
                    m1c_beyond[d][0] += 1
                marker = "✅" if beyond else "—"
                lines.append(
                    f"| {d} | {N} | {feat} | {rho_raw:+.3f} | {rho_res:+.3f} | {marker} |"
                )

    lines.append("\n**Method 1 beyond-density summary:**\n")
    for d in dims:
        p_, t_ = m1c_beyond[d]
        lines.append(f"- d={d}: {p_}/{t_} features pass beyond-density")

    # ==================================================================
    # METHOD 2: Density-Matched Sub-sampling
    # ==================================================================
    lines.append("\n---\n")
    lines.append("## Method 2: Density-Matched Sub-sampling\n")
    lines.append("Full [0.1, 1.0] sprinkling, split early/late by median time,")
    lines.append("downsample larger bin to match smaller bin's N. K=5 sub-samples averaged.\n")

    # M2-A: Paired early vs late
    lines.append("### M2-A: Paired Early vs Late (matched N)\n")
    lines.append("| d | N | p | n_matched | early w_max | late w_max | Δw (E−L) | direction |")
    lines.append("|---|---|---|-----------|------------|-----------|----------|-----------|")

    m2_pos = 0
    m2_total = 0

    for d in dims:
        for N in ns:
            for p_frw in p_values:
                early_list = [r for r in matched_rows
                              if r.d == d and r.N == N
                              and abs(r.p_frw - p_frw) < 0.01
                              and r.bin_label == "early"]
                late_list = [r for r in matched_rows
                             if r.d == d and r.N == N
                             and abs(r.p_frw - p_frw) < 0.01
                             and r.bin_label == "late"]

                if not early_list or not late_list:
                    continue

                e_by_rep = {r.rep: r for r in early_list}
                l_by_rep = {r.rep: r for r in late_list}
                common = sorted(set(e_by_rep) & set(l_by_rep))

                if not common:
                    continue

                dw_list = []
                for rep in common:
                    ew = e_by_rep[rep].w_max_ratio
                    lw = l_by_rep[rep].w_max_ratio
                    if not (np.isnan(ew) or np.isnan(lw)):
                        dw_list.append(ew - lw)

                if not dw_list:
                    continue

                dw_mean = np.mean(dw_list)
                n_pos = sum(1 for dw in dw_list if dw > 0)
                n_pairs = len(dw_list)
                m2_total += n_pairs
                m2_pos += n_pos

                mean_e_w = np.mean([e_by_rep[rep].w_max_ratio for rep in common])
                mean_l_w = np.mean([l_by_rep[rep].w_max_ratio for rep in common])
                n_matched = early_list[0].n_matched if early_list else 0

                direction = "early>late ✅" if dw_mean > 0 else "late>early ❌"
                lines.append(
                    f"| {d} | {N} | {p_frw:.2f} | {n_matched} | "
                    f"{mean_e_w:.4f} | {mean_l_w:.4f} | {dw_mean:+.4f} ({n_pos}/{n_pairs}) | {direction} |"
                )

    if m2_total > 0:
        lines.append(f"\n**Method 2 overall: {m2_pos}/{m2_total} "
                     f"({100*m2_pos/m2_total:.1f}%) early > late**\n")

    # M2-B: Density-residualized comparison
    lines.append("### M2-B: Density-Residualized Matched Features vs H_local\n")
    lines.append("| d | Feature | ρ_raw(H_local, feat) | ρ_resid | Beyond density? |")
    lines.append("|---|---------|---------------------|---------|-----------------|")

    m2b_beyond = {d: [0, 0] for d in dims}
    local_feats = ["w_max_ratio", "n_layers", "layer_ratio", "mean_layer_width"]

    for d in dims:
        sub = [r for r in matched_rows if r.d == d]
        if len(sub) < 10:
            continue
        h_arr = np.array([r.H_local for r in sub])
        dens_arr = np.array([r.n_causal_pairs for r in sub], dtype=float)

        for feat in local_feats:
            fa = np.array([getattr(r, feat) for r in sub])
            mask = ~np.isnan(fa) & ~np.isnan(h_arr)
            if mask.sum() < 10:
                continue

            rho_raw, _ = sp_stats.spearmanr(h_arr[mask], fa[mask])
            resid = ols_residualize(fa[mask], dens_arr[mask])
            rho_res, pv_res = sp_stats.spearmanr(h_arr[mask], resid)
            beyond = abs(rho_res) > 0.3 and pv_res < 0.05
            m2b_beyond[d][1] += 1
            if beyond:
                m2b_beyond[d][0] += 1
            marker = "✅" if beyond else "—"
            lines.append(
                f"| {d} | {feat} | {rho_raw:+.3f} | {rho_res:+.3f} | {marker} |"
            )

    lines.append("\n**Method 2 beyond-density summary:**\n")
    for d in dims:
        p_, t_ = m2b_beyond[d]
        lines.append(f"- d={d}: {p_}/{t_} features pass beyond-density")

    # ==================================================================
    # VERDICT
    # ==================================================================
    lines.append("\n---\n")
    lines.append("## Verdict\n")

    # Method 1 summary
    m1_total_pos = sum(m1_summary[d]["positive"] for d in dims)
    m1_total_sig = sum(m1_summary[d]["significant"] for d in dims)
    m1_total_n = sum(m1_summary[d]["total"] for d in dims)
    lines.append(f"### Method 1 (Independent Patches): "
                 f"{m1_total_pos}/{m1_total_n} positive direction, "
                 f"{m1_total_sig}/{m1_total_n} significant\n")

    if m1_total_n > 0:
        pos_frac = m1_total_pos / m1_total_n
        if pos_frac > 0.7:
            lines.append("**Strong support for local H(t) tracking.** The majority of (d,N,p) cells")
            lines.append("show the expected positive correlation between H_local and w_max_ratio")
            lines.append("in independently sprinkled patches — no density confound possible.\n")
        elif pos_frac > 0.5:
            lines.append("**Moderate support for local H(t) tracking.** More than half of cells")
            lines.append("show the expected direction, but the signal is not overwhelming.\n")
        else:
            lines.append("**Weak or no support for local H(t) tracking.** The expected direction")
            lines.append("does not dominate, suggesting the observable may not respond to local")
            lines.append("curvature in power-law FRW (at least at these N values).\n")

    # Method 2 summary
    if m2_total > 0:
        m2_frac = m2_pos / m2_total
        lines.append(f"### Method 2 (Density-Matched): "
                     f"{m2_pos}/{m2_total} ({100*m2_frac:.1f}%) early > late\n")

        if m2_frac > 0.6:
            lines.append("**Density-matched sub-sampling REVERSES the §4.1.34 confound.**")
            lines.append("After equalizing element counts, the early (high-H) bin shows wider")
            lines.append("antichains, confirming that the §4.1.34 null was indeed methodological.\n")
        elif m2_frac > 0.4:
            lines.append("**Density-matched result is ambiguous.** Approximately 50/50, suggesting")
            lines.append("the §4.1.34 confound was real but the underlying signal may be weak.\n")
        else:
            lines.append("**Even with density matching, late > early persists.** This would suggest")
            lines.append("the §4.1.34 result was not purely a density artifact.\n")

    # Overall
    lines.append("### Overall Assessment\n")

    # Combine evidence
    m1c_total_beyond = sum(m1c_beyond[d][0] for d in dims)
    m1c_total_tested = sum(m1c_beyond[d][1] for d in dims)
    m2b_total_beyond = sum(m2b_beyond[d][0] for d in dims)
    m2b_total_tested = sum(m2b_beyond[d][1] for d in dims)

    lines.append(f"- Method 1 beyond-density: {m1c_total_beyond}/{m1c_total_tested}")
    lines.append(f"- Method 2 beyond-density: {m2b_total_beyond}/{m2b_total_tested}")
    if m1b_total_n > 0:
        lines.append(f"- Method 1 paired direction: {m1b_total_pos}/{m1b_total_n} "
                     f"({100*m1b_total_pos/m1b_total_n:.1f}%)")
    if m2_total > 0:
        lines.append(f"- Method 2 paired direction: {m2_pos}/{m2_total} "
                     f"({100*m2_pos/m2_total:.1f}%)")

    lines.append("")
    lines.append("### Comparison with §4.1.34 Phase B\n")
    lines.append("| Metric | §4.1.34 Phase B | §4.1.35 Method 1 | §4.1.35 Method 2 |")
    lines.append("|--------|----------------|------------------|------------------|")
    lines.append(f"| N balance | Unequal (confounded) | Equal (independent) | Matched (subsampled) |")
    if m1b_total_n > 0 and m2_total > 0:
        lines.append(f"| early>late fraction | 0/54 (0%) | {m1b_total_pos}/{m1b_total_n} "
                     f"({100*m1b_total_pos/m1b_total_n:.1f}%) | {m2_pos}/{m2_total} "
                     f"({100*m2_pos/m2_total:.1f}%) |")
    lines.append(f"| Beyond density | N/A (confounded) | {m1c_total_beyond}/{m1c_total_tested} | "
                 f"{m2b_total_beyond}/{m2b_total_tested} |")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="§4.1.35: Density-Matched Local H(t) Tracking"
    )
    ap.add_argument("--dims", nargs="*", type=int, default=[2, 3, 4])
    ap.add_argument("--ns", nargs="*", type=int, default=[128, 256])
    ap.add_argument("--p-values", nargs="*", type=float,
                    default=[0.5, 1.0, 1.5, 2.0])
    ap.add_argument("--t-centers", nargs="*", type=float,
                    default=[0.2, 0.4, 0.6, 0.8])
    ap.add_argument("--t-half-width", type=float, default=0.08)
    ap.add_argument("--reps", type=int, default=8)
    ap.add_argument("--seed", type=int, default=2035)
    ap.add_argument("--out-patch-csv",
                    default="outputs_unified_functional/conjecture_e_local_h_patches.csv")
    ap.add_argument("--out-matched-csv",
                    default="outputs_unified_functional/conjecture_e_local_h_matched.csv")
    ap.add_argument("--out-report",
                    default="outputs_unified_functional/conjecture_e_local_h_tracking.md")
    args = ap.parse_args()

    dims = args.dims
    ns = args.ns
    p_values = args.p_values
    t_centers = args.t_centers
    t_hw = args.t_half_width
    reps = args.reps
    seed = args.seed

    Path(args.out_patch_csv).parent.mkdir(parents=True, exist_ok=True)

    # =============================================
    # Method 1: Independent Patch Sprinklings
    # =============================================
    print("=" * 60)
    print("Method 1: Independent Patch Sprinklings")
    print("=" * 60)

    patch_rows: list[PatchRow] = []
    total_runs = len(dims) * len(ns) * len(p_values) * len(t_centers) * reps
    count = 0
    for d in dims:
        for N in ns:
            for p_frw in p_values:
                for t_c in t_centers:
                    for rep in range(reps):
                        count += 1
                        if count % 30 == 0 or count == total_runs:
                            print(f"  [{count}/{total_runs}] d={d} N={N} p={p_frw:.2f} "
                                  f"t_c={t_c:.2f} rep={rep}")
                        try:
                            row = run_patch(d, N, p_frw, rep, t_c, t_hw, seed)
                            patch_rows.append(row)
                        except Exception as e:
                            print(f"  WARNING: d={d} N={N} p={p_frw} t_c={t_c} "
                                  f"rep={rep} failed: {e}")

    # Save patch CSV
    with open(args.out_patch_csv, "w", newline="", encoding="utf-8") as f:
        fieldnames = [fld.name for fld in fields(PatchRow)]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in patch_rows:
            w.writerow({fld.name: getattr(row, fld.name) for fld in fields(PatchRow)})
    print(f"\nPatch CSV saved: {args.out_patch_csv}")
    print(f"  Total rows: {len(patch_rows)}")

    # =============================================
    # Method 2: Density-Matched Sub-sampling
    # =============================================
    print("\n" + "=" * 60)
    print("Method 2: Density-Matched Sub-sampling")
    print("=" * 60)

    matched_rows: list[MatchedRow] = []
    total_runs2 = len(dims) * len(ns) * len(p_values) * reps
    count = 0
    for d in dims:
        for N in ns:
            for p_frw in p_values:
                for rep in range(reps):
                    count += 1
                    if count % 20 == 0 or count == total_runs2:
                        print(f"  [{count}/{total_runs2}] d={d} N={N} p={p_frw:.2f} rep={rep}")
                    try:
                        rows = run_matched_bins(d, N, p_frw, rep, seed)
                        matched_rows.extend(rows)
                    except Exception as e:
                        print(f"  WARNING: d={d} N={N} p={p_frw} rep={rep} failed: {e}")

    # Save matched CSV
    with open(args.out_matched_csv, "w", newline="", encoding="utf-8") as f:
        fieldnames = [fld.name for fld in fields(MatchedRow)]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in matched_rows:
            w.writerow({fld.name: getattr(row, fld.name) for fld in fields(MatchedRow)})
    print(f"\nMatched CSV saved: {args.out_matched_csv}")
    print(f"  Total rows: {len(matched_rows)}")

    # =============================================
    # Generate report
    # =============================================
    print("\nGenerating report...")
    report = generate_report(patch_rows, matched_rows, dims, ns, p_values, t_centers)
    with open(args.out_report, "w", encoding="utf-8") as f:
        f.write(report)
    print(f"Report saved: {args.out_report}")

    # Quick summary
    print("\n" + "=" * 60)
    print("Quick Summary")
    print("=" * 60)
    print(f"Method 1: {len(patch_rows)} patch rows "
          f"(d={dims}, N={ns}, p={p_values}, t_c={t_centers})")
    print(f"Method 2: {len(matched_rows)} matched rows")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
