"""Conjecture E — Shape-Ratio Family Experiment (§4.1.24).

Following the Two-Family Decomposition (§4.1.23), this experiment focuses
EXCLUSIVELY on C_k ratio/shape observables — metrics that capture how the
interval size distribution RESHAPES under curvature, rather than how much
total connectivity there is.

Design principle: all candidates must be RATIO or SHAPE metrics that are
insensitive to overall density scaling. If C_k all scale by the same factor,
a good shape metric should not change.

Candidate shape-ratio metrics:
  S1: C1/C0                      — first-order ratio (bdg_d2c's core signal)
  S2: C2/C1                      — second-order ratio
  S3: C2/C0                      — two-step ratio
  S4: (C1+C2)/(C0+C1+C2)        — non-link fraction among low-order
  S5: C1²/(C0·C2)               — curvature of log-distribution at k=1
  S6: log(C1/C0) - log(C2/C1)   — log-concavity (second difference of log C_k)
  S7: Σ(k·C_k)/Σ(C_k) / (C0/total)  — mean_interval / link_fraction
  S8: (N·C1 - 2·C0·C1) / (C0² + C1²)  — bdg_d2c rewritten as shape
  S9: gini(C_k distribution)     — Gini coefficient of {C_k}
  S10: KL(C_k || C_k^flat)       — KL divergence from flat baseline
  S11: C0/total (link fraction)  — reference (known ↑ with H)
  S12: bdg_d2c / N               — intensive bdg_d2c (reference)

Also compute flat-baseline ratios by averaging H=0 runs per (d,N) slice
and computing deviation metrics relative to flat.

Design: d=2/3/4, N=128/256/512, H=0–2, 8 reps = 360 realizations
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, fields
from pathlib import Path
from collections import defaultdict

import numpy as np
from scipy import stats as sp_stats

from curvature_layer2_de_sitter_scan import sprinkle_de_sitter_like_diamond
from curvature_layer2_recovery_fast import (
    build_causal_matrix_de_sitter,
    compute_tau_and_intervals,
)
from curvature_layer2_baseline import fit_curvature_proxy
from bd_action import IntervalCounts, bd_ratio_metric, bdg_action_d2_corrected


def interval_counts_from_causal(causal: np.ndarray) -> IntervalCounts:
    """Compute IntervalCounts from boolean causal matrix."""
    c_int = causal.astype(np.int32)
    sizes = c_int @ c_int
    ks = sizes[causal].astype(np.int64)
    total = int(ks.size)
    if total == 0:
        return IntervalCounts(counts={}, total_relations=0)
    bc = np.bincount(ks)
    nonzero = np.nonzero(bc)[0]
    out = {int(k): int(bc[k]) for k in nonzero}
    return IntervalCounts(counts=out, total_relations=total)


# ── Shape-ratio metrics ──

def safe_ratio(a: float, b: float) -> float:
    return a / b if b > 0 else float("nan")


def safe_log(x: float) -> float:
    return math.log(x) if x > 0 else float("nan")


def gini_coefficient(counts_dict: dict[int, int], total: int) -> float:
    """Gini coefficient of the C_k distribution."""
    if total <= 0:
        return float("nan")
    vals = sorted(counts_dict.values())
    n = len(vals)
    if n <= 1:
        return 0.0
    cum = np.cumsum(vals)
    s = float(cum[-1])
    if s <= 0:
        return float("nan")
    # Gini = 1 - 2 * area under Lorenz curve
    lorenz_sum = sum((2 * (i + 1) - n - 1) * v for i, v in enumerate(vals))
    return float(lorenz_sum) / (n * s)


def compute_shape_metrics(counts: IntervalCounts, N: int) -> dict[str, float]:
    """Compute all shape-ratio metrics from interval counts."""
    total = float(counts.total_relations) if counts.total_relations > 0 else 0.0
    c0 = float(counts.get(0))
    c1 = float(counts.get(1))
    c2 = float(counts.get(2))
    c3 = float(counts.get(3))

    # Mean interval size (for S7)
    mean_k = sum(float(k) * float(v) for k, v in counts.counts.items()) / total if total > 0 else 0.0
    link_frac = c0 / total if total > 0 else 1.0

    metrics: dict[str, float] = {}

    # S1: C1/C0
    metrics["S1_c1_over_c0"] = safe_ratio(c1, c0)

    # S2: C2/C1
    metrics["S2_c2_over_c1"] = safe_ratio(c2, c1)

    # S3: C2/C0
    metrics["S3_c2_over_c0"] = safe_ratio(c2, c0)

    # S4: (C1+C2)/(C0+C1+C2) — non-link fraction among first 3 orders
    low_total = c0 + c1 + c2
    metrics["S4_nonlink_low"] = safe_ratio(c1 + c2, low_total)

    # S5: C1²/(C0·C2) — curvature of log-distribution at k=1
    metrics["S5_log_curvature"] = safe_ratio(c1 * c1, c0 * c2)

    # S6: log(C1/C0) - log(C2/C1) = 2·log(C1) - log(C0) - log(C2)
    lc0, lc1, lc2 = safe_log(c0), safe_log(c1), safe_log(c2)
    if not (math.isnan(lc0) or math.isnan(lc1) or math.isnan(lc2)):
        metrics["S6_log_concavity"] = 2 * lc1 - lc0 - lc2
    else:
        metrics["S6_log_concavity"] = float("nan")

    # S7: mean_interval / link_fraction — shape-adjusted depth
    metrics["S7_depth_over_linkfrac"] = safe_ratio(mean_k, link_frac) if link_frac > 0 else float("nan")

    # S8: bdg_d2c rewritten: (N - 2C0 + 2C1) / N = 1 - 2(C0-C1)/N
    metrics["S8_bdg_d2c_shape"] = 1.0 - 2.0 * (c0 - c1) / N if N > 0 else float("nan")

    # S9: Gini coefficient of {C_k}
    metrics["S9_gini"] = gini_coefficient(counts.counts, counts.total_relations)

    # S10: placeholder — will compute KL vs flat baseline in post-processing
    # For now store the raw distribution fingerprint
    metrics["S10_kl_vs_flat"] = float("nan")  # computed in post-processing

    # S11: link fraction (reference, known ↑ with H)
    metrics["S11_link_frac"] = link_frac

    # S12: bdg_d2c / N (intensive reference)
    d2c = bdg_action_d2_corrected(counts, N, normalized=True)
    metrics["S12_bdg_d2c_intensive"] = d2c / N if N > 0 else float("nan")

    # Bonus: consecutive ratio chain
    # S13: C3/C2
    metrics["S13_c3_over_c2"] = safe_ratio(c3, c2)

    # S14: geometric mean ratio = (C1/C0 · C2/C1 · C3/C2)^(1/3) = (C3/C0)^(1/3)
    metrics["S14_geo_mean_ratio"] = safe_ratio(c3, c0) ** (1.0 / 3.0) if c0 > 0 and c3 > 0 else float("nan")

    # S15: ratio of ratios = (C2/C1) / (C1/C0) = C0·C2/C1²
    # This is 1/S5 — measures how fast the ratio changes
    metrics["S15_ratio_of_ratios"] = safe_ratio(c0 * c2, c1 * c1)

    return metrics


@dataclass(frozen=True)
class Row:
    d: int
    N: int
    hubble: float
    rep: int
    R_dS: float
    R_hat: float
    n_causal_pairs: int
    bd_ratio: float
    occupancy_R: float
    # Shape metrics
    S1_c1_over_c0: float
    S2_c2_over_c1: float
    S3_c2_over_c0: float
    S4_nonlink_low: float
    S5_log_curvature: float
    S6_log_concavity: float
    S7_depth_over_linkfrac: float
    S8_bdg_d2c_shape: float
    S9_gini: float
    S10_kl_vs_flat: float
    S11_link_frac: float
    S12_bdg_d2c_intensive: float
    S13_c3_over_c2: float
    S14_geo_mean_ratio: float
    S15_ratio_of_ratios: float
    # Raw counts for post-processing
    C0: int
    C1: int
    C2: int
    C3: int
    C4: int
    total_pairs: int


def run_single(
    d: int, N: int, hubble: float, rep: int, base_seed: int,
    tau_quantile: float = 0.75,
) -> Row:
    seed = base_seed + d * 100000 + N * 1000 + rep * 10 + int(hubble * 100)

    points = sprinkle_de_sitter_like_diamond(N, d - 1, hubble=hubble, seed=seed)
    causal = build_causal_matrix_de_sitter(points, hubble)

    tau, k = compute_tau_and_intervals(points, hubble, causal)
    n_pairs = int(tau.size)

    nan = float("nan")
    if n_pairs < 10:
        R_dS = d * (d - 1) * hubble ** 2
        return Row(
            d=d, N=N, hubble=hubble, rep=rep, R_dS=R_dS, R_hat=nan,
            n_causal_pairs=n_pairs, bd_ratio=nan, occupancy_R=nan,
            S1_c1_over_c0=nan, S2_c2_over_c1=nan, S3_c2_over_c0=nan,
            S4_nonlink_low=nan, S5_log_curvature=nan, S6_log_concavity=nan,
            S7_depth_over_linkfrac=nan, S8_bdg_d2c_shape=nan, S9_gini=nan,
            S10_kl_vs_flat=nan, S11_link_frac=nan, S12_bdg_d2c_intensive=nan,
            S13_c3_over_c2=nan, S14_geo_mean_ratio=nan, S15_ratio_of_ratios=nan,
            C0=0, C1=0, C2=0, C3=0, C4=0, total_pairs=0,
        )

    y = k / max(N - 2, 1)
    x = np.power(tau, d)
    ratio_arr = y / np.clip(x, 1e-12, None)
    fit_mask = tau <= np.quantile(tau, tau_quantile)
    if fit_mask.sum() < 10:
        fit_mask = np.ones_like(tau, dtype=bool)
    _, r_hat = fit_curvature_proxy(tau[fit_mask], ratio_arr[fit_mask], d)

    counts = interval_counts_from_causal(causal)
    br = bd_ratio_metric(counts)
    total_rel = counts.total_relations
    c0 = counts.get(0)
    f_lk = c0 / total_rel if total_rel > 0 else 1.0
    occ_R = 1.0 - f_lk
    R_dS = d * (d - 1) * hubble ** 2

    sm = compute_shape_metrics(counts, N)

    return Row(
        d=d, N=N, hubble=hubble, rep=rep, R_dS=R_dS, R_hat=r_hat,
        n_causal_pairs=n_pairs, bd_ratio=br, occupancy_R=occ_R,
        S1_c1_over_c0=sm["S1_c1_over_c0"],
        S2_c2_over_c1=sm["S2_c2_over_c1"],
        S3_c2_over_c0=sm["S3_c2_over_c0"],
        S4_nonlink_low=sm["S4_nonlink_low"],
        S5_log_curvature=sm["S5_log_curvature"],
        S6_log_concavity=sm["S6_log_concavity"],
        S7_depth_over_linkfrac=sm["S7_depth_over_linkfrac"],
        S8_bdg_d2c_shape=sm["S8_bdg_d2c_shape"],
        S9_gini=sm["S9_gini"],
        S10_kl_vs_flat=sm["S10_kl_vs_flat"],
        S11_link_frac=sm["S11_link_frac"],
        S12_bdg_d2c_intensive=sm["S12_bdg_d2c_intensive"],
        S13_c3_over_c2=sm["S13_c3_over_c2"],
        S14_geo_mean_ratio=sm["S14_geo_mean_ratio"],
        S15_ratio_of_ratios=sm["S15_ratio_of_ratios"],
        C0=counts.get(0), C1=counts.get(1), C2=counts.get(2),
        C3=counts.get(3), C4=counts.get(4),
        total_pairs=total_rel,
    )


# ── Post-processing: compute KL vs flat baseline ──

def compute_kl_vs_flat(rows: list[Row], dims: list[int], ns: list[int]) -> None:
    """Compute S10 = KL(p_k || p_k^flat) for each row using H=0 average as baseline."""
    flat_dists: dict[tuple[int, int], dict[int, float]] = {}

    for d in dims:
        for N in ns:
            flat_rows = [r for r in rows if r.d == d and r.N == N and r.hubble == 0.0]
            if not flat_rows:
                continue
            # Average the C_k distribution
            max_k = 20
            avg_dist: dict[int, float] = {}
            for kk in range(max_k + 1):
                vals = []
                for r in flat_rows:
                    if kk == 0: vals.append(float(r.C0))
                    elif kk == 1: vals.append(float(r.C1))
                    elif kk == 2: vals.append(float(r.C2))
                    elif kk == 3: vals.append(float(r.C3))
                    elif kk == 4: vals.append(float(r.C4))
                if vals:
                    avg_dist[kk] = float(np.mean(vals))
            flat_dists[(d, N)] = avg_dist

    # Now compute KL for each row
    new_rows = []
    for r in rows:
        flat_d = flat_dists.get((r.d, r.N))
        if flat_d is None or r.total_pairs <= 0:
            new_rows.append(r)
            continue

        # Build distributions
        total_flat = sum(flat_d.values())
        total_cur = float(r.total_pairs)
        if total_flat <= 0 or total_cur <= 0:
            new_rows.append(r)
            continue

        kl = 0.0
        for kk in range(5):  # k = 0..4
            if kk == 0: ck = float(r.C0)
            elif kk == 1: ck = float(r.C1)
            elif kk == 2: ck = float(r.C2)
            elif kk == 3: ck = float(r.C3)
            elif kk == 4: ck = float(r.C4)
            else: ck = 0.0

            p = ck / total_cur
            q = flat_d.get(kk, 0.0) / total_flat

            if p > 1e-10 and q > 1e-10:
                kl += p * math.log(p / q)

        # Replace the row with updated S10
        d_row = {f.name: getattr(r, f.name) for f in fields(r)}
        d_row["S10_kl_vs_flat"] = kl
        new_rows.append(Row(**d_row))

    rows.clear()
    rows.extend(new_rows)


# ── Analysis ──

SHAPE_METRICS = [
    ("S1_c1_over_c0", "S1: C₁/C₀"),
    ("S2_c2_over_c1", "S2: C₂/C₁"),
    ("S3_c2_over_c0", "S3: C₂/C₀"),
    ("S4_nonlink_low", "S4: (C₁+C₂)/(C₀+C₁+C₂)"),
    ("S5_log_curvature", "S5: C₁²/(C₀·C₂)"),
    ("S6_log_concavity", "S6: log-concavity"),
    ("S7_depth_over_linkfrac", "S7: mean_k / f_link"),
    ("S8_bdg_d2c_shape", "S8: 1-2(C₀-C₁)/N"),
    ("S9_gini", "S9: Gini(C_k)"),
    ("S10_kl_vs_flat", "S10: KL(p||p_flat)"),
    ("S11_link_frac", "S11: f_link (ref ↑)"),
    ("S12_bdg_d2c_intensive", "S12: bdg_d2c/N (ref ↑)"),
    ("S13_c3_over_c2", "S13: C₃/C₂"),
    ("S14_geo_mean_ratio", "S14: (C₃/C₀)^(1/3)"),
    ("S15_ratio_of_ratios", "S15: C₀C₂/C₁²"),
]


def analyze(rows: list[Row], dims: list[int], ns: list[int]) -> list[str]:
    lines: list[str] = []
    lines.append("# Conjecture E — Shape-Ratio Family (§4.1.24)\n")
    lines.append(f"Total realizations: {len(rows)}\n")
    lines.append("Focus: ONLY ratio/shape metrics from C_k distribution.\n")
    lines.append("Goal: find metrics that positively track H (curvature) via distribution shape change.\n")

    # ── Monotonicity table ──
    lines.append("\n## Spearman(H², metric) — Monotonicity Tournament\n")

    results: dict[str, list[tuple[int, int, float, float]]] = {}
    for attr, _ in SHAPE_METRICS:
        results[attr] = []

    for d in dims:
        for N in ns:
            for attr, _ in SHAPE_METRICS:
                sub = [r for r in rows if r.d == d and r.N == N
                       and not math.isnan(getattr(r, attr))]
                if len(sub) < 5:
                    results[attr].append((d, N, float("nan"), float("nan")))
                    continue
                h2 = np.array([r.hubble ** 2 for r in sub])
                vals = np.array([getattr(r, attr) for r in sub])
                rho, p = sp_stats.spearmanr(h2, vals)
                results[attr].append((d, N, rho, p))

    # Focus on ↑ metrics (positive correlation with H²)
    lines.append("\n### Positive-direction candidates (want ↑ with H)\n")
    lines.append("| Metric | d | N | ρ(H²) | p | sig? |")
    lines.append("|--------|---|---|--------|---|------|")

    for attr, label in SHAPE_METRICS:
        for d, N, rho, p in results[attr]:
            if math.isnan(rho):
                continue
            if rho > 0.2:  # only show positive candidates
                sig = "✅" if p < 0.05 else "❌"
                lines.append(f"| {label[:22]} | {d} | {N} | **{rho:+.3f}** | {p:.1e} | {sig} |")

    lines.append("\n### Negative-direction (density camp, for reference)\n")
    lines.append("| Metric | d | N | ρ(H²) | p | sig? |")
    lines.append("|--------|---|---|--------|---|------|")

    for attr, label in SHAPE_METRICS:
        for d, N, rho, p in results[attr]:
            if math.isnan(rho):
                continue
            if rho < -0.5 and p < 0.05:
                lines.append(f"| {label[:22]} | {d} | {N} | {rho:+.3f} | {p:.1e} | ✅ |")

    # ── Winner summary ──
    lines.append("\n## Winner Summary (ranked by positive-direction performance)\n")
    lines.append("| Rank | Metric | # ↑sig | # total | mean ρ(↑) | d=2? | d=3? | d=4? |")
    lines.append("|------|--------|--------|---------|-----------|------|------|------|")

    scored = []
    for attr, label in SHAPE_METRICS:
        up_sig = 0
        up_rhos = []
        total_valid = 0
        per_d: dict[int, list[float]] = {2: [], 3: [], 4: []}
        for d, N, rho, p in results[attr]:
            if math.isnan(rho):
                continue
            total_valid += 1
            per_d[d].append(rho)
            if rho > 0 and p < 0.05:
                up_sig += 1
                up_rhos.append(rho)
        mean_up = float(np.mean(up_rhos)) if up_rhos else 0.0

        d2_ok = any(r > 0.3 for r in per_d[2])
        d3_ok = any(r > 0.3 for r in per_d[3])
        d4_ok = any(r > 0.3 for r in per_d[4])

        scored.append((up_sig, mean_up, label, total_valid, d2_ok, d3_ok, d4_ok, attr))

    scored.sort(key=lambda x: (x[0], x[1]), reverse=True)
    for rank, (us, mu, label, tv, d2, d3, d4, attr) in enumerate(scored, 1):
        d2s = "✅" if d2 else "❌"
        d3s = "✅" if d3 else "❌"
        d4s = "✅" if d4 else "❌"
        lines.append(f"| {rank} | {label[:25]} | {us}/{tv} | {tv} | {mu:.3f} | {d2s} | {d3s} | {d4s} |")

    # ── Correlation with R_hat ──
    lines.append("\n## Correlation with R_hat (positive direction only)\n")
    lines.append("| Metric | d | N | ρ(metric, R_hat) | p |")
    lines.append("|--------|---|---|-----------------|---|")

    for attr, label in SHAPE_METRICS:
        for d in dims:
            for N in ns:
                sub = [r for r in rows if r.d == d and r.N == N
                       and not math.isnan(r.R_hat) and not math.isnan(getattr(r, attr))]
                if len(sub) < 5:
                    continue
                rhat = np.array([r.R_hat for r in sub])
                vals = np.array([getattr(r, attr) for r in sub])
                rho, p = sp_stats.spearmanr(vals, rhat)
                if rho > 0.3 and p < 0.05:
                    lines.append(f"| {label[:22]} | {d} | {N} | **{rho:+.3f}** | {p:.1e} |")

    # ── Per-metric detail for top candidates ──
    top_attrs = [attr for _, _, _, _, _, _, _, attr in scored[:5]]
    lines.append("\n## Top 5 Candidates — Full Detail\n")
    for attr in top_attrs:
        label = dict(SHAPE_METRICS).get(attr, attr)
        lines.append(f"\n### {label}\n")
        lines.append("| d | N | ρ(H²) | p | mean(H=0) | mean(H=2) | Δ% |")
        lines.append("|---|---|--------|---|-----------|-----------|-----|")
        for d in dims:
            for N in ns:
                sub = [r for r in rows if r.d == d and r.N == N
                       and not math.isnan(getattr(r, attr))]
                if len(sub) < 5:
                    continue
                h2 = np.array([r.hubble ** 2 for r in sub])
                vals = np.array([getattr(r, attr) for r in sub])
                rho, p = sp_stats.spearmanr(h2, vals)

                flat_vals = [getattr(r, attr) for r in sub if r.hubble == 0.0]
                high_vals = [getattr(r, attr) for r in sub if r.hubble == 2.0]
                m_flat = float(np.mean(flat_vals)) if flat_vals else float("nan")
                m_high = float(np.mean(high_vals)) if high_vals else float("nan")
                if m_flat != 0 and not math.isnan(m_flat) and not math.isnan(m_high):
                    delta_pct = 100.0 * (m_high - m_flat) / abs(m_flat)
                else:
                    delta_pct = float("nan")

                lines.append(
                    f"| {d} | {N} | {rho:+.3f} | {p:.1e} | {m_flat:.4f} | {m_high:.4f} | {delta_pct:+.1f}% |"
                )

    return lines


def main() -> int:
    ap = argparse.ArgumentParser(description="Shape-ratio family experiment")
    ap.add_argument("--dims", nargs="*", type=int, default=[2, 3, 4])
    ap.add_argument("--ns", nargs="*", type=int, default=[128, 256, 512])
    ap.add_argument("--hubbles", nargs="*", type=float, default=[0.0, 0.25, 0.5, 1.0, 2.0])
    ap.add_argument("--reps", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20260324)
    ap.add_argument("--out", default="outputs_unified_functional/conjecture_e_shape_ratio_family.csv")
    ap.add_argument("--report", default="outputs_unified_functional/conjecture_e_shape_ratio_family.md")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[Row] = []
    total = len(args.dims) * len(args.ns) * len(args.hubbles) * args.reps
    done = 0

    for d in args.dims:
        for N in args.ns:
            for H in args.hubbles:
                for rep in range(args.reps):
                    row = run_single(d, N, H, rep, args.seed)
                    rows.append(row)
                    done += 1
                    if done % 30 == 0 or done == total:
                        print(f"  [{done}/{total}] d={d} N={N} H={H:.2f} "
                              f"S1={row.S1_c1_over_c0:.4f} "
                              f"S5={row.S5_log_curvature:.4f} "
                              f"S10={row.S10_kl_vs_flat:.4f}")

    # Post-process: compute KL vs flat
    print("Computing KL divergence vs flat baseline...")
    compute_kl_vs_flat(rows, args.dims, args.ns)

    # Save CSV
    fn = [f.name for f in fields(Row)]
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fn)
        w.writeheader()
        for r in rows:
            d_row = {}
            for fk in fn:
                v = getattr(r, fk)
                d_row[fk] = f"{v:.6f}" if isinstance(v, float) else str(v)
            w.writerow(d_row)
    print(f"\nSaved: {out_path}")

    # Analysis
    report_lines = analyze(rows, args.dims, args.ns)
    report_text = "\n".join(report_lines) + "\n"
    report_path = Path(args.report)
    report_path.write_text(report_text, encoding="utf-8")
    print(f"Saved: {report_path}")

    # Console summary
    print("\n" + "=" * 60)
    for line in report_lines[-30:]:
        print(line)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
