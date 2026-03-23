"""Conjecture E — Curvature-Sensitive BD Metrics.

The standard BDG d=4 action has almost no signal on de Sitter sprinklings
(§4.1.22: 9/9 slices not significant). This script designs and tests
ALTERNATIVE coordinate-free BD-derived quantities that respond monotonically
to curvature.

Key insight: curvature changes the SHAPE of the interval size distribution
{C_k}. In flat space, most causal pairs are links (C_0 dominates). As
curvature increases (de Sitter H↑), causal pairs become sparser and the
few surviving non-link pairs have smaller intervals. This shifts the
distribution toward C_0.

Candidate metrics (all coordinate-free):
  M1: mean_interval_size = Σ(k·C_k) / Σ(C_k)  — average interval depth
  M2: interval_entropy = -Σ p_k log p_k where p_k = C_k/total
  M3: C1_ratio = C_1 / (C_0 + C_1)  — fraction of order-1 intervals
  M4: tail_weight = Σ_{k≥2} C_k / total  — fraction of deep intervals
  M5: log_interval_ratio = log(C_1/C_0) if C_1>0  — log ratio of first two
  M6: bdg_d2c_shifted = bdg_d2c - bdg_d2c(flat baseline)  — curvature shift
  M7: spectral_gap = (C_1 - C_0·expected_ratio) / C_0  — deviation from flat

Design: same grid as §4.1.22 (d=2/3/4, N=128/256/512, H=0–2, 8 reps)
but now compute all 7 candidate metrics and test monotonicity with H².
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

from curvature_layer2_de_sitter_scan import (
    sprinkle_de_sitter_like_diamond,
)
from curvature_layer2_recovery_fast import (
    build_causal_matrix_de_sitter,
    compute_tau_and_intervals,
)
from curvature_layer2_baseline import fit_curvature_proxy, fit_flat_volume_law
from bd_action import (
    IntervalCounts,
    bd_ratio_metric,
    bdg_action_d2_corrected,
    bdg_action_d4_standard,
)


def interval_counts_from_causal(causal: np.ndarray) -> IntervalCounts:
    """Compute IntervalCounts directly from a boolean causal matrix."""
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


# ── Candidate metrics ──

def mean_interval_size(counts: IntervalCounts) -> float:
    """M1: average interval depth."""
    if counts.total_relations <= 0:
        return 0.0
    total = float(counts.total_relations)
    s = sum(float(k) * float(v) for k, v in counts.counts.items())
    return s / total


def interval_entropy(counts: IntervalCounts) -> float:
    """M2: Shannon entropy of the interval size distribution."""
    if counts.total_relations <= 0:
        return 0.0
    total = float(counts.total_relations)
    ent = 0.0
    for v in counts.counts.values():
        p = float(v) / total
        if p > 0:
            ent -= p * math.log(p)
    return ent


def c1_ratio(counts: IntervalCounts) -> float:
    """M3: C_1 / (C_0 + C_1)."""
    c0 = float(counts.get(0))
    c1 = float(counts.get(1))
    denom = c0 + c1
    return c1 / denom if denom > 0 else 0.0


def tail_weight(counts: IntervalCounts) -> float:
    """M4: fraction of pairs with k >= 2."""
    if counts.total_relations <= 0:
        return 0.0
    total = float(counts.total_relations)
    c0 = float(counts.get(0))
    c1 = float(counts.get(1))
    return 1.0 - (c0 + c1) / total


def log_interval_ratio(counts: IntervalCounts) -> float:
    """M5: log(C_1 / C_0), or nan if C_1 = 0."""
    c0 = float(counts.get(0))
    c1 = float(counts.get(1))
    if c0 <= 0 or c1 <= 0:
        return float("nan")
    return math.log(c1 / c0)


def weighted_depth_ratio(counts: IntervalCounts) -> float:
    """M6: Σ(k² · C_k) / Σ(k · C_k) — depth-weighted mean interval size.
    Higher values mean longer tails in the distribution. More sensitive
    to curvature than simple mean because k² weights emphasize deep intervals.
    """
    if counts.total_relations <= 0:
        return 0.0
    num = sum(float(k * k) * float(v) for k, v in counts.counts.items())
    denom = sum(float(k) * float(v) for k, v in counts.counts.items())
    return num / denom if denom > 0 else 0.0


def normalized_bd_shift(counts: IntervalCounts, n: int) -> float:
    """M7: (bdg_d2c - expected_flat) / N — shift from flat baseline.
    In flat space bdg_d2c → a known value; curvature shifts it.
    We normalize by N to make it intensive.
    """
    d2c = bdg_action_d2_corrected(counts, n, normalized=True)
    return d2c


# ── Row dataclass ──

@dataclass(frozen=True)
class Row:
    d: int
    N: int
    hubble: float
    rep: int
    R_dS: float
    R_hat: float
    n_causal_pairs: int
    # standard metrics for reference
    bd_ratio: float
    bdg_d4: float
    bdg_d2c: float
    occupancy_R: float
    # NEW candidate metrics
    M1_mean_interval: float
    M2_interval_entropy: float
    M3_c1_ratio: float
    M4_tail_weight: float
    M5_log_interval_ratio: float
    M6_weighted_depth_ratio: float
    M7_norm_bd_shift: float
    # distribution summary
    C0: int
    C1: int
    C2: int
    C3: int
    total_pairs: int


def run_single(
    d: int, N: int, hubble: float, rep: int, base_seed: int,
    tau_quantile: float = 0.75,
) -> Row:
    seed = base_seed + d * 100000 + N * 1000 + rep * 10 + int(hubble * 100)

    points = sprinkle_de_sitter_like_diamond(
        N, d - 1, hubble=hubble, seed=seed,
    )
    causal = build_causal_matrix_de_sitter(points, hubble)

    # Layer 2: curvature proxy
    tau, k = compute_tau_and_intervals(points, hubble, causal)
    n_pairs = int(tau.size)

    if n_pairs < 10:
        R_dS = d * (d - 1) * hubble ** 2
        nan = float("nan")
        return Row(
            d=d, N=N, hubble=hubble, rep=rep, R_dS=R_dS, R_hat=nan,
            n_causal_pairs=n_pairs,
            bd_ratio=nan, bdg_d4=nan, bdg_d2c=nan, occupancy_R=nan,
            M1_mean_interval=nan, M2_interval_entropy=nan,
            M3_c1_ratio=nan, M4_tail_weight=nan,
            M5_log_interval_ratio=nan, M6_weighted_depth_ratio=nan,
            M7_norm_bd_shift=nan,
            C0=0, C1=0, C2=0, C3=0, total_pairs=0,
        )

    y = k / max(N - 2, 1)
    x = np.power(tau, d)
    ratio_arr = y / np.clip(x, 1e-12, None)
    fit_mask = tau <= np.quantile(tau, tau_quantile)
    if fit_mask.sum() < 10:
        fit_mask = np.ones_like(tau, dtype=bool)
    _, r_hat = fit_curvature_proxy(tau[fit_mask], ratio_arr[fit_mask], d)

    # Layer 1: BD actions + new metrics
    counts = interval_counts_from_causal(causal)
    br = bd_ratio_metric(counts)
    d4 = bdg_action_d4_standard(counts, N, normalized=True)
    d2c = bdg_action_d2_corrected(counts, N, normalized=True)

    total_rel = counts.total_relations
    c0 = counts.get(0)
    f_lk = c0 / total_rel if total_rel > 0 else 1.0
    occ_R = 1.0 - f_lk

    R_dS = d * (d - 1) * hubble ** 2

    # New candidate metrics
    m1 = mean_interval_size(counts)
    m2 = interval_entropy(counts)
    m3 = c1_ratio(counts)
    m4 = tail_weight(counts)
    m5 = log_interval_ratio(counts)
    m6 = weighted_depth_ratio(counts)
    m7 = normalized_bd_shift(counts, N)

    return Row(
        d=d, N=N, hubble=hubble, rep=rep, R_dS=R_dS, R_hat=r_hat,
        n_causal_pairs=n_pairs,
        bd_ratio=br, bdg_d4=d4, bdg_d2c=d2c, occupancy_R=occ_R,
        M1_mean_interval=m1, M2_interval_entropy=m2,
        M3_c1_ratio=m3, M4_tail_weight=m4,
        M5_log_interval_ratio=m5, M6_weighted_depth_ratio=m6,
        M7_norm_bd_shift=m7,
        C0=counts.get(0), C1=counts.get(1),
        C2=counts.get(2), C3=counts.get(3),
        total_pairs=total_rel,
    )


# ── Analysis ──

ALL_METRICS = [
    ("bd_ratio", "bd_ratio (ref)"),
    ("bdg_d4", "bdg_d4 (ref)"),
    ("bdg_d2c", "bdg_d2c (ref)"),
    ("occupancy_R", "occupancy R (ref)"),
    ("M1_mean_interval", "M1: mean interval size"),
    ("M2_interval_entropy", "M2: interval entropy"),
    ("M3_c1_ratio", "M3: C1/(C0+C1)"),
    ("M4_tail_weight", "M4: tail weight (k≥2)"),
    ("M5_log_interval_ratio", "M5: log(C1/C0)"),
    ("M6_weighted_depth_ratio", "M6: weighted depth ratio"),
    ("M7_norm_bd_shift", "M7: norm BD d2c shift"),
]


def analyze(rows: list[Row], dims: list[int], ns: list[int]) -> list[str]:
    lines: list[str] = []
    lines.append("# Conjecture E — Curvature-Sensitive BD Metrics\n")
    lines.append(f"Total realizations: {len(rows)}\n")

    # ── Monotonicity tournament ──
    lines.append("\n## Monotonicity Tournament: Spearman(H², metric) per (d, N)\n")
    lines.append("Goal: find metrics with **consistent** sign and **high |ρ|** across all slices.\n")

    # Build results table
    results: dict[str, list[tuple[int, int, float, float]]] = {}
    for attr, _ in ALL_METRICS:
        results[attr] = []

    for d in dims:
        for N in ns:
            for attr, _ in ALL_METRICS:
                sub = [r for r in rows if r.d == d and r.N == N
                       and not math.isnan(getattr(r, attr))]
                if len(sub) < 5:
                    results[attr].append((d, N, float("nan"), float("nan")))
                    continue
                h2 = np.array([r.hubble ** 2 for r in sub])
                vals = np.array([getattr(r, attr) for r in sub])
                rho, p = sp_stats.spearmanr(h2, vals)
                results[attr].append((d, N, rho, p))

    # Print per-metric summary
    for attr, label in ALL_METRICS:
        lines.append(f"\n### {label}\n")
        lines.append("| d | N | ρ(H², metric) | p-value | significant? |")
        lines.append("|---|---|--------------|---------|-------------|")

        pass_count = 0
        total_count = 0
        sign_consistent = True
        first_sign = None

        for d, N, rho, p in results[attr]:
            if math.isnan(rho):
                lines.append(f"| {d} | {N} | nan | nan | — |")
                continue
            sig = "✅" if p < 0.05 else "❌"
            if p < 0.05:
                pass_count += 1
            total_count += 1
            if first_sign is None:
                first_sign = 1 if rho > 0 else -1
            elif (rho > 0 and first_sign < 0) or (rho < 0 and first_sign > 0):
                sign_consistent = False
            lines.append(f"| {d} | {N} | {rho:+.3f} | {p:.2e} | {sig} |")

        sign_str = "✅ consistent" if sign_consistent and total_count > 0 else "❌ mixed"
        lines.append(f"\nPassed: **{pass_count}/{total_count}** | Sign: {sign_str}")

    # ── Correlation with R_hat ──
    lines.append("\n## Correlation with R_hat (curvature proxy)\n")
    lines.append("Spearman(metric, R_hat) within each (d, N):\n")

    lines.append("| Metric | d | N | ρ(metric, R_hat) | p | direction |")
    lines.append("|--------|---|---|-----------------|---|-----------|")

    for attr, label in ALL_METRICS:
        for d in dims:
            for N in ns:
                sub = [r for r in rows if r.d == d and r.N == N
                       and not math.isnan(r.R_hat) and not math.isnan(getattr(r, attr))]
                if len(sub) < 5:
                    continue
                rhat = np.array([r.R_hat for r in sub])
                vals = np.array([getattr(r, attr) for r in sub])
                rho, p = sp_stats.spearmanr(vals, rhat)
                direction = "+" if rho > 0.3 else ("-" if rho < -0.3 else "~")
                if p < 0.05 and abs(rho) > 0.3:
                    lines.append(
                        f"| {label[:20]} | {d} | {N} | **{rho:+.3f}** | {p:.1e} | {direction} |"
                    )

    # ── Winner summary ──
    lines.append("\n## Winner Summary\n")
    lines.append("Ranking by (# significant slices, sign consistency, mean |ρ|):\n")
    lines.append("| Rank | Metric | # sig | sign | mean |ρ| | direction |")
    lines.append("|------|--------|-------|------|---------|-----------|")

    scored = []
    for attr, label in ALL_METRICS:
        rhos = [rho for _, _, rho, p in results[attr]
                if not math.isnan(rho)]
        sigs = sum(1 for _, _, rho, p in results[attr]
                   if not math.isnan(rho) and p < 0.05)
        n_total = len(rhos)
        if n_total == 0:
            continue
        signs = [1 if r > 0 else -1 for r in rhos if abs(r) > 0.1]
        sign_ok = len(set(signs)) <= 1 if signs else False
        mean_abs_rho = float(np.mean([abs(r) for r in rhos]))
        dominant_dir = "↓" if sum(signs) < 0 else ("↑" if sum(signs) > 0 else "~")
        scored.append((sigs, sign_ok, mean_abs_rho, label, dominant_dir, n_total))

    scored.sort(key=lambda x: (x[0], x[1], x[2]), reverse=True)
    for rank, (sigs, sign_ok, mar, label, ddir, nt) in enumerate(scored, 1):
        sign_str = "✅" if sign_ok else "❌"
        lines.append(f"| {rank} | {label[:25]} | {sigs}/{nt} | {sign_str} | {mar:.3f} | {ddir} |")

    return lines


def main() -> int:
    ap = argparse.ArgumentParser(description="Curvature-sensitive BD metrics")
    ap.add_argument("--dims", nargs="*", type=int, default=[2, 3, 4])
    ap.add_argument("--ns", nargs="*", type=int, default=[128, 256, 512])
    ap.add_argument("--hubbles", nargs="*", type=float, default=[0.0, 0.25, 0.5, 1.0, 2.0])
    ap.add_argument("--reps", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20260324)
    ap.add_argument("--out", default="outputs_unified_functional/conjecture_e_curvature_sensitive_bd.csv")
    ap.add_argument("--report", default="outputs_unified_functional/conjecture_e_curvature_sensitive_bd.md")
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
                              f"M1={row.M1_mean_interval:.3f} "
                              f"M2={row.M2_interval_entropy:.3f} "
                              f"M3={row.M3_c1_ratio:.3f}")

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
    for line in report_lines[-20:]:
        print(line)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
