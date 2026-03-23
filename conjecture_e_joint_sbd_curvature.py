"""Conjecture E — Joint S_BD × Curvature Experiment.

On the SAME de Sitter sprinklings, compute both:
  1. S_BD (BD action variants: bd_ratio, bdg_d4, bdg_d2c)  — Layer 1/3
  2. R_hat (discrete curvature proxy from interval-volume law) — Layer 2

Then test the core EH recovery prediction:
  S_BD ∝ ∫ R √(-g) d⁴x  ∝  R_dS · V(H,d)

Three key questions:
  Q1: Does S_BD respond monotonically to H (curvature parameter)?
  Q2: Does S_BD correlate with R_hat within each (d, N) slice?
  Q3: Does S_BD / R_dS stabilize as N → ∞ (calibration)?

Design: d ∈ {2,3,4}, N ∈ {128,256,512}, H ∈ {0,0.25,0.5,1.0,2.0}, 8 reps
Total: 3 × 3 × 5 × 8 = 360 realizations
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

# ── Layer 2 imports ──
from curvature_layer2_baseline import fit_curvature_proxy, fit_flat_volume_law
from curvature_layer2_de_sitter_scan import (
    de_sitter_scale_factor,
    sprinkle_de_sitter_like_diamond,
)
from curvature_layer2_recovery_fast import (
    build_causal_matrix_de_sitter,
    compute_tau_and_intervals,
)

# ── Layer 1 imports ──
from bd_action import (
    IntervalCounts,
    bd_ratio_metric,
    bdg_action_d2_corrected,
    bdg_action_d4_standard,
)
from generators import Poset


# ─────────────────────────────────────────────────────────────────────
# Helpers: build Poset + IntervalCounts from de Sitter sprinkling
# ─────────────────────────────────────────────────────────────────────

def poset_from_causal_matrix(causal: np.ndarray) -> Poset:
    """Wrap a boolean causal matrix into a Poset (already transitive)."""
    return Poset(causal)


def interval_counts_from_causal(causal: np.ndarray) -> IntervalCounts:
    """Compute IntervalCounts directly from a boolean causal matrix."""
    c_int = causal.astype(np.int32)
    sizes = c_int @ c_int  # sizes[i,j] = #{k: i≺k≺j}
    ks = sizes[causal].astype(np.int64)
    total = int(ks.size)
    if total == 0:
        return IntervalCounts(counts={}, total_relations=0)
    bc = np.bincount(ks)
    nonzero = np.nonzero(bc)[0]
    out = {int(k): int(bc[k]) for k in nonzero}
    return IntervalCounts(counts=out, total_relations=total)


# ─────────────────────────────────────────────────────────────────────
# Row dataclass
# ─────────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class Row:
    d: int
    N: int
    hubble: float
    rep: int
    seed: int
    # theoretical
    R_dS: float         # d(d-1)H²
    # Layer 2: curvature proxy
    R_hat: float
    flat_prefactor: float
    flat_r2: float
    n_causal_pairs: int
    # Layer 1: BD actions
    bd_ratio: float
    bdg_d4: float
    bdg_d2c: float
    occupancy_R: float  # 1 - f_link = 1 - C0/total
    f_link: float
    # derived
    wall_sigma: float   # σ((R - Rc)/w) with Rc=0.25, w=0.015


def sigmoid(x: float) -> float:
    if x > 500:
        return 1.0
    if x < -500:
        return 0.0
    return 1.0 / (1.0 + math.exp(-x))


def run_single(
    d: int, N: int, hubble: float, rep: int, base_seed: int,
    tau_quantile: float = 0.75,
) -> Row:
    """Run one realization: sprinkle, compute Layer 2 + Layer 1 observables."""
    seed = base_seed + d * 100000 + N * 1000 + rep * 10 + int(hubble * 100)

    # Sprinkle
    points = sprinkle_de_sitter_like_diamond(
        N, d - 1, hubble=hubble, seed=seed,
    )

    # Build causal matrix (coordinate-based, already transitive)
    causal = build_causal_matrix_de_sitter(points, hubble)

    # ── Layer 2: curvature proxy ──
    tau, k = compute_tau_and_intervals(points, hubble, causal)
    n_pairs = int(tau.size)

    if n_pairs < 10:
        R_dS = d * (d - 1) * hubble ** 2
        return Row(
            d=d, N=N, hubble=hubble, rep=rep, seed=seed,
            R_dS=R_dS,
            R_hat=float("nan"), flat_prefactor=float("nan"),
            flat_r2=float("nan"), n_causal_pairs=n_pairs,
            bd_ratio=float("nan"), bdg_d4=float("nan"),
            bdg_d2c=float("nan"), occupancy_R=float("nan"),
            f_link=float("nan"), wall_sigma=float("nan"),
        )

    y = k / max(N - 2, 1)
    x = np.power(tau, d)
    ratio_arr = y / np.clip(x, 1e-12, None)

    fit_mask = tau <= np.quantile(tau, tau_quantile)
    if fit_mask.sum() < 10:
        fit_mask = np.ones_like(tau, dtype=bool)

    flat_pref, _, flat_r2 = fit_flat_volume_law(tau[fit_mask], y[fit_mask], d)
    _, r_hat = fit_curvature_proxy(tau[fit_mask], ratio_arr[fit_mask], d)

    # ── Layer 1: BD actions ──
    counts = interval_counts_from_causal(causal)
    br = bd_ratio_metric(counts)
    d4 = bdg_action_d4_standard(counts, N, normalized=True)
    d2c = bdg_action_d2_corrected(counts, N, normalized=True)

    total_rel = counts.total_relations
    c0 = counts.get(0)
    f_lk = c0 / total_rel if total_rel > 0 else 1.0
    occ_R = 1.0 - f_lk

    R_dS = d * (d - 1) * hubble ** 2
    wall = sigmoid((occ_R - 0.25) / 0.015)

    return Row(
        d=d, N=N, hubble=hubble, rep=rep, seed=seed,
        R_dS=R_dS,
        R_hat=r_hat, flat_prefactor=flat_pref, flat_r2=flat_r2,
        n_causal_pairs=n_pairs,
        bd_ratio=br, bdg_d4=d4, bdg_d2c=d2c,
        occupancy_R=occ_R, f_link=f_lk,
        wall_sigma=wall,
    )


# ─────────────────────────────────────────────────────────────────────
# Analysis
# ─────────────────────────────────────────────────────────────────────

def analyze(rows: list[Row], dims: list[int], ns: list[int],
            hubbles: list[float]) -> list[str]:
    lines: list[str] = []
    lines.append("# Conjecture E — Joint S_BD × Curvature Experiment\n")
    lines.append(f"Total realizations: {len(rows)}\n")

    # ── Q1: S_BD monotonicity with H ──
    lines.append("\n## Q1: Does S_BD respond monotonically to H?\n")
    lines.append("Spearman(H², metric) per (d, N):\n")

    for metric_name in ["bd_ratio", "bdg_d4", "bdg_d2c", "occupancy_R"]:
        lines.append(f"\n### Metric: {metric_name}\n")
        lines.append("| d | N | ρ(H², metric) | p-value | mean@H=0 | mean@H=2 | direction |")
        lines.append("|---|---|--------------|---------|----------|----------|-----------|")

        for d in dims:
            for N in ns:
                sub = [r for r in rows if r.d == d and r.N == N
                       and not math.isnan(getattr(r, metric_name))]
                if len(sub) < 5:
                    continue
                h2 = np.array([r.hubble ** 2 for r in sub])
                vals = np.array([getattr(r, metric_name) for r in sub])
                rho, p = sp_stats.spearmanr(h2, vals)
                m0 = np.mean([getattr(r, metric_name) for r in sub if r.hubble == 0])
                m2 = np.mean([getattr(r, metric_name) for r in sub
                              if abs(r.hubble - 2.0) < 0.01])
                direction = "↑" if rho > 0.3 else ("↓" if rho < -0.3 else "~")
                lines.append(
                    f"| {d} | {N} | {rho:+.3f} | {p:.2e} | {m0:.4f} | {m2:.4f} | {direction} |"
                )

    # ── Q2: S_BD vs R_hat correlation within (d, N) ──
    lines.append("\n## Q2: Does S_BD correlate with R_hat within each (d, N)?\n")
    lines.append("| d | N | ρ(bd_ratio, R_hat) | p | ρ(bdg_d4, R_hat) | p | ρ(occ_R, R_hat) | p | n |")
    lines.append("|---|---|-------------------|---|------------------|---|-----------------|---|---|")

    for d in dims:
        for N in ns:
            sub = [r for r in rows if r.d == d and r.N == N
                   and not math.isnan(r.R_hat) and not math.isnan(r.bd_ratio)]
            if len(sub) < 5:
                continue
            rhat = np.array([r.R_hat for r in sub])
            bd = np.array([r.bd_ratio for r in sub])
            d4 = np.array([r.bdg_d4 for r in sub])
            oR = np.array([r.occupancy_R for r in sub])

            rho_bd, p_bd = sp_stats.spearmanr(bd, rhat)
            rho_d4, p_d4 = sp_stats.spearmanr(d4, rhat)
            rho_oR, p_oR = sp_stats.spearmanr(oR, rhat)

            lines.append(
                f"| {d} | {N} | {rho_bd:+.3f} | {p_bd:.1e} | "
                f"{rho_d4:+.3f} | {p_d4:.1e} | "
                f"{rho_oR:+.3f} | {p_oR:.1e} | {len(sub)} |"
            )

    # ── Q3: S_BD / R_dS calibration ──
    lines.append("\n## Q3: Does S_BD / R_dS stabilize as N grows?\n")
    lines.append("(Only H > 0 rows used; ratio = metric / R_dS)\n")

    for metric_name in ["bd_ratio", "bdg_d4", "bdg_d2c"]:
        lines.append(f"\n### Metric: {metric_name}\n")
        lines.append("| d | N | mean(metric/R_dS) | std | CV | n |")
        lines.append("|---|---|-------------------|-----|-----|---|")

        for d in dims:
            for N in ns:
                sub = [r for r in rows if r.d == d and r.N == N
                       and r.hubble > 0
                       and not math.isnan(getattr(r, metric_name))
                       and r.R_dS > 0]
                if len(sub) < 3:
                    continue
                ratios = np.array([getattr(r, metric_name) / r.R_dS for r in sub])
                m = float(np.mean(ratios))
                s = float(np.std(ratios))
                cv = s / abs(m) if abs(m) > 1e-12 else float("inf")
                lines.append(f"| {d} | {N} | {m:.4f} | {s:.4f} | {cv:.2f} | {len(sub)} |")

    # ── Q2b: Pooled correlation (all d, N combined) ──
    lines.append("\n## Q2b: Pooled Correlation (all data)\n")
    valid = [r for r in rows if not math.isnan(r.R_hat) and not math.isnan(r.bd_ratio)]
    if len(valid) > 10:
        rhat = np.array([r.R_hat for r in valid])
        bd = np.array([r.bd_ratio for r in valid])
        d4 = np.array([r.bdg_d4 for r in valid])
        oR = np.array([r.occupancy_R for r in valid])
        wall = np.array([r.wall_sigma for r in valid])

        for name, arr in [("bd_ratio", bd), ("bdg_d4", d4),
                          ("occupancy_R", oR), ("wall_sigma", wall)]:
            rho, p = sp_stats.spearmanr(arr, rhat)
            lines.append(f"- ρ({name}, R_hat) = {rho:+.3f} (p={p:.2e}, n={len(valid)})")

    # ── Q4: Joint scaling: does bd_ratio track R_dS · V(H)? ──
    lines.append("\n## Q4: bd_ratio vs R_dS (direct proportionality test)\n")
    lines.append("Linear fit: bd_ratio = a · R_dS + b, per (d, N)\n")
    lines.append("| d | N | a (slope) | b (intercept) | R² | Spearman ρ | p |")
    lines.append("|---|---|-----------|---------------|-----|------------|---|")

    for d in dims:
        for N in ns:
            sub = [r for r in rows if r.d == d and r.N == N
                   and not math.isnan(r.bd_ratio)]
            if len(sub) < 5:
                continue
            x = np.array([r.R_dS for r in sub])
            y = np.array([r.bd_ratio for r in sub])
            # include H=0 points
            rho_sp, p_sp = sp_stats.spearmanr(x, y)
            # linear fit
            if np.std(x) > 1e-12:
                slope, intercept, r_val, _, _ = sp_stats.linregress(x, y)
                r2 = r_val ** 2
            else:
                slope, intercept, r2 = float("nan"), float("nan"), float("nan")
            lines.append(
                f"| {d} | {N} | {slope:+.4f} | {intercept:+.4f} | {r2:.3f} | "
                f"{rho_sp:+.3f} | {p_sp:.2e} |"
            )

    # ── Summary ──
    lines.append("\n## Summary\n")
    lines.append("### Key findings:\n")

    # Auto-summarize Q1
    q1_pass = 0
    q1_total = 0
    for d in dims:
        for N in ns:
            sub = [r for r in rows if r.d == d and r.N == N
                   and not math.isnan(r.bd_ratio)]
            if len(sub) < 5:
                continue
            h2 = np.array([r.hubble ** 2 for r in sub])
            bd = np.array([r.bd_ratio for r in sub])
            rho, p = sp_stats.spearmanr(h2, bd)
            q1_total += 1
            if p < 0.05:
                q1_pass += 1
    lines.append(f"- **Q1** (S_BD monotone with H): {q1_pass}/{q1_total} (d,N) slices significant")

    # Auto-summarize Q2
    q2_pass = 0
    q2_total = 0
    for d in dims:
        for N in ns:
            sub = [r for r in rows if r.d == d and r.N == N
                   and not math.isnan(r.R_hat) and not math.isnan(r.bd_ratio)]
            if len(sub) < 5:
                continue
            rhat = np.array([r.R_hat for r in sub])
            bd = np.array([r.bd_ratio for r in sub])
            rho, p = sp_stats.spearmanr(bd, rhat)
            q2_total += 1
            if abs(rho) > 0.3 and p < 0.05:
                q2_pass += 1
    lines.append(f"- **Q2** (S_BD correlates with R_hat): {q2_pass}/{q2_total} slices with |ρ|>0.3 & p<0.05")

    return lines


# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser(description="Joint S_BD × Curvature experiment")
    ap.add_argument("--dims", nargs="*", type=int, default=[2, 3, 4])
    ap.add_argument("--ns", nargs="*", type=int, default=[128, 256, 512])
    ap.add_argument("--hubbles", nargs="*", type=float, default=[0.0, 0.25, 0.5, 1.0, 2.0])
    ap.add_argument("--reps", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20260323)
    ap.add_argument("--tau-quantile", type=float, default=0.75)
    ap.add_argument("--out", default="outputs_unified_functional/conjecture_e_joint_sbd_curvature.csv")
    ap.add_argument("--report", default="outputs_unified_functional/conjecture_e_joint_sbd_curvature.md")
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
                    row = run_single(d, N, H, rep, args.seed, args.tau_quantile)
                    rows.append(row)
                    done += 1
                    if done % 20 == 0 or done == total:
                        print(f"  [{done}/{total}] d={d} N={N} H={H:.2f} rep={rep} "
                              f"R_hat={row.R_hat:.2f} bd_ratio={row.bd_ratio:.4f} "
                              f"occ_R={row.occupancy_R:.3f}")

    # Save CSV
    fn = [f.name for f in fields(Row)]
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fn)
        w.writeheader()
        for r in rows:
            d_row = {}
            for k in fn:
                v = getattr(r, k)
                d_row[k] = f"{v:.6f}" if isinstance(v, float) else str(v)
            w.writerow(d_row)
    print(f"\nSaved: {out_path}")

    # Analysis
    report_lines = analyze(rows, args.dims, args.ns, args.hubbles)

    report_text = "\n".join(report_lines) + "\n"
    report_path = Path(args.report)
    report_path.write_text(report_text, encoding="utf-8")
    print(f"Saved: {report_path}")

    # Console summary
    print("\n" + "=" * 60)
    print("JOINT S_BD × CURVATURE: SUMMARY")
    print("=" * 60)
    for line in report_lines[-10:]:
        print(line)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
