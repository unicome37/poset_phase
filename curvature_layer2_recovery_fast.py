"""Conjecture E — Layer 2: Curvature Recovery (Fast, Large-N version).

Optimized for N >= 2048 by skipping transitive closure.
For coordinate-aware sprinklings, the causal order from coordinates IS
already transitive — no Warshall needed. Interval sizes are computed
directly via matrix multiplication on the adjacency matrix.

This script is functionally identical to curvature_layer2_recovery.py
but avoids the O(N^3) transitive closure bottleneck.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np

from curvature_layer2_baseline import fit_curvature_proxy, fit_flat_volume_law
from curvature_layer2_de_sitter_scan import (
    de_sitter_scale_factor,
    sprinkle_de_sitter_like_diamond,
)


def build_causal_matrix_de_sitter(points: np.ndarray, hubble: float) -> np.ndarray:
    """Build boolean causal relation matrix directly from coordinates.
    
    For coordinate sprinklings the causal order is already transitive,
    so this IS the transitive closure — no Warshall needed.
    Returns: bool array of shape (N, N) where [i,j] = True iff i ≺ j.
    """
    t = points[:, 0]
    spatial = points[:, 1:]
    dt = t[None, :] - t[:, None]  # dt[i,j] = t_j - t_i

    if hubble == 0.0:
        horizon = np.clip(dt, 0.0, None)
    else:
        ti = t[:, None]
        tj = t[None, :]
        horizon = (np.exp(-hubble * ti) - np.exp(-hubble * tj)) / hubble
        horizon = np.clip(horizon, 0.0, None)

    spatial_d2 = np.sum((spatial[None, :, :] - spatial[:, None, :]) ** 2, axis=2)
    causal = (dt > 0.0) & (spatial_d2 <= horizon * horizon)
    return causal


def compute_tau_and_intervals(
    points: np.ndarray, hubble: float, causal: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Compute proper-time proxy and interval sizes for all causal pairs."""
    t = points[:, 0]
    spatial = points[:, 1:]
    N = len(t)

    # Proper-time / conformal-distance proxy
    dt = t[None, :] - t[:, None]
    if hubble == 0.0:
        tau_matrix = np.sqrt(np.clip(dt * dt - np.sum(
            (spatial[None, :, :] - spatial[:, None, :]) ** 2, axis=2
        ), 0.0, None))
    else:
        ti = t[:, None]
        tj = t[None, :]
        chi = (np.exp(-hubble * ti) - np.exp(-hubble * tj)) / hubble
        tau_matrix = np.clip(chi, 0.0, None)

    # Interval sizes: |I(x,y)| = number of z with x ≺ z ≺ y
    # = sum_z causal[x,z] & causal[z,y]
    # = (causal @ causal)[x,y]  (treating bool as int)
    c_int = causal.astype(np.int32)
    interval_matrix = c_int @ c_int

    # Extract values for causal pairs only
    tau = tau_matrix[causal]
    k = interval_matrix[causal].astype(float)
    return tau, k


@dataclass(frozen=True)
class Row:
    d: int
    N: int
    hubble: float
    rep: int
    R_dS: float
    R_hat: float
    ratio: float
    flat_prefactor: float
    flat_r2: float
    n_causal_pairs: int
    intercept: float


def run_single(
    d: int, N: int, hubble: float, rep: int, seed: int, tau_quantile: float
) -> Row:
    """Run one realization."""
    points = sprinkle_de_sitter_like_diamond(
        N, d - 1, hubble=hubble,
        seed=seed + d * 100000 + N * 100 + rep + int(hubble * 1000),
    )
    causal = build_causal_matrix_de_sitter(points, hubble)
    tau, k = compute_tau_and_intervals(points, hubble, causal)

    if tau.size < 10:
        R_dS = d * (d - 1) * hubble ** 2
        return Row(d=d, N=N, hubble=hubble, rep=rep, R_dS=R_dS,
                   R_hat=float('nan'), ratio=float('nan'),
                   flat_prefactor=float('nan'), flat_r2=float('nan'),
                   n_causal_pairs=int(tau.size), intercept=float('nan'))

    y = k / max(N - 2, 1)
    x = np.power(tau, d)
    ratio_arr = y / np.clip(x, 1e-12, None)

    fit_mask = tau <= np.quantile(tau, tau_quantile)
    if fit_mask.sum() < 10:
        fit_mask = np.ones_like(tau, dtype=bool)

    flat_pref, _, flat_r2 = fit_flat_volume_law(tau[fit_mask], y[fit_mask], d)
    intercept, r_hat = fit_curvature_proxy(tau[fit_mask], ratio_arr[fit_mask], d)

    R_dS = d * (d - 1) * hubble ** 2
    r_ratio = r_hat / R_dS if R_dS > 0 else float('nan')

    return Row(
        d=d, N=N, hubble=hubble, rep=rep,
        R_dS=R_dS, R_hat=r_hat, ratio=r_ratio,
        flat_prefactor=flat_pref, flat_r2=flat_r2,
        n_causal_pairs=int(tau.size), intercept=intercept,
    )


def main() -> int:
    ap = argparse.ArgumentParser(description="Layer 2: Curvature Recovery (Fast)")
    ap.add_argument("--dims", nargs="*", type=int, default=[2, 3, 4])
    ap.add_argument("--ns", nargs="*", type=int, default=[2048])
    ap.add_argument("--hubbles", nargs="*", type=float, default=[0.0, 0.25, 0.5, 1.0, 2.0])
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--tau-quantile", type=float, default=0.75)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out", default="outputs_unified_functional/curvature_layer2_recovery_n2048.csv")
    ap.add_argument("--report", default="outputs_unified_functional/curvature_layer2_recovery_n2048.md")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[Row] = []
    total = len(args.dims) * len(args.ns) * len(args.hubbles) * args.reps
    done = 0

    for d in args.dims:
        for N in args.ns:
            for hubble in args.hubbles:
                for rep in range(args.reps):
                    row = run_single(d, N, hubble, rep, args.seed, args.tau_quantile)
                    rows.append(row)
                    done += 1
                    if done % 5 == 0 or done == total:
                        print(f"  [{done}/{total}] d={d} N={N} H={hubble:.1f} "
                              f"R_hat={row.R_hat:+.1f} pairs={row.n_causal_pairs}")

    # Save raw CSV
    fieldnames = [f.name for f in fields(Row)]
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow({fn: getattr(row, fn) for fn in fieldnames})
    print(f"\nSaved raw: {out_path}")

    # Generate report
    from scipy import stats as sp_stats
    report_path = Path(args.report)
    lines: list[str] = []
    lines.append("# Layer 2: N=2048 Curvature Recovery\n")
    lines.append(f"Total: {len(rows)} realizations\n")

    # R_hat table
    lines.append("\n## R_hat vs H\n")
    lines.append("| d | N | H | H^2 | R_dS | mean R_hat | std | pairs |")
    lines.append("|---|---|---|-----|------|-----------|-----|-------|")
    for d in args.dims:
        for N in args.ns:
            for hubble in args.hubbles:
                subset = [r for r in rows if r.d == d and r.N == N and r.hubble == hubble]
                rh = [r.R_hat for r in subset if not math.isnan(r.R_hat)]
                np_val = [r.n_causal_pairs for r in subset]
                if not rh:
                    continue
                R_dS = d * (d - 1) * hubble ** 2
                lines.append(
                    f"| {d} | {N} | {hubble:.2f} | {hubble**2:.3f} | {R_dS:.2f} | "
                    f"{np.mean(rh):+.2f} | {np.std(rh):.2f} | {int(np.mean(np_val))} |"
                )

    # Monotonicity
    lines.append("\n## Monotonicity\n")
    lines.append("| d | N | Spearman | p-value | slope | mono? |")
    lines.append("|---|---|---------|---------|-------|-------|")
    slope_table = {}
    for d in args.dims:
        for N in args.ns:
            h2_vals, rh_vals = [], []
            for hubble in args.hubbles:
                subset = [r for r in rows if r.d == d and r.N == N and r.hubble == hubble]
                rh = [r.R_hat for r in subset if not math.isnan(r.R_hat)]
                for v in rh:
                    h2_vals.append(hubble ** 2)
                    rh_vals.append(v)
            if len(h2_vals) < 5:
                continue
            rho_s, p_s = sp_stats.spearmanr(h2_vals, rh_vals)
            sl, ic, _, _, _ = sp_stats.linregress(h2_vals, rh_vals)
            slope_table[(d, N)] = sl
            mono = "YES" if rho_s > 0.5 and p_s < 0.05 else "NO"
            lines.append(f"| {d} | {N} | {rho_s:+.3f} | {p_s:.2e} | {sl:+.1f} | {mono} |")

    # Cross-N slope comparison (include prior results)
    lines.append("\n## Slope Stability (all N)\n")
    lines.append("Prior slopes from N=256/512/1024 run:\n")
    prior = {
        (2, 256): 416.2, (2, 512): 405.3, (2, 1024): 425.2,
        (3, 256): 277.2, (3, 512): 213.1, (3, 1024): 264.5,
        (4, 256): -33.5, (4, 512): 88.0, (4, 1024): 51.3,
    }
    for d in args.dims:
        lines.append(f"\n**d={d}**:")
        for N in [256, 512, 1024]:
            sl = prior.get((d, N), float('nan'))
            lines.append(f"  N={N:5d}: slope = {sl:+.1f}")
        for N in args.ns:
            sl = slope_table.get((d, N), float('nan'))
            lines.append(f"  N={N:5d}: slope = {sl:+.1f}  <-- NEW")

    report_text = "\n".join(lines) + "\n"
    report_path.write_text(report_text, encoding="utf-8")
    print(f"Saved report: {report_path}")

    # Console summary
    print("\n" + "=" * 60)
    print("N=2048 SUMMARY")
    print("=" * 60)
    for d in args.dims:
        for N in args.ns:
            sl = slope_table.get((d, N), float('nan'))
            subset_flat = [r for r in rows if r.d == d and r.N == N and r.hubble == 0.0]
            rh_flat = [r.R_hat for r in subset_flat if not math.isnan(r.R_hat)]
            bias = np.mean(rh_flat) if rh_flat else float('nan')
            print(f"  d={d} N={N}: slope={sl:+.1f}  bias(H=0)={bias:+.1f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
