"""Conjecture E — Layer 2: Curvature Recovery Experiment.

Production-scale test of whether the discrete curvature proxy R_hat
responds monotonically to the de Sitter Hubble parameter H, and whether
the response slope stabilises as N → ∞.

The raw R_hat carries an unknown calibration constant that depends on
coordinate choice and finite-size effects.  Therefore the correct test
is NOT "R_hat = R_dS" but:
  (a) R_hat is monotonically increasing in H² at each N;
  (b) the slope dR_hat / d(H²) stabilises as N grows;
  (c) R_hat(H=0) → 0 as N → ∞ (flat baseline control).

Experiment design:
  - Dimensions: d = 2, 3, 4
  - N values: 256, 512, 1024 (finite-size scaling)
  - H values: 0.0, 0.25, 0.5, 1.0, 2.0 (curvature sweep)
  - Reps: 10 per (d, N, H) for error bars
  - Output: per-rep R_hat, plus summary with monotonicity/slope analysis

Key metrics:
  1. Spearman(R_hat, H²) at each (d, N) — should be positive & significant
  2. Linear slope dR_hat/d(H²) at each (d, N) — should stabilise with N
  3. R_hat(H=0) bias — should shrink with N
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
    de_sitter_interval_proxy,
    poset_from_de_sitter_points,
    sprinkle_de_sitter_like_diamond,
)


@dataclass(frozen=True)
class Row:
    d: int
    N: int
    hubble: float
    rep: int
    R_dS: float          # analytic: d(d-1)H²
    R_hat: float          # fitted curvature proxy
    ratio: float          # R_hat / R_dS (nan if R_dS=0)
    flat_prefactor: float
    flat_r2: float
    n_causal_pairs: int
    intercept: float


def run_single(
    d: int, N: int, hubble: float, rep: int, seed: int, tau_quantile: float
) -> Row:
    """Run one realization and return results."""
    points = sprinkle_de_sitter_like_diamond(
        N, d - 1, hubble=hubble,
        seed=seed + d * 100000 + N * 100 + rep + int(hubble * 1000),
    )
    poset = poset_from_de_sitter_points(points, hubble)

    # Get causal distance proxy
    t = points[:, 0]
    spatial = points[:, 1:]
    dt = t[None, :] - t[:, None]
    if hubble == 0.0:
        tau_all = np.sqrt(np.clip(dt * dt - np.sum(
            (spatial[None, :, :] - spatial[:, None, :]) ** 2, axis=2
        ), 0.0, None))
    else:
        ti = t[:, None]
        tj = t[None, :]
        chi = (np.exp(-hubble * ti) - np.exp(-hubble * tj)) / hubble
        tau_all = np.clip(chi, 0.0, None)

    causal = poset.closure
    tau = tau_all[causal]
    # Interval sizes
    c = poset.closure.astype(np.int32)
    k = (c @ c)[causal].astype(float)

    if tau.size < 10:
        return Row(d=d, N=N, hubble=hubble, rep=rep,
                   R_dS=d * (d - 1) * hubble ** 2,
                   R_hat=float('nan'), ratio=float('nan'),
                   flat_prefactor=float('nan'), flat_r2=float('nan'),
                   n_causal_pairs=int(tau.size), intercept=float('nan'))

    y = k / max(N - 2, 1)
    x = np.power(tau, d)
    ratio_arr = y / np.clip(x, 1e-12, None)

    # Fit in lower tau quantile
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
    ap = argparse.ArgumentParser(description="Layer 2: Curvature Recovery")
    ap.add_argument("--dims", nargs="*", type=int, default=[2, 3, 4])
    ap.add_argument("--ns", nargs="*", type=int, default=[256, 512, 1024])
    ap.add_argument("--hubbles", nargs="*", type=float, default=[0.0, 0.25, 0.5, 1.0, 2.0])
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--tau-quantile", type=float, default=0.75)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out", default="outputs_unified_functional/curvature_layer2_recovery.csv")
    ap.add_argument("--report", default="outputs_unified_functional/curvature_layer2_recovery.md")
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
                    if done % 50 == 0 or done == total:
                        print(f"  [{done}/{total}] d={d} N={N} H={hubble:.1f} rep={rep}")

    # Save raw CSV
    fieldnames = [f.name for f in fields(Row)]
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow({fn: getattr(row, fn) for fn in fieldnames})
    print(f"\nSaved raw: {out_path}")

    # Generate summary report
    report_path = Path(args.report)
    lines: list[str] = []
    lines.append("# Conjecture E — Layer 2: Curvature Recovery Report\n")
    lines.append(f"Total realizations: {len(rows)}\n")

    # === Flat baseline (H=0) ===
    lines.append("\n## Flat Baseline (H=0)\n")
    lines.append("| d | N | mean R_hat | std R_hat | flat_r2 | n_pairs |")
    lines.append("|---|---|-----------|-----------|---------|---------|")
    for d in args.dims:
        for N in args.ns:
            subset = [r for r in rows if r.d == d and r.N == N and r.hubble == 0.0]
            if not subset:
                continue
            rh = [r.R_hat for r in subset if not math.isnan(r.R_hat)]
            fr2 = [r.flat_r2 for r in subset if not math.isnan(r.flat_r2)]
            np_val = [r.n_causal_pairs for r in subset]
            lines.append(
                f"| {d} | {N} | {np.mean(rh):+.3f} | {np.std(rh):.3f} | "
                f"{np.mean(fr2):.3f} | {int(np.mean(np_val))} |"
            )

    # === Full R_hat table ===
    lines.append("\n## R_hat vs H at each (d, N)\n")
    lines.append("| d | N | H | H² | R_dS | mean R_hat | std R_hat | n_pairs |")
    lines.append("|---|---|---|-----|------|-----------|-----------|---------|")
    for d in args.dims:
        for N in args.ns:
            for hubble in args.hubbles:
                subset = [r for r in rows if r.d == d and r.N == N and r.hubble == hubble]
                if not subset:
                    continue
                R_dS = d * (d - 1) * hubble ** 2
                rh = [r.R_hat for r in subset if not math.isnan(r.R_hat)]
                np_val = [r.n_causal_pairs for r in subset]
                if not rh:
                    continue
                lines.append(
                    f"| {d} | {N} | {hubble:.2f} | {hubble**2:.3f} | {R_dS:.2f} | "
                    f"{np.mean(rh):+.2f} | {np.std(rh):.2f} | {int(np.mean(np_val))} |"
                )

    # === Monotonicity & slope analysis ===
    from scipy import stats as sp_stats
    lines.append("\n## Monotonicity Test: Spearman(R_hat, H²) at each (d, N)\n")
    lines.append("| d | N | Spearman ρ | p-value | slope dR_hat/d(H²) | intercept | monotone? |")
    lines.append("|---|---|-----------|---------|--------------------|-----------|-----------|")  
    slope_table = {}  # (d, N) -> slope
    for d in args.dims:
        for N in args.ns:
            h2_vals = []
            rh_vals = []
            for hubble in args.hubbles:
                subset = [r for r in rows if r.d == d and r.N == N and r.hubble == hubble]
                rh = [r.R_hat for r in subset if not math.isnan(r.R_hat)]
                if rh:
                    for v in rh:
                        h2_vals.append(hubble ** 2)
                        rh_vals.append(v)
            if len(h2_vals) < 5:
                continue
            h2_arr = np.array(h2_vals)
            rh_arr = np.array(rh_vals)
            rho_s, p_s = sp_stats.spearmanr(h2_arr, rh_arr)
            # Linear fit
            sl, ic, _, _, _ = sp_stats.linregress(h2_arr, rh_arr)
            slope_table[(d, N)] = sl
            mono = "✅" if rho_s > 0.5 and p_s < 0.05 else "❌"
            lines.append(
                f"| {d} | {N} | {rho_s:+.3f} | {p_s:.2e} | {sl:+.2f} | {ic:+.2f} | {mono} |"
            )

    # === Slope stability ===
    lines.append("\n## Slope Stability: dR_hat/d(H²) vs N\n")
    lines.append("Does the slope stabilise as N grows?\n")
    for d in args.dims:
        slopes = [(N, slope_table.get((d, N), float('nan'))) for N in args.ns
                  if (d, N) in slope_table]
        if len(slopes) < 2:
            continue
        lines.append(f"\n**d={d}**:")
        for N, sl in slopes:
            lines.append(f"  N={N:5d}: slope = {sl:+.2f}")
        # Check if slope is becoming more stable
        if len(slopes) >= 3:
            diffs = [abs(slopes[i+1][1] - slopes[i][1]) for i in range(len(slopes)-1)]
            if diffs[-1] < diffs[0]:
                lines.append(f"  → slope **stabilising** (Δ: {diffs[0]:.2f} → {diffs[-1]:.2f})")
            else:
                lines.append(f"  → slope not yet stable (Δ: {diffs[0]:.2f} → {diffs[-1]:.2f})")

    # === Flat baseline bias ===
    lines.append("\n## Flat Baseline Bias: R_hat(H=0) vs N\n")
    lines.append("Does R_hat(H=0) → 0 as N grows?\n")
    for d in args.dims:
        biases = []
        for N in args.ns:
            subset = [r for r in rows if r.d == d and r.N == N and r.hubble == 0.0]
            rh = [r.R_hat for r in subset if not math.isnan(r.R_hat)]
            if rh:
                biases.append((N, np.mean(rh), np.std(rh)))
        if biases:
            lines.append(f"\n**d={d}**:")
            for N, m, s in biases:
                lines.append(f"  N={N:5d}: R_hat(H=0) = {m:+.2f} ± {s:.2f}")
            if len(biases) >= 2:
                if abs(biases[-1][1]) < abs(biases[0][1]):
                    lines.append(f"  → bias **shrinking**")
                else:
                    lines.append(f"  → bias not shrinking")

    # === Summary ===
    lines.append("\n## Summary\n")
    n_mono = 0
    n_total = 0
    for d in args.dims:
        for N in args.ns:
            h2_vals = []
            rh_vals = []
            for hubble in args.hubbles:
                subset = [r for r in rows if r.d == d and r.N == N and r.hubble == hubble]
                rh = [r.R_hat for r in subset if not math.isnan(r.R_hat)]
                if rh:
                    h2_vals.extend([hubble**2] * len(rh))
                    rh_vals.extend(rh)
            if len(h2_vals) >= 5:
                n_total += 1
                rho_s, p_s = sp_stats.spearmanr(h2_vals, rh_vals)
                if rho_s > 0.5 and p_s < 0.05:
                    n_mono += 1
    lines.append(f"- Monotonicity (Spearman > 0.5, p < 0.05): **{n_mono}/{n_total}**")
    lines.append(f"- Slope values: see table above")
    lines.append(f"- Flat bias: see table above")

    report_text = "\n".join(lines) + "\n"
    report_path.write_text(report_text, encoding="utf-8")
    print(f"Saved report: {report_path}")

    # Print summary to console
    print("\n" + "=" * 70)
    print("CURVATURE RECOVERY SUMMARY")
    print("=" * 70)
    print(f"Monotonicity: {n_mono}/{n_total}")
    for d in args.dims:
        print(f"\n  d={d}:")
        for N in args.ns:
            sl = slope_table.get((d, N), float('nan'))
            subset_flat = [r for r in rows if r.d == d and r.N == N and r.hubble == 0.0]
            bias = np.mean([r.R_hat for r in subset_flat if not math.isnan(r.R_hat)]) if subset_flat else float('nan')
            print(f"    N={N}: slope={sl:+.1f}  bias(H=0)={bias:+.1f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
