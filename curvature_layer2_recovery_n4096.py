"""Conjecture E — Layer 2: Curvature Recovery for N=4096.

Optimized to avoid the O(N^3) full matrix multiply bottleneck.
Instead of computing interval_matrix = causal @ causal (N^3 ops),
we SAMPLE causal pairs and count intervals per pair using vectorized
column indexing: |I(x,y)| = sum(causal[x,:] & causal[:,y]).

For N=4096, d=4, we typically have ~500K causal pairs.
We sample up to max_pairs (default 50000) for the regression fit,
which is more than enough for stable R_hat estimation.
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
    sprinkle_de_sitter_like_diamond,
)


def build_causal_matrix_de_sitter(points: np.ndarray, hubble: float) -> np.ndarray:
    """Build boolean causal relation matrix directly from coordinates."""
    t = points[:, 0]
    spatial = points[:, 1:]
    dt = t[None, :] - t[:, None]

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


def compute_tau_and_intervals_sampled(
    points: np.ndarray, hubble: float, causal: np.ndarray,
    max_pairs: int = 50000, seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Compute proper-time proxy and interval sizes for SAMPLED causal pairs.
    
    Returns (tau_array, k_array, total_causal_pairs).
    """
    t = points[:, 0]
    spatial = points[:, 1:]
    N = len(t)

    # Find all causal pairs
    ii, jj = np.where(causal)
    n_total = len(ii)

    if n_total == 0:
        return np.array([]), np.array([]), 0

    # Sample if too many
    rng = np.random.RandomState(seed)
    if n_total > max_pairs:
        idx = rng.choice(n_total, max_pairs, replace=False)
        ii_s = ii[idx]
        jj_s = jj[idx]
    else:
        ii_s = ii
        jj_s = jj

    # Proper-time proxy for sampled pairs
    if hubble == 0.0:
        dt_s = t[jj_s] - t[ii_s]
        sp_diff = spatial[jj_s] - spatial[ii_s]
        sp_d2 = np.sum(sp_diff ** 2, axis=1)
        tau = np.sqrt(np.clip(dt_s * dt_s - sp_d2, 0.0, None))
    else:
        chi = (np.exp(-hubble * t[ii_s]) - np.exp(-hubble * t[jj_s])) / hubble
        tau = np.clip(chi, 0.0, None)

    # Interval sizes for sampled pairs: |I(x,y)| = sum_z (causal[x,z] & causal[z,y])
    # Vectorized batch computation using einsum for the diagonal.
    batch_size = 2000
    k = np.zeros(len(ii_s), dtype=np.float64)
    
    causal_uint8 = causal.astype(np.uint8)  # for faster arithmetic
    
    for start in range(0, len(ii_s), batch_size):
        end = min(start + batch_size, len(ii_s))
        batch_i = ii_s[start:end]
        batch_j = jj_s[start:end]
        
        # rows[p, z] = causal[i_p, z],  cols[z, p] = causal[z, j_p]
        rows = causal_uint8[batch_i]       # (batch, N)
        cols = causal_uint8[:, batch_j]    # (N, batch)
        
        # k[p] = sum_z rows[p,z] * cols[z,p] = einsum('pz,zp->p', rows, cols)
        k[start:end] = np.einsum('pz,zp->p', rows, cols).astype(np.float64)

    return tau, k, n_total


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
    n_sampled_pairs: int
    intercept: float


def run_single(
    d: int, N: int, hubble: float, rep: int, seed: int,
    tau_quantile: float, max_pairs: int,
) -> Row:
    """Run one realization."""
    print(f"    Starting d={d} N={N} H={hubble:.1f} rep={rep}...", flush=True)
    
    points = sprinkle_de_sitter_like_diamond(
        N, d - 1, hubble=hubble,
        seed=seed + d * 100000 + N * 100 + rep + int(hubble * 1000),
    )
    
    print(f"    Building causal matrix...", flush=True)
    causal = build_causal_matrix_de_sitter(points, hubble)

    print(f"    Computing intervals (sampled)...", flush=True)
    pair_seed = seed + d * 1000 + N + rep * 7 + int(hubble * 100)
    tau, k, n_total = compute_tau_and_intervals_sampled(
        points, hubble, causal, max_pairs=max_pairs, seed=pair_seed,
    )
    n_sampled = len(tau)

    if tau.size < 10:
        R_dS = d * (d - 1) * hubble ** 2
        return Row(d=d, N=N, hubble=hubble, rep=rep, R_dS=R_dS,
                   R_hat=float('nan'), ratio=float('nan'),
                   flat_prefactor=float('nan'), flat_r2=float('nan'),
                   n_causal_pairs=n_total, n_sampled_pairs=n_sampled,
                   intercept=float('nan'))

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
        n_causal_pairs=n_total, n_sampled_pairs=n_sampled,
        intercept=intercept,
    )


def main() -> int:
    ap = argparse.ArgumentParser(description="Layer 2: Curvature Recovery N=4096 (Sampled)")
    ap.add_argument("--dims", nargs="*", type=int, default=[4])
    ap.add_argument("--ns", nargs="*", type=int, default=[4096])
    ap.add_argument("--hubbles", nargs="*", type=float, default=[0.0, 0.25, 0.5, 1.0, 2.0])
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--max-pairs", type=int, default=50000)
    ap.add_argument("--tau-quantile", type=float, default=0.75)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out", default="outputs_unified_functional/curvature_layer2_recovery_n4096.csv")
    ap.add_argument("--report", default="outputs_unified_functional/curvature_layer2_recovery_n4096.md")
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
                    row = run_single(d, N, hubble, rep, args.seed,
                                     args.tau_quantile, args.max_pairs)
                    rows.append(row)
                    done += 1
                    print(f"  [{done}/{total}] d={d} N={N} H={hubble:.1f} rep={rep} "
                          f"R_hat={row.R_hat:+.1f} total_pairs={row.n_causal_pairs} "
                          f"sampled={row.n_sampled_pairs}", flush=True)

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
    lines.append("# Layer 2: N=4096 Curvature Recovery (d=4)\n")
    lines.append(f"Total: {len(rows)} realizations, max_pairs={args.max_pairs}\n")

    # R_hat table
    lines.append("\n## R_hat vs H\n")
    lines.append("| d | N | H | R_dS | mean R_hat | std | total_pairs | sampled |")
    lines.append("|---|---|---|------|-----------|-----|-------------|---------|")
    for d in args.dims:
        for N in args.ns:
            for hubble in args.hubbles:
                subset = [r for r in rows if r.d == d and r.N == N and r.hubble == hubble]
                rh = [r.R_hat for r in subset if not math.isnan(r.R_hat)]
                if not rh:
                    continue
                R_dS = d * (d - 1) * hubble ** 2
                tp = [r.n_causal_pairs for r in subset]
                sp = [r.n_sampled_pairs for r in subset]
                lines.append(
                    f"| {d} | {N} | {hubble:.2f} | {R_dS:.2f} | "
                    f"{np.mean(rh):+.2f} | {np.std(rh):.2f} | "
                    f"{int(np.mean(tp))} | {int(np.mean(sp))} |"
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

    # Calibration
    lines.append("\n## Calibration\n")
    lines.append("Prior c_eff values: d=4 N=256→-2, N=512→7, N=1024→4, N=2048→5.7\n")
    for d in args.dims:
        for N in args.ns:
            subset_h0 = [r for r in rows if r.d == d and r.N == N and r.hubble == 0.0]
            rh_flat = [r.R_hat for r in subset_h0 if not math.isnan(r.R_hat)]
            bias = float(np.mean(rh_flat)) if rh_flat else 0.0
            bias_std = float(np.std(rh_flat)) if rh_flat else 0.0

            subset_h = [r for r in rows if r.d == d and r.N == N and r.hubble > 0]
            if len(subset_h) < 4:
                continue
            R_dS_arr = np.array([r.R_dS for r in subset_h])
            R_hat_corr = np.array([r.R_hat - bias for r in subset_h
                                   if not math.isnan(r.R_hat)])
            R_dS_arr = R_dS_arr[:len(R_hat_corr)]

            if np.sum(R_dS_arr ** 2) > 0:
                c_eff = float(np.sum(R_dS_arr * R_hat_corr) / np.sum(R_dS_arr ** 2))
                pred = c_eff * R_dS_arr
                ss_res = float(np.sum((R_hat_corr - pred) ** 2))
                ss_tot = float(np.sum((R_hat_corr - np.mean(R_hat_corr)) ** 2))
                r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
                n_obs = len(R_hat_corr)
                se = math.sqrt(ss_res / max(n_obs - 1, 1) / np.sum(R_dS_arr ** 2))

                lines.append(f"- d={d}, N={N}: bias={bias:+.1f}±{bias_std:.1f}, "
                             f"**c_eff={c_eff:.1f}** ± {se:.1f}, R²={r2:.3f}")
                lines.append(f"- c_eff trajectory: -2 → 7 → 4 → 5.7 → **{c_eff:.1f}**")

    # Slope comparison
    lines.append("\n## Slope Stability (all N)\n")
    prior = {
        (4, 256): -33.5, (4, 512): 88.0, (4, 1024): 51.3, (4, 2048): 74.1,
    }
    for d in args.dims:
        lines.append(f"\n**d={d}**:")
        for N in [256, 512, 1024, 2048]:
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
    print("N=4096 SUMMARY")
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
