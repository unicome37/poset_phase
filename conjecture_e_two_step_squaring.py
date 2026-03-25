"""Conjecture E — T3: Two-Step Analytical Squaring Experiment.

Tests the simplest possible H→R bridge: since b1_std tracks H (first-order),
calibrate Ĥ from b1_std, then compute R̂ = d(d-1)·Ĥ² analytically.

The idea:
  Step 1 (calibration): Within each (d, N) slice, fit b1_std → Ĥ using
          a monotone calibration (OLS on H or median-based).
  Step 2 (squaring):    R̂ = d(d-1) · Ĥ²

If this two-step procedure achieves R²(R̂, H²) > R²(b1_std, H²), then
the bridge is not intrinsically stiff — it just requires explicit squaring
of a first-order observable. If it doesn't help, the stiffness is deeper.

Strategies:
  S1: Global OLS calibration
      Fit Ĥ = a·|b1_std| + b on training data, then R̂ = d(d-1)·Ĥ²
  S2: Leave-one-out cross-validated OLS
      For each realization, fit calibration on all OTHER realizations at same (d,N)
  S3: Rank-preserving (Spearman-optimal)
      Map b1_std to Ĥ via rank-matching (preserves Spearman ρ perfectly)
  S4: Group-mean calibration
      Average b1_std per H level, fit on group means, apply to individuals
  S5: Antichain channel (w_max_ratio → Ĥ → R̂)
      Same two-step but using w_max_ratio as the first-order observable

For each strategy, we measure:
  - Spearman ρ(R̂, H²) — should be same as ρ(observable, H) if monotone
  - R²(R̂, H²) vs R²(raw_observable, H²) — the KEY test
  - R²(R̂, R_dS) — direct match to target
  - α_eff from power-law scan on R̂
  - Calibration stability: std of fitted slope across LOO folds

Design:
  - d = 4 (critical dimension)
  - N = 128, 256, 512, 1024
  - H = 0, 0.25, 0.5, 1.0, 2.0
  - 16 reps per cell → 320 total realizations
  - Reuses T1/T5 sprinkling infrastructure

References:
  - T1 (§4.1.38): raw b1_std is the unique spectral bulk observable
  - T5 (§4.1.37): α_eff monotonically increases with N
  - §4.1.32: b1_std α≈1.00 at d=4 (tracks H, not H²)
"""

from __future__ import annotations

import argparse
import csv
import time
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats
from scipy.optimize import curve_fit

from curvature_layer2_de_sitter_scan import sprinkle_de_sitter_like_diamond
from curvature_layer2_recovery_fast import build_causal_matrix_de_sitter
from conjecture_e_sorkin_dalembertian import build_dalembertian_matrix
from conjecture_e_antichain_structure import compute_antichain_features


# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class T3Row:
    d: int
    N: int
    hubble: float
    rep: int
    R_dS: float
    H2: float
    n_causal_pairs: int
    b1_mean: float
    b1_std: float
    w_max_ratio: float
    mean_layer_width: float
    layer_width_std: float


def run_single(d: int, N: int, hubble: float, rep: int, seed: int) -> T3Row:
    points = sprinkle_de_sitter_like_diamond(
        N, d - 1, hubble=hubble,
        seed=seed + d * 100000 + N * 100 + rep + int(hubble * 1000),
    )
    causal = build_causal_matrix_de_sitter(points, hubble)
    n_causal_pairs = int(np.sum(causal))
    R_dS = d * (d - 1) * hubble ** 2
    H2 = hubble ** 2

    # B_ℓ spectral features
    B = build_dalembertian_matrix(causal, d, rho=1.0)
    ones = np.ones(N)
    b1 = B @ ones
    b1_mean = float(np.mean(b1))
    b1_std = float(np.std(b1))

    # Antichain features
    ac_feats = compute_antichain_features(causal, compute_dilworth=(N <= 600))

    return T3Row(
        d=d, N=N, hubble=hubble, rep=rep,
        R_dS=R_dS, H2=H2,
        n_causal_pairs=n_causal_pairs,
        b1_mean=b1_mean, b1_std=b1_std,
        w_max_ratio=ac_feats.get("w_max_ratio", float("nan")),
        mean_layer_width=ac_feats.get("mean_layer_width", float("nan")),
        layer_width_std=ac_feats.get("layer_width_std", float("nan")),
    )


# ---------------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------------
def ols_r2(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.std(x) < 1e-15 or np.std(y) < 1e-15:
        return 0.0
    coeffs = np.polyfit(x, y, 1)
    pred = np.polyval(coeffs, x)
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return max(0.0, 1.0 - ss_res / ss_tot) if ss_tot > 1e-15 else 0.0


def ols_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Return (slope, intercept) of OLS fit y = slope*x + intercept."""
    if len(x) < 3 or np.std(x) < 1e-15:
        return (0.0, np.mean(y))
    coeffs = np.polyfit(x, y, 1)
    return (float(coeffs[0]), float(coeffs[1]))


def alpha_grid(H_arr: np.ndarray, Y_arr: np.ndarray,
               alpha_range: np.ndarray | None = None) -> dict:
    if alpha_range is None:
        alpha_range = np.arange(0.25, 8.05, 0.25)
    results = {}
    for alpha in alpha_range:
        x = np.power(np.abs(H_arr), alpha)
        results[round(alpha, 2)] = ols_r2(x, Y_arr)
    return results


def best_alpha(ag: dict) -> tuple[float, float]:
    sorted_ag = sorted(ag.items(), key=lambda x: -x[1])
    return sorted_ag[0]


# ---------------------------------------------------------------------------
# Two-step squaring strategies
# ---------------------------------------------------------------------------
def strategy_s1_global_ols(obs: np.ndarray, H: np.ndarray, d: int) -> np.ndarray:
    """S1: Global OLS calibration. Fit |obs| → H, then R̂ = d(d-1)·Ĥ²."""
    slope, intercept = ols_fit(np.abs(obs), H)
    H_hat = slope * np.abs(obs) + intercept
    H_hat = np.maximum(H_hat, 0.0)  # Ĥ must be non-negative
    R_hat = d * (d - 1) * H_hat ** 2
    return R_hat


def strategy_s2_loo_ols(obs: np.ndarray, H: np.ndarray, d: int) -> tuple[np.ndarray, float]:
    """S2: Leave-one-out cross-validated OLS.
    Returns (R_hat_array, std_of_slopes_across_folds)."""
    n = len(obs)
    R_hat = np.zeros(n)
    slopes = []
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        slope, intercept = ols_fit(np.abs(obs[mask]), H[mask])
        slopes.append(slope)
        H_hat_i = slope * np.abs(obs[i]) + intercept
        H_hat_i = max(H_hat_i, 0.0)
        R_hat[i] = d * (d - 1) * H_hat_i ** 2
    return R_hat, float(np.std(slopes))


def strategy_s3_rank_preserving(obs: np.ndarray, H: np.ndarray, d: int) -> np.ndarray:
    """S3: Rank-preserving map. Assign Ĥ = H[rank(obs)] then square.
    This is the theoretical best case: if rank correlation is perfect,
    R̂ = d(d-1)·H²[rank] gives perfect R̂."""
    # Map obs to H via rank matching
    obs_ranks = sp_stats.rankdata(np.abs(obs))
    h_sorted = np.sort(H)
    # obs is anti-correlated with H for b1_std, so reverse
    rho_sign, _ = sp_stats.spearmanr(obs, H)
    if rho_sign < 0:
        h_sorted = h_sorted[::-1]
    # Map each rank to corresponding H value
    n = len(obs)
    H_hat = np.zeros(n)
    for i in range(n):
        idx = int(obs_ranks[i]) - 1  # 0-indexed
        idx = min(idx, n - 1)
        H_hat[i] = h_sorted[idx]
    R_hat = d * (d - 1) * H_hat ** 2
    return R_hat


def strategy_s4_group_mean(obs: np.ndarray, H: np.ndarray, d: int,
                            H_levels: np.ndarray) -> np.ndarray:
    """S4: Group-mean calibration. Average obs per H level, fit on group
    means, apply to individuals."""
    # Compute group means
    group_obs = []
    group_h = []
    for h_val in H_levels:
        mask = np.isclose(H, h_val, atol=1e-6)
        if np.sum(mask) >= 2:
            group_obs.append(np.mean(np.abs(obs[mask])))
            group_h.append(h_val)
    group_obs = np.array(group_obs)
    group_h = np.array(group_h)

    if len(group_obs) < 2:
        return np.full(len(obs), np.nan)

    slope, intercept = ols_fit(group_obs, group_h)
    H_hat = slope * np.abs(obs) + intercept
    H_hat = np.maximum(H_hat, 0.0)
    R_hat = d * (d - 1) * H_hat ** 2
    return R_hat


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=16)
    parser.add_argument("--seed", type=int, default=20260325)
    args = parser.parse_args()

    dims = [4]
    ns_by_d = {4: [128, 256, 512, 1024]}
    hubbles = [0.0, 0.25, 0.5, 1.0, 2.0]
    hubbles_nonzero = [h for h in hubbles if h > 0]
    reps = args.reps

    total = sum(len(ns_by_d[d]) * len(hubbles) * reps for d in dims)

    out_dir = Path(__file__).parent / "outputs_unified_functional"
    out_dir.mkdir(exist_ok=True)

    print(f"T3: Two-Step Analytical Squaring")
    print(f"  Dims: {dims}")
    print(f"  H values: {hubbles}")
    print(f"  Reps: {reps}")
    print(f"  Total: {total} realizations")
    print()

    rows: list[T3Row] = []
    t0 = time.time()
    count = 0

    for d in dims:
        ns = ns_by_d[d]
        for N in ns:
            t_n = time.time()
            for h in hubbles:
                for rep in range(reps):
                    row = run_single(d, N, h, rep, args.seed)
                    rows.append(row)
                    count += 1
                    if count % 20 == 0:
                        elapsed = time.time() - t0
                        eta = elapsed / count * (total - count)
                        print(f"  [{count}/{total}] d={d} N={N} H={h} rep={rep} "
                              f"({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining)")
            dt_n = time.time() - t_n
            print(f"  d={d} N={N} complete ({dt_n:.1f}s)")

    total_time = time.time() - t0
    print(f"\nAll realizations complete ({total_time:.1f}s)")

    # Save CSV
    csv_path = out_dir / "conjecture_e_two_step_squaring.csv"
    field_names = [f.name for f in fields(T3Row)]
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=field_names)
        writer.writeheader()
        for r in rows:
            writer.writerow({fn: getattr(r, fn) for fn in field_names})
    print(f"CSV saved: {csv_path}")

    # ---------------------------------------------------------------------------
    # Analysis
    # ---------------------------------------------------------------------------
    import pandas as pd

    df = pd.DataFrame([{fn: getattr(r, fn) for fn in field_names} for r in rows])

    md_lines: list[str] = []
    md = md_lines.append

    md("# T3: Two-Step Analytical Squaring — Results\n")
    md("")
    md("## Design\n")
    md(f"- Dimensions: {dims}")
    md(f"- N values: { {d: ns_by_d[d] for d in dims} }")
    md(f"- H values: {hubbles}")
    md(f"- Reps per cell: {reps}")
    md(f"- Total realizations: {len(df)}")
    md(f"- Runtime: {total_time:.1f}s")
    md("")
    md("## Method\n")
    md("Two-step bridge: observable → Ĥ (calibration) → R̂ = d(d-1)·Ĥ² (analytical squaring)")
    md("")
    md("Strategies:")
    md("- **S1**: Global OLS: fit |obs| → H on all data, then square")
    md("- **S2**: LOO cross-validated OLS: fit on N-1 points, predict 1, then square")
    md("- **S3**: Rank-preserving: map obs ranks to H values, then square (theoretical upper bound)")
    md("- **S4**: Group-mean: fit on per-H-level averages, apply to individuals, then square")
    md("- **S5**: Same as S1/S2 but using w_max_ratio (antichain channel)")
    md("")

    # -----------------------------------------------------------------------
    # Section 1: Baseline — raw observable correlations
    # -----------------------------------------------------------------------
    md("## 1. Baseline: Raw Observable Performance\n")
    md("| d | N | Observable | ρ(obs, H) | ρ(obs, H²) | R²(obs, H) | R²(obs, H²) | R²(H²)/R²(H) | α_best |")
    md("|---|---|------------|-----------|-------------|------------|--------------|---------------|--------|")

    d = 4
    for N in ns_by_d[d]:
        mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] > 0)
        sub = df[mask]
        if len(sub) < 5:
            continue
        h1 = sub["hubble"].values
        h2 = sub["H2"].values
        R_dS = sub["R_dS"].values

        for obs_name in ["b1_std", "w_max_ratio", "mean_layer_width"]:
            y = sub[obs_name].values
            if np.std(y) < 1e-15 or np.any(np.isnan(y)):
                continue
            rho_h, _ = sp_stats.spearmanr(y, h1)
            rho_h2, _ = sp_stats.spearmanr(y, h2)
            r2_h = ols_r2(h1, y)
            r2_h2 = ols_r2(h2, y)
            ratio = r2_h2 / r2_h if r2_h > 1e-6 else float("inf")
            ag = alpha_grid(h1, y)
            a_best, _ = best_alpha(ag)
            md(f"| {d} | {N} | {obs_name} | {rho_h:+.4f} | {rho_h2:+.4f} | "
               f"{r2_h:.4f} | {r2_h2:.4f} | {ratio:.3f} | {a_best:.2f} |")

    md("")

    # -----------------------------------------------------------------------
    # Section 2: Two-step results — b1_std channel
    # -----------------------------------------------------------------------
    md("## 2. Two-Step Squaring: b1_std Channel\n")
    md("| d | N | Strategy | ρ(R̂, H²) | R²(R̂, H²) | R²(R̂, R_dS) | R²_raw(obs,H²) | Δ R² | α_eff(R̂) | Calibration slope | LOO slope std |")
    md("|---|---|----------|-----------|------------|-------------|-----------------|------|-----------|-------------------|---------------|")

    for N in ns_by_d[d]:
        mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] > 0)
        sub = df[mask]
        if len(sub) < 5:
            continue

        h1 = sub["hubble"].values
        h2 = sub["H2"].values
        R_dS = sub["R_dS"].values
        obs = sub["b1_std"].values

        # Raw baseline
        r2_raw_h2 = ols_r2(h2, obs)

        # S1: Global OLS
        R_hat_s1 = strategy_s1_global_ols(obs, h1, d)
        slope_s1, _ = ols_fit(np.abs(obs), h1)

        # S2: LOO OLS
        R_hat_s2, loo_slope_std = strategy_s2_loo_ols(obs, h1, d)

        # S3: Rank-preserving
        R_hat_s3 = strategy_s3_rank_preserving(obs, h1, d)

        # S4: Group-mean
        R_hat_s4 = strategy_s4_group_mean(obs, h1, d, np.array(hubbles_nonzero))

        strategies = [
            ("S1: Global OLS", R_hat_s1, f"{slope_s1:.4f}", "—"),
            ("S2: LOO OLS", R_hat_s2, "—", f"{loo_slope_std:.4f}"),
            ("S3: Rank-preserving", R_hat_s3, "—", "—"),
            ("S4: Group-mean", R_hat_s4, "—", "—"),
        ]

        for sname, R_hat, slope_str, loo_str in strategies:
            if np.any(np.isnan(R_hat)):
                md(f"| {d} | {N} | {sname} | — | — | — | {r2_raw_h2:.4f} | — | — | {slope_str} | {loo_str} |")
                continue
            rho_s, _ = sp_stats.spearmanr(R_hat, h2)
            r2_rhat_h2 = ols_r2(h2, R_hat)
            r2_rhat_rds = ols_r2(R_dS, R_hat)
            delta_r2 = r2_rhat_h2 - r2_raw_h2
            ag = alpha_grid(h1, R_hat)
            a_best, _ = best_alpha(ag)
            md(f"| {d} | {N} | {sname} | {rho_s:+.4f} | {r2_rhat_h2:.4f} | "
               f"{r2_rhat_rds:.4f} | {r2_raw_h2:.4f} | {delta_r2:+.4f} | "
               f"{a_best:.2f} | {slope_str} | {loo_str} |")

    md("")

    # -----------------------------------------------------------------------
    # Section 3: Two-step results — antichain channel (w_max_ratio)
    # -----------------------------------------------------------------------
    md("## 3. Two-Step Squaring: Antichain Channel (w_max_ratio)\n")
    md("| d | N | Strategy | ρ(R̂, H²) | R²(R̂, H²) | R²(R̂, R_dS) | R²_raw(obs,H²) | Δ R² | α_eff(R̂) |")
    md("|---|---|----------|-----------|------------|-------------|-----------------|------|-----------|")

    for N in ns_by_d[d]:
        mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] > 0)
        sub = df[mask]
        if len(sub) < 5:
            continue

        h1 = sub["hubble"].values
        h2 = sub["H2"].values
        R_dS = sub["R_dS"].values
        obs = sub["w_max_ratio"].values

        if np.any(np.isnan(obs)):
            md(f"| {d} | {N} | S5: Antichain | — | — | — | — | — | — |")
            continue

        r2_raw_h2 = ols_r2(h2, obs)

        # S5a: Global OLS on w_max_ratio
        R_hat_s5a = strategy_s1_global_ols(obs, h1, d)

        # S5b: LOO
        R_hat_s5b, _ = strategy_s2_loo_ols(obs, h1, d)

        # S5c: Group-mean
        R_hat_s5c = strategy_s4_group_mean(obs, h1, d, np.array(hubbles_nonzero))

        for sname, R_hat in [("S5a: Global OLS", R_hat_s5a),
                              ("S5b: LOO OLS", R_hat_s5b),
                              ("S5c: Group-mean", R_hat_s5c)]:
            if np.any(np.isnan(R_hat)):
                md(f"| {d} | {N} | {sname} | — | — | — | {r2_raw_h2:.4f} | — | — |")
                continue
            rho_s, _ = sp_stats.spearmanr(R_hat, h2)
            r2_rhat_h2 = ols_r2(h2, R_hat)
            r2_rhat_rds = ols_r2(R_dS, R_hat)
            delta_r2 = r2_rhat_h2 - r2_raw_h2
            ag = alpha_grid(h1, R_hat)
            a_best, _ = best_alpha(ag)
            md(f"| {d} | {N} | {sname} | {rho_s:+.4f} | {r2_rhat_h2:.4f} | "
               f"{r2_rhat_rds:.4f} | {r2_raw_h2:.4f} | {delta_r2:+.4f} | "
               f"{a_best:.2f} |")

    md("")

    # -----------------------------------------------------------------------
    # Section 4: N-scaling of R²(R̂, H²) improvement
    # -----------------------------------------------------------------------
    md("## 4. N-Scaling: Does Two-Step Squaring Improve with N?\n")
    md("| Observable | Strategy | N=128 R²(R̂,H²) | N=256 | N=512 | N=1024 | ρ(N, R²) | Trend |")
    md("|------------|----------|-----------------|-------|-------|--------|----------|-------|")

    for obs_name in ["b1_std", "w_max_ratio"]:
        for sname_label, strat_func in [("Global OLS", "s1"), ("LOO OLS", "s2"), ("Group-mean", "s4")]:
            r2_by_n = {}
            for N in ns_by_d[d]:
                mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] > 0)
                sub = df[mask]
                if len(sub) < 5:
                    r2_by_n[N] = float("nan")
                    continue

                h1 = sub["hubble"].values
                h2 = sub["H2"].values
                obs = sub[obs_name].values

                if np.any(np.isnan(obs)) or np.std(obs) < 1e-15:
                    r2_by_n[N] = float("nan")
                    continue

                if strat_func == "s1":
                    R_hat = strategy_s1_global_ols(obs, h1, d)
                elif strat_func == "s2":
                    R_hat, _ = strategy_s2_loo_ols(obs, h1, d)
                elif strat_func == "s4":
                    R_hat = strategy_s4_group_mean(obs, h1, d, np.array(hubbles_nonzero))
                else:
                    R_hat = np.full(len(obs), np.nan)

                if np.any(np.isnan(R_hat)):
                    r2_by_n[N] = float("nan")
                    continue

                r2_by_n[N] = ols_r2(h2, R_hat)

            ns_valid = [N for N in ns_by_d[d] if not np.isnan(r2_by_n.get(N, float("nan")))]
            vals_valid = [r2_by_n[N] for N in ns_valid]

            if len(ns_valid) >= 3:
                rho_trend, _ = sp_stats.spearmanr(ns_valid, vals_valid)
                trend_str = "↑" if rho_trend > 0.5 else ("↓" if rho_trend < -0.5 else "→")
            else:
                rho_trend = float("nan")
                trend_str = "N/A"

            cols = []
            for N in [128, 256, 512, 1024]:
                v = r2_by_n.get(N, float("nan"))
                cols.append(f"{v:.4f}" if not np.isnan(v) else "N/A")

            rho_str = f"{rho_trend:+.2f}" if not np.isnan(rho_trend) else "N/A"
            md(f"| {obs_name} | {sname_label} | {cols[0]} | {cols[1]} | {cols[2]} | {cols[3]} | "
               f"{rho_str} | {trend_str} |")

    md("")

    # -----------------------------------------------------------------------
    # Section 5: Δ R² improvement over raw — the critical test
    # -----------------------------------------------------------------------
    md("## 5. Critical Test: Δ R² = R²(R̂, H²) − R²(raw, H²)\n")
    md("Positive Δ R² means the two-step squaring HELPS.\n")
    md("| Observable | Strategy | N=128 Δ R² | N=256 | N=512 | N=1024 | Mean Δ R² | Verdict |")
    md("|------------|----------|-----------|-------|-------|--------|-----------|---------|")

    for obs_name in ["b1_std", "w_max_ratio"]:
        for sname_label, strat_func in [("Global OLS", "s1"), ("LOO OLS", "s2"), ("Group-mean", "s4")]:
            delta_by_n = {}
            for N in ns_by_d[d]:
                mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] > 0)
                sub = df[mask]
                if len(sub) < 5:
                    delta_by_n[N] = float("nan")
                    continue

                h1 = sub["hubble"].values
                h2 = sub["H2"].values
                obs = sub[obs_name].values

                if np.any(np.isnan(obs)) or np.std(obs) < 1e-15:
                    delta_by_n[N] = float("nan")
                    continue

                r2_raw = ols_r2(h2, obs)

                if strat_func == "s1":
                    R_hat = strategy_s1_global_ols(obs, h1, d)
                elif strat_func == "s2":
                    R_hat, _ = strategy_s2_loo_ols(obs, h1, d)
                elif strat_func == "s4":
                    R_hat = strategy_s4_group_mean(obs, h1, d, np.array(hubbles_nonzero))
                else:
                    R_hat = np.full(len(obs), np.nan)

                if np.any(np.isnan(R_hat)):
                    delta_by_n[N] = float("nan")
                    continue

                r2_rhat = ols_r2(h2, R_hat)
                delta_by_n[N] = r2_rhat - r2_raw

            cols = []
            for N in [128, 256, 512, 1024]:
                v = delta_by_n.get(N, float("nan"))
                cols.append(f"{v:+.4f}" if not np.isnan(v) else "N/A")

            valid_deltas = [delta_by_n[N] for N in ns_by_d[d] if not np.isnan(delta_by_n.get(N, float("nan")))]
            mean_delta = np.mean(valid_deltas) if valid_deltas else float("nan")
            if np.isnan(mean_delta):
                verdict = "N/A"
            elif mean_delta > 0.01:
                verdict = "✅ HELPS"
            elif mean_delta > -0.01:
                verdict = "→ NEUTRAL"
            else:
                verdict = "❌ HURTS"

            mean_str = f"{mean_delta:+.4f}" if not np.isnan(mean_delta) else "N/A"
            md(f"| {obs_name} | {sname_label} | {cols[0]} | {cols[1]} | {cols[2]} | {cols[3]} | "
               f"{mean_str} | {verdict} |")

    md("")

    # -----------------------------------------------------------------------
    # Section 6: R̂/R_dS calibration accuracy
    # -----------------------------------------------------------------------
    md("## 6. Calibration Accuracy: R̂/R_dS Ratio\n")
    md("Does R̂ → R_dS = d(d-1)H²? If ratio → 1, calibration is accurate.\n")
    md("| N | H | R_dS | mean(R̂_S1) | mean(R̂_S2) | ratio_S1 | std_S1 | ratio_S2 | std_S2 |")
    md("|---|---|------|-----------|-----------|----------|--------|----------|--------|")

    for N in ns_by_d[d]:
        # Need full slice for calibration
        full_mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] > 0)
        full_sub = df[full_mask]
        if len(full_sub) < 5:
            continue
        h1_full = full_sub["hubble"].values
        obs_full = full_sub["b1_std"].values

        for h_val in hubbles_nonzero:
            R_dS = d * (d - 1) * h_val ** 2
            h_mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] == h_val)
            h_sub = df[h_mask]
            if len(h_sub) < 2:
                continue

            # S1: Global calibration on full slice
            R_hat_s1_full = strategy_s1_global_ols(obs_full, h1_full, d)
            # Extract indices corresponding to this H
            h_indices = full_sub.index[full_sub["hubble"] == h_val]
            # Recompute: simpler to just apply to this H's obs
            slope_s1, intercept_s1 = ols_fit(np.abs(obs_full), h1_full)
            obs_h = h_sub["b1_std"].values
            H_hat_s1 = slope_s1 * np.abs(obs_h) + intercept_s1
            H_hat_s1 = np.maximum(H_hat_s1, 0.0)
            R_hat_s1_h = d * (d - 1) * H_hat_s1 ** 2

            # S2: LOO on full slice, extract this H's predictions
            R_hat_s2_full, _ = strategy_s2_loo_ols(obs_full, h1_full, d)
            # Find positions of this H in full_sub
            pos_mask = (full_sub["hubble"].values == h_val)
            R_hat_s2_h = R_hat_s2_full[pos_mask]

            if R_dS > 1e-10:
                ratio_s1 = R_hat_s1_h / R_dS
                ratio_s2 = R_hat_s2_h / R_dS
                md(f"| {N} | {h_val} | {R_dS:.4f} | {np.mean(R_hat_s1_h):.4f} | "
                   f"{np.mean(R_hat_s2_h):.4f} | {np.mean(ratio_s1):.4f} | "
                   f"{np.std(ratio_s1):.4f} | {np.mean(ratio_s2):.4f} | {np.std(ratio_s2):.4f} |")

    md("")

    # -----------------------------------------------------------------------
    # Section 7: Verdict
    # -----------------------------------------------------------------------
    md("## 7. Verdict\n")
    md("### Key Questions:\n")
    md("1. Does two-step squaring (obs → Ĥ → R̂ = d(d-1)Ĥ²) improve R² against H²?")
    md("2. Does Δ R² increase with N?")
    md("3. Does R̂/R_dS converge to a constant (ideally 1)?")
    md("4. Is the improvement from LOO vs global calibration significant?")
    md("5. Which channel benefits more: spectral (b1_std) or antichain (w_max_ratio)?")
    md("")
    md("*Analysis to be completed after reviewing numerical results.*")
    md("")
    md("---\n")
    md(f"*Generated by conjecture_e_two_step_squaring.py*")
    md(f"*{len(df)} realizations (d=4, N=128/256/512/1024, H=0–2, {reps} reps), runtime {total_time:.1f}s*")

    # Save report
    md_path = out_dir / "conjecture_e_two_step_squaring.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(md_lines))
    print(f"Report saved: {md_path}")


if __name__ == "__main__":
    main()
