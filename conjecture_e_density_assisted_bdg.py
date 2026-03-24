"""Conjecture E — T1: Density-Assisted BDG Calibration Experiment.

Tests whether including the ρ^{2/d} prefactor and calibrating on flat (H=0)
sprinklings improves the BDG action's correlation with R = d(d-1)H².

The BDG theorem states:
    lim_{ρ→∞} ρ^{-2/d} ⟨S_BDG⟩ / N → -½ ∫ R √g d^d x / V

So the properly normalized per-element curvature estimator is:
    R_hat(x) = -2 · ρ^{2/d} · (B1)(x)

where ρ = N/V is the sprinkling density.

Problem: at finite N, b1_mean carries a large density-dependent offset that
overwhelms the curvature signal. T1 tests three correction strategies:

  Strategy A: Raw ρ^{2/d} scaling
    R_hat_A = -2 · ρ_eff^{2/d} · b1_mean

  Strategy B: Flat-baseline subtraction
    R_hat_B = -2 · ρ_eff^{2/d} · (b1_mean - b1_mean_flat)
    where b1_mean_flat = mean(b1_mean) at H=0 for same N

  Strategy C: OLS density residualization + ρ^{2/d} scaling
    R_hat_C = -2 · ρ_eff^{2/d} · resid(b1_mean | n_causal_pairs)

  Strategy D: b1_std (already known to track H) with ρ^{2/d} scaling
    R_hat_D = ρ_eff^{2/d} · b1_std

For each strategy, we measure:
  - Spearman ρ(R_hat, H²) at each N
  - R²(H²) vs R²(H) ratio (does corrected observable prefer H² over H?)
  - N-scaling: does R²(H²) improve with N?
  - α_eff from power-law grid scan

Design:
  - d = 2, 3, 4
  - N = 128, 256, 512 (+ 1024 for d=4 only)
  - H = 0, 0.25, 0.5, 1.0, 2.0
  - 16 reps per cell
  - ρ_eff estimated as n_causal_pairs / N (proxy for sprinkling density)

References:
  - §4.1.27: b1_std beyond-density at d=4 (ρ_resid = -0.703)
  - §4.1.33: No algebraic bridge (0/9)
  - §4.1.37: α_eff(N) monotonically increases but very slowly
  - Benincasa-Dowker (2010): BDG action continuum limit
"""

from __future__ import annotations

import argparse
import csv
import time
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

from curvature_layer2_de_sitter_scan import sprinkle_de_sitter_like_diamond
from curvature_layer2_recovery_fast import build_causal_matrix_de_sitter
from conjecture_e_sorkin_dalembertian import build_dalembertian_matrix


# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class T1Row:
    d: int
    N: int
    hubble: float
    rep: int
    R_dS: float
    H2: float
    n_causal_pairs: int
    rho_eff: float          # n_causal_pairs / N as density proxy
    b1_mean: float
    b1_std: float
    b1_median: float
    # BDG action = N * b1_mean (sum of B1 over all elements)
    s_bdg: float
    # ρ^{2/d} scaled versions
    b1_mean_scaled: float   # ρ_eff^{2/d} * b1_mean
    b1_std_scaled: float    # ρ_eff^{2/d} * b1_std
    s_bdg_scaled: float     # ρ_eff^{2/d} * s_bdg


def run_single(d: int, N: int, hubble: float, rep: int, seed: int) -> T1Row:
    points = sprinkle_de_sitter_like_diamond(
        N, d - 1, hubble=hubble,
        seed=seed + d * 100000 + N * 100 + rep + int(hubble * 1000),
    )
    causal = build_causal_matrix_de_sitter(points, hubble)
    n_causal_pairs = int(np.sum(causal))
    R_dS = d * (d - 1) * hubble ** 2
    H2 = hubble ** 2

    # Density proxy: n_causal_pairs / N
    rho_eff = n_causal_pairs / N if N > 0 else 1.0

    # Build B_ℓ matrix and apply to constant field
    B = build_dalembertian_matrix(causal, d, rho=1.0)
    ones = np.ones(N)
    b1 = B @ ones
    b1_mean = float(np.mean(b1))
    b1_std = float(np.std(b1))
    b1_median = float(np.median(b1))

    # BDG action = sum of B1
    s_bdg = float(np.sum(b1))

    # ρ^{2/d} scaled versions
    rho_factor = rho_eff ** (2.0 / d) if rho_eff > 0 else 1.0
    b1_mean_scaled = rho_factor * b1_mean
    b1_std_scaled = rho_factor * b1_std
    s_bdg_scaled = rho_factor * s_bdg

    return T1Row(
        d=d, N=N, hubble=hubble, rep=rep,
        R_dS=R_dS, H2=H2,
        n_causal_pairs=n_causal_pairs,
        rho_eff=rho_eff,
        b1_mean=b1_mean, b1_std=b1_std, b1_median=b1_median,
        s_bdg=s_bdg,
        b1_mean_scaled=b1_mean_scaled,
        b1_std_scaled=b1_std_scaled,
        s_bdg_scaled=s_bdg_scaled,
    )


# ---------------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------------
def ols_residualize(feat: np.ndarray, dens: np.ndarray) -> np.ndarray:
    if len(feat) < 3 or np.std(dens) < 1e-15:
        return feat - np.mean(feat)
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
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=16)
    parser.add_argument("--seed", type=int, default=20260324)
    parser.add_argument("--d4-only", action="store_true",
                        help="Only run d=4 for speed")
    args = parser.parse_args()

    dims = [4] if args.d4_only else [2, 3, 4]
    ns_by_d = {2: [128, 256, 512], 3: [128, 256, 512], 4: [128, 256, 512, 1024]}
    hubbles = [0.0, 0.25, 0.5, 1.0, 2.0]
    reps = args.reps

    total = sum(len(ns_by_d[d]) * len(hubbles) * reps for d in dims)

    out_dir = Path(__file__).parent / "outputs_unified_functional"
    out_dir.mkdir(exist_ok=True)

    print(f"T1: Density-Assisted BDG Calibration")
    print(f"  Dims: {dims}")
    print(f"  H values: {hubbles}")
    print(f"  Reps: {reps}")
    print(f"  Total: {total} realizations")
    print()

    rows: list[T1Row] = []
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
    csv_path = out_dir / "conjecture_e_density_assisted_bdg.csv"
    field_names = [f.name for f in fields(T1Row)]
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

    md("# T1: Density-Assisted BDG Calibration — Results\n")
    md("")
    md("## Design\n")
    md(f"- Dimensions: {dims}")
    md(f"- N values: { {d: ns_by_d[d] for d in dims} }")
    md(f"- H values: {hubbles}")
    md(f"- Reps per cell: {reps}")
    md(f"- Total realizations: {len(df)}")
    md(f"- Runtime: {total_time:.1f}s")
    md("")

    # -----------------------------------------------------------------------
    # Section 1: Raw correlations (baseline)
    # -----------------------------------------------------------------------
    md("## 1. Baseline: Raw Feature Correlations with H²\n")
    md("| d | N | Feature | Spearman ρ(feat, H²) | p-value | R²(H²) | R²(H) | R²(H²)/R²(H) |")
    md("|---|---|---------|----------------------|---------|---------|-------|---------------|")

    raw_features = ["b1_mean", "b1_std", "b1_mean_scaled", "b1_std_scaled", "s_bdg", "s_bdg_scaled"]

    for d in dims:
        for N in ns_by_d[d]:
            mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] > 0)
            sub = df[mask]
            if len(sub) < 5:
                continue
            h2 = sub["H2"].values
            h1 = sub["hubble"].values
            for feat_name in raw_features:
                y = sub[feat_name].values
                if np.std(y) < 1e-15:
                    continue
                rho_s, p_s = sp_stats.spearmanr(y, h2)
                r2_h2 = ols_r2(h2, y)
                r2_h1 = ols_r2(h1, y)
                ratio = r2_h2 / r2_h1 if r2_h1 > 1e-6 else float("inf")
                md(f"| {d} | {N} | {feat_name} | {rho_s:+.4f} | {p_s:.2e} | "
                   f"{r2_h2:.4f} | {r2_h1:.4f} | {ratio:.3f} |")

    md("")

    # -----------------------------------------------------------------------
    # Section 2: Flat-baseline subtraction (Strategy B)
    # -----------------------------------------------------------------------
    md("## 2. Strategy B: Flat-Baseline Subtraction\n")
    md("Subtract mean(b1_mean) at H=0 for same (d, N), then apply ρ^{2/d}.\n")
    md("| d | N | R_hat_B = -2·ρ^{2/d}·(b1_mean - b1_flat) | Spearman ρ(R_hat_B, H²) | p | R²(H²) | R²(H) | ratio | α_best |")
    md("|---|---|-------------------------------------------|-------------------------|---|---------|-------|-------|--------|")

    for d in dims:
        for N in ns_by_d[d]:
            # Get flat baseline
            flat_mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] == 0.0)
            b1_flat = df.loc[flat_mask, "b1_mean"].mean()

            # Get non-zero H data
            mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] > 0)
            sub = df[mask]
            if len(sub) < 5:
                continue

            h2 = sub["H2"].values
            h1 = sub["hubble"].values
            rho_eff = sub["rho_eff"].values
            b1_mean = sub["b1_mean"].values

            # Strategy B: baseline-corrected, ρ-scaled
            rho_factor = rho_eff ** (2.0 / d)
            R_hat_B = -2.0 * rho_factor * (b1_mean - b1_flat)

            if np.std(R_hat_B) < 1e-15:
                continue

            rho_s, p_s = sp_stats.spearmanr(R_hat_B, h2)
            r2_h2 = ols_r2(h2, R_hat_B)
            r2_h1 = ols_r2(h1, R_hat_B)
            ratio = r2_h2 / r2_h1 if r2_h1 > 1e-6 else float("inf")
            ag = alpha_grid(h1, R_hat_B)
            a_best, r2_best = best_alpha(ag)

            md(f"| {d} | {N} | b1_flat={b1_flat:.4f} | {rho_s:+.4f} | {p_s:.2e} | "
               f"{r2_h2:.4f} | {r2_h1:.4f} | {ratio:.3f} | {a_best:.2f} |")

    md("")

    # -----------------------------------------------------------------------
    # Section 3: OLS density residualization + ρ^{2/d} (Strategy C)
    # -----------------------------------------------------------------------
    md("## 3. Strategy C: OLS Density Residualization + ρ^{2/d}\n")
    md("OLS remove n_causal_pairs from b1_mean, then scale by ρ^{2/d}.\n")
    md("| d | N | Spearman ρ(R_hat_C, H²) | p | R²(H²) | R²(H) | ratio | α_best |")
    md("|---|---|-------------------------|---|---------|-------|-------|--------|")

    for d in dims:
        for N in ns_by_d[d]:
            mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] > 0)
            sub = df[mask]
            if len(sub) < 5:
                continue

            h2 = sub["H2"].values
            h1 = sub["hubble"].values
            rho_eff = sub["rho_eff"].values
            b1_mean = sub["b1_mean"].values
            n_cp = sub["n_causal_pairs"].values.astype(float)

            # Strategy C: OLS residualize, then ρ-scale
            resid = ols_residualize(b1_mean, n_cp)
            rho_factor = rho_eff ** (2.0 / d)
            R_hat_C = -2.0 * rho_factor * resid

            if np.std(R_hat_C) < 1e-15:
                continue

            rho_s, p_s = sp_stats.spearmanr(R_hat_C, h2)
            r2_h2 = ols_r2(h2, R_hat_C)
            r2_h1 = ols_r2(h1, R_hat_C)
            ratio = r2_h2 / r2_h1 if r2_h1 > 1e-6 else float("inf")
            ag = alpha_grid(h1, R_hat_C)
            a_best, r2_best = best_alpha(ag)

            md(f"| {d} | {N} | {rho_s:+.4f} | {p_s:.2e} | "
               f"{r2_h2:.4f} | {r2_h1:.4f} | {ratio:.3f} | {a_best:.2f} |")

    md("")

    # -----------------------------------------------------------------------
    # Section 4: b1_std with ρ^{2/d} (Strategy D) — known best channel
    # -----------------------------------------------------------------------
    md("## 4. Strategy D: b1_std with ρ^{2/d} Scaling\n")
    md("| d | N | Spearman ρ(R_hat_D, H²) | p | R²(H²) | R²(H) | ratio | α_best |")
    md("|---|---|-------------------------|---|---------|-------|-------|--------|")

    for d in dims:
        for N in ns_by_d[d]:
            mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] > 0)
            sub = df[mask]
            if len(sub) < 5:
                continue

            h2 = sub["H2"].values
            h1 = sub["hubble"].values
            R_hat_D = sub["b1_std_scaled"].values

            if np.std(R_hat_D) < 1e-15:
                continue

            rho_s, p_s = sp_stats.spearmanr(R_hat_D, h2)
            r2_h2 = ols_r2(h2, R_hat_D)
            r2_h1 = ols_r2(h1, R_hat_D)
            ratio = r2_h2 / r2_h1 if r2_h1 > 1e-6 else float("inf")
            ag = alpha_grid(h1, R_hat_D)
            a_best, r2_best = best_alpha(ag)

            md(f"| {d} | {N} | {rho_s:+.4f} | {p_s:.2e} | "
               f"{r2_h2:.4f} | {r2_h1:.4f} | {ratio:.3f} | {a_best:.2f} |")

    md("")

    # -----------------------------------------------------------------------
    # Section 5: Strategy E — b1_std density-residualized + ρ^{2/d}
    # -----------------------------------------------------------------------
    md("## 5. Strategy E: b1_std Density-Residualized + ρ^{2/d}\n")
    md("| d | N | Spearman ρ(R_hat_E, H²) | p | R²(H²) | R²(H) | ratio | α_best |")
    md("|---|---|-------------------------|---|---------|-------|-------|--------|")

    for d in dims:
        for N in ns_by_d[d]:
            mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] > 0)
            sub = df[mask]
            if len(sub) < 5:
                continue

            h2 = sub["H2"].values
            h1 = sub["hubble"].values
            rho_eff = sub["rho_eff"].values
            b1_std_vals = sub["b1_std"].values
            n_cp = sub["n_causal_pairs"].values.astype(float)

            resid = ols_residualize(b1_std_vals, n_cp)
            rho_factor = rho_eff ** (2.0 / d)
            R_hat_E = rho_factor * resid

            if np.std(R_hat_E) < 1e-15:
                continue

            rho_s, p_s = sp_stats.spearmanr(R_hat_E, h2)
            r2_h2 = ols_r2(h2, R_hat_E)
            r2_h1 = ols_r2(h1, R_hat_E)
            ratio = r2_h2 / r2_h1 if r2_h1 > 1e-6 else float("inf")
            ag = alpha_grid(h1, R_hat_E)
            a_best, r2_best = best_alpha(ag)

            md(f"| {d} | {N} | {rho_s:+.4f} | {p_s:.2e} | "
               f"{r2_h2:.4f} | {r2_h1:.4f} | {ratio:.3f} | {a_best:.2f} |")

    md("")

    # -----------------------------------------------------------------------
    # Section 6: Cross-strategy comparison at d=4
    # -----------------------------------------------------------------------
    md("## 6. Cross-Strategy Comparison at d=4\n")
    md("| N | Strategy | Feature | ρ(H²) | R²(H²) | R²(H) | R²(H²)/R²(H) | α_best |")
    md("|---|----------|---------|--------|---------|-------|---------------|--------|")

    d = 4
    for N in ns_by_d.get(d, []):
        flat_mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] == 0.0)
        b1_flat = df.loc[flat_mask, "b1_mean"].mean()

        mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] > 0)
        sub = df[mask]
        if len(sub) < 5:
            continue

        h2 = sub["H2"].values
        h1 = sub["hubble"].values
        rho_eff_vals = sub["rho_eff"].values
        rho_factor = rho_eff_vals ** (2.0 / d)
        n_cp = sub["n_causal_pairs"].values.astype(float)

        strategies = {}

        # A: raw b1_mean (no correction)
        strategies["A: raw b1_mean"] = sub["b1_mean"].values

        # A2: raw b1_mean_scaled
        strategies["A2: ρ^{2/d}·b1_mean"] = sub["b1_mean_scaled"].values

        # B: flat-baseline + ρ scaling
        strategies["B: baseline+ρ"] = -2.0 * rho_factor * (sub["b1_mean"].values - b1_flat)

        # C: OLS resid + ρ scaling
        resid_c = ols_residualize(sub["b1_mean"].values, n_cp)
        strategies["C: OLS resid+ρ"] = -2.0 * rho_factor * resid_c

        # D: b1_std scaled
        strategies["D: ρ^{2/d}·b1_std"] = sub["b1_std_scaled"].values

        # E: b1_std resid + ρ
        resid_e = ols_residualize(sub["b1_std"].values, n_cp)
        strategies["E: b1_std resid+ρ"] = rho_factor * resid_e

        # F: b1_std raw (no ρ correction — baseline from T5)
        strategies["F: raw b1_std"] = sub["b1_std"].values

        for name, y_vals in strategies.items():
            if np.std(y_vals) < 1e-15:
                md(f"| {N} | {name} | — | — | — | — | — | — |")
                continue
            rho_s, _ = sp_stats.spearmanr(y_vals, h2)
            r2_h2 = ols_r2(h2, y_vals)
            r2_h1 = ols_r2(h1, y_vals)
            ratio = r2_h2 / r2_h1 if r2_h1 > 1e-6 else float("inf")
            ag = alpha_grid(h1, y_vals)
            a_best_val, _ = best_alpha(ag)
            md(f"| {N} | {name} | — | {rho_s:+.4f} | {r2_h2:.4f} | "
               f"{r2_h1:.4f} | {ratio:.3f} | {a_best_val:.2f} |")

    md("")

    # -----------------------------------------------------------------------
    # Section 7: N-scaling of R²(H²)/R²(H) ratio per strategy
    # -----------------------------------------------------------------------
    md("## 7. N-Scaling: R²(H²)/R²(H) Ratio Trend (d=4)\n")
    md("Does any strategy show faster convergence toward H² preference?\n")
    md("| Strategy | N=128 | N=256 | N=512 | N=1024 | ρ(N, ratio) | Trend |")
    md("|----------|-------|-------|-------|--------|-------------|-------|")

    d = 4
    strategy_names_d4 = ["A: raw b1_mean", "B: baseline+ρ", "C: OLS resid+ρ",
                         "D: ρ^{2/d}·b1_std", "E: b1_std resid+ρ", "F: raw b1_std"]

    for strat_name in strategy_names_d4:
        ratios_by_n = {}
        for N in ns_by_d.get(d, []):
            flat_mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] == 0.0)
            b1_flat = df.loc[flat_mask, "b1_mean"].mean()

            mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] > 0)
            sub = df[mask]
            if len(sub) < 5:
                ratios_by_n[N] = float("nan")
                continue

            h2 = sub["H2"].values
            h1 = sub["hubble"].values
            rho_eff_vals = sub["rho_eff"].values
            rho_factor = rho_eff_vals ** (2.0 / d)
            n_cp = sub["n_causal_pairs"].values.astype(float)

            if strat_name == "A: raw b1_mean":
                y = sub["b1_mean"].values
            elif strat_name == "B: baseline+ρ":
                y = -2.0 * rho_factor * (sub["b1_mean"].values - b1_flat)
            elif strat_name == "C: OLS resid+ρ":
                y = -2.0 * rho_factor * ols_residualize(sub["b1_mean"].values, n_cp)
            elif strat_name == "D: ρ^{2/d}·b1_std":
                y = sub["b1_std_scaled"].values
            elif strat_name == "E: b1_std resid+ρ":
                y = rho_factor * ols_residualize(sub["b1_std"].values, n_cp)
            elif strat_name == "F: raw b1_std":
                y = sub["b1_std"].values
            else:
                y = sub["b1_mean"].values

            if np.std(y) < 1e-15:
                ratios_by_n[N] = float("nan")
                continue

            r2_h2 = ols_r2(h2, y)
            r2_h1 = ols_r2(h1, y)
            ratios_by_n[N] = r2_h2 / r2_h1 if r2_h1 > 1e-6 else float("nan")

        ns_valid = [N for N in ns_by_d.get(d, []) if not np.isnan(ratios_by_n.get(N, float("nan")))]
        vals_valid = [ratios_by_n[N] for N in ns_valid]

        if len(ns_valid) >= 3:
            rho_trend, _ = sp_stats.spearmanr(ns_valid, vals_valid)
            trend_str = "↑ increasing" if rho_trend > 0.5 else ("↓ decreasing" if rho_trend < -0.5 else "→ flat")
        else:
            rho_trend = float("nan")
            trend_str = "N/A"

        cols = []
        for N in [128, 256, 512, 1024]:
            v = ratios_by_n.get(N, float("nan"))
            cols.append(f"{v:.3f}" if not np.isnan(v) else "N/A")

        md(f"| {strat_name} | {cols[0]} | {cols[1]} | {cols[2]} | {cols[3]} | "
           f"{rho_trend:+.2f}" if not np.isnan(rho_trend) else f"| {strat_name} | {cols[0]} | {cols[1]} | {cols[2]} | {cols[3]} | N/A"
           + f" | {trend_str} |")

    md("")

    # -----------------------------------------------------------------------
    # Section 8: R_hat / R_dS convergence check
    # -----------------------------------------------------------------------
    md("## 8. R_hat / R_dS Convergence (Strategy B, d=4)\n")
    md("If R_hat_B → R_dS, then R_hat_B / R_dS → constant.\n")
    md("| N | H | mean(R_hat_B) | R_dS | ratio | std(ratio) |")
    md("|---|---|---------------|------|-------|------------|")

    d = 4
    for N in ns_by_d.get(d, []):
        flat_mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] == 0.0)
        b1_flat = df.loc[flat_mask, "b1_mean"].mean()

        for h in [h for h in hubbles if h > 0]:
            mask = (df["d"] == d) & (df["N"] == N) & (df["hubble"] == h)
            sub = df[mask]
            if len(sub) < 2:
                continue
            rho_eff_vals = sub["rho_eff"].values
            rho_factor = rho_eff_vals ** (2.0 / d)
            R_hat_B = -2.0 * rho_factor * (sub["b1_mean"].values - b1_flat)
            R_dS = d * (d - 1) * h ** 2
            mean_rhat = np.mean(R_hat_B)
            if abs(R_dS) > 1e-10:
                ratio_vals = R_hat_B / R_dS
                md(f"| {N} | {h} | {mean_rhat:.4f} | {R_dS:.4f} | "
                   f"{np.mean(ratio_vals):.4f} | {np.std(ratio_vals):.4f} |")
            else:
                md(f"| {N} | {h} | {mean_rhat:.4f} | {R_dS:.4f} | N/A | N/A |")

    md("")

    # -----------------------------------------------------------------------
    # Section 9: Verdict
    # -----------------------------------------------------------------------
    md("## 9. Verdict\n")
    md("### Key Questions:\n")
    md("1. Does ρ^{2/d} scaling improve H² preference over raw features?")
    md("2. Does flat-baseline subtraction help extract curvature signal?")
    md("3. Does any strategy achieve α ≈ 2 at finite N?")
    md("4. Does R_hat_B / R_dS converge to a constant with N?")
    md("5. Which strategy shows fastest N-scaling of H² preference?")
    md("")
    md("*Analysis to be completed after reviewing numerical results.*")
    md("")
    md("---\n")
    md(f"*Generated by conjecture_e_density_assisted_bdg.py*")
    md(f"*{len(df)} realizations, runtime {total_time:.1f}s*")

    # Save report
    md_path = out_dir / "conjecture_e_density_assisted_bdg.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(md_lines))
    print(f"Report saved: {md_path}")


if __name__ == "__main__":
    main()
