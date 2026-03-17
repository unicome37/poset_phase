#!/usr/bin/env python3
"""Prediction C — N-Scaling Law for d(log_H)/dk.

How does the entropy-layer slope scale with poset size N?

Analytical prediction for complete layered posets:
    log_H = k * log((N/k)!)
    d(log_H)/dk ≈ -(N/k + 1/2) ← for fixed k, linear in N

This predicts: at fixed k, |slope| ∝ N (linear scaling).

We test this with:
1. Exact computation at N = 10, 12, 14, 16, 18 (maybe 20)
2. Power-law fit: slope = a * N^alpha
3. Comparison with analytical prediction
4. Extrapolation to large N
"""
from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from generators import generate_random_layered
from entropy_exact import count_linear_extensions_exact
from matched_residual_freedom_check import layer_index_by_minima

OUT_DIR = Path("outputs_exploratory/prediction_c_n_scaling")
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_VALUES = [10, 12, 14, 16, 18, 20]  # 20 may be slow
K_VALUES = [2, 3, 4, 5, 6]
REPS_PER_CELL = 50
ADJACENT_P = 0.35
SKIP_P = 0.08
SEED_BASE = 20260318


def log_H(poset) -> float:
    return math.log(count_linear_extensions_exact(poset, prefer_c=False))


def actual_layer_count(poset) -> int:
    idx = layer_index_by_minima(poset.closure)
    return int(idx.max() + 1) if len(idx) > 0 else 0


def log_H_complete(N: int, k: int) -> float:
    if N % k != 0:
        m = N // k
        r = N % k
        sizes = [m + 1] * r + [m] * (k - r)
        return sum(math.lgamma(s + 1) for s in sizes)
    m = N // k
    return k * math.lgamma(m + 1)


def run():
    rows = []
    t0 = time.time()

    for N in N_VALUES:
        for k in K_VALUES:
            if k > N:
                continue
            for rep in range(REPS_PER_CELL):
                seed = SEED_BASE + N * 10000 + k * 1000 + rep
                poset = generate_random_layered(
                    n=N, n_layers=k,
                    imbalance_mode="uniform",
                    adjacent_p=ADJACENT_P,
                    skip_p=SKIP_P,
                    seed=seed,
                )
                lh = log_H(poset)
                ak = actual_layer_count(poset)
                rows.append({
                    "N": N, "k_nominal": k, "k_actual": ak,
                    "log_H": lh, "seed": seed,
                })
        elapsed = time.time() - t0
        print(f"  N={N}: done ({elapsed:.0f}s)")

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "n_scaling_raw.csv", index=False)
    print(f"\nTotal samples: {len(df)}")

    # ── Analysis ────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("N-SCALING LAW: d(log_H)/dk as function of N")
    print("=" * 70)

    # 1. Per-N slope on actual_k
    slope_rows = []
    for N in sorted(df["N"].unique()):
        sub = df[df["N"] == N]
        slope, intercept, r, p_val, se = stats.linregress(sub["k_actual"], sub["log_H"])
        slope_rows.append({
            "N": N, "slope": slope, "intercept": intercept,
            "r": r, "p": p_val, "se": se,
            "n_obs": len(sub),
        })
        print(f"  N={N}: slope={slope:.4f}, r={r:.3f}, se={se:.4f}")

    slope_df = pd.DataFrame(slope_rows)
    slope_df.to_csv(OUT_DIR / "n_scaling_slopes.csv", index=False)

    # 2. Power-law fit: |slope| = a * N^alpha
    Ns = np.array([r["N"] for r in slope_rows])
    abs_slopes = np.array([abs(r["slope"]) for r in slope_rows])
    log_Ns = np.log(Ns)
    log_slopes = np.log(abs_slopes)

    alpha_slope, alpha_intercept, alpha_r, alpha_p, alpha_se = stats.linregress(log_Ns, log_slopes)
    a_coef = math.exp(alpha_intercept)
    print(f"\n  Power-law fit: |slope| = {a_coef:.4f} * N^{alpha_slope:.3f}")
    print(f"    r^2 = {alpha_r**2:.4f}, p = {alpha_p:.2e}")

    # 3. Analytical prediction for complete layered
    print("\n--- Analytical slopes for complete layered (p=1) ---")
    ana_slope_rows = []
    for N in N_VALUES:
        # Compute analytical slope using finite differences
        # d(log_H)/dk at k=4 (midpoint)
        k_center = 4
        h = 0.5
        lh_plus = log_H_complete(N, k_center + 1)
        lh_minus = log_H_complete(N, k_center - 1)
        ana_slope = (lh_plus - lh_minus) / 2
        ana_slope_rows.append({"N": N, "ana_slope": ana_slope})
        print(f"  N={N}: analytical slope at k=4: {ana_slope:.4f}")

    # 4. Log-log regression for analytical slopes
    ana_Ns = np.array([r["N"] for r in ana_slope_rows])
    ana_abs_slopes = np.array([abs(r["ana_slope"]) for r in ana_slope_rows])
    ana_alpha, ana_int, _, _, _ = stats.linregress(np.log(ana_Ns), np.log(ana_abs_slopes))
    print(f"\n  Analytical power-law: |slope| ∝ N^{ana_alpha:.3f}")

    # 5. Extrapolation
    print("\n--- Extrapolation ---")
    for N_ext in [30, 50, 100, 200]:
        predicted_slope = a_coef * N_ext ** alpha_slope
        ana_slope_ext = -(N_ext / 4 + 0.5)  # Stirling approximation at k=4
        print(f"  N={N_ext}: predicted |slope| = {predicted_slope:.2f} "
              f"(analytical bound = {abs(ana_slope_ext):.2f})")

    # 6. Plot N-scaling
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Linear plot
    ax1.plot(Ns, abs_slopes, "o-", color="C0", label="Empirical (sparse)")
    ana_line = np.array([abs(r["ana_slope"]) for r in ana_slope_rows])
    ax1.plot(Ns, ana_line, "s--", color="C1", label="Analytical (complete)")
    ax1.set_xlabel("N (poset size)")
    ax1.set_ylabel("|slope| = |d(log_H)/dk|")
    ax1.set_title("N-Scaling of Entropy-Layer Slope")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Log-log plot
    ax2.plot(np.log(Ns), np.log(abs_slopes), "o", color="C0")
    x_fit = np.linspace(np.log(Ns.min()), np.log(Ns.max()) * 1.3, 50)
    ax2.plot(x_fit, alpha_slope * x_fit + alpha_intercept, "--", color="C0",
             label=f"fit: $\\alpha$ = {alpha_slope:.2f}")
    ax2.plot(np.log(Ns), np.log(ana_line), "s", color="C1")
    ax2.set_xlabel("log(N)")
    ax2.set_ylabel("log(|slope|)")
    ax2.set_title("Log-Log: Power Law Fit")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(OUT_DIR / "n_scaling_law.png", dpi=150)
    print(f"\n  Saved: n_scaling_law.png")

    # 7. Summary
    print("\n" + "=" * 70)
    print(f"  N-SCALING LAW: |slope| = {a_coef:.4f} × N^{alpha_slope:.3f}")
    print(f"    (R² = {alpha_r**2:.4f})")
    if abs(alpha_slope - 1.0) < 0.3:
        print(f"  ✓ Approximately LINEAR scaling (α ≈ {alpha_slope:.2f} ≈ 1)")
        print("    Consistent with analytical prediction: slope ∝ N/k")
    elif alpha_slope > 1.0:
        print(f"  ↑ Super-linear (α = {alpha_slope:.2f} > 1)")
    else:
        print(f"  ↓ Sub-linear (α = {alpha_slope:.2f} < 1)")
    print("=" * 70)


if __name__ == "__main__":
    run()
