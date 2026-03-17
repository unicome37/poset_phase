#!/usr/bin/env python3
"""Reanalyze dose-response data using actual layer count instead of nominal k."""
import pandas as pd
from scipy import stats
import numpy as np
from pathlib import Path

OUT_DIR = Path("outputs_exploratory/prediction_c_dose_response")
df = pd.read_csv(OUT_DIR / "dose_response_raw.csv")

print("=== Actual-k Dose-Response Reanalysis ===\n")

# 1. Regression on actual_k
print("--- Per-N regression: log_H ~ k_actual ---")
reg_rows = []
for N in sorted(df["N"].unique()):
    sub = df[df["N"] == N]
    slope, intercept, r, p, se = stats.linregress(sub["k_actual"], sub["log_H"])
    k_min, k_max = int(sub.k_actual.min()), int(sub.k_actual.max())
    print(f"  N={N}: slope={slope:.4f}, r={r:.3f}, p={p:.2e}, k=[{k_min},{k_max}]")
    reg_rows.append({"N": N, "predictor": "k_actual", "slope": slope, "r": r, "p": p})

print("\n--- Per-N regression: log_H ~ k_nominal (for comparison) ---")
for N in sorted(df["N"].unique()):
    sub = df[df["N"] == N]
    slope, intercept, r, p, se = stats.linregress(sub["k_nominal"], sub["log_H"])
    print(f"  N={N}: slope={slope:.4f}, r={r:.3f}, p={p:.2e}")
    reg_rows.append({"N": N, "predictor": "k_nominal", "slope": slope, "r": r, "p": p})

pd.DataFrame(reg_rows).to_csv(OUT_DIR / "dose_response_actual_k_regression.csv", index=False)

# 2. Means by actual_k
print("\n--- Mean log_H by actual_k ---")
for N in sorted(df["N"].unique()):
    sub = df[df["N"] == N]
    print(f"\n  N={N}:")
    grp = sub.groupby("k_actual").agg(
        mean_log_H=("log_H", "mean"),
        std_log_H=("log_H", "std"),
        n=("log_H", "count"),
    ).reset_index()
    for _, row in grp.iterrows():
        ka = int(row.k_actual)
        print(f"    k_actual={ka:2d}: mean_log_H={row.mean_log_H:.3f} +/- {row.std_log_H:.3f} (n={int(row.n)})")

# 3. Spearman on cell means
print("\n--- Spearman rho on cell means (actual_k -> mean_log_H) ---")
for N in sorted(df["N"].unique()):
    sub = df[df["N"] == N]
    grp = sub.groupby("k_actual")["log_H"].mean().reset_index()
    rho, p = stats.spearmanr(grp["k_actual"], grp["log_H"])
    print(f"  N={N}: Spearman rho={rho:.3f}, p={p:.4f}, n_cells={len(grp)}")

# 4. Key finding
print("\n" + "=" * 60)
all_actual = [r for r in reg_rows if r["predictor"] == "k_actual"]
all_nominal = [r for r in reg_rows if r["predictor"] == "k_nominal"]
print("Key: Using actual_k vs nominal_k")
for a, n in zip(all_actual, all_nominal):
    ratio = a["r"] / n["r"] if n["r"] != 0 else float("nan")
    print(f"  N={a['N']}: r(actual)={a['r']:.3f} vs r(nominal)={n['r']:.3f}  (ratio={ratio:.2f}x)")
print("=" * 60)
