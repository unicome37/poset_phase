"""
Prediction B — Extended γ-range analysis
Extracts γ_c from extended scan (γ up to 5.0) and compares with original [0, 2.0] results.
Also tests near-wall N=48,52,56 using SIS.
"""
import numpy as np
import pandas as pd
from pathlib import Path

OUT_DIR = Path("outputs_exploratory/prediction_b_gamma_extended")

raw = pd.read_csv(OUT_DIR / "raw_samples.csv")
summary = pd.read_csv(OUT_DIR / "summary.csv")

print("=== Extended γ scan: lor2d vs KR_like (A2 exact) ===\n")

results = []
for n in sorted(summary["n"].unique()):
    sub = summary[summary["n"] == n]
    gammas = sorted(sub["gamma"].unique())

    for g in gammas:
        slice_g = sub[sub["gamma"] == g]
        lor = slice_g[slice_g["family"] == "lorentzian_like_2d"]["mean_score_norm"].values
        kr = slice_g[slice_g["family"] == "KR_like"]["mean_score_norm"].values
        if len(lor) > 0 and len(kr) > 0:
            results.append({"n": n, "gamma": g, "lor2d_score": lor[0], "kr_score": kr[0],
                           "diff_lor_minus_kr": lor[0] - kr[0],
                           "winner": "lor2d" if lor[0] > kr[0] else "KR"})

rdf = pd.DataFrame(results)

# Find γ_c for each N (crossover point)
print(f"{'N':>4} | {'γ_c (original)':>14} | {'γ_c (extended)':>14} | {'max γ tested':>12} | Winner at γ=5")
print("-" * 75)

for n in sorted(rdf["n"].unique()):
    sub = rdf[rdf["n"] == n].sort_values("gamma")
    diffs = sub["diff_lor_minus_kr"].values
    gammas = sub["gamma"].values

    # Find first sign change (lor2d starts winning at low γ, KR takes over at high γ)
    gc = None
    for i in range(len(diffs) - 1):
        if diffs[i] > 0 and diffs[i+1] <= 0:
            # Linear interpolation
            gc = gammas[i] + (gammas[i+1] - gammas[i]) * diffs[i] / (diffs[i] - diffs[i+1])
            break

    winner_5 = sub[sub["gamma"] == 5.0]["winner"].values[0] if 5.0 in gammas else "?"
    gc_str = f"{gc:.4f}" if gc is not None else "no cross"
    print(f"{n:4d} | {'':>14s} | {gc_str:>14s} | {gammas.max():>12.1f} | {winner_5}")

print()

# Detailed crossover table
rdf.to_csv(OUT_DIR / "gamma_extended_crossover.csv", index=False)

# --- Visualization ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for ax, n in zip(axes.flat, sorted(rdf["n"].unique())):
    sub = rdf[rdf["n"] == n].sort_values("gamma")
    ax.plot(sub["gamma"], sub["lor2d_score"], 'b-o', markersize=4, label="lor2d")
    ax.plot(sub["gamma"], sub["kr_score"], 'r-s', markersize=4, label="KR_like")
    ax.axhline(0, color='gray', linestyle=':', alpha=0.5)

    # Mark crossover
    diffs = sub["diff_lor_minus_kr"].values
    gammas = sub["gamma"].values
    for i in range(len(diffs) - 1):
        if diffs[i] > 0 and diffs[i+1] <= 0:
            gc = gammas[i] + (gammas[i+1] - gammas[i]) * diffs[i] / (diffs[i] - diffs[i+1])
            ax.axvline(gc, color='green', linestyle='--', alpha=0.7, label=f"γ_c={gc:.3f}")
            break

    ax.set_xlabel("γ")
    ax.set_ylabel("Normalized score")
    ax.set_title(f"N = {n}")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

plt.suptitle("Prediction B: Extended γ scan [0, 5.0]", fontsize=14)
plt.tight_layout()
fig.savefig(OUT_DIR / "gamma_extended_crossover.png", dpi=150)
print(f"Saved to {OUT_DIR}")
