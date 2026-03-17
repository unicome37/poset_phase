"""
Prediction A — Full BDG Action vs Link-Proxy Comparison

GPT / literature review identified a critical distinction:
- Our S_BD_d2 = N - 2*C_0 is the "link action" / d=2 BD action
- The standard d=4 BDG action from Machet-Wang (2020) / Carlip (2024) is:
    S^(4) = N - C_0 + 9*C_1 - 16*C_2 + 8*C_3
  where C_k = intervals with k interior elements

This script re-analyzes our existing data with the CORRECT BDG coefficients
to determine whether the 4D selection window survives under the full action.

Data source: outputs_exploratory/prediction_a_bd_extended/base_observables.csv
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd

DATA_DIR = Path("outputs_exploratory/prediction_a_bd_extended")
OUT_DIR = Path("outputs_exploratory/prediction_a_bdg_full")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Load existing data
base = pd.read_csv(DATA_DIR / "base_observables.csv")

N_VALUES = sorted(base["n"].unique())
BETA = 1.0

# ===== Define Action Variants =====
# Each variant computes S/N from the interval counts

def s_link_d2(row) -> float:
    """Link action / d=2 BD: S = N - 2*C_0"""
    return (row["n"] - 2 * row["C0_links"]) / row["n"]

def s_bdg_d4_standard(row) -> float:
    """Standard BDG d=4 from Machet-Wang / Carlip:
    S^(4) = N - C_0 + 9*C_1 - 16*C_2 + 8*C_3
    """
    s = row["n"] - row["C0_links"] + 9*row["C1"] - 16*row["C2"] + 8*row["C3"]
    return s / row["n"]

def s_bdg_d4_variant_a(row) -> float:
    """Alternative: drop the N term (pure interval contribution):
    S = -C_0 + 9*C_1 - 16*C_2 + 8*C_3
    """
    s = -row["C0_links"] + 9*row["C1"] - 16*row["C2"] + 8*row["C3"]
    return s / row["n"]

def s_bdg_d2_intervals(row) -> float:
    """d=2 BDG with interval correction:
    S^(2,full) = N - 2*C_0 + 2*C_1
    (first correction term from BDG expansion)
    """
    s = row["n"] - 2*row["C0_links"] + 2*row["C1"]
    return s / row["n"]

ACTION_VARIANTS = {
    "link_d2":          s_link_d2,           # Our original: N - 2*C_0
    "BDG_d4_standard":  s_bdg_d4_standard,   # Literature standard: N - C_0 + 9C_1 - 16C_2 + 8C_3
    "BDG_d4_pure":      s_bdg_d4_variant_a,  # Pure interval: -C_0 + 9C_1 - 16C_2 + 8C_3
    "BDG_d2_corrected": s_bdg_d2_intervals,  # d=2 with first correction: N - 2C_0 + 2C_1
}

# ===== Compute Normalized Action for Each Variant =====
print("=" * 76)
print("FULL BDG ACTION vs LINK-PROXY COMPARISON")
print("=" * 76)

# Add action columns
for vname, func in ACTION_VARIANTS.items():
    base[f"S_{vname}_norm"] = base.apply(func, axis=1)

# Print action profiles per dimension
lor_families = [f for f in base["family"].unique() if f.startswith("lorentzian_like_")]

print("\n--- Normalized Action Values by Dimension & Variant ---")
for n in N_VALUES:
    print(f"\n  N={n}:")
    for fam in sorted(lor_families) + ["KR_like"]:
        sub = base[(base["n"] == n) & (base["family"] == fam)]
        label = fam.replace("lorentzian_like_", "").replace("KR_like", "KR")
        parts = [f"{label:>4}:"]
        for vname in ACTION_VARIANTS:
            val = sub[f"S_{vname}_norm"].mean()
            parts.append(f"{vname}={val:+.3f}")
        parts.append(f"logH={sub['log_H'].mean():.2f}")
        print(f"    {' | '.join(parts)}")

# ===== Scoring & Winner Analysis =====
LAMBDA_VALUES = [0.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0]

# Filter Lorentzian only
lor = base[base["family"].str.startswith("lorentzian_like_")]

score_rows = []
for _, row in lor.iterrows():
    for vname in ACTION_VARIANTS:
        s_norm = row[f"S_{vname}_norm"]
        for lbd in LAMBDA_VALUES:
            score = -BETA * row["log_H"] + lbd * s_norm
            score_rows.append({
                "n": row["n"], "family": row["family"],
                "sample_id": row["sample_id"],
                "variant": vname, "lambda": lbd,
                "score": score,
            })

scores = pd.DataFrame(score_rows)

summary = scores.groupby(["n", "lambda", "variant", "family"]).agg(
    mean_score=("score", "mean"),
).reset_index()

# Winner per (n, λ, variant)
winner_rows = []
for (n, lbd, vname), grp in summary.groupby(["n", "lambda", "variant"]):
    ordered = grp.sort_values("mean_score")
    best = ordered.iloc[0]
    second = ordered.iloc[1]
    margin = float(second["mean_score"] - best["mean_score"])
    winner_rows.append({
        "n": int(n), "lambda": lbd, "variant": vname,
        "winner": best["family"].replace("lorentzian_like_", ""),
        "margin": margin,
        "runner_up": second["family"].replace("lorentzian_like_", ""),
    })

winners = pd.DataFrame(winner_rows)
winners.to_csv(OUT_DIR / "winners_bdg_comparison.csv", index=False)

# ===== Print Winner Tables =====
print("\n" + "=" * 76)
print("WINNER TABLES: link_d2 vs BDG_d4_standard")
print("=" * 76)

for vname in ["link_d2", "BDG_d4_standard"]:
    print(f"\n  {vname}:")
    w = winners[winners["variant"] == vname]
    print(f"  {'λ':>5}", end="")
    for n in N_VALUES:
        print(f"  {n:>4}", end="")
    print()
    for lbd in LAMBDA_VALUES:
        print(f"  {lbd:5.1f}", end="")
        for n in N_VALUES:
            row = w[(w["n"] == n) & (w["lambda"] == lbd)]
            if not row.empty:
                print(f"  {row.iloc[0]['winner']:>4}", end="")
            else:
                print(f"  {'?':>4}", end="")
        print()

# All variants side by side
print("\n" + "=" * 76)
print("4D WIN COUNT BY (VARIANT, λ)")
print("=" * 76)

print(f"\n  {'λ':>5}", end="")
for vname in ACTION_VARIANTS:
    print(f"  {vname:>16}", end="")
print()

for lbd in LAMBDA_VALUES:
    print(f"  {lbd:5.1f}", end="")
    for vname in ACTION_VARIANTS:
        w = winners[(winners["variant"] == vname) & (winners["lambda"] == lbd)]
        n4d = (w["winner"] == "4d").sum()
        total = len(w)
        label = f"{n4d}/{total}"
        if n4d == total:
            label += "★"
        print(f"  {label:>16}", end="")
    print()

# ===== Critical Comparison: link_d2 vs BDG_d4_standard =====
print("\n" + "=" * 76)
print("CRITICAL COMPARISON: Does 4D plateau survive under full BDG d=4?")
print("=" * 76)

for lbd in [6.0, 7.0, 8.0]:
    print(f"\n  λ = {lbd}:")
    for vname in ["link_d2", "BDG_d4_standard"]:
        w = winners[(winners["variant"] == vname) & (winners["lambda"] == lbd)]
        n4d = (w["winner"] == "4d").sum()
        print(f"    {vname:>20}: 4D wins {n4d}/{len(w)}")
        for _, row in w.sort_values("n").iterrows():
            marker = "★" if row["winner"] == "4d" else " "
            print(f"      N={int(row['n']):3d}: {row['winner']:>4} beats {row['runner_up']:>4} by {row['margin']:.3f} {marker}")

# ===== BDG_d4_standard: find best λ for 4D =====
print("\n--- BDG_d4_standard: Best λ for unanimous 4D ---")
w_bdg = winners[winners["variant"] == "BDG_d4_standard"]
for lbd in LAMBDA_VALUES:
    w = w_bdg[w_bdg["lambda"] == lbd]
    n4d = (w["winner"] == "4d").sum()
    dominant = w["winner"].value_counts().index[0]
    avg_m = w[w["winner"] == "4d"]["margin"].mean() if n4d > 0 else 0
    print(f"  λ={lbd:5.1f}: 4D wins {n4d}/{len(w)}, dominant={dominant}, avg_margin={avg_m:.3f}")

# ===== Score details at key λ for BDG_d4_standard =====
print("\n--- BDG_d4_standard: Full rankings at critical λ ---")
for lbd in [5.0, 6.0, 7.0, 8.0, 10.0]:
    print(f"\n  λ = {lbd}")
    sub = summary[(summary["variant"] == "BDG_d4_standard") & (summary["lambda"] == lbd)]
    for n in N_VALUES:
        s_n = sub[sub["n"] == n].sort_values("mean_score")
        ranking = " > ".join(
            f"{r['family'].replace('lorentzian_like_', '')}({r['mean_score']:.1f})"
            for _, r in s_n.iterrows()
        )
        print(f"    N={n}: {ranking}")

# ===== Key diagnostic: WHY does BDG_d4 differ from link_d2? =====
print("\n" + "=" * 76)
print("DIAGNOSTIC: Action component breakdown")
print("=" * 76)

for n in [44, 52, 68]:
    print(f"\n  N={n}:")
    for fam in sorted(lor_families):
        sub = lor[(lor["n"] == n) & (lor["family"] == fam)]
        label = fam.replace("lorentzian_like_", "")
        c0 = sub["C0_links"].mean()
        c1 = sub["C1"].mean()
        c2 = sub["C2"].mean()
        c3 = sub["C3"].mean()
        # BDG_d4 contributions
        link_term = -c0
        c1_term = 9 * c1
        c2_term = -16 * c2
        c3_term = 8 * c3
        print(f"    {label:>4}: C0={c0:6.1f} C1={c1:6.1f} C2={c2:5.1f} C3={c3:5.1f}"
              f" | -C0={link_term:+8.1f}  +9C1={c1_term:+8.1f}"
              f"  -16C2={c2_term:+8.1f}  +8C3={c3_term:+7.1f}"
              f" | net={link_term + c1_term + c2_term + c3_term:+8.1f}")

print(f"\nAll outputs saved to {OUT_DIR}")
