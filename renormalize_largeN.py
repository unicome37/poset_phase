"""Fair renormalization for large-N experiment."""
import pandas as pd, numpy as np

df = pd.read_csv("outputs_info_largeN/raw_samples.csv")

def robust_zscore(group):
    s = group["score"]
    med = s.median()
    mad = (s - med).abs().median() * 1.4826
    if mad < 1e-12:
        mad = 1.0
    group = group.copy()
    group["z_fair"] = (s - med) / mad
    return group

df = df.groupby(["n", "gamma", "action_mode"], group_keys=False).apply(robust_zscore)

print("Large-N fair normalization [n, gamma, mode]:")
print("=" * 90)
for n_val in [60, 80, 100]:
    print(f"\n  n = {n_val}")
    for g in [0.0, 0.2, 0.5, 1.0]:
        sub = df[(df["n"] == n_val) & (df["gamma"] == g) & (df["action_mode"] == "A4")]
        if sub.empty:
            continue
        grp = sub.groupby("family")["z_fair"].mean().sort_values()
        ranking = list(grp.index)
        lor_rank = ranking.index("lorentzian_like_4d") + 1
        lor_z = grp["lorentzian_like_4d"]
        top5 = []
        for i in range(min(5, len(ranking))):
            fam = ranking[i]
            tag = " **" if "lorentzian" in fam and "4d" in fam else ""
            top5.append(f"{fam}({grp.iloc[i]:.3f}){tag}")
        print(f"    γ={g:.1f}: Lor4D #{lor_rank} (z={lor_z:+.3f})")
        print(f"           Top 5: {', '.join(top5)}")

# Trend table for gamma=1.0
print(f"\n{'='*90}")
print(f"Lor4D trend at γ=1.0 (fair norm):")
print(f"{'n':>5} {'rank':>6} {'z_fair':>8} {'gap to #1':>10} {'#1 family':>30}")
for n_val in [60, 80, 100]:
    sub = df[(df["n"] == n_val) & (df["gamma"] == 1.0) & (df["action_mode"] == "A4")]
    if sub.empty:
        continue
    grp = sub.groupby("family")["z_fair"].mean().sort_values()
    ranking = list(grp.index)
    lor_rank = ranking.index("lorentzian_like_4d") + 1
    lor_z = grp["lorentzian_like_4d"]
    top_z = grp.iloc[0]
    top_fam = ranking[0]
    print(f"{n_val:5d} #{lor_rank:5d} {lor_z:+8.3f} {lor_z - top_z:+10.3f} {top_fam:>30}")
