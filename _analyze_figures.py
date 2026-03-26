"""Analyze control group figures data in detail."""
import pandas as pd
import numpy as np

df = pd.read_csv("outputs/raw_samples.csv")

# A2, N=80
sub = df[(df["action_mode"] == "A2") & (df["n"] == 80)]
agg = sub.groupby(["gamma", "family"]).agg(
    mean_score_norm=("score_norm", "mean"),
    std_score_norm=("score_norm", "std"),
    mean_logH=("log_H_mean", "mean"),
    mean_geo=("penalty_geometric", "mean"),
    mean_neutral=("penalty_neutral", "mean"),
).reset_index()

CTRL = {"KR_2layer", "KR_4layer"}
LOR = {"lorentzian_like_2d", "lorentzian_like_3d", "lorentzian_like_4d", "lorentzian_like_5d"}

print("=" * 90)
print("TABLE 1: A2, N=80 — Score ranking per gamma")
print("=" * 90)
for g in sorted(agg["gamma"].unique()):
    sg = agg[agg["gamma"] == g].sort_values("mean_score_norm", ascending=False)
    print(f"\n--- gamma = {g} ---")
    for rank, (_, r) in enumerate(sg.iterrows(), 1):
        fam = r["family"]
        tag = " <<<CTRL>>>" if fam in CTRL else (" [Lor]" if fam in LOR else "")
        print(f"  {rank:2d}. {fam:35s}  S_norm={r['mean_score_norm']:+.3f}  "
              f"logH={r['mean_logH']:.1f}  geo={r['mean_geo']:.1f}{tag}")

print("\n" + "=" * 90)
print("TABLE 2: Cross-N comparison for controls and Lor2D (A2, gamma=0.2)")
print("=" * 90)
focus = ["KR_2layer", "KR_4layer", "KR_like", "lorentzian_like_2d", "lorentzian_like_4d"]
sub2 = df[(df["action_mode"] == "A2") & (df["gamma"] == 0.2)]
agg2 = sub2.groupby(["n", "family"]).agg(
    mean_score_norm=("score_norm", "mean"),
    mean_logH=("log_H_mean", "mean"),
    mean_geo=("penalty_geometric", "mean"),
).reset_index()
agg2 = agg2[agg2["family"].isin(focus)]
pivot = agg2.pivot(index="family", columns="n", values="mean_score_norm")
print(pivot.to_string(float_format=lambda x: f"{x:+.3f}"))

print("\n" + "=" * 90)
print("TABLE 3: KR_2layer logH vs geometric penalty across N (gamma=0)")
print("=" * 90)
sub3 = df[(df["action_mode"] == "A2") & (df["gamma"] == 0.0) & (df["family"] == "KR_2layer")]
agg3 = sub3.groupby("n").agg(
    mean_logH=("log_H_mean", "mean"),
    std_logH=("log_H_mean", "std"),
    mean_geo=("penalty_geometric", "mean"),
    std_geo=("penalty_geometric", "std"),
    mean_neutral=("penalty_neutral", "mean"),
).reset_index()
print(agg3.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

print("\n" + "=" * 90)
print("TABLE 4: gamma sensitivity — Rank of KR_2layer across gamma values")
print("=" * 90)
for n_val in [20, 40, 60, 80]:
    sub4 = df[(df["action_mode"] == "A2") & (df["n"] == n_val)]
    agg4 = sub4.groupby(["gamma", "family"])["score_norm"].mean().reset_index()
    ranks = []
    for g in sorted(agg4["gamma"].unique()):
        sg = agg4[agg4["gamma"] == g].sort_values("score_norm", ascending=False).reset_index(drop=True)
        r = sg[sg["family"] == "KR_2layer"].index[0] + 1
        ranks.append(f"g={g}: #{r}")
    print(f"  N={n_val:3d}  " + "  ".join(ranks))
