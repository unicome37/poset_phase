"""
Re-analyze full experiment with per-(n,gamma,mode) normalization
to check whether the flagship n=40 gamma=1.0 result was a normalization artifact.
"""
import pandas as pd, numpy as np

# --- 1. Load full experiment raw data ---
df = pd.read_csv("outputs_info_full/raw_samples.csv")
print(f"Full experiment: {len(df)} rows, columns: {list(df.columns)}")

# --- 2. Re-normalize by (n, gamma, action_mode) instead of [n] ---
def robust_zscore(group):
    s = group["score"]
    med = s.median()
    mad = (s - med).abs().median() * 1.4826
    if mad < 1e-12:
        mad = 1.0
    group = group.copy()
    group["z_fair"] = (s - med) / mad
    return group

df_fair = df.groupby(["n", "gamma", "action_mode"], group_keys=False).apply(robust_zscore)

# --- 3. Summary under fair normalization ---
print("\n" + "=" * 90)
print("Re-normalized rankings: group_by=[n, gamma, action_mode]")
print("=" * 90)

for mode in ["A4"]:
    for n_val in [10, 20, 40, 60]:
        print(f"\n{'─'*90}\nn = {n_val}, mode = {mode}")
        for g in [0.0, 0.2, 0.5, 1.0]:
            sub = df_fair[(df_fair["n"] == n_val) & (df_fair["gamma"] == g) & (df_fair["action_mode"] == mode)]
            if sub.empty:
                continue
            grp = sub.groupby("family")["z_fair"].agg(["mean", "std"]).sort_values("mean")
            ranking = list(grp.index)
            lor_rank = ranking.index("lorentzian_like_4d") + 1 if "lorentzian_like_4d" in ranking else -1
            lor_z = grp.loc["lorentzian_like_4d", "mean"] if "lorentzian_like_4d" in grp.index else np.nan
            kr_z = grp.loc["KR_like", "mean"] if "KR_like" in grp.index else np.nan
            top3 = ", ".join([f"{ranking[i]}({grp.iloc[i]['mean']:.3f})" for i in range(min(3, len(ranking)))])
            print(f"  γ={g:.1f}: Lor4D rank #{lor_rank} (z={lor_z:+.3f}), "
                  f"KR_like (z={kr_z:+.3f}), gap={lor_z-kr_z:+.4f}")
            print(f"         Top 3: {top3}")

# --- 4. Compare with original normalization ---
print("\n\n" + "=" * 90)
print("Comparison: original [n]-only normalization vs fair [n,g,mode] normalization")
print("=" * 90)

print(f"\n{'n':>3} {'gamma':>5} {'mode':>4} | {'Orig rank':>9} {'Orig z':>8} | {'Fair rank':>9} {'Fair z':>8}")
print("─" * 65)

for n_val in [10, 20, 40, 60]:
    for g in [0.0, 0.5, 1.0]:
        sub_orig = df[(df["n"] == n_val) & (df["gamma"] == g) & (df["action_mode"] == "A4")]
        sub_fair = df_fair[(df_fair["n"] == n_val) & (df_fair["gamma"] == g) & (df_fair["action_mode"] == "A4")]
        
        if sub_orig.empty:
            continue
            
        # Original
        if "score_norm" in sub_orig.columns:
            grp_o = sub_orig.groupby("family")["score_norm"].mean().sort_values()
            lor_rank_o = list(grp_o.index).index("lorentzian_like_4d") + 1
            lor_z_o = grp_o["lorentzian_like_4d"]
        else:
            lor_rank_o = -1
            lor_z_o = np.nan
        
        # Fair
        grp_f = sub_fair.groupby("family")["z_fair"].mean().sort_values()
        lor_rank_f = list(grp_f.index).index("lorentzian_like_4d") + 1
        lor_z_f = grp_f["lorentzian_like_4d"]
        
        print(f"{n_val:3d} {g:5.1f} {'A4':>4} | #{lor_rank_o:8d} {lor_z_o:+8.3f} | #{lor_rank_f:8d} {lor_z_f:+8.3f}")
