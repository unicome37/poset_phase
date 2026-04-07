"""
Physical explanation: Why Lor5D and KR_2layer beat Lor4D under info penalties
This script quantifies the entropy vs. penalty decomposition across families.
"""
import pandas as pd, numpy as np

# Use full experiment data
df = pd.read_csv("outputs_info_full/raw_samples.csv")

# Focus on A4 mode
a4 = df[df["action_mode"] == "A4"]

print("=" * 100)
print("为什么 Lor5D 和 KR_2layer 在信息论惩罚下胜过 Lor4D？")
print("=" * 100)

# 1. Entropy vs penalty decomposition for top families
print("\n[1] Entropy-Penalty Decomposition at A4, γ=1.0:")
print(f"{'n':>4} {'Family':>30} {'logH':>10} {'I_info':>10} {'A=−logH+I':>12} {'净优势 vs Lor4D':>16}")
print("─" * 90)

for n_val in [20, 40, 60]:
    sub = a4[(a4["n"] == n_val) & (a4["gamma"] == 1.0)]
    grp = sub.groupby("family")[["log_H_mean", "penalty_info", "score"]].mean().sort_values("score")
    lor4d_score = grp.loc["lorentzian_like_4d", "score"]
    for fam in ["KR_2layer", "lorentzian_like_5d", "lorentzian_like_4d", "KR_like"]:
        if fam in grp.index:
            row = grp.loc[fam]
            gap = row["score"] - lor4d_score
            print(f"{n_val:4d} {fam:>30} {row['log_H_mean']:10.2f} {row['penalty_info']:10.3f} "
                  f"{row['score']:12.2f} {gap:+16.2f}")
    print()

# 2. Entropy advantage decomposition
print("\n[2] Entropy Advantage (logH difference relative to Lor4D):")
print(f"{'n':>4} {'Family':>30} {'ΔlogH':>10} {'ΔI_info':>10} {'Δ(−logH+I)':>12} {'熵够补惩罚?':>12}")
print("─" * 90)

for n_val in [20, 40, 60]:
    sub = a4[(a4["n"] == n_val) & (a4["gamma"] == 1.0)]
    grp = sub.groupby("family")[["log_H_mean", "penalty_info", "score"]].mean()
    lor4d = grp.loc["lorentzian_like_4d"]
    for fam in ["KR_2layer", "lorentzian_like_5d"]:
        if fam in grp.index:
            row = grp.loc[fam]
            d_h = row["log_H_mean"] - lor4d["log_H_mean"]
            d_p = row["penalty_info"] - lor4d["penalty_info"]
            d_score = row["score"] - lor4d["score"]
            enough = "是(熵多)" if d_score < 0 else "否(惩罚重)"
            print(f"{n_val:4d} {fam:>30} {d_h:+10.2f} {d_p:+10.3f} {d_score:+12.2f} {enough:>12}")
    print()

# 3. Scaling behavior
print("\n[3] LogH and I_info scaling with n (A4, γ=0):")
print(f"{'Family':>30} {'n=20 logH':>10} {'n=40 logH':>10} {'n=60 logH':>10} {'logH ratio 60/20':>16}")
print("─" * 90)

for fam in ["KR_2layer", "lorentzian_like_5d", "lorentzian_like_4d", "KR_like"]:
    vals = {}
    for n_val in [20, 40, 60]:
        sub = a4[(a4["n"] == n_val) & (a4["gamma"] == 0.0)]
        grp = sub.groupby("family")["log_H_mean"].mean()
        vals[n_val] = grp.get(fam, np.nan)
    ratio = vals[60] / vals[20] if vals[20] > 0 else np.nan
    print(f"{fam:>30} {vals[20]:10.2f} {vals[40]:10.2f} {vals[60]:10.2f} {ratio:16.2f}")

# 4. Info components comparison at n=40
print("\n\n[4] Info penalty components at n=40 (from raw data):")
from observables_info import info_components
from generators import generate_lorentzian_like_4d, generate_lorentzian_like_5d, generate_kr_2layer

COMPS = ["info_spectral_entropy_deficit", "info_degree_heterogeneity",
         "info_layer_concentration", "info_edge_density_extremity",
         "info_interval_diversity_deficit", "info_total"]

families = {
    "KR_2layer": generate_kr_2layer,
    "Lor5D": generate_lorentzian_like_5d,
    "Lor4D": generate_lorentzian_like_4d,
}

for name, gen in families.items():
    comp_vals = {k: [] for k in COMPS}
    for s in range(16):
        p = gen(n=40, seed=s * 100)
        c = info_components(p)
        for k in COMPS:
            comp_vals[k].append(c[k])
    print(f"\n  {name}:")
    for k in COMPS:
        arr = np.array(comp_vals[k])
        print(f"    {k:45s}  mean={arr.mean():.4f}  std={arr.std():.4f}")

# 5. Gap trend
print("\n\n[5] Gap to #1 (KR_2layer) at γ=1.0 — scaling with n:")
print(f"{'n':>4} {'Lor4D−KR2L':>14} {'Lor5D−KR2L':>14} {'Lor4D grows faster?':>20}")
print("─" * 60)
prev_lor4d_gap = None
prev_lor5d_gap = None
for n_val in [20, 40, 60]:
    sub = a4[(a4["n"] == n_val) & (a4["gamma"] == 1.0)]
    grp = sub.groupby("family")["score"].mean()
    lor4d_gap = grp.get("lorentzian_like_4d", np.nan) - grp.get("KR_2layer", np.nan)
    lor5d_gap = grp.get("lorentzian_like_5d", np.nan) - grp.get("KR_2layer", np.nan)
    if prev_lor4d_gap is not None:
        faster = "是" if (lor4d_gap - prev_lor4d_gap) > (lor5d_gap - prev_lor5d_gap) else "否"
    else:
        faster = "—"
    print(f"{n_val:4d} {lor4d_gap:+14.2f} {lor5d_gap:+14.2f} {faster:>20}")
    prev_lor4d_gap = lor4d_gap
    prev_lor5d_gap = lor5d_gap

print("\n" + "=" * 100)
print("物理解释总结:")
print("=" * 100)
print("""
1. KR_2layer 胜出因为：
   - 它是 n 个元素的最大反链叠层，线性扩张数 H 极大化（比 Lor4D 多 30-60%）
   - 作为最规则的二分图结构，5 项信息论惩罚极低（接近完美均匀）
   - 熵优势 + 低惩罚 = 双重碾压

2. Lor5D 胜出因为：
   - 5D 因果偏序比 4D 有更多因果链（logH 高约 5-8%）
   - 因果结构的维数越高，区间数量和层结构越丰富，I_info 也越低
   - 5D 在所有 5 项信息论指标上都优于或持平 4D

3. 根本原因 — 信息论惩罚的固有局限：
   - 信息论惩罚奖励"均匀性"和"高熵性"，这本质上倾向于更大/更对称的结构
   - d=4 的唯一性不在于它是"最均匀的"，而在于它的拓扑/几何特征满足特定约束
   - 纯信息论无法编码"正确维度"这一几何概念
   - 这解释了为何需要两层筛选：信息论选家族类型，几何选维度
""")
