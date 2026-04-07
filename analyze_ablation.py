"""消融实验深入分析"""
import pandas as pd
import numpy as np

df = pd.read_csv('outputs_info_ablation/ablation_raw.csv')

# 1) interval_diversity 对各族惩罚的贡献
print("=" * 80)
print("§1  interval_diversity_deficit 对 penalty_info 的贡献 (γ=1.0)")
print("=" * 80)

base = df[(df['dropped_term']=='none')&(df['gamma']==1.0)].groupby(['family','n'])['penalty_info'].mean()
drop = df[(df['dropped_term']=='info_interval_diversity_deficit')&(df['gamma']==1.0)].groupby(['family','n'])['penalty_info'].mean()
contrib = (base - drop).unstack('n')
base_un = base.unstack('n')

for n_val in sorted(contrib.columns):
    print(f"\n--- n={n_val} ---")
    ranked = contrib[n_val].sort_values(ascending=False)
    for fam, val in ranked.items():
        total = base_un.loc[fam, n_val]
        pct = val / total * 100 if total > 0 else 0
        tag = " <<< Lor4D" if "4d" in str(fam) else ""
        print(f"  {str(fam):40s} interval_contrib={val:7.3f}  total={total:7.3f}  pct={pct:5.1f}%{tag}")

# 2) 各信息论分量的相对杀伤力（对 Lor4D）
print("\n" + "=" * 80)
print("§2  各信息论分量对 Lor4D 的杀伤力 (Δrank at γ=1.0)")
print("=" * 80)

summary = pd.read_csv('outputs_info_ablation/ablation_summary.csv')
baseline_ranks = summary[(summary['gamma']==1.0)&(summary['dropped_term']=='none')].set_index('n')['lor4d_rank']

for term in sorted(summary['dropped_term'].unique()):
    if term == 'none':
        continue
    term_ranks = summary[(summary['gamma']==1.0)&(summary['dropped_term']==term)].set_index('n')['lor4d_rank']
    changes = baseline_ranks - term_ranks
    parts = [f"n={n}: {changes.get(n, 0):+.0f}" for n in [10, 20, 40]]
    print(f"  Drop {term:45s}  " + "  ".join(parts))

# 3) γ 临界阈值分析：Lor4D 何时从 top-3 跌出
print("\n" + "=" * 80)
print("§3  Lor4D rank vs γ (baseline, all terms)")
print("=" * 80)

for n_val in [10, 20, 40]:
    ranks = summary[(summary['dropped_term']=='none')&(summary['n']==n_val)].sort_values('gamma')
    print(f"\n  n={n_val}:")
    for _, row in ranks.iterrows():
        bar = "#" * int(17 - row['lor4d_rank'])
        print(f"    γ={row['gamma']:.1f}: rank #{row['lor4d_rank']:2.0f}/17  {bar}  top1={row['top1_family']}")

# 4) 互相关：各分量的冗余度
print("\n" + "=" * 80)
print("§4  分量间冗余 (penalty_info correlation across drops)")
print("=" * 80)

terms = [t for t in sorted(summary['dropped_term'].unique()) if t != 'none']
# Use n=20 gamma=1.0 raw data for correlation
data_sets = {}
for term in terms:
    sub = df[(df['dropped_term']==term)&(df['n']==20)&(df['gamma']==1.0)]
    data_sets[term] = sub.groupby('family')['penalty_info'].mean()

corr_df = pd.DataFrame(data_sets).corr()
print(corr_df.to_string(float_format=lambda x: f"{x:.3f}"))

# 5) 作用量分解：entropy vs penalty at key conditions
print("\n" + "=" * 80)
print("§5  Lor4D 作用量分解 (n=20)")
print("=" * 80)

lor4d_base = df[(df['dropped_term']=='none')&(df['family']=='lorentzian_like_4d')&(df['n']==20)]
for gamma_val in [0.0, 0.2, 0.5, 1.0]:
    sub = lor4d_base[lor4d_base['gamma']==gamma_val]
    logH = sub['log_H'].mean()
    penalty = sub['penalty_info'].mean()
    score = sub['score'].mean()
    print(f"  γ={gamma_val:.1f}: logH={logH:.3f}  penalty_info={penalty:.3f}  score={score:.3f}  (penalty/|logH|)={penalty/abs(logH):.2%}")
