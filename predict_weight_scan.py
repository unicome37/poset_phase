"""基于消融数据推算权重扫描的理论预期"""
import pandas as pd
import numpy as np

# Load ablation raw data
df = pd.read_csv('outputs_info_ablation/ablation_raw.csv')

# Strategy: from ablation data, we can reconstruct each family's 
# per-component penalty contribution, then compute what happens
# under different weight configurations

# For each family at n=20, compute the implied per-component penalties:
#   penalty(all) - penalty(drop_X) = contribution of X
# Then total penalty under custom weights = sum(weight_i * contribution_i / original_weight_i)

terms = [
    'info_spectral_entropy_deficit',
    'info_degree_heterogeneity', 
    'info_layer_concentration',
    'info_edge_density_extremity',
    'info_interval_diversity_deficit',
]

for n_val in [20, 40]:
    print(f"\n{'='*80}")
    print(f"n={n_val}: Reconstructed per-component penalties")
    print(f"{'='*80}")
    
    # Get baseline penalty_info for each family
    base = df[(df['dropped_term']=='none') & (df['n']==n_val) & (df['gamma']==1.0)]
    base_means = base.groupby('family')['penalty_info'].mean()
    
    # Get penalty when each term is dropped
    contribs = {}
    for term in terms:
        dropped = df[(df['dropped_term']==term) & (df['n']==n_val) & (df['gamma']==1.0)]
        drop_means = dropped.groupby('family')['penalty_info'].mean()
        contribs[term] = (base_means - drop_means)
    
    contrib_df = pd.DataFrame(contribs)
    
    # Per-unit contribution (divide by weight=5.0)
    per_unit = contrib_df / 5.0
    
    print("\nPer-unit contributions (penalty / weight=5.0):")
    for fam in ['lorentzian_like_4d', 'KR_2layer', 'lorentzian_like_5d', 'lorentzian_like_3d', 'KR_like']:
        if fam in per_unit.index:
            vals = per_unit.loc[fam]
            parts = [f"{t.replace('info_','')[:8]}={v:.4f}" for t,v in vals.items()]
            print(f"  {fam:35s}  {', '.join(parts)}")
    
    # Simulate custom weight configurations
    configs = {
        'baseline_5x5': {t: 5.0 for t in terms},
        'interval_0': {**{t: 5.0 for t in terms}, 'info_interval_diversity_deficit': 0.0},
        'interval_1': {**{t: 5.0 for t in terms}, 'info_interval_diversity_deficit': 1.0},
        'interval_2': {**{t: 5.0 for t in terms}, 'info_interval_diversity_deficit': 2.0},
        'no_int_boost_2x': {**{t: 10.0 for t in terms}, 'info_interval_diversity_deficit': 0.0},
        'int_1_rest_1.5x': {**{t: 7.5 for t in terms}, 'info_interval_diversity_deficit': 1.0},
        'int_0.5_rest_2x': {**{t: 10.0 for t in terms}, 'info_interval_diversity_deficit': 0.5},
    }
    
    print(f"\nSimulated rankings under different weight configs (n={n_val}, gamma=1.0):")
    
    # Get logH for each family
    logH = base.groupby('family')['log_H_mean'].mean()
    
    for config_name, weights in configs.items():
        # Reconstruct penalty under new weights
        new_penalty = sum(per_unit[t] * weights[t] for t in terms)
        # Score = -logH + penalty (at gamma=1.0, beta=1.0)
        new_score = -logH + new_penalty
        ranked = new_score.sort_values()
        families_ranked = list(ranked.index)
        lor4d_rank = families_ranked.index('lorentzian_like_4d') + 1
        lor4d_score = ranked['lorentzian_like_4d']
        top3 = [(f, ranked[f]) for f in families_ranked[:3]]
        
        top3_str = ', '.join([f"{f.replace('lorentzian_like_','lor')}({s:.1f})" for f,s in top3])
        print(f"  {config_name:25s}: Lor4D=#{lor4d_rank:2d}  top3=[{top3_str}]")

# Also check gamma=0.5
print(f"\n{'='*80}")
print(f"Simulated rankings at gamma=0.5 (n=20)")
print(f"{'='*80}")

base05 = df[(df['dropped_term']=='none') & (df['n']==20) & (df['gamma']==0.5)]
base05_means = base05.groupby('family')['penalty_info'].mean()
logH_20 = df[(df['dropped_term']=='none') & (df['n']==20) & (df['gamma']==0.0)].groupby('family')['log_H_mean'].mean()

# At gamma=0.5: score = -logH + 0.5*penalty
for config_name, weights in configs.items():
    # Approximate: scale baseline penalty by weight ratio
    # penalty_new = sum(per_unit_i * new_weight_i) for n=20
    per_unit_20 = pd.DataFrame({t: (
        df[(df['dropped_term']=='none')&(df['n']==20)&(df['gamma']==1.0)].groupby('family')['penalty_info'].mean() -
        df[(df['dropped_term']==t)&(df['n']==20)&(df['gamma']==1.0)].groupby('family')['penalty_info'].mean()
    ) / 5.0 for t in terms})
    
    new_penalty = sum(per_unit_20[t] * weights[t] for t in terms)
    new_score = -logH_20 + 0.5 * new_penalty
    ranked = new_score.sort_values()
    families_ranked = list(ranked.index)
    lor4d_rank = families_ranked.index('lorentzian_like_4d') + 1
    
    top3 = list(families_ranked[:3])
    top3_short = [f.replace('lorentzian_like_', 'lor') for f in top3]
    print(f"  {config_name:25s}: Lor4D=#{lor4d_rank:2d}  top3={top3_short}")
