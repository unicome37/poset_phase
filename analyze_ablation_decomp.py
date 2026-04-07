"""Lor4D vs KR_2layer 作用量分解对比"""
import pandas as pd

df = pd.read_csv('outputs_info_ablation/ablation_raw.csv')

for fam_name in ['lorentzian_like_4d', 'KR_2layer', 'lorentzian_like_5d']:
    fam = df[(df['dropped_term'] == 'none') & (df['family'] == fam_name)]
    by_ng = fam.groupby(['n', 'gamma']).agg({
        'score': 'mean',
        'penalty_info': 'mean',
        'log_H_mean': 'mean',
    }).reset_index()

    print(f"\n=== {fam_name} ===")
    for _, row in by_ng.iterrows():
        logH = row['log_H_mean']
        p = row['penalty_info']
        s = row['score']
        g = row['gamma']
        ratio = g * p / abs(logH) if logH != 0 else 0
        n = int(row['n'])
        print(f"  n={n:3d} γ={g:.1f}: logH={logH:9.3f}  I_info={p:7.3f}  A={s:10.3f}  γI/|logH|={ratio:.1%}")

# Direct comparison at γ=1.0
print("\n=== γ=1.0 direct competition ===")
baseline = df[(df['dropped_term'] == 'none') & (df['gamma'] == 1.0)]
by_fn = baseline.groupby(['family', 'n']).agg({
    'score': 'mean',
    'penalty_info': 'mean',
    'log_H_mean': 'mean',
}).reset_index()

for n_val in [10, 20, 40]:
    sub = by_fn[by_fn['n'] == n_val].sort_values('score')
    print(f"\n  n={n_val} (sorted by A = -logH + I_info):")
    for _, row in sub.head(5).iterrows():
        fam = row['family']
        tag = " <<" if "4d" in fam else ""
        print(f"    {fam:40s} logH={row['log_H_mean']:9.3f}  I_info={row['penalty_info']:7.3f}  A={row['score']:10.3f}{tag}")
    # Where is Lor4D?
    lor4d_row = sub[sub['family'] == 'lorentzian_like_4d']
    if not lor4d_row.empty:
        idx = list(sub['family']).index('lorentzian_like_4d') + 1
        r = lor4d_row.iloc[0]
        print(f"    ... (Lor4D at #{idx})")
        print(f"    {'lorentzian_like_4d':40s} logH={r['log_H_mean']:9.3f}  I_info={r['penalty_info']:7.3f}  A={r['score']:10.3f} <<")
