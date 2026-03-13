"""Simpson's Paradox 深度分析"""
import pandas as pd
import numpy as np

df = pd.read_csv('outputs_exploratory/prediction_c_comprehensive/tier1_all_family_raw.csv')

# 1. N as confound
print('=== Within-family: N vs HII and N vs log_H ===')
for fam, sub in df.groupby('family'):
    r_n_hii = np.corrcoef(sub['n'], sub['hierarchy_integration_index'])[0,1]
    r_n_logh = np.corrcoef(sub['n'], sub['log_H'])[0,1]
    r_hii_logh = np.corrcoef(sub['hierarchy_integration_index'], sub['log_H'])[0,1]
    print(f'{fam:30s}  r(N,HII)={r_n_hii:+.4f}  r(N,logH)={r_n_logh:+.4f}  r(HII,logH)={r_hii_logh:+.4f}')

# 2. Fixed-N analysis
print('\n=== Fixed-N: cross-family HII vs log_H ===')
for n_val in sorted(df['n'].unique()):
    subn = df[df['n'] == n_val]
    means = subn.groupby('family').agg(
        hii=('hierarchy_integration_index', 'mean'),
        logh=('log_H', 'mean')
    ).reset_index()
    r_cross = np.corrcoef(means['hii'], means['logh'])[0, 1]
    r_all = np.corrcoef(subn['hierarchy_integration_index'], subn['log_H'])[0, 1]
    print(f'N={int(n_val):2d}: cross-family_means r={r_cross:+.4f}, all_samples r={r_all:+.4f} (n={len(subn)})')

# 3. Per-family fixed-N
print('\n=== Within N=14, per-family HII vs log_H ===')
sub14 = df[df['n'] == 14]
for fam, s in sub14.groupby('family'):
    if len(s) < 5:
        continue
    r = np.corrcoef(s['hierarchy_integration_index'], s['log_H'])[0, 1]
    hii_m = s['hierarchy_integration_index'].mean()
    lh_m = s['log_H'].mean()
    print(f'{fam:30s}  n={len(s):3d}  r={r:+.4f}  hii={hii_m:+.3f}  log_H={lh_m:.2f}')

# 4. Partial correlation controlling for N and family
print('\n=== Partial corr: HII vs log_H | N + family dummies ===')
dummies = pd.get_dummies(df['family'], prefix='fam', drop_first=True).astype(float)
X = np.column_stack([np.ones(len(df)), df['n'].values, dummies.values])
y = df['log_H'].values
x_hii = df['hierarchy_integration_index'].values

# Residualize
from numpy.linalg import lstsq
b_y, *_ = lstsq(X, y, rcond=None)
b_x, *_ = lstsq(X, x_hii, rcond=None)
y_res = y - X @ b_y
x_res = x_hii - X @ b_x
r_partial = np.corrcoef(x_res, y_res)[0, 1]
print(f'partial_r(HII, log_H | N + family) = {r_partial:+.4f}')

# 5. Component decomposition at fixed N
print('\n=== Component-level at fixed N (pooled N=14) ===')
sub14 = df[df['n'] == 14].copy()
cols = ['layer_count', 'mean_layer_gap', 'long_edge_fraction',
        'adjacent_edge_fraction', 'reduction_edge_density',
        'cover_density', 'layer_signature_redundancy']
for c in cols:
    r = np.corrcoef(sub14[c], sub14['log_H'])[0, 1]
    print(f'{c:35s}  r(vs log_H) = {r:+.4f}')

# 6. The key insight: layer_count correlation with N
print('\n=== Layer count scaling with N by family ===')
for fam, sub in df.groupby('family'):
    means = sub.groupby('n')['layer_count'].mean()
    vals = [f'N={int(n)}:{v:.1f}' for n, v in means.items()]
    print(f'{fam:30s}  {", ".join(vals)}')

# 7. HII z-score decomposition
print('\n=== HII decomposition: mean z-scores by family ===')
z_cols = ['z_layer_count', 'z_mean_layer_gap', 'z_long_edge_fraction',
          'z_adjacent_edge_fraction', 'z_reduction_edge_density']
for fam, sub in df.groupby('family'):
    vals = {c: sub[c].mean() for c in z_cols}
    parts = [f'{c.replace("z_",""):20s}={v:+.3f}' for c, v in vals.items()]
    print(f'{fam:30s}  HII={sub["hierarchy_integration_index"].mean():+.3f}  | {" | ".join(parts)}')
