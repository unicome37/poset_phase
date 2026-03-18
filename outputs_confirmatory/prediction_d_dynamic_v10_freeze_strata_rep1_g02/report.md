# Prediction D v8 Refine

- v6 strata CSV: `outputs_exploratory/prediction_d_dynamic_v6/cg_blockperm_delta_rank_stratified.csv`
- zeta-rank CSV: `outputs_confirmatory/prediction_d_dynamic_v9_confirm/cg_zeta_scan_rankings_variants.csv`
- require_variants: `['full', 'switch', 'no_switch']`
- alpha (v6 selection): `0.05`
- require_positive: `False`
- explicit keys: `['30:0.6:0.2']`
- selected (n,keep_ratio,gamma) strata: `1`
- refined n_perm: `200000`

Outputs:
- `cg_blockperm_delta_rank_stratified_refined.csv`

Best (lowest p) strata (refined):
- `full` n=30 keep=0.60 gamma=0.2: rho=0.553, p=2e-05
- `switch` n=30 keep=0.60 gamma=0.2: rho=0.486, p=0.00018
- `no_switch` n=30 keep=0.60 gamma=0.2: rho=0.262, p=0.05183

Worst (highest p) strata (refined):
- `no_switch` n=30 keep=0.60 gamma=0.2: rho=0.262, p=0.05183
- `switch` n=30 keep=0.60 gamma=0.2: rho=0.486, p=0.00018
- `full` n=30 keep=0.60 gamma=0.2: rho=0.553, p=2e-05
