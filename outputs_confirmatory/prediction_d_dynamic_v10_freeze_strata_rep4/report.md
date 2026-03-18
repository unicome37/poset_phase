# Prediction D v8 Refine

- v6 strata CSV: `outputs_exploratory/prediction_d_dynamic_v6/cg_blockperm_delta_rank_stratified.csv`
- zeta-rank CSV: `outputs_confirmatory/prediction_d_dynamic_v9_confirm_rep4/cg_zeta_scan_rankings_variants.csv`
- require_variants: `['full', 'switch', 'no_switch']`
- alpha (v6 selection): `0.05`
- require_positive: `False`
- explicit keys: `['30:0.6:0.8']`
- selected (n,keep_ratio,gamma) strata: `1`
- refined n_perm: `200000`

Outputs:
- `cg_blockperm_delta_rank_stratified_refined.csv`

Best (lowest p) strata (refined):
- `full` n=30 keep=0.60 gamma=0.8: rho=0.299, p=0.02566
- `no_switch` n=30 keep=0.60 gamma=0.8: rho=-0.003, p=0.9693
- `switch` n=30 keep=0.60 gamma=0.8: rho=0.003, p=0.9704

Worst (highest p) strata (refined):
- `switch` n=30 keep=0.60 gamma=0.8: rho=0.003, p=0.9704
- `no_switch` n=30 keep=0.60 gamma=0.8: rho=-0.003, p=0.9693
- `full` n=30 keep=0.60 gamma=0.8: rho=0.299, p=0.02566
