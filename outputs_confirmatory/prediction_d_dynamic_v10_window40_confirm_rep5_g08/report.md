# Prediction D v8 Refine

- v6 strata CSV: `outputs_exploratory/prediction_d_dynamic_v6/cg_blockperm_delta_rank_stratified.csv`
- zeta-rank CSV: `outputs_confirmatory/prediction_d_dynamic_v9_confirm_rep5/cg_zeta_scan_rankings_variants.csv`
- require_variants: `['full', 'switch', 'no_switch']`
- alpha (v6 selection): `0.05`
- require_positive: `False`
- explicit keys: `['40:0.6:0.8']`
- selected (n,keep_ratio,gamma) strata: `1`
- refined n_perm: `200000`

Outputs:
- `cg_blockperm_delta_rank_stratified_refined.csv`

Best (lowest p) strata (refined):
- `full` n=40 keep=0.60 gamma=0.8: rho=0.626, p=5e-06
- `no_switch` n=40 keep=0.60 gamma=0.8: rho=0.818, p=5e-06
- `switch` n=40 keep=0.60 gamma=0.8: rho=0.688, p=5e-06

Worst (highest p) strata (refined):
- `full` n=40 keep=0.60 gamma=0.8: rho=0.626, p=5e-06
- `no_switch` n=40 keep=0.60 gamma=0.8: rho=0.818, p=5e-06
- `switch` n=40 keep=0.60 gamma=0.8: rho=0.688, p=5e-06
