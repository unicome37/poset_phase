# Prediction D v8 Refine

- v6 strata CSV: `outputs_exploratory/prediction_d_dynamic_v6/cg_blockperm_delta_rank_stratified.csv`
- zeta-rank CSV: `outputs_confirmatory/prediction_d_dynamic_v9_confirm/cg_zeta_scan_rankings_variants.csv`
- require_variants: `['full', 'switch', 'no_switch']`
- alpha (v6 selection): `0.05`
- require_positive: `False`
- explicit keys: `['30:0.6:0.2', '30:0.6:0.8', '40:0.6:0.8', '52:0.6:0.8']`
- selected (n,keep_ratio,gamma) strata: `4`
- refined n_perm: `100000`

Outputs:
- `cg_blockperm_delta_rank_stratified_refined.csv`

Best (lowest p) strata (refined):
- `full` n=30 keep=0.60 gamma=0.2: rho=0.553, p=3e-05
- `full` n=30 keep=0.60 gamma=0.8: rho=0.506, p=6e-05
- `switch` n=30 keep=0.60 gamma=0.2: rho=0.486, p=0.0002
- `no_switch` n=30 keep=0.60 gamma=0.8: rho=0.475, p=0.00024
- `switch` n=30 keep=0.60 gamma=0.8: rho=0.444, p=0.00085
- `full` n=52 keep=0.60 gamma=0.8: rho=0.361, p=0.00687
- `full` n=40 keep=0.60 gamma=0.8: rho=0.361, p=0.00704
- `switch` n=40 keep=0.60 gamma=0.8: rho=0.340, p=0.01128
- `no_switch` n=30 keep=0.60 gamma=0.2: rho=0.262, p=0.05454
- `switch` n=52 keep=0.60 gamma=0.8: rho=0.262, p=0.05458

Worst (highest p) strata (refined):
- `no_switch` n=52 keep=0.60 gamma=0.8: rho=0.247, p=0.06989
- `no_switch` n=40 keep=0.60 gamma=0.8: rho=0.252, p=0.05961
- `switch` n=52 keep=0.60 gamma=0.8: rho=0.262, p=0.05458
- `no_switch` n=30 keep=0.60 gamma=0.2: rho=0.262, p=0.05454
- `switch` n=40 keep=0.60 gamma=0.8: rho=0.340, p=0.01128
- `full` n=40 keep=0.60 gamma=0.8: rho=0.361, p=0.00704
- `full` n=52 keep=0.60 gamma=0.8: rho=0.361, p=0.00687
- `switch` n=30 keep=0.60 gamma=0.8: rho=0.444, p=0.00085
- `no_switch` n=30 keep=0.60 gamma=0.8: rho=0.475, p=0.00024
- `switch` n=30 keep=0.60 gamma=0.2: rho=0.486, p=0.0002
