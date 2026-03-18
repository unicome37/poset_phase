# Prediction D v8 Refine

- v6 strata CSV: `outputs_exploratory/prediction_d_dynamic_v6/cg_blockperm_delta_rank_stratified.csv`
- zeta-rank CSV: `outputs_exploratory/prediction_d_dynamic_v7_blockperm/cg_zeta_scan_rankings_variants.csv`
- require_variants: `['full', 'switch', 'no_switch']`
- alpha (v6 selection): `0.05`
- require_positive: `True`
- selected (n,keep_ratio,gamma) strata: `14`
- refined n_perm: `20000`

Outputs:
- `selected_strata_v6.csv`
- `cg_blockperm_delta_rank_stratified_refined.csv`

Best (lowest p) strata (refined):
- `full` n=16 keep=0.60 gamma=0.2: rho=0.621, p=5e-05
- `full` n=16 keep=0.80 gamma=0.2: rho=0.756, p=5e-05
- `full` n=30 keep=0.60 gamma=0.2: rho=0.792, p=5e-05
- `full` n=30 keep=0.60 gamma=0.8: rho=0.610, p=5e-05
- `full` n=30 keep=0.70 gamma=0.8: rho=0.636, p=5e-05
- `full` n=30 keep=0.90 gamma=0.8: rho=0.725, p=5e-05
- `full` n=40 keep=0.60 gamma=0.8: rho=0.657, p=5e-05
- `full` n=52 keep=0.70 gamma=0.8: rho=0.548, p=5e-05
- `no_switch` n=16 keep=0.60 gamma=0.2: rho=0.761, p=5e-05
- `full` n=52 keep=0.60 gamma=0.8: rho=0.600, p=5e-05

Worst (highest p) strata (refined):
- `switch` n=52 keep=0.70 gamma=0.8: rho=0.288, p=0.03425
- `no_switch` n=16 keep=0.80 gamma=0.2: rho=0.288, p=0.03195
- `no_switch` n=52 keep=0.80 gamma=0.8: rho=0.314, p=0.01885
- `no_switch` n=52 keep=0.60 gamma=0.2: rho=0.319, p=0.01655
- `switch` n=30 keep=0.80 gamma=0.2: rho=0.325, p=0.01595
- `no_switch` n=30 keep=0.80 gamma=0.2: rho=0.340, p=0.0116
- `full` n=30 keep=0.80 gamma=0.2: rho=0.340, p=0.01085
- `no_switch` n=52 keep=0.90 gamma=0.8: rho=0.366, p=0.00755
- `switch` n=52 keep=0.90 gamma=0.8: rho=0.361, p=0.0071
- `switch` n=40 keep=0.80 gamma=0.2: rho=0.361, p=0.00605
