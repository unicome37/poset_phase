# Prediction D: Independent-Target Quasi-Intervention Check

- input: `outputs_exploratory/prediction_d_perturbation/perturbation_sample_cg.csv`
- target: `Y = Δscore_local` (no CG term)
- predictor: `X = Δpenalty_cg`
- p-values: permutation test, `n_perm=5000` (two-sided, +1 correction)

## Pooled (All N + All Families)

| perturb   |   n_obs |   rho_pooled |   p_pooled |
|:----------|--------:|-------------:|-----------:|
| p05       |     144 |    -0.212415 | 0.010198   |
| p10       |     144 |    -0.156442 | 0.0523895  |
| p20       |     144 |    -0.305884 | 0.00059988 |

## Stratified (by N × family; mean Spearman across strata)

| perturb   |   n_blocks |   obs_mean_spearman_stratified |   p_stratified |
|:----------|-----------:|-------------------------------:|---------------:|
| p05       |         18 |                     -0.141534  |       0.107379 |
| p10       |         18 |                     -0.0449735 |       0.611678 |
| p20       |         18 |                     -0.0608466 |       0.497101 |

## Pooled by N

| perturb   |   n |   n_obs |   rho_pooled |   p_pooled |
|:----------|----:|--------:|-------------:|-----------:|
| p05       |  30 |      48 |   -0.313287  | 0.0253949  |
| p05       |  40 |      48 |   -0.0870604 | 0.552689   |
| p05       |  52 |      48 |   -0.388949  | 0.00679864 |
| p10       |  30 |      48 |   -0.355515  | 0.0129974  |
| p10       |  40 |      48 |   -0.142206  | 0.331334   |
| p10       |  52 |      48 |   -0.0650239 | 0.65107    |
| p20       |  30 |      48 |   -0.527573  | 0.00039992 |
| p20       |  40 |      48 |   -0.241424  | 0.0983803  |
| p20       |  52 |      48 |   -0.181502  | 0.207958   |

## Stratified by N (within N, stratify by family)

| perturb   |   n |   n_blocks |   obs_mean_spearman_stratified |   p_stratified |
|:----------|----:|-----------:|-------------------------------:|---------------:|
| p05       |  30 |          6 |                     -0.166667  |      0.295941  |
| p05       |  40 |          6 |                      0.0833333 |      0.597081  |
| p05       |  52 |          6 |                     -0.34127   |      0.0267946 |
| p10       |  30 |          6 |                     -0.0436508 |      0.774445  |
| p10       |  40 |          6 |                     -0.138889  |      0.370126  |
| p10       |  52 |          6 |                      0.047619  |      0.735453  |
| p20       |  30 |          6 |                     -0.0634921 |      0.69966   |
| p20       |  40 |          6 |                     -0.103175  |      0.514297  |
| p20       |  52 |          6 |                     -0.015873  |      0.922615  |