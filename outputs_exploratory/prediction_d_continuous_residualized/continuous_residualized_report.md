# Prediction D: Continuous-Y + Rich-Residualized Within-Stratum Test

## Design
- **Samples**: 32/family × 6 families × 3 N = 576 per perturbation level
- **X**: Δpenalty_cg (cover-removal perturbation)
- **Y targets**: delta_score_local, delta_log_H, delta_penalty_local, delta_geo_total, delta_penalty_neutral
- **Residualization confounders**: baseline sig_comp, sig_d_eff, sig_height_ratio, sig_width_ratio, sig_degree_var, layer_count, mean_layer_gap, log_H, penalty_local, geo_total, mean_penalty_cg
- **Boundary samples**: |Δpenalty_cg| in [25%, 75%] quantile per stratum
- **Permutations**: 100,000

## Results

| perturb   | y_target              | subset   |   n_obs |   rho_pooled |    p_pooled |    rho_strat |     p_strat |   n_blocks |   mean_within_rho |   neg_frac |
|:----------|:----------------------|:---------|--------:|-------------:|------------:|-------------:|------------:|-----------:|------------------:|-----------:|
| p05       | delta_score_local     | all      |     576 |  -0.1311     | 0.00153998  | -0.0998086   | 0.0186698   |         18 |      -0.0998086   |   0.611111 |
| p05       | delta_score_local     | boundary |     288 |  -0.0943861  | 0.108869    | -0.111928    | 0.0655493   |         18 |      -0.111928    |   0.666667 |
| p05       | delta_log_H           | all      |     576 |   0.120985   | 0.00363996  |  0.0988311   | 0.0191798   |         18 |       0.0988311   |   0.333333 |
| p05       | delta_log_H           | boundary |     288 |   0.0895268  | 0.129519    |  0.101144    | 0.097959    |         18 |       0.101144    |   0.333333 |
| p05       | delta_penalty_local   | all      |     576 |  -0.0519404  | 0.213688    | -0.0394876   | 0.350836    |         18 |      -0.0394876   |   0.611111 |
| p05       | delta_penalty_local   | boundary |     288 |  -0.0391614  | 0.506045    | -0.000326797 | 0.99542     |         18 |      -0.000326797 |   0.444444 |
| p05       | delta_geo_total       | all      |     576 |  -0.0510664  | 0.220588    | -0.0392229   | 0.353966    |         18 |      -0.0392229   |   0.611111 |
| p05       | delta_geo_total       | boundary |     288 |  -0.0384953  | 0.516045    | -0.000163399 | 0.99693     |         18 |      -0.000163399 |   0.444444 |
| p05       | delta_penalty_neutral | all      |     576 |  -0.130095   | 0.00181998  | -0.138645    | 0.00123999  |         18 |      -0.138645    |   0.777778 |
| p05       | delta_penalty_neutral | boundary |     288 |  -0.12948    | 0.0282597   | -0.097549    | 0.108769    |         18 |      -0.097549    |   0.611111 |
| p10       | delta_score_local     | all      |     576 |  -0.124584   | 0.00303997  | -0.075839    | 0.0719093   |         18 |      -0.075839    |   0.611111 |
| p10       | delta_score_local     | boundary |     288 |  -0.167665   | 0.00414996  | -0.111928    | 0.0658893   |         18 |      -0.111928    |   0.722222 |
| p10       | delta_log_H           | all      |     576 |   0.135636   | 0.000909991 |  0.0826002   | 0.0509595   |         18 |       0.0826002   |   0.388889 |
| p10       | delta_log_H           | boundary |     288 |   0.182097   | 0.00192998  |  0.124673    | 0.0409196   |         18 |       0.124673    |   0.277778 |
| p10       | delta_penalty_local   | all      |     576 |   0.0992775  | 0.0168798   |  0.0736193   | 0.0826392   |         18 |       0.0736193   |   0.333333 |
| p10       | delta_penalty_local   | boundary |     288 |   0.147127   | 0.0126999   |  0.10719     | 0.0770192   |         18 |       0.10719     |   0.388889 |
| p10       | delta_geo_total       | all      |     576 |   0.0993852  | 0.0169698   |  0.0713587   | 0.0930291   |         18 |       0.0713587   |   0.333333 |
| p10       | delta_geo_total       | boundary |     288 |   0.145563   | 0.0135599   |  0.105719    | 0.0821392   |         18 |       0.105719    |   0.388889 |
| p10       | delta_penalty_neutral | all      |     576 |  -0.0784797  | 0.0597194   | -0.0783643   | 0.0639994   |         18 |      -0.0783643   |   0.611111 |
| p10       | delta_penalty_neutral | boundary |     288 |  -0.00173915 | 0.97618     | -0.074183    | 0.222608    |         18 |      -0.074183    |   0.555556 |
| p20       | delta_score_local     | all      |     576 |  -0.227898   | 9.9999e-06  | -0.158765    | 0.000199998 |         18 |      -0.158765    |   0.777778 |
| p20       | delta_score_local     | boundary |     288 |  -0.230999   | 9.9999e-05  | -0.179575    | 0.00294997  |         18 |      -0.179575    |   0.722222 |
| p20       | delta_log_H           | all      |     576 |   0.262805   | 9.9999e-06  |  0.184954    | 3.99996e-05 |         18 |       0.184954    |   0.222222 |
| p20       | delta_log_H           | boundary |     288 |   0.256934   | 2.99997e-05 |  0.184314    | 0.00233998  |         18 |       0.184314    |   0.333333 |
| p20       | delta_penalty_local   | all      |     576 |   0.238839   | 9.9999e-06  |  0.174629    | 6.99993e-05 |         18 |       0.174629    |   0.222222 |
| p20       | delta_penalty_local   | boundary |     288 |   0.230842   | 5.99994e-05 |  0.152614    | 0.0117499   |         18 |       0.152614    |   0.277778 |
| p20       | delta_geo_total       | all      |     576 |   0.239962   | 9.9999e-06  |  0.175587    | 6.99993e-05 |         18 |       0.175587    |   0.222222 |
| p20       | delta_geo_total       | boundary |     288 |   0.231248   | 0.000129999 |  0.150817    | 0.0134799   |         18 |       0.150817    |   0.333333 |
| p20       | delta_penalty_neutral | all      |     576 |  -0.162239   | 0.000109999 | -0.186767    | 1.99998e-05 |         18 |      -0.186767    |   0.833333 |
| p20       | delta_penalty_neutral | boundary |     288 |  -0.107621   | 0.0695393   | -0.133824    | 0.0274597   |         18 |      -0.133824    |   0.777778 |

## Interpretation

**14 test(s) reached within-stratum significance.**
This suggests a genuine within-stratum mechanism channel worth pursuing.