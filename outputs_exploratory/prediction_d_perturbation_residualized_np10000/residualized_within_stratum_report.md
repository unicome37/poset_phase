# Prediction D: Residualized Within-Stratum Test

- input: `D:/Kiro/理论体系/poset_phase/outputs_exploratory/prediction_d_perturbation_n32/perturbation_sample_cg_n32.csv`
- target: `delta_score_local` (continuous, mechanism-independent)
- predictor: `delta_penalty_cg`
- transforms:
  - `raw`: original deltas
  - `demeaned`: within-(perturb,n,family) demeaning
  - `residualized`: within-(perturb,n,family) residualization on `pen0`, `score0`, `icg0`

| perturb   | transform    |   n_obs |   n_blocks |   rho_pooled |   p_pooled |   rho_stratified |   p_stratified |   mean_within_rho |   neg_fraction |
|:----------|:-------------|--------:|-----------:|-------------:|-----------:|-----------------:|---------------:|------------------:|---------------:|
| p05       | raw          |     576 |         18 |  -0.109976   | 0.00679932 |       -0.0581215 |      0.167583  |        -0.0581215 |       0.555556 |
| p05       | demeaned     |     576 |         18 |  -0.0928688  | 0.0228977  |       -0.0581215 |      0.167583  |        -0.0581215 |       0.555556 |
| p05       | residualized |     576 |         18 |  -0.135944   | 0.00149985 |       -0.0860419 |      0.0378962 |        -0.0860419 |       0.555556 |
| p10       | raw          |     576 |         18 |  -0.120639   | 0.00379962 |        0.0111193 |      0.791621  |         0.0111193 |       0.444444 |
| p10       | demeaned     |     576 |         18 |   0.00079258 | 0.983902   |        0.0111193 |      0.791621  |         0.0111193 |       0.444444 |
| p10       | residualized |     576 |         18 |  -0.136879   | 0.00129987 |       -0.100481  |      0.0172983 |        -0.100481  |       0.611111 |
| p20       | raw          |     576 |         18 |  -0.250107   | 9.999e-05  |       -0.0211999 |      0.616638  |        -0.0211999 |       0.555556 |
| p20       | demeaned     |     576 |         18 |  -0.0415941  | 0.319868   |       -0.0211999 |      0.616638  |        -0.0211999 |       0.555556 |
| p20       | residualized |     576 |         18 |  -0.258665   | 9.999e-05  |       -0.191838  |      9.999e-05 |        -0.191838  |       0.722222 |

Interpretation rule:
- if signal reappears after residualization, that supports a within-stratum continuous channel worth pursuing;
- if it stays null, the current perturbation design still supports only structural covariation.