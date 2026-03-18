# Prediction D: Layer-Controlled Tier-3 Refinement

- input: `outputs_confirmatory/prediction_d_dynamic_v9_confirm_rep5/cg_zeta_scan_rankings_variants.csv`
- keys: `52:0.6:0.2`
- variants: `full,switch,no_switch`
- n_perm per block: `20000`

## Stratum Results (Mean Across Blocks)

|   n |   keep_ratio |   gamma | icg_variant   |   n_blocks |   n_perm |   obs_mean_spearman_raw |   p_abs_mean_spearman_raw |   obs_mean_spearman_partial_lc_gap |   p_abs_mean_spearman_partial_lc_gap |   obs_mean_spearman_partial_lc |   p_abs_mean_spearman_partial_lc |   obs_mean_spearman_partial_gap |   p_abs_mean_spearman_partial_gap |   mean_spearman_icg_layer_count |   mean_spearman_icg_mean_layer_gap |
|----:|-------------:|--------:|:--------------|-----------:|---------:|------------------------:|--------------------------:|-----------------------------------:|-------------------------------------:|-------------------------------:|---------------------------------:|--------------------------------:|----------------------------------:|--------------------------------:|-----------------------------------:|
|  52 |          0.6 |     0.2 | full          |         11 |    20000 |                0.615024 |               4.99975e-05 |                           0.702694 |                          4.99975e-05 |                       0.733065 |                      4.99975e-05 |                        0.733065 |                       4.99975e-05 |                      -0.0285714 |                         -0.0285714 |
|  52 |          0.6 |     0.2 | no_switch     |         11 |    20000 |                0.553192 |               4.99975e-05 |                           0.536659 |                          0.000149993 |                       0.597399 |                      4.99975e-05 |                        0.597399 |                       4.99975e-05 |                      -0.0285714 |                         -0.0285714 |
|  52 |          0.6 |     0.2 | switch        |         11 |    20000 |                0.564377 |               4.99975e-05 |                           0.58769  |                          4.99975e-05 |                       0.64843  |                      4.99975e-05 |                        0.64843  |                       4.99975e-05 |                      -0.0285714 |                         -0.0285714 |
