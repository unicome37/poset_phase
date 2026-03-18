# Prediction D: Layer-Controlled Tier-3 Refinement

- input: `outputs_confirmatory/prediction_d_dynamic_v9_confirm_rep3/cg_zeta_scan_rankings_variants.csv`
- keys: `30:0.6:0.2`
- variants: `full,switch,no_switch`
- n_perm per block: `20000`

## Stratum Results (Mean Across Blocks)

|   n |   keep_ratio |   gamma | icg_variant   |   n_blocks |   n_perm |   obs_mean_spearman_raw |   p_abs_mean_spearman_raw |   obs_mean_spearman_partial_lc_gap |   p_abs_mean_spearman_partial_lc_gap |   obs_mean_spearman_partial_lc |   p_abs_mean_spearman_partial_lc |   obs_mean_spearman_partial_gap |   p_abs_mean_spearman_partial_gap |   mean_spearman_icg_layer_count |   mean_spearman_icg_mean_layer_gap |
|----:|-------------:|--------:|:--------------|-----------:|---------:|------------------------:|--------------------------:|-----------------------------------:|-------------------------------------:|-------------------------------:|---------------------------------:|--------------------------------:|----------------------------------:|--------------------------------:|-----------------------------------:|
|  30 |          0.6 |     0.2 | full          |         11 |    20000 |                0.78077  |               4.99975e-05 |                           0.829409 |                          4.99975e-05 |                       0.836762 |                      4.99975e-05 |                        0.836762 |                       4.99975e-05 |                       -0.314286 |                          -0.314286 |
|  30 |          0.6 |     0.2 | no_switch     |         11 |    20000 |                0.716226 |               4.99975e-05 |                           0.708734 |                          4.99975e-05 |                       0.714371 |                      4.99975e-05 |                        0.714371 |                       4.99975e-05 |                       -0.142857 |                          -0.142857 |
|  30 |          0.6 |     0.2 | switch        |         11 |    20000 |                0.810413 |               4.99975e-05 |                           0.825193 |                          4.99975e-05 |                       0.831254 |                      4.99975e-05 |                        0.831254 |                       4.99975e-05 |                       -0.2      |                          -0.2      |
