# Conjecture E S_BD mapping check

- Input raw features: `outputs_unified_functional/raw_features.csv`
- Input BD table: `outputs_unified_functional/bd_actions.csv`
- Output CSV: `outputs_unified_functional/conjecture_e_sbd_mapping_check.csv`

## Mapping claim

The current empirical bridge term `bd_ratio` is the same interval-richness quantity used as `S_BD` in the theory note:

`S_BD(X) = (1 - C0 / total_relations) * (1 + mean_k)`

So the operational map is simply:

`S_BD_proxy(X) = bd_ratio(X)`

## Identity check

Max |S_BD_exact - bd_ratio_stored| = 0.000000e+00
Mean |S_BD_exact - bd_ratio_stored| = 0.000000e+00
Spearman(S_BD_exact, bd_ratio_stored) = +1.0000

## Residual bridge

Baseline model: `F5 ~ N + family` with R² = 0.9692
Extended model: `F5 ~ N + family + S_BD_proxy` with R² = 0.9869
ΔR² = +0.0178
Corr(S_BD_proxy, residual F5) = -0.2786

## Interpretation

This is the cleanest possible bridge statement at the current stage: the empirical `bd_ratio` is not just inspired by `S_BD`; it is the same normalized interval-richness observable. So the next theoretical step is not to redefine the metric, but to justify why this interval-richness term should enter the discrete effective action with a positive penalty sign.
