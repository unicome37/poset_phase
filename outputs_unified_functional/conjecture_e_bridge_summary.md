# Conjecture E bridge summary

- bd_ratio scan: `outputs_unified_functional/conjecture_e_bd_ratio_scan.csv`
- candidate compare: `outputs_unified_functional/conjecture_e_bridge_candidate_compare.csv`
- output CSV: `outputs_unified_functional/conjecture_e_bridge_summary.csv`

## Current recommendation

The best-supported third-layer bridge proxy remains `bd_ratio`: it has the largest residual gain beyond `F5 ~ N + family`, while the best finite `alpha` in the ordering scans still stays at `0`.

## Compact evidence table

| candidate | ΔR² | coef | residual corr | best alpha | best win rate |
|---|---:|---:|---:|---:|---:|
| bd_ratio | +0.0178 | -18.5279 | -0.2786 | +0.00 | 1.000 |
| bdg_d2_corrected_norm | +0.0013 | -1.3998 | -0.1079 | +0.00 | 1.000 |

## Next step

With the empirical ordering now stabilized, the next useful move is to write down the explicit theoretical map from `bd_ratio` to the intended `S_BD` correction term and test that map against the same data.
