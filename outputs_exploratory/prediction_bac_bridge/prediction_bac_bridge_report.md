# Prediction B/A/C Bridge Analysis

## Core findings

### Prediction B ↔ Prediction C

For `Lor2D` vs `KR_like`, the hierarchy-depth gap tracks both sides of the bounded transition: `corr(delta_hii, delta_penalty) = -0.581` and `corr(delta_hii, delta_log_H) = -0.680`.
Interpretation: deeper Lorentzian-like structures pay an entropy cost (`delta_log_H < 0`) but gain a larger geometric penalty advantage (`delta_penalty < 0`), which is exactly the tradeoff needed for a finite `gamma_c`.
Across the confirmatory `N=20..44` line, the estimated crossing window for `Lor2D` vs `KR_like` stays finite with mean `gamma_cross ≈ 1.032`.

### Prediction A ↔ Prediction C

The A/C relation is not a simple extension of the B/C mechanism. `Lor4D` has higher HII than `Lor2D` in `0/14` tested `N` values, and higher HII than `Lor3D` in `0/14` tested `N` values.
So the 4D wins under consistency actions are not being selected by maximal hierarchy depth. They are being selected on a different axis: higher entropy plus non-target-anchored dimensional consistency.

### Resulting interpretation

Prediction C currently explains a local refinement mechanism inside the geometric window opened by Prediction B: deeper hierarchy suppresses combinatorial entropy among matched quasi-geometric competitors.
Prediction A is adjacent to that mechanism but not reducible to it: once the action is made dimension-agnostic, 4D dominance persists even though 4D is not the deepest family on the C-style hierarchy axis.

## Consistency-action win rates

- `A2_replace_dim_with_consistency`: `lorentzian_like_4d` vs `KR_like` win rate = `1.000`, mean `delta_hii = +0.940`
- `A2_replace_dim_with_consistency`: `lorentzian_like_4d` vs `lorentzian_like_2d` win rate = `1.000`, mean `delta_hii = -1.533`
- `A2_replace_dim_with_consistency`: `lorentzian_like_4d` vs `lorentzian_like_3d` win rate = `0.939`, mean `delta_hii = -0.486`
- `A2_replace_dim_with_multi_consistency`: `lorentzian_like_4d` vs `KR_like` win rate = `1.000`, mean `delta_hii = +0.940`
- `A2_replace_dim_with_multi_consistency`: `lorentzian_like_4d` vs `lorentzian_like_2d` win rate = `1.000`, mean `delta_hii = -1.533`
- `A2_replace_dim_with_multi_consistency`: `lorentzian_like_4d` vs `lorentzian_like_3d` win rate = `0.908`, mean `delta_hii = -0.486`
