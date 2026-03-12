# Prediction A Pilot Summary

## Scope

- Families: `lorentzian_like_2d`, `lorentzian_like_3d`, `lorentzian_like_4d`, `KR_like`
- Sizes: `N = 20, 24, 28, 32, 36`
- Samples per family: `4`
- Action variants:
  - `A2_full`
  - `A2_no_comp_window`
  - `width_height_dim_consistency`
  - `width_height_dim_multi_consistency`

## Main findings

1. Under `A2_full`, `Lor4D` does **not** become the dominant family in the tested range.
   - Winner counts: `Lor2D = 16`, `Lor4D = 9`, `Lor3D = 5`, `KR = 5`.
   - This is not substantially changed by removing `comparability_window`.

2. Under the reduced non-targeted backbones, `Lor4D` becomes the dominant winner.
   - `width_height_dim_consistency`: `Lor4D` wins `32/35` tested `(N, gamma)` cells.
   - `width_height_dim_multi_consistency`: `Lor4D` wins `32/35` tested `(N, gamma)` cells.

3. Pairwise comparison shows a sharp asymmetry:
   - In the reduced backbones, `Lor4D` beats `Lor2D` in `35/35` tested cells.
   - In the reduced backbones, `Lor4D` beats `Lor3D` in `32/35` tested cells.

4. The low-dimensional preference of `A2_full` is therefore not explained solely by `comparability_window`.
   - `A2_no_comp_window` remains qualitatively similar to `A2_full`.
   - The bias appears to be carried by the broader `A2_full` geometric bundle rather than by that single term alone.

## Current interpretation

The first pilot does not yet support the claim that `3+1` dimensions are selected by the current full geometric action. However, it strongly suggests that once the action is reduced to a more neutral backbone built from `width_height_balance` plus dimensional self-consistency, a `Lor4D` preference becomes visible already at modest `N`.

This means Prediction A is now in a tractable state:

- `A2_full` currently looks low-dimension biased.
- The reduced non-targeted backbones look much more favorable to a `3+1` window.
- The next step should be to identify which remaining `A2_full` subterms push the system back toward `Lor2D`.
