# Prediction D: N=52 Denoise Frozen Protocol

> Last updated: 2026-03-18  
> Purpose: define a *frozen* confirmatory protocol for the N=52 window where the third-layer signal is variance-sensitive under the standard settings.

---

## Frozen Target Window

- `N = 52`
- `keep_ratio = 0.6`
- `gamma = 0.2`
- `action_mode = A2`
- `zeta scan`: use the default list in `prediction_d_dynamic_validation.py`
- `eta = 0.0` (do not mix global-consistency penalty into the confirmatory claim)

## Frozen “Denoise” Settings

These are frozen to reduce variance at N=52:

- `sis_runs = 128`
- `samples_per_family = 24`
- `coarse_grain.repeats = 30`
- families: same as `config_prediction_d_grid.yaml` (Lor2D, MLR, KR, interval_order, random_layered_k6_uniform, transitive_percolation)
- independence: enforce via `experiment.seed_offset` (each confirmatory replication uses a new seed_offset)

## Keep-Ratio Side-Window Attempts (Standard Protocol)

To test whether a *standard-protocol* alternative to the denoise window exists at `N=52`, we tested nearby keep ratios at `gamma=0.2` with:

- `sis_runs=128`, `samples_per_family=8`, `coarse_grain.repeats=10`, `action_mode=A2`

Results (v9 stratified quick-check):

- `keep_ratio=0.65`: not a viable window.
  - Rep1 (`seed_offset=181000`): `full rho≈0.03 (p≈0.84)`, `no_switch rho≈0.09 (p≈0.54)`.
  - Rep2 (`seed_offset=182000`): `full rho≈-0.18 (p≈0.17)`, `no_switch rho≈-0.37 (p≈0.007)` (sign flip).
- `keep_ratio=0.75`: not a viable window.
  - Rep1 (`seed_offset=183000`): `full rho≈0.25 (p≈0.068)`, `no_switch rho≈0.06 (p≈0.67)`.
  - Rep2 (`seed_offset=184000`): `full rho≈0.065 (p≈0.63)`, `no_switch rho≈0.034 (p≈0.82)`.

Therefore, at least for `keep_ratio in {0.65,0.75}`, we do not observe a stable standard-protocol window that improves over the existing `keep_ratio=0.6` (which itself is variance-sensitive at N=52).

## What Counts as Passing (Layer 3)

On this frozen window, we require (at minimum):

- `full/switch/no_switch` all have `obs_mean_spearman > 0`
- and are significant under stratified block permutation after refinement:
  - `n_perm = 200000` (so p-resolution floor is `~5e-06`)

Outputs used as evidence:

- v9: `outputs_confirmatory/prediction_d_dynamic_v9_window52_denoise_confirm_v*/`
- v10: `outputs_confirmatory/prediction_d_dynamic_v10_window52_denoise_confirm_v*_g02/`
- summary: `outputs_confirmatory/prediction_d_dynamic_v10_window52_denoise_g02_summary.csv`

## Repro Commands

1. Generate the CG data:

```powershell
Set-Location 'd:\Kiro\理论体系\poset_phase'
python .\experiment_cg.py --config .\config_prediction_d_freeze_window52_denoise_confirm_v1.yaml
python .\experiment_cg.py --config .\config_prediction_d_freeze_window52_denoise_confirm_v2.yaml
```

2. Run block-permutation outputs:

```powershell
python .\prediction_d_dynamic_validation.py --only-blockperm --cg-source-contains prediction_d_freeze_window52_denoise_confirm_v1 --n-perm-blockperm 20000 --n-perm-stratified 1000 --out outputs_confirmatory/prediction_d_dynamic_v9_window52_denoise_confirm_v1
python .\prediction_d_dynamic_validation.py --only-blockperm --cg-source-contains prediction_d_freeze_window52_denoise_confirm_v2 --n-perm-blockperm 20000 --n-perm-stratified 1000 --out outputs_confirmatory/prediction_d_dynamic_v9_window52_denoise_confirm_v2
```

3. Refine the frozen stratum:

```powershell
python .\prediction_d_strata_refine.py --zeta-rank-csv outputs_confirmatory/prediction_d_dynamic_v9_window52_denoise_confirm_v1/cg_zeta_scan_rankings_variants.csv --keys 52:0.6:0.2 --require-variants full,switch,no_switch --n-perm 200000 --out outputs_confirmatory/prediction_d_dynamic_v10_window52_denoise_confirm_v1_g02
python .\prediction_d_strata_refine.py --zeta-rank-csv outputs_confirmatory/prediction_d_dynamic_v9_window52_denoise_confirm_v2/cg_zeta_scan_rankings_variants.csv --keys 52:0.6:0.2 --require-variants full,switch,no_switch --n-perm 200000 --out outputs_confirmatory/prediction_d_dynamic_v10_window52_denoise_confirm_v2_g02
```
