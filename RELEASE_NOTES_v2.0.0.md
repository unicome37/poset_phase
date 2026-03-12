# poset_phase v2.0.0

## Prediction A: Dimension-Agnostic Geometric Dominance in Finite Causal Posets

### Highlights

- **New 4D Lorentzian-like family** (`lorentzian_like_4d`): Four-dimensional causal diamond posets added to the ensemble
- **Dimension scan experiments** (N=20–72): Systematic comparison of 2D, 3D, 4D Lorentzian-like and KR posets under three action variants
- **Unconditional 4D dominance**: Under consistency-based action (g_dim → g_con), 4D posets win 92/98 configurations with growing margin
- **Seed robustness confirmed**: 100% win rate for 4D at N=72 across 3 independent seeds (21/21) under consistency variants
- **Growing margin of victory**: +7 at N=20 → +57 at N=72 (monotonic increase)
- **Full manuscript**: Two LaTeX versions (standard article + MDPI Entropy template) with figures

### New Files

#### Prediction A Manuscripts (`preA/`)
- `preA/manuscript.tex` — Standard LaTeX article format
- `preA/manuscript_entropy.tex` — MDPI Entropy template format
- `preA/MANUSCRIPT_DRAFT.md` — Markdown draft for review
- `preA/README.md` — Prediction A overview
- `preA/cover_letter_entropy_mdpi.md` — Cover letter and submission guide
- `preA/manuscript_figures/` — Publication figures (margin of victory, winner phase comparison)
- `preA/Definitions/` — MDPI LaTeX template files

#### Prediction A Code
- `prediction_a_dimension_scan.py` — Main 4-family dimension scan driver
- `prediction_a_geometric_ablation.py` — Geometric ablation for extended ensemble
- `prediction_a_seed_sensitivity.py` — Seed sensitivity analysis
- `prediction_a_margin_plot.py` — Margin of victory visualization
- `prediction_a_winner_phase_plot.py` — Winner phase comparison plots
- `prediction_a_seed_sensitivity_plot.py` — Seed sensitivity plots
- `seed_sensitivity_c_threshold.py` — γ_c threshold seed sensitivity
- `pairwise_blind_identity_scan.py` — Blind identity pairwise scan
- `pairwise_identity_robustness_plot.py` — Identity robustness visualization
- `pairwise_switch_catchup_analysis.py` — Switch catchup analysis
- `pairwise_switch_catchup_plot.py` — Switch catchup plots
- `plot_c_threshold_summary.py` — γ_c threshold summary plots
- `width_height_degeneration_check.py` — Width-height degeneration validation
- `multiview_dim_consistency_gamma_c.py` — Multi-view dimension consistency

#### Prediction A Configs (`config_prediction_a_*.yaml`)
- Dimension scan configs for N=44–72 (pilot, robust, mixed variants)
- Seed sensitivity configs (N=64, 68, 72)
- Geometric ablation and dimension replacement configs
- C-threshold sensitivity configs (various N ranges)

#### Prediction A Results (`outputs_exploratory/prediction_a_*/`)
- Dimension scan results for N=52–72
- Margin of victory summary across all N
- Winner phase statistics
- Seed sensitivity analysis at N=68, 72
- Geometric ablation results
- Dimension replacement experiments

### Modified Files
- `entropy_exact.py` — Extended for larger N support
- `observables_geo.py` — Added consistency penalty variants
- `runtime_utils.py` — Performance improvements
- `.zenodo.json` — Updated metadata for v2.0.0 with Prediction A
- `README.md` — Added Prediction A documentation

### Performance
- `entropy_dp.c` / `entropy_exact_c.py` — C-accelerated entropy computation backend (new)

### Reproducibility

```bash
# Prediction A pilot scan
python prediction_a_dimension_scan.py --config config_prediction_a_pilot_fast.yaml

# Full Prediction A scan (N=20–72)
python prediction_a_dimension_scan.py --config config_prediction_a_robust.yaml

# Seed sensitivity at N=72
python prediction_a_seed_sensitivity.py --config config_prediction_a_seed_sensitivity_n72.yaml

# Generate figures
python prediction_a_margin_plot.py
python prediction_a_winner_phase_plot.py
```

### Citation

If you use this code, please cite:
> Zhang, G. (2026). Dimension-Agnostic Geometric Dominance in Finite Causal Posets: Higher-Dimensional Lorentzian Structures Emerge Under Consistency-Based Actions. https://github.com/unicome37/poset_phase
