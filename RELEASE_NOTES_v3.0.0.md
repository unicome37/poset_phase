# poset_phase v3.0.0

## Prediction C: Hierarchy Depth Observables Predict Combinatorial Entropy

### Highlights

- **Hierarchy Integration Index (HII)**: Pre-registered composite z-score of five structural observables (layer_count, mean_layer_gap, long_edge_fraction, adjacent_edge_fraction, reduction_edge_density)
- **Three-tier validation design**: (i) all-family partial correlation (8 families, N=10–16, 320 samples), (ii) matched-pair Δ-analysis (46 Lor2D–MLR pairs, N=30–56), (iii) coarse-graining identity stability (92 samples)
- **Strong negative correlation**: Tier 2 r = −0.834, stable across three filter stringencies (variation < 0.005)
- **Simpson's Paradox diagnosed**: Naïve r = +0.336 → r = −0.578 after controlling for N
- **Component decomposition**: `layer_count` and `mean_layer_gap` carry most signal; full HII never exceeds best constituent
- **Full manuscript**: 788 lines, 8,720 words, 14 tables, three rounds of review completed

### Key Quantitative Results

| Tier | Design | Metric | Value |
|------|--------|--------|-------|
| 1 | All-family partial correlation | r(HII, log H \| N) | −0.578 |
| 1 | Fixed-N within-family | r range | −0.53 to −0.86 |
| 2 | Matched-pair (P5–P95, 46 pairs) | r(ΔHII, Δlog H) | −0.834 |
| 2 | Matched-pair (P10–P90, 34 pairs) | r | −0.836 |
| 2 | Matched-pair (P0–P100, 50 pairs) | r | −0.839 |
| 3 | CG stability (92 samples) | layer_count → σ_CG | −0.874 |
| 3 | CG stability (92 samples) | HII → σ_CG | −0.820 |

### New Files

#### Prediction C Manuscripts (`preC/`)
- `preC/MANUSCRIPT_PredictionC_Full.md` — Final merged manuscript (788 lines)
- `preC/MANUSCRIPT_Section1_Abstract_Introduction.md` — Abstract + Introduction
- `preC/MANUSCRIPT_Section2_Methods.md` — Methods (ensemble, HII, matching, stats)
- `preC/MANUSCRIPT_Section3_Results.md` — Results (all 3 tiers)
- `preC/MANUSCRIPT_Section4_SimpsonsParadox.md` — Simpson's Paradox diagnosis
- `preC/MANUSCRIPT_Section5_ComponentDecomposition.md` — Component decomposition
- `preC/MANUSCRIPT_Section6_Discussion.md` — Discussion, limitations, conclusion
- `preC/merge_script.py` — Automated section-to-full merger

#### Prediction C Code
- `prediction_c_comprehensive.py` — Three-tier comprehensive validation driver
- `prediction_c_pairwise_validation.py` — Tier 2 matched-pair validation
- `prediction_c_pooled_regression.py` — Pooled regression analysis
- `prediction_c_bronze_matched_validation.py` — Bronze-level matched validation
- `prediction_c_switch_enhancement_scan.py` — Switch enhancement scan (Tier 3)
- `_simpson_analysis.py` — Simpson's Paradox analysis script
- `augment_prediction_c_features.py` — Feature augmentation for C analysis

#### Prediction C Configs
- `config_prediction_c_comprehensive.yaml` — Main 3-tier experiment
- `config_prediction_c_pairwise_validation*.yaml` — Pairwise validation (multiple stringencies)
- `config_prediction_c_tier3_extended.yaml` — Extended Tier 3 (N=52/56)
- `config_prediction_c_pooled_regression.yaml` — Pooled regression
- `config_prediction_c_switch_enhancement_scan*.yaml` — Switch enhancement variants
- `config_frozen_exploratory_submodel_duel_{44,48,52,56}_*.yaml` — Per-N duel configs

#### Prediction C Results (`outputs_exploratory/`)
- `prediction_c_comprehensive/` — Original 3-tier results
- `prediction_c_pairwise_validation_nearwall_moderate/` — Definitive Tier 2 (46 pairs)
- `prediction_c_tier3_extended/` — Extended Tier 3 (92 samples)
- `frozen_exploratory_submodel_{52,56}_moderate/` — Per-N duel data
- `matched_residual_freedom_check_nearwall_moderate/` — Residual metrics
- `mlr_survivor_matched_lor2d_nearwall_moderate/` — MLR survivors for matching

### Breaking Changes

None. All Prediction A and B code and results remain unchanged.

### Full Changelog

See git log from v2.1.0 to v3.0.0.
