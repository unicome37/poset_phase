# poset_phase — Action-Weighted Poset Ensemble Analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18980657.svg)](https://doi.org/10.5281/zenodo.18980657)

Exact numerical framework for studying geometric phase transitions in finite partially ordered sets (posets). Investigates whether Lorentzian-like structures can emerge as competitive phases against high-entropy non-geometric structures in action-weighted poset ensembles.

## Preprints

| Paper | DOI | Status |
|---|---|---|
| Prediction B | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19048146.svg)](https://doi.org/10.5281/zenodo.19048146) | Preprint |
| Prediction A | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19048324.svg)](https://doi.org/10.5281/zenodo.19048324) | Preprint |
| Prediction C | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19048405.svg)](https://doi.org/10.5281/zenodo.19048405) | Preprint |

## Papers

This repository accompanies three manuscripts:

1. **Prediction B** (Companion Paper): "Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy"
   - Manuscripts: `manuscript.tex` / `manuscript_foundphys.tex` / `mdpi_template/entropy_manuscript.tex`
   - Preprint: [DOI:10.5281/zenodo.19048146](https://doi.org/10.5281/zenodo.19048146)
   - Establishes bounded γ_c across N=10–44; ablation-verified non-circularity

2. **Prediction A**: "Dimension-Agnostic Geometric Dominance in Finite Causal Posets: Higher-Dimensional Lorentzian Structures Emerge Under Consistency-Based Actions"
   - Manuscripts: `preA/manuscript.tex` / `preA/manuscript_entropy.tex`
   - Preprint: [DOI:10.5281/zenodo.19048324](https://doi.org/10.5281/zenodo.19048324)
   - 4D Lorentzian unconditional dominance under consistency actions (N=20–72); seed-robust 100% win rate

3. **Prediction C**: "Hierarchy Depth Observables Predict Combinatorial Entropy in Finite Causal Posets: A Three-Tier Correlational Study"
   - Manuscript: `preC/MANUSCRIPT_PredictionC_Full.md`
   - Preprint: [DOI:10.5281/zenodo.19048405](https://doi.org/10.5281/zenodo.19048405)
   - Three-tier validation: HII–entropy correlation r = −0.834 (Tier 2, 46 matched pairs); Simpson's Paradox diagnosed and resolved

## Key Results

### Companion Paper (Prediction B)

1. **Bounded phase transition threshold**: Under action path A2 (neutral + geometric penalty), the critical coupling γ_c for Lorentzian-like 2D vs Kleitman–Rothschild posets remains O(1) across N = 10–44, computed with exact linear extension counts.

2. **Ablation-verified non-circularity**: The phase transition window is primarily driven by two structural constraints — global layer-shape balance (`width_height_balance`) and a dimension-consistency penalty (`dim_consistency_penalty`). The latter does not require anchoring to a target dimension d = 2; a non-target-anchored version fully preserves γ_c.

3. **Combinatorial compressibility window**: A two-sieve framework (antichain compressibility + coarse-grained identity retention) filters 7 candidate families down to Lorentzian-like 2D as the dominant geometric phase.

### Prediction A — Dimension-Agnostic Geometric Dominance

4. **4D unconditional dominance**: Under consistency-based action (g_dim → g_con), 4D Lorentzian-like posets achieve unconditional dominance across all tested N=20–72 and γ values (92/98 configurations).

5. **Growing margin of victory**: The margin of 4D over 2D/3D grows monotonically with N (+7 at N=20 → +57 at N=72).

6. **Seed robustness**: At N=72 across 3 independent generator seeds, 4D posets win 100% under both consistency variants (21/21), vs 43% under original target-anchored action.

### Prediction C — Hierarchy Depth Predicts Entropy

7. **Negative correlation at fixed N**: The Hierarchy Integration Index (HII) — a 5-component composite z-score of depth observables — correlates negatively with log H at fixed N: partial r = −0.578 (Tier 1, 8 families, 320 samples).

8. **Strong matched-pair signal**: Lor2D–MLR matched-pair r(ΔHII, Δlog H) = −0.834 (46 pairs, N=30–56), stable across three filter stringencies (variation < 0.005).

9. **Simpson's Paradox diagnosed**: Naïve cross-N r = +0.336 reverses to −0.578 after controlling for N, establishing N as the dominant sign-determining confound.

10. **Coarse-graining linkage**: layer_count predicts CG switch rate at r = −0.874 (Tier 3, 92 samples), extending the association chain to identity stability (classifier-contingent).

11. **Component decomposition**: `layer_count` and `mean_layer_gap` carry most of the signal; full HII never exceeds its best constituent.

## Installation

```bash
pip install -r requirements.txt
```

Requires Python ≥ 3.10. Dependencies: `numpy`, `pandas`, `matplotlib`, `pyyaml`, `scipy`.

## Quick Start

```bash
# Smoke test
python experiment.py --config config_smoke.yaml

# Confirmatory exact mainline (N=10–16)
python experiment.py --config config_frozen_exact.yaml

# Extended exact mainline (N=20–44)
python experiment.py --config config_confirmatory_medium_exact_scan.yaml
```

## Reproducing Key Results

### Confirmatory γ_c curve (Table 1 in manuscript)

```bash
python experiment.py --config config_frozen_exact.yaml
python experiment.py --config config_confirmatory_medium_exact_scan.yaml
python gamma_c_report.py
```

### Geometric sub-term ablation (Table 2)

```bash
python geometric_ablation_gamma_c.py --config config_geometric_ablation_gamma_c.yaml
```

### Non-cyclic dim_consistency replacement (Table 3)

```bash
python noncyclic_dim_replacement_gamma_c.py --config config_noncyclic_dim_replacement_gamma_c.yaml
```

### Exact timing benchmark

```bash
python benchmark_exact_timing.py --config config_exact_timing.yaml
```

### Prediction C — Three-tier HII–entropy analysis

```bash
# Tier 1: All-family partial correlation (N=10–16)
python prediction_c_comprehensive.py --config config_prediction_c_comprehensive.yaml

# Tier 2: Matched-pair validation (Lor2D vs MLR, N=30–56)
python prediction_c_pairwise_validation.py --config config_prediction_c_pairwise_validation_nearwall_moderate.yaml

# Tier 3: Extended CG stability (N=30–56)
python prediction_c_switch_enhancement_scan.py --config config_prediction_c_tier3_extended.yaml

# Simpson's Paradox analysis
python _simpson_analysis.py
```

### Prediction B/A/C bridge analysis

```bash
python prediction_bac_bridge.py --config config_prediction_bac_bridge.yaml
```

## Project Structure

### Core Modules

| File | Description |
|------|-------------|
| `generators.py` | Poset data structure and 7 family generators |
| `observables.py` | Neutral observables: comparable fraction, degree stats, layer distribution |
| `observables_geo.py` | Geometric penalty: 7 sub-terms |
| `entropy_exact.py` | Exact linear extension count via DP over antichains |
| `entropy_sis.py` | SIS (Sequential Importance Sampling) entropy estimator |
| `action.py` | Action value computation: A1 / A2 / A3 paths |
| `experiment.py` | Main experiment driver with parameter sweep |
| `bootstrap.py` | Poset-level bootstrap with optional SIS error propagation |
| `normalization.py` | Cross-N normalization |

### Analysis Scripts

| File | Description |
|------|-------------|
| `geometric_ablation_gamma_c.py` | Single and combined sub-term ablation |
| `noncyclic_dim_replacement_gamma_c.py` | dim_proxy → dim_consistency replacement |
| `width_height_consistency_scale_scan.py` | Weight scan for consistency-only backbone |
| `pairwise_compressibility_duel.py` | Two-sieve Lor2D vs MLR analysis |
| `mlr_survivor_profile.py` | MLR window survivor structural profiling |
| `pairwise_locality_delta_validation.py` | Paired locality delta ↔ ΔlogH validation |
| `gamma_c_report.py` | γ_c extraction from summary CSV |

### Prediction C Scripts

| File | Description |
|------|-------------|
| `prediction_c_comprehensive.py` | Three-tier HII–entropy validation driver |
| `prediction_c_pairwise_validation.py` | Tier 2 matched-pair Δ-analysis |
| `prediction_c_pooled_regression.py` | Pooled regression analysis |
| `prediction_c_bronze_matched_validation.py` | Bronze-level matched validation |
| `prediction_c_switch_enhancement_scan.py` | Tier 3 CG switch enhancement scan |
| `_simpson_analysis.py` | Simpson's Paradox diagnosis and resolution |
| `augment_prediction_c_features.py` | Feature augmentation for C analysis |
| `prediction_bac_bridge.py` | Bridge analysis: B↔C mechanism and A↔C relation test |

### Prediction A Scripts

| File | Description |
|------|-------------|
| `prediction_a_dimension_scan.py` | 4-family dimension scan (2D/3D/4D/KR) across N=20–72 |
| `prediction_a_geometric_ablation.py` | Geometric ablation for Prediction A families |
| `prediction_a_seed_sensitivity.py` | Seed sensitivity analysis at N=64/68/72 |
| `prediction_a_margin_plot.py` | Margin of victory visualization |
| `prediction_a_winner_phase_plot.py` | Winner phase comparison plot |
| `prediction_a_seed_sensitivity_plot.py` | Seed sensitivity visualization |
| `seed_sensitivity_c_threshold.py` | γ_c threshold seed sensitivity |

### Output Organization

Results are stratified into two evidence tiers:

- **`outputs_confirmatory/`** — Frozen mainline results; citable as primary evidence
  - `frozen_exact/` — Small-N exact (N=10–16)
  - `medium_exact/` — Extended exact (N=20–40)
  - `medium_exact_scan/` — Full exact scan (N=20–44, incl. Lor3D)
  - `exact_timing/` — Exact computation benchmark
- **`outputs_exploratory/`** — Discovery, diagnostics, and method extensions
  - `pairwise_compressibility_duel/` — Two-sieve analysis
  - `pairwise_locality_delta_validation/` — Paired validation
  - `noncyclic_dim_replacement_gamma_c/` — Ablation & replacement
  - `geometric_ablation_gamma_c/` — Sub-term ablation
  - `prediction_a_n{52..72}_mixed/` — Prediction A dimension scans (N=52–72)
  - `prediction_a_margin_summary/` — Margin of victory across all N
  - `prediction_a_phase_summary/` — Winner phase statistics
  - `prediction_a_seed_sensitivity_n{68,72}/` — Seed robustness analysis
  - `prediction_a_geometric_ablation/` — Prediction A ablation results
  - `prediction_a_dim_replacement*/` — Prediction A dimension replacement experiments
  - `prediction_bac_bridge/` — Bridge analysis across Prediction B, A, and C
  - `prediction_c_comprehensive/` — Prediction C Tier 1 all-family results
  - `prediction_c_pairwise_validation_nearwall_moderate/` — Prediction C Tier 2 (46 matched pairs)
  - `prediction_c_tier3_extended/` — Prediction C Tier 3 CG stability (92 samples)
  - `frozen_exploratory_submodel_{52,56}_moderate/` — Near-wall duel data
  - `matched_residual_freedom_check_nearwall_moderate/` — Residual metrics
  - `mlr_survivor_matched_lor2d_nearwall_moderate/` — MLR survivors for matching

See `RESULTS_INDEX.md` for full index.

### Configuration

All experiments are driven by YAML config files (`config_*.yaml`). Key frozen configs:

| Config | Purpose |
|--------|---------|
| `config_frozen_exact.yaml` | Confirmatory small-N exact mainline |
| `config_confirmatory_medium_exact_scan.yaml` | Confirmatory N=20–44 exact |
| `config_geometric_ablation_gamma_c.yaml` | Ablation experiment |
| `config_noncyclic_dim_replacement_gamma_c.yaml` | Non-cyclic replacement |
| `config_prediction_c_comprehensive.yaml` | Prediction C three-tier validation |
| `config_prediction_c_pairwise_validation_nearwall_moderate.yaml` | Prediction C Tier 2 (P5–P95) |
| `config_prediction_c_tier3_extended.yaml` | Prediction C Tier 3 extended |

## Candidate Poset Families

| Family | Description |
|--------|-------------|
| `lorentzian_like_2d` | 2D causal-diamond–like layered posets |
| `lorentzian_like_3d` | 3D analogue |
| `lorentzian_like_4d` | 4D Lorentzian-like posets (new in v2.0.0) |
| `KR_like` | Kleitman–Rothschild 3-layer maximum-entropy posets |
| `multi_layer_random` | Random layered posets with variable layer count |
| `transitive_percolation` | Transitive percolation on a chain |
| `interval_order` | Interval orders from random intervals |
| `absolute_layered` | Fully connected layered posets |

## Action Paths

Three action paths separate neutral and geometric contributions:

- **A1** = −βH + γ · I_neutral (neutral penalty only)
- **A2** = −βH + γ · (I_neutral + I_geometric) (neutral + geometric)
- **A3** = −βH + γ · I_geometric (geometric only)

where H = log(linear extensions) and I = structure penalty.

## Citation

If you use this code, please cite the software and the relevant paper(s):

```bibtex
@software{poset_phase_2026,
  title   = {poset\_phase: Action-Weighted Poset Ensemble Analysis},
  author  = {Zhang, Gang},
  year    = {2026},
  url     = {https://github.com/unicome37/poset_phase},
  doi     = {10.5281/zenodo.18980657}
}

@article{zhang2026predC,
  title   = {Hierarchy Depth Observables Predict Combinatorial Entropy in Finite
             Causal Posets: A Three-Tier Correlational Study},
  author  = {Zhang, Gang},
  journal = {Entropy},
  year    = {2026},
  note    = {Submitted}
}
```

## License

MIT License. See [LICENSE](LICENSE).
