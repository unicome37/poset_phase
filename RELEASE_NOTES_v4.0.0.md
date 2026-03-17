## What's New in v4.0.0

### Prediction C Major Upgrade: From Correlational Evidence to Quasi-Causal Intervention

Upgraded from correlational-only evidence to a **9-experiment quasi-causal evidence tower**:

| # | Experiment | Cohen's d / Effect | p-value |
|---|-----------|-------------------|---------|
| 1 | Fisher z stratified regression | \|r\| ~ 0.35–0.54 | CI excludes 0 |
| 2 | Single-edge intervention | d = 1.05 | < 10⁻³² |
| 3 | Placebo-controlled (Lor 2D/3D/4D) | d = 1.40–1.83 | < 10⁻¹³³ |
| 4 | Reverse intervention (merge) | d = 2.68, 100% directional | < 10⁻⁴ |
| 5 | Dose-response (k = 2→8) | All slopes < 0 | < 10⁻²⁹ |
| 6 | Edge density universality | No phase transition | All p < 10⁻⁷ |
| 7 | N-scaling power law | \|slope\| ∝ N⁰·⁷⁰ | Growing with N |
| 8 | Analytical bound (complete layered) | Strict monotone decrease | 10⁻¹⁵ precision |
| 9 | Large-N SIS (N ≤ 36) | d = 0.69–1.29 | < 10⁻⁵ all N |

### Manuscript
- **Title**: *Hierarchy Depth Predicts Combinatorial Entropy in Finite Causal Posets: From Correlational Evidence to Quasi-Causal Intervention*
- **File**: `preC/MANUSCRIPT_PredictionC_Full.md` (1,020 lines / 11,915 words)
- **PDF**: `preC/MANUSCRIPT_PredictionC_v2.pdf` (43 pages, 816 KB)

### New Experiment Scripts (13 total, 2,796 lines)
- `prediction_c_stratified_regression.py` — Fisher z within-N correction
- `prediction_c_intervention.py` — Single-edge split intervention
- `prediction_c_dose_response.py` — k-family dose-response
- `prediction_c_placebo_intervention.py` — Same-layer vs cross-layer placebo
- `prediction_c_analytical_bound.py` — Complete layered poset theorem
- `prediction_c_reverse_intervention.py` — Edge removal → layer merge → entropy increase
- `prediction_c_edge_density_phase.py` — p = 0.05→1.0 density sweep
- `prediction_c_n_scaling.py` — N⁰·⁷⁰ power-law scaling
- `prediction_c_large_n_intervention.py` — SIS paired-seed design (N = 16–36)

### Output Data
All results in `outputs_exploratory/prediction_c_*/` with CSV tables and figures.

### Zenodo
`.zenodo.json` updated to v4.0.0 for DOI archival.
