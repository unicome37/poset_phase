# poset_phase v1.0.1

## Action-Weighted Poset Ensemble Analysis: Phase Transitions in Discrete Causal Orders

### Highlights

- **Exact linear-extension entropy** for 7 poset families at N=10–44 using bitmask dynamic programming
- **Bounded phase transition thresholds** γ_c ∈ O(1) across all tested sizes (N=10–44), confirmed via frozen confirmatory mainline
- **Systematic ablation study**: 7 geometric penalty sub-terms reduced to 2-term skeleton (width-height balance + dimension proxy)
- **Non-circularity verification**: dimension proxy replaced by target-free dimension consistency penalty with γ_c preserved at all N
- **Two-sieve framework**: combinatorial compressibility + coarse-graining identity stability for structure selection
- **MLR survivor analysis**: multi-layer random survivors within compressibility window show lor2d_nearest_rate = 1.00

### What's New in v1.0.1

- Updated `.gitignore` to exclude LaTeX build artifacts and compiled PDFs
- All prior content from v1.0.0 intact (code, configs, data, manuscripts)

### Contents

- `generators.py` — 8 poset family generators (Lorentzian 2D/3D/4D, KR-like, transitive percolation, multi-layer random, interval order, absolute layered)
- `entropy_exact.py` — Exact linear extension count via bitmask DP
- `observables_geo.py` — 7-component geometric penalty with `dimension_consistency_penalty` alternative
- `experiment.py` — Main experiment loop with config-driven parameter scanning
- `geometric_ablation_gamma_c.py` — Systematic single-removal and combined ablation
- `noncyclic_dim_replacement_gamma_c.py` — Non-circular replacement experiments
- `manuscript.tex` / `manuscript_foundphys.tex` — LaTeX manuscripts (standard article / Foundations of Physics format)
- `manuscript_figures/` — Publication-quality figures (Fig 1–3, PNG 300dpi + PDF)
- `outputs_confirmatory/` — Frozen mainline results (exact γ_c, timing benchmarks)
- `outputs_exploratory/` — Two-sieve analysis, pairwise duels, survivor profiling

### Reproducibility

```bash
pip install -r requirements.txt
python experiment.py --config config_frozen_exact.yaml   # Confirmatory exact mainline
python geometric_ablation_gamma_c.py                     # Ablation study
python noncyclic_dim_replacement_gamma_c.py              # Non-circular replacement
```

### Citation

If you use this code, please cite:
> Action-weighted partition functions on finite poset ensembles: bounded phase transitions and ablation-verified geometric selection. (2026). https://github.com/unicome37/poset_phase
