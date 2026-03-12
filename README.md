# poset_phase — Action-Weighted Poset Ensemble Analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18963421.svg)](https://doi.org/10.5281/zenodo.18963421)

Exact numerical framework for studying geometric phase transitions in finite partially ordered sets (posets). Investigates whether Lorentzian-like structures can emerge as competitive phases against high-entropy non-geometric structures in action-weighted poset ensembles.

## Key Results

1. **Bounded phase transition threshold**: Under action path A2 (neutral + geometric penalty), the critical coupling γ_c for Lorentzian-like 2D vs Kleitman–Rothschild posets remains O(1) across N = 10–44, computed with exact linear extension counts.

2. **Ablation-verified non-circularity**: The phase transition window is primarily driven by two structural constraints — global layer-shape balance (`width_height_balance`) and a dimension-consistency penalty (`dim_consistency_penalty`). The latter does not require anchoring to a target dimension d = 2; a non-target-anchored version fully preserves γ_c.

3. **Combinatorial compressibility window**: A two-sieve framework (antichain compressibility + coarse-grained identity retention) filters 7 candidate families down to Lorentzian-like 2D as the dominant geometric phase.

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

See `RESULTS_INDEX.md` for full index.

### Configuration

All experiments are driven by YAML config files (`config_*.yaml`). Key frozen configs:

| Config | Purpose |
|--------|---------|
| `config_frozen_exact.yaml` | Confirmatory small-N exact mainline |
| `config_confirmatory_medium_exact_scan.yaml` | Confirmatory N=20–44 exact |
| `config_geometric_ablation_gamma_c.yaml` | Ablation experiment |
| `config_noncyclic_dim_replacement_gamma_c.yaml` | Non-cyclic replacement |

## Candidate Poset Families

| Family | Description |
|--------|-------------|
| `lorentzian_like_2d` | 2D causal-diamond–like layered posets |
| `lorentzian_like_3d` | 3D analogue |
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

If you use this code, please cite:

```
@software{poset_phase_2026,
  title = {poset\_phase: Action-Weighted Poset Ensemble Analysis},
  year = {2026},
  url = {https://github.com/unicome37/poset_phase},
  doi = {10.5281/zenodo.18963421}
}
```

## License

MIT License. See [LICENSE](LICENSE).
