# Release Notes — v5.0.0

**Physica A Submission Release**

This release accompanies the manuscript:

> *Constrained structural selection on finite causal orders: historical loading, dimensional admissibility, and Lorentzian windows*
> Gang Zhang, submitted to Physica A: Statistical Mechanics and its Applications (2026)

## What's new in v5.0.0

### Unified structural cost functional (F7)
- Interval-richness admissibility wall with occupancy-identity derivation
- Sigmoid threshold with finite-size N-scaling (q = -0.5)
- Pareto parameter search identifying 915 viable configurations

### Anti-overfitting test battery (`_anti_illusion_battery.py`)
Seven code-level diagnostics:
1. Permutation test (p < 0.001)
2. Random-functional Monte Carlo (0.02% joint match)
3. Ablation with noise replacement (each term non-redundant)
4. 4-fold cross-validation (ΔA = +0.009)
5. Sigmoid-wall N-breakdown (ρ = +0.846, informative)
6. New-seed generalization (ΔA = -0.062)
7. Fake-family discrimination (66.7%)

### Spectral recovery experiments (`_d_recovery_experiment.py`)
- W₁ (interval-spectrum Wasserstein distance): ρ_partial = -0.645 (p < 10⁻⁴)
- ΔH_int (spectral-entropy drift): p = 0.029 after controlling F and N

### Constrained dynamics
- Microcanonical swap dynamics (N = 16–36)
- Cluster moves (N = 36)
- Simulated annealing (N = 64, hit rate 80%)

## Previous versions
- v4.0.0: Prediction C quasi-causal evidence tower (9 experiments)
- v3.0.0: Prediction C correlational
- v2.1.0: Prediction B bounded γ_c
- v2.0.0: Prediction A 4D Lorentzian dominance
