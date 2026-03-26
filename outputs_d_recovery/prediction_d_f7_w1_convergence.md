# F7–W₁ Convergence at Large N: Physical Analysis

**Data**: outputs_d_recovery/prediction_d_large_n_r15.csv (300 samples)


## A1: Comparable Fraction as Common Driver

comp_frac (cf) = fraction of comparable pairs = 2·|{(i,j): i<j}| / N(N-1)

| N | ρ(cf,F7) | ρ(cf,W₁) | ρ(cf,Pb) | ρ(cf,ΔH) |
|---|----------|----------|----------|----------|
| 36 | -0.538 | +0.954 | -0.689 | +0.526 |
| 52 | -0.828 | +0.946 | -0.789 | +0.451 |
| 72 | -0.874 | +0.971 | -0.716 | +0.481 |
| 100 | -0.868 | +0.976 | -0.813 | +0.370 |

## A2: Family Membership Explains Variance (η²)

η² = SS_between / SS_total — how much of each variable is family-determined

| N | η²(F7) | η²(W₁) | η²(cf) | η²(Pb) | η²(ΔH) |
|---|--------|--------|--------|--------|--------|
| 36 | 0.926 | 0.915 | 0.948 | 0.856 | 0.645 |
| 52 | 0.964 | 0.950 | 0.978 | 0.819 | 0.658 |
| 72 | 0.986 | 0.972 | 0.974 | 0.865 | 0.783 |
| 100 | 0.986 | 0.983 | 0.986 | 0.891 | 0.879 |

## A3: Within-Family ρ(F7, W₁)

If convergence is purely family-driven, within-family ρ should be weak

| N | Lor2D | Lor3D | Lor4D | Lor5D | KR | pooled |
|---|-------|-------|-------|-------|-----|--------|
| 36 | +0.436 | -0.332 | -0.464 | +0.361 | -0.023 | -0.475 |
| 52 | -0.029 | +0.154 | -0.389 | +0.229 | -0.286 | -0.812 |
| 72 | -0.286 | -0.350 | -0.414 | -0.350 | -0.264 | -0.870 |
| 100 | -0.096 | -0.368 | -0.786★ | -0.304 | -0.196 | -0.866 |

## A4: Family Mean Profiles

| N | family | mean_F7 | mean_W₁ | mean_cf | F7_rank | W₁_rank | rank_match |
|---|--------|---------|---------|---------|---------|---------|------------|
| 36 | KR | 77.9 | 0.5266 | 0.390 | 4 | 4 | ✓ |
| 36 | Lor2D | 57.5 | 1.0257 | 0.495 | 1 | 5 |  |
| 36 | Lor3D | 75.9 | 0.3043 | 0.289 | 3 | 3 | ✓ |
| 36 | Lor4D | 74.5 | 0.1130 | 0.166 | 2 | 2 | ✓ |
| 36 | Lor5D | 78.8 | 0.0347 | 0.102 | 5 | 1 |  |
| 52 | KR | 120.0 | 0.7277 | 0.380 | 3 | 4 |  |
| 52 | Lor2D | 86.4 | 1.5808 | 0.500 | 1 | 5 |  |
| 52 | Lor3D | 114.8 | 0.4623 | 0.286 | 2 | 3 |  |
| 52 | Lor4D | 123.0 | 0.1396 | 0.172 | 4 | 2 |  |
| 52 | Lor5D | 128.0 | 0.0515 | 0.104 | 5 | 1 |  |
| 72 | KR | 179.3 | 1.0169 | 0.380 | 3 | 4 |  |
| 72 | Lor2D | 123.3 | 2.0685 | 0.496 | 1 | 5 |  |
| 72 | Lor3D | 167.2 | 0.6177 | 0.285 | 2 | 3 |  |
| 72 | Lor4D | 188.1 | 0.1649 | 0.170 | 4 | 2 |  |
| 72 | Lor5D | 196.3 | 0.0502 | 0.103 | 5 | 1 |  |
| 100 | KR | 270.3 | 1.3409 | 0.380 | 3 | 4 |  |
| 100 | Lor2D | 184.3 | 2.8860 | 0.497 | 1 | 5 |  |
| 100 | Lor3D | 245.1 | 0.8824 | 0.287 | 2 | 3 |  |
| 100 | Lor4D | 279.0 | 0.2686 | 0.171 | 4 | 2 |  |
| 100 | Lor5D | 292.0 | 0.0870 | 0.105 | 5 | 1 |  |

## A5: Family Rank Concordance

Kendall's tau between family rankings by F7 vs W₁:

- N=36: Kendall τ = -0.400 (p=0.483)
- N=52: Kendall τ = -0.800 (p=0.083)
- N=72: Kendall τ = -0.800 (p=0.083)
- N=100: Kendall τ = -0.800 (p=0.083)

## A6: Within-Family Residual ρ(F7, W₁)

Remove family means, compute ρ on residuals:

- N=36: ρ_resid(F7, W₁) = -0.087 (p=0.4599)
- N=52: ρ_resid(F7, W₁) = -0.062 (p=0.6000)
- N=72: ρ_resid(F7, W₁) = -0.362 (p=0.0014)
- N=100: ρ_resid(F7, W₁) = -0.327 (p=0.0042)

## A7: Between-Family vs Within-Family Variance Ratio

| N | F7 between/within | W₁ between/within |
|---|-------------------|-------------------|
| 36 | 12.5 | 10.7 |
| 52 | 26.8 | 18.9 |
| 72 | 69.0 | 35.1 |
| 100 | 72.6 | 57.7 |

## Synthesis


### Physical Interpretation

**F7** measures "Lorentzian-ness" — a composite of:
- logH (entropy of linear extensions — captures order complexity)
- Π_geo (geometric functional — path structure)  
- Σ_hist (history functional — information content)
- Ξ_d (dimensional estimator)
- sigmoid wall (density-based admissibility filter)

**W₁** measures "spectral stability under coarse-graining" — the Wasserstein distance
between the interval spectrum before and after 30% node deletion.

Both quantities are fundamentally about **how structured/ordered the poset is**.
A poset with high F7 (Lorentzian-like) has a rigid, self-similar interval spectrum
that is robust to coarse-graining (low W₁). A poset with low F7 (KR/random-like)
has an irregular spectrum that changes significantly under CG (high W₁).

The convergence at large N reflects **spectral rigidity**: as N grows, the interval
spectrum of each family type becomes increasingly determined by the family's geometric
properties alone (less finite-size noise). Both F7 and W₁ become better estimators of
the underlying family type, hence more correlated with each other.

This is analogous to how in statistical mechanics, different thermodynamic observables
become more correlated as system size increases — they all converge to measuring the
same underlying phase/state.
