# Prediction A — F7 Dimension Selection at Large N

**Data**: 36 samples, N ∈ {20, 36, 52}

**F7 model**: §5.10.7 definitive (α₀=16, q=−0.5, λ=10, η=0.6, Rc=0.25, w=0.015)


## Q1: Mean F7 by Family and N

| N | Lor2D | Lor3D | Lor4D | Lor5D | Winner | 4D_rank |
|---|-------|-------|-------|-------|--------|---------|
| 20 | 34.83 | 32.43 | 31.96 | 35.04 | **Lor4D** | 1 |
| 36 | 58.63 | 78.54 | 74.68 | 78.96 | **Lor2D** | 2 |
| 52 | 85.43 | 115.19 | 124.38 | 129.23 | **Lor2D** | 3 |

## Q2: Pairwise Win Rates (Lor4D vs Others)


### Lor4D vs Lor2D

| N | wins | total | win% | mean_ΔF7 | p_value | sig |
|---|------|-------|------|----------|---------|-----|
| 20 | 3 | 3 | 100% | -2.87 | nan |  |
| 36 | 0 | 3 | 0% | +16.05 | nan |  |
| 52 | 0 | 3 | 0% | +38.96 | nan |  |

**Overall**: 3/9 (33.3%)

### Lor4D vs Lor3D

| N | wins | total | win% | mean_ΔF7 | p_value | sig |
|---|------|-------|------|----------|---------|-----|
| 20 | 1 | 3 | 33% | -0.47 | nan |  |
| 36 | 3 | 3 | 100% | -3.85 | nan |  |
| 52 | 0 | 3 | 0% | +9.19 | nan |  |

**Overall**: 4/9 (44.4%)

### Lor4D vs Lor5D

| N | wins | total | win% | mean_ΔF7 | p_value | sig |
|---|------|-------|------|----------|---------|-----|
| 20 | 3 | 3 | 100% | -3.08 | nan |  |
| 36 | 3 | 3 | 100% | -4.28 | nan |  |
| 52 | 3 | 3 | 100% | -4.84 | nan |  |

**Overall**: 9/9 (100.0%)


## Q3: F7 Component Decomposition (Mean per Family per N)

| N | family | logH | −λΣ_hist | ηΞ_d | wall | F7 |
|---|--------|------|----------|------|------|----|
| 20 | Lor2D | 22.31 | -4.25 | +0.77 | 16.00 | 34.83 |
| 20 | Lor3D | 30.13 | -2.75 | +0.00 | 5.04 | 32.43 |
| 20 | Lor4D | 33.95 | -2.00 | +0.00 | 0.00 | 31.96 |
| 20 | Lor5D | 36.77 | -1.75 | +0.00 | 0.00 | 35.04 |
| 36 | Lor2D | 49.10 | -3.19 | +0.80 | 11.93 | 58.63 |
| 36 | Lor3D | 67.50 | -1.67 | +0.78 | 11.92 | 78.54 |
| 36 | Lor4D | 75.54 | -1.39 | +0.00 | 0.52 | 74.68 |
| 36 | Lor5D | 80.19 | -1.25 | +0.00 | 0.00 | 78.96 |
| 52 | Lor2D | 78.16 | -3.17 | +0.52 | 9.92 | 85.43 |
| 52 | Lor3D | 106.34 | -1.63 | +0.55 | 9.92 | 115.19 |
| 52 | Lor4D | 119.31 | -1.25 | +0.45 | 5.86 | 124.38 |
| 52 | Lor5D | 130.07 | -0.87 | +0.00 | 0.00 | 129.23 |

## Q4: N-Scaling of 4D Margin

Margin = mean F7(competitor) − mean F7(Lor4D)

| N | Δ(2D−4D) | Δ(3D−4D) | Δ(5D−4D) | min_margin |
|---|----------|----------|----------|------------|
| 20 | +2.87 | +0.47 | +3.08 | +0.47 |
| 36 | -16.05 | +3.85 | +4.28 | -16.05 |
| 52 | -38.96 | -9.19 | +4.84 | -38.96 |
- **vs 2D** N-trend: ρ=-1.000 (p=0.000) — shrinking
- **vs 3D** N-trend: ρ=-0.500 (p=0.667) — stable
- **vs 5D** N-trend: ρ=+1.000 (p=0.000) — growing

## Q5: Verdict

**PARTIAL**: Lor4D does NOT win at N=[36, 52].

### Component Attribution

Which F7 terms favor 4D?

- N=20: term_logH → Δ(4D−5D) = -2.82 (favors 4D)
- N=20: term_sigma → Δ(4D−5D) = -0.25 (favors 4D)
- N=20: term_xi → Δ(4D−5D) = +0.00 (favors 5D)
- N=20: term_wall → Δ(4D−5D) = -0.00 (favors 4D)

- N=36: term_logH → Δ(4D−5D) = -4.65 (favors 4D)
- N=36: term_sigma → Δ(4D−5D) = -0.14 (favors 4D)
- N=36: term_xi → Δ(4D−5D) = +0.00 (favors 5D)
- N=36: term_wall → Δ(4D−5D) = +0.52 (favors 5D)

- N=52: term_logH → Δ(4D−5D) = -10.76 (favors 4D)
- N=52: term_sigma → Δ(4D−5D) = -0.38 (favors 4D)
- N=52: term_xi → Δ(4D−5D) = +0.45 (favors 5D)
- N=52: term_wall → Δ(4D−5D) = +5.86 (favors 5D)

### Mechanism I: Low-Dimension Rejection

- N=20: wall(2D)=16.00 vs wall(4D)=0.00, R(2D)=0.601 vs R(4D)=0.045
- N=36: wall(2D)=11.93 vs wall(4D)=0.52, R(2D)=0.736 vs R(4D)=0.157
- N=52: wall(2D)=9.92 vs wall(4D)=5.86, R(2D)=0.802 vs R(4D)=0.252

### Mechanism II: High-Dimension Barrier

- N=20: ΔlogH(5D−4D)=+2.82, ΔΞ(5D−4D)=+0.00, ΔΣ(5D−4D)=+0.25
- N=36: ΔlogH(5D−4D)=+4.65, ΔΞ(5D−4D)=-0.00, ΔΣ(5D−4D)=+0.14
- N=52: ΔlogH(5D−4D)=+10.76, ΔΞ(5D−4D)=-0.45, ΔΣ(5D−4D)=+0.38