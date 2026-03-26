# Prediction A — F9 Refinement

F9 = logH − μ·N·R + γ·N·(R−R*)² − λ·Σ_hist + η·ξ_d

## Phase 1: Fine Grid Search

Goal: find (μ, R*, γ, λ, η) where 4D is rank 1 at ALL N and ALL pairwise >80%.

### Top 20 Parameter Sets

| # | μ | R* | γ | λ | η | ranks | min(4D<2D) | min(4D<3D) | min(4D<5D) | score |
|---|---|-----|---|---|---|-------|------------|------------|------------|-------|

## Phase 2: Best Parameter Set Analysis

**F9** = logH − 0.00·N·R + 3·N·(R − 0.100)² − 0·Σ_hist + 0.6·ξ_d

### Per-N Dimension Ordering

| N | F9(2D) | F9(3D) | F9(4D) | F9(5D) | ordering | 4D_rank |
|---|--------|--------|--------|--------|----------|---------|
| 20 | 38.6 | 32.7 | 33.4 | 35.9 | 3D < 4D < 5D < 2D | 2/4 |
| 36 | 93.6 | 79.3 | 75.6 | 80.7 | 4D < 3D < 5D < 2D | 1/4 |
| 52 | 155.6 | 133.1 | 124.7 | 129.4 | 4D < 5D < 3D < 2D | 1/4 |
| 72 | 238.3 | 210.7 | 191.6 | 197.4 | 4D < 5D < 3D < 2D | 1/4 |
| 100 | 359.9 | 327.1 | 294.9 | 296.2 | 4D < 5D < 3D < 2D | 1/4 |

### Per-N Pairwise Win Rates

| N | 4D<2D | 4D<3D | 4D<5D |
|---|----------|----------|----------|
| 20 | 92% (p=0.0000) ★★★ | 30% (p=0.9967)  | 85% (p=0.0000) ★★★ |
| 36 | 100% (p=0.0000) ★★★ | 88% (p=0.0000) ★★★ | 95% (p=0.0000) ★★★ |
| 52 | 100% (p=0.0000) ★★★ | 98% (p=0.0000) ★★★ | 95% (p=0.0000) ★★★ |
| 72 | 100% (p=0.0000) ★★★ | 100% (p=0.0000) ★★★ | 85% (p=0.0000) ★★★ |
| 100 | 100% (p=0.0000) ★★★ | 100% (p=0.0000) ★★★ | 60% (p=0.0897)  |

### Component Decomposition

| N | family | logH | −μNR | γN(R−R*)² | −λΣ | ηξ | F9 |
|---|--------|------|------|-----------|-----|-----|-----|
| 20 | Lor2D | 21.5 | -0.0 | 16.4 | -0.0 | 0.7 | 38.6 |
| 20 | Lor3D | 28.6 | -0.0 | 3.6 | -0.0 | 0.5 | 32.7 |
| 20 | Lor4D | 32.9 | -0.0 | 0.5 | -0.0 | 0.0 | 33.4 |
| 20 | Lor5D | 35.5 | -0.0 | 0.4 | -0.0 | 0.0 | 35.9 |
| 36 | Lor2D | 48.7 | -0.0 | 44.3 | -0.0 | 0.6 | 93.6 |
| 36 | Lor3D | 65.8 | -0.0 | 12.6 | -0.0 | 0.8 | 79.3 |
| 36 | Lor4D | 74.4 | -0.0 | 1.1 | -0.0 | 0.1 | 75.6 |
| 36 | Lor5D | 80.4 | -0.0 | 0.4 | -0.0 | 0.0 | 80.7 |
| 52 | Lor2D | 79.4 | -0.0 | 75.7 | -0.0 | 0.6 | 155.6 |
| 52 | Lor3D | 106.1 | -0.0 | 26.2 | -0.0 | 0.8 | 133.1 |
| 52 | Lor4D | 119.7 | -0.0 | 4.6 | -0.0 | 0.4 | 124.7 |
| 52 | Lor5D | 129.2 | -0.0 | 0.3 | -0.0 | 0.0 | 129.4 |
| 72 | Lor2D | 119.7 | -0.0 | 118.1 | -0.0 | 0.5 | 238.3 |
| 72 | Lor3D | 160.5 | -0.0 | 49.6 | -0.0 | 0.7 | 210.7 |
| 72 | Lor4D | 183.6 | -0.0 | 7.7 | -0.0 | 0.3 | 191.6 |
| 72 | Lor5D | 197.1 | -0.0 | 0.2 | -0.0 | 0.0 | 197.4 |
| 100 | Lor2D | 180.2 | -0.0 | 179.4 | -0.0 | 0.4 | 359.9 |
| 100 | Lor3D | 241.1 | -0.0 | 85.4 | -0.0 | 0.6 | 327.1 |
| 100 | Lor4D | 273.8 | -0.0 | 20.6 | -0.0 | 0.4 | 294.9 |
| 100 | Lor5D | 295.2 | -0.0 | 1.0 | -0.0 | 0.1 | 296.2 |

## Phase 3: Physical Interpretation

The quadratic R-well γ·N·(R−R*)² has a clear physical meaning:

- R (occupancy/density) is a **dimension proxy**: R(2D)≈0.87 > R(3D)≈0.63 > R(4D)≈0.36 > R(5D)≈0.15
- The quadratic well creates a **preferred density band** around R*
- Deviations from R* are penalized proportionally to N (grows O(N))
- If R* ≈ R(4D), then 4D sits at the well minimum while 2D/3D (too dense) and 5D (too sparse) are penalized
- The −μ·N·R linear term shifts the well center and slope

**Physical correspondence**:
- In causal set theory, R = 1 − f_link relates to the causal connectivity density
- R(d) = 1 − E[exp(−ρ·c_d·τ^d)] where c_d is a dimension-dependent volume factor
- The quadratic well R(d)=R* effectively selects d* by matching causal connectivity
- This is analogous to the Einstein-Hilbert action selecting d=4 through the balance
  between R·√g (curvature prefers low d) and the path integral measure (entropy prefers high d)

## Phase 4: F9 Cross-Check with Prediction B

Does F9(Lor) < F9(KR) hold?

| N | pair | win% | sig |
|---|------|------|-----|
| 20 | Lor2D<KR | 0% |  |
| 20 | Lor3D<KR | 70% | ★★ |
| 20 | Lor4D<KR | 48% |  |
| 20 | Lor5D<KR | 5% |  |
| 36 | Lor2D<KR | 0% |  |
| 36 | Lor3D<KR | 2% |  |
| 36 | Lor4D<KR | 12% |  |
| 36 | Lor5D<KR | 0% |  |
| 52 | Lor2D<KR | 0% |  |
| 52 | Lor3D<KR | 0% |  |
| 52 | Lor4D<KR | 2% |  |
| 52 | Lor5D<KR | 0% |  |
| 72 | Lor2D<KR | 0% |  |
| 72 | Lor3D<KR | 0% |  |
| 72 | Lor4D<KR | 0% |  |
| 72 | Lor5D<KR | 0% |  |
| 100 | Lor2D<KR | 0% |  |
| 100 | Lor3D<KR | 0% |  |
| 100 | Lor4D<KR | 2% |  |
| 100 | Lor5D<KR | 0% |  |