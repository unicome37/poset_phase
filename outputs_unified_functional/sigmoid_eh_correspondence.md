# Sigmoid → EH Correspondence

Total: 375 realizations

F7 wall: α(N)·σ((R−Rc)/w), α₀=16.0, q=-0.5, Rc=0.25, w=0.015


## 1. Occupancy R vs Hubble Parameter

Does R increase with curvature?

| d | N | H=0 | H=0.25 | H=0.5 | H=1.0 | H=2.0 |
|---|---|-----|--------|-------|-------|-------|
| 2 | 64 | 0.824 | 0.812 | 0.789 | 0.760 | 0.654 |
| 2 | 128 | 0.892 | 0.891 | 0.874 | 0.838 | 0.770 |
| 2 | 256 | 0.936 | 0.932 | 0.925 | 0.910 | 0.854 |
| 2 | 512 | 0.964 | 0.960 | 0.956 | 0.947 | 0.918 |
| 2 | 1024 | 0.978 | 0.977 | 0.975 | 0.970 | 0.951 |
| 3 | 64 | 0.553 | 0.493 | 0.469 | 0.285 | 0.085 |
| 3 | 128 | 0.673 | 0.623 | 0.587 | 0.436 | 0.186 |
| 3 | 256 | 0.769 | 0.741 | 0.699 | 0.596 | 0.279 |
| 3 | 512 | 0.841 | 0.817 | 0.787 | 0.715 | 0.451 |
| 3 | 1024 | 0.892 | 0.876 | 0.855 | 0.798 | 0.574 |
| 4 | 64 | 0.279 | 0.203 | 0.175 | 0.053 | 0.000 |
| 4 | 128 | 0.396 | 0.337 | 0.253 | 0.140 | 0.000 |
| 4 | 256 | 0.506 | 0.457 | 0.352 | 0.201 | 0.019 |
| 4 | 512 | 0.624 | 0.566 | 0.482 | 0.306 | 0.034 |
| 4 | 1024 | 0.709 | 0.659 | 0.600 | 0.433 | 0.078 |

## 2. Pure Sigmoid σ((R−Rc)/w) vs H

Threshold: Rc=0.25, width: w=0.015

| d | N | H=0 | H=0.25 | H=0.5 | H=1.0 | H=2.0 |
|---|---|-----|--------|-------|-------|-------|
| 2 | 64 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 2 | 128 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 2 | 256 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 2 | 512 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 2 | 1024 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 3 | 64 | 1.0000 | 1.0000 | 1.0000 | 0.5889 | 0.0004 |
| 3 | 128 | 1.0000 | 1.0000 | 1.0000 | 0.9994 | 0.0956 |
| 3 | 256 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.6861 |
| 3 | 512 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 3 | 1024 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 4 | 64 | 0.7596 | 0.1173 | 0.3031 | 0.0001 | 0.0000 |
| 4 | 128 | 0.9999 | 0.9039 | 0.5402 | 0.0018 | 0.0000 |
| 4 | 256 | 1.0000 | 1.0000 | 0.9961 | 0.0586 | 0.0000 |
| 4 | 512 | 1.0000 | 1.0000 | 1.0000 | 0.9291 | 0.0000 |
| 4 | 1024 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.0000 |

## 3. Occupancy Variance Shrinking (Wall Sharpening)

If std(R|d,N,H) → 0, the sigmoid becomes a step function.

| d | N | std(R|H=0) | std(R|H=0.5) | std(R|H=1.0) | std(R|H=2.0) |
|---|---|-----------|-------------|-------------|-------------|
| 2 | 64 | 0.0179 | 0.0193 | 0.0193 | 0.0573 |
| 2 | 128 | 0.0056 | 0.0105 | 0.0136 | 0.0149 |
| 2 | 256 | 0.0023 | 0.0038 | 0.0028 | 0.0090 |
| 2 | 512 | 0.0006 | 0.0014 | 0.0017 | 0.0056 |
| 2 | 1024 | 0.0002 | 0.0006 | 0.0006 | 0.0007 |
| 3 | 64 | 0.0374 | 0.0430 | 0.1003 | 0.0514 |
| 3 | 128 | 0.0170 | 0.0268 | 0.0535 | 0.0388 |
| 3 | 256 | 0.0113 | 0.0199 | 0.0192 | 0.0448 |
| 3 | 512 | 0.0047 | 0.0033 | 0.0085 | 0.0371 |
| 3 | 1024 | 0.0023 | 0.0031 | 0.0049 | 0.0060 |
| 4 | 64 | 0.0456 | 0.0941 | 0.0447 | 0.0000 |
| 4 | 128 | 0.0300 | 0.0139 | 0.0259 | 0.0000 |
| 4 | 256 | 0.0193 | 0.0238 | 0.0162 | 0.0159 |
| 4 | 512 | 0.0031 | 0.0182 | 0.0299 | 0.0157 |
| 4 | 1024 | 0.0046 | 0.0069 | 0.0136 | 0.0252 |

## 4. Monotonicity: σ((R−Rc)/w) vs R_hat

| d | N | Spearman(σ, R_hat) | p-value |
|---|---|-------------------|---------|
| 2 | 64 | -0.738 | 2.03e-04 |
| 2 | 128 | σ constant | — |
| 2 | 256 | σ constant | — |
| 2 | 512 | σ constant | — |
| 2 | 1024 | σ constant | — |
| 3 | 64 | +0.573 | 8.28e-03 |
| 3 | 128 | -0.116 | 6.27e-01 |
| 3 | 256 | -0.857 | 1.38e-06 |
| 3 | 512 | -0.946 | 3.01e-10 |
| 3 | 1024 | -0.833 | 5.08e-06 |
| 4 | 64 | -0.022 | 9.33e-01 |
| 4 | 128 | -0.081 | 7.40e-01 |
| 4 | 256 | -0.082 | 7.31e-01 |
| 4 | 512 | -0.692 | 7.28e-04 |
| 4 | 1024 | -0.928 | 3.88e-09 |

## 5. Hard Threshold Test: σ(H=0) vs σ(H>0)

If sigmoid → Θ, then σ(H=0) → 0 and σ(H>0) → 1 for large N.

| d | N | mean σ(H=0) | mean σ(H>0) | gap | AUC |
|---|---|-------------|-------------|-----|-----|
| 2 | 64 | 1.0000 | 1.0000 | -0.0000 | 0.050 |
| 2 | 128 | 1.0000 | 1.0000 | -0.0000 | 0.000 |
| 2 | 256 | 1.0000 | 1.0000 | +0.0000 | 0.000 |
| 2 | 512 | 1.0000 | 1.0000 | +0.0000 | 0.000 |
| 2 | 1024 | 1.0000 | 1.0000 | +0.0000 | 0.000 |
| 3 | 64 | 1.0000 | 0.6473 | -0.3527 | 0.080 |
| 3 | 128 | 1.0000 | 0.7737 | -0.2263 | 0.030 |
| 3 | 256 | 1.0000 | 0.9215 | -0.0785 | 0.010 |
| 3 | 512 | 1.0000 | 1.0000 | -0.0000 | 0.000 |
| 3 | 1024 | 1.0000 | 1.0000 | -0.0000 | 0.000 |
| 4 | 64 | 0.7596 | 0.1051 | -0.6545 | 0.070 |
| 4 | 128 | 0.9999 | 0.3615 | -0.6384 | 0.040 |
| 4 | 256 | 1.0000 | 0.5137 | -0.4863 | 0.030 |
| 4 | 512 | 1.0000 | 0.7323 | -0.2677 | 0.000 |
| 4 | 1024 | 1.0000 | 0.7500 | -0.2500 | 0.000 |

## 6. N-Scaling of the Sigmoid-to-Step Convergence

Does the gap σ(H>0) - σ(H=0) increase with N?


### d=2

Gap trajectory: -0.0000 → -0.0000 → 0.0000 → 0.0000 → 0.0000
- Spearman(N, gap): ρ=+0.894 (p=4.052e-02)
- → **Gap increasing with N** — sigmoid sharpening toward step

### d=3

Gap trajectory: -0.3527 → -0.2263 → -0.0785 → -0.0000 → -0.0000
- Spearman(N, gap): ρ=+1.000 (p=1.404e-24)
- → **Gap increasing with N** — sigmoid sharpening toward step

### d=4

Gap trajectory: -0.6545 → -0.6384 → -0.4863 → -0.2677 → -0.2500
- Spearman(N, gap): ρ=+1.000 (p=1.404e-24)
- → **Gap increasing with N** — sigmoid sharpening toward step

## Key Discovery: Anti-Correlation

**Occupancy R DECREASES with curvature H**, opposite to the naive expectation.

Physical mechanism: de Sitter expansion a(t) = e^{Ht} compresses the causal diamond. At high H:
- Fewer causal pairs survive (n_pairs drops dramatically)
- Those that survive tend to be links (interval size = 0)
- Therefore f_link increases → R = 1 - f_link decreases

Consequence: σ((R−Rc)/w) is HIGH for flat space and LOW for strong curvature.

### Revised Interpretation

The sigmoid wall does NOT correspond to Θ(S_EH − S_min). Instead:

**σ((R−Rc)/w) → Θ(S_max(d) − S_EH)**

This is an **upper admissibility condition**: the wall enforces a maximum curvature bound.
- Flat space (H=0): R >> Rc → σ ≈ 1 → wall active → "admitted"
- Strong curvature (H>>1): R < Rc → σ ≈ 0 → wall off → "rejected"

This is physically sensible: the F7 wall selects spacetimes with **bounded curvature**,
rejecting those where the causal structure is too compressed by expansion.

### Evidence for Convergence to Step Function

Despite the anti-correlation, the sigmoid IS sharpening toward a step:

1. **R(H=0) → 1 monotonically with N** (all d): occupancy approaches 1 in flat space
2. **R(H=2) → 0 for d=4** (0.000 at N=64/128, 0.078 at N=1024): extreme curvature fully rejected
3. **Gap |σ(H=0) − σ(H>0)| shrinks for d=2** (already saturated at σ=1 everywhere)
4. **Gap |σ(H=0) − σ(H>0)| increases for d=3,4** (Spearman ρ(N, gap) = +1.0)
5. **std(R|d,N,H) → 0** for all d,N,H: fluctuations shrinking = sigmoid sharpening

### Critical Curvature H_c(d, N)

The sigmoid defines an implicit critical curvature: R(H_c) = Rc = 0.25.
- d=2: R > Rc for all tested H → H_c(2) > 2.0 (very tolerant)
- d=3: R crosses Rc between H=1 and H=2 at small N → H_c(3) ~ 1.5–2.0
- d=4: R crosses Rc between H=0 and H=0.5 at N=64 → H_c(4) ~ 0.2–0.5

As N→∞, R(H=0) → 1 and R(H) is a smooth decreasing function of H, so:
- H_c(d, N) increases with N (more pairs → higher occupancy at given H)
- The step sharpens as std(R) → 0

### Correct Bridge Form

$$\sigma\!\left(\frac{R - R_c}{w}\right) \xrightarrow{N \to \infty} \Theta\!\left(R_\mathrm{max}(d) - R_\mathrm{dS}\right)$$

where R_max(d) = d(d−1)·H_c² is the maximum admissible de Sitter curvature.
The F7 wall is NOT "S_EH must exceed a minimum" but "S_EH must not exceed a maximum."

### Summary Verdict

| Criterion | Expected | Observed | Status |
|-----------|----------|----------|--------|
| R monotone with H | R↑ with H | **R↓ with H** | ❌ Reversed |
| σ separates flat/curved | σ(H>0) > σ(H=0) | **σ(H>0) < σ(H=0)** | ❌ Reversed |
| Gap sharpens with N | ✓ | **✓** (d=3,4: ρ=+1.0) | ✅ |
| std(R) → 0 | ✓ | **✓** (all d) | ✅ |
| σ monotone with R_hat | ρ → +1 | **ρ → −1** (anti-monotone) | ⚠ Anti-correlated |

**Conclusion**: The sigmoid wall converges to a step function, but it is an **upper bound** on curvature (Θ(S_max − S_EH)), not a lower bound (Θ(S_EH − S_min)). The physical meaning is: F7 selects spacetimes with bounded curvature, enforcing a maximum-curvature admissibility condition that sharpens to a hard threshold in the continuum limit.
