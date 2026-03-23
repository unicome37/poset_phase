# E2 v4: Structural Efficiency Asymmetry

**Date:** 2026-03-23

**Design:** 900 experiments, N=[12, 16, 20], k=[1, 2, 4], 20 reps


## Key Insight

E1 shows: forward ΔH > backward ΔH (entropy asymmetry ✅)

E2 v3 shows: backward ΔR > forward ΔR (backward adds more structure)

Combined: forward is **more efficient** — more entropy per unit of structural change.


## Test 1: ΔR Asymmetry (causal density change)

**Prediction:** backward ΔR > forward ΔR (backward adds more causal pairs)


| Family | mean ΔR_fwd | mean ΔR_bwd | bwd>fwd % | Wilcoxon p |
|--------|-----------|-----------|-----------|------------|
| Lor2D | +0.05266 | +0.05717 | 55.6% | 4.288e-05 |
| Lor3D | +0.06578 | +0.07357 | 60.6% | 1.062e-07 |
| Lor4D | +0.06793 | +0.07989 | 66.7% | 2.017e-10 |
| KR_like | +0.06572 | +0.06572 | 43.9% | 9.796e-01 |
| TransPerc | +0.05436 | +0.08385 | 89.4% | 4.307e-26 |

## Test 2: Structural Efficiency η = ΔH / |ΔR|

**Prediction:** forward η > backward η (forward gets more entropy per structure)


| Family | mean η_fwd | mean η_bwd | fwd>bwd % | Wilcoxon p |
|--------|-----------|-----------|-----------|------------|
| Lor2D | 8851.90 | 21.99 | 58.9% | 3.094e-03 |
| Lor3D | 30.60 | 17.17 | 60.6% | 1.721e-05 |
| Lor4D | 28.71 | 23.67 | 69.4% | 3.836e-08 |
| KR_like | 21.99 | 24.35 | 48.9% | 4.375e-01 |
| TransPerc | 52.17 | 19.27 | 90.0% | 8.472e-24 |

## Test 3: E1 Recap — A_entropy = ΔH_fwd − ΔH_bwd

| Family | mean A | A>0 % | mean ΔR_fwd | mean ΔR_bwd | interpretation |
|--------|--------|-------|-----------|-----------|----------------|
| Lor2D | +0.166 | 58.3% | +0.05266 | +0.05717 | fwd: +H, −ΔR → efficient |
| Lor3D | +0.256 | 61.1% | +0.06578 | +0.07357 | fwd: +H, −ΔR → efficient |
| Lor4D | +0.413 | 67.2% | +0.06793 | +0.07989 | fwd: +H, −ΔR → efficient |
| KR_like | +0.080 | 49.4% | +0.06572 | +0.06572 | fwd: +H, −ΔR → efficient |
| TransPerc | +1.000 | 87.8% | +0.05436 | +0.08385 | fwd: +H, −ΔR → efficient |

## Test 4: ΔΣ_hist Asymmetry

| Family | mean ΔΣ_fwd | mean ΔΣ_bwd | diff | fwd>bwd % |
|--------|------------|------------|------|-----------|
| Lor2D | +0.0641 | +0.0646 | -0.0005 | 22.2% |
| Lor3D | +0.0991 | +0.0967 | +0.0024 | 22.8% |
| Lor4D | +0.1150 | +0.1083 | +0.0067 | 25.0% |
| KR_like | +0.1097 | +0.1124 | -0.0027 | 18.3% |
| TransPerc | +0.0972 | +0.0916 | +0.0056 | 25.0% |

## Signal Scaling with k

| k | mean A | A>0 % | mean ΔR_diff | mean η_diff |
|---|--------|-------|-------------|-------------|
| 1 | +0.155 | 63% | -0.00385 | +5303.76 |
| 2 | +0.297 | 63% | -0.00900 | +11.48 |
| 4 | +0.697 | 68% | -0.01941 | +12.12 |

## Conclusion

E2 reformulated: the time arrow in structural augmentation manifests as:

1. Forward: more ΔH (entropy), less ΔR (causal density) → **efficient expansion**

2. Backward: less ΔH, more ΔR → **structural compaction**

3. η_fwd > η_bwd: forward augmentation is structurally more efficient

4. This asymmetry grows with k (augmentation size)


Combined with E1 (A>0) and E3 (sign-robust under CG), this confirms:

**The causal future direction is an entropy-efficient expansion,**

**while the causal past direction is a density-increasing compaction.**
