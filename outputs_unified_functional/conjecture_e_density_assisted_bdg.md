# T1: Density-Assisted BDG Calibration — Results


## Design

- Dimensions: [4]
- N values: {4: [128, 256, 512, 1024]}
- H values: [0.0, 0.25, 0.5, 1.0, 2.0]
- Reps per cell: 16
- Total realizations: 320
- Runtime: 337.1s

## 1. Baseline: Raw Feature Correlations with H²

| d | N | Feature | Spearman ρ(feat, H²) | p-value | R²(H²) | R²(H) | R²(H²)/R²(H) |
|---|---|---------|----------------------|---------|---------|-------|---------------|
| 4 | 128 | b1_mean | -0.0670 | 5.99e-01 | 0.0001 | 0.0044 | 0.012 |
| 4 | 128 | b1_std | -0.9487 | 1.06e-32 | 0.6935 | 0.8290 | 0.837 |
| 4 | 128 | b1_mean_scaled | -0.2792 | 2.55e-02 | 0.0590 | 0.0899 | 0.657 |
| 4 | 128 | b1_std_scaled | -0.9623 | 9.08e-37 | 0.5672 | 0.7318 | 0.775 |
| 4 | 128 | s_bdg | -0.0670 | 5.99e-01 | 0.0001 | 0.0044 | 0.012 |
| 4 | 128 | s_bdg_scaled | -0.2792 | 2.55e-02 | 0.0590 | 0.0899 | 0.657 |
| 4 | 256 | b1_mean | +0.0976 | 4.43e-01 | 0.0009 | 0.0015 | 0.595 |
| 4 | 256 | b1_std | -0.9600 | 5.38e-36 | 0.7859 | 0.9036 | 0.870 |
| 4 | 256 | b1_mean_scaled | -0.1649 | 1.93e-01 | 0.0447 | 0.0462 | 0.966 |
| 4 | 256 | b1_std_scaled | -0.9661 | 3.62e-38 | 0.6216 | 0.7825 | 0.794 |
| 4 | 256 | s_bdg | +0.0976 | 4.43e-01 | 0.0009 | 0.0015 | 0.595 |
| 4 | 256 | s_bdg_scaled | -0.1649 | 1.93e-01 | 0.0447 | 0.0462 | 0.966 |
| 4 | 512 | b1_mean | +0.1089 | 3.91e-01 | 0.0058 | 0.0075 | 0.770 |
| 4 | 512 | b1_std | -0.9600 | 5.38e-36 | 0.8727 | 0.9616 | 0.908 |
| 4 | 512 | b1_mean_scaled | -0.1067 | 4.01e-01 | 0.0065 | 0.0053 | 1.235 |
| 4 | 512 | b1_std_scaled | -0.9669 | 1.82e-38 | 0.6964 | 0.8504 | 0.819 |
| 4 | 512 | s_bdg | +0.1089 | 3.91e-01 | 0.0058 | 0.0075 | 0.770 |
| 4 | 512 | s_bdg_scaled | -0.1067 | 4.01e-01 | 0.0065 | 0.0053 | 1.235 |
| 4 | 1024 | b1_mean | -0.0961 | 4.50e-01 | 0.0000 | 0.0001 | 0.400 |
| 4 | 1024 | b1_std | -0.9661 | 3.62e-38 | 0.9044 | 0.9786 | 0.924 |
| 4 | 1024 | b1_mean_scaled | -0.3631 | 3.19e-03 | 0.0395 | 0.0483 | 0.818 |
| 4 | 1024 | b1_std_scaled | -0.9684 | 4.37e-39 | 0.7185 | 0.8704 | 0.826 |
| 4 | 1024 | s_bdg | -0.0961 | 4.50e-01 | 0.0000 | 0.0001 | 0.400 |
| 4 | 1024 | s_bdg_scaled | -0.3631 | 3.19e-03 | 0.0395 | 0.0483 | 0.818 |

## 2. Strategy B: Flat-Baseline Subtraction

Subtract mean(b1_mean) at H=0 for same (d, N), then apply ρ^{2/d}.

| d | N | R_hat_B = -2·ρ^{2/d}·(b1_mean - b1_flat) | Spearman ρ(R_hat_B, H²) | p | R²(H²) | R²(H) | ratio | α_best |
|---|---|-------------------------------------------|-------------------------|---|---------|-------|-------|--------|
| 4 | 128 | b1_flat=0.4082 | +0.1649 | 1.93e-01 | 0.0218 | 0.0402 | 0.543 | 0.25 |
| 4 | 256 | b1_flat=1.9980 | -0.4153 | 6.43e-04 | 0.1351 | 0.1569 | 0.861 | 0.25 |
| 4 | 512 | b1_flat=0.3027 | -0.0121 | 9.24e-01 | 0.0000 | 0.0001 | 0.385 | 0.25 |
| 4 | 1024 | b1_flat=1.1835 | -0.2708 | 3.04e-02 | 0.0306 | 0.0306 | 1.000 | 1.50 |

## 3. Strategy C: OLS Density Residualization + ρ^{2/d}

OLS remove n_causal_pairs from b1_mean, then scale by ρ^{2/d}.

| d | N | Spearman ρ(R_hat_C, H²) | p | R²(H²) | R²(H) | ratio | α_best |
|---|---|-------------------------|---|---------|-------|-------|--------|
| 4 | 128 | +0.0272 | 8.31e-01 | 0.0002 | 0.0003 | 0.547 | 0.25 |
| 4 | 256 | -0.0658 | 6.05e-01 | 0.0000 | 0.0000 | 0.000 | 0.25 |
| 4 | 512 | -0.0416 | 7.44e-01 | 0.0000 | 0.0000 | 3.171 | 2.25 |
| 4 | 1024 | +0.1127 | 3.75e-01 | 0.0000 | 0.0000 | inf | 8.00 |

## 4. Strategy D: b1_std with ρ^{2/d} Scaling

| d | N | Spearman ρ(R_hat_D, H²) | p | R²(H²) | R²(H) | ratio | α_best |
|---|---|-------------------------|---|---------|-------|-------|--------|
| 4 | 128 | -0.9623 | 9.08e-37 | 0.5672 | 0.7318 | 0.775 | 0.25 |
| 4 | 256 | -0.9661 | 3.62e-38 | 0.6216 | 0.7825 | 0.794 | 0.25 |
| 4 | 512 | -0.9669 | 1.82e-38 | 0.6964 | 0.8504 | 0.819 | 0.25 |
| 4 | 1024 | -0.9684 | 4.37e-39 | 0.7185 | 0.8704 | 0.826 | 0.25 |

## 5. Strategy E: b1_std Density-Residualized + ρ^{2/d}

| d | N | Spearman ρ(R_hat_E, H²) | p | R²(H²) | R²(H) | ratio | α_best |
|---|---|-------------------------|---|---------|-------|-------|--------|
| 4 | 128 | +0.0734 | 5.64e-01 | 0.0003 | 0.0006 | 0.431 | 0.25 |
| 4 | 256 | +0.0711 | 5.77e-01 | 0.0021 | 0.0019 | 1.069 | 0.25 |
| 4 | 512 | +0.1263 | 3.20e-01 | 0.0035 | 0.0065 | 0.534 | 0.25 |
| 4 | 1024 | +0.1816 | 1.51e-01 | 0.0057 | 0.0052 | 1.095 | 0.25 |

## 6. Cross-Strategy Comparison at d=4

| N | Strategy | Feature | ρ(H²) | R²(H²) | R²(H) | R²(H²)/R²(H) | α_best |
|---|----------|---------|--------|---------|-------|---------------|--------|
| 128 | A: raw b1_mean | — | -0.0670 | 0.0001 | 0.0044 | 0.012 | 0.25 |
| 128 | A2: ρ^{2/d}·b1_mean | — | -0.2792 | 0.0590 | 0.0899 | 0.657 | 0.25 |
| 128 | B: baseline+ρ | — | +0.1649 | 0.0218 | 0.0402 | 0.543 | 0.25 |
| 128 | C: OLS resid+ρ | — | +0.0272 | 0.0002 | 0.0003 | 0.547 | 0.25 |
| 128 | D: ρ^{2/d}·b1_std | — | -0.9623 | 0.5672 | 0.7318 | 0.775 | 0.25 |
| 128 | E: b1_std resid+ρ | — | +0.0734 | 0.0003 | 0.0006 | 0.431 | 0.25 |
| 128 | F: raw b1_std | — | -0.9487 | 0.6935 | 0.8290 | 0.837 | 0.25 |
| 256 | A: raw b1_mean | — | +0.0976 | 0.0009 | 0.0015 | 0.595 | 0.25 |
| 256 | A2: ρ^{2/d}·b1_mean | — | -0.1649 | 0.0447 | 0.0462 | 0.966 | 1.25 |
| 256 | B: baseline+ρ | — | -0.4153 | 0.1351 | 0.1569 | 0.861 | 0.25 |
| 256 | C: OLS resid+ρ | — | -0.0658 | 0.0000 | 0.0000 | 0.000 | 0.25 |
| 256 | D: ρ^{2/d}·b1_std | — | -0.9661 | 0.6216 | 0.7825 | 0.794 | 0.25 |
| 256 | E: b1_std resid+ρ | — | +0.0711 | 0.0021 | 0.0019 | 1.069 | 0.25 |
| 256 | F: raw b1_std | — | -0.9600 | 0.7859 | 0.9036 | 0.870 | 0.25 |
| 512 | A: raw b1_mean | — | +0.1089 | 0.0058 | 0.0075 | 0.770 | 0.25 |
| 512 | A2: ρ^{2/d}·b1_mean | — | -0.1067 | 0.0065 | 0.0053 | 1.235 | 3.00 |
| 512 | B: baseline+ρ | — | -0.0121 | 0.0000 | 0.0001 | 0.385 | 0.25 |
| 512 | C: OLS resid+ρ | — | -0.0416 | 0.0000 | 0.0000 | 3.171 | 2.25 |
| 512 | D: ρ^{2/d}·b1_std | — | -0.9669 | 0.6964 | 0.8504 | 0.819 | 0.25 |
| 512 | E: b1_std resid+ρ | — | +0.1263 | 0.0035 | 0.0065 | 0.534 | 0.25 |
| 512 | F: raw b1_std | — | -0.9600 | 0.8727 | 0.9616 | 0.908 | 0.50 |
| 1024 | A: raw b1_mean | — | -0.0961 | 0.0000 | 0.0001 | 0.400 | 0.25 |
| 1024 | A2: ρ^{2/d}·b1_mean | — | -0.3631 | 0.0395 | 0.0483 | 0.818 | 0.25 |
| 1024 | B: baseline+ρ | — | -0.2708 | 0.0306 | 0.0306 | 1.000 | 1.50 |
| 1024 | C: OLS resid+ρ | — | +0.1127 | 0.0000 | 0.0000 | inf | 8.00 |
| 1024 | D: ρ^{2/d}·b1_std | — | -0.9684 | 0.7185 | 0.8704 | 0.826 | 0.25 |
| 1024 | E: b1_std resid+ρ | — | +0.1816 | 0.0057 | 0.0052 | 1.095 | 0.25 |
| 1024 | F: raw b1_std | — | -0.9661 | 0.9044 | 0.9786 | 0.924 | 0.75 |

## 7. N-Scaling: R²(H²)/R²(H) Ratio Trend (d=4)

Does any strategy show faster convergence toward H² preference?

| Strategy | N=128 | N=256 | N=512 | N=1024 | ρ(N, ratio) | Trend |
|----------|-------|-------|-------|--------|-------------|-------|
| A: raw b1_mean | 0.012 | 0.595 | 0.770 | 0.400 | +0.40
| B: baseline+ρ | 0.543 | 0.861 | 0.385 | 1.000 | +0.40
| C: OLS resid+ρ | 0.547 | 0.000 | 3.171 | N/A | +0.50
| D: ρ^{2/d}·b1_std | 0.775 | 0.794 | 0.819 | 0.826 | +1.00
| E: b1_std resid+ρ | 0.431 | 1.069 | 0.534 | 1.095 | +0.80
| F: raw b1_std | 0.837 | 0.870 | 0.908 | 0.924 | +1.00

## 8. R_hat / R_dS Convergence (Strategy B, d=4)

If R_hat_B → R_dS, then R_hat_B / R_dS → constant.

| N | H | mean(R_hat_B) | R_dS | ratio | std(ratio) |
|---|---|---------------|------|-------|------------|
| 128 | 0.25 | -5.9700 | 0.7500 | -7.9600 | 14.7801 |
| 128 | 0.5 | 0.0377 | 3.0000 | 0.0126 | 2.1664 |
| 128 | 1.0 | 0.3058 | 12.0000 | 0.0255 | 0.1502 |
| 128 | 2.0 | -0.3300 | 48.0000 | -0.0069 | 0.0016 |
| 256 | 0.25 | 11.3844 | 0.7500 | 15.1792 | 16.8874 |
| 256 | 0.5 | 7.6898 | 3.0000 | 2.5633 | 3.3707 |
| 256 | 1.0 | 4.8200 | 12.0000 | 0.4017 | 0.2063 |
| 256 | 2.0 | 1.2489 | 48.0000 | 0.0260 | 0.0054 |
| 512 | 0.25 | 1.2713 | 0.7500 | 1.6951 | 27.7221 |
| 512 | 0.5 | -2.8573 | 3.0000 | -0.9524 | 5.7609 |
| 512 | 1.0 | -1.2169 | 12.0000 | -0.1014 | 0.6170 |
| 512 | 2.0 | -0.5464 | 48.0000 | -0.0114 | 0.0027 |
| 1024 | 0.25 | 7.8784 | 0.7500 | 10.5045 | 33.2122 |
| 1024 | 0.5 | 8.3970 | 3.0000 | 2.7990 | 6.1668 |
| 1024 | 1.0 | 5.5230 | 12.0000 | 0.4602 | 0.3554 |
| 1024 | 2.0 | 1.2021 | 48.0000 | 0.0250 | 0.0081 |

## 9. Verdict

### Answers to Key Questions:

**Q1: Does ρ^{2/d} scaling improve H² preference?**
**NO.** Adding ρ^{2/d} to b1_std (Strategy D vs F) *reduces* R²(H²) at all N:
- N=128: 0.694→0.567, N=256: 0.786→0.622, N=512: 0.873→0.696, N=1024: 0.904→0.719
- The ρ^{2/d} factor introduces density-dependent noise that dilutes the curvature signal.
- For b1_mean, ρ^{2/d} scaling (A2) slightly improves correlation but remains weak (|ρ| < 0.37).

**Q2: Does flat-baseline subtraction help?**
**NO.** Strategy B (baseline + ρ) gives erratic results:
- ρ(R_hat_B, H²) flips sign across N values (+0.16, −0.42, −0.01, −0.27)
- R² never exceeds 0.14 — far worse than raw b1_std (R² > 0.69)
- The flat baseline b1_flat itself is unstable (0.30–2.00 across N), revealing that
  b1_mean at H=0 has large finite-N fluctuations that swamp the correction.

**Q3: Does any strategy achieve α ≈ 2?**
**NO.** All strategies give α_best ≤ 1.50 for b1_mean-based corrections. The best
α_best from Strategy C at N=512 (2.25) and N=1024 (8.00) are artifacts of near-zero
R² (< 0.001) — noise fitting, not signal. Raw b1_std shows α_best climbing:
0.25→0.25→0.50→0.75, consistent with T5's finding of very slow α growth.

**Q4: Does R_hat_B / R_dS converge?**
**NO.** §8 shows wildly unstable ratios:
- At H=0.25: ratio oscillates −7.96→+15.18→+1.70→+10.50 across N — no convergence
- At H=2.0: ratio is small and negative at some N (−0.007, −0.011) — wrong sign
- std(ratio) >> |mean(ratio)| at all N and H — SNR < 1

**Q5: Which strategy shows fastest N-scaling?**
**F (raw b1_std) and D (ρ^{2/d}·b1_std)** both show R²(H²)/R²(H) ratio
monotonically increasing with ρ(N, ratio) = +1.00:
- F: 0.837→0.870→0.908→0.924
- D: 0.775→0.794→0.819→0.826
But F is strictly better than D at every N — ρ^{2/d} scaling hurts, not helps.

### Overall Assessment

| Strategy | Signal? | α ≈ 2? | Faster than T5? | Verdict |
|----------|---------|--------|-----------------|---------|
| A: raw b1_mean | ❌ No (|ρ|<0.11) | ❌ | — | **FAILED** |
| A2: ρ^{2/d}·b1_mean | ❌ Weak (|ρ|<0.37) | ❌ | — | **FAILED** |
| B: baseline+ρ | ❌ Erratic (sign flips) | ❌ | — | **FAILED** |
| C: OLS resid+ρ | ❌ No (|ρ|<0.12) | ❌ | — | **FAILED** |
| D: ρ^{2/d}·b1_std | ✅ Strong (|ρ|≈0.97) | ❌ | ❌ Worse than F | **Not helpful** |
| E: b1_std resid+ρ | ❌ Destroyed (|ρ|<0.18) | ❌ | — | **FAILED** |
| F: raw b1_std | ✅ Best (|ρ|≈0.97) | ❌ | = Baseline | **Baseline winner** |

### Physical Interpretation — Four-Point Analysis

**1. Strategy A/B/C fail not because of tuning, but because the object is wrong.**
`b1_mean`'s correlation with H² is nearly zero; flat-baseline subtraction and OLS
residualization cannot rescue it. At current N it is noise-dominated, not a usable
curvature proxy. Combined with the completely non-convergent R̂/R_dS calibration (§8):

> **b1_mean cannot be calibrated into an R estimator at these scales.**

**2. The truly stable quantity is raw `b1_std`.**
Its R²(H) is already 0.83→0.98 and continues improving with N; crucially, the
R²(H²)/R²(H) ratio rises monotonically to 0.924, showing that while the dominant
target remains first-order H, second-order sensitivity is slowly strengthening.
This is fully consistent with the established picture: finite-N bulk channels measure
H-like / extrinsic-curvature-trace, and the bridge to R ~ H² is a very slow
asymptotic flow. The proper identity of `b1_std` is therefore not "R estimator" but:

> **The best first-order spectral bulk observable.**

**3. Adding ρ^{2/d} scaling makes things worse — density correction is overcorrection.**
This is significant: it shows that `b1_std`'s bulk information is NOT a simple
"density-contaminated signal waiting to be cleaned." At current sample sizes,
raw `b1_std` is already the optimal effective observable; forcing a density prefactor
only corrupts the primary signal.

**4. OLS density residualization washes out the signal entirely (|ρ| < 0.18).**
This means `b1_std` and density are NOT in a "linearly separable signal + contamination"
relationship. Together with earlier results, the emerging picture is:
- Density is the dominant background
- `b1_std` is a bulk-sensitive effective mode that develops ON this background
- But at finite N, the two are not simply decomposable into "signal + pollution"

Therefore the standard DDT-style linear density removal (which works for raw {C_k})
cannot be applied to `b1_std`. **The spectral bulk mode lives in a nonlinearly
entangled regime with density at finite N.**

### Formal Verdict

> **Verdict.** At accessible finite N, the `b1_mean` route should be closed: neither
> raw regression, flat-baseline subtraction, nor OLS residualization yields a stable
> H²-tracking observable, and the attempted R̂/R_dS calibration is non-convergent.
> By contrast, raw `b1_std` emerges as the unique robust spectral bulk observable:
> it tracks the first-order geometric target H with very high explanatory power,
> and its relative sensitivity to H² increases monotonically with N, consistent
> with an asymptotically real but numerically stiff H→R bridge. Density-prefactor
> scaling and residualization both degrade this signal, indicating that the finite-N
> bulk mode encoded in `b1_std` is not recoverable through simple linear density removal.

> **判决：在当前 finite-N 范围内，`b1_mean` 路线应正式关闭；无论原始拟合、平直基线
> 扣除还是 OLS 去密度，都不能把它变成稳定的 H² 追踪量，R̂/R_dS 标定也不收敛。
> 相反，raw `b1_std` 已成为唯一可靠的谱通道 bulk 观测量：它高精度追踪一阶几何
> 靶标 H，同时对 H² 的相对敏感性随 N 单调增强，符合"H→R bridge 渐近真实但数值
> 僵硬"的总体图景。加入密度前因子或做线性去密度都会削弱这一信号，说明 finite-N
> 上 `b1_std` 承载的 bulk 模态并不能通过简单的线性密度剥离来提纯。**

### Impact on Framework

This result converges the spectral channel from "multiple candidates in parallel" to:

> **Unique main line = raw `b1_std`**

| Observable | Status | Role |
|------------|--------|------|
| `b1_mean` (= S_BDG/N) | **CLOSED** | DDT-trapped, noise-dominated at finite N |
| `b1_std` | **MAIN LINE** | Unique robust first-order spectral bulk observable |
| eigenvalue features | Auxiliary | Beyond-density at d=4 (§4.1.27) but weaker than antichain |

**Spectral bulk is carried by `b1_std`, not `b1_mean`.**

### Confidence Update

| Component | Before T1 | After T1 | Change |
|-----------|-----------|----------|--------|
| E-bulk-second-order | 78–83% | **78–83%** | ±0% |

No confidence change: T1 is a negative result that **eliminates one candidate acceleration
path** but does not change the underlying convergence picture established by T5.

### Updated Priority Path

1. ~~T4 (spectral ratios)~~ — EXCLUDED
2. ~~T5 (N-scaling)~~ — CONFIRMED (α monotonically increases)
3. ~~T1 (density-assisted BDG)~~ — **NEGATIVE** (ρ^{2/d} calibration doesn't help)
4. **T3 (two-step analytical squaring)** — now the last candidate: compute H_hat from
   best first-order channel, then square analytically to get R_hat = d(d-1)·H_hat²
5. **T2 (cross-scale variance)** — test Var(H_hat) ~ H² as an alternative path

---

*Generated by conjecture_e_density_assisted_bdg.py*
*320 realizations (d=4, N=128/256/512/1024, H=0–2, 16 reps), runtime 337.1s*