# N=16 C Residual Check v2

**Date:** 2026-03-23

**Purpose:** Test C with TWO definitions under F7 and F8a_v3


## F7

### Definition 1: ALL families at each N (original weight scan)

| N | n | ρ(Σ_hist, F) | p-value | direction |
|---|---|-------------|---------|-----------|
| 16 | 50 | +0.3262 | 2.078e-02 | ❌ |
| 20 | 50 | +0.1276 | 3.773e-01 | ❌ |
| 28 | 50 | -0.5495 | 3.587e-05 | ✅ |
| 36 | 50 | -0.6941 | 2.278e-08 | ✅ |

**C (Def1) = 2/4 = 0.50**

### Definition 2: Lor2D only at each N (f8a_acd style)

| N | n | ρ(Σ_hist, F) | p-value | direction |
|---|---|-------------|---------|-----------|
| 16 | 10 | -0.8200 | 3.675e-03 | ✅ |
| 20 | 10 | -0.3596 | 3.075e-01 | ✅ |
| 28 | 10 | -0.7490 | 1.266e-02 | ✅ |
| 36 | 10 | -0.6039 | 6.450e-02 | ✅ |

**C (Def2) = 4/4 = 1.00**

### F7: N=16 Component Analysis (all families)

| Family | mean(R) | mean(logH) | mean(Σ_hist) | mean(F) |
|--------|---------|-----------|-------------|---------|
| Lor2D | 0.5706 | 15.70 | 0.5250 | 29.1072 |
| Lor3D | 0.2622 | 21.04 | 0.3094 | 29.5827 |
| Lor4D | 0.0793 | 24.22 | 0.2437 | 21.8311 |
| Lor5D | 0.0536 | 26.17 | 0.2062 | 26.6646 |
| KR_like | 0.3178 | 19.07 | 0.2812 | 36.5380 |

**corr(Σ_hist, wall) = +0.9145**, corr(Σ_hist, logH) = -0.8270

**std(wall) = 8.4192**, std(logH) = 3.8603

## F8a_v3

### Definition 1: ALL families at each N (original weight scan)

| N | n | ρ(Σ_hist, F) | p-value | direction |
|---|---|-------------|---------|-----------|
| 16 | 50 | -0.8154 | 5.612e-13 | ✅ |
| 20 | 50 | -0.7592 | 1.659e-10 | ✅ |
| 28 | 50 | -0.7619 | 1.306e-10 | ✅ |
| 36 | 50 | -0.7741 | 4.313e-11 | ✅ |

**C (Def1) = 4/4 = 1.00**

### Definition 2: Lor2D only at each N (f8a_acd style)

| N | n | ρ(Σ_hist, F) | p-value | direction |
|---|---|-------------|---------|-----------|
| 16 | 10 | -0.9535 | 1.940e-05 | ✅ |
| 20 | 10 | -0.7416 | 1.408e-02 | ✅ |
| 28 | 10 | -0.8739 | 9.485e-04 | ✅ |
| 36 | 10 | -0.9217 | 1.496e-04 | ✅ |

**C (Def2) = 4/4 = 1.00**

### F8a_v3: N=16 Component Analysis (all families)

| Family | mean(R) | mean(logH) | mean(Σ_hist) | mean(F) |
|--------|---------|-----------|-------------|---------|
| Lor2D | 0.5706 | 15.70 | 0.5250 | 13.8154 |
| Lor3D | 0.2622 | 21.04 | 0.3094 | 34.2005 |
| Lor4D | 0.0793 | 24.22 | 0.2437 | 102.8454 |
| Lor5D | 0.0536 | 26.17 | 0.2062 | 130.3949 |
| KR_like | 0.3178 | 19.07 | 0.2812 | 20.7193 |


## Higher Power Test (50 reps, all families)

### F7 (50 reps, Def1 = all families)

| N | n | ρ(Σ_hist, F) | p-value | direction |
|---|---|-------------|---------|-----------|
| 16 | 250 | +0.3028 | 1.064e-06 | ❌ |
| 20 | 250 | -0.0406 | 5.232e-01 | ✅ |
| 28 | 250 | -0.4998 | 3.329e-17 | ✅ |
| 36 | 250 | -0.6448 | 8.937e-31 | ✅ |

**C = 3/4 = 0.75**

### F8a_v3 (50 reps, Def1 = all families)

| N | n | ρ(Σ_hist, F) | p-value | direction |
|---|---|-------------|---------|-----------|
| 16 | 250 | -0.8181 | 1.642e-61 | ✅ |
| 20 | 250 | -0.8086 | 4.768e-59 | ✅ |
| 28 | 250 | -0.6903 | 1.019e-36 | ✅ |
| 36 | 250 | -0.7118 | 6.318e-40 | ✅ |

**C = 4/4 = 1.00**


## Conclusion

**Root cause of F7 C16 failure (Definition 1):**

At N=16, the F7 wall term std = 8.42 >> logH std = 3.86, and corr(Σ_hist, wall) = **+0.91**.
The wall's positive correlation with Σ_hist (deeper-layered Lor2D has higher R → higher wall)
overwhelms the negative logH contribution, flipping the net ρ(Σ_hist, F) to +0.33.

This is a **Definition 1 artifact**: it mixes all families at each N, where the wall term
varies by 0–18 across families (Lor2D R≈0.57 → wall≈17.9; Lor4D R≈0.08 → wall≈0).
Definition 2 (Lor2D only) never fails because within Lor2D, wall is nearly constant.

**F8a_v3 resolution:**

F8a divides logH by NlogN, reducing the absolute scale. More importantly, F8a_v3's
two-sided wall + one-sided Ξ_{d,+} redistribute the penalty landscape so that Σ_hist
dominates the within-N ranking regardless of family mixture.

**Final verdict:**
- F7: C = 0.50 (Def1) / 1.00 (Def2) — N=16 AND N=20 fail under Def1
- F8a_v3: **C = 1.00 (both definitions)** — all N pass, including N=16 (ρ = −0.82, p < 10⁻⁶¹)
- **Open problem (6) is RESOLVED**: the N=16 C failure was specific to F7's wall scaling, not a fundamental finite-size limitation
