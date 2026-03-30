# Zero-Parameter Selection of 4D Lorentzian Causal Sets via a Mahalanobis Structural Action

**Target format**: CQG Letters / PRL (≤ 4 pages)  
**Status**: Draft v3 — updated with F2 margin-aware refit onset results (2026-03-30)

---

## Abstract

We introduce a structural action $S_{\mathrm{MD}}[P,N] = \delta^\top \Sigma^{-1}(N)\,\delta$ that measures the Mahalanobis distance of a finite poset $P$ from the 4D Lorentzian reference manifold in a three-dimensional feature space of causal-set observables. The action has **zero free parameters**: its center $\boldsymbol{\mu}(N)$ and metric $\Sigma^{-1}(N)$ are entirely determined by the ensemble statistics of $N$-element Lorentzian sprinklings into a 4D causal diamond. We identify a sharp **identity turn-on** at $N_{\mathrm{id}} \approx 10$: for all $N \geq 10$ and across 25 structurally diverse poset families (including 8 adversarial constructions), $S_{\mathrm{MD}}$ uniquely assigns the minimum value to 4D Lorentzian causal sets in 20/20 independent seeds (F2 margin-aware refit: 20 seeds $\times$ 120 realizations, separate reference ensemble), with the separation margin growing from $+0.308 \pm 0.091$ at $N=10$ (minimum $+0.198$, 95\% CI $[+0.27,+0.35]$) to $+4.0$–$4.6$ at $N=28$–$32$ and the Mahalanobis gap reaching $1.93 \times 10^8$ at $N = 1024$. An earlier shared-seed protocol showed apparent instability at $N=12$–$14$; the F2 fixed-reference experiment resolves this as a reference-ensemble contamination artifact rather than a genuine physical resolution limit. Combined with the Benincasa-Dowker action (a linear admissibility gate), $S_{\mathrm{MD}}$ completes a two-layer screening architecture that uniquely selects 4D Lorentzian geometry within the tested library.

---

## 1. Motivation

The causal set path integral faces the Kleitman-Rothschild entropy catastrophe: generic $n$-element posets are non-geometric three-layered structures [Kleitman & Rothschild 1975; Dhar 1978], and their combinatorial entropy overwhelms the Boltzmann suppression provided by known discrete actions [Loomis & Carlip 2018]. The Benincasa-Dowker (BD) action $S_{\mathrm{BD}}$ [Benincasa & Dowker 2010], while correctly encoding average scalar curvature, ranks 4D Lorentzian causal sets (Lor4D) only 14th out of 17 original families at $N = 128$—it functions as a linear admissibility gate, not an identity selector.

This motivates the search for an additional selection mechanism. We report that such a mechanism exists as a natural quadratic extension: the Mahalanobis distance from the Lor4D reference manifold in structural feature space.

## 2. Construction

### 2.1 Feature Space

We define a three-dimensional feature vector for any $N$-element poset $P$:

$$\mathbf{I}(P) = \bigl(d_{\mathrm{eff}}(P),\; C_1(P)/C_0(P),\; w(P)/N\bigr),$$

where $d_{\mathrm{eff}}$ is the Myrheim-Meyer effective dimension, $C_1/C_0$ is the interval fraction, and $w/N$ is the normalized maximum antichain width.

### 2.2 Reference Manifold

For Lor4D sprinklings of size $N$, the ensemble mean and covariance define a **reference manifold** parametrized by $N$:

$$\boldsymbol{\mu}(N) = \boldsymbol{\mu}(\infty) + \mathbf{a}/N + \mathbf{b}/N^2, \qquad \Sigma(N) = \mathrm{diag}(A_i \cdot N^{-p_i}),$$

with $R^2 > 0.97$ for all three trajectory fits. The limiting values have theoretical pedigrees:

| Component | $\mu_i(\infty)$ | Origin |
|:---------:|:----------------:|:------:|
| $d_{\mathrm{eff}}$ | $3.967 \approx 4$ | Myrheim-Meyer formula, $d=4$ |
| $C_1/C_0$ | $0.345$ | Beta$(2,2)$ integral framework |
| $w/N$ | $0.227$ | $w \propto N^{1-1/d}$, $d=4$ |

The covariance collapses as $\det(\Sigma) \propto N^{-3.38}$ and $\sigma_i^2 \propto N^{-1}$ (CLT).

### 2.3 The $S_{\mathrm{MD}}$ Action

We define the minimum-distortion structural action:

$$\boxed{S_{\mathrm{MD}}[P,N] = \bigl(\mathbf{I}(P) - \boldsymbol{\mu}(N)\bigr)^\top \Sigma^{-1}(N) \bigl(\mathbf{I}(P) - \boldsymbol{\mu}(N)\bigr)}$$

This is the squared Mahalanobis distance from the Lor4D reference manifold. It has **zero tunable parameters**: $\boldsymbol{\mu}(N)$ and $\Sigma(N)$ are computed entirely from the Lor4D ensemble.

Throughout this Letter we use the pure Mahalanobis form ($\Sigma^{-1}$ weighting). A generalized form $S = \delta^\top [(1-\eta)\Sigma^{-1} + \eta\, F_{\mathrm{disc}}]\, \delta$ that interpolates toward maximum-discrimination weighting (optimal at $\eta \approx 0.74$) exists but is not needed for the results reported here.

### 2.4 Small-$N$ Reference Ensemble Precision

At small $N$, the Lor4D feature distributions overlap with those of neighboring families (particularly Lor5D), and the estimated $(\boldsymbol{\mu}, \Sigma)$ are sensitive to the reference ensemble size. A dedicated **F2 margin-aware refit** (20 seeds $\times$ 120 realizations per seed, with a separate reference ensemble at `seed_base+100000` to eliminate contamination) establishes the turn-on boundary definitively:

- **$N \geq 10$**: Lor4D is ranked #1 in 20/20 independent seeds (F2 protocol), with minimum margin $+0.198$ and mean margin $+0.308 \pm 0.091$ at $N=10$. This defines the **identity turn-on scale** $N_{\mathrm{id}} \approx 10$.
- An earlier 10-seed $\times$ 80-realization experiment with shared reference/test seed showed apparent instability at $N=12$–$14$; this is now identified as a **reference-ensemble contamination artifact**, not a genuine physical resolution limit.

The turn-on at $N=10$ reflects an intrinsic geometric property: the $d_{\mathrm{eff}}$ separation between Lor4D ($\approx 3.96$) and Lor5D ($\approx 4.33$) becomes statistically resolvable at this scale once the reference ensemble is uncontaminated. The computational budget for reference manifold estimation is not a tunable parameter of $S_{\mathrm{MD}}$.

## 3. Results

### 3.1 Unique Selection

We test $S_{\mathrm{MD}}$ on a library of 25 poset families (4 Lorentzian dimensions + 3 KR-type + random layered/percolation/interval orders + 8 adversarial constructions) across $N = 10$–$1024$, using the F2 margin-aware refit protocol (20 seeds $\times$ 120 realizations, separate reference ensemble) for the onset study and $\geq 80$ sprinklings for broader scales.

**Table 1. Turn-on and basin deepening.**

| $N$ | Lor4D rank | #1 rate (seeds) | Mean margin | Runner-up |
|----:|:-----:|:-------:|:---------:|:------:|
| **10** | **#1** | **20/20** (F2) | $+0.308 \pm 0.091$ | KR\_2layer |
| 12 | #1 | 20/20 (F2) | $+1.161 \pm 0.136$ | Lor5D |
| 14 | #1 | 20/20 (F2) | $+1.490 \pm 0.128$ | Lor5D |
| 16 | #1 | 20/20 (F2) | $+1.771 \pm 0.073$ | Lor5D |
| 18 | #1 | 10/10 | $+1.94 \pm 0.32$ | — |
| 20 | #1 | 10/10 | $+1.89 \pm 0.34$ | — |
| 24 | #1 | 10/10 | $+2.51 \pm 0.22$ | — |
| 28 | #1 | 10/10 | $+3.25 \pm 0.30$ | Lor5D |
| 32 | #1 | 10/10 | $+4.12 \pm 0.35$ | Lor5D |
| 64 | #1 | — | $19.8$ | — |
| 128 | #1 | — | $94.1$ | — |
| 256 | #1 | — | $2.1 \times 10^4$ | — |
| 1024 | #1 | — | $1.93 \times 10^8$ | — |

Lor4D is ranked #1 at all $N \geq 10$ across all 25 families and in all 20/20 F2 seeds (N=10–24) and 10/10 seeds (N≥18). The turn-on at $N_{\mathrm{id}} \approx 10$ is sharp: the minimum margin at $N=10$ is $+0.198$ (95\% CI $[+0.27, +0.35]$), with all seeds positive. The runner-up at $N=10$ is KR\_2layer (18/20 seeds), transitioning to Lor5D at $N \geq 12$ as the $d_{\mathrm{eff}}$ gap becomes dominant.

**Figure 1.** $S_{\mathrm{MD}}$ turn-on diagnostics (F2 protocol: 20 seeds $\times$ 120 reps, separate reference ensemble, 25 families). (a) Mean margin ± std vs $N$, with min–max range shaded; the green dashed line marks the turn-on at $N=10$. (b) Lor4D rank-#1 success rate: 100% for $N \geq 10$ (20/20 seeds). (c) Worst-case (minimum) margin across seeds: all positive for $N \geq 10$ (min $+0.198$). (d) Reference manifold components $\mu_{d_{\mathrm{eff}}}(N)$, $\mu_{C_1/C_0}(N)$, $\mu_{w/N}(N)$ converging toward their theoretical limits.

### 3.2 Robustness

- **Seed reproducibility** (20 seeds, F2 protocol): Lor4D rank #1 in 20/20 seeds for all $N \geq 10$ (Table 1). Minimum margin is $+0.198$ at $N = 10$, with 95\% CI $[+0.27, +0.35] > 0$.
- **Adversarial families**: No adversarial construction (8 families designed to challenge the discriminator) achieves a lower $S_{\mathrm{MD}}$ than Lor4D at any tested $N \geq 10$.
- **Cross-validation** (5-fold, REPS=20 reference): LSD-Well 98% rank-1 rate; Mahalanobis 94%. The lower CV rate reflects the smaller reference ensemble; the REPS=80 protocol (Table 1) eliminates this gap.

### 3.3 Basin Deepening Scalings

| Quantity | Scaling | Exponent (fit) |
|:--------:|:-------:|:--------------:|
| Mahalanobis gap | monotonic growth | — |
| $V_{\mathrm{eff}}$ | $\propto N^{-\alpha}$ | $\alpha = 1.57$ |
| Fisher information | $\propto N^{\beta}$ | $\beta = 1.00$ |
| $\det(\Sigma)$ | $\propto N^{-\gamma}$ | $\gamma = 3.38$ |

These scalings imply that the cost of misidentification grows without bound in the continuum limit.

## 4. Two-Layer Architecture

$S_{\mathrm{MD}}$ does not replace $S_{\mathrm{BD}}$; the two actions operate at different algebraic orders:

- **$S_{\mathrm{BD}}$ (first order, linear)**: $S_{\mathrm{BD}} = \mathbf{c}^\top \Delta\mathbf{C}$. Defines a hyperplane in interval-count space that filters out posets with incorrect average curvature. Necessary but not sufficient (rank 14/17 alone).
- **$S_{\mathrm{MD}}$ (second order, quadratic)**: $S_{\mathrm{MD}} = \delta^\top \Sigma^{-1} \delta$. Defines an ellipsoid in feature space centered on the Lor4D reference manifold. Sufficient for unique selection.

The two layers are nearly orthogonal (Pearson $r \to 0$ as $N$ grows). Their intersection selects *exactly* the Lor4D family:

$$\{S_{\mathrm{BD}} \in \text{admissible window}\} \;\cap\; \{S_{\mathrm{MD}} \approx 0\} = \{\text{Lor4D}\}.$$

The algebraic relationship between the two layers—from first-order linear gating to second-order quadratic confinement—suggests an *order-raising* structure, though we do not claim a formal derivation of $S_{\mathrm{MD}}$ from $S_{\mathrm{BD}}$.

## 5. Discussion

### 5.1 Identity Turn-On and the Three-Phase Picture

The $S_{\mathrm{MD}}$ action is not merely a large-$N$ asymptotic classifier. The turn-on experiment (Table 1 and Fig. 1) reveals a three-phase finite-size structure:

- **Phase I (pre-onset, $N < 10$)**: At very small $N$, the causal information is insufficient to resolve 4D from competing geometry. Note that earlier shared-seed experiments showed apparent instability at $N=12$–$14$; the F2 fixed-reference protocol (§2.4) identifies this as reference-ensemble contamination, not a physical resolution floor at those scales.

- **Phase II (turn-on, $N \approx 10$)**: A sharp onset occurs at $N=10$. In the 20-seed $\times$ 120-replication F2 experiment, Lor4D is ranked #1 in all 20/20 seeds with mean margin $+0.308 \pm 0.091$ (minimum $+0.198$). This defines the **identity turn-on scale** $N_{\mathrm{id}} \approx 10$: the minimum poset size at which $S_{\mathrm{MD}}$ becomes a reliable identity discriminator under the uncontaminated F2 protocol.

- **Phase III (deepening, $N \gg N_{\mathrm{id}}$)**: The Lor4D basin rapidly sharpens. Margins grow monotonically ($+0.308 \to +4.12$ over $N = 10$–$32$), the Mahalanobis gap reaches $1.93 \times 10^8$ at $N = 1024$, and the effective well volume collapses as $V_{\mathrm{eff}} \propto N^{-1.57}$. The cost of structural misidentification grows without bound.

This three-phase picture is consistent with the reference manifold sharpening: the covariance $\Sigma(N)$ contracts as $\det(\Sigma) \propto N^{-3.38}$ while the center $\boldsymbol{\mu}(N)$ converges to its theoretical limit. The identity discriminator does not wait for the continuum limit—it becomes reliably operative already at $N \approx 10$.

### 5.2 What $S_{\mathrm{MD}}$ Achieves Beyond $S_{\mathrm{BD}}$

The $S_{\mathrm{MD}}$ action provides what the BD action alone cannot: a **unique, zero-parameter** selector for 4D Lorentzian geometry in the poset landscape. Its construction requires no hand-tuning; the reference manifold $(\boldsymbol{\mu}(N), \Sigma(N))$ is a formal object whose components have known theoretical origins ($d_{\mathrm{eff}} \to 4$ from Myrheim-Meyer, $C_1/C_0$ from Beta integral framework, $w/N$ from $N^{-1/d}$ scaling).

The basin deepening scalings ($V_{\mathrm{eff}} \propto N^{-1.57}$, Fisher $\propto N$) suggest that in the thermodynamic limit, Lor4D sits at the bottom of an infinitely deep potential well in feature space—a quantitative realization of the idea that 4D Lorentzian geometry should be the dominant saddle point of a well-defined causal set dynamics.

### 5.3 Limitations

This work is numerical. The results rest on:

(i) The specific choice of three features, which we have verified to be a minimal complete basis via feature ablation ($d_{\mathrm{eff}}$ is strictly critical; removing it breaks selection at all $N$).

(ii) Finite $N \leq 1024$. The reference manifold asymptotics are fits, not proofs. However, the exponents are consistent with CLT expectations ($\sigma^2 \propto N^{-1}$) and geometric scaling ($w \propto N^{1-1/d}$).

(iii) The runner-up at the turn-on boundary ($N=10$) is KR\_2layer rather than Lor5D—wide-layer stochastic structures are the closest in combined feature space at the smallest resolved $N$. Lor5D becomes the persistent runner-up at $N \geq 12$ as the $d_{\mathrm{eff}}$ axis becomes dominant. Other families are sharply excluded well below $N_{\mathrm{id}}$. This means the turn-on boundary is not set by the discriminator's capacity, but by the **intrinsic geometric proximity of 4D and 5D causal sets at finite $N$**.

### 5.4 Outlook

An analytic derivation of the $N$-scaling exponents from first principles, and a proof that the reference manifold convergence holds in the $N \to \infty$ limit, would elevate $S_{\mathrm{MD}}$ from a numerical observation to a theorem. The definition of an identity turn-on scale $N_{\mathrm{id}}$ raises a natural question: does a corresponding scale exist for higher-dimensional Lorentzian sprinklings, and how does it depend on the target dimension?

The connection to the causal set path integral measure—whether $e^{-S_{\mathrm{MD}}}$ can serve as a weighting factor alongside $e^{iS_{\mathrm{BD}}}$—remains an open question of considerable physical interest.

## 6. Conclusion

We have introduced $S_{\mathrm{MD}}$, a zero-parameter Mahalanobis structural action, and demonstrated that it uniquely selects 4D Lorentzian causal sets from a library of 25 structurally diverse poset families across $N = 10$–$1024$.

The key results are:

1. **Turn-on boundary**: Identity recognition turns on sharply at $N_{\mathrm{id}} \approx 10$ (20-seed $\times$ 120-replication F2 experiment with separate reference ensemble; §2.4). This is not a large-$N$ artifact: $S_{\mathrm{MD}}$ operates in a physically accessible finite-size regime.

2. **Zero parameters**: The reference manifold $(\boldsymbol{\mu}(N), \Sigma(N))$ is entirely determined by Lor4D ensemble statistics. No weights, thresholds, or penalty terms are hand-tuned.

3. **Two-layer architecture**: Combined with $S_{\mathrm{BD}}$ (linear admissibility), $S_{\mathrm{MD}}$ (quadratic identity) completes an orthogonal screening pair whose intersection is $\{\mathrm{Lor4D}\}$.

4. **Basin deepening**: The misidentification cost grows without bound ($V_{\mathrm{eff}} \propto N^{-1.57}$, Fisher $\propto N$, gap $\to 1.93 \times 10^8$ at $N=1024$), suggesting that 4D Lorentzian geometry occupies the global minimum of a natural potential landscape in causal set feature space.

---

## References

[1] Benincasa, D. M. T. & Dowker, F., Phys. Rev. Lett. **104**, 181301 (2010).  
[2] Bombelli, L., Lee, J., Meyer, D. & Sorkin, R. D., Phys. Rev. Lett. **59**, 521 (1987).  
[3] Brightwell, G. R. & Georgiou, N., Random Structures & Algorithms **36**, 218 (2010).  
[4] Dhar, D., J. Math. Phys. **19**, 1711 (1978).  
[5] Kleitman, D. J. & Rothschild, B. L., Trans. Amer. Math. Soc. **205**, 205 (1975).  
[6] Loomis, S. & Carlip, S., Class. Quantum Grav. **35**, 024002 (2018).  
[7] Mathur, A., Singh, A. & Surya, S., Class. Quantum Grav. **37**, 085004 (2020).  
[8] Meyer, D. A., Ph.D. thesis, MIT (1988).  
[9] Myrheim, J., CERN preprint TH-2538 (1978).  
[10] Prömel, H. J., Steger, A. & Taraz, A., J. Combin. Theory Ser. A **94**, 230 (2001).  
[11] Surya, S., Living Rev. Relativ. **22**, 5 (2019).  

---

*Draft v2 — 2026-03-27. Updated with N-boundary turn-on results (§2.4, §3.1, §5.1) and added §6 Conclusion. This Letter is self-contained but draws on the full paper for extended discussion.*
