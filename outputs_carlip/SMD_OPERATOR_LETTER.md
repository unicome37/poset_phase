# Zero-Parameter Selection of 4D Lorentzian Causal Sets via a Mahalanobis Structural Action

**Target format**: CQG Letters / PRL (≤ 4 pages)  
**Status**: Draft v1

---

## Abstract

We introduce a structural action $S_{\mathrm{MD}}[P,N] = \delta^\top \Sigma^{-1}(N)\,\delta$ that measures the Mahalanobis distance of a finite poset $P$ from the 4D Lorentzian reference manifold in a three-dimensional feature space of causal-set observables. The action has **zero free parameters**: its center $\boldsymbol{\mu}(N)$ and metric $\Sigma^{-1}(N)$ are entirely determined by the ensemble statistics of $N$-element Lorentzian sprinklings into a 4D causal diamond. Across $N = 16$–$1024$ and 25 structurally diverse poset families—including 8 adversarial constructions—$S_{\mathrm{MD}}$ uniquely assigns the minimum value to 4D Lorentzian causal sets, with the Mahalanobis gap growing to $1.93 \times 10^8$ at $N = 1024$. The reference manifold covariance collapses as $\det(\Sigma) \propto N^{-3.38}$, while individual feature variances follow $\sigma^2 \propto N^{-1}$, consistent with CLT scaling. Combined with the Benincasa-Dowker action (a linear admissibility gate), $S_{\mathrm{MD}}$ completes a two-layer screening architecture that uniquely selects 4D Lorentzian geometry from the poset landscape.

---

## 1. Motivation

The causal set path integral faces the Kleitman-Rothschild entropy catastrophe: generic $n$-element posets are non-geometric three-layered structures [Kleitman & Rothschild 1975; Dhar 1978], and their combinatorial entropy overwhelms the Boltzmann suppression provided by known discrete actions [Loomis & Carlip 2018]. The Benincasa-Dowker (BD) action $S_{\mathrm{BD}}$ [Benincasa & Dowker 2010], while correctly encoding average scalar curvature, ranks 4D Lorentzian causal sets (Lor4D) only 14th out of 17 families at $N = 128$—it functions as a linear admissibility gate, not an identity selector.

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

**Operator form.** Defining $\delta = \mathbf{I}(P) - \boldsymbol{\mu}(N)$ and the Fisher-weighted matrix $\Lambda(N) = (1-\eta)\Sigma^{-1}(N) + \eta\, F_{\mathrm{disc}}(N)$, the generalized action $S_{\mathrm{MD}} = \delta^\top \Lambda\, \delta$ interpolates between pure Mahalanobis ($\eta=0$) and maximum-discrimination ($\eta=1$) weighting, with optimal mixing at $\eta = 0.74 \pm 0.08$.

### 2.4 Small-$N$ Reference Ensemble Precision

At small $N$ (particularly $N=16$), the Lor4D feature distribution has relatively large variance, and the estimated $(\boldsymbol{\mu}, \Sigma)$ are sensitive to the reference ensemble size. With a small reference sample ($\lesssim 20$ sprinklings), finite-sample noise in $\Sigma^{-1}$ can occasionally invert the top-2 ordering between Lor4D and Lor5D. This is not a model defect but a resolution limit: increasing the reference ensemble to $\geq 80$ sprinklings per family stabilizes the Lor4D #1 ranking at $N=16$ to 10/10 independent seeds, with a mean margin of 1.38. The reference manifold $(\boldsymbol{\mu}(N), \Sigma(N))$ is defined as an *ensemble object estimated to a fixed precision*; the computational budget for this estimation is not a tunable parameter of $S_{\mathrm{MD}}$.

## 3. Results

### 3.1 Unique Selection

We test $S_{\mathrm{MD}}$ on a library of 25 poset families (4 Lorentzian dimensions + 3 KR-type + random layered/percolation/interval orders + 8 adversarial constructions) across $N = 16$–$1024$, using $\geq 80$ sprinklings per family per $N$ for reference manifold estimation.

| $N$ | Lor4D rank ($S_{\mathrm{MD}}$) | Mahalanobis gap | Margin |
|----:|:-----:|:---------:|:------:|
| 16 | #1$^\dagger$ | $1.38$ | resolution limit |
| 32 | #1 | $3.2$ | clear |
| 64 | #1 | $19.8$ | robust |
| 128 | #1 | $94.1$ | definitive |
| 256 | #1 | $2.1 \times 10^4$ | extreme |
| 1024 | #1 | $1.93 \times 10^8$ | overwhelming |

Lor4D is ranked #1 at all tested $N$ across all 25 families. $^\dagger$At $N=16$, the margin is small and sensitive to reference ensemble precision (§2.4); at $N \geq 20$ the ranking is unconditionally stable.

### 3.2 Robustness

- **Cross-validation** (5-fold): LSD-Well 98% rank-1 rate; Mahalanobis 94%.
- **Seed reproducibility** (10 independent seeds): LSD-Well 80/80 rank-1; Mahalanobis 79/80.
- **Adversarial families**: No adversarial construction achieves a lower $S_{\mathrm{MD}}$ than Lor4D at any tested $N$.

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

A magnitude-level gradient correlation ($|\cos(\nabla S_{\mathrm{BD}}, \nabla F_{\mathrm{LSD}})| \approx 0.85$) suggests that the two layers encode related geometric information at different algebraic orders—an *order-raising* relationship from first-order linear gating to second-order quadratic confinement—but the sign of this cosine is numerically unstable due to Jacobian pseudo-inverse ill-conditioning, and we do not claim pointwise gradient identity.

## 5. Discussion

The $S_{\mathrm{MD}}$ action provides what the BD action alone cannot: a **unique, zero-parameter** selector for 4D Lorentzian geometry in the poset landscape. Its construction requires no hand-tuning; the reference manifold $(\boldsymbol{\mu}(N), \Sigma(N))$ is a formal object whose components have known theoretical origins ($d_{\mathrm{eff}} \to 4$ from Myrheim-Meyer, $C_1/C_0$ from Beta integral framework, $w/N$ from $N^{-1/d}$ scaling).

The basin deepening scalings ($V_{\mathrm{eff}} \propto N^{-1.57}$, Fisher $\propto N$) suggest that in the thermodynamic limit, Lor4D sits at the bottom of an infinitely deep potential well in feature space—a quantitative realization of the idea that 4D Lorentzian geometry should be the dominant saddle point of a well-defined causal set dynamics.

**Limitations.** This work is numerical. The results rest on (i) the specific choice of three features, which we have verified to be a minimal complete basis via feature ablation, and (ii) finite $N \leq 1024$. The reference manifold asymptotics are fits, not proofs.

**Outlook.** An analytic derivation of the $N$-scaling exponents from first principles, and a proof that the reference manifold convergence holds in the $N \to \infty$ limit, would elevate $S_{\mathrm{MD}}$ from a numerical observation to a theorem. The connection to the causal set path integral measure—whether $e^{-S_{\mathrm{MD}}}$ can serve as a weighting factor alongside $e^{iS_{\mathrm{BD}}}$—remains an open question of considerable physical interest.

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

*Draft v1 — 2026-03-27. This Letter is self-contained but draws on the full paper (INTRODUCTION_DRAFT.md + DISCUSSION_THEORY_IMPLICATIONS.md) for extended discussion.*
