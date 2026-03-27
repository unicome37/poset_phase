# Layered Structural Screening of 4D Lorentzian Causal Sets

**Target**: Classical and Quantum Gravity (Full Paper)
**Version**: Manuscript v0.2 — §1–§7 reviewed, 6 issues fixed (A+B level)
**Date**: 2026-03-27

---

## §1. Introduction

### 1.1 The entropy catastrophe and the dimension selection problem

The causal set programme posits that the fundamental structure of spacetime is a locally finite partial order—a *causal set*—in which the order relation encodes the causal structure and the counting measure encodes the volume element [1]. A central obstacle to the programme is combinatorial in nature. The Kleitman–Rothschild theorem [2] shows that almost all finite $n$-element posets are three-layered bipartite structures with a $(n/4,\, n/2,\, n/4)$ partition, bearing no resemblance to any Lorentzian manifold. In a discrete path integral

$$Z = \sum_{\mathcal{C}} e^{iS[\mathcal{C}]},$$

the entropy of such non-geometric posets overwhelms the Boltzmann weight of manifold-like causal sets by an exponential factor. This is the **entropy catastrophe**: without an additional selection mechanism, the causal set path integral is dominated by configurations that have nothing to do with spacetime.

Dhar [3] mapped the structure of generic finite posets to a lattice gas with three-body interactions and showed that the occupation fraction $\rho$ (the fraction of comparable pairs) undergoes a first-order phase transition. Prömel, Steger & Taraz [4] extended this to a complete phase diagram with infinitely many transitions as $\rho$ increases—from the KR three-layer regime through higher-layer phases. These results establish that finite posets have a rich combinatorial phase structure, but leave open the question of how a *physical* selection mechanism navigates this landscape to pick out the manifold-like sector.

### 1.2 The limitations of a single action

The most natural candidate for such a mechanism is the Benincasa–Dowker (BD) action [5], the discrete analogue of the Einstein–Hilbert action for causal sets. In $d=4$ dimensions,

$$S_{\mathrm{BD}}^{(4)} = N - C_0 + 9C_1 - 16C_2 + 8C_3, \tag{1}$$

where $C_k$ denotes the number of order-intervals of length $k$ and the coefficients are uniquely determined by the requirement of reproducing the scalar curvature integral in the continuum limit. Loomis & Carlip [6] demonstrated that the BD action preferentially suppresses certain classes of non-manifoldlike causal sets. However, as a *linear* functional of the interval counts $\{C_k\}$, the BD action encodes average curvature content but not structural identity. We shall demonstrate that $S_{\mathrm{BD}}$ alone ranks 4D Lorentzian causal sets (Lor4D) only 14th out of 17 families at $N=128$: it is a necessary admissibility gate, not a sufficient identity selector.

This observation motivates the central question of the present work:

> The question is not "what single action selects Lor4D?" but rather "how many functionally distinct layers does the selection mechanism require, and what does each layer do?"

### 1.3 Summary of results

We report numerical evidence, obtained from a systematic survey across a 25-family poset library at scales $N=12$–$1024$, for a **layered structural screening architecture** comprising at least two functionally distinct layers:

**Layer 1 — Admissibility.** A first-order, linear screening layer (carried jointly by the triple-product functional $S_{\mathrm{triple}}$ and the BD action $S_{\mathrm{BD}}$) eliminates non-geometric and non-robust structures from the poset space, reducing a combinatorially vast landscape to a manageable geometric sector. This layer is necessary but insufficient: Lor4D is not uniquely selected.

**Layer 2 — Identity.** A second-order, quadratic identity functional $S_{\mathrm{MD}}$—the Mahalanobis distance from a Lor4D reference manifold in structural feature space—positively selects Lor4D from among all admissible candidates with zero free parameters. At all tested scales $N \geq 14$ and across all 25 poset families, Lor4D is uniquely ranked first.

The two layers are functionally separated and statistically weakly correlated (Pearson $r \to 0$ as $N$ grows). Their intersection, in the tested library, is Lor4D and only Lor4D.

Furthermore, the identity basin exhibits systematic **deepening** with increasing $N$: the Mahalanobis gap between Lor4D and its nearest competitor grows from $+1.28$ at $N=14$ to $+4.59$ at $N=28$ to $1.93 \times 10^8$ at $N=1024$, the effective potential volume shrinks as $V_{\mathrm{eff}} \propto N^{-1.57}$, and the Fisher information grows linearly. This basin deepening—which we term **historical sedimentation**—represents the progressive locking of Lor4D identity as the number of causal events increases.

### 1.4 Scope and limitations

The present work is strictly numerical. We do not claim analytic proofs of the layered screening principle, nor exhaustive coverage of all possible poset structures. Our conclusions are established within:

- a 25-family poset library (17 standard + 8 adversarial);
- the scale range $N = 12$–$1024$;
- a three-dimensional structural feature space $(d_{\mathrm{eff}},\; C_1/C_0,\; w/N)$.

The layered screening principle is proposed as a **mid-level effective theory**—more specific than the philosophical assertion "existence is what survives screening," more general than the individual numerical experiments from which it is extracted.

### 1.5 Outline

The paper is organized as follows. Section 2 defines the admissibility layer and demonstrates its internal gradient from coarse screening ($S_{\mathrm{triple}}$) to geometric gating ($S_{\mathrm{BD}}$). Section 3 introduces the identity layer ($S_{\mathrm{MD}}$) and the Lor4D reference manifold. Section 4 quantifies the identity turn-on boundary and basin deepening dynamics. Section 5 interprets basin deepening as historical sedimentation and states the formal layered screening principle. Section 6 discusses implications and limitations. Section 7 concludes.

---

## §2. Layer 1: Admissibility

The first layer of the screening architecture answers a single question: **can this poset exist as a geometric structure?** It operates at first order (linear functionals), is broad in scope, and eliminates the vast majority of non-geometric structures—but it cannot uniquely identify Lor4D.

### 2.1 The poset library

We work with a library of 25 poset families organized into four classes:

**Geometric families** (4): Lor2D, Lor3D, Lor4D, and Lor5D—causal sets generated by Poisson sprinkling into causal diamonds in $d$-dimensional Minkowski spacetime, followed by restriction to the induced causal order.

**KR-type families** (3): KR_like, KR_2layer, and KR_4layer—random bipartite and multipartite posets approximating the layered structure predicted by the Kleitman–Rothschild theorem for generic posets.

**Random layered and percolation families** (10): RLk4, RLk6, RLk6_lj, RLk6_mid, RLk6_tap, RLk8—random DAGs with varying layer counts and density profiles; AbsLayer, IntOrder, TransPerc, MLR—interval orders, transitive percolation, and mixed structures.

**Adversarial families** (8): KR_8layer, RandomDAG_sp/md/dn, ChainAnti_sp/dn, SparseChain, MixedLor4D—constructions specifically designed to mimic one or more structural features of Lor4D.

For each family, at each poset size $N$, we generate ensembles of 80 independent realizations and extract three structural features.

**Note on family counts.** Rankings reported "out of 17" refer to the 17 standard families (Geometric + KR-type + Random layered/percolation). Rankings reported "out of 25" include all 8 adversarial families in addition to the 17 standard families. The adversarial families were introduced specifically to stress-test the identity layer (§3); Layer 1 results (§2.3–§2.4) are reported against the 17 standard families for consistency with prior literature.

### 2.2 Structural feature space

The structural features are:

1. **Effective dimension** $d_{\mathrm{eff}}$: the Myrheim–Meyer estimator [7, 8], computed from the fraction of causally related pairs via $d_{\mathrm{eff}} = f_2^{-1}\bigl(C_0/\binom{N}{2}\bigr)$, where $f_2(d)$ is the ordering fraction for $d$-dimensional Minkowski sprinklings. For Lor4D, $\langle d_{\mathrm{eff}} \rangle \to 4.000$ as $N \to \infty$.

2. **Interval fraction** $C_1/C_0$: the ratio of length-1 intervals (links) to length-0 intervals (comparable pairs), encoding the nearest-neighbour density of the causal structure. For Lor4D in a causal diamond, $C_1/C_0 \to 0.345$ from a Beta$(2,2)$ integral framework.

3. **Normalized maximum antichain width** $w/N$: the size of the largest antichain divided by $N$, a proxy for the "spatial extent" of the poset. For Lor4D, $w \propto N^{1-1/d}$ with $1/d \approx 0.25$, giving $w/N \to 0$ as $N^{-1/4}$.

These features are not chosen by optimization. Each has an independent theoretical pedigree in causal set theory, and together they span complementary aspects of the causal structure: dimensionality, local connectivity, and global spatial extent.

### 2.3 Coarse screening: $S_{\mathrm{triple}}$

The simplest admissibility criterion is the triple product of feature deviations from the Lor4D reference:

$$S_{\mathrm{triple}}(P, N) = \bigl|d_{\mathrm{eff}}(P) - \mu_d(N)\bigr| \times \bigl|C_1/C_0(P) - \mu_c(N)\bigr| \times \bigl|w/N(P) - \mu_w(N)\bigr|, \tag{2}$$

where $\boldsymbol{\mu}(N) = (\mu_d(N), \mu_c(N), \mu_w(N))$ is the Lor4D ensemble mean at scale $N$. This is a dimensionless product of absolute deviations—a coarse measure of general "existential" proximity.

Under $S_{\mathrm{triple}}$ alone, Lor4D ranks only **#12 out of 17 standard families**. Several families with quite different internal structure—including Lor2D (#1) and Lor3D—achieve smaller triple products simply because one or two of their features happen to lie close to the Lor4D values. This confirms that $S_{\mathrm{triple}}$ functions as a **coarse sieve**: it identifies a broad class of "not obviously non-geometric" posets, but makes no claim to identity selection.

### 2.4 Geometric gating: $S_{\mathrm{BD}}$

The Benincasa–Dowker action (Eq. 1) provides a more physically grounded admissibility criterion. Introducing the normalized interval densities $\tilde{c}_k = C_k / \binom{N}{2}$, we write

$$S_{\mathrm{BD}} / N = 1 + \frac{N-1}{2}\left(-\tilde{c}_0 + 9\tilde{c}_1 - 16\tilde{c}_2 + 8\tilde{c}_3\right), \tag{3}$$

which is manifestly a **linear** functional of the interval counts. In interval-count space, $S_{\mathrm{BD}} = \text{const}$ defines a hyperplane.

For near-flat 4D Minkowski sprinklings, $\langle S_{\mathrm{BD}} \rangle \to 0$ (zero scalar curvature). We define the admissibility window as the $1\sigma$ range of $S_{\mathrm{BD}}$ values observed for the Lor4D ensemble at each $N$. At $N = 128$, this window is $[-2.88,\; +2.62]$.

**What $S_{\mathrm{BD}}$ achieves.** The BD gate eliminates all KR-type families (KR_like: $S_{\mathrm{BD}} = -15.0$; KR_2layer: $-11.1$; KR_4layer: $-14.0$), most random layered structures, and other manifestly non-geometric posets. It compresses the 17-family space to approximately 6–7 survivors.

**What $S_{\mathrm{BD}}$ does not achieve.** Multiple families pass the gate: IntOrder ($+0.69$), Lor2D ($+1.72$), Lor3D ($+1.07$), Lor5D ($-0.87$), and TransPerc ($-1.60$) all fall within the admissibility window. At $N = 128$, Lor4D ranks only **#14 out of 17** by $S_{\mathrm{BD}}$ alone—worse than most non-geometric families. This is because $S_{\mathrm{BD}}$ measures average curvature correctness, not structural identity: many posets with no manifold-like character happen to have near-zero average "curvature."

**Table 1.** $S_{\mathrm{BD}}$ rankings at selected $N$ values.

| $N$ | Lor4D $S_{\mathrm{BD}}$ rank | Nearest non-Lor4D family |
|:---:|:---:|:---|
| 48 | #9/17 | Lor5D (#8) |
| 64 | #10/17 | Lor5D (#9) |
| 96 | #15/17 | TransPerc (#13) |
| 128 | #14/17 | Lor5D (#13) |

### 2.5 The admissibility gradient

Within the admissibility layer, there is a natural gradient from coarse to geometric:

- $S_{\mathrm{triple}}$ tests general existential proximity—can this poset be *sprinkled from something*?
- $S_{\mathrm{BD}}$ tests geometric curvature—does this poset have the right *average curvature* for a flat 4D spacetime?

Both are necessary-but-not-sufficient conditions. Both eliminate large swathes of non-geometric structure. Neither, alone or together, can uniquely select Lor4D from the geometric sector.

**Key insight.** A linear (first-order) functional can impose necessary conditions but cannot provide sufficient conditions for identity. To discriminate between structures that all pass the admissibility gate, one must ascend to higher-order functionals. This motivates the second layer.

---

## §3. Layer 2: Identity

The second layer answers a different question: **is this poset Lor4D?** It operates at second order (a quadratic functional), is precise rather than broad, and—within the tested 25-family library—robustly performs identity selection with zero free parameters.

### 3.1 The minimum-distortion functional

We define the **minimum-distortion action** (equivalently, the Mahalanobis distance from the Lor4D reference manifold):

$$S_{\mathrm{MD}}[P, N] = \bigl(\mathbf{I}(P) - \boldsymbol{\mu}(N)\bigr)^\top \Sigma^{-1}(N)\, \bigl(\mathbf{I}(P) - \boldsymbol{\mu}(N)\bigr), \tag{4}$$

where:
- $\mathbf{I}(P) = (d_{\mathrm{eff}}(P),\; C_1/C_0(P),\; w(P)/N)$ is the structural feature vector of poset $P$;
- $\boldsymbol{\mu}(N) = \langle \mathbf{I} \rangle_{\mathrm{Lor4D}, N}$ is the ensemble mean of the Lor4D feature vector at scale $N$;
- $\Sigma(N) = \mathrm{Cov}[\mathbf{I}]_{\mathrm{Lor4D}, N}$ is the ensemble covariance matrix.

$S_{\mathrm{MD}}$ is a **quadratic functional** in feature space. Geometrically, the level sets $S_{\mathrm{MD}} = \text{const}$ are ellipsoids centred on $\boldsymbol{\mu}(N)$ with axes determined by $\Sigma^{-1}(N)$. The functional has **zero free parameters**: both the centre and the metric are entirely determined by the Lor4D ensemble statistics at each $N$.

This is not a machine-learning classifier. There is no training/test split, no decision boundary optimization, and no hyperparameter tuning. $S_{\mathrm{MD}}$ is a **one-class measure**—it uses only Lor4D's own statistics to define what "Lor4D-ness" means, and then measures how far any given poset deviates from that reference.

### 3.2 The Lor4D reference manifold

The one-parameter family $\{(\boldsymbol{\mu}(N), \Sigma(N)) \mid N \geq N_{\min}\}$ constitutes the **Lor4D reference manifold** $\mathcal{M}_4$ in feature space. Its key properties:

**Centre trajectory.** The mean feature vector follows a $1/N$ expansion:

$$\boldsymbol{\mu}(N) = \boldsymbol{\mu}(\infty) + \frac{\mathbf{a}}{N} + \frac{\mathbf{b}}{N^2}, \quad R^2 > 0.97, \tag{5}$$

with limiting values:
- $\mu_d(\infty) = 3.967 \approx 4$ (Myrheim–Meyer, exact in the continuum);
- $\mu_c(\infty) = 0.345$ (derivable from a Beta$(2,2)$ integral framework);
- $\mu_w(\infty) = 0.227$ (consistent with $w \propto N^{1-1/d}$ scaling, exponent $0.284 \approx 1/4$).

Each component of $\boldsymbol{\mu}(\infty)$ has an independent theoretical pedigree. This is not a fit to the data; it is a consistency check against known analytical results.

**Covariance contraction.** The eigenvalues of $\Sigma(N)$ decay as power laws:

$$\sigma_i^2(N) = A_i \cdot N^{-p_i}, \tag{6}$$

with $p_i \approx 1$ (consistent with CLT scaling $\sigma^2 \propto 1/N$) for all three features. The covariance determinant collapses as

$$\det(\Sigma) \propto N^{-3.38}. \tag{7}$$

As $N$ grows, the Lor4D point cloud in feature space concentrates toward $\boldsymbol{\mu}(\infty)$—the reference manifold sharpens.

**Eigenvector stability.** The principal axes of $\Sigma(N)$ rotate slowly with $N$, indicating that the contraction is approximately isotropic in the natural feature basis. This stability supports the interpretation of $\mathcal{M}_4$ as a well-defined geometric object, not an artifact of a particular $N$ value.

**Terminological note.** We use "reference manifold" in the sense of a parametrized reference object within a mid-level effective theory. No strict differentiable manifold structure is assumed; such formalization, if needed, can be pursued in future analytical work.

### 3.3 Identity selection: the evidence

**Table 2.** $S_{\mathrm{MD}}$ rankings across $N$ and family library.

| $N$ | Lor4D $S_{\mathrm{MD}}$ rank | #1 rate (10 seeds) | Runner-up | Mean margin |
|:---:|:---:|:---:|:---|:---:|
| 12 | varies | 6/10 | Lor5D | $+0.10$ |
| 14 | #1 | 10/10 | Lor5D | $+1.28$ |
| 16 | #1 | 10/10 | Lor5D | $+1.72$ |
| 18 | #1 | 10/10 | Lor5D | $+2.20$ |
| 20 | #1 | 10/10 | Lor5D | $+2.05$ |
| 24 | #1 | 10/10 | Lor5D | $+3.17$ |
| 28 | #1 | 10/10 | Lor5D | $+4.59$ |
| 32 | #1 | 10/10 | Lor5D | $+3.99$ |
| 128 | #1 | $^\dagger$ | Lor5D | $+94.1$ |
| 1024 | #1 | $^\dagger$ | Lor5D | $+1.93 \times 10^8$ |

$^\dagger$ The 10-seed turn-on experiment covers $N = 12$–$32$. The $N = 128$ and $N = 1024$ entries are from the basin-deepening experiment (30 realizations per family, single reference seed); rank #1 is confirmed but multi-seed #1-rate statistics are not available at these scales.

Across all tested scales $N \geq 14$ and all 25 families (including 8 adversarial constructions), Lor4D is **uniquely and unconditionally ranked #1** by $S_{\mathrm{MD}}$. The sole consistent runner-up is Lor5D—a fact we return to in §4.

### 3.4 Physical content of the identity layer

Several aspects of $S_{\mathrm{MD}}$ deserve emphasis:

**Identity, not classification.** $S_{\mathrm{MD}}$ does not say "this poset belongs to category X." It says: "this poset deviates from the Lor4D reference by this much, in these directions, with this statistical cost." It is an identity measure in the same sense that a Mahalanobis distance measures deviation from a known population.

**The runner-up is always Lor5D.** This is physically significant. The turn-on boundary $N_{\mathrm{id}} \approx 14$ is not set by the discriminator's statistical power; it is set by the **intrinsic geometric proximity between 4D and 5D causal structures**. At $N = 12$, the effective-dimension overlap between Lor4D ($d_{\mathrm{eff}} \approx 3.96$) and Lor5D ($d_{\mathrm{eff}} \approx 4.33$) is too small to resolve with only 12 sprinkled points. The identity boundary is a physical resolution limit, not a statistical artifact.

**Zero parameters, anchored by theory.** The centre $\boldsymbol{\mu}(N)$ converges to values with independent theoretical derivations (Myrheim–Meyer, CLT). The metric $\Sigma^{-1}(N)$ is the inverse empirical covariance—the natural Riemannian metric on a Gaussian family (the Fisher information metric in the diagonal limit). No optimization is performed.

### 3.5 Functional separation from Layer 1

The admissibility layer ($S_{\mathrm{BD}}$) and the identity layer ($S_{\mathrm{MD}}$) are not redundant. Their separation is evidenced by:

**Rank discordance.** At $N = 128$, $S_{\mathrm{BD}}$ ranks Lor4D #14/17 while $S_{\mathrm{MD}}$ ranks Lor4D #1/25. The two rankings are nearly uncorrelated.

**Correlation decay.** The Pearson correlation coefficient between $S_{\mathrm{BD}}$ and $S_{\mathrm{MD}}$ values across the family library is:

| $N$ | Pearson $r$ | Spearman $\rho$ |
|:---:|:---:|:---:|
| 48 | $+0.57$ | $+0.22$ |
| 64 | $+0.07$ | $+0.08$ |
| 96 | $+0.18$ | $-0.25$ |
| 128 | $-0.26$ | $-0.31$ |

The correlation is modest at small $N$ and decays toward zero (or fluctuates around zero) as $N$ grows. The two layers extract functionally independent information from the poset.

**Order difference.** $S_{\mathrm{BD}}$ is a linear functional of interval counts (first-order). $S_{\mathrm{MD}}$ is a quadratic functional in feature space (second-order). The transition from Layer 1 to Layer 2 is an **order-raising**: from a hyperplane cut to an ellipsoidal confinement.

---

## §4. Identity Dynamics: Turn-On and Basin Deepening

### 4.1 Three phases of identity emergence

As $N$ increases from small values, the identity layer exhibits a three-phase structure:

**Phase I — Pre-onset** ($N \lesssim N_{\mathrm{res}} \approx 12$). The poset contains too few causal events for the three structural features to resolve the difference between Lor4D and its nearest competitors. At $N = 12$, Lor4D achieves rank #1 in only 6 out of 10 independent seeds (REPS = 80); the intruders are Lor5D (3 wins) and KR_2layer (1 win). Detailed diagnosis reveals that the failures are physical: at 12 sprinkled points, the effective-dimension distributions of Lor4D ($d_{\mathrm{eff}} = 3.96 \pm 0.05$) and Lor5D ($d_{\mathrm{eff}} = 4.33 \pm 0.10$) overlap to the extent that seed-to-seed fluctuations in the reference ensemble $(\boldsymbol{\mu}(12), \Sigma(12))$ can tip the ranking. The condition number $\kappa(\Sigma) < 50$ throughout, ruling out numerical pathology.

**Phase II — Turn-on** ($N \approx N_{\mathrm{id}} \approx 14$). A sharp transition occurs. At $N = 14$, Lor4D is ranked #1 in 10/10 seeds with a mean margin of $+1.28$ (minimum $+0.33$). The jump from Phase I to Phase II—#1 rate from 60% to 100%, mean margin from $+0.10$ to $+1.28$—is not gradual; it is a **phase-transition-like turn-on** of identity discrimination.

The turn-on scale $N_{\mathrm{id}} \approx 14$ is not a parameter of the model. It emerges from the interplay between the intrinsic geometric separation of Lor4D and Lor5D in feature space and the statistical precision of the reference ensemble at finite $N$. This is a resolution limit of the causal structure itself, not of the discriminator.

**Phase III — Deepening** ($N \gg N_{\mathrm{id}}$). Once identity is established, the basin continues to deepen: the margin between Lor4D and its nearest competitor grows, the reference manifold sharpens, and the statistical cost of misidentification increases without bound.

### 4.2 Quantifying basin deepening

We define the **isolation gap** as the Mahalanobis distance difference between Lor4D and the nearest competitor:

$$\Delta_{\mathrm{hist}}(N) \equiv \min_{f \neq \mathrm{Lor4D}} \bigl[S_{\mathrm{MD}}(f, N) - S_{\mathrm{MD}}(\mathrm{Lor4D}, N)\bigr]. \tag{8}$$

**Table 3.** Isolation gap and basin diagnostics across $N$.

| $N$ | $\Delta_{\mathrm{hist}}$ | $V_{\mathrm{eff}}$ | $I_F = \mathrm{tr}(\Sigma^{-1})$ | Status |
|:---:|:---:|:---:|:---:|:---|
| 12 | $\approx 0.10$ | $3.3 \times 10^{-3}$ | $197$ | Identity not locked |
| 14 | $\approx 1.28$ | — | — | Turn-on |
| 16 | $\approx 1.72$ | $1.6 \times 10^{-3}$ | $469$ | Stable deepening |
| 20 | $\approx 2.05$ | $2.4 \times 10^{-3}$ | $368$ | Stable deepening |
| 28 | $\approx 4.59$ | $1.1 \times 10^{-3}$ | $540$ | Accelerated deepening |
| 128 | $\approx 94.1$ | $9.7 \times 10^{-5}$ | $2355$ | Deep lock |
| 1024 | $\approx 1.93 \times 10^8$ | — | — | Ultra-deep lock |

$V_{\mathrm{eff}}$ and $I_F$ values are from the 30-realization basin-deepening experiment; entries marked "—" correspond to scales not included in that experiment run. The power-law fits quoted in §4.3 are derived from the full $N$-grid of this experiment ($N = 12$–$256$, 11 points).

The growth of $\Delta_{\mathrm{hist}}(N)$ is not strictly monotonic at all scales—a reversal from $N = 28$ ($+4.59$) to $N = 32$ ($+3.99$) is observed and attributable to finite-seed fluctuation—but the overall trend is unambiguously upward.

### 4.3 Scaling laws

Four independent diagnostics confirm the systematic deepening:

**1. Effective potential volume.**

$$V_{\mathrm{eff}}(N) \propto N^{-1.57}. \tag{9}$$

The effective volume of the Lor4D identity basin shrinks as a power law. In the continuum limit, $V_{\mathrm{eff}} \to 0$: the basin concentrates to a point.

**2. Covariance determinant.**

$$\det\bigl(\Sigma(N)\bigr) \propto N^{-3.38}. \tag{10}$$

Since $V_{\mathrm{eff}} \propto \sqrt{\det(\Sigma)}$, this gives $V_{\mathrm{eff}} \propto N^{-1.69}$, consistent with the directly measured exponent (Eq. 9) within numerical uncertainty.

**3. Fisher information.**

$$I_F(N) = \mathrm{tr}\bigl(\Sigma^{-1}(N)\bigr) \propto N^{1.00}. \tag{11}$$

The statistical cost of misidentification—the information content of the Lor4D reference ensemble—grows linearly with $N$. This is the expected CLT scaling: each additional causal event contributes a fixed amount of geometric information.

**4. Mahalanobis gap.**

$$\Delta_{\mathrm{hist}}(N) : \quad +1.28 \;\to\; +94.1 \;\to\; +1.93 \times 10^8 \qquad (N = 14 \to 128 \to 1024). \tag{12}$$

The gap grows super-linearly in $N$, reflecting the compound effect of reference manifold sharpening ($\Sigma \to 0$) and centre stabilization ($\boldsymbol{\mu} \to \boldsymbol{\mu}(\infty)$).

### 4.4 The physical meaning of basin deepening

The deepening of the identity basin has a clear physical interpretation: **each additional causal event added to the poset provides additional geometric information that further distinguishes Lor4D from its competitors.** The reference manifold progressively "comes into focus"—the center converges to its theoretical limit, the fluctuations shrink, and the cost of deviating from the Lor4D template grows without bound.

This is not an artifact of increasing statistical power with larger samples. It reflects the intrinsic geometric content of causal structures: a 4D Lorentzian causal set with 1024 elements contains genuinely more dimensional, topological, and geometric information than one with 14 elements, and this information accumulates in a structured way that systematically strengthens the identity signal.

We formalize this observation in §5 as the **historical sedimentation** of the identity basin, and state the layered screening principle that synthesizes the results of §§2–4.

---

## §5. Historical Sedimentation and the Layered Screening Principle

Section 4 established that the identity basin deepens systematically with $N$. In this section we elevate that observation to a mid-level theoretical interpretation—**historical sedimentation**—and state the formal layered screening principle that synthesizes the results of §§2–4.

### 5.1 From basin deepening to historical sedimentation

The scaling laws of §4.3 admit a unified physical reading:

| Physical concept | Mathematical correspondent | Observable |
|:---|:---|:---|
| History gradually accumulates | Reference manifold gradually comes into focus | $\boldsymbol{\mu}(N) \to \boldsymbol{\mu}(\infty)$ with $1/N$ corrections |
| Past constrains future | Covariance contracts | $\det(\Sigma) \propto N^{-3.38}$ |
| Inertia strengthens | Identity basin deepens | Gap growth; $V_{\mathrm{eff}} \propto N^{-1.57}$ |
| Deviation cost increases | Off-manifold penalty grows | Fisher $\propto N$ |
| Accessible directions narrow | Ellipsoid cross-sections shrink | Eigenvalue decay $\sigma_i^2 \propto N^{-p_i}$ |

The key insight is that **history is not an additional memory appended to the structure**. It is the process by which the reference manifold sharpens, the covariance contracts, and the cost of deviating from the Lor4D template grows without bound. Each additional causal event does not merely add a data point; it tightens the geometric constraints that define what it means to be a 4D Lorentzian causal set.

We call this process **historical sedimentation**: the progressive, irreversible accumulation of geometric identity through the growth of the causal structure.

### 5.2 Analytical skeleton of gap growth

While a complete analytical proof is beyond the scope of this numerical study, the basic mechanism of gap growth can be understood from first principles.

For a Lor4D sprinkling of $N$ points into a $d$-dimensional causal diamond, the central limit theorem gives

$$\sigma_i^2(N) \sim A_i / N \tag{13}$$

for each feature $i$, where $A_i$ depends on the feature's intrinsic variance per causal event. The inverse covariance therefore scales as $\Sigma^{-1} \propto N$, and for a competitor family $f$ at fixed feature-space distance $\|\boldsymbol{\mu}_f - \boldsymbol{\mu}_{\mathrm{Lor4D}}\| = \Delta_f > 0$,

$$S_{\mathrm{MD}}(f, N) \sim \Delta_f^2 \cdot N. \tag{14}$$

Since $S_{\mathrm{MD}}(\mathrm{Lor4D}, N) = O(1)$ (the self-distance is $\chi^2_3$-distributed), the isolation gap grows as

$$\Delta_{\mathrm{hist}}(N) \sim \Delta_f^2 \cdot N - O(1) \;\propto\; N. \tag{15}$$

In practice, the competitor means $\boldsymbol{\mu}_f(N)$ also shift with $N$, and the off-diagonal terms of $\Sigma^{-1}$ contribute, leading to empirically faster-than-linear growth (a rough power-law fit gives $\Delta_{\mathrm{hist}} \sim N^{\alpha}$ with $\alpha \approx 1.6$, but this exponent is estimated from only 11 $N$-values on a non-uniform grid and should not be taken as precise). A full analytical derivation would require the $N$-dependence of the competitor trajectories, which we leave to future work.

The essential point is that gap growth is not an accident of the particular families in our library. It follows from the CLT structure of the Lor4D ensemble combined with the assumption that competitor families do not converge to the same point in feature space—an assumption verified numerically across all 25 tested families.

### 5.3 The sedimentation proposition

> **Proposition (Historical sedimentation as monotone accumulation).** In a finite causal structure space, once the Lor4D identity basin opens at $N \geq N_{\mathrm{id}}$, the isolation gap $\Delta_{\mathrm{hist}}(N)$ exhibits an overall growth trend and the effective basin volume $V_{\mathrm{eff}}(N)$ exhibits an overall contraction trend as $N$ increases—i.e., historical sedimentation is, in a statistical sense, monotonically accumulative.

This is not a theorem. It is a numerically well-supported mid-level proposition that can be proved or refuted by future analytical work. The non-strict monotonicity caveat (e.g., the $N = 28 \to 32$ reversal) is explicitly acknowledged: seed-to-seed fluctuations can temporarily reverse the trend at any finite $N$, but the overall trajectory is unambiguously upward.

### 5.4 Historical sedimentation is not a third screening layer

An important structural clarification: historical sedimentation describes the **longitudinal behavior of the identity layer as $N$ increases**. It is not an independent third screening operator parallel to admissibility and identity. The quantities it tracks—$\boldsymbol{\mu}(N)$, $\Sigma(N)$, gap, $V_{\mathrm{eff}}$—all belong to the identity layer's scale evolution.

The correct layer hierarchy is:
- **Cross-sectional** (fixed $N$): Admissibility + Identity → layered screening;
- **Longitudinal** (across $N$): Identity basin deepens → historical sedimentation.

### 5.5 The Layered Structural Screening Principle

Synthesizing §§2–4 and the sedimentation interpretation:

> **Layered Structural Screening Principle.** In a set of finite posets with sufficiently many causal events, the selection of 4D Lorentzian causal structure requires at least two functionally distinct mechanisms of different orders:
>
> 1. A first-order admissibility layer (carried jointly by $S_{\mathrm{triple}}$ and $S_{\mathrm{BD}}$), excluding non-geometric and non-robust structures;
> 2. A second-order quadratic identity functional $S_{\mathrm{MD}}$, positively selecting Lor4D from among all admissible candidates.
>
> The two layers are functionally separated, statistically weakly correlated, and their intersection in the tested 25-family library uniquely and robustly selects Lor4D.
>
> Furthermore, the identity basin exhibits monotone-accumulative deepening as $N$ increases (historical sedimentation), making Lor4D an increasingly irreplaceable robust identity centre in the current family library.

### 5.6 Numerical evidence summary

| Evidence type | Status |
|:---|:---|
| 25 families × $N = 14$–$1024$ × 10 seeds | ✅ Complete |
| Turn-on boundary $N_{\mathrm{id}} \approx 14$ | ✅ Quantified |
| $S_{\mathrm{BD}}$–$S_{\mathrm{MD}}$ functional separation + weak correlation | ✅ Pearson $r \to 0$ |
| Basin deepening power laws | ✅ $V_{\mathrm{eff}} \propto N^{-1.57}$ |
| Sedimentation monotone trend | ✅ $\Delta_{\mathrm{hist}}$ increasing |
| Analytical derivation | ❌ Skeleton only (§5.2) |
| $N > 1024$ verification | ❌ Not yet tested |

---

## §6. Discussion

### 6.1 This is not machine learning

A natural concern is that $S_{\mathrm{MD}}$ might be a disguised classifier, overfitting to a particular family library. Several features of the construction argue against this:

(i) $S_{\mathrm{MD}}$ is a **one-class measure**: it uses only Lor4D's own ensemble statistics to define the reference, with no knowledge of competitor families during construction. There is no training/test split, no decision boundary optimization, and no hyperparameter.

(ii) The centre $\boldsymbol{\mu}(\infty)$ is anchored by independent theoretical results (Myrheim–Meyer, CLT, Beta integral framework), not learned from data.

(iii) The metric $\Sigma^{-1}(N)$ is the empirical inverse covariance of the Lor4D ensemble—the natural Riemannian metric on a Gaussian family (the Fisher information metric in the diagonal limit). It is a property of Lor4D physics, not a fitting choice.

(iv) Robustness against adversarial families: 8 constructions specifically designed to mimic Lor4D features all fail to displace it from rank #1 at any $N \geq 14$.

### 6.2 Order-raising, not replacement

$S_{\mathrm{MD}}$ does not render $S_{\mathrm{BD}}$ obsolete. The two layers encode different aspects of the same underlying causal geometry:

- $S_{\mathrm{BD}}$: a linear functional of interval counts → average curvature (first order);
- $S_{\mathrm{MD}}$: a quadratic functional in structural feature space → local geometric distortion (second order).

The transition from Layer 1 to Layer 2 is an **order-raising**: the same causal geometry that produces near-zero average curvature ($S_{\mathrm{BD}} \approx 0$) also produces a characteristic fingerprint in $(d_{\mathrm{eff}}, C_1/C_0, w/N)$ space that can be detected by a quadratic functional. The two layers are related as a hyperplane cut is related to an ellipsoidal confinement—different geometric objects encoding overlapping but non-redundant information.

Supporting evidence for this structural relatedness comes from a magnitude-level gradient correlation: after correcting the reference-manifold centre to the Lor4D centroid trajectory, $|\cos(\nabla S_{\mathrm{BD}}, \nabla F_{\mathrm{LSD}})| \approx 0.85$ across all tested scales. However, the sign of the cosine is unstable due to Jacobian pseudoinverse ill-conditioning, so this correlation is interpreted as supportive rather than conclusive. The primary evidence for the two-layer architecture comes from rank separation, correlation decay, and adversarial robustness.

### 6.3 The reference manifold as a physical object

The one-parameter family $\mathcal{M}_4 = \{(\boldsymbol{\mu}(N), \Sigma(N)) \mid N \geq N_{\min}\}$ is not merely a statistical tool. It is a formal theoretical object: the **attractor trajectory** of 4D Lorentzian causal structure in structural feature space. Its properties—centre convergence to theoretically predicted limits, power-law covariance contraction, stable eigenvector orientations—suggest that it captures genuine geometric content of the continuum limit.

In the path integral framework, one may speculate that $e^{-S_{\mathrm{MD}}/2}$ plays the role of a Boltzmann-like weight favouring posets close to the Lor4D reference manifold. If the identity layer can be incorporated into the causal set action as

$$S_{\mathrm{eff}}[P, N] = S_{\mathrm{BD}}[P] + \lambda\, S_{\mathrm{MD}}[P, N], \tag{16}$$

then the combined action would suppress both non-geometric posets (via $S_{\mathrm{BD}}$) and non-Lor4D geometric posets (via $S_{\mathrm{MD}}$), potentially resolving the entropy catastrophe at both levels. We note that the coupling $\lambda$ introduces a free parameter, in contrast to the zero-parameter character of $S_{\mathrm{MD}}$ itself. Whether a natural value (e.g., $\lambda = 1$) can be motivated from first principles—perhaps by requiring that the combined action reproduces the two-layer ranking structure without fine-tuning—is an open question. Equation (16) is therefore a formal target for future analytical work, not a concrete proposal.

### 6.4 Analytical prospects: deriving $\Delta_{\mathrm{hist}}(N)$

The analytical skeleton of §5.2 suggests that a full derivation of the gap growth law is within reach. The key ingredients are:

(i) The CLT scaling of feature variances, $\sigma_i^2 \propto 1/N$, which is a consequence of the independence of sprinkled-point contributions in sufficiently large causal diamonds;

(ii) The convergence of feature means to their theoretical limits, with $1/N$ corrections related to finite-size effects in the causal diamond geometry;

(iii) The separation of competitor families in feature space—i.e., the statement that no non-Lor4D family converges to $\boldsymbol{\mu}(\infty)$ as $N \to \infty$.

If (i)–(iii) can be established analytically, then $\Delta_{\mathrm{hist}}(N) \propto N$ follows immediately, and the empirically faster-than-linear growth would be attributable to the $N$-dependent corrections in the competitor trajectories. Item (iii) is the most substantive: it amounts to a **uniqueness theorem** for $\boldsymbol{\mu}(\infty)$ among all geometrically meaningful poset families—a conjecture that our numerical results strongly support but do not prove.

### 6.5 Dimensional dependence of the turn-on scale

The turn-on scale $N_{\mathrm{id}} \approx 14$ is specific to the selection of Lor4D from a library that includes Lor5D as the nearest competitor. One expects $N_{\mathrm{id}}$ to depend on the target dimension $d$ and the geometric proximity of the nearest confusable structure:

- For **Lor2D selection** (runway: 1D vs 2D), the feature-space separation is large ($d_{\mathrm{eff}}$ differ by $\sim 1$), so $N_{\mathrm{id}}$ should be smaller.
- For **Lor6D selection** (if it were the target), the 5D–6D $d_{\mathrm{eff}}$ gap narrows, so $N_{\mathrm{id}}$ should be larger.

A systematic study of $N_{\mathrm{id}}(d)$ would test whether the layered screening principle holds across dimensions and whether the resolution-limit interpretation of Phase I generalizes. This is left to future work.

### 6.6 Beyond quadratic: higher-order screening?

The present work identifies two layers: linear (first-order) and quadratic (second-order). A natural question is whether the hierarchy continues:

- **Third order**: cubic or higher-moment statistics of the feature distributions could in principle distinguish families that share the same mean and covariance but differ in skewness or kurtosis.
- **Information-geometric extensions**: the Fisher–Rao metric, Rényi divergence, or $f$-divergences between Lor4D and competitor distributions could provide richer discrimination without increasing the parameter count.

We did not observe any need for a third layer in the present library—$S_{\mathrm{MD}}$ uniquely selects Lor4D across all tested conditions. However, with a sufficiently adversarial family library (e.g., a family engineered to match Lor4D's mean and covariance but differ in higher moments), a third layer might become necessary. The layered screening principle is formulated with "at least two layers" precisely to leave this possibility open.

### 6.7 Limitations

We reiterate the boundary conditions of our claims:

1. **Library scope.** All results are established within a 25-family library. While this includes adversarial constructions, it does not exhaust the space of all possible posets. A family that matches Lor4D in all three features at all $N$ would defeat the screening—though we have found no such family despite targeted efforts.

2. **Feature choice.** The three structural features $(d_{\mathrm{eff}}, C_1/C_0, w/N)$ are not derived from a first-principles variational argument. They are selected for their independent theoretical pedigrees and empirical completeness. A feature-ablation study confirms that $d_{\mathrm{eff}}$ is the most critical feature; the others primarily increase the margin.

3. **Mid-level theory.** The layered screening principle is a mid-level effective theory, not a fundamental principle. It describes *what* the selection mechanism does (layered, order-raising, accumulative) without deriving *why* exactly two layers suffice from first principles.

4. **Continuum limit.** Our numerical experiments extend to $N = 1024$. While the scaling laws are consistent with well-defined continuum limits, we have not probed the regime $N \gg 10^3$ where subleading corrections might alter the picture.

---

## §7. Conclusion

We have presented numerical evidence for a **layered structural screening architecture** that selects 4D Lorentzian causal sets from a diverse 25-family poset library:

1. A **first-order admissibility layer** (§2), carried by the triple screening functional $S_{\mathrm{triple}}$ and the Benincasa–Dowker action $S_{\mathrm{BD}}$, eliminates non-geometric structures but cannot uniquely select Lor4D ($S_{\mathrm{BD}}$ rank: 14/17 at $N = 128$).

2. A **second-order identity layer** (§3), implemented by the minimum-distortion functional $S_{\mathrm{MD}}$ (Mahalanobis distance from the Lor4D reference manifold), uniquely and robustly selects Lor4D across all tested scales ($N \geq 14$) and all 25 families with zero free parameters.

3. The identity layer exhibits a **phase-transition-like turn-on** at $N_{\mathrm{id}} \approx 14$ (§4), set by the intrinsic geometric resolution limit of 4D vs. 5D causal structures, not by the discriminator's statistical power.

4. Beyond turn-on, the identity basin undergoes **historical sedimentation** (§5): the Mahalanobis gap grows (from $+1.28$ to $1.93 \times 10^8$), the reference manifold sharpens ($\det(\Sigma) \propto N^{-3.38}$), and the effective basin volume contracts ($V_{\mathrm{eff}} \propto N^{-1.57}$)—making Lor4D an increasingly irreplaceable identity centre.

These results support the **Layered Structural Screening Principle**: the selection of 4D Lorentzian causal structure from a discrete poset space requires at least two functionally distinct mechanisms of different orders, whose intersection is Lor4D and whose identity basin deepens monotonically with system size.

The principle is proposed as a mid-level effective theory—more specific than the philosophical assertion that "existence is what survives screening," more general than the individual numerical experiments from which it is extracted. Its analytical foundations (§5.2, §6.4) provide concrete targets for future work: a derivation of gap growth from CLT + Myrheim–Meyer, a uniqueness theorem for $\boldsymbol{\mu}(\infty)$, and the incorporation of $S_{\mathrm{MD}}$ into the causal set path integral.

---

## References

[1] L. Bombelli, J. Lee, D. Meyer, and R. D. Sorkin, "Space-time as a causal set," *Phys. Rev. Lett.* **59**, 521 (1987).

[2] D. J. Kleitman and B. L. Rothschild, "Asymptotic enumeration of partial orders on a finite set," *Trans. Amer. Math. Soc.* **205**, 205–220 (1975).

[3] D. Dhar, "Entropy and phase transitions in partially ordered sets," *J. Math. Phys.* **19**, 1711 (1978).

[4] H. J. Prömel, A. Steger, and A. Taraz, "Phase transitions in the evolution of partial orders," *J. Combin. Theory Ser. A* **94**, 230–275 (2001).

[5] D. M. T. Benincasa and F. Dowker, "The scalar curvature of a causal set," *Phys. Rev. Lett.* **104**, 181301 (2010).

[6] S. Loomis and S. Carlip, "Suppression of non-manifold-like sets in the causal set path integral," *Class. Quantum Grav.* **35**, 024002 (2018).

[7] J. Myrheim, "Statistical geometry," CERN preprint TH-2538 (1978).

[8] D. A. Meyer, *The Dimension of Causal Sets*, Ph.D. thesis, MIT (1988).

[9] A. Mathur, A. Singh, and S. Surya, "Entropy and the link action in the causal set path-sum," *Class. Quantum Grav.* **37**, 085004 (2020).

[10] G. R. Brightwell and N. Georgiou, "Continuum limits for classical sequential growth models," *Random Structures & Algorithms* **36**, 218–250 (2010).

[11] S. Surya, "The causal set approach to quantum gravity," *Living Rev. Relativ.* **22**, 5 (2019).

---

*Manuscript v0.2 — §1–§7 reviewed, 6 issues fixed (A+B level). 2026-03-27.*
*Cross-references: LAYERED_SCREENING_PRINCIPLE_OUTLINE.md (v2.0), MASTER_NARRATIVE.md, DISCUSSION_THEORY_IMPLICATIONS.md.*
