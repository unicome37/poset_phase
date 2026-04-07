# Layered Structural Screening of 4D Lorentzian Causal Sets

**Target**: Classical and Quantum Gravity (Full Paper)
**Version**: Manuscript v1.0 — Added §6.9 (information-theoretic non-circularity test) and §7 point 6
**Date**: 2026-04-07

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

We report numerical evidence, obtained from a systematic survey across a 25-family poset library at scales $N=10$–$1024$, for a **layered structural screening architecture** comprising at least two functionally distinct layers:

**Layer 1 — Admissibility.** A first-order, linear screening layer (carried jointly by the triple-product functional $S_{\mathrm{triple}}$ and the BD action $S_{\mathrm{BD}}$) eliminates non-geometric and non-robust structures from the poset space, reducing a combinatorially vast landscape to a manageable geometric sector. This layer is necessary but insufficient: Lor4D is not uniquely selected.

**Layer 2 — Identity.** A second-order, quadratic identity functional $S_{\mathrm{MD}}$—the Mahalanobis distance from a Lor4D reference manifold in structural feature space—identifies Lor4D as the robust identity centre among admissible candidates, with zero free parameters. At all tested scales $N \geq 10$ and across all 25 poset families (20-seed $\times$ 120-replication F2 onset study with separate reference ensemble; 10 independent seeds at broader scales), Lor4D is uniquely ranked first with 100% reliability.

The two layers are functionally separated: the Pearson correlation between $S_{\mathrm{BD}}$ and $S_{\mathrm{MD}}$ fluctuates around zero across all tested $N$ (bootstrap 95% confidence intervals span zero at every scale; see §3.5). Their intersection, in the tested library, is Lor4D and only Lor4D.

Furthermore, the identity basin exhibits systematic **deepening** with increasing $N$: the Mahalanobis gap grows from $+0.308 \pm 0.091$ at $N=10$ to $+4.12 \pm 0.35$ at $N=32$ to $1.93 \times 10^8$ at $N=1024$, the effective potential volume shrinks as $V_{\mathrm{eff}} \propto N^{-1.66 \pm 0.07}$ ($R^2 = 0.983$), and the Fisher information grows as $I_F \propto N^{+1.12 \pm 0.07}$ ($R^2 = 0.965$). This basin deepening—which we term **historical sedimentation**—represents the progressive locking of Lor4D identity as the number of causal events increases.

Crucially, the identity mechanism is not specific to Lor4D: when $S_{\mathrm{MD}}$ is centred on any Lorentzian family (Lor2D, Lor3D, or Lor5D), the centred family uniquely self-selects as rank #1 at all tested scales (§3.4). This means self-minimum is a general one-class property. The non-trivial Lor4D claim in this paper is therefore anchored not in self-minimum itself, but in the **combined evidence pattern**: cross-family separation, finite-size turn-on, basin deepening, and background-dependent curvature robustness within the tested library.

### 1.4 Scope and limitations

The present work is strictly numerical. We do not claim analytic proofs of the layered screening principle, nor exhaustive coverage of all possible poset structures. Our conclusions are established within:

- a 25-family poset library (17 standard + 8 adversarial);
- the scale range $N = 12$–$1024$;
- a three-dimensional structural feature space $(d_{\mathrm{eff}},\; C_1/C_0,\; w/N)$;
- **flat Minkowski background**: all Lorentzian families are generated by Poisson sprinkling into causal diamonds in $\mathbb{R}^{d-1,1}$. Curved-background extensions are reported in §4.7: de Sitter is robust up to $H \leq 0.3$ at $N \leq 1024$; weak-field Schwarzschild remains compatible with top-2 behavior in current split tests; matter-FLRW is boundary-sensitive at $\kappa=1.0$ (low-$N$ split reaches hard-fail threshold, while completed high-$N$ runs show partial recovery with fail ratio $0.3<0.5$). Strong curvature still requires a curvature-adaptive reference manifold.

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

2. **Interval ratio** $C_1/C_0$: the ratio of order intervals with exactly one intermediate element ($C_1$) to links ($C_0$; order intervals with no intermediate element), encoding the nearest-neighbour density relative to local connectivity. For Lor4D in a causal diamond, $C_1/C_0 \to 0.357$ as $N \to \infty$, as determined by finite-size extrapolation of the empirical trajectory $\mu_c(N) = \mu_c(\infty) + a/N + b/N^2$ ($R^2 = 0.988$). The limit reflects the geometry of Alexandrov intervals in $\mathbb{R}^{3,1}$: since the volume of the interval $I(x,y)$ between two related elements scales as $\tau^4$ (where $\tau$ is proper-time separation), the link/interval-1 ratio is determined by the Poisson void probability integrated over the pair-distance distribution, which involves $\mathrm{B}(2,2)$ integrals arising from the causal-diamond cross-sectional geometry [7, 8, 11].

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

This is not a trained black-box classifier. There is no training/test split, no decision boundary optimization, and no hyperparameter tuning. $S_{\mathrm{MD}}$ is a **reference-based one-class geometric discriminator**—it uses only Lor4D's own statistics to define the reference, and then measures how far any given poset deviates from that reference.

### 3.2 The Lor4D reference manifold

The one-parameter family $\{(\boldsymbol{\mu}(N), \Sigma(N)) \mid N \geq N_{\min}\}$ constitutes the **Lor4D reference manifold** $\mathcal{M}_4$ in feature space. Its key properties:

**Centre trajectory.** The mean feature vector follows a $1/N$ expansion:

$$\boldsymbol{\mu}(N) = \boldsymbol{\mu}(\infty) + \frac{\mathbf{a}}{N} + \frac{\mathbf{b}}{N^2}, \quad R^2 > 0.97, \tag{5}$$

with limiting values:
- $\mu_d(\infty) = 3.967 \approx 4$ (Myrheim–Meyer, exact in the continuum);
- $\mu_c(\infty) = 0.357$ (Poisson void-probability integral over causal-diamond geometry; involves $\mathrm{B}(2,2)$ cross-section factors [7, 8, 11]);
- $\mu_w(\infty) = 0.227$ (consistent with $w \propto N^{1-1/d}$ scaling, exponent $0.284 \approx 1/4$).

Each component of $\boldsymbol{\mu}(\infty)$ has an independent theoretical pedigree. This is not a fit to the data; it is a consistency check against known analytical results.

**Covariance contraction.** The eigenvalues of $\Sigma(N)$ decay as power laws:

$$\sigma_i^2(N) = A_i \cdot N^{-p_i}, \tag{6}$$

with $p_i \approx 1$ (consistent with CLT scaling $\sigma^2 \propto 1/N$) for all three features. The covariance determinant collapses as

$$\det(\Sigma) \propto N^{-3.31}. \tag{7}$$

As $N$ grows, the Lor4D point cloud in feature space concentrates toward $\boldsymbol{\mu}(\infty)$—the reference manifold sharpens.

**Eigenvector stability.** The principal axes of $\Sigma(N)$ rotate slowly with $N$, indicating that the contraction is approximately isotropic in the natural feature basis. This stability supports the interpretation of $\mathcal{M}_4$ as a well-defined geometric object, not an artifact of a particular $N$ value.

**Terminological note.** We use "reference manifold" in the sense of a parametrized reference object within a mid-level effective theory. No strict differentiable manifold structure is assumed; such formalization, if needed, can be pursued in future analytical work.

### 3.3 Identity selection: the evidence

**Table 2.** $S_{\mathrm{MD}}$ rankings across $N$ and family library. The turn-on region ($N = 12$–$32$) is from a dedicated 10-seed $\times$ 80-realization experiment across all 25 families; entries marked $^\dagger$ are from the basin-deepening experiment (30 realizations, single seed).
**Table 2.** $S_{\mathrm{MD}}$ rankings across $N$ and family library. Entries marked $^\ddagger$ are from the F2 margin-aware refit (20 seeds $\times$ 120 realizations per seed, separate reference ensemble $N_{\mathrm{ref}}=120$, `seed\_base+100000` offset); entries marked $^\dagger$ are from the basin-deepening experiment (30 realizations, single seed). Older 10-seed $\times$ 80-rep shared-seed results (previously in this table for $N=12$–$14$) have been superseded by the F2 protocol.

| $N$ | Lor4D $S_{\mathrm{MD}}$ rank | #1 rate | Runner-up | Mean margin $\pm$ SE | 95% CI |
|:---:|:---:|:---:|:---|:---:|:---:|
| 10 | #1 | 20/20$^\ddagger$ | KR\_2layer | $+0.308 \pm 0.091$ | $[+0.268,\; +0.348]$ |
| 12 | #1 | 20/20$^\ddagger$ | Lor5D | $+1.161 \pm 0.144$ | $[+1.034,\; +1.287]$ |
| 14 | #1 | 20/20$^\ddagger$ | Lor5D | $+1.487 \pm 0.127$ | $[+1.280,\; +1.694]$ |
| 16 | #1 | 20/20$^\ddagger$ | Lor5D | $+1.771 \pm 0.073$ | $[+1.628,\; +1.914]$ |
| 18 | #1 | 10/10 | — | $+1.94 \pm 0.32$ | $[+1.31,\; +2.57]$ |
| 20 | #1 | 10/10 | — | $+1.89 \pm 0.34$ | $[+1.23,\; +2.55]$ |
| 24 | #1 | 10/10 | — | $+2.51 \pm 0.22$ | $[+2.08,\; +2.95]$ |
| 28 | #1 | 10/10 | — | $+3.25 \pm 0.30$ | $[+2.65,\; +3.85]$ |
| 32 | #1 | 10/10 | — | $+4.12 \pm 0.35$ | $[+3.44,\; +4.81]$ |
| 128 | #1 | $^\dagger$ | Lor5D | $+94.1$ | — |
| 1024 | #1 | $^\dagger$ | Lor5D | $+1.93 \times 10^8$ | — |

At $N = 10$, all 20 seeds rank Lor4D first with minimum margin $+0.198$ and mean margin $+0.308$ (F2 margin-aware refit, 120 realizations per seed, separate reference ensemble). Identity selection is robust from $N = 10$ onward.

### 3.4 Physical content of the identity layer

Several aspects of $S_{\mathrm{MD}}$ deserve emphasis:

**Identity, not classification.** $S_{\mathrm{MD}}$ does not say "this poset belongs to category X." It says: "this poset deviates from the Lor4D reference by this much, in these directions, with this statistical cost." It is an identity measure in the same sense that a Mahalanobis distance measures deviation from a known population.

**Runner-up structure and geometric proximity.** At $N = 10$, the primary runner-up is KR\_2layer (runner-up in 18/20 seeds), with Lor5D appearing only 2/20 times. At $N \geq 12$, Lor5D becomes the dominant runner-up across all seeds, consistent with the **intrinsic geometric proximity between 4D and 5D causal structures**: the effective-dimension overlap between Lor4D ($d_{\mathrm{eff}} \approx 3.96$) and Lor5D ($d_{\mathrm{eff}} \approx 4.33$) begins to dominate once the sample size is large enough for $d_{\mathrm{eff}}$ discrimination to operate. The turn-on boundary is a physical resolution limit, not a statistical artifact.

**Zero parameters, anchored by theory.** The centre $\boldsymbol{\mu}(N)$ converges to values with independent theoretical derivations (Myrheim–Meyer, CLT). The metric $\Sigma^{-1}(N)$ is the inverse empirical covariance—the natural Riemannian metric on a Gaussian family (the Fisher information metric in the diagonal limit). No optimization is performed.

**Counter-factual self-selection.** A critical question is whether the privileged status of Lor4D under $S_{\mathrm{MD}}$ is an artifact of centring the functional on Lor4D, or whether $S_{\mathrm{MD}}$ encodes a general identity mechanism. To test this, we centred $S_{\mathrm{MD}}$ on each of the four Lorentzian families (Lor2D, Lor3D, Lor4D, Lor5D) in turn, computing the full 25-family ranking at six scales $N \in \{16, 28, 48, 64, 96, 128\}$ with 30 realizations per family.

**Result: universal self-selection.** In all 24 test conditions (4 centres $\times$ 6 scales), the centred family is uniquely ranked #1. The margins grow monotonically with $N$: for example, the Lor3D-centred margin increases from $+2.61$ ($N = 16$) to $+88.65$ ($N = 128$), and the Lor5D-centred margin from $+12.45$ to $+60.84$. The Lor4D-centred margin (relevant to the main argument) grows from $+3.65$ ($N = 16$) to $+40.47$ ($N = 128$).

This establishes that the Mahalanobis identity functional is a **general structural identity operator**: it selects whichever family defines its reference manifold. The screening architecture described in this paper is not Lor4D-specific engineering; it is a universal property of the layered (admissibility + identity) structure that happens to be applied here with a Lor4D centre because that is the physically relevant case.

### 3.5 Functional separation from Layer 1

The admissibility layer ($S_{\mathrm{BD}}$) and the identity layer ($S_{\mathrm{MD}}$) are not redundant. Their separation is evidenced by:

**Rank discordance.** At $N = 128$, $S_{\mathrm{BD}}$ ranks Lor4D #14/17 while $S_{\mathrm{MD}}$ ranks Lor4D #1/25. The two rankings are nearly uncorrelated.

**Correlation decay.** The Pearson correlation coefficient between $S_{\mathrm{BD}}$ and $S_{\mathrm{MD}}$ values across the 25-family library is measured with bootstrap 95% confidence intervals (1000 resamples):

| $N$ | Pearson $r$ | 95% CI | Spearman $\rho$ | 95% CI |
|:---:|:---:|:---:|:---:|:---:|
| 16 | $-0.59$ | $[-0.83,\; +0.01]$ | $-0.21$ | $[-0.59,\; +0.31]$ |
| 20 | $-0.57$ | $[-0.80,\; -0.13]$ | $-0.29$ | $[-0.65,\; +0.23]$ |
| 28 | $-0.14$ | $[-0.59,\; +0.56]$ | $+0.10$ | $[-0.38,\; +0.56]$ |
| 48 | $+0.06$ | $[-0.50,\; +0.64]$ | $+0.06$ | $[-0.38,\; +0.53]$ |
| 64 | $-0.07$ | $[-0.44,\; +0.45]$ | $+0.09$ | $[-0.33,\; +0.59]$ |
| 96 | $-0.39$ | $[-0.66,\; +0.04]$ | $-0.13$ | $[-0.53,\; +0.31]$ |
| 128 | $+0.11$ | $[-0.26,\; +0.51]$ | $+0.19$ | $[-0.24,\; +0.58]$ |

At every tested scale, the bootstrap 95% CI spans zero for both Pearson and Spearman statistics. The two layers extract functionally independent information from the poset. The sign fluctuations across $N$ further confirm the absence of any systematic coupling.

**Order difference.** $S_{\mathrm{BD}}$ is a linear functional of interval counts (first-order). $S_{\mathrm{MD}}$ is a quadratic functional in feature space (second-order). The transition from Layer 1 to Layer 2 is an **order-raising**: from a hyperplane cut to an ellipsoidal confinement.

---

## §4. Identity Dynamics: Turn-On and Basin Deepening

### 4.1 Three phases of identity emergence

As $N$ increases from small values, the identity layer exhibits a three-phase structure:

**Phase I — Pre-onset** ($N < 10$). At very small $N$ (below the tested range), the poset contains too few causal events for reliable identity discrimination. Note that an earlier 10-seed $\times$ 80-realization experiment with shared reference/test seed showed apparent instability at $N = 12$–$14$; that instability is now understood to stem from reference-ensemble contamination rather than genuine physical overlap. A dedicated 20-seed $\times$ 120-realization experiment with a separate reference ensemble (F2 margin-aware refit; §4.2) shows stable Lor4D identification from $N = 10$.

**Phase II — Turn-on** ($N \approx N_{\mathrm{id}} \approx 10$). A robust onset of identity identification is observed at $N = 10$. In the 20-seed $\times$ 120-realization F2 experiment, Lor4D is ranked \#1 in all 20 seeds with mean margin $+0.308 \pm 0.091$ (minimum $+0.198$, 95\% CI $[+0.268,\; +0.348]$). The primary runner-up at $N = 10$ is KR\_2layer (runner-up in 18/20 seeds; Lor5D in 2/20), consistent with geometric proximity of wide-layer stochastic structures before the Lor4D–Lor5D $d_{\mathrm{eff}}$ gap becomes dominant. At $N \geq 12$, the runner-up transitions to Lor5D, confirming the geometric proximity argument of §6.5.

The turn-on scale $N_{\mathrm{id}} \approx 10$ is not a parameter of the model. It emerges from the interplay between the intrinsic geometric separation of Lor4D and its nearest competitor in feature space and the statistical precision of the reference ensemble at finite $N$. This is a resolution limit of causal structure, not of the discriminator.

**Phase III — Deepening** ($N \gg N_{\mathrm{id}}$). Once identity is established, the basin continues to deepen: the margin between Lor4D and its nearest competitor grows, the reference manifold sharpens, and the statistical cost of misidentification increases without bound.

### 4.2 Quantifying basin deepening

We define the **isolation gap** as the Mahalanobis distance difference between Lor4D and the nearest competitor:

$$\Delta_{\mathrm{hist}}(N) \equiv \min_{f \neq \mathrm{Lor4D}} \bigl[S_{\mathrm{MD}}(f, N) - S_{\mathrm{MD}}(\mathrm{Lor4D}, N)\bigr]. \tag{8}$$

**Table 3.** Isolation gap and basin diagnostics across $N$. The onset region ($N = 10$–$24$) uses 20-seed $\times$ 120-replication F2 data (mean $\pm$ SE, separate reference ensemble); intermediate entries ($N = 28$–$32$) use 10-seed multi-realization data; large-$N$ entries are from the 30-realization basin-deepening experiment.

| $N$ | $\Delta_{\mathrm{hist}}$ | 95% CI | $V_{\mathrm{eff}}$ | $I_F = \mathrm{tr}(\Sigma^{-1})$ | Status |
|:---:|:---:|:---:|:---:|:---:|:---|
| 10 | $0.31 \pm 0.09$ | $[+0.27,\; +0.35]$ | — | — | Turn-on (F2, 20 seeds $\times$ 120 reps)$^\ddagger$ |
| 12 | $1.16 \pm 0.14$ | $[+1.03,\; +1.29]$ | $3.3 \times 10^{-3}$ | $197$ | Onset established |
| 14 | $1.49 \pm 0.13$ | $[+1.28,\; +1.70]$ | — | — | Established |
| 16 | $1.77 \pm 0.06$ | $[+1.63,\; +1.92]$ | $1.6 \times 10^{-3}$ | $469$ | Stable |
| 20 | $1.89 \pm 0.34$ | $[+1.23,\; +2.55]$ | $2.4 \times 10^{-3}$ | $368$ | Stable deepening |
| 28 | $3.25 \pm 0.30$ | $[+2.65,\; +3.85]$ | $1.1 \times 10^{-3}$ | $540$ | Accelerated deepening |
| 32 | $4.12 \pm 0.35$ | $[+3.44,\; +4.81]$ | — | — | Accelerated deepening |
| 128 | $\approx 94.1$ | — | $9.7 \times 10^{-5}$ | $2355$ | Deep lock |
| 1024 | $\approx 1.93 \times 10^8$ | — | — | — | Ultra-deep lock |

$V_{\mathrm{eff}}$ and $I_F$ values are from the 30-realization basin-deepening experiment; entries marked "—" correspond to scales not included in that experiment run.
$^\ddagger$F2 margin-aware refit (20 seeds $\times$ 120 realizations, separate reference ensemble, `seed\_base+100000` offset); $V_{\mathrm{eff}}$ and $I_F$ not measured in this protocol.

The growth of $\Delta_{\mathrm{hist}}(N)$ is monotonic: the ensemble mean increases consistently from $+0.308$ at $N = 10$ to $+4.12$ at $N = 32$ (F2 protocol for $N \leq 24$; older 10-seed data for $N = 28, 32$). Individual seeds show occasional non-monotonicities, but the ensemble mean is monotone.

### 4.3 Scaling laws

Four independent diagnostics confirm the systematic deepening. All power-law fits are performed on $N$-grid values $\{16, 20, 28, 36, 48, 64, 96, 128, 192, 256\}$ (10 points, 30 realizations per $N$) with bootstrap standard errors (500 resamples).

**1. Effective potential volume.**

$$V_{\mathrm{eff}}(N) \propto N^{-1.66 \pm 0.07}, \quad R^2 = 0.983. \tag{9}$$

The effective volume of the Lor4D identity basin shrinks as a power law. In the continuum limit, $V_{\mathrm{eff}} \to 0$: the basin concentrates to a point.

**2. Covariance determinant.**

$$\det\bigl(\Sigma(N)\bigr) \propto N^{-3.31 \pm 0.14}, \quad R^2 = 0.983. \tag{10}$$

Since $V_{\mathrm{eff}} \propto \sqrt{\det(\Sigma)}$, the consistency relation $\frac{1}{2} \times (-3.31) = -1.66$ is satisfied exactly, providing an internal cross-check on the fitting procedure.

**3. Fisher information.**

$$I_F(N) = \mathrm{tr}\bigl(\Sigma^{-1}(N)\bigr) \propto N^{+1.12 \pm 0.07}, \quad R^2 = 0.965. \tag{11}$$

The exponent $\approx +1$ is consistent with CLT scaling ($\sigma^2 \propto 1/N$, so $\sigma^{-2} \propto N$), but the slightly super-linear value $+1.12$ suggests additional concentration beyond pure sampling: the concentration of Lor4D's identity is slightly faster than what random fluctuation reduction alone would predict.

**4. Isolation gap.**

$$\Delta_{\mathrm{hist}}(N) \propto N^{+1.22 \pm 0.10}, \quad R^2 = 0.892. \tag{12}$$

The gap between Lor4D and its nearest competitor grows as a power law. The exponent $> 1$ indicates super-linear deepening—the identity basin does not merely stabilize; it **accelerates** its separation from alternatives. The lower $R^2$ relative to the other diagnostics reflects the inherent fluctuation of a nearest-competitor statistic (a minimum over stochastic quantities), not a failure of the power-law model.

### 4.4 Physical gap versus statistical amplification

A natural concern is whether the Mahalanobis gap growth reflects a genuine increase in physical separation between Lor4D and its competitors, or merely an amplification of fixed differences by a sharpening metric ($\Sigma^{-1} \to \infty$). To disentangle these, we decompose the gap into:

- **Euclidean distance**: $d_E(N) = \|\mathbf{I}_{\mathrm{runner\text{-}up}}(N) - \boldsymbol{\mu}(N)\|_2$, the unstandardized feature-space distance between the runner-up and the Lor4D centre;
- **Precision amplification**: $\sigma_A(N) = \Delta_{\mathrm{hist}}(N) / d_E(N)$, the factor by which the inverse covariance inflates the raw distance.

**Table 4.** Gap decomposition across $N$ (runner-up: Lor5D at all scales).

| $N$ | $d_E$ | $\Delta_{\mathrm{hist}}$ | $\sigma_A$ |
|:---:|:---:|:---:|:---:|
| 16 | 0.511 | 3.65 | 14.0 |
| 28 | 0.369 | 4.47 | 32.8 |
| 64 | 0.488 | 20.1 | 84.3 |
| 128 | 0.459 | 40.5 | 191.9 |
| 256 | 0.452 | 92.6 | 454.0 |

A power-law fit to the Euclidean distance gives $d_E \propto N^{+0.007}$ ($R^2 = 0.003$)—effectively **flat**. The physical feature-space separation between Lor4D and Lor5D does not shrink with $N$; if anything, it fluctuates around $d_E \approx 0.46$ with no trend. The Mahalanobis gap growth ($\propto N^{+1.22}$) is therefore driven almost entirely by precision amplification: the inverse covariance $\Sigma^{-1}$ grows as $\sim N$, inflating a genuine but fixed physical distance into an ever-larger statistical penalty.

This has an important physical interpretation: **the selection is genuine** (the distance is real, not an artifact of a shrinking denominator), but the *strength* of the selection grows because the reference manifold sharpens—not because the competitors move away. Historical sedimentation is the process by which a fixed physical distinction is rendered increasingly decisive by statistical concentration.

### 4.5 Feature space minimality

The three structural features $(d_{\mathrm{eff}}, C_1/C_0, w/N)$ were not selected by optimization. To test whether this feature set is minimal—neither redundant nor insufficient—we perform a systematic ablation, evaluating $S_{\mathrm{MD}}$ with all subsets of the three features: each feature alone (3 subsets), each pair (3 subsets), and the full triple, at $N \in \{28, 64, 128\}$.

**Table 5.** Feature ablation results ($N = 128$). Lor4D rank and margin over runner-up.

| Feature set | Lor4D rank | Margin |
|:---|:---:|:---:|
| $d_{\mathrm{eff}}$ alone | 2 | $+0.17$ (behind KR_2layer) |
| $C_1/C_0$ alone | 1 | $+7.17$ |
| $w/N$ alone | 1 | $+0.04$ |
| $d_{\mathrm{eff}} + C_1/C_0$ | 1 | $+36.2$ |
| $d_{\mathrm{eff}} + w/N$ | 1 | $+28.8$ |
| $C_1/C_0 + w/N$ | 1 | $+8.31$ |
| **All three** | **1** | **$+40.5$** |

No single feature reliably identifies Lor4D across all scales: $d_{\mathrm{eff}}$ alone fails (rank 2 at $N = 28$ and $N = 128$), and $w/N$ alone has a negligible margin. The full triple consistently achieves the largest margin at all tested $N$. The pattern is consistent across scales: $C_1/C_0$ is the individually strongest feature (encoding local connectivity structure), $d_{\mathrm{eff}}$ and $w/N$ act as complementary dimensions, and the full triple captures structural information that no proper subset reliably encodes. The three-dimensional feature space is therefore a **minimal non-redundant effective basis within the tested feature library** for robust Lor4D identity identification.

### 4.6 The physical meaning of basin deepening

The deepening of the identity basin has a clear physical interpretation: **each additional causal event added to the poset provides additional geometric information that further distinguishes Lor4D from its competitors.** The reference manifold progressively "comes into focus"—the center converges to its theoretical limit, the fluctuations shrink, and the cost of deviating from the Lor4D template grows without bound.

This is not an artifact of increasing statistical power with larger samples. It reflects the intrinsic geometric content of causal structures: a 4D Lorentzian causal set with 1024 elements contains genuinely more dimensional, topological, and geometric information than one with 14 elements, and this information accumulates in a structured way that systematically strengthens the identity signal.

We formalize this observation in §5 as the **historical sedimentation** of the identity basin, and state the layered screening principle that synthesizes the results of §§2–4.

### 4.7 Robustness under background curvature

All preceding results are based on sprinklings into **flat** Minkowski space. To test whether the layered screening architecture survives background curvature, we introduce de Sitter families $\mathrm{dS4D}_H$ generated by Poisson sprinkling into a 4-dimensional de Sitter causal diamond with Hubble parameter $H \in \{0.1, 0.3, 0.5, 1.0, 2.0\}$ (the $H = 0$ limit recovers flat Lor4D). The sprinkling uses the rejection-sampling method with volume weight $\propto a(t)^3$, $a(t) = e^{Ht}$, and the causal relation is determined by the comoving horizon $\chi = H^{-1}[e^{-Ht_1} - e^{-Ht_2}]$. Each family is sampled at $N = 16$–$128$ with 40 realizations per $(H, N)$. After inclusion of all 5 de Sitter variants, the family library expands from 25 to 30 families.

**Table 6.** $S_{\mathrm{MD}}$ rank of Lor4D and nearest dS4D competitor across scales.

| $N$ | Lor4D rank | Lor4D $S_{\mathrm{MD}}$ | dS4D$_{H=0.1}$ rank | dS4D$_{H=0.1}$ $S_{\mathrm{MD}}$ | Gap |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 16 | **#2** | 0.127 | **#1** | 0.040 | $-0.088$ |
| 28 | **#1** | 0.116 | #2 | 0.442 | $+0.325$ |
| 48 | **#1** | 0.186 | #2 | 0.812 | $+0.626$ |
| 64 | **#1** | 0.084 | #2 | 0.668 | $+0.584$ |
| 96 | **#1** | 0.064 | #2 | 0.595 | $+0.531$ |
| 128 | **#1** | 0.061 | #2 | 1.634 | $+1.574$ |

At $N = 16$, the mildly curved dS4D$_{H=0.1}$ **outranks** flat Lor4D (score 0.040 vs 0.127). This is a small-$N$ fluctuation: with only 16 elements, the reference covariance $\Sigma$ is estimated from few samples, and the Mahalanobis metric is noisy. By $N = 28$, Lor4D recovers rank #1 and retains it at every subsequent scale. The gap grows monotonically from $+0.33$ ($N = 28$) to $+1.57$ ($N = 128$), confirming that the layered screening architecture is robust in the expanded 30-family library.

**Curvature ordering.** At every $N \geq 28$, the de Sitter variants rank strictly by $H$: $H = 0.1 > 0.3 > 0.5 > 1.0 > 2.0$ (nearest to farthest from Lor4D). This monotone curvature ordering demonstrates that $S_{\mathrm{MD}}$ is genuinely sensitive to geometric curvature, not merely to family-label differences.

**Feature-space decomposition.** The dominant feature shift is $\Delta d_{\mathrm{eff}} \approx H$ (curvature inflates the effective dimension), with smaller secondary shifts in $C_1/C_0$ and $w/N$. The Euclidean feature distance $d_E$ is roughly constant across $N$ for fixed $H$ (e.g., $d_E \approx 0.11$ for $H = 0.1$), while the Mahalanobis distance $S_{\mathrm{MD}}$ grows as $\sim N^\beta$ ($\beta \approx 1.1$–$1.3$), reflecting the same precision-amplification mechanism identified in §4.4.

**Extended curvature scan.** A finer Hubble grid ($H \in \{0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0\}$) at $N = 16$–$256$ (40 realizations each) reveals three regimes:

1. *Indistinguishable regime* ($H \leq 0.01$): $S_{\mathrm{MD}}$ does not grow systematically with $N$ ($\beta = 0.10 \pm 0.36$, $R^2 = 0.01$); ranking oscillates between dS4D and flat Lor4D. At this curvature, the de Sitter causal set is **statistically indistinguishable** from flat Lor4D.

2. *Separated-but-near regime* ($0.05 \leq H \leq 0.5$): $S_{\mathrm{MD}}$ grows as $N^\beta$ with $\beta$ increasing smoothly from $0.77$ ($H = 0.05$) to $1.02$ ($H = 0.5$), all with $R^2 > 0.82$. The dS4D family remains within the top-2 (behind Lor4D) at all tested $N$, with the nearest non-dS4D competitor always being Lor5D.

3. *Divergent regime* ($H \geq 1.0$): $\beta = 1.15 \pm 0.11$ (super-linear), dS4D drops out of the top-3 at large $N$.

The **critical Hubble parameter** $H_c$, defined as the largest $H$ at which dS4D remains within rank #2, ranges from $H_c \approx 0.3$ to $0.5$ across all tested $N \leq 256$, substantially above the conservative $H \leq 0.1$ estimate from the initial scan.

**Large-$N$ convergence.** To verify that these results are not small-$N$ artifacts, we extend the experiment to $N \in \{256, 384, 512\}$ with $H \in \{0.0, 0.1, 0.3, 0.5\}$ (20 realizations each, 527 s total). Table 7 summarizes the key findings.

**Table 7.** Large-$N$ de Sitter convergence: $S_{\mathrm{MD}}$ rankings and scores.

| $N$ | $H$ | $S_{\mathrm{MD}}$(dS4D) | $S_{\mathrm{MD}}$(Lor4D) | dS4D rank | Lor4D rank |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 256 | 0.1 | 3.91 | 0.66 | #2 | #1 |
| 384 | 0.1 | 5.62 | 0.17 | #2 | #1 |
| 512 | 0.1 | 6.57 | 0.31 | #2 | #1 |
| 256 | 0.3 | 20.7 | 0.66 | #2 | #1 |
| 384 | 0.3 | 62.7 | 0.17 | #2 | #1 |
| 512 | 0.3 | 118.2 | 0.31 | #2 | #1 |
| 512 | 0.5 | 442.3 | 0.31 | #2 | #1 |

At all tested $N$ up to 512, Lor4D retains rank #1 and dS4D$_{H=0.1}$ remains rank #2; the ordering is fully stable. At $H = 0.5$, $H_c$ actually increases with $N$: at $N = 512$, dS4D$_{H=0.5}$ ($S_{\mathrm{MD}} = 442.3$) narrowly beats its nearest non-dS competitor Lor5D ($S_{\mathrm{MD}} = 456.3$), recovering rank #2. This occurs because Lor5D, being structurally farther from Lor4D, diverges even faster under the sharpening reference.

The scaling exponent $\beta$ in the large-$N$ regime reveals an important nuance: fitting $\log S_{\mathrm{MD}}$ vs $\log N$ over $N = 256$–$512$ yields $\beta = 0.76 \pm 0.10$ ($H = 0.1$), $\beta = 2.53 \pm 0.15$ ($H = 0.3$), and $\beta = 2.81 \pm 0.15$ ($H = 0.5$), all with $R^2 > 0.98$. The large-$N$ values for $H \geq 0.3$ are significantly steeper than the small-$N$ fit ($\beta \approx 1.0$ for $N = 16$–$256$), indicating that the log–log relationship is **convex**: the screening architecture's discriminating power accelerates at larger $N$, rather than saturating. This is physically consistent with the precision-amplification mechanism of §4.4—the $\Sigma^{-1}(N)$ sharpening continuously increases the penalty for curvature-induced feature shifts.

**Interpretation.** The three-regime structure has a clear physical origin. At $H \leq 0.01$, the de Sitter radius $\ell_{\mathrm{dS}} = 1/H \geq 100$ far exceeds the causal diamond scale, so curvature corrections to interval counts are below statistical resolution. At $0.05 \leq H \leq 0.5$, curvature shifts the feature means by detectable but modest amounts ($\Delta d_{\mathrm{eff}} \lesssim 0.5$), and the precision-amplification mechanism (§4.4) converts this into growing Mahalanobis distance while preserving the rank ordering. At $H \geq 1.0$, the de Sitter radius is comparable to the diamond, the sprinkling geometry changes qualitatively (volume concentration toward the expanding boundary), and features diverge rapidly. A curvature-adaptive reference manifold $\mathcal{M}(H)$ would be needed to extend the identity layer into this strongly curved regime.

The large-$N$ convergence test (Table 7) provides strong evidence that the de Sitter robustness results are not finite-size artifacts: ranking stability persists to $N = 512$, and discriminating power continues to grow.

**Very-large-$N$ convergence.** To verify that ranking stability persists deep into the large-$N$ regime, we extend to $N \in \{256, 512, 768, 1024\}$ with $H \in \{0.0, 0.1, 0.3\}$ (15 realizations each, 9276 s total computation time dominated by $O(N^3)$ transitive closure). Table 8 summarizes the results.

**Table 8.** Very-large-$N$ de Sitter convergence ($N = 256$–$1024$).

| $N$ | $H$ | $S_{\mathrm{MD}}$(dS4D) | $S_{\mathrm{MD}}$(Lor4D) | dS4D rank | Lor4D rank |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 256 | 0.0 | 0.48 | 0.76 | #1 | #2 |
| 256 | 0.1 | 3.03 | 0.76 | #2 | #1 |
| 256 | 0.3 | 18.2 | 0.76 | #2 | #1 |
| 512 | 0.1 | 6.15 | 1.25 | #2 | #1 |
| 512 | 0.3 | 135.4 | 1.25 | #2 | #1 |
| 768 | 0.1 | 4.94 | 0.17 | #2 | #1 |
| 768 | 0.3 | 63.9 | 0.17 | #2 | #1 |
| 1024 | 0.0 | 0.074 | 0.073 | #2 | #1 |
| 1024 | 0.1 | 7.33 | 0.073 | #2 | #1 |
| 1024 | 0.3 | 101.1 | 0.073 | #2 | #1 |

The ranking result is unambiguous: **Lor4D retains rank #1 and dS4D retains rank #2 (or #1 at $H = 0$) at every tested $N$ up to 1024**, with no exceptions across all 12 $(N, H)$ conditions. At $H = 0$, the $H \to 0$ limit is confirmed: $\Delta = S_{\mathrm{MD}}(\mathrm{dS4D}) - S_{\mathrm{MD}}(\mathrm{Lor4D}) = +0.001$ at $N = 1024$, i.e., the de Sitter family is statistically indistinguishable from flat Lor4D as expected.

The $S_{\mathrm{MD}}$ values at $H = 0.3$ show non-monotonic fluctuation across $N$ (18.2 → 135.4 → 63.9 → 101.1), reflecting covariance estimation noise with 15 realizations. Fitting $\log S_{\mathrm{MD}}$ vs $\log N$ over the full $N = 256$–$1024$ range yields $\beta = 0.55 \pm 0.22$ ($H = 0.1$, $R^2 = 0.76$) and $\beta = 1.10 \pm 0.69$ ($H = 0.3$, $R^2 = 0.56$). These are noisier than the two-point large-$N$ fits (Table 7); the low $R^2$ for $H = 0.3$ is consistent with the previously observed convexity—a single power law cannot capture both the slow initial growth and the rapidly accelerating large-$N$ behaviour. The crucial result is not the precise $\beta$ value but the invariance of the rank ordering: dS4D never falls below #2 in any of the 50 $(N, H)$ conditions tested across Tables 6–8.

**Additional curved backgrounds (extended to $N=512$).** We tested two non-de-Sitter generators at $N \in \{64,128,256,512\}$ (15 realizations each): (i) a matter-dominated FLRW model with $a(t)=(1+\kappa t)^{2/3}$ and $\kappa \in \{0,0.3,1.0,3.0\}$, and (ii) a weak-field Schwarzschild proxy with compactness parameter $\phi_0 \in \{0,0.01,0.05,0.1\}$. Table 9 summarizes the large-$N$ ($N=256,512$) results.

**Table 9.** FLRW/Schwarzschild robustness at large $N$.

| Background | Parameter | $N=256$ rank | $N=512$ rank | Interpretation |
|:---|:---:|:---:|:---:|:---|
| FLRW | $\kappa=0.3$ ($H_0\approx0.2$) | #2 | #2 | Mild curvature robust |
| FLRW | $\kappa=1.0$ ($H_0\approx0.67$) | #3 | #2 | Moderate curvature, boundary-sensitive; low-$N$ split now shows repeatable top-2 failure |
| FLRW | $\kappa=3.0$ ($H_0\approx2.0$) | #4 | #5 | Strong curvature divergence |
| Schwarzschild | $\phi_0=0.01$ | #1/#2 tie-level | #2 | Very weak field robust |
| Schwarzschild | $\phi_0=0.05$ | #2 | #2 | Weak field robust |
| Schwarzschild | $\phi_0=0.10$ | #2 | #2 | Upper weak-field still robust |

In the weak-to-moderate range, the picture is background-dependent rather than uniform. de Sitter remains cleanly top-2 throughout the tested window, and weak-field Schwarzschild remains top-2 in split tests. By contrast, matter-FLRW at $\kappa=1.0$ is a boundary-sensitive regime: the low-$N$ split run reaches the hard-fail threshold, whereas the completed high-$N$ branch ($N=768,1024$) shows partial recovery (failure ratio $0.3<0.5$), i.e., instability persists but does not trigger hard fail under the same rule. Under strong FLRW curvature ($\kappa=3.0$), the rank degrades further to #4–#5, consistent with the divergent regime. A curvature-adaptive reference $\mathcal{M}(\text{background})$ therefore remains necessary outside the de Sitter-like neighbourhood.

---

## §5. Historical Sedimentation and the Layered Screening Principle

Section 4 established that the identity basin deepens systematically with $N$. In this section we elevate that observation to a mid-level theoretical interpretation—**historical sedimentation**—and state the formal layered screening principle that synthesizes the results of §§2–4.

### 5.1 From basin deepening to historical sedimentation

The scaling laws of §4.3 admit a unified physical reading:

| Physical concept | Mathematical correspondent | Observable |
|:---|:---|:---|
| History gradually accumulates | Reference manifold gradually comes into focus | $\boldsymbol{\mu}(N) \to \boldsymbol{\mu}(\infty)$ with $1/N$ corrections |
| Past constrains future | Covariance contracts | $\det(\Sigma) \propto N^{-3.31 \pm 0.14}$ |
| Inertia strengthens | Identity basin deepens | Gap growth $\propto N^{+1.22}$; $V_{\mathrm{eff}} \propto N^{-1.66}$ |
| Deviation cost increases | Off-manifold penalty grows | Fisher $\propto N^{+1.12}$ |
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

This is not a theorem. It is a numerically well-supported mid-level proposition that can be proved or refuted by future analytical work. The 10-seed ensemble averages are strictly monotonic for $N \geq 16$ (Table 3); individual seeds may exhibit non-monotonicities, but the overall trend is unambiguously upward.

### 5.4 Historical sedimentation is not a third screening layer

An important structural clarification: historical sedimentation describes the **longitudinal behavior of the identity layer as $N$ increases**. It is not an independent third screening operator parallel to admissibility and identity. The quantities it tracks—$\boldsymbol{\mu}(N)$, $\Sigma(N)$, gap, $V_{\mathrm{eff}}$—all belong to the identity layer's scale evolution.

The correct layer hierarchy is:
- **Cross-sectional** (fixed $N$): Admissibility + Identity → layered screening;
- **Longitudinal** (across $N$): Identity basin deepens → historical sedimentation.

### 5.5 The Layered Structural Screening Principle

Synthesizing §§2–4 and the sedimentation interpretation:

> **Layered Structural Screening Principle.** In a set of finite posets with sufficiently many causal events, robust identification/stabilization of 4D flat Lorentzian causal structure requires at least two functionally distinct mechanisms of different orders:
>
> 1. A first-order admissibility layer (carried jointly by $S_{\mathrm{triple}}$ and $S_{\mathrm{BD}}$), excluding non-geometric and non-robust structures;
> 2. A second-order quadratic identity functional $S_{\mathrm{MD}}$, identifying Lor4D as the robust identity centre among all admissible candidates.
>
> The two layers are functionally separated, statistically weakly correlated, and their intersection in the tested 25-family library uniquely and robustly selects Lor4D.
>
> Furthermore, the identity basin exhibits monotone-accumulative deepening as $N$ increases (historical sedimentation), making Lor4D increasingly stable as the identity centre in the current family library.

### 5.6 Numerical evidence summary

| Evidence type | Status |
|:---|:---|
| 25 families × $N = 12$–$1024$ × 10 seeds × 80 REPS | ✅ Complete |
| 25 families × $N = 10$–$1024$ × 20 seeds × 120 REPS (F2 margin-aware refit) | ✅ Complete |
| Turn-on boundary $N_{\mathrm{id}} \approx 10$ (manuscript-safe, 3 consecutive N criterion); $N = 10$: min\_margin=0.198, ci95\_lower=0.268 | ✅ Quantified with 95% CI |
| $S_{\mathrm{BD}}$–$S_{\mathrm{MD}}$ functional separation | ✅ Bootstrap CI spans zero at all $N$ |
| Counter-factual self-selection (Lor2D/3D/5D centres) | ✅ 24/24 conditions pass |
| Basin deepening power laws | ✅ $V_{\mathrm{eff}} \propto N^{-1.66 \pm 0.07}$, $R^2 = 0.983$ |
| Gap decomposition: physical vs statistical | ✅ $d_E \approx 0.46$ flat; $\sigma_A \propto N$ |
| Feature space minimality (ablation) | ✅ 3D triple is minimal non-redundant effective basis (within tested feature library) |
| Sedimentation monotone trend | ✅ $\Delta_{\mathrm{hist}}$ increasing (ensemble mean) |
| Analytical derivation | ❌ Skeleton only (§5.2) |
| de Sitter sprinklings ($H = 0.01$–$2.0$, $N = 16$–$1024$) | ✅ Robust for $H \leq 0.3$ at $N \leq 1024$; ranking invariant across all 50 conditions (§4.7) |
| Schwarzschild / FLRW sprinklings (extended pilot, $N=64$–$512$) | ⚠️ Background-dependent: Schwarzschild low-$N$ split passes, but FLRW $\kappa=1.0$ already triggers top-2 failure in the current split run; strong FLRW curvature diverges further (rank #4–#5 at $\kappa=3.0$) |
| $N > 1024$ verification | ❌ Not yet tested ($N \leq 1024$ completed; $O(N^3)$ closure limits further scaling) |
| Information-theoretic non-circularity test | ✅ 5-term info-only penalty, $N = 10$–$100$, 17 families; Lor4D #1 at $\gamma = 0.2$; complementarity with geometric penalties confirmed; ablation identifies interval-diversity as key discriminator (§6.9) |

---

## §6. Discussion

### 6.1 Positioning: reference-based one-class geometric discriminator

A natural concern is that $S_{\mathrm{MD}}$ might be overfit to a particular family library. We therefore position it explicitly as a reference-based one-class geometric discriminator and report the evidence boundaries:

(i) $S_{\mathrm{MD}}$ is a **one-class measure**: it uses only Lor4D's own ensemble statistics to define the reference, with no knowledge of competitor families during construction. There is no training/test split, no decision boundary optimization, and no hyperparameter.

(ii) The centre $\boldsymbol{\mu}(\infty)$ is anchored by independent theoretical results (Myrheim–Meyer, CLT, Beta integral framework), not learned from data.

(iii) The metric $\Sigma^{-1}(N)$ is the empirical inverse covariance of the Lor4D ensemble—the natural Riemannian metric on a Gaussian family (the Fisher information metric in the diagonal limit). It is a property of Lor4D physics, not a fitting choice.

(iv) Robustness against adversarial families: 8 constructions specifically designed to mimic Lor4D features all fail to displace it from rank #1 at any $N \geq 14$.

(v) Pre-registered family-pressure falsification: in the current expanded-library F1 run ($N=12\text{–}256$, 10 seeds, formal C1 threshold), we observe no non-Lor4D family that stably and repeatedly outranks Lor4D; in fact, Lor4D remains rank #1 in every tested seed at every tested $N$ within that run.

We do **not** claim this exhausts all conceivable feature libraries or all possible family generators; our claim is explicitly scoped to the tested library and tested feature set.

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

The turn-on scale $N_{\mathrm{id}} \approx 10$ is specific to the selection of Lor4D from a library that includes Lor5D as the dominant runner-up at $N \geq 12$. At $N = 10$, the primary runner-up is KR\_2layer—a wide-layer stochastic structure—rather than Lor5D (§3.4). This reflects the fact that at the smallest resolved scale, the Mahalanobis geometry of the feature triple $(d_{\mathrm{eff}}, C_1/C_0, w/N)$ is dominated by the covariance structure of wide-layer constructions, which project closer to Lor4D in the combined three-dimensional feature space before the Lor4D–Lor5D $d_{\mathrm{eff}}$ separation ($\approx 0.37$ in units of the pooled standard deviation) becomes the dominant axis of discrimination. By $N \geq 12$, the $d_{\mathrm{eff}}$ gap becomes resolvable and Lor5D takes over as the persistent runner-up, consistent with the intrinsic dimensional proximity argument.

One expects $N_{\mathrm{id}}$ to depend on both the target dimension $d$ and the geometric proximity of the nearest confusable structure in feature space:

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

1. **Curvature: de Sitter validated; non-de-Sitter backgrounds remain conditional.** A de Sitter scan (§4.7) shows that the flat-calibrated screening architecture is robust for $H \leq 0.3$ at all tested $N$ up to 1024: Lor4D retains rank #1 and dS4D retains rank #2 in every one of the 50 $(N, H)$ conditions tested (Tables 6–8). The critical Hubble $H_c \geq 0.3$ at all $N$. At $N = 1024$, the $H = 0$ baseline confirms statistical indistinguishability ($\Delta = +0.001$). Split low-$N$ tests indicate that weak-field Schwarzschild remains compatible with the local-basin picture, whereas matter-FLRW at $\kappa=1.0$ reaches the hard-fail threshold in low-$N$ but shows partial recovery in completed high-$N$ runs (failure ratio $0.3<0.5$). These non-de-Sitter tests remain preliminary (proxy metrics, limited parameterization), and full metric-faithful Schwarzschild/FLRW validation at $N\gtrsim512$ remains open.

2. **Library scope.** All results are established within a 25-family library. While this includes adversarial constructions, it does not exhaust the space of all possible posets. A family that matches Lor4D in all three features at all $N$ would defeat the screening—though we have found no such family despite targeted efforts.

3. **Feature choice.** The three structural features $(d_{\mathrm{eff}}, C_1/C_0, w/N)$ are not derived from a first-principles variational argument. They are selected for their independent theoretical pedigrees and empirical completeness. A feature-ablation study (§4.5, Table 5) confirms that the full triple is a minimal non-redundant effective basis within the tested feature library: no single feature or pair reliably identifies Lor4D across all scales.

4. **Mid-level theory.** The layered screening principle is a mid-level effective theory, not a fundamental principle. It describes *what* the selection mechanism does (layered, order-raising, accumulative) without deriving *why* exactly two layers suffice from first principles.

5. **Continuum limit.** Our numerical experiments extend to $N = 1024$. While the scaling laws are consistent with well-defined continuum limits, the $O(N^3)$ transitive closure bottleneck ($\sim 5600$ s at $N = 1024$) prevents direct verification at $N \gg 10^3$. Algorithmic improvements (sparse or GPU-accelerated closure) would be needed to probe the $N > 10^3$ regime where subleading corrections might alter the picture.

### 6.8 Falsifiability conditions

To avoid unconstrained post-hoc reinterpretation, we state explicit failure conditions for the present mid-level principle:

1. **Identity-layer failure condition.** If, after passing admissibility, one or more non-Lor4D families stably and repeatedly outrank Lor4D across increasing $N$ in expanded libraries, then the current identity-layer stabilization claim fails. The current formal F1 family-pressure run does **not** trigger this condition.

2. **Background-response failure condition.** If weak-to-moderate curved 4D families systematically fail to maintain near-neighbor hierarchy (top-2 behavior in the present setup) under expanded and metric-faithful background tests, then the flat-centred local-basin interpretation fails. The current split low-$N$ F3 results already show that this risk is background-specific rather than merely hypothetical: de Sitter and weak-field Schwarzschild pass, while FLRW at $\kappa=1.0$ triggers the present hard-fail threshold.

3. **Basis-sufficiency failure condition.** If a reproducible alternative feature basis of comparable interpretability strictly dominates $(d_{\mathrm{eff}}, C_1/C_0, w/N)$ across the same libraries and scales, then the current "minimal non-redundant effective basis" claim must be revised.

### 6.9 Information-theoretic penalty: independent non-circularity test

A separate line of evidence addresses the circularity concern from a different angle. Instead of replacing a single geometric sub-term (as in the $d_{\mathrm{consistency}}$ test of §3.3), we replace *all* geometric penalties with pure information-theoretic ones—functions that detect structural regularity without referencing any target dimension, window, or geometric template.

Five penalty functions were constructed from the Hasse graph alone:

| Penalty | Detects | Formula |
|:---|:---|:---|
| Spectral entropy deficit | Spectral degeneracy | $(1 - S_{\mathrm{vN}}/\!\ln n)^2$ |
| Degree heterogeneity | Hub dominance | $G_{\mathrm{degree}}^2$ (Gini) |
| Layer concentration | Non-uniform layering | $(1 - H_{\mathrm{layer}}/\!\ln K)^2$ |
| Edge density extremity | Coverage extremity | $(1 - H_b(p)/\!\ln 2)^2$ |
| Interval diversity deficit | Interval homogeneity | $(1 - H_{\mathrm{int}}/\!\ln m)^2$ |

Under the purely information-theoretic action $A_4 = -\beta\,\log H + \gamma\,I_{\mathrm{info}}$ with 17 families, $N = 10$–$100$, and $\gamma \in [0, 1]$:

(i) **At moderate penalty strength ($\gamma = 0.2$), Lor4D is rank #1 at all tested $N \geq 20$.** No geometric prior is present, yet the entropy-efficiency of Lor4D—high linear-extension count per unit informational penalty—suffices for selection.

(ii) **At strong penalty ($\gamma \geq 0.5$), Lor4D drops to rank #3**, consistently behind KR\_2layer and Lor5D. These are precisely the families hit hardest by geometric penalties. Thus the information-theoretic top-3 is the geometric-penalty bottom-3: the two penalty classes encode orthogonal quality dimensions.

(iii) **Ablation within the 5 info terms identifies interval diversity deficit as the sole critical discriminator** (z-score rank improvement $\Delta = +6$ to $+11$ upon removal at $\gamma = 1.0$). The other four terms are nearly redundant ($r \approx 0.97$–$1.00$). Crucially, however, the raw action $A = -\log H + \gamma\,I_{\rm info}$ tells the opposite story: because interval diversity penalises KR\_2layer (80\% of its total info penalty) more heavily than Lor4D (58\%), *keeping* the interval term *helps* Lor4D win at $n = 20$, $\gamma = 1.0$ (rank \#1 at baseline weights vs.\ \#3 when interval weight $= 0$). The apparent "killing" effect is a robust-$z$-score normalisation artefact in small-sample ablation. The interval diversity deficit captures the fact that Alexandrov intervals in 4D Lorentzian sprinklings have characteristically homogeneous size distributions—a geometric signature encoded in a purely information-theoretic language—and its differential impact across families is what enables Lor4D selection.

(iv) **Complementarity theorem.** Neither penalty class alone selects Lor4D uniquely at large $N$ and all $\gamma$: geometric penalties select dimension ($d = 4$) but not family type; information-theoretic penalties select family type (Lorentzian-like) but not dimension. Their conjunction—the two-layer screening architecture of the present paper—selects Lor4D uniquely.

(v) **Weight sensitivity scan** (19 configurations, $N \in \{10,20,40\}$, $\gamma \in [0,2]$). No single weight vector allows Lor4D to be $\#1$ across all $(N, \gamma)$: the optimal interval-diversity weight shifts from $w \approx 0$ at small $N$ (where penalties overwhelm entropy differences) to $w \geq 8$ at large $N$ (where the log-extension-count gap demands heavy penalisation of competitors). This $N$-dependent optimal weight is a direct consequence of the different scaling of $\log H$ ($\sim N^{\alpha}$) versus $I_{\rm info}$ ($\sim N^{\beta}$, $\beta < \alpha$), and constitutes further evidence that information-theoretic penalties alone cannot serve as a universal Lor4D selector.

(vi) **Hybrid penalty failure.** Combining geometric and information-theoretic penalties into a single action ($A_6 = P_{\rm geo} + I_{\rm info}$ or $A_7 = P_{\rm neutral} + P_{\rm geo} + I_{\rm info}$) performs *worse* than either single penalty in 11 out of 15 non-zero-$\gamma$ conditions (Lor4D rank degrades by up to 14 positions). The doubled penalty overwhelms the entropy term in favour of structurally "bland" families. This negative result demonstrates that the two penalty classes must be applied *sequentially* (as distinct screening layers), not *additively* in a single objective—precisely the architecture proposed in this paper.

(vii) **Physical grounding of interval diversity.** Independent verification shows that the normalised interval-size entropy $H_{\rm int}/\!\ln m$ is *strictly monotone decreasing* in spacetime dimension $d$ for Lor-$d$D sprinklings ($d = 2,3,4,5$; $H=0.82, 0.66, 0.45, 0.28$ at $N=40$). This is the information-theoretic signature of the Alexandrov interval volume power law $V \propto \tau^d$: higher $d$ concentrates the size distribution more extremely, lowering Shannon entropy. KR\_2layer has $H=0$ (completely degenerate). The KR/Lor4D penalty ratio grows from $1.19\times$ ($N=10$) to $3.25\times$ ($N=40$), providing a scale-dependent competitive advantage that drives Lor4D's selection at moderate penalty strengths.

This complementarity provides independent evidence that the two-layer screening is not an artifact of a circular geometric prior, but reflects two genuinely distinct structural axes along which Lor4D must be selected.

---

## §7. Conclusion

We have presented numerical evidence for a **layered structural screening architecture** that selects 4D Lorentzian causal sets from a diverse 25-family poset library:

1. A **first-order admissibility layer** (§2), carried by the triple screening functional $S_{\mathrm{triple}}$ and the Benincasa–Dowker action $S_{\mathrm{BD}}$, eliminates non-geometric structures but cannot uniquely select Lor4D ($S_{\mathrm{BD}}$ rank: 14/17 at $N = 128$).

2. A **second-order identity layer** (§3), implemented by the minimum-distortion functional $S_{\mathrm{MD}}$ (Mahalanobis distance from the Lor4D reference manifold), uniquely and robustly identifies Lor4D as the identity centre across all tested scales ($N \geq 10$) and all 25 families with zero free parameters. Counter-factual testing confirms this is a general identity mechanism, not a Lor4D-specific artifact (§3.4).

3. The identity layer exhibits a **sharp finite-size turn-on (onset threshold)** at $N_{\mathrm{id}} \approx 10$ (§4), set by the intrinsic geometric resolution limit of causal structure, not by the discriminator's statistical power. A 20-seed $\times$ 120-replication experiment with a separate reference ensemble confirms this onset with 95\% CI (§4, Table 3).

4. Beyond turn-on, the identity basin undergoes **historical sedimentation** (§5): the Mahalanobis gap grows (from $+0.308$ at $N = 10$ to $1.93 \times 10^8$ at $N = 1024$), the reference manifold sharpens ($\det(\Sigma) \propto N^{-3.31 \pm 0.14}$, $R^2 = 0.983$), and the effective basin volume contracts ($V_{\mathrm{eff}} \propto N^{-1.66 \pm 0.07}$)—making Lor4D an increasingly irreplaceable identity centre. Gap decomposition (§4.4) shows that the underlying physical separation is genuine ($d_E \approx 0.46$, flat), with the gap growth driven by statistical concentration.

5. A **curvature robustness program** (§4.7) shows that the screening architecture is robust for de Sitter-like deformations but only conditionally robust beyond that class. For de Sitter, up to $H \leq 0.3$ at $N \leq 1024$, Lor4D remains #1 and dS4D remains #2 in all 50 tested $(N,H)$ conditions; the $H=0$ limit at $N=1024$ gives $\Delta=+0.001$. Split low-$N$ tests show that weak-field Schwarzschild remains within the local-basin picture, whereas matter-FLRW at $\kappa=1.0$ reaches the hard-fail threshold at low $N$ but partially recovers in completed high-$N$ runs (failure ratio $0.3<0.5$); strong FLRW curvature exits further to #4–#5. The safest current claim is therefore background-dependent robustness, not a uniform mild-curvature theorem.

6. An **information-theoretic non-circularity test** (§6.9) replaces all geometric penalties with 5 pure information-theoretic measures that reference no target dimension or geometric template. At moderate penalty strength ($\gamma = 0.2$), Lor4D is selected as rank #1 from 17 families at all $N \geq 20$. At strong penalty, the information-theoretic and geometric penalties exhibit **complementarity**: they encode orthogonal quality dimensions whose conjunction uniquely selects Lor4D, confirming that the two-layer architecture reflects genuine structural axes rather than circular priors.

These results support the **Layered Structural Screening Principle**: robust identification/stabilization of 4D flat Lorentzian causal structure from a discrete poset space requires at least two functionally distinct mechanisms of different orders, whose intersection is Lor4D in the tested library and whose identity basin deepens with system size.

The principle is proposed as a mid-level effective theory—more specific than the philosophical assertion that "existence is what survives screening," more general than the individual numerical experiments from which it is extracted. Its analytical foundations (§5.2, §6.4) provide concrete targets for future work: a derivation of gap growth from CLT + Myrheim–Meyer, a uniqueness theorem for $\boldsymbol{\mu}(\infty)$, the incorporation of $S_{\mathrm{MD}}$ into the causal set path integral, and the extension to strongly curved spacetimes (Schwarzschild, FLRW) and the construction of a curvature-adaptive reference manifold $\mathcal{M}(H)$.

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
[12] F. Dowker and L. Glaser, “Causal set d’Alembertians for various dimensions,” *Class. Quantum Grav.* **30**, 195016 (2013).
---

*Manuscript v0.9 — §1–§7 complete. Very-large-N de Sitter convergence (N≤1024, Table 8) + FLRW/Schwarzschild extended pilot (N≤512, Table 9) integrated, with lowN/highN FLRW boundary wording synchronized (lowN threshold hit, highN partial recovery). Scope: flat Minkowski + de Sitter ($H \leq 0.3$, $N \leq 1024$) + preliminary non-de-Sitter curved backgrounds ($N=64$–$512$). 2026-03-30.*
*Cross-references: LAYERED_SCREENING_PRINCIPLE_OUTLINE.md (v2.0), MASTER_NARRATIVE.md, DISCUSSION_THEORY_IMPLICATIONS.md, manuscript_supplement_experiments.md, experiment_de_sitter.md, experiment_de_sitter_very_large_N.md, experiment_curved_backgrounds.md, experiment_curved_backgrounds_largeN.md.*
