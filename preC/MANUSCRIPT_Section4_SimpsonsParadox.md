# Prediction C — Manuscript Draft: Section 4

## 4. Simpson's Paradox: Poset Size as the Dominant Sign-Determining Confound

### 4.1 The Anomaly

The most counterintuitive result in Tier 1 is the *sign* of the naïve partial correlation. Controlling for antichain width, comparable fraction, and geometric dimension — the standard controls from Prediction B — gives

$$r_\text{partial}(\text{HII},\;\log H \mid \text{aw}, \text{cf}, \text{geo\_dim}) \;=\; +0.336 \qquad (p = 0.0005).$$

This is positive: higher HII appears to be associated with *higher* entropy, in direct contradiction to the prediction. If accepted at face value, this result would falsify Prediction C. This section shows that the positive sign is a Simpson's Paradox driven predominantly by the poset size $N$ — the dominant confound in Tier 1 — and that the true structural relationship is negative.

### 4.2 Stepwise Confound Identification

To isolate the responsible variable, we systematically vary the control set and observe the sign of the partial correlation.

**Table 9.** Stepwise diagnosis of the Simpson's Paradox ($n = 320$).

| Step | Control variables | $r_\text{partial}$ | Sign | Interpretation |
|------|------------------|---------------------|------|----------------|
| 0 | None (raw Pearson) | $+0.12$ | $+$ | Weak positive; $N$-scaling dominates |
| 1 | aw, cf, geo_dim | $+0.336$ | $+$ | **Simpson's Paradox** — controls absorb wrong variance |
| 2 | **$N$ only** | $\mathbf{-0.578}$ | $\mathbf{-}$ | **Sign flip** — $N$ resolves the paradox |
| 3 | $N$ + family dummies | $-0.250$ | $-$ | Same direction; attenuated by within-family noise |
| 4 | $N$ + aw + cf + geo_dim | $-0.434$ | $-$ | Full controls with $N$: still negative |

The diagnostic pattern is clear: $N$ is the dominant sign-determining confound. Adding $N$ at step 2 produces the largest absolute change in $r$ ($+0.336 \to -0.578$, a swing of $0.914$). No other single control achieves a comparable effect.

### 4.3 Why $N$ Is a Confound

$N$ simultaneously drives both HII and $\log H$, but with different magnitudes:

1. **$N \to \log H$**: The number of linear extensions grows combinatorially with $N$. Across all 8 families, $r(N, \log H) > 0.91$. This is a first-order scaling effect: more elements create exponentially more compatible orderings.

2. **$N \to \text{HII}$**: Larger posets generically admit deeper layer structures, so HII tends to increase with $N$ across families (e.g., Lor2D has $r(N, \text{HII}) \approx +0.65$). This is a second-order structural effect: more elements permit — but do not guarantee — deeper hierarchy.

3. **Cross-$N$ contamination**: Without controlling for $N$, the analysis conflates two effects:

   - A **scaling effect** (positive): larger $N$ $\Rightarrow$ both HII $\uparrow$ and $\log H \uparrow$.
   - A **structural effect** (negative): at fixed $N$, deeper hierarchy $\Rightarrow$ fewer compatible orderings $\Rightarrow$ $\log H \downarrow$.

   Because the scaling effect operates on a larger dynamic range than the structural effect (entropy grows super-linearly in $N$, while HII grows sub-linearly), the aggregate relationship is dominated by scaling and appears positive.

### 4.4 An Ecological Correlation

This is a textbook instance of an *ecological fallacy*: the aggregate-level trend (positive) is the reverse of the within-group trend (negative). The "groups" here are the fixed-$N$ slices.

The fixed-$N$ Pearson correlations (Table 6) confirm this directly: at every $N$ from 10 to 16, the cross-family $r(\text{HII}, \log H)$ is negative ($-0.86$ to $-0.53$). The aggregate positive trend exists *only* because families with higher HII (e.g., Lor2D, Interval Order) also have higher *baseline* entropy at every $N$, and larger $N$ amplifies both quantities.

### 4.5 Physical Analogy: Crystallization at Fixed Temperature

An intuitive analogy clarifies the structure of the paradox:

- **Temperature** (analogue of $N$): A hotter melt has more thermal kinetic energy (higher entropy baseline) and more molecular mobility (higher structural complexity baseline).
- **At fixed temperature**: Crystalline solids have lower entropy than liquids because their internal structure constrains the microstate count.
- **Across temperatures**: Comparing a cold crystal to a hot gas shows the crystal having *lower* entropy with *more* order — but this correlation disguises the fact that temperature (system size) is the true driver.

The Prediction C analogy:
- **$N$** plays the role of temperature.
- **HII** plays the role of crystalline order.
- **$\log H$** plays the role of entropy.

Studying the structure–entropy relationship requires *fixing the temperature* (controlling for $N$).

### 4.6 Verification: Family Controls Alone Are Insufficient

A natural question: can family dummies substitute for $N$ control? The answer is no:

$$r_\text{partial}(\text{HII},\;\log H \mid \text{family dummies only}) \;=\; +0.148 \qquad (p = 0.043).$$

The sign remains positive, because within each family, $N$ still acts as a confound. Family membership absorbs between-family variance but not between-$N$ variance. Only explicit $N$ control (or the matched-pair design of Tier 2) resolves the paradox.

### 4.7 Why the Original Controls Fail

The Prediction B controls (aw, cf, geo_dim) were designed for a different purpose: controlling for structural similarity in the action-score comparison. They fail here because:

1. **Collinearity with $N$**: `antichain_width` scales approximately as $N^{1/2}$; `comparable_fraction` weakly depends on $N$; `geo_dim_eff` has family-dependent $N$-scaling. These controls absorb some $N$-variance but not enough to overcome the dominant $N \to \log H$ pathway.

2. **Wrong regression target**: The Prediction B controls were not chosen to purge $N$-scaling from an entropy residual — they were chosen to equalize structural baselines across families. Re-purposing them for the HII–$\log H$ test introduces an inadequate adjustment that paradoxically amplifies the confound (from $r = +0.12$ to $r = +0.336$).

### 4.8 Methodological Implication

The Simpson's Paradox has a general lesson for cross-size poset studies:

> **Raw correlations between structural observables and entropy-like quantities should not be trusted without explicit $N$ control.** The super-linear growth of $\log H$ with $N$ creates confounds that can reverse structural relationships.

This conclusion extends beyond the present analysis: any study comparing poset families across different sizes must either (i) include $N$ as a regression control, (ii) stratify by $N$, or (iii) use a matched-pair design at identical $N$.

The matched-pair Tier 2 design was motivated precisely by this concern, and its results (Section 3.3) confirm the negative direction without any $N$-related ambiguity.
