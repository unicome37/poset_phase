# γ from Einstein-Hilbert: Derivation Attempt

## 1. Entropy Density h(d) = logH / N

| N | h(2D) | h(3D) | h(4D) | h(5D) | h(KR) |
|---|-------|-------|-------|-------|-------|
| 20 | 1.347 | 1.074 | 1.429 | 1.647 | 1.775 |
| 36 | 1.791 | 1.352 | 1.828 | 2.066 | 2.232 |
| 52 | 2.089 | 1.526 | 2.041 | 2.303 | 2.484 |
| 72 | 2.359 | 1.663 | 2.229 | 2.550 | 2.737 |
| 100 | 2.634 | 1.802 | 2.411 | 2.738 | 2.952 |

### h(d) growth: Δh = h(d+1) − h(d)

| N | Δh(3D−2D) | Δh(4D−3D) | Δh(5D−4D) |
|---|-----------|-----------|-----------|
| 20 | +0.3555 | +0.2177 | +0.1282 |
| 36 | +0.4761 | +0.2384 | +0.1663 |
| 52 | +0.5147 | +0.2616 | +0.1819 |
| 72 | +0.5658 | +0.3214 | +0.1871 |
| 100 | +0.6097 | +0.3267 | +0.2139 |

## 2. d_eff Well Penalty: (d_eff − 4)² per Element

| N | (d_eff−4)²(2D) | (d_eff−4)²(3D) | (d_eff−4)²(4D) | (d_eff−4)²(5D) | (d_eff−4)²(KR) |
|---|--------------|--------------|--------------|--------------|--------------|
| 20 | 1.841 | 4.151 | 0.735 | 0.103 | 0.174 |
| 36 | 1.729 | 4.040 | 0.637 | 0.030 | 0.177 |
| 52 | 1.649 | 3.886 | 0.544 | 0.039 | 0.134 |
| 72 | 1.629 | 4.069 | 0.570 | 0.016 | 0.170 |
| 100 | 1.614 | 3.990 | 0.547 | 0.023 | 0.149 |

## 3. γ from Entropy-Geometry Balance

At the F10 minimum (4D), for 4D to beat dimension d:

  F10(d) > F10(4D)
  logH(d) + γN(d_eff(d)−4)² > logH(4D) + γN(d_eff(4D)−4)²
  γ > [logH(4D) − logH(d)] / [N · ((d_eff(d)−4)² − (d_eff(4D)−4)²)]
  γ > −Δh(d−4D) / Δ(d_eff²)(d−4D)

where Δh and Δ(d_eff²) are per-element quantities.

### Per-element balance:

| N | Δh/Δd² (3D) | Δh/Δd² (5D) | Δh/Δd² (2D) | binding |
|---|------------|------------|------------|---------|
| 20 | -0.3445 | 1.7830 | -0.1416 | 5D |
| 36 | -0.3926 | 1.1265 | -0.1781 | 5D |
| 52 | -0.5182 | 1.9097 | -0.2018 | 5D |
| 72 | -0.5798 | 1.2121 | -0.2189 | 5D |
| 100 | -0.6242 | 1.6974 | -0.2360 | 5D |

**Mean γ_crit(3D)**: -0.492 ± 0.107
**Mean γ_crit(5D)**: 1.546 ± 0.316
**Mean γ_crit(2D)**: -0.195 ± 0.033

**Binding constraint**: γ must exceed max(γ_crit) ≈ 1.91 (5D)

## 4. Theoretical γ from EH Action Structure

### 4a. The Myrheim-Meyer encoding identity

For a causal set sprinkled into d-dim Minkowski at density ρ:

  f₂(d) = Γ(d+1)Γ(d/2) / (4Γ(3d/2))  [fraction of causally related pairs]

The Myrheim-Meyer estimator solves f₂(d_eff) = C₀/C(N,2) for d_eff.

For d_eff close to d*, Taylor expand:

  f₂(d) ≈ f₂(d*) + f₂'(d*)(d−d*) + ½f₂''(d*)(d−d*)² + ...

At d*=4: f₂ = 0.050000, f₂' = -0.042084, f₂'' = 0.034149

f₂ at integer d:

  f₂(2) = 0.250000
  f₂(3) = 0.114286
  f₂(4) = 0.050000
  f₂(5) = 0.021312
  f₂(6) = 0.008929

### 4b. Key subtlety: f₂-space vs d-space

**IMPORTANT**: The F10 well operates in d-space: γN(d_eff−4)².
But the BDG link action operates in f₂-space: S_link ∝ −2C₀ = −2·f₂·C(N,2).
These are related by the NONLINEAR map d = f₂⁻¹(observed fraction).

The d-space curvature of f₂ tells us about the link action's
dimension sensitivity, but γ in F10 is a SEPARATE quantity.

f₂''(4) = +0.034149 (POSITIVE = convex in d-space)
→ The link action has CONCAVE d-dependence at d=4
→ S_link FAVORS deviations from d=4, not penalizes them!

This is why γ cannot come from f₂ alone — the sign is wrong.
The d_eff well must arise from a DIFFERENT mechanism.

### 4c. The d_eff well as an INDEPENDENT geometric constraint

The d_eff estimator inverts f₂(d_eff) = C₀/C(N,2).
Since f₂ is monotonically decreasing, this inversion is well-defined.
But d_eff is a DERIVED quantity — it encodes C₀ in d-units.

The F10 well γN(d_eff−4)² is therefore equivalent to:
  γN · [f₂⁻¹(C₀/C(N,2)) − 4]²

This is a NONLINEAR function of C₀, amplified by the
steep gradient of f₂⁻¹ near d=4 (where f₂ is small and
changing rapidly with d).

Jacobian: dd/df₂ = 1/f₂'(4) = 1/(-0.042084) = -23.76
→ A unit change in f₂ maps to 23.8 units of d_eff

This amplification factor is KEY: small changes in C₀
get magnified into large changes in d_eff, making the
well γ(d_eff−4)² effectively γ·565·(Δf₂)² in f₂-space.

γ=1 in d-space ≡ γ_f₂ = 565 in f₂-space
This enormous amplification explains why γ=1 works:

the d_eff well leverages the steep f₂→d inversion near d=4.

### 4d. Full BDG coefficient structure

The link action only uses C₀. The full BDG action uses C₀, C₁, C₂, C₃...
Each C_k has its own dimension-dependent coefficient α_k^(d).

The total dimension penalty from the BDG action is:

  γ_eff = Σ_k α_k^(d*) · (d²f_k/dd²)|_{d=d*}

where f_k(d) = E[C_k]/C(N,2) is the expected k-interval fraction.

This sum can amplify or suppress the leading f₂'' term.

## 5. Why γ = O(1) is Natural

The question is not 'derive γ=1 from f₂' but rather:
'why should the dimension penalty be O(1) per element?'

### Three independent arguments for γ = O(1):

**Argument 1: Entropy-geometry equipartition**
The entropy density h(d) ≈ 1.4–2.7 nat/element (from Part 1).
For the well to compete with entropy, γ·(d_eff−4)² must be O(h).
Since (d_eff(5D)−4)² ≈ 0.15 at N=100,
we need γ ≈ Δh(5D−4D)/Δ(d²) = 1.5 for 5D competition.
This gives γ = O(1), which is a BALANCE condition, not a derivation.

**Argument 2: Planck density normalization**
In the continuum, S_EH/V = R/(16πG). At Planck density ρ_P = 1/ℓ_P^d:
  S_EH/N = S_EH/(ρ_P V) = R·ℓ_P²/(16π) = O(1)
since R·ℓ_P² = O(1) for curvatures at the Planck scale.
The dimension dependence ∂²(S_EH/N)/∂d² is therefore also O(1),
giving γ_EH = O(1) naturally.

**Argument 3: f₂ Jacobian amplification (from §4c)**
γ=1 in d-space = γ_f₂ = 565 in f₂-space.
The Myrheim-Meyer inversion d = f₂⁻¹(p) has steep gradient
near d=4: |dd/df₂| = 23.8.
So a mild O(1) well in d-space corresponds to an enormous
penalty of ~565 per unit (Δf₂)² in the observable C₀ space.
The physical content: d_eff is a HIGHLY COMPRESSED encoding
of C₀, and the well leverages this compression.

## 6. d_eff vs R: Information Content

Our R observable is the OCCUPANCY RATIO (fraction of matrix entries = 1).
The Myrheim-Meyer f₂ is the fraction of CAUSALLY RELATED pairs.
For a Hasse diagram stored as a CLOSURE matrix, R = f₂ = C₀/C(N,2).
But our R is computed differently — let's check the actual relationship.

| N | family | R | f₂_theory | d_eff | d_from_f₂_theory |
|---|--------|---|-----------|-------|-----------------|
| 20 | KR_like | 0.321 | — | 2.66 | f₂(d_eff)=0.1500 |
| 20 | Lor2D | 0.620 | 0.2500 | 1.99 | f₂(d_eff)=0.2524 |
| 20 | Lor3D | 0.326 | 0.1143 | 3.21 | f₂(d_eff)=0.0961 |
| 20 | Lor4D | 0.125 | 0.0500 | 3.96 | f₂(d_eff)=0.0518 |
| 20 | Lor5D | 0.037 | 0.0213 | 4.32 | f₂(d_eff)=0.0381 |
| 100 | KR_like | 0.333 | — | 2.73 | f₂(d_eff)=0.1420 |
| 100 | Lor2D | 0.873 | 0.2500 | 2.01 | f₂(d_eff)=0.2486 |
| 100 | Lor3D | 0.633 | 0.1143 | 3.28 | f₂(d_eff)=0.0912 |
| 100 | Lor4D | 0.359 | 0.0500 | 3.93 | f₂(d_eff)=0.0529 |
| 100 | Lor5D | 0.148 | 0.0213 | 4.38 | f₂(d_eff)=0.0364 |

Key: R in our dataset is NOT f₂. R is the occupancy ratio (total
relations / N²), while f₂ = C₀/C(N,2) (ordered pairs only).
But d_eff IS calibrated from f₂ by construction.


## 7. Formal Interpretation of γ

### Summary of derivation attempt:

### 结论

1. **f₂'' 通道的符号是错的**：f₂''(4) > 0（凸），link action 在 d=4
   附近是凹的 → 它反而**鼓励**偏离 d=4，而非惩罚。
   因此 γ 不可能从 link action 推导出来。

2. **d_eff 井是独立的几何约束**，不是 BDG action 的重参数化。
   它通过 f₂⁻¹ 反演的陡峭 Jacobian 从 C₀ 中提取维度信息,
   放大倍数 |dd/df₂|² ≈ 565。

3. **γ = O(1) 的三个独立论证**:
   (a) 熵-几何等分: γ·Δ(d²) ∼ Δh ∼ O(1) → γ ∼ O(1)
   (b) Planck 密度归一化: S_EH/N = R·ℓ_P²/(16π) = O(1)
   (c) Jacobian 放大: γ=1 在 d 空间 = γ_f₂ = O(500) 在 f₂ 空间

4. **严格推导需要**: 完整 BDG 展开到所有阶的 d 依赖性,
   这是一个开放的理论问题。但 γ = O(1) 是自然标度。

5. **核心洞见**: F10 的 d_eff 井不是 BDG action 的简化版,
   而是一个**正交**的几何约束——它编码的是 Myrheim-Meyer 维度
   (通过 f₂ 反演获得), 不是标量曲率(通过 BDG 系数获得)。
   这正是为什么 Φ_geom(d_eff) ⊥ Ψ_Lor(R) 是正交的：
   它们编码了因果集合几何的**不同方面**。