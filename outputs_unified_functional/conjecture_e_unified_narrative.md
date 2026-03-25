# Conjecture E: Unified Narrative of §4.1.21–42

# 猜想 E：§4.1.21–42 统一叙事

> **One-sentence summary**: Twenty-two numerical experiments, comprising over 10,000 causal-set sprinklings across dimensions $d=2,3,4$ and element counts $N=64$–$2048$, progressively establish that a finite causal set encodes spacetime curvature at three distinct levels — an admissibility wall, a first-order bulk recovery of expansion rate $H$, and a constructible second-order bridge to scalar curvature $R$ via two-step calibration-then-squaring ($R^2 = 0.987$/$0.9996$) — with the signal triple-confirmed in non-constant-curvature backgrounds (global, local, and N-scaling). DDT Condition C2 (non-uniform curvature) partially answered: 3+1D Schwarzschild yields 4/32 beyond density; the antichain channel requires $R \neq 0$ (primarily responds to scalar curvature, not Weyl).

> **一句话摘要**：二十二项数值实验、超过 10,000 次因果集 sprinkling（覆盖 $d=2,3,4$ 维、$N=64$–$2048$ 元素），逐步建立了有限因果集在三个层次上编码时空曲率的证据——准入壁（wall）、膨胀率 $H$ 的一阶 bulk 恢复、以及通过两步校准-平方法可构造的到标量曲率 $R$ 的二阶桥接（$R^2 = 0.987$/$0.9996$）——且该信号在非常曲率背景中已获三重确认（全局 + 局域 + N-scaling）。DDT 条件 C2（非均匀曲率）部分回答：3+1D Schwarzschild 4/32 超越密度；反链通道需 $R \neq 0$（主要响应标量曲率，非 Weyl）。

---

## Overview: The Five Phases / 概览：五个阶段

The experimental programme divides naturally into five phases:

| Phase | Experiments | Central Question | Outcome |
|-------|------------|-----------------|---------|
| **I. Wall & Density** | §4.1.21–26b | How does a causal set encode curvature at the coarsest level? | Sigmoid wall = curvature upper bound; all interval statistics are density-dominated |
| **II. DDT Escape** | §4.1.27–30 | Can any purely causal observable carry curvature information *beyond* density? | Yes: spectral (B_ℓ eigenvalues) and transverse (antichain structure) channels escape DDT |
| **III. Unification** | §4.1.31–33 | Are the escape channels independent, and what geometric quantity do they track? | Single post-density DoF; target = $H$ at $d=4$; EH bridge = continuum limit |
| **IV. Generalization** | §4.1.34–35 | Does the signal survive outside de Sitter? | Yes: global + local confirmation in power-law FRW (anti-monotone direction) |
| **V. EH Bridge** | §4.1.36–40 | Can the $H \to R \sim H^2$ bridge be constructed at finite $N$? | T3 two-step squaring positive ($R^2=0.987$/$0.9996$); T1/T2/T4 negative; second-order = squared amplitude of first-order deviation |

实验纲领自然分为五个阶段：

| 阶段 | 实验 | 核心问题 | 结论 |
|------|------|---------|------|
| **I. Wall 与密度** | §4.1.21–26b | 因果集如何在最粗层次编码曲率？ | Sigmoid wall = 曲率上界；所有区间统计量均受密度主导 |
| **II. DDT 逃逸** | §4.1.27–30 | 是否有纯因果可观测量能携带密度之外的曲率信息？ | 有：谱通道（$B_\ell$ 特征值）和横向通道（反链结构）逃逸 DDT |
| **III. 统一** | §4.1.31–33 | 逃逸通道是否独立？它们追踪的几何量是什么？ | 单一 post-density DoF；$d=4$ 靶标 = $H$；EH 桥接 = 连续极限 |
| **IV. 推广** | §4.1.34–35 | 信号在 de Sitter 之外是否存活？ | 是：幂律 FRW 中全局 + 局域确认（反单调方向） |
| **V. EH 桥接** | §4.1.36–40 | 有限 $N$ 下能否构造 $H \to R \sim H^2$ 桥接？ | T3 两步平方阳性（$R^2=0.987$/$0.9996$）；T1/T2/T4 阴性；二阶 = 一阶偏移之平方幅度 |

---

## Phase I: Wall Establishment & Density Dominance (§4.1.21–26b)

## 第一阶段：Wall 建立与密度主导（§4.1.21–26b）

### The Sigmoid → EH Correspondence (§4.1.21)

The first experiment tested whether the sigmoid wall term $\sigma((R - R_c)/w)$ in the unified functional $\mathcal{F}_7$ has a physical counterpart in the Einstein–Hilbert action.

**Key discovery**: The occupancy ratio $R$ (fraction of non-link causal pairs) **anti-correlates** with expansion rate $H$ — opposite to naïve expectation. At $d=4$, $N=1024$: $R(H=0) = 0.709$ vs $R(H=2) = 0.078$.

**Physical mechanism**: De Sitter expansion compresses causal diamonds → fewer causal pairs per element → survivors are mostly links → $f_\text{link} \uparrow$ → $R = 1 - f_\text{link} \downarrow$.

**Corrected correspondence**:
$$\sigma\!\left(\frac{R - R_c}{w}\right) \;\longrightarrow\; \Theta(S_\max - S_\text{EH})$$
The wall is a **curvature upper bound** (admissibility gate), not a lower bound. Flat spacetime ($R \gg R_c$) is "admitted"; strongly curved spacetime ($R < R_c$) is "rejected."

**Convergence evidence**: $\text{std}(R)$ shrinks monotonically with $N$ (all $d$, $H$), and $|\sigma(H=0) - \sigma(H>0)|$ increases — the sigmoid sharpens toward a step function in the continuum limit.

> **E-wall status**: Established (90–95% confidence). The mechanism, direction, scaling, and convergence are all confirmed.

### 核心发现
占有率 $R$ 与膨胀率 $H$ **反相关**——这修正了朴素预期。物理原因是 de Sitter 膨胀压缩因果钻石。Sigmoid wall 编码的是**曲率上界**（准入门槛）。

---

### The Density Dominance Theorem (§4.1.22–26b)

A systematic investigation of 26 causal-interval statistics ($\{C_k\}$-based metrics) across 6 experiments established a foundational negative result:

**§4.1.22 (BD ratio scan)**: All density-based metrics (bd_ratio, M1, M2, M4, M6, occupancy $R$) anti-correlate with $H$ at $|\rho| \approx 0.975$–$0.980$. No pure causal-interval metric positively tracks curvature.

**§4.1.23 (Curvature-sensitive BD)**: Two camps identified:
- **↓ Camp** (9/9 significant, $|\rho| \approx 0.98$): All measure causal connectivity density
- **↑ Camp** (6/9, $|\rho| \approx 0.69$): Only `bdg_d2c` ($N - 2C_0 + 2C_1$), works for $d \geq 3$ only

**§4.1.24 (Bridge candidates)**: KL divergence (7/9) strongest positive candidate, but requires flat-space reference distribution.

**§4.1.25 (PCA)**: PC1 = 88.6% (density), PC2+ = 11.4% (apparent "shape").

**§4.1.26 (PC-Space proxy)**: Two-line architecture: Line A = total density → wall, Line B = differential decay rate → bdg_d2c.

**§4.1.26b (Density Degeneracy Theorem)**: The critical correction — after removing $\Sigma C_k$, **all** residual correlations drop to $|\rho| < 0.24$ (none significant). `bdg_d2c`'s positive tracking is not an independent "shape signal" but differential density decay ($C_0/C_1$ ratio shift). The two-line architecture collapses: **both lines are density**, with zero independent shape degrees of freedom within $\{C_k\}$.

> **Conclusion**: The Alexandrov interval statistics $\{C_k\}$ encode curvature **exclusively** through total causal density. Any beyond-density signal must come from **outside** the $\{C_k\}$ family. This is the Density Dominance Theorem (DDT).

### 密度简并定理
经过 6 项实验的系统排查：$\{C_k\}$ 族内部不存在独立于总密度的曲率自由度。bdg_d2c 的正相关源于差分密度衰减（$C_0/C_1$ 比值偏移），而非独立形状信号。**任何 beyond-density 信号必须来自 $\{C_k\}$ 族之外。**

---

## Phase II: DDT Escape — The Bulk Breakthrough (§4.1.27–30)

## 第二阶段：DDT 逃逸——Bulk 的突破（§4.1.27–30）

Phase I established that $\{C_k\}$ (Alexandrov intervals) are a dead end for beyond-density curvature information. Phase II searched for observables that **escape** the DDT.

### Channel I.1: Spectral Escape via B_ℓ (§4.1.27)

The Sorkin d'Alembertian $B_\ell$ is an $N \times N$ matrix built from BDG coefficients applied to interval counts. Its **constant-field projection** (= BDG action per element) is a linear combination of $\{C_k\}$ and therefore density-dominated (0/9 slices pass — confirming DDT). But its **eigenvalue spectrum** encodes topological correlation structure that goes beyond "how many pairs have $k$ elements between them."

**Density-residual analysis** (the critical test):
| Dimension | Beyond-density features | Peak $|\rho_\text{resid}|$ |
|-----------|----------------------|--------------------------|
| $d=2$ | 0/6 | < 0.29 |
| $d=3$ | 2/6 | 0.569 |
| $d=4$ | **4/6** | **0.703** |
| **Total** | **6/18** | |

The spectral channel carries curvature information beyond density, with strength increasing with dimension. **DDT Condition C1 escaped for the first time.**

### 通道 I.1：$B_\ell$ 谱逃逸（§4.1.27）
常数场投影 0/9 确认 DDT；但**特征值谱**在去密度后保留 6/18 的 beyond-density 信号（$d=4$: 4/6，峰值 $|\rho_\text{resid}| = 0.703$）。DDT 条件 C1 首次被逃逸。

---

### Channel I.2: Transverse Escape via Antichain Structure (§4.1.28)

While $\{C_k\}$ and $B_\ell$ projections measure **longitudinal** structure (causal intervals = time direction), antichains measure **transverse** structure (spacelike hypersurfaces = space direction). This is a physically orthogonal channel that DDT does not constrain.

Seven antichain/layer statistics were extracted: Dilworth width $w_\max/N$, number of layers, mean/std/cv of layer widths, max layer width ratio, and layer entropy.

**Density-residual analysis**:
| Dimension | Beyond-density features | Peak $|\rho_\text{resid}|$ |
|-----------|----------------------|--------------------------|
| $d=2$ | **7/7** | **0.817** |
| $d=3$ | **7/7** | 0.755 |
| $d=4$ | **7/7** | 0.740 |
| **Total** | **21/21** | |

**100% pass rate across all dimensions.** The antichain channel is the strongest DDT escape path found. It covers $d=2$ (where $B_\ell$ spectral fails), achieves stronger residuals, and is physically interpretable: de Sitter expansion narrows causal diamonds but **widens** spatial slices — a geometric effect beyond density.

> **Verdict**: DDT C1 thoroughly escaped via the transverse channel. The antichain channel is the primary evidence carrier for E-bulk.

### 通道 I.2：反链结构横向逃逸（§4.1.28）
反链/层结构统计量 **21/21 全部通过去密度检验**（峰值 0.817）。全维度覆盖。物理上，反链测量的是"类空超曲面"（横向），与 $\{C_k\}$ 的"因果区间"（纵向）正交——DDT 不约束横向。

---

### Failed Paths: Schwarzschild (§4.1.29) and Chain Statistics (§4.1.30)

**§4.1.29 (Schwarzschild)**: 1+1D Schwarzschild sprinkling tested DDT Condition C2 (non-uniform curvature). Only 2/27 features passed at $N=512$. The weak result is expected: 1+1D Schwarzschild is Ricci-flat, and $d=2$ has no Weyl transverse degrees of freedom. **Verdict: weak positive, C2 remains open.**

**§4.1.30 (Chain statistics $N_k$)**: $k$-chain counts ($N_2$–$N_5$) belong to the same longitudinal sector as $\{C_k\}$. Result: $d=2,3$: 0/72 (complete density absorption); $d=4$: 6/36 nominally pass but signal vanishes at $N=512$ — finite-size artifact. **Path I.4 unreliable.**

### 失败路径
- **Schwarzschild**（§4.1.29）：1+1D Ricci 平坦背景中仅 2/27 通过——弱阳性。
- **链统计 $N_k$**（§4.1.30）：同属纵向扇区，密度吸收有效——路径不可靠。

---

## Phase III: Unification & Geometric Target (§4.1.31–33)

## 第三阶段：统一与几何靶标识别（§4.1.31–33）

### Dual-Channel Unification (§4.1.31)

The two successful escape channels — spectral ($B_\ell$) and transverse (antichain) — were extracted **simultaneously** from the same poset in 360 de Sitter realizations.

**Cross-channel residual correlations** (after removing density):
| Dimension | $N=128$ | $N=256$ | $N=512$ | $N$-trend |
|-----------|---------|---------|---------|-----------|
| $d=3$ | 0.49 | 0.71 | **0.86** | $\rho(N, |\rho_\text{cross}|) = +1.0$ |
| $d=4$ | 0.48 | 0.66 | **0.85** | $\rho(N, |\rho_\text{cross}|) = +1.0$ |

**Joint $\Delta R^2 \approx 0$**: Adding the second channel to a regression already containing the first provides essentially zero additional explanatory power. The two channels are **not independent** — they are complementary projections of a **single** post-density geometric degree of freedom.

> **Conclusion**: At $d \geq 3$, large $N$, the antichain and spectral channels converge to a single geometric DoF. The cross-channel correlation monotonically strengthens with $N$, indicating convergence to a physical quantity in the continuum limit.

### 双通道统一（§4.1.31）
去密度后的交叉通道残差在 $d=3,4$ 达到 $|\rho| \to 0.86$，且 $N$ 增大时单调增强（$\rho(N, |\rho_\text{cross}|) = +1.0$）。联合 $\Delta R^2 \approx 0$。**两条通道收敛到同一 post-density 几何自由度。**

---

### Geometric Target Identification (§4.1.32)

What continuous geometric quantity does this single post-density DoF correspond to? A Pearson $R^2$ power-law scan ($\alpha \in [0.25, 8.0]$) was performed on pooled group-mean data (24 realizations per $H$ level — the most reliable estimate).

**$d=4$ results** (physical spacetime dimension):
| Feature | Optimal $\alpha$ | Interpretation |
|---------|-----------------|---------------|
| w_max_ratio (antichain) | 1.25 | ≈ $H^1$ |
| b1_std (spectral) | 1.00 | = $H^1$ |

**Both channels converge to $\alpha \approx 1$**: the target is **expansion rate $H$** (equivalently, the trace of extrinsic curvature $K = (d-1)H$), **not** scalar curvature $R = d(d-1)H^2$.

**Dimension dependence**: $d=3$: $\alpha \approx 1.75$–$2.50$; $d=2$: scattered (channels inconsistent). The $d=2$ anomaly is explained: $B_\ell$ has only 1 BDG coefficient layer, and antichain width has minimal dynamic range in 1 spatial dimension.

> **Physical content**: At finite $N$, causal-set observables are genuinely first-order in curvature. The EH action $S_\text{EH} \sim R = d(d-1)H^2$ is a second-order quantity; its recovery is a **continuum-limit** phenomenon, not a finite-$N$ algebraic operation.

### 几何靶标（§4.1.32）
$d=4$ 的两条通道一致收敛到 $\alpha \approx 1$：靶标是**膨胀率 $H$**（外曲率迹），而非标量曲率 $R = d(d-1)H^2$。离散可观测量在有限 $N$ 下是真正一阶的。

---

### The EH Bridge (§4.1.33)

Five candidate constructions were tested for lifting $\alpha$ from 1 ($H$-tracking) to 2 ($R$-tracking):
- C1: Residual squared
- C2: Cross-channel product
- C3: Feature squared then density-removed
- C4: Raw cross-product then density-removed
- C5: Quadratic regression diagnostic

**$d=4$ hit rate: 0/9** — no candidate achieves $\alpha \approx 2$.

**C5 diagnostic** (most informative): Fitting $\text{resid} = a H^2 + b H + c$ shows $H$ dominates at $d=4$ ($t_H = 2.9/3.3$, $t_{H^2} = 1.7/0.2$), while $H^2$ dominates at $d=2,3$ ($t_{H^2} > 2$). This dimension-dependent crossover independently confirms the $\alpha$-decreasing trend from §4.1.32.

> **Conclusion**: The bridge from $H$ to $R$ is not a missing algebraic trick — it is the **continuum limit itself**. The BDG theorem ($S_\text{BD} \to \int R \sqrt{g}\, d^4x$ as $N \to \infty$) completes the squaring. At finite $N$, first-order $H$-tracking is the correct physical content.

### EH 桥接（§4.1.33）
5 种候选构造在 $d=4$ 均未达 $\alpha \approx 2$（0/9）。C5 诊断确认 $d=4$ 残差中 $H$ 主导。**从 $H$ 到 $R$ 的桥接是连续极限本身——BDG 定理在 $N \to \infty$ 时完成平方。**

---

## Phase IV: Beyond de Sitter (§4.1.34–35)

## 第四阶段：超越 de Sitter（§4.1.34–35）

All Phase I–III results were obtained in de Sitter spacetime ($H = \text{const}$). A critical question remained: **is the signal a de Sitter symmetry artifact, or does it generalize to non-constant-curvature backgrounds?**

### Global Extension (§4.1.34 Phase A)

432 power-law FRW sprinklings ($a(t) = (t/t_0)^p$, $p \in \{0.4, 0.5, 0.67, 1.0, 1.5, 2.0\}$, $d=2/3/4$, $N=128/256/512$, 8 reps). Time-dependent expansion rate: $H(t) = p/t$.

**Global beyond-density results** (features correlated with $p^2$ after density removal):
| Dimension | Beyond-density | Total |
|-----------|---------------|-------|
| $d=2$ | 22/42 | Strong |
| $d=3$ | 10/42 | Moderate |
| $d=4$ | 5/42 | Weak |
| **Total** | **37/126** | |

The antichain channel (w_max_ratio) partially generalizes. The $B_\ell$ spectral channel does **not** generalize to power-law FRW. The $\alpha$ grid search does not converge to a clean power-law target.

### 全局推广（§4.1.34 Phase A）
幂律 FRW 中反链通道部分推广（37/126 去密度通过）。$B_\ell$ 谱不推广。

---

### Local H(t) Tracking (§4.1.35) — The Critical Experiment

§4.1.34 Phase B (early/late binning within a single sprinkling) yielded 0/54 in the "expected" direction — but this was diagnosed as a **methodological confound** due to unequal bin sizes from FRW volume weighting. §4.1.35 was designed to resolve this with two clean methods.

**Method 1 (Independent Patches)**: 768 sprinklings — for each $(d, N, p)$ combination, sprinkle independently at 4 different epoch centers $t_c \in \{0.2, 0.4, 0.6, 0.8\}$, each with $N$ elements in a narrow time window $[t_c - \delta, t_c + \delta]$. This guarantees equal element counts at each epoch.

**Result**: w_max_ratio anti-correlates with local $H(t_c) = p/t_c$ in **24/24** $(d, N, p)$ cells:
- All $|\rho| \approx 0.88$–$0.97$
- All $p < 10^{-3}$ (highly significant)
- **0/24 positive correlations** — 100% anti-monotone

**Method 2 (Density-Matched Sub-sampling)**: 384 matched-bin comparisons within single sprinklings, equalizing element counts between early and late bins by random sub-sampling.

**Result**: 0/192 cases where early bin (high $H$) shows higher w_max_ratio than late bin (low $H$). Sign test: $p < 10^{-58}$.

### The Anti-Monotone Insight

The naïve prediction was "higher $H$ → wider antichains." This was wrong. The correct physics:

1. **Higher local $H(t)$** → faster expansion within the patch → causal diamonds compressed → sparser poset
2. In a sparser poset, the maximum antichain is a larger fraction of $N$ (trivially, an empty poset is one giant antichain)
3. But **w_max_ratio = w_max / N** in a dense, well-structured poset (low $H$, late epoch) reflects meaningful transverse slicing — and is **higher** than in a sparse, near-trivial poset (high $H$, early epoch)

**This is exactly the wall mechanism operating locally**: higher curvature → sparser causal structure → anti-monotone observable shift. The sigmoid wall $\sigma((R - R_c)/w)$ is anti-monotone in $R$; the local patch experiment confirms this anti-monotone response at the **local** level.

### Reinterpretation of §4.1.34 Phase B

The §4.1.34 Phase B result (0/54 early > late) was originally diagnosed as "methodological confound." §4.1.35 reveals this diagnosis was **partially wrong**: the direction (late > early for w_max_ratio) is the **correct physical direction** — anti-monotone curvature response, not density artifact. The density imbalance in §4.1.34 amplified a real signal; the signal direction was already correct.

### Beyond-Density in Local Test

| Method | $d=2$ | $d=3$ | $d=4$ | Total |
|--------|-------|-------|-------|-------|
| Method 1 (patches) | 0/10 | 6/10 | 6/10 | **12/30** |
| Method 2 (matched) | 3/4 | 1/4 | 0/4 | 4/12 |

Method 1 (independent patches) is the cleaner test. **12/30 features survive density removal** at $d \geq 3$, dominated by w_max_ratio, mean_layer_width, and layer_width_std. The antichain channel carries genuine **local** curvature information beyond density, even in non-constant-$H$ backgrounds.

### 局域 $H(t)$ 追踪（§4.1.35）
独立 patch sprinkling 在 24/24 cells 中确认 w_max_ratio 与局域 $H(t_c)$ 反单调相关。密度匹配 sub-sampling 确认方向（0/192 early > late）。12/30 beyond-density（$d \geq 3$）。§4.1.34 Phase B 重新解读为正确的反单调信号。**Wall 机制在局域层面运作得到直接验证。**

---

## Synthesis: The Evidence Architecture

## 综合：证据架构

### The Three-Layer Structure of Conjecture E

```
┌─────────────────────────────────────────────────────────────┐
│                    CONJECTURE E                              │
│   "A finite causal set encodes spacetime curvature"          │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  Layer 1: E-wall (ESTABLISHED, 90–95%)                       │
│  ┌─────────────────────────────────────────────┐             │
│  │ σ((R-Rc)/w) → Θ(S_max - S_EH)              │             │
│  │ = curvature upper bound / admissibility gate │             │
│  │ Evidence: §4.1.21 (anti-correlation),        │             │
│  │   §4.1.22–26b (density dominance theorem)    │             │
│  └─────────────────────────────────────────────┘             │
│                                                              │
│  Layer 2: E-bulk-first-order (STRONG SUPPORT, 91–95%)        │
│  ┌─────────────────────────────────────────────┐             │
│  │ Post-density observables track H (d=4)       │             │
│  │ Two channels: spectral (B_ℓ) + transverse    │             │
│  │   (antichain) → single DoF convergence        │             │
│  │ Evidence: §4.1.27 (spectral 6/18),           │             │
│  │   §4.1.28 (antichain 21/21),                 │             │
│  │   §4.1.31 (unification |ρ|→0.86),            │             │
│  │   §4.1.32 (target = H, α≈1),                 │             │
│  │   §4.1.34 (global FRW 37/126),               │             │
│  │   §4.1.35 (local H(t) 24/24 anti-monotone)   │             │
│  └─────────────────────────────────────────────┘             │
│                                                              │
│  Layer 3: E-bulk-second-order (CHARACTERIZED, 75–80%)        │
│  ┌─────────────────────────────────────────────┐             │
│  │ H → R = d(d-1)H² bridge                      │             │
│  │ = continuum limit / BDG theorem               │             │
│  │ Evidence: §4.1.33 (0/9 at d=4; diagnosed     │             │
│  │   as continuum-limit construction)            │             │
│  └─────────────────────────────────────────────┘             │
│                                                              │
└─────────────────────────────────────────────────────────────┘
```

### Information Content Mapping

The complete picture of what a causal set "knows" about curvature, organized by information layer:

| Layer | Content | Independent DoF | Curvature response | Encoding |
|-------|---------|----------------|-------------------|----------|
| **Zeroth-order** | Total density $D = \Sigma C_k$ | 1 | $\rho \approx -0.98$ (strong anti-correlation) | Sigmoid wall (E-wall) |
| **First-order** | Differential decay $C_0/C_1$ | 0 (internal partition of $D$) | Partial positive via bdg_d2c (6/9) | Density-derived, not independent |
| **Second-order+** | Density-removed residuals within $\{C_k\}$ | 0 | $|\rho| < 0.24$ (noise) | None (DDT closed) |
| **Longitudinal path** | $k$-chain counts $N_k$ | 0 | $d=2,3$: 0; $d=4$: finite-size | Density-absorbed (§4.1.30) |
| **Spectral** | $B_\ell$ eigenvalues | **≥1** | $|\rho_\text{resid}| \approx 0.50$–$0.70$ ($d=4$) | **Confirmed (§4.1.27)** |
| **Transverse** | Antichain/layer structure | **≥1** | $|\rho_\text{resid}| \approx 0.56$–$0.82$ (all $d$) | **Confirmed (§4.1.28) — strongest** |
| **Spectral × Transverse** | Cross-channel residual | **1** ($d \geq 3$) | $|\rho_\text{cross}| \to 0.86$, monotonically strengthening with $N$ | **Single DoF convergence (§4.1.31)** |
| **Geometric target** | Shared DoF's continuous identity | — | $d=4$: $\alpha \approx 1$ (both channels) | **= expansion rate $H$ (§4.1.32)** |
| **Non-uniform** | Schwarzschild background | ? | $|\rho_\text{resid}| \approx 0.33$–$0.37$ ($N=512$ only) | Weak (§4.1.29), needs higher $d$ |

### 信息内容映射
因果集对曲率的"知识"完整图谱：零阶 = 总密度（wall），纵向 = 密度吸收（DDT 关闭），谱层 + 横向层 = beyond-density（两条通道收敛到单一 DoF = $H$）。

---

## The DDT Escape Channel Summary

| Path | DDT Condition | Experiment | Beyond density | Status |
|------|--------------|-----------|---------------|--------|
| I.1 $B_\ell$ spectral | C1 | §4.1.27 | 6/18 ($d=4$: 4/6) | ✅ Confirmed |
| I.2 Antichain | C1 | §4.1.28 | **21/21** (all $d$) | ✅ **Strongest** |
| I.3 Schwarzschild | C2 | §4.1.29 | 2/27 ($N=512$ only) | ⚠️ Weak |
| I.4 Chain $N_k$ | C1 | §4.1.30 | 0/72 ($d=2,3$); finite-size at $d=4$ | ❌ Unreliable |

---

## Phase V: EH Bridge — Systematic Construction Testing (§4.1.36–40)

## 第五阶段：EH 桥接——系统构造测试（§4.1.36–40）

Phase III identified the geometric target as the expansion rate $H$ (first-order) and characterized the second-order bridge ($H \to R \sim H^2$) as a continuum-limit construction. Phase V systematically tests **five candidate construction paths** on 2,088 de Sitter realizations ($d=4$, $N=128$–$1024$, $H=0$–$2$) to determine whether the bridge can be realized at finite $N$.

第三阶段识别了几何靶标为膨胀率 $H$（一阶），并将二阶桥接（$H \to R \sim H^2$）定性为连续极限构造。第五阶段在 2,088 次 de Sitter realization（$d=4$, $N=128$–$1024$, $H=0$–$2$）上系统测试**五条候选构造路径**，确定桥接能否在有限 $N$ 实现。

### T4: Spectral Ratio Bridge (§4.1.36) ❌

Eight eigenvalue ratio/gap features tested against $\alpha \approx 2$ target. All saturate at grid ceiling ($\alpha = 8.00$), with weak residual signal and step-function artifacts from only 4 non-zero $H$ levels. `b1_std_over_mean` shows $|\rho_\text{resid}| > 0.3$ at $d=4$ but **sign flips** between $N=256$ (+0.495) and $N=512$ (−0.677). **Conclusion**: eigenvalue algebraic operations cannot produce an $R$-tracking shortcut.

### T5: $N$-Scaling of $\alpha_\text{eff}$ (§4.1.37) ✅

320 realizations. Three independent methods all confirm $\alpha_\text{eff}(N)$ **monotonically increasing** ($\rho = +1.00$): (A) $R^2$ grid scan, (B) log-log slope, (C) $R^2(H^2)/R^2(H)$ ratio. $H^2$ explains 2.2–2.9× more variance than $H$ at per-realization level. Extrapolation convergence rate $\gamma \approx 0.17$–$0.22$ (slow; $N \sim 10^9$ to reach $\alpha = 1.9$). **Significance**: the second-order signal is real and growing, not a finite-size artifact.

### T1: Density-Assisted BDG Calibration (§4.1.38) ❌

320 realizations, 6 strategies. Raw `b1_mean` is DDT-trapped ($|\rho| < 0.11$). Density correction $\rho^{2/d}$ **overcorrects** rather than helps ($R^2$: $0.904 \to 0.719$). OLS residualization destroys signal (b1_std and density are nonlinearly entangled). **b1_mean route permanently closed**; raw `b1_std` confirmed as **unique robust spectral bulk observable**.

### T3: Two-Step Analytical Squaring (§4.1.39) ✅ **POSITIVE**

320 realizations, 5 strategies. **Core discovery**: $H^2$ information is **extraction-sensitive**, not absent.

| Strategy | $R^2(\hat{R}, H^2)$ @ $N\!=\!1024$ | $\alpha_\text{eff}$ |
|----------|--------------------------------------|---------------------|
| S1 Global OLS (b1_std) | **0.987** | 1.50 |
| S3 Rank-preserving | **0.9996** | **2.00** at ALL $N$ |
| S5 Antichain (w_max_ratio) | 0.993 ($N\!=\!512$) | 1.75 |

All 5 strategies × all $N$ × both channels: $\Delta R^2 > 0$ (universal improvement). **Correction to §4.1.33**: the original "bridge = continuum limit" was overly conservative — $H^2$ information is already encoded in the finite-$N$ rank structure; what was missing was the extraction method, not the information itself. Correct convergence target: $\alpha_\text{raw} \to 1$ (then square to 2), estimated $N \sim 10^4$–$10^5$ (far better than T5's $10^9$).

### T2: Cross-Scale Variance (§4.1.40) ❌

768 realizations, three methods all negative. (A) Inter-realization variance anti-correlates $H$ (CLT-driven, $\beta < 0$). (B) Sub-block variance fails (structure destroyed). (C) Per-realization $\delta = (\text{obs} - \text{baseline})^2$ achieves $R^2(H^2) = 0.992$ but is **algebraically equivalent to T3**.

**KEY CONCEPTUAL INSIGHT**: The second-order EH bridge is **NOT an independent variance mode** but the **squared amplitude of the first-order bulk deviation from flat baseline** — obs measures $H$, $(\text{obs} - \text{baseline})^2$ measures $H^2 \sim R$. One natural chain, not two independent channels.

### T1–T5 全景结论 / Panoramic Conclusion

| Construction | § | Result | Key finding |
|-------------|---|--------|------------|
| T4 Spectral ratios | 4.1.36 | ❌ | 0/8 features reach $\alpha \approx 2$ |
| T5 $N$-scaling | 4.1.37 | ✅ | $\alpha_\text{eff}(N)$ monotonically increasing ($\rho = +1.00$) |
| T1 Density calibration | 4.1.38 | ❌ | $\rho^{2/d}$ overcorrects; b1_mean permanently closed |
| T3 Two-step squaring | 4.1.39 | ✅ | $R^2 = 0.987$ (OLS) / $0.9996$ (rank-preserving) |
| T2 Cross-scale variance | 4.1.40 | ❌ | No independent variance channel; $\delta \equiv$ T3 algebraically |

Three paths closed (T1, T2, T4), one path confirmed (T3), convergence direction verified (T5). The second-order bridge is not a separate theoretical construction but the natural consequence of squaring a well-calibrated first-order observable.

三条路径关闭（T1、T2、T4），一条路径确认（T3），收敛方向已验证（T5）。二阶桥接不是独立的理论构造，而是对校准良好的一阶可观测量取平方的自然结果。

---

## Confidence Assessment (Post-§4.1.40, T1–T5 Complete)

| Layer | Confidence | Key evidence |
|-------|-----------|-------------|
| **E-wall** | **90–95%** | Sigmoid → EH as curvature upper bound; direction, scaling, convergence all confirmed |
| **E-bulk-first-order** | **91–95%** | DDT escape (antichain 21/21, spectral 6/18); dual-channel convergence to single DoF ($|\rho| \to 0.86$); geometric target = $H$ at $d=4$; global FRW extension 37/126; **local $H(t)$ tracking 24/24 anti-monotone, 12/30 beyond-density** |
| **E-bulk-second-order** | **85–90%** | T3 two-step squaring $R^2 = 0.987$ (OLS) / $0.9996$ (rank-preserving); T5 confirms $\alpha_\text{eff}(N)$ monotonically increasing; T1/T2/T4 negative but informative; second-order = squared amplitude of first-order deviation |
| **Overall Conjecture E** | **92–95%** | Main theory closed in de Sitter; beyond-dS dual-confirmed; T1–T5 complete — viable bridge path identified (T3), convergence direction confirmed (T5) |

---

## Remaining Open Questions

1. ~~**Second-order EH bridge**~~: **Resolved by T1–T5 (§4.1.36–40)**. T3 two-step squaring demonstrates the bridge is constructible at finite $N$ ($R^2 = 0.987$ OLS / $0.9996$ rank-preserving). The remaining task is convergence: $\alpha_\text{raw} \to 1$ at $N \sim 10^4$–$10^5$.

2. ~~**$d=4$ local signal N-scaling**~~: **Resolved by §4.1.41**. Extending §4.1.35 to $N=512$ confirms: `layer_width_std` $|\rho_\text{resid}|$: $0.527 \to 0.571 \to 0.604$ (Spearman $\rho(N, |\rho|) = +1.0$); total beyond-density 12/30 → **19/42**; d=2 signal emerges at $N=512$. Anti-monotone direction 0/480 early > late.

3. ~~**DDT Condition C2**~~: **Partially answered by §4.1.42**. 3+1D Schwarzschild sprinkling (150 realizations) yields 4/32 beyond density — improvement over 1+1D (2/27) but far weaker than de Sitter (21/21). The antichain channel is density-dominated in Ricci-flat backgrounds; beyond-density signals come from interval distribution shape ($C_1$, $C_2$). **Physical interpretation**: the antichain channel primarily responds to scalar curvature $R$ (de Sitter expansion widens spatial slices), not Weyl curvature. DDT C2 is partially confirmed: non-uniform backgrounds do allow weak escape, but the strongest channel (antichains) requires $R \neq 0$.

4. **English paper**: The only remaining fully open problem from the overall theory system.

### 剩余开放问题
1. ~~**二阶 EH 桥接**~~：**已由 T1–T5（§4.1.36–40）解决**。T3 两步平方证明桥接在有限 $N$ 可构造（$R^2 = 0.987$ OLS / $0.9996$ 秩保持）。剩余任务是收敛：$\alpha_\text{raw} \to 1$，估计 $N \sim 10^4$–$10^5$。
2. ~~**$d=4$ 局域信号 $N$-scaling**~~：**已由 §4.1.41 解决**。将 §4.1.35 扩展至 $N=512$：`layer_width_std` $|\rho_\text{resid}|$: $0.527 \to 0.571 \to 0.604$（$\rho(N, |\rho|) = +1.0$）；beyond-density 总计 12/30 → **19/42**；$d=2$ 新信号涌现。反单调方向 0/480。
3. ~~**DDT 条件 C2**~~：**§4.1.42 部分回答**。3+1D Schwarzschild（150 次 realization）：4/32 超越密度（vs 1+1D 2/27），但远弱于 de Sitter（21/21）。反链通道在 Ricci 平坦背景下被密度吸收；超密度信号来自 $C_1$/$C_2$ 区间分布形状。反链通道主要响应标量曲率 $R$，非 Weyl 曲率。DDT C2 部分确认：非均匀背景允许弱逃逸，但最强通道需 $R \neq 0$。
4. **英文论文**：整个理论体系唯一完全开放的问题。

---

## Experiment Index

| § | Experiment | Realizations | Key result |
|---|-----------|-------------|-----------|
| 4.1.21 | Sigmoid → EH correspondence | 375 | $R$ anti-correlates $H$; wall = curvature upper bound |
| 4.1.22 | BD ratio scan | 360 | All density metrics $|\rho| \approx 0.98$ anti-correlate $H$ |
| 4.1.23 | Curvature-sensitive BD | 360 | Two camps; no positive density metric |
| 4.1.24 | Bridge candidate compare | 360 | KL 7/9 strongest; bdg_d2c 6/9 |
| 4.1.25 | PCA decomposition | 360 | PC1=88.6% density; PC2+=11.4% |
| 4.1.26 | PC-Space proxy design | 360 | Two-line architecture |
| 4.1.26b | Density Degeneracy Theorem | 360 | 0 independent shape DoF within $\{C_k\}$ |
| 4.1.27 | Sorkin d'Alembertian $B_\ell$ | 360 | Spectral escape: 6/18 beyond density |
| 4.1.28 | Antichain structure | 360 | Transverse escape: **21/21** beyond density |
| 4.1.29 | Schwarzschild sprinkling | 180 | Weak C2: 2/27 at $N=512$ |
| 4.1.30 | Chain statistics $N_k$ | 432 | Failed: density-absorbed |
| 4.1.31 | Dual-channel unification | 360 | Single DoF: $|\rho_\text{cross}| \to 0.86$ |
| 4.1.32 | Geometric target identification | 360 | Target = $H$ at $d=4$ ($\alpha \approx 1$) |
| 4.1.33 | EH bridge candidates | 360 | 0/9 at $d=4$; bridge = continuum limit |
| 4.1.34 | Beyond de Sitter (global) | 432 | Phase A: 37/126 beyond density |
| 4.1.35 | Local $H(t)$ tracking | 768+384 | 24/24 anti-monotone; 12/30 beyond density |
| 4.1.36 | T4: Spectral ratio bridge | 360 | 0/8 features reach $\alpha \approx 2$; eigenvalue algebra path closed |
| 4.1.37 | T5: $N$-scaling of $\alpha_\text{eff}$ | 320 | $\alpha_\text{eff}(N)$ monotonically increasing ($\rho = +1.00$); $\gamma \approx 0.2$ |
| 4.1.38 | T1: Density-assisted BDG | 320 | 0/6 strategies improve on raw b1_std; density correction overcorrects |
| 4.1.39 | T3: Two-step squaring | 320 | $R^2 = 0.987$ (OLS) / $0.9996$ (rank-preserving); **POSITIVE** |
| 4.1.40 | T2: Cross-scale variance | 768 | No independent variance channel; $\delta \equiv$ T3 algebraically |
| 4.1.41 | Local $H(t)$ N-scaling | 1152+576 | $d=4$ beyond-density strengthens with $N$; 19/42 total; 0/480 anti-monotone |
| 4.1.42 | 3+1D Schwarzschild | 150 | 4/32 beyond density; antichain density-dominated in Ricci-flat; $C_1$/$C_2$ weak escape |

**Total**: ~10,066+ realizations across 22 experiments.

---

*Document generated: 2026-03-24; updated 2026-03-25 (Phase V: T1–T5; §4.1.41: N-scaling; §4.1.42: 3+1D Schwarzschild)*
*Status: Conjecture E — main theory closed, beyond-dS triple-confirmed, EH bridge T1–T5 complete, DDT C2 partially answered (Weyl weak / antichain requires R≠0), confidence 92–95%*
