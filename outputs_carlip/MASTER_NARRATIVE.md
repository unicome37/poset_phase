# 总叙事页：因果集结构选择的两层筛选理论

**Master Narrative — 供摘要、导论、结论统一调用**

---

## 一句话核心

> A two-layer screening architecture—linear admissibility via Benincasa-Dowker action plus quadratic identity via Mahalanobis reference-manifold distance—uniquely selects 4D Lorentzian causal sets from a 25-family poset library with zero free parameters, with the identity basin deepening as $N^{-1.57}$ and the reference manifold covariance collapsing as $N^{-3.38}$.

中文：

> 一套两层筛选架构——以 Benincasa-Dowker 作用量实现的线性准入 + 以 Mahalanobis 参考流形距离实现的二次认同——在25族偏序集库中以零自由参数唯一选择4D洛伦兹因果集，认同域深度以 $N^{-1.57}$ 加深，参考流形协方差以 $N^{-3.38}$ 塌缩。

---

## 四段叙事（分别对应摘要的四句话）

### ① 问题

在因果集量子引力的路径积分中，Kleitman-Rothschild 定理表明绝大多数有限偏序集是非几何的三层结构。**为什么4维洛伦兹时空被选择？** 这是因果集纲领必须回答的核心问题。

> In the causal set path integral, the Kleitman-Rothschild theorem implies that almost all finite posets are non-geometric three-layered structures. What mechanism selects 4D Lorentzian spacetime?

### ② 方法

我们构建了包含25个偏序集族（4个洛伦兹维度 + 3个 KR 型 + 多种随机分层/渗流/区间序结构 + 8个对抗性族）的完整样本空间，并在 $N=16$–$1024$ 范围内系统比较三个结构特征的判别能力：有效维度 $d_{\mathrm{eff}}$（Myrheim-Meyer）、区间比率 $C_1/C_0$、最大反链宽度比 $w/N$。

> We construct a 25-family poset library spanning Lorentzian causal sets (2D–5D), KR-type generic posets, random layered structures, and 8 adversarial families, and systematically compare the discriminative power of three structural features—effective dimension, interval fraction, and maximum antichain width—across $N=16$ to $1024$.

### ③ 发现

一个两层筛选结构从数值中涌现。**第一层**：$S_{\mathrm{BD}} = \mathbf{c}^\top \Delta\mathbf{C}$（Benincasa-Dowker action）是 interval count 空间中的线性超平面，淘汰曲率不对的结构，但不能唯一选择 Lor4D（单独排名仅 14/17）。**第二层**：$S_{\mathrm{MD}} = \delta^\top \Sigma^{-1} \delta$（Mahalanobis 距离）是特征空间中的二次椭球，以 Lor4D 参考流形 $\boldsymbol{\mu}(N)$ 为中心，零自由参数，在所有 $N$ 和全部25族中唯一将 Lor4D 排名第1。两层的交集恰好且仅仅是 Lor4D。

> A two-layer screening structure emerges. **Layer 1**: The Benincasa-Dowker action $S_{\mathrm{BD}} = \mathbf{c}^\top \Delta\mathbf{C}$, a linear functional of interval counts, filters out posets with incorrect curvature but cannot uniquely select Lor4D (rank 14/17 alone at $N=128$). **Layer 2**: The Mahalanobis distance $S_{\mathrm{MD}} = \delta^\top \Sigma^{-1} \delta$ from the Lor4D reference manifold $\boldsymbol{\mu}(N)$, a zero-parameter quadratic functional, uniquely ranks Lor4D first across all $N$ and all 25 families. Their intersection is Lor4D and only Lor4D.

### ④ 意义

随着 $N$ 增大，身份域持续深化（Mahalanobis gap: $1.38 \to 1.93 \times 10^8$；Fisher 信息 $\propto N$），参考流形协方差塔缩为 $\det(\Sigma) \propto N^{-3.38}$。两层之间存在**升阶关联**（梯度幅度级 $|\cos| \approx 0.85$），但不是逐点梯度同一——这是同一因果几何在不同阶数上的自然编码，从一阶线性准入升阶为二阶二次约束。

> As $N$ grows, the identity basin deepens (Mahalanobis gap: $1.38 \to 1.93 \times 10^8$; Fisher information $\propto N$) and the reference manifold covariance collapses as $\det(\Sigma) \propto N^{-3.38}$. The two layers share an **order-raising correlation** (gradient magnitude $|\cos| \approx 0.85$), not a pointwise gradient identity—they encode the same causal geometry at different orders, from first-order linear gating to second-order quadratic confinement.

---

## 关键数据表（供正文引用）

| 量 | 值 / 行为 | 来源实验 |
|:--:|:--------:|:--------:|
| 样本空间 | 25 族（17标准 + 8对抗） | expanded_family_robustness |
| $N$ 范围 | 16 – 1024 | large_n_extreme_test |
| $S_{\mathrm{BD}}$ 单独排名 | 14/17 (N=128) | smd_sbd_connection |
| $S_{\mathrm{MD}}$ 排名 | 全 N 全族 #1 | mahalanobis_lsd_test |
| 两层交集 | $= \{\text{Lor4D}\}$ | smd_sbd_connection |
| Mahal gap | $1.38$ (N=16, REPS≥80) $\to 94.1$ (N=128) | basin_deepening / n16_stability |
| Margin @ N=1024 | $1.93 \times 10^8$ | large_n_extreme_test |
| $V_{\mathrm{eff}}$ scaling | $\propto N^{-1.57}$ | basin_deepening_experiment |
| Fisher scaling | $\propto N^{1.00}$ | basin_deepening_experiment |
| $\det(\Sigma)$ scaling | $\propto N^{-3.38}$ | mu_trajectory_theory |
| $\sigma^2$ scaling | $\propto N^{-1}$ (CLT) | mu_trajectory_theory |
| $\mu_d(\infty)$ | 3.967 ≈ 4 | cstar_wstar_first_principles |
| $\mu_c(\infty)$ | 0.345 | cstar_wstar_first_principles |
| $\mu_w(\infty)$ | 0.227 | gradient_alignment_v2 |
| $w(N)$ exponent | 0.284 ≈ 1/d | cstar_wstar_first_principles |
| Mixing exponent η | 0.74 ± 0.08 | min_distortion_verify |
| CV #1 rate | LSD 98%, Mahal 94% | prediction_b_cross_validation |
| Seed 再现 | LSD 80/80, Mahal 79/80 | prediction_b_seed_reproducibility |
| 梯度幅度 | $|\cos| \approx 0.85$ | gradient_alignment_v2 |
| 梯度符号稳定性 | 5/11 正, 6/11 负 → 不可靠 | gradient_alignment_v2 |

---

## 论文结构对照

| 论文章节 | 从本页调用 | 补充展开 |
|:--------:|:----------:|:--------:|
| **Abstract** | 四句话叙事 ①②③④ | — |
| **Introduction** | ① 问题 + 文献综述 | Dhar, Prömel, Carlip, Surya |
| **Methods** | ② 方法 | 生成器详述, 三特征定义 |
| **Results** | ③ 发现 + 数据表 | 分 N 详细表格, 图 |
| **Discussion §5.1** | 两层架构 | DISCUSSION §6.1 |
| **Discussion §5.2** | Basin deepening | DISCUSSION §6.2 |
| **Discussion §5.3** | 梯度桥（辅助） | DISCUSSION §6.3 |
| **Discussion §5.4** | Reference manifold | DISCUSSION §6.4 |
| **Conclusion** | ④ 意义 + 一句话核心 | — |

---

## 禁止在正文中出现的说法

| 旧说法 | 风险 | 替代 |
|:------:|:----:|:----:|
| "cos=0.97 proves first-principles bridge" | Jacobian 病态 | "$|\cos| \approx 0.85$ suggests structural relatedness" |
| "gradient alignment confirms derivation" | 符号不稳 | "magnitude-level correlation supports order-raising interpretation" |
| "logH = physical entropy" | Carlip C1 | "residual ordering freedom proxy" |
| "LSD-Well optimally selects Lor4D" | 有参数 | "Mahalanobis LSD achieves zero-parameter selection" |
| “all families fail except Lor4D” | N=16 参考集精度敏感 | “Lor4D is uniquely ranked #1 at all $N$; at $N=16$ the margin is small and requires $\geq 80$ reference sprinklings” |

---

*本页为统一调用接口。详细推导见各专题文档；英文段落见 DISCUSSION_THEORY_IMPLICATIONS.md §6。*
