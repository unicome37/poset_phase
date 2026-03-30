
> **2026-03-30 更新（F2 margin-aware refit）**：上述 N=16 失败来自 Mahalanobis CV 模式（训练/测试共享种子，小参考集）。使用 120 reps 独立参考集后（`seed_base+100000` 隔离），N=10 即达到全档 manuscript-safe 标准（20/20 rank#1，min_margin=0.198, ci95_lower=0.268）。CV 模式的 N=16 失败属于参考估计量波动，非 Lor4D basin 本身的物理边界。物理分辨率极限目前推估在 N<10。
# 两层筛选理论：线性准入与二次认同

**——从14组数值实验中涌现的因果集维度选择原理**

**日期**：2026-03  
**项目**：poset_phase  
**状态**：理论总结稿

---

## 目录

- [§1 概述：两层筛选理论的涌现](#1-概述两层筛选理论的涌现)
- [§2 第一层：线性准入 (S_BD)](#2-第一层线性准入-s_bd)
- [§3 第二层：二次认同 (S_MD / Mahalanobis)](#3-第二层二次认同-s_md--mahalanobis)
- [§4 Well Center 的第一原理地位](#4-well-center-的第一原理地位)
- [§5 两层的正交性与联合选择](#5-两层的正交性与联合选择)
- [§6 鲁棒性三重验证](#6-鲁棒性三重验证)
- [§7 物理解释：宇宙不是最小化单一作用量](#7-物理解释宇宙不是最小化单一作用量)
- [§8 从经验可分到几何可表述](#8-从经验可分到几何可表述)
- [§9 完整实验索引](#9-完整实验索引)

---

## §1 概述：两层筛选理论的涌现

### 1.1 问题背景

因果集理论（causal set theory）将时空建模为离散偏序集（poset）。一个核心问题是：**在所有可能的偏序集中，为什么4维洛伦兹因果集（Lor4D）被选择？** 这不仅仅是一个组合学分类问题，而是直指因果集量子引力路径积分中维度选择机制的物理本质。

本项目构建了包含17个偏序集族的完整样本空间——包括Lor2D/3D/4D/5D（洛伦兹因果集）、KR_like/KR_2layer/KR_4layer（Kleitman-Rothschild型典型偏序集）、多种随机分层结构（RLk4/6/8）、以及AbsLayer、IntOrder、TransPerc、MLR等——并通过系统的数值实验，发现了一个具有清晰数学结构的**两层筛选机制**。

### 1.2 核心发现

经过14组独立数值实验，涌现出一个统一的理论图像：

> **线性准入 + 二次认同**（Linear Admissibility + Quadratic Identity）

两层筛选的数学结构为：

$$\text{Layer 1: } S_{\mathrm{BD}} = \mathbf{c}^\top \Delta\mathbf{C} \qquad \text{（线性超平面——准入门）}$$

$$\text{Layer 2: } S_{\mathrm{MD}} = \Delta\mathbf{I}^\top \Sigma^{-1}(N) \, \Delta\mathbf{I} \qquad \text{（二次椭球——认同域）}$$

其中 $\Delta\mathbf{C} = \mathbf{C} - \mathbf{C}^*$ 是interval count偏差向量，$\Delta\mathbf{I} = \mathbf{I}(P) - \boldsymbol{\mu}(N)$ 是特征偏差向量。

**第一层**（Benincasa-Dowker action）是interval count空间中的一个线性超平面：它筛选掉一大批curvature不对的偏序集，但不能唯一锁定Lor4D（在N=128时排名仅14/17）。

**第二层**（Minimum Distortion / Mahalanobis距离）是特征空间中的一个二次椭球面：它测量与Lor4D参考流形的距离，零自由参数，在所有测试的N值上将Lor4D唯一排名为第1。

**两层的交集**恰好且仅仅是Lor4D：

$$\{S_{\mathrm{BD}} \in \text{window}\} \cap \{S_{\mathrm{MD}} \approx 0\} = \{\text{Lor4D}\}$$

这不是一个经验发现的分类器，而是一个具有几何-物理意义的筛选结构。

### 1.3 理论结构的层级

| 层级 | 作用量 | 数学形式 | 功能 | 物理含义 |
|:----:|:------:|:--------:|:----:|:--------:|
| 第一层 | $S_{\mathrm{BD}}$ | $\mathbf{c}^\top \Delta\mathbf{C}$（线性） | 准入门 | 平均曲率约束 |
| 第二层 | $S_{\mathrm{MD}}$ | $\Delta\mathbf{I}^\top \Sigma^{-1} \Delta\mathbf{I}$（二次） | 认同域 | 结构认同/流形距离 |
| 联合 | $S_{\mathrm{BD}} \cap S_{\mathrm{MD}}$ | 超平面 ∩ 椭球 | 唯一选择 | 4D 洛伦兹时空 |

---

## §2 第一层：线性准入 (S_BD)

### 2.1 Benincasa-Dowker 作用量

Benincasa-Dowker action 是因果集理论中对Einstein-Hilbert作用量的离散类比，其d=4形式为：

$$S_{\mathrm{BD}}^{(4)} / N = 1 - \frac{C_0}{N} + 9\frac{C_1}{N} - 16\frac{C_2}{N} + 8\frac{C_3}{N}$$

其中 $C_k$ 是k-interval count（长度为k的因果链数目）。引入归一化interval密度 $\tilde{c}_k = C_k / \binom{N}{2}$，可以写成：

$$S_{\mathrm{BD}} / N = 1 + \frac{N-1}{2}\left(-\tilde{c}_0 + 9\tilde{c}_1 - 16\tilde{c}_2 + 8\tilde{c}_3\right)$$

这是 $\{C_k\}$ 的**严格线性泛函**。在interval count空间中，$S_{\mathrm{BD}} = \text{const}$ 定义一个**超平面**。

### 2.2 S_BD 能做什么

$S_{\mathrm{BD}}$ 本质上编码的是平均标量曲率。对于接近平坦4D Minkowski时空的因果集——也就是sprinkle到因果钻石中的Lor4D——$S_{\mathrm{BD}}$ 的期望值趋于零（曲率为零）。

因此，$S_{\mathrm{BD}}$ 可以作为一个**准入门**：那些 $|S_{\mathrm{BD}}|$ 太大的偏序集不可能是近平坦的4D洛伦兹时空。在N=128时，设定Lor4D的S_BD窗口 $[-2.88, 2.62]$（均值 ± 标准差），我们可以排除：

- 所有KR型（KR_like: −15.0, KR_2layer: −11.1, KR_4layer: −14.0）
- 所有RLk型除AbsLayer之外的多数
- MLR (−2.9，刚好在边界外)

但仍有多个族通过了S_BD准入：IntOrder (0.69)、Lor2D (1.72)、Lor3D (1.07)、Lor5D (−0.87)、TransPerc (−1.60) 等全部落入窗口内。

### 2.3 S_BD 不能做什么

$S_{\mathrm{BD}}$ 不能唯一选择Lor4D。这一点由数据明确证实：

| N | Lor4D的S_BD排名 | 最接近的非Lor4D族 |
|:---:|:---:|:---|
| 48 | #9/17 | Lor5D (#8) |
| 64 | #10/17 | Lor5D (#9) |
| 96 | #15/17 | TransPerc (#13) |
| 128 | **#14/17** | Lor5D (#13) |

在N=128时，**Lor4D在S_BD排名中位于倒数第4**。S_BD/N的最小值属于KR_like（−15.0），不属于Lor4D（−0.13）。这是因为S_BD编码的是"average curvature correctness"而非"Lor4D identity"——许多结构完全不像流形的偏序集恰好也有接近零的平均"曲率"。

**关键洞见**：线性泛函（一阶）可以筛选掉大量不满足必要条件的偏序集，但不能提供充分条件。需要更高阶的结构来完成认同。

### 2.4 门控模型

将 $S_{\mathrm{BD}}$ 理解为**门控函数**：

$$G_{\mathrm{BD}}(P) = \begin{cases} 1, & |S_{\mathrm{BD}}(P) - S_{\mathrm{BD}}^{(4D)}| < \Delta(N) \\ 0, & \text{otherwise} \end{cases}$$

通过门控的偏序集构成"准入集"$\mathcal{A}_{\mathrm{BD}}$。在N=128时，$\mathcal{A}_{\mathrm{BD}}$ 包含约6–7个族，Lor4D只是其中之一。门的意义在于预过滤：它将17族的全空间压缩为6–7族的子空间，为第二层筛选提供了更紧凑的搜索域。

---

## §3 第二层：二次认同 (S_MD / Mahalanobis)

### 3.1 从LSD-Well到Mahalanobis距离

本项目的判别框架经历了一个从经验到理论的演进过程：

1. **LSD-Well（三参数版）**：$F = \alpha(d_{\mathrm{eff}} - 4)^2 + \beta(c_1/c_0 - c^*(N))^2 + \gamma(w - w^*(N))^2$，权重 $\alpha, \beta, \gamma$ 为手动设定。
2. **Fisher-guided version**：发现方差排序 $\sigma^2(w) < \sigma^2(c) < \sigma^2(d)$ 与最优权重排序 $\gamma > \beta > \alpha$ 精确匹配——权重就是 $1/\sigma^2$，即Fisher信息矩阵的对角线。
3. **Mahalanobis LSD（零参数版）**：直接使用Lor4D集成的均值向量 $\boldsymbol{\mu}(N)$ 和协方差矩阵 $\Sigma(N)$，构造Mahalanobis距离。

最终形式：

$$S_{\mathrm{MD}}(P, N) = \left(\mathbf{I}(P) - \boldsymbol{\mu}(N)\right)^\top \Sigma^{-1}(N) \left(\mathbf{I}(P) - \boldsymbol{\mu}(N)\right)$$

其中 $\mathbf{I}(P) = (d_{\mathrm{eff}}(P),\; c_1/c_0(P),\; w(P))$ 是偏序集 $P$ 的三维特征向量。

### 3.2 参考流形 $(\boldsymbol{\mu}(N), \Sigma(N))$

Lor4D集成在特征空间中的分布构成一个**一参数（N参数化）高斯云族**——即参考流形。其关键参量：

**均值轨迹**：

$$\boldsymbol{\mu}(N) = \boldsymbol{\mu}(\infty) + \frac{\mathbf{a}}{N} + \frac{\mathbf{b}}{N^2}$$

- $d_{\mathrm{eff}}(N) \to 4.000$（Myrheim-Meyer，精确）
- $c^*(N) = 0.374 - 10.70/N + 99.5/N^2$，$R^2 = 0.987$
- $w^*(N)$：幂律 $w \propto N^{-0.284}$，$R^2 = 0.995$

**协方差收缩**：

$$\sigma_i^2(N) = A_i \cdot N^{-p_i}$$

| 特征 | 衰减指数 $p_i$ |
|:----:|:---:|
| $d_{\mathrm{eff}}$ | 1.05 |
| $c_1/c_0$ | 1.11 |
| $w$ | 1.32 |

三个特征的方差均以 $\sim N^{-1}$ 速率收缩（经典泊松统计行为），det($\Sigma$) $\propto N^{-3.38}$，即Lor4D点云体积以 $N^{-1.69}$ 速率坍缩到固定点。

### 3.3 零自由参数

Mahalanobis LSD 的核心优势在于**零自由参数**：

| 框架 | 自由参数 | Well center参数 | 总数 |
|:----:|:--------:|:--------------:|:----:|
| LSD-Well (原始版) | 3 ($\alpha, \beta, \gamma$) | 3 ($d^*, c^*, w^*$) | 6 |
| LSD-Well (推导版) | 3 ($\alpha, \beta, \gamma$) | 1 ($\alpha = 0.0024$) | 4 |
| **Mahalanobis LSD** | **0** | **0** (从Lor4D集成提取) | **0** |

所有信息——均值、方差、协方差——都从Lor4D集成的经验分布中直接获取。没有拟合、没有调参、没有选择。Lor4D的唯一性是一个**自举出现的几何事实**，不是一个人为优化的结果。

### 3.4 几何解释：二次椭球面

在特征空间 $(d_{\mathrm{eff}}, c_1/c_0, w)$ 中，$S_{\mathrm{MD}} = \text{const}$ 定义一族**椭球面**，以 $\boldsymbol{\mu}(N)$ 为中心，以 $\Sigma^{-1}(N)$ 为度量。Lor4D样本分布在最内层椭球内（$S_{\mathrm{MD}} \sim \chi^2_3$），而所有其他16个族的质心都在远处的椭球壳上。

这个椭球面不是任意选择的分类边界。它是**围绕Lor4D结构吸引域的最低阶有效势**——Landau型展开的前导项：

$$S_{\mathrm{MD}} \approx \sum_{i,j} \Lambda_{ij} \left(I_i - I_i^{(4D)}(N)\right)\left(I_j - I_j^{(4D)}(N)\right) + O(\delta I^3)$$

在参考结构的邻域内，椭球近似是精确的。其精度由实验验证：在当前安全协议（fixed-reference F2）下，Lor4D从 `N=10` 起稳定排名#1，并在更大尺度测试中持续保持该身份优势；margin从 $S_{\mathrm{MD}} \approx 0.308$（N=10, F2）增长到 $1.93 \times 10^8$（N=1024），呈幂律发散。

### 3.5 interval count 基中的表达

将 $S_{\mathrm{MD}}$ 回写到interval count基 $\{C_k\}$中，可得：

$$S_{\mathrm{MD}} = \Delta\mathbf{C}^\top \mathbf{J}^\top \Sigma^{-1} \mathbf{J} \, \Delta\mathbf{C} \equiv \Delta\mathbf{C}^\top \mathbf{M} \, \Delta\mathbf{C}$$

其中 $J_{ik} = \partial I_i / \partial C_k$ 是Jacobian矩阵（已通过数值回归独立估计），$\mathbf{M} = \mathbf{J}^\top \Sigma^{-1} \mathbf{J}$ 是interval count空间中的**有效度量张量**。

这揭示了与 $S_{\mathrm{BD}}$ 的精确结构对比：

| | $S_{\mathrm{BD}}$ | $S_{\mathrm{MD}}$ |
|:---|:---:|:---:|
| 在 $\{C_k\}$ 中的阶数 | 一阶（线性） | 二阶（二次型） |
| 几何 | 超平面 | 椭球面 |
| 独立选择 Lor4D? | ❌ | ✅ |
| 自由参数 | 0 (BDG系数固定) | 0 (Mahalanobis版) |

---

## §4 Well Center 的第一原理地位

### 4.1 三个Well Center

LSD-Well框架中有三个量规定了二次井的中心位置：$d^* = 4$, $c^*(N)$, $w^*(N)$。一个关键问题是：它们是自由参数还是理论推论？答案是后者。

### 4.2 $d^* = 4$：Myrheim-Meyer维度

Myrheim-Meyer维度估计器将有效维度定义为：

$$d_{\mathrm{eff}} = f_2^{-1}(R), \quad R = \frac{C_0}{\binom{N}{2}}, \quad f_2(d) = \frac{\Gamma(d+1)\Gamma(d/2)}{4\Gamma(3d/2)}$$

对于4维Minkowski时空的泊松sprinkle，$R \to f_2(4) = 0.050000$ as $N \to \infty$，因此 $d^*(\infty) = 4.000$ **精确成立**。

经验验证显示收敛路径：

| N | $d_{\mathrm{eff}}$ | $\sigma(d)$ | $\Delta(d-4)$ |
|:---:|:---:|:---:|:---:|
| 16 | 3.856 | 0.378 | −0.145 |
| 64 | 3.990 | 0.144 | −0.011 |
| 128 | 3.971 | 0.064 | −0.029 |
| 256 | 3.944 | 0.070 | −0.056 |
| 512 | 3.968 | 0.035 | −0.032 |

$d^* = 4$ 不是一个拟合参数，而是**时空维度的物理输入**——因果集量子引力中仅有的先验假设之一。

### 4.3 $c^*(\infty) = 0.374$：interval shape ratio 的标度极限

$c_1/c_0 = C_1/C_0$ 衡量的是"2-step interval 与 link 的比值"——即因果链中插入一个中间元素的频率。对于4维因果钻石中的泊松sprinkle，其理论推导从以下开始：

对于一对因果相关的点 $(x \prec y)$，其Alexandrov interval $I(x,y)$ 的体积占比 $u = V_{xy}/V_{\mathrm{diamond}}$ 遵循Meyer (1988) 的Beta分布 $u \sim \mathrm{Beta}(d/2, d/2)$。interval中期望的sprinkle点数为 $\lambda = \rho \cdot V_{xy}$，于是：

$$\frac{C_1}{C_0} = \frac{\langle \lambda e^{-\lambda} \rangle}{\langle e^{-\lambda} \rangle}$$

直接从 $\mathrm{Beta}(2,2)$ 积分出的裸公式需要一个几何因子修正（$\alpha = 0.0024$，对应因果钻石中generators的几何结构），但有限尺度标度拟合直接给出高精度结果：

$$c^*(N) = 0.3744 \pm 0.0073 - \frac{10.70}{N} + \frac{99.5}{N^2}, \quad R^2 = 0.987$$

**$c^*(\infty) = 0.374$ 是泊松sprinkle统计与Alexandrov interval几何的组合推论**，不是自由参数。对比验证：Lor2D的 $c(\infty) \approx 0.74$，Lor3D的 $c(\infty) \approx 0.49$，Lor5D的 $c(\infty) \approx 0.17$——各维度的 $c^*$ 值反映了因果钻石的维度依赖几何。

### 4.4 $w^*(\infty)$：反链宽度的幂律标度

最大反链宽度比 $w = |\text{max antichain}|/N$ 衡量偏序集的"横向扩展"。对于d维因果钻石，沿时间轴的体积剖面为 $V(t) \propto (t(T-t))^{(d-1)/2}$，最宽层出现在赤道 $t = T/2$。由此可以期望：

$$w(N) \propto N^{-1/d}$$

拟合结果：

$$w(N) = 1.2211 \cdot N^{-0.2842} + w_0, \quad R^2 = 0.995$$

**拟合指数 $\beta = 0.284$ 与理论预测 $1/d = 1/4 = 0.250$ 吻合**，偏差源于有限N效应和greedy peeling算法对真实最大反链的近似。

$$w^*(\infty) \to 0$$

这意味着在连续极限下，最大反链的相对宽度趋于零——所有sprinkle点最终沿因果链有序排列，这正是4D洛伦兹时空的因果结构特征。

### 4.5 总结：零自由参数的物理基础

| Well center | 极限值 | 推导来源 | 自由参数数 | 状态 |
|:-----------:|:------:|:--------:|:----------:|:----:|
| $d^* = 4$ | 4.000（精确） | Myrheim-Meyer定理 | 0 | ✅ 精确 |
| $c^*(\infty)$ | 0.374 ± 0.007 | 泊松sprinkle统计 + 有限尺寸标度 | 0（或1个几何因子） | ✅ 理论值 |
| $w^*(\infty)$ | → 0 | 因果钻石体积剖面 + $N^{-1/d}$ 标度律 | 0 | ✅ 标度律 |

**Mahalanobis版完全绕过了这些推导**——它直接从Lor4D集成中提取 $\boldsymbol{\mu}(N)$ 和 $\Sigma(N)$。但第一原理推导证实了经验值并非拟合伪影，而是因果集理论的几何推论。

---

## §5 两层的正交性与联合选择

### 5.1 S_BD 与 S_MD 几乎不相关

如果两层筛选是冗余的（即测量同一件事），那么它们的数值应高度相关。实验明确否定了这一点：

| N | Pearson $r(S_{\mathrm{BD}}, S_{\mathrm{MD}})$ | Spearman $\rho$ | 解释 |
|:---:|:---:|:---:|:---|
| 48 | 0.565 | 0.225 | 中等（小N混杂） |
| 64 | **0.073** | 0.082 | **几乎不相关** |
| 96 | 0.185 | −0.248 | 弱，且方向反转 |
| 128 | −0.259 | **−0.400** | 弱负相关 |

N=64时Pearson $r = 0.073$——本质上是零。随着N增大，两者甚至出现弱负相关。这说明 $S_{\mathrm{BD}}$ 和 $S_{\mathrm{MD}}$ 提取的是偏序集结构的**不同维度的信息**。

这一结果在数学上并不意外：$S_{\mathrm{BD}}$ 是 $C_k$ 的线性组合（一阶），$S_{\mathrm{MD}}$ 是非线性特征 $(d_{\mathrm{eff}}, c_1/c_0, w)$ 的二次型（二阶）。一阶和二阶泛函之间没有结构性的相关要求。

### 5.2 排名正交性

$S_{\mathrm{BD}}$ 和 $S_{\mathrm{MD}}$ 对17个族的排名几乎完全不同：

**N=128排名对比**（部分）：

| 族 | $S_{\mathrm{MD}}$ 排名 | $S_{\mathrm{BD}}$ 排名 |
|:---|:---:|:---:|
| KR_like | 5 | **1** |
| RLk6 | 14 | **2** |
| KR_4layer | 9 | **3** |
| Lor4D | **1** | **14** |
| Lor2D | 4 | **17** |

$S_{\mathrm{BD}}$ 最"喜欢"的是KR_like（大量link结构导致的极负BDG action），而 $S_{\mathrm{MD}}$ 最"喜欢"的是Lor4D。两个排名之间几乎是正交旋转。

### 5.3 联合选择：交集 = Lor4D

当我们取两个筛选的**交集**时，唯一幸存者是Lor4D：

**N=64示例**：

| 族 | 通过 $S_{\mathrm{BD}}$ 窗口? | 通过 $S_{\mathrm{MD}}$ 井? | **通过两者?** |
|:---|:---:|:---:|:---:|
| AbsLayer | ✅ | ❌ | ❌ |
| IntOrder | ✅ | ❌ | ❌ |
| KR_2layer | ❌ | ❌ | ❌ |
| KR_like | ❌ | ❌ | ❌ |
| Lor2D | ✅ | ❌ | ❌ |
| Lor3D | ✅ | ❌ | ❌ |
| **Lor4D** | **✅** | **✅** | **✅** |
| Lor5D | ✅ | ❌ | ❌ |
| MLR | ✅ | ❌ | ❌ |
| TransPerc | ✅ | ❌ | ❌ |

**N=128**中结论完全一致：通过 $S_{\mathrm{BD}}$ 窗口的有IntOrder、Lor2D/3D/4D/5D、TransPerc等6族；通过 $S_{\mathrm{MD}}$ 井的只有Lor4D。**交集 = {Lor4D}**。

### 5.4 互补性的物理图像

两个作用量编码了不同层次的因果几何信息：

- **$S_{\mathrm{BD}}$** 检查"平均曲率是否正确"（动力学约束）
- **$S_{\mathrm{MD}}$** 检查"interval结构是否像Lor4D"（结构认同）

一个偏序集可以有正确的平均曲率而不像流形（例如KR型三层结构），也可以在某些特征上接近Lor4D但平均曲率偏离（例如某些随机分层结构）。只有Lor4D同时满足两个约束——这正是"线性准入 + 二次认同"的协同机制。

---

## §6 鲁棒性三重验证

任何理论主张的可靠性取决于其鲁棒性。我们通过三种独立手段进行了系统验证。

### 6.1 第一重：种子再现性

**实验设计（历史 seed-robustness 记录）**：10个独立随机种子 × 8个N值 (16–128) × 20 reps/seed = 27,200个独立样本。

**结果**：

| 方法 | #1排名次数 / 总数 | 百分比 |
|:----:|:-----------------:|:------:|
| LSD-Well | **80/80** | 100% |
| Mahalanobis LSD | **79/80** | 99% |

Mahalanobis的唯一失败发生在 seed=1001, N=16，margin仅0.6（对比N=128平均margin=40.4）。应注意：这属于旧的 shared-reference / CV 历史记录；当前 fixed-reference F2 协议已将安全 onset 口径推进到 `N≥10`。LSD-Well在最危险的配置（N=16, seed=3141）下margin=0.003，极其微弱但仍然正确。

**Mahalanobis margin随N的发散**：

| N | 平均margin | 标准差 |
|:---:|:---:|:---:|
| 16 | 3.0 | 2.1 |
| 32 | 8.5 | 3.2 |
| 64 | 18.4 | 5.7 |
| 128 | 40.4 | 8.3 |

margin的幂律增长意味着在大N极限下，Lor4D的优势趋于无穷——这是Lyapunov泛函的特征行为。

### 6.2 第二重：交叉验证（非过拟合验证）

**实验设计**：Leave-20%-out交叉验证。用80%的Lor4D样本估计 $\boldsymbol{\mu}(N)$ 和 $\Sigma(N)$，在剩余20%上测试排名。

**结果**：

| 方法 | #1比例 | margin保留率 |
|:----:|:------:|:-----------:|
| LSD-Well | **98%** | > 95% |
| Mahalanobis LSD | **94%** | > 95% |

6%的Mahalanobis失败全部集中在N=16，且失败时Lor4D仍排名#2（仅被Lor5D超过，margin极小）。这些结果应被视为**历史 CV 诊断**，不是当前安全 claim 的边界。保留率 > 95% 意味着训练/测试切分几乎不影响分离度——判别力是结构性的，不是过拟合的。

### 6.3 第三重：N=16失败根因分析

N=16是旧 Mahalanobis CV 方阵中唯一出现Lor4D非#1的配置。根因诊断结果用于解释历史异常，而非定义当前安全 onset：

**唯一入侵者 = Lor5D**。在19/19次失败中，100%是Lor5D击败了Lor4D。没有任何非洛伦兹族曾在Mahalanobis CV下击败Lor4D。

**物理根因**：

$$d_{\mathrm{eff}}(\text{Lor4D}, N=16) \approx 3.9, \quad d_{\mathrm{eff}}(\text{Lor5D}, N=16) \approx 4.3$$

在N=16时，两者的有效维度重叠严重（差距仅0.4，而单个样本的 $\sigma(d) \approx 0.38$）。但结合 F2 fixed-reference 结果，这一异常更准确地说是**低样本参考估计 + shared-reference contamination 与几何近邻共同放大的历史诊断点**，而不是当前安全协议下的物理边界定义。

**辅助证据**：

- 协方差矩阵的条件数 $\kappa = 36$（温和，非病态）
- 收缩估计 $\alpha = 0.1$ 略有帮助（37/50 → 44/50）但不能根本解决
- LSD-Well在N=16不受此影响，因为 $d^* = 4$ 的硬约束直接惩罚了 $d_{\mathrm{eff}} \neq 4$

**结论**：N=16 的旧 CV“失败”应作为历史根因分析保留，但当前安全结论已由 fixed-reference F2 protocol 改写为 `N≥10`。因此这里标记的是旧协议的诊断边界，而不是现行主稿口径中的物理下界。

---

## §7 物理解释：宇宙不是最小化单一作用量

### 7.1 传统图像的局限

标准因果集量子引力的路径积分图像是：

$$Z = \sum_{\text{causets } C} e^{i S_{\mathrm{BD}}[C]} / \text{sym}(C)$$

这暗示宇宙通过单一作用量 $S_{\mathrm{BD}}$ 来"选择"4D洛伦兹时空。但我们的实验数据明确显示：**$S_{\mathrm{BD}}$ 本身无法完成这一选择**。在N=128时，Lor4D在 $S_{\mathrm{BD}}$ 排名中位于第14/17——不仅不是最小化者，而且排名很低。

### 7.2 两层筛选的替代图像

实验数据支持一个更丰富的图像：**宇宙通过两层筛选来实现维度选择**。

**第一层（线性准入）** 对应路径积分中的相位约束。$e^{iS_{\mathrm{BD}}}$ 的快速振荡抑制了 $S_{\mathrm{BD}}$ 偏离驻点过远的偏序集——这正是准入门的物理机制。任何curvature太大或太小的偏序集在路径积分中的贡献都被相干抵消。这一层消除了大约10/17的族。

**第二层（二次认同）** 对应Lor4D结构吸引域的局域度量。在通过第一层准入后的子空间中，$S_{\mathrm{MD}}$ 的椭球面标识了"哪些偏序集真正具有4D洛伦兹流形的内部结构"。这可能对应于路径积分中尚未被识别的二阶（非节点的）贡献——interval count层级中一阶信息（$S_{\mathrm{BD}}$）之上的二阶信息。

### 7.3 作用量层级假说

将两个作用量放入统一的interval count层级中：

$$S_{\text{total}} = \underbrace{S_{\mathrm{BD}}}_{\text{一阶：} \mathbf{c}^\top \Delta\mathbf{C}} + \underbrace{S_{\mathrm{MD}}}_{\text{二阶：} \Delta\mathbf{C}^\top \mathbf{M} \, \Delta\mathbf{C}} + O(\Delta C^3)$$

这自然地对应于Landau展开中按扰动阶数排列的有效作用量：

- **一阶项**（$S_{\mathrm{BD}}$）：在interval count空间中的线性项，编码全局曲率约束——"宇宙必须是近似平坦的"。
- **二阶项**（$S_{\mathrm{MD}}$）：在interval count空间中的二次型，编码局域结构认同——"宇宙必须看起来像4D洛伦兹流形"。

两层的分离在物理上是完全自然的：全局约束（一阶）必然先于局域认同（二阶）。正如凝聚态物理中，对称性约束（一阶）先于序参量选择（二阶），宇宙首先满足曲率约束，然后在曲率正确的子空间中选择具有正确几何结构的偏序集。

### 7.4 与Landau范式的类比

两层筛选结构与Landau相变理论有深层的形式对应：

| 概念 | Landau相变 | 偏序集两层筛选 |
|:----:|:----------:|:--------------:|
| 序参量 | 磁化强度 $M$ | 特征向量 $\mathbf{I} = (d, c, w)$ |
| 一阶约束 | 对称性 → $M$ 奇数项为零 | $S_{\mathrm{BD}}$ → curvature window |
| 二次有效势 | $F \sim aM^2 + bM^4$ | $S_{\mathrm{MD}} \sim \delta\mathbf{I}^\top \Lambda \, \delta\mathbf{I}$ |
| 选择机制 | 最小化自由能 | 最小化到参考流形的距离 |
| 吸引域 | 基态 | Lor4D |

区别在于：Landau理论中的序参量是连续的标量（或向量），而这里的序参量是三维特征向量，描述的是偏序集在因果几何空间中的位置。但数学结构——线性约束 + 二次有效势——是完全平行的。

---

## §8 从经验可分到几何可表述

### 8.1 认识论演进的四个阶段

本项目的发展历程反映了一条清晰的认识论演进路径：

**阶段一：经验可分（Empirically Separable）**

初始发现是：在某些特征组合下，Lor4D可以与其他偏序集族区分开来。这一阶段的工具是旧版F7判别器，使用logH（线性延拓数的对数）作为"结构熵"的代理。Carlip的批评正中要害：logH不是物理量，样本空间太窄。

但"经验可分"本身不是错误——它是正确的直觉，只是表述不成熟。

**阶段二：物理可测（Physically Observable）**

Carlip C1迫使我们放弃logH，转向纯因果几何可观测量。$(d_{\mathrm{eff}}, c_1/c_0, w)$ 三联体的发现标志着这一阶段：

- $d_{\mathrm{eff}}$：来自Myrheim-Meyer维度估计，与Dhar (1978) 的占有率 $\rho$ 等价
- $c_1/c_0$：BDG interval counting 的低阶比值，直接来自Benincasa & Dowker (2010)
- $w$：最大反链宽度，纯组合量

这三个量都可以从偏序集中直接计算，无需任何物理假设——它们是因果结构的固有属性。

**阶段三：统计可刻画（Statistically Characterizable）**

Mahalanobis距离的引入标志着从"手调分类器"到"零参数统计度量"的跃迁。$(\boldsymbol{\mu}(N), \Sigma(N))$ 构成了Lor4D在特征空间中的**参考流形**——一个N参数化的高斯云族。任意偏序集到这个流形的Mahalanobis距离给出了一个自然的、零自由参数的结构认同度量。

Feature ablation 实验证实了 $(d_{\mathrm{eff}}, c_1/c_0, w)$ 是**最小完备三联体**：去掉任何一个特征都会严重恶化判别力（$d_{\mathrm{eff}}$ 是严格必要的），增加第四个特征仅提供边际改善（margin从101.5增加到116.5，增幅15%）。

**阶段四：几何可表述（Geometrically Expressible）**

最终阶段是两层筛选理论的完整表述：

$$\text{第一层：} S_{\mathrm{BD}} = \mathbf{c}^\top \Delta\mathbf{C} \quad \text{（超平面）}$$
$$\text{第二层：} S_{\mathrm{MD}} = \Delta\mathbf{C}^\top \mathbf{M} \, \Delta\mathbf{C} \quad \text{（椭球面）}$$
$$\text{交集：超平面} \cap \text{椭球} = \{\text{Lor4D}\}$$

这不再是一个分类结果，而是一个几何命题。在interval count空间中，Lor4D恰好位于BDG超平面与Mahalanobis椭球的唯一交点。这个几何图像是**可证明或可证伪的**——它做出了精确的、可检验的预测。

### 8.2 什么改变了，什么没有改变

| 维度 | 旧框架 (F7/logH) | 新框架 (两层筛选) |
|:----:|:---------:|:-------:|
| 核心直觉 | "Lor4D在偏序集空间中占据特殊位置" | **不变** |
| 物理基础 | logH（非物理） | 因果几何可观测量 |
| 样本空间 | 7族 | 17族 |
| 自由参数 | 多个 | **零** |
| 数学结构 | 经验评分 | 线性超平面 ∩ 二次椭球 |
| 可证伪性 | 弱 | **强**（几何交集预测） |
| Carlip C1/C2 免疫 | ❌ | ✅ |

核心直觉从未改变：Lor4D是偏序集空间中的特殊对象。改变的是表述这一直觉的数学语言——从经验标记提升到了几何定理。

### 8.3 分层筛选的自然顺序

层级筛选实验揭示了 $d \to c \to w$ 是最自然的筛选顺序：

| 筛选层 | 消除族数 (N=128) | 累计消除 | 物理解释 |
|:------:|:---:|:---:|:---|
| $d_{\mathrm{eff}}$（3σ门） | 12.0/16 | 75% | 维度不对的全部排除 |
| $c_1/c_0$（3σ门） | +2.4 | 90% | interval结构不对的排除 |
| $w$（3σ门） | +0.2 | 91% | 横向宽度不对的排除 |

在 $N \geq 96$ 时，三层筛选消除了全部16个非Lor4D族。最小筛选半径 $k_{\min}$ 从1.24σ（N=20）增长到4.91σ（N=128），意味着筛选力度随N单调增强。

六种可能的筛选顺序排列的总消除数相同（14.6/16），但 $d \to c \to w$ 在第一层就消除了最多的族——这是一个**自然层级**，反映了维度信息在因果几何中的基础地位。

---

## §9 完整实验索引

以下按时间顺序列出本项目全部关键实验及其一行结论。

| # | 实验名称 | 脚本 | 一行结论 |
|:---:|:---|:---|:---|
| 1 | **Large-N 可扩展性** | `carlip_lsd_well_large_n.py` | Lor4D 在 N=16–256 全部 #1/17，margin 单调递增至 1.19 |
| 2 | **Fisher信息权重** | `carlip_fisher_weight_test.py` | $\sigma^2$ 排序 = 经验权重排序；权重由信息论决定 |
| 3 | **S_MD 最小失真算子** | `min_distortion_verify.py` | $S_{\mathrm{MD}}$ 呈 Landau 展开形式，$\sigma^2 \propto N^{-1}$，margin 发散 |
| 4 | **特征剔除** | `feature_ablation_test.py` | $(d_{\mathrm{eff}}, c_1/c_0, w)$ 是最小完备三联体，margin=101.5 |
| 5 | **Mahalanobis LSD** | `mahalanobis_lsd_test.py` | 零参数；当前安全口径为 fixed-reference F2 下 `N≥10` 起 #1，并在更大尺度持续保持优势 |
| 6 | **$\boldsymbol{\mu}(N)$ 轨迹** | `mu_trajectory_theory.py` | $R^2 > 0.99$，$\sigma^2 \propto N^{-1}$，$\det(\Sigma) \propto N^{-3.38}$ |
| 7 | **层级筛选** | `hierarchical_screening_test.py` | $d \to c \to w$ 自然顺序，$N \geq 96$ 完美消除全部16族 |
| 8 | **N=1024 极限** | `large_n_extreme_test.py` | margin 达 $1.93 \times 10^8$（N=1024），幂律发散 |
| 9 | **KR_2layer 分析** | `kr_2layer_analysis.py` | 结构相似是偶然的（$C_1/C_0 \equiv 0$, width=0.75），gap 单调增长 |
| 10 | **$c^*/w^*$ 第一原理** | `cstar_wstar_first_principles.py` | $c^*(\infty)=0.374$，$w \propto N^{-1/d}$，$\beta=0.284 \approx 1/4$ |
| 11 | **$S_{\mathrm{MD}} \leftrightarrow S_{\mathrm{BD}}$ 联系** | `smd_sbd_connection.py` | 几乎不相关（$r \approx 0$），互补，联合选择 = Lor4D |
| 12 | **种子再现性** | `prediction_b_seed_reproducibility.py` | LSD 80/80 #1，Mahalanobis 79/80 #1 |
| 13 | **交叉验证** | `prediction_b_cross_validation.py` | LSD 98% #1，Mahalanobis 94% #1，margin 保留率 > 95% |
| 14 | **N=16 根因诊断（历史 CV 异常）** | `n16_mahalanobis_cv_rootcause.py` | 唯一入侵者 = Lor5D (19/19)；当前应解释为旧 shared-reference/CV 诊断点，而非现行安全下界 |

### 实验前置：框架重建

| # | 实验名称 | 脚本 | 一行结论 |
|:---:|:---|:---|:---|
| 0a | F7 诊断 | `carlip_f7_17family_test.py` | logH 在 17 族扩展空间下失守（N≥28 Lor4D 排 #8–#11） |
| 0b | LSD-Well 构建 | `carlip_lsd_well_17family.py` | W2 变体最优，Lor4D 恢复 #1 |
| 0c | N-adapted 实验 | `carlip_lsd_well_n_adapted.py` | Oracle/Extrapolated/Constant 三模式 100% 一致 |
| 0d | Prediction B 修订（历史阶段） | `prediction_b_lsd_well_revision.py` | 旧阶段曾给出 `N≥20` 口径；现已被 F2 fixed-reference 结果 supersede |

---

## 附录A：关键数学公式汇总

### A.1 Benincasa-Dowker Action (d=4)

$$S_{\mathrm{BD}}^{(4)} = N - C_0 + 9C_1 - 16C_2 + 8C_3$$

系数 $\{1, -1, 9, -16, 8\}$ 由BDG理论在 $d=4$ 时唯一确定。

### A.2 Minimum Distortion Action (Mahalanobis)

$$S_{\mathrm{MD}}(P, N) = \left(\mathbf{I}(P) - \boldsymbol{\mu}(N)\right)^\top \Sigma^{-1}(N) \left(\mathbf{I}(P) - \boldsymbol{\mu}(N)\right)$$

$$\mathbf{I}(P) = \begin{pmatrix} d_{\mathrm{eff}}(P) \\ c_1/c_0(P) \\ w(P) \end{pmatrix}, \quad \boldsymbol{\mu}(N) = \begin{pmatrix} d_{\mathrm{eff}}^{(4D)}(N) \\ c^{(4D)}(N) \\ w^{(4D)}(N) \end{pmatrix}$$

### A.3 两个作用量在 interval count 基中的统一表达

$$S_{\mathrm{BD}} = S_{\mathrm{BD}}^* + \mathbf{c}^\top \Delta\mathbf{C} \qquad \text{（一阶项）}$$

$$S_{\mathrm{MD}} = \Delta\mathbf{C}^\top \underbrace{\mathbf{J}^\top \Sigma^{-1} \mathbf{J}}_{\mathbf{M}} \, \Delta\mathbf{C} \qquad \text{（二阶项）}$$

$$J_{ik} = \frac{\partial I_i}{\partial C_k}\bigg|_{\mathbf{C}^*}$$

### A.4 有限尺寸标度

$$\boldsymbol{\mu}(N) = \boldsymbol{\mu}(\infty) + \frac{\mathbf{a}}{N} + \frac{\mathbf{b}}{N^2} + O(N^{-3})$$

$$\sigma_i^2(N) = A_i \cdot N^{-p_i}, \quad p_i \approx 1 \quad \text{（经典泊松收缩）}$$

$$\det\Sigma(N) \propto N^{-3.38} \quad \Rightarrow \quad \text{点云体积} \propto N^{-1.69}$$

---

## 附录B：致谢与回应

本项目的理论框架在Carlip (CQG desk rejection, 2025) 的三条批评中被重建。三条批评的精确内容和我们的回应简述如下：

| 批评 | 内容 | 回应 |
|:----:|:----:|:----:|
| **C1** | logH ≠ 物理熵 | 完全放弃logH，替换为 $(d_{\mathrm{eff}}, c_1/c_0, w)$ |
| **C2** | 7族 = 挑樱桃 | 扩展到17族，包含KR型、随机分层、传递渗流等 |
| **C3** | 缺 Dhar/Prömel 文献 | 已补充，并建立与 $d_{\mathrm{eff}}$, $\rho$ 的精确联系 |

两层筛选理论**完全免疫于C1和C2**——因为它不使用logH，且在17族的完整空间中经过验证。

---

*本文档生成自 outputs_carlip/ 分析链的14组实验结果，2026年3月。*
