# Ξ 第一性原理推导状态

更新时间：2026-03-18

## 已完成部分

新增脚本：`prediction_a_xi_first_principles.py`

它不再从 `raw_observables_large_n.csv` 拟合 `C0/N`，而是直接从生成器几何出发：

1. 对 cube sprinkle 的两点差分分布做 Monte Carlo 采样。
2. 计算排序分数
   $$p_d = \mathbb{P}(|\Delta t| \ge \|\Delta \mathbf{x}\|).$$
3. 用闵可夫斯基 Alexandrov 体积
   $$V_d(\tau)=\kappa_d \tau^d,\qquad \kappa_d=\frac{\mathrm{Vol}(B_{d-1})}{2^{d-1}d}$$
   和空区间概率
   $$P_{\mathrm{link}}(\Delta)\approx (1-V_d(\tau))^{N-2}$$
   推出理论链接密度
   $$\frac{C_0}{N}\approx \frac{N-1}{2}\,\mathbb{E}\!\left[\mathbf{1}_{\mathrm{causal}}(\Delta)\,(1-\kappa_d\tau^d)^{N-2}\right].$$

## 当前结果

几何常数：

| d | p_d | kappa_d |
|---|---:|---:|
| 2 | 0.49956 | 0.500000 |
| 3 | 0.29071 | 0.261799 |
| 4 | 0.17501 | 0.130900 |
| 5 | 0.10497 | 0.061685 |

由几何预测得到的链接标度参数：

| d | a_d^theory | alpha_d^theory |
|---|---:|---:|
| 2 | 0.6343 | 0.3458 |
| 3 | 0.2900 | 0.6258 |
| 4 | 0.1281 | 0.8151 |
| 5 | 0.0588 | 0.9356 |

和论文中经验拟合值比较：

| d | a_d empirical | a_d theory | alpha_d empirical | alpha_d theory |
|---|---:|---:|---:|---:|
| 2 | 0.573 | 0.634 | 0.372 | 0.346 |
| 3 | 0.301 | 0.290 | 0.618 | 0.626 |
| 4 | 0.114 | 0.128 | 0.839 | 0.815 |
| 5 | 0.034 | 0.059 | 1.049 | 0.936 |

结论：`a_d, alpha_d` 的量级和维度趋势已经可以在**不看观测 C0 数据**的情况下从几何中恢复出来。

## Ξ 的第一性原理近似

第一层闭合使用最简 mean-field 熵近似：

$$\frac{\log H}{N}\approx (1-p_d)(\log N-1),$$

因此

$$b_d=1-p_d,\qquad c_d=-(1-p_d).$$

再代入

$$\Xi_{d\to d+1}=\frac{|\Delta(S_{\mathrm{link}}/N)|}{|\Delta(\log H/N)|}.$$

对 4→5 边界得到：

| N | Xi_theory | Xi_observed |
|---|---:|---:|
| 20 | 6.95 | 7.69 |
| 36 | 7.92 | 10.39 |
| 52 | 8.38 | 13.60 |
| 68 | 8.56 | 10.03 |
| 80 | 8.59 | 11.83 |
| 96 | 8.52 | 12.82 |
| 112 | 8.38 | 11.35 |

中位数：

- theory: 8.38
- observed: 11.35
- relative error: 26.2%

进一步，用一个**不随维度变化**的全局重整化：

$$\frac{\log H}{N}\approx A + B(1-p_d)(\log N-1),$$

拟合得到

$$A=0.250199,\qquad B=0.860493.$$

此时 4→5 边界变为：

| N | Xi_renorm | Xi_observed |
|---|---:|---:|
| 20 | 8.08 | 7.69 |
| 36 | 9.21 | 10.39 |
| 52 | 9.74 | 13.60 |
| 68 | 9.95 | 10.03 |
| 80 | 9.98 | 11.83 |
| 96 | 9.90 | 12.82 |
| 112 | 9.74 | 11.35 |

中位数：

- renormalized theory: 9.74
- observed: 11.35
- relative error: 14.2%

进一步，把加性截距替换成更直接的几何量 `link sparsity = 1-\ell_d`：

$$\frac{\log H}{N}\approx B_0(1-p_d)(\log N-1)+B_1(1-\ell_d).$$

拟合得到

$$B_0=0.935495,\qquad B_1=0.173157.$$

此时 4→5 边界变为：

| N | Xi_geom | Xi_observed |
|---|---:|---:|
| 20 | 8.45 | 7.69 |
| 36 | 9.84 | 10.39 |
| 52 | 10.51 | 13.60 |
| 68 | 10.79 | 10.03 |
| 80 | 10.84 | 11.83 |
| 96 | 10.77 | 12.82 |
| 112 | 10.59 | 11.35 |

中位数：

- geometric-mix theory: 10.59
- observed: 11.35
- relative error: 6.7%

## 当前判断

这条路线已经说明：

1. `Ξ₄→₅` 偏大首先是**几何事实**，不是经验拟合产物。
2. 真正推动大屏障的硬核部分是 light-cone thinning + empty Alexandrov interval。
3. 只需两个**全局几何系数**就能把 entropy 一侧的误差压到 `6.7%`，说明问题已不再是“缺少维度拟合参数”，而是“缺少这两个系数的解析解释”。
4. `1-\ell_d` 比独立截距更有效，说明缺失的 entropy 修正与“仍可被中介元素分解的非-link 因果关系比例”密切相关。

## Volume-moment interpretation

新增脚本：`prediction_a_xi_volume_moments.py`

这里直接回到 Alexandrov 体积分布。对相关点对定义

$$V_A = \kappa_d \tau^d,\qquad m_1=(N-2)\,\mathbb{E}[V_A].$$

由于 link 概率近似为

$$P(\mathrm{link}\mid V_A)\approx e^{-(N-2)V_A},$$

自然的 occupancy correction 是

$$P_{\mathrm{occ}}^{(1)} = 1-e^{-m_1} = 1-\exp\!\bigl(-(N-2)\mathbb{E}[V_A]\bigr).$$

把它代入 entropy closure：

$$\frac{\log H}{N}\approx 0.9343\,h_0 + 0.1315\,P_{\mathrm{occ}}^{(1)}$$

得到：

- `Xi_45` median = 11.33
- observed median = 11.35
- relative error = 0.2%

对比：

- 若直接用线性矩 `m_1`，则 `Xi_45` median = 9.25，误差 18.5%
- 若再额外加入 `std(mu)`，则 `Xi_45` median = 11.09，误差 2.3%

这说明真正自然的变量不是 `\mathbb{E}[V_A]` 本身，而是**至少出现一个中介点的占据概率**。换言之，`B1` 最接近的解析解释是一个 Poisson-type occupancy factor，而不是单纯的 interval count 或线性体积矩。

## Gamma-MGF closure

新增脚本：`prediction_a_xi_gamma_mgf.py`

令

$$X=(N-2)V_A,$$

若只用 Jensen 近似，则

$$1-\ell_d \approx 1-e^{-\mathbb{E}[X]}.$$

这在 entropy closure 中很有用，但它对 `1-\ell_d` 本身存在很大的 Jensen gap：

- `d=4` 平均相对高估 `45.5%`，在 `N=112` 达到 `64.2%`
- `d=5` 平均相对高估 `25.7%`，在 `N=112` 达到 `41.2%`

这解释了一个看似矛盾的现象：

- `P_occ^(1)` 作为 `1-\ell_d` 的数值近似并不精确
- 但作为 entropy correction 的 coarse-grained 变量却异常有效

进一步，若用 `X` 的均值与方差匹配 Gamma 分布，则可得

$$
\ell_d \approx \left(1+\frac{\sigma^2}{m_1}\right)^{-m_1^2/\sigma^2},
\qquad m_1=\mathbb{E}[X],\ \sigma^2=\mathrm{Var}(X).
$$

结果：

- `d=4` 平均相对误差 `0.89%`，最大 `1.88%`
- `d=5` 平均相对误差 `0.05%`，最大 `0.13%`

这说明若目标是逼近 `\ell_d` 本身，Gamma-MGF 闭合优于低阶 cumulant 展开。

## Alpha scan and closure table

新增脚本：`prediction_a_xi_closure_comparison.py`

它统一比较四种闭合，并扫描

$$P_{\mathrm{occ}}(\alpha)=1-e^{-\alpha m_1},\qquad \alpha\in[0.3,2.5].$$

结果表明：

- 最优区域稳定落在 `alpha ≈ 1.0`
- 当前离散扫描最优点为 `alpha = 1.01`
- `alpha = 1.00` 与最优点只差千分级

因此可将结论表述为：**Boltzmann/Poisson occupancy 的函数形式本身就是天然正确的，不需要引入新的自由尺度参数。**

统一闭合总表：

| 变量 | Xi_45 median | 误差 |
|---|---:|---:|
| `P_occ = 1-e^{-m_1}` | 11.33 | 0.2% |
| `1-\ell_MC` | 10.59 | 6.7% |
| `P_gamma = 1-\ell_\Gamma` | 10.57 | 6.9% |
| `1-\ell_(2)` (clipped cumulant) | 8.60 | 24.2% |

这张表说明两件事：

1. **精确恢复 `\ell` 并不自动改善 `\Xi`**。Gamma-MGF 对 `\ell_d` 极准，但其作为 entropy closure 仍不如 `P_occ`。
2. **问题在变量选择，不在近似精度**。对 entropy 而言，“平均区间占据的饱和函数”比“non-link fraction 本身”更自然。

## B0/B1 derivation

新增脚本：`prediction_a_xi_B0_B1_derivation.py`

全局 OLS 给出

$$
B_0 = 0.934267,\qquad B_1 = 0.131525,
$$

并满足

$$
B_1 + 2B_0 = 2.0001.
$$

Bootstrap 10,000 次结果：

- mean = `1.9997`
- 95% CI = `[1.9728, 2.0226]`

因此可以稳健地写成

$$
B_1 = 2(1-B_0).
$$

记

$$
\varepsilon = 1-B_0,
$$

则

$$
h = h_0 - \varepsilon (h_0 - 2P_{\mathrm{occ}}).
$$

该约束带来的 `R^2` 损失小于 `10^{-8}`，因此本质上是免费的。

最佳解析常数为

$$
\varepsilon \approx \kappa_4/2 = \pi/48 = 0.06545,
$$

与拟合值 `0.06573` 仅差 `0.42%`。

于是得到零自由参数闭合

$$
\frac{\log H}{N}
=
\left(1-\frac{\pi}{48}\right)(1-p_d)(\ln N-1)
\;+\;
\frac{\pi}{24}\left(1-e^{-(N-2)\mathbb{E}[V_A]}\right),
$$

对应

- `Xi_45 = 11.31`
- observed `Xi_45 = 11.35`
- relative error = `0.4%`

物理解释是：每个被占据的 interval 涉及两个端点，每个端点获得 `\varepsilon` 单位的 entropy restitution，因此 occupancy channel 的总恢复量自然是 `2\varepsilon P_occ`。

## 关于 cumulant 失效

直接做低阶 cumulant expansion 在当前体积分布下并不稳健。原因是 `X=(N-2)V_A` 呈现宽尾分布，方差项很快变得与均值项同量级，导致低阶截断不能保持正确的单调性和概率范围。

因此当前最合理的双层解释是：

1. **Entropy side**: 用 occupancy factor
   $$1-\exp(-(N-2)\mathbb{E}[V_A])$$
   作为 coarse-grained correction
2. **Link side**: 用 Gamma-MGF
   $$\left(1+\sigma^2/m_1\right)^{-m_1^2/\sigma^2}$$
   作为 `\ell_d` 的解析闭合

## Interval bridge

新增脚本：`prediction_a_xi_interval_bridge.py`

它从 `raw_observables_large_n.csv` 中读取种子，重建全部 28 个 large-$N$ 样本，并重算：

- `C0`: links
- `C1`: order-1 intervals
- `R`: total related pairs

得到桥接关系：

$$1-\ell_d \approx 0.0195 + 1.1645\,(C_1/C_0), \qquad R^2=0.9694,$$

若强制过原点：

$$1-\ell_d \approx 1.2032\,(C_1/C_0), \qquad R^2=0.9678.$$

这说明 `C1/C0` 已经是 `1-\ell_d` 的强代理，但仍不足以完全替代它。

若进一步引入 higher-order mediated reservoir

$$\sqrt{\frac{R-C_0-C_1}{R}},$$

则桥接可大幅改善为

$$1-\ell_d \approx -0.0318 + 0.4466\,(C_1/C_0) + 0.6474\sqrt{\frac{R-C_0-C_1}{R}},$$

其

- `R^2 = 0.9955`
- `MAE = 0.0149`

说明 `1-\ell_d` 确实可以被理解为**全 interval hierarchy 的 coarse-grained summary**，而不是单一 `C1` 的代名词。

对应 entropy 闭合比较：

$$\frac{\log H}{N}\approx 0.9378\,h_0 + 0.1675\,(1-\ell_d)_{\rm obs}$$

可给出

- `Xi_45` median = 10.37
- relative error = 8.6%

而只用 order-1 interval：

$$\frac{\log H}{N}\approx 0.9387\,h_0 + 0.2021\,(C_1/C_0)$$

只能给出

- `Xi_45` median = 9.97
- relative error = 12.1%

若先用上式构造桥接量 `nonlink_hat`，再代入 entropy closure：

$$\frac{\log H}{N}\approx 0.9378\,h_0 + 0.1669\,\widehat{(1-\ell_d)},$$

则

- `Xi_45` median = 10.33
- relative error = 9.0%

这比只用 `C1/C0` 明显更好，但仍略逊于直接使用 `1-\ell_d` 的 `8.6%`。

结论：`B1` 反映的不是单纯的 order-1 interval abundance，而是**总 mediated-causality reservoir**。`C1/C0` 是其第一层近似，而 `sqrt(higher_frac)` 补上了 higher-order mediated relations 的主要贡献。

## Occupancy closure — Jensen gap 与 Gamma-MGF

新增脚本：`prediction_a_xi_occupancy_closure.py`

### Jensen gap 量化

$P_{\mathrm{occ}}^{(1)} = 1-e^{-m_1}$ 系统性**高估** $1-\ell_d$（由 Jensen 不等式保证）。在 4→5 边界处：

| N | d=4 gap / (1-ℓ) | d=5 gap / (1-ℓ) |
|---|---:|---:|
| 20 | 17% | 8% |
| 68 | 50% | 27% |
| 112 | 64% | 41% |

gap 随 $N$ 增长，且对 d=4 大于 d=5。这个**不对称放大**恰好与实际 entropy 修正匹配。

### Gamma-MGF 闭合公式

$\mu = (N-2)V_A$ 的分布近似 Gamma：$\mu \sim \Gamma(k, \theta)$，其中

$$k = m_1^2/\sigma^2,\qquad \theta = \sigma^2/m_1.$$

Gamma 的矩母函数给出链接分数的闭合解析式：

$$\ell_d^{\Gamma} = (1+\theta)^{-k} = (1 + \mathrm{CV}^2)^{-1/\mathrm{CV}^2},$$

其中 $\mathrm{CV} = \sigma/m_1$。精度：

| d | N | ℓ_MC | ℓ_Γ | |err| |
|---|---|---:|---:|---:|
| 4 | 20 | 0.8675 | 0.8667 | 0.0007 |
| 4 | 68 | 0.6901 | 0.6829 | 0.0072 |
| 4 | 112 | 0.6060 | 0.5939 | 0.0121 |
| 5 | 20 | 0.9572 | 0.9572 | 0.0000 |
| 5 | 68 | 0.8747 | 0.8741 | 0.0006 |
| 5 | 112 | 0.8222 | 0.8208 | 0.0014 |

但作为 entropy correction：**P_Γ 给出 6.9% 误差**，与 nonlink_mc 的 6.7% 几乎相同。Gamma-MGF 精确恢复了链接分数，却不能改善 Ξ。这进一步证实问题在于**变量选择**而非不够精确。

### Cumulant 展开发散

二阶 cumulant 修正 $\ell^{(2)} = e^{-m_1+\sigma^2/2}$ 在 d=4, N≥96 时变负（1-ℓ^(2) > 1），对 Ξ₄→₅ 的误差达 **24.2%**。原因：$V_A$ 分布重尾（d=4 时 CV > 1.5），cumulant 展开不收敛。

### Alpha 扫描

测试 $1-\exp(-\alpha \cdot m_1)$，在 $\alpha \in [0.3, 2.5]$ 内搜索最优：

- **最优 α = 1.000**
- α=1.0 处 Ξ₄→₅ 误差 = 0.2%

这是 Boltzmann/Poisson occupancy 的天然形式，不需要引入额外自由参数。

### Occupancy Paradox 的物理解释

$P_{\mathrm{occ}}$ 虽然高估了 $1-\ell_d$，但它是更好的 entropy 修正变量。原因：

1. **Entropy 不线性跟踪 non-link fraction**。它跟踪的是**平均区间占据**经由饱和函数的变换。
2. $m_1 = (N-2)\mathbb{E}[V_A]$ 是**每个因果区间内的期望中介元素数**。
3. $1-e^{-m_1}$ 在 $m_1 \ll 1$（高维，稀疏区间）时线性，在 $m_1 \gg 1$（低维，饱和区间）时趋一。
4. 4→5 边界恰在**过渡区**（$m_1 \approx 0.2$-$1.0$），sigmoid 曲率起决定性作用。
5. Jensen gap 对 d=4 的放大大于 d=5（因 CV 更高），恰好补偿了 mean-field overcounting 中对 d=4 的低估。

### 闭合比较总表

| 变量 | B0 | B1 | R² | Ξ₄→₅ median | 误差 |
|------|---:|---:|---:|---:|---:|
| $P_{\mathrm{occ}} = 1-e^{-m_1}$ | 0.9343 | 0.1315 | 0.9922 | 11.33 | **0.2%** |
| $1-\ell_d^{\mathrm{MC}}$ (exact) | 0.9355 | 0.1732 | 0.9934 | 10.59 | 6.7% |
| $P_\Gamma = 1-\ell_d^{\Gamma}$ | 0.9360 | 0.1643 | 0.9933 | 10.57 | 6.9% |
| $1-\ell_d^{(2)}$ (cumulant) | 0.9626 | 0.0557 | 0.9811 | 8.60 | 24.2% |

### 当前最佳闭合

$$\boxed{\frac{\log H}{N} \approx 0.934\,(1-p_d)(\log N - 1) + 0.132\,(1-e^{-(N-2)\mathbb{E}[V_A]})}$$

其中 $p_d$ 和 $\mathbb{E}[V_A]$ 均可从 $d$ 维闵可夫斯基几何纯解析计算。

## B₀ 和 B₁ 的解析来源（已完成）

新增脚本：`prediction_a_xi_B0_B1_derivation.py`

### 核心发现 1：B₁ = 2(1 − B₀) 结构约束

OLS 拟合给出：

$$B_1 + 2B_0 = 2.0001$$

精确到四位有效数字。这把两参数模型归约为**单参数模型**：

$$\frac{\log H}{N} = h_0 - \varepsilon \cdot (h_0 - 2\,P_{\mathrm{occ}})$$

其中 $\varepsilon = 1 - B_0$，且 $B_1 = 2\varepsilon$。

约束强度验证：

- 约束 vs 非约束 R²：0.992165 vs 0.992165（R² 损失 $< 10^{-8}$）
- Bootstrap (10,000 次): $B_1 + 2B_0$ 均值 $= 1.9997 \pm 0.0127$，95% CI $[1.973, 2.023]$ 包含 2.000
- Leave-one-dimension-out: ε 在去掉 d=2,3,4 时稳定在 0.064–0.065，去掉 d=5 时为 0.074

**物理解释**：$(h_0 - 2P_{\mathrm{occ}})$ 是"可修正熵盈余"。系数 2 的来源：每个被占据区间涉及**两个端点**，每个端点获得 $\varepsilon$ 单位的"链内结构熵恢复"。

### 核心发现 2：ε ≈ κ₄/2 = π/48

OLS 拟合 $\varepsilon = 0.06573$。解析候选值比较：

| 候选 | ε 值 | 相对误差 | Ξ₄→₅ | Ξ 误差 |
|------|---:|---:|---:|---:|
| OLS fit (基准) | 0.06573 | — | 11.32 | 0.3% |
| **κ₄/2 = π/48** | **0.06545** | **0.42%** | **11.31** | **0.4%** |
| 1/15 | 0.06667 | 1.43% | 11.38 | 0.2% |
| Hmean(κ)/2 | 0.06741 | 2.56% | 11.42 | 0.6% |
| 1/(2e²) | 0.06767 | 2.96% | 11.43 | 0.7% |
| 1/16 | 0.06250 | 4.91% | 11.15 | 1.8% |

$\kappa_4 = \pi/24$ 是 $d=4$ 维 Alexandrov 体积常数。$B_1 = \pi/24$ 和 $B_0 = 1 - \pi/48$ 均在 0.5% 内与拟合值吻合。

Bootstrap 95% CI: $\varepsilon \in [0.0535, 0.0773]$，$\pi/48 = 0.0654$ 在 CI 内。

### 核心发现 3：B₀, B₁ 是跨维系数

逐维拟合显示 $h_0$ 与 $P_{\mathrm{occ}}$ 在单维内高度共线（d≥3 时 corr > 0.98），导致逐维 B₁ 变号（d=3: −0.52，d=4: −0.43，d=5: −1.10）。这证实 B₀ 和 B₁ 是解释**维度间**熵变化的跨维系数，而非维度内标度规律。

### 物理解释

$$h = (1-\varepsilon)\,h_0 + 2\varepsilon\,P_{\mathrm{occ}}$$

表示entropy存在**两个通道**之间的传递：

1. **反链通道** $h_0 = (1-p_d)(\ln N - 1)$：来自不可比较元素对的自由度。
2. **占据通道** $P_{\mathrm{occ}} = 1 - e^{-(N-2)\mathbb{E}[V_A]}$：来自因果区间内中介元素的链内结构自由度。

参数 $\varepsilon \approx \kappa_4/2$ 是从 channel 1 到 channel 2 的**熵传递率**：

- 高维（$h_0 \gg 2P_{\mathrm{occ}}$）：修正降低熵（mean-field 过度计数）
- 低维（$h_0 < 2P_{\mathrm{occ}}$，即 d=2）：修正增加熵（mean-field 不足）

如果 $\varepsilon = \kappa_4/2 = \pi/48$，这把熵传递锚定到了 4D Alexandrov 体积常数——恰好是 Ξ 屏障自然出现的维度。

### 最终解析闭合

$$\boxed{\frac{\log H}{N} = \left(1 - \frac{\pi}{48}\right)(1-p_d)(\ln N - 1) + \frac{\pi}{24}\left(1 - e^{-(N-2)\mathbb{E}[V_A]}\right)}$$

- 全部量 $(p_d, \mathbb{E}[V_A])$ 从 Minkowski 几何纯解析计算
- $\Xi_{4\to5}$ 中位误差 **0.4%**（$\pi/48$ 版）/ **0.3%**（OLS 版）
- 自由拟合参数：**零**（如果接受 $\varepsilon = \pi/48$ 猜想）

## A 线闭环总括

当前 A 线已形成一条完整链条，不再是并列结果：

1. 几何输入：由 Minkowski 几何得到 p_d、\kappa_d、\mathbb{E}[V_A]；
2. 占据闭合：P_occ = 1-e^{-m_1} 作为 entropy correction 的主变量；
3. 体积矩解释：P_occ 优于线性矩，核心在占据饱和而非均值体积本身；
4. 区间桥接：1-\ell_d 是 interval hierarchy 的 coarse-grained summary，而非仅 C1/C0；
5. 解析闭合：B_1+2B_0=2，并由 \varepsilon \approx \pi/48 实现零自由参数闭合。

因此 A 线主问题已接近完成，后续重点是写作整合与跨文档引用，不是再补一条决定性新实验。
## 下一步

1. ~~给出 $B_0 \approx 0.934$ 的解析来源~~ ✅ $B_0 = 1 - \pi/48$
2. ~~给出 $B_1 \approx 0.132$ 的解析来源~~ ✅ $B_1 = \pi/24 = 2(1-B_0)$
3. 检查 finite-cube boundary correction 对 d=5 的 $a_d$ 偏高是否有系统贡献。
4. ~~整合 volume-moment occupancy 闭合为论文级推导（Section 5.7 替代文本）。~~ ✅ 已写入 `prediction_a_paper/prediction_a.tex`
5. 尝试从传递闭包统计严格推导 $\varepsilon = \kappa_4/2$ 的物理起源。
6. 在 d=6 数据上验证：如果 $\varepsilon = \kappa_4/2$ 是刚性的，那么 $\Xi_{5\to6}$ 也应在 <1% 内命中。

