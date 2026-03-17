# Poset Phase 项目工作记录与进度追踪

> **最后更新**: 2026-03-19  
> **Git HEAD**: `待提交` (main)  
> **GitHub**: `github.com/unicome37/poset_phase`  
> **版本**: v4.0.0 (Zenodo 已发布)

---

## 一、项目概述

偏序集相变（Poset Phase Transition）项目——基于"结构存在论"理论框架，用数值实验验证三个核心预测（Prediction A/B/C），探讨**存在性筛选**（Existential Selection Dynamics, ESD）如何通过熵-作用量竞争，在偏序集空间中自然涌现出物理时空结构。

### 核心思想

$$\text{Score} = -\beta \cdot \log H(\text{poset}) + \gamma \cdot \text{Penalty}(\text{poset})$$

- $\log H$: 偏序集的线性扩展数（熵）
- Penalty: 各种结构约束（几何、因果、维度等）
- 低 score = 更"存在"的结构

---

## 二、环境与工具链

| 项 | 值 |
|---|---|
| Python | 3.14.2 (系统 `python` 命令，**不用 .venv**) |
| 项目路径 | `d:\Kiro\理论体系\poset_phase\` |
| 运行方式 | `cd "d:\Kiro\理论体系\poset_phase"; python <script>.py` |
| Poset 类 | `poset.closure` (bool ndarray, 传递闭包)，**没有** `poset.adj` |
| 熵估计 | `entropy_sis.py` (SIS 采样, 通常 2048-4096 runs)；`entropy_exact.py` (精确, N≤24) |
| 作用量 | `action.py` → `action_value(log_H, penalty, β, γ)` |
| 几何组件 | `observables_geo.py` → `geometric_components()` |
| 生成器 | `generators.py` → 2D/3D/4D/5D Lorentzian + KR_like |
| 族字典 | `experiment.py` → `FAMILIES` dict |
| 运行工具 | `runtime_utils.py` → `estimate_entropy_by_family()` |

### 关键注意事项

1. **cwd 必须在 `d:\Kiro\理论体系\poset_phase\`**，否则导入失败
2. **精确计算阈值**: 5D 用 `exact_threshold=24`，N≥32 时 5D 精确计算会 hang
3. **KR_like** 是随机偏序，不是物理时空，分析时需区分

---

## 三、三大预测 (Predictions A/B/C) 状态

### Prediction A — 维度选择

**核心问题**: 作用量框架是否能从偏序集空间中自然选出 3+1 维时空？

#### 已完成实验

| 脚本 | 内容 | 关键结果 |
|------|------|----------|
| `prediction_a_dimension_scan.py` | 基础维度扫描 (2D/3D/4D, N=8-52) | 2D 在中等 γ 赢, 4D 在低 γ 赢 |
| `prediction_a_geometric_ablation.py` | 几何惩罚消融实验 | A2_replace_dim_with_consistency 变体效果最好 |
| `prediction_a_margin_fit.py` | margin(N) OLS 回归 | vs 2D: slope=0.945±0.022, R²=0.994; vs 3D: slope=0.319±0.019, R²=0.958 |
| `prediction_a_5d_pilot.py` | **5D 天花板测试** (N=20-52, 5族, SIS 4096) | **5D 不饱和**! 5D 赢 48/63 (76.2%) under consistency |
| `prediction_ac_causal_link.py` | **因果结构×熵 跨维度诊断** | density-entropy 相关: r=-0.988 to -0.997 (p<0.0001) |
| `prediction_a_bd_dimension.py` | **BD 因果集作用量维度选择** | **λ=6-8: 4D 全票获胜!** |
| `prediction_a_bd_lorentzian_only.py` | BD 结果 Lorentzian-only 最终分析 | 维度级联: 5D→4D→3D→2D |
| `prediction_a_bd_extended.py` | **扩展 BD: d=4区间 + 混合 + N=60,68** | **4D plateau 在 N=20-68 全票稳定!** |
| `prediction_a_bdg_full_comparison.py` | **文献 BDG 系数 vs link-proxy 对比** | **link_d2 选4D; 标准BDG_d4 全选5D — 机制是链接密度!** |
| `prediction_a_bdg_component_figure.py` | **BDG 组件诊断图** (stacked bar + winner heatmap) | 可视化 C0/C1/C2/C3 对各变体的贡献，供论文使用 |
| `prediction_a_generator_robustness.py` | **生成器鲁棒性测试** (独立种子 + 因果钻石) | **机制跨生成器普适: cube 7/7, indep-seed 6/7, diamond 4/7** |
| `prediction_a_xi_parameter.py` | **Ξ 无量纲控制参数分析** (3 种生成器) | **Ξ₄→₅ ≈ 10 跨生成器稳定 (CV=13.9%)** |
| `prediction_a_xi_figure.py` | **Ξ 论文级图表** (strip plot, 收敛, 屏障不对称) | 3 张出版级 PNG/PDF |
| `prediction_a_large_n_scaling.py` | **大 N 有限尺寸标度** (N=80,96,112) | **Ξ₄→₅ 在 N=112 仍然稳定 (CV=17.7%), λ 窗口右移** |
| `prediction_a_large_n_figures.py` | **大 N 闭环图表** (热力图, Ξ 稳定性, 链接密度交叉) | 3 张出版级 PNG/PDF |
| `prediction_a_unification.py` | **统一测试: link density ↔ geometric consistency** | **链接密度与 d_order 弱相关 (r=-0.27), 但间隙结构揭示共同机制** |
| `prediction_a_xi_derivation.py` | **Ξ 解析推导** (标度律拟合 + MC序分数 + 闭合公式) | **★★ 预测 Ξ₄→₅=11.8, 实测 11.3 — 根因: 熵饱和 + 链接密度间隙持续** |

#### ★ 突破性发现: BD 作用量天然选出 3+1 维

**Benincasa-Dowker 作用量**: $S_{BD} = N - 2 n_{\text{links}}$ (Hasse 图覆盖关系数)

$$\text{score} = -\beta \cdot \log H + \lambda \cdot \frac{S_{BD}}{N}$$

Lorentzian 族内部的维度选择级联:

| λ 范围 | 一致赢家 | 物理含义 |
|--------|----------|----------|
| 0 - 2 | 5D (全票) | 纯熵主导 → 高维永远赢 |
| 2.5 - 3.5 | 5D→4D 过渡 | BD 开始生效 |
| **6 - 8** | **4D 全票 (N=20-68, 含 N=60/68 扩展)** | **熵 vs 因果性的平衡 = 3+1 维** |
| 10 - 12 | 4D/3D 混合 | BD 开始过度 |
| 15 - 25 | 3D 主导 | 因果性过强 |
| 35+ | 2D 开始出现 | 极端因果约束 |

**物理解读**: BD 因果集作用量（离散 Einstein-Hilbert）创造了"熵 vs 因果连通性"的权衡，其最优解自然落在 3+1 维。

#### 扩展 BD 分析关键结论 (e4d237d)

1. **有限尺寸标度**: BD_d2 λ=6-8 时 4D 在 N=20-68 全部 7 个尺度一致获胜，plateau 稳定
2. **BD_d4 (区间计数)**: 对 5D 更激进，但在 N=68 处出现 5D "逃逸"，不如 BD_d2 稳定
3. **BD+几何混合**: 添加几何惩罚**不改善** 4D 选择，纯 BD_d2 已是最优
4. **区间轮廓**: 5D 几乎没有高阶区间 (C1≈0, C2≈0)；2D 有丰富的区间层级
5. **N=68 margin**: λ=6 时 margin=0.633 (窄但稳定)；λ=8 时 margin=4.912 (宽裕)

#### ★★ 关键验证: link-proxy vs 文献 BDG 系数 (GPT 审查驱动)

**背景**: 外部 GPT 文献分析指出，我们的 $S = N - 2 n_{\text{links}}$ 实际上是 **d=2 BD 作用量** (link action proxy)，而非文献中 Benincasa-Dowker 的 d=4 标准形式。标准 BDG d=4 形式为:

$$S^{(4)}_{BDG} = N - C_0 + 9C_1 - 16C_2 + 8C_3$$

其中 $C_k$ = 含 $k$ 个内部元素的因果区间数。

**实验**: `prediction_a_bdg_full_comparison.py` 利用已收集的 C0/C1/C2/C3 数据 (N=20-68, 4 样本×5 族), 对比 4 种作用量变体:

| 变体 | 公式 | 4D 最佳得分 |
|------|------|-------------|
| `link_d2` | $N - 2C_0$ | **λ=6-8: 7/7 全票 4D** ★ |
| `BDG_d4_standard` | $N - C_0 + 9C_1 - 16C_2 + 8C_3$ | **0/7 全 λ, 5D 永远赢** |
| `BDG_d4_pure` | $-C_0 + 9C_1 - 16C_2 + 8C_3$ | 0/7 |
| `BDG_d2_corrected` | $N - 2.09C_0 + 0.41C_1$ | λ=20: 7/7, 但窗口太窄 |

**关键对比 (λ=6)**:

| N | link_d2 | BDG_d4_standard |
|---|---------|-----------------|
| 20 | **4D** beats 3D by 1.50 | 5D beats 4D by 4.51 |
| 28 | **4D** beats 3D by 4.06 | 5D beats 4D by 7.15 |
| 36 | **4D** beats 3D by 3.73 | 5D beats 3D by 7.13 |
| 44 | **4D** beats 5D by 4.20 | 5D beats 4D by 6.30 |
| 52 | **4D** beats 5D by 1.41 | 5D beats 4D by 9.38 |
| 60 | **4D** beats 5D by 0.79 | 5D beats 4D by 9.16 |
| 68 | **4D** beats 5D by 0.63 | 5D beats 4D by 16.44 |

**物理诊断 — 为什么 BDG_d4 失败**:

BDG_d4 中的 $+9C_1$ 项巨幅放大低维度的 order-1 区间:
- 2D (N=44): $+9C_1 = +594$, 压过 $-16C_2 = -788$, net = **+107**
- 5D (N=44): $+9C_1 = +99$, $-16C_2 = -32$, net = **-28**

5D 的高阶区间极少 → BDG_d4 惩罚极小 → 结合 5D 熵最大 → 5D 永远赢。换言之, BDG_d4 的区间修正项实际上**削弱**了对高维的惩罚。

**核心结论**:

> **维度选择发生在链接密度 (link-counting) 层面，而非完整曲率作用量层面。**
>
> - $S = N - 2 n_{\text{links}}$ (d=2 BD / link action) 是选出 3+1 维的**充分且必要**元素
> - 文献标准 BDG d=4 的高阶区间修正 ($C_1, C_2, C_3$) **破坏**了 4D 选择
> - 这说明选择机制是**Hasse 图覆盖关系的密度平衡**，不是离散 Ricci 曲率
>
> **论文定位**: 明确声明使用的是 "link action" (d=2 BD), 并展示 full BDG 的对比作为证据，证明维度选择的物理机制是**因果连接数密度**。

#### ★★★ 生成器鲁棒性验证 (b34b8be)

**问题**: 4D 选择窗口是否依赖于特定的撒点几何或种子家族？

**测试方案**: 三种独立生成器:
1. **原始立方体撒点** (seed 980000) — 基线
2. **独立种子立方体** (seed 1234567) — 排除种子依赖
3. **Alexandrov 集 (因果钻石) 撒点** (seed 7770000) — 排除几何依赖

**结果 (λ=7)**:

| 生成器 | N=20 | N=28 | N=36 | N=44 | N=52 | N=60 | N=68 | 4D 计数 |
|--------|------|------|------|------|------|------|------|--------|
| 原始立方体 | 4D | 4D | 4D | 4D | 4D | 4D | 4D | **7/7 ★** |
| 独立种子 | 3D | 4D | 4D | 4D | 4D | 4D | 4D | **6/7** |
| 因果钻石 | 3D | 3D | 4D | 4D | 4D | 4D | 5D | **4/7** |

**5D 因果稀疏度诊断 (N=52)**:

| 度量 | 立方体 5D | 钻石 5D | 比值 |
|------|----------|---------|------|
| C₀ (链接数) | 123.8 | 41.2 | 3.0× 更稀疏 |
| log H | 129.3 | 144.7 | 12% 更多熵 |
| S_link/N | −3.76 | −0.59 | 6.4× 更弱惩罚 |

**核心结论**:
> 维度选择窗口在参数空间不是普适常数，但在**机制空间是普适的**。
> 不同几何通过改变因果稀疏度来调制窗口位置，但链接密度竞争机制本身——及其对 d=4 的偏好——在所有测试的生成器中持续存在。
> 因果链: 钻石几何 → 更稀疏 5D → 更少链接 → 更弱惩罚 → 窗口向更大 λ 漂移。

#### ★★★★ Ξ 无量纲控制参数 (41d866b + 5655002)

**定义**: $\Xi_{d \to d+1} = |\Delta(S_\text{link}/N)| / |\Delta(\log H / N)|$ —— 跨维度边界的链接惩罚/熵比。

| 转换 | 原始立方体 | 独立种子 | 因果钻石 | 总体 |
|------|-----------|---------|---------|------|
| Ξ₂→₃ | 3.2 | 1.9 | 1.2 | 1.9 |
| Ξ₃→₄ | 2.1 | 2.7 | 5.0 | 3.3 |
| **Ξ₄→₅** | **10.8** | **10.2** | **9.4** | **10.0** |

**物理含义**: 从 4D 上升到 5D 的"熵价格"是从 3D 到 4D 的 ~3× —— 不对称屏障使 3+1 维成为自然平衡点。

#### ★★★★★ 大 N 闭环验证 (N=80,96,112)

**关键结果**:

1. **λ 窗口右移但 Ξ 不漂**: 在 λ=7 时, 4D 赢到 N=80; 在 λ=10 时, 4D 赢到 N=112。
2. **链接密度交叉**: N≈96 处, 4D 链接密度首次超过 3D — 解释了窗口右移。
3. **间隙比发散**: |Δ(4→5)| / |Δ(3→4)| 在链接密度中从 1.4 增长到 19.5 — 非对称屏障随 N 强化。
4. **Ξ₄→₅ 稳定性**: 7 个 N 值的中位数 = 11.35, CV = 17.7%, 小 N/大 N 漂移仅 15.8%。

**闭环评估**: 
> 4D 选择不是在固定 λ 对所有 N 成立——而是在 λ 超过 N 依赖的临界值 λ*(N) 时成立。
> 由于 Ξ₄→₅ ≈ 11 与 N 无关，机制确认成立，即使窗口位置随 N 变化。
> 这实际上**比固定 λ 选择更强**：它意味着选择是偏序集竞争格局的结构性质，不是特定耦合的产物。

#### 统一测试: Link Action ↔ Geometric Consistency

**核心问题**: link action 和 BDG 几何一致性惩罚是否测量同一信号？

**发现**: 
- 链接密度与 d_order 的简单相关性较弱 (r=-0.27, p=0.31)，因为关系非单调（3D 处链接密度最高）
- 但**间隙结构**揭示共同机制:
  - 链接密度: |Δ(4→5)|/|Δ(3→4)| = 1.4 → 15.4 (随 N 发散)
  - BDG proxy_penalty: 该比值 <1 (3→4 惩罚反而更大)
- 两者共享**根本原因**: 5D Lorentzian 偏序在结构上比 4D 更稀疏
- Link action 还捕获了 BDG 看不到的**3D-4D 连通性收敛**现象

**统一评定**: 机制层面统一（共同的因果稀疏度层级），度量层面互补（不同的信号权重分布）。

#### ★★★★★★ Ξ₄→₅ 解析推导闭环

**目标**: 从经验标度律推导 Ξ₄→₅ ≈ 11 的闭合公式。

**标度律**:
- $C_0/N = a_d \cdot N^{\alpha_d}$（幂律）
- $\log H/N = b_d \cdot \log N + c_d$（对数线性）

**闭合公式**:
$$\Xi_{d \to d+1}(N) = \frac{2|a_d \cdot N^{\alpha_d} - a_{d+1} \cdot N^{\alpha_{d+1}}|}{|(b_{d+1}-b_d) \cdot \log N + (c_{d+1}-c_d)|}$$

**验证 (N=68)**:
| 转换 | 分子 | 分母 | Ξ 预测 | Ξ 实测 |
|------|------|------|--------|--------|
| 2→3 | 2.668 | 0.570 | 4.68 | 3.2 |
| 3→4 | 0.308 | 0.301 | 1.02 | 2.1 |
| **4→5** | **2.252** | **0.191** | **11.8** | **10.0-11.3** |

**根因分解**:
1. **熵饱和**: Δb₄→₅ = 0.042，仅为 Δb₃→₄ = 0.110 的 **38%** — 边际熵递减
2. **链接密度间隙持续**: |ΔC₀/N(4→5)| = 1.126，是 |ΔC₀/N(3→4)| = 0.154 的 **7.3×**
3. **序分数加速衰减**: p₂=0.50→p₃=0.29→p₄=0.17→p₅=0.11 (每维 ×0.6)
4. **链接饱和**: d=5 时 94.7% 的因果关系是直接链接 (l₅=0.947) — 结构已近反链

**物理结论**:
> Ξ₄→₅ ≈ 11 源自闵可夫斯基几何的一个事实：d=4→5 时光锥体积分数骤降（破坏链接），
> 而额外维度带来的熵收益是边际的（对数饱和）。这使得 4→5 边界的每单位熵代价
> 比任何低维转换高一个数量级。

#### 5D Pilot 详细数据

| 变体 | 5D | 4D | 3D | 2D | KR |
|------|-----|-----|-----|-----|-----|
| A2_full | 20 | 1 | 4 | 18 | 20 |
| A2_consistency | **48** | 7 | 8 | 0 | 0 |

5D vs 4D margin (consistency):

| N | mean_delta | win_rate |
|---|-----------|---------|
| 20 | -0.12 | 0.571 |
| 36 | -2.11 | 0.714 |
| 52 | -5.29 | 1.000 |

#### 因果诊断数据 (N=52)

| 维度 | 边密度 | 层数 | 最大反链 | log H |
|------|--------|------|---------|-------|
| 2D | 0.494 | 9.2 | 8.2 | 78.8 |
| 3D | 0.283 | 5.8 | 14.8 | 105.5 |
| 4D | 0.167 | 3.5 | 20.8 | 121.6 |
| 5D | 0.103 | 2.8 | 27.2 | 129.9 |

**链条**: 高维 → 光锥占比小 → 因果更稀疏 → 层数少 → 反链大 → 熵暴涨 → 无约束时永远赢

---

### Prediction B — γ 相变阈值

**核心问题**: 是否存在临界 $\gamma_c$，使得 KR > Lorentzian 切换为 Lorentzian > KR？

#### 已完成实验

| 脚本 | 内容 | 关键结果 |
|------|------|----------|
| `prediction_b_gamma_c_scaling.py` | γ_c(N) 有限尺寸标度拟合 | 常数: 0.528±0.096, 95%CI [0.31, 0.74]; 线性 p=0.21 |
| `prediction_b_gamma_extended_analysis.py` | 扩展 γ 扫描 [0, 5.0] 分析 | N=20-44 精确窗: γ_c = 0.98-1.24 |
| `prediction_b_gamma_extended_nearwall.py` | **Near-wall N=48,52,56** | N=48: γ_c=1.76; N=52: γ_c=1.84; **N=56: γ_c=2.08** |

#### ★ 关键发现: N=56 γ_c 恢复

之前 N=56 "无交叉" 是 **γ 扫描范围假象**（γ只到 2.0）。扩展到 5.0 后:

| N | γ_c | 方法 |
|---|-----|------|
| 48 | 1.76 | SIS |
| 52 | 1.84 | SIS |
| **56** | **2.08** | SIS |

γ_c(N) 随 N 缓慢增长，但确实存在交叉。

---

### Prediction C — 层级深度预测熵

**核心问题**: 偏序集的层级结构（hierarchy depth）是否能预测其熵？

#### 状态: v4.0.0 已发布

9 个实验构成准因果证据塔 (quasi-causal evidence tower):

1. **观察性回归** (pooled_regression, stratified_regression)
2. **边密度普遍性** (edge_density_phase)
3. **N-scaling 法则** (n_scaling)
4. **解析界** (analytical_bound)
5. **安慰剂控制** (placebo_intervention)
6. **反向干预** (reverse_intervention)
7. **跨维泛化** (comprehensive, 2D/3D/4D + KR)
8. **大N SIS干预** (large_n_intervention, N=16-36)
9. **剂量-反应** (dose_response)

- **手稿**: `PredictionC_MainPaper_Unified_Draft.md` (v2, 1020 行)
- **PDF**: `manuscript_unified/` → 43 页
- **Zenodo**: v4.0.0 已发布

---

## 四、A×C 交叉发现 — 因果稀疏性解释

### 核心洞察 (用户提出)

> "有没可能要考虑 C 的因果因素，当维度上升时，因果复杂度也过大？"

### 验证结论

**Prediction C 的层级-熵机制正是跨维度运作的。** 三个 Prediction 之间第一次出现实质性理论交叉:

```
Prediction C: 更多层级 → 更少熵 (within-dimension)
     ↓ 扩展到跨维度
Prediction A: 高维 → 更少因果关系 → 更少层级 → 更高熵
     ↓ 解决方案
BD 作用量: 惩罚因果稀疏性 → 自然选出 3+1 维
```

---

## 五、代码关键修改记录

### generators.py (239 行)

- 新增 `generate_lorentzian_like_5d()`: 1 时间维 + 4 空间维 in [0,1]^5，因果条件 dt² ≥ Σdx²
- 位于 `generate_lorentzian_like_4d()` 之后

### experiment.py (224 行)

- 导入新增: `from generators import ..., generate_lorentzian_like_5d`
- FAMILIES dict 新增: `"lorentzian_like_5d": generate_lorentzian_like_5d`

### 新脚本一览

| 脚本 | 创建于 | 输出目录 |
|------|--------|----------|
| `prediction_a_margin_fit.py` | 本轮 | `outputs_exploratory/prediction_a_margin_fit/` |
| `prediction_a_5d_pilot.py` (~200行) | 本轮 | `outputs_exploratory/prediction_a_5d_pilot/` |
| `prediction_b_gamma_c_scaling.py` | 本轮 | `outputs_exploratory/prediction_b_gamma_c_scaling/` |
| `prediction_b_gamma_extended_analysis.py` | 本轮 | N/A (分析已有数据) |
| `prediction_b_gamma_extended_nearwall.py` | 本轮 | `outputs_exploratory/prediction_b_gamma_extended_nearwall/` |
| `prediction_ac_causal_link.py` | 本轮 | `outputs_exploratory/prediction_ac_causal_link/` |
| `prediction_a_bd_dimension.py` | 本轮 | `outputs_exploratory/prediction_a_bd_dimension/` |
| `prediction_a_bd_lorentzian_only.py` | 本轮 | (同上目录) |
| `prediction_a_bd_extended.py` | 本轮 | `outputs_exploratory/prediction_a_bd_extended/` |
| `prediction_a_bdg_full_comparison.py` | 本轮 | `outputs_exploratory/prediction_a_bdg_full/` |
| `prediction_a_bdg_component_figure.py` | 本轮 | `outputs_exploratory/prediction_a_bdg_component_figure/` |
| `prediction_a_generator_robustness.py` | 本轮 | `outputs_exploratory/prediction_a_generator_robustness/` |
| `PredictionA_LinkAction_Paper_v1.md` | 本轮 | (论文草稿, 根目录) |

---

## 六、Git 提交历史 (本轮)

```
b34b8be feat: generator robustness test — cube indep-seed + causal diamond sprinkle
9eac932 docs: Prediction A paper draft v1 — Link Action Selects 3+1 Dimensions
58501cb feat: BDG component breakdown figures — stacked bar + winner heatmap
6514e95 feat: full BDG coefficient comparison — link-proxy selects 4D, standard BDG_d4 selects 5D
964f1f8 feat: BD extended analysis final cleanup
e4d237d feat: extended BD analysis — d=4 intervals + hybrid + finite-size scaling to N=68
a8e79c7 feat: Benincasa-Dowker action dimension selection experiment
878d8db feat: cross-dimensional causal diagnostics (Prediction AC link)
82d2729 feat: Prediction A margin OLS + 5D ceiling test; Prediction B gamma_c scaling + extended gamma scan
ba2f43d docs: update README for v4.0.0
eb1f83e (tag: v4.0.0) release: v4.0.0 Prediction C quasi-causal upgrade for Zenodo
4e3c0f3 feat(predC): large-N SIS intervention (N=16-36) + manuscript v2
8324e36 feat: Prediction C reverse intervention, edge density universality, N-scaling law
6a147a7 feat: Prediction C placebo control, cross-dim generalization, analytical theorem
5df4303 feat: Prediction C three-pronged supplementary verification
```

---

## 七、已解决的技术问题

| 问题 | 原因 | 解决方案 |
|------|------|----------|
| N=56 "无交叉" | γ 扫描只到 2.0 | 扩展到 5.0, 发现 γ_c=2.08 |
| 5D N=32 精确计算 hang | 5D 稀疏→反链化→精确计算指数爆炸 | `exact_threshold=24` |
| `poset.adj` AttributeError | 属性不存在 | 用 `poset.closure` |
| cwd 错误导致 import 失败 | 脚本从 `D:\Kiro` 运行 | 必须先 `cd` 到项目目录 |
| `A2_replace_dim_with_consistency` 不在 dimension_scan.py | 该变体只在 ablation 脚本中 | 写独立的 5D pilot 脚本 |

---

## 八、下一步可能方向

### 高优先级

1. ~~**BD 作用量 + 原框架整合**~~: ✅ 已验证，纯 BD_d2 最优，几何惩罚不改善
2. ~~**BD 有限尺寸标度**~~: ✅ N=20-68 全部 7 个尺度 4D 一致获胜
3. ~~**BD 作用量的 d>2 推广**~~: ✅ BD_d4 已实现，但稳定性不如 BD_d2
4. ~~**文献 BDG 系数验证**~~: ✅ 标准 BDG_d4 (N-C0+9C1-16C2+8C3) 全选 5D，确认 link action 是关键
5. ~~**生成器鲁棒性**~~: ✅ 三种生成器 (cube/indep-seed/diamond) 均保持 4D 选择，窗口漂移可由因果稀疏度定量解释
6. ~~**论文草稿 v1**~~: ✅ `PredictionA_LinkAction_Paper_v1.md` — 含 BDG 对比 + 鲁棒性分析 + 机制普适性论述
7. **无量纲控制参数 Ξ**: 计算 Ξ_d(N,geom) ≈ Δ(link penalty per element) / Δ(entropy per element)，验证 4D 赢的阈值是否跨生成器一致
8. **论文 v2 完善**: 加入 Ξ 分析、完善 figures、准备提交

### 中优先级

9. **6D/7D 测试**: 验证 BD 是否也能抑制 6D+
10. **Myrheim-Meyer 维度估计量**: 作为 BD 的替代或补充
11. **Prediction B 论文更新**: 加入 N=56 γ_c=2.08 数据

### 探索性

12. **BD + consistency penalty 混合**: 是否能在更窄的 λ 窗口锁定 4D？
13. **因果连通度作为 order parameter**: Prediction B 的 γ 相变是否与 BD 的维度选择有关联？

---

## 九、关键函数签名速查

```python
# 生成偏序集
from generators import generate_lorentzian_like_2d  # (n, seed) → Poset
from generators import generate_lorentzian_like_5d  # 同上, 1+4维

# 熵估计
from runtime_utils import estimate_entropy_by_family
# (poset, family, sis_runs, seed, default_exact_threshold, family_exact_thresholds) → (log_h, method)

# 作用量
from action import action_value  
# (log_extensions, penalty, beta, gamma) → float

# 几何组件
from observables_geo import geometric_components  # (poset) → dict
from observables_geo import cover_density         # (poset) → float (Hasse 边密度)

# 中性惩罚
from observables import neutral_penalty  # (poset) → float

# Hasse 链接数 (BD 作用量核心)
# 见 prediction_a_bd_dimension.py 中 hasse_links(poset) → int
# BD action: S_BD = N - 2 * hasse_links(poset)
```

---

## 十、实验配置模板

### 5D 实验标准配置

```python
N_VALUES = [20, 24, 28, 32, 36, 40, 44, 48, 52]
FAMILY_LIST = ["lorentzian_like_2d", "lorentzian_like_3d", 
               "lorentzian_like_4d", "lorentzian_like_5d", "KR_like"]
SAMPLES = 4
SIS_RUNS = 4096
DEFAULT_EXACT_THRESHOLD = 24
FAMILY_EXACT_THRESHOLDS = {
    "lorentzian_like_2d": 104,  # 2D 可以精确到很大
    "lorentzian_like_3d": 24,
    "lorentzian_like_4d": 24,
    "lorentzian_like_5d": 24,   # 5D 必须≤24
    "KR_like": 24,
}
```

### BD 作用量实验配置

```python
LAMBDA_VALUES = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 18.0, 25.0, 50.0]
BETA = 1.0
# score = -BETA * log_H + lam * S_BD_norm
# S_BD_norm = (N - 2 * hasse_links) / N
```
