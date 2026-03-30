# Poset Phase 项目进度总览

> **生成日期**: 2026-03-27  
> **覆盖范围**: 项目全生命周期，截至 v5.0.0 + Carlip 后重建  
> **GitHub**: `github.com/unicome37/poset_phase` | **Branch**: main  
> **Zenodo 代码 DOI**: `10.5281/zenodo.18980657`

## 🔄 状态快照更新（2026-03-30）

- **当前核验基线**：`HEAD = eaea483`（`main` 与 `origin/main` 对齐）
- **关键新增完成项**：
  - F1 family-pressure 正式检验通过（Hard fail = NO）
  - F2 margin-aware turn-on refit 完成（20 seeds × 120 reps），onset 三档统一前推至 `N=10`
  - F3 highN 三路任务（FLRW/Schwarzschild/de Sitter）5/5 完成
  - FLRW `kappa=1.0` 在 highN 为边界敏感（fail_ratio=0.3），未触发 hard fail
  - P0 FLRW metric-faithful PhaseA 完成（10/10 seeds）：metric `fail_ratio = 11/60 = 0.183`，proxy `0/60`；`κ=0.3` 全通，`κ=1.0` 为边界敏感但未触发 hard fail
- **当前推荐总口径**：`background-dependent robustness`（避免“统一弱曲率稳健”）
- **最新总览建议优先参考**：`进展.md`（2026-03-30 重写版）

---

## 一、项目核心

偏序集相变（Poset Phase Transition）项目基于"结构存在论"理论框架，用数值实验验证核心预测（A/B/C/D/E），探讨**存在性筛选**（ESD）如何通过熵-作用量竞争，在偏序集空间中自然涌现出物理时空结构。

### 关键警告：Carlip 转折

**2026-03-17 前后，项目经历了关键转折点。** Steve Carlip 教授的直接反馈击穿了早期框架的三个核心薄弱环节（详见第六节），导致：

1. **Prediction A CQG 投稿 (CQG-115402) 被 desk reject**
2. **v5.0.0 F7 泛函方案被 retire** — logH 在 17 家族空间 N≥28 失败
3. **研究方向全面转向** → "两层筛选"架构（S_BD 线性准入 + S_MD Mahalanobis 身份认同）

当前论文定位：从 logH/F7 体系退出，转向纯因果几何可观测量 $(d_{\text{eff}}, C_1/C_0, \text{width})$ 构建零参数身份选择器。

---

## 二、各预测验证总状态

| 预测 | 核心问题 | 数值验证 | 投稿/发布状态 | 当前定位 |
|------|---------|---------|-------------|---------|
| **A — 维度选择** | Link action 能否自然选出 3+1 维？ | ✅ 数值完成 | ~~CQG-115402~~ **desk reject** | Ξ 参数理论成立，但 logH 基础被质疑 |
| **B — γ 相变阈值** | 是否存在有界 γ_c 实现几何相变？ | ✅ 完成 | Zenodo preprint DOI:10.5281/zenodo.19048146 | 最硬主线，但 logH 物理意义待辩护 |
| **C — 层级预测熵** | 层级深度能否预测组合熵？ | ✅ 完成 (v4.0) | Zenodo preprint DOI:10.5281/zenodo.19048405 | 机制解释方向，9 实验准因果塔 |
| **D — 两层筛选** | Mahalanobis 零参数选择 Lor4D？ | ✅ 完成 | 撰写中 | **Carlip 后核心方向**，25 族 N=10-1024 |
| **E — F7 泛函 (v5.0.0)** | 统一结构代价泛函？ | ⚠️ **已 retire** | ~~Physica A~~ **desk reject** | logH 主导项在大 N 失效 |

---

## 三、Prediction A — 维度选择（数值完成，投稿被拒）

### 3.1 发现链（Carlip 前）

| Phase | 内容 | Commit | 核心结果 |
|-------|------|--------|---------|
| 1 | BD 链接作用量发现 | a8e79c7→e4d237d | Link action λ=6-8 选出 4D, N=20-68 全票 |
| 2 | BDG 对比验证 | 6514e95→58501cb | BDG d=4 全选 5D — 高阶修正破坏选择 |
| 3 | 生成器鲁棒性 | b34b8be | 3 种生成器均保持 4D 选择 |
| 4 | Ξ 参数发现 | 41d866b→5655002 | Ξ₄→₅≈10 跨生成器 (CV=13.9%) |
| 5 | 大 N 闭环 | 36a8f45 | N=80-112: Ξ 不漂移 (CV=17.7%) |
| 6 | 解析推导 | 8c36896 | 闭合公式 Ξ₄→₅=11.8 (实测 11.3) |
| 7 | d≥6 盲预测 | afd6bfa | 5→6: 64→46 (29%), 6→7: 51→38 (26%) |
| 8 | LaTeX + Zenodo | 4c8e06d→f70ff5b | 31pp, 16 fig, DOI:10.5281/zenodo.19068927 |
| 9 | 第一性原理推导 | 9c5481e→4ca1397 | 零参数闭合 ε=π/48, Ξ₄→₅ 误差 0.4% |
| 10 | AI 审阅 + CQG 投稿 | b6379b2→1ee2654 | v2.6, 33pp, ~~CQG-115402~~ **desk reject** |

### 3.2 零自由参数闭合公式（理论仍成立）

$$\frac{\log H}{N} = \left(1-\frac{\pi}{48}\right)(1-p_d)(\ln N-1) + \frac{\pi}{24}\left(1-e^{-(N-2)\mathbb{E}[V_A]}\right)$$

- 输入：闵可夫斯基几何常数 $p_d$, $\kappa_d$, $\mathbb{E}[V_A]$
- 输出：Ξ₄→₅ = 11.31（观测 11.35，误差 0.4%）
- **注意**：此公式基于 logH，而 logH 的物理地位已被 Carlip 质疑（见第六节）

### 3.3 物理解释

- **为什么是 3+1 维**：5D 因果极稀疏（p₅=0.11, 94.7% 链接），轻微耦合即抑制；4D 仍有足够熵击败 3D
- **不对称屏障**：Ξ₄→₅≈10 >> Ξ₃→₄≈3，4D 上方硬屏障，下方无
- **d≥6 更强**：Ξ₅→₆≈46, Ξ₆→₇≈38，4→5 是最弱屏障

### 3.4 CQG 投稿与 Desk Rejection

- **Manuscript ID**: CQG-115402
- **投稿日期**: 2026-03-17
- **结果**: **Desk Reject**
- **直接原因**: Carlip 教授指出三大核心问题（见第六节）
- **Zenodo 存档**: DOI `10.5281/zenodo.19079466` (v2.6, 33p)

### 3.5 Pred A 遗产

虽然论文被拒，Ξ 参数理论和零参数闭合推导本身作为**数学结果**仍然成立。问题不在数值计算，而在 logH 这个基础物理量是否等同于 CST 物理熵。该问题已在 Carlip 后重建中被直面。

---

## 四、Prediction B — 有界几何相变

### 4.1 核心结果

A2 作用量（neutral + geometric penalty）下，临界几何耦合 γ_c 在 N=12-44 保持 O(1)：

| N | 12 | 14 | 16 | 20 | 24 | 28 | 32 | 36 | 40 | 44 |
|---|----|----|----|----|----|----|----|----|----|----|
| γ_c | 0.74 | 0.78 | 0.69 | 0.15 | 0.26 | 1.00 | 0.69 | 0.57 | 0.25 | 0.15 |

### 4.2 消融与非循环性

- **Backbone**（必需）：width-height balance + dimension proxy
- **Shell**（增强）：其余 5 个几何子项
- **非循环替代**：dimension consistency $g_{\text{con}}$ 不指定目标维度，相变保持

### 4.3 候选池冻结

- **Tier-1（5 家族）**：Lor2D/3D/4D, KR_like, transitive_percolation — 评估核心声明
- **Tier-2（7 家族）**：random_layered 系列 — 压力测试, 暴露 action 判别天花板
- $aw/\sqrt{N}$ 标度律清晰分离两 tier (Mann-Whitney p=1.94×10⁻²⁸)

### 4.4 状态

- Zenodo preprint 已发布
- Near-wall N=52-56 mixed 结果存在 near-degeneracy 但无 crossing 确认

---

## 五、Prediction C — 层级深度预测熵 (准因果)

### 5.1 九实验证据塔

| # | 实验 | 关键结果 |
|---|------|---------|
| 1 | Stratified Fisher z | \|r\|~0.35-0.54 |
| 2 | 单边因果干预 | Cohen's d=1.05, p<10⁻³² |
| 3 | 安慰剂控制跨维 | d=1.40-1.83, p<10⁻¹³³ |
| 4 | 反向干预 | 100% 方向一致, d=2.68 |
| 5 | k-family 剂量-反应 | all slopes negative, p<10⁻²⁹ |
| 6 | 边密度普遍性 | 无相变, all slopes < 0 |
| 7 | N-scaling | \|slope\|∝N⁰·⁷⁰ |
| 8 | 解析下界 | logH=k·log(m!), 验证至 10⁻¹⁵ |
| 9 | 大 N SIS | N=16-36, Cohen's d=0.69-1.29 |

### 5.2 B↔A↔C 交叉

$$\text{Pred C: 更多层级→更少熵} \xrightarrow{\text{跨维}} \text{Pred A: 高维→更少因果→更高熵} \xrightarrow{\text{Link action}} \text{自然选出 3+1 维}$$

验证：density-entropy 反相关 r = -0.988 to -0.997

---

## 六、⚡ Carlip 转折点 — 项目分水岭

### 6.1 Carlip 邮件（2026-03-17）

在向 Carlip 教授发送论文草稿后收到直接批评，三条命中要害：

| 编号 | Carlip 批评 | 评估 | 冲击 |
|------|-----------|------|------|
| **C1** | logH（线性扩展数）≠ 物理熵 | **Carlip 正确** | 整个 F7/logH 体系的物理基础被动摇 |
| **C2** | 7 家族 = cherry-picking；缺 2-layer/4-layer | **部分成立** | KR_2layer 在 F7 下与 Lor5D 几乎不可区分 |
| **C3** | 错引文献（Carlip/Surya/Loomis 归因全错）；缺 Dhar(1978)/Prömel(2001) | **完全成立** | AI 幻觉导致的文献错误被直接指出 |

**Carlip 原话（节选）**：

> *"Your definition of entropy in terms of linear extensions has no apparent relation to the quantity relevant to physics."*

> *"If you are not familiar enough with the literature to describe basic results, or to recognize when some AI is producing slop, you're not likely to find much of interest."*

### 6.2 直接后果

1. **CQG-115402 desk reject** — Prediction A 论文
2. **F7 泛函方案 retire** — logH 主导项在 17 家族 N≥28 失败，Lor4D 排名跌至 #11/17
3. **Physica A 投稿 desk reject** — v5.0.0 F7 泛函论文
4. **文献审计启动** — 49 句强动词引用句全部核实，Dhar/Prömel 补入

### 6.3 F7 泛函失败根因（N=20 诊断）

| Family | logH | wall | F7 | 说明 |
|--------|------|------|-----|------|
| Lor5D | 33.5 | 0.00 | 31.5 | wall=0, 纯靠低 Σ_hist |
| KR_2layer | 33.8 | 0.00 | 32.3 | **与 Lor5D 几乎不可区分** |
| Lor4D | 32.6 | 2.15 | 32.4 | 微弱 wall 帮助 |
| KR_like | 27.0 | 15.84 | 43.3 | wall 惩罚成功 |

**根因**：logH 以 O(N) 增长而 wall 以 O(N⁻⁰·⁵) 衰减 → 大 N 时 layered 结构自动胜出。

### 6.4 教训

1. **AI 生成的文献引用必须逐条人工核实** — AI 幻觉是学术致命伤
2. **logH 是组合量，不能直接等同物理熵** — 需要 Rideout-Sorkin 路径积分或 BDG 替代
3. **样本空间必须包含数学上真正占优的结构** — 2-layer/4-layer 是 KR 之后的主导相
4. **内部自洽 ≠ 外部接口通畅** — "齿轮完美但插头不匹配"

---

## 七、Carlip 后重建 — 两层筛选架构（当前方向）

### 7.1 核心转变

| 维度 | Carlip 前 | Carlip 后 |
|------|----------|----------|
| 核心变量 | logH（无物理推导） | $d_{\text{eff}}, C_1/C_0, \text{width}$（因果几何） |
| 样本空间 | 7-17 家族 | 25 家族（17 标准 + 8 对抗） |
| 选择器架构 | F7 泛函（单层） | 两层筛选: S_BD + S_MD |
| 自由参数 | α,β,γ + logH 权重 | **零自由参数**（Mahalanobis） |
| N 范围 | N=10-112 | N=10-1024 |

### 7.2 两层筛选架构

**Layer 1 — 线性准入 (S_BD)**:
$$S_{\text{BD}} = \mathbf{c}^\top \Delta\mathbf{C}$$
Benincasa-Dowker 作用量，interval count 空间中的线性超平面。淘汰曲率不对的结构，但**不能**唯一选择 Lor4D（单独排名仅 14/17 at N=128）。

**Layer 2 — 二次身份认同 (S_MD)**:
$$S_{\text{MD}} = (\mathbf{I}(P) - \boldsymbol{\mu}(N))^\top \boldsymbol{\Sigma}(N)^{-1} (\mathbf{I}(P) - \boldsymbol{\mu}(N))$$
Mahalanobis 距离，以 Lor4D 参考流形 $\boldsymbol{\mu}(N)$ 为中心，零自由参数。当前安全口径是：在 fixed-reference F2 protocol 下，Lor4D 从 `N≥10` 起稳定排名第 1，并在更大尺度和 25 族库中持续保持该身份优势。

**关键**：两层几乎正交（Pearson r≈0），交集 = {Lor4D} only。

### 7.3 高置信结论

| 结论 | 来源 | 关键数据 |
|------|------|---------|
| Mahalanobis 零参数身份选择 | `f2_turnon_margin_refit` + `mahalanobis_lsd_test` | 当前安全口径：`N≥10`（fixed-reference F2）；旧 seed/CV 记录仅作历史诊断 |
| Basin deepening | basin_deepening_experiment + F2 onset | gap: 0.308 (N=10, F2) → 1.93×10⁸ (N=1024); V_eff∝N⁻¹·⁵⁷ |
| 参考流形坍缩 | mu_trajectory_theory | det(Σ)∝N⁻³·³⁸, σ²∝N⁻¹ (CLT) |
| d_eff→4 收敛 | cstar_wstar_first_principles | μ_d(∞)=3.967≈4 |
| 对抗鲁棒性 | expanded_family_robustness | 25 家族既有 tested range 内 Lor4D 持续 #1；主稿外口径以 F2 `N≥10` 为准 |
| LSD-Well N-adapted | lsd_well_n_adapted | 三种模式 100%, N=16-256 |

### 7.4 主要风险

| 风险 | 状态 | 缓解 |
|------|------|------|
| 小N turn-on 边界不稳（旧口径 N≥14） | ✅ 已解决 | F2 margin-aware refit（20 seeds × 120 reps，独立参考集）：N≥10 全部 20/20 rank#1，manuscript_safe 亦在 N=10 成立 |
| 梯度桥过度声明 | ✅ 已降级 | 仅作辅助，不声称逐点梯度同一 |
| "零参数"措辞精确度 | ⚠️ 需谨慎 | 明确 μ(N)/Σ(N) 估计协议 |

---

## 八、v5.0.0 F7 泛函 — 已 Retire

### 8.1 原设计

- Interval-richness 准入墙 + Sigmoid 阈值 + Pareto 参数搜索（915 可行配置）
- 配合论文投稿 Physica A → **desk reject**

### 8.2 Anti-overfitting 测试（7 项全通过，但前提已失效）

排列检验 p<0.001 / 随机泛函 MC 0.02% / 消融非冗余 / 4-fold CV ΔA=+0.009 / 新种子泛化 ΔA=-0.062

### 8.3 Retire 原因

F7 内部逻辑自洽，但核心变量 logH 在 17 家族空间 N≥28 不能区分 Lor4D 和 KR_2layer/layered，被 Carlip C1 批评命中。Anti-overfitting 通过不代表物理基础成立。

### 8.4 遗产继承

- 约束动力学（微正则 swap, cluster moves, 模拟退火）→ 可复用于新架构
- 谱恢复实验（W₁ Wasserstein distance）→ 可能整合入两层筛选

---

## 九、理论框架 — 深层分析

### 9.1 LSD-Well 本质（仍成立）

不是 Action，而是**有效势**——三重筛选（可存在/可延续/可显影）在 Lor4D 附近的局域二阶展开。这一解释在两层筛选架构下得到形式化：S_BD = 准入，S_MD = 身份认同井。

### 9.2 三层统一架构

| 层面 | 表现 | 语言 |
|------|------|------|
| **本体层** | 幸存条件 | 可存在 / 可延续 / 可显影 |
| **中层机制** | 准入势 | wall / admissibility / 门槛 |
| **观察层** | 固定点距离 | 距离 / 参考流形 / 吸引域 |

### 9.3 核心命题

> 宇宙的根本不是"有什么"，而是"什么能长期低失真地幸存并显影成现实"。Lor4D 可能正是这样一种主幸存几何。

---

## 十、证据分层体系

| 层级 | 用途 | 示例 |
|------|------|------|
| **Confirmatory** | 冻结主线，可直接引用 | frozen_exact, medium_exact_scan |
| **Exploratory** | 发现/诊断/方法扩展 | prediction_c_intervention, BAC bridge |
| **Control** | 对照组验证 | 17-family + KR_2layer/4layer |
| **Carlip Response** | Carlip 后重建实验 | F7 17-family, LSD-Well, Mahalanobis, S_BD/S_MD |

---

## 十一、发布资产（含状态更新）

| 资产 | DOI | 类型 | 当前状态 |
|------|-----|------|---------|
| 代码仓库 v4.0.0 | 10.5281/zenodo.18980657 | Software | ✅ 公开 |
| Prediction B | 10.5281/zenodo.19048146 | Preprint | ✅ 最硬主线 |
| Prediction A v1.5 | 10.5281/zenodo.19068927 | Preprint | ⚠️ CQG desk reject, Zenodo 存档有效 |
| Prediction C v2 | 10.5281/zenodo.19048405 | Preprint | ✅ 准因果证据塔 |
| CQG 投稿稿 v2.6 | 10.5281/zenodo.19079466 | Archive | ⚠️ Desk reject, 存档仍可引用 |

**优先权状态**：所有关键文本和代码已有国际通行 Zenodo DOI，优先权已锁定，不依赖期刊。

---

## 十二、待完成工作与 2-3 周工作计划

### 当前主线（两层筛选论文）

| 周次 | 任务 | 交付物 |
|------|------|--------|
| Week 1 | ✅ 锁定声明边界 + 解决小 N 问题 | **N≥10 turn-on 边界确认**（F2: 20×120，三档 onset 一致） |
| Week 1 | ✅ N-boundary turn-on 实验完成 | `mahalanobis_n_boundary_turnon.md` + `fig_margin_vs_n.png` |
| Week 2 | 参考流形从"拟合"变为"理论对象" | Methods 节完整化 |
| Week 2 | 特征极小性和消融故事 | 消融表 + Limitations 节 |
| Week 2 | 梯度桥降级为辅助 | Discussion 子节 |
| Week 3 | Letter 定稿（CQG Letters / PRL 风格） | 1 fig, 1 table, self-contained |
| Week 3 | Submission packaging | arXiv gr-qc 首发目标 |

### 投稿策略（传播先行）

> **传播先行，期刊并行，CQG 可投但不再卡主线。**

1. **arXiv 首发**（最高优先）：主口径 Prediction B + 两层筛选，首选 gr-qc
2. **Zenodo 已锁优先权**：不把公开节奏押在单一审稿流程上
3. **期刊并行**：CQG Letters / PRL / Found. Phys. 可选

### 技术遗留

| 项目 | 优先级 | 说明 |
|------|--------|------|
| ε=κ₄/2 严格物理推导 | 中 | 当前为 numerical identification |
| d=6 零参数闭合验证 | 中 | Ξ₅→₆ 精度待测 |
| 弯曲时空撒点 | 低 | de Sitter / Schwarzschild 下 Ξ 稳定性 |
| Action 判别天花板突破 | 低 | 将 aw/√N 嵌入 action |

### Carlip 指引的三条前进路径

- **A: 重新定义度量** — S_BDG 替代 logH，或 Rideout-Sorkin sequential growth 概率
- **B: ρ 条件对比** — 固定 ρ=f₂(4)≈0.05，该密度下采样所有偏序集与 Lor4D 对比
- **C: 信息几何方法** — (ρ, Σ_hist, d_eff, ...) 参数化偏序集，展示 Lor4D 独特性

---

## 十三、技术栈

| 项 | 值 |
|---|---|
| Python | 3.14.2 |
| 项目路径 | `d:\Kiro\理论体系\poset_phase\` |
| 生成器 | 2D-7D Lorentzian + KR + layered 系列（25 家族含 8 对抗） |
| 熵计算 | 精确 DP (N≤44) + SIS 采样 (大 N) |
| 特征空间 | $d_{\text{eff}}, C_1/C_0, \text{width\_ratio}$ |
| LaTeX | MiKTeX (pdflatex) |
| 版本管理 | Git → GitHub → Zenodo |

---

## 十四、项目时间线概览

```
2025 Q3-Q4: Prediction B 有界 γ_c (exact N=10-44)
2026 Q1 前: Prediction C 准因果证据塔 (v4.0)
2026-01~02: Prediction A 维度选择链完成 (Phase 1-9)
2026-03-17: CQG-115402 投稿 + Carlip 邮件
2026-03-17: ⚡ CQG desk reject + F7 retire
2026-03-18~26: Carlip 后重建 → 两层筛选架构
2026-03-26: LSD-Well/Mahalanobis 25 族实验完成
2026-03-27: 当前 — 2-3 周 Letter 撰写阶段
```

---

*本文档由 9 份核心项目文件 + Carlip 后重建文档整合生成，反映 2026-03-27 真实项目状态。*
