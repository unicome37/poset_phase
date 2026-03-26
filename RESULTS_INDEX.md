# Results Index

`poset_phase` 当前输出按证据等级分为两层：

## Confirmatory

这些结果对应已经冻结的主实验设定，可直接进入当前正文主证据链。

- `outputs_confirmatory/frozen_exact`
  - 来源配置：`config_frozen_exact.yaml`
  - 作用：小 `N` 精确主线复现，验证 `lorentzian_like_2d` 相对 `KR_like` 的临界切换
  - 关键文件：
    - `summary.csv`
    - `bootstrap_summary.csv`
    - `gamma_c_report.csv`

- `outputs_confirmatory/frozen_cg`
  - 来源配置：`config_frozen_cg.yaml`
  - 作用：冻结版粗粒化稳定性复现，检验 `I_cg` 是否提供附加区分力
  - 关键文件：
    - `cg_summary.csv`
    - `cg_rank_summary.csv`

- `outputs_confirmatory/medium_exact`
  - 来源配置：`config_confirmatory_medium_exact.yaml`
  - 作用：将确认性精确主线扩展到 `N=20~40`，检验 `lorentzian_like_2d` 相对 `KR_like` 的临界切换是否在更大 `N` 上继续存在
  - 关键文件：
    - `summary.csv`
    - `bootstrap_summary.csv`
    - `gamma_c_report.csv`

- `outputs_confirmatory/medium_exact_scan`
  - 来源配置：`config_confirmatory_medium_exact_scan.yaml`
  - 作用：对 `Lor2D / Lor3D / KR` 做更完整的确认性精确扫描，将 `lorentzian_like_2d` 相对 `KR_like` 的交点追踪推进到 `N=44`
  - 关键文件：
    - `summary.csv`
    - `bootstrap_summary.csv`
    - `gamma_c_report.csv`

- `outputs_confirmatory/exact_timing`
  - 来源配置：`config_exact_timing.yaml`
  - 作用：基准测试精确线性延拓数计算的性能窗口，验证 Lorentzian-like 2D 在更大 `N` 上仍可直接使用 exact
  - 关键文件：
    - `exact_timing_benchmark.csv`
    - `exact_timing_summary.csv`

## Exploratory

这些结果用于现象发现、诊断和方法扩展，不直接作为主结论依据。

- `outputs_exploratory/smallN_exact`
  - 早期小 `N` 精确主实验结果

- `outputs_exploratory/medium`
  - 中规模参数扫描与相图

- `outputs_exploratory/gamma_scan`
  - 更宽 `gamma` 扫描

- `outputs_exploratory/cg_smoke`
  - 粗粒化稳定性烟雾测试

- `outputs_exploratory/cg_pilot`
  - 粗粒化稳定性 pilot

- `outputs_exploratory/cg_dimension`
  - 维度扫描与维度偏置诊断

- `outputs_exploratory/controlled_dimension`
  - 受控窗口下的维度比较

- `outputs_exploratory/smoke`
  - 最小链路烟雾测试

- `outputs_exploratory/prediction_c_comprehensive`
  - 来源配置：`config_prediction_c_comprehensive.yaml`
  - 来源脚本：`prediction_c_comprehensive.py`
  - 作用：Prediction C 三层综合验证——因果层级整合 (HII) 与 log_H 的偏相关 (Tier 1)、配对Δ分析 (Tier 2)、CG 稳定性关联 (Tier 3)
  - 关键文件：
    - `tier1_overall_summary.csv` — 总体偏相关
    - `tier1_by_family.csv` — 8 族分解
    - `tier1_components.csv` — 组件级偏相关
    - `tier2_pairwise_summary.csv` — 34 配对Δ相关汇总
    - `tier3_cg_linkage_summary.csv` — CG 稳定性关联
  - 精确表述文档：`PredictionC_精确表述与验证报告.md`
  - 核心结论：
    - 跨族 HII→log_H: r=-0.836 (p<0.001, 34 配对)
    - HII→CG_switch_rate: r=-0.835 (p<0.001, 68 样本)
    - Simpson's Paradox 根因为尺度变量 N，控制后方向翻转

- `outputs_exploratory/prediction_bac_bridge`
  - 来源配置：`config_prediction_bac_bridge.yaml`
  - 来源脚本：`prediction_bac_bridge.py`
  - 作用：桥接 Prediction B、A、C，检验 B 的有界相变是否沿 C 的层级深度机制展开，以及 A 是否延续同一机制
  - 关键文件：
    - `prediction_bac_relation_summary.csv`
    - `prediction_bac_bridge_report.md`
    - `prediction_a_family_depth_summary.csv`
    - `prediction_b_family_depth_summary.csv`
  - 核心结论：
    - B↔C：`Lor2D` 相对 `KR_like` 的层级优势同时对应更低惩罚和更低 `log_H`，形成有限 `gamma_c`
    - A↔C：`Lor4D` 在 `N=20..72` 上从未高于 `Lor2D` 或 `Lor3D` 的 HII，因此 A 的 4D 优势并不等同于“更深层级”

## Control Group (added post hoc)

- `outputs/` (full 17-family experiment)
  - 来源配置：`config.yaml`
  - 来源脚本：`experiment.py`
  - 作用：包含全部 17 个 family（含 KR_2layer/KR_4layer 对照组）的完整实验，N=20/40/60/80，16 samples × 5 gammas × 3 actions = 16,320 行
  - 关键文件：
    - `raw_samples.csv` — 完整原始数据
    - `summary.csv` — 分组汇总统计
    - `bootstrap_summary.csv` — Bootstrap CI

- `outputs_control/` (dedicated control experiment)
  - 来源脚本：`control_group_experiment.py`
  - 作用：9 族子集 × N=10-20 × 8 samples × 5 gammas，精确 logH 对比
  - 关键结论：
    - KR_2layer 拥有最高 logH 但最弱几何结构 → γ>0.4 时被有效筛除（排名 16-17/17）
    - KR_4layer 与 KR_like 表现接近 → 不产生误导信号

- `manuscript_figures/fig_control_*.{png,pdf}` (3 figures)
  - 来源脚本：`plot_control_group_figures.py`
  - **Fig A** `fig_control_ranking_heatmap`: 17 族 × 5 γ 排名热图（A2, N=80）
    - KR_2layer 在 γ≤0.2 稳定排名 #17；γ=0.4 跃升至 #4（logH 主导）
    - 证明 geometric penalty 在适当 γ 下有效区分几何良/劣
  - **Fig B** `fig_control_separation`: 关键 family score_norm 轨迹（N=60, 80）
    - Lor2D 稳定在 +1.0（不受 γ 影响）；KR_2layer 从 -1.9 线性增长至 +3.6
    - 分离点 γ≈0.3 处 KR_2layer 超越 KR_like → 标志 logH 主导的伪信号转折
  - **Fig C** `fig_entropy_geometry_scatter`: logH vs Π_geo 散点（N=80, A2）
    - 各 family 类别在熵-几何平面上清晰分离
    - Lor2D 左下角（紧致）、KR_2layer 右上角（松散）、Layered 变体聚集中间

### KR_2layer 大 γ 翻转的解释

KR_2layer 在 γ≥0.4 时排名跃升**不代表**其具有几何优势。机制如下：
- `score = β·logH - γ·penalty`（A2 模式下 penalty = penalty_geometric + penalty_neutral）
- KR_2layer 的 logH=229.7 是所有 17 族最高（因 2-layer 结构几乎不施加偏序约束）
- 当 γ 增大时，`-γ·penalty` 项的绝对值增长受限于 penalty 有限（56.4），但 `logH` 贡献恒定
- **Normalization 机制**：score_norm 按 n 组做 robust z-score。大 γ 下所有 family 的 score 绝对值拉开，KR_2layer 的 logH 优势被放大
- **论文不应在 γ > γ_c 区域下结论**——γ_c 定义为 Lor2D 超越 KR_like 的阈值（~0.2-0.4），超过此点结论稳健；但继续增大 γ 会引入 logH 主导的伪信号

## Carlip Critique Response (2026-03-26)

针对 Carlip (CQG desk rejection) 三条批评的系统回应。

### 文件

- `carlip_f7_17family_test.py` — F7 泛函在全部 17 家族上的测试脚本
- `outputs_carlip/carlip_response_analysis.md` — 完整 8 节回应文档（§1-§8）
- `outputs_carlip/f7_17family_test.md` — F7 17 家族测试详细输出
- `outputs_carlip/诊断总结_F7_17家族.md` — 中文诊断总结

### Carlip 三条批评摘要

| 编号 | 批评内容 | 我方评估 |
|------|---------|---------|
| C1 | logH (线性扩展数) ≠ 物理熵 | **Carlip 是对的。** logH 计数标签排列数，非 CST 路径积分中的物理测度 |
| C2 | 7 家族 = cherry-picking | **部分成立。** KR_2layer 与 Lor5D 在 F7 下几乎不可区分 (R≈0, wall=0) |
| C3 | 缺失 Dhar/Prömel 文献 | **已补充。** Dhar 1978 的 ρ ≡ 我们的 f₂，Prömel 2001 的层结构相变恰是 KR_2layer/4layer |

### F7 17 家族测试关键结果

| N | Lor4D 排名 | 状态 | 说明 |
|---|-----------|------|------|
| 16 | #1 | ✅ | 击败所有非 Lorentzian |
| 20 | #3 | ≈ | 与 KR_2layer 持平 |
| 28 | #8 | ❌ | 输给 6 个 layered 家族 |
| 36 | #11 | ❌ | 输给 9 个非 Lorentzian |

**根因**：logH 以 O(N) 增长且 Lor4D/5D 的 logH 最高，而 wall 项以 O(N^{-0.5}) 衰减。大 N 时 layered 结构（低 logH + 满 wall）自动胜出。

### 核心诊断（N=20 组件分解）

| Family | logH | wall | F7 | logH rank → F7 rank |
|--------|------|------|-----|---------------------|
| Lor5D | 33.5 | 0.00 | 31.5 | 15 → 1 (wall=0, 纯靠低 Σ_hist) |
| KR_2layer | 33.8 | 0.00 | 32.3 | 16 → 2 (**与 Lor5D 几乎不可区分**) |
| Lor4D | 32.6 | 2.15 | 32.4 | 14 → 3 (微弱 wall 帮助) |
| KR_like | 27.0 | 15.84 | 43.3 | 11 → 17 (wall 惩罚成功) |

### 前进路径

- **A: 重新定义度量** — 用 S_BDG 替代 logH，或用 Rideout-Sorkin sequential growth 概率
- **B: ρ 条件对比** — 固定 ρ = f₂(4) ≈ 0.05，在该密度下采样所有偏序集，与 Lor4D 对比
- **C: 信息几何方法** — 用 (ρ, Σ_hist, d_eff, ...) 参数化偏序集，展示 Lor4D 在几何不变量上的独特性

---

## Well Center Physics & LSD-W2 (2026-03-26)

从经验到第一原理：将 LSD-Well selector 的 well center 物理化分析。

### 文件

- `carlip_well_center_physics.py` — d* 推导 + c*(N) 解析公式 + w*(N) 钻石几何
- `carlip_lsd_17family.py` — LSD 全 17 族判别测试
- `carlip_lsd_well_17family.py` — LSD-Well 版本测试
- `carlip_link_action_17family.py` — BDG link action 替代 logH 方案
- `carlip_info_geometry_17family.py` — 信息几何方法
- `outputs_carlip/well_center_physics.md` — well center 物理推导报告
- `outputs_carlip/lsd_17family_test.md` — LSD 判别结果
- `outputs_carlip/lsd_well_17family_test.md` — LSD-Well 测试结果
- `outputs_carlip/link_action_17family_test.md` — link action 测试结果
- `outputs_carlip/info_geometry_17family.md` — 信息几何分析
- `outputs_carlip/综合诊断_两轮实验.md` — 五轮实验累积诊断

### LSD-W2 Well Center 三参数评估

| 参数 | 值 (LSD-W2) | 可从第一原理推导? | N-stable? | 物理来源 |
|------|------------|----------------|-----------|---------|
| d* | 3.93 → 4.0 | ✅ 完全可导 | ✅ (→4) | 时空维数 d=4 |
| c* (C₁/C₀) | 0.213 | ✅ 有解析公式 | ❌ 漂移 63% | Beta(2,2) 积分 (需 finite-N 修正) |
| w* (width) | 0.408 | ⚠️ 部分可导 | ❌ 漂移 34% | 因果钻石空间截面 |

### 关键洞察

1. **d\*=4 是物理输入**，不是拟合参数 — 这是目标时空维数，完全可导
2. **c\*(N) 有封闭形式**：E[Nu·e^{-Nu}]/E[e^{-Nu}], u ~ Beta(2,2) — 是导出函数，非自由参数
3. **w\*(N) 来自因果钻石几何**，但精确公式需要完整 Alexandrov interval 体积积分
4. LSD-W2 的**真正自由参数只有 α, β, γ（三个相对权重）**，well 中心全部可推导

### LSD-W2 vs F10 的认识论质变

| 维度 | F10 | LSD-W2 |
|------|-----|--------|
| 核心变量 | logH (无物理推导) | d_eff, C₁/C₀, width (因果几何) |
| 问题类型 | 核心变量本身站不住 | 核心变量站得住，well center 需 N-adapted 化 |
| 自由参数 | α,β,γ + logH 权重 | 仅 α,β,γ (well center 可推导) |
| Carlip C1 | ❌ 无法回应 | ✅ 不依赖 logH |

### N-adapted 升级方向

正式形式：
```
F_LSD(N) = α·[d_eff − 4]² + β·[C₁/C₀ − c*(N)]² + γ·[w − w*(N)]²
```
- d* = 4 固定
- c*(N) 由 finite-diamond interval-volume law 给出
- w*(N) 由 finite-diamond 横截面几何给出

> **well center 的物理分析表明，LSD-W2 的成功并非纯经验现象。d_eff 的中心值可完全由第一原理给出，并收敛到物理时空维数 4；C₁/C₀ 的中心值原则上也有封闭解析表达，但无限因果钻石的 Beta(2,2) 近似在当前 finite-N 范围内严重失效，必须引入 N-适配的有限体积修正；width 中心则仅部分可导，并同样表现出强烈的 N-依赖。因此，真正的下一代 selector 不应是常数中心井，而应是一个 N-适配的 Lorentzian 结构判别器，其唯一真正自由参数仅是 α,β,γ 三个相对权重。**

### N-adapted 实验结果 (Commit 2822dc1)

- `carlip_lsd_well_n_adapted.py` — 三种模式 N-adapted LSD-Well 实验
- `outputs_carlip/lsd_well_n_adapted.md` — 完整结果报告
- `_lsd_well_mode_analysis.py` — 三种模式一致性深度分析
- `outputs_carlip/lsd_well_n_adapted_why_identical.md` — 为什么三种模式结果一致

**三种模式全部达到 100% 判别率，Lor4D #1/17 at every N**：

| 模式 | 描述 | Beat% | 最差排名 |
|------|------|:-----:|:-------:|
| A (Oracle) | 同N Lor4D 中心值 | 100% | #1/17 |
| B (Extrapolated) | 幂律拟合 N≤36 → 外推 | 100% | #1/17 |
| C (Constant) | 固定 N=48 中心 | 100% | #1/17 |

**最佳权重**：α=0.5 (d_eff), β=1.0 (C₁/C₀), γ=5.0 (width)

**有限尺寸标度**：
```
c*(N) = 0.2485 − 2.33/N  (C₁/C₀ 中心)
w*(N) = 0.3255 + 3.80/N  (width 中心)
```

**为什么三种模式结果一致**：
1. d_eff 项使用固定 d*=4，不受模式影响 → 主导 Lor vs Layered/KR 的类别间分离
2. width 权重最大 (γ=5.0) 但 Lor4D 到 runner-up 的 margin 远大于 w* 变化量
3. c\*/w\* 在三种模式间的最大变化 ~0.10，而最小 margin = 0.043 → 二次项 O(Δx\*²) ≈ 10⁻⁴ 完全无法翻转排名
4. **核心结论**：LSD-Well 的判别力来自 Lor4D 在特征空间中的结构性孤立，而非 well center 的精确调优

---

## Rule

简单规则：

- 冻结配置直接产出的结果，归 `confirmatory`
- 为解释当前失败点或扩展候选空间新增的试验，归 `exploratory`
- 对照组实验及图表，归 `control group`
- Carlip 批评回应及文献补充，归 `carlip response`
- Well center 物理化分析及 LSD 相关实验，归 `well center / LSD`
