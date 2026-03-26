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

## Rule

简单规则：

- 冻结配置直接产出的结果，归 `confirmatory`
- 为解释当前失败点或扩展候选空间新增的试验，归 `exploratory`
- 对照组实验及图表，归 `control group`
