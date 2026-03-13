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

## Rule

简单规则：

- 冻结配置直接产出的结果，归 `confirmatory`
- 为解释当前失败点或扩展候选空间新增的试验，归 `exploratory`
