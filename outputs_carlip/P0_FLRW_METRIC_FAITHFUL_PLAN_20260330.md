# P0 执行细化：FLRW metric-faithful 边界验证（2026-03-30）

## 1. 目标

把 F3 中 `FLRW κ=1.0` 的争议点从“已有 lowN/highN 现象”推进为“可复核的跨尺度边界结论”，并与 C2 预注册 hard-fail 规则一致。

> 当前焦点不是再证明 de Sitter/Schwarzschild，而是把 FLRW 边界做成可审稿防守的主证据链。

---

## 2. 当前代码基线（已核对）

`curvature_backgrounds.py` 当前 FLRW 生成链：
- 体积权重：`a(t)=(1+κt)^(2/3)` 的 rejection sampling；
- 因果判据：`χ=(3/κ)[(1+κt2)^(1/3)-(1+κt1)^(1/3)]` 的 comoving horizon；
- 在 runner 中由 `falsify_c3_metric_runner.py` 调用。

该实现已是“几何驱动”的背景生成，不是纯玩具随机扰动；本轮先按该实现完成跨尺度封口。

---

## 3. 两阶段执行矩阵

## Phase A（边界复核）
- 配置：`config_falsify_c3_flrw_metricfaithful_phaseA.yaml`
- N 网格：`[512, 768, 1024]`
- κ 网格：`[0.3, 1.0]`
- seeds：`10`，reps：`15`
- 目标：复核 `κ=1.0` 在中高 N 区间是否维持“部分失败但未硬失败”。

## Phase A-Control（对照分支）
- 配置：`config_falsify_c3_flrw_proxy_control_phaseA.yaml`
- 模式：`flrw_mode=proxy_flatcausal`（保留 FLRW 抽样，平直化因果锥）
- 目标：分离“抽样权重效应”与“因果锥几何效应”。

## Phase B（跨尺度封口）
- 配置：`config_falsify_c3_flrw_metricfaithful_phaseB.yaml`
- N 网格：`[1024, 1536, 2048]`
- κ 网格：`[0.3, 1.0]`
- seeds：`10`，reps：`12`
- 目标：判断 fail_ratio 随 N 的走势（回落 / 平台 / 再恶化）。

---

## 4. 统一判定规则（与 C2 对齐）

- top-k 要求：`rank <= 2`
- hard-fail 单背景阈值：`fail_ratio >= 0.5`
- 本项目当前关键参数：`family=flrw, κ=1.0`

## 结论标签
- **PASS**：`fail_ratio < 0.5` 且跨 N 无恶化趋势
- **Boundary-sensitive PASS**：`fail_ratio < 0.5`，但局部 N 有重复跌出 top-2
- **HARD FAIL**：`fail_ratio >= 0.5`

---

## 5. 交付物（运行后）

- `outputs_carlip/falsify_c2_background_response_flrw_metricfaithful_phaseA.json`
- `outputs_carlip/falsify_c2_background_response_flrw_metricfaithful_phaseA.md`
- `outputs_carlip/falsify_c2_background_response_flrw_proxy_control_phaseA.json`
- `outputs_carlip/falsify_c2_background_response_flrw_proxy_control_phaseA.md`
- `outputs_carlip/falsify_c2_background_response_flrw_metricfaithful_phaseB.json`
- `outputs_carlip/falsify_c2_background_response_flrw_metricfaithful_phaseB.md`

建议汇总表字段：
- `N, κ, fail_cells, total_cells, fail_ratio, hard_fail, trend_note`

---

## 6. 文稿口径映射（P0 完成后直接替换）

若 `κ=1.0` 仍 `fail_ratio < 0.5`：
> FLRW at κ=1.0 remains boundary-sensitive: lowN threshold hit, while highN-to-extended-N runs show partial recovery without triggering hard fail.

若 `κ=1.0` 升到 `>=0.5`：
> FLRW at κ=1.0 triggers sustained hard fail under the current C2 rule, indicating that the flat-centered local-basin interpretation does not hold in this background window.

---

## 7. 后续升级位点（可选）

若要再上一个层级（严格“对照法”）：
- 在 `falsify_c3_metric_runner.py` 新增 `flrw_mode`（`metric_current` vs `proxy_control`）分支；
- 在同一 N/seed 下输出双轨差分（同 seeds 对照），把“生成器差异”从叙事变成量化项。
