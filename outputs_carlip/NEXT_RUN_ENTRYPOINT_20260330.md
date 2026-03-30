# 下一步执行入口（2026-03-30）

## 1) 先跑 smoke（确认环境与输出路径）

- 配置：`config_falsify_c1_turnon_refit_smoke.yaml`
- 目标：快速验证 `falsify_c1_runner.py` 在当前环境可正常出 JSON/MD。

## 2) 正式跑 P1（F2 turn-on 再估计）

- 配置：`config_falsify_c1_turnon_refit.yaml`
- 目标：得到 `N=10..24` 的稳定开启区间与 `N_id` 后验描述。

## 3) 正式跑 P0-PhaseA（FLRW 边界复核）

- 配置：`config_falsify_c3_flrw_metricfaithful_phaseA.yaml`
- 目标：复核 `κ=1.0` 在 `N=512/768/1024` 的 fail_ratio。

## 3.5) 跑 P0-PhaseA-Control（FLRW 对照分支）

- 配置：`config_falsify_c3_flrw_proxy_control_phaseA.yaml`
- 目标：与 PhaseA 同网格对照，分离抽样权重与因果锥效应。

## 4) 正式跑 P0-PhaseB（跨尺度封口）

- 配置：`config_falsify_c3_flrw_metricfaithful_phaseB.yaml`
- 目标：确认 `N=1024/1536/2048` 的趋势归类（回落/平台/恶化）。

---

## 结果文稿联动

- 口径清单：`outputs_carlip/WRITING_UNIFICATION_CHECKLIST_20260330.md`
- P0 文档：`outputs_carlip/P0_FLRW_METRIC_FAITHFUL_PLAN_20260330.md`
- P1 文档：`outputs_carlip/P1_F2_TURNON_REESTIMATION_PLAN_20260330.md`

---

## 建议先后顺序

`P1 -> P0A -> P0A-Control -> P0B -> 文稿口径一次性同步`

原因：先锁定小N边界，再做 FLRW 主分支+对照分支封口，最后避免反复改稿。