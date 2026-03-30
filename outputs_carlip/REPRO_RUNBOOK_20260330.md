# 可复跑说明书：F1 / F2 / F3 / P0（2026-03-30）

## 1. 目的

这份 runbook 只服务当前 post-Carlip 主线：

- F1 家族压力可证伪
- F2 margin-aware onset refit
- P1 winner-only onset refit
- F3 / P0 背景响应（重点：FLRW `κ=1.0`）

> 历史探索脚本很多，但当前投稿/防守应优先以本页列出的 runner + config + output 映射为准。

---

## 2. 环境前提

- 工作目录：仓库根目录 `d:\Kiro\理论体系\poset_phase`
- Python 环境：安装 `requirements.txt`
- 默认输出目录：`outputs_carlip/`

---

## 3. 推荐执行顺序

`F1 -> F2 -> P1 -> P0-PhaseA -> P0-PhaseA-Control -> (可选) P0-PhaseB`

理由：

1. 先锁定家族压力与小N边界；
2. 再处理 FLRW 边界主分支；
3. 最后决定是否需要更高 N 的封口。

---

## 4. Runner / Config / Output 映射

| 任务 | Runner | Config | 主要输出 | 当前用途 |
|---|---|---|---|---|
| F1 family pressure | `falsify_c1_runner.py` | `config_falsify_c1.yaml` | `falsify_c1_family_pressure.{json,md}` | 家族压力可证伪主结果 |
| F2 margin-aware onset | `f2_turnon_margin_runner.py` | `config_f2_turnon_margin_refit.yaml` | `f2_turnon_margin_refit.{json,md}` | 当前 onset = `N≥10` 的主证据 |
| F2 smoke | `f2_turnon_margin_runner.py` | `config_f2_turnon_margin_refit_smoke.yaml` | `f2_turnon_margin_refit_smoke.{json,md}` | 环境/路径快速确认 |
| P1 winner-only refit | `falsify_c1_runner.py` | `config_falsify_c1_turnon_refit.yaml` | `falsify_c1_turnon_refit.{json,md}` | 与 F2 交叉印证 |
| P1 smoke | `falsify_c1_runner.py` | `config_falsify_c1_turnon_refit_smoke.yaml` | 对应 `*_smoke` 输出 | 快速检查 runner |
| P0 PhaseA (metric) | `falsify_c3_metric_runner.py` | `config_falsify_c3_flrw_metricfaithful_phaseA.yaml` | `falsify_c2_background_response_flrw_metricfaithful_phaseA.{json,md}` | FLRW metric-faithful 主证据 |
| P0 PhaseA control | `falsify_c3_metric_runner.py` | `config_falsify_c3_flrw_proxy_control_phaseA.yaml` | `falsify_c2_background_response_flrw_proxy_control_phaseA.{json,md}` | 区分 sampling vs causal-cone effect |
| P0 PhaseB (optional) | `falsify_c3_metric_runner.py` | `config_falsify_c3_flrw_metricfaithful_phaseB.yaml` | `falsify_c2_background_response_flrw_metricfaithful_phaseB.{json,md}` | 更高 N 封口 |

---

## 5. 标准执行命令

### F1

`python falsify_c1_runner.py --config config_falsify_c1.yaml`

### F2

`python f2_turnon_margin_runner.py --config config_f2_turnon_margin_refit.yaml`

### P1

`python falsify_c1_runner.py --config config_falsify_c1_turnon_refit.yaml`

### P0 PhaseA (metric)

`python falsify_c3_metric_runner.py --config config_falsify_c3_flrw_metricfaithful_phaseA.yaml`

### P0 PhaseA control

`python falsify_c3_metric_runner.py --config config_falsify_c3_flrw_proxy_control_phaseA.yaml`

### P0 PhaseB (optional)

`python falsify_c3_metric_runner.py --config config_falsify_c3_flrw_metricfaithful_phaseB.yaml`

---

## 6. 如何读结果

### F1

- 看 `hard_fail`
- 当前安全结论：`Hard fail = NO`

### F2

- 看三档 onset：`winner_only / operational / manuscript_safe`
- 当前安全结论：三档都在 **`N=10`** 开启

### P1

- 看 Lor4D 是否在每个 `(seed, N)` 都 rank #1
- 当前已知结果：`160/160` 全 rank #1

### P0 / F3

- 看 `fail_ratio`
- 当前安全解释：
  - lowN: threshold hit
  - highN: partial recovery (`0.3 < 0.5`)
  - P0 PhaseA: metric `11/60 = 0.183`, proxy `0/60`

---

## 7. 口径映射

| 结果 | 应转成的文字 |
|---|---|
| F2 onset = 10 | `Lor4D identity turn-on is supported from N≥10 under the fixed-reference protocol.` |
| F1 hard_fail = NO | `The current family-pressure falsification does not trigger hard fail.` |
| FLRW lowN hit + highN partial recovery | `background-dependent robustness`, not `uniform mild-curvature robustness` |
| P0 metric 11/60, proxy 0/60 | `boundary-sensitive degradation in the metric branch, with stable proxy control` |

---

## 8. 相关文档

- 真值源：`进展.md`
- 口径清单：`WRITING_UNIFICATION_CHECKLIST_20260330.md`
- 短执行入口：`NEXT_RUN_ENTRYPOINT_20260330.md`
- FLRW 防守包：`FLRW_KAPPA1_DEFENSE_BRIEF_20260330.md`

---

## 9. 当前建议

若目标是投稿/防守，而不是继续探索：

1. 直接引用本 runbook 对应的 outputs；
2. 所有文稿统一为 `N≥10` + `background-dependent robustness`；
3. 仅在需要更强防守时再追加 P0 PhaseB。