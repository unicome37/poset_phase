# P1 执行细化：F2 turn-on 再估计（2026-03-30）

## 1. 目标

在小 N 边界区间重估 `N_id`（identity turn-on），把“边界叙事”升级为“可重复统计判定”。

---

## 2. 执行配置

- 主配置：`config_falsify_c1_turnon_refit.yaml`
- smoke：`config_falsify_c1_turnon_refit_smoke.yaml`
- runner：`falsify_c1_runner.py`

### 关键参数
- N 网格：`[10,12,14,16,18,20,22,24]`
- seeds：`20`
- reps：`120`
- 规则映射：`min_seed_success=16`（等价于原 10-seed 规则中的 8/10 强度）

---

## 3. 判定表（建议写入最终报告）

| 指标 | 判定条件 | 解释 |
|---|---|---|
| 稳定开启（Strict） | 连续 3 个 N 档位均 `20/20` Lor4D #1 且 min-margin > 0 | 最强封口 |
| 稳定开启（Practical） | 连续 3 个 N 档位均 `>=18/20` Lor4D #1 且下置信界 > 0 | 现实可接受封口 |
| 边界波动区 | 出现 `16~17/20` 的零散档位 | 需继续补样本/补参考 reps |
| 未开启 | `<=15/20` 且无连续改善 | 不应宣称 turn-on 已稳定 |

---

## 4. 推荐输出字段

从 `per_seed_records` 聚合生成：
- `N`
- `lor4d_rank1_hits`
- `rank1_rate`
- `runner_up_family_census`
- `min_margin / mean_margin / CI95`
- `is_in_stable_block`

建议最终追加：
- `N_id_lower, N_id_point, N_id_upper`

---

## 5. 与主线口径的关系

P1 的价值不是“再做一次小实验”，而是修复两层主线最后一块容易被问的区域：

- 如果 `N_id` 在 `14~16` 稳定：可维持当前主线表达；
- 如果 `N_id` 推迟到 `18+`：需同步调整所有文稿的下限口径；
- 如果出现多次反复：需在方法学中强调 reference 估计精度对小 N 的影响。

---

## 6. 一句话写法（跑完可直接用）

> The turn-on boundary is re-estimated on a dense low-N grid (N=10–24, 20 seeds, 120 reps), yielding a reproducible onset block rather than a single-point claim.
