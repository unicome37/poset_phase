# F2 Margin-Aware Refit 方案（2026-03-30）

## 目标

把当前 `winner-only` 的 low-N refit，升级为与旧版 turn-on 表可直接比较的 `margin-aware` 版本，从而回答：

> manuscript 的 turn-on 边界是否可以从 `N>=14` 正式前推？

---

## 为什么需要新 refit

当前已完成结果：
- `falsify_c1_turnon_refit.json`
- `N=10..24`, `20 seeds`, `120 reps`
- Lor4D = `160/160` 全 #1

但当前 runner 不输出：
- `mean margin`
- `min margin`
- `runner-up census`
- `cond(Σ)`

因此它更像是：
- **family-pressure 风格 winner test**

而不是：
- **turn-on stability test**

---

## 需要补的统计量

对每个 `(seed, N)`，除 rank 外新增：

1. `lor4d_score`
2. `runner_up_family`
3. `runner_up_score`
4. `margin = runner_up_score - lor4d_score`
5. `cond_sigma`
6. `mu_d, mu_c, mu_w`

对每个 `N` 聚合：

1. `rank1_rate`
2. `mean_margin`
3. `min_margin`
4. `margin_ci95`
5. `top_intruder_census`
6. `max_cond_sigma`

---

## 推荐判定分级

### A. Manuscript-safe onset
满足以下全部条件：
- 连续 3 个 N 档位 `20/20` 为 #1
- `min_margin > 0`
- `margin_ci95` 下界 > 0
- `max_cond_sigma < 60`

### B. Operational onset
满足：
- 连续 3 个 N 档位 `>=18/20` 为 #1
- `mean_margin > 0`

适合内部状态汇报，不直接替换论文边界。

### C. Winner-only onset
仅满足：
- 连续 3 个 N 档位 `20/20` 为 #1

这是当前 `falsify_c1_turnon_refit` 的等级。

---

## 推荐输出文件

- `outputs_carlip/f2_turnon_margin_refit.json`
- `outputs_carlip/f2_turnon_margin_refit.md`
- `outputs_carlip/f2_turnon_margin_refit.seedlog.jsonl`
- `outputs_carlip/f2_turnon_margin_refit.partial.json`

---

## 最小实现路径

### 方案 1（推荐）
复制 `falsify_c1_runner.py` 为新 runner：
- `f2_turnon_margin_runner.py`

理由：
- 不污染当前 C1 family-pressure runner
- turn-on 统计逻辑与 C1 判定逻辑可分开维护
- 更容易在文稿中解释“这是专用 turn-on refit”

### 方案 2（不推荐）
继续在 `falsify_c1_runner.py` 上堆分支。

问题：
- C1 与 F2 目标开始混杂
- 后续更难解释哪个输出服务于哪个理论命题

---

## manuscript 建议用语（跑完后）

若 margin-aware 仍支持 `N=10/12`：
> A denser margin-aware low-N refit confirms that the onset block can be pushed below the previously quoted `N>=14` boundary.

若 only winner-aware 支持而 margin-aware 不支持：
> Winner-only refits can be pushed to lower N, but the conservative onset boundary in the manuscript is still kept at `N>=14` because margin-stability criteria have not yet been met below that scale.

---

## 当前建议

- **立即可做**：继续保持 P0 运行，不打断；
- **下一轮可做**：新建 `f2_turnon_margin_runner.py`，按上面字段补全输出；
- **对外口径**：在 margin-aware refit 完成前，不改 manuscript 边界。