# FLRW κ=1.0 一页式防守包（2026-03-30）

## 1. 最安全的一句话

> FLRW at `κ=1.0` is currently a **boundary-sensitive** case of the flat-centered local-basin picture: lowN split reaches the hard-fail threshold, highN shows partial recovery (`fail_ratio = 0.3 < 0.5`), and the completed metric-faithful PhaseA run yields metric-branch degradation without triggering the preregistered hard-fail rule (`11/60 = 0.183`).

中文：

> FLRW 的 `κ=1.0` 当前应被定位为 flat-centered local-basin 解释中的**边界敏感区**：lowN split 触及 hard-fail 阈值，highN 出现部分回升（`fail_ratio = 0.3 < 0.5`），而 metric-faithful PhaseA 完整运行显示 metric 分支退化但未触发预注册 hard-fail 规则（`11/60 = 0.183`）。

---

## 2. 三段证据链

| 证据层 | 设计 | 结果 | 安全解释 |
|---|---|---|---|
| lowN split | `N=256,512` | `κ=1.0` 触及 hard-fail threshold | 说明 FLRW 不是“统一稳健”背景 |
| highN completion | `N=768,1024`, 5/5 seeds | `fail_ratio = 0.3 < 0.5` | 出现部分回升，但不能写成 fully recovered |
| P0 metric-faithful PhaseA | `N=512,768,1024`, 10 seeds, metric vs proxy | metric `11/60 = 0.183`, proxy `0/60` | 问题集中在 metric-faithful 分支的边界敏感退化，而非全背景崩塌 |

---

## 3. 必须成对出现的事实

写 FLRW `κ=1.0` 时，必须同时写：

1. `lowN threshold hit`（或 lowN hard-fail signal）
2. `highN partial recovery (fail_ratio 0.3 < 0.5)`
3. 若提到 P0，则补上：metric `11/60 = 0.183`, proxy `0/60`

缺任一项，叙事都会失真。

---

## 4. 推荐口径

### 英文短版（摘要 / cover letter 可用）

The current curved-background evidence supports a **background-dependent robustness** claim rather than a uniform mild-curvature theorem. In matter-FLRW at `κ=1.0`, the low-N split reaches the hard-fail threshold, completed high-N runs show partial recovery (`fail_ratio = 0.3 < 0.5`), and the metric-faithful PhaseA run records boundary-sensitive degradation in the metric branch (`11/60`) while the proxy branch remains top-2 throughout.

### 中文短版（内部统一）

当前曲率结果支持“**分背景稳健**”而非“统一弱曲率稳健”。在 FLRW 的 `κ=1.0` 上，lowN split 触及 hard-fail 阈值，highN 出现部分回升（`0.3 < 0.5`），而 metric-faithful PhaseA 中 metric 分支出现边界敏感退化（`11/60`），proxy 分支则保持全程 top-2。

---

## 5. 禁止写法

- `FLRW fully recovered`
- `all weak-to-moderate curved backgrounds pass`
- `uniform mild-curvature robustness`
- `FLRW failure falsifies the whole two-layer framework`

---

## 6. 审稿人若追问，应怎么答

### Q1. 这是否意味着两层筛选理论失败？

**答：不是。**

当前结果表明失败是**背景特异**的，不是全背景普遍失效。de Sitter-like 与 weak-field Schwarzschild 在当前窗口仍兼容 local-basin 图景；FLRW `κ=1.0` 是边界敏感点。

### Q2. 为什么不直接说 FLRW 不稳健？

**答：因为证据不是单向崩塌。**

lowN 与 highN 给出的不是同一方向的完全失败，而是“threshold hit → partial recovery”的混合格局；P0 进一步表明问题主要出现在 metric-faithful 分支，而 proxy 分支保持稳定。

### Q3. 下一步最有价值的数据是什么？

**答：PhaseB。**

如果需要更强防守，应继续跑 `N=1024/1536/2048` 的 metric-faithful PhaseB，看 `κ=1.0` 的 `fail_ratio` 是继续下降、平台化，还是重新上升。

---

## 7. 文件索引

- 真值源：`进展.md`
- 风险总览：`CURRENT_REVIEW_RISK_STATUS_20260328.md`
- 口径清单：`WRITING_UNIFICATION_CHECKLIST_20260330.md`
- 执行细化：`P0_FLRW_METRIC_FAITHFUL_PLAN_20260330.md`
- 当前主稿：`MANUSCRIPT_SECTIONS_1_4.md`