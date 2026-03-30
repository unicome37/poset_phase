# P1 结果为何强于旧版 turn-on 表：诊断说明（2026-03-30）

## 一句话结论

`falsify_c1_turnon_refit` 给出 `N=10..24` 全 #1，并不等于可以直接用它替换旧版 `N>=14` turn-on 口径。更准确的说法是：

> 在**当前 family-pressure 风格 runner** 与 **当前 reps=120 的集成均值判定**下，Lor4D 在 `N=10..24` 上表现为全档稳定 #1；但这比旧版 turn-on 表更宽松/更实用，不是同一严格度下的逐项可比复现。

---

## 1. 旧版 turn-on 表在测什么

旧文档：`outputs_carlip/mahalanobis_n_boundary_turnon.md`

其核心量包括：
- `#1 rate`
- `mean rank`
- `mean margin`
- `min margin`
- `cond(Σ)`
- `top intruder`

并明确把 `N=12` 判为边界外，原因不是单纯“有没有 #1”，而是：
- 10 seeds 中只有 `6/10` 为 #1；
- 还出现负 margin；
- intruder（Lor5D / KR_2layer）真实压过 Lor4D；
- 文档将其解释为**物理分辨率极限**。

这套口径本质上是：

> 既看 winner，也看 margin 稳定性与 reference 条件数。

---

## 2. 新版 P1 runner 在测什么

当前运行器：`falsify_c1_runner.py`

它的正式输出是：
- 每个 `(seed_run, N)` 下 Lor4D 是否 rank #1；
- 任意非 Lor4D 家族是否构成 C1 型连续压制；
- 不输出 `margin`；
- 不输出 `min margin`；
- 不把“近乎并列但仍 #1”与“强分离 #1”区分开。

因此当前 P1 结果回答的是：

> 在当前配置下，有没有任何 family 真的压过 Lor4D？

而不是：

> Lor4D 是否已经进入“有明显安全边际”的 turn-on 稳定区？

---

## 3. 为什么这次结果会更强

至少有三点实质原因：

### (A) 判定标准更偏“winner-only”
旧表要求的是“#1 + 边际稳定”。
当前 P1 主要看“是不是 #1”。

这会天然把边界往更小 N 推，因为：
- 只要 Lor4D 以很小优势赢一次，仍算成功；
- 不要求 `min margin > 0.3` 之类的安全边际。

### (B) reps 从 80 提升到 120
当前 P1：
- `20 seeds`
- `120 reps`

旧表：
- `10 seeds`
- `80 reps`

更高 reps 会降低 Lor4D reference 估计噪声，使小 N 下的 family 均值更稳定，更不容易被 Lor5D/KR_2layer 偶发穿透。

### (C) 当前 P1 是“family-pressure式再估计”，不是旧 turn-on 脚本原样复跑
因此它更适合回答：
- `C1 是否在小N就已经完全不触发？`

但未必能直接回答：
- `论文里的 onset boundary 是否应从 14 改写成 10？`

---

## 4. 最安全的新口径

### 不建议直接写
- “turn-on boundary 已更新为 N>=10”
- “旧的 N>=14 结论被推翻”

### 建议写法

> Under the current family-pressure style re-estimation (20 seeds, 120 reps), Lor4D remains rank #1 throughout the dense low-N grid down to N=10. However, this result is not strictly identical to the earlier turn-on table, which also tracked margin stability and conditioning; the conservative manuscript boundary may therefore remain phrased as `N>=14` until a margin-aware refit is completed.

中文内部版：

> 在当前 family-pressure 风格的再估计下（20 seeds, 120 reps），Lor4D 在低 N 密集网格上可一直保持 #1 到 N=10；但该结果与旧版同时考察 margin 稳定性的 turn-on 表并非同一判定协议，因此主稿仍宜暂时保守维持 `N>=14`，直到完成 margin-aware refit 为止。

---

## 5. 下一步若要“正式改边界”，还差什么

需要一个真正与旧表同口径的 `margin-aware F2 refit`：
- 输出 `mean margin / min margin / CI95 / cond(Σ)`；
- 在 `N=10..24` 上跑 `20 seeds, 120 reps`；
- 再决定是否把 `N>=14` 改写为 `N>=10` 或 `N>=12`。

---

## 6. 当前建议

- **对外主稿**：继续保守写 `N>=14`。
- **对内状态更新**：记录“winner-only refit 已前推到 N=10”。
- **方法学备注**：说明当前 P1 比旧表更偏 operational / falsification 风格，而非原封不动的 turn-on stability 表。
