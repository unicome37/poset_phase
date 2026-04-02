# PhaseB 阶段报告（草案）

生成时间：2026-04-02  
项目：`falsify_c3_metric_runner.py`（FLRW metric-faithful，PhaseB）

## 一、执行概况

本阶段聚焦 `flrw` 家族在两组参数（κ=0.3, κ=1.0）下的 rank 稳定性与 hard-fail 判定，采用多组配置并行验证：

- `quick_n1024`：reps=6, seed_runs=4
- `ultra_n1024`：reps=4, seed_runs=3
- `ultra_n1536`：reps=4, seed_runs=3
- `quick_n1536`：reps=6, seed_runs=4

主线 `phaseB`（reps=12, N=[1024,1536,2048], seed_runs=10）仍在长期运行。

## 二、已完成实验结果

### 1) quick_n1024

- 结果：`HARD_FAIL = NO`
- 记录：
	- κ=0.3: [2, 2, 2, 2]
	- κ=1.0: [3, 2, 2, 2]

### 2) ultra_n1024

- 结果：`HARD_FAIL = YES`
- 记录：
	- κ=0.3: [2, 1, 2]
	- κ=1.0: [5, 4, 2]

### 3) ultra_n1536

- 结果：`HARD_FAIL = YES`
- 记录：
	- κ=0.3: [3, 2, 2]
	- κ=1.0: [4, 2, 3]

### 4) quick_n1536

- 结果：`HARD_FAIL = YES`
- 记录：
	- κ=0.3: [2, 2, 2, 2]
	- κ=1.0: [3, 2, 3, 2]

## 三、阶段性统计结论

### κ=0.3（全部已完成样本）

- rank 范围：1~3
- 主峰：2
- 解释：弱参数窗口表现稳定，区分性强。

### κ=1.0（全部已完成样本）

- rank 范围：2~5
- 在多个配置中出现 `rank>2` 的情形
- 解释：中等参数窗口更敏感，方差明显更高。

### hard-fail 结果聚合

- YES：3 组
- NO：1 组

## 四、当前风险与解释

1. **判定对配置敏感**：不同 reps/seed_runs 组合下，hard-fail 结论可发生变化。  
2. **脚本汇总口径影响解读**：报告中的 `fail_ratio=1.00 (1/1)` 属于当前 bin 级口径，应与原始 rank 记录联合解读。  
3. **2048 已进入首轮 checkpoint**：30分钟轮询第1轮后，`quick_n2048` 与 `ultra_n2048` 均已写出首个 `partial`。

## 五、在跑任务状态（截至本草案）

- 主线 `phaseB`：运行中（进度慢）
- `quick_n2048`：`partial` 已出现，当前 `completed_seed_runs=1/4`
- `ultra_n2048`：`partial` 已出现，当前 `completed_seed_runs=1/3`

## 六、下一步建议

1. 对 2048 两路采用**30分钟低频轮询**，确认首个 `partial` 出现时间。  
2. 主线 `phaseB` 维持运行，作为最终高可信结论来源。  
3. 在 2048 产出首个 seed 后，补充一次“跨 N 对比”中期报告。

---

## 方法学附录：为何出现 `fail_ratio=1.00 (1/1)`

### A. 现象

在 `ultra_n1024 / ultra_n1536 / quick_n1536` 的最终 `json` 中，`hard_fail_items` 出现：

- `family=flrw`
- `param=1.0`
- `fail_ratio=1.00 (1/1)`

这与读者直觉中的“按 seeds 统计失败比例”不同（例如 4 个 seed 不应显示为 1/1）。

### B. 原因（当前脚本口径）

当前 runner 的 hard-fail 汇总并非直接按 `seed_runs` 计数，而是将失败事件聚合到更高层的 **bin（参数-网格统计桶）** 再计算比率。对于当前配置，这个口径常对应单一桶：

- `total_bins = 1`
- 若该桶触发失败条件，则 `fail_bins = 1`
- 从而显示为 `fail_ratio = 1.00 (1/1)`

### C. 正确解读方式

应采用“双通道解读”而不是只看 `hard_fail_items`：

1. **判定通道**：看 `hard_fail` 与 `hard_fail_items`（脚本当前官方口径）
2. **证据通道**：回到 `records`，直接检查各 seed 的 rank 序列

本报告所有核心结论已按第2通道做了交叉核对（逐 seed rank 已列出）。

### D. 对后续报告的建议措辞

在正文中建议统一写法：

> “当前 hard-fail 采用 bin-level 汇总口径，`fail_ratio=1.00(1/1)` 表示该参数统计桶触发失败；详细 seed 级证据见 records 表。”

这样可同时保持与脚本输出一致，并避免误读为“只有 1 个 seed”。

---

## 轮询附记（30分钟第1轮）

轮询时间：2026-04-02 20:30 左右

- `quick_n2048.partial.json`：已出现，进度 `1/4`
  - seed0 ranks：κ=0.3 → 2，κ=1.0 → 2
- `ultra_n2048.partial.json`：已出现，进度 `1/3`
  - seed0 ranks：κ=0.3 → 2，κ=1.0 → 7

这表明 2048 在高参数（κ=1.0）下的 rank 波动可能进一步扩大，需等待更多 seed 收敛后再下定论。

---

## 附录B：2048 风险预警（并入版）

### 当前事实（首轮checkpoint）

- `quick_n2048`：`completed_seed_runs=1/4`
  - seed0：κ=0.3 → rank=2，κ=1.0 → rank=2
- `ultra_n2048`：`completed_seed_runs=1/3`
  - seed0：κ=0.3 → rank=2，κ=1.0 → rank=7

### 风险分级

1. **结论不稳定风险：高**  
   同一 N=2048 下，κ=1.0 在 quick/ultra 首seed已出现 2 vs 7 的分叉。

2. **误判风险（假阴性/假阳性）：中高**  
   仅凭首seed容易得出方向性错误结论。

3. **进度延迟风险：中**  
   2048 单seed成本高，收敛慢，需控制并行度并采用低频轮询。

### 执行建议（已生效）

- 已切换 30 分钟低频轮询策略（按指令可暂停/恢复）。
- 下一关键观察点：两路从 `1/x` 推进到 `2/x` 后再做趋势定性。

