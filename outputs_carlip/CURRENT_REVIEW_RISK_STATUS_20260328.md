# 当前审稿风险状态摘要（2026-03-28）

## 一句话结论

当前最主要的审稿风险已从“F1 家族压力是否击穿 identity layer”转移为“F3 中 FLRW 背景是否实质性击穿 local-basin 解释”。

## 当前已锁定结论

### 1. F1 风险：当前基本解除

- 正式 family-pressure run：`N=12–256`, `10 seeds`
- 结果：**Hard fail = NO**
- Lor4D 在全部 tested N / 全部 seeds 上均保持 rank #1
- 当前扩展库下 competitor pressure = 0

因此，C1 当前不是主要审稿攻击面。

### 2. F3 风险：已出现实质性背景依赖

当前 lowN split（`N=256,512`, `seed_runs=5`, `reps=10`）结果：

- `de Sitter`：Pass
- `Schwarzschild`：Pass
- `FLRW (kappa=1.0)`：**Hard Fail**

这意味着：

- 不能再把“weak/mild curvature robustness”写成统一成立
- 当前最安全表述应为：**background-dependent robustness**
- 当前 local-basin 解释并未被全面推翻，但已被明确限定：它在不同背景上的成立性不均匀

### 3. 小N turn-on 风险：当前下降，但 manuscript 口径暂不前推

- 新的 winner-only refit（`N=10..24`, `20 seeds`, `120 reps`）给出：Lor4D 在全部 160 个 `(seed, N)` 单元中均为 #1。
- 这说明在当前 family-pressure 风格判定下，小N边界风险显著下降。
- 但该 refit 没有复刻旧版 turn-on 表中的 `min margin / cond(Σ)` 维度，因此**不宜直接把主稿边界从 `N>=14` 改写成 `N>=10`**。
- 当前最稳策略：
   - **内部状态**：记录 lowN winner-only 全胜已前推到 `N=10`；
   - **外部 manuscript**：继续保守维持 `N>=14`，直到完成 margin-aware refit。

## 风险排序（从高到低）

1. **最高风险：FLRW highN 是否继续 fail**
   - 如果 `kappa=1.0` 在 `N=768/1024` 继续 fail，主稿需要进一步降级；
   - 如果 fail 扩大，则 C2 风险从“局部边界现象”升级为“更系统的背景响应失败”。

2. **中等风险：文稿口径是否及时跟上结果**
   - 这一项已经开始处理；主稿/提纲已完成首轮降级。

3. **较低风险：de Sitter 主结论是否动摇**
   - 当前没有这个迹象；de Sitter 仍是最稳的一条线。

4. **较低风险：Schwarzschild 是否成为第二个失败源**
   - 当前 lowN 没有迹象，但 highN 尚未验证。

5. **较低风险：是否过早前推 turn-on 边界**
   - 当前已有更强内部结果，但若不区分判定协议，容易造成文稿自相矛盾；
   - 该风险可通过“双层口径”（内部前推 / 对外保守）控制。

## 当前已完成的风控动作

- 主稿已从“统一弱曲率稳健”改为“分背景成立”
- 提纲已同步改写
- 可证伪计划与执行清单已写入 lowN split 结果
- 已单独启动 `FLRW highN` 继续追踪风险是否持续

## 下一步最关键问题

> `FLRW (kappa=1.0)` 在 highN（768/1024）是继续 fail，还是恢复到 top-2？

这将决定：

- 我们只是需要把 FLRW 写成“边界敏感”；
- 还是必须进一步改写为“当前 local-basin 解释对 FLRW 不成立”。
