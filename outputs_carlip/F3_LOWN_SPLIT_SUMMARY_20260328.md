# F3 lowN 分项结果摘要（2026-03-28）

## 范围

本摘要汇总 F3 在 lowN 窗口 `N = [256, 512]`、`seed_runs = 5`、`reps = 10` 下的分项结果。

对应目标：定位 C2（background-response failure）在 lowN 条件下究竟由哪类背景触发。

## 总结论

### 1. de Sitter（lowN）

- `N=256`：**Hard fail = NO**
- `N=512`：**Hard fail = NO**
- 观测：
  - `H=0.1` 与 `H=0.3` 基本保持 top-2
  - 个别 seed 在 `N=256, H=0.1` 出现 rank #1，但不影响“近邻层级稳定”结论

结论：在当前 lowN 分项测试下，**de Sitter 未触发 C2 hard fail**。

### 2. FLRW（lowN）

- `N=[256,512]`：**Hard fail = YES**
- 触发项：`flrw, kappa=1.0`
- 判定细节：
  - `fail_ratio = 0.50 (1/2)`
  - 满足当前 lowN 配置下的 hard-fail 门槛
- 具体表现：
  - `kappa=0.3` 始终保持 rank #2
  - `kappa=1.0` 在 `N=256` 多次跌出 top-2（rank #3），并在 `N=512` 个别 seed 跌至 rank #4

结论：在当前 lowN 分项测试下，**FLRW 是当前 C2 hard fail 的来源**。

### 3. Schwarzschild（lowN）

- `N=[256,512]`：**Hard fail = NO**
- 观测：
  - `phi0=0.01` 多数为 rank #2，个别 seed 在 `N=256` 达到 rank #1
  - `phi0=0.05` 稳定为 rank #2

结论：在当前 lowN 分项测试下，**Schwarzschild 未触发 C2 hard fail**。

## 当前解释边界

在 lowN 窗口内，C2 的风险并不是“所有弱/中曲率背景都失败”，而是更具体地表现为：

- `de Sitter`：通过
- `Schwarzschild`：通过
- `FLRW(kappa=1.0)`：触发 hard fail

因此，当前最合理的文稿口径不是笼统宣称“弱/中曲率全部稳健”，而应改成：

> lowN 分项结果表明，de Sitter 与 weak-field Schwarzschild 在当前设置下维持近邻层级，而 matter-FLRW 在 `kappa=1.0` 处已出现可重复的 top-2 失守信号，提示 background-response 稳健性具有背景依赖性，而非统一成立。

## 相关文件

- `outputs_carlip/falsify_c2_background_response_desitter_n256.json`
- `outputs_carlip/falsify_c2_background_response_desitter_n512.json`
- `outputs_carlip/falsify_c2_background_response_flrw_lown.json`
- `outputs_carlip/falsify_c2_background_response_schwarz_lown.json`

## 建议下一步

1. 跑 `FLRW highN`，验证 `kappa=1.0` 的 fail 是否持续扩大；
2. 跑 `de Sitter highN`，检查其 lowN 通过是否延续到 `N=768/1024`；
3. 将主稿中“weak-to-moderate curvature robust”改写为分背景口径，至少在 FLRW 上加入条件限制。
