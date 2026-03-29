# F1 一页式结论摘要（2026-03-28）

## 核心结论

**C1（Identity-layer failure）当前未被触发。**

在预注册正式阈值下，F1 全网格结果给出：

- **Hard fail = NO**
- 判定规则：任一非 Lor4D 家族需在 **至少 3 个连续 N 档位**、且每档 **至少 8/10 seeds** 下压过 Lor4D，才构成 C1-HF
- 实验网格：`N = [12, 14, 16, 20, 28, 48, 64, 96, 128, 256]`
- seeds：10

## 最关键的经验事实

### 1. Lor4D 在全部正式档位均为稳定 rank #1

所有 N 档位都满足：

- `mean rank = 1.00`
- `min rank = 1`
- `max rank = 1`
- `winner census = Lor4D:10`

换言之，在当前扩展库与当前判定设置下，Lor4D 在 10 个正式 N 档位、每档 10 个 seed 中，**全数获得 rank #1**。

### 2. 扩展家族对 Lor4D 的“压制压力”为零

`outrank_counts` 汇总显示：

- 当前 top competitors 的 `max hits at one N = 0`
- `total hits across N = 0`

这意味着：在当前扩展 family 库中，没有任何一个非 Lor4D 家族在任何一个正式 N 档位、任何一次正式统计中，形成过对 Lor4D 的稳定穿透。

## 对理论口径的含义

### 可以更有把握保留的陈述

以下口径在 **tested library** 范围内目前是安全的：

- `S_MD` 的 identity-layer 稳健识别主张，当前**未被扩展 family pressure 否证**
- Lor4D 在当前 family 库下呈现 **robust identity centre**
- C1 的最直接红队攻击路径（“存在近邻家族稳定压过 Lor4D”）当前**没有证据支持**

### 仍然不能越界升级的陈述

即便 F1 结果极强，仍不应升级为：

- 对所有 conceivable families 的普适定理
- 不依赖 tested library 的底层唯一性原理
- 与 feature basis / curved background 无关的最终结论

更稳妥表述仍应保持为：

> Within the tested family library and under the pre-registered F1 criterion, no non-Lor4D family exhibits stable repeated outranking of Lor4D; hence the current identity-layer stabilization claim survives the family-pressure falsification test.

## 建议用于主稿/摘要的短句

### 中文版

在扩展 family-pressure 检验中，我们未观察到任何非 Lor4D 家族在预注册阈值下对 Lor4D 形成稳定重复压制；因此，当前 identity-layer stabilization 主张在 tested library 范围内通过了 F1 可证伪检验。

### 英文版

Under the pre-registered F1 family-pressure test, we find no non-Lor4D family that stably and repeatedly outranks Lor4D. The current identity-layer stabilization claim therefore survives C1 falsification within the tested library.

## 相关产物

- 正式结果：`outputs_carlip/falsify_c1_family_pressure.json`
- 正式报告：`outputs_carlip/falsify_c1_family_pressure.md`
- 汇总报告：`outputs_carlip/f1_fullgrid_summary.md`
- 图：`outputs_carlip/f1_fullgrid_lor4d_rank.png`

## 当前限制

- 结论仍受限于当前 family registry
- 尚未包含 F5 合成对手家族压力测试
- 尚未与 F3（background-response）与 F4（basis challenge）联合更新最终文稿口径
