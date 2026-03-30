# 主稿/摘要/风险表口径统一清单（2026-03-30）

## A. 统一总原则（必须同时满足）

1. 不再使用“uniform mild-curvature robustness”。
2. 统一使用：`background-dependent robustness`。
3. FLRW `κ=1.0` 必须成对出现：
   - `lowN threshold hit`（或 hard-fail signal）；
   - `highN partial recovery (fail_ratio 0.3 < 0.5)`。
4. 任何版本不得写成“FLRW fully recovered”。

---

## B. 文件级检查点

| 文件 | 必须出现 | 禁止出现 |
|---|---|---|
| `进展.md` | `background-dependent` + `lowN hit + highN partial recovery` | `uniformly robust` |
| `PROJECT_OVERVIEW_2026Q1.md` | 顶部状态快照 + FLRW边界敏感提示 | “已全面恢复”式句子 |
| `PROJECT_PROGRESS.md` | 历史档标注 + 迁移到`进展.md` | 把旧口径当当前事实 |
| `outputs_carlip/DISCUSSION_THEORY_IMPLICATIONS.md` | FLRW边界段（英/中） | “all weak-to-moderate curvature pass” |
| `outputs_carlip/MANUSCRIPT_SECTIONS_1_4.md` | v0.9 口径（已完成） | 旧版“mild-curvature theorem”措辞 |

---

## C. 摘要可直接粘贴版本

### 英文 2 句版
The layered screening architecture is robust on de Sitter-like backgrounds (up to $H\le 0.3$, $N\le 1024$) but is not uniformly robust across all weak-to-moderate curved backgrounds. In matter-FLRW at $\kappa=1.0$, the low-$N$ split reaches the hard-fail threshold, while completed high-$N$ branches show partial recovery (failure ratio $0.3<0.5$), indicating a boundary-sensitive regime rather than a universal mild-curvature theorem.

### 中文内部版
两层筛选架构在 de Sitter 类背景上表现稳健，但在弱到中等曲率背景上并非统一稳健。对于 FLRW 的 $\kappa=1.0$，lowN 触发阈值信号，而 highN 显示部分回升（failure ratio $0.3<0.5$），因此应表述为边界敏感区而非统一曲率定理。

---

## D. 风险表统一标签（建议）

| 风险项 | 标签 | 当前状态 |
|---|---|---|
| C1 家族压力 | LOW | 通过（Hard fail=NO） |
| C2 背景响应（FLRW κ=1.0） | MED-HIGH (controlled) | lowN hit + highN partial recovery |
| C3 特征基充分性 | MED | 待 F4 挑战赛封口 |

---

## E. 发布前最终自检（5 项）

- [ ] 摘要/结论是否仍含“uniform”字样
- [ ] 是否同时写出 lowN 与 highN 两端事实
- [ ] 是否把 `0.3<0.5` 与 hard-fail 规则对应写清
- [ ] 是否在 limitations 标注“background-dependent”边界
- [ ] 是否避免把 FLRW 边界写成主理论失败
