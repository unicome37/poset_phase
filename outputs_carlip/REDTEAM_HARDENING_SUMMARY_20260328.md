# Red-Team Hardening Summary (2026-03-28)

## 本轮目标

针对“表述过强、对象混淆、审稿人降格为 one-class novelty detector”的风险，对核心文稿做降强度与可证伪化修订。

## 已修改文件

1. `outputs_carlip/MANUSCRIPT_SECTIONS_1_4.md`
2. `outputs_carlip/LAYERED_SCREENING_PRINCIPLE_OUTLINE.md`

## 关键修订点

### 1) 强词降半步（已执行）

- `selects / positively selects` → `identifies / identity centre`
- `phase-transition-like turn-on` → `sharp finite-size turn-on / onset threshold`
- `minimum sufficient` → `minimal non-redundant effective basis within the tested feature library`
- `This is not machine learning` → `reference-based one-class geometric discriminator (not trained black-box)`

### 2) 防“one-class降格攻击”（已执行）

在主稿中新增明确锚点：

- self-minimum 是通用 one-class 属性；
- Lor4D 非平凡性不锚定 self-minimum，锚定 **组合证据**：
  - cross-family separation
  - finite-size turn-on
  - basin deepening
  - mild-curvature robustness（在 tested library 内）

### 3) 增加可证伪条件（已执行）

新增专节：`§6.8 Falsifiability conditions`（主稿）与 `8.5 可证伪条件`（提纲），包含 3 条失败判据：

1. Identity-layer failure
2. Background-response failure
3. Basis-sufficiency failure

### 4) 统计措辞降风险（已执行）

将“近弱相关（r→0）”改为“弱耦合/无稳定单调依赖”口径，避免“正交化”过度解释。

### 5) 曲率稳健性口径再降级（已执行）

随着 F3 lowN split 结果落地，曲率相关表述已进一步从“weak/mild curvature 基本稳健”降级为“**background-dependent robustness**”：

- `de Sitter`：通过
- `weak-field Schwarzschild`：通过
- `FLRW (\kappa=1.0)`：已触发当前 C2 hard-fail 门槛

因此，主稿与提纲中所有“统一弱曲率稳健”的说法都已改写为“分背景成立、且 FLRW 边界已出现失守信号”。

## 仍建议后续处理（下一轮）

1. 对 `SMD_OPERATOR_LETTER.md` 做同口径降强度（仍有 “Phase I/II + phase-like” 语气）。
2. 对 `PROJECT_OVERVIEW_2026Q1.md`、`WORKPLAN_2_3_WEEKS.md` 等对外摘要做术语统一，避免旧措辞回流。
3. 若准备投稿版，建议在摘要与结论首段都重复“mid-level effective theory + tested library scope”。

## 本轮结论

红队风险从“容易被降格/反写”降到“可辩护且边界清晰”。
当前版本更像：

- **有边界的 reference-based effective geometric discriminator theory**，
而非
- 过度宣称的“底层普适选择原理”。

补充：在曲率问题上，风险结构已进一步清晰化——真正的当前审稿攻击点不是 de Sitter，而是 `FLRW (\kappa=1.0)` 对 local-basin 解释的边界挑战。
