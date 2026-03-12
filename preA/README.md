# Prediction A — Why 3+1?

## Paper

**Title**: Dimensional Selection Without Dimensional Priors: 4D Lorentzian Dominance in Finite Causal Poset Ensembles Under Consistency-Based Actions

**Author**: Gang Zhang (Independent Researcher)

**Target Journal**: MDPI *Entropy*

**Companion Paper**: "Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy" (Prediction B, in `../mdpi_template/` and `../manuscript.tex`)

## Core Claim (Diagnose → Ablate → Replace → Confirm)

1. **Diagnose**: Under the original target-anchored action (A2 full), lor4d and KR are nearly tied (36:35 wins). 4D does *not* win.
2. **Ablate**: Systematic ablation identifies a single term (g_dim) as the dominant source of low-dimensional bias.
3. **Replace**: Replacing g_dim with non-target-anchored g_con yields systematic 4D dominance across the 14×7 base grid (92/98 wins).
4. **Confirm**: Margin shows strong upward trend (+7 at N=20 → +57 at N=72) with minor finite-size dips at N=44 and N=64. Seed-robust: 100% win rate at N=68 and N=72 across 3 seeds each (42/42).

## Files

| File | Description |
|------|-------------|
| `manuscript_entropy.tex` | MDPI Entropy 模板格式论文（推荐投稿用） |
| `manuscript.tex` | 标准 LaTeX 格式论文（通用版） |
| `MANUSCRIPT_DRAFT.md` | Markdown 草稿（便于阅读和编辑） |
| `cover_letter_entropy_mdpi.md` | Cover Letter + 投稿操作指南 |
| `Definitions/` | MDPI LaTeX 模板定义文件 |
| `manuscript_figures/` | 论文图片（fig1_margin_of_victory.png, fig2_winner_phase_comparison.png） |

## Data Sources

实验数据位于 `../outputs_exploratory/` 下：

- `prediction_a_dim_replacement_sp8/` — N=20~40 基础网格（K=8, exact computation）
- `prediction_a_dim_replacement_n44_n48/` — N=44, N=48 扩展数据
- `prediction_a_n{52,56,60,64,68,72}_mixed/` — 各 N 值实验数据（K=4, SIS regime）
- `prediction_a_geometric_ablation/` — 系统化消融实验（N=20~36）
- `prediction_a_seed_sensitivity_n68/` — N=68 种子敏感性分析
- `prediction_a_seed_sensitivity_n72/` — N=72 种子敏感性分析
- `prediction_a_margin_summary/` — 全 N 值 margin 汇总
- `prediction_a_phase_summary/` — winner phase 统计

## 编译

```bash
# MDPI 模板版
cd preA
pdflatex manuscript_entropy.tex

# 标准版
pdflatex manuscript.tex
```

## 投稿顺序

1. 先投 Prediction B（companion paper）
2. 获得稿件编号后，在 Prediction A 中引用
3. 或两篇同时投稿，注明 companion 关系

## 维护说明

**Source of Truth**: `manuscript.tex`（标准版）为内容主稿，`manuscript_entropy.tex`（MDPI版）通过同步更新保持一致。修改内容时应先改 `manuscript.tex`，再同步到 `manuscript_entropy.tex`。
