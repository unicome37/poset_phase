# Prediction A — Dimension-Agnostic Geometric Dominance

## Paper

**Title**: Dimension-Agnostic Geometric Dominance in Finite Causal Posets: Higher-Dimensional Lorentzian Structures Emerge Under Consistency-Based Actions

**Author**: Gang Zhang (Independent Researcher)

**Target Journal**: MDPI *Entropy*

**Companion Paper**: "Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy" (Prediction B, in `../mdpi_template/` and `../manuscript.tex`)

## Core Claim

If the geometric action rewards dimensional *consistency* (not proximity to d=2), then higher-dimensional Lorentzian-like posets should dominate lower-dimensional ones. This is confirmed:

- **A2 full** (target d=2): lor4d vs KR nearly tied (36:35 wins)
- **A2 consistency**: lor4d achieves unconditional dominance (92/98 wins)
- **Margin grows monotonically**: +7 at N=20 → +57 at N=72
- **Seed-robust**: 100% win rate at N=72 across 3 seeds (21/21)

## Files

| File | Description |
|------|-------------|
| `manuscript_entropy.tex` | MDPI Entropy 模板格式论文（推荐投稿用） |
| `manuscript.tex` | 标准 LaTeX 格式论文（通用版） |
| `MANUSCRIPT_DRAFT.md` | Markdown 草稿（便于阅读和编辑） |
| `cover_letter_entropy_mdpi.md` | Cover Letter + 投稿操作指南 |
| `Definitions/` | MDPI LaTeX 模板定义文件 |
| `manuscript_figures/` | 论文图片（待生成） |

## Data Sources

实验数据位于 `../outputs_exploratory/` 下：

- `prediction_a_n72_mixed/` — N=72 主实验数据
- `prediction_a_n{52,56,60,64,68}_mixed/` — 各 N 值实验数据
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
