# Prediction C — Hierarchy Depth Mechanism

## Paper

**Title**: Hierarchy Depth Observables Predict Combinatorial Entropy in Finite Causal Posets: A Three-Tier Correlational Study

**Author**: Gang Zhang (Independent Researcher)

**Target Journal**: MDPI *Entropy*

**Companion Papers**:
- "Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy" (Prediction B)
- "Dimensional Selection Without Dimensional Priors" (Prediction A, in `../preA/`)

## Core Claim (Three-Tier Correlational Study)

1. **Tier 1**: All-family partial correlation (8 families, N=10–16, 320 samples) — r_partial(HII, log H | N) = −0.578
2. **Tier 2**: Matched-pair Δ-analysis (Lor2D vs MLR, 46 pairs, N=30–56) — r = −0.834, stable across 3 filter stringencies (variation < 0.005)
3. **Tier 3**: CG identity-stability linkage (92 samples) — layer_count predicts switch rate at r = −0.874
4. **Simpson's Paradox**: Naïve cross-N correlation is positive (+0.336); controlling for N reverses sign to −0.578
5. **Component Decomposition**: layer_count is the single strongest predictor; HII composite never exceeds its best constituent

## Files

| File | Description |
|------|-------------|
| `manuscript_predictionC.tex` | MDPI Entropy 模板格式论文（投稿用） |
| `MANUSCRIPT_PredictionC_Full.md` | Markdown 完整稿（便于阅读和编辑） |
| `cover_letter_entropy_mdpi.md` | Cover Letter + 投稿操作指南 |
| `CoverLetter_PredictionC.md` | Cover Letter 早期版本 |
| `SUBMISSION_CHECKLIST.md` | 投稿准备清单 |
| `Definitions/` | MDPI LaTeX 模板定义文件 |
| `manuscript_predictionC.pdf` | 编译好的 PDF（投稿用） |
| `manuscript_figures/` | 论文图表（fig1–fig3，png+pdf 格式） |
| `generate_figures.py` | 图表生成脚本（基于实际数据，可复现） |
| `MANUSCRIPT_Section*.md` | 分段 Markdown 稿件 |
| `PredictionC_Section*.md` | 分段 Markdown 稿件（早期版本） |

## Relation to Other Papers

| Paper | Core claim | Relation to C |
|-------|-----------|---------------|
| **Prediction A** | g_con selects higher-d Lorentzian; 4D wins under consistency | C provides the mechanism: deeper hierarchy → lower entropy |
| **Prediction B** | Bounded γ_c for Lor2D vs KR; non-circular dim_consistency replacement | C explains *why* the structural advantage exists |

### Logical Dependency

- **Logically independent**: C's HII–log_H correlation can be verified without A or B being true.
- **Semantically layered**: If B+C both hold, C elevates B from "numerical coincidence" to "structurally supported regularity."
- **Three-prediction tower**: A (dimension selection) → B (self-consistency) → C (structural mechanism). Each layer adds explanatory depth.

## Data Sources

实验数据位于 `../outputs_exploratory/` 下：

- `prediction_c_comprehensive/tier1_*.csv` — Tier 1 全家族精确熵数据
- `prediction_c_comprehensive/tier2_*.csv` — Tier 2 匹配对数据
- `prediction_c_comprehensive/tier3_*.csv` — Tier 3 粗粒化稳定性数据
- Full validation report: `PredictionC_精确表述与验证报告.md`

## 编译

```bash
cd preC
python generate_figures.py          # 生成 manuscript_figures/fig1–fig3
pdflatex manuscript_predictionC.tex   # 首次编译
pdflatex manuscript_predictionC.tex   # 二次编译（解析引用）
```

## 投稿顺序

1. 先投 Prediction B（基础论文）
2. 再投 Prediction A（维度选择）
3. 最后投 Prediction C（本文，结构机制）
4. 或三篇同时投稿，注明 companion 关系

## 维护说明

**Source of Truth**: `MANUSCRIPT_PredictionC_Full.md` 为内容主稿，`manuscript_predictionC.tex`（MDPI版）通过同步更新保持一致。修改内容时应先改 Markdown 稿，再同步到 tex 文件。
