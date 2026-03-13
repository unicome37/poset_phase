# Prediction C — Entropy (MDPI) 投稿准备清单

## 投稿通道
- 期刊：MDPI *Entropy* (ISSN 1099-4300)
- 系统：SuSy https://susy.mdpi.com/
- 影响因子：2.0 | CiteScore 5.2 | 首次决定 22 天
- 首次投稿接受 Free Format（无需套模板）

## 清单

### 1. 稿件文件
- [x] `preC/MANUSCRIPT_PredictionC_Full.md` — 完整稿件 (788 lines / 8,720 words / 14 tables)
- [x] 标题："Hierarchy Depth Observables Predict Combinatorial Entropy in Finite Causal Posets: A Three-Tier Correlational Study"
- [x] 摘要 ≈ 190 words（Entropy 要求 ~200 words）
- [x] 关键词：causal posets, hierarchy integration, combinatorial entropy, linear extensions, coarse-graining stability, Simpson's Paradox, causal set theory

### 2. Cover Letter
- [x] `preC/CoverLetter_PredictionC.md`
- [ ] 填写日期 [INSERT DATE]
- [ ] 确认推荐审稿人信息

### 3. 参考文献
- [x] [1] Prediction B — Zhang, G. *Entropy* 2026, submitted.
- [x] [2] Prediction A — Zhang, G. *Entropy* 2026, submitted.
- [x] [3] Surya 2012 — 格式符合 MDPI 样式
- [x] [4] Brightwell & Winkler 1991 — 格式符合 MDPI 样式
- [x] [5] Pratt & Gibbons 1981 — 格式符合 MDPI 样式

### 4. Data Availability Statement
粘贴至投稿系统：

> **Data Availability Statement:** All data generated during this study are available in the GitHub repository https://github.com/unicome37/poset_phase, archived at Zenodo (DOI: [10.5281/zenodo.18980657](https://doi.org/10.5281/zenodo.18980657)). Raw output CSV files are in the `outputs_exploratory/` directory; configuration files for full reproduction are provided as `config_prediction_c_*.yaml`.

### 5. Author Contributions (CRediT)
粘贴至投稿系统：

> **Author Contributions:** Conceptualization, G.Z.; methodology, G.Z.; software, G.Z.; validation, G.Z.; formal analysis, G.Z.; investigation, G.Z.; data curation, G.Z.; writing—original draft preparation, G.Z.; writing—review and editing, G.Z.

### 6. Funding Statement
> **Funding:** This research received no external funding.

### 7. Conflicts of Interest
> **Conflicts of Interest:** The author declares no conflicts of interest.

### 8. AI Disclosure
在稿件中已声明。投稿时同时在系统中勾选并声明：

> Large language model assistants (Claude, Anthropic) contributed to code implementation, data analysis, and manuscript drafting. All scientific hypotheses, experimental design decisions, and result interpretation were solely the responsibility of the human author.

### 9. GitHub + Zenodo

#### 当前状态
- GitHub 仓库：https://github.com/unicome37/poset_phase
- 现有 DOI（v2.1.0）：10.5281/zenodo.18980657（概念 DOI，自动更新）
- Zenodo–GitHub 联动已启用

#### 发布步骤
1. `git add` 所有 Prediction C 新文件
2. `git commit -m "feat: add Prediction C — Hierarchy Depth Predicts Entropy (v3.0.0)"`
3. `git tag v3.0.0`
4. `git push origin main --tags`
5. GitHub → Releases → New Release → Tag: v3.0.0 → 粘贴 RELEASE_NOTES_v3.0.0.md
6. Zenodo 自动归档 → 新版本 DOI 自动挂载在概念 DOI 下
7. 确认 DOI badge 链接正确

#### DOI 证明
- 概念 DOI（始终指向最新版）：https://doi.org/10.5281/zenodo.18980657
- v3.0.0 版本 DOI：发布后自动生成，格式 10.5281/zenodo.XXXXXXX
- 投稿中使用概念 DOI 即可

### 10. 三篇论文投稿顺序建议
1. **Prediction B** 先投（基础：建立竞争窗口和非循环性）
2. **Prediction A** 次投（延伸：维度选择）
3. **Prediction C** 最后投（机制：解释 "为什么"）

或三篇同时投稿（不同编辑/审稿人），在 cover letter 中互相引用为 "submitted"。

---

## 投稿后检查清单
- [ ] 收到 SuSy 确认邮件
- [ ] 记录 Manuscript ID
- [ ] 三篇论文的 Manuscript ID 互相补充到引用中（如果同时投稿）
- [ ] Zenodo v3.0.0 DOI 生成后更新稿件中的 DOI 链接
