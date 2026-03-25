# Cover Letter — Foundations of Physics 投稿

> 请将下方英文 Cover Letter 粘贴到 Springer 在线投稿系统的 Cover Letter 栏中。  
> 标记 `[...]` 处需替换为真实信息。

---

**Date:** March 12, 2026

**To:** The Editors, *Foundations of Physics*

**Re:** Submission of manuscript "Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy"

---

Dear Editors,

We are pleased to submit the above manuscript for consideration in *Foundations of Physics*.

**Summary.** A longstanding obstacle in causal set theory (CST) is the *entropy catastrophe*: the Kleitman–Rothschild theorem implies that almost all finite posets are 3-layered and non-geometric, dominating any entropy-weighted partition function and preventing geometric phases from emerging. While Carlip (2017) reviewed the evidence for dimensional reduction in quantum gravity approaches and Loomis & Carlip (2018) analytically demonstrated the suppression of non-manifold-like sets in the causal set path integral, no prior work has established an explicit, finite-size phase transition threshold computed with exact entropy across a systematic range of poset sizes.

**Key contributions.** In this work we:

1. Construct a 7-family poset ensemble spanning structurally diverse causal orders and compute *exact* linear extension counts (via antichain-based dynamic programming) for N = 10–44 — far beyond the reach of brute-force enumeration.

2. Demonstrate that under a decomposable action with neutral and geometric penalty terms, the critical geometric coupling γ\_c for Lorentzian-like 2D posets to overtake KR posets remains bounded and O(1) across all tested sizes — the first such finite-size evidence.

3. Perform systematic ablation of seven geometric sub-terms, identifying a minimal backbone of two constraints: global layer-shape balance and dimensional scale consistency.

4. Show that the dimension-dependent term can be replaced by a *non-target-anchored* alternative (penalizing local–global dimension inconsistency rather than deviation from d = 2), with γ\_c preserved at the same order of magnitude — directly addressing circularity concerns.

**Relevance to Foundations of Physics.** The paper addresses a core foundational question in discrete quantum gravity: whether geometric order can emerge from structural self-consistency rather than being imposed as a prior. This aligns with the journal's focus on conceptual and mathematical foundations of physics, particularly in quantum gravity and spacetime structure.

**Declarations.**
- This manuscript has not been published elsewhere and is not under consideration by another journal.
- Code and data are publicly available at: https://github.com/unicome37/poset_phase (archived at Zenodo: https://doi.org/10.5281/zenodo.18963421).
- AI assistants (Claude, Anthropic) contributed to code implementation, data analysis, and manuscript drafting. All scientific decisions, theoretical framework design, and result interpretation were made by the human author. This is disclosed in the Acknowledgments section.

**Suggested reviewers** (optional — fill in or remove):
1. Sumati Surya, Raman Research Institute, ssurya@rri.res.in — expertise in causal set theory
2. Graham Brightwell, London School of Economics, g.r.brightwell@lse.ac.uk — expertise in combinatorics / poset enumeration
3. Lisa Glaser, University of Vienna, lisa.glaser@univie.ac.at — expertise in discrete quantum gravity

Thank you for your consideration. We look forward to your response.

Sincerely,

Gang Zhang  
Independent Researcher  
unicome@gmail.com

---

# 投稿操作步骤（Foundations of Physics, Springer）

## 1. 在线投稿入口
- 网址：https://www.editorialmanager.com/fophy/
- 如无账号，先注册（Register）

## 2. 投稿流程
1. **Login** → **Submit New Manuscript**
2. **Article Type**: 选 "Original Research"
3. **Title**: `Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy`
4. **Abstract**: 直接粘贴手稿中的 abstract
5. **Keywords**: `causal sets, posets, phase transition, linear extensions, Kleitman–Rothschild, discrete quantum gravity`
6. **Authors**: 填写真实姓名、机构、邮箱
7. **Cover Letter**: 粘贴上方英文内容
8. **Manuscript File**: 上传 `manuscript_foundphys.tex`（或编译好的 PDF）
9. **Figures**: 上传 `manuscript_figures/` 中的三张 PDF 图
10. **Supplementary**: 可选上传代码仓库说明
11. **Review** → **Submit**

## 3. arXiv 同步（可选，建议优先）
- 先投 arXiv 拿到预印本编号，投稿时在 Cover Letter 中补充 arXiv ID
- arXiv 分类建议：`gr-qc`（主）、`hep-th`（副）、`math-ph`（副）
- arXiv 提交页：https://arxiv.org/submit

---

# 投稿邮件模板（备用，如需邮件联系编辑）

**Subject:** Manuscript Submission: Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy

Dear Editorial Office,

Please find attached our manuscript entitled "Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy" for consideration in Foundations of Physics.

The paper presents the first finite-size evidence of a bounded geometric phase transition threshold in causal poset ensembles, computed with exact linear extension entropy for N = 10–44. Through systematic ablation, we identify a minimal non-circular backbone of structural constraints sufficient for geometric phase emergence — directly addressing the Kleitman–Rothschild entropy catastrophe in causal set theory.

We believe this work is well suited for Foundations of Physics given its focus on foundational questions in discrete quantum gravity and spacetime emergence.

A detailed cover letter is included with the submission. Code and data are publicly available on GitHub and archived on Zenodo.

We look forward to your consideration.

Best regards,  
Gang Zhang  
unicome@gmail.com
