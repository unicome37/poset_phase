# Cover Letter — MDPI *Entropy* 投稿 (Prediction C)

> MDPI *Entropy* 投稿通过在线系统 SuSy 完成：https://susy.mdpi.com/  
> 下方 Cover Letter 粘贴至投稿系统的 Cover Letter 栏。  

---

**Date:** March 13, 2026

**To:** The Academic Editor, *Entropy* (MDPI)

**Re:** Submission of manuscript entitled "Hierarchy Depth Observables Predict Combinatorial Entropy in Finite Causal Posets: A Three-Tier Correlational Study"

---

Dear Editor,

I am pleased to submit the enclosed manuscript for consideration in *Entropy*. This paper is the third in a series of three companion papers: Prediction B ("Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy") and Prediction A ("Dimensional Selection Without Dimensional Priors"), both concurrently submitted to *Entropy*. The present paper (Prediction C) addresses the structural mechanism question: *why* do certain poset families achieve systematically lower combinatorial entropy?

**Motivation and scope fit.** The central objects are combinatorial entropy $H(P) = \log |L(P)|$ of finite causal posets and a pre-registered Hierarchy Integration Index (HII) — a composite z-score of five structural depth observables. The interplay between structural order and combinatorial counting connects directly to *Entropy's* scope on information-theoretic and statistical-mechanical foundations of discrete systems.

**Summary of contributions.**

1. **Three-tier validation design.** We test the HII–entropy hypothesis across three independent statistical tiers: (i) all-family partial correlation (8 families, $N = 10$–$16$, 320 samples with exact entropy), (ii) matched-pair $\Delta$-analysis between Lor2D and MLR families (46 pairs, $N = 30$–$56$), and (iii) coarse-graining identity-stability linkage (92 samples). Each tier addresses different potential confounds.

2. **Strong, stable correlation.** At fixed $N$, deeper hierarchy — principally `layer_count` and `mean_layer_gap` — correlates negatively with $\log H$: matched-pair $r = -0.834$, stable across three levels of filter stringency (variation $< 0.005$). The association extends to coarse-graining stability ($r = -0.874$, classifier-contingent).

3. **Simpson's Paradox diagnosis.** The naïve cross-$N$ correlation is positive ($r = +0.336$) and reverses sign only after controlling for $N$ ($r = -0.578$). We provide a complete decomposition of this paradox — a methodological finding of independent value for researchers working with size-heterogeneous discrete structure data.

4. **Component decomposition.** The five-component HII composite never exceeds its best single constituent (`layer_count`), revealing that the signal is carried by a two-component depth pair.

5. **Hard limitations stated.** Family-specific HII reversals (KR-like, absolute-layered), classifier-dependent Tier 3, and weak near-wall statistical power at $N \geq 52$ are documented as explicit boundaries.

**Novelty.** No previous study of causal set ensembles has: (a) defined a pre-registered composite structural index for hierarchy depth; (b) validated it across three independent statistical designs spanning two disjoint $N$ ranges; (c) diagnosed a Simpson's Paradox in the structural–entropy relationship; or (d) extended the association chain from hierarchy to coarse-graining stability.

**Declarations.**
- This manuscript is original, has not been published elsewhere, and is not under consideration by another journal.
- Code and data: https://github.com/unicome37/poset_phase (Zenodo archive: https://doi.org/10.5281/zenodo.18980657).
- AI Assistance: Large language model assistants (Claude, Anthropic) contributed to code implementation, data analysis, and manuscript drafting. All scientific decisions were made by the human author. This is fully disclosed in the manuscript.
- No external funding was received.
- The author declares no conflicts of interest.

Thank you for considering this work. I look forward to your editorial feedback.

Sincerely,

Gang Zhang  
Independent Researcher  
unicome@gmail.com

---

# MDPI *Entropy* 投稿操作步骤 (Prediction C)

## 1. 在线系统
- **投稿入口**：https://susy.mdpi.com/
- 注册/登录后，选择 **Submit to Journal** → **Entropy**

## 2. 投稿流程
1. **Journal**: Entropy
2. **Section**: 选 **"Complexity"** 或 **"Statistical Physics"**（推荐前者）
3. **Special Issue**: 若无匹配的特刊，选 **"Regular Submission"**
4. **Article Type**: **Article**
5. **Title**: `Hierarchy Depth Observables Predict Combinatorial Entropy in Finite Causal Posets: A Three-Tier Correlational Study`
6. **Authors**: Gang Zhang（填 Independent Researcher，邮箱 unicome@gmail.com）
7. **Abstract**: 粘贴手稿中的 abstract
8. **Keywords**: `causal posets; hierarchy integration; combinatorial entropy; linear extensions; coarse-graining stability; Simpson's Paradox; causal set theory`
9. **Cover Letter**: 粘贴上方英文内容
10. **Manuscript**: 上传编译好的 **PDF**
11. **Supplementary Materials**: 可选上传代码说明文档
12. **Related Manuscript**: 注明与 Prediction A 和 Prediction B companion papers 的关系
13. **Suggest Reviewers**（可选）:
    - Sumati Surya (Raman Research Institute) — causal set quantum gravity
    - Graham Brightwell (London School of Economics) — linear extension counting, poset combinatorics
    - Seth Major (Hamilton College) — discrete quantum gravity
    - Lisa Glaser (University of Vienna) — causal set numerical simulation
14. **Review** → **Confirm** → **Submit**

## 3. 关于 APC
- *Entropy* 是 **开放获取** 期刊，录用后需支付 APC（约 **2600 CHF / ≈2800 USD**）
- 可申请 **fee waiver**

## 4. 投稿顺序建议
- **推荐**：先投 Prediction B（基础论文），再投 Prediction A，最后投 Prediction C
- 或三篇同时投稿，在各自 cover letter 中说明是 companion papers

## 5. 三篇论文的关系
- **Prediction B** = "Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy"（基础）
- **Prediction A** = "Dimensional Selection Without Dimensional Priors"（维度选择）
- **Prediction C**（本文）= "Hierarchy Depth Observables Predict Combinatorial Entropy"（结构机制）
- 三篇可以独立阅读，但 C 引用 A 和 B 作为基础，提供结构层面的解释
