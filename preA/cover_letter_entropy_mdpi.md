# Cover Letter — MDPI *Entropy* 投稿 (Prediction A)

> MDPI *Entropy* 投稿通过在线系统 SuSy 完成：https://susy.mdpi.com/  
> 下方 Cover Letter 粘贴至投稿系统的 Cover Letter 栏。  

---

**Date:** March 12, 2026

**To:** The Academic Editor, *Entropy* (MDPI)

**Re:** Submission of manuscript entitled "Dimensional Selection Without Dimensional Priors: 4D Lorentzian Dominance in Finite Causal Poset Ensembles Under Consistency-Based Actions"

---

Dear Editor,

I am pleased to submit the enclosed manuscript for consideration in *Entropy*. This paper is a companion to a concurrently submitted work "Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy" (Prediction B), which establishes the existence of a bounded geometric phase transition in finite causal posets and introduces a non-target-anchored consistency penalty as a replacement for dimension-specific priors.

**Motivation and scope fit.** This paper asks a second-order question arising from the companion work: among Lorentzian-like structures of different dimensionality, *which dimension does a consistency-based action prefer?* The central objects are combinatorial entropy $H(P) = \log |L(P)|$ of finite posets and the competition between entropy maximization and structural consistency penalties — a free energy competition that falls squarely within *Entropy's* scope on information-theoretic and statistical-mechanical foundations.

**Summary of contributions.** The paper follows a "diagnose → ablate → replace → confirm" narrative:

1. **Diagnosis.** Under the original target-anchored action ($A_2$ full), a 4D Lorentzian-like family and KR posets are locked in a near-perfect tie (36 vs 35 wins out of 98 configurations). The 4D family does *not* achieve stable dominance — the negative result is presented honestly.

2. **Ablation.** Systematic removal of action components identifies a single term ($g_{\text{dim}}$, the target-anchored dimension proxy) as the dominant source of low-dimensional bias. This parallels the companion paper's finding that the same term is the critical backbone for the 2D phase transition.

3. **Replacement.** Replacing $g_{\text{dim}}$ with the non-target-anchored consistency penalty $g_{\text{con}}$ yields systematic 4D dominance across the tested grid (92/98 configurations on a $14\,N \times 7\,\gamma$ grid), with the 4D signal surviving inside a full geometric action package.

4. **Scaling and robustness.** The margin of victory shows a strong upward trend from +7 ($N = 20$) to +57 ($N = 72$), with finite-size dips at $N = 44$ and $N = 64$. Three independent generator seeds at $N = 68$ and $N = 72$ confirm 100% win rate (42/42 configurations) under consistency variants, versus 43% under the original action.

**Novelty.** This is the first finite-size numerical evidence that dimensional consistency constraints, without dimension-specific priors, naturally favor higher-dimensional Lorentzian structures in discrete causal order ensembles — a nontrivial step toward turning "why 3+1?" from a metaphysical presupposition into a quantitatively testable structural question.

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

# MDPI *Entropy* 投稿操作步骤 (Prediction A)

## 1. 在线系统
- **投稿入口**：https://susy.mdpi.com/
- 注册/登录后，选择 **Submit to Journal** → **Entropy**

## 2. 投稿流程
1. **Journal**: Entropy
2. **Section**: 选 **"Complexity"** 或 **"Statistical Physics"**（推荐前者）
3. **Special Issue**: 若无匹配的特刊，选 **"Regular Submission"**
4. **Article Type**: **Article**
5. **Title**: `Dimensional Selection Without Dimensional Priors: 4D Lorentzian Dominance in Finite Causal Poset Ensembles Under Consistency-Based Actions`
6. **Authors**: Gang Zhang（填 Independent Researcher，邮箱 unicome@gmail.com）
7. **Abstract**: 粘贴手稿中的 abstract
8. **Keywords**: `causal sets; posets; phase transition; dimensional consistency; dimension selection; higher-dimensional geometry; linear extensions; discrete quantum gravity`
9. **Cover Letter**: 粘贴上方英文内容
10. **Manuscript**: 上传编译好的 **PDF**
11. **Supplementary Materials**: 可选上传代码说明文档
12. **Related Manuscript**: 注明与 Prediction B companion paper 的关系
13. **Suggest Reviewers**（可选）:
    - Sumati Surya (Raman Research Institute) — causal set theory
    - Lisa Glaser (University of Vienna) — causal set numerical simulation
    - William Cunningham — discrete quantum gravity
    - Graham Brightwell (London School of Economics) — combinatorics of posets
14. **Review** → **Confirm** → **Submit**

## 3. 关于 APC
- *Entropy* 是 **开放获取** 期刊，录用后需支付 APC（约 **2600 CHF / ≈2800 USD**）
- 可申请 **fee waiver**

## 4. 投稿顺序建议
- **推荐**：先投 Prediction B（companion paper），获得稿件编号后在 Prediction A 的 cover letter 和正文中引用
- 或两篇同时投稿，在 cover letter 中说明是 companion papers

## 5. 与 Prediction B 的关系
- Prediction B = "Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy"
- Prediction A（本文）= 基于 Prediction B 的非循环替换结果做出的预测性检验
- 两篇可以独立阅读，但 Prediction A 引用 Prediction B 作为基础
