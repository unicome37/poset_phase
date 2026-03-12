# Cover Letter — MDPI *Entropy* 投稿

> MDPI *Entropy* 投稿通过在线系统 SuSy 完成：https://susy.mdpi.com/  
> 下方 Cover Letter 粘贴至投稿系统的 Cover Letter 栏。  
> 标记 `[...]` 处替换为真实信息。

---

**Date:** March 12, 2026

**To:** The Academic Editor, *Entropy* (MDPI)

**Re:** Submission of manuscript entitled "Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy"

---

Dear Editor,

We are pleased to submit the enclosed manuscript for consideration in *Entropy*.

**Motivation and scope fit.** This paper studies how combinatorial entropy — quantified by exact linear extension counts of finite partially ordered sets — competes with structural order constraints to determine which causal configurations emerge as dominant phases. The central mathematical object is $H(P) = \log |L(P)|$, the logarithm of the number of linear extensions of a poset $P$, which serves as the discrete analogue of thermodynamic entropy in this setting. The interplay between entropy maximization and geometric penalty terms directly mirrors the free energy competition familiar in statistical mechanics. We believe this falls squarely within *Entropy's* scope on information-theoretic and statistical-mechanical foundations.

**Summary of contributions.**

1. **Exact entropy computation.** We construct a 7-family poset ensemble and compute exact linear extension counts via antichain-based dynamic programming for posets with $N = 10$–$44$ elements. This exact approach avoids the systematic biases inherent in Monte Carlo approximations and reveals that computational cost itself varies dramatically across structural types — a structural fact about the posets, not an algorithmic artifact.

2. **Bounded phase transition threshold.** Under a decomposable action combining entropic and geometric terms, we show that the critical coupling $\gamma_c$ — the minimum geometric penalty strength needed for Lorentzian-like posets to overtake the entropically dominant Kleitman–Rothschild (KR) structures — remains bounded and $O(1)$ across all tested sizes. This constitutes the first finite-size evidence of such a bounded threshold.

3. **Systematic ablation.** Ablation of seven individual geometric sub-terms identifies a minimal two-component backbone (global layer-shape balance + dimensional scale information) necessary and sufficient for the phase transition.

4. **Non-circular replacement.** The dimension-dependent backbone component can be replaced by a non-target-anchored alternative — penalizing local–global dimension *inconsistency* rather than deviation from a specific target $d = 2$ — with the phase transition preserved at the same order of magnitude. This directly addresses the concern that geometric penalties might circularly "reward looking geometric."

**Novelty.** While prior work (Carlip 2017; Loomis & Carlip 2018) has argued qualitatively and via MCMC for KR suppression in causal set models, no previous study has provided an explicit $\gamma_c(N)$ curve with exact entropy across a systematic range of sizes, combined with full ablation and non-circular replacement analysis.

**Declarations.**
- This manuscript is original, has not been published elsewhere, and is not under consideration by another journal.
- Code and data: https://github.com/unicome37/poset_phase (Zenodo archive: https://doi.org/10.5281/zenodo.18963421).
- AI Assistance: Large language model assistants (Claude, Anthropic) contributed to code implementation, data analysis, and manuscript drafting. All scientific decisions were made by the human author. This is fully disclosed in the manuscript.
- No external funding was received.
- The author declares no conflicts of interest.

Thank you for considering this work. We look forward to your editorial feedback.

Sincerely,

Gang Zhang  
Independent Researcher  
unicome@gmail.com

---

# MDPI *Entropy* 投稿操作步骤

## 1. 在线系统
- **投稿入口**：https://susy.mdpi.com/
- 注册/登录后，选择 **Submit to Journal** → **Entropy**

## 2. 投稿流程
1. **Journal**: Entropy
2. **Section**: 选 **"Complexity"** 或 **"Statistical Physics"**（推荐前者——因为论文涉及因果结构的涌现、相变、组合熵竞争）
3. **Special Issue**: 若无匹配的特刊，选 **"Regular Submission"**
4. **Article Type**: **Article**
5. **Title**: `Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy`
- **Authors**: Gang Zhang（填 Independent Researcher，邮箱 unicome@gmail.com）
7. **Abstract**: 粘贴手稿中的 abstract
8. **Keywords**: `causal sets; posets; phase transition; linear extensions; Kleitman–Rothschild; discrete quantum gravity; Lorentzian geometry; dimensional consistency`
9. **Cover Letter**: 粘贴上方英文内容
10. **Manuscript**: 上传编译好的 **PDF** 或 **Word**（MDPI 接受 PDF 初投，录用后再转 LaTeX 模板）
    - 推荐上传从 `mdpi_template/manuscript_entropy.tex` 编译的 PDF
11. **Supplementary Materials**: 可选上传代码说明文档
12. **Suggest Reviewers**（可选）:
    - Sumati Surya (Raman Research Institute) — causal set theory
    - Lisa Glaser (University of Vienna) — causal set numerical simulation
    - William Cunningham — discrete quantum gravity
    - Graham Brightwell (London School of Economics) — combinatorics of posets
13. **Review** → **Confirm** → **Submit**

## 3. 关于 MDPI 的 APC（文章处理费）
- *Entropy* 是 **开放获取** 期刊，录用后需支付 APC（目前约 **2600 CHF / ≈2800 USD**）
- 可在投稿时申请 **fee waiver**（根据资金情况），MDPI 对来自发展中国家或无资助的独立研究者有减免政策
- 在投稿过程最后一步会看到 APC 相关声明，确认即可

## 4. 审稿周期
- MDPI 以快速审稿著称，通常 **2-4 周** 内返回第一轮审稿意见
- 修改稿通常 **1-2 周** 内返回

## 5. 关于 arXiv
- 虽然暂时无法投 arXiv，但 MDPI 许可作者在 **录用后** 将 preprint 版本上传到 arXiv 或其他预印本服务器
- 也可考虑使用 **Research Square**、**SSRN** 或 **Zenodo** 作为预印本替代平台（这些不需要推荐资格）

---

# 备选：免推荐的预印本平台

| 平台 | 审核 | 费用 | DOI | 说明 |
|------|------|------|-----|------|
| **Zenodo** | 无同行评审 | 免费 | ✅ | 你已有 Zenodo 存档，可直接关联 |
| **Research Square** | 轻度审核 | 免费 | ✅ | Springer Nature 旗下，接受度好 |
| **SSRN** | 无同行评审 | 免费 | ✅ | Elsevier 旗下，偏社科但也收物理 |
| **OSF Preprints** | 轻度审核 | 免费 | ✅ | COS 运营，学术界认可 |
| **Preprints.org** | MDPI 旗下 | 免费 | ✅ | 与 MDPI 投稿流程天然兼容 |

**推荐**：先在 **Preprints.org**（MDPI 旗下）上传预印本 → 获取 DOI → 再通过 SuSy 正式投稿 *Entropy*。这样既有预印本可引用，又与 MDPI 生态无缝衔接。
