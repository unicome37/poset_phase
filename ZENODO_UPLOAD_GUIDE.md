# Zenodo 预印本上传指南 — 三篇论文时间戳

## 已编译 PDF 文件位置

| 论文 | PDF 路径 | 页数 |
|---|---|---|
| Prediction B | `D:\Kiro\理论体系\poset_phase\manuscript_foundphys.pdf` | 7pp |
| Prediction A | `D:\Kiro\理论体系\poset_phase\preA\manuscript.pdf` | 9pp |
| Prediction C | `D:\Kiro\理论体系\poset_phase\preC\manuscript_predictionC.pdf` | 14pp |

---

## 操作步骤

### 1. 登录 Zenodo
- 打开 https://zenodo.org/uploads/new
- 用之前上传代码仓库的同一账号登录

### 2. 上传 Prediction B（第一篇）

**Files**: 上传 `manuscript_foundphys.pdf`

**Metadata**:
- **Resource type**: `Preprint`
- **Title**: `Bounded Geometric Phase Transition in Finite Causal Posets with Exact Entropy`
- **Authors**: Gang Zhang (Independent Researcher, unicome@gmail.com, ORCID 如有)
- **Description** (粘贴以下):

```
We construct a 7-family poset ensemble and compute exact linear extension counts for N = 10–44 elements to study action-weighted competition between Lorentzian-like geometric structures and Kleitman–Rothschild (KR) high-entropy posets. Under a decomposable action combining neutral and geometric penalty terms, the critical geometric coupling γ_c remains O(1) across all tested sizes — the first finite-size evidence of a bounded phase transition threshold in this setting. Systematic ablation identifies a minimal backbone of two constraints: global layer-shape balance and dimensional scale consistency, the latter replaceable by a non-target-anchored variant with γ_c preserved.
```

- **Keywords**: `causal sets`, `posets`, `phase transition`, `linear extensions`, `Kleitman-Rothschild`, `discrete quantum gravity`
- **License**: Creative Commons Attribution 4.0 International (CC BY 4.0)
- **Related identifiers**: 
  - `is supplemented by` → `https://doi.org/10.5281/zenodo.18980657` (代码仓库)

### 3. 上传 Prediction A（第二篇）

**Files**: 上传 `manuscript.pdf` (来自 preA/)

**Metadata**:
- **Resource type**: `Preprint`
- **Title**: `Dimensional Selection Without Dimensional Priors: 4D Lorentzian Dominance in Finite Causal Poset Ensembles Under Consistency-Based Actions`
- **Authors**: Gang Zhang (Independent Researcher)
- **Description** (粘贴以下):

```
We extend a finite causal poset ensemble to include 2D, 3D, and 4D Lorentzian-like families and compare them under original target-anchored and non-target-anchored consistency-based actions across N = 20–72. Systematic ablation localizes low-dimensional bias to a single term (g_dim). Replacing it with a consistency penalty yields systematic 4D dominance (92/98 configurations), with growing margin and 100% seed-robust win rate at N = 68–72 (42/42). These results demonstrate that dimensional selection without dimensional priors is achievable in discrete causal order ensembles.
```

- **Keywords**: `causal sets`, `posets`, `dimensional consistency`, `dimension selection`, `linear extensions`, `discrete quantum gravity`
- **License**: CC BY 4.0
- **Related identifiers**: 
  - `is supplemented by` → `https://doi.org/10.5281/zenodo.18980657`
  - `references` → Prediction B 的 Zenodo DOI（上传后获得）

### 4. 上传 Prediction C（第三篇）

**Files**: 上传 `manuscript_predictionC.pdf` (来自 preC/)

**Metadata**:
- **Resource type**: `Preprint`
- **Title**: `Hierarchy Depth Observables Predict Combinatorial Entropy in Finite Causal Posets: A Three-Tier Correlational Study`
- **Authors**: Gang Zhang (Independent Researcher)
- **Description** (粘贴以下):

```
We define a Hierarchy Integration Index (HII) — a composite z-score of five structural depth observables — and test whether hierarchy depth predicts combinatorial entropy log H across three independent tiers: (i) all-family exact computation (8 families, N = 10–16, 320 samples), (ii) matched-pair Δ-analysis (Lor2D vs MLR, 46 pairs, N = 30–56), and (iii) coarse-graining identity-stability linkage (92 samples). Matched-pair r = −0.834, stable across three filter stringencies. Simpson's Paradox (naïve r = +0.336) is diagnosed and resolved through N-control.
```

- **Keywords**: `causal posets`, `hierarchy integration`, `combinatorial entropy`, `linear extensions`, `Simpson's Paradox`, `coarse-graining stability`
- **License**: CC BY 4.0
- **Related identifiers**: 
  - `is supplemented by` → `https://doi.org/10.5281/zenodo.18980657`
  - `references` → Prediction B 和 A 的 Zenodo DOI

---

## 上传后操作

1. **三个 DOI 已获取** ✅:
   - Pred B DOI: **10.5281/zenodo.19048146**
   - Pred A DOI: **10.5281/zenodo.19048324**
   - Pred C DOI: **10.5281/zenodo.19048405**

2. **更新相互引用**: Prediction A 引用 B 的 DOI，C 引用 A 和 B 的 DOI

3. **更新代码仓库 README**: 添加三篇预印本的 DOI 链接

4. **注意**: 三篇分开上传（各自独立 DOI），不要合并为一个 upload

---

## Zenodo 版本管理说明

- 上传后获得两个 DOI:
  - **Concept DOI**: 永久指向所有版本的通用 DOI
  - **Version DOI**: 指向当前具体版本
- 后续修改（如合并版）可作为新版本上传，Concept DOI 自动指向最新版
- 在论文中引用时使用 Version DOI（指向具体版本）
