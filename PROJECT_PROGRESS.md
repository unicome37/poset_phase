# Poset Phase 项目工作记录与进度追踪

> **最后更新**: 2026-03-18  
> **Git HEAD**: `1ee2654` (main)  
> **GitHub**: `github.com/unicome37/poset_phase`  
> **版本**: v4.0.0 (Zenodo 代码仓库)  
> **论文版本**: Prediction A v2.6 (**已被 CQG desk rejection, Manuscript ID: CQG-115402**; Carlip 批评 logH 非主流物理量), Prediction C v2 (准因果版)  
> **Zenodo DOI**: `10.5281/zenodo.19079466`

## ⚠️ 文档时效说明（2026-03-30）

本文件主体为历史工作记录，含大量有效过程信息，但**状态口径已部分过时**（如 F3 highN、falsification 执行完成度等）。

> **维护策略（当前）**：本文件作为“历史档案”保留，不再承担滚动状态面板职责。涉及当前结论（尤其 turn-on 边界、FLRW 风险口径、可证伪执行状态）请以 `进展.md` 与 `PROJECT_OVERVIEW_2026Q1.md` 为准。

- **当前仓库基线**：`HEAD = eaea483`
- **最新状态总览**：请优先参考 `进展.md`（2026-03-30 重写版）
- **本文件建议用途**：用于追溯历史实验链路与提交时间线，不作为唯一“当前状态”依据

---

## 一、项目概述

偏序集相变（Poset Phase Transition）项目——基于"结构存在论"理论框架，用数值实验验证三个核心预测（Prediction A/B/C），探讨**存在性筛选**（Existential Selection Dynamics, ESD）如何通过熵-作用量竞争，在偏序集空间中自然涌现出物理时空结构。

### 核心思想

$$\text{Score} = -\beta \cdot \log H(\text{poset}) + \lambda \cdot \text{Action}(\text{poset})/N$$

- $\log H$: 偏序集的线性扩展数（组合熵）
- Action: 因果动力学项（link action 或 BDG action）
- 低 score = 更"存在"的结构

### 核心发现

**Link action $S = N - 2C_0$ 在 λ≈6-8 自然选出 3+1 维时空。** 关键不是 λ 的具体值（因生成器而异），而是无量纲控制参数 **Ξ₄→₅ ≈ 10-11**——从 4D 跃迁到 5D 的链接惩罚/熵比，跨生成器、跨尺寸稳定。这创造了一个**不对称屏障**：抑制 5D 只需很小的耦合，而 4D 仍能在熵上击败 3D。

---

## 二、环境与工具链

| 项 | 值 |
|---|---|
| Python | 3.14.2 (系统 `python` 命令，**不用 .venv**) |
| 项目路径 | `d:\Kiro\理论体系\poset_phase\` |
| 运行方式 | `Set-Location "d:\Kiro\理论体系\poset_phase"; python <script>.py` |
| 编码 | 脚本含中文时必须 `$env:PYTHONIOENCODING='utf-8'` |
| Poset 类 | `poset.closure` (bool ndarray, 传递闭包)，**没有** `poset.adj` |
| 熵估计 | `entropy_sis.py` (SIS 采样, 通常 4096 runs)；`entropy_exact.py` (精确, N≤24) |
| 作用量 | `action.py` → `action_value(log_H, penalty, β, γ)` |
| 几何组件 | `observables_geo.py` → `geometric_components()` |
| 生成器 | `generators.py` → 2D/3D/4D/5D/6D/7D Lorentzian + KR_like |
| LaTeX | MiKTeX (pdflatex)，论文目录 `prediction_a_paper/` |

### 关键注意事项

1. **cwd 必须在 `d:\Kiro\理论体系\poset_phase\`**，否则导入失败
2. **精确计算阈值**: 5D 用 `exact_threshold=24`，N≥32 时 5D 精确计算会 hang
3. **背景终端** 的 cwd 总是 `D:\Kiro`，必须在命令中先 `Set-Location`
4. **KR_like** 是随机偏序（Kleitman-Rothschild），不是物理时空

---

## 三、Zenodo 发布状态

| 资产 | DOI | 类型 | 说明 |
|------|-----|------|------|
| 代码仓库 v4.0.0 | `10.5281/zenodo.18980657` | Software | A/B/C 三预测全部代码 |
| Prediction B | `10.5281/zenodo.19048146` | Preprint | 有界几何相变 (7pp) |
| **Prediction A v1.5** | **`10.5281/zenodo.19068927`** | **Preprint** | **不对称维度屏障 (31pp, 16 fig, LaTeX)** |
| Prediction C v2 | `10.5281/zenodo.19048405` | Preprint | 层级深度预测熵 (43pp, 9 实验) |

---

## 四、三大预测 (Predictions A/B/C) 完整状态

### Prediction A — 维度选择 ★★★★★★ (核心突破)

**核心问题**: 作用量框架是否能从偏序集空间中自然选出 3+1 维时空？

**答案**: 是。Link action 通过不对称的链接密度屏障，使 4D 成为熵-连通性竞争的自然平衡点。

#### 完整实验清单 (21 个脚本)

| 脚本 | 内容 | 关键结果 | Commit |
|------|------|----------|--------|
| `prediction_a_dimension_scan.py` | 基础维度扫描 (2D/3D/4D) | 初步维度选择信号 | 早期 |
| `prediction_a_geometric_ablation.py` | 几何惩罚消融 | A2_consistency 最优 | 早期 |
| `prediction_a_margin_fit.py` | margin(N) OLS 回归 | slope=0.945 (vs2D), R²=0.994 | 82d2729 |
| `prediction_a_margin_plot.py` | margin 可视化 | 出版级图 | 82d2729 |
| `prediction_a_5d_pilot.py` | 5D 天花板测试 | 5D 不饱和! 76.2% 胜率 | 82d2729 |
| `prediction_ac_causal_link.py` | 因果结构×熵跨维诊断 | r = -0.997 (density vs entropy) | 878d8db |
| `prediction_a_bd_dimension.py` | BD 链接作用量维度选择 | **λ=6-8: 4D 全票** | a8e79c7 |
| `prediction_a_bd_lorentzian_only.py` | BD Lorentzian-only 分析 | 维度级联: 5D→4D→3D→2D | a8e79c7 |
| `prediction_a_bd_extended.py` | 扩展 BD (N=60,68, 混合) | 4D plateau N=20-68 全票 | e4d237d |
| `prediction_a_bd_bridge.py` | BD 4D 桥接扫描 | 扫描的 `α_BD` 范围内未发现 Lor2D/Lor3D/Lor4D/Lor5D 的排序翻转；BD-inspired 修正未改写维度次序 |
| `prediction_a_bdg_full_comparison.py` | 文献 BDG 系数 vs link-proxy | link_d2→4D; BDG_d4→5D | 6514e95 |
| `prediction_a_bdg_component_figure.py` | BDG 组件诊断图 | C₁ 项破坏 4D 选择 | 58501cb |
| `prediction_a_generator_robustness.py` | 三种生成器鲁棒性 | 机制跨生成器普适 | b34b8be |
| `prediction_a_seed_sensitivity.py` | 种子敏感性分析 | 稳健 | 早期 |
| `prediction_a_seed_sensitivity_plot.py` | 种子敏感性可视化 | — | 早期 |
| `prediction_a_winner_phase_plot.py` | 赢家相图 | 出版级图 | 早期 |
| `prediction_a_xi_parameter.py` | **Ξ 无量纲控制参数** | **Ξ₄→₅≈10, CV=13.9%** | 41d866b |
| `prediction_a_xi_figure.py` | Ξ 出版级三图 | strip/convergence/asymmetry | 41d866b |
| `prediction_a_large_n_scaling.py` | **大 N 实验 (N=80-112)** | **Ξ₄→₅ 稳定, CV=17.7%** | 36a8f45 |
| `prediction_a_large_n_figures.py` | 大 N 三图 | heatmap/stability/crossover | 36a8f45 |
| `prediction_a_unification.py` | 统一测试: link ↔ geometric | 机制层面统一 | 36a8f45 |
| `prediction_a_xi_derivation.py` | **Ξ 解析推导** | **闭合公式 Ξ₄→₅=11.8** | 8c36896 |
| `prediction_a_d6_extrapolation.py` | **d≥6 盲预测+验证** | **Ξ₅→₆=64/46, Ξ₆→₇=51/38** | afd6bfa |
| `prediction_a_xi_first_principles.py` | **第一性原理 p_d, κ_d 推导** | **从几何恢复 a_d, α_d** | 9c5481e |
| `prediction_a_xi_interval_bridge.py` | Interval hierarchy 桥接 | 1-ℓ_d ≈ interval 汇总, R²=0.9955 | 9c5481e |
| `prediction_a_xi_volume_moments.py` | **体积矩占据闭合** | **P_occ: Ξ₄→₅=11.33, 误差 0.2%** | 9c5481e |
| `prediction_a_xi_gamma_mgf.py` | Gamma-MGF ℓ_d 闭合 | ℓ_d 误差 <1.2%, 但 Ξ 仍 6.9% | 9c5481e |
| `prediction_a_xi_occupancy_closure.py` | Jensen gap + α 扫描验证 | **α*=1.000, 无需自由参数** | 9c5481e |
| `prediction_a_xi_closure_comparison.py` | 统一闭合总表 | 四种闭合 Ξ 误差对比 | 4ca1397 |
| `prediction_a_xi_B0_B1_derivation.py` | **B₀/B₁ 解析推导** | **B₁+2B₀=2, ε=π/48, 零参数** | 4ca1397 |
| `prediction_a_xi_analytic_entropy.py` | 解析熵预测 | Interval hierarchy 预测验证 | 4ca1397 |

#### 发现时间线

**Phase 1: BD 链接作用量发现 (a8e79c7 → e4d237d)**
- Link action $S = N - 2C_0$ 在 λ=6-8 选出 4D，N=20-68 全票稳定
- 级联: 5D(低λ) → 4D(中λ) → 3D(高λ) → 2D(极高λ)

**Phase 2: BDG 对比验证 (6514e95 → 58501cb)**
- 标准 BDG d=4 ($N-C_0+9C_1-16C_2+8C_3$) **全选 5D** — 高阶区间修正破坏维度选择
- 物理诊断: $+9C_1$ 项巨幅放大低维 order-1 区间, 5D 高阶区间极少 → BDG_d4 惩罚极小
- 确认维度选择机制在**链接密度层面**，非曲率层面

**Phase 3: 生成器鲁棒性 (b34b8be)**
- 三种生成器 (cube/indep-seed/diamond) 均保持 4D 选择
- 窗口位置因几何而异 (cube λ≈6-8, diamond λ≈7-10)，机制不变
- 因果链: 不同几何 → 不同因果稀疏度 → 窗口位置漂移 → 但链接密度竞争机制普适

**Phase 4: Ξ 参数发现 (41d866b → 5655002)**
- 定义 $\Xi_{d \to d+1} = |\Delta(S_\text{link}/N)| / |\Delta(\log H/N)|$
- Ξ₄→₅ ≈ 10 跨三种生成器稳定 (CV=13.9%)
- 不对称屏障: $\lambda^*_{4\to5} \approx 0.1 \ll \lambda^*_{3\to4} \approx 0.3$
- 物理含义: 从 4D→5D 的每单位熵代价比 3D→4D 高 ~3×

**Phase 5: 大 N 闭环 (36a8f45)**
- N=80,96,112: λ 窗口右移但 Ξ 不漂 (CV=17.7%)
- 链接密度交叉 N≈96: 4D 密度首次低于 3D — 解释窗口右移
- 间隙比发散: |Δ(4→5)|/|Δ(3→4)| 随 N 增长 (1.4→15.4)
- 统一测试: link density ↔ geometric consistency 机制层面统一

**Phase 6: 解析推导 (8c36896)**
- 标度律: $C_0/N = a_d N^{\alpha_d}$, $\log H/N = b_d \log N + c_d$
- 闭合公式预测 Ξ₄→₅ = 11.8 (实测 11.3，误差 <17%)
- 根因分解:
  - **熵饱和**: Δb₄→₅ = 0.042，仅为 Δb₃→₄ = 0.110 的 38%
  - **链接密度间隙持续**: |ΔC₀/N(4→5)| = 1.126，是 |ΔC₀/N(3→4)| = 0.154 的 7.3×
  - **序分数衰减**: p₂=0.50→p₃=0.29→p₄=0.17→p₅=0.11 (每维 ×0.6)
  - **链接饱和**: d=5 时 94.7% 因果关系是直接链接

**Phase 7: d≥6 盲预测验证 (afd6bfa)**
- 新增 6D/7D 生成器到 `generators.py`
- 盲预测 → 数值验证:

| 转换 | 盲预测 | 实测 | 误差 |
|------|--------|------|------|
| 5→6 | 64.4 | 45.8 | 29% |
| 6→7 | 51.3 | 37.8 | 26% |

- 三项定性预测全部确认: 巨幅放大、非单调平台、量级正确
- **结论**: 4→5 屏障(Ξ≈11)是 d≥3 以上**最弱**的! d≥5 屏障 3-4× 更强

**Phase 8: LaTeX 投稿版 + Zenodo 更新 (4c8e06d → f70ff5b)**
- 完整 LaTeX 论文: `prediction_a_paper/prediction_a.tex` (31pp)
- 16 张出版级图片
- 编译零错误
- Zenodo v1.5 发布: DOI `10.5281/zenodo.19068927`
- 所有文件中的 DOI 引用已同步更新 (10 个文件)

**Phase 9: 第一性原理推导完成 (9c5481e → 62bd420)**
- **闵可夫斯基几何输入**: 从 cube-sprinkle 差分分布推出 $p_d$、$\kappa_d$、$\mathbb{E}[V_A]$，不依赖任何经验拟合
- **四层逐步逼近**:
  1. 最简 mean-field: $h \approx (1-p_d)(\ln N-1)$，Ξ₄→₅ 误差 26%
  2. 全局重整化: $h \approx A + B(1-p_d)(\ln N-1)$，误差 14.2%
  3. 几何混合: $h \approx B_0 h_0 + B_1(1-\ell_d)$，误差 6.7%
  4. **占据闭合**: $h \approx 0.934 h_0 + 0.132 P_{\mathrm{occ}}$，**误差 0.2%**
- **B₁+2B₀=2 结构约束**: 将两参数归约为单参数 $\varepsilon=1-B_0$，R² 损失 $<10^{-8}$
- **ε=κ₄/2=π/48 解析识别**: 与 OLS 拟合仅差 0.42%，bootstrap CI 包含
- **零自由参数闭合**:
  $$\frac{\log H}{N} = \left(1-\frac{\pi}{48}\right)(1-p_d)(\ln N-1) + \frac{\pi}{24}\left(1-e^{-(N-2)\mathbb{E}[V_A]}\right)$$
  给出 Ξ₄→₅=11.31，观测 11.35，误差 **0.4%**
- **论文整合**: Section 5.7 已改写为第一性原理推导，Discussion、Summary 同步更新
- **辅助验证**: α 扫描确认 α*=1.000，Gamma-MGF 闭合 ℓ_d 误差 <1.2%，interval bridge R²=0.9955
- **A 线收口说明**: occupancy closure + volume moments + interval bridge + B₀/B₁ 已形成同一闭环链路，当前以写作整合与跨文档引用为主，不再依赖新增决定性实验

**Phase 10: AI 多轮审阅 + CQG 投稿 (b6379b2 → 1ee2654)**

四轮 AI 审阅（DeepSeek / Gemini / 通义千问 / GPT×4），论文从 v1.5 升级到 v2.6：

| 轮次 | 审阅者 | 修改数 | 版本 | Commit |
|------|--------|--------|------|--------|
| 1 | DeepSeek | 4 项 | — | b6379b2 |
| 2 | Gemini | 3 项 | — | e0db163 |
| 3 | 通义千问 | 3 项 | — | 83a1530 |
| 4 | GPT (v2.0) | 15 项结构修订 | v2.0 | 62a0d3d |
| 5 | GPT (v2.1) | 8 项精修 | v2.1 | 27c490a |
| 6 | GPT (v2.2) | 5 项措辞校准 | v2.2 | fde13f9 |
| 7 | CQG 合规 | Abstract 388→197 词, keywords, refs 9→16, Ack/Funding/Conflict/Data | v2.3 | 802957c |
| 8 | De-AI pass | 去除公式化 AI 痕迹, 34→33 页 | v2.4 | 012d7a1 |
| 9 | Zenodo DOI | Data Availability 写入 DOI | v2.5 | 819c6fe |
| 10 | GPT 最终校准 | 3 处过强措辞软化 | v2.6 | 3915f4b |

- **Cover Letter**: `prediction_a_paper/cover_letter.tex/.pdf` (1 页, commit 4895fe6)
- **投稿包**: `prediction_a_paper/CQG_submission/`
  - `Zhang_CQG_2026_manuscript.pdf` — 完整论文 (33 页, v2.6)
  - `Zhang_CQG_2026_cover_letter.pdf` — 投稿信
  - `source_files.zip` — LaTeX 源文件 + 15 张图
  - `ScholarOne_metadata.txt` — ScholarOne 表单逐步复制粘贴指南
  - `submission_confirmation.txt` — 投稿确认记录
- **Zenodo**: DOI `10.5281/zenodo.19079466` (v2.4 archive, 5.13 MB)

**❌ 已被 CQG desk rejection**（Carlip 批评 logH/线性延拓数不被认为是主流物理量，非物理熵的合法代理）
- **Manuscript ID**: CQG-115402
- **Journal**: Classical and Quantum Gravity
- **Date Submitted**: 17-Mar-2026
- **Title**: Dimensional Selection in Finite Causal-Set Link Dynamics: Evidence for an Asymmetric 4→5 Barrier
- **Authors**: Zhang, Gang
- **查稿网址**: https://mc04.manuscriptcentral.com/cqg-iop

#### 核心数据表

**Ξ 跨生成器 (Table 8 of paper)**:

| 转换 | 原始立方体 | 独立种子 | 因果钻石 | 总体 |
|------|-----------|---------|---------|------|
| Ξ₂→₃ | 3.2 | 1.9 | 1.2 | 1.9 |
| Ξ₃→₄ | 2.1 | 2.7 | 5.0 | 3.3 |
| **Ξ₄→₅** | **10.8** | **10.2** | **9.4** | **10.0** |

**Ξ₄→₅ 跨尺寸 (N=20-112, cube sprinkle)**:

| N | 20 | 36 | 52 | 68 | 80 | 96 | 112 |
|---|----|----|----|----|----|----|-----|
| Ξ₄→₅ | 4.8 | 9.8 | 10.6 | 10.8 | 13.3 | 11.1 | 11.8 |

Overall CV = 17.7%, Small-N median = 10.21, Large-N median = 11.83

**d≥6 盲预测 vs 实测**:

| 转换 | 盲预测 | 实测 | 误差 |
|------|--------|------|------|
| 5→6 | 64.4 | 45.8 | 29% |
| 6→7 | 51.3 | 37.8 | 26% |

**Link action vs BDG d=4 赢家对比 (N=68)**:

| λ | Link ($N-2C_0$) | BDG d=4 |
|---|-----------------|---------|
| 0 | 5D | 5D |
| 5 | 4D | 5D |
| 6 | **4D** | 5D |
| 8 | **4D** | 5D |
| 10 | 3D | 5D |
| 20 | 2D | 5D |

#### 论文状态

- **Markdown**: `PredictionA_LinkAction_Paper_v1.md` (806 行, v1.5)
- **LaTeX**: `prediction_a_paper/prediction_a.tex` (~1673 行, Section 5.7 已改写为零自由参数第一性原理推导)
- **PDF**: `prediction_a_paper/prediction_a.pdf`
- **Zenodo 上传副本**: `prediction_a_paper/PredictionA_LinkAction_v1.5.pdf`
- **Zenodo DOI**: `10.5281/zenodo.19068927`
- **结构**: 8 节 + 3 附录, 22 表, 16 图, 9 参考文献
- **标题**: "An Asymmetric Dimensional Barrier in Causal-Set Link Dynamics Selects 3+1 Dimensions"
- **编译**: pdflatex 零错误、零未定义引用 (commit 62bd420 验证)

#### 出版级图表清单 (16 张, 在 `prediction_a_paper/figures/`)

| 编号 | 文件 | 内容 | 来源目录 |
|------|------|------|----------|
| Fig 1 | `fig01_winner_heatmap_comparison.png` | Link vs BDG 赢家热力图 | prediction_a_bdg_component_figure |
| Fig 2 | `fig02_bdg_component_breakdown.png` | BDG 组件分解图 | prediction_a_bdg_component_figure |
| Fig 3 | `fig03_margin_of_victory.png` | 4D 胜利边际 | prediction_a_margin_summary |
| Fig 4 | `fig04_cascade_phase_diagram.png` | 维度级联相图 | prediction_a_phase_summary |
| Fig 6 | `fig06_xi_core_figure.png` | **Ξ 核心strip图 (3转换×3生成器)** | prediction_a_xi_parameter |
| Fig 7 | `fig07_xi_45_convergence.png` | Ξ₄→₅ 收敛图 | prediction_a_xi_parameter |
| Fig 8 | `fig08_xi_barrier_asymmetry.png` | 屏障不对称性 | prediction_a_xi_parameter |
| Fig 9 | `fig09_large_n_winner_heatmap.png` | 大 N 赢家热力图 | prediction_a_large_n_scaling |
| Fig 10 | `fig10_xi_stability_large_n.png` | Ξ 大 N 稳定性 | prediction_a_large_n_scaling |
| Fig 11 | `fig11_link_density_crossover.png` | 链接密度交叉 | prediction_a_large_n_scaling |
| Fig 12 | `fig12_xi_master_derivation.png` | **解析推导主图 (4面板)** | prediction_a_xi_derivation |
| Fig 13 | `fig13_xi_decomposition.png` | Ξ 分解图 | prediction_a_xi_derivation |
| Fig 14 | `fig14_ordering_fraction_geometry.png` | 序分数几何根 | prediction_a_xi_derivation |
| Fig 15 | `fig15_parameter_extrapolation.png` | 参数外推到 d=6,7 | prediction_a_d6_extrapolation |
| Fig 16 | `fig16_xi_staircase.png` | **Ξ 阶梯图: 盲预测 vs 实测** | prediction_a_d6_extrapolation |
| Fig 17 | `fig17_d567_comparison.png` | d=5,6,7 可观测量对比 | prediction_a_d6_extrapolation |

---

### Prediction B — γ 相变阈值

**核心问题**: 是否存在临界 $\gamma_c$，使得 KR > Lorentzian 切换为 Lorentzian > KR？

#### 已完成实验

| 脚本 | 内容 | 关键结果 |
|------|------|----------|
| `prediction_b_gamma_c_scaling.py` | γ_c(N) 有限尺寸标度 | 常数: 0.528±0.096, 95%CI [0.31, 0.74] |
| `prediction_b_gamma_extended_analysis.py` | 扩展 γ 扫描 [0, 5.0] | N=20-44 精确窗: γ_c = 0.98-1.24 |
| `prediction_b_gamma_extended_nearwall.py` | Near-wall N=48,52,56 | **N=56: γ_c=2.08** (扩展范围后恢复) |
| `prediction_b_candidate_expansion.py` | 候选族扩张（layered 近邻族） | `A2_full`: Lor2D 5/42, KR_like 20/42, layered 17/42；`A2_replace_dim_with_consistency`: Lor2D 0/42 |
| `expanded_family_robustness.py` | 25 家族广义 DAG / interval-order 扩张 | `Lor4D` 在 `N=16,28,48,64,96` 的 LSD-Well 与 Mahalanobis 均为 #1，新增 8 对手未击穿 |
| `prediction_b_seed_reproducibility.py` | 跨 seed 鲁棒性 | `Lor4D` 在 10/10 seed base、8 个 N 上 LSD-Well 与 Mahalanobis 均保持 #1（N=16 也全中） |
| `prediction_b_cross_validation.py` | train/test cross-validation | LSD-Well `24/24` 全中，Mahalanobis `116/120` 为 #1 |
| `prediction_bac_bridge.py` | B/A/C 桥接分析 | B↔C: `corr(delta_hii, delta_penalty) = -0.581`, `corr(delta_hii, delta_log_H) = -0.680`; A↔C: `Lor4D` 在 HII 上不高于 `Lor2D/Lor3D`，但仍保持 A 线赢家地位 |

#### 状态
- v2.1.0 已发布
- γ_c(N) 随 N 缓慢增长但有界
- 候选族扩张表明：Lor2D 对一般 layered 近邻族没有稳健优势，当前 B 的更强版本应收缩到 `KR suppression + local Lorentzian competition`
- 25 家族扩张表明：Lor4D 在当前广义 DAG / percolation / interval-order 基线上仍保持全局赢家地位，但结论仍限定于当前 library 与评分函数
- seed reproducibility 与 CV 表明：Lor4D 的赢家地位并非单一 seed 或训练集自洽造成，而是跨复现稳定
- B/A/C 桥接表明：B 的几何门槛与 C 的历史沉积之间存在稳定负相关，而 A 的 4D 胜出并不依赖“层级更深”这一单轴解释
- **Zenodo DOI**: `10.5281/zenodo.19048146`

---

### Prediction C — 层级深度预测熵

**核心问题**: 偏序集的层级结构（hierarchy depth）是否能预测其组合熵？

#### 状态: v4.0.0 已发布 (准因果证据塔)

**注**: 下列 C 线条目中，
`nearwall` / `LOCO` / `seed 复现` 属于当前双层筛选主线；`F7 large-N` 仅作为旧 `Prediction C` 的辅助验证，不作为主线结论依据。

9 个实验:
1. Stratified Fisher z 回归 (|r|~0.35-0.54)
2. 单边因果干预 (Cohen's d=1.05)
3. 安慰剂控制跨维 (d=1.40-1.83, p<1e-133)
4. 反向干预 100% 方向一致 (d=2.68)
5. k-family 剂量-反应 (all slopes negative)
6. 边密度普遍性 (无相变)
7. N-scaling 法则 (|slope|∝N^0.70)
8. 解析下界 (验证至 1e-15)
9. 大 N SIS 干预 N=16-36

- **手稿**: `PredictionC_MainPaper_Unified_Draft.md` → 43 页 PDF
- **Zenodo DOI**: `10.5281/zenodo.19048405`

- **2026-03-31 C 线 nearwall 补强**:
  - mlr_survivor_matched_lor2d_nearwall_expanded 与 matched_residual_freedom_check_nearwall_expanded 已补齐；
  - prediction_c_pairwise_validation_nearwall_expanded：n_pairs=35，locality_dominance -> ΔlogH Pearson +0.850，hierarchy_integration -> ΔlogH Pearson -0.848；
  - prediction_c_switch_enhancement_scan_nearwall_expanded：baseline mean zeta_cross=4.307，best mean 约 -0.365；
  - prediction_c_switch_enhancement_scan_nearwall_p5p95：baseline mean zeta_cross=4.942，best mean 约 -0.197；
  - prediction_c_pairwise_validation_nearwall_rescue：n_pairs=50，locality_dominance -> ΔlogH Pearson +0.840，hierarchy_integration -> ΔlogH Pearson -0.839；
  - prediction_c_switch_enhancement_scan_nearwall_rescue：baseline mean zeta_cross=5.101，best mean 约 -0.084；
  - 结论：nearwall 的 expanded / p5p95 / rescue 三档已形成同向闭环，主线鲁棒性补强完成。
  - 结论：C 线在 nearwall 多口径下保持同向，已从 local pilot 推进到鲁棒性补强阶段。
- **2026-03-31 C 线效应量张力诊断**:
  - prediction_c_effect_size_tension.py 已完成，输出 outputs_exploratory/effect_size_tension；
  - matched-pair (50 对) 下 layer_count/mean_layer_gap/HII 对 ΔlogH 的相关分别为 -0.8268/-0.8175/-0.8389，bootstrap 95% CI 全为负；
  - leave-k-out 去除 15 对极端样本后仍保持 
≈-0.73~-0.75，信号并非单点驱动；
  - 结论：matched-pair 提供强局部信号，但总体效应量应按中等强度解释（与 Fisher-z 口径一致）。
- **2026-03-31 C 线 LOCO 混杂剥离**:
  - 三组 leave-one-confounder-out 已完成，汇总于 outputs_exploratory/prediction_c_pooled_regression_loo_summary.csv；
  - 在 controls_plus_n_family_fe 下，drop antichain_width 会将 HII 从 +0.086 拉到 -0.495，并显著放大 layer_count/mean_layer_gap 的负相关；
  - drop comparable_fraction 与 drop geo_dim_eff 与 full controls 几乎同结果，提示二者在当前样本里高度冗余。

- **2026-03-31 C 线 F7 large-N fresh（辅助验证，非双层筛选主线）**:
  - prediction_c_f7_large_n.py --fresh --reps 8 --ns 20 36 52 72 100 已完成；
  - within-cell 方向正确率 16/19=84%；
  - per-N 跨族 pooled ρ(Σ_hist,logH) 分别为 -0.805/-0.682/-0.769/-0.743/-0.742，均为显著负相关；
  - 结论：大 N fresh 数据继续支持旧 Prediction C 的方向判断，但仅作为辅助验证，不纳入双层筛选主线。

- **2026-03-31 D 线辅助收口（large-N / weakening / convergence）**:
  - prediction_d_large_n_verification.md 已给出大 N 验证：W1 与 ΔH_int 在控制 F7 后仍显著，控制 F7+N 后仍保持正相关且显著；
  - prediction_d_weakening_diagnosis.md 已给出根因诊断：P_basin 方差随 N 降低，W1 的 CV 基本稳定，weakening 更像 basin variance compression 而非观测量退化；
  - prediction_d_f7_w1_convergence.md 已给出收敛解释：F7 与 W1 在 large-N 下共享 family-level rigidity，within-family residual 相关在大 N 仍显著；
  - 结论：D 的主体实验链已闭合，剩余 TODO 主要是写作升级与是否向“无条件准因果”推进的设计问题。

- **2026-03-31 C 线独立 seed 复现**:
  - prediction_c_comprehensive.py 以 seed=20260331 和 seed=20260401 各完成一次全链重跑；
  - Tier2 核心相关系数保持一致：hii_delta=-0.8363，layer_count_delta=-0.8229，mean_layer_gap_delta=-0.8144；
  - 仅置换 p 值出现千分位波动，方向与量级稳定。

---

### A×C 交叉发现 — 因果稀疏性解释

三个 Prediction 之间的理论链接:

```
Prediction C: 更多层级 → 更少熵 (within-dimension)
     ↓ 扩展到跨维度
Prediction A: 高维 → 更少因果关系 → 更少层级 → 更高熵
     ↓ 解决方案  
Link action: 惩罚因果稀疏性 → 自然选出 3+1 维
```

验证: density-entropy 反相关 r = -0.988 to -0.997 (跨 N=20-52)

---

## 五、关键物理结论总结

### 为什么是 3+1 维？

1. **Link action** ($S = N - 2C_0$) 惩罚因果稀疏性 (链接少 = 惩罚高)
2. **5D** 因果极稀疏 (p₅=0.11, 94.7% 关系是直接链接) → 被轻微耦合即可抑制
3. **4D** 仍有足够熵击败 3D → 在中等耦合下胜出
4. **不对称屏障**: Ξ₄→₅ ≈ 10 >> Ξ₃→₄ ≈ 3 → 4D 上方有硬屏障，下方无
5. **结论**: 4D 不是"特殊"的——4D 上方的屏障是特殊的

### 为什么 BDG d=4 失败？

- $+9C_1$ 项巨幅放大低维 order-1 区间
- 5D 高阶区间极少 → BDG 惩罚极小 → 5D 永远赢
- **双层图像**: 维度选择在**预几何连通性层面** (links)；曲率恢复在**后几何连续极限** (full BDG)

---

## 六、代码结构概览

### 核心基础设施

| 文件 | 功能 |
|------|------|
| `generators.py` | 2D-7D Lorentzian + KR_like 偏序集生成器 |
| `entropy_sis.py` | SIS 熵估计 |
| `entropy_exact.py` | 精确线性扩展计数 |
| `observables.py` | 中性惩罚等基础观测量 |
| `observables_geo.py` | 几何组件 (维度估计等) |
| `action.py` | 作用量计算 |
| `experiment.py` | FAMILIES 字典 + 实验框架 |
| `runtime_utils.py` | 辅助函数 |

### 输出目录 (`outputs_exploratory/prediction_a_*`)

29 个子目录，包含所有 Prediction A 实验的原始数据和图表。

### 论文目录 (`prediction_a_paper/`)

```
prediction_a_paper/
├── prediction_a.tex              # LaTeX 主文件 (v2.6, 33pp)
├── prediction_a.pdf              # 编译 PDF (33pp)
├── cover_letter.tex              # CQG 投稿信
├── cover_letter.pdf              # 投稿信 PDF
├── PredictionA_LinkAction_v1.5.pdf  # Zenodo 上传副本
├── CQG_submission/               # 投稿包 (已提交)
│   ├── Zhang_CQG_2026_manuscript.pdf
│   ├── Zhang_CQG_2026_cover_letter.pdf
│   ├── source_files.zip
│   ├── ScholarOne_metadata.txt
│   └── submission_confirmation.txt
└── figures/                      # 15 张出版级 PNG
```

---

## 七、Git 提交历史 (Prediction A 相关)

```
1ee2654 record CQG submission confirmation: CQG-115402, 17-Mar-2026
4895fe6 add CQG cover letter
3915f4b paper v2.6: final phrase calibration
819c6fe paper v2.5: add Zenodo DOI to Data and Code Availability
012d7a1 paper v2.4: de-AI editorial pass  reduce formulaic patterns
802957c paper v2.3: CQG submission compliance
fde13f9 paper v2.2: final tone calibration for CQG submission
27c490a paper v2.1: GPT final polish  8 refinements
62a0d3d paper v2.0: CQG submission revision  GPT review
83a1530 paper: address Qianwen review
e0db163 paper: address Gemini review
b6379b2 paper: address DeepSeek review
4eb622a docs: comprehensive work records update
62bd420 paper: remove duplicate first-principles section, fix undefined ref
4ca1397 feat: analytical derivation of B0=1-pi/48, B1=pi/24 (zero free parameters)
9c5481e feat: first-principles Xi derivation  occupancy closure (0.2% error)
46a69c5 refactor: revise Prediction A paper per DeepSeek review
f162978 docs: comprehensive PROJECT_PROGRESS.md update
f70ff5b docs: update Prediction A DOI to v1.5 (10.5281/zenodo.19068927)
4c8e06d feat: LaTeX submission version of Prediction A paper (31pp, 16 figures, v1.5)
afd6bfa feat: d>=6 blind extrapolation -- Xi_5->6=64(blind)/46(meas), Xi_6->7=51/38
8c36896 feat: analytical derivation of Xi_4->5 = 11.8 (numerical: 11.3)
36a8f45 feat: large-N closure (N=112) + unification test + paper v1.3
5655002 feat: paper v1.2 — asymmetric barrier framing + Xi as central variable
41d866b feat: Xi dimensionless control parameter analysis + paper v1.1
b34b8be feat: generator robustness test — cube indep-seed + causal diamond
9eac932 docs: Prediction A paper draft v1
58501cb feat: BDG component breakdown figures
6514e95 feat: full BDG coefficient comparison
e4d237d feat: extended BD analysis — d=4 intervals + hybrid + N=68
a8e79c7 feat: Benincasa-Dowker action dimension selection
878d8db feat: cross-dimensional causal diagnostics (Prediction AC link)
82d2729 feat: Prediction A margin OLS + 5D ceiling test
```

---

## 八、已解决的技术问题

| 问题 | 原因 | 解决方案 |
|------|------|----------|
| N=56 "无交叉" | γ 扫描只到 2.0 | 扩展到 5.0, 发现 γ_c=2.08 |
| 5D N=32 精确计算 hang | 5D 稀疏→反链化→指数爆炸 | `exact_threshold=24` |
| `poset.adj` AttributeError | 属性不存在 | 用 `poset.closure` |
| cwd 错误导致 import 失败 | 从 `D:\Kiro` 运行 | 必须 `Set-Location` 到项目目录 |
| GBK 编码错误 | Windows 默认编码 | `$env:PYTHONIOENCODING='utf-8'` |
| BDG d=4 全选 5D | 高阶区间修正破坏选择 | 确认 link action 才是机制 |
| 大 N λ 窗口右移 | 链接密度交叉 N≈96 | Ξ 不漂移，确认机制不变 |
| LaTeX exit code 1 | pdflatex 返回 1 但 PDF 正确 | 多次编译解析引用即可 |

---

## 九、下一步可能方向

### 已完成 ✅

1. ~~BD 作用量天然选出 3+1 维~~ ✅ (a8e79c7)
2. ~~有限尺寸标度 N=20-112~~ ✅ (36a8f45)
3. ~~BDG 文献系数验证~~ ✅ (6514e95)
4. ~~生成器鲁棒性 (3 种)~~ ✅ (b34b8be)
5. ~~Ξ 无量纲控制参数~~ ✅ (41d866b)
6. ~~Ξ 解析推导 (闭合公式)~~ ✅ (8c36896)
7. ~~d≥6 盲预测验证~~ ✅ (afd6bfa)
8. ~~论文 v1.5 LaTeX 投稿版~~ ✅ (4c8e06d)
9. ~~Zenodo v1.5 发布~~ ✅ (f70ff5b)
10. ~~所有 DOI 引用更新~~ ✅ (f70ff5b)
11. ~~Ξ 第一性原理推导 — 占据闭合 (P_occ 0.2% 误差)~~ ✅ (9c5481e)
12. ~~B₀/B₁ 解析推导 — B₁+2B₀=2, ε=π/48 零自由参数~~ ✅ (4ca1397)
13. ~~论文 Section 5.7 改写 — 第一性原理推导整合~~ ✅ (62bd420)
14. ~~四轮 AI 审阅 (DeepSeek/Gemini/通义千问/GPT×4)~~ ✅ (b6379b2 → 3915f4b)
15. ~~CQG 合规 (Abstract/Keywords/Refs/Ack/Funding/Conflict/Data)~~ ✅ (802957c)
16. ~~De-AI 编辑 + 措辞终校~~ ✅ (012d7a1 → 3915f4b)
17. ~~Cover Letter~~ ✅ (4895fe6)
18. ~~Zenodo DOI 写入 Data Availability~~ ✅ (819c6fe)
19. ~~**CQG 投稿 — Manuscript ID: CQG-115402**~~ ✅ (1ee2654, 17-Mar-2026)

### 待做 (高优先级)

1. **Ξ 第一性原理推导（已完成 ✅）**: 详见 `FIRST_PRINCIPLES_XI_STATUS.md`。
   - 推导链完整: Minkowski 几何 → $p_d, \kappa_d, \mathbb{E}[V_A]$ → $C_0/N$ → $P_{\mathrm{occ}}$ → 熵闭合 ($\varepsilon=\pi/48$) → Ξ₄→₅ 0.4% 误差
   - 零自由参数闭合:
     $$\frac{\log H}{N} = \left(1-\frac{\pi}{48}\right)(1-p_d)(\ln N-1) + \frac{\pi}{24}\left(1-e^{-(N-2)\mathbb{E}[V_A]}\right)$$
   - 论文已同步: Section 5.7 改写, R5 + Discussion 更新, pdflatex 编译零错误 (62bd420)
   - **剩余开放问题**:
     - [ ] finite-cube boundary correction 对 d=5 的 $a_d$ 偏高
     - [ ] $\varepsilon = \kappa_4/2$ 的严格物理推导 (当前为 numerical identification)
     - [ ] d=6 上验证零参数闭合的 Ξ₅→₆ 精度
2. ~~**论文投稿**: 选择期刊 (CQG)，调整格式提交~~ ✅ **已投稿 CQG-115402 (17-Mar-2026)**
   - ✅ 四轮 AI 审阅 (DeepSeek / Gemini / 通义千问 / GPT×4)，论文 v1.5 → v2.6
   - ✅ CQG 合规: Abstract 197 词, 8 keywords, 16 refs, Ack/Funding/Conflict/Data 声明
   - ✅ De-AI 编辑 + 最终措辞校准 (33 页)
   - ✅ Cover Letter 编写 (1 页)
   - ✅ Zenodo DOI: `10.5281/zenodo.19079466`
   - ✅ ScholarOne 投稿包: `prediction_a_paper/CQG_submission/`
   - ✅ **已提交**: Manuscript ID **CQG-115402**, 17-Mar-2026
   - 下一动作: 等待编辑初审 (1-2 周) → 同行评审 (6-12 周)
3. **代码仓库 v5.0.0**: 包含 d≥6 生成器 + LaTeX 论文 → Zenodo 更新

### 待做 (中优先级)

4. **弯曲时空撒点**: de Sitter / Schwarzschild 几何下测试 Ξ 稳定性
5. **随机因果三角化**: 非撒点类型生成器
6. **Bootstrap 置信区间**: 加强统计推断力度
7. **Prediction B 论文更新**: 加入 N=56 γ_c=2.08
8. **推论D（动态过程）确认性窗口冻结**: 从 v8 给出的候选冻结窗口里选 2-3 个，对应 strata 做更高 `n_perm` 与跨 source 复现
   - ✅ 已完成确认性冻结（已更新）：冻结窗口收敛为 `N=30, keep_ratio=0.6, gamma=0.2`
   - 说明：`gamma=0.8` 在后续确认性复现（`sis_runs=128`）下出现不稳定（`no_switch` 乃至 `switch` 可能坍缩），因此不宜冻结为确认性主窗口
   - 主证据：在 `sis_runs=128` 的四次独立复现（`seed_offset=172000/174000/175000/176000`）上，对 `full/switch/no_switch` 三口径以 `n_perm=200000` 做第三层定点加密，均保持 `rho>0` 且强显著（见 `outputs_confirmatory/prediction_d_dynamic_v10_freeze_strata_rep3_g02/` 等）

9. **推论D 论文/章节撰写**: 三层证据 + 冻结窗口 + 负结果边界 → 论文级叙事整理
   - 草稿：`PREDICTION_D_WRITEUP_DRAFT.md`
10. **更多 N 的窗口扩展**: 在确认性复现族上寻找并定点加密 N=40/52 的可重复窗口
   - ✅ N=40 次冻结候选：`N=40, keep_ratio=0.6, gamma=0.8` 在 rep3–rep6 上定点加密后 `full/switch/no_switch` 均显著（见 `outputs_confirmatory/prediction_d_dynamic_v10_window40_g08_summary.csv`）
   - ⚠️ N=52 暂不稳定：`N=52, keep_ratio=0.6, gamma=0.2` 在 rep3/4/6 很强，但在 rep5 出现 `no_switch` 失效与 `switch` 边界显著（见 `outputs_confirmatory/prediction_d_dynamic_v10_window52_g02_summary.csv`），当前不冻结
   - ✅ N=52 去噪补充：专用去噪 run（更大 `samples_per_family`/`repeats`）可让 `N=52, keep=0.6, gamma=0.2` 的 `full/switch/no_switch` 三口径恢复强显著（见 `outputs_confirmatory/prediction_d_dynamic_v10_window52_denoise_confirm_v1_g02/`）
   - ✅ N=52 去噪复现2：`outputs_confirmatory/prediction_d_dynamic_v10_window52_denoise_confirm_v2_g02/`（并汇总于 `outputs_confirmatory/prediction_d_dynamic_v10_window52_denoise_g02_summary.csv`）
   - ⚠️ keep_ratio=0.7 次窗口（标准规程）不更稳：`N=52, keep_ratio=0.7, gamma=0.2` 两次独立复现结果不一致（见 `outputs_confirmatory/prediction_d_dynamic_v9_window52_keep07_standard_confirm_v1/` 与 `..._v2/`）
   - ❌ keep_ratio=0.65 次窗口（标准规程）失败：`N=52, keep_ratio=0.65, gamma=0.2` 两次独立复现 v9 显示 `full` 接近 0 且不显著，`no_switch` 甚至出现符号翻转（见 `outputs_confirmatory/prediction_d_dynamic_v9_window52_keep065_standard_confirm_v1/` 与 `..._v2/`）
   - ❌ keep_ratio=0.75 次窗口（标准规程）失败：`N=52, keep_ratio=0.75, gamma=0.2` 两次独立复现 v9 未出现 `full/switch/no_switch` 三口径稳定显著（见 `outputs_confirmatory/prediction_d_dynamic_v9_window52_keep075_standard_confirm_v1/` 与 `..._v2/`）
   - ✅ 跨推论最小交叉检验（C → D 控制变量）：在冻结窗口上把 `layer_count/mean_layer_gap` 加入第三层作为控制，`I_cg → improve_rank` 仍强显著，说明第三层不只是“层深代理”。产物：`outputs_confirmatory/prediction_d_dynamic_v11_layercontrol_rep3_n30_k060_g02_np20000/`（方法见 `结构存在论_推论D验证方案.md`）   - ✅ **v12 tie-aware Spearman 修复**：v10 使用 `np.argsort(np.argsort())` 作为 rank 近似，对 `improve_rank`（整数型、高 N 大量 ties）引入噪声并压低 ρ。v12 替换为平均秩实现（`_tieaware_rank()`），**彻底解决 N=52 rep5 的"不稳"**：`no_switch` 从 ρ=0.096/p=0.48（不显著）恢复为 ρ=0.553/p=5e-6（强显著）。N=30 冻结窗口 ρ 全面上升（full 0.70→0.78, no_switch 0.50→0.72, switch 0.74→0.81）。产物：`outputs_confirmatory/prediction_d_dynamic_v12_tieaware_rep{3,4,5}_g02/`
   - ✅ **N=52 升格为次级确认性窗口**：在 v12 下，rep4/rep5 的 N=52 三口径均 ρ>0.55, p=5e-6，跨复现一致。之前的"不稳定"结论撤销（系 tie 处理伪影）
   - ✅ **剂量-响应分析（dose-response）**：`prediction_d_dose_response.py` 对 rep3–rep6 共 2160 行做 5 项分析：(1) 块内 OLS 斜率全正，r ∈ [0.52, 0.95]；(2) I_cg 三等分 → improve_rank 单调递增 **9/9**；(3) 斜率未随 N 递增（N=30 ≈ 3.3, N=40/52 ≈ 1.5），但 Spearman ρ 跨 N 稳定（0.65–0.86）；(4) 池化 p 全部 < 1e-30；(5) 成分分解：retain_identity / penalty_cg（ρ≈±0.80）为主驱动，结构形状特征几乎无信号。产物：`outputs_exploratory/prediction_d_dose_response/`
   - ✅ **CG 扰动准介入实验**：`prediction_d_perturbation_intervention.py` — 随机移除 5%/10%/20% Hasse 覆盖关系，重建传递闭包后跑完整 CG 管线。核心发现：(1) 剂量强度确认——mean|ΔI_cg| 随扰动比例单调递增（3/3 N）；(2) 样本级 ΔI_cg ↔ Δpenalty 近完美锁定（ρ ≈ -0.98，族内同样成立）；(3) 刚性检验混合——N=40 低扰动下 CG 稳定样本更鲁棒（ρ=-0.35, p=0.015），但 N=52/20% 出现符号翻转。产物：`outputs_exploratory/prediction_d_perturbation/`
   - ✅ **准介入证据边界（机制独立 target）**：用 `Y=Δscore_local`（不含 CG 项）检验 `Δpenalty_cg` 的预测力：pooled 下存在剂量依赖信号（p20 ρ≈−0.306, p≈6e-4），但在 (N,family) 分层置换后失显（18 strata×8 样本均不显著）→ 更像“结构共变”而非层内准因果。产物：`outputs_exploratory/prediction_d_perturbation/perturbation_independent_target_report.md`
   - ✅ **分层置换准因果检验**：`prediction_d_perturbation_stratified_test.py` — 用独立目标 Y=Δscore_local（不含 CG 项）测试扰动诱导的 Δpenalty_cg 是否预测竞争力变化。核心发现：(1) 池化信号确认——p20: ρ=-0.31, p=0.0002***，信号随扰动剂量递增（dose-response）；(2) **分层后失显**——按 (N, family) 分层置换，ρ̄ 全部不显著（p>0.11）；(3) N=30 最敏感——p20 下 ρ=-0.53, p=0.0002***。结论：Δpenalty→Δscore_local 属"结构共变"而非族内准因果效应，D 的证据水平从"准因果"下修为"结构共变"。产物：`outputs_exploratory/prediction_d_perturbation/stratified_intervention_results.csv`
   - ✅ **高功效分层重测 (n=32/层)**：`prediction_d_stratified_highpower.py` — 将每层样本量从 8 提升至 32（576 对观测），以排除 §6.5 零结果是功效不足假象的可能。决定性结论：(1) 池化效应加强——三级扰动均达 p<0.01（p20: ρ=-0.25, p=0.00001***）；(2) **分层 p 值进一步远离显著**——p05: 0.113→0.170, p10: 0.613→0.793, p20: 0.497→0.618；(3) 层内 ρ̄ 衰减至 ≈0（p10 甚至变为 +0.011 正值）；(4) %negative ≈ 50%，完全符合零假设。**最终裁定：结构共变确认，非功效不足假象。** 产物：`outputs_exploratory/prediction_d_perturbation_n32/perturbation_sample_cg_n32.csv`

   - ✅ **连续Y+富残差化层内准因果检验**：`prediction_d_continuous_residualized.py` — 针对"结构共变"裁定的设计改进：(1) 5 个连续 Y 目标（Δscore_local, Δlog_H, Δpenalty_local, Δgeo_total, Δpenalty_neutral）拆分信号通道；(2) 11 维混杂变量层内 OLS 残差化（结构签名 + 基线适应度）；(3) 边界样本子分析。**突破性发现**：残差化后层内信号恢复！p05 下 Δscore_local/Δlog_H/Δpenalty_neutral 均层内显著（p=0.019/0.019/0.001）；**p20 下全部 5 个 Y 目标均达显著（含 Δpenalty_neutral all/boundary：ρ̄=-0.187,p=0.00002；ρ̄=-0.134,p=0.027）**。剂量递增确认（ρ̄ 和显著性均随扰动递增）。边界样本 p20 同样全面显著。**直接逆转 §6.5–6.6 的"层内零"结论——信号被层内混杂掩盖，非不存在。D 证据水平上修为"条件准因果"。** 产物：`outputs_exploratory/prediction_d_continuous_residualized/`
   - ✅ **最小复现锚点（n=32, mechanism-independent 连续Y）**：`prediction_d_residualized_within_stratum_test.py` — 不依赖 rich 11 维特征重建，仅用现成 `perturbation_sample_cg_n32.csv` 做层内去均值 + 基线残差化（`pen0/score0/icg0`）。`n_perm=10000` 下 stratified 结果：p05 `ρ̄≈-0.086,p≈0.0379`；p10 `ρ̄≈-0.100,p≈0.0173`；p20 `ρ̄≈-0.192,p≈1e-4`。说明“层内信号回潮”不只存在于 rich 模型，也可由轻量残差化独立复现。产物：`outputs_exploratory/prediction_d_perturbation_residualized_np10000/`

   - ⏭️ **TODO（进一步巩固"条件准因果"）**：前 3 项（独立 seed、partial ρ 上限、leave-one-confounder-out）已通过 D 线辅助收口报告部分覆盖；当前仍开放的是 (4) 是否可将条件准因果升级为无条件准因果（需更多 family 和更大 N 的设计）。
### 探索性

8. **连续极限**: N→∞ 下 Ξ 是否收敛到精确常数？
9. **Myrheim-Meyer 维度估计量**: 替代维度选择机制？
10. **BD + consistency 混合**: 更窄 λ 窗口锁定 4D？
11. **统一 B+A 论文** (`manuscript_unified/manuscript_cqg.tex`)

---

## 十、关键函数签名速查

```python
# 生成偏序集
from generators import generate_lorentzian_like_2d   # (n, seed) → Poset
from generators import generate_lorentzian_like_5d   # (n, seed) → Poset
from generators import generate_lorentzian_like_6d   # (n, seed) → Poset
from generators import generate_lorentzian_like_7d   # (n, seed) → Poset

# 熵估计  
from entropy_sis import estimate_entropy
# estimate_entropy(poset, sis_runs, seed, exact_threshold) → (log_h, method)

# 观测量
from observables import neutral_penalty              # (poset) → float
from observables_geo import geometric_components     # (poset) → dict
from observables_geo import cover_density            # (poset) → float

# 区间计数 (在各脚本中定义)
# hasse_links(poset) → int (C₀)
# count_intervals(poset, k) → int (C_k)

# Actions:
# Link action: S = N - 2*C₀
# BDG d=4: S = N - C₀ + 9*C₁ - 16*C₂ + 8*C₃
# Score: score = -β*log_H + λ*S/N
```

---

## 十一、实验配置模板

### 标准维度选择实验

```python
N_VALUES = [20, 28, 36, 44, 52, 60, 68]
FAMILIES = ["lorentzian_like_2d", "lorentzian_like_3d",
            "lorentzian_like_4d", "lorentzian_like_5d"]
SAMPLES_PER_FAMILY = 4
SIS_RUNS = 4096
FAMILY_EXACT_THRESHOLDS = {
    "lorentzian_like_2d": 104,
    "lorentzian_like_3d": 24,
    "lorentzian_like_4d": 24,
    "lorentzian_like_5d": 24,
}
LAMBDA_VALUES = [0, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20]
BETA = 1.0
```

### 大 N 实验配置

```python
N_VALUES = [20, 36, 52, 68, 80, 96, 112]
LAMBDA_VALUES = [5, 6, 7, 8, 10]
SEED_BASE = 980000
```

### d≥6 外推实验配置

```python
FAMILIES_HIGH_D = ["lorentzian_like_5d", "lorentzian_like_6d", "lorentzian_like_7d"]
N_VALUES = [20, 36, 52, 68]
SIS_RUNS = 4096
```

---

## 十二、论文核心方程速查

**Link action**:
$$S_\text{link} = N - 2C_0$$

**Score function**:
$$\text{Score}(G;\lambda) = -\beta \cdot H(G) + \lambda \cdot S(G)/N$$

**Ξ 定义**:
$$\Xi_{d \to d+1} = \frac{|\overline{S}_{d+1} - \overline{S}_d|}{|\overline{h}_{d+1} - \overline{h}_d|}$$

**Ξ 闭合公式**:
$$\Xi_{d \to d+1}(N) = \frac{2|a_d \cdot N^{\alpha_d} - a_{d+1} \cdot N^{\alpha_{d+1}}|}{|(b_{d+1}-b_d) \cdot \log N + (c_{d+1}-c_d)|}$$

**不对称屏障**:
$$\lambda^*_{4 \to 5} \approx 0.1 \;\ll\; \lambda^*_{3 \to 4} \approx 0.3$$

**BDG d=4 (对比用)**:
$$S^{(4)} = N - C_0 + 9C_1 - 16C_2 + 8C_3$$

**标度律参数 (Table 9 of paper)**:

| $d$ | $a_d$ | $\alpha_d$ | $b_d$ | $c_d$ |
|-----|--------|-----------|--------|--------|
| 2 | 0.573 | 0.372 | 0.487 | −0.399 |
| 3 | 0.301 | 0.618 | 0.621 | −0.393 |
| 4 | 0.114 | 0.839 | 0.731 | −0.555 |
| 5 | 0.034 | 1.049 | 0.773 | −0.541 |








