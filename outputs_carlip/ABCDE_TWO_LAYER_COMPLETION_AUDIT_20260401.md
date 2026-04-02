# ABCDE 五线在两层筛选新主线下的完成度审计

> **日期**: 2026-04-01
> **基线**: `HEAD = eaea483`, `outputs_carlip/` 全部文件 + `进展.md` 重写版
> **口径**: 以两层筛选（S_BD 准入 + S_MD 身份）为当前投稿主线
> **背景**: Prediction A 旧论文（CQG-115402）已被 desk rejection（Carlip 批评 logH 非物理量）

---

## 总览

| 线 | 新主线角色 | 实验完成度 | 写作完成度 | 收口差距 |
|---|---|---|---|---|
| A（维度选择） | Layer 1 物理根因 + Layer 2 唯一选择 | **95%** | **70%** | 需整合进主稿 §2-§3 |
| B（有界阈值/家族选择） | 两层筛选核心主张 | **95%** | **75%** | 主稿 §5-§6 未完，Letter 未转 LaTeX |
| C（层级-熵/历史沉积） | Basin deepening 机制解释 | **90%** | **50%** | 缺独立证据总结文档 |
| D（粗粒化稳定性） | 跨尺度身份保持 | **85%** | **40%** | WRITEUP_DRAFT 未完成 |
| E（曲率编码/时间箭头） | Layer 1 物理合法性支撑 | **90%** | **80%** | 仅需 Discussion 桥接段 |

---

## A 线：维度选择（新主线下已基本完成）

### 新主线下 A 线的重新定位

旧 A 线问的是"link action 能否选出 4D"，答案是 Ξ₄→₅ ≈ 10 的不对称屏障。
新 A 线问的是"为什么两层筛选的交集恰好是 Lor4D"，答案分两层：

- **Layer 1（S_BD 准入）**: S_BD 单独排名 Lor4D 仅 #14/17（N=128），但它淘汰了全部 KR 型和大部分 layered 族，把 17 族压缩到 ~6 族。
- **Layer 2（S_MD 身份）**: Mahalanobis 零参数，从 N≥10 起稳定 #1/25，margin 从 +0.308 发散到 1.93×10⁸。
- **交集**: {S_BD ∈ window} ∩ {S_MD ≈ 0} = {Lor4D}，在 N=64 和 N=128 均已验证。

### 已完成的实验证据（全部在 outputs_carlip/）

#### Layer 1 相关

| 实验 | 文件 | 关键结论 |
|---|---|---|
| S_BD 单独排名 | `smd_sbd_connection.md` | N=128: Lor4D #14/17，KR_like #1 |
| S_BD vs S_MD 正交性 | `smd_sbd_connection.md` | Pearson r→0（bootstrap 95% CI 跨零），排名几乎正交 |
| 联合选择 | `smd_sbd_connection.md` | N=64/128: 交集 = {Lor4D} |
| Link action 17 族测试 | `link_action_17family_test.md` | S_link 单独不能选出 Lor4D（全 N 均 LOSES） |
| 层级筛选原理 | `hierarchical_screening.md` | d→c→w 自然顺序，N≥96 完美消除全部 16 族 |
| 三重筛选定量化 | `triple_screening_results.txt` | S_triple 排 Lor4D #12/17，仅为准入门 |

#### Layer 2 相关

| 实验 | 文件 | 关键结论 |
|---|---|---|
| Mahalanobis 零参数 | `mahalanobis_lsd.md` | 全 N≥10 #1/25，零自由参数 |
| F2 onset refit | `f2_turnon_margin_refit.md` | 20 seeds×120 reps，三档 onset 均为 N=10 |
| Basin deepening | `basin_deepening_results.txt` | V_eff ∝ N⁻¹·⁵⁷，Fisher ∝ N，gap 单调增长 |
| 参考流形轨迹 | `mu_trajectory.md` | σ² ∝ N⁻¹，det(Σ) ∝ N⁻³·³⁸，μ(∞) 有理论推导 |
| Well center 第一原理 | `cstar_wstar_first_principles.md` | d*=4 精确，c*(∞)=0.374，w∝N⁻¹/⁴ |
| 特征剔除 | `feature_ablation_test.md` | (d_eff, C₁/C₀, w) 最小完备三联体 |
| N=1024 极限 | `large_n_extreme.md` | margin 达 1.93×10⁸，幂律发散 |
| 25 族鲁棒性 | `expanded_family_results.txt` | +8 对抗性族，Lor4D 全 N #1 |
| Counter-factual self-selection | `manuscript_supplement_experiments.md` Exp 1 | 4 个 Lor 中心×6 个 N，24/24 自选 |
| Seed 再现性 | `prediction_b_seed_reproducibility.md` | LSD 80/80 #1，Mahal 79/80 #1 |
| 交叉验证 | `prediction_b_cross_validation.md` | LSD 98%，Mahal 94%，margin 保留率 >95% |
| Bootstrap 置信区间 | `bootstrap_confidence.md` | N≥20 P(#1)≥97%，margin CI 全正 |
| KR_2layer 深度分析 | `kr_2layer_analysis.md` | 结构相似是偶然的，C₁/C₀≡0，gap 单调增长 |
| Lor5D 维度编码 | `lor5d_dimension_encoding.md` | MM 分辨率 ~2.4σ at N=16，临界 N≈20 |
| Fisher 权重假说 | `fisher_weight_hypothesis.md` | Σ⁻¹ 排序 = 经验权重排序 |
| 最小失真验证 | `min_distortion_verification.md` | Landau 展开形式，η=0.74±0.08 |

#### 曲率鲁棒性（F3）

| 实验 | 文件 | 关键结论 |
|---|---|---|
| de Sitter | `falsify_c2_background_response_desitter_*.md` | Pass（lowN + highN） |
| Schwarzschild | `falsify_c2_background_response_schwarz_*.md` | Pass（lowN + highN） |
| FLRW κ=1.0 | `falsify_c2_background_response_flrw_*.md` | 边界敏感（lowN hit，highN fail_ratio=0.3<0.5） |
| FLRW metric-faithful PhaseA | `falsify_c2_background_response_flrw_metricfaithful_phaseA.md` | metric 分支 11/60=0.183，proxy 0/60 |
| FLRW 防御简报 | `FLRW_KAPPA1_DEFENSE_BRIEF_20260330.md` | background-dependent robustness 口径 |

#### 可证伪链路（F1）

| 实验 | 文件 | 关键结论 |
|---|---|---|
| 家族压力 | `falsify_c1_family_pressure.md` | Hard fail = NO，全 N 全 seed Lor4D #1 |
| F1 总结 | `F1_ONEPAGE_SUMMARY_20260328.md` | C1 风险基本解除 |

### A 线结论

**新 A 线在实验层面已基本完成。** 核心证据链：
1. S_BD 单独不能选 Lor4D（#14/17）→ 需要第二层
2. S_MD 零参数从 N≥10 起唯一选出 Lor4D → 身份层有效
3. 两层正交（r→0）且交集 = {Lor4D} → 两层互补
4. Basin deepening 随 N 发散 → 身份锁定越来越强
5. 25 族 + 8 对抗族 + 20 seed + CV 全部通过 → 鲁棒性充分
6. 曲率响应 = background-dependent robustness → 边界已标定

**缺口**: 写作整合。旧 A 线的 Ξ 屏障分析可作为"为什么 S_BD 窗口在 d=4 附近有效"的物理解释，整合进主稿 §2。

---

## B 线：有界阈值 / 家族选择（新主线核心，接近收口）

### 新主线下 B 线的重新定位

旧 B 线问的是"Lor 是否在 F7 下赢 KR"，已被 Carlip 击破。
新 B 线 = 两层筛选主线本身：Lor4D 在 S_BD + S_MD 下被唯一选出。

### 已完成的实验

与 A 线 Layer 2 完全重合（见上表）。B 线的独立贡献在于：
- 从 7 族扩展到 17 族再到 25 族的完整样本空间演进
- LSD-Well → Mahalanobis 的方法论演进
- Prediction B 修订版（`prediction_b_lsd_well_revision.md`）

### B 线结论

**实验完成度 ~95%。** 所有核心实验已完成。

**写作缺口**:
1. `SMD_OPERATOR_LETTER.md` v3 已就绪，需转 LaTeX
2. `MANUSCRIPT_SECTIONS_1_4.md` 已写到 §4.3，§5（historical sedimentation）和 §6（discussion + limitations）未完
3. `manuscript_cqg_revised.tex` 已存在但需与最新口径同步
4. 曲率口径需统一为 background-dependent robustness

---

## C 线：层级-熵 / 历史沉积（实验链已闭合，缺写作收束）

### 新主线下 C 线的重新定位

旧 C 线问的是"层深是否预测 logH"。
新 C 线问的是"为什么 Lor4D 的 basin 会随 N 深化"——答案是层级整合压低组合自由度。

### 已完成的实验（来自推论C验证方案）

| 实验层 | 状态 | 关键结果 |
|---|---|---|
| 配对检验（实验 B） | ✅ 完成 | LDI→ΔlogH: r=+0.759; HII→ΔlogH: r=-0.757; layer_count→ΔlogH: r=-0.783 |
| Nearwall 三档闭环 | ✅ 完成 | expanded/p5p95/rescue 三档方向一致，Pearson ≈ -0.82~-0.85 |
| 全样本受控回归（实验 A） | ✅ 完成 | layer_count partial r=-0.450 (p=0.00025)，控制 width/comp/dim_eff 后仍显著 |
| LOCO 混杂剥离 | ✅ 完成 | antichain_width 是主混杂通道，drop 后 HII 从 +0.086 拉到 -0.495 |
| 筛子增强扫描（实验 C） | ✅ 完成 | baseline zeta_cross≈4~5 降到 ≈-0.1~-0.5，层级项有筛选增强潜力 |
| 效应量张力诊断 | ✅ 完成 | matched-pair 强局部信号，leave-k-out 去极端后仍 ≈-0.73~-0.75 |
| 独立 seed 复现 | ✅ 完成 | 双 seed 重跑，Tier2 相关系数完全一致（hii_delta=-0.8363） |
| F7 large-N fresh | ✅ 完成 | within-cell 方向正确率 84%，per-N pooled ρ 全部 -0.68~-0.81 |

### C 线结论

**实验完成度 ~90%。** 三层实验（配对/回归/增强）全部通过，nearwall 三档闭环，独立 seed 复现稳定。

**缺口**:
1. 缺一份独立的 `PREDICTION_C_EVIDENCE_SUMMARY.md` 收束文档
2. 需要在主稿 §5（historical sedimentation）中整合 C 线结论作为 basin deepening 的机制解释
3. C 线的理论定位需明确写出：不再是"logH 被层深预测"，而是"层级整合是 basin deepening 的结构根因"

---

## D 线：粗粒化稳定性（证据链已闭合，写作最弱）

### 新主线下 D 线的重新定位

D 线问的是"Lor4D 的身份在粗粒化后是否保持"，对应总论稿中"换了尺度还能不能认出来"。

### 已完成的实验（来自 PROJECT_PROGRESS.md）

| 实验 | 状态 | 关键结果 |
|---|---|---|
| 确认性冻结窗口 | ✅ | N=30, keep=0.6, γ=0.2，四次独立复现均强显著 |
| N=40 次冻结 | ✅ | rep3-rep6 定点加密，full/switch/no_switch 均显著 |
| N=52 升格 | ✅ | v12 tie-aware 修复后 rep4/rep5 三口径均 ρ>0.55, p=5e-6 |
| v12 tie-aware Spearman | ✅ | 彻底解决 N=52 rep5 "不稳"（系 tie 处理伪影） |
| 剂量-响应分析 | ✅ | I_cg 三等分→improve_rank 单调递增 9/9，池化 p<1e-30 |
| CG 扰动准介入 | ✅ | ΔI_cg↔Δpenalty 近完美锁定 ρ≈-0.98 |
| 分层置换检验 | ✅ | 池化信号确认，分层后失显→"结构共变" |
| 高功效分层重测 | ✅ | n=32/层，分层 p 进一步远离显著→确认非功效不足 |
| 连续Y+残差化 | ✅ | **突破**：残差化后层内信号恢复，p20 全 5 个 Y 显著 |
| 最小复现锚点 | ✅ | 轻量残差化独立复现，p20 ρ̄≈-0.192, p≈1e-4 |
| 跨推论交叉检验 | ✅ | C→D 控制变量后 I_cg→improve_rank 仍强显著 |

### D 线结论

**实验完成度 ~85%。** 证据水平已从"结构共变"上修为"条件准因果"。

**缺口**:
1. `PREDICTION_D_WRITEUP_DRAFT.md` 存在但未完成——这是最大缺口
2. 需要把三层证据（冻结窗口 + 剂量-响应 + 残差化层内信号恢复）整理成论文级叙事
3. 在主稿中的位置：§4 或 §5 的子节，回答"身份在尺度变换中是否保持"

---

## E 线：曲率编码 / 时间箭头（最独立，基本收口）

### 新主线下 E 线的重新定位

E 线完全不依赖 logH，不依赖家族选择，是整个体系中最坚实的部分。
23 项实验的结论（因果集编码曲率）是纯物理结果。

### 与新主线的联系

- E 线的曲率编码结果为 S_BD 作为 Layer 1 准入层提供了物理合法性支撑
- S_BD 本身就是 BDG action 的 d=4 形式，而 E 线验证了 BDG action 确实编码标量曲率
- 时间箭头假说（结构增长不对称）与 basin deepening（身份随 N 锁定）方向一致

### E 线结论

**实验完成度 ~90%，写作完成度 ~80%。**

**缺口**: 仅需在主稿 Discussion 中加一段 E→S_BD 的桥接，说明 E 线为什么支持 Layer 1 的物理合法性。

---

## 收口优先级排序

### 第一优先：B 线论文收口（距投稿最近）

1. 把 `SMD_OPERATOR_LETTER.md` v3 转成 LaTeX
2. 补完 `MANUSCRIPT_SECTIONS_1_4.md` 的 §5（historical sedimentation）和 §6（discussion + limitations）
3. 统一曲率口径为 background-dependent robustness
4. 整合 `FLRW_KAPPA1_DEFENSE_BRIEF` 进主稿

### 第二优先：C 线证据总结

5. 写 `PREDICTION_C_EVIDENCE_SUMMARY.md`，收束配对+回归+LOCO+nearwall+seed 复现
6. 整合进主稿 §5 作为 basin deepening 机制解释

### 第三优先：A 线写作整合

7. 把旧 A 线 Ξ 屏障分析整合进主稿 §2（Layer 1 物理根因）
8. 不需要新实验，只需写作

### 第四优先：D 线写作收口

9. 完成 `PREDICTION_D_WRITEUP_DRAFT.md`
10. 整合进主稿作为跨尺度身份保持的证据

### 第五优先：E 线桥接

11. 在 Discussion 中加 E→S_BD 桥接段

---

## 一句话结论

> 五条线在新主线下的实验工作已基本完成（A 95%、B 95%、C 90%、D 85%、E 90%），当前瓶颈全部在写作整合。最值得立刻做的是把 B 线的两层筛选论文写完并投稿，同时把 C/D 的证据总结作为支撑材料准备好。
