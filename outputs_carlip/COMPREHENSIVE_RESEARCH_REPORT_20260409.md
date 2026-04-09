# poset_phase 综合研究报告（截至 2026-04-09）

> 项目：`d:\Kiro\理论体系\poset_phase`  
> 主题：Layered Structural Screening of 4D Lorentzian Causal Sets  
> 报告范围：方向 α / β / γ / E（含最新 μ(N,κ) 强统计复跑与投稿准备）

---

## 1. 执行摘要（Executive Summary）

本轮研究已形成一个**可复现、可投稿、可扩展**的完整闭环：

1. **两层筛选架构被定量证实**：
   - 第一层（admissibility）负责几何可行性筛选；
   - 第二层（identity, `S_MD`）负责身份中心识别；
   - 两层不可简单加和，需序贯作用。

2. **方向 α（纯信息论替换）结论完成修正**：
   - 成功反驳“纯几何先验唯一决定结果”的强循环批评；
   - 但信息惩罚单独不足以全域唯一选出 Lor4D（补性定理成立）。

3. **方向 β 三分支全部闭环**：
   - β-1 序贯筛选：Lor4D 在 `N={10,20,40}` 均可达 #1（含 `N=40` 突破）；
   - β-2 `S_MD` 路径积分整合：`smd_only` 全 n #1（确认身份锚效应，附循环性边界）；
   - β-3 曲率评估：FLRW κ=1.0 在 `N=1024→1536` 出局率上升后于 `N=2048` 回落，判定为**有限尺度过渡现象**而非单调退出。

4. **方向 γ（工程化收口）完成**：
   - LaTeX 主稿与 submission 同步，编译通过；
   - 新增 arXiv metadata 草案与投稿 checklist；
   - 结论口径与报告已一致化。

5. **方向 E（最新扩展）完成**：
   - 新增并强化 `μ(N,κ)` 自适应参考原型；
   - 强统计复跑（`ref_reps=20, eval_reps=8`）确认：
     - flat 参考下 Lor4D 稳定 #1，FLRW #2；
     - adaptive 参考下 FLRW #1，Lor4D #2；
     - `||μ(N,κ)-μ(N,0)||₂` 从 `0.528` 增至 `0.558`。

---

## 2. 核心理论与方法框架

### 2.1 两层筛选（Layered Structural Screening）

- **Layer 1（Admissibility）**：
  以 `S_triple` + `S_BD` 进行一阶可行性筛选，剔除明显非几何族。
- **Layer 2（Identity）**：
  以 Mahalanobis 作用量
  $$
  S_{\mathrm{MD}}=(\mathbf I-\boldsymbol\mu)^\top\Sigma^{-1}(\mathbf I-\boldsymbol\mu)
  $$
  进行二阶身份识别。

### 2.2 关键判据

- Lor4D 排名（全家族排序）
- runner-up gap / identity basin deepening
- 曲率背景鲁棒性（de Sitter / Schwarzschild / FLRW）
- 信息论惩罚可解释性与非循环性

---

## 3. 方向 α：纯信息论替换实验（含修正）

### 3.1 主要结论

1. 用 5 项纯信息论惩罚替换几何惩罚后，Lor4D 在中等条件下可居前列，证明“非几何先验也能驱动 Lor4D 竞争力”。
2. **标准化伪影已定位并修正**：`group_by=[n]` 改为公平的 `group_by=[n,γ,mode]`。
3. 修正后：信息惩罚单独不能在全 `(N,γ)` 空间统一 #1，确认补性而非替代性。

### 3.2 关键技术发现

- `interval_diversity_deficit` 为主导判别项；
- 混合惩罚（geo+info）并非增益，反而常显著劣化（见 β 负结果）。

---

## 4. 方向 β：三分支完成情况

### 4.1 β-1 序贯筛选（正结果）

协议：`A2(γ1)` 预筛 → TOP-K → `A4(γ2)` 精选（`γ1`,`γ2` 独立）。

| N | 最优序贯结果 | 参数 | 对 A4-only 提升 |
|---|---|---|---|
| 10 | #1 | γ1=0, γ2=0.5, K=3 | +1 |
| 20 | #1 | γ1=0, γ2=0.5, K=3 | 0 |
| 40 | #1 | γ1=0.5, γ2=2.0, K=5 | +1 |

> 亮点：`N=40` 在单模态失败背景下实现突破，证明“软准入 + 硬精选”机制有效。

### 4.2 β-2 `S_MD` 路径积分整合（正结果，带边界）

作用量：
$$
A=-\log H+\gamma\,P+\alpha\,S_{\mathrm{MD}}
$$

- `smd_only` 在 `N={10,20,40}` 均可 #1；
- `info+smd` 可修复高 γ 下 info-only 的退化；
- 但 `S_MD` 参考 Lor4D 定义，存在先验循环性边界（定位为二层身份锚而非一层替代）。

### 4.3 β-3 曲率背景评估（完成）

- de Sitter / Schwarzschild：总体通过；
- FLRW κ=1.0：
  - `N=1024`: 43% fail
  - `N=1536`: 57% fail
  - `N=2048`: **14% fail**（7 seeds 合并）

> 结论更新：并非单调退出，而是**有限尺度过渡现象** + 尾部离群风险。

---

## 5. 方向 γ：文稿与工程收口

### 5.1 文稿状态

- `two_layer_paper/two_layer_screening.tex`：已集成 β 结果与最新口径；
- `two_layer_paper/two_layer_submission/two_layer_screening.tex`：同步更新；
- `outputs_carlip/MANUSCRIPT_SECTIONS_1_4.md`：完成到 `§6.11` 与结论更新。

### 5.2 编译与提交准备

- submission 目录可独立编译（15 页）;
- 图资源路径已修复（本地 `figures/`）；
- 新增：
  - `arxiv_metadata_draft.md`
  - `arxiv_submission_checklist.md`

---

## 6. 方向 E：μ(N,κ) 自适应参考原型（最新）

### 6.1 原型脚本与实验

- 脚本：`prototype_mu_kappa_adaptive.py`
- strong 复跑：`outputs_mu_kappa_prototype_strong/`
  - `mu_kappa_raw.csv`
  - `mu_kappa_summary.csv`
  - `mu_kappa_report.md`

### 6.2 强统计结果（κ=1.0）

| N | FLRW rank(flat) | FLRW rank(adapt) | Lor4D rank(flat) | Lor4D rank(adapt) | Δμ范数 |
|---|---:|---:|---:|---:|---:|
| 256 | 2 | 1 | 1 | 2 | 0.528 |
| 512 | 2 | 1 | 1 | 2 | 0.545 |
| 1024 | 2 | 1 | 1 | 2 | 0.558 |

**解释**：flat 与 adaptive 参考中心发生可重复交换，支持“参考流形族”扩展路线。

---

## 7. 证据链与可复现资产

### 7.1 主报告与主文稿

- `REPORT_INFO_PENALTY_EXPERIMENT_20260407.md`
- `outputs_carlip/MANUSCRIPT_SECTIONS_1_4.md`

### 7.2 关键实验脚本（摘）

- `experiment_info_ablation.py`
- `experiment_info_weight_scan.py`
- `experiment_info_hybrid.py`
- `analyze_sequential_enhanced.py`
- `experiment_smd_integration.py`
- `prototype_mu_kappa_adaptive.py`

### 7.3 关键数据产物（摘）

- `outputs_info_hybrid/sequential_screening_enhanced.csv`
- `outputs_smd_integration/smd_integration_summary.csv`
- `outputs_mu_kappa_prototype_medium/*`
- `outputs_mu_kappa_prototype_strong/*`

---

## 8. 里程碑时间线（近程）

- `66b933b`：Direction α 核心实验 + manuscript v1.0
- `e081132`：权重扫描 + z-score/raw reversal
- `37ec998`：混合惩罚负结果 + interval 物理验证
- `cb89ab0`：Report 审计修正 + Abstract
- `5004af2`：Direction β 三方向整合
- `2f74788`：LaTeX 集成 β 结果并编译
- `64aefea`：γ-1 FLRW 诊断
- `204f339`：γ-2 N=2048 全量收口
- `fb90a7f`：submission 一致性 + μ(N,κ) prototype
- `3dd7afd`：D1+D2 medium-grid + §6.11
- `08e5c65`：E1-E3（强复跑 + arXiv 准备 + tex 同步）

---

## 9. 综合结论（Final Synthesis）

1. **“两层筛选”从经验现象升级为可操作框架**：
   - 线性层做可行性筛；
   - 二次层做身份锚定；
   - 序贯优于加和。

2. **非循环性问题得到结构化回应**：
   - 纯信息论路径可提供独立证据；
   - 但全域唯一性仍依赖层间互补，而非单层万能。

3. **曲率问题从“失败/成功二元”升级为“尺度-背景依赖”表述**：
   - FLRW κ=1.0 在中尺度出现波动，超大 N 回落；
   - 更适合以“有限尺度过渡 + 尾部风险”描述。

4. **扩展方向已具工程抓手**：
   - `μ(N,κ)` 原型已可运行、可复跑、可回灌文稿；
   - 投稿资产已进入可执行清单阶段。

---

## 10. 下一步建议（优先级）

### P0（投稿前必做）
- 完成 arXiv 表单最终文案（plain text abstract / comments / categories）。
- 固定 submission 快照（tag + 文件清单锁定）。

### P1（论文增强）
- 在 `two_layer_submission.tex` 增补 `§6.11` 数值表（若期刊篇幅允许）。
- 对 overfull 行做最小排版修整（不改实质结论）。

### P2（研究继续）
- 扩展 κ 扫描：`κ ∈ {0.1,0.3,0.5,1.0,1.5}` + 固定强统计配置；
- 检验 `||Δμ||` 的 κ 标度律（是否近似线性或幂律）。

---

## 11. 当前状态判定

- **科学状态**：结论链条闭环，核心命题可 defend。
- **工程状态**：脚本-数据-文稿已连通，可重复。
- **投稿状态**：接近就绪（metadata 仍需最终定稿）。

> 结论一句话：
> **该项目已从“实验探索阶段”进入“可投稿与可扩展并行阶段”。**
