# Falsifiability Expansion Plan (v1)

**Date:** 2026-03-28  
**Scope:** 将主稿 §6.8 的 3 条可证伪条件扩展为可执行实验计划（含阈值、统计标准、失败判定、交付物）。

---

## 0. 目标与原则

本计划不追求“更多漂亮结果”，而追求：

1. 若理论错误，实验能把它打脸；
2. 若理论暂时成立，结论边界清楚、可重复、可审计；
3. 所有判断有预注册阈值，避免事后解释弹性。

---

## 1. 对应三条可证伪条件

### C1 — Identity-layer failure

> 若通过准入后，存在非 Lor4D 家族在扩展库中对 Lor4D 稳定压制（随 N 持续出现），则认同层稳定化主张失败。

### C2 — Background-response failure

> 若弱/中曲率 4D 背景在 metric-faithful 测试下系统性失去近邻层级（当前基线为 top-2），则 flat-centered local basin 解释失败。

### C3 — Basis-sufficiency failure

> 若存在同等可解释性的替代特征基，在同库同尺度下可重复且显著优于 $(d_{\mathrm{eff}}, C_1/C_0, w/N)$，则“tested library 内最小非冗余有效基”主张失败。

---

## 2. 预注册判定标准（必须先锁定）

### 2.1 显著失败（Hard Fail）

任一条件满足即判定对应主张失败：

- **C1-HF**：在扩展库中，存在某非 Lor4D 家族在至少 **3 个连续 N 档位**、且每档 **≥8/10 seeds** 下 rank #1 压过 Lor4D。
- **C2-HF**：在弱/中曲率窗口（预注册：de Sitter $H\le0.3$，FLRW $\kappa\le1.0$，Schwarzschild $\phi\le0.05$）中，4D curved family 在 **≥50% N 档位** 跌出 top-2，且跨 seed 可重复。
- **C3-HF**：替代基在同库同 N 上相对基线（当前 3 特征）实现：
  - 平均 margin 提升 **>20%**，且
  - rank #1 稳定性提升 **>10%**（或把边界 N 从 16 显著降至 ≤12），且
  - bootstrap 95% CI 不跨 0。

### 2.2 软失败（Soft Fail / 需降级措辞）

- 指标单次越线但不可重复（跨 seed / 跨配置不稳）。
- 仅在单一背景生成器下失败，跨模型不一致。

对应动作：不宣告理论失败，但必须在摘要/结论降级为“条件性成立”。

---

## 3. 实验矩阵（可执行）

## F1 — 扩展家族库抗压（对应 C1）

### 设计

- 基线库：25 families。
- 扩展目标：新增 12–20 家“近邻伪装”家族（优先围绕 Lor4D/Lor5D 边界）。
- N 网格：`N = 12,14,16,20,28,48,64,96,128,256`。
- Seeds：10；每 seed REPS=80（边界 N 可增至 120）。

### 指标

- Lor4D rank #1 率；
- 最近竞争者类型分布；
- 连续 N 压制段长度（对 C1-HF 核心）。

### 失败判定

按 C1-HF。

### 交付物

- `outputs_carlip/falsify_c1_family_pressure.md`
- `outputs_carlip/falsify_c1_family_pressure.json`

---

## F2 — Identity turn-on 稳定性再估计（C1 辅助）

### 设计

- 聚焦边界：`N = 10..24`（步长 2）。
- Seeds：20（由 10 提升）；REPS=120。
- 同步记录 min-margin、CI、runner-up census。

### 指标

- $N_{\mathrm{id}}$ 的后验区间；
- “连续稳定开启”判据（例如 3 档连续 10/10 且 min-margin>0）。

### 交付物

- `outputs_carlip/falsify_c1_turnon_refit.md`

---

## F3 — Metric-faithful 曲率测试（对应 C2）

### 设计

- de Sitter：保留现有管线，N 扩展到 `1536, 2048`（资源允许时 3072）。
- FLRW：从 proxy 向 metric-faithful 版本迁移（保留旧版对照）。
- Schwarzschild：弱场近似 + 更严格采样几何（明确坐标系与截断）。
- 每背景 3 强度（弱/中/强），每档 10 seeds。

### 指标

- top-2 维持率（按背景与 N 分层）；
- 与 Lor5D 的 rank-swap 频率；
- 曲率-分离斜率（是否仅精度放大而非结构重排）。

### 失败判定

按 C2-HF。

### 当前进展（2026-03-28 更新）

- `de Sitter` lowN split（`N=256,512`）当前结果：**Pass**。
- `Schwarzschild` lowN split（`N=256,512`）当前结果：**Pass**。
- `FLRW` lowN split（`N=256,512`）当前结果：**Hard Fail**，由 `\kappa=1.0` 触发（`fail_ratio = 0.50`）。

这意味着当前 C2 风险已从“理论上的可证伪条件”转为“已在特定背景中被触发的背景依赖性风险”，后续 highN 的首要任务是判定该信号是持续、恶化还是回落。

### 交付物

- `outputs_carlip/falsify_c2_background_response.md`
- `outputs_carlip/falsify_c2_background_response.csv`

---

## F4 — 替代特征基挑战赛（对应 C3）

### 候选基（预注册）

- B0（基线）：$(d_{\mathrm{eff}}, C_1/C_0, w/N)$
- B1：B0 + 局域区间尾部统计（如 $C_2/C_0$）
- B2：B0 + order-fraction 层级差分
- B3：B0 去掉一维并加一项独立可解释统计

### 规则

- 只允许“可解释且可复现实验定义”的特征；
- 禁止黑箱嵌入特征；
- 所有基在同一库、同一 N、同一 seeds 上横评。

### 指标

- rank #1 稳定率；
- margin 分布；
- $N_{\mathrm{id}}$ 变化；
- 复杂度惩罚（维数 + 计算成本）。

### 失败判定

按 C3-HF。

### 交付物

- `outputs_carlip/falsify_c3_basis_challenge.md`

---

## F5 — 反证导向合成对手（C1/C3 交叉）

### 设计

构造“匹配均值+协方差、但高阶矩不同”的合成家族，测试是否能长期穿透当前 identity layer。

### 目的

- 若可穿透：证明二阶层不完备，提示三阶层必要；
- 若不可穿透：增强“二阶层在当前库内充分”的可信度。

### 交付物

- `outputs_carlip/falsify_f5_synthetic_adversary.md`

---

## 4. 统计与报告协议

### 4.1 统计

- 所有主结论使用 bootstrap 95% CI；
- 关键比较同时给 effect size（不只 p-value）；
- 对多重比较做 FDR 控制（Benjamini-Hochberg）。

### 4.2 报告模板（每个 F* 必须一致）

1. 预注册阈值（复制粘贴）
2. 实验配置（N/seeds/REPS/版本哈希）
3. 原始结果（表 + 图）
4. 判定（Pass / Soft Fail / Hard Fail）
5. 结论动作（保留/降级/推翻哪条陈述）

---

## 5. 资源预算（初版）

- 预计总运行：~2–4 天（CPU 并行，按当前 N 上限）
- 优先级：`F1 > F3 > F4 > F2 > F5`
- 快速止损：若 F1 或 F3 出现 Hard Fail，立即暂停扩展并回到主稿降级措辞。

---

## 6. 与论文文本的联动规则

- 任一 Hard Fail：摘要与结论必须在同次提交中同步降级。
- 任一 Soft Fail：在 Limitations 增加条件性说明。
- 全部 Pass：允许保留当前“tested library 内稳健识别”口径，但不得升级为“普适底层原理”。
- 当前已执行动作：因 `FLRW(\kappa=1.0)` 在 lowN split 中触发 C2 Hard Fail，主稿与提纲中的曲率稳健性表述已降级为 **background-dependent robustness** 口径。

---

## 7. 立即执行顺序（建议）

1. 先跑 F1（最快暴露 C1 风险）
2. 并行准备 F3 的 metric-faithful 生成器
3. F4 在 F1/F3 后执行（避免无效特征搜索）

---

## 8. 当前状态标记

- 本文档为 `v1`（计划版）
- 下一步：生成 `config_falsify_c1.yaml` 与 `falsify_c1_runner.py` 的执行骨架
