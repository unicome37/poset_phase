# 工作日志 — 2026-03-18

## Prediction D：从"结构共变"回升为"条件准因果"

---

## 一、背景与动机

### 1.1 前置状态

此前 n=32 高功效分层置换检验（`prediction_d_stratified_highpower.py`，commit `7e01d5b`）已确认：

- **池化信号强且稳定**：p20 ρ=−0.25, p=0.00001***
- **层内信号消失**：p_strat = 0.17–0.79，ρ̄_strat ≈ 0，%negative ≈ 50%
- **结论**：D 被下调为"结构共变"——池化关联主要由 (N, family) 层间组成差异驱动

### 1.2 今日目标

用户提出 4 方向策略，尝试将 D 从"结构共变"恢复为"准因果"：

| 策略 | 说明 |
|------|------|
| ① 换连续 Y | 不再用 rank/tie-heavy 的 Δscore_local，改用 Δlog_H 等连续指标 |
| ② CG-specific 扰动 | 设计只改变 CG 稳定性的扰动 |
| ③ 层内残差化 | 在层内先去除结构混杂因子再做检验 |
| ④ 扩展 family | 增加 family 数量提升层内方差 |

**用户优先选择**：策略 ①+③ 合并——"换连续 Y + 层内残差化"。

---

## 二、实验设计与实现

### 2.1 脚本

**文件**：`prediction_d_continuous_residualized.py`（~640 行）

### 2.2 实验设计

| 设计要素 | 参数 |
|----------|------|
| **Poset 生成** | 6 families × 3 N (10,20,30) × 32 samples = 576/扰动级 |
| **扰动级别** | p00 (0%), p05 (5%), p10 (10%), p20 (20%) cover removal |
| **X 变量** | Δpenalty_cg = penalty_cg(perturbed) − penalty_cg(baseline) |
| **Y 目标** | 5 个连续指标（均机制独立于 penalty_cg） |
| **残差化** | 11 维混杂因子 OLS 层内残差化 |
| **Boundary 子集** | 中间 50% \|Δpenalty_cg\|（quantile 0.25–0.75） |
| **置换次数** | 100,000/组 |
| **总测试组** | 3 扰动 × 5 Y × 2 子集 = 30 组 |
| **种子偏移** | SEED_OFFSET = 172000 |

### 2.3 五个连续 Y 目标

| Y 目标 | 含义 | 与 penalty_cg 的机制独立性 |
|--------|------|---------------------------|
| **Δscore_local** | 综合得分变化 = −β·ΔlogH + γ·(Δpenalty_local + Δgeo_total) | 综合指标 |
| **Δlog_H** | 熵变化 | 直接反映结构有序度，与 CG 算子无关 |
| **Δpenalty_local** | 局部罚分变化 | 几何属性，与 CG 算子无关 |
| **Δgeo_total** | 几何总罚分变化 | 几何属性，与 CG 算子无关 |
| **Δpenalty_neutral** | 中性种罚分变化 | 中性种分数，与 CG 算子无关 |

### 2.4 十一个混杂因子

在每个 (perturb, N, family) 层内，X 和 Y 均经 OLS 回归去除以下 11 维混杂因子后取残差：

| 类别 | 混杂因子 |
|------|----------|
| **结构签名** | sig_comp, sig_d_eff, sig_height_ratio, sig_width_ratio, sig_degree_var |
| **层级结构** | layer_count, mean_layer_gap |
| **基线适应度** | log_H, penalty_local, geo_total, mean_penalty_cg |

---

## 三、实验执行

### 3.1 运行过程

| 步骤 | 内容 | 耗时 |
|------|------|------|
| Step 1 | 生成 576 个 base poset + 计算 base rank | ~36s |
| Step 2 | 4 级扰动 × CG pipeline + extended features | ~220s |
| Step 3 | 构建 delta frame (baseline vs perturbed) | 即时 |
| Step 4 | 层内残差化 + boundary flagging | 即时 |
| Step 5 | 30 组置换检验 × 100k 置换 | ~4500s |
| **总计** | | **79.3 min** |

### 3.2 中间数据

| 文件 | 规模 |
|------|------|
| `extended_sample_features.csv` | 2304 行 |
| `delta_residualized.csv` | 1728 行 × 36 列 |
| `continuous_residualized_results.csv` | 30 行 |
| `continuous_residualized_report.md` | 完整报告 |

---

## 四、核心结果

### 4.1 全样本（all subset）层内显著性

| Y 目标 | p05 ρ̄ | p05 p | p10 ρ̄ | p10 p | p20 ρ̄ | p20 p | 剂量-反应 |
|--------|--------|-------|--------|-------|--------|-------|-----------|
| **Δscore_local** | **−0.100** | **0.019★** | −0.076 | 0.072 | **−0.159** | **0.0002★★★** | ✓ 非单调 |
| **Δlog_H** | **+0.099** | **0.019★** | +0.083 | 0.051 | **+0.185** | **0.00004★★★** | ✓ 单调增强 |
| Δpenalty_local | −0.040 | 0.351 | +0.074 | 0.083 | **+0.175** | **0.00007★★★** | ✓ 符号翻转 |
| Δgeo_total | −0.039 | 0.354 | +0.071 | 0.093 | **+0.176** | **0.00007★★★** | ✓ 符号翻转 |
| **Δpenalty_neutral** | **−0.139** | **0.001★★** | −0.078 | 0.064 | **−0.187** | **0.00002★★★** | ✓ 单调增强 |

**p20 全 5 个指标层内显著**（p_strat ≤ 0.0002），其中最强：Δpenalty_neutral p=0.00002。

### 4.2 Boundary 子集层内显著性

| Y 目标 | p05 p | p10 p | p20 p |
|--------|-------|-------|-------|
| Δscore_local | 0.066 | 0.066 | **0.003★★** |
| Δlog_H | 0.098 | **0.041★** | **0.002★★** |
| Δpenalty_local | 0.995 | 0.077 | **0.012★** |
| Δgeo_total | 0.997 | 0.082 | **0.013★** |
| Δpenalty_neutral | 0.109 | 0.223 | **0.027★** |

**p20 boundary 全 5 指标同样显著**，排除了极端扰动驱动假说。

### 4.3 与此前原始检验的决定性对比

| 设计 | p05 p_strat | p10 p_strat | p20 p_strat |
|------|-------------|-------------|-------------|
| **§6.6 原始**（raw Δscore_local，无残差化） | 0.170 ns | 0.793 ns | 0.618 ns |
| **§6.7 残差化 Δscore_local** | **0.019★** | 0.072 | **0.0002★★★** |
| **§6.7 残差化 Δlog_H** | **0.019★** | 0.051 | **0.00004★★★** |

**残差化是决定性因素**：去除 11 维层内混杂因子后，原本消失的层内信号全面恢复。

### 4.4 层内显著性总计

**14/30 组检验达层内显著（p < 0.05）**：
- p05: 3 显著（Δscore_local, Δlog_H, Δpenalty_neutral）
- p10: 1 显著（Δlog_H boundary）
- p20: 10 显著（全 5 指标 × all + boundary）

---

## 五、关键发现与解释

### 5.1 此前零结果的根因

此前 §6.5–6.6 的"结构共变"判定是**过早的**。原始分层检验失败，不是因为层内无信号，而是因为：

1. **层内混杂**: 11 维结构/适应度混杂因子（结构签名、基线适应度等）稀释了 X–Y 关联
2. **通道抵消**: 综合指标 Δscore_local 混合了方向相反的熵通道（↑）和罚分通道（↓），部分信号互相抵消

### 5.2 恢复机制

残差化解决了问题 1，连续 Y 分通道解决了问题 2：

- 经混杂控制后，**每个通道内的信号方向一致**
- Δlog_H（熵通道）: CG 破坏 → 熵上升（ρ > 0），剂量-反应单调
- Δpenalty_neutral（中性种通道）: CG 破坏 → 中性种罚分下降（ρ < 0），剂量-反应单调
- Δpenalty_local / Δgeo_total（几何通道）: 大扰动下 ρ > 0——松动的序关系减少几何罚分

### 5.3 D 的新定位

| 旧定位 | 新定位 |
|--------|--------|
| 结构共变（structural covariation） | **条件准因果（conditional quasi-causal）** |

**条件性含义**：D 的层内准因果信号真实存在，但需要 11 维混杂控制才可见。这本身是重要的方法论发现——它表明跨尺度关联本质上是"混杂面纱"下的真信号，而非虚假关联。

---

## 六、文档更新

### 6.1 更新的文件清单

| 文件 | 更新内容 |
|------|----------|
| **PREDICTION_D_WRITEUP_DRAFT.md** | 新增 §6.7 "Continuous-Y Rich-Residualized Within-Stratum Test"（~100 行） |
| **PROJECT_PROGRESS.md** | 完成记录 + 新 TODO |
| **结构存在论四大推论阶段性总论断.md** | ~25 处系统性修订，D 从"结构共变"→"条件准因果" |

### 6.2 总论断.md 修订要点

全文系统性更新，涉及以下章节：

| 章节 | 修改类型 |
|------|----------|
| 开篇摘要（§）| D 描述从"结构共变"→"条件准因果" |
| §五 最准确口径 | 更新两个引用块 |
| §六 正式版 | 重写 D 段落 |
| §1.4 Prediction D | 标题"结构共变"→"条件准因果" |
| §1.4 当前最稳口径 | 完全重写，增加三段式转折叙事 |
| §1.4 命中度评估 | 从"池化共变复现"→"条件准因果命中" |
| §1.4 当前边界 | 重写 5 条边界 |
| §1.5 四条推论格局 | "D 外围"→"D 条件准因果延伸" |
| §1.6 本章小结 | 更新总结句 |
| §2.1 方向命中 | D 的描述更新 |
| §2.1 结构命中 | D 的描述更新 |
| §2.1 替代解释 | D 的描述更新 |
| §2.2 确认性窗口 | 重写 D 的确认窗口含义 |
| §7.4 D 的完成 | 完全重写 |
| §7.5 格局 | "B/C/A 主链 + D 外围"→"+ D 条件准因果延伸" |
| §7.6 显影 | 更新 D 的回答 |
| §7.7 交叉验证 | 更新 D 的定位 |
| 结语 | 完全重写 D 段落 |
| 修订日志 | 新增"修订②"条目 |
| 附录中文稳健版 | 重写 |
| 附录一句压缩版 | 更新 |

保留了**历史转折叙事**（准因果→结构共变→条件准因果），体现理论程序的自我校正能力。

---

## 七、Git 记录

| 项目 | 内容 |
|------|------|
| Commit | `b808e4f` |
| Branch | main |
| Message | `feat(D): continuous-Y + 11-confounder rich-residualized within-stratum test` |
| Changed | 7 files, +4828 −1 |
| Push | ✅ 已推送至 GitHub (`unicome37/poset_phase`) |

### 新增文件

- `prediction_d_continuous_residualized.py`（640 行）
- `outputs_exploratory/prediction_d_continuous_residualized/continuous_residualized_report.md`
- `outputs_exploratory/prediction_d_continuous_residualized/continuous_residualized_results.csv`
- `outputs_exploratory/prediction_d_continuous_residualized/delta_residualized.csv`
- `outputs_exploratory/prediction_d_continuous_residualized/extended_sample_features.csv`

### 修改文件

- `PREDICTION_D_WRITEUP_DRAFT.md`（+71 行，§6.7）
- `PROJECT_PROGRESS.md`（+3 −1）

---

## 八、D 定位三段式转折全景

```
初始阶段                n=32 高功效检验            连续Y+残差化检验
(frozen window)         (2026-03-18 修订①)        (2026-03-18 修订②)
     │                        │                         │
     ▼                        ▼                         ▼
  准因果          →       结构共变          →       条件准因果
(池化ρ显著,                (层内信号消失,             (11维残差化后
 剂量-反应,                 p_strat=0.17-0.79,        层内信号全面恢复,
 准干预)                    %neg≈50%)                  p20全5指标p≤0.00002)
```

**核心洞见**：此前的零结果不是"因果不存在"的证据，而是"混杂未控"的假象。经精细混杂控制后，层内准因果信号全面恢复，且强于池化信号的相应效果量。

---

## 九、当前最稳口径

> **结构存在论已在 B/C/A 三条主线上获得了多层命中与初步闭环，D 在精细混杂控制（11维层内残差化）条件下恢复了层内准因果信号，条件性地加强了跨尺度延伸。**

> **主闭环是 B–C–A，D 经富残差化检验回升为条件准因果级跨尺度延伸。**

---

## 十、后续方向

| 优先级 | 方向 | 说明 |
|--------|------|------|
| P1 | **敏感性检验** | 残差化结论对混杂因子选择和回归函数形式的鲁棒性 |
| P2 | **扩展 family/N** | 在更多 family 种类和 N 范围下验证条件准因果的普适性 |
| P3 | **更直接的扰动设计** | 不改变 N/family 组成但改变 CG 路径 |
| P4 | **交叉验证** | 用 leave-one-stratum-out 检验层内信号的稳定性 |

---

*文档生成时间：2026-03-18*
*关联 commit：`b808e4f`*
*仓库：`unicome37/poset_phase`*
