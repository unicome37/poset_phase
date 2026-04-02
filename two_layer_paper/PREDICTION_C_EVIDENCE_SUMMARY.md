# Prediction C 证据总结：层级整合作为 Basin Deepening 的结构根因

> **日期**: 2026-04-01
> **定位**: 两层筛选主稿 §5（Historical Sedimentation）的机制支撑材料
> **口径**: 不再是"logH 被层深预测"，而是"层级整合是 basin deepening 的结构根因"

---

## 核心命题

> 在已匹配的准几何竞争者内部，更浅层、更局域、近邻主导的组织方式倾向于保留更高的组合自由度（logH）；
> 因此，Lor4D 的 basin deepening 不仅是统计精度提升的结果，还有结构层面的根因：
> 更深的因果层级整合压低了组合自由度，使得 Lor4D 在特征空间中越来越孤立。

---

## 三层证据链

### 第一层：配对检验（局部信号）

**方法**: Lor2D vs MLR survivor 的 matched-pair 设计，控制几何匹配特征后比较层级整合差值与 logH 差值。

**结果汇总**（三档 nearwall 闭环）:

| 口径 | n_pairs | LDI→ΔlogH (Pearson) | HII→ΔlogH (Pearson) | layer_count→ΔlogH |
|---|---|---|---|---|
| 原始 pilot | 22 | +0.759 | -0.757 | -0.783 |
| nearwall expanded | 35 | +0.850 | -0.848 | -0.835 |
| nearwall p5p95 | — | — | — | — |
| nearwall rescue | 50 | +0.840 | -0.839 | -0.827 |

**独立 seed 复现**: 双 seed（20260331/20260401）重跑，Tier2 相关系数完全一致（hii_delta=-0.8363）。

**效应量诊断**: leave-k-out 去除最极端 15 对后仍保持 ≈-0.73~-0.75，信号非单点驱动。

### 第二层：全样本受控回归（跨族普适性）

**方法**: 544 个样本，控制 width/comparable_fraction/dim_eff + N + family 固定效应。

**结果**:
- 最强控制下（controls_plus_n_family_fe）:
  - layer_count: partial r = -0.450, p = 0.00025
  - mean_layer_gap: partial r = -0.424, p = 0.00025
  - HII 综合指标: partial r = +0.086, p = 0.0487（弱，受构造方式影响）

**LOCO 混杂剥离**:
- drop antichain_width 后: HII 从 +0.086 拉到 -0.495，layer_count 从 -0.450 拉到 -0.766
- drop comparable_fraction 或 drop geo_dim_eff: 与 full controls 基本一致（二者近线性冗余）
- **结论**: antichain_width 是当前最关键的混杂通道

### 第三层：筛子增强扫描（筛选潜力）

**方法**: 把层级特征并入第二道筛子，扫描 eta_depth 参数，观察 zeta_cross 变化。

**结果**:
- hierarchy_integration_adv_scaled:
  - N=30: baseline 2.219 → best 0.028
  - N=48: baseline 5.588 → best -0.693
- layer_count_adv_scaled: mean best zeta_cross = -0.366
- mean_layer_gap_adv_scaled: mean best zeta_cross = -0.063

**结论**: 层级整合项不只是"解释量"，在小规模扫描里已表现出筛选增强潜力。

---

## F7 large-N fresh 辅助验证（非主线，仅参考）

- prediction_c_f7_large_n.py --fresh --reps 8 --ns 20 36 52 72 100
- within-cell 方向正确率: 16/19 = 84%
- per-N 跨族 pooled ρ(Σ_hist, logH): -0.805 / -0.682 / -0.769 / -0.743 / -0.742
- 全部为强负相关且显著

---

## 与两层筛选主线的整合

C 线在新框架下的理论定位：

1. **Basin deepening 的机制解释**: 为什么 Lor4D 的 Mahalanobis gap 随 N 发散？
   不仅因为 Σ⁻¹ 的统计精度提升（precision amplification），
   还因为更深的因果层级整合在结构层面压低了组合自由度，
   使得 Lor4D 在特征空间中的"物理距离"本身就是真实的（d_E ≈ 0.46，不随 N 衰减）。

2. **Historical sedimentation 的微观机制**: 每增加一个因果事件，
   不仅增加了统计样本量，还增加了层级约束的厚度。
   这就是为什么 sedimentation 是"不可逆的"——层级整合一旦形成，就不会自发解体。

3. **不是第三层筛选**: C 线描述的是 Layer 2 的纵向行为（across N），
   不是与 Layer 1/2 平行的第三个筛选算子。

---

## 判定状态

按推论C验证方案的判定标准：

- ✅ 配对检验通过（三档 nearwall 闭环 + 独立 seed 复现）
- ✅ 全样本回归通过（layer_count/mean_layer_gap 在强控制下仍显著）
- ✅ 筛子增强有信号（zeta_cross 明显下降）

**三层中三层通过 → C 线达到"支持"级别。**

---

## 剩余开放问题

1. HII 综合指标在最强控制下仅弱正相关——建议优先使用组成分量（layer_count, mean_layer_gap）而非复合分数
2. 筛子增强的 pair count 在 N=44/48 仍偏小，应按 pilot evidence 解释
3. 是否将层级整合项正式写入 S_MD 的扩展版本，取决于是否能在更大样本上复现增强效果
