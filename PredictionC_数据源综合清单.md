# Prediction C 数据源综合清单

> 生成时间: 2025-07  
> 目的: 盘点 `poset_phase` 项目中所有可用于扩展 Prediction C 样本量的 CSV 数据源

---

## 一、Prediction C 分析背景

**核心假设**: 更深的因果层级整合系统性地降低组合熵 log_H。

**关键观测量**:
- `layer_count` — 层级数
- `mean_layer_gap` — 平均层间距
- `long_edge_fraction` — 长边比例
- `reduction_edge_density` — 约化边密度
- `cover_density` — 覆盖密度
- `layer_signature_redundancy` — 层签名冗余度
- 综合指标: `hierarchy_integration_index`

**响应变量**: `log_H` (组合熵对数)

---

## 二、数据分级

### GOLD 级: 完整包含全部 6 个层级特征 + log_H

| 数据集 | 行数 | N 值 | 族 | 文件路径 |
|--------|------|------|-----|---------|
| prediction_c_pairwise_validation | 22 | 30,40,44,48 | mlr vs lor2d (配对 Δ 格式) | `outputs_exploratory/prediction_c_pairwise_validation/prediction_c_pairwise_validation_raw.csv` |
| hierarchy_integration_validation | 44 | 30,40,44,48 | lor2d, mlr | `outputs_exploratory/hierarchy_integration_validation/hierarchy_integration_validation_raw.csv` |
| locality_dominance_validation | 44 | 30,40,44,48 | lor2d, mlr | `outputs_exploratory/locality_dominance_validation/locality_dominance_validation_raw.csv` |
| matched_residual_freedom_check | 44 | 30,40,44,48 | lor2d, mlr | `outputs_exploratory/matched_residual_freedom_check/matched_residual_freedom_raw.csv` |

**小计**: ~132 行独立样本（去重后可能更少），仅覆盖 **2 族 × 4 个 N 值**。

**列结构** (hierarchy_integration_validation 为例):
```
family, n, seed, log_H, antichain_width, comparable_fraction, geo_dim_eff,
layer_count, mean_layer_size, layer_size_std, cover_density,
reduction_edge_count, reduction_edge_density, layer_signature_redundancy,
long_edge_fraction, mean_layer_gap,
z_layer_count, z_mean_layer_gap, z_long_edge_fraction, z_reduction_edge_density,
hierarchy_integration_index
```

**prediction_c_pairwise_validation 列结构** (配对 Δ 格式):
```
n, mlr_seed, lor_seed,
log_H_delta_mlr_minus_lor2d, score_A2_gamma_delta_mlr_minus_lor2d,
layer_count_mlr, mean_layer_gap_mlr, ..._lor2d,
layer_count_delta_mlr_minus_lor2d, ...,
locality_dominance_delta_index, hierarchy_integration_delta_index
```

---

### SILVER 级: 部分 Prediction C 特征 + log_H

| 数据集 | 行数 | N 值 | 族 | 含有的 pred_c 特征 |
|--------|------|------|-----|-------------------|
| mlr_survivor_profile | 90 | 30,40,44 | mlr | layer_count, comparable_fraction, antichain_width, geo_dim_eff |
| controlled_compressibility_wide | 72 | 30,40 | abs_lay, int_ord, lor2d, mlr | comparable_fraction, antichain_width |
| compressibility_validation | 12 | 30,40,44 | 4 families | comparable_fraction, antichain_width |

**小计**: ~174 行，缺失部分核心层级特征。

---

### BRONZE 级: 有 log_H + 几何代理特征 (geo_dim_eff, geo_cover_density, geo_layer_smoothness)

需对原始 poset 追加计算层级特征后才能用于 Prediction C。

#### prediction_a 系列

| 数据集 | 行数 | N 值 | 族 | 文件路径 |
|--------|------|------|-----|---------|
| prediction_a_dim_replacement_sp8 | 8,064 | 20-40 | KR, lor2d, lor3d, lor4d | `outputs_exploratory/.../prediction_a_ablation_raw.csv` |
| prediction_a_geometric_ablation | 5,600 | 20-36 | KR, lor2d, lor3d, lor4d | 同上 |
| prediction_a_dim_replacement | 2,240 | 20-36 | KR, lor2d, lor3d, lor4d | 同上 |
| prediction_a_pilot_fast | 2,240 | 20-36 | KR, lor2d, lor3d, lor4d | `prediction_a_raw.csv` |
| prediction_a_dim_replacement_n44_n48 | 1,344 | 44,48 | KR, lor2d, lor3d, lor4d | |
| prediction_a_n52_mixed | 336 | 52 | KR, lor2d, lor3d, lor4d | |
| prediction_a_n56_mixed | 336 | 56 | KR, lor2d, lor3d, lor4d | |
| prediction_a_n60_mixed | 336 | 60 | KR, lor2d, lor3d, lor4d | |
| prediction_a_n64_mixed | 336 | 64 | KR, lor2d, lor3d, lor4d | |
| prediction_a_n68_mixed | 336 | 68 | KR, lor2d, lor3d, lor4d | (top-level outputs_exploratory) |
| prediction_a_n72_mixed | 336 | 72 | KR, lor2d, lor3d, lor4d | |
| prediction_a_seed_sensitivity_n64 | 1,008 | 64 | KR, lor2d, lor3d, lor4d | |
| prediction_a_seed_sensitivity_n68 | 1,008 | 68 | KR, lor2d, lor3d, lor4d | |
| prediction_a_seed_sensitivity_n72 | 1,008 | 72 | KR, lor2d, lor3d, lor4d | |

**pred_a 小计**: ~24,528 行, N=20-72, 4 族

**列结构**:
```
n, family, sample_id, seed, gamma, beta, variant, log_H, entropy_method,
penalty_neutral, penalty_geometric, penalty_effective,
score, score_per_node, score_per_nlogn,
geo_dim_eff, geo_cover_density, geo_layer_smoothness, geo_width_height, geo_dim_proxy_penalty, ...
```

#### width_height_consistency_scale_scan 系列

| 数据集 | 行数 | N 值 | 族 | 位置 |
|--------|------|------|-----|------|
| scale_scan (base) | 2,016 | 20-44 | KR, lor2d | top-level outputs_exploratory |
| scale_scan_low | 3,528 | 20-44 | KR, lor2d | top-level outputs_exploratory |
| scale_scan_mid | 2,268 | 36-44 | KR, lor2d | top-level outputs_exploratory |
| scale_scan_mid_sp4 | 4,752 | 36-44 | KR, lor2d | top-level outputs_exploratory |
| scale_scan_n44_fine | 1,152 | 44 | KR, lor2d | top-level outputs_exploratory |
| scale_scan_n48_fine | 576 | 48 | KR, lor2d | poset_phase/outputs_exploratory |
| scale_scan_n48_sp2 | 756 | 48 | KR, lor2d | poset_phase/outputs_exploratory |
| scale_scan_n48_n52_mixed | 648 | 48,52 | KR, lor2d | poset_phase/outputs_exploratory |
| scale_scan_n52_fine_mixed | 504 | 52 | KR, lor2d | poset_phase/outputs_exploratory |
| scale_scan_n52_lower_mixed | 504 | 52 | KR, lor2d | poset_phase/outputs_exploratory |
| scale_scan_n52_ultrafine_mixed | 396 | 52 | KR, lor2d | poset_phase/outputs_exploratory |

**scale_scan 小计**: ~17,100 行, N=20-52, 2 族

#### 其他 BRONZE 数据集

| 数据集 | 行数 | N 值 | 族 |
|--------|------|------|-----|
| geometric_ablation_gamma_c | 3,276 | 20-44 | KR, lor2d |
| noncyclic_dim_replacement_gamma_c | 1,512 | 20-44 | KR, lor2d |
| seed_sensitivity_c_threshold (多变体) | ~4,000+ | 44-52 | KR, lor2d |

**BRONZE 级总计**: ~50,000+ 行

---

### 基础级: 仅有 log_H（无任何层级/几何特征）

| 数据集 | 行数 | N 值 | 族数 | 备注 |
|--------|------|------|------|------|
| outputs_frozen_exact | 11,520 | 10-16 | 8 | 最大精确熵单源 |
| outputs_smallN_exact | 10,080 | 10-16 | 7 | frozen_exact 子集 |
| medium_exact | 864 | 20-40 | 3 | |
| medium_exact_scan | 1,134 | 20-44 | 3 | |
| outputs_medium | 1,512 | 12-20 | 7 | |
| outputs_gamma_scan | 1,176 | 20 | 7 | |
| outputs_smoke | 84 | 10 | 7 | 测试用 |
| CG 系列 (6 datasets) | ~1,020 | 16-30 | 5 | 粗粒化 |

**基础级总计**: ~27,000+ 行

---

## 三、数据层次拓扑图

```
N 值:     10  12  14  16  20  24  28  32  36  40  44  48  52  56  60  64  68  72
─────────────────────────────────────────────────────────────────────────────────
GOLD (2族)                              ████████████████████
SILVER                                  ████████████████
BRONZE-predA                            ████████████████████████████████████████
  (4族)
BRONZE-scan                             ████████████████████████████████
  (2族)
基础-frozen                 ████████████
  (8族)
基础-medium                     ████████████████████████████
  (3族)
```

---

## 四、扩样路径建议

### 路径 A: 重算层级特征（推荐、性价比最高）

对 BRONZE 级数据中已保存的 poset 实例重新计算 `layer_count`, `mean_layer_gap`, `long_edge_fraction`, `reduction_edge_density`, `cover_density`, `layer_signature_redundancy`。

- **潜在收益**: 样本量从 ~132 → ~50,000+
- **所需工作**: 编写特征计算管道，遍历已有 poset pickle/adjacency 文件
- **前提条件**: BRONZE 级数据的原始 poset 结构是否已保存（需确认）

### 路径 B: 生成新实验数据

针对缺失族 (abs_lay, int_ord, tran_perc, lor3d, lor4d) 和更多 N 值运行新的生成 + 熵计算 + 特征提取全流程。

- **潜在收益**: 补全族覆盖 + N 范围
- **所需工作**: 计算成本高（精确熵计算是瓶颈）

### 路径 C: 混合策略

先用路径 A 快速扩充 BRONZE → GOLD，再用路径 B 定向补充缺失族。

---

## 五、注意事项

1. **重复数据**: `outputs_confirmatory/frozen_exact/raw_samples.csv` 与 `outputs_frozen_exact/raw_samples.csv` 完全相同（11,520 行），合并分析时需去重。
2. **族覆盖不均**: GOLD 级仅含 lor2d + mlr；BRONZE 级最多 4 族；只有基础级覆盖全部 8 族。
3. **配对格式**: prediction_c_pairwise_validation 采用 Δ(mlr - lor2d) 格式，需特殊处理。
4. **文件位置分布**: 数据分散在 `d:\Kiro\outputs_exploratory\` 和 `d:\Kiro\理论体系\poset_phase\outputs_exploratory\` 两个目录树中。

---

## 六、文件位置索引

### 顶层目录: `d:\Kiro\outputs_exploratory\`
包含: prediction_c_pairwise_validation, hierarchy_integration_validation, locality_dominance_validation, matched_residual_freedom_check, mlr_survivor_profile, mlr_survivor_nearest_family_check, controlled_compressibility_wide, compressibility_validation, pairwise_compressibility_duel, pairwise_locality_delta_validation, width_height_consistency_scale_scan (base/low/mid/mid_sp4/n44_fine), prediction_a_n68_mixed, prediction_a_seed_sensitivity_n64

### poset_phase 目录: `d:\Kiro\理论体系\poset_phase\outputs_exploratory\`
包含: prediction_a (多个变体), scale_scan (n48-n52 系列), geometric_ablation_gamma_c, noncyclic_dim_replacement_gamma_c, seed_sensitivity variants, prediction_a_seed_sensitivity (n68, n72)

### 冻结/确认: `d:\Kiro\理论体系\poset_phase\`
包含: outputs_frozen_exact, outputs_frozen_cg, outputs_confirmatory, outputs_smallN_exact, outputs_medium, outputs_gamma_scan, outputs_smoke, outputs_cg_*

---

## 七、BRONZE 增强试点结果 (prediction_c_bronze_augmented_pilot)

### 7.1 增强管道

- **脚本**: `augment_prediction_c_features.py` (184 行)
- **配置**: `config_augment_prediction_c_features_bronze_pilot.yaml` (25 个源 CSV)
- **原理**: 通过 `(family, n, seed)` 三元组重建 poset，调用 `residual_metrics()` 计算全部 6 项层级特征，合并回原 BRONZE 数据
- **输出目录**: `outputs_exploratory/prediction_c_bronze_augmented_pilot/`

### 7.2 增强数据概况

| 指标 | 值 |
|------|-----|
| 总行数 | 41,628 |
| 总列数 | 48 |
| 唯一 poset (family×n×seed) | 544 |
| 族数 | 4 (KR_like, lor2d, lor3d, lor4d) |
| N 范围 | 20–72 |
| Prediction C 六特征覆盖率 | **100%** (无缺失) |

**唯一 poset 族分布**: KR_like: 148, lor2d: 148, lor3d: 124, lor4d: 124

**Family × N 分布矩阵** (单元格 = 行数):

| N | KR_like | lor2d | lor3d | lor4d |
|---|---------|-------|-------|-------|
| 20 | 1,482 | 1,482 | 1,482 | 1,482 |
| 24 | 1,482 | 1,482 | 1,482 | 1,482 |
| 28 | 1,482 | 1,482 | 1,482 | 1,482 |
| 32 | 1,482 | 1,482 | 1,482 | 1,482 |
| 36 | 2,406 | 2,406 | 504 | 504 |
| 40 | 756 | 756 | 504 | 504 |
| 44 | 336 | 336 | 336 | 336 |
| 48 | 336 | 336 | 336 | 336 |
| 52 | 84 | 84 | 84 | 84 |
| 56 | 84 | 84 | 84 | 84 |
| 60 | 84 | 84 | 84 | 84 |
| 64 | 252 | 252 | 252 | 252 |
| 68 | 252 | 252 | 252 | 252 |
| 72 | 252 | 252 | 252 | 252 |

### 7.3 输出文件清单

| 文件 | 说明 |
|------|------|
| `prediction_c_augmented_combined.csv` | 合并后完整数据 (41,628 行 × 48 列) |
| `prediction_c_feature_cache.csv` | 特征缓存 (544 唯一 poset 的层级指标) |
| `prediction_c_augmented_source_summary.csv` | 按源文件汇总 |
| `prediction_c_augmented_family_summary.csv` | 按族汇总 |
| `prediction_c_augmented_n_summary.csv` | 按 N 值汇总 |

---

## 八、Simpson's Paradox 诊断

### 8.1 三层相关性分解

对增强数据中 `layer_count ↔ log_H` 关系进行三级分解:

| 分析层级 | Pearson r | 说明 |
|----------|-----------|------|
| **整体** (raw pooled) | -0.049 | 接近零——N 主效应淹没层级信号 |
| **控制 N** (within-N pooled) | **-0.892** | 非常强——但含族间差异 |
| **控制 N + 族** (within-family-within-N) | -0.31 ~ -0.51 | **真实效应**，中等强度 |

### 8.2 Simpson's Paradox 机制

**族间均值差异 (N=32 为例)**:

| 族 | 平均 log_H | 平均 layer_count |
|-----|-----------|-----------------|
| lor2d | 41.69 | 8.04 |
| lor3d | 52.35 | 4.21 |
| lor4d | 62.72 | 3.55 |
| KR_like | 43.20 | 3.00 (常量) |

lor2d **同时**具有最低 log_H 和最高 layer_count → 控制 N 但不控制族时，族间差异膨胀 r 值至 -0.89。

**KR_like 退化问题**: KR_like 在所有 N 值下 `layer_count = 3.0`（零变异），无法提供任何 within-family 信号。

### 8.3 真实 within-family-within-N 相关性

| 族 | layer_count | mean_layer_gap | long_edge_fraction |
|----|------------|----------------|-------------------|
| lor2d | **-0.508** | -0.476 | -0.357 |
| lor3d | **-0.354** | -0.446 | -0.495 |
| lor4d | **-0.305** | -0.501 | -0.468 |
| KR_like | N/A (零变异) | +0.643 (反向!) | +0.643 (反向!) |

**结论**: 三个 Lorentzian 族均支持 Prediction C 假设（deeper hierarchy → lower log_H），效应量中等 (|r| ≈ 0.3–0.5)。KR_like 是退化案例，不适合参与 within-family 分析。

### 8.4 reduction_edge_density 符号反转

| 分析层级 | r(reduction_edge_density, log_H) |
|----------|--------------------------------|
| 控制 N (pooled) | **+0.36** |
| within-lor3d-within-N | -0.25 |
| within-lor4d-within-N | **-0.64** |

**解读**: 池化方向为正，但 within-family 方向为负——又一个 Simpson's Paradox 实例。真实效应方向与 Prediction C 一致（更密的约化图 → 更低的 log_H）。

### 8.5 分析建议

1. **必须按族分层**: 任何 Prediction C 回归/检验都必须在族内进行，或加入族固定效应
2. **排除 KR_like 的层级回归**: KR_like 无 layer_count 变异，仅可用于 cross-family 配对设计
3. **推荐模型**: `log_H ~ layer_count + mean_layer_gap + long_edge_fraction + N | family` (分族回归或混合效应模型)
4. **配对扩展**: 544 唯一 poset 可构建远多于现有 22 对的跨族配对样本

---

## 九、数据层级更新总结

| 层级 | 更新前行数 | 更新后行数 | 族覆盖 | 用途 |
|------|-----------|-----------|--------|------|
| GOLD | ~132 | ~132 | 2 (lor2d, mlr) | 参考基准 |
| SILVER | ~174 | ~174 | 4 | 辅助 |
| **BRONZE → 增强 GOLD** | 0 | **41,628** | **4** (KR, lor2d, lor3d, lor4d) | **主力分析** |
| 基础 | ~27,000 | ~27,000 | 8 | 特征扩展候选 |
| **有效分析总量** | ~132 | **~41,760** | — | ×316 倍扩增 |

---

## 十、跨族配对验证深度诊断 (prediction_c_bronze_matched_validation)

### 10.1 管道概况

- **脚本**: `prediction_c_bronze_matched_validation.py` (277 行)
- **配置**: `config_prediction_c_bronze_matched_validation.yaml`
- **方法**: 以 lor2d 为参考族，对 KR_like / lor3d / lor4d 做 per-N greedy nearest-neighbor 配对 (协变量: antichain_width, comparable_fraction, geo_dim_eff, geo_interval_shape)
- **输出**: 396 对 (KR: 148, lor3d: 124, lor4d: 124), N=20–72

### 10.2 GPT 报告的表面结论

| 配对组 | hierarchy_integration_delta ↔ log_H_delta (Pearson r) | 置换 p |
|--------|------------------------------------------------------|--------|
| KR_like vs lor2d | -0.929 | 0.0002 |
| lor3d vs lor2d | -0.773 | 0.0002 |
| lor4d vs lor2d | -0.863 | 0.0002 |

### 10.3 深度诊断发现的问题

#### 问题 A: 匹配质量灾难性失败 (Critical)

| 配对组 | 匹配距离中位数 | antichain_width SMD | comparable_fraction SMD |
|--------|---------------|--------------------|-----------------------|
| KR vs lor2d | 3.87 | **+2.31** | **-3.96** |
| lor3d vs lor2d | 3.51 | +1.92 | **-4.92** |
| lor4d vs lor2d | 3.49 | +2.63 | **-8.59** |

> 标准匹配质量阈值: SMD < 0.25 为良好, < 0.5 为可接受。本验证 SMD = 2–9, 远超阈值。

**逐 N 的 antichain_width 差异** (KR vs lor2d):
N=20: +3.9, N=32: +7.6, N=44: +12.4, N=64: +20.1, N=72: +24.0 → 差异随 N 系统性增长，匹配本质上失败。

#### 问题 B: KR_like 结果完全由 lor2d 驱动 (Critical)

- KR_like 的 `layer_count` 在所有 N、所有 seed 下恒为 **3.0** (std = 0.000000)
- 因此 `layer_count_delta = 3.0 - lor2d_layer_count`
- delta 的变异 **100% 来自 lor2d 一侧**
- 报告的 r = -0.92 实质等价于 `r(lor2d_layer, lor2d_logH)` 的池化值，并非 KR 自身层级效应的证据

#### 问题 C: 二阶 Simpson's Paradox — 跨 N 池化膨胀 (Major)

**KR_like Δ(layer_count) vs Δ(log_H)**:

| 分析层级 | r |
|----------|---|
| 池化所有 N | **-0.916** |
| per-N 中位数 | **-0.477** |
| per-N 均值 | -0.479 |

**膨胀机制**: 随 N 增长，|Δ(layer)| 从 3.5 → 9.1，|Δ(logH)| 从 5.5 → 49.1，两个 delta 的绝对值同步放大，跨 N 池化时产生伪相关。

**Per-N 逐族 r(layer_delta, logH_delta)**:

| N | KR_like | lor3d | lor4d |
|---|---------|-------|-------|
| 20 | -0.45 | -0.61 | -0.54 |
| 24 | -0.56 | -0.61 | -0.61 |
| 28 | -0.77 | -0.92 | -0.82 |
| 32 | -0.45 | -0.62 | -0.19 |
| 36 | -0.65 | -0.56 | -0.05 |
| 40 | -0.50 | -0.50 | -0.20 |
| 56 | -0.10 | **+0.41** | **+0.25** |
| 64 | **+0.11** | +0.05 | -0.01 |

> 大 N 区间 (N≥56) 信号不稳定甚至反转，但这些区间每族仅 4–12 对。

#### 问题 D: 特征冗余与多重共线性 (Moderate)

1. **cover_density ≡ reduction_edge_density** — 所有族、所有样本完全恒等，两个名字指向同一个量
2. **adjacent_edge_fraction + long_edge_fraction = 1** (恒等) — delta 也精确互补为 0
3. **hierarchy_integration_delta_index** 是这些冗余特征的线性组合，该复合指标的高 r 值不可与真实自由度对应

#### 问题 E: long_edge_fraction 方向因族而异 (Moderate)

| 配对组 | r(long_edge_delta, logH_delta) |
|--------|-------------------------------|
| KR vs lor2d | **-0.81** |
| lor3d vs lor2d | **+0.31** |
| lor4d vs lor2d | **+0.44** |

> Prediction C 预期方向: 负。仅 KR 配对支持；lor3d、lor4d 方向相反。

### 10.4 什么结论可以成立

尽管存在上述方法学弱点，以下发现仍然成立：

1. **方向一致性**: `layer_count_delta` 和 `mean_layer_gap_delta` 在绝大多数 per-N 分析中为负相关，与 Prediction C 一致
2. **真实效应量**: within-N 中位 |r| ≈ 0.3–0.6，属于中等效应
3. **三族覆盖**: lor2d / lor3d / lor4d 在 N=20–52 区间信号方向稳健
4. **大 N 不确定**: N≥56 区间样本过少 (4–12 对)，信号不可靠

### 10.5 修正建议

| 问题 | 修正方案 |
|------|---------|
| 匹配质量 | 设置 caliper (如 d < 1.0)，拒绝低质量配对；或改用 propensity score matching |
| KR 退化 | 从 within-family 分析中排除 KR_like；仅保留跨族均值比较 |
| 跨 N 膨胀 | 逐 N 报告 per-N 相关性；用 Fisher z-transform 加权聚合 |
| 特征冗余 | 删除 cover_density (与 reduction_edge_density 相同)；删除 adjacent_edge_fraction (= 1 - long_edge_fraction) |
| 复合指标 | 分拆 hierarchy_integration_delta_index 为独立成分，分别报告 |
