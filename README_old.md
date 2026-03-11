# poset_phase

最小可执行的偏序相图实验骨架。

当前版本只覆盖文档《结构存在论_更大规模稳健数值实验方案》中的最小子集：

- 7 类偏序结构生成器
- 中立惩罚 `I_neutral`
- 几何惩罚代理 `I_geometric`
- 小 `N` 精确线性延拓数
- 中大 `N` 的 SIS 近似熵估计
- `A_1 / A_2 / A_3` 三套作用量扫描
- 基于 poset 样本的外层 Bootstrap 汇总
- 跨 `N` 的基础归一化列生成
- 基于 SIS 误差估计的两层 Bootstrap 置信区间

尚未实现：

- `A_4` 的 BD 作用量附加对照
- `A_5` 的自动编码器探索分析
- 更丰富的标度敏感性分析模块
- 更严格的几何惩罚定义与物理论证

## 文件说明

- `generators.py`: 结构族生成器与 `Poset` 数据结构
- `observables.py`: 可比率、度统计、层分布、中立惩罚
- `observables_geo.py`: 几何惩罚代理
- `entropy_exact.py`: `N<=20` 可用的精确线性延拓数
- `entropy_sis.py`: SIS 近似估计 `log N_prec`
- `action.py`: 最小作用量接口
- `experiment.py`: 参数扫描与结果汇总
- `calibration.py`: 小 `N` 精确熵对 SIS 的系统校准
- `plots.py`: 从 CSV 结果生成相图与置信区间图
- `normalization.py`: 跨 `N` 归一化
- `bootstrap.py`: 以 poset 为单位的 Bootstrap，可选传播 SIS 误差
- `config.yaml`: 实验参数配置

## 结果分层

为避免把“冻结主线结果”和“后续探索性试验”混写，当前目录下的输出建议按两类理解：

- `outputs_confirmatory/`: 冻结主线、可作为当前正文主证据链引用的结果
- `outputs_exploratory/`: 诊断性、扩展性、赛后试探性结果

对应说明见：

- `RESULTS_INDEX.md`

## 运行方式

```bash
python experiment.py
```

也可以显式指定配置：

```bash
python experiment.py --config config.yaml
```

快速烟雾测试：

```bash
python experiment.py --config config_smoke.yaml
```

做 SIS 校准：

```bash
python calibration.py --config config_smoke.yaml
```

生成图：

```bash
python plots.py --config config_smoke.yaml
```

默认会：

1. 读取 `config.yaml`
2. 生成原始样本表
3. 追加归一化列
4. 生成 `score_norm_std_est`
5. 生成 Bootstrap 摘要
6. 写入 `outputs/`

如果再运行 `calibration.py` 和 `plots.py`，还会额外产出：

- `calibration_raw.csv`
- `calibration_summary.csv`
- `plots/*.png`

## 注意

当前结果只应视为实验框架检查，不应直接拿来支持最终理论结论。严肃比较前还需要：

1. 在小 `N` 上完成 SIS 对精确熵的系统校准。
2. 对当前几何惩罚代理做更严格的物理论证或替换。
3. 对 `H` 和 `I` 的标度归一化方法做敏感性分析，而不只依赖单一方案。
4. 对当前“高斯近似的 SIS 误差传播”做更严格的层级化替代。

补充：

- 当前最适合作为确认性结果引用的是 `outputs_frozen_exact` 与 `outputs_frozen_cg`
- 当前新增的确认性精确扩展结果位于 `outputs_confirmatory/medium_exact`
- 维度扫描、窗口约束、生成器重标定等结果默认应归入探索性分析
