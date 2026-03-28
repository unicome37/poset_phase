# RIS Curved Backgrounds Dataset v1

公开数据集：`RIS` 在曲率背景（FLRW / Schwarzschild）下的大规模鲁棒性验证。

## 数据范围

- 网络规模：`N = 512, 1024, 2048`
- 背景配置：`6` 组（FLRW: κ=0.5/1.0/3.0；Schwarzschild: φ=0.01/0.03/0.05）
- 重复次数：每配置 `10` 次
- 总试验：`180`

## 主要结论

- 本批公开数据中，所有配置成功率均为 `100%`（top node 命中 seeds）。
- 成本侧可行：`N=2048` 约 `41.6MB` 内存级别（按当前估算模型）。

## 目录结构

- `data/large_scale_robustness_results.json`：主结果数据
- `data/cost_analysis.json`：成本估算
- `reports/large_scale_summary.md`：结果摘要
- `reports/cost_analysis.md`：成本摘要
- `code/README.md`：复现实验脚本索引
- `CITATION.cff`：引用信息

## 复现实验

请直接使用仓库根目录脚本：

- `experiment_curved_backgrounds_largeN.py`
- `experiment_curved_backgrounds.py`
- `curvature_backgrounds.py`

详见 `code/README.md`。

## 许可证

本数据集沿用仓库许可证：见根目录 `LICENSE`。
