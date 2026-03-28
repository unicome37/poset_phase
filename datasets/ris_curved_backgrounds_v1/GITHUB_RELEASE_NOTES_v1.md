# Release: RIS Curved Backgrounds Dataset v1

## Highlights

- 首次公开 `RIS` 在曲率背景（FLRW / Schwarzschild）下的大规模鲁棒性数据集。
- 实验覆盖：`N = 512, 1024, 2048`；`6` 组背景配置；每配置 `10` 次重复；总计 `180` 次试验。
- 本数据集中：所有配置 `success_rate = 100%`（top node 命中 seed 集）。
- 质量巡检：`19/19` 检查通过（结构、JSON解析、关键语义、CITATION字段）。

## Assets

- `ris_curved_backgrounds_v1.zip`
- `ris_curved_backgrounds_v1.zip.sha256.txt`
- `release_manifest.json`
- `RELEASE_ARTIFACT_SUMMARY.md`

## Integrity

- SHA256 (`ris_curved_backgrounds_v1.zip`):

`7f910be71aa008c10c88dc8a2f8faf456e705950eef8473562b500f4a49736de`

## Dataset Structure

- `README.md`
- `CITATION.cff`
- `data/large_scale_robustness_results.json`
- `data/cost_analysis.json`
- `reports/large_scale_summary.md`
- `reports/cost_analysis.md`
- `reports/quality_check_report.md`
- `reports/quality_check_report.json`
- `code/README.md`

## Reproducibility

核心脚本位于仓库根目录：

- `experiment_curved_backgrounds_largeN.py`
- `experiment_curved_backgrounds.py`
- `curvature_backgrounds.py`

## Notes

- 该数据集目录位于：`datasets/ris_curved_backgrounds_v1/`
- 建议发布标签：`dataset-v1.0.0`
- 建议 release 标题：`Dataset v1.0.0 — RIS Curved Backgrounds (N=512..2048)`
