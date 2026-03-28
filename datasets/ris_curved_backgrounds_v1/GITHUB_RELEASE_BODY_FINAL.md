# Dataset v1.0.0 — RIS Curved Backgrounds (N=512..2048)

Tag: `dataset-v1.0.0`  
Commit: `0069670`

---

## 中文

本次发布公开数据集 **RIS Curved Backgrounds Dataset v1**，用于验证 RIS 在曲率背景下的鲁棒性表现。

### 核心内容

- 网络规模：`N = 512, 1024, 2048`
- 背景配置：`6` 组（FLRW: κ=0.5/1.0/3.0；Schwarzschild: φ=0.01/0.03/0.05）
- 重复次数：每配置 `10` 次
- 总试验数：`180`
- 结果概览：本数据集中所有配置 `success_rate = 100%`

### 质量状态

- 数据集巡检：`19/19 PASS`
- 覆盖：文件结构、JSON 解析、关键语义一致性、CITATION 核心字段

### 资产完整性

- Archive: `ris_curved_backgrounds_v1.zip`
- SHA256: `7f910be71aa008c10c88dc8a2f8faf456e705950eef8473562b500f4a49736de`

### 发布资产（建议随 Release 上传）

- `ris_curved_backgrounds_v1.zip`
- `ris_curved_backgrounds_v1.zip.sha256.txt`
- `release_manifest.json`
- `RELEASE_ARTIFACT_SUMMARY.md`

### 数据集目录

`datasets/ris_curved_backgrounds_v1/`

包含：
- `README.md`
- `CITATION.cff`
- `data/large_scale_robustness_results.json`
- `data/cost_analysis.json`
- `reports/large_scale_summary.md`
- `reports/cost_analysis.md`
- `reports/quality_check_report.md`
- `reports/quality_check_report.json`
- `code/README.md`

### 复现脚本（仓库根目录）

- `experiment_curved_backgrounds_largeN.py`
- `experiment_curved_backgrounds.py`
- `curvature_backgrounds.py`

---

## English

This release publishes **RIS Curved Backgrounds Dataset v1**, an open dataset for robustness verification of RIS under curved backgrounds.

### Scope

- Network sizes: `N = 512, 1024, 2048`
- Background configurations: `6` total (FLRW: κ=0.5/1.0/3.0; Schwarzschild: φ=0.01/0.03/0.05)
- Repeats: `10` per configuration
- Total trials: `180`
- Headline result: `success_rate = 100%` for all configurations in this dataset

### Quality status

- Dataset QA: `19/19 PASS`
- Checks include: structure, JSON parsing, semantic consistency, citation metadata sanity

### Integrity

- Archive: `ris_curved_backgrounds_v1.zip`
- SHA256: `7f910be71aa008c10c88dc8a2f8faf456e705950eef8473562b500f4a49736de`

### Recommended release assets

- `ris_curved_backgrounds_v1.zip`
- `ris_curved_backgrounds_v1.zip.sha256.txt`
- `release_manifest.json`
- `RELEASE_ARTIFACT_SUMMARY.md`

### Dataset path

`datasets/ris_curved_backgrounds_v1/`

### Reproducibility scripts (repo root)

- `experiment_curved_backgrounds_largeN.py`
- `experiment_curved_backgrounds.py`
- `curvature_backgrounds.py`

---

## Citation

Please use `datasets/ris_curved_backgrounds_v1/CITATION.cff` for citation metadata.
