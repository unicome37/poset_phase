# Code Index for Reproducibility

本目录只存索引，实际代码位于仓库根目录。

## 核心脚本

- `../../experiment_curved_backgrounds_largeN.py`
  - 大规模曲率背景实验（含 N>512 扩展）
- `../../experiment_curved_backgrounds.py`
  - 基础曲率背景实验
- `../../curvature_backgrounds.py`
  - FLRW / Schwarzschild 背景实现

## 运行建议

1. 优先在虚拟环境运行。
2. 先执行小规模 smoke，再跑大规模。
3. 每个分区完成后及时落盘 JSON，避免长任务中断损失。
