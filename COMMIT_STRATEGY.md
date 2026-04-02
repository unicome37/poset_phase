# Atomic Commit Strategy (config / outputs / paper)

本仓库后续提交采用“三类原子提交”策略，避免一次提交混合多种资产造成审阅困难。

## 1) config 类（实验配置/脚本参数）

**范围**
- `config_*.yaml`
- 与配置强关联的轻量脚本改动（仅参数读取/调度，不含结果叙事）

**提交信息模板**
- `feat(config): add <experiment-scope> configs`
- `refactor(config): normalize <scope> parameter layout`

**提交前检查**
- 配置文件命名与实验目标一致
- smoke 配置与 full 配置成对存在（若适用）
- 不混入结果文件（`outputs*/`）

---

## 2) outputs 类（实验产物/审计报告）

**范围**
- `outputs*/**/*.json`
- `outputs*/**/*.md`
- 阶段性执行摘要（如 `phaseB_*`）

**提交信息模板**
- `chore(outputs): sync <experiment-scope> results`
- `docs(outputs): add <scope> audit notes`

**提交前检查**
- 结果文件带清晰时间戳或 scope 标识
- 不混入配置结构变更（`config_*.yaml`）
- 不混入论文正文改写（`*.tex`）

---

## 3) paper 类（论文与投稿资产）

**范围**
- `two_layer_paper/**/*.tex`
- `two_layer_paper/**/*.md`
- `two_layer_paper/figures/*`
- 其他投稿目录（cover letter / submission package）

**提交信息模板**
- `docs(paper): update <paper-name> draft`
- `docs(paper): refresh figures for <section/claim>`

**提交前检查**
- 图文编号一致（figure/table/citation）
- 参考文献与正文口径一致
- 不混入运行配置与原始实验输出

---

## 推荐节奏

1. 先提交 `config`
2. 再提交 `outputs`
3. 最后提交 `paper`

如必须跨类改动，拆成多提交并在提交信息中显式说明依赖顺序。

## 快速分批（建议流程）

- 第一步：仅暂存 config 相关文件并提交
- 第二步：仅暂存 outputs 相关文件并提交
- 第三步：仅暂存 paper 相关文件并提交

保持每次提交可单独审阅、可回滚、可复现实验意图。

## 可复制的 commit message 清单

### config
- `feat(config): add <experiment-scope> configs`
- `fix(config): correct <scope> parameter bounds`
- `refactor(config): normalize <scope> parameter layout`

### outputs
- `chore(outputs): sync <experiment-scope> results`
- `docs(outputs): add <scope> audit notes`
- `docs(outputs): update <scope> stage report`

### paper
- `docs(paper): update <paper-name> draft`
- `docs(paper): refresh figures for <section/claim>`
- `docs(paper): align references and captions for <paper-name>`

### repo / workflow
- `chore(repo): enforce line endings and atomic commit templates`
- `docs(repo): add contributing guide referencing commit strategy`

> 建议：同一次提交只选一个 scope；如需跨 scope，拆成多提交并按依赖顺序推送。
