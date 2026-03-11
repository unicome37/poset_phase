# 获取 Zenodo DOI 操作指引

## 方式一：GitHub Release → Zenodo 自动归档（推荐）

### 第一步：关联 Zenodo 与 GitHub

1. 访问 https://zenodo.org/login/ ，用 GitHub 账号登录
2. 进入 https://zenodo.org/account/settings/github/
3. 在仓库列表中找到 `unicome37/poset_phase`
4. 点击开关启用（Toggle ON）

### 第二步：在 GitHub 创建 Release

1. 访问 https://github.com/unicome37/poset_phase/releases/new
2. 填写以下信息：
   - **Tag**: 选择 `v1.0.1`（如果不存在会自动创建）
   - **Release title**: `v1.0.1 — Action-Weighted Poset Ensemble Analysis`
   - **Description**: 复制 `RELEASE_NOTES_v1.0.1.md` 的内容
   - **Attach binaries**: 上传 `zenodo_release/poset_phase_v1.0.zip`（可选，Zenodo 会自动获取源码包）
3. 点击 **Publish release**

### 第三步：获取 DOI

1. Zenodo 会在几分钟内自动归档
2. 访问 https://zenodo.org/account/settings/github/ 查看归档状态
3. 点击归档条目获取 DOI（格式如 `10.5281/zenodo.XXXXXXX`）
4. DOI badge 会显示在仓库中（如果添加了 badge）

---

## 方式二：直接上传到 Zenodo

如果不想关联 GitHub，可以手动上传：

1. 访问 https://zenodo.org/deposit/new
2. 上传 `zenodo_release/poset_phase_v1.0.zip`
3. 填写元数据（已在 `.zenodo.json` 中定义）：
   - **Title**: Poset Phase: Action-Weighted Partition Functions on Finite Poset Ensembles
   - **Upload type**: Software
   - **Description**: [从 RELEASE_NOTES 复制]
   - **License**: MIT
   - **Keywords**: causal set theory, phase transition, linear extensions, Kleitman-Rothschild, quantum gravity
4. 点击 **Publish**
5. 获得 DOI

---

## 获取 DOI 后的更新清单

拿到 DOI 后需要更新以下文件中的占位符：

1. `MANUSCRIPT_DRAFT.md` — 搜索 `[Zenodo DOI]` 替换为实际 DOI
2. `manuscript.tex` — 搜索 `[Zenodo DOI]` 替换为实际 DOI  
3. `manuscript_foundphys.tex` — 同上
4. `mdpi_template/entropy_manuscript.tex` — 同上
5. `README.md` — 添加 DOI badge:
   ```markdown
   [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
   ```
