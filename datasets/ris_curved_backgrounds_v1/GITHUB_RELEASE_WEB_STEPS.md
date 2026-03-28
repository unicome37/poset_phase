# GitHub Release 网页填写步骤（dataset-v1.0.1）

适用仓库：`unicome37/poset_phase`  
推荐标签：`dataset-v1.0.1`（已包含正式正文文件）

---

## 1) 进入 Release 页面

1. 打开仓库主页：`https://github.com/unicome37/poset_phase`
2. 点击右侧 **Releases**
3. 点击 **Draft a new release**

---

## 2) 基础字段填写

### Tag version

填写：`dataset-v1.0.1`

> 若页面未自动识别，可手工输入并选择已存在 tag。

### Release title

建议：`Dataset v1.0.1 — RIS Curved Backgrounds (N=512..2048)`

### Target

保持 `main`（默认即可）

---

## 3) Release 正文

将下列文件内容完整粘贴到正文框：

- `datasets/ris_curved_backgrounds_v1/GITHUB_RELEASE_BODY_FINAL.md`

> 该正文已含：中英双语说明、数据范围、QA 状态、SHA256、复现入口。

---

## 4) 上传发布资产（Attach binaries）

从本地目录上传以下 4 个文件：

目录：`datasets/ris_curved_backgrounds_v1/dist/`

- `ris_curved_backgrounds_v1.zip`
- `ris_curved_backgrounds_v1.zip.sha256.txt`
- `release_manifest.json`
- `RELEASE_ARTIFACT_SUMMARY.md`

---

## 5) 发布前检查清单

- [ ] Tag = `dataset-v1.0.1`
- [ ] 标题无误
- [ ] 正文来自 `GITHUB_RELEASE_BODY_FINAL.md`
- [ ] 4 个资产全部上传完成
- [ ] SHA256 与正文一致：
  `7f910be71aa008c10c88dc8a2f8faf456e705950eef8473562b500f4a49736de`

---

## 6) 发布动作

点击 **Publish release**。

发布后建议立刻做两件事：

1. 在仓库 `README` 或项目总索引中加入该 Release 链接。  
2. 用 `ZENODO_METADATA_DATASET_V1.json` 作为模板完成 Zenodo 记录创建。

---

## 可选：若坚持使用 dataset-v1.0.0

也可发布 `dataset-v1.0.0`，但正文文件不在该标签快照中；需要从 `main` 分支手动拷贝正文。  
因此仍建议优先发布 `dataset-v1.0.1`。
