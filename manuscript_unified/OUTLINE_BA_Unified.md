# B+A 合并论文框架设计

## 投稿目标

**期刊**: Classical and Quantum Gravity (IOP Publishing)
- CST 论文天然家园（Surya, Loomis, Glaser, Cunningham 均发此刊）
- 接受数值方法 + 新框架
- 影响因子 ~3.6, 同行评审
- 格式：LaTeX (iopart class), 无字数限制

**备选**: Foundations of Physics (Springer) — 已有 manuscript_foundphys.tex 模板

---

## 标题

**主标题（推荐）**:
> Geometric Phase Transitions and Emergent Dimensional Selection in Finite Causal Poset Ensembles

**备选标题**:
> From Entropy Catastrophe to Dimensional Emergence: Bounded Phase Transitions and Consistency-Driven Selection in Finite Causal Posets

---

## 合并叙事逻辑

### 核心故事线（一条线贯穿）

1. **问题提出**: CST 的熵灾难 — KR 定理意味着泛型 poset 是三层的、最大熵的，几何结构如何涌现？
2. **框架构建**: 7族+4族 集成集合，精确线性扩展计数，可分解作用量
3. **第一个结果**: 有界相变 — γ_c 在 N=10-44 保持 O(1)（原 Prediction B 核心）
4. **消融揭示机制**: 7个几何子项中，只有两个构成骨架（g_wh + g_dim）
5. **关键转折**: g_dim 可被非锚定 g_con 替换 → 相变保持
6. **自然预测**: 如果作用量奖励一致性而非特定维度，高维应该胜出
7. **第二个结果**: 4D 无条件支配 — 92/98 配置，线性增长的优势边际（原 Prediction A 核心）
8. **鲁棒性确认**: 种子敏感性 100% 胜率, 权重跨4倍稳定

### 为什么合并比分拆更好？

- 消除"companion paper"交叉引用问题（两篇互引未发表稿是审稿大忌）
- 叙事链完整：问题→框架→相变→机制→预测→确认 = 一个自洽的科学故事
- CQG 接受长文（20-30页正常），不需要压缩

---

## 论文结构

### Abstract (~250 words)
统一叙事："我们解决两个问题——几何相变的存在性和维度的选择性"

### 1. Introduction (~2 pages)
- 1.1 熵灾难与因果集理论
- 1.2 已有工作：BD action, MCMC, 有限尺度标度
- 1.3 两个开放问题：(a) 相变阈值是否有界？(b) 维度是先验的还是涌现的？
- 1.4 本文方法：精确熵 + 可分解作用量 + 系统消融

### 2. Framework (~3 pages)
- 2.1 Poset ensemble (8族: lor2d/3d/4d, kr, mlr, tp, int, abs)
- 2.2 精确线性扩展计数与 SIS 交叉验证
- 2.3 可分解作用量定义（A₁, A₂, A₃）
- 2.4 几何惩罚子项：7 个细目 + 归一化
- 2.5 临界耦合 γ_c 定义

### 3. Bounded Phase Transition (原 Pred B) (~3 pages)
- 3.1 γ_c(N) 曲线：N=10-44，全部精确
- 3.2 对照实验：A₁（纯中性）无交叉
- 3.3 时序探针：lor2d 精确计算可达 N=104

### 4. Ablation and Minimal Backbone (~2 pages)
- 4.1 单项移除消融热图
- 4.2 骨架-外壳结构：g_wh + g_dim 为骨架
- 4.3 非循环替换：g_dim → g_con
- 4.4 最小非循环骨架 (g_wh + g_con)

### 5. Dimensional Selection (原 Pred A) (~4 pages)
- 5.1 扩展至 4D 家族
- 5.2 原始作用量下的 36:35 平手
- 5.3 消融定位：g_dim 是唯一抑制源
- 5.4 一致性替换：92/98 4D 支配
- 5.5 优势边际的线性增长 (+7 → +57)
- 5.6 种子鲁棒性 (42/42, 100%)
- 5.7 一致性权重敏感性 (w_con ∈ [0.3, 1.2])

### 6. Discussion (~2 pages)
- 6.1 与先前工作的关系（Carlip, Surya, Loomis, Glaser, Cunningham）
- 6.2 与 BD action 的关系（明确区分 + 未来联结方向）
- 6.3 局限性（有限 N, 有限家族, SIS 估计, 生成器依赖）
- 6.4 展望：5D/6D 扩展, 更快精确算法, 场论约束

### 7. Conclusion (~0.5 page)

### Appendices
- A. 精确线性扩展算法细节与计时基准
- B. g_con 权重校准方法
- C. SIS 交叉验证表

### References (~25-30 refs)

---

## 从现有材料的映射

| 合并论文章节 | B 论文来源 | A 论文来源 |
|---|---|---|
| §1 Introduction | §1 全部 | §1 2nd-order question 部分 |
| §2 Framework | §2 全部 | §2.1 4D 族定义 |
| §3 Bounded Phase | §3.1 全部 | — |
| §4 Ablation | §3.2-3.3 全部 | §3 消融表的扩展版 |
| §5 Dim Selection | — | §3-5 全部 |
| §6 Discussion | §4 全部 | §5 全部（合并去重） |

---

## Cover Letter 策略

**关键要素**：
1. 明确定位为 CST 数值方法论文（不是哲学性的）
2. 强调与 CQG 已发表工作的对话（Surya 2012, Loomis 2018, Glaser 2018, Cunningham 2020）
3. 突出"精确熵"作为方法学贡献（vs MCMC）
4. 不提"companion paper" — 这是一篇自包含论文
5. AI 辅助声明简洁化（一句即可）

**开头示例**：
> "We present the first finite-size computation of geometric phase transition thresholds in causal poset ensembles using exact linear extension counts, and demonstrate that dimensional consistency constraints — without dimension-specific priors — naturally select higher-dimensional Lorentzian structures."

---

## Zenodo 上传规划

**上传物**：
1. `manuscript_unified_CQG.tex` — 合并版 LaTeX
2. `manuscript_unified_CQG.pdf` — 编译 PDF
3. Figures 目录
4. Supplementary data tables

**元数据**：
- Title: (同论文标题)
- Type: Preprint
- License: CC BY 4.0
- Related identifiers: 链接到代码仓库 DOI

---

## 时间线

- Day 1-2: 创建 CQG LaTeX 模板，合并 Framework 和 Bounded Phase 部分
- Day 3-4: 合并 Ablation 和 Dimensional Selection 部分
- Day 5: 撰写统一 Introduction、Discussion、Abstract
- Day 6: 编译检查, 图表集成
- Day 7: Zenodo 上传 + Cover letter
- Day 8: 提交 CQG
