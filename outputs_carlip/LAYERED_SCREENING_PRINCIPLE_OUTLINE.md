# 分层结构筛选原理：从准入、认同到历史沉积

**核心立论**：结构不是被单一机制选择出来的，而是经由准入和认同**分层确认**，并在历史沉积中**纵向锁定**的。

**地位**：这是 causal set 数值研究中提取的中层有效理论（mid-level effective principle），介于底层数值实验与顶层"存在即稳定"哲学母句之间。

---

## 提纲

### 一、引言：为什么需要分层

- Kleitman-Rothschild 熵灾难：generic poset ≠ 几何
- BD action 的局限：rank 14/17，线性准入但不足以选身份
- 关键问题的重新提法：
  - 不是"什么 action 最好？"
  - 而是"结构选择的机制本身有几层？每层做什么？"

### 二、第一层——准入（Admissibility）

#### 2.1 定义
- 准入层回答：这个 poset **能不能活**？
- 低阶、宽容、大范围淘汰

#### 2.2 准入层内部的粗筛与几何筛分工

在准入层内部，存在从一般到几何的筛选梯度：

- **$S_{\mathrm{triple}}$**（三重积分）：刻画一般存在资格——可显影、可延续的粗筛
- **$S_{\mathrm{BD}} = \mathbf{c}^\top \Delta\mathbf{C}$**（Benincasa-Dowker 线性超平面）：更接近 Lorentzian 几何的曲率型准入门

两者均属准入级：排除明显非几何结构（KR-type、random DAG），但都无法从几何簇中正选 Lor4D。

#### 2.3 关键特征
- 必要不充分
- $S_{\mathrm{BD}}$ 下 Lor4D 仅排 14/17 original families——不是认同器，是看门人
- $S_{\mathrm{triple}}$ 下 Lor4D 仅排 rank #12——证明三重筛选=准入门非身份选择器
- 对 N 不太敏感（在所有测试 N 下均有效）

### 三、第二层——认同（Identity）

#### 3.1 定义
- 认同层回答：这个 poset **是不是 Lor4D**？
- 高阶、精确、在当前 25 族库中稳健执行身份锁定

#### 3.2 数学载体
- $S_{\mathrm{MD}}[P,N] = \delta^\top \Sigma^{-1}(N)\,\delta$（二次椭球）
- 零自由参数：$(\boldsymbol{\mu}(N), \Sigma(N))$ 完全由 Lor4D 集成统计量决定
- $(\boldsymbol{\mu}(N), \Sigma(N))$ 不仅是统计工具，更是 Lor4D 吸引域的 reference object（详见§六）

#### 3.3 物理含义
- 身份认同不是大 N 极限的远方产物，而是在有限尺度 $N \approx 14$ 就已可操作地存在
- $N \geq 14$ 时唯一 runner-up 始终是 Lor5D——turn-on boundary 不由判别器能力设定，而由 **4D/5D 因果结构的内禀几何邻近性** 决定

### 四、认同层的尺度动力学：Turn-on 与 Basin Deepening

#### 4.1 三阶段图景

**Phase I — 预开启**（$N \lesssim N_{\mathrm{res}} \approx 12$）
- 信息量不足，4D/5D 因果几何无法分辨
- 参考流形尚未"显影"到支撑可靠认同
- N=12 时 Lor5D 和 KR_2layer 均可作为 intruder（$d_{\mathrm{eff}}$ 重叠）

**Phase II — 开启**（$N \approx N_{\mathrm{id}} \approx 14$）
- 10/10 seeds rank #1，margin 从 +0.10 跳到 +1.28
- 身份判别器从"偶尔有效"突变为"系统可靠"
- **这是一个相变式 turn-on**

**Phase III — 加深**（$N \gg N_{\mathrm{id}}$）
- Margin 持续增长：$+1.28 \to +4.59$（N=14→28），非严格单调但趋势明确
- Gap 增长至 $1.93 \times 10^8$（N=1024）
- $V_{\mathrm{eff}} \propto N^{-1.57}$，Fisher $\propto N$
- 身份井越来越深、越来越窄——偏离代价无界增长

#### 4.2 Basin deepening 的定量刻画

$$\Delta_{\mathrm{hist}}(N) \equiv \min_{f \neq \mathrm{Lor4D}} \bigl[S_{\mathrm{MD}}(f, N) - S_{\mathrm{MD}}(\mathrm{Lor4D}, N)\bigr]$$

即 Lor4D 与最近竞争者之间的 Mahalanobis 隔离间隙。

| $N$ | $\Delta_{\mathrm{hist}}$ | 状态 |
|---|---|---|
| 12 | $\approx 0.10$ | 沉积不足，身份未锁 |
| 14 | $\approx 1.28$ | turn-on |
| 16 | $\approx 1.72$ | 稳定加深 |
| 20 | $\approx 2.05$ | 稳定加深 |
| 28 | $\approx 4.59$ | 深度加深 |
| 1024 | $\to 10^8$ | 深度锁定 |

### 五、纵向维度——历史沉积（Historical Sedimentation）

§四建立了 basin deepening 的数值事实。本节将其提升为中层理论解释：**身份井加深的过程就是"历史沉积"的数学形式。**

#### 5.1 从哲学到数学

原有表述："过去保留在结构中，参与塑形未来。"

新的 basin 解释：**历史不是额外添加的记忆，而是参考流形逐渐显影、协方差逐渐收缩、偏离代价逐渐增大的 basin 过程。**

| 哲学概念 | 数学对应 | 可测量 |
|---|---|---|
| 历史逐渐沉积 | 参考流形逐渐显影 | $\boldsymbol{\mu}(N)$ 轨迹收敛 |
| 过去改写未来 | 协方差逐渐收缩 | $\det(\Sigma) \propto N^{-3.38}$ |
| 惯性增强 | 身份井加深 | gap 增长、$V_{\mathrm{eff}}$ 收缩 |
| 偏离代价增大 | off-manifold cost 增长 | Fisher $\propto N$ |
| 方向逐渐锁定 | 可达方向减少 | 椭球截面 shrink |

#### 5.2 核心命题

> **Proposition (历史沉积的单向累积)**：
> 在有限因果结构空间中，一旦 Lor4D 身份井在 $N \geq N_{\mathrm{id}}$ 开启，
> 随后随着 $N$ 的增大，隔离间隙 $\Delta_{\mathrm{hist}}(N)$ 呈总体增长趋势、
> 身份井体积 $V_{\mathrm{eff}}(N)$ 呈总体收缩趋势——
> 即历史沉积在统计意义上表现为单向累积。

这不是定理（目前），但它是数值强支持的、可以被未来的解析推导证明或证伪的**中层命题**。

#### 5.3 区分：历史沉积不是独立的第三层筛选

历史沉积描述的是**认同层在 $N$ 增大下的纵向锁定行为**，不是与准入/认同并列的第三个独立筛选算子。它依赖的对象——$\boldsymbol{\mu}(N)$、$\Sigma(N)$、gap、$V_{\mathrm{eff}}$——全属认同层的尺度演化量。

正确的层级关系：
- **横截面**（fixed $N$）：准入层 + 认同层 → 分层筛选
- **纵向**（across $N$）：认同层的 basin 随 $N$ 持续加深 → 历史沉积

### 六、参考流形作为理论对象

#### 6.1 地位升格
- 不再是"Mahalanobis 的工具组件"
- 而是 Lor4D 在有限尺度特征空间中的**吸引域几何本身**

#### 6.2 形式定义

**Lor4D 参考流形** $\mathcal{M}_4$ 是一条带收缩横截面的参数化参考轨迹：

$$\mathcal{M}_4 = \bigl\{(\boldsymbol{\mu}(N),\, \Sigma(N)) \;\big|\; N \geq N_{\min}\bigr\}$$

其中：
- $\boldsymbol{\mu}(N) = \boldsymbol{\mu}(\infty) + \mathbf{a}/N + \mathbf{b}/N^2$（中心轨迹）
- $\Sigma(N) = \mathrm{diag}(A_i \cdot N^{-p_i})$（局域几何 / 横截面）
- 极限 $\boldsymbol{\mu}(\infty)$ 各分量有独立的理论渊源（$d_{\mathrm{eff}} \to 4$ by Myrheim-Meyer）
- $N_{\min}$ 是统计量可靠估计的最小尺度，不与 $N_{\mathrm{id}}$ 绑定

> **术语说明**：这里的"流形"是中层有效理论意义下的参数化参考对象，不预设严格微分流形结构。严格数学化——若有需要——可以在未来的解析推导中完成。

#### 6.3 三个子结构

| 子结构 | 数学 | 物理 |
|---|---|---|
| 中心轨迹 | $\boldsymbol{\mu}(N) \to \boldsymbol{\mu}(\infty)$ | Lor4D 结构中心在有限 N 下的偏移 |
| 横截面 | $\Sigma(N) \to 0$ | 集成涨落收缩，身份锐化 |
| 身份井 | $S_{\mathrm{MD}} = \delta^\top \Sigma^{-1} \delta$ | 距参考流形的二次代价 |

### 七、分层结构筛选原理

本节汇总§二—§六的结果，给出全文中心命题。

#### 7.1 正式表述

> **分层结构筛选原理（Layered Structural Screening Principle）**：
>
> 在具有足够多因果事件的有限 poset 集合中，4D Lorentzian 因果结构的选择
> 至少需要两个不同阶次的独立机制：
>
> 1. 一阶准入层（由 $S_{\mathrm{triple}}$ 与 $S_{\mathrm{BD}}$ 共同承载），排除非几何与非稳健结构；
> 2. 二阶二次身份函数 $S_{\mathrm{MD}}$，将 Lor4D 从所有准入者中正选出来。
>
> 两层在功能上分离，在统计上近弱相关（$r \to 0$），其交集在当前 25 族库中唯一稳健选出 Lor4D。
>
> 此外，认同层的身份井在 $N$ 增大时表现出单向累积的 basin deepening（历史沉积），使 Lor4D 在当前家族库中成为越来越不可替代的唯一稳健身份中心。

#### 7.2 数值证据强度

| 证据类型 | 现状 |
|---|---|
| 25 族 × N=14-1024 × 10 seeds | ✅ 完成 |
| Turn-on boundary $N_{\mathrm{id}} \approx 14$ | ✅ 定量 |
| $S_{\mathrm{BD}}$ 与 $S_{\mathrm{MD}}$ 功能分离 + 弱相关 | ✅ Pearson $r \to 0$ |
| Basin deepening 幂律 | ✅ $V_{\mathrm{eff}} \propto N^{-1.57}$ |
| 历史沉积单向性 | ✅ $\Delta_{\mathrm{hist}}$ 递增趋势 |
| 解析推导 | ❌ 缺失 |
| $N > 1024$ 验证 | ❌ 缺失 |

#### 7.3 与已有理论的关系

| 已有概念 | 本框架中的定位 |
|---|---|
| BD action | 第一层载体 |
| KR 熵灾 | 第一层存在的动机 |
| Mahalanobis 距离 | 第二层的统计形式（非分类器） |
| CLT | $\sigma^2 \propto N^{-1}$ 的理论来源 |
| Myrheim-Meyer | $d_{\mathrm{eff}} \to 4$ 的解析保证 |
| 历史沉积论 | 纵向维度（basin deepening）的哲学对应 |

### 八、讨论

#### 8.1 这不是机器学习
- 零自由参数
- 参考流形由理论极限（Myrheim-Meyer, CLT）锚定，不是从数据"学"到
- $S_{\mathrm{MD}}$ 是**单类判定**（one-class）——只依赖 Lor4D 自身统计量，无训练/测试拆分、无决策边界优化
- 分层不是超参数选择，而是阶次递增的物理机制

#### 8.2 升阶而非替代
- $S_{\mathrm{MD}}$ 不替代 $S_{\mathrm{BD}}$，而是升阶
- 从一阶线性 → 二阶二次 是自然的分析力学扩展
- 暗示：可能存在更高阶的层次——信息几何、Rényi 距离？

#### 8.3 与"存在即稳定"的关系
- 母句："存在即筛选后的残留"
- 本文给出的子句："Lor4D 是准入层与认同层交集中的唯一稳健胜出者"
- 母句的纵向子句——历史沉积——现在有了可量化的 $\Delta_{\mathrm{hist}}(N)$

#### 8.4 边界与未完成项

本文主动声明以下边界：

- **这是中层有效理论**，不是终极 action、不是基本原理推导
- **这是在当前 25 族扩展库下的稳健结论**，不是全结构空间的穷尽证明
- **历史沉积命题是数值强支持**，不是解析定理；严格单调性尚未在所有 N 下排除反例（N=28→32 有波动）
- **Reference manifold 是有效对象**，不预设严格微分流形结构
- **Gradient bridge（$\cos(\nabla S_{\mathrm{BD}}, \nabla F_{\mathrm{LSD}}) = 0.97$）只是辅助支持**，不承担主论证

### 九、下一步

1. **解析推导**：从 Myrheim-Meyer 公式和 CLT 导出 $\Delta_{\mathrm{hist}}(N)$ 的增长律
2. **更高维**：$N_{\mathrm{id}}$ 如何随目标维度 $d$ 变化？
3. **可能的更高阶层**：beyond 二次——信息几何（Fisher-Rao）、Rényi divergence？
4. **路径积分连接**：$e^{-S_{\mathrm{MD}}}$ 能否作为路径积分中的权重因子？

---

## 关键术语表

| 术语 | 定义 | 首次出现 |
|---|---|---|
| **准入层** | 一阶线性筛选（含粗筛 $S_{\mathrm{triple}}$ 与几何筛 $S_{\mathrm{BD}}$），排除非几何结构 | §二 |
| **认同层** | 二阶二次筛选 $S_{\mathrm{MD}}$，在当前族库中稳健正选 Lor4D | §三 |
| **参考流形** $\mathcal{M}_4$ | $\{(\boldsymbol{\mu}(N), \Sigma(N))\}$，中层有效理论意义下的参数化参考对象 | §六 |
| **Turn-on scale** $N_{\mathrm{id}}$ | 认同层可靠开启的最小 N（$\approx 14$） | §四.1 |
| **Resolution floor** $N_{\mathrm{res}}$ | 信息量不足以支撑认同的尺度（$\approx 12$） | §四.1 |
| **隔离间隙** $\Delta_{\mathrm{hist}}(N)$ | Lor4D 与最近竞争者的 Mahalanobis score 差 | §四.2 |
| **历史沉积** | 认同层 basin 随 N 增大而持续加深的纵向锁定行为；不是独立的第三层筛选 | §五 |
| **分层结构筛选原理** | 选择至少需要两个不同阶次的独立机制，其功能分离且弱相关 | §七.1 |

---

*提纲版本 v2.0 — 2026-03-27*
*基于 poset_phase 项目 N-boundary turn-on 实验结果及 Carlip 后两层筛选架构*
*v1.0→v2.0 修订：历史沉积从"第三层"改为"纵向维度"；核心句/原理表述/术语精度全面修正；增加§8.4边界声明*
