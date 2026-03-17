# Occupancy Derivation Note for Prediction A

更新时间：2026-03-17

## 目标

把 `Prediction A` 中 entropy correction 的几何来源，从经验闭合推进到可写入论文的解析近似：

$$
\frac{\log H}{N}
\approx
\beta_0 (1-p_d)(\log N-1)
\;+\;
\beta_1\,P_{\mathrm{occ}}^{(1)},
$$

其中

$$
P_{\mathrm{occ}}^{(1)}
\equiv
1-\exp\!\bigl(-(N-2)\mathbb{E}[V_A]\bigr).
$$

这里 `V_A` 是两相关点之间 Alexandrov interval 的体积。

---

## 1. 几何起点

对 cube sprinkle 中一对相关点，设其 proper time 为 `tau`，则平直 Minkowski 空间中的 Alexandrov interval 体积为

$$
V_A = \kappa_d \tau^d,
\qquad
\kappa_d = \frac{\mathrm{Vol}(B_{d-1})}{2^{d-1}d}.
$$

在当前生成器里，`tau` 的分布由 cube 差分分布与闵氏因果约束共同决定。

---

## 2. Link 概率与空区间

若在该 interval 中没有第三个点，则该相关对是 link。对给定 interval 体积 `V_A`，其 link 概率近似为

$$
P(\mathrm{link}\mid V_A)
\approx
(1-V_A)^{N-2}.
$$

当 `V_A` 不太大时，

$$
(1-V_A)^{N-2} \approx e^{-(N-2)V_A}.
$$

因此非 link 概率近似为

$$
P(\mathrm{nonlink}\mid V_A)
\approx
1-e^{-(N-2)V_A}.
$$

对相关点对取平均：

$$
1-\ell_d
\approx
1-\mathbb{E}\!\left[e^{-(N-2)V_A}\right].
$$

这一步说明：`1-\ell_d` 的自然解释不是“区间数本身”，而是“至少出现一个中介点”的 occupancy probability。

---

## 3. 一阶 occupancy 近似

若只保留体积的一阶矩，可定义

$$
m_1=(N-2)\mathbb{E}[V_A].
$$

于是得到最简单的 occupancy closure：

$$
P_{\mathrm{occ}}^{(1)}
\equiv
1-e^{-m_1}
=
1-\exp\!\bigl(-(N-2)\mathbb{E}[V_A]\bigr).
$$

这不是对 `1-\ell_d` 的严格恒等式，而是把

$$
1-\mathbb{E}[e^{-X}]
$$

替换为

$$
1-e^{-\mathbb{E}[X]}
$$

的第一层 mean-field closure，其中 `X=(N-2)V_A`。

---

## 4. 为什么不是线性矩 `m_1`

如果直接用 `m_1` 作为修正项，即

$$
\frac{\log H}{N} \sim h_0 + c\,m_1,
\qquad
h_0=(1-p_d)(\log N-1),
$$

则会系统高估大体积区间的贡献，因为：

1. `m_1` 把体积影响当成线性可加。
2. 真正决定 non-link / entropy correction 的是“有无至少一个中介点”，这是一个 occupancy 事件。
3. occupancy 事件天然带指数饱和：当 `m_1` 增大时，修正项不应继续线性增长，而应逐渐饱和到 1。

因此

$$
1-e^{-m_1}
$$

比 `m_1` 本身更符合物理结构。

---

## 5. 数值结果

来自 `prediction_a_xi_volume_moments.py`：

### 5.1 用 occupancy 近似的 entropy closure

$$
\frac{\log H}{N}
\approx
0.9343\,h_0 + 0.1315\,P_{\mathrm{occ}}^{(1)}.
$$

对应：

- `Xi_45` median = `11.33`
- observed median = `11.35`
- relative error = `0.2%`

### 5.2 对比

若直接使用线性矩 `m_1`：

- `Xi_45` median = `9.25`
- relative error = `18.5%`

若在 occupancy 近似外再加宽度修正 `std(mu)`：

- `Xi_45` median = `11.09`
- relative error = `2.3%`

结论：主导修正已经由一阶 occupancy 因子给出，宽度修正只提供次级改进。

---

## 6. 与 interval bridge 的关系

此前得到：

$$
1-\ell_d
\approx
0.020 + 1.164(C_1/C_0)
$$

以及

$$
1-\ell_d
\approx
-0.032 + 0.447(C_1/C_0) + 0.647\sqrt{\frac{R-C_0-C_1}{R}}.
$$

这些 bridge 说明 `1-\ell_d` 可以被 interval hierarchy 近似重建。

而本说明进一步表明：

- interval hierarchy 是 `1-\ell_d` 的**微观分解**
- occupancy factor `1-e^{-m_1}` 是 `1-\ell_d` 的**解析 coarse-grained 近似**

两者并不冲突，而是对应“微观分解”与“宏观闭合”的两层描述。

---

## 7. 对 `Xi_{4→5}` 的物理解释

把 occupancy 近似代入 entropy closure 后，`4→5` 屏障可以解释为：

1. 随维度升高，相关点对比例 `p_d` 持续下降。
2. 同时，平均 interval 体积 `E[V_A]` 也下降。
3. 但真正控制 entropy correction 的不是 `E[V_A]` 线性量，而是
   $$1-\exp(-(N-2)E[V_A]).$$
4. 到 `d=5` 时，interval occupancy 已显著收缩，意味着绝大多数相关对都几乎没有“中介深度”。
5. 于是 5D 既在 link 结构上趋于饱和，又在 entropy correction 上失去 mediated-causality reservoir，最终把 `Xi_{4→5}` 推高到异常大的量级。

---

## 8. 当前最可用的论文表述

如果要写进正文，当前最稳妥的说法是：

> The entropy correction is naturally controlled by an interval-occupancy probability rather than by a linear interval-volume moment. For a related pair with Alexandrov volume \(V_A\), the non-link probability is approximately \(1-e^{-(N-2)V_A}\). Averaging over related pairs yields the first-moment closure \(P_{\mathrm{occ}}^{(1)}=1-\exp(-(N-2)\mathbb{E}[V_A])\), which reproduces the median \(\Xi_{4\to5}\) to within \(0.2\%\).

---

## 9. 剩余工作

要把这条链从“强解析近似”推进到“严格推导”，还需要：

1. 明确从
   $$1-\mathbb{E}[e^{-X}]$$
   到
   $$1-e^{-\mathbb{E}[X]}$$
   的误差项控制。
2. 解释 cube boundary correction 对 `E[V_A]` 和 `p_d` 的影响。
3. 判断是否可用 cumulant expansion 写出
   $$1-\ell_d \approx 1-\exp(-\kappa_1 + \kappa_2/2 - \cdots)$$
   的系统修正。
