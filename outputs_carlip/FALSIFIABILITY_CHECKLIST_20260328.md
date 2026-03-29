# Falsifiability Execution Checklist (v1)

## 预注册锁定（执行前必须完成）

- [ ] 锁定 C1/C2/C3 Hard Fail 阈值（不得事后改）
- [ ] 锁定 N 网格与 seed 数
- [ ] 锁定 family 扩展名单（F1）
- [ ] 锁定背景参数窗口（F3）
- [ ] 锁定候选特征基集合（F4）

## F1: 扩展家族库抗压

- [x] 生成扩展 family 清单（12–20）
- [x] 跑 N=12..256, seeds=10, REPS=80
- [x] 汇总连续 N 压制段统计
- [x] 给出 C1 Pass/SoftFail/HardFail（当前结果：Pass / Hard fail = NO）

## F2: Turn-on 再估计

- [ ] 跑 N=10..24（步长2）, seeds=20, REPS=120
- [ ] 输出 N_id 后验区间
- [ ] 记录 min-margin 连续稳定段

## F3: 背景响应可证伪

- [ ] de Sitter 扩展到 N=1536/2048
- [ ] FLRW metric-faithful 版本与 proxy 对照
- [ ] Schwarzschild 弱场严格采样版本
- [~] 给出 top-2 维持率与 C2 判定（当前 lowN split：de Sitter Pass，Schwarzschild Pass，FLRW `\kappa=1.0` Hard Fail）

## F4: 特征基挑战

- [ ] 在同库同配置下评估 B0/B1/B2/B3
- [ ] 输出 margin、rank稳定率、N_id 变化
- [ ] 给出 C3 判定

## F5: 合成对手

- [ ] 构造匹配均值+协方差的对手家族
- [ ] 测试是否穿透二阶 identity layer
- [ ] 给出“是否需要三阶层”结论

## 文稿联动

- [x] 若任一 Hard Fail，立即降级摘要与结论表述（已因 FLRW lowN split 结果降级曲率稳健性口径）
- [ ] 若 Soft Fail，补充 limitations 条件说明
- [ ] 全 Pass 时仍保持“tested library 范围内”口径
