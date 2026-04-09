# μ(N,κ) Adaptive Reference Prototype

- kappa: 1.0
- n_values: [64]
- ref_reps: 3
- eval_reps: 2

## Summary

|   n |   kappa |   flrw_rank_flat_ref |   flrw_rank_adapt_ref |   lor4d_rank_flat_ref |   lor4d_rank_adapt_ref |   delta_mu_norm |
|----:|--------:|---------------------:|----------------------:|----------------------:|-----------------------:|----------------:|
|  64 |       1 |                    1 |                     2 |                     2 |                      1 |        0.477681 |

## Interpretation

- 若 `flrw_rank_adapt_ref < flrw_rank_flat_ref`，说明曲率自适应参考对 FLRW κ 分支更友好。
- 若 `lor4d_rank_flat_ref` 维持低位，说明平坦参考仍能稳定识别 Lor4D。
- `delta_mu_norm` 衡量 μ(N,κ) 与 μ(N,0) 的几何偏移规模。