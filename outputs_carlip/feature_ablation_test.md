# Feature Ablation + Minimal Complete Basis Test


## 1. Lor4D Rank by Configuration

| Config | N=16 | N=20 | N=28 | N=36 | N=48 | N=64 | N=96 | N=128 | All#1? |
|--------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:-----:|
| Full (d,c,w) | #1 | #1 | #1 | #1 | #1 | #1 | #1 | #1 | ✅ |
| Drop d: (c,w) | #2 | #1 | #1 | #1 | #1 | #1 | #1 | #1 | ❌ |
| Drop c: (d,w) | #1 | #1 | #1 | #1 | #1 | #1 | #1 | #1 | ✅ |
| Drop w: (d,c) | #1 | #1 | #1 | #1 | #1 | #1 | #1 | #1 | ✅ |
| d only | #2 | #2 | #2 | #3 | #2 | #2 | #2 | #1 | ❌ |
| c only | #2 | #4 | #1 | #1 | #1 | #1 | #1 | #1 | ❌ |
| w only | #1 | #1 | #2 | #2 | #2 | #2 | #1 | #3 | ❌ |
| +height: (d,c,w,h) | #1 | #1 | #1 | #1 | #1 | #1 | #1 | #1 | ✅ |
| +order_frac: (d,c,w,R) | #1 | #1 | #1 | #1 | #1 | #1 | #1 | #1 | ✅ |


## 2. Margin (to nearest non-Lor) by Configuration

| Config | N=16 | N=20 | N=28 | N=36 | N=48 | N=64 | N=96 | N=128 | Mean |
|--------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:----:|
| Full (d,c,w) | 9.6 | 7.3 | 29.9 | 36.0 | 63.7 | 153.4 | 169.9 | 342.1 | 101.5 |
| Drop d: (c,w) | -0.2 | 0.6 | 0.8 | 2.8 | 13.9 | 21.2 | 10.0 | 5.8 | 6.9 |
| Drop c: (d,w) | 2.5 | 3.1 | 4.6 | 4.7 | 9.4 | 33.7 | 142.7 | 280.4 | 60.1 |
| Drop w: (d,c) | 1.0 | 0.6 | 4.0 | 6.0 | 19.1 | 33.7 | 56.2 | 96.6 | 27.2 |
| d only | -0.5 | -0.3 | -0.3 | -0.7 | -0.7 | -0.6 | -0.0 | 0.4 | -0.3 |
| c only | -0.3 | -0.7 | 0.1 | 1.6 | 12.1 | 21.3 | 3.9 | 5.4 | 5.4 |
| w only | 0.1 | 0.5 | -0.8 | -0.0 | -0.7 | -0.5 | 0.3 | -0.5 | -0.2 |
| +height: (d,c,w,h) | 9.7 | 6.4 | 30.2 | 37.4 | 85.2 | 159.9 | 169.4 | 353.4 | 106.5 |
| +order_frac: (d,c,w,R) | 9.0 | 6.7 | 32.5 | 37.3 | 64.9 | 177.0 | 214.7 | 394.4 | 117.0 |


## 3. Ablation Impact Analysis

### Dropping d_eff

- **Failures** (Lor4D not #1): N = [16]
- **Mean margin loss**: 94.6
- **Impact**: CRITICAL

### Dropping C₁/C₀

- **No failures**: Lor4D still #1 at all N
- **Mean margin loss**: 41.3
- **Impact**: Moderate

### Dropping width

- **No failures**: Lor4D still #1 at all N
- **Mean margin loss**: 74.3
- **Impact**: Moderate


## 4. Fourth Feature Benefit

- **height_ratio**: mean margin gain = 5.0
- **order_frac**: mean margin gain = 15.5

**Conclusion**: If gains are small and no new failures are resolved, the triple is the minimal complete basis.


## 5. Verdict

🟡 **Critical features**: ['d_eff']
Other features are helpful but not strictly necessary for rank #1.