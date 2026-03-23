# logH ↔ ∫√(-g)d⁴x Correspondence


## 1. logH Scaling with N

| d | H | slope(logH vs NlogN) | R² | logH/NlogN at N=128 |
|---|---|---------------------|----|--------------------|
| 2 | 0.0 | 0.3941 | 0.9999 | 0.3913 |
| 2 | 0.2 | 0.4036 | 0.9997 | 0.3985 |
| 2 | 0.5 | 0.4213 | 0.9999 | 0.4174 |
| 2 | 1.0 | 0.4555 | 0.9994 | 0.4493 |
| 2 | 2.0 | 0.5084 | 0.9992 | 0.5125 |
| 3 | 0.0 | 0.5251 | 0.9999 | 0.5250 |
| 3 | 0.2 | 0.5481 | 1.0000 | 0.5445 |
| 3 | 0.5 | 0.5690 | 0.9999 | 0.5640 |
| 3 | 1.0 | 0.6175 | 0.9998 | 0.6123 |
| 3 | 2.0 | 0.7202 | 1.0000 | 0.7167 |
| 4 | 0.0 | 0.6107 | 0.9999 | 0.6028 |
| 4 | 0.2 | 0.6235 | 1.0000 | 0.6175 |
| 4 | 0.5 | 0.6431 | 0.9999 | 0.6367 |
| 4 | 1.0 | 0.7169 | 0.9998 | 0.7055 |
| 4 | 2.0 | 0.8009 | 1.0000 | 0.7877 |

## 2. logH/N vs Hubble Parameter

logH/N = combinatorial entropy density. If logH ∝ V = N/ρ,
then logH/N ≈ const. Deviations reveal curvature effects.

| d | N | logH/N (H=0) | logH/N (H=1) | logH/N (H=2) | ρ(H, logH/N) |
|---|---|-------------|-------------|-------------|-------------|
| 2 | 32 | 1.342 | 1.506 | 1.832 | +0.816 |
| 2 | 48 | 1.464 | 1.729 | 2.121 | +0.941 |
| 2 | 64 | 1.601 | 1.825 | 2.189 | +0.914 |
| 2 | 96 | 1.780 | 2.094 | 2.414 | +0.969 |
| 2 | 128 | 1.899 | 2.180 | 2.487 | +0.941 |
| 3 | 32 | 1.778 | 2.085 | 2.429 | +0.875 |
| 3 | 48 | 2.037 | 2.349 | 2.738 | +0.906 |
| 3 | 64 | 2.186 | 2.492 | 2.934 | +0.926 |
| 3 | 96 | 2.375 | 2.821 | 3.253 | +0.965 |
| 3 | 128 | 2.548 | 2.971 | 3.477 | +0.969 |
| 4 | 32 | 2.009 | 2.321 | 2.523 | +0.934 |
| 4 | 48 | 2.238 | 2.656 | 2.910 | +0.969 |
| 4 | 64 | 2.472 | 2.838 | 3.175 | +0.937 |
| 4 | 96 | 2.753 | 3.242 | 3.552 | +0.961 |
| 4 | 128 | 2.925 | 3.423 | 3.822 | +0.973 |

## 3. Order Fraction r vs H

r = n_causal / C(N,2). This is a proxy for the causal structure density.

| d | N | r(H=0) | r(H=1) | r(H=2) | ρ(H, r) |
|---|---|--------|--------|--------|---------|
| 2 | 32 | 0.4798 | 0.3254 | 0.1762 | -0.879 |
| 2 | 48 | 0.5126 | 0.3250 | 0.1629 | -0.969 |
| 2 | 64 | 0.5072 | 0.3342 | 0.1868 | -0.957 |
| 2 | 96 | 0.4946 | 0.2998 | 0.1594 | -0.981 |
| 2 | 128 | 0.5187 | 0.3134 | 0.1901 | -0.973 |
| 3 | 32 | 0.2565 | 0.1097 | 0.0169 | -0.895 |
| 3 | 48 | 0.2686 | 0.1021 | 0.0177 | -0.934 |
| 3 | 64 | 0.2640 | 0.1096 | 0.0275 | -0.945 |
| 3 | 96 | 0.2923 | 0.0972 | 0.0253 | -0.977 |
| 3 | 128 | 0.2988 | 0.1136 | 0.0193 | -0.977 |
| 4 | 32 | 0.1536 | 0.0343 | 0.0028 | -0.911 |
| 4 | 48 | 0.1810 | 0.0371 | 0.0014 | -0.962 |
| 4 | 64 | 0.1686 | 0.0393 | 0.0021 | -0.949 |
| 4 | 96 | 0.1638 | 0.0300 | 0.0020 | -0.945 |
| 4 | 128 | 0.1779 | 0.0318 | 0.0018 | -0.981 |

## 4. logH vs n_causal (Causal Pair Count)

If logH ~ f(n_causal), then the volume term is mediated by causal structure.

| d | ρ(log_H, n_causal) | ρ(log_H, log(n_causal)) | ρ(logH/N, r) |
|---|-------------------|------------------------|--------------|
| 2 | +0.836 | +0.836 | -0.608 |
| 3 | +0.601 | +0.601 | -0.546 |
| 4 | +0.408 | +0.408 | -0.528 |

## 5. Theoretical Decomposition

For a poset with N elements and order fraction r:
- Antichain (r=0): logH = log(N!) = N·log(N) - N + O(log(N))
- Total order (r≈0.5): logH = 0
- General: logH ≈ N·log(N)·(1 - g(r)) where g captures ordering constraints

The key question: does logH/N decrease monotonically with H?
If yes, then logH encodes an effective volume-like quantity:
fewer causal constraints (lower H → more pairs) → lower combinatorial
entropy per element → the structure is MORE constrained.

### Empirical direction check

- d=2, N=128: H=0.0: 243.0 → H=0.2: 247.5 → H=0.5: 259.2 → H=1.0: 279.0 → H=2.0: 318.3
- d=3, N=128: H=0.0: 326.1 → H=0.2: 338.2 → H=0.5: 350.3 → H=1.0: 380.3 → H=2.0: 445.1
- d=4, N=128: H=0.0: 374.4 → H=0.2: 383.5 → H=0.5: 395.4 → H=1.0: 438.1 → H=2.0: 489.2


## 6. Connection to Einstein-Hilbert

### The Volume Connection

In CST, N = ρ·V where V is the spacetime volume and ρ the sprinkling density.
For fixed N (our experiment), the sprinkled region has FIXED volume V = N/ρ.
But the CAUSAL structure within that volume changes with H.

The EH action has two parts:
  S_EH = (1/16πG) ∫ R·√(-g) d⁴x = curvature_part + cosmological_part
  = (1/16πG)(R_dS·V - 2Λ·V)

In our F7 functional:
- logH encodes the COMBINATORIAL volume/entropy of the causal set
- wall encodes the CURVATURE admissibility (sigmoid → step function)
- Together: F7 ≈ S_combinatorial + wall(curvature)

The correspondence is therefore:
  logH ↔ cosmological term (volume-extensive entropy)
  wall ↔ curvature threshold (admissibility)
  F7 ↔ S_EH (with wall as finite-size rounding of the curvature cut)
