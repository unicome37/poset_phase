# Mahalanobis Small-N Prior Test — Anchored d*=4

Date: 2026-03-27 11:01
Families: 25
N values: [16, 20, 28, 48, 64, 96]
Reps: 20

## 1. Definition

- Baseline (full): `S = (I-mu)^T Sigma^{-1} (I-mu)` with mu,Sigma from Lor4D ensemble at each N.
- Anchored (full): same Sigma, but `mu_anchor = (4.0, mu_c, mu_w)` (physical d*=4 prior).
- Anchored (diag): assume feature independence, `Sigma^{-1} -> diag(1/var_i)` from Lor4D.

## 2. Lor4D Rank Summary

| N | Baseline winner | Lor4D rank | OAS winner | Lor4D rank | Anchored(full) winner | Lor4D rank | Anchored(diag) winner | Lor4D rank |
|---|---|---:|---|---:|---|---:|---|---:|
| 16 | Lor4D | 1 | Lor4D | 1 | Lor4D | 1 | Lor4D | 1 |
| 20 | Lor4D | 1 | Lor4D | 1 | Lor4D | 1 | Lor4D | 1 |
| 28 | Lor4D | 1 | Lor4D | 1 | Lor4D | 1 | Lor5D | 2 |
| 48 | Lor4D | 1 | Lor4D | 1 | Lor4D | 1 | Lor4D | 1 |
| 64 | Lor4D | 1 | Lor4D | 1 | Lor4D | 1 | Lor4D | 1 |
| 96 | Lor4D | 1 | Lor4D | 1 | Lor4D | 1 | Lor4D | 1 |

## 3. N=16 Top-5 (Baseline vs Anchored)

### Baseline Top 5 (N=16)

| Rank | Family | Score |
|---:|---|---:|
| 1 | Lor4D | 2.8500 |
| 2 | KR_2layer | 4.4958 |
| 3 | Lor5D | 5.4829 |
| 4 | KR_like | 18.3104 |
| 5 | Lor3D | 19.9695 |

### Anchored Top 5 (N=16)

| Rank | Family | Score |
|---:|---|---:|
| 1 | Lor4D | 2.8718 |
| 2 | KR_2layer | 4.7550 |
| 3 | Lor5D | 5.1876 |
| 4 | KR_like | 19.4501 |
| 5 | Lor3D | 20.3864 |

### Anchored (diag) Top 5 (N=16)

| Rank | Family | Score |
|---:|---|---:|
| 1 | Lor4D | 2.8702 |
| 2 | KR_2layer | 4.9740 |
| 3 | Lor5D | 5.5333 |
| 4 | KR_like | 19.9884 |
| 5 | Lor3D | 24.2535 |

## 4. mu_d(N) Drift (Lor4D empirical centroid)

| N | mu_d | |mu_d-4| | mu_c | mu_w |
|---|---:|---:|---:|---:|
| 16 | 3.9504 | 0.0496 | 0.0596 | 0.5531 |
| 20 | 3.9175 | 0.0825 | 0.0776 | 0.5300 |
| 28 | 3.9268 | 0.0732 | 0.1297 | 0.4893 |
| 48 | 3.9767 | 0.0233 | 0.1806 | 0.4250 |
| 64 | 3.9618 | 0.0382 | 0.2266 | 0.3680 |
| 96 | 3.9653 | 0.0347 | 0.2680 | 0.3229 |

Runtime: 19.7s
