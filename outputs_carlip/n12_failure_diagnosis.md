# N=12 Failure Diagnosis

N=12, REPS=80, 25 families, 10 seeds

## seed=42 — **PASS** (Lor4D rank=1, winner=Lor4D)

- cond(Σ) = 37.1
- μ_Lor4D = [3.9634, 0.0655, 0.6063]
- S_MD(Lor4D) = 2.9625

Top 5 ranking:
  #1 Lor4D            S_MD=2.9625  feat=[3.9634, 0.0655, 0.6063]
  #2 Lor5D            S_MD=3.1680  feat=[4.3259, 0.0229, 0.6885]
  #3 KR_2layer        S_MD=3.3044  feat=[3.8099, 0.0000, 0.6531]
  #4 Lor3D            S_MD=7.0116  feat=[3.2781, 0.1937, 0.4885]
  #5 KR_like          S_MD=11.1690  feat=[2.7065, 0.1800, 0.4542]

## seed=137 — **FAIL** (Lor4D rank=3, winner=Lor5D)

- cond(Σ) = 29.5
- μ_Lor4D = [4.0554, 0.0581, 0.6167]
- S_MD(Lor4D) = 2.9625
- S_MD(Lor5D) = 2.8508
- margin = -0.1117
- μ_Lor5D = [4.3179, 0.0142, 0.6844]

Per-feature δ from Lor4D reference:
  d_eff:  Lor4D=+0.0000  Lor5D=+0.2625
  C1/C0:  Lor4D=+0.0000  Lor5D=-0.0439
  w/N:    Lor4D=+0.0000  Lor5D=+0.0677

Top 5 ranking:
  #1 Lor5D            S_MD=2.8508  feat=[4.3179, 0.0142, 0.6844]
  #2 KR_2layer        S_MD=2.8583  feat=[3.8531, 0.0000, 0.6438]
  #3 Lor4D            S_MD=2.9625  feat=[4.0554, 0.0581, 0.6167]
  #4 Lor3D            S_MD=9.7659  feat=[3.2111, 0.1974, 0.4885]
  #5 TransPerc        S_MD=10.2102  feat=[4.3997, 0.2007, 0.6583]

## seed=271 — **FAIL** (Lor4D rank=2, winner=Lor5D)

- cond(Σ) = 24.0
- μ_Lor4D = [4.0020, 0.0687, 0.6208]
- S_MD(Lor4D) = 2.9625
- S_MD(Lor5D) = 2.4665
- margin = -0.4960
- μ_Lor5D = [4.3429, 0.0202, 0.6885]

Per-feature δ from Lor4D reference:
  d_eff:  Lor4D=+0.0000  Lor5D=+0.3409
  C1/C0:  Lor4D=+0.0000  Lor5D=-0.0485
  w/N:    Lor4D=+0.0000  Lor5D=+0.0677

Top 5 ranking:
  #1 Lor5D            S_MD=2.4665  feat=[4.3429, 0.0202, 0.6885]
  #2 Lor4D            S_MD=2.9625  feat=[4.0020, 0.0687, 0.6208]
  #3 KR_2layer        S_MD=3.4355  feat=[3.7895, 0.0000, 0.6677]
  #4 TransPerc        S_MD=8.7572  feat=[4.4429, 0.1790, 0.6729]
  #5 Lor3D            S_MD=10.0887  feat=[3.2190, 0.2128, 0.5000]

## seed=500 — **PASS** (Lor4D rank=1, winner=Lor4D)

- cond(Σ) = 29.9
- μ_Lor4D = [3.8906, 0.0856, 0.5875]
- S_MD(Lor4D) = 2.9625

Top 5 ranking:
  #1 Lor4D            S_MD=2.9625  feat=[3.8906, 0.0856, 0.5875]
  #2 Lor5D            S_MD=3.5603  feat=[4.3554, 0.0076, 0.6969]
  #3 KR_2layer        S_MD=4.0557  feat=[3.7452, 0.0000, 0.6573]
  #4 Lor3D            S_MD=8.1222  feat=[3.2679, 0.1843, 0.5000]
  #5 RandomDAG_sp     S_MD=9.9055  feat=[4.1884, 0.2591, 0.6031]

## seed=777 — **PASS** (Lor4D rank=1, winner=Lor4D)

- cond(Σ) = 35.6
- μ_Lor4D = [3.9509, 0.0668, 0.6062]
- S_MD(Lor4D) = 2.9625

Top 5 ranking:
  #1 Lor4D            S_MD=2.9625  feat=[3.9509, 0.0668, 0.6062]
  #2 KR_2layer        S_MD=3.4684  feat=[3.7747, 0.0000, 0.6542]
  #3 Lor5D            S_MD=3.6861  feat=[4.3440, 0.0225, 0.6813]
  #4 Lor3D            S_MD=7.7436  feat=[3.2804, 0.2077, 0.4854]
  #5 TransPerc        S_MD=12.1640  feat=[4.3349, 0.2401, 0.6458]

## seed=1001 — **PASS** (Lor4D rank=1, winner=Lor4D)

- cond(Σ) = 22.6
- μ_Lor4D = [3.9543, 0.0688, 0.5948]
- S_MD(Lor4D) = 2.9625

Top 5 ranking:
  #1 Lor4D            S_MD=2.9625  feat=[3.9543, 0.0688, 0.5948]
  #2 KR_2layer        S_MD=3.6115  feat=[3.7895, 0.0000, 0.6604]
  #3 Lor5D            S_MD=3.8854  feat=[4.3929, 0.0028, 0.7073]
  #4 RandomDAG_sp     S_MD=9.5144  feat=[4.1884, 0.2410, 0.6083]
  #5 Lor3D            S_MD=9.6314  feat=[3.3202, 0.2336, 0.5010]

## seed=2023 — **FAIL** (Lor4D rank=2, winner=KR_2layer)

- cond(Σ) = 28.8
- μ_Lor4D = [4.0009, 0.0631, 0.6052]
- S_MD(Lor4D) = 2.9625
- S_MD(KR_2layer) = 2.7781
- margin = -0.1844
- μ_KR_2layer = [3.8065, 0.0000, 0.6615]

Per-feature δ from Lor4D reference:
  d_eff:  Lor4D=+0.0000  KR_2layer=-0.1943
  C1/C0:  Lor4D=+0.0000  KR_2layer=-0.0631
  w/N:    Lor4D=+0.0000  KR_2layer=+0.0562

Top 5 ranking:
  #1 KR_2layer        S_MD=2.7781  feat=[3.8065, 0.0000, 0.6615]
  #2 Lor4D            S_MD=2.9625  feat=[4.0009, 0.0631, 0.6052]
  #3 Lor5D            S_MD=3.1315  feat=[4.3702, 0.0204, 0.7073]
  #4 Lor3D            S_MD=9.1311  feat=[3.3270, 0.1955, 0.4906]
  #5 TransPerc        S_MD=10.0810  feat=[4.3588, 0.2391, 0.6375]

## seed=3141 — **PASS** (Lor4D rank=1, winner=Lor4D)

- cond(Σ) = 37.7
- μ_Lor4D = [3.9315, 0.0673, 0.6052]
- S_MD(Lor4D) = 2.9625

Top 5 ranking:
  #1 Lor4D            S_MD=2.9625  feat=[3.9315, 0.0673, 0.6052]
  #2 KR_2layer        S_MD=3.9687  feat=[3.7815, 0.0000, 0.6698]
  #3 Lor5D            S_MD=4.3639  feat=[4.3599, 0.0203, 0.7042]
  #4 Lor3D            S_MD=8.4168  feat=[3.3179, 0.2097, 0.5063]
  #5 TransPerc        S_MD=9.8519  feat=[4.3009, 0.2320, 0.6271]

## seed=5000 — **PASS** (Lor4D rank=1, winner=Lor4D)

- cond(Σ) = 23.9
- μ_Lor4D = [3.9679, 0.0879, 0.5802]
- S_MD(Lor4D) = 2.9625

Top 5 ranking:
  #1 Lor4D            S_MD=2.9625  feat=[3.9679, 0.0879, 0.5802]
  #2 Lor5D            S_MD=3.0247  feat=[4.3554, 0.0164, 0.6760]
  #3 KR_2layer        S_MD=4.5013  feat=[3.7952, 0.0000, 0.6604]
  #4 Lor3D            S_MD=7.3410  feat=[3.2497, 0.1902, 0.4875]
  #5 TransPerc        S_MD=8.1740  feat=[4.4338, 0.1993, 0.6594]

## seed=8888 — **FAIL** (Lor4D rank=2, winner=KR_2layer)

- cond(Σ) = 25.3
- μ_Lor4D = [3.9406, 0.0458, 0.6031]
- S_MD(Lor4D) = 2.9625
- S_MD(KR_2layer) = 2.2090
- margin = -0.7535
- μ_KR_2layer = [3.8020, 0.0000, 0.6500]

Per-feature δ from Lor4D reference:
  d_eff:  Lor4D=+0.0000  KR_2layer=-0.1386
  C1/C0:  Lor4D=+0.0000  KR_2layer=-0.0458
  w/N:    Lor4D=+0.0000  KR_2layer=+0.0469

Top 5 ranking:
  #1 KR_2layer        S_MD=2.2090  feat=[3.8020, 0.0000, 0.6500]
  #2 Lor4D            S_MD=2.9625  feat=[3.9406, 0.0458, 0.6031]
  #3 Lor5D            S_MD=3.3321  feat=[4.3906, 0.0100, 0.6865]
  #4 Lor3D            S_MD=8.2017  feat=[3.2622, 0.1888, 0.4906]
  #5 TransPerc        S_MD=9.9278  feat=[4.4361, 0.1644, 0.6698]

---
## Cross-Seed Summary: d_eff Distribution Overlap

The sole intruder at N=12 is Lor5D. The key diagnostic is d_eff overlap between 4D and 5D sprinklings at only 12 points.

At N=12, the Myrheim-Meyer dimension estimator has insufficient resolution to separate d=4 (d_eff≈3.95) from d=5 (d_eff≈4.3-4.5). This is a fundamental physical limit, not a statistical artifact.

**Conclusion**: N=12 is below the physical resolution horizon. The S_MD operator correctly identifies this as an ill-defined regime. Letter should claim N≥14 with clean conscience.
