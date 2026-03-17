"""
Prediction A — Extended BD Action Analysis

Three advances in one experiment:

1. GENERALIZED BD ACTION (d>2):
   S_BD^(d) uses interval counting (Alexandrov intervals),
   not just links. For d=2: S = N - 2*n_links.
   For general d: S = Σ_k α_k(d) * C_k
   where C_k = number of intervals of cardinality k.
   
   We use the d=4 form from Benincasa-Dowker (2010):
     S_BD^(4) = (4/√6) * [ -N + (9/√(16π)) * C_0
                              - (16/√(6π)) * C_1
                              + (64/(3√(4π))) * C_2
                              - ... ]
   
   Simplified practical form: use the ratio of order-2 intervals
   to links as a dimension-sensitive correction.

2. BD + GEOMETRIC PENALTY HYBRID:
   score = -β*logH + λ_BD * S_BD/N + γ_geo * P_geo
   Tests whether combining BD with geometric structure penalty
   gives sharper 4D selection.

3. FINITE-SIZE SCALING:
   Extend N to 60, 68 to check 4D plateau stability.
"""
from __future__ import annotations
from pathlib import Path
from math import comb

import numpy as np
import pandas as pd

from experiment import FAMILIES
from runtime_utils import estimate_entropy_by_family
from observables import neutral_penalty
from observables_geo import geometric_components

OUT_DIR = Path("outputs_exploratory/prediction_a_bd_extended")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# --- Configuration ---
N_VALUES = [20, 28, 36, 44, 52, 60, 68]
FAMILY_LIST = [
    "lorentzian_like_2d",
    "lorentzian_like_3d",
    "lorentzian_like_4d",
    "lorentzian_like_5d",
    "KR_like",
]
SAMPLES = 4
SIS_RUNS = 4096
BETA = 1.0
DEFAULT_EXACT_THRESHOLD = 24
FAMILY_EXACT_THRESHOLDS = {
    "lorentzian_like_2d": 104,
    "lorentzian_like_3d": 24,
    "lorentzian_like_4d": 24,
    "lorentzian_like_5d": 24,
    "KR_like": 24,
}
CONSISTENCY_WEIGHT = 0.682882


# ===== BD Observable Functions =====

def hasse_links(poset) -> int:
    """Count covering relations (Hasse diagram edges) = C_0 intervals."""
    c = poset.closure.astype(np.uint8, copy=False)
    has_intermediate = (c @ c).astype(bool, copy=False)
    cover = poset.closure & ~has_intermediate
    np.fill_diagonal(cover, False)
    return int(cover.sum())


def interval_counts(poset, max_k: int = 3) -> dict[int, int]:
    """Count Alexandrov intervals by interior cardinality.

    C_k = number of pairs (x, y) with x < y and exactly k elements
    strictly between x and y.

    C_0 = links (covering relations)
    C_1 = intervals with 1 interior element
    C_2 = intervals with 2 interior elements
    ...
    """
    c = poset.closure
    n = poset.n
    counts = {k: 0 for k in range(max_k + 1)}

    # For each comparable pair (i, j), count interior
    pairs = np.argwhere(c)
    for i, j in pairs:
        if i == j:
            continue
        # Interior: elements strictly between i and j
        interior = c[i] & c[:, j]
        interior[i] = False
        interior[j] = False
        k = int(interior.sum())
        if k <= max_k:
            counts[k] += 1

    return counts


def bd_action_d2(n: int, c0: int) -> float:
    """BD action (d=2): S = N - 2*C_0."""
    return float(n - 2 * c0)


def bd_action_d4(n: int, c0: int, c1: int, c2: int) -> float:
    """BD action (d=4) — Benincasa-Dowker generalized form.

    S_BD^(4) ∝ N - α₁·C₀ + α₂·C₁ - α₃·C₂
    
    From Benincasa & Dowker (2010), for d=4 Minkowski:
    The coefficients alternate in sign and grow.
    We use the normalized form with standard coefficients.
    
    Note: exact coefficients depend on the discretization,
    but the KEY physics is: C_0, C_1, C_2 all penalize
    causal richness, with higher intervals being MORE
    dimension-sensitive.
    """
    # Benincasa-Dowker d=4 coefficients (from the literature):
    # S = (4/√6) * [N - (9/(4√π))·C₀ + (16/(3√(6π)))·C₁
    #                  - (64/(9√(4π)))·C₂]
    # Let's use a simplified proportional form:
    alpha1 = 9.0 / (4.0 * np.sqrt(np.pi))       # ≈ 1.270
    alpha2 = 16.0 / (3.0 * np.sqrt(6.0 * np.pi)) # ≈ 1.229
    alpha3 = 64.0 / (9.0 * np.sqrt(4.0 * np.pi)) # ≈ 2.006

    prefactor = 4.0 / np.sqrt(6.0)  # ≈ 1.633
    return float(prefactor * (n - alpha1 * c0 + alpha2 * c1 - alpha3 * c2))


def geometric_penalty_consistency(geo: dict[str, float]) -> float:
    """A2_replace_dim_with_consistency variant."""
    width_height = 2.0 * float(geo["geo_width_height"])
    comp_window = 6.0 * float(geo["geo_comparability_window"])
    cover = 3.0 * float(geo["geo_cover_density"])
    interval_profile = 4.0 * float(geo["geo_interval_profile"])
    interval_shape = 4.0 * float(geo["geo_interval_shape"])
    layer_smoothness = 2.0 * float(geo["geo_layer_smoothness"])
    consistency = CONSISTENCY_WEIGHT * float(geo["geo_dim_consistency"])
    return (width_height + consistency + comp_window + cover
            + interval_profile + interval_shape + layer_smoothness)


# ===== Data Collection =====

print("=" * 72)
print("EXTENDED BD ACTION ANALYSIS")
print("  1. Generalized BD (d=2 and d=4 interval counting)")
print("  2. BD + geometric penalty hybrid")
print("  3. Finite-size scaling to N=60, 68")
print("=" * 72)

print(f"\nN: {N_VALUES}, Families: {len(FAMILY_LIST)}, Samples: {SAMPLES}")
print("Collecting base observables...\n")

base_rows = []
for n in N_VALUES:
    print(f"  N={n}", end="", flush=True)
    for family in FAMILY_LIST:
        generator = FAMILIES[family]
        for sid in range(SAMPLES):
            seed = 960000 + 1000 * n + sid
            poset = generator(n=n, seed=seed)

            # Entropy
            log_h, method = estimate_entropy_by_family(
                poset, family=family, sis_runs=SIS_RUNS, seed=seed,
                default_exact_threshold=DEFAULT_EXACT_THRESHOLD,
                family_exact_thresholds=FAMILY_EXACT_THRESHOLDS,
            )

            # Interval counts
            ic = interval_counts(poset, max_k=3)
            c0, c1, c2, c3 = ic[0], ic[1], ic[2], ic[3]

            # BD actions
            s_d2 = bd_action_d2(n, c0)
            s_d4 = bd_action_d4(n, c0, c1, c2)

            # Geometric penalty (consistency variant)
            geo = geometric_components(poset)
            p_geo = geometric_penalty_consistency(geo)
            p_neutral = neutral_penalty(poset)

            n_relations = int(poset.closure.sum())
            rel_density = n_relations / (n * (n - 1) / 2) if n > 1 else 0.0

            base_rows.append({
                "n": n, "family": family, "sample_id": sid, "seed": seed,
                "log_H": float(log_h), "method": method,
                "C0_links": c0, "C1": c1, "C2": c2, "C3": c3,
                "n_relations": n_relations, "rel_density": rel_density,
                "S_BD_d2": s_d2, "S_BD_d2_norm": s_d2 / max(n, 1),
                "S_BD_d4": s_d4, "S_BD_d4_norm": s_d4 / max(n, 1),
                "P_geo_consistency": p_geo,
                "P_neutral": p_neutral,
            })
        print(f" [{family[:5]}]", end="", flush=True)
    print()

base = pd.DataFrame(base_rows)
base.to_csv(OUT_DIR / "base_observables.csv", index=False)

# ===== Analysis 1: Interval Count Profiles =====
print("\n" + "=" * 72)
print("ANALYSIS 1: INTERVAL COUNT PROFILES BY DIMENSION")
print("=" * 72)

lor_families = [f for f in FAMILY_LIST if f.startswith("lorentzian_like_")]
for n in N_VALUES:
    print(f"\n  N={n}:")
    for fam in lor_families + ["KR_like"]:
        sub = base[(base["n"] == n) & (base["family"] == fam)]
        label = fam.replace("lorentzian_like_", "").replace("KR_like", "KR")
        print(f"    {label:>4}: C0={sub['C0_links'].mean():7.1f}  "
              f"C1={sub['C1'].mean():7.1f}  C2={sub['C2'].mean():7.1f}  "
              f"C3={sub['C3'].mean():7.1f}  "
              f"S_d2/N={sub['S_BD_d2_norm'].mean():+.3f}  "
              f"S_d4/N={sub['S_BD_d4_norm'].mean():+.3f}  "
              f"logH={sub['log_H'].mean():7.2f}")

# ===== Analysis 2: Scoring — Multiple Action Variants =====
print("\n" + "=" * 72)
print("ANALYSIS 2: DIMENSION SELECTION UNDER DIFFERENT ACTIONS")
print("=" * 72)

# Lambda scans for different action combinations
LAMBDA_BD = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 15.0]
GAMMA_GEO = [0.0, 0.5, 1.0, 2.0]

ACTION_VARIANTS = {
    "pure_entropy": lambda row, lbd, ggeo: -BETA * row["log_H"],
    "BD_d2": lambda row, lbd, ggeo: -BETA * row["log_H"] + lbd * row["S_BD_d2_norm"],
    "BD_d4": lambda row, lbd, ggeo: -BETA * row["log_H"] + lbd * row["S_BD_d4_norm"],
    "BD_d2+geo": lambda row, lbd, ggeo: -BETA * row["log_H"] + lbd * row["S_BD_d2_norm"] + ggeo * row["P_geo_consistency"],
    "BD_d4+geo": lambda row, lbd, ggeo: -BETA * row["log_H"] + lbd * row["S_BD_d4_norm"] + ggeo * row["P_geo_consistency"],
}

score_rows = []
for _, row in base.iterrows():
    for lbd in LAMBDA_BD:
        for ggeo in GAMMA_GEO:
            for vname, scorer in ACTION_VARIANTS.items():
                s = scorer(row, lbd, ggeo)
                score_rows.append({
                    "n": row["n"], "family": row["family"],
                    "sample_id": row["sample_id"],
                    "lambda_bd": lbd, "gamma_geo": ggeo,
                    "variant": vname, "score": s,
                })

scores = pd.DataFrame(score_rows)

# Winner analysis (Lorentzian only)
lor_scores = scores[scores["family"].str.startswith("lorentzian_like_")]
summary = lor_scores.groupby(["n", "lambda_bd", "gamma_geo", "variant", "family"]).agg(
    mean_score=("score", "mean"),
).reset_index()

winner_rows = []
for (n, lbd, ggeo, vname), grp in summary.groupby(["n", "lambda_bd", "gamma_geo", "variant"]):
    ordered = grp.sort_values("mean_score")
    best = ordered.iloc[0]
    second = ordered.iloc[1] if len(ordered) > 1 else ordered.iloc[0]
    margin = float(second["mean_score"] - best["mean_score"])
    winner_rows.append({
        "n": int(n), "lambda_bd": lbd, "gamma_geo": ggeo,
        "variant": vname,
        "winner": best["family"].replace("lorentzian_like_", ""),
        "margin": margin,
    })

winners = pd.DataFrame(winner_rows)
winners.to_csv(OUT_DIR / "winners_extended.csv", index=False)

# ===== Print Winner Tables =====

# 2a. BD_d2 vs BD_d4 comparison (γ_geo=0)
print("\n--- BD_d2 vs BD_d4: Winner by (λ, N) [γ_geo=0, Lorentzian only] ---")
for vname in ["BD_d2", "BD_d4"]:
    print(f"\n  {vname}:")
    w = winners[(winners["variant"] == vname) & (winners["gamma_geo"] == 0.0)]
    print(f"  {'λ':>5}", end="")
    for n in N_VALUES:
        print(f"  {n:>4}", end="")
    print()
    for lbd in LAMBDA_BD:
        print(f"  {lbd:5.1f}", end="")
        for n in N_VALUES:
            row = w[(w["n"] == n) & (w["lambda_bd"] == lbd)]
            if not row.empty:
                print(f"  {row.iloc[0]['winner']:>4}", end="")
            else:
                print(f"  {'?':>4}", end="")
        print()

# 2b. BD_d2+geo hybrid: find λ, γ_geo where 4D wins most
print("\n--- BD_d2+geo hybrid: 4D win count by (λ_BD, γ_geo) ---")
w_hybrid = winners[winners["variant"] == "BD_d2+geo"]
print(f"  {'λ\\γ':>5}", end="")
for ggeo in GAMMA_GEO:
    print(f"  γ={ggeo:.1f}", end="")
print()
for lbd in LAMBDA_BD:
    print(f"  {lbd:5.1f}", end="")
    for ggeo in GAMMA_GEO:
        w = w_hybrid[(w_hybrid["lambda_bd"] == lbd) & (w_hybrid["gamma_geo"] == ggeo)]
        n4d = (w["winner"] == "4d").sum()
        total = len(w)
        label = f"{n4d}/{total}"
        if n4d == total:
            label += "★"
        print(f"  {label:>6}", end="")
    print()

# 2c. BD_d4+geo hybrid
print("\n--- BD_d4+geo hybrid: 4D win count by (λ_BD, γ_geo) ---")
w_hybrid4 = winners[winners["variant"] == "BD_d4+geo"]
print(f"  {'λ\\γ':>5}", end="")
for ggeo in GAMMA_GEO:
    print(f"  γ={ggeo:.1f}", end="")
print()
for lbd in LAMBDA_BD:
    print(f"  {lbd:5.1f}", end="")
    for ggeo in GAMMA_GEO:
        w = w_hybrid4[(w_hybrid4["lambda_bd"] == lbd) & (w_hybrid4["gamma_geo"] == ggeo)]
        n4d = (w["winner"] == "4d").sum()
        total = len(w)
        label = f"{n4d}/{total}"
        if n4d == total:
            label += "★"
        print(f"  {label:>6}", end="")
    print()

# 2d. Best (λ, γ) for 4D selection across all variants
print("\n--- Best parameter combination for unanimous 4D selection ---")
for vname in ["BD_d2", "BD_d4", "BD_d2+geo", "BD_d4+geo"]:
    w = winners[winners["variant"] == vname]
    best_count = 0
    best_params = None
    for (lbd, ggeo), grp in w.groupby(["lambda_bd", "gamma_geo"]):
        n4d = (grp["winner"] == "4d").sum()
        avg_margin = grp[grp["winner"] == "4d"]["margin"].mean() if n4d > 0 else 0
        if n4d > best_count or (n4d == best_count and avg_margin > (best_margin if best_params else 0)):
            best_count = n4d
            best_params = (lbd, ggeo)
            best_margin = avg_margin
    if best_params:
        print(f"  {vname:>12}: λ={best_params[0]:.1f}, γ_geo={best_params[1]:.1f} → "
              f"4D wins {best_count}/{len(N_VALUES)}, avg_margin={best_margin:.3f}")

# ===== Analysis 3: Finite-Size Scaling (N=60, 68 check) =====
print("\n" + "=" * 72)
print("ANALYSIS 3: FINITE-SIZE SCALING — Does 4D plateau hold at larger N?")
print("=" * 72)

# At λ=6 (where BD_d2 gave unanimous 4D in original experiment)
for lbd_check in [6.0, 8.0]:
    print(f"\n  λ_BD = {lbd_check}, γ_geo = 0 (BD_d2):")
    w = winners[(winners["variant"] == "BD_d2") &
                (winners["lambda_bd"] == lbd_check) &
                (winners["gamma_geo"] == 0.0)]
    for _, row in w.sort_values("n").iterrows():
        marker = " ★" if row["winner"] == "4d" else ""
        print(f"    N={int(row['n']):3d}: winner={row['winner']:>4}, margin={row['margin']:.3f}{marker}")

# Dimension order at critical λ
print("\n--- Full ranking at λ=6, γ_geo=0, BD_d2 ---")
for n in N_VALUES:
    sub = summary[(summary["variant"] == "BD_d2") &
                  (summary["lambda_bd"] == 6.0) &
                  (summary["gamma_geo"] == 0.0) &
                  (summary["n"] == n)].sort_values("mean_score")
    ranking = " > ".join(
        f"{r['family'].replace('lorentzian_like_', '')}({r['mean_score']:.1f})"
        for _, r in sub.iterrows()
    )
    print(f"  N={n}: {ranking}")

# ===== Analysis 4: BD_d2 vs BD_d4 — Which is more selective? =====
print("\n" + "=" * 72)
print("ANALYSIS 4: BD_d2 vs BD_d4 SELECTIVITY COMPARISON")
print("=" * 72)

for vname in ["BD_d2", "BD_d4"]:
    w = winners[(winners["variant"] == vname) & (winners["gamma_geo"] == 0.0)]
    print(f"\n  {vname} — λ range where 4D wins ALL N:")
    for lbd in LAMBDA_BD:
        w_l = w[w["lambda_bd"] == lbd]
        n4d = (w_l["winner"] == "4d").sum()
        avg_m = w_l[w_l["winner"] == "4d"]["margin"].mean() if n4d > 0 else 0
        dominant = w_l["winner"].value_counts().index[0]
        print(f"    λ={lbd:5.1f}: 4D wins {n4d}/{len(N_VALUES)}, "
              f"dominant={dominant}, avg_margin={avg_m:.3f}")

print(f"\nAll outputs saved to {OUT_DIR}")
