"""
Prediction A — Benincasa-Dowker Action Dimension Selection

Tests whether the BD causal set action naturally selects a preferred
spacetime dimension by penalizing causally sparse (high-D) structures.

BD action (d=2 form): S_BD = N - 2 * n_links
where n_links = Hasse covering relations.

Score = -β * log_H + λ * S_BD / N

At λ=0: pure entropy → 5D wins (known).
As λ↑: BD penalty kicks in → selects lower-D.
Key question: is there a λ window where 3+1 (4D) or 2+1 (3D)
is the natural winner?
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd

from experiment import FAMILIES
from runtime_utils import estimate_entropy_by_family
from observables_geo import cover_density

OUT_DIR = Path("outputs_exploratory/prediction_a_bd_dimension")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# --- Configuration ---
N_VALUES = [20, 28, 36, 44, 52]
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

# λ scan: from entropy-dominated to BD-dominated
LAMBDA_VALUES = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 18.0, 25.0,
                 35.0, 50.0, 75.0, 100.0]


def hasse_links(poset) -> int:
    """Count covering relations (Hasse diagram edges)."""
    c = poset.closure.astype(np.uint8, copy=False)
    has_intermediate = (c @ c).astype(bool, copy=False)
    cover = poset.closure & ~has_intermediate
    np.fill_diagonal(cover, False)
    return int(cover.sum())


def bd_action_d2(poset) -> float:
    """BD action (d=2): S = N - 2 * n_links."""
    return float(poset.n - 2 * hasse_links(poset))


def bd_action_normalized(poset) -> float:
    """BD action per element: S/N = 1 - 2*n_links/N."""
    return bd_action_d2(poset) / max(poset.n, 1)


# --- Data collection ---
print("=" * 70)
print("BENINCASA-DOWKER ACTION DIMENSION SELECTION")
print("=" * 70)
print(f"\nN values: {N_VALUES}")
print(f"Families: {FAMILY_LIST}")
print(f"λ range:  {LAMBDA_VALUES[0]} .. {LAMBDA_VALUES[-1]}")
print(f"Samples:  {SAMPLES}, SIS runs: {SIS_RUNS}\n")

# Phase 1: Collect base observables (logH, BD action) per sample
print("Phase 1: Collecting observables...")
base_rows = []
for n in N_VALUES:
    print(f"  N={n}", end="", flush=True)
    for family in FAMILY_LIST:
        generator = FAMILIES[family]
        for sid in range(SAMPLES):
            seed = 950000 + 1000 * n + sid
            poset = generator(n=n, seed=seed)
            log_h, method = estimate_entropy_by_family(
                poset, family=family, sis_runs=SIS_RUNS, seed=seed,
                default_exact_threshold=DEFAULT_EXACT_THRESHOLD,
                family_exact_thresholds=FAMILY_EXACT_THRESHOLDS,
            )
            n_links = hasse_links(poset)
            s_bd = bd_action_d2(poset)
            s_bd_norm = s_bd / max(n, 1)
            cd = cover_density(poset)
            n_relations = int(poset.closure.sum())
            rel_density = n_relations / (n * (n - 1) / 2) if n > 1 else 0.0

            base_rows.append({
                "n": n, "family": family, "sample_id": sid, "seed": seed,
                "log_H": float(log_h), "method": method,
                "n_links": n_links, "n_relations": n_relations,
                "rel_density": rel_density, "cover_density": cd,
                "S_BD": s_bd, "S_BD_norm": s_bd_norm,
            })
        print(f" [{family[:5]}]", end="", flush=True)
    print()

base = pd.DataFrame(base_rows)
base.to_csv(OUT_DIR / "base_observables.csv", index=False)

# Print BD action summary
print("\n--- BD Action & Entropy per family ---")
for n in N_VALUES:
    print(f"\n  N={n}:")
    for fam in FAMILY_LIST:
        sub = base[(base["n"] == n) & (base["family"] == fam)]
        print(f"    {fam:24s}: links={sub['n_links'].mean():6.1f}  "
              f"S_BD/N={sub['S_BD_norm'].mean():+.3f}  "
              f"logH={sub['log_H'].mean():7.2f}  "
              f"rel_dens={sub['rel_density'].mean():.3f}")

# Phase 2: Score under λ scan
print("\nPhase 2: Scoring under λ scan...")
score_rows = []
for _, row in base.iterrows():
    for lam in LAMBDA_VALUES:
        # Pure BD action: score = -β*logH + λ*S_BD/N
        score = -BETA * row["log_H"] + lam * row["S_BD_norm"]
        score_rows.append({
            "n": row["n"], "family": row["family"],
            "sample_id": row["sample_id"],
            "lambda": lam, "score": score,
            "log_H": row["log_H"], "S_BD_norm": row["S_BD_norm"],
        })

scores = pd.DataFrame(score_rows)

# Phase 3: Winner analysis
print("Phase 3: Winner analysis...")
summary = scores.groupby(["n", "lambda", "family"]).agg(
    mean_score=("score", "mean"),
    std_score=("score", "std"),
).reset_index()

winner_rows = []
for (n, lam), grp in summary.groupby(["n", "lambda"]):
    ordered = grp.sort_values("mean_score")
    best = ordered.iloc[0]
    second = ordered.iloc[1] if len(ordered) > 1 else ordered.iloc[0]
    margin = float(second["mean_score"] - best["mean_score"])
    winner_rows.append({
        "n": int(n), "lambda": lam,
        "winner": best["family"],
        "winner_score": best["mean_score"],
        "margin": margin,
        "runner_up": second["family"],
    })

winners = pd.DataFrame(winner_rows)
winners.to_csv(OUT_DIR / "winners_by_lambda.csv", index=False)

# Phase 4: Print results
print("\n" + "=" * 70)
print("BD ACTION DIMENSION SELECTION — RESULTS")
print("=" * 70)

# Winner table
print("\n--- Winner by (N, λ) ---")
pivot = winners.pivot(index="lambda", columns="n", values="winner")
print("\n  λ\\N  ", end="")
for n in N_VALUES:
    print(f"  {n:>12}", end="")
print()
for lam in LAMBDA_VALUES:
    print(f"  {lam:5.1f}", end="")
    for n in N_VALUES:
        row = winners[(winners["n"] == n) & (winners["lambda"] == lam)]
        if not row.empty:
            w = str(row.iloc[0]["winner"]).replace("lorentzian_like_", "").replace("KR_like", "KR")
            print(f"  {w:>12}", end="")
        else:
            print(f"  {'?':>12}", end="")
    print()

# Find λ transition points
print("\n--- λ Transition Points ---")
for n in N_VALUES:
    w_n = winners[winners["n"] == n].sort_values("lambda")
    transitions = []
    prev_winner = None
    for _, row in w_n.iterrows():
        curr = row["winner"]
        if prev_winner is not None and curr != prev_winner:
            transitions.append(f"λ={row['lambda']:.1f}: {prev_winner.replace('lorentzian_like_', '')} → {curr.replace('lorentzian_like_', '')}")
        prev_winner = curr
    if transitions:
        print(f"  N={n}: {', '.join(transitions)}")
    else:
        print(f"  N={n}: no transitions (always {w_n.iloc[0]['winner'].replace('lorentzian_like_', '')})")

# Margin analysis: at each λ, average margin between winner and runner-up
print("\n--- Mean margin by λ ---")
for lam in LAMBDA_VALUES:
    w_lam = winners[winners["lambda"] == lam]
    print(f"  λ={lam:5.1f}: mean_margin={w_lam['margin'].mean():.3f}  "
          f"winners={dict(w_lam['winner'].str.replace('lorentzian_like_', '').str.replace('KR_like', 'KR').value_counts())}")

# Key question: is there a λ* where 3D or 4D wins across all N?
print("\n--- Dimensional preference sweep ---")
for target_dim in ["lorentzian_like_3d", "lorentzian_like_4d"]:
    dim_label = target_dim.replace("lorentzian_like_", "")
    for lam in LAMBDA_VALUES:
        w_lam = winners[winners["lambda"] == lam]
        n_wins = (w_lam["winner"] == target_dim).sum()
        if n_wins > 0:
            winning_ns = list(w_lam[w_lam["winner"] == target_dim]["n"])
            print(f"  {dim_label} wins at λ={lam:5.1f}: {n_wins}/{len(N_VALUES)} N values ({winning_ns})")

# Phase 5: Score breakdown table for key λ values
print("\n--- Score breakdown at key λ values ---")
key_lambdas = [0.0, 5.0, 12.0, 25.0, 50.0]
for lam in key_lambdas:
    print(f"\n  λ = {lam}")
    sub = summary[summary["lambda"] == lam]
    for n in N_VALUES:
        print(f"    N={n}:")
        s_n = sub[sub["n"] == n].sort_values("mean_score")
        for _, row in s_n.iterrows():
            fam = row["family"].replace("lorentzian_like_", "").replace("KR_like", "KR")
            marker = " ★" if row["mean_score"] == s_n.iloc[0]["mean_score"] else ""
            print(f"      {fam:>5}: score={row['mean_score']:9.2f}{marker}")

print(f"\nAll outputs saved to {OUT_DIR}")
