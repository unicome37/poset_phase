"""
Post-analysis: BD action dimension selection — Lorentzian families only.

Reads the base_observables.csv from the BD dimension experiment
and re-analyzes winners excluding KR_like (which is a random partial order,
not a physical spacetime).

Key question: Among physical (Lorentzian) sprinklings, does the BD action
select a preferred dimension?
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np

DATA_DIR = Path("outputs_exploratory/prediction_a_bd_dimension")
OUT_DIR = DATA_DIR  # append to same directory

base = pd.read_csv(DATA_DIR / "base_observables.csv")

# Filter to Lorentzian families only
lor = base[base["family"].str.startswith("lorentzian_like_")]

LAMBDA_VALUES = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0,
                 6.0, 8.0, 10.0, 12.0, 15.0, 18.0, 25.0, 35.0, 50.0]
N_VALUES = sorted(lor["n"].unique())
BETA = 1.0

# Score
score_rows = []
for _, row in lor.iterrows():
    for lam in LAMBDA_VALUES:
        score = -BETA * row["log_H"] + lam * row["S_BD_norm"]
        score_rows.append({
            "n": row["n"], "family": row["family"],
            "sample_id": row["sample_id"],
            "lambda": lam, "score": score,
        })

scores = pd.DataFrame(score_rows)

# Summary
summary = scores.groupby(["n", "lambda", "family"]).agg(
    mean_score=("score", "mean"),
).reset_index()

# Winner analysis
winner_rows = []
for (n, lam), grp in summary.groupby(["n", "lambda"]):
    ordered = grp.sort_values("mean_score")
    best = ordered.iloc[0]
    second = ordered.iloc[1]
    margin = float(second["mean_score"] - best["mean_score"])
    winner_rows.append({
        "n": int(n), "lambda": lam,
        "winner": best["family"].replace("lorentzian_like_", ""),
        "winner_score": best["mean_score"],
        "margin": margin,
        "runner_up": second["family"].replace("lorentzian_like_", ""),
    })

winners = pd.DataFrame(winner_rows)
winners.to_csv(OUT_DIR / "winners_lorentzian_only.csv", index=False)

# Print results
print("=" * 72)
print("BD ACTION — LORENTZIAN FAMILIES ONLY (Excluding KR)")
print("=" * 72)

print("\n--- Winner by (N, λ) ---")
print(f"\n  {'λ':>5}", end="")
for n in N_VALUES:
    print(f"  {n:>6}", end="")
print()
print("  " + "-" * (5 + 8 * len(N_VALUES)))

for lam in LAMBDA_VALUES:
    print(f"  {lam:5.1f}", end="")
    for n in N_VALUES:
        row = winners[(winners["n"] == n) & (winners["lambda"] == lam)]
        if not row.empty:
            w = row.iloc[0]["winner"]
            print(f"  {w:>6}", end="")
        else:
            print(f"  {'?':>6}", end="")
    print()

# Transition points
print("\n--- λ Transition Points (Lorentzian only) ---")
for n in N_VALUES:
    w_n = winners[winners["n"] == n].sort_values("lambda")
    transitions = []
    prev_winner = None
    for _, row in w_n.iterrows():
        curr = row["winner"]
        if prev_winner is not None and curr != prev_winner:
            transitions.append(f"λ={row['lambda']:.1f}: {prev_winner}→{curr}")
        prev_winner = curr
    if transitions:
        print(f"  N={n}: {', '.join(transitions)}")
    else:
        print(f"  N={n}: always {w_n.iloc[0]['winner']}")

# Which λ gives consistent 3D or 4D wins?
print("\n--- Dimensional consistency across N ---")
for lam in LAMBDA_VALUES:
    w_lam = winners[winners["lambda"] == lam]
    counts = dict(w_lam["winner"].value_counts())
    dom = w_lam["winner"].value_counts().index[0]
    dom_count = counts[dom]
    print(f"  λ={lam:5.1f}: {counts}  {'← UNANIMOUS '+dom if dom_count==len(N_VALUES) else ''}")

# Detailed score table for critical λ range
print("\n--- Score table at critical λ values ---")
critical_lambdas = [2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0]
for lam in critical_lambdas:
    print(f"\n  λ = {lam}")
    sub = summary[summary["lambda"] == lam]
    for n in N_VALUES:
        s_n = sub[sub["n"] == n].sort_values("mean_score")
        ranking = " > ".join(
            f"{r['family'].replace('lorentzian_like_', '')}({r['mean_score']:.1f})"
            for _, r in s_n.iterrows()
        )
        print(f"    N={n}: {ranking}")

# Cross-N consistency: for each λ, which dim is best at most N values?
print("\n--- Best λ for each target dimension ---")
for target in ["2d", "3d", "4d", "5d"]:
    best_lambda = None
    best_count = 0
    for lam in LAMBDA_VALUES:
        w_lam = winners[winners["lambda"] == lam]
        count = (w_lam["winner"] == target).sum()
        if count > best_count:
            best_count = count
            best_lambda = lam
    if best_count > 0:
        w_at_best = winners[(winners["lambda"] == best_lambda) & (winners["winner"] == target)]
        winning_ns = list(w_at_best["n"].values)
        print(f"  {target}: best at λ={best_lambda:.1f}, wins {best_count}/{len(N_VALUES)} N values ({winning_ns})")
    else:
        print(f"  {target}: never wins at any λ")

# Margin analysis: 4D vs 3D and 4D vs 5D
print("\n--- Pairwise margins (positive = left wins) ---")
for lam in [2.5, 3.0, 3.5, 4.0, 5.0]:
    sub = summary[summary["lambda"] == lam]
    print(f"\n  λ = {lam}")
    for n in N_VALUES:
        s_n = sub[sub["n"] == n].set_index("family")
        s4 = s_n.loc["lorentzian_like_4d", "mean_score"]
        s3 = s_n.loc["lorentzian_like_3d", "mean_score"]
        s5 = s_n.loc["lorentzian_like_5d", "mean_score"]
        s2 = s_n.loc["lorentzian_like_2d", "mean_score"]
        print(f"    N={n}: 4D-3D={s3-s4:+.2f}  4D-5D={s5-s4:+.2f}  "
              f"4D-2D={s2-s4:+.2f}  [4D score={s4:.2f}]")

print(f"\nSaved to {OUT_DIR / 'winners_lorentzian_only.csv'}")
