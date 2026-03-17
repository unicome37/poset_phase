"""
Prediction B — Near-wall γ-extended scan with SIS
Tests N=48,52,56 with γ up to 5.0 using mixed method
(lor2d exact + KR SIS).
"""
import numpy as np
import pandas as pd
from pathlib import Path

from generators import generate_lorentzian_like_2d, generate_kr_like
from entropy_exact import log_linear_extensions_exact
from entropy_sis import estimate_log_linear_extensions_sis
from observables import neutral_penalty
from observables_geo import geometric_penalty
from action import action_value, get_action_penalty

OUT_DIR = Path("outputs_exploratory/prediction_b_gamma_extended_nearwall")
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_VALUES = [48, 52, 56]
GAMMAS = [0.0, 0.2, 0.4, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
SAMPLES = 4
SIS_RUNS = 4096
BETA = 1.0

rows = []

for n in N_VALUES:
    print(f"\n=== N = {n} ===")
    for family_name, generator, use_exact in [
        ("lorentzian_like_2d", generate_lorentzian_like_2d, True),
        ("KR_like", generate_kr_like, False),
    ]:
        for sid in range(SAMPLES):
            seed = 1000 * n + sid
            poset = generator(n=n, seed=seed)

            if use_exact:
                log_h = log_linear_extensions_exact(poset)
                log_h_std = 0.0
                method = "exact"
            else:
                log_h, log_h_std = estimate_log_linear_extensions_sis(
                    poset, n_runs=SIS_RUNS, seed=seed)
                method = "sis"

            penalty = get_action_penalty(poset, "A2")

            for g in GAMMAS:
                score = action_value(log_extensions=log_h, penalty=penalty,
                                     beta=BETA, gamma=g)
                rows.append({
                    "n": n, "family": family_name, "sample_id": sid,
                    "gamma": g, "log_H": log_h, "log_H_std": log_h_std,
                    "penalty": penalty, "score": score, "method": method,
                })

            print(f"  {family_name} sid={sid} ({method}): logH={log_h:.2f}, penalty={penalty:.4f}")

df = pd.DataFrame(rows)
df.to_csv(OUT_DIR / "raw_samples.csv", index=False)

# Summarize and find γ_c
summary = df.groupby(["n", "gamma", "family"]).agg(
    mean_score=("score", "mean"),
    std_score=("score", "std"),
    count=("score", "count"),
).reset_index()

summary.to_csv(OUT_DIR / "summary.csv", index=False)

print("\n\n=== Near-wall γ_c extraction ===")
print(f"{'N':>4} | {'γ_c':>8} | {'Winner at γ=5':>15}")
print("-" * 40)

gc_results = []

for n in N_VALUES:
    sub = summary[summary["n"] == n]
    gammas_arr = sorted(sub["gamma"].unique())
    diffs = []
    for g in gammas_arr:
        sg = sub[sub["gamma"] == g]
        lor = sg[sg["family"] == "lorentzian_like_2d"]["mean_score"].values
        kr = sg[sg["family"] == "KR_like"]["mean_score"].values
        if len(lor) > 0 and len(kr) > 0:
            diffs.append(lor[0] - kr[0])
        else:
            diffs.append(np.nan)

    gc = None
    gammas_arr = np.array(gammas_arr)
    diffs = np.array(diffs)
    for i in range(len(diffs) - 1):
        if not np.isnan(diffs[i]) and not np.isnan(diffs[i+1]):
            if diffs[i] > 0 and diffs[i+1] <= 0:
                gc = gammas_arr[i] + (gammas_arr[i+1] - gammas_arr[i]) * diffs[i] / (diffs[i] - diffs[i+1])
                break

    winner_5 = "?"
    for g_val in [5.0]:
        sg = sub[sub["gamma"] == g_val]
        lor = sg[sg["family"] == "lorentzian_like_2d"]["mean_score"].values
        kr = sg[sg["family"] == "KR_like"]["mean_score"].values
        if len(lor) > 0 and len(kr) > 0:
            winner_5 = "lor2d" if lor[0] > kr[0] else "KR"

    gc_str = f"{gc:.4f}" if gc is not None else "no cross"
    print(f"{n:4d} | {gc_str:>8s} | {winner_5:>15s}")
    gc_results.append({"n": n, "gamma_c": gc, "winner_at_5": winner_5})

pd.DataFrame(gc_results).to_csv(OUT_DIR / "nearwall_gamma_c.csv", index=False)
print(f"\nSaved to {OUT_DIR}")
