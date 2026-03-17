"""
Prediction A — 5D Dimension Ceiling Test (Pilot)

Tests whether 5D Lorentzian-like posets continue to dominate over 4D,
or if there is a dimensional saturation effect.

Uses the same framework as prediction_a_geometric_ablation.py but with
5D added and extended pairwise comparisons.
"""
from __future__ import annotations
from pathlib import Path

import pandas as pd

from action import action_value
from experiment import FAMILIES
from observables import neutral_penalty
from observables_geo import geometric_components
from runtime_utils import estimate_entropy_by_family

OUT_DIR = Path("outputs_exploratory/prediction_a_5d_pilot")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# --- Configuration ---
N_VALUES = [20, 24, 28, 32, 36, 40, 44, 48, 52]
FAMILY_LIST = [
    "lorentzian_like_2d",
    "lorentzian_like_3d",
    "lorentzian_like_4d",
    "lorentzian_like_5d",
    "KR_like",
]
GAMMAS = [0.0, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0]
VARIANTS = ["A2_full", "A2_replace_dim_with_consistency"]
SAMPLES = 4
SIS_RUNS = 4096
BETA = 1.0
CONSISTENCY_WEIGHT = 0.682882
MULTI_CONSISTENCY_WEIGHT = 0.682882
DEFAULT_EXACT_THRESHOLD = 24
FAMILY_EXACT_THRESHOLDS = {
    "lorentzian_like_2d": 104,
    "lorentzian_like_3d": 24,
    "lorentzian_like_4d": 24,
    "lorentzian_like_5d": 24,
    "KR_like": 24,
}


def geometric_penalty_from_components(
    geo: dict[str, float], variant: str, consistency_weight: float
) -> float:
    width_height = 2.0 * float(geo["geo_width_height"])
    dim_proxy = 8.0 * float(geo["geo_dim_proxy_penalty"])
    comp_window = 6.0 * float(geo["geo_comparability_window"])
    cover = 3.0 * float(geo["geo_cover_density"])
    interval_profile = 4.0 * float(geo["geo_interval_profile"])
    interval_shape = 4.0 * float(geo["geo_interval_shape"])
    layer_smoothness = 2.0 * float(geo["geo_layer_smoothness"])

    if variant == "A2_full":
        return (width_height + dim_proxy + comp_window + cover
                + interval_profile + interval_shape + layer_smoothness)
    if variant == "A2_replace_dim_with_consistency":
        return (width_height
                + consistency_weight * float(geo["geo_dim_consistency"])
                + comp_window + cover + interval_profile
                + interval_shape + layer_smoothness)
    raise ValueError(f"Unknown variant: {variant}")


# --- Data collection ---
print("Collecting data...")
rows = []
for n in N_VALUES:
    print(f"  N={n}", end="", flush=True)
    for family in FAMILY_LIST:
        generator = FAMILIES[family]
        for sid in range(SAMPLES):
            seed = 900000 + 1000 * n + sid
            poset = generator(n=n, seed=seed)
            log_h, method = estimate_entropy_by_family(
                poset, family=family, sis_runs=SIS_RUNS, seed=seed,
                default_exact_threshold=DEFAULT_EXACT_THRESHOLD,
                family_exact_thresholds=FAMILY_EXACT_THRESHOLDS,
            )
            neutral = neutral_penalty(poset)
            geo = geometric_components(poset)
            for gamma in GAMMAS:
                for variant in VARIANTS:
                    pen_geo = geometric_penalty_from_components(
                        geo, variant, CONSISTENCY_WEIGHT)
                    pen = neutral + pen_geo
                    score = action_value(float(log_h), float(pen),
                                         beta=BETA, gamma=gamma)
                    rows.append({
                        "n": n, "family": family, "sample_id": sid,
                        "seed": seed, "gamma": gamma, "variant": variant,
                        "log_H": float(log_h), "method": method,
                        "penalty": float(pen), "score": score,
                    })
        print(f" [{family[:6]}]", end="", flush=True)
    print()

raw = pd.DataFrame(rows)
raw.to_csv(OUT_DIR / "raw.csv", index=False)

# --- Summary ---
summary = raw.groupby(["n", "gamma", "variant", "family"]).agg(
    mean_score=("score", "mean"),
    std_score=("score", "std"),
    count=("score", "count"),
).reset_index()

summary.to_csv(OUT_DIR / "summary.csv", index=False)

# --- Winner table ---
winner_rows = []
for (n, gamma, variant), grp in summary.groupby(["n", "gamma", "variant"]):
    ordered = grp.sort_values("mean_score")
    best = ordered.iloc[0]
    winner_rows.append({
        "n": int(n), "gamma": gamma, "variant": variant,
        "winner": best["family"],
        "winner_score": best["mean_score"],
    })

winners = pd.DataFrame(winner_rows)
winners.to_csv(OUT_DIR / "winners.csv", index=False)

# --- Pairwise: 5D vs {2D, 3D, 4D} ---
PAIRWISE_TARGETS = [
    ("lorentzian_like_5d", "lorentzian_like_2d"),
    ("lorentzian_like_5d", "lorentzian_like_3d"),
    ("lorentzian_like_5d", "lorentzian_like_4d"),
    ("lorentzian_like_4d", "lorentzian_like_2d"),
    ("lorentzian_like_4d", "lorentzian_like_3d"),
]

pair_rows = []
for (n, gamma, variant), grp in summary.groupby(["n", "gamma", "variant"]):
    lookup = {str(r["family"]): r for _, r in grp.iterrows()}
    for left, right in PAIRWISE_TARGETS:
        if left not in lookup or right not in lookup:
            continue
        delta = float(lookup[left]["mean_score"] - lookup[right]["mean_score"])
        pair_rows.append({
            "n": int(n), "gamma": gamma, "variant": variant,
            "left": left, "right": right,
            "delta_score": delta,
            "left_wins": delta < 0,
        })

pairwise = pd.DataFrame(pair_rows)
pairwise.to_csv(OUT_DIR / "pairwise.csv", index=False)

# --- Margin summary across N (for each variant) ---
print("\n=== DIMENSION CEILING ANALYSIS ===\n")

for variant in VARIANTS:
    print(f"--- {variant} ---")
    pw = pairwise[pairwise["variant"] == variant]

    # 5D vs 4D margin per N (averaged over γ)
    m54 = pw[(pw["left"] == "lorentzian_like_5d") &
             (pw["right"] == "lorentzian_like_4d")]
    if not m54.empty:
        margin_by_n = m54.groupby("n").agg(
            mean_delta=("delta_score", "mean"),
            wins=("left_wins", "sum"),
            count=("left_wins", "count"),
        ).reset_index()
        margin_by_n["win_rate"] = margin_by_n["wins"] / margin_by_n["count"]
        print("\n  5D vs 4D:")
        print(f"  {'N':>4} | {'mean_delta':>10} | {'win_rate':>8} | {'wins/count':>10}")
        for _, row in margin_by_n.iterrows():
            delta_sign = "+" if row["mean_delta"] > 0 else ""
            print(f"  {int(row['n']):4d} | {delta_sign}{row['mean_delta']:9.4f} | {row['win_rate']:8.3f} | {int(row['wins'])}/{int(row['count'])}")

        overall_wins = int(margin_by_n["wins"].sum())
        overall_count = int(margin_by_n["count"].sum())
        print(f"\n  Overall 5D vs 4D: {overall_wins}/{overall_count} = {overall_wins/overall_count:.3f}")

    # 4D vs 2D comparison (reference)
    m42 = pw[(pw["left"] == "lorentzian_like_4d") &
             (pw["right"] == "lorentzian_like_2d")]
    if not m42.empty:
        margin_by_n_42 = m42.groupby("n").agg(
            mean_delta=("delta_score", "mean"),
            wins=("left_wins", "sum"),
            count=("left_wins", "count"),
        ).reset_index()
        overall_wins_42 = int(margin_by_n_42["wins"].sum())
        overall_count_42 = int(margin_by_n_42["count"].sum())
        print(f"  Reference 4D vs 2D: {overall_wins_42}/{overall_count_42} = {overall_wins_42/overall_count_42:.3f}")

    print()

# --- Winner count summary ---
print("=== WINNER COUNTS ===\n")
for variant in VARIANTS:
    w = winners[winners["variant"] == variant]
    counts = w["winner"].value_counts()
    total = len(w)
    print(f"{variant}:")
    for fam, cnt in counts.items():
        print(f"  {fam}: {cnt}/{total} ({cnt/total*100:.1f}%)")
    print()

print(f"All outputs saved to {OUT_DIR}")
