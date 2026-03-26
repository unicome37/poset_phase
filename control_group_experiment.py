"""
对照组实验：KR_2layer / KR_4layer vs 原始 7 个家族

目的：回应 Carlip 批评——验证 2-layer 和 4-layer poset 在结构泛函 F[X] 下的
行为，证明原有结论不是对 7 个家族的 cherry-picking。

输出：outputs_control/control_comparison.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from generators import (
    generate_kr_like,
    generate_kr_2layer,
    generate_kr_4layer,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_transitive_percolation,
    generate_interval_order,
    generate_absolute_layered,
)
from entropy_exact import log_linear_extensions_exact
from observables import neutral_penalty
from observables_geo import geometric_penalty
from action import action_value, get_action_penalty

FAMILIES = {
    # Original 7
    "KR_like": generate_kr_like,
    "lorentzian_like_2d": generate_lorentzian_like_2d,
    "lorentzian_like_3d": generate_lorentzian_like_3d,
    "lorentzian_like_4d": generate_lorentzian_like_4d,
    "transitive_percolation": generate_transitive_percolation,
    "interval_order": generate_interval_order,
    "absolute_layered": generate_absolute_layered,
    # New controls
    "KR_2layer": generate_kr_2layer,
    "KR_4layer": generate_kr_4layer,
}

N_VALUES = [10, 12, 14, 16, 18, 20]
SAMPLES = 8
BETA = 1.0
GAMMAS = [0.0, 0.2, 0.4, 0.8, 1.6]


def run():
    rows = []
    total = len(FAMILIES) * len(N_VALUES) * SAMPLES
    done = 0

    for fname, gen_fn in FAMILIES.items():
        for n in N_VALUES:
            for s in range(SAMPLES):
                seed = 1000 * hash(fname) % (2**31) + 100 * n + s
                poset = gen_fn(n, seed=seed % (2**31))
                log_h = log_linear_extensions_exact(poset)
                np_val = neutral_penalty(poset)
                gp_val = geometric_penalty(poset)
                penalty_a2 = get_action_penalty(poset, "A2")
                penalty_a3 = get_action_penalty(poset, "A3")

                for gamma in GAMMAS:
                    score_a2 = action_value(log_h, penalty_a2, BETA, gamma)
                    score_a3 = action_value(log_h, penalty_a3, BETA, gamma)
                    rows.append({
                        "family": fname,
                        "n": n,
                        "sample": s,
                        "gamma": gamma,
                        "log_H": log_h,
                        "neutral_penalty": np_val,
                        "geometric_penalty": gp_val,
                        "score_A2": score_a2,
                        "score_A3": score_a3,
                    })

                done += 1
                if done % 20 == 0:
                    print(f"  [{done}/{total}] {fname} n={n} s={s}", flush=True)

    df = pd.DataFrame(rows)
    out_dir = Path("outputs_control")
    out_dir.mkdir(exist_ok=True)
    df.to_csv(out_dir / "control_comparison.csv", index=False)

    # Summary table
    summary = (
        df.groupby(["family", "n", "gamma"])
        .agg(
            log_H_mean=("log_H", "mean"),
            log_H_std=("log_H", "std"),
            geo_penalty_mean=("geometric_penalty", "mean"),
            score_A3_mean=("score_A3", "mean"),
            score_A3_std=("score_A3", "std"),
        )
        .reset_index()
    )
    summary.to_csv(out_dir / "control_summary.csv", index=False)

    # Print key comparisons
    print("\n" + "=" * 80)
    print("CONTROL GROUP COMPARISON: KR_2layer / KR_4layer vs Original Families")
    print("=" * 80)

    for gamma in [0.4, 0.8]:
        print(f"\n--- gamma = {gamma} ---")
        sub = summary[summary["gamma"] == gamma].copy()
        for n in N_VALUES:
            sn = sub[sub["n"] == n].sort_values("score_A3_mean")
            print(f"\n  N = {n}:")
            for _, row in sn.iterrows():
                tag = " ★ CONTROL" if row["family"] in ("KR_2layer", "KR_4layer") else ""
                print(
                    f"    {row['family']:30s}  score_A3={row['score_A3_mean']:+8.3f} "
                    f"±{row['score_A3_std']:.3f}  logH={row['log_H_mean']:.2f}{tag}"
                )

    print(f"\nResults saved to {out_dir}/")


if __name__ == "__main__":
    run()
