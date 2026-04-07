"""信息论惩罚消融实验 (Ablation Study)

逐项剔除 observables_info 的 5 项信息论惩罚，观察哪项对 Lor4D 选择最关键。
每轮实验保留 4/5 项（将一项权重置零），运行 A4 模态排名竞赛。

输出: ablation_results.csv  — 每行一个 (dropped_term, n, γ, family, rank, score)
"""

from __future__ import annotations

import itertools
from copy import deepcopy
from pathlib import Path

import pandas as pd
import numpy as np

from experiment import run_experiment, FAMILIES
from normalization import add_normalized_columns, add_size_scaled_columns
from observables_info import DEFAULT_INFO_WEIGHTS


def run_ablation(
    n_values: tuple[int, ...] = (10, 20, 40),
    gammas: tuple[float, ...] = (0.0, 0.2, 0.5, 1.0),
    beta: float = 1.0,
    samples_per_family: int = 8,
    sis_runs: int = 32,
    output_dir: str = "outputs_info_ablation",
) -> pd.DataFrame:
    """Run ablation study: one experiment per dropped term + baseline (all terms)."""

    terms = list(DEFAULT_INFO_WEIGHTS.keys())
    all_results: list[pd.DataFrame] = []

    # Baseline: all terms active
    print("=== Baseline (all terms) ===")
    df_base = run_experiment(
        n_values=n_values,
        gammas=gammas,
        beta=beta,
        samples_per_family=samples_per_family,
        sis_runs=sis_runs,
        action_modes=("A4",),
    )
    df_base = add_size_scaled_columns(df_base)
    df_base = add_normalized_columns(df_base, method="robust_zscore", group_cols=("n", "gamma", "action_mode"))
    df_base["dropped_term"] = "none"
    all_results.append(df_base)

    # Ablation: drop one term at a time
    for drop_term in terms:
        print(f"=== Ablation: drop {drop_term} ===")
        ablated_weights = deepcopy(DEFAULT_INFO_WEIGHTS)
        ablated_weights[drop_term] = 0.0

        # We need to override the weight in experiment — but get_action_penalty
        # doesn't pass info_weights yet in the experiment loop. Instead, we
        # temporarily monkey-patch DEFAULT_INFO_WEIGHTS.
        import observables_info
        original = observables_info.DEFAULT_INFO_WEIGHTS.copy()
        observables_info.DEFAULT_INFO_WEIGHTS[drop_term] = 0.0

        try:
            df_abl = run_experiment(
                n_values=n_values,
                gammas=gammas,
                beta=beta,
                samples_per_family=samples_per_family,
                sis_runs=sis_runs,
                action_modes=("A4",),
            )
            df_abl = add_size_scaled_columns(df_abl)
            df_abl = add_normalized_columns(df_abl, method="robust_zscore", group_cols=("n", "gamma", "action_mode"))
            df_abl["dropped_term"] = drop_term
            all_results.append(df_abl)
        finally:
            observables_info.DEFAULT_INFO_WEIGHTS.update(original)

    combined = pd.concat(all_results, ignore_index=True)

    # Summary: Lor4D rank per (dropped_term, n, γ)
    summary_rows = []
    for (drop, n, gamma), group in combined.groupby(["dropped_term", "n", "gamma"]):
        ranked = group.groupby("family")["score_norm"].mean().sort_values()
        families_ranked = list(ranked.index)
        lor4d_rank = families_ranked.index("lorentzian_like_4d") + 1 if "lorentzian_like_4d" in families_ranked else -1
        lor4d_score = ranked.get("lorentzian_like_4d", np.nan)
        top1 = families_ranked[0]
        top1_score = ranked.iloc[0]
        summary_rows.append({
            "dropped_term": drop,
            "n": n,
            "gamma": gamma,
            "lor4d_rank": lor4d_rank,
            "lor4d_score": lor4d_score,
            "top1_family": top1,
            "top1_score": top1_score,
            "total_families": len(families_ranked),
        })

    summary = pd.DataFrame(summary_rows)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out_dir / "ablation_raw.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(out_dir / "ablation_summary.csv", index=False, encoding="utf-8-sig")

    print("\n=== Ablation Summary: Lor4D rank by dropped term (n=20, γ=1.0) ===")
    pivot = summary[(summary["n"] == n_values[-1]) & (summary["gamma"] == gammas[-1])]
    print(pivot[["dropped_term", "lor4d_rank", "lor4d_score", "top1_family"]].to_string(index=False))

    return summary


if __name__ == "__main__":
    run_ablation()
