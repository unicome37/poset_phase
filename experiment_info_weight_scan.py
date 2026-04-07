"""信息论惩罚权重敏感性扫描实验

系统扫描不同权重配比下 Lor4D 的排名变化：
1. interval_diversity 权重从 0 到 10 递增，其余固定 5.0
2. 所有项统一缩放 (0.5x, 1x, 2x, 3x)
3. "去掉 interval + 放大其余" vs "降低 interval"
4. 找到使 Lor4D 在宽 γ 范围内保持 top-3 的最优配比

输出: outputs_info_weight_scan/weight_scan_summary.csv
"""
from __future__ import annotations

import itertools
from copy import deepcopy
from pathlib import Path

import pandas as pd
import numpy as np

from experiment import run_experiment, FAMILIES
from normalization import add_normalized_columns, add_size_scaled_columns
import observables_info


BASELINE_WEIGHTS = {
    "info_spectral_entropy_deficit": 5.0,
    "info_degree_heterogeneity": 5.0,
    "info_layer_concentration": 5.0,
    "info_edge_density_extremity": 5.0,
    "info_interval_diversity_deficit": 5.0,
}

# --- Scan configurations ---

def make_scan_configs() -> list[tuple[str, dict[str, float]]]:
    """Generate named weight configurations to scan."""
    configs = []

    # 1) Baseline
    configs.append(("baseline_5x5", BASELINE_WEIGHTS.copy()))

    # 2) Interval weight sweep: 0, 0.5, 1, 2, 3, 5, 8, 10
    for iw in [0.0, 0.5, 1.0, 2.0, 3.0, 8.0, 10.0]:
        w = BASELINE_WEIGHTS.copy()
        w["info_interval_diversity_deficit"] = iw
        configs.append((f"interval_{iw:.1f}", w))

    # 3) Uniform scaling: all weights multiplied by factor
    for factor in [0.5, 2.0, 3.0]:
        w = {k: v * factor for k, v in BASELINE_WEIGHTS.items()}
        configs.append((f"uniform_{factor:.1f}x", w))

    # 4) Remove interval + boost others
    for boost in [1.0, 1.5, 2.0, 3.0]:
        w = {k: v * boost for k, v in BASELINE_WEIGHTS.items()}
        w["info_interval_diversity_deficit"] = 0.0
        configs.append((f"no_interval_boost_{boost:.1f}x", w))

    # 5) Optimal candidates: low interval + moderate boost
    for iw, boost in [(1.0, 1.5), (2.0, 1.5), (1.0, 2.0), (0.5, 2.0)]:
        w = {k: v * boost for k, v in BASELINE_WEIGHTS.items()}
        w["info_interval_diversity_deficit"] = iw
        configs.append((f"interval_{iw:.1f}_rest_{boost:.1f}x", w))

    return configs


def run_weight_scan(
    n_values: tuple[int, ...] = (10, 20, 40),
    gammas: tuple[float, ...] = (0.0, 0.2, 0.5, 1.0, 2.0),
    beta: float = 1.0,
    samples_per_family: int = 10,
    sis_runs: int = 48,
    output_dir: str = "outputs_info_weight_scan",
) -> pd.DataFrame:
    """Run weight sensitivity scan."""

    configs = make_scan_configs()
    all_summaries = []

    for config_name, weights in configs:
        print(f"\n{'='*60}")
        print(f"Config: {config_name}")
        print(f"  Weights: { {k.replace('info_',''): v for k,v in weights.items()} }")
        print(f"{'='*60}")

        # Monkey-patch weights
        original = observables_info.DEFAULT_INFO_WEIGHTS.copy()
        observables_info.DEFAULT_INFO_WEIGHTS.update(weights)

        try:
            df = run_experiment(
                n_values=n_values,
                gammas=gammas,
                beta=beta,
                samples_per_family=samples_per_family,
                sis_runs=sis_runs,
                action_modes=("A4",),
            )
            df = add_size_scaled_columns(df)
            df = add_normalized_columns(
                df, method="robust_zscore",
                group_cols=("n", "gamma", "action_mode"),
            )

            # Extract Lor4D rank per (n, gamma)
            for (n, gamma), group in df.groupby(["n", "gamma"]):
                ranked = group.groupby("family")["score"].mean().sort_values()
                families_ranked = list(ranked.index)
                lor4d_rank = (
                    families_ranked.index("lorentzian_like_4d") + 1
                    if "lorentzian_like_4d" in families_ranked
                    else -1
                )
                lor4d_score = ranked.get("lorentzian_like_4d", np.nan)
                top1 = families_ranked[0]
                top1_score = ranked.iloc[0]
                gap = lor4d_score - top1_score if lor4d_rank > 1 else 0.0

                all_summaries.append({
                    "config": config_name,
                    "n": n,
                    "gamma": gamma,
                    "lor4d_rank": lor4d_rank,
                    "lor4d_raw_score": lor4d_score,
                    "top1_family": top1,
                    "top1_score": top1_score,
                    "gap_to_top1": gap,
                    "interval_weight": weights["info_interval_diversity_deficit"],
                    "total_families": len(families_ranked),
                })

        finally:
            observables_info.DEFAULT_INFO_WEIGHTS.update(original)

    summary = pd.DataFrame(all_summaries)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_dir / "weight_scan_summary.csv", index=False, encoding="utf-8-sig")

    # Print key results
    print("\n" + "=" * 80)
    print("WEIGHT SCAN RESULTS: Lor4D rank at gamma=0.5 (critical regime)")
    print("=" * 80)

    crit = summary[summary["gamma"] == 0.5]
    pivot = crit.pivot_table(
        index="config", columns="n",
        values="lor4d_rank", aggfunc="first",
    )
    pivot["mean_rank"] = pivot.mean(axis=1)
    pivot = pivot.sort_values("mean_rank")
    print(pivot.to_string(float_format=lambda x: f"{x:.1f}"))

    print("\n" + "=" * 80)
    print("WEIGHT SCAN RESULTS: Lor4D rank at gamma=1.0 (strong penalty)")
    print("=" * 80)

    strong = summary[summary["gamma"] == 1.0]
    pivot2 = strong.pivot_table(
        index="config", columns="n",
        values="lor4d_rank", aggfunc="first",
    )
    pivot2["mean_rank"] = pivot2.mean(axis=1)
    pivot2 = pivot2.sort_values("mean_rank")
    print(pivot2.to_string(float_format=lambda x: f"{x:.1f}"))

    # Best config: lowest mean rank across all (n, gamma) combinations
    best = summary.groupby("config")["lor4d_rank"].mean().sort_values()
    print(f"\n=== Overall best configs (mean Lor4D rank across all conditions) ===")
    for config_name, mean_rank in best.head(10).items():
        print(f"  {config_name:40s}  mean_rank={mean_rank:.2f}")

    return summary


if __name__ == "__main__":
    run_weight_scan()
