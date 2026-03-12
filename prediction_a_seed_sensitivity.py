from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from action import action_value
from experiment import FAMILIES
from observables import neutral_penalty
from observables_geo import geometric_components
from prediction_a_geometric_ablation import geometric_penalty_from_components
from runtime_utils import estimate_entropy_by_family


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["seed_label", "n", "gamma", "variant", "family"], sort=True)
        .agg(
            mean_score=("score", "mean"),
            std_score=("score", "std"),
            mean_log_H=("log_H", "mean"),
            mean_penalty=("penalty_effective", "mean"),
            count=("score", "count"),
        )
        .reset_index()
        .sort_values(["seed_label", "variant", "n", "gamma", "mean_score"])
    )


def winner_summary(summary_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (seed_label, n, gamma, variant), group in summary_df.groupby(
        ["seed_label", "n", "gamma", "variant"], sort=True
    ):
        ordered = group.sort_values("mean_score")
        best = ordered.iloc[0]
        second = ordered.iloc[1] if len(ordered) > 1 else None
        rows.append(
            {
                "seed_label": str(seed_label),
                "n": int(n),
                "gamma": float(gamma),
                "variant": str(variant),
                "winner_family": str(best["family"]),
                "winner_score": float(best["mean_score"]),
                "runnerup_family": None if second is None else str(second["family"]),
                "runnerup_score": None if second is None else float(second["mean_score"]),
                "winner_margin": None if second is None else float(second["mean_score"] - best["mean_score"]),
            }
        )
    return pd.DataFrame(rows)


def pairwise_summary(summary_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    targets = [("lorentzian_like_4d", "lorentzian_like_2d"), ("lorentzian_like_4d", "lorentzian_like_3d")]
    for (seed_label, n, gamma, variant), group in summary_df.groupby(
        ["seed_label", "n", "gamma", "variant"], sort=True
    ):
        lookup = {str(row["family"]): row for _, row in group.iterrows()}
        for left, right in targets:
            if left not in lookup or right not in lookup:
                continue
            rows.append(
                {
                    "seed_label": str(seed_label),
                    "n": int(n),
                    "gamma": float(gamma),
                    "variant": str(variant),
                    "left_family": left,
                    "right_family": right,
                    "delta_score": float(lookup[left]["mean_score"] - lookup[right]["mean_score"]),
                    "delta_log_H": float(lookup[left]["mean_log_H"] - lookup[right]["mean_log_H"]),
                    "delta_penalty": float(lookup[left]["mean_penalty"] - lookup[right]["mean_penalty"]),
                    "left_wins": bool(float(lookup[left]["mean_score"]) < float(lookup[right]["mean_score"])),
                }
            )
    return pd.DataFrame(rows)


def seed_sensitivity_summary(winners_df: pd.DataFrame, pairwise_df: pd.DataFrame) -> pd.DataFrame:
    winner_counts = (
        winners_df.groupby(["seed_label", "variant", "winner_family"], sort=True)
        .size()
        .rename("wins")
        .reset_index()
    )
    pairwise_stats = (
        pairwise_df.groupby(["variant", "right_family"], sort=True)
        .agg(
            mean_margin=("delta_score", lambda s: float((-s).mean())),
            min_margin=("delta_score", lambda s: float((-s).min())),
            max_margin=("delta_score", lambda s: float((-s).max())),
            lor4d_win_rate=("left_wins", "mean"),
        )
        .reset_index()
    )
    return winner_counts, pairwise_stats


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prediction A seed sensitivity under mixed near-wall settings.")
    parser.add_argument(
        "--config",
        default="config_prediction_a_seed_sensitivity_n64.yaml",
        help="Path to YAML config.",
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    exp_cfg = config["experiment"]
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    n_values = [int(x) for x in exp_cfg["n_values"]]
    families = [str(x) for x in exp_cfg["families"]]
    gammas = [float(x) for x in exp_cfg["gammas"]]
    variants = [str(x) for x in exp_cfg["variants"]]
    samples_per_family = int(exp_cfg["samples_per_family"])
    beta = float(exp_cfg.get("beta", 1.0))
    sis_runs = int(exp_cfg.get("sis_runs", 4096))
    default_exact_threshold = int(exp_cfg.get("exact_threshold", 32))
    consistency_weight = float(exp_cfg.get("consistency_weight", 0.682882))
    multi_consistency_weight = float(exp_cfg.get("multi_consistency_weight", 0.682882))
    seed_base = int(exp_cfg.get("seed_base", 800000))
    seed_offsets = [int(x) for x in exp_cfg.get("seed_offsets", [0])]
    family_exact_thresholds = {
        str(k): int(v) for k, v in (exp_cfg.get("family_exact_thresholds", {}) or {}).items()
    }

    rows = []
    for seed_offset in seed_offsets:
        seed_label = f"offset_{seed_offset}"
        for n in n_values:
            for family in families:
                generator = FAMILIES[family]
                for sample_id in range(samples_per_family):
                    seed = seed_base + seed_offset + 1000 * n + sample_id
                    poset = generator(n=n, seed=seed)
                    log_h, entropy_method = estimate_entropy_by_family(
                        poset,
                        family=family,
                        sis_runs=sis_runs,
                        seed=seed,
                        default_exact_threshold=default_exact_threshold,
                        family_exact_thresholds=family_exact_thresholds,
                    )
                    neutral = neutral_penalty(poset)
                    geo = geometric_components(poset)
                    for gamma in gammas:
                        for variant in variants:
                            penalty_effective = neutral + geometric_penalty_from_components(
                                geo,
                                variant,
                                consistency_weight=consistency_weight,
                                multi_consistency_weight=multi_consistency_weight,
                            )
                            rows.append(
                                {
                                    "seed_label": seed_label,
                                    "seed_offset": seed_offset,
                                    "n": n,
                                    "family": family,
                                    "sample_id": sample_id,
                                    "seed": seed,
                                    "gamma": gamma,
                                    "beta": beta,
                                    "variant": variant,
                                    "log_H": float(log_h),
                                    "entropy_method": entropy_method,
                                    "penalty_neutral": float(neutral),
                                    "penalty_effective": float(penalty_effective),
                                    "score": action_value(float(log_h), float(penalty_effective), beta=beta, gamma=gamma),
                                    **geo,
                                }
                            )

    raw_df = pd.DataFrame(rows)
    summary_df = summarize(raw_df)
    winners_df = winner_summary(summary_df)
    pairwise_df = pairwise_summary(summary_df)
    winner_counts_df, pairwise_stats_df = seed_sensitivity_summary(winners_df, pairwise_df)

    raw_df.to_csv(out_dir / "prediction_a_seed_sensitivity_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "prediction_a_seed_sensitivity_summary.csv", index=False, encoding="utf-8-sig")
    winners_df.to_csv(out_dir / "prediction_a_seed_sensitivity_winners.csv", index=False, encoding="utf-8-sig")
    pairwise_df.to_csv(out_dir / "prediction_a_seed_sensitivity_pairwise.csv", index=False, encoding="utf-8-sig")
    winner_counts_df.to_csv(out_dir / "prediction_a_seed_sensitivity_winner_counts.csv", index=False, encoding="utf-8-sig")
    pairwise_stats_df.to_csv(out_dir / "prediction_a_seed_sensitivity_pairwise_stats.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "prediction_a_seed_sensitivity_winner_counts.csv").as_posix())
    print()
    print(winner_counts_df.to_string(index=False))
