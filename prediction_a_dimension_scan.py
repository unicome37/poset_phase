from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from action import action_value
from experiment import FAMILIES
from observables import neutral_penalty
from observables_geo import geometric_components
from runtime_utils import estimate_entropy_by_family


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def effective_penalty(geo: dict[str, float], variant: str, consistency_weight: float) -> float:
    neutral = float(geo["penalty_neutral"])
    if variant == "A2_full":
        return neutral + float(geo["geo_total"])
    if variant == "width_height_dim_consistency":
        return neutral + 2.0 * float(geo["geo_width_height"]) + consistency_weight * float(geo["geo_dim_consistency"])
    if variant == "width_height_dim_multi_consistency":
        return neutral + 2.0 * float(geo["geo_width_height"]) + consistency_weight * float(geo["geo_dim_multi_consistency"])
    if variant == "A2_no_comp_window":
        return neutral + (
            float(geo["geo_total"]) - 6.0 * float(geo["geo_comparability_window"])
        )
    raise ValueError(f"Unknown variant: {variant}")


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["n", "gamma", "variant", "family"], sort=True)
        .agg(
            mean_score=("score", "mean"),
            std_score=("score", "std"),
            mean_log_H=("log_H", "mean"),
            mean_penalty=("penalty_effective", "mean"),
            mean_dim_order=("geo_d_order", "mean"),
            mean_dim_chain=("geo_d_chain", "mean"),
            mean_dim_width=("geo_d_width", "mean"),
            mean_dim_consistency=("geo_dim_consistency", "mean"),
            mean_dim_multi_consistency=("geo_dim_multi_consistency", "mean"),
            mean_comp_window=("geo_comparability_window", "mean"),
            count=("score", "count"),
        )
        .reset_index()
        .sort_values(["variant", "n", "gamma", "mean_score"])
    )


def winner_summary(summary_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (n, gamma, variant), group in summary_df.groupby(["n", "gamma", "variant"], sort=True):
        ordered = group.sort_values("mean_score")
        best = ordered.iloc[0]
        second = ordered.iloc[1] if len(ordered) > 1 else None
        rows.append(
            {
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
    for (n, gamma, variant), group in summary_df.groupby(["n", "gamma", "variant"], sort=True):
        lookup = {str(row["family"]): row for _, row in group.iterrows()}
        for left, right in targets:
            if left not in lookup or right not in lookup:
                continue
            rows.append(
                {
                    "n": int(n),
                    "gamma": float(gamma),
                    "variant": str(variant),
                    "left_family": left,
                    "right_family": right,
                    "delta_score": float(lookup[left]["mean_score"] - lookup[right]["mean_score"]),
                    "delta_log_H": float(lookup[left]["mean_log_H"] - lookup[right]["mean_log_H"]),
                    "delta_penalty": float(lookup[left]["mean_penalty"] - lookup[right]["mean_penalty"]),
                    "delta_dim_consistency": float(
                        lookup[left]["mean_dim_consistency"] - lookup[right]["mean_dim_consistency"]
                    ),
                    "delta_comp_window": float(
                        lookup[left]["mean_comp_window"] - lookup[right]["mean_comp_window"]
                    ),
                    "left_wins": bool(float(lookup[left]["mean_score"]) < float(lookup[right]["mean_score"])),
                }
            )
    return pd.DataFrame(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prediction A pilot: cross-dimensional Lorentzian scan.")
    parser.add_argument("--config", default="config_prediction_a_pilot.yaml", help="Path to YAML config.")
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
    sis_runs = int(exp_cfg.get("sis_runs", 1024))
    default_exact_threshold = int(exp_cfg.get("exact_threshold", 40))
    family_exact_thresholds = {
        str(k): int(v) for k, v in (exp_cfg.get("family_exact_thresholds", {}) or {}).items()
    }
    consistency_weight = float(exp_cfg.get("consistency_weight", 1.0))

    rows = []
    for n in n_values:
        for family in families:
            generator = FAMILIES[family]
            for sample_id in range(samples_per_family):
                seed = 700000 + 1000 * n + sample_id
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
                geo["penalty_neutral"] = neutral
                for gamma in gammas:
                    for variant in variants:
                        penalty_effective = effective_penalty(geo, variant, consistency_weight=consistency_weight)
                        rows.append(
                            {
                                "n": n,
                                "family": family,
                                "sample_id": sample_id,
                                "seed": seed,
                                "gamma": gamma,
                                "beta": beta,
                                "variant": variant,
                                "log_H": float(log_h),
                                "entropy_method": entropy_method,
                                "penalty_effective": float(penalty_effective),
                                "score": action_value(float(log_h), float(penalty_effective), beta=beta, gamma=gamma),
                                **geo,
                            }
                        )

    raw_df = pd.DataFrame(rows)
    summary_df = summarize(raw_df)
    winners_df = winner_summary(summary_df)
    pairwise_df = pairwise_summary(summary_df)

    raw_df.to_csv(out_dir / "prediction_a_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "prediction_a_summary.csv", index=False, encoding="utf-8-sig")
    winners_df.to_csv(out_dir / "prediction_a_winners.csv", index=False, encoding="utf-8-sig")
    pairwise_df.to_csv(out_dir / "prediction_a_pairwise.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "prediction_a_winners.csv").as_posix())
    print()
    print(winners_df.head(24).to_string(index=False))
