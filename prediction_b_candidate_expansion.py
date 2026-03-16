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


def geometric_penalty_from_variant(
    geo: dict[str, float],
    variant: str,
    consistency_weight: float,
    multi_consistency_weight: float,
) -> float:
    width_height = 2.0 * float(geo["geo_width_height"])
    dim_proxy = 8.0 * float(geo["geo_dim_proxy_penalty"])
    comp_window = 6.0 * float(geo["geo_comparability_window"])
    cover = 3.0 * float(geo["geo_cover_density"])
    interval_profile = 5.0 * float(geo["geo_interval_profile"])
    interval_shape = 5.0 * float(geo["geo_interval_shape"])
    layer_smoothness = 2.0 * float(geo["geo_layer_smoothness"])

    if variant == "A2_full":
        return width_height + dim_proxy + comp_window + cover + interval_profile + interval_shape + layer_smoothness
    if variant == "A2_replace_dim_with_consistency":
        return width_height + consistency_weight * float(geo["geo_dim_consistency"]) + comp_window + cover + interval_profile + interval_shape + layer_smoothness
    if variant == "A2_replace_dim_with_multi_consistency":
        return width_height + multi_consistency_weight * float(geo["geo_dim_multi_consistency"]) + comp_window + cover + interval_profile + interval_shape + layer_smoothness
    raise ValueError(f"Unknown variant: {variant}")


def summarize(raw_df: pd.DataFrame) -> pd.DataFrame:
    return (
        raw_df.groupby(["n", "gamma", "variant", "family"], sort=True)
        .agg(
            mean_score=("score", "mean"),
            std_score=("score", "std"),
            mean_log_H=("log_H", "mean"),
            mean_penalty=("penalty_effective", "mean"),
            count=("score", "count"),
        )
        .reset_index()
        .sort_values(["variant", "n", "gamma", "mean_score", "family"])
    )


def winner_summary(summary_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (n, gamma, variant), group in summary_df.groupby(["n", "gamma", "variant"], sort=True):
        ordered = group.sort_values(["mean_score", "family"])
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


def pairwise_summary(summary_df: pd.DataFrame, pairwise_targets: list[tuple[str, str]]) -> pd.DataFrame:
    rows = []
    for (n, gamma, variant), group in summary_df.groupby(["n", "gamma", "variant"], sort=True):
        lookup = {str(row["family"]): row for _, row in group.iterrows()}
        for left, right in pairwise_targets:
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
                    "left_wins": bool(float(lookup[left]["mean_score"]) < float(lookup[right]["mean_score"])),
                }
            )
    return pd.DataFrame(rows)


def estimate_gamma_crossings(pairwise_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    grouped = pairwise_df.groupby(["variant", "left_family", "right_family", "n"], sort=True)
    for (variant, left, right, n), group in grouped:
        ordered = group.sort_values("gamma")
        deltas = ordered["delta_score"].tolist()
        gamma_c = None
        status = "no_cross_left_never_wins"

        if all(d <= 0 for d in deltas):
            status = "no_cross_left_always_wins"
        elif any(d < 0 for d in deltas):
            status = "no_cross_left_sometimes_wins"

        for i in range(len(ordered) - 1):
            d1 = float(ordered.iloc[i]["delta_score"])
            d2 = float(ordered.iloc[i + 1]["delta_score"])
            g1 = float(ordered.iloc[i]["gamma"])
            g2 = float(ordered.iloc[i + 1]["gamma"])
            if d1 == 0.0:
                gamma_c = g1
                status = "crossing"
                break
            if d1 * d2 < 0:
                gamma_c = g1 + (0.0 - d1) * (g2 - g1) / (d2 - d1)
                status = "crossing"
                break

        rows.append(
            {
                "variant": str(variant),
                "left_family": str(left),
                "right_family": str(right),
                "n": int(n),
                "gamma_c_est": gamma_c,
                "status": status,
                "min_delta_score": float(ordered["delta_score"].min()),
                "max_delta_score": float(ordered["delta_score"].max()),
                "delta_score_at_min_gamma": float(ordered.iloc[0]["delta_score"]),
                "delta_score_at_max_gamma": float(ordered.iloc[-1]["delta_score"]),
            }
        )
    return pd.DataFrame(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prediction B candidate-family expansion scan.")
    parser.add_argument(
        "--config",
        default="config_prediction_b_candidate_expansion_layered.yaml",
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
    samples_per_family = int(exp_cfg.get("samples_per_family", 4))
    beta = float(exp_cfg.get("beta", 1.0))
    sis_runs = int(exp_cfg.get("sis_runs", 2048))
    default_exact_threshold = int(exp_cfg.get("exact_threshold", 32))
    consistency_weight = float(exp_cfg.get("consistency_weight", 0.682882))
    multi_consistency_weight = float(exp_cfg.get("multi_consistency_weight", consistency_weight))
    family_exact_thresholds = {
        str(k): int(v) for k, v in (exp_cfg.get("family_exact_thresholds", {}) or {}).items()
    }
    pairwise_targets = [
        (str(item["left"]), str(item["right"])) for item in exp_cfg.get("pairwise_targets", [])
    ]

    rows = []
    total_jobs = len(n_values) * len(families) * samples_per_family
    completed_jobs = 0
    for n in n_values:
        print(f"[prediction_b_candidate_expansion] start n={n}", flush=True)
        for family in families:
            generator = FAMILIES[family]
            for sample_id in range(samples_per_family):
                seed = 900000 + 1000 * n + sample_id
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
                        penalty_effective = neutral + geometric_penalty_from_variant(
                            geo,
                            variant=variant,
                            consistency_weight=consistency_weight,
                            multi_consistency_weight=multi_consistency_weight,
                        )
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
                                "penalty_neutral": float(neutral),
                                "penalty_effective": float(penalty_effective),
                                "score": action_value(float(log_h), float(penalty_effective), beta=beta, gamma=gamma),
                                **geo,
                            }
                        )
                completed_jobs += 1
                print(
                    f"[prediction_b_candidate_expansion] done n={n} family={family} sample={sample_id} "
                    f"entropy={entropy_method} progress={completed_jobs}/{total_jobs}",
                    flush=True,
                )

        partial_raw_df = pd.DataFrame(rows)
        partial_summary_df = summarize(partial_raw_df)
        partial_winners_df = winner_summary(partial_summary_df)
        partial_pairwise_df = pairwise_summary(partial_summary_df, pairwise_targets=pairwise_targets)
        partial_gamma_df = estimate_gamma_crossings(partial_pairwise_df) if not partial_pairwise_df.empty else pd.DataFrame()

        partial_raw_df.to_csv(out_dir / "prediction_b_candidate_raw.partial.csv", index=False, encoding="utf-8-sig")
        partial_summary_df.to_csv(out_dir / "prediction_b_candidate_summary.partial.csv", index=False, encoding="utf-8-sig")
        partial_winners_df.to_csv(out_dir / "prediction_b_candidate_winners.partial.csv", index=False, encoding="utf-8-sig")
        partial_pairwise_df.to_csv(out_dir / "prediction_b_candidate_pairwise.partial.csv", index=False, encoding="utf-8-sig")
        partial_gamma_df.to_csv(out_dir / "prediction_b_candidate_gamma_c.partial.csv", index=False, encoding="utf-8-sig")
        print(f"[prediction_b_candidate_expansion] wrote partial outputs for n={n}", flush=True)

    raw_df = pd.DataFrame(rows)
    summary_df = summarize(raw_df)
    winners_df = winner_summary(summary_df)
    pairwise_df = pairwise_summary(summary_df, pairwise_targets=pairwise_targets)
    gamma_df = estimate_gamma_crossings(pairwise_df) if not pairwise_df.empty else pd.DataFrame()

    raw_df.to_csv(out_dir / "prediction_b_candidate_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "prediction_b_candidate_summary.csv", index=False, encoding="utf-8-sig")
    winners_df.to_csv(out_dir / "prediction_b_candidate_winners.csv", index=False, encoding="utf-8-sig")
    pairwise_df.to_csv(out_dir / "prediction_b_candidate_pairwise.csv", index=False, encoding="utf-8-sig")
    gamma_df.to_csv(out_dir / "prediction_b_candidate_gamma_c.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "prediction_b_candidate_winners.csv").as_posix())
    print()
    print(winners_df.head(24).to_string(index=False))
