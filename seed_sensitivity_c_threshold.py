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


def estimate_gamma_crossing(summary_df: pd.DataFrame, left: str, right: str) -> pd.DataFrame:
    rows = []
    for (seed_label, n), group in summary_df.groupby(["seed_label", "n"], sort=True):
        left_df = group[group["family"] == left][["scale_factor", "gamma", "mean_score"]].rename(
            columns={"mean_score": "left_score"}
        )
        right_df = group[group["family"] == right][["scale_factor", "gamma", "mean_score"]].rename(
            columns={"mean_score": "right_score"}
        )
        merged = left_df.merge(right_df, on=["scale_factor", "gamma"], how="inner").sort_values(["scale_factor", "gamma"])
        if merged.empty:
            continue
        rows_for_case = []
        for scale_factor, scale_group in merged.groupby("scale_factor", sort=True):
            delta = scale_group["left_score"].to_numpy() - scale_group["right_score"].to_numpy()
            gammas = scale_group["gamma"].to_numpy(dtype=float)
            gamma_c = None
            status = "no_cross_left_never_wins"
            if all(float(d) <= 0.0 for d in delta):
                status = "no_cross_left_always_wins"
            elif any(float(d) < 0.0 for d in delta):
                status = "no_cross_left_sometimes_wins"
            for i in range(len(delta) - 1):
                d1 = float(delta[i])
                d2 = float(delta[i + 1])
                g1 = float(gammas[i])
                g2 = float(gammas[i + 1])
                if d1 == 0.0:
                    gamma_c = g1
                    status = "crossing"
                    break
                if d1 * d2 < 0:
                    gamma_c = g1 + (0.0 - d1) * (g2 - g1) / (d2 - d1)
                    status = "crossing"
                    break
            rows_for_case.append(
                {
                    "seed_label": str(seed_label),
                    "n": int(n),
                    "scale_factor": float(scale_factor),
                    "gamma_c_est": gamma_c,
                    "status": status,
                    "min_delta_score": float(delta.min()),
                    "max_delta_score": float(delta.max()),
                }
            )
        case_df = pd.DataFrame(rows_for_case).sort_values("scale_factor")
        threshold = None
        threshold_status = "no_cross_anywhere"
        if any(case_df["status"] == "crossing"):
            threshold_status = "crossing_recovered"
            threshold = float(case_df.loc[case_df["status"] == "crossing", "scale_factor"].min())
        rows.append(
            {
                "seed_label": str(seed_label),
                "n": int(n),
                "c_threshold_est": threshold,
                "threshold_status": threshold_status,
                "scale_factor_min": float(case_df["scale_factor"].min()),
                "scale_factor_max": float(case_df["scale_factor"].max()),
            }
        )
    return pd.DataFrame(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Seed sensitivity scan for c_N thresholds.")
    parser.add_argument("--config", default="config_seed_sensitivity_c_threshold.yaml", help="Path to YAML config.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    exp_cfg = config["experiment"]
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    n_values = [int(x) for x in exp_cfg["n_values"]]
    families = tuple(exp_cfg["families"])
    gammas = tuple(float(g) for g in exp_cfg["gammas"])
    samples_per_family = int(exp_cfg.get("samples_per_family", 2))
    beta = float(exp_cfg.get("beta", 1.0))
    fixed_matched_weight = float(exp_cfg["fixed_matched_weight"])
    sis_runs = int(exp_cfg.get("sis_runs", 1024))
    default_exact_threshold = int(exp_cfg.get("exact_threshold", 44))
    family_exact_thresholds = {
        str(k): int(v) for k, v in (exp_cfg.get("family_exact_thresholds", {}) or {}).items()
    }
    seed_offsets = [int(x) for x in exp_cfg["seed_offsets"]]
    per_n_scale_factors = {
        int(k): [float(x) for x in v] for k, v in exp_cfg["per_n_scale_factors"].items()
    }

    base_rows = []
    for seed_offset in seed_offsets:
        seed_label = f"offset_{seed_offset}"
        for n in n_values:
            for family in families:
                generator = FAMILIES[family]
                for sample_id in range(samples_per_family):
                    seed = seed_offset + 1000 * n + sample_id
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
                    base_rows.append(
                        {
                            "seed_label": seed_label,
                            "seed_offset": seed_offset,
                            "n": n,
                            "family": family,
                            "sample_id": sample_id,
                            "seed": seed,
                            "beta": beta,
                            "log_H": log_h,
                            "entropy_method": entropy_method,
                            "penalty_neutral": neutral,
                            **geo,
                        }
                    )

    base_df = pd.DataFrame(base_rows)

    rows = []
    for row in base_df.to_dict(orient="records"):
        for scale_factor in per_n_scale_factors[int(row["n"])]:
            consistency_weight = fixed_matched_weight * scale_factor
            penalty = (
                row["penalty_neutral"]
                + 2.0 * row["geo_width_height"]
                + consistency_weight * row["geo_dim_consistency"]
            )
            for gamma in gammas:
                rows.append(
                    {
                        **row,
                        "scale_factor": scale_factor,
                        "matched_weight": fixed_matched_weight,
                        "consistency_weight": consistency_weight,
                        "gamma": gamma,
                        "penalty_effective": penalty,
                        "score": action_value(row["log_H"], penalty, beta=beta, gamma=gamma),
                    }
                )

    raw_df = pd.DataFrame(rows)
    summary_df = (
        raw_df.groupby(["seed_label", "n", "scale_factor", "gamma", "family"])
        .agg(
            mean_score=("score", "mean"),
            std_score=("score", "std"),
            mean_penalty=("penalty_effective", "mean"),
            mean_log_H=("log_H", "mean"),
            entropy_method=("entropy_method", lambda s: ",".join(sorted(set(str(x) for x in s)))),
            count=("score", "count"),
        )
        .reset_index()
        .sort_values(["seed_label", "n", "scale_factor", "gamma", "family"])
    )

    threshold_df = estimate_gamma_crossing(summary_df, "lorentzian_like_2d", "KR_like")
    sensitivity_df = (
        threshold_df.groupby("n")
        .agg(
            n_runs=("c_threshold_est", "count"),
            mean_c_threshold=("c_threshold_est", "mean"),
            min_c_threshold=("c_threshold_est", "min"),
            max_c_threshold=("c_threshold_est", "max"),
            std_c_threshold=("c_threshold_est", "std"),
        )
        .reset_index()
    )

    raw_df.to_csv(out_dir / "seed_sensitivity_c_threshold_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "seed_sensitivity_c_threshold_summary.csv", index=False, encoding="utf-8-sig")
    threshold_df.to_csv(out_dir / "seed_sensitivity_c_thresholds.csv", index=False, encoding="utf-8-sig")
    sensitivity_df.to_csv(out_dir / "seed_sensitivity_c_threshold_sensitivity.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "seed_sensitivity_c_threshold_sensitivity.csv").as_posix())
    print()
    print(sensitivity_df.to_string(index=False))
