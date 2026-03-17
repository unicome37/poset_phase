from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from action import action_value, get_action_penalty
from bootstrap import bootstrap_group_summary
from entropy_exact import log_linear_extensions_exact
from entropy_sis import estimate_log_linear_extensions_sis
from generators import (
    generate_absolute_layered,
    generate_interval_order,
    generate_kr_like,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
    generate_multi_layer_random,
    generate_random_layered_k4_uniform,
    generate_random_layered_k6_middle_heavy,
    generate_random_layered_k6_tapered,
    generate_random_layered_k6_longjump,
    generate_random_layered_k6_uniform,
    generate_random_layered_k8_uniform,
    generate_transitive_percolation,
)
from normalization import add_normalized_columns, add_size_scaled_columns
from observables import neutral_penalty
from observables_geo import geometric_penalty


FAMILIES = {
    "absolute_layered": generate_absolute_layered,
    "KR_like": generate_kr_like,
    "lorentzian_like_2d": generate_lorentzian_like_2d,
    "lorentzian_like_3d": generate_lorentzian_like_3d,
    "lorentzian_like_4d": generate_lorentzian_like_4d,
    "lorentzian_like_5d": generate_lorentzian_like_5d,
    "transitive_percolation": generate_transitive_percolation,
    "interval_order": generate_interval_order,
    "multi_layer_random": generate_multi_layer_random,
    "random_layered_k4_uniform": generate_random_layered_k4_uniform,
    "random_layered_k6_uniform": generate_random_layered_k6_uniform,
    "random_layered_k8_uniform": generate_random_layered_k8_uniform,
    "random_layered_k6_tapered": generate_random_layered_k6_tapered,
    "random_layered_k6_middle_heavy": generate_random_layered_k6_middle_heavy,
    "random_layered_k6_longjump": generate_random_layered_k6_longjump,
}


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def run_experiment(
    n_values: tuple[int, ...] = (20, 40, 60, 80),
    families: tuple[str, ...] | None = None,
    gammas: tuple[float, ...] = (0.0, 0.1, 0.2, 0.4, 0.8),
    beta: float = 1.0,
    samples_per_family: int = 16,
    sis_runs: int = 64,
    action_modes: tuple[str, ...] = ("A1", "A2", "A3"),
    exact_threshold: int = 0,
) -> pd.DataFrame:
    rows = []
    selected_families = families or tuple(FAMILIES.keys())

    for n in n_values:
        for family in selected_families:
            generator = FAMILIES[family]
            for sample_id in range(samples_per_family):
                seed = 1000 * n + sample_id
                poset = generator(n=n, seed=seed)
                if n <= exact_threshold:
                    log_h_mean = log_linear_extensions_exact(poset)
                    log_h_std = 0.0
                    entropy_method = "exact"
                else:
                    log_h_mean, log_h_std = estimate_log_linear_extensions_sis(
                        poset,
                        n_runs=sis_runs,
                        seed=seed,
                    )
                    entropy_method = "sis"
                penalty_neutral = neutral_penalty(poset)
                penalty_geometric = geometric_penalty(poset)

                for gamma in gammas:
                    for action_mode in action_modes:
                        penalty_effective = get_action_penalty(poset, action_mode)
                        rows.append(
                            {
                                "n": n,
                                "family": family,
                                "sample_id": sample_id,
                                "gamma": gamma,
                                "beta": beta,
                                "action_mode": action_mode,
                                "log_H_mean": log_h_mean,
                                "log_H_std": log_h_std,
                                "entropy_method": entropy_method,
                                "penalty_neutral": penalty_neutral,
                                "penalty_geometric": penalty_geometric,
                                "penalty_effective": penalty_effective,
                                "score": action_value(
                                    log_extensions=log_h_mean,
                                    penalty=penalty_effective,
                                    beta=beta,
                                    gamma=gamma,
                                ),
                            }
                        )

    return pd.DataFrame(rows)


def summarize_scores(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["n", "gamma", "family", "action_mode"])
        .agg(
            mean_score=("score", "mean"),
            std_score=("score", "std"),
            mean_log_H=("log_H_mean", "mean"),
            mean_penalty=("penalty_effective", "mean"),
            count=("score", "count"),
        )
        .reset_index()
        .sort_values(["action_mode", "n", "gamma", "mean_score"])
    )


def summarize_normalized_scores(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["n", "gamma", "family", "action_mode"])
        .agg(
            mean_score_norm=("score_norm", "mean"),
            std_score_norm=("score_norm", "std"),
            mean_score_norm_std_est=("score_norm_std_est", "mean"),
            mean_log_H_norm=("log_H_norm", "mean"),
            mean_penalty_norm=("penalty_norm", "mean"),
            count=("score_norm", "count"),
        )
        .reset_index()
        .sort_values(["action_mode", "n", "gamma", "mean_score_norm"])
    )


def write_outputs(
    raw_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    bootstrap_df: pd.DataFrame | None,
    config: dict,
) -> None:
    output_cfg = config["output"]
    out_dir = Path(output_cfg["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    raw_df.to_csv(out_dir / output_cfg["raw_samples_csv"], index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / output_cfg["summary_csv"], index=False, encoding="utf-8-sig")

    if bootstrap_df is not None:
        bootstrap_df.to_csv(out_dir / output_cfg["bootstrap_csv"], index=False, encoding="utf-8-sig")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run minimal poset phase experiment.")
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config file.",
    )
    return parser


if __name__ == "__main__":
    parser = build_arg_parser()
    args = parser.parse_args()

    config = load_config(args.config)
    exp_cfg = config["experiment"]
    norm_cfg = config["normalization"]
    boot_cfg = config["bootstrap"]

    df = run_experiment(
        n_values=tuple(exp_cfg["n_values"]),
        families=tuple(exp_cfg.get("families", list(FAMILIES.keys()))),
        gammas=tuple(exp_cfg["gammas"]),
        beta=float(exp_cfg["beta"]),
        samples_per_family=int(exp_cfg["samples_per_family"]),
        sis_runs=int(exp_cfg["sis_runs"]),
        action_modes=tuple(exp_cfg.get("action_modes", ["A1", "A2", "A3"])),
        exact_threshold=int(exp_cfg.get("exact_threshold", 0)),
    )
    df = add_size_scaled_columns(df)
    df = add_normalized_columns(
        df,
        method=str(norm_cfg["method"]),
        group_cols=tuple(norm_cfg["group_by"]),
    )

    summary = summarize_normalized_scores(df)

    bootstrap_df = None
    if boot_cfg.get("enabled", False):
        bootstrap_df = bootstrap_group_summary(
            df,
            value_col=str(boot_cfg["value_col"]),
            group_cols=("n", "gamma", "family", "action_mode"),
            n_bootstrap=int(boot_cfg["n_bootstrap"]),
            seed=int(boot_cfg["seed"]),
            value_std_col=boot_cfg.get("value_std_col"),
        )

    write_outputs(df, summary, bootstrap_df, config)

    print(summary.head(40).to_string(index=False))
    if bootstrap_df is not None:
        print("\nBootstrap summary:")
        print(bootstrap_df.head(20).to_string(index=False))
