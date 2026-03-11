from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from entropy_exact import log_linear_extensions_exact
from entropy_sis import estimate_log_linear_extensions_sis
from experiment import FAMILIES, load_config


def run_calibration(
    n_values: tuple[int, ...],
    families: tuple[str, ...],
    samples_per_family: int,
    sis_runs: int,
) -> pd.DataFrame:
    rows = []
    for n in n_values:
        for family in families:
            generator = FAMILIES[family]
            for sample_id in range(samples_per_family):
                seed = 5000 * n + sample_id
                poset = generator(n=n, seed=seed)
                log_exact = log_linear_extensions_exact(poset)
                log_sis_mean, log_sis_std = estimate_log_linear_extensions_sis(
                    poset,
                    n_runs=sis_runs,
                    seed=seed,
                )
                rows.append(
                    {
                        "n": n,
                        "family": family,
                        "sample_id": sample_id,
                        "log_H_exact": log_exact,
                        "log_H_sis_mean": log_sis_mean,
                        "log_H_sis_std": log_sis_std,
                        "abs_error": abs(log_sis_mean - log_exact),
                        "signed_error": log_sis_mean - log_exact,
                    }
                )
    return pd.DataFrame(rows)


def summarize_calibration(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["n", "family"])
        .agg(
            mean_abs_error=("abs_error", "mean"),
            max_abs_error=("abs_error", "max"),
            mean_signed_error=("signed_error", "mean"),
            mean_sis_std=("log_H_sis_std", "mean"),
            count=("abs_error", "count"),
        )
        .reset_index()
        .sort_values(["n", "mean_abs_error", "family"])
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run SIS vs exact entropy calibration.")
    parser.add_argument("--config", default="config.yaml", help="Path to YAML config file.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    calib_cfg = config["calibration"]
    out_cfg = config["output"]

    families = tuple(calib_cfg.get("families", list(FAMILIES.keys())))
    df = run_calibration(
        n_values=tuple(calib_cfg["n_values"]),
        families=families,
        samples_per_family=int(calib_cfg["samples_per_family"]),
        sis_runs=int(calib_cfg["sis_runs"]),
    )
    summary = summarize_calibration(df)

    out_dir = Path(out_cfg["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / out_cfg["calibration_raw_csv"], index=False, encoding="utf-8-sig")
    summary.to_csv(out_dir / out_cfg["calibration_summary_csv"], index=False, encoding="utf-8-sig")

    print(summary.to_string(index=False))
