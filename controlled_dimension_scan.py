from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from action import action_value
from experiment import FAMILIES
from observables import comparable_fraction, neutral_penalty
from observables_geo import height_ratio, width_ratio
from runtime_utils import estimate_entropy


WINDOW_METRICS = ["comp", "height_ratio", "width_ratio"]
DIMENSION_FAMILIES = ["lorentzian_like_2d", "lorentzian_like_3d", "lorentzian_like_4d"]


def load_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def sample_metrics(poset) -> dict[str, float]:
    return {
        "comp": comparable_fraction(poset),
        "height_ratio": height_ratio(poset),
        "width_ratio": width_ratio(poset),
        "neutral_penalty": neutral_penalty(poset),
    }


def build_reference_windows(
    *,
    n: int,
    reference_family: str,
    pool_size: int,
    quantile_low: float,
    quantile_high: float,
) -> dict[str, tuple[float, float]]:
    generator = FAMILIES[reference_family]
    rows = []
    for sample_id in range(pool_size):
        poset = generator(n=n, seed=10_000 * n + sample_id)
        rows.append(sample_metrics(poset))
    frame = pd.DataFrame(rows)
    windows = {}
    for metric in WINDOW_METRICS:
        low = float(frame[metric].quantile(quantile_low))
        high = float(frame[metric].quantile(quantile_high))
        windows[metric] = (low, high)
    return windows


def in_windows(
    metrics: dict[str, float],
    windows: dict[str, tuple[float, float]],
    active_metrics: list[str],
) -> bool:
    for metric in active_metrics:
        low, high = windows[metric]
        value = float(metrics[metric])
        if value < low or value > high:
            return False
    return True


def collect_controlled_samples(
    *,
    n: int,
    family: str,
    windows: dict[str, tuple[float, float]],
    active_metrics: list[str],
    accepted_target: int,
    max_attempts: int,
    sis_runs: int,
    exact_threshold: int,
    beta: float,
    gamma: float,
) -> tuple[list[dict], dict[str, float]]:
    generator = FAMILIES[family]
    accepted = []
    attempts = 0
    while len(accepted) < accepted_target and attempts < max_attempts:
        seed = 50_000 * n + attempts
        poset = generator(n=n, seed=seed)
        metrics = sample_metrics(poset)
        attempts += 1
        if not in_windows(metrics, windows, active_metrics):
            continue
        log_h, entropy_method = estimate_entropy(
            poset,
            sis_runs=sis_runs,
            seed=seed,
            exact_threshold=exact_threshold,
        )
        score_local = action_value(
            log_extensions=log_h,
            penalty=metrics["neutral_penalty"],
            beta=beta,
            gamma=gamma,
        )
        accepted.append(
            {
                "n": n,
                "family": family,
                "seed": seed,
                "entropy_method": entropy_method,
                "log_H": log_h,
                "score_local": score_local,
                **metrics,
            }
        )

    stats = {
        "n": n,
        "family": family,
        "attempts": attempts,
        "accepted": len(accepted),
        "acceptance_rate": float(len(accepted) / attempts) if attempts else 0.0,
    }
    return accepted, stats


def summarize_results(raw_df: pd.DataFrame, attempts_df: pd.DataFrame) -> pd.DataFrame:
    if raw_df.empty:
        return attempts_df.copy()
    summary = (
        raw_df.groupby(["regime", "n", "family"])
        .agg(
            mean_score_local=("score_local", "mean"),
            mean_log_H=("log_H", "mean"),
            mean_neutral_penalty=("neutral_penalty", "mean"),
            mean_comp=("comp", "mean"),
            mean_height_ratio=("height_ratio", "mean"),
            mean_width_ratio=("width_ratio", "mean"),
            count=("score_local", "count"),
        )
        .reset_index()
    )
    summary = summary.merge(attempts_df, on=["regime", "n", "family"], how="right")
    summary["rank_local"] = (
        summary.groupby("regime")["mean_score_local"].rank(method="dense", ascending=True).astype("Int64")
    )
    return summary.sort_values(["regime", "n", "rank_local", "family"], na_position="last")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Controlled dimension scan under matched structural windows.")
    parser.add_argument("--config", default="config_controlled_dimension.yaml", help="Path to YAML config file.")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    config = load_config(args.config)

    exp_cfg = config["experiment"]
    ctl_cfg = config["controlled_scan"]
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    stats_rows = []
    window_rows = []
    regimes = config["controlled_scan"]["regimes"]
    for n in exp_cfg["n_values"]:
        windows = build_reference_windows(
            n=int(n),
            reference_family=str(ctl_cfg["reference_family"]),
            pool_size=int(ctl_cfg["reference_pool_size"]),
            quantile_low=float(ctl_cfg["quantile_low"]),
            quantile_high=float(ctl_cfg["quantile_high"]),
        )
        for metric, bounds in windows.items():
            window_rows.append({"n": int(n), "metric": metric, "lower": bounds[0], "upper": bounds[1]})
        for regime_name, regime_metrics in regimes.items():
            for family in ctl_cfg["families"]:
                accepted, stats = collect_controlled_samples(
                    n=int(n),
                    family=str(family),
                    windows=windows,
                    active_metrics=list(regime_metrics),
                    accepted_target=int(ctl_cfg["accepted_target"]),
                    max_attempts=int(ctl_cfg["max_attempts"]),
                    sis_runs=int(exp_cfg["sis_runs"]),
                    exact_threshold=int(exp_cfg["exact_threshold"]),
                    beta=float(exp_cfg["beta"]),
                    gamma=float(exp_cfg["gamma"]),
                )
                for row in accepted:
                    row["regime"] = str(regime_name)
                stats["regime"] = str(regime_name)
                rows.extend(accepted)
                stats_rows.append(stats)

    raw_df = pd.DataFrame(rows)
    attempts_df = pd.DataFrame(stats_rows)
    windows_df = pd.DataFrame(window_rows)
    summary_df = summarize_results(raw_df, attempts_df)

    raw_df.to_csv(out_dir / "controlled_dimension_raw.csv", index=False, encoding="utf-8-sig")
    attempts_df.to_csv(out_dir / "controlled_dimension_attempts.csv", index=False, encoding="utf-8-sig")
    windows_df.to_csv(out_dir / "controlled_dimension_windows.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "controlled_dimension_summary.csv", index=False, encoding="utf-8-sig")

    print("Reference windows:")
    print(windows_df.to_string(index=False))
    print()
    print("Controlled dimension summary:")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
