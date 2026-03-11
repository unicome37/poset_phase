from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd
import yaml

from action import action_value
from entropy_exact import count_linear_extensions_exact
from experiment import FAMILIES
from observables import antichain_width, antichain_width_ratio, comparable_fraction, neutral_penalty
from observables_geo import geometric_penalty


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def quantile_window(values: pd.Series, lower_q: float, upper_q: float) -> tuple[float, float]:
    return float(values.quantile(lower_q)), float(values.quantile(upper_q))


def build_reference_samples(
    family: str,
    n: int,
    n_samples: int,
    seed_base: int,
) -> pd.DataFrame:
    generator = FAMILIES[family]
    rows: list[dict] = []
    for sample_id in range(n_samples):
        seed = seed_base + sample_id
        poset = generator(n=n, seed=seed)
        rows.append(
            {
                "family": family,
                "n": n,
                "seed": seed,
                "antichain_width": antichain_width(poset),
                "antichain_width_ratio": antichain_width_ratio(poset),
                "comparable_fraction": comparable_fraction(poset),
            }
        )
    return pd.DataFrame(rows)


def within_window(row: dict, criterion: str, width_window: tuple[float, float], comp_window: tuple[float, float]) -> bool:
    width_ok = width_window[0] <= row["antichain_width"] <= width_window[1]
    comp_ok = comp_window[0] <= row["comparable_fraction"] <= comp_window[1]
    if criterion == "width_only":
        return width_ok
    if criterion == "comp_only":
        return comp_ok
    if criterion == "width_and_comp":
        return width_ok and comp_ok
    raise ValueError(f"Unsupported criterion: {criterion}")


def scan_family(
    family: str,
    n: int,
    accepted_target: int,
    max_attempts: int,
    criterion: str,
    width_window: tuple[float, float],
    comp_window: tuple[float, float],
    beta: float,
    gamma: float,
    seed_base: int,
) -> tuple[pd.DataFrame, dict]:
    generator = FAMILIES[family]
    rows: list[dict] = []
    accepted = 0
    attempts = 0

    while attempts < max_attempts and accepted < accepted_target:
        seed = seed_base + attempts
        poset = generator(n=n, seed=seed)
        row = {
            "family": family,
            "n": n,
            "seed": seed,
            "criterion": criterion,
            "antichain_width": antichain_width(poset),
            "antichain_width_ratio": antichain_width_ratio(poset),
            "comparable_fraction": comparable_fraction(poset),
        }
        attempts += 1
        if not within_window(row, criterion, width_window, comp_window):
            continue

        count = count_linear_extensions_exact(poset)
        log_h = math.log(count)
        pen_neutral = neutral_penalty(poset)
        pen_geo = geometric_penalty(poset)

        row.update(
            {
                "accepted": 1,
                "log_H": log_h,
                "penalty_neutral": pen_neutral,
                "penalty_geometric": pen_geo,
                "score_A1_gamma": action_value(log_h, pen_neutral, beta=beta, gamma=gamma),
                "score_A2_gamma": action_value(log_h, pen_neutral + pen_geo, beta=beta, gamma=gamma),
            }
        )
        rows.append(row)
        accepted += 1

    summary = {
        "family": family,
        "n": n,
        "criterion": criterion,
        "attempts": attempts,
        "accepted": accepted,
        "acceptance_rate": float(accepted / attempts) if attempts else 0.0,
    }
    return pd.DataFrame(rows), summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Controlled compressibility scan with exact A2 comparisons.")
    parser.add_argument("--config", default="config_controlled_compressibility.yaml", help="Path to YAML config file.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)

    exp_cfg = config["experiment"]
    ref_cfg = config["reference"]
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    beta = float(exp_cfg["beta"])
    gamma = float(exp_cfg["gamma"])
    n_values = tuple(int(v) for v in exp_cfg["n_values"])
    families = tuple(str(v) for v in exp_cfg["families"])
    criteria = tuple(str(v) for v in exp_cfg["criteria"])
    accepted_target = int(exp_cfg["accepted_target"])
    max_attempts = int(exp_cfg["max_attempts"])

    ref_family = str(ref_cfg["family"])
    ref_samples = int(ref_cfg["samples"])
    lower_q = float(ref_cfg["lower_quantile"])
    upper_q = float(ref_cfg["upper_quantile"])
    ref_seed_base = int(ref_cfg.get("seed_base", 100000))
    candidate_seed_base = int(exp_cfg.get("seed_base", 200000))

    reference_rows: list[pd.DataFrame] = []
    accepted_rows: list[pd.DataFrame] = []
    summary_rows: list[dict] = []
    window_rows: list[dict] = []

    for n in n_values:
        ref_df = build_reference_samples(
            family=ref_family,
            n=n,
            n_samples=ref_samples,
            seed_base=ref_seed_base + 1000 * n,
        )
        reference_rows.append(ref_df)
        width_window = quantile_window(ref_df["antichain_width"], lower_q, upper_q)
        comp_window = quantile_window(ref_df["comparable_fraction"], lower_q, upper_q)
        window_rows.append(
            {
                "n": n,
                "reference_family": ref_family,
                "width_lower": width_window[0],
                "width_upper": width_window[1],
                "comp_lower": comp_window[0],
                "comp_upper": comp_window[1],
            }
        )

        for criterion in criteria:
            for family in families:
                df_acc, summary = scan_family(
                    family=family,
                    n=n,
                    accepted_target=accepted_target,
                    max_attempts=max_attempts,
                    criterion=criterion,
                    width_window=width_window,
                    comp_window=comp_window,
                    beta=beta,
                    gamma=gamma,
                    seed_base=candidate_seed_base + 100000 * n + 1000 * criteria.index(criterion) + 100 * families.index(family),
                )
                if not df_acc.empty:
                    accepted_rows.append(df_acc)
                summary_rows.append(summary)
                print(
                    f"{criterion:15s} n={n:<3d} {family:22s} "
                    f"accepted={summary['accepted']:<2d}/{summary['attempts']:<3d} "
                    f"rate={summary['acceptance_rate']:.3f}"
                )

    reference_all = pd.concat(reference_rows, ignore_index=True)
    accepted_all = pd.concat(accepted_rows, ignore_index=True) if accepted_rows else pd.DataFrame()
    summary_df = pd.DataFrame(summary_rows)
    windows_df = pd.DataFrame(window_rows)

    ranking_df = pd.DataFrame()
    if not accepted_all.empty:
        ranking_df = (
            accepted_all.groupby(["criterion", "n", "family"])
            .agg(
                mean_score_A2=("score_A2_gamma", "mean"),
                std_score_A2=("score_A2_gamma", "std"),
                mean_log_H=("log_H", "mean"),
                mean_antichain_width=("antichain_width", "mean"),
                mean_comp=("comparable_fraction", "mean"),
                count=("score_A2_gamma", "count"),
            )
            .reset_index()
            .sort_values(["criterion", "n", "mean_score_A2"])
        )

    reference_all.to_csv(out_dir / "controlled_compressibility_reference.csv", index=False, encoding="utf-8-sig")
    accepted_all.to_csv(out_dir / "controlled_compressibility_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "controlled_compressibility_summary.csv", index=False, encoding="utf-8-sig")
    windows_df.to_csv(out_dir / "controlled_compressibility_windows.csv", index=False, encoding="utf-8-sig")
    ranking_df.to_csv(out_dir / "controlled_compressibility_rankings.csv", index=False, encoding="utf-8-sig")

    if not ranking_df.empty:
        print()
        print(ranking_df.to_string(index=False))
