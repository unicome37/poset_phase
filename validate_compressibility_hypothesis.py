from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

import pandas as pd
import yaml

from action import action_value
from entropy_exact import count_linear_extensions_exact
from experiment import FAMILIES
from observables import antichain_width, antichain_width_ratio, comparable_fraction, neutral_penalty
from observables_geo import geometric_penalty, height_ratio, width_ratio


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def run_cases(cases: list[dict], beta: float, gamma: float) -> pd.DataFrame:
    rows: list[dict] = []
    for case in cases:
        family = str(case["family"])
        n = int(case["n"])
        seed = int(case.get("seed", 1000 * n))
        generator = FAMILIES[family]
        poset = generator(n=n, seed=seed)

        start = time.perf_counter()
        count = count_linear_extensions_exact(poset)
        elapsed = time.perf_counter() - start
        log_h = math.log(count)
        pen_neutral = neutral_penalty(poset)
        pen_geo = geometric_penalty(poset)

        rows.append(
            {
                "family": family,
                "n": n,
                "seed": seed,
                "elapsed_seconds": elapsed,
                "log_elapsed": math.log(max(elapsed, 1e-9)),
                "log_H": log_h,
                "antichain_width": antichain_width(poset),
                "antichain_width_ratio": antichain_width_ratio(poset),
                "comparable_fraction": comparable_fraction(poset),
                "height_ratio": height_ratio(poset),
                "layer_width_ratio": width_ratio(poset),
                "penalty_neutral": pen_neutral,
                "penalty_geometric": pen_geo,
                "score_A1_gamma": action_value(log_h, pen_neutral, beta=beta, gamma=gamma),
                "score_A2_gamma": action_value(log_h, pen_neutral + pen_geo, beta=beta, gamma=gamma),
            }
        )
        print(
            f"{family:22s} n={n:<3d} time={elapsed:8.3f}s "
            f"aw={rows[-1]['antichain_width']:<3d} "
            f"logH={log_h:8.3f} "
            f"A1={rows[-1]['score_A1_gamma']:8.3f} "
            f"A2={rows[-1]['score_A2_gamma']:8.3f}"
        )
    return pd.DataFrame(rows)


def correlation_table(df: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "elapsed_seconds",
        "log_elapsed",
        "log_H",
        "antichain_width",
        "antichain_width_ratio",
        "comparable_fraction",
        "penalty_neutral",
        "penalty_geometric",
        "score_A1_gamma",
        "score_A2_gamma",
    ]
    corr = df[cols].corr(method="spearman")
    return corr.reset_index().rename(columns={"index": "metric"})


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate compressibility hypothesis with exact metrics.")
    parser.add_argument("--config", default="config_compressibility_validation.yaml", help="Path to YAML config file.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    beta = float(config["experiment"]["beta"])
    gamma = float(config["experiment"]["gamma"])
    df = run_cases(config["cases"], beta=beta, gamma=gamma)
    corr = correlation_table(df)

    df.to_csv(out_dir / "compressibility_validation_raw.csv", index=False, encoding="utf-8-sig")
    corr.to_csv(out_dir / "compressibility_validation_corr.csv", index=False, encoding="utf-8-sig")

    print()
    print("Correlation table:")
    print(corr.to_string(index=False))
