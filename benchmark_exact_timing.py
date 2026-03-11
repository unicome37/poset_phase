from __future__ import annotations

import argparse
import time
from pathlib import Path

import pandas as pd
import yaml

from entropy_exact import count_linear_extensions_exact
from experiment import FAMILIES
from observables import antichain_width, antichain_width_ratio, comparable_fraction
from observables_geo import height_ratio, width_ratio


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def run_benchmark(cases: list[dict]) -> pd.DataFrame:
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

        rows.append(
            {
                "family": family,
                "n": n,
                "seed": seed,
                "elapsed_seconds": elapsed,
                "log10_count": len(str(count)) - 1 if count > 0 else 0,
                "comparable_fraction": comparable_fraction(poset),
                "antichain_width": antichain_width(poset),
                "antichain_width_ratio": antichain_width_ratio(poset),
                "height_ratio": height_ratio(poset),
                "width_ratio": width_ratio(poset),
            }
        )
        print(
            f"{family:22s} n={n:<3d} "
            f"time={elapsed:8.3f}s "
            f"comp={rows[-1]['comparable_fraction']:.3f} "
            f"aw={rows[-1]['antichain_width']:<3d} "
            f"h={rows[-1]['height_ratio']:.3f} "
            f"w={rows[-1]['width_ratio']:.3f}"
        )
    return pd.DataFrame(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Benchmark exact linear-extension timing.")
    parser.add_argument("--config", default="config_exact_timing.yaml", help="Path to YAML config file.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    df = run_benchmark(config["cases"])
    out_path = out_dir / "exact_timing_benchmark.csv"
    df.to_csv(out_path, index=False, encoding="utf-8-sig")

    summary = (
        df.pivot(index="n", columns="family", values="elapsed_seconds")
        .sort_index()
        .reset_index()
    )
    summary_path = out_dir / "exact_timing_summary.csv"
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")

    width_summary = (
        df.pivot(index="n", columns="family", values="antichain_width")
        .sort_index()
        .reset_index()
    )
    width_summary_path = out_dir / "exact_width_summary.csv"
    width_summary.to_csv(width_summary_path, index=False, encoding="utf-8-sig")

    print()
    print("Timing summary:")
    print(summary.to_string(index=False))
    print()
    print("Width summary:")
    print(width_summary.to_string(index=False))
