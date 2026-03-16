"""Fault-tolerant version of benchmark_exact_timing.py.

Saves results incrementally after each case, so a crash on one case
does not lose previously computed results.
"""
from __future__ import annotations

import argparse
import csv
import gc
import time
from pathlib import Path

import yaml

from entropy_exact import count_linear_extensions_exact
from experiment import FAMILIES
from observables import antichain_width, antichain_width_ratio, comparable_fraction
from observables_geo import height_ratio, width_ratio


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


FIELDNAMES = [
    "family", "n", "seed", "elapsed_seconds", "log10_count",
    "comparable_fraction", "antichain_width", "antichain_width_ratio",
    "height_ratio", "width_ratio", "status",
]


def run_single(family: str, n: int, seed: int) -> dict:
    generator = FAMILIES[family]
    poset = generator(n=n, seed=seed)

    start = time.perf_counter()
    count = count_linear_extensions_exact(poset)
    elapsed = time.perf_counter() - start

    return {
        "family": family,
        "n": n,
        "seed": seed,
        "elapsed_seconds": round(elapsed, 3),
        "log10_count": len(str(count)) - 1 if count > 0 else 0,
        "comparable_fraction": round(comparable_fraction(poset), 4),
        "antichain_width": antichain_width(poset),
        "antichain_width_ratio": round(antichain_width_ratio(poset), 4),
        "height_ratio": round(height_ratio(poset), 4),
        "width_ratio": round(width_ratio(poset), 4),
        "status": "ok",
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Benchmark exact timing (fault-tolerant).")
    parser.add_argument("--config", default="config_exact_timing.yaml")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "exact_timing_benchmark.csv"

    # Write header
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()

    results = []
    for case in config["cases"]:
        family = str(case["family"])
        n = int(case["n"])
        seed = int(case.get("seed", 1000 * n))

        try:
            row = run_single(family, n, seed)
            print(
                f"{family:22s} n={n:<3d} "
                f"time={row['elapsed_seconds']:8.3f}s "
                f"comp={row['comparable_fraction']:.3f} "
                f"aw={row['antichain_width']:<3d} "
                f"h={row['height_ratio']:.3f} "
                f"w={row['width_ratio']:.3f}",
                flush=True,
            )
        except Exception as e:
            row = {
                "family": family, "n": n, "seed": seed,
                "elapsed_seconds": -1, "log10_count": -1,
                "comparable_fraction": -1, "antichain_width": -1,
                "antichain_width_ratio": -1, "height_ratio": -1,
                "width_ratio": -1, "status": f"FAIL:{type(e).__name__}",
            }
            print(f"{family:22s} n={n:<3d}  FAILED: {type(e).__name__}: {e}", flush=True)

        results.append(row)

        # Append to CSV immediately
        with open(csv_path, "a", newline="", encoding="utf-8-sig") as f:
            writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
            writer.writerow(row)

        # Force garbage collection to free lru_cache memory from previous cases
        gc.collect()

    # Write summary tables
    import pandas as pd
    df = pd.DataFrame(results)

    ok = df[df["status"] == "ok"]
    if not ok.empty:
        summary = ok.pivot(index="n", columns="family", values="elapsed_seconds").sort_index().reset_index()
        summary.to_csv(out_dir / "exact_timing_summary.csv", index=False, encoding="utf-8-sig")

        width_summary = ok.pivot(index="n", columns="family", values="antichain_width").sort_index().reset_index()
        width_summary.to_csv(out_dir / "exact_width_summary.csv", index=False, encoding="utf-8-sig")

        print("\nTiming summary:")
        print(summary.to_string(index=False))
        print("\nWidth summary:")
        print(width_summary.to_string(index=False))

    failed = df[df["status"] != "ok"]
    if not failed.empty:
        print(f"\n{len(failed)} cases failed:")
        for _, r in failed.iterrows():
            print(f"  {r['family']} n={r['n']}: {r['status']}")
