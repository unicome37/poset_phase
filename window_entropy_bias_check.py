from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd
import yaml

from entropy_exact import count_linear_extensions_exact
from experiment import FAMILIES
from observables import antichain_width, comparable_fraction
from pairwise_compressibility_duel import reference_window


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def in_window(width: int, comp: float, width_window: tuple[float, float], comp_window: tuple[float, float]) -> bool:
    return width_window[0] <= width <= width_window[1] and comp_window[0] <= comp <= comp_window[1]


def sample_family(
    family: str,
    n: int,
    n_samples: int,
    width_window: tuple[float, float],
    comp_window: tuple[float, float],
    seed_base: int,
) -> pd.DataFrame:
    generator = FAMILIES[family]
    rows: list[dict] = []
    for sample_id in range(n_samples):
        seed = seed_base + 1000 * n + sample_id
        poset = generator(n=n, seed=seed)
        width = antichain_width(poset)
        comp = comparable_fraction(poset)
        count = count_linear_extensions_exact(poset)
        rows.append(
            {
                "family": family,
                "n": n,
                "seed": seed,
                "sample_id": sample_id,
                "antichain_width": width,
                "comparable_fraction": comp,
                "log_H": math.log(count),
                "accepted_by_window": int(in_window(width, comp, width_window, comp_window)),
            }
        )
    return pd.DataFrame(rows)


def summarize_bias(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for (n, family), sub in df.groupby(["n", "family"]):
        accepted = sub[sub["accepted_by_window"] == 1]
        rejected = sub[sub["accepted_by_window"] == 0]
        rows.append(
            {
                "n": int(n),
                "family": str(family),
                "count_total": int(len(sub)),
                "count_accepted": int(len(accepted)),
                "acceptance_rate": float(len(accepted) / max(len(sub), 1)),
                "mean_log_H_all": float(sub["log_H"].mean()),
                "mean_log_H_accepted": float(accepted["log_H"].mean()) if len(accepted) else float("nan"),
                "mean_log_H_rejected": float(rejected["log_H"].mean()) if len(rejected) else float("nan"),
                "accepted_minus_all": float(accepted["log_H"].mean() - sub["log_H"].mean()) if len(accepted) else float("nan"),
                "accepted_minus_rejected": float(accepted["log_H"].mean() - rejected["log_H"].mean()) if len(accepted) and len(rejected) else float("nan"),
                "max_log_H_all": float(sub["log_H"].max()),
                "max_log_H_accepted": float(accepted["log_H"].max()) if len(accepted) else float("nan"),
            }
        )
    return pd.DataFrame(rows).sort_values(["n", "family"])


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Check whether the first-sieve window preferentially keeps high-entropy MLR samples.")
    parser.add_argument("--config", default="config_window_entropy_bias_check.yaml", help="Path to YAML config file.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_cfg = config["reference"]
    exp_cfg = config["experiment"]

    rows = []
    window_rows = []
    for n in exp_cfg["n_values"]:
        width_window, comp_window, _ = reference_window(
            family=str(ref_cfg["family"]),
            n=int(n),
            samples=int(ref_cfg["samples"]),
            lower_q=float(ref_cfg["lower_quantile"]),
            upper_q=float(ref_cfg["upper_quantile"]),
            seed_base=int(ref_cfg["seed_base"]),
        )
        window_rows.append(
            {
                "n": int(n),
                "width_lower": width_window[0],
                "width_upper": width_window[1],
                "comp_lower": comp_window[0],
                "comp_upper": comp_window[1],
            }
        )
        for family, n_samples in exp_cfg["families"].items():
            df_family = sample_family(
                family=str(family),
                n=int(n),
                n_samples=int(n_samples),
                width_window=width_window,
                comp_window=comp_window,
                seed_base=int(exp_cfg["seed_base"]) + 100000 * list(exp_cfg["families"].keys()).index(family),
            )
            rows.append(df_family)
            print(
                f"n={n:<3d} {family:20s} accepted={int(df_family['accepted_by_window'].sum()):<2d}/{len(df_family):<2d} "
                f"mean_logH_all={df_family['log_H'].mean():.3f}"
            )

    raw_df = pd.concat(rows, ignore_index=True)
    summary_df = summarize_bias(raw_df)
    windows_df = pd.DataFrame(window_rows)

    raw_df.to_csv(out_dir / "window_entropy_bias_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "window_entropy_bias_summary.csv", index=False, encoding="utf-8-sig")
    windows_df.to_csv(out_dir / "window_entropy_bias_windows.csv", index=False, encoding="utf-8-sig")

    print()
    print(summary_df.to_string(index=False))
