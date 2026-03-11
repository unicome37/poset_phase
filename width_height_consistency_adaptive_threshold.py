from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from entropy_exact import log_linear_extensions_exact
from experiment import FAMILIES
from observables import neutral_penalty
from observables_geo import geometric_components


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_base_df(n_values: tuple[int, ...], samples_per_family: int) -> pd.DataFrame:
    rows = []
    for n in n_values:
        for family in ("KR_like", "lorentzian_like_2d"):
            generator = FAMILIES[family]
            for sample_id in range(samples_per_family):
                seed = 1000 * n + sample_id
                poset = generator(n=n, seed=seed)
                log_h = log_linear_extensions_exact(poset)
                neutral = neutral_penalty(poset)
                geo = geometric_components(poset)
                rows.append(
                    {
                        "n": n,
                        "family": family,
                        "sample_id": sample_id,
                        "seed": seed,
                        "log_H": log_h,
                        "penalty_neutral": neutral,
                        "geo_width_height": geo["geo_width_height"],
                        "geo_dim_proxy_penalty": geo["geo_dim_proxy_penalty"],
                        "geo_dim_consistency": geo["geo_dim_consistency"],
                    }
                )
    return pd.DataFrame(rows)


def mean_delta_score(base_df: pd.DataFrame, n: int, factor: float, matched_weight: float) -> tuple[float, float, float]:
    sub = base_df[base_df["n"] == n].copy()
    sub["penalty_effective"] = (
        sub["penalty_neutral"] + 2.0 * sub["geo_width_height"] + (matched_weight * factor) * sub["geo_dim_consistency"]
    )
    agg = (
        sub.groupby("family")
        .agg(mean_log_h=("log_H", "mean"), mean_penalty=("penalty_effective", "mean"))
        .reset_index()
    )
    lor = agg[agg["family"] == "lorentzian_like_2d"].iloc[0]
    kr = agg[agg["family"] == "KR_like"].iloc[0]
    h_lor = float(lor["mean_log_h"])
    h_kr = float(kr["mean_log_h"])
    p_lor = float(lor["mean_penalty"])
    p_kr = float(kr["mean_penalty"])
    delta_h = h_lor - h_kr
    delta_p = p_lor - p_kr
    gamma_c = (delta_h / delta_p) if (delta_p != 0) else float("nan")
    return delta_h, delta_p, gamma_c


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Adaptive threshold search for width_height + dim_consistency factor c_N.")
    parser.add_argument(
        "--config",
        default="config_width_height_consistency_adaptive_threshold.yaml",
        help="Path to YAML config.",
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    exp_cfg = config["experiment"]
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    n_values = tuple(exp_cfg["n_values"])
    samples_per_family = int(exp_cfg.get("samples_per_family", 2))
    factor_min = float(exp_cfg.get("factor_min", 0.2))
    factor_max = float(exp_cfg.get("factor_max", 3.0))
    coarse_step = float(exp_cfg.get("coarse_step", 0.1))
    refine_step = float(exp_cfg.get("refine_step", 0.01))

    base_df = build_base_df(n_values, samples_per_family)
    mean_dim_proxy = float(base_df["geo_dim_proxy_penalty"].mean())
    mean_dim_consistency = float(base_df["geo_dim_consistency"].mean())
    matched_weight = 8.0 * mean_dim_proxy / mean_dim_consistency if mean_dim_consistency > 0 else 8.0

    threshold_rows = []
    scan_rows = []
    for n in n_values:
        coarse_factors = []
        f = factor_min
        while f <= factor_max + 1e-12:
            coarse_factors.append(round(f, 10))
            f += coarse_step

        found_low = None
        found_high = None
        for factor in coarse_factors:
            delta_h, delta_p, gamma_c = mean_delta_score(base_df, n, factor, matched_weight)
            crossing = (delta_h * delta_p > 0) and (delta_h / delta_p > 0)
            scan_rows.append(
                {
                    "stage": "coarse",
                    "n": n,
                    "factor": factor,
                    "delta_h": delta_h,
                    "delta_p": delta_p,
                    "gamma_c_est": gamma_c,
                    "crossing_possible": crossing,
                }
            )
            if crossing:
                found_high = factor
                break
            found_low = factor

        min_factor = None
        min_gamma_c = None
        if found_high is not None:
            start = max(factor_min, found_high - coarse_step)
            refine_factors = []
            f = start
            while f <= found_high + 1e-12:
                refine_factors.append(round(f, 10))
                f += refine_step
            for factor in refine_factors:
                delta_h, delta_p, gamma_c = mean_delta_score(base_df, n, factor, matched_weight)
                crossing = (delta_h * delta_p > 0) and (delta_h / delta_p > 0)
                scan_rows.append(
                    {
                        "stage": "refine",
                        "n": n,
                        "factor": factor,
                        "delta_h": delta_h,
                        "delta_p": delta_p,
                        "gamma_c_est": gamma_c,
                        "crossing_possible": crossing,
                    }
                )
                if crossing:
                    min_factor = factor
                    min_gamma_c = gamma_c
                    break

        threshold_rows.append(
            {
                "n": n,
                "matched_weight": matched_weight,
                "min_restoring_factor": min_factor,
                "gamma_c_at_min_factor": min_gamma_c,
            }
        )

    scan_df = pd.DataFrame(scan_rows)
    threshold_df = pd.DataFrame(threshold_rows)
    base_df.to_csv(out_dir / "width_height_consistency_adaptive_base.csv", index=False, encoding="utf-8-sig")
    scan_df.to_csv(out_dir / "width_height_consistency_adaptive_scan.csv", index=False, encoding="utf-8-sig")
    threshold_df.to_csv(out_dir / "width_height_consistency_adaptive_thresholds.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "width_height_consistency_adaptive_thresholds.csv").as_posix())
    print()
    print(f"matched_weight={matched_weight:.6f}")
    print()
    print(threshold_df.to_string(index=False))
