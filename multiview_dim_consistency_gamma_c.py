from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from action import action_value
from entropy_exact import log_linear_extensions_exact
from experiment import FAMILIES
from observables import neutral_penalty
from observables_geo import geometric_components


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def estimate_gamma_crossing(summary_df: pd.DataFrame, variant: str, left: str, right: str) -> pd.DataFrame:
    rows = []
    sub = summary_df[summary_df["variant"] == variant].copy()
    for n, group in sub.groupby("n", sort=True):
        left_df = group[group["family"] == left][["gamma", "mean_score"]].rename(columns={"mean_score": "left_score"})
        right_df = group[group["family"] == right][["gamma", "mean_score"]].rename(columns={"mean_score": "right_score"})
        merged = left_df.merge(right_df, on="gamma", how="inner").sort_values("gamma")
        if merged.empty:
            continue
        merged["delta"] = merged["left_score"] - merged["right_score"]
        gamma_c = None
        status = "no_cross_left_never_wins"
        if all(float(d) <= 0.0 for d in merged["delta"]):
            status = "no_cross_left_always_wins"
        elif any(float(d) < 0.0 for d in merged["delta"]):
            status = "no_cross_left_sometimes_wins"
        for i in range(len(merged) - 1):
            d1 = float(merged.iloc[i]["delta"])
            d2 = float(merged.iloc[i + 1]["delta"])
            g1 = float(merged.iloc[i]["gamma"])
            g2 = float(merged.iloc[i + 1]["gamma"])
            if d1 == 0.0:
                gamma_c = g1
                status = "crossing"
                break
            if d1 * d2 < 0:
                gamma_c = g1 + (0.0 - d1) * (g2 - g1) / (d2 - d1)
                status = "crossing"
                break
        rows.append(
            {
                "variant": variant,
                "n": int(n),
                "gamma_c_est": gamma_c,
                "status": status,
                "min_delta_score": float(merged["delta"].min()),
                "max_delta_score": float(merged["delta"].max()),
            }
        )
    return pd.DataFrame(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Test multi-estimator dimension consistency as a replacement for dim_proxy."
    )
    parser.add_argument("--config", default="config_multiview_dim_consistency_gamma_c.yaml", help="Path to YAML config.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    exp_cfg = config["experiment"]
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    n_values = tuple(exp_cfg["n_values"])
    gammas = tuple(float(g) for g in exp_cfg["gammas"])
    beta = float(exp_cfg.get("beta", 1.0))
    samples_per_family = int(exp_cfg.get("samples_per_family", 2))
    families = tuple(exp_cfg["families"])

    base_rows = []
    for n in n_values:
        for family in families:
            generator = FAMILIES[family]
            for sample_id in range(samples_per_family):
                seed = 1000 * n + sample_id
                poset = generator(n=n, seed=seed)
                log_h = log_linear_extensions_exact(poset)
                neutral = neutral_penalty(poset)
                geo = geometric_components(poset)
                base_rows.append(
                    {
                        "n": n,
                        "family": family,
                        "sample_id": sample_id,
                        "seed": seed,
                        "beta": beta,
                        "log_H": log_h,
                        "penalty_neutral": neutral,
                        **geo,
                    }
                )

    base_df = pd.DataFrame(base_rows)
    mean_dim_proxy = float(base_df["geo_dim_proxy_penalty"].mean())
    mean_dim_consistency = float(base_df["geo_dim_consistency"].mean())
    mean_dim_multi = float(base_df["geo_dim_multi_consistency"].mean())
    scale_consistency = 8.0 * mean_dim_proxy / mean_dim_consistency if mean_dim_consistency > 0 else 8.0
    scale_multi = 8.0 * mean_dim_proxy / mean_dim_multi if mean_dim_multi > 0 else 8.0

    variants = {
        "A2_full": lambda row: row["penalty_neutral"] + row["geo_total"],
        "replace_dim_with_consistency_scale_matched": lambda row: (
            row["penalty_neutral"]
            + row["geo_total"]
            - 8.0 * row["geo_dim_proxy_penalty"]
            + scale_consistency * row["geo_dim_consistency"]
        ),
        "replace_dim_with_multi_consistency_scale_matched": lambda row: (
            row["penalty_neutral"]
            + row["geo_total"]
            - 8.0 * row["geo_dim_proxy_penalty"]
            + scale_multi * row["geo_dim_multi_consistency"]
        ),
        "width_height_plus_multi_consistency_only": lambda row: (
            row["penalty_neutral"]
            + 2.0 * row["geo_width_height"]
            + scale_multi * row["geo_dim_multi_consistency"]
        ),
    }

    rows = []
    for row in base_df.to_dict(orient="records"):
        for gamma in gammas:
            for variant_name, penalty_fn in variants.items():
                penalty = float(penalty_fn(row))
                rows.append(
                    {
                        **row,
                        "gamma": gamma,
                        "variant": variant_name,
                        "scale_consistency": scale_consistency,
                        "scale_multi": scale_multi,
                        "score": action_value(row["log_H"], penalty, beta=beta, gamma=gamma),
                        "penalty_effective": penalty,
                    }
                )

    raw_df = pd.DataFrame(rows)
    summary_df = (
        raw_df.groupby(["variant", "n", "gamma", "family"])
        .agg(
            mean_score=("score", "mean"),
            std_score=("score", "std"),
            mean_log_H=("log_H", "mean"),
            mean_penalty=("penalty_effective", "mean"),
            mean_dim_consistency=("geo_dim_consistency", "mean"),
            mean_dim_multi_consistency=("geo_dim_multi_consistency", "mean"),
            mean_d_order=("geo_d_order", "mean"),
            mean_d_chain=("geo_d_chain", "mean"),
            mean_d_width=("geo_d_width", "mean"),
            count=("score", "count"),
        )
        .reset_index()
        .sort_values(["variant", "n", "gamma", "mean_score"])
    )

    report_df = pd.concat(
        [estimate_gamma_crossing(summary_df, variant_name, "lorentzian_like_2d", "KR_like") for variant_name in variants],
        ignore_index=True,
    )

    family_views_df = (
        base_df.groupby(["n", "family"])
        .agg(
            mean_dim_consistency=("geo_dim_consistency", "mean"),
            mean_dim_multi_consistency=("geo_dim_multi_consistency", "mean"),
            mean_d_order=("geo_d_order", "mean"),
            mean_d_chain=("geo_d_chain", "mean"),
            mean_d_width=("geo_d_width", "mean"),
            count=("geo_d_order", "count"),
        )
        .reset_index()
        .sort_values(["n", "family"])
    )

    raw_df.to_csv(out_dir / "multiview_dim_consistency_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "multiview_dim_consistency_summary.csv", index=False, encoding="utf-8-sig")
    report_df.to_csv(out_dir / "multiview_dim_consistency_gamma_c_report.csv", index=False, encoding="utf-8-sig")
    family_views_df.to_csv(out_dir / "multiview_dim_consistency_family_views.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "multiview_dim_consistency_gamma_c_report.csv").as_posix())
    print()
    print(f"scale_consistency={scale_consistency:.6f}")
    print(f"scale_multi={scale_multi:.6f}")
    print()
    print(report_df.to_string(index=False))
