from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from action import action_value
from experiment import FAMILIES
from observables import neutral_penalty
from observables_geo import geometric_components
from runtime_utils import estimate_entropy_by_family


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def estimate_gamma_crossing(summary_df: pd.DataFrame, scale_label: str, left: str, right: str) -> pd.DataFrame:
    rows = []
    sub = summary_df[summary_df["scale_label"] == scale_label].copy()
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
                "scale_label": scale_label,
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
        description="Scan consistency weight for width_height + dim_consistency skeleton."
    )
    parser.add_argument("--config", default="config_width_height_consistency_scale_scan.yaml", help="Path to YAML config.")
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
    scan_factors = tuple(float(x) for x in exp_cfg["scale_factors"])
    sis_runs = int(exp_cfg.get("sis_runs", 512))
    default_exact_threshold = int(exp_cfg.get("exact_threshold", 44))
    family_exact_thresholds = {
        str(k): int(v) for k, v in (exp_cfg.get("family_exact_thresholds", {}) or {}).items()
    }

    base_rows = []
    for n in n_values:
        for family in families:
            generator = FAMILIES[family]
            for sample_id in range(samples_per_family):
                seed = 1000 * n + sample_id
                poset = generator(n=n, seed=seed)
                log_h, entropy_method = estimate_entropy_by_family(
                    poset,
                    family=family,
                    sis_runs=sis_runs,
                    seed=seed,
                    default_exact_threshold=default_exact_threshold,
                    family_exact_thresholds=family_exact_thresholds,
                )
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
                        "entropy_method": entropy_method,
                        "penalty_neutral": neutral,
                        **geo,
                    }
                )

    base_df = pd.DataFrame(base_rows)
    fixed_weight = exp_cfg.get("fixed_matched_weight")
    if fixed_weight is not None:
        matched_weight = float(fixed_weight)
    else:
        mean_dim_proxy = float(base_df["geo_dim_proxy_penalty"].mean())
        mean_dim_consistency = float(base_df["geo_dim_consistency"].mean())
        matched_weight = 8.0 * mean_dim_proxy / mean_dim_consistency if mean_dim_consistency > 0 else 8.0

    rows = []
    for row in base_df.to_dict(orient="records"):
        for factor in scan_factors:
            scale = matched_weight * factor
            scale_label = f"{factor:.3f}"
            penalty = row["penalty_neutral"] + 2.0 * row["geo_width_height"] + scale * row["geo_dim_consistency"]
            for gamma in gammas:
                rows.append(
                    {
                        **row,
                        "gamma": gamma,
                        "scale_factor": factor,
                        "scale_label": scale_label,
                        "matched_weight": matched_weight,
                        "consistency_weight": scale,
                        "penalty_effective": penalty,
                        "score": action_value(row["log_H"], penalty, beta=beta, gamma=gamma),
                    }
                )

    raw_df = pd.DataFrame(rows)
    summary_df = (
        raw_df.groupby(["scale_label", "scale_factor", "n", "gamma", "family"])
        .agg(
            mean_score=("score", "mean"),
            std_score=("score", "std"),
            mean_penalty=("penalty_effective", "mean"),
            mean_log_H=("log_H", "mean"),
            entropy_method=("entropy_method", lambda s: ",".join(sorted(set(str(x) for x in s)))),
            count=("score", "count"),
        )
        .reset_index()
        .sort_values(["scale_factor", "n", "gamma", "mean_score"])
    )

    report_df = pd.concat(
        [estimate_gamma_crossing(summary_df, label, "lorentzian_like_2d", "KR_like") for label in summary_df["scale_label"].drop_duplicates()],
        ignore_index=True,
    )

    raw_df.to_csv(out_dir / "width_height_consistency_scale_scan_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "width_height_consistency_scale_scan_summary.csv", index=False, encoding="utf-8-sig")
    report_df.to_csv(out_dir / "width_height_consistency_scale_scan_report.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "width_height_consistency_scale_scan_report.csv").as_posix())
    print()
    print(f"matched_weight={matched_weight:.6f}")
    print()
    print(report_df.to_string(index=False))
