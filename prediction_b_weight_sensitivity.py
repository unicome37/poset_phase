from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from action import action_value
from entropy_exact import log_linear_extensions_exact
from experiment import FAMILIES
from observables import neutral_penalty
from observables_geo import DEFAULT_GEOMETRIC_WEIGHTS, geometric_components, geometric_penalty_from_components


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def estimate_gamma_crossing(summary_df: pd.DataFrame, left: str, right: str) -> pd.DataFrame:
    rows = []
    grouped = summary_df.groupby(["scan_name", "weight_name", "weight_value", "n"], sort=True)
    for (scan_name, weight_name, weight_value, n), group in grouped:
        left_df = group[group["family"] == left][["gamma", "mean_score"]].rename(columns={"mean_score": "left_score"})
        right_df = group[group["family"] == right][["gamma", "mean_score"]].rename(columns={"mean_score": "right_score"})
        merged = left_df.merge(right_df, on="gamma", how="inner").sort_values("gamma")
        if merged.empty:
            continue

        merged["delta_score_left_minus_right"] = merged["left_score"] - merged["right_score"]
        gamma_c = None
        status = "no_cross_left_never_wins"
        deltas = merged["delta_score_left_minus_right"].tolist()

        if all(d <= 0 for d in deltas):
            status = "no_cross_left_always_wins"
        elif any(d < 0 for d in deltas):
            status = "no_cross_left_sometimes_wins"

        for i in range(len(merged) - 1):
            d1 = float(merged.iloc[i]["delta_score_left_minus_right"])
            d2 = float(merged.iloc[i + 1]["delta_score_left_minus_right"])
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
                "scan_name": scan_name,
                "weight_name": weight_name,
                "weight_value": float(weight_value),
                "n": int(n),
                "left_family": left,
                "right_family": right,
                "gamma_c_est": gamma_c,
                "status": status,
                "min_delta_score": float(merged["delta_score_left_minus_right"].min()),
                "max_delta_score": float(merged["delta_score_left_minus_right"].max()),
            }
        )
    return pd.DataFrame(rows)


def build_scan_specs(config: dict) -> list[dict[str, object]]:
    scans = []
    scan_cfg = config["experiment"]["weight_scans"]
    for weight_name, values in scan_cfg.items():
        for value in values:
            weights = dict(DEFAULT_GEOMETRIC_WEIGHTS)
            weights[weight_name] = float(value)
            scans.append(
                {
                    "scan_name": f"{weight_name}={float(value):g}",
                    "weight_name": weight_name,
                    "weight_value": float(value),
                    "weights": weights,
                }
            )
    return scans


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prediction B weight sensitivity scan for gamma_c.")
    parser.add_argument(
        "--config",
        default="config_prediction_b_weight_sensitivity.yaml",
        help="Path to YAML config.",
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    exp_cfg = config["experiment"]
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    n_values = tuple(int(n) for n in exp_cfg["n_values"])
    gammas = tuple(float(g) for g in exp_cfg["gammas"])
    beta = float(exp_cfg.get("beta", 1.0))
    samples_per_family = int(exp_cfg.get("samples_per_family", 1))
    families = tuple(exp_cfg["families"])
    scans = build_scan_specs(config)

    rows = []
    for n in n_values:
        for family in families:
            generator = FAMILIES[family]
            for sample_id in range(samples_per_family):
                seed = 1000 * n + sample_id
                poset = generator(n=n, seed=seed)
                log_h = log_linear_extensions_exact(poset)
                neutral = neutral_penalty(poset)
                geo = geometric_components(poset)

                for scan in scans:
                    penalty = neutral + geometric_penalty_from_components(geo, weights=scan["weights"])
                    for gamma in gammas:
                        rows.append(
                            {
                                "n": n,
                                "family": family,
                                "sample_id": sample_id,
                                "seed": seed,
                                "gamma": gamma,
                                "beta": beta,
                                "log_H": log_h,
                                "penalty_neutral": neutral,
                                "penalty_geometric_full": float(geo["geo_total"]),
                                "penalty_effective": penalty,
                                "scan_name": scan["scan_name"],
                                "weight_name": scan["weight_name"],
                                "weight_value": scan["weight_value"],
                                "score": action_value(log_h, penalty, beta=beta, gamma=gamma),
                                **geo,
                            }
                        )

    raw_df = pd.DataFrame(rows)
    summary_df = (
        raw_df.groupby(["scan_name", "weight_name", "weight_value", "n", "gamma", "family"])
        .agg(
            mean_score=("score", "mean"),
            std_score=("score", "std"),
            mean_log_H=("log_H", "mean"),
            mean_penalty=("penalty_effective", "mean"),
            count=("score", "count"),
        )
        .reset_index()
        .sort_values(["weight_name", "weight_value", "n", "gamma", "mean_score"])
    )

    report_df = estimate_gamma_crossing(summary_df, "lorentzian_like_2d", "KR_like")

    raw_df.to_csv(out_dir / "prediction_b_weight_sensitivity_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "prediction_b_weight_sensitivity_summary.csv", index=False, encoding="utf-8-sig")
    report_df.to_csv(out_dir / "prediction_b_weight_sensitivity_report.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "prediction_b_weight_sensitivity_report.csv").as_posix())
    print()
    print(report_df.to_string(index=False))
