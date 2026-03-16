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


VARIANT_DROPS = {
    "A2_full": [],
    "A1_neutral_only": list(DEFAULT_GEOMETRIC_WEIGHTS.keys()),
    "drop_geo_width_height": ["geo_width_height"],
    "drop_geo_dim_proxy_penalty": ["geo_dim_proxy_penalty"],
    "drop_geo_comparability_window": ["geo_comparability_window"],
    "drop_geo_cover_density": ["geo_cover_density"],
    "drop_geo_interval_profile": ["geo_interval_profile"],
    "drop_geo_interval_shape": ["geo_interval_shape"],
    "drop_geo_layer_smoothness": ["geo_layer_smoothness"],
    "drop_width_height_plus_dim": ["geo_width_height", "geo_dim_proxy_penalty"],
    "drop_all_interval_terms": ["geo_interval_profile", "geo_interval_shape"],
    "interval_terms_only": [
        "geo_width_height",
        "geo_dim_proxy_penalty",
        "geo_comparability_window",
        "geo_cover_density",
        "geo_layer_smoothness",
    ],
    "non_dim_non_width_only": [
        "geo_width_height",
        "geo_dim_proxy_penalty",
    ],
}


def resolve_variants(config: dict) -> dict[str, dict[str, object]]:
    variant_cfg = config["experiment"].get("variants")
    if not variant_cfg:
        return {
            name: {"drops": list(drops), "weight_overrides": {}}
            for name, drops in VARIANT_DROPS.items()
        }

    resolved: dict[str, dict[str, object]] = {}
    for name, spec in variant_cfg.items():
        resolved[name] = {
            "drops": list(spec.get("drops", [])),
            "weight_overrides": dict(spec.get("weight_overrides", {})),
        }
    return resolved


def build_variant_weights(drops: list[str], weight_overrides: dict[str, float]) -> dict[str, float]:
    weights = dict(DEFAULT_GEOMETRIC_WEIGHTS)
    for key in drops:
        weights[key] = 0.0
    for key, value in weight_overrides.items():
        weights[key] = float(value)
    return weights


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
        merged["delta_score_left_minus_right"] = merged["left_score"] - merged["right_score"]
        gamma_c = None
        status = "no_cross_left_never_wins"
        deltas = merged["delta_score_left_minus_right"].tolist()
        gammas = merged["gamma"].tolist()

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
                "variant": variant,
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


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run A2 geometric-term ablation and estimate gamma_c.")
    parser.add_argument("--config", default="config_geometric_ablation_gamma_c.yaml", help="Path to YAML config.")
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
    variant_specs = resolve_variants(config)
    variants = list(variant_specs.keys())

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
                geo_total = float(geo["geo_total"])

                variant_penalties = {}
                for variant_name, variant_spec in variant_specs.items():
                    weights = build_variant_weights(
                        drops=variant_spec["drops"],
                        weight_overrides=variant_spec["weight_overrides"],
                    )
                    variant_penalties[variant_name] = neutral + geometric_penalty_from_components(geo, weights=weights)

                for gamma in gammas:
                    for variant_name, penalty in variant_penalties.items():
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
                                "penalty_geometric_full": geo_total,
                                "penalty_effective": penalty,
                                "variant": variant_name,
                                "score": action_value(log_h, penalty, beta=beta, gamma=gamma),
                                **geo,
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
            count=("score", "count"),
        )
        .reset_index()
        .sort_values(["variant", "n", "gamma", "mean_score"])
    )

    reports = [
        report
        for variant_name in variants
        for report in [estimate_gamma_crossing(summary_df, variant_name, "lorentzian_like_2d", "KR_like")]
        if not report.empty
    ]
    report_df = pd.concat(reports, ignore_index=True) if reports else pd.DataFrame()

    raw_df.to_csv(out_dir / "geometric_ablation_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "geometric_ablation_summary.csv", index=False, encoding="utf-8-sig")
    report_df.to_csv(out_dir / "geometric_ablation_gamma_c_report.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "geometric_ablation_gamma_c_report.csv").as_posix())
    print()
    print(report_df.to_string(index=False))
