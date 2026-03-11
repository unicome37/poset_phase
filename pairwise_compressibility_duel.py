from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd
import yaml

from action import action_value
from coarse_grain import coarse_grain_delete_nodes
from entropy_exact import count_linear_extensions_exact
from experiment import FAMILIES
from observables import antichain_width, comparable_fraction, neutral_penalty
from observables_geo import geometric_components, geometric_penalty
from stability import (
    coarse_grain_penalty,
    family_centroids,
    nearest_family,
    self_drift,
    signature_dict,
)


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def quantile_window(values: pd.Series, lower_q: float, upper_q: float) -> tuple[float, float]:
    return float(values.quantile(lower_q)), float(values.quantile(upper_q))


def reference_window(
    family: str,
    n: int,
    samples: int,
    lower_q: float,
    upper_q: float,
    seed_base: int,
) -> tuple[tuple[float, float], tuple[float, float], pd.DataFrame]:
    generator = FAMILIES[family]
    rows = []
    for sample_id in range(samples):
        seed = seed_base + 1000 * n + sample_id
        poset = generator(n=n, seed=seed)
        rows.append(
            {
                "family": family,
                "n": n,
                "seed": seed,
                "antichain_width": antichain_width(poset),
                "comparable_fraction": comparable_fraction(poset),
            }
        )
    df = pd.DataFrame(rows)
    return (
        quantile_window(df["antichain_width"], lower_q, upper_q),
        quantile_window(df["comparable_fraction"], lower_q, upper_q),
        df,
    )


def accept(poset, width_window: tuple[float, float], comp_window: tuple[float, float]) -> bool:
    width = antichain_width(poset)
    comp = comparable_fraction(poset)
    return width_window[0] <= width <= width_window[1] and comp_window[0] <= comp <= comp_window[1]


def build_centroids(n: int, families: tuple[str, ...], samples_per_family: int, seed_base: int) -> dict[tuple[int, str], object]:
    rows = []
    for family in families:
        generator = FAMILIES[family]
        for sample_id in range(samples_per_family):
            seed = seed_base + 10000 * n + 100 * families.index(family) + sample_id
            poset = generator(n=n, seed=seed)
            row = {"n": n, "family": family}
            row.update(signature_dict(poset))
            rows.append(row)
    df = pd.DataFrame(rows)
    return family_centroids(df)


def coarse_grain_stats(
    poset,
    family: str,
    n: int,
    centroids,
    keep_ratios: tuple[float, ...],
    repeats: int,
) -> dict[str, float]:
    before = signature_dict(poset)
    drift_values = []
    family_switches = []
    penalty_values = []

    for keep_ratio in keep_ratios:
        for repeat in range(repeats):
            seed = 700000 + 10000 * n + 1000 * int(round(keep_ratio * 100)) + repeat
            cg_poset = coarse_grain_delete_nodes(poset, keep_ratio=keep_ratio, seed=seed)
            after = signature_dict(cg_poset)
            drift = self_drift(before, after)
            nearest, _ = nearest_family(after, centroids, n_value=n)
            switch = 0.0 if nearest == family else 1.0
            penalty = coarse_grain_penalty(
                drift_value=drift,
                family_switch_penalty=switch,
                rank_shift_penalty=0.0,
            )
            drift_values.append(drift)
            family_switches.append(switch)
            penalty_values.append(penalty)

    return {
        "cg_mean_self_drift": float(sum(drift_values) / max(len(drift_values), 1)),
        "cg_family_switch_rate": float(sum(family_switches) / max(len(family_switches), 1)),
        "cg_mean_penalty": float(sum(penalty_values) / max(len(penalty_values), 1)),
    }


def collect_family_samples(
    family: str,
    n: int,
    accepted_target: int,
    max_attempts: int,
    width_window: tuple[float, float],
    comp_window: tuple[float, float],
    beta: float,
    gamma: float,
    centroids,
    keep_ratios: tuple[float, ...],
    cg_repeats: int,
    seed_base: int,
) -> tuple[pd.DataFrame, dict]:
    generator = FAMILIES[family]
    rows = []
    accepted = 0
    attempts = 0

    while attempts < max_attempts and accepted < accepted_target:
        seed = seed_base + 10000 * n + attempts
        poset = generator(n=n, seed=seed)
        attempts += 1
        if not accept(poset, width_window, comp_window):
            continue

        count = count_linear_extensions_exact(poset)
        log_h = math.log(count)
        pen_neutral = neutral_penalty(poset)
        pen_geo = geometric_penalty(poset)
        geo = geometric_components(poset)
        cg = coarse_grain_stats(
            poset,
            family=family,
            n=n,
            centroids=centroids,
            keep_ratios=keep_ratios,
            repeats=cg_repeats,
        )
        row = {
            "family": family,
            "n": n,
            "seed": seed,
            "log_H": log_h,
            "antichain_width": antichain_width(poset),
            "comparable_fraction": comparable_fraction(poset),
            "penalty_neutral": pen_neutral,
            "penalty_geometric": pen_geo,
            "score_A1_gamma": action_value(log_h, pen_neutral, beta=beta, gamma=gamma),
            "score_A2_gamma": action_value(log_h, pen_neutral + pen_geo, beta=beta, gamma=gamma),
        }
        row.update(geo)
        row.update(cg)
        rows.append(row)
        accepted += 1

    summary = {
        "family": family,
        "n": n,
        "attempts": attempts,
        "accepted": accepted,
        "acceptance_rate": float(accepted / attempts) if attempts else 0.0,
    }
    return pd.DataFrame(rows), summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Pairwise Lor2D vs multi-layer compressibility duel.")
    parser.add_argument("--config", default="config_pairwise_compressibility_duel.yaml", help="Path to YAML config file.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)

    exp_cfg = config["experiment"]
    ref_cfg = config["reference"]
    cg_cfg = config["coarse_grain"]
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    beta = float(exp_cfg["beta"])
    gamma = float(exp_cfg["gamma"])
    n_values = tuple(int(v) for v in exp_cfg["n_values"])
    families = tuple(str(v) for v in exp_cfg["families"])
    accepted_target = int(exp_cfg["accepted_target"])
    max_attempts = int(exp_cfg["max_attempts"])
    width_q = (float(ref_cfg["lower_quantile"]), float(ref_cfg["upper_quantile"]))
    ref_samples = int(ref_cfg["samples"])
    ref_family = str(ref_cfg["family"])
    ref_seed_base = int(ref_cfg.get("seed_base", 100000))
    seed_base = int(exp_cfg.get("seed_base", 200000))

    keep_ratios = tuple(float(v) for v in cg_cfg["keep_ratios"])
    cg_repeats = int(cg_cfg["repeats"])
    centroid_samples = int(cg_cfg.get("centroid_samples", 8))
    centroid_seed_base = int(cg_cfg.get("centroid_seed_base", 900000))

    raw_rows = []
    summary_rows = []
    windows_rows = []

    for n in n_values:
        width_window, comp_window, ref_df = reference_window(
            family=ref_family,
            n=n,
            samples=ref_samples,
            lower_q=width_q[0],
            upper_q=width_q[1],
            seed_base=ref_seed_base,
        )
        windows_rows.append(
            {
                "n": n,
                "width_lower": width_window[0],
                "width_upper": width_window[1],
                "comp_lower": comp_window[0],
                "comp_upper": comp_window[1],
            }
        )
        centroids = build_centroids(
            n=n,
            families=families,
            samples_per_family=centroid_samples,
            seed_base=centroid_seed_base,
        )
        for family in families:
            df_family, summary = collect_family_samples(
                family=family,
                n=n,
                accepted_target=accepted_target,
                max_attempts=max_attempts,
                width_window=width_window,
                comp_window=comp_window,
                beta=beta,
                gamma=gamma,
                centroids=centroids,
                keep_ratios=keep_ratios,
                cg_repeats=cg_repeats,
                seed_base=seed_base + 1000 * families.index(family),
            )
            if not df_family.empty:
                raw_rows.append(df_family)
            summary_rows.append(summary)
            print(
                f"n={n:<3d} {family:20s} accepted={summary['accepted']:<2d}/{summary['attempts']:<3d} "
                f"rate={summary['acceptance_rate']:.3f}"
            )

    raw_df = pd.concat(raw_rows, ignore_index=True) if raw_rows else pd.DataFrame()
    summary_df = pd.DataFrame(summary_rows)
    windows_df = pd.DataFrame(windows_rows)
    duel_df = pd.DataFrame()
    if not raw_df.empty:
        duel_df = (
            raw_df.groupby(["n", "family"])
            .agg(
                mean_score_A2=("score_A2_gamma", "mean"),
                mean_log_H=("log_H", "mean"),
                mean_geo_total=("geo_total", "mean"),
                mean_geo_interval_shape=("geo_interval_shape", "mean"),
                mean_geo_interval_profile=("geo_interval_profile", "mean"),
                mean_geo_dim_proxy=("geo_dim_proxy_penalty", "mean"),
                mean_cg_penalty=("cg_mean_penalty", "mean"),
                mean_cg_self_drift=("cg_mean_self_drift", "mean"),
                mean_cg_switch_rate=("cg_family_switch_rate", "mean"),
                mean_width=("antichain_width", "mean"),
                mean_comp=("comparable_fraction", "mean"),
                count=("score_A2_gamma", "count"),
            )
            .reset_index()
            .sort_values(["n", "mean_score_A2"])
        )

    raw_df.to_csv(out_dir / "pairwise_duel_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "pairwise_duel_summary.csv", index=False, encoding="utf-8-sig")
    windows_df.to_csv(out_dir / "pairwise_duel_windows.csv", index=False, encoding="utf-8-sig")
    duel_df.to_csv(out_dir / "pairwise_duel_rankings.csv", index=False, encoding="utf-8-sig")

    if not duel_df.empty:
        print()
        print(duel_df.to_string(index=False))
