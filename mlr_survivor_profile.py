from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from experiment import FAMILIES
from observables import (
    antichain_width,
    comparable_fraction,
    degree_stats,
    layer_imbalance,
    layer_profile,
    normalized_degree_variance,
)
from observables_geo import (
    geometric_components,
    height_ratio,
    interval_empty_fraction,
    interval_shape_penalty,
    mean_interval_size_ratio,
    width_ratio,
)
from pairwise_compressibility_duel import reference_window


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def accept(width: int, comp: float, width_window: tuple[float, float], comp_window: tuple[float, float]) -> bool:
    return width_window[0] <= width <= width_window[1] and comp_window[0] <= comp <= comp_window[1]


def window_distance(width: int, comp: float, width_window: tuple[float, float], comp_window: tuple[float, float]) -> float:
    width_center = 0.5 * (width_window[0] + width_window[1])
    comp_center = 0.5 * (comp_window[0] + comp_window[1])
    width_scale = max(width_window[1] - width_window[0], 1.0)
    comp_scale = max(comp_window[1] - comp_window[0], 1e-6)
    return float(abs(width - width_center) / width_scale + abs(comp - comp_center) / comp_scale)


def profile_metrics(poset, width_window: tuple[float, float], comp_window: tuple[float, float]) -> dict[str, float]:
    layers = layer_profile(poset)
    stats = degree_stats(poset)
    width = antichain_width(poset)
    comp = comparable_fraction(poset)
    geo = geometric_components(poset)
    return {
        "antichain_width": float(width),
        "comparable_fraction": float(comp),
        "window_distance": window_distance(width, comp, width_window, comp_window),
        "layer_count": float(len(layers)),
        "max_layer_size": float(layers.max()),
        "mean_layer_size": float(layers.mean()),
        "layer_cv": float(layers.std(ddof=0) / max(layers.mean(), 1e-12)),
        "layer_imbalance": float(layer_imbalance(poset)),
        "height_ratio": float(height_ratio(poset)),
        "width_ratio": float(width_ratio(poset)),
        "degree_var_norm": float(normalized_degree_variance(poset)),
        "in_mean": float(stats["in_mean"]),
        "out_mean": float(stats["out_mean"]),
        "interval_empty_fraction": float(interval_empty_fraction(poset)),
        "interval_size_ratio": float(mean_interval_size_ratio(poset)),
        "interval_shape_penalty": float(interval_shape_penalty(poset)),
        "geo_dim_eff": float(geo["geo_dim_eff"]),
        "geo_width_height": float(geo["geo_width_height"]),
        "geo_dim_proxy_penalty": float(geo["geo_dim_proxy_penalty"]),
        "geo_comparability_window": float(geo["geo_comparability_window"]),
        "geo_cover_density": float(geo["geo_cover_density"]),
        "geo_interval_profile": float(geo["geo_interval_profile"]),
        "geo_interval_shape": float(geo["geo_interval_shape"]),
        "geo_layer_smoothness": float(geo["geo_layer_smoothness"]),
        "geo_total": float(geo["geo_total"]),
    }


def collect_profiles(
    n: int,
    accepted_target: int,
    rejected_target: int,
    max_attempts: int,
    width_window: tuple[float, float],
    comp_window: tuple[float, float],
    seed_base: int,
) -> pd.DataFrame:
    generator = FAMILIES["multi_layer_random"]
    rows: list[dict] = []
    accepted = 0
    rejected = 0
    attempts = 0

    while attempts < max_attempts and (accepted < accepted_target or rejected < rejected_target):
        seed = seed_base + 10000 * n + attempts
        poset = generator(n=n, seed=seed)
        attempts += 1
        width = antichain_width(poset)
        comp = comparable_fraction(poset)
        is_accepted = accept(width, comp, width_window, comp_window)
        if is_accepted and accepted >= accepted_target:
            continue
        if (not is_accepted) and rejected >= rejected_target:
            continue

        row = {
            "n": int(n),
            "seed": int(seed),
            "group": "accepted" if is_accepted else "rejected",
        }
        row.update(profile_metrics(poset, width_window, comp_window))
        rows.append(row)

        if is_accepted:
            accepted += 1
        else:
            rejected += 1

    return pd.DataFrame(rows)


def summarize_profiles(df: pd.DataFrame) -> pd.DataFrame:
    metric_cols = [c for c in df.columns if c not in {"n", "seed", "group"}]
    summary = (
        df.groupby(["n", "group"])[metric_cols]
        .mean()
        .reset_index()
    )
    return summary


def contrast_profiles(df: pd.DataFrame) -> pd.DataFrame:
    metric_cols = [c for c in df.columns if c not in {"n", "seed", "group"}]
    rows: list[dict] = []
    for n, sub in df.groupby("n"):
        acc = sub[sub["group"] == "accepted"]
        rej = sub[sub["group"] == "rejected"]
        if acc.empty or rej.empty:
            continue
        row = {"n": int(n), "count_accepted": int(len(acc)), "count_rejected": int(len(rej))}
        for col in metric_cols:
            row[f"{col}_accepted_mean"] = float(acc[col].mean())
            row[f"{col}_rejected_mean"] = float(rej[col].mean())
            row[f"{col}_delta"] = float(acc[col].mean() - rej[col].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def top_deltas(contrast_df: pd.DataFrame, top_k: int = 8) -> pd.DataFrame:
    rows: list[dict] = []
    for row in contrast_df.itertuples(index=False):
        n = int(row.n)
        deltas = []
        for col, value in row._asdict().items():
            if col.endswith("_delta"):
                deltas.append((col[:-6], abs(float(value)), float(value)))
        deltas.sort(key=lambda x: x[1], reverse=True)
        for metric, abs_delta, signed_delta in deltas[:top_k]:
            rows.append({"n": n, "metric": metric, "abs_delta": abs_delta, "signed_delta": signed_delta})
    return pd.DataFrame(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Profile MLR survivors vs ordinary rejected samples.")
    parser.add_argument("--config", default="config_mlr_survivor_profile.yaml", help="Path to YAML config file.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_cfg = config["reference"]
    exp_cfg = config["experiment"]

    raw_parts = []
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
        df_n = collect_profiles(
            n=int(n),
            accepted_target=int(exp_cfg["accepted_target"]),
            rejected_target=int(exp_cfg["rejected_target"]),
            max_attempts=int(exp_cfg["max_attempts"]),
            width_window=width_window,
            comp_window=comp_window,
            seed_base=int(exp_cfg["seed_base"]),
        )
        raw_parts.append(df_n)
        print(
            f"n={n:<3d} accepted={int((df_n['group']=='accepted').sum()):<2d} "
            f"rejected={int((df_n['group']=='rejected').sum()):<2d}"
        )

    raw_df = pd.concat(raw_parts, ignore_index=True)
    summary_df = summarize_profiles(raw_df)
    contrast_df = contrast_profiles(raw_df)
    top_df = top_deltas(contrast_df, top_k=int(config["output"].get("top_k", 8)))
    windows_df = pd.DataFrame(window_rows)

    raw_df.to_csv(out_dir / "mlr_survivor_profile_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "mlr_survivor_profile_summary.csv", index=False, encoding="utf-8-sig")
    contrast_df.to_csv(out_dir / "mlr_survivor_profile_contrast.csv", index=False, encoding="utf-8-sig")
    top_df.to_csv(out_dir / "mlr_survivor_profile_top_deltas.csv", index=False, encoding="utf-8-sig")
    windows_df.to_csv(out_dir / "mlr_survivor_profile_windows.csv", index=False, encoding="utf-8-sig")

    print()
    print(top_df.to_string(index=False))
