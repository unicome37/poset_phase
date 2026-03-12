from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from action import action_value
from coarse_grain import coarse_grain_delete_nodes
from entropy_exact import count_linear_extensions_exact
from experiment import FAMILIES
from observables import antichain_width, comparable_fraction, neutral_penalty
from observables_geo import geometric_components, geometric_penalty
from stability import (
    family_centroids,
    nearest_family,
    signature_dict,
    signature_vector_from_dict,
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
) -> tuple[tuple[float, float], tuple[float, float]]:
    generator = FAMILIES[family]
    rows = []
    for sample_id in range(samples):
        seed = seed_base + 1000 * n + sample_id
        poset = generator(n=n, seed=seed)
        rows.append(
            {
                "antichain_width": antichain_width(poset),
                "comparable_fraction": comparable_fraction(poset),
            }
        )
    df = pd.DataFrame(rows)
    return (
        quantile_window(df["antichain_width"], lower_q, upper_q),
        quantile_window(df["comparable_fraction"], lower_q, upper_q),
    )


def accept(poset, width_window: tuple[float, float], comp_window: tuple[float, float]) -> bool:
    width = antichain_width(poset)
    comp = comparable_fraction(poset)
    return width_window[0] <= width <= width_window[1] and comp_window[0] <= comp <= comp_window[1]


def build_centroids(n: int, families: tuple[str, ...], samples_per_family: int, seed_base: int) -> dict[tuple[int, str], np.ndarray]:
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


def collect_family_samples(
    family: str,
    n: int,
    accepted_target: int,
    max_attempts: int,
    width_window: tuple[float, float],
    comp_window: tuple[float, float],
    beta: float,
    gamma: float,
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
        sig = signature_dict(poset)
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
            "poset": poset,
            "sig_dict": sig,
        }
        row.update(geo)
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


def blind_identity_stats(
    sample_row: pd.Series,
    all_original: pd.DataFrame,
    centroids: dict[tuple[int, str], np.ndarray],
    keep_ratios: tuple[float, ...],
    repeats: int,
) -> dict[str, float]:
    poset = sample_row["poset"]
    family = str(sample_row["family"])
    n = int(sample_row["n"])
    seed = int(sample_row["seed"])
    before_sig = sample_row["sig_dict"]
    before_vec = signature_vector_from_dict(before_sig)

    original_pool = []
    for row in all_original.itertuples(index=False):
        if int(row.n) != n:
            continue
        original_pool.append(
            {
                "family": str(row.family),
                "seed": int(row.seed),
                "vec": signature_vector_from_dict(row.sig_dict),
            }
        )

    drift_values = []
    family_switches = []
    blind_hits = []
    blind_knn3_hits = []
    blind_knn5_hits = []

    for keep_ratio in keep_ratios:
        for repeat in range(repeats):
            cg_seed = 700000 + 10000 * n + 1000 * int(round(keep_ratio * 100)) + repeat
            cg_poset = coarse_grain_delete_nodes(poset, keep_ratio=keep_ratio, seed=cg_seed)
            after_sig = signature_dict(cg_poset)
            after_vec = signature_vector_from_dict(after_sig)

            drift = float(np.mean((after_vec - before_vec) ** 2))
            nearest_family_name, _ = nearest_family(after_sig, centroids, n_value=n)
            family_switch = 0.0 if nearest_family_name == family else 1.0

            best_seed = None
            best_dist = float("inf")
            dist_rows = []
            for cand in original_pool:
                dist = float(np.linalg.norm(after_vec - cand["vec"]))
                dist_rows.append((dist, cand["family"], cand["seed"]))
                if dist < best_dist:
                    best_dist = dist
                    best_seed = cand["seed"]

            blind_hit = 1.0 if best_seed == seed else 0.0
            dist_rows.sort(key=lambda x: x[0])
            knn3 = [fam for _dist, fam, _seed in dist_rows[:3]]
            knn5 = [fam for _dist, fam, _seed in dist_rows[:5]]
            blind_knn3 = 1.0 if knn3.count(family) >= 2 else 0.0
            blind_knn5 = 1.0 if knn5.count(family) >= 3 else 0.0

            drift_values.append(drift)
            family_switches.append(family_switch)
            blind_hits.append(blind_hit)
            blind_knn3_hits.append(blind_knn3)
            blind_knn5_hits.append(blind_knn5)

    return {
        "cg_mean_self_drift": float(np.mean(drift_values)),
        "cg_family_switch_rate": float(np.mean(family_switches)),
        "cg_blind_self_hit_rate": float(np.mean(blind_hits)),
        "cg_blind_knn3_family_hit_rate": float(np.mean(blind_knn3_hits)),
        "cg_blind_knn5_family_hit_rate": float(np.mean(blind_knn5_hits)),
    }


def zscore_by_n(df: pd.DataFrame, value_col: str) -> pd.Series:
    out = []
    for _n, sub in df.groupby("n", sort=False):
        mean = float(sub[value_col].mean())
        std = float(sub[value_col].std(ddof=0))
        std = std if std > 1e-12 else 1.0
        out.extend(((sub[value_col] - mean) / std).tolist())
    return pd.Series(out, index=df.index)


def build_scan(raw_df: pd.DataFrame, zetas: list[float], variant: str) -> pd.DataFrame:
    df = raw_df.copy()
    if variant == "switch_zscore":
        df["identity_variant"] = zscore_by_n(df, "cg_family_switch_rate")
    elif variant == "blind_self_hit_zscore":
        df["identity_variant"] = -zscore_by_n(df, "cg_blind_self_hit_rate")
    elif variant == "blind_knn3_zscore":
        df["identity_variant"] = -zscore_by_n(df, "cg_blind_knn3_family_hit_rate")
    elif variant == "blind_knn5_zscore":
        df["identity_variant"] = -zscore_by_n(df, "cg_blind_knn5_family_hit_rate")
    elif variant == "blind_self_hit_centered":
        df["identity_variant"] = -(df.groupby("n")["cg_blind_self_hit_rate"].transform(lambda s: s - float(s.mean())))
    else:
        raise ValueError(f"Unsupported variant: {variant}")

    rows = []
    for zeta in zetas:
        sub = df.copy()
        sub["score_eval"] = sub["score_A2_gamma"] + zeta * sub["identity_variant"]
        summary = (
            sub.groupby(["n", "family"])
            .agg(
                mean_score=("score_eval", "mean"),
                mean_identity_variant=("identity_variant", "mean"),
                mean_score_A2=("score_A2_gamma", "mean"),
                count=("score_A2_gamma", "count"),
            )
            .reset_index()
        )
        summary["zeta"] = float(zeta)
        summary["variant"] = variant
        rows.append(summary)
    return pd.concat(rows, ignore_index=True)


def crossing_report(scan_df: pd.DataFrame) -> pd.DataFrame:
    pivot = scan_df.pivot_table(
        index=["variant", "n", "zeta"],
        columns="family",
        values="mean_score",
    ).reset_index()
    pivot["delta_mlr_minus_lor2d"] = pivot["multi_layer_random"] - pivot["lorentzian_like_2d"]

    rows = []
    for (variant, n), sub in pivot.groupby(["variant", "n"]):
        sub = sub.sort_values("zeta").reset_index(drop=True)
        zeta_cross = None
        for i in range(1, len(sub)):
            y0 = float(sub.loc[i - 1, "delta_mlr_minus_lor2d"])
            y1 = float(sub.loc[i, "delta_mlr_minus_lor2d"])
            if y0 < 0.0 <= y1:
                x0 = float(sub.loc[i - 1, "zeta"])
                x1 = float(sub.loc[i, "zeta"])
                frac = (-y0) / max(y1 - y0, 1e-12)
                zeta_cross = x0 + frac * (x1 - x0)
                break
        rows.append(
            {
                "variant": str(variant),
                "n": int(n),
                "zeta_cross": zeta_cross,
                "delta_at_min_zeta": float(sub.iloc[0]["delta_mlr_minus_lor2d"]),
                "delta_at_max_zeta": float(sub.iloc[-1]["delta_mlr_minus_lor2d"]),
            }
        )
    return pd.DataFrame(rows)


def plot_scan(scan_df: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(
        nrows=len(sorted(scan_df["n"].unique())),
        ncols=1,
        figsize=(8, 3.4 * len(sorted(scan_df["n"].unique()))),
        squeeze=False,
    )
    for ax, (n, sub) in zip(axes[:, 0], scan_df.groupby("n")):
        for (variant, family), fam_df in sub.groupby(["variant", "family"]):
            fam_df = fam_df.sort_values("zeta")
            ax.plot(fam_df["zeta"], fam_df["mean_score"], marker="o", label=f"{family} [{variant}]")
        ax.set_title(f"Blind identity scan, N={n}")
        ax.set_xlabel("zeta")
        ax.set_ylabel("mean score")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Blind self-identity scan without family centroids.")
    parser.add_argument("--config", default="config_pairwise_blind_identity_scan.yaml", help="Path to YAML config.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    exp_cfg = config["experiment"]
    ref_cfg = config["reference"]
    cg_cfg = config["coarse_grain"]
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    summaries = []

    for n in exp_cfg["n_values"]:
        width_window, comp_window = reference_window(
            family=ref_cfg["family"],
            n=int(n),
            samples=int(ref_cfg["samples"]),
            lower_q=float(ref_cfg["lower_quantile"]),
            upper_q=float(ref_cfg["upper_quantile"]),
            seed_base=int(ref_cfg["seed_base"]),
        )
        centroids = build_centroids(
            n=int(n),
            families=tuple(exp_cfg["families"]),
            samples_per_family=int(cg_cfg["centroid_samples"]),
            seed_base=int(cg_cfg["centroid_seed_base"]),
        )

        family_frames = []
        for family in exp_cfg["families"]:
            frame, summary = collect_family_samples(
                family=family,
                n=int(n),
                accepted_target=int(exp_cfg["accepted_target"]),
                max_attempts=int(exp_cfg["max_attempts"]),
                width_window=width_window,
                comp_window=comp_window,
                beta=float(exp_cfg["beta"]),
                gamma=float(exp_cfg["gamma"]),
                seed_base=int(exp_cfg["seed_base"]),
            )
            family_frames.append(frame)
            summaries.append(summary)

        all_original = pd.concat(family_frames, ignore_index=True)
        for _, sample_row in all_original.iterrows():
            cg = blind_identity_stats(
                sample_row,
                all_original=all_original,
                centroids=centroids,
                keep_ratios=tuple(cg_cfg["keep_ratios"]),
                repeats=int(cg_cfg["repeats"]),
            )
            row = {k: v for k, v in sample_row.items() if k not in {"poset", "sig_dict"}}
            row.update(cg)
            rows.append(row)

    raw_df = pd.DataFrame(rows)
    summary_df = pd.DataFrame(summaries)

    zetas = [float(v) for v in config["experiment"]["zetas"]]
    variants = [str(v) for v in config["experiment"]["variants"]]
    scan_df = pd.concat([build_scan(raw_df, zetas, variant=v) for v in variants], ignore_index=True)
    cross_df = crossing_report(scan_df)

    raw_df.to_csv(out_dir / "pairwise_blind_identity_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "pairwise_blind_identity_acceptance.csv", index=False, encoding="utf-8-sig")
    scan_df.to_csv(out_dir / "pairwise_blind_identity_scan.csv", index=False, encoding="utf-8-sig")
    cross_df.to_csv(out_dir / "pairwise_blind_identity_crossings.csv", index=False, encoding="utf-8-sig")
    plot_scan(scan_df, out_dir / "pairwise_blind_identity_scan.png")

    print((out_dir / "pairwise_blind_identity_raw.csv").as_posix())
    print((out_dir / "pairwise_blind_identity_crossings.csv").as_posix())
    print((out_dir / "pairwise_blind_identity_scan.png").as_posix())
    print()
    print(cross_df.to_string(index=False))
