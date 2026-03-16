from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from experiment import FAMILIES
from mlr_survivor_profile import accept, profile_metrics
from pairwise_compressibility_duel import reference_window
from stability import SIGNATURE_COLUMNS, family_centroids, nearest_family, signature_dict


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def select_winner_families(winners_csv: Path, family_prefix: str, min_wins: int) -> pd.DataFrame:
    df = pd.read_csv(winners_csv)
    winner_counts = (
        df[df["winner_family"].astype(str).str.startswith(family_prefix)]
        .groupby("winner_family")
        .size()
        .reset_index(name="winner_count")
        .sort_values(["winner_count", "winner_family"], ascending=[False, True])
    )
    winner_counts = winner_counts[winner_counts["winner_count"] >= min_wins].reset_index(drop=True)
    return winner_counts


def build_reference_centroids(
    n: int,
    families: tuple[str, ...],
    samples: int,
    seed_base: int,
) -> dict[tuple[int, str], np.ndarray]:
    rows = []
    for family_idx, family in enumerate(families):
        generator = FAMILIES[family]
        for sample_id in range(samples):
            seed = seed_base + 10000 * n + 100 * family_idx + sample_id
            poset = generator(n=n, seed=seed)
            row = {"n": int(n), "family": family}
            row.update(signature_dict(poset))
            rows.append(row)
    return family_centroids(pd.DataFrame(rows))


def collect_family_samples(
    family: str,
    n: int,
    samples_per_family: int,
    seed_base: int,
    width_window: tuple[float, float],
    comp_window: tuple[float, float],
    centroids: dict[tuple[int, str], np.ndarray],
) -> pd.DataFrame:
    generator = FAMILIES[family]
    rows: list[dict[str, float | int | str]] = []
    for sample_id in range(samples_per_family):
        seed = seed_base + 100000 * n + 1000 * sorted(FAMILIES.keys()).index(family) + sample_id
        poset = generator(n=n, seed=seed)
        metrics = profile_metrics(poset, width_window, comp_window)
        sig = signature_dict(poset)
        nearest, nearest_dist = nearest_family(sig, centroids, n_value=int(n))
        dist_to_lor2d = float(
            np.linalg.norm(
                np.asarray([sig[c] for c in SIGNATURE_COLUMNS], dtype=float)
                - centroids[(int(n), "lorentzian_like_2d")]
            )
        )
        row = {
            "family": family,
            "n": int(n),
            "sample_id": int(sample_id),
            "seed": int(seed),
            "accepted_lor2d_window": int(
                accept(
                    int(metrics["antichain_width"]),
                    float(metrics["comparable_fraction"]),
                    width_window,
                    comp_window,
                )
            ),
            "nearest_family": nearest,
            "nearest_dist": float(nearest_dist),
            "dist_to_lor2d": dist_to_lor2d,
        }
        row.update(metrics)
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_acceptance(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["family", "n"])
        .agg(
            samples=("seed", "count"),
            accepted=("accepted_lor2d_window", "sum"),
            acceptance_rate=("accepted_lor2d_window", "mean"),
            mean_window_distance=("window_distance", "mean"),
            accepted_mean_window_distance=(
                "window_distance",
                lambda s: float(s[df.loc[s.index, "accepted_lor2d_window"] == 1].mean())
                if (df.loc[s.index, "accepted_lor2d_window"] == 1).any()
                else float("nan"),
            ),
        )
        .reset_index()
    )


def summarize_nearest(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for (family, n), sub in df.groupby(["family", "n"]):
        for group_name, group_df in {
            "all": sub,
            "accepted": sub[sub["accepted_lor2d_window"] == 1],
            "rejected": sub[sub["accepted_lor2d_window"] == 0],
        }.items():
            if group_df.empty:
                continue
            row: dict[str, float | int | str] = {
                "family": family,
                "n": int(n),
                "group": group_name,
                "count": int(len(group_df)),
                "mean_dist_to_lor2d": float(group_df["dist_to_lor2d"].mean()),
            }
            nearest_rate = (
                group_df["nearest_family"].value_counts(normalize=True).sort_index().to_dict()
            )
            for nearest_family_name, rate in nearest_rate.items():
                row[f"nearest_rate__{nearest_family_name}"] = float(rate)
            rows.append(row)
    return pd.DataFrame(rows)


def contrast_accepted_vs_rejected(df: pd.DataFrame) -> pd.DataFrame:
    metric_cols = [
        c
        for c in df.columns
        if c
        not in {
            "family",
            "n",
            "sample_id",
            "seed",
            "accepted_lor2d_window",
            "nearest_family",
            "nearest_dist",
        }
    ]
    rows: list[dict[str, float | int | str]] = []
    for (family, n), sub in df.groupby(["family", "n"]):
        acc = sub[sub["accepted_lor2d_window"] == 1]
        rej = sub[sub["accepted_lor2d_window"] == 0]
        if acc.empty or rej.empty:
            continue
        row: dict[str, float | int | str] = {
            "family": family,
            "n": int(n),
            "count_accepted": int(len(acc)),
            "count_rejected": int(len(rej)),
        }
        for col in metric_cols:
            row[f"{col}_accepted_mean"] = float(acc[col].mean())
            row[f"{col}_rejected_mean"] = float(rej[col].mean())
            row[f"{col}_delta"] = float(acc[col].mean() - rej[col].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Check whether layered winners form Lorentzian-like anomalous subgroups."
    )
    parser.add_argument(
        "--config",
        default="config_prediction_b_layered_winner_subgroup_check.yaml",
        help="Path to YAML config file.",
    )
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)

    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    winner_counts = select_winner_families(
        winners_csv=Path(config["input"]["winners_csv"]),
        family_prefix=str(config["selection"]["winner_family_prefix"]),
        min_wins=int(config["selection"]["min_wins"]),
    )
    active_families = tuple(winner_counts["winner_family"].tolist())
    reference_families = tuple(config["selection"].get("reference_families", []))
    centroid_families = tuple(dict.fromkeys(reference_families + active_families))

    ref_cfg = config["reference"]
    exp_cfg = config["experiment"]

    all_rows = []
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
                "width_lower": float(width_window[0]),
                "width_upper": float(width_window[1]),
                "comp_lower": float(comp_window[0]),
                "comp_upper": float(comp_window[1]),
            }
        )
        centroids = build_reference_centroids(
            n=int(n),
            families=centroid_families,
            samples=int(exp_cfg["centroid_samples"]),
            seed_base=int(exp_cfg["centroid_seed_base"]),
        )
        for family in active_families:
            df_family = collect_family_samples(
                family=family,
                n=int(n),
                samples_per_family=int(exp_cfg["samples_per_family"]),
                seed_base=int(exp_cfg["seed_base"]),
                width_window=width_window,
                comp_window=comp_window,
                centroids=centroids,
            )
            all_rows.append(df_family)
            accepted = int(df_family["accepted_lor2d_window"].sum())
            print(
                f"n={n:<3d} family={family:<30s} accepted={accepted:<3d} "
                f"rate={accepted / max(len(df_family), 1):.3f}",
                flush=True,
            )

    raw_df = pd.concat(all_rows, ignore_index=True)
    acceptance_df = summarize_acceptance(raw_df)
    nearest_df = summarize_nearest(raw_df)
    contrast_df = contrast_accepted_vs_rejected(raw_df)
    windows_df = pd.DataFrame(window_rows)

    winner_counts.to_csv(
        out_dir / "layered_winner_family_counts.csv",
        index=False,
        encoding="utf-8-sig",
    )
    raw_df.to_csv(
        out_dir / "layered_winner_subgroup_raw.csv",
        index=False,
        encoding="utf-8-sig",
    )
    acceptance_df.to_csv(
        out_dir / "layered_winner_subgroup_acceptance.csv",
        index=False,
        encoding="utf-8-sig",
    )
    nearest_df.to_csv(
        out_dir / "layered_winner_subgroup_nearest_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    contrast_df.to_csv(
        out_dir / "layered_winner_subgroup_contrast.csv",
        index=False,
        encoding="utf-8-sig",
    )
    windows_df.to_csv(
        out_dir / "layered_winner_subgroup_windows.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print()
    print((out_dir / "layered_winner_subgroup_acceptance.csv").as_posix())
    print((out_dir / "layered_winner_subgroup_nearest_summary.csv").as_posix())
