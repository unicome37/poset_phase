from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from experiment import FAMILIES
from pairwise_compressibility_duel import reference_window
from stability import SIGNATURE_COLUMNS, family_centroids, nearest_family, signature_dict


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_reference_centroids(
    n: int,
    families: tuple[str, ...],
    samples: int,
    seed_base: int,
) -> dict[tuple[int, str], np.ndarray]:
    rows = []
    for family in families:
        generator = FAMILIES[family]
        for sample_id in range(samples):
            seed = seed_base + 10000 * n + 100 * families.index(family) + sample_id
            poset = generator(n=n, seed=seed)
            row = {"n": n, "family": family}
            row.update(signature_dict(poset))
            rows.append(row)
    return family_centroids(pd.DataFrame(rows))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Check whether MLR survivors are closer to Lor2D in signature space.")
    parser.add_argument("--config", default="config_mlr_survivor_nearest_family_check.yaml", help="Path to YAML config.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    profile_csv = Path(config["input"]["profile_raw_csv"])
    profile_df = pd.read_csv(profile_csv)

    exp_cfg = config["experiment"]
    families = tuple(str(v) for v in exp_cfg["families"])
    samples = int(exp_cfg["centroid_samples"])
    seed_base = int(exp_cfg["centroid_seed_base"])

    rows = []
    for n, sub in profile_df.groupby("n"):
        centroids = build_reference_centroids(
            n=int(n),
            families=families,
            samples=samples,
            seed_base=seed_base,
        )
        lor_centroid = centroids[(int(n), "lorentzian_like_2d")]
        mlr_centroid = centroids[(int(n), "multi_layer_random")]
        for row in sub.itertuples(index=False):
            sig = {
                "sig_comp": float(row.comparable_fraction),
                "sig_d_eff": float(row.geo_dim_eff),
                "sig_height_ratio": float(row.height_ratio),
                "sig_width_ratio": float(row.width_ratio),
                "sig_degree_var": float(row.degree_var_norm),
                "sig_interval_empty": float(row.interval_empty_fraction),
                "sig_interval_size_ratio": float(row.interval_size_ratio),
            }
            vec = np.asarray([sig[c] for c in SIGNATURE_COLUMNS], dtype=float)
            dist_lor = float(np.linalg.norm(vec - lor_centroid))
            dist_mlr = float(np.linalg.norm(vec - mlr_centroid))
            nearest, nearest_dist = nearest_family(sig, centroids, n_value=int(n))
            rows.append(
                {
                    "n": int(n),
                    "seed": int(row.seed),
                    "group": str(row.group),
                    "dist_to_lor2d": dist_lor,
                    "dist_to_mlr": dist_mlr,
                    "dist_delta_mlr_minus_lor2d": dist_mlr - dist_lor,
                    "nearest_family": nearest,
                    "nearest_dist": nearest_dist,
                }
            )

    out_df = pd.DataFrame(rows)
    summary_df = (
        out_df.groupby(["n", "group"])
        .agg(
            mean_dist_to_lor2d=("dist_to_lor2d", "mean"),
            mean_dist_to_mlr=("dist_to_mlr", "mean"),
            mean_dist_delta=("dist_delta_mlr_minus_lor2d", "mean"),
            lor2d_nearest_rate=("nearest_family", lambda s: float((s == "lorentzian_like_2d").mean())),
            mlr_nearest_rate=("nearest_family", lambda s: float((s == "multi_layer_random").mean())),
            count=("seed", "count"),
        )
        .reset_index()
    )

    out_df.to_csv(out_dir / "mlr_survivor_nearest_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "mlr_survivor_nearest_summary.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "mlr_survivor_nearest_raw.csv").as_posix())
    print((out_dir / "mlr_survivor_nearest_summary.csv").as_posix())
    print()
    print(summary_df.to_string(index=False))
