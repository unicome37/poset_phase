from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from experiment import FAMILIES
from observables import layer_profile
from observables_geo import cover_density


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def layer_index_by_minima(closure: np.ndarray) -> np.ndarray:
    indeg = closure.sum(axis=0).astype(int)
    remaining = set(range(closure.shape[0]))
    layer_idx = np.full(closure.shape[0], -1, dtype=int)
    current_layer = 0
    while remaining:
        mins = [i for i in remaining if indeg[i] == 0]
        if not mins:
            break
        for u in mins:
            layer_idx[u] = current_layer
        for u in mins:
            remaining.remove(u)
            for v in np.where(closure[u])[0]:
                indeg[v] -= 1
        current_layer += 1
    return layer_idx


def transitive_reduction_edge_count(closure: np.ndarray) -> int:
    n = closure.shape[0]
    count = 0
    for i in range(n):
        succ = np.where(closure[i])[0]
        for j in succ:
            mid = np.where(closure[i] & closure[:, j])[0]
            if len(mid) == 0:
                count += 1
    return count


def layer_connection_metrics(closure: np.ndarray, layer_idx: np.ndarray) -> dict[str, float]:
    adjacent_edges = 0
    long_edges = 0
    total_edges = int(closure.sum())
    if total_edges == 0:
        return {
            "adjacent_edge_fraction": 0.0,
            "long_edge_fraction": 0.0,
            "mean_layer_gap": 0.0,
        }
    gaps = []
    for i, j in np.argwhere(closure):
        gap = int(layer_idx[j] - layer_idx[i])
        if gap == 1:
            adjacent_edges += 1
        elif gap > 1:
            long_edges += 1
        gaps.append(gap)
    return {
        "adjacent_edge_fraction": float(adjacent_edges / total_edges),
        "long_edge_fraction": float(long_edges / total_edges),
        "mean_layer_gap": float(np.mean(gaps)),
    }


def layer_signature_redundancy(closure: np.ndarray, layer_idx: np.ndarray) -> float:
    signatures = []
    for layer in sorted(set(layer_idx.tolist())):
        nodes = np.where(layer_idx == layer)[0]
        if len(nodes) <= 1:
            continue
        succ_sign = [tuple(np.where(closure[u])[0].tolist()) for u in nodes]
        pred_sign = [tuple(np.where(closure[:, u])[0].tolist()) for u in nodes]
        sigs = list(zip(pred_sign, succ_sign))
        unique = len(set(sigs))
        redundancy = 1.0 - unique / len(sigs)
        signatures.append(redundancy)
    if not signatures:
        return 0.0
    return float(np.mean(signatures))


def residual_metrics(poset) -> dict[str, float]:
    closure = poset.closure
    layers = layer_profile(poset)
    layer_idx = layer_index_by_minima(closure)
    red_edges = transitive_reduction_edge_count(closure)
    total_pairs = poset.n * (poset.n - 1) / 2
    conn = layer_connection_metrics(closure, layer_idx)
    return {
        "layer_count": float(len(layers)),
        "mean_layer_size": float(layers.mean()),
        "layer_size_std": float(layers.std(ddof=0)),
        "cover_density": float(cover_density(poset)),
        "reduction_edge_count": float(red_edges),
        "reduction_edge_density": float(red_edges / total_pairs) if total_pairs else 0.0,
        "layer_signature_redundancy": float(layer_signature_redundancy(closure, layer_idx)),
        **conn,
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compare residual freedom metrics on matched MLR vs Lor2D pairs.")
    parser.add_argument("--config", default="config_matched_residual_freedom_check.yaml", help="Path to YAML config.")
    return parser


if __name__ == "__main__":
    args = build_arg_parser().parse_args()
    config = load_config(args.config)
    out_dir = Path(config["output"]["directory"])
    out_dir.mkdir(parents=True, exist_ok=True)

    pairs_df = pd.read_csv(config["input"]["pairs_csv"])
    rows = []
    for row in pairs_df.itertuples(index=False):
        n = int(row.n)
        for family, seed in [("multi_layer_random", int(row.mlr_seed)), ("lorentzian_like_2d", int(row.lor_seed))]:
            poset = FAMILIES[family](n=n, seed=seed)
            metrics = residual_metrics(poset)
            out = {"n": n, "pair_mlr_seed": int(row.mlr_seed), "pair_lor_seed": int(row.lor_seed), "family": family, "seed": seed}
            out.update(metrics)
            rows.append(out)

    raw_df = pd.DataFrame(rows)
    summary_df = (
        raw_df.groupby(["n", "family"])
        .agg(
            mean_cover_density=("cover_density", "mean"),
            mean_reduction_edge_density=("reduction_edge_density", "mean"),
            mean_adjacent_edge_fraction=("adjacent_edge_fraction", "mean"),
            mean_long_edge_fraction=("long_edge_fraction", "mean"),
            mean_layer_gap=("mean_layer_gap", "mean"),
            mean_signature_redundancy=("layer_signature_redundancy", "mean"),
            mean_layer_count=("layer_count", "mean"),
            mean_layer_size_std=("layer_size_std", "mean"),
            count=("seed", "count"),
        )
        .reset_index()
    )

    contrast_rows = []
    for n, sub in summary_df.groupby("n"):
        mlr = sub[sub["family"] == "multi_layer_random"].iloc[0]
        lor = sub[sub["family"] == "lorentzian_like_2d"].iloc[0]
        contrast_rows.append(
            {
                "n": int(n),
                "cover_density_delta": float(mlr["mean_cover_density"] - lor["mean_cover_density"]),
                "reduction_edge_density_delta": float(mlr["mean_reduction_edge_density"] - lor["mean_reduction_edge_density"]),
                "adjacent_edge_fraction_delta": float(mlr["mean_adjacent_edge_fraction"] - lor["mean_adjacent_edge_fraction"]),
                "long_edge_fraction_delta": float(mlr["mean_long_edge_fraction"] - lor["mean_long_edge_fraction"]),
                "mean_layer_gap_delta": float(mlr["mean_layer_gap"] - lor["mean_layer_gap"]),
                "signature_redundancy_delta": float(mlr["mean_signature_redundancy"] - lor["mean_signature_redundancy"]),
                "layer_count_delta": float(mlr["mean_layer_count"] - lor["mean_layer_count"]),
                "layer_size_std_delta": float(mlr["mean_layer_size_std"] - lor["mean_layer_size_std"]),
            }
        )
    contrast_df = pd.DataFrame(contrast_rows)

    raw_df.to_csv(out_dir / "matched_residual_freedom_raw.csv", index=False, encoding="utf-8-sig")
    summary_df.to_csv(out_dir / "matched_residual_freedom_summary.csv", index=False, encoding="utf-8-sig")
    contrast_df.to_csv(out_dir / "matched_residual_freedom_contrast.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "matched_residual_freedom_summary.csv").as_posix())
    print((out_dir / "matched_residual_freedom_contrast.csv").as_posix())
    print()
    print(contrast_df.to_string(index=False))
