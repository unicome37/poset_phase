from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from experiment import FAMILIES
from matched_residual_freedom_check import residual_metrics


BLOCK_COLS = [
    "source",
    "n",
    "n_cg",
    "keep_ratio",
    "gamma",
    "action_mode",
    "zeta",
    "icg_variant",
]

STRATUM_COLS = ["n", "keep_ratio", "gamma", "icg_variant"]


def _rank(x: np.ndarray) -> np.ndarray:
    # Fast average ranks, 1..n, like Spearman (handles ties by averaging).
    x = np.asarray(x, dtype=float)
    n = len(x)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1, dtype=float)
    xs = x[order]
    i = 0
    while i < n:
        j = i + 1
        while j < n and xs[j] == xs[i]:
            j += 1
        if j - i > 1:
            avg = 0.5 * ((i + 1) + j)  # ranks are 1-based
            ranks[order[i:j]] = avg
        i = j
    return ranks


def _residualize(y: np.ndarray, covs: list[np.ndarray]) -> np.ndarray:
    if not covs:
        return y - float(np.mean(y))
    X = np.column_stack([np.ones(len(y), dtype=float)] + covs)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    x0 = x - float(np.mean(x))
    y0 = y - float(np.mean(y))
    den = float(np.sqrt(np.sum(x0 * x0) * np.sum(y0 * y0)))
    if den <= 0:
        # If either vector is constant after residualization, treat as "no association" for this block.
        return 0.0
    return float(np.sum(x0 * y0) / den)


def _residualizer_matrix(n: int, covs: list[np.ndarray]) -> np.ndarray:
    if not covs:
        v = np.ones((n, 1), dtype=float)
        return np.eye(n, dtype=float) - (v @ v.T) / float(n)
    X = np.column_stack([np.ones(n, dtype=float)] + covs)
    H = X @ np.linalg.pinv(X)
    return np.eye(n, dtype=float) - H


def _corr_residualized(x: np.ndarray, y_res: np.ndarray, R_x: np.ndarray) -> float:
    xr = R_x @ x
    den = float(np.sqrt(np.sum(xr * xr) * np.sum(y_res * y_res)))
    if den <= 0:
        return 0.0
    return float(np.sum(xr * y_res) / den)


def partial_spearman(x: np.ndarray, y: np.ndarray, covs: list[np.ndarray]) -> float:
    xr = _rank(x)
    yr = _rank(y)
    covr = [_rank(c) for c in covs]
    xr_res = _residualize(xr, covr)
    yr_res = _residualize(yr, covr)
    return _pearson(xr_res, yr_res)


@dataclass(frozen=True)
class BlockStat:
    rho_raw: float
    rho_partial: float
    rho_partial_lc: float
    rho_partial_gap: float
    spearman_icg_lc: float
    spearman_icg_gap: float


def _block_stats(block_df: pd.DataFrame) -> BlockStat:
    x = block_df["mean_I_cg"].to_numpy(dtype=float)
    y = block_df["improve_rank"].to_numpy(dtype=float)
    lc = block_df["layer_count"].to_numpy(dtype=float)
    gap = block_df["mean_layer_gap"].to_numpy(dtype=float)

    rho_raw = partial_spearman(x, y, covs=[])
    rho_partial = partial_spearman(x, y, covs=[lc, gap])
    rho_partial_lc = partial_spearman(x, y, covs=[lc])
    rho_partial_gap = partial_spearman(x, y, covs=[gap])

    spearman_icg_lc = partial_spearman(x, lc, covs=[])
    spearman_icg_gap = partial_spearman(x, gap, covs=[])

    return BlockStat(
        rho_raw=rho_raw,
        rho_partial=rho_partial,
        rho_partial_lc=rho_partial_lc,
        rho_partial_gap=rho_partial_gap,
        spearman_icg_lc=spearman_icg_lc,
        spearman_icg_gap=spearman_icg_gap,
    )


def _perm_p_value(obs: float, perm_stats: np.ndarray) -> float:
    perm_stats = perm_stats[np.isfinite(perm_stats)]
    if len(perm_stats) == 0 or not np.isfinite(obs):
        return float("nan")
    # Two-sided permutation p-value with +1 correction.
    return float((np.sum(np.abs(perm_stats) >= abs(obs)) + 1.0) / (len(perm_stats) + 1.0))


def _block_perm_test(
    block_df: pd.DataFrame,
    n_perm: int,
    rng: np.random.Generator,
) -> tuple[BlockStat, dict[str, np.ndarray]]:
    obs = _block_stats(block_df)

    x = block_df["mean_I_cg"].to_numpy(dtype=float)
    y = block_df["improve_rank"].to_numpy(dtype=float)
    lc = block_df["layer_count"].to_numpy(dtype=float)
    gap = block_df["mean_layer_gap"].to_numpy(dtype=float)

    n = len(block_df)
    # Precompute rank-space covariates and residualizer matrices (since covariates + y are fixed under permutations).
    yr = _rank(y)
    lc_r = _rank(lc)
    gap_r = _rank(gap)

    R_none = _residualizer_matrix(n, covs=[])
    R_both = _residualizer_matrix(n, covs=[lc_r, gap_r])
    R_lc = _residualizer_matrix(n, covs=[lc_r])
    R_gap = _residualizer_matrix(n, covs=[gap_r])

    yr_res_none = R_none @ yr
    yr_res_both = R_both @ yr
    yr_res_lc = R_lc @ yr
    yr_res_gap = R_gap @ yr

    perm_raw = np.empty(n_perm, dtype=float)
    perm_partial = np.empty(n_perm, dtype=float)
    perm_partial_lc = np.empty(n_perm, dtype=float)
    perm_partial_gap = np.empty(n_perm, dtype=float)

    for i in range(n_perm):
        xp = rng.permutation(x)
        xpr = _rank(xp)
        perm_raw[i] = _corr_residualized(xpr, yr_res_none, R_none)
        perm_partial[i] = _corr_residualized(xpr, yr_res_both, R_both)
        perm_partial_lc[i] = _corr_residualized(xpr, yr_res_lc, R_lc)
        perm_partial_gap[i] = _corr_residualized(xpr, yr_res_gap, R_gap)

    return obs, {
        "raw": perm_raw,
        "partial": perm_partial,
        "partial_lc": perm_partial_lc,
        "partial_gap": perm_partial_gap,
    }


def parse_keys(keys: list[str]) -> set[tuple[int, float, float]]:
    out: set[tuple[int, float, float]] = set()
    for k in keys:
        parts = k.split(":")
        if len(parts) != 3:
            raise ValueError(f"Invalid key: {k!r} (expected 'n:keep:gamma')")
        n = int(parts[0])
        keep = float(parts[1])
        gamma = float(parts[2])
        out.add((n, keep, gamma))
    return out


def compute_layer_metrics_for_source(source_dir: Path, n_values: set[int] | None = None) -> pd.DataFrame:
    original_csv = source_dir / "cg_original_samples.csv"
    if not original_csv.exists():
        raise FileNotFoundError(f"Missing cg_original_samples.csv in source dir: {source_dir}")

    df = pd.read_csv(original_csv)
    if n_values is not None:
        df = df[df["n"].isin(sorted(n_values))].copy()

    rows = []
    for row in df[["n", "family", "sample_id", "seed_offset"]].drop_duplicates().itertuples(index=False):
        n = int(row.n)
        family = str(row.family)
        sample_id = int(row.sample_id)
        seed_offset = int(row.seed_offset)
        seed_base = seed_offset * 10_000_000
        seed = seed_base + 1000 * n + sample_id

        poset = FAMILIES[family](n=n, seed=seed)
        metrics = residual_metrics(poset)
        rows.append(
            {
                "source": source_dir.as_posix(),
                "n": n,
                "family": family,
                "sample_id": sample_id,
                "layer_count": float(metrics["layer_count"]),
                "mean_layer_gap": float(metrics["mean_layer_gap"]),
            }
        )

    raw = pd.DataFrame(rows)
    agg = (
        raw.groupby(["source", "n", "family"])
        .agg(
            layer_count=("layer_count", "mean"),
            mean_layer_gap=("mean_layer_gap", "mean"),
            n_samples=("sample_id", "count"),
        )
        .reset_index()
    )
    return raw, agg


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Prediction D: refine the block-permutation delta-rank test with layer_count/mean_layer_gap controls."
    )
    p.add_argument("--zeta-rank-csv", required=True, help="Path to cg_zeta_scan_rankings_variants.csv")
    p.add_argument(
        "--keys",
        required=True,
        nargs="+",
        help="Target strata as 'n:keep:gamma' (e.g. 30:0.6:0.2).",
    )
    p.add_argument(
        "--require-variants",
        default="full,switch,no_switch",
        help="Comma-separated icg_variants to include (others are dropped).",
    )
    p.add_argument("--n-perm", type=int, default=20000, help="Permutations per block (default: 20000).")
    p.add_argument("--seed", type=int, default=0, help="RNG seed for permutations.")
    p.add_argument("--out", required=True, help="Output directory.")
    return p


def main() -> None:
    args = build_arg_parser().parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    keys = parse_keys(args.keys)
    n_values = {n for (n, _, _) in keys}
    require_variants = [s.strip() for s in str(args.require_variants).split(",") if s.strip()]

    df = pd.read_csv(args.zeta_rank_csv)
    df = df[df["icg_variant"].isin(require_variants)].copy()
    df["n"] = df["n"].astype(int)
    df["n_cg"] = df["n_cg"].astype(int)
    df["keep_ratio"] = df["keep_ratio"].astype(float)
    df["gamma"] = df["gamma"].astype(float)

    df["improve_rank"] = -(df["rank_eval"].to_numpy(dtype=float) - df["rank_local_eval"].to_numpy(dtype=float))

    df = df[df.apply(lambda r: (int(r["n"]), float(r["keep_ratio"]), float(r["gamma"])) in keys, axis=1)].copy()
    if df.empty:
        raise SystemExit("No rows matched the requested --keys.")

    # Compute layer metrics for each unique source directory appearing in the rankings table.
    unique_sources = sorted(set(df["source"].astype(str).tolist()))
    raw_rows = []
    agg_rows = []
    for src in unique_sources:
        src_dir = Path(src)
        raw, agg = compute_layer_metrics_for_source(src_dir, n_values=n_values)
        raw_rows.append(raw)
        agg_rows.append(agg)

    raw_df = pd.concat(raw_rows, ignore_index=True)
    agg_df = pd.concat(agg_rows, ignore_index=True)
    raw_df.to_csv(out_dir / "layer_metrics_raw.csv", index=False, encoding="utf-8-sig")
    agg_df.to_csv(out_dir / "layer_metrics_by_family.csv", index=False, encoding="utf-8-sig")

    df = df.merge(agg_df[["source", "n", "family", "layer_count", "mean_layer_gap"]], on=["source", "n", "family"], how="left")
    if df[["layer_count", "mean_layer_gap"]].isna().any().any():
        raise SystemExit("Missing layer_count/mean_layer_gap after merge; check cg_original_samples.csv availability.")

    rng = np.random.default_rng(int(args.seed))

    # Compute within-block permutation stats and aggregate to stratum level by averaging block rho.
    block_results = []
    for _, block in df.groupby(BLOCK_COLS):
        if len(block) < 4:
            continue
        obs, perm = _block_perm_test(block, n_perm=int(args.n_perm), rng=rng)
        key = {col: block.iloc[0][col] for col in BLOCK_COLS}
        block_results.append(
            {
                **key,
                "n_perm": float(args.n_perm),
                "n_families": float(len(block)),
                "rho_raw": obs.rho_raw,
                "rho_partial": obs.rho_partial,
                "rho_partial_lc": obs.rho_partial_lc,
                "rho_partial_gap": obs.rho_partial_gap,
                "icg_lc_rho": obs.spearman_icg_lc,
                "icg_gap_rho": obs.spearman_icg_gap,
                "perm_raw": perm["raw"],
                "perm_partial": perm["partial"],
                "perm_partial_lc": perm["partial_lc"],
                "perm_partial_gap": perm["partial_gap"],
            }
        )

    if not block_results:
        raise SystemExit("No blocks produced results. Check input table and required variants.")

    # Convert to per-stratum aggregated stats.
    strat_rows = []
    block_df = pd.DataFrame([{k: v for (k, v) in r.items() if not k.startswith("perm_")} for r in block_results])
    for (n, keep, gamma, var), sub in block_df.groupby(STRATUM_COLS):
        # Find the corresponding perm arrays from block_results to build the aggregated permutation distribution.
        idxs = sub.index.to_list()
        perm_raw = np.vstack([block_results[i]["perm_raw"] for i in idxs]).mean(axis=0)
        perm_partial = np.vstack([block_results[i]["perm_partial"] for i in idxs]).mean(axis=0)
        perm_partial_lc = np.vstack([block_results[i]["perm_partial_lc"] for i in idxs]).mean(axis=0)
        perm_partial_gap = np.vstack([block_results[i]["perm_partial_gap"] for i in idxs]).mean(axis=0)

        obs_raw = float(sub["rho_raw"].mean())
        obs_partial = float(sub["rho_partial"].mean())
        obs_partial_lc = float(sub["rho_partial_lc"].mean())
        obs_partial_gap = float(sub["rho_partial_gap"].mean())

        strat_rows.append(
            {
                "n": int(n),
                "keep_ratio": float(keep),
                "gamma": float(gamma),
                "icg_variant": str(var),
                "n_blocks": float(len(sub)),
                "n_perm": float(args.n_perm),
                "obs_mean_spearman_raw": obs_raw,
                "p_abs_mean_spearman_raw": _perm_p_value(obs_raw, perm_raw),
                "obs_mean_spearman_partial_lc_gap": obs_partial,
                "p_abs_mean_spearman_partial_lc_gap": _perm_p_value(obs_partial, perm_partial),
                "obs_mean_spearman_partial_lc": obs_partial_lc,
                "p_abs_mean_spearman_partial_lc": _perm_p_value(obs_partial_lc, perm_partial_lc),
                "obs_mean_spearman_partial_gap": obs_partial_gap,
                "p_abs_mean_spearman_partial_gap": _perm_p_value(obs_partial_gap, perm_partial_gap),
                "mean_spearman_icg_layer_count": float(sub["icg_lc_rho"].mean()),
                "mean_spearman_icg_mean_layer_gap": float(sub["icg_gap_rho"].mean()),
            }
        )

    strat_df = pd.DataFrame(strat_rows).sort_values(["n", "keep_ratio", "gamma", "icg_variant"]).reset_index(drop=True)
    strat_df.to_csv(out_dir / "cg_blockperm_delta_rank_stratified_layercontrolled.csv", index=False, encoding="utf-8-sig")

    report_lines = []
    report_lines.append("# Prediction D: Layer-Controlled Tier-3 Refinement\n")
    report_lines.append(f"- input: `{Path(args.zeta_rank_csv).as_posix()}`")
    report_lines.append(f"- keys: `{', '.join(args.keys)}`")
    report_lines.append(f"- variants: `{','.join(require_variants)}`")
    report_lines.append(f"- n_perm per block: `{int(args.n_perm)}`")
    report_lines.append("")
    report_lines.append("## Stratum Results (Mean Across Blocks)")
    report_lines.append("")
    report_lines.append(strat_df.to_markdown(index=False))
    report_lines.append("")
    (out_dir / "report.md").write_text("\n".join(report_lines), encoding="utf-8")

    print((out_dir / "cg_blockperm_delta_rank_stratified_layercontrolled.csv").as_posix())
    print((out_dir / "report.md").as_posix())


if __name__ == "__main__":
    main()
