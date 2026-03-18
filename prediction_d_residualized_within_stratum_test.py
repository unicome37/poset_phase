from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def tieaware_rank(x: np.ndarray) -> np.ndarray:
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
            avg = 0.5 * ((i + 1) + j)
            ranks[order[i:j]] = avg
        i = j
    return ranks


def spearman_tieaware(x: np.ndarray, y: np.ndarray) -> float:
    rx = tieaware_rank(x)
    ry = tieaware_rank(y)
    rx = rx - float(rx.mean())
    ry = ry - float(ry.mean())
    den = float(np.sqrt(np.sum(rx * rx) * np.sum(ry * ry)))
    if den <= 1e-12:
        return 0.0
    return float(np.sum(rx * ry) / den)


def residualize(y: np.ndarray, design: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    x = np.asarray(design, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    x = np.column_stack([np.ones(len(y), dtype=float), x])
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ beta


def perm_pvalue_abs_spearman(x: np.ndarray, y: np.ndarray, n_perm: int, seed: int) -> float:
    rng = np.random.default_rng(int(seed))
    obs = abs(spearman_tieaware(x, y))
    count = 1
    for _ in range(int(n_perm)):
        xp = rng.permutation(x)
        if abs(spearman_tieaware(xp, y)) >= obs:
            count += 1
    return float(count / (int(n_perm) + 1))


def stratified_perm_abs_mean_spearman(
    df: pd.DataFrame,
    *,
    strata_cols: list[str],
    x_col: str,
    y_col: str,
    n_perm: int,
    seed: int,
) -> tuple[float, float, int]:
    rng = np.random.default_rng(int(seed))
    blocks: list[tuple[np.ndarray, np.ndarray]] = []
    for _, sub in df.groupby(strata_cols, sort=False):
        x = sub[x_col].to_numpy(float)
        y = sub[y_col].to_numpy(float)
        if len(x) < 4:
            continue
        blocks.append((x, y))
    if not blocks:
        return float("nan"), float("nan"), 0
    obs = float(np.mean([spearman_tieaware(x, y) for x, y in blocks]))
    count = 1
    for _ in range(int(n_perm)):
        vals = []
        for x, y in blocks:
            vals.append(spearman_tieaware(rng.permutation(x), y))
        if abs(float(np.mean(vals))) >= abs(obs):
            count += 1
    return obs, float(count / (int(n_perm) + 1)), len(blocks)


def prepare_delta_frame(sample_cg_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(sample_cg_csv)
    needed = [
        "n",
        "family",
        "sample_id",
        "perturb",
        "mean_penalty_cg",
        "mean_score_local_orig",
        "I_cg_sample",
    ]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing columns in {sample_cg_csv.as_posix()}: {missing}")

    base = (
        df[df["perturb"] == "p00"][
            ["n", "family", "sample_id", "mean_penalty_cg", "mean_score_local_orig", "I_cg_sample"]
        ]
        .copy()
        .rename(
            columns={
                "mean_penalty_cg": "pen0",
                "mean_score_local_orig": "score0",
                "I_cg_sample": "icg0",
            }
        )
    )
    parts: list[pd.DataFrame] = []
    for pert in ["p05", "p10", "p20"]:
        sub = (
            df[df["perturb"] == pert][
                ["n", "family", "sample_id", "mean_penalty_cg", "mean_score_local_orig", "I_cg_sample"]
            ]
            .copy()
            .rename(
                columns={
                    "mean_penalty_cg": "pen1",
                    "mean_score_local_orig": "score1",
                    "I_cg_sample": "icg1",
                }
            )
        )
        merged = sub.merge(base, on=["n", "family", "sample_id"], how="inner")
        merged["perturb"] = pert
        merged["delta_penalty_cg"] = merged["pen1"] - merged["pen0"]
        merged["delta_score_local"] = merged["score1"] - merged["score0"]
        merged["delta_icg"] = merged["icg1"] - merged["icg0"]
        parts.append(merged)
    return pd.concat(parts, ignore_index=True)


def add_within_stratum_transforms(df: pd.DataFrame) -> pd.DataFrame:
    out_frames: list[pd.DataFrame] = []
    for _, sub in df.groupby(["perturb", "n", "family"], sort=False):
        sub = sub.copy()
        sub["x_demeaned"] = sub["delta_penalty_cg"] - float(sub["delta_penalty_cg"].mean())
        sub["y_demeaned"] = sub["delta_score_local"] - float(sub["delta_score_local"].mean())
        design = sub[["pen0", "score0", "icg0"]].to_numpy(float)
        sub["x_resid"] = residualize(sub["delta_penalty_cg"].to_numpy(float), design)
        sub["y_resid"] = residualize(sub["delta_score_local"].to_numpy(float), design)
        out_frames.append(sub)
    return pd.concat(out_frames, ignore_index=True)


def summarize_for_perturb(sub: pd.DataFrame, n_perm: int, seed: int) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    specs = [
        ("raw", "delta_penalty_cg", "delta_score_local"),
        ("demeaned", "x_demeaned", "y_demeaned"),
        ("residualized", "x_resid", "y_resid"),
    ]
    for label, x_col, y_col in specs:
        x = sub[x_col].to_numpy(float)
        y = sub[y_col].to_numpy(float)
        rho_pool = spearman_tieaware(x, y)
        p_pool = perm_pvalue_abs_spearman(x, y, n_perm=n_perm, seed=seed)
        rho_strat, p_strat, n_blocks = stratified_perm_abs_mean_spearman(
            sub,
            strata_cols=["n", "family"],
            x_col=x_col,
            y_col=y_col,
            n_perm=n_perm,
            seed=seed,
        )
        within = []
        for _, blk in sub.groupby(["n", "family"], sort=False):
            within.append(spearman_tieaware(blk[x_col].to_numpy(float), blk[y_col].to_numpy(float)))
        neg_frac = float(np.mean([1.0 if r < 0 else 0.0 for r in within])) if within else float("nan")
        rows.append(
            {
                "transform": label,
                "n_obs": int(len(sub)),
                "n_blocks": int(n_blocks),
                "rho_pooled": float(rho_pool),
                "p_pooled": float(p_pool),
                "rho_stratified": float(rho_strat),
                "p_stratified": float(p_strat),
                "mean_within_rho": float(np.mean(within)) if within else float("nan"),
                "neg_fraction": neg_frac,
            }
        )
    return rows


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Prediction D: within-(N,family) residualized test on perturbation data using a continuous mechanism-independent target."
    )
    p.add_argument(
        "--sample-cg-csv",
        default="outputs_exploratory/prediction_d_perturbation_n32/perturbation_sample_cg_n32.csv",
        help="Path to perturbation sample-level CG CSV.",
    )
    p.add_argument(
        "--out-dir",
        default="outputs_exploratory/prediction_d_perturbation_residualized",
        help="Output directory.",
    )
    p.add_argument("--n-perm", type=int, default=10000, help="Permutation count.")
    p.add_argument("--seed", type=int, default=0, help="RNG seed.")
    return p


def main() -> None:
    args = build_arg_parser().parse_args()
    in_csv = Path(args.sample_cg_csv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    base = prepare_delta_frame(in_csv)
    full = add_within_stratum_transforms(base)
    full.to_csv(out_dir / "residualized_delta_frame.csv", index=False, encoding="utf-8-sig")

    rows = []
    for pert, sub in full.groupby("perturb", sort=False):
        for row in summarize_for_perturb(sub, n_perm=int(args.n_perm), seed=int(args.seed) + hash(pert) % 10000):
            row["perturb"] = pert
            rows.append(row)
    summary = pd.DataFrame(rows)[
        [
            "perturb",
            "transform",
            "n_obs",
            "n_blocks",
            "rho_pooled",
            "p_pooled",
            "rho_stratified",
            "p_stratified",
            "mean_within_rho",
            "neg_fraction",
        ]
    ]
    summary.to_csv(out_dir / "residualized_within_stratum_summary.csv", index=False, encoding="utf-8-sig")

    report_lines = [
        "# Prediction D: Residualized Within-Stratum Test",
        "",
        f"- input: `{in_csv.as_posix()}`",
        "- target: `delta_score_local` (continuous, mechanism-independent)",
        "- predictor: `delta_penalty_cg`",
        "- transforms:",
        "  - `raw`: original deltas",
        "  - `demeaned`: within-(perturb,n,family) demeaning",
        "  - `residualized`: within-(perturb,n,family) residualization on `pen0`, `score0`, `icg0`",
        "",
        summary.to_markdown(index=False),
        "",
        "Interpretation rule:",
        "- if signal reappears after residualization, that supports a within-stratum continuous channel worth pursuing;",
        "- if it stays null, the current perturbation design still supports only structural covariation.",
    ]
    (out_dir / "residualized_within_stratum_report.md").write_text("\n".join(report_lines), encoding="utf-8")

    print((out_dir / "residualized_within_stratum_summary.csv").as_posix())
    print((out_dir / "residualized_within_stratum_report.md").as_posix())


if __name__ == "__main__":
    main()
