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


def perm_pvalue_abs_spearman(x: np.ndarray, y: np.ndarray, n_perm: int, seed: int) -> float:
    rng = np.random.default_rng(int(seed))
    obs = abs(spearman_tieaware(x, y))
    count = 1
    for _ in range(int(n_perm)):
        xp = rng.permutation(x)
        if abs(spearman_tieaware(xp, y)) >= obs:
            count += 1
    return float(count / (int(n_perm) + 1))


def stratified_perm_pvalue_abs_mean_spearman(
    df: pd.DataFrame,
    *,
    strata_cols: list[str],
    x_col: str,
    y_col: str,
    n_perm: int,
    seed: int,
) -> tuple[float, float, int]:
    rng = np.random.default_rng(int(seed))
    blocks = []
    for _k, sub in df.groupby(strata_cols, sort=False):
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
        stats = []
        for x, y in blocks:
            stats.append(spearman_tieaware(rng.permutation(x), y))
        if abs(float(np.mean(stats))) >= abs(obs):
            count += 1
    p = float(count / (int(n_perm) + 1))
    return obs, p, len(blocks)


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Prediction D perturbation experiment: test whether Δpenalty_cg predicts Δscore_local (mechanism-independent target)."
    )
    p.add_argument(
        "--sample-cg-csv",
        default="outputs_exploratory/prediction_d_perturbation/perturbation_sample_cg.csv",
        help="Path to perturbation_sample_cg.csv",
    )
    p.add_argument("--out", default="outputs_exploratory/prediction_d_perturbation", help="Output directory.")
    p.add_argument("--n-perm", type=int, default=5000, help="Permutation count (default: 5000).")
    p.add_argument("--seed", type=int, default=0, help="RNG seed.")
    return p


def main() -> None:
    args = build_arg_parser().parse_args()
    in_csv = Path(args.sample_cg_csv)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(in_csv)
    needed = ["n", "family", "sample_id", "perturb", "mean_penalty_cg", "mean_score_local_orig"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing columns in {in_csv.as_posix()}: {missing}")

    base = df[df["perturb"] == "p00"][["n", "family", "sample_id", "mean_penalty_cg", "mean_score_local_orig"]].copy()
    base = base.rename(columns={"mean_penalty_cg": "pen0", "mean_score_local_orig": "score0"})

    rows_pooled = []
    rows_by_n = []
    rows_strat = []
    rows_strat_by_n = []

    for pert in ["p05", "p10", "p20"]:
        sub = df[df["perturb"] == pert][["n", "family", "sample_id", "mean_penalty_cg", "mean_score_local_orig"]].copy()
        sub = sub.rename(columns={"mean_penalty_cg": "pen1", "mean_score_local_orig": "score1"})
        m = sub.merge(base, on=["n", "family", "sample_id"], how="inner")
        m["delta_penalty_cg"] = m["pen1"].astype(float) - m["pen0"].astype(float)
        m["delta_score_local"] = m["score1"].astype(float) - m["score0"].astype(float)

        x = m["delta_penalty_cg"].to_numpy(float)
        y = m["delta_score_local"].to_numpy(float)
        rho = spearman_tieaware(x, y)
        p = perm_pvalue_abs_spearman(x, y, n_perm=int(args.n_perm), seed=int(args.seed))
        rows_pooled.append({"perturb": pert, "n_obs": int(len(m)), "rho_pooled": rho, "p_pooled": p})

        obs_mean, p_strat, n_blocks = stratified_perm_pvalue_abs_mean_spearman(
            m,
            strata_cols=["n", "family"],
            x_col="delta_penalty_cg",
            y_col="delta_score_local",
            n_perm=int(args.n_perm),
            seed=int(args.seed),
        )
        rows_strat.append(
            {
                "perturb": pert,
                "n_blocks": int(n_blocks),
                "obs_mean_spearman_stratified": obs_mean,
                "p_stratified": p_strat,
            }
        )

        for n_val, mn in m.groupby("n", sort=False):
            x_n = mn["delta_penalty_cg"].to_numpy(float)
            y_n = mn["delta_score_local"].to_numpy(float)
            rho_n = spearman_tieaware(x_n, y_n)
            p_n = perm_pvalue_abs_spearman(x_n, y_n, n_perm=int(args.n_perm), seed=int(args.seed))
            rows_by_n.append(
                {
                    "perturb": pert,
                    "n": int(n_val),
                    "n_obs": int(len(mn)),
                    "rho_pooled": rho_n,
                    "p_pooled": p_n,
                }
            )

            obs_n, p_n_strat, nb_n = stratified_perm_pvalue_abs_mean_spearman(
                mn,
                strata_cols=["family"],
                x_col="delta_penalty_cg",
                y_col="delta_score_local",
                n_perm=int(args.n_perm),
                seed=int(args.seed),
            )
            rows_strat_by_n.append(
                {
                    "perturb": pert,
                    "n": int(n_val),
                    "n_blocks": int(nb_n),
                    "obs_mean_spearman_stratified": obs_n,
                    "p_stratified": p_n_strat,
                }
            )

    pooled_df = pd.DataFrame(rows_pooled)
    by_n_df = pd.DataFrame(rows_by_n).sort_values(["perturb", "n"]).reset_index(drop=True)
    strat_df = pd.DataFrame(rows_strat)
    strat_by_n_df = pd.DataFrame(rows_strat_by_n).sort_values(["perturb", "n"]).reset_index(drop=True)

    pooled_df.to_csv(out_dir / "perturbation_independent_target_pooled.csv", index=False, encoding="utf-8-sig")
    by_n_df.to_csv(out_dir / "perturbation_independent_target_by_n.csv", index=False, encoding="utf-8-sig")
    strat_df.to_csv(out_dir / "perturbation_independent_target_stratified.csv", index=False, encoding="utf-8-sig")
    strat_by_n_df.to_csv(out_dir / "perturbation_independent_target_stratified_by_n.csv", index=False, encoding="utf-8-sig")

    report_lines = []
    report_lines.append("# Prediction D: Independent-Target Quasi-Intervention Check\n")
    report_lines.append(f"- input: `{in_csv.as_posix()}`")
    report_lines.append(f"- target: `Y = Δscore_local` (no CG term)")
    report_lines.append(f"- predictor: `X = Δpenalty_cg`")
    report_lines.append(f"- p-values: permutation test, `n_perm={int(args.n_perm)}` (two-sided, +1 correction)")
    report_lines.append("")
    report_lines.append("## Pooled (All N + All Families)\n")
    report_lines.append(pooled_df.to_markdown(index=False))
    report_lines.append("\n## Stratified (by N × family; mean Spearman across strata)\n")
    report_lines.append(strat_df.to_markdown(index=False))
    report_lines.append("\n## Pooled by N\n")
    report_lines.append(by_n_df.to_markdown(index=False))
    report_lines.append("\n## Stratified by N (within N, stratify by family)\n")
    report_lines.append(strat_by_n_df.to_markdown(index=False))
    (out_dir / "perturbation_independent_target_report.md").write_text("\n".join(report_lines), encoding="utf-8")

    print((out_dir / "perturbation_independent_target_pooled.csv").as_posix())
    print((out_dir / "perturbation_independent_target_stratified.csv").as_posix())
    print((out_dir / "perturbation_independent_target_report.md").as_posix())


if __name__ == "__main__":
    main()

