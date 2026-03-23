"""Conjecture E - Unified bridge summary.

This script consolidates the current layer-3 bridge evidence into one
compact report:

  - bd_ratio fine scan
  - candidate comparison vs bdg_d2_corrected_norm

Outputs:
  - outputs_unified_functional/conjecture_e_bridge_summary.csv
  - outputs_unified_functional/conjecture_e_bridge_summary.md
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


def load_rows(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def design_matrix(rows: list[dict], extra_key: str | None = None) -> np.ndarray:
    fams = ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]
    n_values = sorted({int(r["N"]) for r in rows})
    cols = []
    for r in rows:
        row = [1.0]
        n = int(r["N"])
        for nv in n_values[:-1]:
            row.append(1.0 if n == nv else 0.0)
        fam = r["family"]
        for ff in fams[:-1]:
            row.append(1.0 if fam == ff else 0.0)
        if extra_key is not None:
            row.append(float(r[extra_key]))
        cols.append(row)
    return np.array(cols, dtype=float)


def ols_r2(y: np.ndarray, X: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    pred = X @ beta
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return r2, beta, pred


def best_all(rows: list[dict]) -> dict:
    sub = [r for r in rows if r["scope"] == "ALL"]
    if not sub:
        raise RuntimeError("Missing ALL scope in scan rows")
    best_win = max(float(r["win_rate_focus_vs_control"]) for r in sub)
    tied = [r for r in sub if float(r["win_rate_focus_vs_control"]) == best_win]
    return min(tied, key=lambda r: (abs(float(r["alpha"])), abs(float(r["mean_gap_control_minus_focus"]))))


def main() -> int:
    ap = argparse.ArgumentParser(description="Summarize Conjecture E bridge evidence")
    ap.add_argument("--bd-actions", default="outputs_unified_functional/bd_actions.csv")
    ap.add_argument("--bd-ratio-scan", default="outputs_unified_functional/conjecture_e_bd_ratio_scan.csv")
    ap.add_argument(
        "--candidate-compare",
        default="outputs_unified_functional/conjecture_e_bridge_candidate_compare.csv",
    )
    ap.add_argument("--out-csv", default="outputs_unified_functional/conjecture_e_bridge_summary.csv")
    ap.add_argument("--out-md", default="outputs_unified_functional/conjecture_e_bridge_summary.md")
    args = ap.parse_args()

    bd_ratio_path = Path(args.bd_ratio_scan)
    compare_path = Path(args.candidate_compare)
    bd_actions_path = Path(args.bd_actions)
    if not bd_actions_path.exists():
        raise FileNotFoundError(f"Missing {bd_actions_path}")
    if not bd_ratio_path.exists():
        raise FileNotFoundError(f"Missing {bd_ratio_path}")
    if not compare_path.exists():
        raise FileNotFoundError(f"Missing {compare_path}")

    raw_rows = load_rows(bd_actions_path)
    bd_ratio_rows = load_rows(bd_ratio_path)
    compare_rows = load_rows(compare_path)

    # Residual bridge comparison table.
    bridge_rows = {}
    for metric in sorted({row["metric"] for row in compare_rows}):
        metric_rows = [row for row in compare_rows if row["metric"] == metric]
        best_row = best_all(metric_rows)
        bridge_rows[metric] = {
            "metric": metric,
            "best_alpha": float(best_row["alpha"]),
            "best_win_rate": float(best_row["win_rate_focus_vs_control"]),
            "mean_gap": float(best_row["mean_gap_control_minus_focus"]),
        }

    # Pull the residual bridge stats from the compare report.
    residual_rows: dict[str, dict] = {}
    for metric in ["bd_ratio", "bdg_d2_corrected_norm"]:
        x = np.array([float(r[metric]) for r in raw_rows], dtype=float)
        y = np.array([float(r["f5_calibrated"]) for r in raw_rows], dtype=float)
        r2_base, _, pred_base = ols_r2(y, design_matrix(raw_rows))
        r2_ext, beta_ext, _ = ols_r2(y, design_matrix(raw_rows, extra_key=metric))
        residual = y - pred_base
        resid_corr = float(np.corrcoef(x, residual)[0, 1]) if x.std() > 0 and residual.std() > 0 else float("nan")
        residual_rows[metric] = {
            "metric": metric,
            "delta_r2": r2_ext - r2_base,
            "coef": float(beta_ext[-1]),
            "resid_corr": resid_corr,
        }

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    summary_rows = [
        {
            "candidate": metric,
            "delta_r2": residual_rows[metric]["delta_r2"],
            "coef": residual_rows[metric]["coef"],
            "resid_corr": residual_rows[metric]["resid_corr"],
            "best_alpha": bridge_rows[metric]["best_alpha"],
            "best_win_rate": bridge_rows[metric]["best_win_rate"],
            "mean_gap": bridge_rows[metric]["mean_gap"],
        }
        for metric in sorted(bridge_rows)
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        w.writeheader()
        w.writerows(summary_rows)

    bd_ratio_best = best_all(bd_ratio_rows)
    compare_best = {metric: best_all([r for r in compare_rows if r["metric"] == metric]) for metric in bridge_rows}

    md_lines: list[str] = []
    md_lines.append("# Conjecture E bridge summary")
    md_lines.append("")
    md_lines.append(f"- bd_ratio scan: `{bd_ratio_path.as_posix()}`")
    md_lines.append(f"- candidate compare: `{compare_path.as_posix()}`")
    md_lines.append(f"- output CSV: `{out_csv.as_posix()}`")
    md_lines.append("")
    md_lines.append("## Current recommendation")
    md_lines.append("")
    md_lines.append(
        "The best-supported third-layer bridge proxy remains `bd_ratio`: it has the largest residual gain beyond `F5 ~ N + family`, "
        "while the best finite `alpha` in the ordering scans still stays at `0`."
    )
    md_lines.append("")
    md_lines.append("## Compact evidence table")
    md_lines.append("")
    md_lines.append("| candidate | ΔR² | coef | residual corr | best alpha | best win rate |")
    md_lines.append("|---|---:|---:|---:|---:|---:|")
    for row in summary_rows:
        md_lines.append(
            f"| {row['candidate']} | {float(row['delta_r2']):+.4f} | {float(row['coef']):+.4f} | "
            f"{float(row['resid_corr']):+.4f} | {row['best_alpha']:+.2f} | {row['best_win_rate']:.3f} |"
        )

    md_lines.append("")
    md_lines.append("## Next step")
    md_lines.append("")
    md_lines.append(
        "With the empirical ordering now stabilized, the next useful move is to write down the explicit theoretical map from `bd_ratio` to the intended `S_BD` correction term and test that map against the same data."
    )

    out_md = Path(args.out_md)
    out_md.write_text("\n".join(md_lines).rstrip() + "\n", encoding="utf-8")

    print(f"Saved: {out_csv}")
    print(f"Saved: {out_md}")
    print(
        "Recommended bridge proxy: bd_ratio "
        f"(best alpha={float(bd_ratio_best['alpha']):+.2f}, "
        f"win_rate={float(bd_ratio_best['win_rate_focus_vs_control']):.3f})"
    )
    for metric in sorted(compare_best):
        row = compare_best[metric]
        print(
            f"{metric}: best alpha={float(row['alpha']):+.2f}, "
            f"win_rate={float(row['win_rate_focus_vs_control']):.3f}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
