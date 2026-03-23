"""Conjecture E - Fine calibration around bd_ratio.

This script focuses on the strongest residual bridge candidate from the
first bridge fit:

  F6 = F5_calibrated + alpha * bd_ratio

It answers two related questions:

  1. Does a fine scan in alpha improve Lorentzian-vs-KR ordering?
  2. What is the fitted residual coefficient of bd_ratio after controlling
     for N + family in the current F5 baseline?

Inputs:
  - outputs_unified_functional/bd_actions.csv

Outputs:
  - outputs_unified_functional/conjecture_e_bd_ratio_scan.csv
  - outputs_unified_functional/conjecture_e_bd_ratio_scan.md
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


def load_rows(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def hit_rate(focus_scores: np.ndarray, control_scores: np.ndarray) -> float:
    m = min(focus_scores.size, control_scores.size)
    if m == 0:
        return float("nan")
    return float(np.mean(focus_scores[:m] < control_scores[:m]))


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


def scan_alpha(
    rows: list[dict],
    focus_family: str,
    control_family: str,
    alphas: np.ndarray,
    entropy_sign: float,
) -> list[dict]:
    out: list[dict] = []
    by_n: dict[int, dict[str, list[dict]]] = {}
    for row in rows:
        if row["family"] not in {focus_family, control_family}:
            continue
        n = int(row["N"])
        by_n.setdefault(n, {}).setdefault(row["family"], []).append(row)

    for n in sorted(by_n):
        focus = sorted(by_n[n].get(focus_family, []), key=lambda r: int(r["rep"]))
        control = sorted(by_n[n].get(control_family, []), key=lambda r: int(r["rep"]))
        m = min(len(focus), len(control))
        if m == 0:
            continue
        focus = focus[:m]
        control = control[:m]
        f5_focus = np.array([float(r["f5_calibrated"]) for r in focus], dtype=float)
        f5_control = np.array([float(r["f5_calibrated"]) for r in control], dtype=float)
        bd_focus = np.array([float(r["bd_ratio"]) for r in focus], dtype=float)
        bd_control = np.array([float(r["bd_ratio"]) for r in control], dtype=float)

        for alpha in alphas:
            score_focus = entropy_sign * f5_focus + float(alpha) * bd_focus
            score_control = entropy_sign * f5_control + float(alpha) * bd_control
            out.append(
                {
                    "scope": str(n),
                    "alpha": float(alpha),
                    "focus_family": focus_family,
                    "control_family": control_family,
                    "win_rate_focus_vs_control": hit_rate(score_focus, score_control),
                    "mean_score_focus": float(np.mean(score_focus)),
                    "mean_score_control": float(np.mean(score_control)),
                    "mean_gap_control_minus_focus": float(np.mean(score_control - score_focus)),
                }
            )

    # A pooled view across all matched samples.
    focus_all = sorted(
        [r for r in rows if r["family"] == focus_family],
        key=lambda r: (int(r["N"]), int(r["rep"])),
    )
    control_all = sorted(
        [r for r in rows if r["family"] == control_family],
        key=lambda r: (int(r["N"]), int(r["rep"])),
    )
    by_pair: dict[tuple[int, int], dict[str, dict]] = {}
    for row in focus_all:
        by_pair.setdefault((int(row["N"]), int(row["rep"])), {})[focus_family] = row
    for row in control_all:
        by_pair.setdefault((int(row["N"]), int(row["rep"])), {})[control_family] = row

    pooled_focus_f5: list[float] = []
    pooled_control_f5: list[float] = []
    pooled_focus_bd: list[float] = []
    pooled_control_bd: list[float] = []
    for key in sorted(by_pair):
        pair = by_pair[key]
        if focus_family not in pair or control_family not in pair:
            continue
        pooled_focus_f5.append(float(pair[focus_family]["f5_calibrated"]))
        pooled_control_f5.append(float(pair[control_family]["f5_calibrated"]))
        pooled_focus_bd.append(float(pair[focus_family]["bd_ratio"]))
        pooled_control_bd.append(float(pair[control_family]["bd_ratio"]))
    if pooled_focus_f5 and pooled_control_f5:
        f5_focus = np.array(pooled_focus_f5, dtype=float)
        f5_control = np.array(pooled_control_f5, dtype=float)
        bd_focus = np.array(pooled_focus_bd, dtype=float)
        bd_control = np.array(pooled_control_bd, dtype=float)
        for alpha in alphas:
            score_focus = entropy_sign * f5_focus + float(alpha) * bd_focus
            score_control = entropy_sign * f5_control + float(alpha) * bd_control
            out.append(
                {
                    "scope": "ALL",
                    "alpha": float(alpha),
                    "focus_family": focus_family,
                    "control_family": control_family,
                    "win_rate_focus_vs_control": hit_rate(score_focus, score_control),
                    "mean_score_focus": float(np.mean(score_focus)),
                    "mean_score_control": float(np.mean(score_control)),
                    "mean_gap_control_minus_focus": float(np.mean(score_control - score_focus)),
                }
            )

    return out


def best_row(rows: list[dict]) -> dict:
    return max(rows, key=lambda r: (r["win_rate_focus_vs_control"], -abs(r["alpha"])))


def main() -> int:
    ap = argparse.ArgumentParser(description="Fine calibration around bd_ratio bridge term")
    ap.add_argument("--bd-actions", default="outputs_unified_functional/bd_actions.csv")
    ap.add_argument("--out-csv", default="outputs_unified_functional/conjecture_e_bd_ratio_scan.csv")
    ap.add_argument("--out-md", default="outputs_unified_functional/conjecture_e_bd_ratio_scan.md")
    ap.add_argument("--focus-family", default="Lor4D")
    ap.add_argument("--control-family", default="KR_like")
    ap.add_argument("--alpha-min", type=float, default=-10.0)
    ap.add_argument("--alpha-max", type=float, default=10.0)
    ap.add_argument("--alpha-step", type=float, default=0.25)
    ap.add_argument("--entropy-sign", type=float, default=-1.0)
    args = ap.parse_args()

    path = Path(args.bd_actions)
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}. Run conjecture_e_layer1_bdg_report.py first.")

    rows = load_rows(path)
    rows = [r for r in rows if r["family"] in {args.focus_family, args.control_family, "Lor2D", "Lor3D", "Lor5D"}]
    if not rows:
        raise RuntimeError("No rows loaded from bd_actions.csv")

    alphas = np.arange(args.alpha_min, args.alpha_max + 0.5 * args.alpha_step, args.alpha_step, dtype=float)
    scan_rows = scan_alpha(rows, args.focus_family, args.control_family, alphas, args.entropy_sign)
    if not scan_rows:
        raise RuntimeError("No scan rows generated.")

    y_f5 = np.array([float(r["f5_calibrated"]) for r in rows], dtype=float)
    x_bd = np.array([float(r["bd_ratio"]) for r in rows], dtype=float)
    x_bd_mean = float(np.mean(x_bd))
    x_bd_std = float(np.std(x_bd))
    x_bd_z = (x_bd - x_bd_mean) / x_bd_std if x_bd_std > 0 else x_bd * 0.0

    X_base = design_matrix(rows)
    r2_base, _, pred_base = ols_r2(y_f5, X_base)
    r2_ext, beta_ext, _ = ols_r2(y_f5, design_matrix(rows, extra_key="bd_ratio"))
    residual = y_f5 - pred_base
    resid_corr = float(np.corrcoef(x_bd, residual)[0, 1]) if x_bd.std() > 0 and residual.std() > 0 else float("nan")
    resid_corr_z = float(np.corrcoef(x_bd_z, residual)[0, 1]) if x_bd.std() > 0 and residual.std() > 0 else float("nan")

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(scan_rows[0].keys()))
        w.writeheader()
        w.writerows(scan_rows)

    grouped: dict[str, list[dict]] = {}
    for row in scan_rows:
        grouped.setdefault(row["scope"], []).append(row)

    md_lines: list[str] = []
    md_lines.append("# Conjecture E bd_ratio fine scan")
    md_lines.append("")
    md_lines.append(f"- Input: `{path.as_posix()}`")
    md_lines.append(f"- Output CSV: `{out_csv.as_posix()}`")
    md_lines.append("")
    md_lines.append("## Residual bridge")
    md_lines.append("")
    md_lines.append(f"Baseline model: `F5 ~ N + family` with R² = {r2_base:.4f}")
    md_lines.append(f"Extended model: `F5 ~ N + family + bd_ratio` with R² = {r2_ext:.4f}")
    md_lines.append(f"ΔR² = {r2_ext - r2_base:+.4f}")
    md_lines.append(f"Fitted bd_ratio coefficient (raw units) = {beta_ext[-1]:+.4f}")
    md_lines.append(f"Corr(bd_ratio, residual F5) = {resid_corr:+.4f}")
    md_lines.append(f"Corr(z(bd_ratio), residual F5) = {resid_corr_z:+.4f}")
    md_lines.append("")
    md_lines.append("## Alpha scan")
    md_lines.append("")
    md_lines.append(f"Scoring convention: `score = {args.entropy_sign:+.1f}·F5 + alpha·bd_ratio`")
    md_lines.append("")
    md_lines.append("| scope | best alpha | best win rate | mean gap (control - focus) |")
    md_lines.append("|---|---:|---:|---:|")
    for scope in sorted(grouped, key=lambda s: (s != "ALL", int(s) if s.isdigit() else 0)):
        best = best_row(grouped[scope])
        md_lines.append(
            f"| {scope} | {best['alpha']:+.2f} | {best['win_rate_focus_vs_control']:.3f} | {best['mean_gap_control_minus_focus']:+.4f} |"
        )
    md_lines.append("")
    md_lines.append("## Interpretation")
    md_lines.append("")
    md_lines.append(
        "This scan treats `bd_ratio` as the strongest bridge proxy and checks whether a finite `alpha` can actually improve the current Lor4D-vs-KR_like ordering. "
        "If the best alpha stays near zero while the residual R² gain stays positive, that is a good sign that `bd_ratio` is a residual-corrective bridge rather than a rank-flipping replacement for F5."
    )

    out_md = Path(args.out_md)
    out_md.write_text("\n".join(md_lines).rstrip() + "\n", encoding="utf-8")

    print(f"Saved: {out_csv}")
    print(f"Saved: {out_md}")
    for scope in sorted(grouped, key=lambda s: (s != "ALL", int(s) if s.isdigit() else 0)):
        best = best_row(grouped[scope])
        print(
            f"{scope}: best alpha={best['alpha']:+.2f}, win_rate={best['win_rate_focus_vs_control']:.3f}, "
            f"gap={best['mean_gap_control_minus_focus']:+.4f}"
        )
    print(f"Residual baseline R2(F5 ~ N + family) = {r2_base:.4f}")
    print(f"Residual extended  R2(F5 ~ N + family + bd_ratio) = {r2_ext:.4f}")
    print(f"bd_ratio coef (raw units) = {beta_ext[-1]:+.4f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
