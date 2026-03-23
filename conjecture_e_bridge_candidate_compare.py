"""Conjecture E - Compare bridge candidates side by side.

This script compares the two leading bridge candidates:

  - bd_ratio
  - bdg_d2_corrected_norm

It reports:
  1. Residual bridge power beyond F5 ~ N + family
  2. Fine alpha scans for F6 = -F5 + alpha * candidate

Outputs:
  - outputs_unified_functional/conjecture_e_bridge_candidate_compare.csv
  - outputs_unified_functional/conjecture_e_bridge_candidate_compare.md
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
    metric: str,
    focus_family: str,
    control_family: str,
    alpha_values: np.ndarray,
    entropy_sign: float,
) -> list[dict]:
    out: list[dict] = []
    by_n: dict[int, dict[str, list[dict]]] = {}
    for row in rows:
        if row["family"] not in {focus_family, control_family}:
            continue
        by_n.setdefault(int(row["N"]), {}).setdefault(row["family"], []).append(row)

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
        x_focus = np.array([float(r[metric]) for r in focus], dtype=float)
        x_control = np.array([float(r[metric]) for r in control], dtype=float)

        for alpha in alpha_values:
            score_focus = entropy_sign * f5_focus + float(alpha) * x_focus
            score_control = entropy_sign * f5_control + float(alpha) * x_control
            out.append(
                {
                    "metric": metric,
                    "scope": str(n),
                    "alpha": float(alpha),
                    "win_rate_focus_vs_control": hit_rate(score_focus, score_control),
                    "mean_gap_control_minus_focus": float(np.mean(score_control - score_focus)),
                }
            )

    # pooled scan
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

    f5_focus_vals: list[float] = []
    f5_control_vals: list[float] = []
    x_focus_vals: list[float] = []
    x_control_vals: list[float] = []
    for key in sorted(by_pair):
        pair = by_pair[key]
        if focus_family not in pair or control_family not in pair:
            continue
        f5_focus_vals.append(float(pair[focus_family]["f5_calibrated"]))
        f5_control_vals.append(float(pair[control_family]["f5_calibrated"]))
        x_focus_vals.append(float(pair[focus_family][metric]))
        x_control_vals.append(float(pair[control_family][metric]))
    if f5_focus_vals and f5_control_vals:
        f5_focus = np.array(f5_focus_vals, dtype=float)
        f5_control = np.array(f5_control_vals, dtype=float)
        x_focus = np.array(x_focus_vals, dtype=float)
        x_control = np.array(x_control_vals, dtype=float)
        for alpha in alpha_values:
            score_focus = entropy_sign * f5_focus + float(alpha) * x_focus
            score_control = entropy_sign * f5_control + float(alpha) * x_control
            out.append(
                {
                    "metric": metric,
                    "scope": "ALL",
                    "alpha": float(alpha),
                    "win_rate_focus_vs_control": hit_rate(score_focus, score_control),
                    "mean_gap_control_minus_focus": float(np.mean(score_control - score_focus)),
                }
            )

    return out


def best_row(rows: list[dict]) -> dict:
    best_win = max(float(r["win_rate_focus_vs_control"]) for r in rows)
    tied = [r for r in rows if float(r["win_rate_focus_vs_control"]) == best_win]
    return min(tied, key=lambda r: (abs(float(r["alpha"])), abs(float(r["mean_gap_control_minus_focus"]))))


def main() -> int:
    ap = argparse.ArgumentParser(description="Compare bridge candidates around F5")
    ap.add_argument("--bd-actions", default="outputs_unified_functional/bd_actions.csv")
    ap.add_argument("--out-csv", default="outputs_unified_functional/conjecture_e_bridge_candidate_compare.csv")
    ap.add_argument("--out-md", default="outputs_unified_functional/conjecture_e_bridge_candidate_compare.md")
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

    alpha_values = np.arange(args.alpha_min, args.alpha_max + 0.5 * args.alpha_step, args.alpha_step, dtype=float)
    metrics = ["bd_ratio", "bdg_d2_corrected_norm"]
    scan_rows: list[dict] = []
    for metric in metrics:
        scan_rows.extend(scan_alpha(rows, metric, args.focus_family, args.control_family, alpha_values, args.entropy_sign))

    y_f5 = np.array([float(r["f5_calibrated"]) for r in rows], dtype=float)
    X_base = design_matrix(rows)
    r2_base, _, pred_base = ols_r2(y_f5, X_base)
    residual = y_f5 - pred_base

    bridge_rows: list[dict] = []
    for metric in metrics:
        x = np.array([float(r[metric]) for r in rows], dtype=float)
        r2_ext, beta_ext, _ = ols_r2(y_f5, design_matrix(rows, extra_key=metric))
        resid_corr = float(np.corrcoef(x, residual)[0, 1]) if x.std() > 0 and residual.std() > 0 else float("nan")
        bridge_rows.append(
            {
                "metric": metric,
                "r2_ext": r2_ext,
                "delta_r2": r2_ext - r2_base,
                "coef": float(beta_ext[-1]),
                "resid_corr": resid_corr,
            }
        )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(scan_rows[0].keys()))
        w.writeheader()
        w.writerows(scan_rows)

    grouped: dict[str, list[dict]] = {}
    for row in scan_rows:
        grouped.setdefault(row["metric"], []).append(row)

    md_lines: list[str] = []
    md_lines.append("# Conjecture E bridge candidate compare")
    md_lines.append("")
    md_lines.append(f"- Input: `{path.as_posix()}`")
    md_lines.append(f"- Output CSV: `{out_csv.as_posix()}`")
    md_lines.append("")
    md_lines.append("## Residual bridge comparison")
    md_lines.append("")
    md_lines.append(f"Baseline model: `F5 ~ N + family` with R² = {r2_base:.4f}")
    md_lines.append("")
    md_lines.append("| metric | R² ext | ΔR² | coef | corr(metric, residual F5) |")
    md_lines.append("|---|---:|---:|---:|---:|")
    for row in bridge_rows:
        md_lines.append(
            f"| {row['metric']} | {row['r2_ext']:.4f} | {row['delta_r2']:+.4f} | {row['coef']:+.4f} | {row['resid_corr']:+.4f} |"
        )

    md_lines.append("")
    md_lines.append("## Alpha scan")
    md_lines.append("")
    md_lines.append(f"Scoring convention: `score = {args.entropy_sign:+.1f}·F5 + alpha·metric`")
    md_lines.append("")
    md_lines.append("| metric | scope | best alpha | best win rate | mean gap (control - focus) |")
    md_lines.append("|---|---|---:|---:|---:|")
    for metric in metrics:
        metric_rows = grouped[metric]
        for scope in ["ALL"] + sorted({r["scope"] for r in metric_rows if r["scope"] != "ALL"}, key=lambda s: int(s)):
            sub = [r for r in metric_rows if r["scope"] == scope]
            best = best_row(sub)
            md_lines.append(
                f"| {metric} | {scope} | {best['alpha']:+.2f} | {best['win_rate_focus_vs_control']:.3f} | {best['mean_gap_control_minus_focus']:+.4f} |"
            )

    md_lines.append("")
    md_lines.append("## Takeaway")
    md_lines.append("")
    md_lines.append(
        "The main question is whether any candidate both preserves a nonzero residual bridge and can improve the Lor4D-vs-KR_like ordering at finite alpha. "
        "If one metric keeps the larger ΔR² while both best alphas stay pinned near zero, it is the better bridge proxy for the current stage."
    )

    out_md = Path(args.out_md)
    out_md.write_text("\n".join(md_lines).rstrip() + "\n", encoding="utf-8")

    print(f"Saved: {out_csv}")
    print(f"Saved: {out_md}")
    for row in bridge_rows:
        print(
            f"{row['metric']}: delta_R2={row['delta_r2']:+.4f}, coef={row['coef']:+.4f}, resid_corr={row['resid_corr']:+.4f}"
        )
    for metric in metrics:
        best = best_row([r for r in grouped[metric] if r["scope"] == "ALL"])
        print(
            f"{metric}: best alpha={best['alpha']:+.2f}, win_rate={best['win_rate_focus_vs_control']:.3f}, "
            f"gap={best['mean_gap_control_minus_focus']:+.4f}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
