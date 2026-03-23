"""Aggregate and visualize layer-2 curvature diagnostics.

Inputs are the CSV outputs produced by:
  - curvature_layer2_scaling_scan.py
  - curvature_layer2_toy_curvature_scan.py
  - curvature_layer2_de_sitter_scan.py

The script renders:
  - a markdown summary table
  - a PNG figure with the main trend lines
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_rows(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def group_mean(rows: list[dict], x_key: str, y_key: str, group_keys: tuple[str, ...]) -> dict:
    grouped: dict[tuple, list[tuple[float, float]]] = defaultdict(list)
    for r in rows:
        key = tuple(r[g] for g in group_keys)
        grouped[key].append((float(r[x_key]), float(r[y_key])))
    out = {}
    for key, values in grouped.items():
        xs = np.array([v[0] for v in values], dtype=float)
        ys = np.array([v[1] for v in values], dtype=float)
        out[key] = (xs, ys, float(np.mean(ys)), float(np.std(ys)))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Render curvature layer-2 summary report")
    ap.add_argument("--scaling", default="outputs_unified_functional/curvature_layer2_scaling_scan.csv")
    ap.add_argument("--toy", default="outputs_unified_functional/curvature_layer2_toy_curvature_scan.csv")
    ap.add_argument("--de-sitter", default="outputs_unified_functional/curvature_layer2_de_sitter_scan_dense.csv")
    ap.add_argument("--out-md", default="outputs_unified_functional/curvature_layer2_report.md")
    ap.add_argument("--out-fig", default="outputs_unified_functional/curvature_layer2_report.png")
    args = ap.parse_args()

    scaling_path = Path(args.scaling)
    toy_path = Path(args.toy)
    ds_path = Path(args.de_sitter)
    out_md = Path(args.out_md)
    out_fig = Path(args.out_fig)
    out_md.parent.mkdir(parents=True, exist_ok=True)

    scaling_rows = load_rows(scaling_path)
    toy_rows = load_rows(toy_path)
    ds_rows = load_rows(ds_path)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

    # Panel 1: N scaling
    for d in sorted({int(r["d"]) for r in scaling_rows}):
        sub = [r for r in scaling_rows if int(r["d"]) == d]
        ns = np.array([int(r["N"]) for r in sub], dtype=float)
        pref = np.array([float(r["flat_prefactor_fit"]) for r in sub], dtype=float)
        rhat = np.array([float(r["curvature_r_hat"]) for r in sub], dtype=float)
        axes[0].plot(ns, pref, marker="o", label=f"d={d} pref")
        axes[0].plot(ns, rhat / np.nanmax(np.abs(rhat)) if np.nanmax(np.abs(rhat)) > 0 else rhat, marker="x", linestyle="--", label=f"d={d} R_hat (scaled)")
    axes[0].set_title("Scaling scan")
    axes[0].set_xlabel("N")
    axes[0].set_ylabel("prefactor / scaled R_hat")
    axes[0].legend(fontsize=7)

    # Panel 2: toy kappa scan
    for d in sorted({int(r["d"]) for r in toy_rows}):
        sub = [r for r in toy_rows if int(r["d"]) == d]
        kappas = np.array([float(r["kappa"]) for r in sub], dtype=float)
        rhat = np.array([float(r["curvature_r_hat"]) for r in sub], dtype=float)
        axes[1].plot(kappas, rhat, marker="o", label=f"d={d}")
    axes[1].set_title("Toy curvature scan")
    axes[1].set_xlabel("kappa")
    axes[1].set_ylabel("R_hat")
    axes[1].legend(fontsize=8)

    # Panel 3: de Sitter-like H scan
    for d in sorted({int(r["d"]) for r in ds_rows}):
        sub = [r for r in ds_rows if int(r["d"]) == d]
        hs = np.array([float(r["hubble"]) for r in sub], dtype=float)
        rhat = np.array([float(r["curvature_r_hat"]) for r in sub], dtype=float)
        axes[2].plot(hs, rhat, marker="o", label=f"d={d}")
    axes[2].set_title("de Sitter-like scan")
    axes[2].set_xlabel("H")
    axes[2].set_ylabel("R_hat")
    axes[2].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(out_fig, dpi=200)

    md_lines: list[str] = []
    md_lines.append("# Layer-2 curvature report")
    md_lines.append("")
    md_lines.append(f"- Scaling CSV: `{scaling_path.as_posix()}`")
    md_lines.append(f"- Toy CSV: `{toy_path.as_posix()}`")
    md_lines.append(f"- de Sitter-like CSV: `{ds_path.as_posix()}`")
    md_lines.append(f"- Figure: `{out_fig.as_posix()}`")
    md_lines.append("")
    md_lines.append("## Main observations")
    md_lines.append("")

    for title, rows, key in [
        ("Scaling", scaling_rows, "N"),
        ("Toy curvature", toy_rows, "kappa"),
        ("de Sitter-like", ds_rows, "hubble"),
    ]:
        md_lines.append(f"### {title}")
        for d in sorted({int(r["d"]) for r in rows}):
            sub = [r for r in rows if int(r["d"]) == d]
            xs = np.array([float(r[key]) for r in sub], dtype=float)
            rhat = np.array([float(r["curvature_r_hat"]) for r in sub], dtype=float)
            pref = np.array([float(r["flat_prefactor_fit"]) for r in sub], dtype=float)
            if xs.size > 1:
                slope = float(np.polyfit(xs, rhat, 1)[0])
            else:
                slope = float("nan")
            md_lines.append(
                f"- d={d}: `prefactor` mean={np.mean(pref):.4f}, `R_hat` slope vs `{key}`={slope:+.3f}"
            )
        md_lines.append("")

    out_md.write_text("\n".join(md_lines).rstrip() + "\n", encoding="utf-8")
    print(f"Saved: {out_md}")
    print(f"Saved: {out_fig}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

