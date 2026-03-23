"""Render a Markdown snippet from layer-1 consistency scan CSV.

Usage:
  python conjecture_e_layer1_render_summary.py \
    --in outputs_unified_functional/conjecture_e_layer1_consistency_actionSign-1.csv \
    --out outputs_unified_functional/conjecture_e_layer1_summary.md
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description="Render conjecture E layer-1 scan summary as Markdown")
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--anchors", nargs="*", type=float, default=[10.0, 20.0])
    args = ap.parse_args()

    in_path = Path(args.inp)
    out_path = Path(args.out)
    if not in_path.exists():
        raise FileNotFoundError(str(in_path))

    rows = list(csv.DictReader(in_path.open("r", encoding="utf-8")))
    if not rows:
        raise RuntimeError("Empty input CSV")

    anchors = [float(x) for x in args.anchors]
    anchors_int = [int(x) if float(int(x)) == x else None for x in anchors]

    by = defaultdict(dict)
    for r in rows:
        key = (r["action_field"], int(r["N"]))
        lam = float(r["lambda"])
        if "lambda_cross" not in by[key]:
            by[key]["lambda_cross"] = float(r["lambda_cross_mean_focus_vs_kr"])
            by[key]["win_f5"] = float(r["win_rate_f5_focus_vs_kr"])
        if lam in anchors:
            li = int(lam) if float(int(lam)) == lam else None
            by[key][f"win_{li}"] = float(r["win_rate_focus_vs_kr"])
            by[key][f"agree_{li}"] = float(r["agreement_rate_score_vs_f5"])

    acts = sorted({k[0] for k in by})
    lines: list[str] = []
    lines.append(f"Input: `{in_path.as_posix()}`")
    lines.append("")

    for act in acts:
        lines.append(f"### Layer-1 scan: {act}")
        headers = ["N", "lambda_cross (mean)", "F5 win"]
        for a in anchors:
            li = int(a) if float(int(a)) == a else None
            headers += [f"win@lambda={li}", f"agree@lambda={li}"]

        lines.append("| " + " | ".join(headers) + " |")
        lines.append("|" + "|".join(["---:" for _ in headers]) + "|")

        ns = sorted({k[1] for k in by if k[0] == act})
        for n in ns:
            d = by[(act, n)]

            def fmt(v):
                return "" if v is None else f"{v:.2f}"

            row = [str(n), f"{d['lambda_cross']:.2f}", f"{d['win_f5']:.2f}"]
            for a in anchors:
                li = int(a) if float(int(a)) == a else None
                row += [fmt(d.get(f"win_{li}")), fmt(d.get(f"agree_{li}"))]
            lines.append("| " + " | ".join(row) + " |")
        lines.append("")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    print(f"Saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

