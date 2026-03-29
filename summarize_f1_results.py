from __future__ import annotations

import json
from pathlib import Path
from collections import defaultdict, Counter

def main() -> None:
    base = Path(__file__).resolve().parent
    in_json = base / "outputs_carlip" / "falsify_c1_family_pressure.json"
    out_md = base / "outputs_carlip" / "f1_fullgrid_summary.md"
    out_png = base / "outputs_carlip" / "f1_fullgrid_lor4d_rank.png"

    if not in_json.exists():
        raise FileNotFoundError(f"Missing input: {in_json}")

    data = json.loads(in_json.read_text(encoding="utf-8"))
    records = data.get("per_seed_records", [])
    n_grid = data.get("n_grid", [])
    outrank_counts = data.get("outrank_counts", {})

    # Aggregate Lor4D rank stats by N
    ranks_by_n: dict[int, list[int]] = defaultdict(list)
    winners_by_n: dict[int, Counter] = defaultdict(Counter)

    for r in records:
        n = int(r["N"])
        ranks_by_n[n].append(int(r["lor4d_rank"]))
        winners_by_n[n][str(r["winner"])] += 1

    avg_rank = {n: (sum(v) / len(v) if v else float("nan")) for n, v in ranks_by_n.items()}
    min_rank = {n: (min(v) if v else None) for n, v in ranks_by_n.items()}
    max_rank = {n: (max(v) if v else None) for n, v in ranks_by_n.items()}

    # Top competitor pressure (max outrank count across N)
    pressure = []
    for fam, m in outrank_counts.items():
        max_hit = max(int(v) for v in m.values()) if m else 0
        total_hit = sum(int(v) for v in m.values()) if m else 0
        pressure.append((fam, max_hit, total_hit))
    pressure.sort(key=lambda x: (x[1], x[2], x[0]), reverse=True)

    # Plot Lor4D rank mean/min/max across N (optional)
    xs = sorted(int(n) for n in n_grid)
    y_avg = [float(avg_rank.get(n, float("nan"))) for n in xs]
    y_min: list[float] = []
    y_max: list[float] = []
    for n in xs:
        mn = min_rank.get(n)
        mx = max_rank.get(n)
        y_min.append(float("nan") if mn is None else float(mn))
        y_max.append(float("nan") if mx is None else float(mx))

    plot_ok = False
    try:
        import matplotlib.pyplot as plt  # type: ignore

        plt.figure(figsize=(8, 4.8))
        plt.plot(xs, y_avg, marker="o", label="mean rank")
        plt.plot(xs, y_min, marker="s", linestyle="--", label="min rank")
        plt.plot(xs, y_max, marker="^", linestyle="--", label="max rank")
        plt.ylim(0.8, max(max(y_max), 2) + 0.2)
        plt.xlabel("N")
        plt.ylabel("Lor4D rank")
        plt.title("F1 Full-grid: Lor4D rank vs N")
        plt.grid(alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(out_png, dpi=180)
        plt.close()
        plot_ok = True
    except Exception:
        plot_ok = False

    lines: list[str] = []
    lines.append("# F1 Full-grid Summary")
    lines.append("")
    lines.append(f"- Input: `{in_json.as_posix()}`")
    lines.append(f"- Hard fail: **{'YES' if data.get('hard_fail') else 'NO'}**")
    lines.append(f"- Rule: consecutive_N >= {data['rule']['consecutive_n']}, seed_hits >= {data['rule']['min_seed_success']}/{data['rule']['seed_runs']}")
    lines.append(f"- N grid: {n_grid}")
    lines.append("")
    lines.append("## Lor4D rank stability")
    lines.append("")
    lines.append("| N | mean rank | min rank | max rank | winner census |")
    lines.append("|---:|---:|---:|---:|---|")

    for n in xs:
        wc = winners_by_n.get(n, Counter())
        wc_text = ", ".join([f"{k}:{v}" for k, v in wc.items()]) if wc else "-"
        lines.append(f"| {n} | {avg_rank.get(n, float('nan')):.2f} | {min_rank.get(n, '-') } | {max_rank.get(n, '-') } | {wc_text} |")

    lines.append("")
    lines.append("## Competitor pressure (top 10 by outrank count)")
    lines.append("")
    lines.append("| Family | max hits at one N | total hits across N |")
    lines.append("|---|---:|---:|")
    for fam, mh, th in pressure[:10]:
        lines.append(f"| {fam} | {mh} | {th} |")

    lines.append("")
    lines.append("## Figure")
    lines.append("")
    if plot_ok:
        lines.append(f"- `outputs_carlip/{out_png.name}`")
    else:
        lines.append("- Plot skipped (matplotlib unavailable or interrupted).")

    out_md.write_text("\n".join(lines), encoding="utf-8")
    print(out_md.as_posix())
    if plot_ok:
        print(out_png.as_posix())


if __name__ == "__main__":
    main()
