from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


OUT_DIR = Path("outputs_exploratory/c_threshold_summary")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def build_summary() -> pd.DataFrame:
    rows = [
        {
            "n": 44,
            "point_estimate": 1.41,
            "estimate_label": "exact/mainline",
            "seed_min": 1.38,
            "seed_max": 1.48,
            "seed_note": "one additional offset did not cross by current upper bound 1.50",
            "fully_locked": "no",
        },
        {
            "n": 48,
            "point_estimate": 1.69,
            "estimate_label": "exact",
            "seed_min": 1.66,
            "seed_max": 1.70,
            "seed_note": "three tested offsets all crossed",
            "fully_locked": "yes",
        },
        {
            "n": 52,
            "point_estimate": 1.900,
            "estimate_label": "mixed near-wall (sp4 seed mean)",
            "seed_min": 1.720,
            "seed_max": 2.200,
            "seed_note": "sp4 seeds crossed at 1.72, 1.78, and 2.20 after tail scan; ultrafine baseline was 1.754",
            "fully_locked": "yes",
        },
    ]
    return pd.DataFrame(rows)


def plot_summary(df: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.5))

    x = df["n"].to_numpy(dtype=float)
    y = df["point_estimate"].to_numpy(dtype=float)
    ax.plot(x, y, "-o", color="#2563eb", linewidth=1.4, markersize=6, label=r"point estimate $c_N$")

    for _, row in df.iterrows():
        n = float(row["n"])
        seed_min = float(row["seed_min"])
        seed_max = float(row["seed_max"])
        locked = str(row["fully_locked"]).lower() == "yes"
        if locked:
            ax.vlines(n, seed_min, seed_max, color="#16a34a", linewidth=3.0, alpha=0.85)
        else:
            ax.vlines(n, seed_min, seed_max, color="#d97706", linewidth=3.0, alpha=0.85)
            ax.scatter([n], [seed_max], marker="_", s=260, color="#d97706", linewidths=2.0)
        ax.text(n, row["point_estimate"] + 0.045, f"{row['point_estimate']:.2f}", ha="center", va="bottom", fontsize=8)

    ax.set_xlabel(r"Poset size $N$")
    ax.set_ylabel(r"Restoration threshold $c_N$")
    ax.set_title(r"Current $c_N$ sequence with seed-sensitivity range")
    ax.grid(True, alpha=0.25, linewidth=0.5)
    ax.set_xticks(df["n"].tolist())
    ax.set_ylim(1.25, 2.30)

    legend_handles = [
        plt.Line2D([0], [0], color="#2563eb", marker="o", linewidth=1.4, label="point estimate"),
        plt.Line2D([0], [0], color="#16a34a", linewidth=3.0, label="seed range fully resolved"),
        plt.Line2D([0], [0], color="#d97706", linewidth=3.0, label="seed range truncated by current upper bound"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", fontsize=8, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    df = build_summary()
    csv_path = OUT_DIR / "c_threshold_summary.csv"
    png_path = OUT_DIR / "c_threshold_summary.png"
    md_path = OUT_DIR / "c_threshold_summary.md"

    df.to_csv(csv_path, index=False, encoding="utf-8-sig")
    plot_summary(df, png_path)

    lines = [
        "# c_N Summary",
        "",
        "| N | point estimate | mode | seed range | note |",
        "|---|---:|---|---|---|",
    ]
    for _, row in df.iterrows():
        lines.append(
            f"| {int(row['n'])} | {row['point_estimate']:.3f} | {row['estimate_label']} | "
            f"{row['seed_min']:.3f}–{row['seed_max']:.3f} | {row['seed_note']} |"
        )
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(csv_path.as_posix())
    print(png_path.as_posix())
    print(md_path.as_posix())


if __name__ == "__main__":
    main()
