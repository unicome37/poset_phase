from __future__ import annotations

from pathlib import Path

import pandas as pd

from experiment import FAMILIES
from observables import layer_profile
from observables_geo import width_height_balance_penalty


def summarize_family(family: str, n: int, samples: int = 5) -> dict:
    rows = []
    for sample_id in range(samples):
        poset = FAMILIES[family](n=n, seed=1000 * n + sample_id)
        profile = layer_profile(poset)
        rows.append(
            {
                "height_layers": float(len(profile)),
                "max_layer_size": float(profile.max()),
                "ratio_w_over_h": float(profile.max() / max(len(profile), 1)),
                "width_height_penalty": float(width_height_balance_penalty(poset)),
            }
        )
    df = pd.DataFrame(rows)
    return {
        "family": family,
        "n": n,
        "mean_height_layers": float(df["height_layers"].mean()),
        "mean_max_layer_size": float(df["max_layer_size"].mean()),
        "mean_ratio_w_over_h": float(df["ratio_w_over_h"].mean()),
        "mean_width_height_penalty": float(df["width_height_penalty"].mean()),
        "std_width_height_penalty": float(df["width_height_penalty"].std(ddof=0)),
        "count": int(len(df)),
    }


def main() -> None:
    out_dir = Path(r"D:\Kiro\outputs_exploratory\width_height_degeneration_check")
    out_dir.mkdir(parents=True, exist_ok=True)

    families = [
        "absolute_layered",
        "KR_like",
        "transitive_percolation",
        "multi_layer_random",
        "lorentzian_like_2d",
    ]
    rows = []
    for n in (20, 40):
        for family in families:
            rows.append(summarize_family(family, n=n, samples=5))

    summary = pd.DataFrame(rows)
    summary.to_csv(out_dir / "width_height_degeneration_summary.csv", index=False, encoding="utf-8-sig")

    theoretical = pd.DataFrame(
        [
            {
                "case": "perfect_chain",
                "n": n,
                "height_layers": float(n),
                "max_layer_size": 1.0,
                "ratio_w_over_h": float(1.0 / n),
                "width_height_penalty": float((__import__('math').log(1.0 / n)) ** 2),
            }
            for n in (20, 40)
        ]
        + [
            {
                "case": "perfect_antichain",
                "n": n,
                "height_layers": 1.0,
                "max_layer_size": float(n),
                "ratio_w_over_h": float(n),
                "width_height_penalty": float((__import__('math').log(n)) ** 2),
            }
            for n in (20, 40)
        ]
        + [
            {
                "case": "balanced_sqrt_like",
                "n": n,
                "height_layers": float(n ** 0.5),
                "max_layer_size": float(n ** 0.5),
                "ratio_w_over_h": 1.0,
                "width_height_penalty": 0.0,
            }
            for n in (20, 40)
        ]
    )
    theoretical.to_csv(out_dir / "width_height_theoretical_extremes.csv", index=False, encoding="utf-8-sig")

    print((out_dir / "width_height_degeneration_summary.csv").as_posix())
    print((out_dir / "width_height_theoretical_extremes.csv").as_posix())
    print()
    print(summary.to_string(index=False))
    print()
    print(theoretical.to_string(index=False))


if __name__ == "__main__":
    main()
