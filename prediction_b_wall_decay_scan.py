"""Prediction B — Wall decay exponent scan for Lor4D fix.

Scans q values in α(N) = α₀·(N₀/N)^|q| to find the crossover where
Lor4D consistently beats KR at large N.

Two approaches tested:
  1. q-scan: vary |q| from 0.0 to 0.5
  2. Dimension-dependent wall: α(N,d) = α₀·(N₀/N)^|q| · f(d_eff)
     where f(d) penalizes high-d structures more

Uses cached features from prediction_b_f7_scaling.csv to avoid recomputation.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np


def sigmoid(x: float) -> float:
    x = max(-30.0, min(30.0, x))
    return 1.0 / (1.0 + math.exp(-x))


def load_data(path: str = "outputs_unified_functional/prediction_b_f7_scaling.csv") -> list[dict]:
    rows = []
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append({
                "family": row["family"],
                "N": int(row["N"]),
                "rep": int(row["rep"]),
                "log_H": float(row["log_H"]),
                "pi_geo": float(row["pi_geo"]),
                "sigma_hist": float(row["sigma_hist"]),
                "xi_dim": float(row["xi_dim"]),
                "d_eff": float(row["d_eff"]),
                "R": float(row["R"]),
            })
    return rows


def recompute_F7(row: dict, *, alpha0: float, q: float, lam: float = 10.0,
                 eta: float = 0.6, Rc: float = 0.25, w: float = 0.015,
                 N0: float = 20.0, dim_boost: float = 0.0) -> float:
    """Recompute F7 with given q and optional dim_boost."""
    N = row["N"]
    alpha_N = alpha0 * (N0 / max(N, 1)) ** abs(q)
    # Optional: dimension-dependent boost to wall
    if dim_boost != 0.0:
        d_eff = row["d_eff"]
        # Low d_eff → high R → wall already high → no extra needed
        # High d_eff → low R → wall near 0 → boost to compensate logH growth
        # f(d) = 1 + dim_boost * max(0, d_eff - 3) — penalize d>3 more
        alpha_N *= (1.0 + dim_boost * max(0.0, d_eff - 3.0))

    wall = alpha_N * sigmoid((row["R"] - Rc) / w)
    F7 = (row["log_H"]
          + 0.0004 * row["pi_geo"]
          - lam * row["sigma_hist"]
          + eta * row["xi_dim"]
          + wall)
    return F7


def pairwise_winrate(data: list[dict], left_fam: str, right_fam: str,
                     **f7_kwargs) -> dict:
    """Compute win rate for left < right at each N."""
    from collections import defaultdict
    by_nf: dict[tuple[int, str], list[dict]] = defaultdict(list)
    for r in data:
        by_nf[(r["N"], r["family"])].append(r)

    n_values = sorted(set(r["N"] for r in data))
    results = {}
    total_wins = 0
    total_pairs = 0

    for N in n_values:
        left_rows = by_nf.get((N, left_fam), [])
        right_rows = by_nf.get((N, right_fam), [])
        if not left_rows or not right_rows:
            continue
        n_pairs = min(len(left_rows), len(right_rows))
        wins = 0
        for i in range(n_pairs):
            f_left = recompute_F7(left_rows[i], **f7_kwargs)
            f_right = recompute_F7(right_rows[i], **f7_kwargs)
            if f_left < f_right:
                wins += 1
        results[N] = (wins, n_pairs)
        total_wins += wins
        total_pairs += n_pairs

    return {"per_N": results, "total": (total_wins, total_pairs)}


def main():
    data = load_data()
    print("=" * 80)
    print("Prediction B — Wall Decay Exponent Scan")
    print("=" * 80)

    pairs = [("Lor2D", "KR_like"), ("Lor3D", "KR_like"), ("Lor4D", "KR_like")]
    n_values = sorted(set(r["N"] for r in data))

    # ── Strategy 1: q-scan ────────────────────────────────────────────────
    print("\n## Strategy 1: q-scan (α₀=16, vary |q|)\n")
    q_values = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]

    header = f"{'q':>5s} | {'Lor2D':>6s} | {'Lor3D':>6s} | {'Lor4D':>6s} | "
    header += " | ".join(f"4D@{N}" for N in n_values)
    print(header)
    print("-" * len(header))

    q_scan_results = []
    for q in q_values:
        row_str = f"{q:5.2f} | "
        pair_results = {}
        for left, right in pairs:
            res = pairwise_winrate(data, left, right, alpha0=16.0, q=q)
            tw, tp = res["total"]
            pct = 100 * tw / tp if tp > 0 else 0
            pair_results[left] = pct
            row_str += f"{pct:5.1f}% | "

        # Per-N detail for Lor4D
        res_4d = pairwise_winrate(data, "Lor4D", "KR_like", alpha0=16.0, q=q)
        for N in n_values:
            if N in res_4d["per_N"]:
                w, t = res_4d["per_N"][N]
                row_str += f" {100*w/t:3.0f}% |"
            else:
                row_str += "  N/A |"

        print(row_str)
        q_scan_results.append({
            "q": q,
            "Lor2D_pct": pair_results.get("Lor2D", 0),
            "Lor3D_pct": pair_results.get("Lor3D", 0),
            "Lor4D_pct": pair_results.get("Lor4D", 0),
        })

    # ── Strategy 2: α₀ scan at fixed q ──────────────────────────────────
    print("\n\n## Strategy 2: α₀ scan (q=0.0 i.e. no decay, vary α₀)\n")
    alpha_values = [4, 8, 12, 16, 20, 24, 32, 48, 64]

    header2 = f"{'α₀':>4s} | {'Lor2D':>6s} | {'Lor3D':>6s} | {'Lor4D':>6s} | "
    header2 += " | ".join(f"4D@{N}" for N in n_values)
    print(header2)
    print("-" * len(header2))

    for a0 in alpha_values:
        row_str = f"{a0:4d} | "
        for left, right in pairs:
            res = pairwise_winrate(data, left, right, alpha0=float(a0), q=0.0)
            tw, tp = res["total"]
            pct = 100 * tw / tp if tp > 0 else 0
            row_str += f"{pct:5.1f}% | "

        res_4d = pairwise_winrate(data, "Lor4D", "KR_like", alpha0=float(a0), q=0.0)
        for N in n_values:
            if N in res_4d["per_N"]:
                w, t = res_4d["per_N"][N]
                row_str += f" {100*w/t:3.0f}% |"
            else:
                row_str += "  N/A |"
        print(row_str)

    # ── Strategy 3: dimension-dependent boost ────────────────────────────
    print("\n\n## Strategy 3: dim-dependent boost (q=0.5, α₀=16, vary dim_boost)\n")
    boost_values = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0]

    header3 = f"{'boost':>6s} | {'Lor2D':>6s} | {'Lor3D':>6s} | {'Lor4D':>6s} | "
    header3 += " | ".join(f"4D@{N}" for N in n_values)
    print(header3)
    print("-" * len(header3))

    for boost in boost_values:
        row_str = f"{boost:6.1f} | "
        for left, right in pairs:
            res = pairwise_winrate(data, left, right, alpha0=16.0, q=0.5, dim_boost=boost)
            tw, tp = res["total"]
            pct = 100 * tw / tp if tp > 0 else 0
            row_str += f"{pct:5.1f}% | "

        res_4d = pairwise_winrate(data, "Lor4D", "KR_like", alpha0=16.0, q=0.5, dim_boost=boost)
        for N in n_values:
            if N in res_4d["per_N"]:
                w, t = res_4d["per_N"][N]
                row_str += f" {100*w/t:3.0f}% |"
            else:
                row_str += "  N/A |"
        print(row_str)

    # ── Strategy 4: combined best-fit search ─────────────────────────────
    print("\n\n## Strategy 4: Grid search (α₀, q) — maximize min(Lor2D%, Lor3D%, Lor4D%)\n")
    best_score = 0
    best_params = None
    grid_results = []

    for a0 in [12, 16, 20, 24, 32]:
        for q in [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]:
            pcts = {}
            for left, right in pairs:
                res = pairwise_winrate(data, left, right, alpha0=float(a0), q=q)
                tw, tp = res["total"]
                pcts[left] = 100 * tw / tp if tp > 0 else 0

            min_pct = min(pcts.values())
            mean_pct = sum(pcts.values()) / len(pcts)
            grid_results.append((a0, q, pcts["Lor2D"], pcts["Lor3D"], pcts["Lor4D"], min_pct, mean_pct))

            if min_pct > best_score or (min_pct == best_score and mean_pct > (best_params[6] if best_params else 0)):
                best_score = min_pct
                best_params = (a0, q, pcts["Lor2D"], pcts["Lor3D"], pcts["Lor4D"], min_pct, mean_pct)

    # Sort by min_pct desc, then mean_pct desc
    grid_results.sort(key=lambda x: (x[5], x[6]), reverse=True)

    print(f"{'α₀':>4s} {'q':>5s} | {'Lor2D':>6s} {'Lor3D':>6s} {'Lor4D':>6s} | {'min':>5s} {'mean':>6s}")
    print("-" * 55)
    for a0, q, l2, l3, l4, mn, avg in grid_results[:15]:
        print(f"{a0:4d} {q:5.2f} | {l2:5.1f}% {l3:5.1f}% {l4:5.1f}% | {mn:4.1f}% {avg:5.1f}%")

    if best_params:
        a0, q, l2, l3, l4, mn, avg = best_params
        print(f"\n★ Best: α₀={a0}, q={q:.2f} → Lor2D={l2:.1f}%, Lor3D={l3:.1f}%, Lor4D={l4:.1f}%, min={mn:.1f}%")

        # Show per-N detail for best
        print(f"\nPer-N detail for best (α₀={a0}, q={q:.2f}):")
        for left, right in pairs:
            res = pairwise_winrate(data, left, right, alpha0=float(a0), q=q)
            line = f"  {left} vs {right}: "
            for N in n_values:
                if N in res["per_N"]:
                    w, t = res["per_N"][N]
                    line += f"N={N}:{100*w/t:.0f}% "
            print(line)


if __name__ == "__main__":
    main()
