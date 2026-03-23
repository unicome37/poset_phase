"""F8a with two-sided wall: penalize both R > Rc (too-low-d) and R < Rfloor (too-high-d).

The current sigmoid wall only penalizes R > Rc (structures with too many causal
pairs = too-low dimension). But Lor5D escapes with R ≈ 0.1 << Rc = 0.25,
paying zero wall penalty while having similar entropy density to Lor4D.

Fix: Add a floor penalty for R < Rfloor that penalizes structures with too FEW
causal pairs (too-high dimension / too sparse causality).

Physical meaning: admissible structures must have R in a WINDOW [Rfloor, Rc],
not just below Rc. This corresponds to the two-sided EH admissibility:
structures need enough causal richness (not too sparse) AND not too much (not
too dense).

F8a_v2 = logH/(NlogN) + γΠ_geo - λΣ_hist + ηΞ_d + α(N)·σ_upper + α_f(N)·σ_lower

where σ_upper = σ((R-Rc)/w)   [penalizes R > Rc, original wall]
      σ_lower = σ((Rfloor-R)/w_f) [penalizes R < Rfloor, new floor]
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


def sigmoid(x: float) -> float:
    x = max(-30.0, min(30.0, x))
    return 1.0 / (1.0 + math.exp(-x))


def load_data() -> list[dict]:
    path = Path("outputs_unified_functional/prediction_b_f7_scaling.csv")
    rows = []
    with path.open(encoding="utf-8") as f:
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


def compute_F8a_v2(row: dict, *,
                   eta: float = 0.6, lam: float = 10.0,
                   alpha0: float = 16.0, q: float = 0.5,
                   Rc: float = 0.25, w: float = 0.015,
                   Rfloor: float = 0.12, w_f: float = 0.02,
                   alpha0_f: float = 16.0,
                   N0: float = 20.0) -> float:
    N = row["N"]
    nlogn = N * math.log(N) if N > 1 else 1.0
    log_H_norm = row["log_H"] / nlogn
    R = row["R"]

    alpha_N = alpha0 * (N0 / max(N, 1)) ** q
    wall_upper = alpha_N * sigmoid((R - Rc) / w)

    alpha_N_f = alpha0_f * (N0 / max(N, 1)) ** q
    wall_lower = alpha_N_f * sigmoid((Rfloor - R) / w_f)

    return (log_H_norm
            + 0.0004 * row["pi_geo"]
            - lam * row["sigma_hist"]
            + eta * row["xi_dim"]
            + wall_upper
            + wall_lower)


def evaluate(data: list[dict], model_fn, label: str = "") -> dict:
    by_nf: dict[tuple[int, str], list[dict]] = defaultdict(list)
    for r in data:
        by_nf[(r["N"], r["family"])].append(r)
    n_values = sorted(set(r["N"] for r in data))

    # B: Lor < KR
    b_results = {}
    for lor_fam in ["Lor2D", "Lor3D", "Lor4D"]:
        wins = total = 0
        for N in n_values:
            lr = by_nf.get((N, lor_fam), [])
            kr = by_nf.get((N, "KR_like"), [])
            if not lr or not kr:
                continue
            n_p = min(len(lr), len(kr))
            for i in range(n_p):
                if model_fn(lr[i]) < model_fn(kr[i]):
                    wins += 1
                total += 1
        b_results[lor_fam] = 100 * wins / max(total, 1)

    # A: Lor4D < others
    a_wins = a_total = 0
    for other in ["Lor2D", "Lor3D", "Lor5D"]:
        for N in n_values:
            l4 = by_nf.get((N, "Lor4D"), [])
            oth = by_nf.get((N, other), [])
            if not l4 or not oth:
                continue
            n_p = min(len(l4), len(oth))
            for i in range(n_p):
                if model_fn(l4[i]) < model_fn(oth[i]):
                    a_wins += 1
                a_total += 1
    a_rate = 100 * a_wins / max(a_total, 1)

    # A_4vs5 at N>=52
    a5_wins = a5_total = 0
    for N in n_values:
        if N < 52:
            continue
        l4 = by_nf.get((N, "Lor4D"), [])
        l5 = by_nf.get((N, "Lor5D"), [])
        if not l4 or not l5:
            continue
        n_p = min(len(l4), len(l5))
        for i in range(n_p):
            if model_fn(l4[i]) < model_fn(l5[i]):
                a5_wins += 1
            a5_total += 1
    a5_rate = 100 * a5_wins / max(a5_total, 1)

    # A_4vs5 at ALL N
    a5a_wins = a5a_total = 0
    for N in n_values:
        l4 = by_nf.get((N, "Lor4D"), [])
        l5 = by_nf.get((N, "Lor5D"), [])
        if not l4 or not l5:
            continue
        n_p = min(len(l4), len(l5))
        for i in range(n_p):
            if model_fn(l4[i]) < model_fn(l5[i]):
                a5a_wins += 1
            a5a_total += 1
    a5a_rate = 100 * a5a_wins / max(a5a_total, 1)

    return {
        "label": label,
        "B2": b_results.get("Lor2D", 0),
        "B3": b_results.get("Lor3D", 0),
        "B4": b_results.get("Lor4D", 0),
        "A": a_rate,
        "A5_N52": a5_rate,
        "A5_all": a5a_rate,
    }


def main():
    data = load_data()
    print(f"Loaded {len(data)} rows\n")

    lines = ["# F8a Two-Sided Wall Scan\n"]

    # Baseline: F8a without floor
    def f8a_base(row):
        return compute_F8a_v2(row, alpha0_f=0.0)

    base = evaluate(data, f8a_base, "F8a (no floor)")
    print(f"Baseline: B4D={base['B4']:.0f}% A={base['A']:.1f}% "
          f"A5_all={base['A5_all']:.0f}% A5_N52={base['A5_N52']:.0f}%")

    # Scan Rfloor and alpha0_f
    lines.append("## Rfloor × α₀_f scan (w_f=0.02)\n")
    lines.append("| Rfloor | α₀_f | B_Lor2D | B_Lor3D | B_Lor4D | A | 4vs5_all | 4vs5_N≥52 |")
    lines.append("|--------|------|---------|---------|---------|---|----------|-----------|")

    best_score = -999
    best_params = {}

    for Rfloor in [0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20]:
        for alpha0_f in [2, 4, 8, 12, 16, 20, 30]:
            fn = lambda row, rf=Rfloor, af=alpha0_f: compute_F8a_v2(
                row, Rfloor=rf, alpha0_f=af)
            res = evaluate(data, fn)
            lines.append(
                f"| {Rfloor:.2f} | {alpha0_f} | {res['B2']:.0f}% | {res['B3']:.0f}% | "
                f"{res['B4']:.0f}% | {res['A']:.1f}% | {res['A5_all']:.0f}% | "
                f"{res['A5_N52']:.0f}% |")

            # Score: B4D must be 100%, then maximize A + A5
            if res["B4"] >= 99 and res["B2"] >= 95 and res["B3"] >= 95:
                score = res["A"] + res["A5_all"] * 2 + res["A5_N52"] * 3
                if score > best_score:
                    best_score = score
                    best_params = {"Rfloor": Rfloor, "alpha0_f": alpha0_f}

    lines.append(f"\n**Best: Rfloor={best_params.get('Rfloor', '?')}, "
                 f"α₀_f={best_params.get('alpha0_f', '?')}**\n")
    print(f"\nBest: {best_params}")

    # Detailed per-N results for best params
    if best_params:
        lines.append("## Detailed per-N for best parameters\n")
        Rf = best_params["Rfloor"]
        Af = best_params["alpha0_f"]
        fn = lambda row: compute_F8a_v2(row, Rfloor=Rf, alpha0_f=Af)

        by_nf: dict[tuple[int, str], list[dict]] = defaultdict(list)
        for r in data:
            by_nf[(r["N"], r["family"])].append(r)
        n_values = sorted(set(r["N"] for r in data))

        lines.append("| N | Lor2D | Lor3D | Lor4D | Lor5D | KR | 4D rank | 4D<5D? |")
        lines.append("|---|-------|-------|-------|-------|-----|---------|--------|")

        for N in n_values:
            mfs = {}
            for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]:
                grp = by_nf.get((N, fam), [])
                if grp:
                    mfs[fam] = np.mean([fn(r) for r in grp])

            lor_fs = {k: v for k, v in mfs.items() if k.startswith("Lor")}
            sorted_lor = sorted(lor_fs.items(), key=lambda x: x[1])
            rank4d = [i+1 for i, (k, _) in enumerate(sorted_lor) if k == "Lor4D"]
            rank4d = rank4d[0] if rank4d else "?"
            ok_5 = "✅" if mfs.get("Lor4D", 0) < mfs.get("Lor5D", 0) else "❌"

            lines.append(
                f"| {N} | {mfs.get('Lor2D',0):.2f} | {mfs.get('Lor3D',0):.2f} | "
                f"{mfs.get('Lor4D',0):.2f} | {mfs.get('Lor5D',0):.2f} | "
                f"{mfs.get('KR_like',0):.2f} | {rank4d}/4 | {ok_5} |")

    # Also test w_f scan
    lines.append("\n## w_f scan (Rfloor, α₀_f from best)\n")
    if best_params:
        Rf = best_params["Rfloor"]
        Af = best_params["alpha0_f"]
        lines.append(f"Fixed: Rfloor={Rf}, α₀_f={Af}\n")
        lines.append("| w_f | B4D | A | 4vs5_all | 4vs5_N≥52 |")
        lines.append("|-----|-----|---|----------|-----------|")
        for wf in [0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.08]:
            fn = lambda row, w=wf: compute_F8a_v2(row, Rfloor=Rf, alpha0_f=Af, w_f=w)
            res = evaluate(data, fn)
            lines.append(f"| {wf:.3f} | {res['B4']:.0f}% | {res['A']:.1f}% | "
                         f"{res['A5_all']:.0f}% | {res['A5_N52']:.0f}% |")

    report = "\n".join(lines) + "\n"
    out_path = Path("outputs_unified_functional/f8a_twosided_wall.md")
    out_path.write_text(report, encoding="utf-8")
    print(f"\nReport: {out_path}")


if __name__ == "__main__":
    main()
