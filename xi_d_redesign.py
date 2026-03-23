"""Ξ_d Redesign: fix Lor5D escaping the dimensional barrier.

Current Ξ_d = var(local_d) + (mean_local - global_d)² measures self-consistency
only. Both Lor4D (d≈3.9) and Lor5D (d≈4.3) are self-consistent → Ξ_d ≈ 0.

Candidates:
  Ξ_A: Asymmetric target penalty — (d_eff - 4)² with steeper slope for d > 4
  Ξ_B: One-sided upper barrier — σ((d_eff - d_max) / w_d) penalizes d > d_max
  Ξ_C: Consistency + target hybrid — original Ξ_d + κ·max(0, d_eff - 4)²
  Ξ_D: R-based dimension proxy — use R directly as dimension indicator
        (Lor4D R≈0.29, Lor5D R≈0.10 at N=72 — already captured by floor wall)

Key insight: the two-sided wall already penalizes R<Rfloor (catching Lor5D).
So maybe we DON'T need to redesign Ξ_d at all — we just need to verify the
two-sided wall is sufficient combined with the existing Ξ_d.

But the data shows 4vs5 only 46% at N≥52 with the two-sided wall. So either:
1. The floor wall parameters need further tuning, OR
2. Ξ_d needs a target-aware component

This script tests all candidates systematically.
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


# ---------------------------------------------------------------------------
# Model builders
# ---------------------------------------------------------------------------

def make_model(*, eta: float = 0.6, lam: float = 10.0,
               alpha0: float = 16.0, q: float = 0.5,
               Rc: float = 0.25, w: float = 0.015,
               Rfloor: float = 0.0, alpha0_f: float = 0.0, w_f: float = 0.02,
               xi_mode: str = "original",  # original, asym, upper, hybrid
               d_target: float = 4.0,
               asym_ratio: float = 3.0,  # how much steeper for d > d_target
               kappa_target: float = 0.0,  # weight for target-based term
               d_max: float = 4.2, w_d: float = 0.2,  # for upper barrier
               N0: float = 20.0):
    """Return a model function with given parameters."""

    def model(row: dict) -> float:
        N = row["N"]
        nlogn = N * math.log(N) if N > 1 else 1.0
        log_H_norm = row["log_H"] / nlogn
        R = row["R"]

        alpha_N = alpha0 * (N0 / max(N, 1)) ** q
        wall_upper = alpha_N * sigmoid((R - Rc) / w)

        wall_lower = 0.0
        if alpha0_f > 0:
            alpha_N_f = alpha0_f * (N0 / max(N, 1)) ** q
            wall_lower = alpha_N_f * sigmoid((Rfloor - R) / w_f)

        # Ξ_d calculation based on mode
        d = row["d_eff"]
        xi_orig = row["xi_dim"]

        if xi_mode == "original":
            xi = xi_orig
        elif xi_mode == "asym":
            # Asymmetric target: steeper penalty for d > d_target
            delta = d - d_target
            if delta > 0:
                xi = xi_orig + kappa_target * (asym_ratio * delta) ** 2
            else:
                xi = xi_orig + kappa_target * delta ** 2
        elif xi_mode == "upper":
            # One-sided upper barrier: sigmoid penalizes d > d_max
            xi = xi_orig + kappa_target * sigmoid((d - d_max) / w_d)
        elif xi_mode == "hybrid":
            # Consistency + target hybrid
            xi = xi_orig + kappa_target * max(0, d - d_target) ** 2
        elif xi_mode == "target_only":
            # Pure target-anchored
            xi = kappa_target * (d - d_target) ** 2
        else:
            xi = xi_orig

        return (log_H_norm
                + 0.0004 * row["pi_geo"]
                - lam * row["sigma_hist"]
                + eta * xi
                + wall_upper
                + wall_lower)

    return model


def evaluate(data: list[dict], model_fn) -> dict:
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
            for i in range(min(len(lr), len(kr))):
                if model_fn(lr[i]) < model_fn(kr[i]):
                    wins += 1
                total += 1
        b_results[lor_fam] = 100 * wins / max(total, 1)

    # A: Lor4D < others pairwise
    a_wins = a_total = 0
    for other in ["Lor2D", "Lor3D", "Lor5D"]:
        for N in n_values:
            l4 = by_nf.get((N, "Lor4D"), [])
            oth = by_nf.get((N, other), [])
            if not l4 or not oth:
                continue
            for i in range(min(len(l4), len(oth))):
                if model_fn(l4[i]) < model_fn(oth[i]):
                    a_wins += 1
                a_total += 1
    a_rate = 100 * a_wins / max(a_total, 1)

    # 4vs5 all N
    a5_wins = a5_total = 0
    for N in n_values:
        l4 = by_nf.get((N, "Lor4D"), [])
        l5 = by_nf.get((N, "Lor5D"), [])
        if not l4 or not l5:
            continue
        for i in range(min(len(l4), len(l5))):
            if model_fn(l4[i]) < model_fn(l5[i]):
                a5_wins += 1
            a5_total += 1
    a5_rate = 100 * a5_wins / max(a5_total, 1)

    # 4vs5 N>=52
    a5h_wins = a5h_total = 0
    for N in n_values:
        if N < 52:
            continue
        l4 = by_nf.get((N, "Lor4D"), [])
        l5 = by_nf.get((N, "Lor5D"), [])
        if not l4 or not l5:
            continue
        for i in range(min(len(l4), len(l5))):
            if model_fn(l4[i]) < model_fn(l5[i]):
                a5h_wins += 1
            a5h_total += 1
    a5h_rate = 100 * a5h_wins / max(a5h_total, 1)

    # 4vs2 (should be high — Lor4D should beat Lor2D)
    a2_wins = a2_total = 0
    for N in n_values:
        l4 = by_nf.get((N, "Lor4D"), [])
        l2 = by_nf.get((N, "Lor2D"), [])
        if not l4 or not l2:
            continue
        for i in range(min(len(l4), len(l2))):
            if model_fn(l4[i]) < model_fn(l2[i]):
                a2_wins += 1
            a2_total += 1
    a2_rate = 100 * a2_wins / max(a2_total, 1)

    return {
        "B2": b_results.get("Lor2D", 0),
        "B3": b_results.get("Lor3D", 0),
        "B4": b_results.get("Lor4D", 0),
        "A": a_rate,
        "A5": a5_rate,
        "A5h": a5h_rate,
        "A2": a2_rate,
    }


def main():
    data = load_data()
    print(f"Loaded {len(data)} rows\n")

    # First, print d_eff statistics per family at key N values
    print("=== d_eff statistics ===")
    by_nf: dict[tuple[int, str], list[float]] = defaultdict(list)
    for r in data:
        by_nf[(r["N"], r["family"])].append(r["d_eff"])
    for N in [20, 36, 52, 72]:
        print(f"N={N}:")
        for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]:
            ds = by_nf.get((N, fam), [])
            if ds:
                print(f"  {fam:8s}: d_eff = {np.mean(ds):.3f} ± {np.std(ds):.3f}")

    lines = ["# Ξ_d Redesign: Dimensional Barrier Fix\n"]

    # ── Baseline: F8a with two-sided wall, original Ξ_d ──
    lines.append("## 1. Baseline (F8a + two-sided wall, original Ξ_d)\n")
    base = make_model(Rfloor=0.18, alpha0_f=12)
    res = evaluate(data, base)
    lines.append(f"B: L2={res['B2']:.0f}% L3={res['B3']:.0f}% L4={res['B4']:.0f}%")
    lines.append(f"A={res['A']:.1f}% A(4v5)={res['A5']:.0f}% A(4v5 N≥52)={res['A5h']:.0f}%")
    lines.append(f"A(4v2)={res['A2']:.0f}%\n")
    print(f"\nBaseline: B4={res['B4']:.0f}% A={res['A']:.1f}% 4v5={res['A5']:.0f}%/{res['A5h']:.0f}%")

    # ── Candidate A: Asymmetric target ──
    lines.append("## 2. Candidate A: Asymmetric target penalty\n")
    lines.append("| κ | asym_ratio | B4 | A | 4v5_all | 4v5_N≥52 | 4v2 |")
    lines.append("|---|-----------|-----|---|---------|----------|-----|")

    for kappa in [0.5, 1.0, 2.0, 3.0, 5.0, 8.0]:
        for ar in [2.0, 3.0, 5.0]:
            m = make_model(Rfloor=0.18, alpha0_f=12,
                          xi_mode="asym", kappa_target=kappa, asym_ratio=ar)
            res = evaluate(data, m)
            lines.append(f"| {kappa:.1f} | {ar:.0f} | {res['B4']:.0f}% | {res['A']:.1f}% | "
                        f"{res['A5']:.0f}% | {res['A5h']:.0f}% | {res['A2']:.0f}% |")
            if kappa in (2.0, 5.0) and ar == 3.0:
                print(f"  Asym κ={kappa} ar={ar}: B4={res['B4']:.0f}% A={res['A']:.1f}% "
                      f"4v5={res['A5']:.0f}%/{res['A5h']:.0f}%")

    # ── Candidate B: One-sided upper barrier ──
    lines.append("\n## 3. Candidate B: One-sided upper barrier σ(d-d_max)\n")
    lines.append("| κ | d_max | w_d | B4 | A | 4v5_all | 4v5_N≥52 | 4v2 |")
    lines.append("|---|-------|-----|-----|---|---------|----------|-----|")

    for kappa in [1.0, 2.0, 4.0, 8.0, 12.0]:
        for dmax in [4.0, 4.1, 4.2, 4.3]:
            for wd in [0.1, 0.2]:
                m = make_model(Rfloor=0.18, alpha0_f=12,
                              xi_mode="upper", kappa_target=kappa, d_max=dmax, w_d=wd)
                res = evaluate(data, m)
                lines.append(f"| {kappa:.0f} | {dmax:.1f} | {wd:.1f} | {res['B4']:.0f}% | "
                            f"{res['A']:.1f}% | {res['A5']:.0f}% | {res['A5h']:.0f}% | "
                            f"{res['A2']:.0f}% |")

    # ── Candidate C: Hybrid (consistency + target) ──
    lines.append("\n## 4. Candidate C: Hybrid consistency + max(0,d-4)²\n")
    lines.append("| κ | B4 | A | 4v5_all | 4v5_N≥52 | 4v2 |")
    lines.append("|---|-----|---|---------|----------|-----|")

    for kappa in [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0]:
        m = make_model(Rfloor=0.18, alpha0_f=12,
                      xi_mode="hybrid", kappa_target=kappa)
        res = evaluate(data, m)
        lines.append(f"| {kappa:.1f} | {res['B4']:.0f}% | {res['A']:.1f}% | "
                    f"{res['A5']:.0f}% | {res['A5h']:.0f}% | {res['A2']:.0f}% |")
        if kappa in (2.0, 5.0, 8.0):
            print(f"  Hybrid κ={kappa}: B4={res['B4']:.0f}% A={res['A']:.1f}% "
                  f"4v5={res['A5']:.0f}%/{res['A5h']:.0f}%")

    # ── Candidate D: Pure target-anchored ──
    lines.append("\n## 5. Candidate D: Pure target (d-4)²\n")
    lines.append("| κ | η | B4 | A | 4v5_all | 4v5_N≥52 | 4v2 |")
    lines.append("|---|---|-----|---|---------|----------|-----|")

    for kappa in [0.5, 1.0, 2.0, 4.0, 8.0]:
        for eta in [0.3, 0.6, 1.0]:
            m = make_model(Rfloor=0.18, alpha0_f=12, eta=eta,
                          xi_mode="target_only", kappa_target=kappa)
            res = evaluate(data, m)
            lines.append(f"| {kappa:.1f} | {eta:.1f} | {res['B4']:.0f}% | {res['A']:.1f}% | "
                        f"{res['A5']:.0f}% | {res['A5h']:.0f}% | {res['A2']:.0f}% |")

    # ── Find overall best ──
    lines.append("\n## 6. Best Model Search\n")

    best_score = -999
    best_label = ""
    best_res = {}

    configs = []
    # Scan promising candidates
    for kappa in [1.0, 2.0, 3.0, 5.0, 8.0]:
        for ar in [2.0, 3.0, 5.0]:
            configs.append(("Asym", dict(xi_mode="asym", kappa_target=kappa, asym_ratio=ar)))
    for kappa in [2.0, 4.0, 8.0, 12.0]:
        for dmax in [4.0, 4.1, 4.2]:
            configs.append(("Upper", dict(xi_mode="upper", kappa_target=kappa, d_max=dmax, w_d=0.15)))
    for kappa in [1.0, 2.0, 3.0, 5.0, 8.0, 12.0]:
        configs.append(("Hybrid", dict(xi_mode="hybrid", kappa_target=kappa)))
    for kappa in [1.0, 2.0, 4.0, 8.0]:
        configs.append(("Target", dict(xi_mode="target_only", kappa_target=kappa, eta=0.6)))

    # Also scan floor wall params with each Ξ_d mode
    for rf in [0.15, 0.18, 0.20]:
        for af in [8, 12, 16]:
            for kappa in [2.0, 5.0, 8.0]:
                configs.append(("Hybrid+Wall",
                               dict(xi_mode="hybrid", kappa_target=kappa,
                                    Rfloor=rf, alpha0_f=af)))

    for label, extra in configs:
        params = dict(Rfloor=0.18, alpha0_f=12)
        params.update(extra)
        m = make_model(**params)
        res = evaluate(data, m)

        # Score: B4 must be ≥98, then maximize A + A5 + A5h
        if res["B4"] >= 98 and res["B2"] >= 95 and res["B3"] >= 95:
            score = res["A"] + res["A5"] * 1.5 + res["A5h"] * 2.0
            if score > best_score:
                best_score = score
                best_label = f"{label}({', '.join(f'{k}={v}' for k, v in extra.items())})"
                best_res = res

    if best_label:
        lines.append(f"**Best: {best_label}**")
        lines.append(f"- B: L2={best_res['B2']:.0f}% L3={best_res['B3']:.0f}% L4={best_res['B4']:.0f}%")
        lines.append(f"- A={best_res['A']:.1f}%, 4v5={best_res['A5']:.0f}%/{best_res['A5h']:.0f}%")
        lines.append(f"- 4v2={best_res['A2']:.0f}%\n")
        print(f"\n*** BEST: {best_label}")
        print(f"    B4={best_res['B4']:.0f}% A={best_res['A']:.1f}% "
              f"4v5={best_res['A5']:.0f}%/{best_res['A5h']:.0f}% 4v2={best_res['A2']:.0f}%")

    # Per-N detail for best
    if best_label and best_res:
        lines.append("### Per-N detail for best model\n")
        # Re-extract the best params to build per-N table
        # Use the best hybrid+wall config
        for kappa in [2.0, 5.0, 8.0]:
            for rf in [0.15, 0.18, 0.20]:
                for af in [8, 12, 16]:
                    params = dict(Rfloor=rf, alpha0_f=af,
                                 xi_mode="hybrid", kappa_target=kappa)
                    m = make_model(**params)
                    res = evaluate(data, m)
                    if (res["B4"] >= 98 and res["A5h"] >= 60):
                        lines.append(f"\n**κ={kappa}, Rfloor={rf}, α₀_f={af}:**")
                        lines.append("| N | Lor2D | Lor3D | Lor4D | Lor5D | KR | 4D<5D? |")
                        lines.append("|---|-------|-------|-------|-------|-----|--------|")

                        by_nf2: dict[tuple[int, str], list[dict]] = defaultdict(list)
                        for r in data:
                            by_nf2[(r["N"], r["family"])].append(r)
                        for N in sorted(set(r["N"] for r in data)):
                            mfs = {}
                            for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]:
                                grp = by_nf2.get((N, fam), [])
                                if grp:
                                    mfs[fam] = np.mean([m(r) for r in grp])
                            ok = "✅" if mfs.get("Lor4D", 0) < mfs.get("Lor5D", 0) else "❌"
                            lines.append(f"| {N} | {mfs.get('Lor2D',0):.2f} | "
                                        f"{mfs.get('Lor3D',0):.2f} | {mfs.get('Lor4D',0):.2f} | "
                                        f"{mfs.get('Lor5D',0):.2f} | {mfs.get('KR_like',0):.2f} | {ok} |")

    report = "\n".join(lines) + "\n"
    out_path = Path("outputs_unified_functional/xi_d_redesign.md")
    out_path.write_text(report, encoding="utf-8")
    print(f"\nReport: {out_path}")


if __name__ == "__main__":
    main()
