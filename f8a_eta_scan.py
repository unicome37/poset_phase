"""F8a η-scan: find Ξ_d weight that restores 4D locking at all N.

The F8a normalization makes all terms O(1), so the original η=0.6 may be
insufficient to prevent Lor5D from having the lowest F at large N.

Scans η and also λ (Σ_hist weight) to find the Pareto-optimal region.
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


def compute_F8a(row: dict, *, eta: float = 0.6, lam: float = 10.0,
                alpha0: float = 16.0, q: float = 0.5,
                Rc: float = 0.25, w: float = 0.015, N0: float = 20.0) -> float:
    N = row["N"]
    nlogn = N * math.log(N) if N > 1 else 1.0
    log_H_norm = row["log_H"] / nlogn
    alpha_N = alpha0 * (N0 / max(N, 1)) ** q
    wall = alpha_N * sigmoid((row["R"] - Rc) / w)
    return log_H_norm + 0.0004 * row["pi_geo"] - lam * row["sigma_hist"] + eta * row["xi_dim"] + wall


def evaluate(data: list[dict], eta: float, lam: float) -> dict:
    """Evaluate all key metrics for given (η, λ)."""
    by_nf: dict[tuple[int, str], list[dict]] = defaultdict(list)
    for r in data:
        by_nf[(r["N"], r["family"])].append(r)
    
    n_values = sorted(set(r["N"] for r in data))
    
    # B: pairwise Lor < KR
    b_results = {}
    for lor_fam in ["Lor2D", "Lor3D", "Lor4D"]:
        wins = 0
        total = 0
        for N in n_values:
            l_rows = by_nf.get((N, lor_fam), [])
            k_rows = by_nf.get((N, "KR_like"), [])
            if not l_rows or not k_rows:
                continue
            n_p = min(len(l_rows), len(k_rows))
            for i in range(n_p):
                fl = compute_F8a(l_rows[i], eta=eta, lam=lam)
                fk = compute_F8a(k_rows[i], eta=eta, lam=lam)
                if fl < fk:
                    wins += 1
                total += 1
        b_results[lor_fam] = 100 * wins / max(total, 1)
    
    # A: Lor4D < Lor2D/3D/5D pairwise
    a_wins = 0
    a_total = 0
    for other in ["Lor2D", "Lor3D", "Lor5D"]:
        for N in n_values:
            l4_rows = by_nf.get((N, "Lor4D"), [])
            oth_rows = by_nf.get((N, other), [])
            if not l4_rows or not oth_rows:
                continue
            n_p = min(len(l4_rows), len(oth_rows))
            for i in range(n_p):
                f4 = compute_F8a(l4_rows[i], eta=eta, lam=lam)
                fo = compute_F8a(oth_rows[i], eta=eta, lam=lam)
                if f4 < fo:
                    a_wins += 1
                a_total += 1
    a_rate = 100 * a_wins / max(a_total, 1)
    
    # A large-N: Lor4D < Lor5D at N≥52
    a5_wins = 0
    a5_total = 0
    for N in n_values:
        if N < 52:
            continue
        l4_rows = by_nf.get((N, "Lor4D"), [])
        l5_rows = by_nf.get((N, "Lor5D"), [])
        if not l4_rows or not l5_rows:
            continue
        n_p = min(len(l4_rows), len(l5_rows))
        for i in range(n_p):
            f4 = compute_F8a(l4_rows[i], eta=eta, lam=lam)
            f5 = compute_F8a(l5_rows[i], eta=eta, lam=lam)
            if f4 < f5:
                a5_wins += 1
            a5_total += 1
    a5_rate = 100 * a5_wins / max(a5_total, 1)
    
    # C: corr(Σ_hist, F) < 0 within Lor2D (use all N pooled)
    lor2d = [r for r in data if r["family"] == "Lor2D"]
    sh = np.array([r["sigma_hist"] for r in lor2d])
    fs = np.array([compute_F8a(r, eta=eta, lam=lam) for r in lor2d])
    from scipy import stats as sp_stats
    rho_c, _ = sp_stats.spearmanr(sh, fs)
    
    return {
        "eta": eta,
        "lam": lam,
        "B_Lor2D": b_results.get("Lor2D", 0),
        "B_Lor3D": b_results.get("Lor3D", 0),
        "B_Lor4D": b_results.get("Lor4D", 0),
        "A_overall": a_rate,
        "A_4vs5_N52": a5_rate,
        "C_rho": rho_c,
    }


def main():
    data = load_data()
    print(f"Loaded {len(data)} rows\n")
    
    lines = ["# F8a η-scan: Ξ_d Weight Optimization\n"]
    
    # ── 1. η scan with fixed λ=10 ──
    lines.append("## 1. η-scan (λ=10 fixed)\n")
    lines.append("| η | B_Lor2D | B_Lor3D | B_Lor4D | A_overall | A_4vs5_N≥52 | C_ρ |")
    lines.append("|---|---------|---------|---------|-----------|-------------|------|")
    
    best_eta = 0.6
    best_score = -999
    
    for eta in [0.0, 0.3, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]:
        res = evaluate(data, eta=eta, lam=10.0)
        lines.append(f"| {eta:.1f} | {res['B_Lor2D']:.0f}% | {res['B_Lor3D']:.0f}% | "
                     f"{res['B_Lor4D']:.0f}% | {res['A_overall']:.1f}% | "
                     f"{res['A_4vs5_N52']:.0f}% | {res['C_rho']:+.3f} |")
        
        # Score: B_Lor4D (must be 100) + A + A_4vs5
        if res["B_Lor4D"] >= 99:
            score = res["A_overall"] + res["A_4vs5_N52"] * 0.5 + (1 if res["C_rho"] < -0.3 else 0) * 10
            if score > best_score:
                best_score = score
                best_eta = eta
        
        print(f"  η={eta:.1f}: B4D={res['B_Lor4D']:.0f}% A={res['A_overall']:.1f}% "
              f"A(4v5)={res['A_4vs5_N52']:.0f}% C_ρ={res['C_rho']:+.3f}")
    
    lines.append(f"\n**Best η = {best_eta:.1f}** (B_Lor4D=100% constraint)\n")
    
    # ── 2. (η, λ) grid scan ──
    lines.append("## 2. (η, λ) grid scan\n")
    lines.append("| η | λ | B_Lor4D | A_overall | A_4vs5_N≥52 | C_ρ |")
    lines.append("|---|---|---------|-----------|-------------|------|")
    
    best_pair = (0.6, 10.0)
    best_pair_score = -999
    
    for eta in [1.0, 2.0, 3.0, 4.0, 5.0]:
        for lam in [5.0, 8.0, 10.0, 12.0, 15.0]:
            res = evaluate(data, eta=eta, lam=lam)
            lines.append(f"| {eta:.0f} | {lam:.0f} | {res['B_Lor4D']:.0f}% | "
                         f"{res['A_overall']:.1f}% | {res['A_4vs5_N52']:.0f}% | "
                         f"{res['C_rho']:+.3f} |")
            
            if res["B_Lor4D"] >= 99:
                score = res["A_overall"] + res["A_4vs5_N52"] * 0.5 + (1 if res["C_rho"] < -0.3 else 0) * 10
                if score > best_score:
                    best_pair_score = score
                    best_pair = (eta, lam)
    
    lines.append(f"\n**Best (η, λ) = ({best_pair[0]:.0f}, {best_pair[1]:.0f})**\n")
    
    # ── 3. Detailed results for best parameters ──
    lines.append("## 3. Detailed Results for Best Parameters\n")
    best_res = evaluate(data, eta=best_eta, lam=10.0)
    lines.append(f"**η={best_eta}, λ=10:**")
    lines.append(f"- B: Lor2D={best_res['B_Lor2D']:.0f}%, Lor3D={best_res['B_Lor3D']:.0f}%, Lor4D={best_res['B_Lor4D']:.0f}%")
    lines.append(f"- A: overall={best_res['A_overall']:.1f}%, 4vs5 N≥52={best_res['A_4vs5_N52']:.0f}%")
    lines.append(f"- C: ρ={best_res['C_rho']:+.3f}\n")
    
    report = "\n".join(lines) + "\n"
    out_path = Path("outputs_unified_functional/f8a_eta_scan.md")
    out_path.write_text(report, encoding="utf-8")
    print(f"\nReport: {out_path}")


if __name__ == "__main__":
    main()
