"""Prediction B — Lor4D Large-N Fix.

Root cause: ΔlogH(Lor4D - KR) ~ +0.18N (linear growth), while
Δwall shrinks due to α(N) ~ N^{-0.5}. Any bounded wall is overwhelmed.

Three candidate fixes:
  F8a: Normalize logH by expected entropy for given order fraction r
       logH_norm = logH / (N * log(N) * h(r))  where h(r) = 1 - r/r_max
  F8b: N-linear wall: wall = β·N·σ((R-Rc)/w) with β tuned
  F8c: Entropy density: replace logH with logH/N (per-element entropy)

Uses cached features from prediction_b_f7_scaling.csv.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats


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
                "F7": float(row["F7"]),
                "wall": float(row["wall"]),
            })
    return rows


def compute_F8a(row: dict, *, lam: float = 10.0, eta: float = 0.6,
                alpha0: float = 16.0, q: float = -0.5,
                Rc: float = 0.25, w: float = 0.015, N0: float = 20.0) -> float:
    """F8a: Normalize logH by order-fraction-dependent expected entropy.
    
    logH ~ N*logN * c(r), so we normalize:
    logH_norm = logH / max(N*log(N), 1)
    
    This makes the entropy term O(1) instead of O(N*logN), 
    eliminating the linear growth of ΔlogH.
    """
    N = row["N"]
    log_H = row["log_H"]
    nlogn = N * math.log(N) if N > 1 else 1.0
    log_H_norm = log_H / nlogn
    
    alpha_N = alpha0 * (N0 / max(N, 1)) ** abs(q)
    wall = alpha_N * sigmoid((row["R"] - Rc) / w)
    
    # Scale factor to keep F8a in similar range
    return log_H_norm + 0.0004 * row["pi_geo"] - lam * row["sigma_hist"] + eta * row["xi_dim"] + wall


def compute_F8b(row: dict, *, beta: float = 0.2, lam: float = 10.0,
                eta: float = 0.6, Rc: float = 0.25, w: float = 0.015) -> float:
    """F8b: N-linear wall to match logH growth.
    
    wall = β·N·σ((R-Rc)/w) — grows linearly with N.
    """
    N = row["N"]
    wall = beta * N * sigmoid((row["R"] - Rc) / w)
    return row["log_H"] + 0.0004 * row["pi_geo"] - lam * row["sigma_hist"] + eta * row["xi_dim"] + wall


def compute_F8c(row: dict, *, lam: float = 10.0, eta: float = 0.6,
                alpha0: float = 16.0, q: float = -0.5,
                Rc: float = 0.25, w: float = 0.015, N0: float = 20.0) -> float:
    """F8c: Entropy density (logH/N) instead of logH.
    
    Makes the dominant term O(logN) instead of O(N*logN).
    """
    N = row["N"]
    log_H = row["log_H"]
    log_H_density = log_H / max(N, 1)
    
    alpha_N = alpha0 * (N0 / max(N, 1)) ** abs(q)
    wall = alpha_N * sigmoid((row["R"] - Rc) / w)
    
    return log_H_density + 0.0004 * row["pi_geo"] - lam * row["sigma_hist"] + eta * row["xi_dim"] + wall


def compute_F8d(row: dict, *, lam: float = 10.0, eta: float = 0.6,
                alpha0: float = 16.0, q: float = -0.5,
                Rc: float = 0.25, w: float = 0.015, N0: float = 20.0,
                mu: float = 0.18) -> float:
    """F8d: logH with linear R-correction.
    
    F8d = logH - μ·N·(1 - R) + corrections
    
    Since ΔlogH ~ +0.18N comes from Lor4D having lower R,
    subtracting μ·N·(1-R) compensates exactly.
    The physical meaning: penalize structures for their "disorder surplus"
    relative to a fully connected (R=1) baseline.
    """
    N = row["N"]
    log_H = row["log_H"]
    R = row["R"]
    
    # Corrected logH: subtract expected excess entropy from low R
    log_H_corr = log_H - mu * N * (1.0 - R)
    
    alpha_N = alpha0 * (N0 / max(N, 1)) ** abs(q)
    wall = alpha_N * sigmoid((R - Rc) / w)
    
    return log_H_corr + 0.0004 * row["pi_geo"] - lam * row["sigma_hist"] + eta * row["xi_dim"] + wall


MODELS = {
    "F7 (baseline)": lambda row: row["F7"],
    "F8a (logH/NlogN)": compute_F8a,
    "F8b (N-linear wall, β=0.15)": lambda row: compute_F8b(row, beta=0.15),
    "F8b (N-linear wall, β=0.20)": lambda row: compute_F8b(row, beta=0.20),
    "F8b (N-linear wall, β=0.25)": lambda row: compute_F8b(row, beta=0.25),
    "F8c (logH/N density)": compute_F8c,
    "F8d (logH - μN(1-R), μ=0.15)": lambda row: compute_F8d(row, mu=0.15),
    "F8d (logH - μN(1-R), μ=0.18)": lambda row: compute_F8d(row, mu=0.18),
    "F8d (logH - μN(1-R), μ=0.22)": lambda row: compute_F8d(row, mu=0.22),
}


def pairwise_test(data: list[dict], left: str, right: str,
                  model_fn) -> dict:
    """Compute win rate for left < right at each N."""
    by_nf: dict[tuple[int, str], list[dict]] = defaultdict(list)
    for r in data:
        by_nf[(r["N"], r["family"])].append(r)
    
    n_values = sorted(set(r["N"] for r in data))
    results = {}
    total_wins = 0
    total_pairs = 0
    
    for N in n_values:
        left_rows = by_nf.get((N, left), [])
        right_rows = by_nf.get((N, right), [])
        if not left_rows or not right_rows:
            continue
        n_pairs = min(len(left_rows), len(right_rows))
        wins = 0
        deltas = []
        for i in range(n_pairs):
            fl = model_fn(left_rows[i])
            fr = model_fn(right_rows[i])
            deltas.append(fl - fr)
            if fl < fr:
                wins += 1
        results[N] = {"wins": wins, "total": n_pairs, 
                      "mean_delta": float(np.mean(deltas))}
        total_wins += wins
        total_pairs += n_pairs
    
    return {"per_N": results, "total_wins": total_wins, "total_pairs": total_pairs}


def main():
    data = load_data()
    print(f"Loaded {len(data)} rows")
    
    pairs = [
        ("Lor2D", "KR_like"),
        ("Lor3D", "KR_like"),
        ("Lor4D", "KR_like"),
    ]
    
    lines: list[str] = []
    lines.append("# Prediction B — Lor4D Large-N Fix\n")
    
    # ── 1. Diagnosis recap ──
    lines.append("## 1. Root Cause Diagnosis\n")
    lines.append("| N | Lor4D logH | KR logH | ΔlogH | Lor4D wall | KR wall | Δwall | ΔF7 |")
    lines.append("|---|-----------|---------|-------|-----------|---------|-------|-----|")
    
    by_nf: dict[tuple[int, str], list[dict]] = defaultdict(list)
    for r in data:
        by_nf[(r["N"], r["family"])].append(r)
    
    n_values = sorted(set(r["N"] for r in data))
    for N in n_values:
        l4 = by_nf.get((N, "Lor4D"), [])
        kr = by_nf.get((N, "KR_like"), [])
        if not l4 or not kr:
            continue
        ml4h = np.mean([r["log_H"] for r in l4])
        mkrh = np.mean([r["log_H"] for r in kr])
        ml4w = np.mean([r["wall"] for r in l4])
        mkrw = np.mean([r["wall"] for r in kr])
        ml4f = np.mean([r["F7"] for r in l4])
        mkrf = np.mean([r["F7"] for r in kr])
        lines.append(f"| {N} | {ml4h:.1f} | {mkrh:.1f} | {ml4h-mkrh:+.1f} | "
                     f"{ml4w:.1f} | {mkrw:.1f} | {ml4w-mkrw:+.1f} | {ml4f-mkrf:+.1f} |")
    
    lines.append("\n**ΔlogH ~ +0.18N (linear), |Δwall| shrinks (α ~ N^{-0.5}).**\n")
    
    # ── 2. Model comparison ──
    lines.append("## 2. Model Comparison\n")
    lines.append("| Model | Lor2D<KR | Lor3D<KR | Lor4D<KR | Lor4D N≥52 |")
    lines.append("|-------|----------|----------|----------|------------|")
    
    best_model = None
    best_score = -1
    
    for name, fn in MODELS.items():
        results_by_pair = {}
        for left, right in pairs:
            res = pairwise_test(data, left, right, fn)
            pct = 100 * res["total_wins"] / max(res["total_pairs"], 1)
            results_by_pair[(left, right)] = pct
            
            # Large-N Lor4D detail
            if left == "Lor4D":
                large_n_wins = sum(
                    res["per_N"][N]["wins"]
                    for N in n_values if N >= 52 and N in res["per_N"]
                )
                large_n_total = sum(
                    res["per_N"][N]["total"]
                    for N in n_values if N >= 52 and N in res["per_N"]
                )
                large_n_pct = 100 * large_n_wins / max(large_n_total, 1)
        
        l2 = results_by_pair.get(("Lor2D", "KR_like"), 0)
        l3 = results_by_pair.get(("Lor3D", "KR_like"), 0)
        l4 = results_by_pair.get(("Lor4D", "KR_like"), 0)
        
        lines.append(f"| {name} | {l2:.0f}% | {l3:.0f}% | {l4:.0f}% | {large_n_pct:.0f}% |")
        
        # Score: Lor4D win rate (primary) + Lor2D + Lor3D (secondary)
        score = l4 * 3 + l2 + l3 + large_n_pct * 2
        if score > best_score:
            best_score = score
            best_model = name
        
        print(f"  {name}: Lor2D={l2:.0f}% Lor3D={l3:.0f}% Lor4D={l4:.0f}% (N≥52: {large_n_pct:.0f}%)")
    
    lines.append(f"\n**Best model: {best_model}**\n")
    
    # ── 3. Detailed results for best model + F8d variants ──
    lines.append("## 3. Detailed Win Rates per N\n")
    
    detail_models = ["F7 (baseline)", best_model]
    # Add F8d variants
    for name in MODELS:
        if "F8d" in name and name not in detail_models:
            detail_models.append(name)
    
    for name in detail_models:
        fn = MODELS[name]
        lines.append(f"### {name}\n")
        lines.append("| N | Lor2D<KR ΔF | Lor3D<KR ΔF | Lor4D<KR ΔF | Lor4D win% |")
        lines.append("|---|------------|------------|------------|------------|")
        
        for N in n_values:
            entries = []
            for left, right in pairs:
                l_rows = by_nf.get((N, left), [])
                r_rows = by_nf.get((N, right), [])
                if not l_rows or not r_rows:
                    entries.append(("–", "–"))
                    continue
                n_p = min(len(l_rows), len(r_rows))
                deltas = []
                wins = 0
                for i in range(n_p):
                    d = fn(l_rows[i]) - fn(r_rows[i])
                    deltas.append(d)
                    if d < 0:
                        wins += 1
                md = np.mean(deltas)
                wp = 100 * wins / n_p
                entries.append((f"{md:+.1f}", f"{wp:.0f}%"))
            
            lines.append(f"| {N} | {entries[0][0]} | {entries[1][0]} | {entries[2][0]} | {entries[2][1]} |")
    
    # ── 4. β scan for F8b ──
    lines.append("\n## 4. β-scan for F8b (N-linear wall)\n")
    lines.append("| β | Lor2D<KR | Lor3D<KR | Lor4D<KR | Lor4D N≥52 |")
    lines.append("|---|----------|----------|----------|------------|")
    
    for beta in [0.05, 0.10, 0.12, 0.15, 0.18, 0.20, 0.22, 0.25, 0.30, 0.35]:
        fn = lambda row, b=beta: compute_F8b(row, beta=b)
        res_pairs = {}
        for left, right in pairs:
            res = pairwise_test(data, left, right, fn)
            res_pairs[(left, right)] = res
        
        l2 = 100 * res_pairs[("Lor2D", "KR_like")]["total_wins"] / max(res_pairs[("Lor2D", "KR_like")]["total_pairs"], 1)
        l3 = 100 * res_pairs[("Lor3D", "KR_like")]["total_wins"] / max(res_pairs[("Lor3D", "KR_like")]["total_pairs"], 1)
        l4 = 100 * res_pairs[("Lor4D", "KR_like")]["total_wins"] / max(res_pairs[("Lor4D", "KR_like")]["total_pairs"], 1)
        
        res4 = res_pairs[("Lor4D", "KR_like")]
        ln_w = sum(res4["per_N"][N]["wins"] for N in n_values if N >= 52 and N in res4["per_N"])
        ln_t = sum(res4["per_N"][N]["total"] for N in n_values if N >= 52 and N in res4["per_N"])
        ln_p = 100 * ln_w / max(ln_t, 1)
        
        lines.append(f"| {beta:.2f} | {l2:.0f}% | {l3:.0f}% | {l4:.0f}% | {ln_p:.0f}% |")
    
    # ── 5. μ scan for F8d ──
    lines.append("\n## 5. μ-scan for F8d (logH - μN(1-R))\n")
    lines.append("| μ | Lor2D<KR | Lor3D<KR | Lor4D<KR | Lor4D N≥52 |")
    lines.append("|---|----------|----------|----------|------------|")
    
    for mu in [0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20, 0.22, 0.25, 0.30]:
        fn = lambda row, m=mu: compute_F8d(row, mu=m)
        res_pairs = {}
        for left, right in pairs:
            res = pairwise_test(data, left, right, fn)
            res_pairs[(left, right)] = res
        
        l2 = 100 * res_pairs[("Lor2D", "KR_like")]["total_wins"] / max(res_pairs[("Lor2D", "KR_like")]["total_pairs"], 1)
        l3 = 100 * res_pairs[("Lor3D", "KR_like")]["total_wins"] / max(res_pairs[("Lor3D", "KR_like")]["total_pairs"], 1)
        l4 = 100 * res_pairs[("Lor4D", "KR_like")]["total_wins"] / max(res_pairs[("Lor4D", "KR_like")]["total_pairs"], 1)
        
        res4 = res_pairs[("Lor4D", "KR_like")]
        ln_w = sum(res4["per_N"][N]["wins"] for N in n_values if N >= 52 and N in res4["per_N"])
        ln_t = sum(res4["per_N"][N]["total"] for N in n_values if N >= 52 and N in res4["per_N"])
        ln_p = 100 * ln_w / max(ln_t, 1)
        
        lines.append(f"| {mu:.2f} | {l2:.0f}% | {l3:.0f}% | {l4:.0f}% | {ln_p:.0f}% |")
    
    # ── 6. Physical interpretation ──
    lines.append("\n## 6. Physical Interpretation\n")
    lines.append("### Why logH penalizes the right structure\n")
    lines.append("- Lor4D has **low R** (few causal pairs relative to total) → weak ordering → high logH")
    lines.append("- KR has **moderate R** → more ordering → lower logH")
    lines.append("- The gap ΔlogH ~ +0.18N grows linearly, overwhelming any bounded wall\n")
    lines.append("### F8d interpretation: disorder surplus correction\n")
    lines.append("- logH - μN(1-R) subtracts the 'expected excess entropy' from low R")
    lines.append("- Physical meaning: structures shouldn't be penalized for having")
    lines.append("  the *right* causal structure (few pairs = high-d Lorentzian)")
    lines.append("- μ captures the rate at which entropy grows per unit of 'missing' ordering\n")
    lines.append("### F8b interpretation: extensive wall\n")
    lines.append("- β·N·σ((R-Rc)/w) is an N-extensive admissibility term")
    lines.append("- In the continuum limit, this corresponds to an action density")
    lines.append("  (energy per volume) rather than a finite-size correction")
    lines.append("- Natural form: S_EH is extensive (∝ N), so the wall should also be extensive\n")
    
    report = "\n".join(lines) + "\n"
    out_path = Path("outputs_unified_functional/prediction_b_lor4d_fix.md")
    out_path.write_text(report, encoding="utf-8")
    print(f"\nReport: {out_path}")


if __name__ == "__main__":
    main()
