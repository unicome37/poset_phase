"""
Prediction A — F8 Series Dimension Selection Reanalysis
========================================================
Reanalyze cached prediction_a_f7_dimension_full.csv with F8 variants
to test whether any F8 model fixes the low-dimension rejection failure.

F8 variants:
  F8a: logH/(N·logN) normalization — makes entropy O(1)
  F8b: N-linear wall β·N·σ((R-Rc)/w) — wall grows with N
  F8c: logH/N density — makes entropy O(logN)
  F8d: logH − μ·N·(1−R) correction — penalizes disorder surplus
  F8a_v2: F8a + lower wall (two-sided)

For each model, test:
  - 4D vs 2D/3D/5D pairwise win rates at each N
  - Overall Prediction A pass rate
  - Whether B (Lor < KR) still holds (NOT tested here, only Lorentzian families)

Usage:
    python prediction_a_f8_series.py
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


# ══════════════════════════════════════════════════════════════════════════
# F8 model definitions
# ══════════════════════════════════════════════════════════════════════════

def F7_baseline(r: dict) -> float:
    return r["F7"]


def F8a(r: dict, *, lam=10.0, eta=0.6, alpha0=16.0, q=0.5,
        Rc=0.25, w=0.015, N0=20.0) -> float:
    """F8a: logH/(N·logN) normalization."""
    N = r["N"]
    nlogn = N * math.log(N) if N > 1 else 1.0
    log_H_norm = r["log_H"] / nlogn
    alpha_N = alpha0 * (N0 / max(N, 1)) ** q
    wall = alpha_N * sigmoid((r["R"] - Rc) / w)
    return (log_H_norm + 0.0004 * r["pi_geo"]
            - lam * r["sigma_hist"] + eta * r["xi_dim"] + wall)


def F8b(r: dict, *, beta=0.20, lam=10.0, eta=0.6,
        Rc=0.25, w=0.015) -> float:
    """F8b: N-linear wall."""
    N = r["N"]
    wall = beta * N * sigmoid((r["R"] - Rc) / w)
    return (r["log_H"] + 0.0004 * r["pi_geo"]
            - lam * r["sigma_hist"] + eta * r["xi_dim"] + wall)


def F8c(r: dict, *, lam=10.0, eta=0.6, alpha0=16.0, q=0.5,
        Rc=0.25, w=0.015, N0=20.0) -> float:
    """F8c: logH/N entropy density."""
    N = r["N"]
    log_H_density = r["log_H"] / max(N, 1)
    alpha_N = alpha0 * (N0 / max(N, 1)) ** q
    wall = alpha_N * sigmoid((r["R"] - Rc) / w)
    return (log_H_density + 0.0004 * r["pi_geo"]
            - lam * r["sigma_hist"] + eta * r["xi_dim"] + wall)


def F8d(r: dict, *, mu=0.18, lam=10.0, eta=0.6, alpha0=16.0, q=0.5,
        Rc=0.25, w=0.015, N0=20.0) -> float:
    """F8d: logH − μ·N·(1−R) correction."""
    N = r["N"]
    log_H_corr = r["log_H"] - mu * N * (1.0 - r["R"])
    alpha_N = alpha0 * (N0 / max(N, 1)) ** q
    wall = alpha_N * sigmoid((r["R"] - Rc) / w)
    return (log_H_corr + 0.0004 * r["pi_geo"]
            - lam * r["sigma_hist"] + eta * r["xi_dim"] + wall)


def F8a_v2(r: dict, *, eta=0.6, lam=10.0, alpha0=16.0, q=0.5,
           Rc=0.25, w=0.015, Rfloor=0.12, w_f=0.02,
           alpha0_f=16.0, N0=20.0) -> float:
    """F8a_v2: logH/(N·logN) + two-sided wall."""
    N = r["N"]
    nlogn = N * math.log(N) if N > 1 else 1.0
    log_H_norm = r["log_H"] / nlogn
    R = r["R"]
    alpha_N = alpha0 * (N0 / max(N, 1)) ** q
    wall_upper = alpha_N * sigmoid((R - Rc) / w)
    alpha_N_f = alpha0_f * (N0 / max(N, 1)) ** q
    wall_lower = alpha_N_f * sigmoid((Rfloor - R) / w_f)
    return (log_H_norm + 0.0004 * r["pi_geo"]
            - lam * r["sigma_hist"] + eta * r["xi_dim"]
            + wall_upper + wall_lower)


# ══════════════════════════════════════════════════════════════════════════
# Model registry — includes parameter sweeps for key models
# ══════════════════════════════════════════════════════════════════════════

MODELS = {}

# Baseline
MODELS["F7 (baseline)"] = F7_baseline

# F8a: η sweep (key parameter for Pred A)
for eta_val in [0.6, 1.0, 2.0, 3.0, 5.0, 8.0]:
    label = f"F8a (η={eta_val})"
    MODELS[label] = (lambda r, e=eta_val: F8a(r, eta=e))

# F8b: β sweep (N-linear wall strength)
for beta_val in [0.10, 0.15, 0.20, 0.25, 0.30, 0.50]:
    label = f"F8b (β={beta_val})"
    MODELS[label] = (lambda r, b=beta_val: F8b(r, beta=b))

# F8c: η sweep
for eta_val in [0.6, 2.0, 5.0]:
    label = f"F8c (η={eta_val})"
    MODELS[label] = (lambda r, e=eta_val: F8c(r, eta=e))

# F8d: μ sweep
for mu_val in [0.10, 0.15, 0.18, 0.22, 0.30]:
    label = f"F8d (μ={mu_val})"
    MODELS[label] = (lambda r, m=mu_val: F8d(r, mu=m))

# F8a_v2: two-sided wall, η sweep
for eta_val in [0.6, 2.0, 5.0]:
    label = f"F8a_v2 (η={eta_val})"
    MODELS[label] = (lambda r, e=eta_val: F8a_v2(r, eta=e))

# F8b + high η: combined N-linear wall with boosted Ξ_d
for beta_val in [0.15, 0.20, 0.30]:
    for eta_val in [2.0, 5.0, 8.0]:
        label = f"F8b (β={beta_val},η={eta_val})"
        MODELS[label] = (lambda r, b=beta_val, e=eta_val:
                         F8b(r, beta=b, eta=e))


# ══════════════════════════════════════════════════════════════════════════
# Evaluation
# ══════════════════════════════════════════════════════════════════════════

def load_csv(path: str) -> list[dict]:
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
                "F7": float(row["F7"]),
                "wall": float(row["wall"]),
            })
    return rows


def evaluate_model(data: list[dict], model_fn, n_values: list[int]) -> dict:
    """Evaluate a single model for Prediction A dimension selection."""
    by_nf = defaultdict(list)
    for r in data:
        by_nf[(r["N"], r["family"])].append(r)

    comparisons = [("Lor4D", "Lor2D"), ("Lor4D", "Lor3D"), ("Lor4D", "Lor5D")]
    result = {"per_pair": {}, "per_n": {}}

    for left, right in comparisons:
        pair_key = f"{left}_vs_{right}"
        total_wins = 0
        total_pairs = 0
        per_n_wins = {}

        for N in n_values:
            l_rows = by_nf.get((N, left), [])
            r_rows = by_nf.get((N, right), [])
            n_pairs = min(len(l_rows), len(r_rows))

            wins = 0
            for i in range(n_pairs):
                if model_fn(l_rows[i]) < model_fn(r_rows[i]):
                    wins += 1

            total_wins += wins
            total_pairs += n_pairs
            pct = 100 * wins / n_pairs if n_pairs > 0 else 0
            per_n_wins[N] = pct

        result["per_pair"][pair_key] = {
            "total_win_pct": 100 * total_wins / max(total_pairs, 1),
            "per_n": per_n_wins,
        }

    # Overall Prediction A: 4D < all others at all N
    all_n_4d_lowest = True
    for N in n_values:
        means = {}
        for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]:
            vals = by_nf.get((N, fam), [])
            if vals:
                means[fam] = np.mean([model_fn(r) for r in vals])
        if means and min(means, key=means.get) != "Lor4D":
            all_n_4d_lowest = False
        result["per_n"][N] = {
            "winner": min(means, key=means.get) if means else "N/A",
            "rank_4d": sorted(means, key=means.get).index("Lor4D") + 1 if "Lor4D" in means else 0,
        }

    result["all_n_4d_wins"] = all_n_4d_lowest
    return result


def generate_report(data, n_values):
    lines = ["# Prediction A — F8 Series Dimension Selection\n"]
    lines.append(f"**Data**: {len(data)} cached samples from prediction_a_f7_dimension_full.csv\n")
    lines.append(f"**N values**: {n_values}\n")
    lines.append(f"**Models tested**: {len(MODELS)}\n")

    # ── Summary table ──
    lines.append("## Summary Table\n")
    lines.append("| Model | 4D<2D | 4D<3D | 4D<5D | All-N winner? | "
                 "best_pair | worst_pair |")
    lines.append("|-------|-------|-------|-------|---------------|"
                 "-----------|------------|")

    best_model = None
    best_score = -999

    model_results = {}
    for label, fn in MODELS.items():
        try:
            res = evaluate_model(data, fn, n_values)
        except Exception as e:
            lines.append(f"| {label} | ERROR: {e} |")
            continue

        model_results[label] = res

        w2 = res["per_pair"]["Lor4D_vs_Lor2D"]["total_win_pct"]
        w3 = res["per_pair"]["Lor4D_vs_Lor3D"]["total_win_pct"]
        w5 = res["per_pair"]["Lor4D_vs_Lor5D"]["total_win_pct"]
        allw = "✅" if res["all_n_4d_wins"] else "❌"
        best_p = max(w2, w3, w5)
        worst_p = min(w2, w3, w5)

        lines.append(f"| {label} | {w2:.0f}% | {w3:.0f}% | {w5:.0f}% | "
                     f"{allw} | {best_p:.0f}% | {worst_p:.0f}% |")

        # Score: geometric mean of three win rates
        score = (w2 * w3 * w5) ** (1/3) if min(w2, w3, w5) > 0 else min(w2, w3, w5)
        if score > best_score:
            best_score = score
            best_model = label

    lines.append(f"\n**Best model**: {best_model} (geo-mean score = {best_score:.1f}%)\n")

    # ── Detailed per-N for top 5 models ──
    lines.append("## Top Models — Per-N Detail\n")

    # Rank models
    ranked = sorted(model_results.items(),
                    key=lambda x: min(
                        x[1]["per_pair"]["Lor4D_vs_Lor2D"]["total_win_pct"],
                        x[1]["per_pair"]["Lor4D_vs_Lor3D"]["total_win_pct"],
                        x[1]["per_pair"]["Lor4D_vs_Lor5D"]["total_win_pct"]
                    ), reverse=True)

    for label, res in ranked[:8]:
        lines.append(f"### {label}\n")
        lines.append("| N | winner | 4D_rank | 4D<2D | 4D<3D | 4D<5D |")
        lines.append("|---|--------|---------|-------|-------|-------|")

        for N in n_values:
            ninfo = res["per_n"].get(N, {})
            w2 = res["per_pair"]["Lor4D_vs_Lor2D"]["per_n"].get(N, 0)
            w3 = res["per_pair"]["Lor4D_vs_Lor3D"]["per_n"].get(N, 0)
            w5 = res["per_pair"]["Lor4D_vs_Lor5D"]["per_n"].get(N, 0)
            lines.append(f"| {N} | {ninfo.get('winner', '?')} | "
                         f"{ninfo.get('rank_4d', '?')} | "
                         f"{w2:.0f}% | {w3:.0f}% | {w5:.0f}% |")
        lines.append("")

    # ── Physics analysis: what changed? ──
    lines.append("## Physics Analysis: Why F8 Helps (or Not)\n")

    if best_model and best_model in model_results:
        bres = model_results[best_model]
        fres = model_results.get("F7 (baseline)", {})
        if fres:
            lines.append("### F7 → Best F8 Comparison\n")
            for pair_key in ["Lor4D_vs_Lor2D", "Lor4D_vs_Lor3D", "Lor4D_vs_Lor5D"]:
                f7_pct = fres["per_pair"][pair_key]["total_win_pct"]
                f8_pct = bres["per_pair"][pair_key]["total_win_pct"]
                delta = f8_pct - f7_pct
                lines.append(f"- **{pair_key}**: F7={f7_pct:.0f}% → "
                             f"F8={f8_pct:.0f}% (Δ={delta:+.0f}%)")

    # ── Key insight: is low-d rejection fixed? ──
    lines.append("\n### Low-Dimension Rejection Status\n")
    # Find best model for 4D vs 2D
    best_vs2d = max(model_results.items(),
                    key=lambda x: x[1]["per_pair"]["Lor4D_vs_Lor2D"]["total_win_pct"])
    best_vs3d = max(model_results.items(),
                    key=lambda x: x[1]["per_pair"]["Lor4D_vs_Lor3D"]["total_win_pct"])

    lines.append(f"- Best model for 4D<2D: **{best_vs2d[0]}** "
                 f"({best_vs2d[1]['per_pair']['Lor4D_vs_Lor2D']['total_win_pct']:.0f}%)")
    lines.append(f"- Best model for 4D<3D: **{best_vs3d[0]}** "
                 f"({best_vs3d[1]['per_pair']['Lor4D_vs_Lor3D']['total_win_pct']:.0f}%)")

    # Check if ANY model achieves >80% for all three pairs
    full_pass = [(l, r) for l, r in model_results.items()
                 if (r["per_pair"]["Lor4D_vs_Lor2D"]["total_win_pct"] > 80
                     and r["per_pair"]["Lor4D_vs_Lor3D"]["total_win_pct"] > 80
                     and r["per_pair"]["Lor4D_vs_Lor5D"]["total_win_pct"] > 80)]

    if full_pass:
        lines.append(f"\n**{len(full_pass)} model(s) achieve >80% on ALL three pairs:**")
        for l, r in full_pass:
            w2 = r["per_pair"]["Lor4D_vs_Lor2D"]["total_win_pct"]
            w3 = r["per_pair"]["Lor4D_vs_Lor3D"]["total_win_pct"]
            w5 = r["per_pair"]["Lor4D_vs_Lor5D"]["total_win_pct"]
            lines.append(f"- {l}: {w2:.0f}%/{w3:.0f}%/{w5:.0f}%")
    else:
        lines.append("\n**No model achieves >80% on all three pairs.**")

    return "\n".join(lines)


def main():
    csv_path = "outputs_d_recovery/prediction_a_f7_dimension_full.csv"
    print(f"Loading: {csv_path}")
    data = load_csv(csv_path)
    n_values = sorted(set(r["N"] for r in data))
    print(f"Loaded {len(data)} rows, N = {n_values}")
    print(f"Testing {len(MODELS)} models...\n")

    report = generate_report(data, n_values)

    outdir = Path("outputs_d_recovery")
    md_path = outdir / "prediction_a_f8_series.md"
    md_path.write_text(report, encoding="utf-8")
    print(f"Saved: {md_path}")
    print("\n" + report)


if __name__ == "__main__":
    main()
