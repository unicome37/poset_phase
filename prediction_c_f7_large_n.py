"""
Prediction C — F7 Large-N Verification
=======================================
Core question: Does higher Σ_hist predict lower F7 within families at large N?

Prediction C states: history deposition (measured by Σ_hist) compresses
combinatorial entropy. In the unified functional:
    F7 = logH + 0.0004·Π_geo − λ·Σ_hist + η·Ξ_d + wall
the −λ·Σ_hist term means higher Σ_hist → lower F7 (other things equal).

BUT this is tautological if we just look at F7's formula. The real test is:
  Q1. Within a fixed (family, N) cell, does Σ_hist correlate with logH?
      (i.e., does the structural claim hold: more history → less entropy?)
  Q2. Does the within-cell ρ(Σ_hist, logH) strengthen with N?
  Q3. What fraction of logH variance is explained by Σ_hist within-cell?
  Q4. Are the existing C-evidence statistics (layer_count, mean_layer_gap)
      captured by Σ_hist?
  Q5. Does F7 ranking reflect Σ_hist ordering within families?

This script can:
  (a) Reanalyze cached prediction_a_f7_dimension_full.csv (200 samples)
  (b) Generate NEW data with more reps for statistical power

Usage:
    python prediction_c_f7_large_n.py [--fresh --reps 15]
"""
from __future__ import annotations

import argparse
import csv
import math
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

from generators import (
    Poset,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
    generate_kr_like,
)
from unified_functional import (
    compute_log_H,
    compute_pi_geo,
    compute_sigma_hist,
    compute_xi_dim,
)
from prediction_a_bd_bridge import count_intervals_fast


# ══════════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════════

def sigmoid(x: float) -> float:
    x = max(-30.0, min(30.0, x))
    return 1.0 / (1.0 + math.exp(-x))


def compute_R(poset: Poset) -> float:
    counts = count_intervals_fast(poset)
    total = sum(counts.values())
    if total == 0:
        return 0.0
    C0 = counts.get(0, 0)
    return 1.0 - C0 / total


def adaptive_sis_runs(N: int) -> int:
    if N <= 36:
        return 512
    elif N <= 52:
        return 256
    elif N <= 72:
        return 128
    else:
        return 64


def compute_layer_count(poset: Poset) -> int:
    """Layer count = length of longest chain + 1 = number of antichains in a chain decomposition."""
    n = poset.n
    adj = poset.closure
    # BFS topological order to find longest path
    in_deg = [0] * n
    children = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and adj[i][j]:
                # check if direct (not transitive)
                children[i].append(j)
                in_deg[j] += 1

    # longest path via DP
    from collections import deque
    topo = []
    q = deque([i for i in range(n) if in_deg[i] == 0])
    temp_in = in_deg[:]
    while q:
        u = q.popleft()
        topo.append(u)
        for v in children[u]:
            temp_in[v] -= 1
            if temp_in[v] == 0:
                q.append(v)

    dist = [0] * n
    for u in topo:
        for v in children[u]:
            if dist[v] < dist[u] + 1:
                dist[v] = dist[u] + 1
    return max(dist) + 1 if dist else 1


def compute_all_features(poset: Poset) -> dict:
    N = poset.n
    sis = adaptive_sis_runs(N)
    log_H = compute_log_H(poset, n_runs=sis)
    pi_geo = compute_pi_geo(poset)
    sigma_hist = compute_sigma_hist(poset)
    xi_dim_val, d_eff = compute_xi_dim(poset)
    R = compute_R(poset)

    alpha0, q, lam, eta = 16.0, -0.5, 10.0, 0.6
    Rc, w, N0 = 0.25, 0.015, 20.0
    alpha_N = alpha0 * (N0 / max(N, 1)) ** abs(q)
    wall = alpha_N * sigmoid((R - Rc) / w)
    F7 = log_H + 0.0004 * pi_geo - lam * sigma_hist + eta * xi_dim_val + wall

    layer_count = compute_layer_count(poset)

    return {
        "F7": F7,
        "log_H": log_H,
        "pi_geo": pi_geo,
        "sigma_hist": sigma_hist,
        "xi_dim": xi_dim_val,
        "d_eff": d_eff,
        "R": R,
        "wall": wall,
        "layer_count": layer_count,
    }


FAMILIES = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
    "KR_like": generate_kr_like,
}


# ══════════════════════════════════════════════════════════════════════════
# Data collection
# ══════════════════════════════════════════════════════════════════════════

def collect_fresh_data(n_values, reps, seed_base=77):
    rows = []
    total = len(n_values) * len(FAMILIES) * reps
    count = 0
    t0 = time.time()

    for N in n_values:
        for fam_name, gen in FAMILIES.items():
            for rep in range(reps):
                count += 1
                seed = seed_base + N * 1000 + hash(fam_name) % 100 * 10 + rep
                poset = gen(N, seed=seed)
                comp = compute_all_features(poset)
                row = {"family": fam_name, "N": N, "rep": rep, "seed": seed, **comp}
                rows.append(row)

                if count % 25 == 0 or count == total:
                    elapsed = time.time() - t0
                    eta_s = (elapsed / count) * (total - count)
                    print(f"  [{count}/{total}] {fam_name} N={N}: "
                          f"Σ_hist={comp['sigma_hist']:.3f} logH={comp['log_H']:.1f} "
                          f"(elapsed {elapsed:.0f}s, ETA {eta_s:.0f}s)")

    return rows


def load_cached_data(path: str) -> list[dict]:
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


# ══════════════════════════════════════════════════════════════════════════
# Analysis
# ══════════════════════════════════════════════════════════════════════════

def generate_report(rows, n_values, has_layer_count=False):
    by_nf = defaultdict(list)
    for r in rows:
        by_nf[(r["N"], r["family"])].append(r)

    families = sorted(set(r["family"] for r in rows))
    report = []
    report.append("# Prediction C — F7 Large-N Verification\n")
    report.append(f"**Data**: {len(rows)} samples, N ∈ {{{', '.join(str(n) for n in n_values)}}}\n")
    report.append(f"**Families**: {', '.join(families)}\n")

    # ── Q1: Within-cell ρ(Σ_hist, logH) ──
    report.append("\n## Q1: Within-Cell Correlation ρ(Σ_hist, logH)\n")
    report.append("The core Prediction C claim: higher Σ_hist → lower logH within fixed (family, N).\n")
    report.append("| family | N | n | ρ(Σ,logH) | p_value | sig | direction |")
    report.append("|--------|---|---|-----------|---------|-----|-----------|")

    q1_data = []  # for Q2 scaling
    n_sig = 0
    n_total = 0
    n_correct_dir = 0

    for N in n_values:
        for fam in families:
            vals = by_nf.get((N, fam), [])
            if len(vals) < 5:
                continue
            sh = np.array([r["sigma_hist"] for r in vals])
            lh = np.array([r["log_H"] for r in vals])
            if np.std(sh) < 1e-12 or np.std(lh) < 1e-12:
                report.append(f"| {fam} | {N} | {len(vals)} | — | — | — | degenerate |")
                continue

            rho, p = sp_stats.spearmanr(sh, lh)
            sig = "★★★" if p < 0.001 else "★★" if p < 0.01 else "★" if p < 0.05 else ""
            direction = "✅ neg" if rho < 0 else "❌ pos"
            report.append(f"| {fam} | {N} | {len(vals)} | {rho:+.3f} | {p:.4f} | {sig} | {direction} |")

            n_total += 1
            if p < 0.05:
                n_sig += 1
            if rho < 0:
                n_correct_dir += 1
            q1_data.append({"family": fam, "N": N, "rho": rho, "p": p, "n": len(vals)})

    report.append(f"\n**Summary**: {n_correct_dir}/{n_total} cells have ρ < 0 (correct direction), "
                  f"{n_sig}/{n_total} significant at p<0.05\n")

    # ── Q2: N-scaling of within-cell ρ ──
    report.append("\n## Q2: N-Scaling of Within-Cell ρ(Σ_hist, logH)\n")
    report.append("Does the correlation strengthen with N?\n")

    for fam in families:
        fam_data = [d for d in q1_data if d["family"] == fam]
        if len(fam_data) < 3:
            continue
        ns = [d["N"] for d in fam_data]
        rhos = [d["rho"] for d in fam_data]
        if len(set(ns)) < 3:
            continue
        trend_rho, trend_p = sp_stats.spearmanr(ns, rhos)
        trend = "strengthening" if trend_rho < -0.5 else "weakening" if trend_rho > 0.5 else "stable"
        report.append(f"- **{fam}**: ρ(N, ρ_cell) = {trend_rho:+.3f} (p={trend_p:.3f}) — {trend}")
        for d in fam_data:
            report.append(f"  - N={d['N']}: ρ={d['rho']:+.3f} (p={d['p']:.4f})")

    # ── Q3: Variance explained by Σ_hist within-cell ──
    report.append("\n## Q3: Within-Cell R²(logH ~ Σ_hist)\n")
    report.append("| family | N | R² | Pearson_r |")
    report.append("|--------|---|----| ---------|")

    for N in n_values:
        for fam in families:
            vals = by_nf.get((N, fam), [])
            if len(vals) < 5:
                continue
            sh = np.array([r["sigma_hist"] for r in vals])
            lh = np.array([r["log_H"] for r in vals])
            if np.std(sh) < 1e-12:
                continue
            r_pearson, _ = sp_stats.pearsonr(sh, lh)
            r2 = r_pearson ** 2
            report.append(f"| {fam} | {N} | {r2:.3f} | {r_pearson:+.3f} |")

    # ── Q4: Σ_hist vs layer_count ──
    if has_layer_count:
        report.append("\n## Q4: Σ_hist vs layer_count Relationship\n")
        report.append("| family | N | ρ(Σ_hist, layer_count) | p | direction |")
        report.append("|--------|---|------------------------|---|-----------|")

        for N in n_values:
            for fam in families:
                vals = by_nf.get((N, fam), [])
                if len(vals) < 5 or "layer_count" not in vals[0]:
                    continue
                sh = np.array([r["sigma_hist"] for r in vals])
                lc = np.array([r["layer_count"] for r in vals])
                if np.std(sh) < 1e-12 or np.std(lc) < 1e-12:
                    report.append(f"| {fam} | {N} | — | — | degenerate |")
                    continue
                rho, p = sp_stats.spearmanr(sh, lc)
                dir_s = "positive" if rho > 0 else "negative"
                report.append(f"| {fam} | {N} | {rho:+.3f} | {p:.4f} | {dir_s} |")

    # ── Q5: F7 ranking reflects Σ_hist within-cell ──
    report.append("\n## Q5: Does Higher Σ_hist → Lower F7 Within Cells?\n")
    report.append("Since F7 = logH − λ·Σ_hist + ..., partial correlation "
                  "ρ(F7, Σ_hist | other F7 terms) should be strongly negative.\n")
    report.append("| family | N | ρ(F7, Σ_hist) | ρ(logH_resid, Σ_hist) |")
    report.append("|--------|---|---------------|----------------------|")

    for N in n_values:
        for fam in families:
            vals = by_nf.get((N, fam), [])
            if len(vals) < 5:
                continue
            sh = np.array([r["sigma_hist"] for r in vals])
            f7 = np.array([r["F7"] for r in vals])
            lh = np.array([r["log_H"] for r in vals])
            wall = np.array([r["wall"] for r in vals])
            xi = np.array([r["xi_dim"] for r in vals])

            if np.std(sh) < 1e-12:
                continue

            rho_f7, _ = sp_stats.spearmanr(sh, f7)

            # logH_resid = logH - (wall + 0.6*xi) → remove non-Σ terms
            logH_resid = lh - wall - 0.6 * xi
            rho_resid, _ = sp_stats.spearmanr(sh, logH_resid)

            report.append(f"| {fam} | {N} | {rho_f7:+.3f} | {rho_resid:+.3f} |")

    # ── Cross-family pooled analysis ──
    report.append("\n## Q6: Cross-Family Pooled Analysis (per N)\n")
    report.append("Pooled across families but fixed N — tests whether Σ_hist "
                  "carries cross-family information.\n")
    report.append("| N | n | ρ(Σ_hist, logH) | ρ(Σ_hist, F7) | sig |")
    report.append("|---|---|-----------------|---------------|-----|")

    for N in n_values:
        all_vals = []
        for fam in families:
            all_vals.extend(by_nf.get((N, fam), []))
        if len(all_vals) < 10:
            continue
        sh = np.array([r["sigma_hist"] for r in all_vals])
        lh = np.array([r["log_H"] for r in all_vals])
        f7 = np.array([r["F7"] for r in all_vals])

        rho_lh, p_lh = sp_stats.spearmanr(sh, lh)
        rho_f7, p_f7 = sp_stats.spearmanr(sh, f7)
        sig = "★★★" if p_lh < 0.001 else "★★" if p_lh < 0.01 else "★" if p_lh < 0.05 else ""
        report.append(f"| {N} | {len(all_vals)} | {rho_lh:+.3f} | {rho_f7:+.3f} | {sig} |")

    # ── Verdict ──
    report.append("\n## Verdict\n")
    pct_neg = 100 * n_correct_dir / max(n_total, 1)
    pct_sig = 100 * n_sig / max(n_total, 1)

    if pct_neg >= 80 and pct_sig >= 50:
        report.append(f"**STRONG**: {pct_neg:.0f}% cells show correct direction, "
                      f"{pct_sig:.0f}% significant — Prediction C well-supported at large N.")
    elif pct_neg >= 60:
        report.append(f"**MODERATE**: {pct_neg:.0f}% cells show correct direction, "
                      f"{pct_sig:.0f}% significant — partial support.")
    else:
        report.append(f"**WEAK**: Only {pct_neg:.0f}% cells show correct direction.")

    return "\n".join(report)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fresh", action="store_true",
                        help="Generate fresh data with KR + layer_count")
    parser.add_argument("--reps", type=int, default=15)
    parser.add_argument("--ns", type=int, nargs="+", default=None)
    args = parser.parse_args()

    outdir = Path("outputs_d_recovery")
    outdir.mkdir(exist_ok=True)

    if args.fresh:
        n_values = args.ns or [20, 36, 52, 72, 100]
        print(f"=== Prediction C — Fresh Data Generation ===")
        print(f"N = {n_values}, Families = {list(FAMILIES.keys())}, Reps = {args.reps}")
        rows = collect_fresh_data(n_values, args.reps)
        has_layer_count = True

        csv_path = outdir / "prediction_c_f7_large_n.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"Saved: {csv_path}")
    else:
        # Use cached Prediction A data (no KR, no layer_count)
        cached = "outputs_d_recovery/prediction_a_f7_dimension_full.csv"
        print(f"=== Prediction C — Reanalysis of cached data ===")
        print(f"Loading: {cached}")
        rows = load_cached_data(cached)
        n_values = sorted(set(r["N"] for r in rows))
        has_layer_count = False
        print(f"Loaded {len(rows)} rows, N = {n_values}")

    report = generate_report(rows, n_values, has_layer_count)

    md_path = outdir / "prediction_c_f7_large_n.md"
    md_path.write_text(report, encoding="utf-8")
    print(f"Saved: {md_path}")
    print("\n" + report)


if __name__ == "__main__":
    main()
