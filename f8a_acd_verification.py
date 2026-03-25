"""F8a Impact Verification on Predictions A, C, D.

Uses cached data from prediction_b_f7_scaling.csv (which has all 5 families
at N=20-72) to recompute F8a and check:
  A: Lor4D is local minimum (vs Lor2D/3D/5D ordering)
  B: Already supported at 100% pass rate
  C: corr(Σ_hist, F8a) < 0 within Lor2D (deeper layers → lower cost)
  D: Needs coarse-graining data — will generate fresh if needed

Also generates fresh data at multiple N values for A/C tests.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats

from bd_action import count_intervals_fast
from generators import (
    Poset,
    generate_kr_like,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
)
from unified_functional import (
    compute_log_H,
    compute_pi_geo,
    compute_sigma_hist,
    compute_xi_dim,
)
from coarse_grain import coarse_grain_delete_nodes


def sigmoid(x: float) -> float:
    x = max(-30.0, min(30.0, x))
    return 1.0 / (1.0 + math.exp(-x))


def compute_R(poset: Poset) -> float:
    counts = count_intervals_fast(poset, k_max=3)
    if counts.total_relations <= 0:
        return 0.0
    return 1.0 - float(counts.get(0)) / float(counts.total_relations)


def compute_F7(row: dict) -> float:
    N = row["N"]
    alpha_N = 16.0 * (20.0 / max(N, 1)) ** 0.5
    wall = alpha_N * sigmoid((row["R"] - 0.25) / 0.015)
    return row["log_H"] + 0.0004 * row["pi_geo"] - 10.0 * row["sigma_hist"] + 0.6 * row["xi_dim"] + wall


def compute_F8a(row: dict) -> float:
    N = row["N"]
    nlogn = N * math.log(N) if N > 1 else 1.0
    log_H_norm = row["log_H"] / nlogn
    alpha_N = 16.0 * (20.0 / max(N, 1)) ** 0.5
    wall = alpha_N * sigmoid((row["R"] - 0.25) / 0.015)
    return log_H_norm + 0.0004 * row["pi_geo"] - 10.0 * row["sigma_hist"] + 0.6 * row["xi_dim"] + wall


def load_cached() -> list[dict]:
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


def generate_fresh(n_values: list[int], reps: int, seed: int) -> list[dict]:
    """Generate fresh data with all components for A/C/D tests."""
    families = {
        "Lor2D": generate_lorentzian_like_2d,
        "Lor3D": generate_lorentzian_like_3d,
        "Lor4D": generate_lorentzian_like_4d,
        "Lor5D": generate_lorentzian_like_5d,
        "KR_like": generate_kr_like,
    }
    rows = []
    total = len(n_values) * len(families) * reps
    count = 0
    for N in n_values:
        for fam_name, gen in families.items():
            for rep in range(reps):
                count += 1
                s = seed + rep * 1000 + N * 100
                poset = gen(N, seed=s)
                log_H = compute_log_H(poset, n_runs=512)
                pi_geo = compute_pi_geo(poset)
                sigma_hist = compute_sigma_hist(poset)
                xi_dim_val, d_eff = compute_xi_dim(poset)
                R = compute_R(poset)
                
                rows.append({
                    "family": fam_name,
                    "N": N,
                    "rep": rep,
                    "log_H": log_H,
                    "pi_geo": pi_geo,
                    "sigma_hist": sigma_hist,
                    "xi_dim": xi_dim_val,
                    "d_eff": d_eff,
                    "R": R,
                })
                if count % 50 == 0 or count == total:
                    print(f"  [{count}/{total}] {fam_name} N={N} rep={rep}")
    return rows


def test_prediction_A(data: list[dict], model_fn, model_name: str) -> list[str]:
    """Test Prediction A: dimensional selection.
    
    A is tested as pairwise win rates:
    - Lor4D < Lor2D (vs 2D)
    - Lor4D < Lor3D (vs 3D)  
    - Lor4D < Lor5D (vs 5D)
    
    These are the A sub-metrics from the manuscript.
    """
    lines = [f"### Prediction A — {model_name}\n"]
    
    by_nf: dict[tuple[int, str], list[dict]] = defaultdict(list)
    for r in data:
        by_nf[(r["N"], r["family"])].append(r)
    
    n_values = sorted(set(r["N"] for r in data))
    dim_families = ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]
    
    # Mean F per (N, family)
    lines.append("| N | Lor2D | Lor3D | Lor4D | Lor5D | KR | Lor4D rank |")
    lines.append("|---|-------|-------|-------|-------|-----|------------|")
    
    a_wins_total = 0
    a_total = 0
    
    for N in n_values:
        mean_fs = {}
        for fam in dim_families + ["KR_like"]:
            grp = by_nf.get((N, fam), [])
            if grp:
                mean_fs[fam] = np.mean([model_fn(r) for r in grp])
        
        if len(mean_fs) < 4:
            continue
        
        # Lor4D rank among Lorentzian families
        lor_fs = {k: v for k, v in mean_fs.items() if k.startswith("Lor")}
        sorted_lor = sorted(lor_fs.items(), key=lambda x: x[1])
        lor4d_rank = [i for i, (k, _) in enumerate(sorted_lor) if k == "Lor4D"][0] + 1
        
        lines.append(f"| {N} | {mean_fs.get('Lor2D', 0):.2f} | {mean_fs.get('Lor3D', 0):.2f} | "
                     f"{mean_fs.get('Lor4D', 0):.2f} | {mean_fs.get('Lor5D', 0):.2f} | "
                     f"{mean_fs.get('KR_like', 0):.2f} | {lor4d_rank}/4 |")
        
        # A pairwise: Lor4D vs Lor2D, vs Lor3D, vs Lor5D
        for other in ["Lor2D", "Lor3D", "Lor5D"]:
            l4_rows = by_nf.get((N, "Lor4D"), [])
            oth_rows = by_nf.get((N, other), [])
            if l4_rows and oth_rows:
                n_pairs = min(len(l4_rows), len(oth_rows))
                wins = sum(1 for i in range(n_pairs) 
                          if model_fn(l4_rows[i]) < model_fn(oth_rows[i]))
                a_wins_total += wins
                a_total += n_pairs
    
    a_rate = 100 * a_wins_total / max(a_total, 1)
    lines.append(f"\n**A overall (Lor4D < others): {a_wins_total}/{a_total} ({a_rate:.1f}%)**\n")
    
    # Barrier asymmetry check per N
    lines.append("| N | ΔF(2→3) | ΔF(3→4) | ΔF(4→5) | Barrier↑ at 4? |")
    lines.append("|---|---------|---------|---------|----------------|")
    
    for N in n_values:
        mean_fs = {}
        for fam in dim_families:
            grp = by_nf.get((N, fam), [])
            if grp:
                mean_fs[fam] = np.mean([model_fn(r) for r in grp])
        
        if len(mean_fs) < 4:
            continue
        
        d23 = mean_fs["Lor3D"] - mean_fs["Lor2D"]
        d34 = mean_fs["Lor4D"] - mean_fs["Lor3D"]
        d45 = mean_fs["Lor5D"] - mean_fs["Lor4D"]
        barrier = "✅" if d45 > d34 else "❌"
        lines.append(f"| {N} | {d23:+.2f} | {d34:+.2f} | {d45:+.2f} | {barrier} |")
    
    return lines


def test_prediction_C(data: list[dict], model_fn, model_name: str) -> list[str]:
    """Test Prediction C: deeper hierarchy → lower F (within Lor2D).
    
    C = corr(Σ_hist, F) should be negative.
    """
    lines = [f"### Prediction C — {model_name}\n"]
    
    by_n: dict[int, list[dict]] = defaultdict(list)
    for r in data:
        if r["family"] == "Lor2D":
            by_n[r["N"]].append(r)
    
    lines.append("| N | n_samples | corr(Σ_hist, F) | p-value | C direction |")
    lines.append("|---|-----------|-----------------|---------|-------------|")
    
    all_c_ok = 0
    all_c_total = 0
    
    for N in sorted(by_n.keys()):
        grp = by_n[N]
        if len(grp) < 5:
            continue
        
        sh = np.array([r["sigma_hist"] for r in grp])
        fs = np.array([model_fn(r) for r in grp])
        
        if sh.std() < 1e-10:
            lines.append(f"| {N} | {len(grp)} | – | – | no variance |")
            continue
        
        rho, pval = sp_stats.spearmanr(sh, fs)
        direction = "✅ (↑Σ → ↓F)" if rho < 0 else "❌"
        lines.append(f"| {N} | {len(grp)} | {rho:+.3f} | {pval:.3e} | {direction} |")
        
        all_c_total += 1
        if rho < 0:
            all_c_ok += 1
    
    c_rate = all_c_ok / max(all_c_total, 1)
    lines.append(f"\n**C pass rate: {all_c_ok}/{all_c_total} ({100*c_rate:.0f}%)**\n")
    return lines


def test_prediction_D(n_values: list[int], reps: int, model_fn, model_name: str, seed: int) -> list[str]:
    """Test Prediction D: coarse-graining stability.
    
    |ΔF| after 30% node deletion should be smaller for Lor than KR.
    """
    lines = [f"### Prediction D — {model_name}\n"]
    
    families = {
        "Lor2D": generate_lorentzian_like_2d,
        "Lor4D": generate_lorentzian_like_4d,
        "KR_like": generate_kr_like,
    }
    
    results: dict[tuple[int, str], list[float]] = defaultdict(list)
    total = len(n_values) * len(families) * reps
    count = 0
    
    for N in n_values:
        for fam_name, gen in families.items():
            for rep in range(reps):
                count += 1
                s = seed + rep * 1000 + N * 100
                poset = gen(N, seed=s)
                
                # Compute components for original
                log_H = compute_log_H(poset, n_runs=256)
                pi_geo = compute_pi_geo(poset)
                sigma_hist = compute_sigma_hist(poset)
                xi_dim_val, d_eff = compute_xi_dim(poset)
                R = compute_R(poset)
                
                row_orig = {"N": N, "log_H": log_H, "pi_geo": pi_geo,
                           "sigma_hist": sigma_hist, "xi_dim": xi_dim_val,
                           "d_eff": d_eff, "R": R, "family": fam_name}
                f_orig = model_fn(row_orig)
                
                # Coarse-grain
                cg_poset = coarse_grain_delete_nodes(poset, keep_ratio=0.7, seed=s + 99)
                N_cg = cg_poset.n
                log_H_cg = compute_log_H(cg_poset, n_runs=256)
                pi_geo_cg = compute_pi_geo(cg_poset)
                sigma_hist_cg = compute_sigma_hist(cg_poset)
                xi_dim_cg, d_eff_cg = compute_xi_dim(cg_poset)
                R_cg = compute_R(cg_poset)
                
                row_cg = {"N": N_cg, "log_H": log_H_cg, "pi_geo": pi_geo_cg,
                         "sigma_hist": sigma_hist_cg, "xi_dim": xi_dim_cg,
                         "d_eff": d_eff_cg, "R": R_cg, "family": fam_name}
                f_cg = model_fn(row_cg)
                
                delta = abs(f_cg - f_orig)
                results[(N, fam_name)].append(delta)
                
                if count % 30 == 0 or count == total:
                    print(f"  D [{count}/{total}] {fam_name} N={N}")
    
    lines.append("| N | Lor2D |ΔF| | Lor4D |ΔF| | KR |ΔF| | Lor2D<KR? | Lor4D<KR? |")
    lines.append("|---|---------|---------|---------|-----------|-----------|")
    
    d_pass = 0
    d_total = 0
    
    for N in n_values:
        l2 = np.mean(results.get((N, "Lor2D"), [0]))
        l4 = np.mean(results.get((N, "Lor4D"), [0]))
        kr = np.mean(results.get((N, "KR_like"), [0]))
        l2_ok = "✅" if l2 < kr else "❌"
        l4_ok = "✅" if l4 < kr else "❌"
        lines.append(f"| {N} | {l2:.3f} | {l4:.3f} | {kr:.3f} | {l2_ok} | {l4_ok} |")
        d_total += 2
        if l2 < kr: d_pass += 1
        if l4 < kr: d_pass += 1
    
    lines.append(f"\n**D pass rate: {d_pass}/{d_total} ({100*d_pass/max(d_total,1):.0f}%)**\n")
    return lines


def main():
    print("=" * 70)
    print("F8a Impact Verification on Predictions A, C, D")
    print("=" * 70)
    
    # Load cached data (has N=20-72, 5 families, 10 reps each)
    cached = load_cached()
    print(f"Loaded {len(cached)} cached rows")
    
    # Also generate fresh data at smaller N for broader coverage
    print("\nGenerating fresh data at N=16,20,28,36...")
    fresh = generate_fresh([16, 20, 28, 36], reps=10, seed=42)
    
    all_data = cached + fresh
    print(f"Total data: {len(all_data)} rows")
    
    report_lines: list[str] = ["# F8a Impact Verification: Predictions A, C, D\n"]
    
    models = {
        "F7 (baseline)": compute_F7,
        "F8a (logH/NlogN)": compute_F8a,
    }
    
    # ── Prediction A ──
    report_lines.append("## Prediction A — Dimensional Selection\n")
    for name, fn in models.items():
        report_lines.extend(test_prediction_A(all_data, fn, name))
    
    # ── Prediction C ──
    report_lines.append("\n## Prediction C — Historical Sedimentation\n")
    for name, fn in models.items():
        report_lines.extend(test_prediction_C(all_data, fn, name))
    
    # ── Prediction D (fresh generation needed) ──
    report_lines.append("\n## Prediction D — Coarse-Graining Stability\n")
    d_ns = [16, 20, 28]
    d_reps = 10
    for name, fn in models.items():
        print(f"\nRunning Prediction D for {name}...")
        report_lines.extend(test_prediction_D(d_ns, d_reps, fn, name, seed=77777))
    
    # ── Summary ──
    report_lines.append("\n## Summary\n")
    report_lines.append("| Prediction | F7 | F8a | Change |")
    report_lines.append("|------------|-----|------|--------|")
    
    report = "\n".join(report_lines) + "\n"
    out_path = Path("outputs_unified_functional/f8a_acd_verification.md")
    out_path.write_text(report, encoding="utf-8")
    print(f"\nReport: {out_path}")


if __name__ == "__main__":
    main()
