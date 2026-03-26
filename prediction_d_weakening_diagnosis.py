"""
Prediction D: Diagnose Per-N W₁ Weakening Root Cause
=====================================================
Three hypotheses:
  H1: Recovery steps insufficient at large N (40 vs 120)
  H2: P_basin variance compression (floor/ceiling effects)
  H3: W₁ distribution changes with N (scale effect)

Analysis plan:
  A1: P_basin distribution per N per family — check variance compression
  A2: W₁ distribution per N — check if range/variance changes
  A3: Scatter W₁ vs P_basin per N — visual check
  A4: Recovery step budget analysis — are 40 steps enough at N=100?
  A5: Controlled experiment — rerun N=100 with MORE recovery steps (120)
"""

import csv
import numpy as np
from scipy import stats
from pathlib import Path
import argparse

def load_csv(path):
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = []
        for r in reader:
            row = {}
            for k, v in r.items():
                try:
                    row[k] = float(v)
                except (ValueError, TypeError):
                    row[k] = v
            rows.append(row)
    return rows

def partial_corr(x, y, z):
    """Spearman partial correlation of x,y controlling for z."""
    if len(set(x)) < 3 or len(set(y)) < 3:
        return float('nan'), 1.0
    rx = stats.spearmanr(x, z)[0]
    ry = stats.spearmanr(y, z)[0]
    res_x = np.array(x) - rx * np.array(z)
    res_y = np.array(y) - ry * np.array(z)
    # Rank-based partial
    from scipy.stats import spearmanr
    rho, p = spearmanr(x, z)
    residuals_x = stats.rankdata(x) - rho * stats.rankdata(z)
    rho2, p2 = spearmanr(y, z)
    residuals_y = stats.rankdata(y) - rho2 * stats.rankdata(z)
    return spearmanr(residuals_x, residuals_y)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=str, default="outputs_d_recovery/prediction_d_large_n_r15.csv")
    args = parser.parse_args()

    data = load_csv(args.csv)
    ns = sorted(set(int(r["N"]) for r in data))
    families = sorted(set(r["family"] for r in data))

    report = []
    report.append("# Prediction D: Per-N W₁ Weakening — Root Cause Diagnosis\n")
    report.append(f"**Data**: {args.csv} ({len(data)} samples)\n")

    # ═══════════════════════════════════════════════════════════════
    # A1: P_basin distribution per N
    # ═══════════════════════════════════════════════════════════════
    report.append("\n## A1: P_basin Distribution per N\n")
    report.append("| N | n | mean | std | min | max | n_zero | n_one | unique_vals |")
    report.append("|---|---|------|-----|-----|-----|--------|-------|-------------|")
    
    pb_stats = {}
    for n_val in ns:
        subset = [r for r in data if int(r["N"]) == n_val]
        pbs = [r["P_basin"] for r in subset]
        n_zero = sum(1 for p in pbs if p == 0.0)
        n_one = sum(1 for p in pbs if p >= 1.0)
        unique = len(set(f"{p:.4f}" for p in pbs))
        pb_stats[n_val] = {
            "mean": np.mean(pbs), "std": np.std(pbs),
            "min": min(pbs), "max": max(pbs),
            "n_zero": n_zero, "n_one": n_one, "unique": unique,
            "values": pbs
        }
        report.append(f"| {n_val} | {len(subset)} | {np.mean(pbs):.3f} | {np.std(pbs):.3f} | "
                      f"{min(pbs):.3f} | {max(pbs):.3f} | {n_zero} | {n_one} | {unique} |")

    # P_basin per family per N
    report.append("\n### P_basin per Family per N\n")
    report.append("| N | family | n | mean_Pb | std_Pb | n_zero | n_nonzero |")
    report.append("|---|--------|---|---------|--------|--------|-----------|")
    for n_val in ns:
        for fam in families:
            subset = [r for r in data if int(r["N"]) == n_val and r["family"] == fam]
            pbs = [r["P_basin"] for r in subset]
            n_zero = sum(1 for p in pbs if p == 0.0)
            report.append(f"| {n_val} | {fam} | {len(subset)} | {np.mean(pbs):.3f} | "
                          f"{np.std(pbs):.3f} | {n_zero} | {len(subset)-n_zero} |")

    # ═══════════════════════════════════════════════════════════════
    # A2: W₁ distribution per N
    # ═══════════════════════════════════════════════════════════════
    report.append("\n## A2: W₁ Distribution per N\n")
    report.append("| N | mean_W1 | std_W1 | min_W1 | max_W1 | range | CV |")
    report.append("|---|---------|--------|--------|--------|-------|-----|")
    for n_val in ns:
        subset = [r for r in data if int(r["N"]) == n_val]
        w1s = [r["wasserstein"] for r in subset]
        mn = np.mean(w1s)
        sd = np.std(w1s)
        report.append(f"| {n_val} | {mn:.4f} | {sd:.4f} | {min(w1s):.4f} | {max(w1s):.4f} | "
                      f"{max(w1s)-min(w1s):.4f} | {sd/mn:.3f} |")

    # ═══════════════════════════════════════════════════════════════
    # A3: Variance decomposition — what fraction of P_basin variance
    #     is explained by family alone vs W₁?
    # ═══════════════════════════════════════════════════════════════
    report.append("\n## A3: Variance Decomposition per N\n")
    report.append("Check: does W₁ add information beyond family membership?\n")
    report.append("| N | ρ(W1,Pb) | ρ(W1,Pb\\|F7) | ρ(W1,Pb) within-Lor5D | Pb_var | Pb_var_nonLor5D |")
    report.append("|---|----------|-------------|----------------------|--------|-----------------|")
    
    for n_val in ns:
        subset = [r for r in data if int(r["N"]) == n_val]
        w1 = [r["wasserstein"] for r in subset]
        pb = [r["P_basin"] for r in subset]
        f7 = [r["F7"] for r in subset]
        
        rho_raw, p_raw = stats.spearmanr(w1, pb)
        rho_part, p_part = partial_corr(w1, pb, f7)
        
        # Within Lor5D only
        lor5d = [r for r in subset if r["family"] == "Lor5D"]
        if len(lor5d) > 3:
            w1_5d = [r["wasserstein"] for r in lor5d]
            pb_5d = [r["P_basin"] for r in lor5d]
            rho_5d, p_5d = stats.spearmanr(w1_5d, pb_5d)
            rho_5d_str = f"{rho_5d:+.3f} (p={p_5d:.3f})"
        else:
            rho_5d_str = "N/A"
        
        # P_basin variance for non-Lor5D
        non5d = [r for r in subset if r["family"] != "Lor5D"]
        pb_non5d = [r["P_basin"] for r in non5d]
        
        report.append(f"| {n_val} | {rho_raw:+.3f} | {rho_part:+.3f} | {rho_5d_str} | "
                      f"{np.var(pb):.4f} | {np.var(pb_non5d):.4f} |")

    # ═══════════════════════════════════════════════════════════════
    # A4: P_basin floor/ceiling analysis
    # ═══════════════════════════════════════════════════════════════
    report.append("\n## A4: P_basin Floor/Ceiling Effect\n")
    report.append("Fraction of samples with P_basin = 0 or P_basin ≥ 0.5 per N:\n")
    report.append("| N | frac_zero | frac_low(<0.17) | frac_mid | frac_high(≥0.5) |")
    report.append("|---|-----------|-----------------|----------|-----------------|")
    
    for n_val in ns:
        subset = [r for r in data if int(r["N"]) == n_val]
        pbs = [r["P_basin"] for r in subset]
        n = len(pbs)
        frac_zero = sum(1 for p in pbs if p == 0.0) / n
        frac_low = sum(1 for p in pbs if 0 < p < 0.17) / n
        frac_high = sum(1 for p in pbs if p >= 0.5) / n
        frac_mid = 1 - frac_zero - frac_low - frac_high
        report.append(f"| {n_val} | {frac_zero:.1%} | {frac_low:.1%} | {frac_mid:.1%} | {frac_high:.1%} |")

    # ═══════════════════════════════════════════════════════════════
    # A5: Information content — mutual information proxy
    # ═══════════════════════════════════════════════════════════════
    report.append("\n## A5: Per-N Spearman Correlations (All Pairs)\n")
    report.append("| N | ρ(W1,Pb) | ρ(ΔH,Pb) | ρ(F7,Pb) | ρ(W1,F7) | ρ(W1,ΔH) |")
    report.append("|---|----------|----------|----------|----------|----------|")
    
    for n_val in ns:
        subset = [r for r in data if int(r["N"]) == n_val]
        w1 = [r["wasserstein"] for r in subset]
        dh = [r["delta_H_int"] for r in subset]
        pb = [r["P_basin"] for r in subset]
        f7 = [r["F7"] for r in subset]
        
        r1, _ = stats.spearmanr(w1, pb)
        r2, _ = stats.spearmanr(dh, pb)
        r3, _ = stats.spearmanr(f7, pb)
        r4, _ = stats.spearmanr(w1, f7)
        r5, _ = stats.spearmanr(w1, dh)
        
        report.append(f"| {n_val} | {r1:+.3f} | {r2:+.3f} | {r3:+.3f} | {r4:+.3f} | {r5:+.3f} |")

    # ═══════════════════════════════════════════════════════════════
    # A6: Recovery step analysis — compare actual recovery outcomes
    # ═══════════════════════════════════════════════════════════════
    report.append("\n## A6: Recovery Step Budget vs Outcome\n")
    report.append("Adaptive recovery steps: N=36→120, N=52→80, N=72→60, N=100→40\n")
    report.append("| N | steps | Lor4D_Pb | Lor5D_Pb | Lor4D_spread | Lor5D_spread |")
    report.append("|---|-------|----------|----------|-------------|-------------|")
    
    step_map = {36: 120, 52: 80, 72: 60, 100: 40}
    for n_val in ns:
        steps = step_map.get(n_val, "?")
        for fam in ["Lor4D", "Lor5D"]:
            subset = [r for r in data if int(r["N"]) == n_val and r["family"] == fam]
            pbs = [r["P_basin"] for r in subset]
        lor4d = [r for r in data if int(r["N"]) == n_val and r["family"] == "Lor4D"]
        lor5d = [r for r in data if int(r["N"]) == n_val and r["family"] == "Lor5D"]
        pb4 = [r["P_basin"] for r in lor4d]
        pb5 = [r["P_basin"] for r in lor5d]
        report.append(f"| {n_val} | {steps} | {np.mean(pb4):.1%} | {np.mean(pb5):.1%} | "
                      f"{np.std(pb4):.3f} | {np.std(pb5):.3f} |")

    # ═══════════════════════════════════════════════════════════════
    # DIAGNOSIS: Which hypothesis is supported?
    # ═══════════════════════════════════════════════════════════════
    report.append("\n## Diagnosis Summary\n")
    
    # Check H2: variance compression
    pb_vars = [pb_stats[n]["std"] for n in ns]
    rho_var, p_var = stats.spearmanr(ns, pb_vars)
    report.append(f"**H2 (variance compression)**: Pb_std trend with N: ρ={rho_var:+.3f} (p={p_var:.3f})")
    if rho_var < -0.5:
        report.append("  → **SUPPORTED**: P_basin variance decreases with N\n")
    else:
        report.append("  → NOT supported: P_basin variance does not decrease with N\n")
    
    # Check H3: W₁ distribution change
    w1_cvs = []
    for n_val in ns:
        subset = [r for r in data if int(r["N"]) == n_val]
        w1s = [r["wasserstein"] for r in subset]
        w1_cvs.append(np.std(w1s) / np.mean(w1s))
    rho_cv, p_cv = stats.spearmanr(ns, w1_cvs)
    report.append(f"**H3 (W₁ scale effect)**: W1_CV trend with N: ρ={rho_cv:+.3f} (p={p_cv:.3f})")
    if abs(rho_cv) > 0.5:
        report.append(f"  → W₁ distribution {'widens' if rho_cv > 0 else 'narrows'} with N\n")
    else:
        report.append("  → W₁ CV stable across N\n")

    # Check: fraction of P_basin=0
    frac_zeros = []
    for n_val in ns:
        subset = [r for r in data if int(r["N"]) == n_val]
        pbs = [r["P_basin"] for r in subset]
        frac_zeros.append(sum(1 for p in pbs if p == 0.0) / len(pbs))
    rho_fz, p_fz = stats.spearmanr(ns, frac_zeros)
    report.append(f"**Floor effect**: frac(Pb=0) trend with N: ρ={rho_fz:+.3f} (p={p_fz:.3f})")
    report.append(f"  Values: {[f'{fz:.1%}' for fz in frac_zeros]}")
    if rho_fz > 0.5:
        report.append("  → **SUPPORTED**: More samples hit P_basin=0 floor at large N\n")
    else:
        report.append("  → Floor effect does not worsen with N\n")

    # Print and save
    report_text = "\n".join(report)
    print(report_text)
    
    outdir = Path("outputs_d_recovery")
    outdir.mkdir(exist_ok=True)
    with open(outdir / "prediction_d_weakening_diagnosis.md", "w", encoding="utf-8") as f:
        f.write(report_text)
    print(f"\nSaved to {outdir / 'prediction_d_weakening_diagnosis.md'}")

if __name__ == "__main__":
    main()
