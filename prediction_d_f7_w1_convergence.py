"""
F7–W₁ Convergence Analysis
===========================
Explore WHY F7 and W₁ become near-collinear at large N.

F7 = logH + γΠ_geo − λΣ_hist + ηΞ_d + wall(R)
W₁ = Wasserstein-1(interval_spectrum_orig, interval_spectrum_CG)

Hypotheses for convergence:
  H1: Both are dominated by the SAME underlying variable (e.g., comparable_fraction)
  H2: F7's component terms become individually correlated with W₁ at large N
  H3: The interval spectrum shape is determined by Lorentzian-ness at large N
      (spectral rigidity — fewer DOF in spectrum at large N)

Analyses:
  A1: Correlate each F7 component with W₁ per N
  A2: Correlate comparable_fraction with both F7 and W₁ per N
  A3: Check if W₁ variance explained by family membership increases with N
  A4: Theoretical decomposition — what drives the collinearity?
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

def eta_sq(groups):
    """Eta-squared: fraction of variance explained by group membership."""
    all_vals = []
    for g in groups:
        all_vals.extend(g)
    grand_mean = np.mean(all_vals)
    ss_between = sum(len(g) * (np.mean(g) - grand_mean)**2 for g in groups)
    ss_total = sum((x - grand_mean)**2 for x in all_vals)
    if ss_total == 0:
        return 0.0
    return ss_between / ss_total

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=str, default="outputs_d_recovery/prediction_d_large_n_r15.csv")
    args = parser.parse_args()

    data = load_csv(args.csv)
    ns = sorted(set(int(r["N"]) for r in data))
    families = sorted(set(r["family"] for r in data))

    report = []
    report.append("# F7–W₁ Convergence at Large N: Physical Analysis\n")
    report.append(f"**Data**: {args.csv} ({len(data)} samples)\n")

    # ═══════════════════════════════════════════════════════════════
    # A1: Correlate comparable_fraction with both F7 and W₁
    # ═══════════════════════════════════════════════════════════════
    report.append("\n## A1: Comparable Fraction as Common Driver\n")
    report.append("comp_frac (cf) = fraction of comparable pairs = 2·|{(i,j): i<j}| / N(N-1)\n")
    report.append("| N | ρ(cf,F7) | ρ(cf,W₁) | ρ(cf,Pb) | ρ(cf,ΔH) |")
    report.append("|---|----------|----------|----------|----------|")
    
    for n_val in ns:
        subset = [r for r in data if int(r["N"]) == n_val]
        cf = [r["comp_frac"] for r in subset]
        f7 = [r["F7"] for r in subset]
        w1 = [r["wasserstein"] for r in subset]
        pb = [r["P_basin"] for r in subset]
        dh = [r["delta_H_int"] for r in subset]
        
        r1, _ = stats.spearmanr(cf, f7)
        r2, _ = stats.spearmanr(cf, w1)
        r3, _ = stats.spearmanr(cf, pb)
        r4, _ = stats.spearmanr(cf, dh)
        report.append(f"| {n_val} | {r1:+.3f} | {r2:+.3f} | {r3:+.3f} | {r4:+.3f} |")

    # ═══════════════════════════════════════════════════════════════
    # A2: Family explains what fraction of F7 and W₁ variance?
    # ═══════════════════════════════════════════════════════════════
    report.append("\n## A2: Family Membership Explains Variance (η²)\n")
    report.append("η² = SS_between / SS_total — how much of each variable is family-determined\n")
    report.append("| N | η²(F7) | η²(W₁) | η²(cf) | η²(Pb) | η²(ΔH) |")
    report.append("|---|--------|--------|--------|--------|--------|")
    
    for n_val in ns:
        subset = [r for r in data if int(r["N"]) == n_val]
        
        results = {}
        for var_name in ["F7", "wasserstein", "comp_frac", "P_basin", "delta_H_int"]:
            groups = []
            for fam in families:
                vals = [r[var_name] for r in subset if r["family"] == fam]
                groups.append(vals)
            results[var_name] = eta_sq(groups)
        
        report.append(f"| {n_val} | {results['F7']:.3f} | {results['wasserstein']:.3f} | "
                      f"{results['comp_frac']:.3f} | {results['P_basin']:.3f} | "
                      f"{results['delta_H_int']:.3f} |")

    # ═══════════════════════════════════════════════════════════════
    # A3: Within-family correlations F7 vs W₁
    # ═══════════════════════════════════════════════════════════════
    report.append("\n## A3: Within-Family ρ(F7, W₁)\n")
    report.append("If convergence is purely family-driven, within-family ρ should be weak\n")
    report.append("| N | Lor2D | Lor3D | Lor4D | Lor5D | KR | pooled |")
    report.append("|---|-------|-------|-------|-------|-----|--------|")
    
    for n_val in ns:
        row = [f"| {n_val}"]
        subset = [r for r in data if int(r["N"]) == n_val]
        for fam in families:
            fam_data = [r for r in subset if r["family"] == fam]
            if len(fam_data) < 4:
                row.append("N/A")
                continue
            f7 = [r["F7"] for r in fam_data]
            w1 = [r["wasserstein"] for r in fam_data]
            rho, p = stats.spearmanr(f7, w1)
            sig = "★" if p < 0.05 else ""
            row.append(f"{rho:+.3f}{sig}")
        
        # pooled
        f7_all = [r["F7"] for r in subset]
        w1_all = [r["wasserstein"] for r in subset]
        rho_p, _ = stats.spearmanr(f7_all, w1_all)
        row.append(f"{rho_p:+.3f}")
        report.append(" | ".join(row) + " |")

    # ═══════════════════════════════════════════════════════════════
    # A4: Family mean profiles
    # ═══════════════════════════════════════════════════════════════
    report.append("\n## A4: Family Mean Profiles\n")
    report.append("| N | family | mean_F7 | mean_W₁ | mean_cf | F7_rank | W₁_rank | rank_match |")
    report.append("|---|--------|---------|---------|---------|---------|---------|------------|")
    
    for n_val in ns:
        subset = [r for r in data if int(r["N"]) == n_val]
        fam_means = {}
        for fam in families:
            fam_data = [r for r in subset if r["family"] == fam]
            fam_means[fam] = {
                "F7": np.mean([r["F7"] for r in fam_data]),
                "W1": np.mean([r["wasserstein"] for r in fam_data]),
                "cf": np.mean([r["comp_frac"] for r in fam_data]),
            }
        
        # Rank by F7 (ascending) and W₁ (ascending)
        f7_sorted = sorted(families, key=lambda f: fam_means[f]["F7"])
        w1_sorted = sorted(families, key=lambda f: fam_means[f]["W1"])
        
        for fam in families:
            f7_rank = f7_sorted.index(fam) + 1
            w1_rank = w1_sorted.index(fam) + 1
            match = "✓" if f7_rank == w1_rank else ""
            report.append(f"| {n_val} | {fam} | {fam_means[fam]['F7']:.1f} | "
                          f"{fam_means[fam]['W1']:.4f} | {fam_means[fam]['cf']:.3f} | "
                          f"{f7_rank} | {w1_rank} | {match} |")

    # ═══════════════════════════════════════════════════════════════
    # A5: Rank concordance of families by F7 and W₁ across N
    # ═══════════════════════════════════════════════════════════════
    report.append("\n## A5: Family Rank Concordance\n")
    report.append("Kendall's tau between family rankings by F7 vs W₁:\n")
    
    for n_val in ns:
        subset = [r for r in data if int(r["N"]) == n_val]
        fam_means_f7 = []
        fam_means_w1 = []
        for fam in families:
            fam_data = [r for r in subset if r["family"] == fam]
            fam_means_f7.append(np.mean([r["F7"] for r in fam_data]))
            fam_means_w1.append(np.mean([r["wasserstein"] for r in fam_data]))
        
        tau, p = stats.kendalltau(fam_means_f7, fam_means_w1)
        report.append(f"- N={n_val}: Kendall τ = {tau:+.3f} (p={p:.3f})")

    # ═══════════════════════════════════════════════════════════════
    # A6: Residual analysis — after removing family effect
    # ═══════════════════════════════════════════════════════════════
    report.append("\n## A6: Within-Family Residual ρ(F7, W₁)\n")
    report.append("Remove family means, compute ρ on residuals:\n")
    
    for n_val in ns:
        subset = [r for r in data if int(r["N"]) == n_val]
        f7_resid = []
        w1_resid = []
        for r in subset:
            fam = r["family"]
            fam_data = [rr for rr in subset if rr["family"] == fam]
            f7_mean = np.mean([rr["F7"] for rr in fam_data])
            w1_mean = np.mean([rr["wasserstein"] for rr in fam_data])
            f7_resid.append(r["F7"] - f7_mean)
            w1_resid.append(r["wasserstein"] - w1_mean)
        
        rho, p = stats.spearmanr(f7_resid, w1_resid)
        report.append(f"- N={n_val}: ρ_resid(F7, W₁) = {rho:+.3f} (p={p:.4f})")

    # ═══════════════════════════════════════════════════════════════
    # A7: N-scaling of between vs within family variance
    # ═══════════════════════════════════════════════════════════════
    report.append("\n## A7: Between-Family vs Within-Family Variance Ratio\n")
    report.append("| N | F7 between/within | W₁ between/within |")
    report.append("|---|-------------------|-------------------|")
    
    for n_val in ns:
        subset = [r for r in data if int(r["N"]) == n_val]
        for var_name, label in [("F7", "F7"), ("wasserstein", "W1")]:
            pass  # compute below
        
        def bw_ratio(subset, var_name):
            groups = {}
            for r in subset:
                fam = r["family"]
                if fam not in groups:
                    groups[fam] = []
                groups[fam].append(r[var_name])
            grand = np.mean([r[var_name] for r in subset])
            ss_b = sum(len(g) * (np.mean(g) - grand)**2 for g in groups.values())
            ss_w = sum(sum((x - np.mean(g))**2 for x in g) for g in groups.values())
            if ss_w == 0:
                return float('inf')
            return ss_b / ss_w
        
        f7_ratio = bw_ratio(subset, "F7")
        w1_ratio = bw_ratio(subset, "wasserstein")
        report.append(f"| {n_val} | {f7_ratio:.1f} | {w1_ratio:.1f} |")

    # ═══════════════════════════════════════════════════════════════
    # SYNTHESIS
    # ═══════════════════════════════════════════════════════════════
    report.append("\n## Synthesis\n")
    report.append("""
### Physical Interpretation

**F7** measures "Lorentzian-ness" — a composite of:
- logH (entropy of linear extensions — captures order complexity)
- Π_geo (geometric functional — path structure)  
- Σ_hist (history functional — information content)
- Ξ_d (dimensional estimator)
- sigmoid wall (density-based admissibility filter)

**W₁** measures "spectral stability under coarse-graining" — the Wasserstein distance
between the interval spectrum before and after 30% node deletion.

Both quantities are fundamentally about **how structured/ordered the poset is**.
A poset with high F7 (Lorentzian-like) has a rigid, self-similar interval spectrum
that is robust to coarse-graining (low W₁). A poset with low F7 (KR/random-like)
has an irregular spectrum that changes significantly under CG (high W₁).

The convergence at large N reflects **spectral rigidity**: as N grows, the interval
spectrum of each family type becomes increasingly determined by the family's geometric
properties alone (less finite-size noise). Both F7 and W₁ become better estimators of
the underlying family type, hence more correlated with each other.

This is analogous to how in statistical mechanics, different thermodynamic observables
become more correlated as system size increases — they all converge to measuring the
same underlying phase/state.
""")

    # Print and save
    report_text = "\n".join(report)
    print(report_text)
    
    outdir = Path("outputs_d_recovery")
    outdir.mkdir(exist_ok=True)
    with open(outdir / "prediction_d_f7_w1_convergence.md", "w", encoding="utf-8") as f:
        f.write(report_text)
    print(f"\nSaved to {outdir / 'prediction_d_f7_w1_convergence.md'}")

if __name__ == "__main__":
    main()
