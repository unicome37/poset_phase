"""
Prediction C — Lor5D Weakness Diagnosis
========================================
Diagnose why Lor5D reverses the ρ(Σ_hist, logH) correlation at N=52/72.

Hypotheses:
  H1: Σ_hist has near-zero variance in Lor5D (too few distinct values)
  H2: logH variance is dominated by non-layer structural noise
  H3: Lor5D generators produce near-degenerate layer structures
  H4: The reversal is a finite-N artifact that disappears at larger N

Tests:
  T1: Distribution of Σ_hist per (family, N) — unique values, variance, range
  T2: Distribution of layer_count per (family, N)
  T3: logH conditional on layer_count — does the ordering hold?
  T4: Component decomposition: what drives the reversal in Lor5D?
"""
from __future__ import annotations

import csv
import math
from collections import defaultdict, Counter
from pathlib import Path

import numpy as np
from scipy import stats as sp_stats


def load_csv(path: str) -> list[dict]:
    rows = []
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            d = {}
            for k, v in row.items():
                try:
                    d[k] = float(v)
                except ValueError:
                    d[k] = v
            d["N"] = int(d["N"])
            d["rep"] = int(d["rep"])
            if "layer_count" in d:
                d["layer_count"] = int(d["layer_count"])
            rows.append(d)
    return rows


def main():
    csv_path = "outputs_d_recovery/prediction_c_f7_large_n.csv"
    print(f"Loading: {csv_path}")
    data = load_csv(csv_path)
    n_values = sorted(set(r["N"] for r in data))
    families = sorted(set(r["family"] for r in data))
    print(f"Loaded {len(data)} rows, N={n_values}, families={families}\n")

    report = []
    report.append("# Prediction C — Lor5D Weakness Diagnosis\n")
    report.append(f"**Data**: {len(data)} samples from prediction_c_f7_large_n.csv\n")

    # ── T1: Σ_hist distribution ──
    report.append("## T1: Σ_hist Distribution per (family, N)\n")
    report.append("| family | N | n | unique_vals | min | max | range | std | CV |")
    report.append("|--------|---|---|-------------|-----|-----|-------|-----|-----|")

    by_nf = defaultdict(list)
    for r in data:
        by_nf[(r["N"], r["family"])].append(r)

    for N in n_values:
        for fam in families:
            vals = by_nf.get((N, fam), [])
            if not vals:
                continue
            sh = np.array([r["sigma_hist"] for r in vals])
            unique = len(set(sh))
            mn, mx = sh.min(), sh.max()
            rng = mx - mn
            std = sh.std()
            cv = std / abs(sh.mean()) if abs(sh.mean()) > 1e-12 else float("inf")
            report.append(f"| {fam} | {N} | {len(vals)} | {unique} | {mn:.4f} | {mx:.4f} | "
                          f"{rng:.4f} | {std:.4f} | {cv:.2f} |")

    # ── T2: layer_count distribution ──
    has_lc = "layer_count" in data[0]
    if has_lc:
        report.append("\n## T2: layer_count Distribution\n")
        report.append("| family | N | unique_lc | distribution | entropy_bits |")
        report.append("|--------|---|-----------|-------------|--------------|")

        for N in n_values:
            for fam in families:
                vals = by_nf.get((N, fam), [])
                if not vals or "layer_count" not in vals[0]:
                    continue
                lc = [r["layer_count"] for r in vals]
                counter = Counter(lc)
                unique = len(counter)
                dist_str = ", ".join(f"{k}:{v}" for k, v in sorted(counter.items()))
                # entropy
                total = len(lc)
                probs = np.array([v / total for v in counter.values()])
                entropy = -np.sum(probs * np.log2(probs + 1e-30))
                report.append(f"| {fam} | {N} | {unique} | {dist_str} | {entropy:.2f} |")

    # ── T3: logH conditional on layer_count ──
    if has_lc:
        report.append("\n## T3: Mean logH by layer_count (Lor5D vs Lor2D control)\n")
        report.append("Does more layers → less logH within Lor5D?\n")

        for fam in ["Lor2D", "Lor5D"]:
            report.append(f"\n### {fam}\n")
            report.append("| N | layer_count | n | mean_logH | std_logH |")
            report.append("|---|------------|---|-----------|----------|")

            for N in n_values:
                vals = by_nf.get((N, fam), [])
                if not vals:
                    continue
                by_lc = defaultdict(list)
                for r in vals:
                    by_lc[r["layer_count"]].append(r["log_H"])

                for lc in sorted(by_lc.keys()):
                    lh = np.array(by_lc[lc])
                    report.append(f"| {N} | {lc} | {len(lh)} | {lh.mean():.2f} | {lh.std():.2f} |")

    # ── T4: Component decomposition for reversal cells ──
    report.append("\n## T4: Component Decomposition at Reversal Cells\n")
    report.append("Why does Lor5D@N=52/72 reverse?\n")

    for N in [52, 72]:
        vals = by_nf.get((N, "Lor5D"), [])
        if not vals:
            continue

        report.append(f"\n### Lor5D @ N={N} (n={len(vals)})\n")

        sh = np.array([r["sigma_hist"] for r in vals])
        lh = np.array([r["log_H"] for r in vals])
        R = np.array([r["R"] for r in vals])
        wall = np.array([r["wall"] for r in vals])
        xi = np.array([r["xi_dim"] for r in vals])
        pi = np.array([r["pi_geo"] for r in vals])
        f7 = np.array([r["F7"] for r in vals])

        report.append("**Distributions:**")
        report.append(f"- Σ_hist: {sh.min():.4f}–{sh.max():.4f}, std={sh.std():.4f}, "
                      f"unique={len(set(sh))}")
        report.append(f"- logH: {lh.min():.1f}–{lh.max():.1f}, std={lh.std():.1f}")
        report.append(f"- R: {R.min():.4f}–{R.max():.4f}, std={R.std():.4f}")
        report.append(f"- wall: {wall.min():.2f}–{wall.max():.2f}, std={wall.std():.2f}")
        report.append(f"- xi_dim: {xi.min():.3f}–{xi.max():.3f}, std={xi.std():.3f}")

        if sh.std() > 1e-12 and lh.std() > 1e-12:
            rho_sh_lh, p_sh_lh = sp_stats.spearmanr(sh, lh)
            rho_sh_R, _ = sp_stats.spearmanr(sh, R)
            rho_sh_wall, _ = sp_stats.spearmanr(sh, wall)
            rho_sh_xi, _ = sp_stats.spearmanr(sh, xi)
            rho_R_lh, _ = sp_stats.spearmanr(R, lh)

            report.append(f"\n**Correlations with Σ_hist:**")
            report.append(f"- ρ(Σ_hist, logH) = {rho_sh_lh:+.3f} (p={p_sh_lh:.4f}) {'❌ REVERSED' if rho_sh_lh > 0 else '✅'}")
            report.append(f"- ρ(Σ_hist, R) = {rho_sh_R:+.3f}")
            report.append(f"- ρ(Σ_hist, wall) = {rho_sh_wall:+.3f}")
            report.append(f"- ρ(Σ_hist, xi_dim) = {rho_sh_xi:+.3f}")
            report.append(f"- ρ(R, logH) = {rho_R_lh:+.3f}")

            # Check if the reversal comes from R confounding
            # Partial correlation: ρ(Σ_hist, logH | R)
            # Compute rank-based partial
            from scipy.stats import rankdata
            r_sh = rankdata(sh)
            r_lh = rankdata(lh)
            r_R = rankdata(R)

            # Residualize
            from numpy.polynomial.polynomial import polyfit
            c_sh = np.polyfit(r_R, r_sh, 1)
            c_lh = np.polyfit(r_R, r_lh, 1)
            res_sh = r_sh - np.polyval(c_sh, r_R)
            res_lh = r_lh - np.polyval(c_lh, r_R)
            rho_partial, _ = sp_stats.spearmanr(res_sh, res_lh)
            report.append(f"- ρ_partial(Σ_hist, logH | R) = {rho_partial:+.3f}")
        else:
            report.append("\n**Σ_hist is degenerate (std≈0) — cannot compute correlations**")

    # ── T5: Cross-family comparison at same N ──
    report.append("\n## T5: Family Comparison of Σ_hist Variance\n")
    report.append("Ratio of Lor5D variance to Lor2D variance:\n")
    report.append("| N | std(Σ_hist) Lor2D | std(Σ_hist) Lor5D | ratio 5D/2D |")
    report.append("|---|-------------------|-------------------|-------------|")

    for N in n_values:
        v2d = by_nf.get((N, "Lor2D"), [])
        v5d = by_nf.get((N, "Lor5D"), [])
        if not v2d or not v5d:
            continue
        std2d = np.std([r["sigma_hist"] for r in v2d])
        std5d = np.std([r["sigma_hist"] for r in v5d])
        ratio = std5d / std2d if std2d > 1e-12 else float("inf")
        report.append(f"| {N} | {std2d:.4f} | {std5d:.4f} | {ratio:.3f} |")

    # ── Verdict ──
    report.append("\n## Diagnosis Verdict\n")

    out = "\n".join(report)
    outdir = Path("outputs_d_recovery")
    md_path = outdir / "prediction_c_lor5d_diagnosis.md"
    md_path.write_text(out, encoding="utf-8")
    print(f"Saved: {md_path}")
    print("\n" + out)


if __name__ == "__main__":
    main()
