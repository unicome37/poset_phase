"""Enhanced E2 v3: Direction-agnostic structural complexity metrics.

The original E2 (layer count ratio) was weak because:
  1. Layer count is integer and nearly symmetric (fwd ≈ bwd).
  2. Mean element depth has backward bias (adding to bottom shifts all layers up).

This version uses normalized, direction-AGNOSTIC metrics:
  1. Height ratio h = n_layers / N_total — structural depth per element
  2. Width entropy S_w = -Σ p_ℓ log p_ℓ where p_ℓ = |layer_ℓ|/N — layer complexity
  3. Order fraction R = comparable_pairs / all_pairs — causal density
  4. Σ_hist = historical sedimentation — the functional's own metric

The E2 prediction is reformulated as:
  "Forward augmentation PRESERVES structural complexity better than backward."
  Equivalently: Δh_fwd > Δh_bwd, ΔS_w_fwd > ΔS_w_bwd, etc.

Physical motivation: CST sprinkling into the future should extend the causal
structure coherently (preserving depth-per-element), while sprinkling into the
past should DILUTE the existing structure (many shallow predecessors compress
the layer structure).
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, asdict, fields
from pathlib import Path
import csv

import numpy as np
from scipy import stats as sp_stats

from generators import (
    Poset,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_kr_like,
    generate_transitive_percolation,
)
from observables import layer_profile
from prediction_e_experiment import (
    augment_forward,
    augment_backward,
    _layer_assignment,
)
from unified_functional import compute_sigma_hist

OUT_DIR = Path("outputs_unified_functional")
OUT_DIR.mkdir(parents=True, exist_ok=True)

FAMILY_GENS = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "KR_like": generate_kr_like,
    "TransPerc": lambda n, seed=None: generate_transitive_percolation(n, p=0.08, seed=seed),
}


def compute_metrics(poset: Poset) -> dict:
    """Compute direction-agnostic structural complexity metrics."""
    N = poset.n
    lp = layer_profile(poset)
    n_layers = len(lp)

    # 1. Height ratio: structural depth per element
    height_ratio = n_layers / max(N, 1)

    # 2. Width entropy: Shannon entropy of normalized layer profile
    p = lp / max(lp.sum(), 1)
    p = p[p > 0]
    width_entropy = -float(np.sum(p * np.log(p)))

    # 3. Order fraction R
    c = poset.closure
    comp = int(c.sum()) - N  # exclude diagonal
    total_pairs = N * (N - 1)
    R = comp / max(total_pairs, 1)

    # 4. Sigma_hist
    sigma_hist = compute_sigma_hist(poset)

    # 5. Max width / N (layer concentration)
    max_width_frac = float(lp.max()) / max(N, 1)

    return {
        "N": N,
        "n_layers": n_layers,
        "height_ratio": height_ratio,
        "width_entropy": width_entropy,
        "R": R,
        "sigma_hist": sigma_hist,
        "max_width_frac": max_width_frac,
    }


def main():
    FAMILIES = ["Lor2D", "Lor3D", "Lor4D", "KR_like", "TransPerc"]
    N_VALUES = [16, 24, 36]
    T_VALUES = [4, 8, 16]
    REPS = 30
    BASE_SEED = 20260323

    total = len(FAMILIES) * len(N_VALUES) * len(T_VALUES) * REPS
    print("=" * 70)
    print("E2 v3: Direction-Agnostic Structural Complexity")
    print("=" * 70)
    print(f"Total experiments: {total}")

    t0 = time.time()
    rows = []
    count = 0

    for fam in FAMILIES:
        gen = FAMILY_GENS[fam]
        for n0 in N_VALUES:
            for T in T_VALUES:
                for rep in range(REPS):
                    count += 1
                    s = BASE_SEED + rep * 1000 + n0 * 100 + T * 10 + hash(fam) % 10000

                    poset = gen(n0, seed=s)
                    base_m = compute_metrics(poset)

                    # Forward T steps
                    p_fwd = poset
                    for t in range(T):
                        p_fwd = augment_forward(p_fwd, k=1, seed=s + t * 7)
                    fwd_m = compute_metrics(p_fwd)

                    # Backward T steps
                    p_bwd = poset
                    for t in range(T):
                        p_bwd = augment_backward(p_bwd, k=1, seed=s + t * 7)
                    bwd_m = compute_metrics(p_bwd)

                    row = {
                        "family": fam, "N0": n0, "T": T, "rep": rep, "seed": s,
                    }
                    for key in ["height_ratio", "width_entropy", "R", "sigma_hist", "max_width_frac"]:
                        row[f"base_{key}"] = base_m[key]
                        row[f"fwd_{key}"] = fwd_m[key]
                        row[f"bwd_{key}"] = bwd_m[key]
                        row[f"delta_fwd_{key}"] = fwd_m[key] - base_m[key]
                        row[f"delta_bwd_{key}"] = bwd_m[key] - base_m[key]
                        row[f"diff_{key}"] = (fwd_m[key] - base_m[key]) - (bwd_m[key] - base_m[key])

                    rows.append(row)

                    if count % 200 == 0:
                        dh = row["diff_height_ratio"]
                        ds = row["diff_sigma_hist"]
                        dw = row["diff_width_entropy"]
                        print(f"  [{count}/{total}] {fam} N={n0} T={T}: "
                              f"Δh={dh:+.4f} ΔS_w={dw:+.4f} ΔΣ_hist={ds:+.4f}")

    elapsed = time.time() - t0
    print(f"\nTotal time: {elapsed:.1f}s")

    # Save CSV
    csv_path = OUT_DIR / "prediction_e2_v3.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"Saved: {csv_path}")

    # ══════════════════════════════════════════════════════════════════
    # Analysis
    # ══════════════════════════════════════════════════════════════════
    METRICS = [
        ("height_ratio", "Height Ratio (n_layers / N)", "fwd deepens more"),
        ("width_entropy", "Width Entropy S_w", "fwd produces richer layer structure"),
        ("R", "Order Fraction R", "fwd preserves causal density better"),
        ("sigma_hist", "Σ_hist (Historical Sedimentation)", "fwd produces deeper sedimentation"),
        ("max_width_frac", "Max Width Fraction", "lower = deeper structure"),
    ]

    lines = ["# E2 v3: Direction-Agnostic Structural Complexity\n"]
    lines.append(f"**Date:** 2026-03-23\n")
    lines.append(f"**Design:** {len(FAMILIES)} families × {len(N_VALUES)} N × {len(T_VALUES)} T × {REPS} reps = {total}\n\n")

    lines.append("## Reformulated Prediction\n")
    lines.append("E2 (original): \"Forward augmentation deepens layer structure faster\" — WEAK (layer count is symmetric)\n")
    lines.append("E2 (v3): \"Forward augmentation preserves/enhances structural complexity better than backward\"\n")
    lines.append("- Height ratio h = n_layers/N: forward should maintain depth-per-element\n")
    lines.append("- Width entropy S_w: forward should produce richer, less concentrated layer profiles\n")
    lines.append("- Order fraction R: forward should preserve causal density\n")
    lines.append("- Σ_hist: forward should increase historical sedimentation\n\n")

    for metric_key, metric_name, prediction in METRICS:
        lines.append(f"## {metric_name}\n")
        lines.append(f"**Prediction:** {prediction}\n")

        # Per-family table
        lines.append("### Per-family (all T pooled)\n")
        lines.append("| Family | mean Δ_fwd | mean Δ_bwd | mean diff | fwd>bwd % | Wilcoxon p | verdict |")
        lines.append("|--------|-----------|-----------|-----------|-----------|------------|---------|")

        print(f"\n{'='*70}")
        print(f"METRIC: {metric_name}")
        print(f"{'='*70}")

        for fam in FAMILIES:
            sub = [r for r in rows if r["family"] == fam]
            d_f = np.array([r[f"delta_fwd_{metric_key}"] for r in sub])
            d_b = np.array([r[f"delta_bwd_{metric_key}"] for r in sub])
            diff = d_f - d_b
            pct = 100.0 * np.mean(diff > 0)
            nz = diff[np.abs(diff) > 1e-10]
            if len(nz) >= 10:
                _, pval = sp_stats.wilcoxon(nz)
            else:
                pval = 1.0
            verdict = "✅" if pct > 50 and pval < 0.05 else ("⚠️" if pct > 50 else "❌")
            lines.append(f"| {fam} | {d_f.mean():+.4f} | {d_b.mean():+.4f} | {diff.mean():+.4f} | {pct:.1f}% | {pval:.3e} | {verdict} |")
            print(f"  {fam:>10s}: diff={diff.mean():+.4f}, fwd>bwd={pct:.1f}%, p={pval:.3e} {verdict}")

        # Per-T
        lines.append(f"\n### Per-T (all families pooled)\n")
        lines.append("| T | mean diff | fwd>bwd % | Wilcoxon p |")
        lines.append("|---|-----------|-----------|------------|")
        for T in T_VALUES:
            sub = [r for r in rows if r["T"] == T]
            diff = np.array([r[f"diff_{metric_key}"] for r in sub])
            pct = 100.0 * np.mean(diff > 0)
            nz = diff[np.abs(diff) > 1e-10]
            if len(nz) >= 10:
                _, pval = sp_stats.wilcoxon(nz)
            else:
                pval = 1.0
            lines.append(f"| {T} | {diff.mean():+.4f} | {pct:.1f}% | {pval:.3e} |")

        # Family × T
        lines.append(f"\n### Family × T (fwd>bwd rate)\n")
        hdr = "| Family |"
        for T in T_VALUES:
            hdr += f" T={T} |"
        lines.append(hdr)
        lines.append("|--------|" + "---:|" * len(T_VALUES))
        for fam in FAMILIES:
            row_str = f"| {fam} |"
            for T in T_VALUES:
                sub = [r for r in rows if r["family"] == fam and r["T"] == T]
                diff = np.array([r[f"diff_{metric_key}"] for r in sub])
                pct = 100.0 * np.mean(diff > 0)
                row_str += f" {pct:.0f}% |"
            lines.append(row_str)

        lines.append("")

    # ── Grand Summary ──
    lines.append("## Grand Summary\n")
    lines.append("| Metric | Lorentzian fwd>bwd % | KR fwd>bwd % | TransPerc fwd>bwd % | Signal? |")
    lines.append("|--------|---------------------|-------------|--------------------|---------|\n")

    lor_fams = ["Lor2D", "Lor3D", "Lor4D"]
    for metric_key, metric_name, _ in METRICS:
        lor_sub = [r for r in rows if r["family"] in lor_fams]
        kr_sub = [r for r in rows if r["family"] == "KR_like"]
        tp_sub = [r for r in rows if r["family"] == "TransPerc"]

        pct_lor = 100.0 * np.mean([r[f"diff_{metric_key}"] > 0 for r in lor_sub])
        pct_kr = 100.0 * np.mean([r[f"diff_{metric_key}"] > 0 for r in kr_sub])
        pct_tp = 100.0 * np.mean([r[f"diff_{metric_key}"] > 0 for r in tp_sub])

        signal = "✅" if pct_lor > 55 else ("⚠️" if pct_lor > 50 else "❌")
        lines.append(f"| {metric_name} | {pct_lor:.1f}% | {pct_kr:.1f}% | {pct_tp:.1f}% | {signal} |")

    lines.append("\n## Interpretation\n")
    lines.append("If forward augmentation consistently produces HIGHER structural complexity\n")
    lines.append("(height ratio, width entropy, Σ_hist) than backward, this confirms E2:\n")
    lines.append("the causal future direction is structurally richer than the past direction.\n")
    lines.append("This is complementary to E1 (entropy asymmetry) and strengthens the\n")
    lines.append("\"time arrow from structural asymmetry\" thesis.\n")

    md_path = OUT_DIR / "prediction_e2_v3.md"
    md_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nReport: {md_path}")
    print("=" * 70)
    print("DONE")


if __name__ == "__main__":
    main()
