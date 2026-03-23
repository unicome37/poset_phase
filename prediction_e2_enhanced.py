"""Enhanced E2: Causal Depth Directional Growth — Stronger verification.

Original E2 weaknesses:
  1. Layer count is integer → depth_gain ∈ {0,1,2}, R_depth extremely noisy
  2. Small N (12–20) → base depth only 4–6 layers
  3. Small T (2–8) → insufficient augmentation to accumulate signal

Enhancements:
  1. Larger N (16, 24, 36, 48) for deeper base structures
  2. Larger T (4, 8, 16, 32) for stronger directional accumulation
  3. Multiple depth metrics:
     - n_layers: antichain layer count (original)
     - longest_chain: length of longest chain (= height of poset)
     - mean_depth: mean per-element depth (continuous-valued)
     - width_depth_ratio: max_layer_width / n_layers
  4. Paired analysis: Wilcoxon signed-rank test on Δdepth_fwd - Δdepth_bwd
  5. More reps (30) for better statistics
"""
from __future__ import annotations

import csv
import math
import time
from dataclasses import dataclass, asdict
from pathlib import Path

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

OUT_DIR = Path("outputs_unified_functional")
OUT_DIR.mkdir(parents=True, exist_ok=True)

FAMILY_GENS = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "KR_like": generate_kr_like,
    "TransPerc": lambda n, seed=None: generate_transitive_percolation(n, p=0.08, seed=seed),
}


def longest_chain_length(poset: Poset) -> int:
    """Compute the longest chain (= height) of the poset via DP on layers."""
    layers = _layer_assignment(poset.closure)
    return int(layers.max()) + 1 if len(layers) > 0 else 0


def mean_element_depth(poset: Poset) -> float:
    """Mean per-element layer depth (continuous-valued)."""
    layers = _layer_assignment(poset.closure)
    return float(layers.mean()) if len(layers) > 0 else 0.0


def max_layer_width(poset: Poset) -> int:
    """Width of the widest antichain layer."""
    lp = layer_profile(poset)
    return int(lp.max()) if len(lp) > 0 else 0


def compute_depth_metrics(poset: Poset) -> dict:
    """Compute all depth metrics for a poset."""
    lp = layer_profile(poset)
    layers = _layer_assignment(poset.closure)
    n_layers = len(lp)
    lc = int(layers.max()) + 1 if len(layers) > 0 else 0
    md = float(layers.mean()) if len(layers) > 0 else 0.0
    mw = int(lp.max()) if len(lp) > 0 else 0
    wdr = mw / max(n_layers, 1)
    return {
        "n_layers": n_layers,
        "longest_chain": lc,
        "mean_depth": md,
        "max_width": mw,
        "width_depth_ratio": wdr,
    }


@dataclass
class E2EnhRow:
    family: str
    N0: int
    T: int
    rep: int
    seed: int
    # Base metrics
    base_n_layers: int
    base_longest_chain: int
    base_mean_depth: float
    base_max_width: int
    # Forward metrics
    fwd_n_layers: int
    fwd_longest_chain: int
    fwd_mean_depth: float
    fwd_max_width: int
    # Backward metrics
    bwd_n_layers: int
    bwd_longest_chain: int
    bwd_mean_depth: float
    bwd_max_width: int
    # Deltas
    delta_layers_fwd: int
    delta_layers_bwd: int
    delta_chain_fwd: int
    delta_chain_bwd: int
    delta_mean_depth_fwd: float
    delta_mean_depth_bwd: float
    # Paired difference (fwd - bwd)
    diff_layers: int
    diff_chain: int
    diff_mean_depth: float


def main():
    FAMILIES = ["Lor2D", "Lor3D", "Lor4D", "KR_like", "TransPerc"]
    N_VALUES = [16, 24, 36, 48]
    T_VALUES = [4, 8, 16, 32]
    REPS = 30
    BASE_SEED = 20260323

    print("=" * 70)
    print("Enhanced E2: Causal Depth Directional Growth")
    print("=" * 70)
    print(f"Families: {FAMILIES}")
    print(f"N values: {N_VALUES}")
    print(f"T values: {T_VALUES}")
    print(f"Reps: {REPS}")

    total = len(FAMILIES) * len(N_VALUES) * len(T_VALUES) * REPS
    print(f"Total experiments: {total}")

    t0 = time.time()
    rows: list[E2EnhRow] = []
    count = 0

    for fam in FAMILIES:
        gen = FAMILY_GENS[fam]
        for n0 in N_VALUES:
            for T in T_VALUES:
                for rep in range(REPS):
                    count += 1
                    s = BASE_SEED + rep * 1000 + n0 * 100 + T * 10 + hash(fam) % 10000

                    poset = gen(n0, seed=s)
                    base_m = compute_depth_metrics(poset)

                    # Forward T steps
                    p_fwd = poset
                    for t in range(T):
                        p_fwd = augment_forward(p_fwd, k=1, seed=s + t * 7)
                    fwd_m = compute_depth_metrics(p_fwd)

                    # Backward T steps
                    p_bwd = poset
                    for t in range(T):
                        p_bwd = augment_backward(p_bwd, k=1, seed=s + t * 7)
                    bwd_m = compute_depth_metrics(p_bwd)

                    dl_f = fwd_m["n_layers"] - base_m["n_layers"]
                    dl_b = bwd_m["n_layers"] - base_m["n_layers"]
                    dc_f = fwd_m["longest_chain"] - base_m["longest_chain"]
                    dc_b = bwd_m["longest_chain"] - base_m["longest_chain"]
                    dm_f = fwd_m["mean_depth"] - base_m["mean_depth"]
                    dm_b = bwd_m["mean_depth"] - base_m["mean_depth"]

                    rows.append(E2EnhRow(
                        family=fam, N0=n0, T=T, rep=rep, seed=s,
                        base_n_layers=base_m["n_layers"],
                        base_longest_chain=base_m["longest_chain"],
                        base_mean_depth=base_m["mean_depth"],
                        base_max_width=base_m["max_width"],
                        fwd_n_layers=fwd_m["n_layers"],
                        fwd_longest_chain=fwd_m["longest_chain"],
                        fwd_mean_depth=fwd_m["mean_depth"],
                        fwd_max_width=fwd_m["max_width"],
                        bwd_n_layers=bwd_m["n_layers"],
                        bwd_longest_chain=bwd_m["longest_chain"],
                        bwd_mean_depth=bwd_m["mean_depth"],
                        bwd_max_width=bwd_m["max_width"],
                        delta_layers_fwd=dl_f, delta_layers_bwd=dl_b,
                        delta_chain_fwd=dc_f, delta_chain_bwd=dc_b,
                        delta_mean_depth_fwd=dm_f, delta_mean_depth_bwd=dm_b,
                        diff_layers=dl_f - dl_b,
                        diff_chain=dc_f - dc_b,
                        diff_mean_depth=dm_f - dm_b,
                    ))

                    if count % 200 == 0:
                        print(f"  [{count}/{total}] {fam} N={n0} T={T}: "
                              f"Δlayers fwd={dl_f} bwd={dl_b}, "
                              f"Δmean_depth fwd={dm_f:+.2f} bwd={dm_b:+.2f}")

    elapsed = time.time() - t0
    print(f"\nTotal time: {elapsed:.1f}s")

    # Save CSV
    csv_path = OUT_DIR / "prediction_e2_enhanced.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()))
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))
    print(f"Saved: {csv_path}")

    # ═══════════════════════════════════════════════════════════════════
    # Analysis
    # ═══════════════════════════════════════════════════════════════════
    data = [asdict(r) for r in rows]

    lines = ["# Enhanced E2: Causal Depth Directional Growth\n"]
    lines.append(f"**Date:** 2026-03-23\n")
    lines.append(f"**Design:** {len(FAMILIES)} families × {len(N_VALUES)} N × {len(T_VALUES)} T × {REPS} reps = {total} experiments\n")
    lines.append(f"**Families:** {', '.join(FAMILIES)}\n")
    lines.append(f"**N values:** {N_VALUES}\n")
    lines.append(f"**T values:** {T_VALUES}\n\n")

    # ── Metric 1: Layer count (original) ──
    lines.append("## Metric 1: Layer Count (original E2 metric)\n")
    lines.append("### Per-family (all T pooled)\n")
    lines.append("| Family | mean Δ_fwd | mean Δ_bwd | mean diff | fwd>bwd % | Wilcoxon p | n |")
    lines.append("|--------|-----------|-----------|-----------|-----------|------------|---|")

    print(f"\n{'='*70}")
    print("METRIC 1: LAYER COUNT")
    print(f"{'='*70}")

    for fam in FAMILIES:
        sub = [d for d in data if d["family"] == fam]
        dl_f = np.array([d["delta_layers_fwd"] for d in sub])
        dl_b = np.array([d["delta_layers_bwd"] for d in sub])
        diff = dl_f - dl_b
        pct = 100.0 * np.mean(diff > 0)
        nonzero = diff[diff != 0]
        if len(nonzero) > 0:
            stat, pval = sp_stats.wilcoxon(nonzero)
        else:
            stat, pval = 0, 1.0
        lines.append(f"| {fam} | {dl_f.mean():.2f} | {dl_b.mean():.2f} | {diff.mean():+.3f} | {pct:.1f}% | {pval:.3e} | {len(sub)} |")
        print(f"  {fam:>10s}: mean diff={diff.mean():+.3f}, fwd>bwd={pct:.1f}%, p={pval:.3e}")

    # ── Metric 2: Mean element depth (continuous!) ──
    lines.append("\n## Metric 2: Mean Element Depth (continuous-valued)\n")
    lines.append("### Per-family (all T pooled)\n")
    lines.append("| Family | mean Δ_fwd | mean Δ_bwd | mean diff | fwd>bwd % | Wilcoxon p | n |")
    lines.append("|--------|-----------|-----------|-----------|-----------|------------|---|")

    print(f"\n{'='*70}")
    print("METRIC 2: MEAN ELEMENT DEPTH")
    print(f"{'='*70}")

    for fam in FAMILIES:
        sub = [d for d in data if d["family"] == fam]
        dm_f = np.array([d["delta_mean_depth_fwd"] for d in sub])
        dm_b = np.array([d["delta_mean_depth_bwd"] for d in sub])
        diff = dm_f - dm_b
        pct = 100.0 * np.mean(diff > 0)
        nonzero = diff[np.abs(diff) > 1e-10]
        if len(nonzero) > 0:
            stat, pval = sp_stats.wilcoxon(nonzero)
        else:
            stat, pval = 0, 1.0
        lines.append(f"| {fam} | {dm_f.mean():.3f} | {dm_b.mean():.3f} | {diff.mean():+.4f} | {pct:.1f}% | {pval:.3e} | {len(sub)} |")
        print(f"  {fam:>10s}: mean diff={diff.mean():+.4f}, fwd>bwd={pct:.1f}%, p={pval:.3e}")

    # ── Per-T breakdown (mean_depth, pooled families) ──
    lines.append("\n### Per-T breakdown (mean element depth, all families pooled)\n")
    lines.append("| T | mean Δ_fwd | mean Δ_bwd | mean diff | fwd>bwd % | Wilcoxon p | n |")
    lines.append("|---|-----------|-----------|-----------|-----------|------------|---|")

    print(f"\nPer-T (mean_depth):")
    for T in T_VALUES:
        sub = [d for d in data if d["T"] == T]
        dm_f = np.array([d["delta_mean_depth_fwd"] for d in sub])
        dm_b = np.array([d["delta_mean_depth_bwd"] for d in sub])
        diff = dm_f - dm_b
        pct = 100.0 * np.mean(diff > 0)
        nonzero = diff[np.abs(diff) > 1e-10]
        if len(nonzero) > 0:
            stat, pval = sp_stats.wilcoxon(nonzero)
        else:
            stat, pval = 0, 1.0
        lines.append(f"| {T} | {dm_f.mean():.3f} | {dm_b.mean():.3f} | {diff.mean():+.4f} | {pct:.1f}% | {pval:.3e} | {len(sub)} |")
        print(f"  T={T:>3d}: diff={diff.mean():+.4f}, fwd>bwd={pct:.1f}%, p={pval:.3e}")

    # ── Per-N breakdown (mean_depth) ──
    lines.append("\n### Per-N breakdown (mean element depth, all families pooled)\n")
    lines.append("| N | mean diff | fwd>bwd % | Wilcoxon p | n |")
    lines.append("|---|-----------|-----------|------------|---|")

    print(f"\nPer-N (mean_depth):")
    for N in N_VALUES:
        sub = [d for d in data if d["N0"] == N]
        dm_f = np.array([d["delta_mean_depth_fwd"] for d in sub])
        dm_b = np.array([d["delta_mean_depth_bwd"] for d in sub])
        diff = dm_f - dm_b
        pct = 100.0 * np.mean(diff > 0)
        nonzero = diff[np.abs(diff) > 1e-10]
        if len(nonzero) > 0:
            stat, pval = sp_stats.wilcoxon(nonzero)
        else:
            stat, pval = 0, 1.0
        lines.append(f"| {N} | {diff.mean():+.4f} | {pct:.1f}% | {pval:.3e} | {len(sub)} |")
        print(f"  N={N:>3d}: diff={diff.mean():+.4f}, fwd>bwd={pct:.1f}%, p={pval:.3e}")

    # ── Family × T heatmap (mean_depth diff) ──
    lines.append("\n### Family × T (mean diff in mean_depth)\n")
    hdr = "| Family |"
    for T in T_VALUES:
        hdr += f" T={T} |"
    lines.append(hdr)
    lines.append("|--------|" + "---:|" * len(T_VALUES))

    print(f"\nFamily × T (mean_depth diff):")
    for fam in FAMILIES:
        row = f"| {fam} |"
        for T in T_VALUES:
            sub = [d for d in data if d["family"] == fam and d["T"] == T]
            diff = np.array([d["delta_mean_depth_fwd"] - d["delta_mean_depth_bwd"] for d in sub])
            row += f" {diff.mean():+.3f} |"
        lines.append(row)
        print(f"  {fam}")

    # ── Family × T heatmap (fwd>bwd %) ──
    lines.append("\n### Family × T (fwd > bwd rate for mean_depth)\n")
    hdr = "| Family |"
    for T in T_VALUES:
        hdr += f" T={T} |"
    lines.append(hdr)
    lines.append("|--------|" + "---:|" * len(T_VALUES))

    for fam in FAMILIES:
        row = f"| {fam} |"
        for T in T_VALUES:
            sub = [d for d in data if d["family"] == fam and d["T"] == T]
            diff = np.array([d["delta_mean_depth_fwd"] - d["delta_mean_depth_bwd"] for d in sub])
            pct = 100.0 * np.mean(diff > 0)
            row += f" {pct:.0f}% |"
        lines.append(row)

    # ── Global summary ──
    lines.append("\n## Global Summary\n")

    # All Lorentzian pooled
    lor_fams = ["Lor2D", "Lor3D", "Lor4D"]
    lor_sub = [d for d in data if d["family"] in lor_fams]
    diff_lor = np.array([d["delta_mean_depth_fwd"] - d["delta_mean_depth_bwd"] for d in lor_sub])
    pct_lor = 100.0 * np.mean(diff_lor > 0)
    nz_lor = diff_lor[np.abs(diff_lor) > 1e-10]
    _, p_lor = sp_stats.wilcoxon(nz_lor) if len(nz_lor) > 0 else (0, 1.0)

    # KR_like only
    kr_sub = [d for d in data if d["family"] == "KR_like"]
    diff_kr = np.array([d["delta_mean_depth_fwd"] - d["delta_mean_depth_bwd"] for d in kr_sub])
    pct_kr = 100.0 * np.mean(diff_kr > 0)
    nz_kr = diff_kr[np.abs(diff_kr) > 1e-10]
    _, p_kr = sp_stats.wilcoxon(nz_kr) if len(nz_kr) > 0 else (0, 1.0)

    # TransPerc only
    tp_sub = [d for d in data if d["family"] == "TransPerc"]
    diff_tp = np.array([d["delta_mean_depth_fwd"] - d["delta_mean_depth_bwd"] for d in tp_sub])
    pct_tp = 100.0 * np.mean(diff_tp > 0)
    nz_tp = diff_tp[np.abs(diff_tp) > 1e-10]
    _, p_tp = sp_stats.wilcoxon(nz_tp) if len(nz_tp) > 0 else (0, 1.0)

    lines.append("| Group | mean diff | fwd>bwd % | Wilcoxon p | n |")
    lines.append("|-------|-----------|-----------|------------|---|")
    lines.append(f"| Lorentzian (2D/3D/4D) | {diff_lor.mean():+.4f} | {pct_lor:.1f}% | {p_lor:.3e} | {len(lor_sub)} |")
    lines.append(f"| KR_like (control) | {diff_kr.mean():+.4f} | {pct_kr:.1f}% | {p_kr:.3e} | {len(kr_sub)} |")
    lines.append(f"| TransPerc | {diff_tp.mean():+.4f} | {pct_tp:.1f}% | {p_tp:.3e} | {len(tp_sub)} |")

    print(f"\nGlobal summary (mean_depth):")
    print(f"  Lorentzian: diff={diff_lor.mean():+.4f}, fwd>bwd={pct_lor:.1f}%, p={p_lor:.3e}")
    print(f"  KR_like:    diff={diff_kr.mean():+.4f}, fwd>bwd={pct_kr:.1f}%, p={p_kr:.3e}")
    print(f"  TransPerc:  diff={diff_tp.mean():+.4f}, fwd>bwd={pct_tp:.1f}%, p={p_tp:.3e}")

    # ── Scaling: does signal grow with T? ──
    lines.append("\n## Signal Scaling with T\n")
    lines.append("Does the directional asymmetry grow with augmentation steps?\n")
    lines.append("| Family | Spearman(T, mean_diff) | p-value |")
    lines.append("|--------|----------------------|---------|")

    print(f"\nScaling with T:")
    for fam in FAMILIES:
        per_T_diffs = []
        for T in T_VALUES:
            sub = [d for d in data if d["family"] == fam and d["T"] == T]
            diff = np.mean([d["delta_mean_depth_fwd"] - d["delta_mean_depth_bwd"] for d in sub])
            per_T_diffs.append((T, diff))
        Ts = [x[0] for x in per_T_diffs]
        diffs = [x[1] for x in per_T_diffs]
        rho, pval = sp_stats.spearmanr(Ts, diffs)
        lines.append(f"| {fam} | {rho:+.3f} | {pval:.3e} |")
        print(f"  {fam:>10s}: rho(T, diff)={rho:+.3f}, p={pval:.3e}")

    # ── Conclusion ──
    lines.append("\n## Conclusion\n")
    lines.append("(Filled based on results above)\n")

    # Write report
    md_path = OUT_DIR / "prediction_e2_enhanced.md"
    md_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nReport: {md_path}")
    print("=" * 70)
    print("DONE")


if __name__ == "__main__":
    main()
