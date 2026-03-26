"""
LSD-Well: Lorentzian Structural Discriminator with WELL potentials
====================================================================
Key insight from LSD1-4: monotone rewards fail because
  - C₁/C₀: layered >> Lor4D >> KR_2layer
  - width_ratio: KR_2layer >> Lor4D >> layered
  → No monotone function of either can select Lor4D

Solution: WELL potentials — penalize deviation from Lor4D's characteristic values.
This is the same architecture as F10's d_eff-well, but now using the
features that ACTUALLY separate families (not d_eff/R which conflate KR with Lor).

Candidates:
  LSD-W1: α·N·(C₁/C₀ - c*)² + β·N·(width - w*)²
           → pure interval-shape + width well, O(N) penalty
  LSD-W2: LSD-W1 + wall(R)
           → adds sparse-structure penalty
  LSD-W3: α·N·(C₁/C₀ - c*)² + β·N·(width - w*)² + γ·N·(d_eff - d*)²
           → triple well: shape + width + dimension
  LSD-W4: Mahalanobis-like distance from Lor4D centroid in feature space
           → uses full covariance structure

Scan over (c*, w*, d*, α, β, γ).
"""
from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path
from itertools import product

import numpy as np

from bd_action import count_intervals_fast
from generators import (
    Poset,
    generate_kr_like,
    generate_kr_2layer,
    generate_kr_4layer,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
    generate_transitive_percolation,
    generate_interval_order,
    generate_absolute_layered,
    generate_multi_layer_random,
    generate_random_layered_k4_uniform,
    generate_random_layered_k6_uniform,
    generate_random_layered_k8_uniform,
    generate_random_layered_k6_tapered,
    generate_random_layered_k6_middle_heavy,
    generate_random_layered_k6_longjump,
)
from unified_functional import compute_sigma_hist, compute_xi_dim


def longest_chain_length(poset: Poset) -> int:
    n = poset.n
    if n <= 1:
        return n
    c = poset.closure
    dp = [1] * n
    for i in range(n):
        for j in range(i + 1, n):
            if c[i, j]:
                dp[j] = max(dp[j], dp[i] + 1)
    return max(dp)


def max_antichain_width(poset: Poset) -> int:
    n = poset.n
    if n <= 1:
        return n
    c = poset.closure
    remaining = set(range(n))
    max_w = 0
    while remaining:
        minimals = []
        for i in remaining:
            is_min = True
            for j in remaining:
                if j != i and c[j, i] and not c[i, j]:
                    is_min = False
                    break
            if is_min:
                minimals.append(i)
        if not minimals:
            break
        max_w = max(max_w, len(minimals))
        for m in minimals:
            remaining.discard(m)
    return max_w


def compute_raw_features(poset: Poset, N: int) -> dict:
    counts = count_intervals_fast(poset, k_max=5)
    C0 = counts.get(0)
    C1 = counts.get(1)
    C2 = counts.get(2)
    C3 = counts.get(3)
    total_rel = counts.total_relations

    R = 1.0 - float(C0) / float(total_rel) if total_rel > 0 else 0.0
    c1_c0 = C1 / max(1, C0)
    interval_richness = (C1 + C2 + C3) / max(1, C0)

    xi_val, d_eff = compute_xi_dim(poset)
    sigma_hist = compute_sigma_hist(poset)

    lc = longest_chain_length(poset)
    aw = max_antichain_width(poset)
    chain_ratio = lc / max(1, N)
    width_ratio = aw / max(1, N)

    return {
        "R": R,
        "c1_c0": c1_c0,
        "interval_richness": interval_richness,
        "d_eff": d_eff,
        "sigma_hist": sigma_hist,
        "chain_ratio": chain_ratio,
        "width_ratio": width_ratio,
    }


FAMILIES = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
    "KR_like": generate_kr_like,
    "KR_2layer": generate_kr_2layer,
    "KR_4layer": generate_kr_4layer,
    "AbsLayer": generate_absolute_layered,
    "MLR": generate_multi_layer_random,
    "RLk4": generate_random_layered_k4_uniform,
    "RLk6": generate_random_layered_k6_uniform,
    "RLk8": generate_random_layered_k8_uniform,
    "RLk6_tap": generate_random_layered_k6_tapered,
    "RLk6_mid": generate_random_layered_k6_middle_heavy,
    "RLk6_lj": generate_random_layered_k6_longjump,
    "TransPerc": generate_transitive_percolation,
    "IntOrder": generate_interval_order,
}

CATEGORY = {}
for f in ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]:
    CATEGORY[f] = "Lorentzian"
for f in ["KR_like", "KR_2layer", "KR_4layer"]:
    CATEGORY[f] = "KR-family"
for f in ["AbsLayer", "MLR", "RLk4", "RLk6", "RLk8", "RLk6_tap", "RLk6_mid", "RLk6_lj"]:
    CATEGORY[f] = "Layered"
for f in ["TransPerc", "IntOrder"]:
    CATEGORY[f] = "Other"


def evaluate(by_nf, N_VALUES):
    """Return (lor4d_rank_dict, beaten_frac)."""
    total_beaten = 0
    total_possible = 0
    ranks = {}
    for N in N_VALUES:
        means = {}
        for fam in FAMILIES:
            rows = by_nf.get((N, fam), [])
            if rows:
                means[fam] = np.mean([r["F"] for r in rows])
        ranked = sorted(means, key=means.get)
        r4d = ranked.index("Lor4D") + 1 if "Lor4D" in ranked else 99
        ranks[N] = r4d
        lor4d_mean = means.get("Lor4D", float("inf"))
        non_lor = [f for f in FAMILIES if CATEGORY[f] != "Lorentzian"]
        beaten = sum(1 for f in non_lor if f in means and lor4d_mean < means[f])
        total_beaten += beaten
        total_possible += len(non_lor)
    frac = total_beaten / max(1, total_possible)
    return ranks, frac


def main():
    N_VALUES = [16, 20, 28, 36, 48, 64]
    REPS = 15
    SEED_BASE = 42

    print("=" * 80)
    print("LSD-WELL: Lorentzian Structural Discriminator with Well Potentials")
    print("=" * 80)

    # === Phase 1: Generate all raw features ===
    all_feats = []
    total = len(FAMILIES) * len(N_VALUES) * REPS
    done = 0

    for fam_name, gen_fn in FAMILIES.items():
        for N in N_VALUES:
            for rep in range(REPS):
                seed = SEED_BASE + hash(fam_name) % 10000 + N * 100 + rep
                seed = seed % (2**31)
                try:
                    poset = gen_fn(N, seed=seed)
                    feat = compute_raw_features(poset, N)
                    feat["family"] = fam_name
                    feat["category"] = CATEGORY[fam_name]
                    feat["N"] = N
                    feat["rep"] = rep
                    all_feats.append(feat)
                except Exception as e:
                    print(f"  ERROR: {fam_name} N={N} rep={rep}: {e}")
                done += 1
                if done % 200 == 0:
                    print(f"  [{done}/{total}]", flush=True)

    print(f"  Generated {len(all_feats)} samples.")

    # === Phase 1b: Compute Lor4D centroids per N (for adaptive well centers) ===
    lor4d_centroids = {}
    for N in N_VALUES:
        lor_rows = [f for f in all_feats if f["family"] == "Lor4D" and f["N"] == N]
        if lor_rows:
            lor4d_centroids[N] = {
                "c1_c0": np.mean([r["c1_c0"] for r in lor_rows]),
                "width_ratio": np.mean([r["width_ratio"] for r in lor_rows]),
                "d_eff": np.mean([r["d_eff"] for r in lor_rows]),
                "chain_ratio": np.mean([r["chain_ratio"] for r in lor_rows]),
                "R": np.mean([r["R"] for r in lor_rows]),
            }

    print("\nLor4D centroids:")
    for N in N_VALUES:
        c = lor4d_centroids.get(N, {})
        print(f"  N={N}: c1/c0={c.get('c1_c0',0):.4f}, "
              f"width={c.get('width_ratio',0):.4f}, "
              f"d_eff={c.get('d_eff',0):.2f}, "
              f"chain={c.get('chain_ratio',0):.4f}, "
              f"R={c.get('R',0):.4f}")

    report = []
    report.append("# LSD-Well: Lorentzian Structural Discriminator with Well Potentials\n")
    report.append("**No logH. Wells in interval-shape and transverse-structure space.**\n")
    report.append(f"N = {N_VALUES}, reps = {REPS}\n")

    # === Phase 2: LSD-W1 — (C₁/C₀) well + width well ===
    report.append("\n## LSD-W1: α·(C₁/C₀ − c*)² + β·(width − w*)²\n")
    report.append("**Fixed well centers from Lor4D large-N centroid.**\n")

    # Use N=48 centroid as reference (stable enough)
    c48 = lor4d_centroids.get(48, lor4d_centroids.get(36, {}))
    c_star_ref = c48.get("c1_c0", 0.20)
    w_star_ref = c48.get("width_ratio", 0.41)
    d_star_ref = c48.get("d_eff", 3.93)

    report.append(f"Reference (N=48 centroid): c*={c_star_ref:.4f}, w*={w_star_ref:.4f}, d*={d_star_ref:.2f}\n")

    best_w1 = None
    best_w1_score = -1

    # Scan c*, w* around Lor4D centroid + α, β
    c_stars = [c_star_ref * f for f in [0.8, 1.0, 1.2]]
    w_stars = [w_star_ref * f for f in [0.8, 1.0, 1.2]]
    alphas = [0.5, 1.0, 2.0, 5.0, 10.0]
    betas = [0.5, 1.0, 2.0, 5.0, 10.0]

    for cs, ws, a, b in product(c_stars, w_stars, alphas, betas):
        by_nf = defaultdict(list)
        for feat in all_feats:
            F = a * (feat["c1_c0"] - cs)**2 + b * (feat["width_ratio"] - ws)**2
            fc = dict(feat)
            fc["F"] = F
            by_nf[(feat["N"], feat["family"])].append(fc)
        ranks, frac = evaluate(by_nf, N_VALUES)
        if frac > best_w1_score:
            best_w1_score = frac
            best_w1 = (cs, ws, a, b, ranks)

    cs, ws, a, b, ranks = best_w1
    report.append(f"**Best**: c*={cs:.4f}, w*={ws:.4f}, α={a}, β={b}")
    report.append(f"**Lor4D beats {best_w1_score*100:.1f}% of non-Lor across all N**\n")
    report.append("| N | Lor4D Rank |")
    report.append("|---|:----------:|")
    for N in N_VALUES:
        report.append(f"| {N} | #{ranks[N]}/17 |")

    # Full ranking at each N for best W1
    by_nf_w1 = defaultdict(list)
    for feat in all_feats:
        F = a * (feat["c1_c0"] - cs)**2 + b * (feat["width_ratio"] - ws)**2
        fc = dict(feat)
        fc["F"] = F
        by_nf_w1[(feat["N"], feat["family"])].append(fc)

    for N in N_VALUES:
        means = {}
        for fam in FAMILIES:
            rows = by_nf_w1.get((N, fam), [])
            if rows:
                means[fam] = np.mean([r["F"] for r in rows])
        ranked = sorted(means, key=means.get)
        report.append(f"\n### LSD-W1 N={N}\n")
        report.append("| Rank | Family | Category | F |")
        report.append("|------|--------|----------|:-:|")
        for rank, fam in enumerate(ranked, 1):
            cat = CATEGORY[fam]
            tag = " ◆" if cat == "Lorentzian" else (" ●" if cat == "KR-family" else "")
            report.append(f"| {rank} | {fam}{tag} | {cat} | {means[fam]:.4f} |")

    # === Phase 3: LSD-W2 — triple well (shape + width + d_eff) ===
    report.append("\n\n## LSD-W2: α·(C₁/C₀−c*)² + β·(w−w*)² + γ·(d_eff−d*)²\n")
    report.append("**Triple well: interval shape + width + dimension.**\n")

    best_w2 = None
    best_w2_score = -1

    gammas = [0.5, 1.0, 2.0, 5.0]
    d_stars = [d_star_ref * f for f in [0.95, 1.0, 1.05]]

    for cs, ws, ds, a, b, g in product(
        [c_star_ref], [w_star_ref],
        d_stars,
        [1.0, 2.0, 5.0, 10.0],  # α
        [1.0, 2.0, 5.0],         # β
        gammas                     # γ
    ):
        by_nf = defaultdict(list)
        for feat in all_feats:
            F = (a * (feat["c1_c0"] - cs)**2
                 + b * (feat["width_ratio"] - ws)**2
                 + g * (feat["d_eff"] - ds)**2)
            fc = dict(feat)
            fc["F"] = F
            by_nf[(feat["N"], feat["family"])].append(fc)
        ranks, frac = evaluate(by_nf, N_VALUES)
        if frac > best_w2_score:
            best_w2_score = frac
            best_w2 = (cs, ws, ds, a, b, g, ranks)

    cs, ws, ds, a, b, g, ranks = best_w2
    report.append(f"**Best**: c*={cs:.4f}, w*={ws:.4f}, d*={ds:.2f}, α={a}, β={b}, γ={g}")
    report.append(f"**Lor4D beats {best_w2_score*100:.1f}% of non-Lor across all N**\n")
    report.append("| N | Lor4D Rank |")
    report.append("|---|:----------:|")
    for N in N_VALUES:
        report.append(f"| {N} | #{ranks[N]}/17 |")

    by_nf_w2 = defaultdict(list)
    for feat in all_feats:
        F = (a * (feat["c1_c0"] - cs)**2
             + b * (feat["width_ratio"] - ws)**2
             + g * (feat["d_eff"] - ds)**2)
        fc = dict(feat)
        fc["F"] = F
        by_nf_w2[(feat["N"], feat["family"])].append(fc)

    for N in N_VALUES:
        means = {}
        for fam in FAMILIES:
            rows = by_nf_w2.get((N, fam), [])
            if rows:
                means[fam] = np.mean([r["F"] for r in rows])
        ranked = sorted(means, key=means.get)
        report.append(f"\n### LSD-W2 N={N}\n")
        report.append("| Rank | Family | Category | F |")
        report.append("|------|--------|----------|:-:|")
        for rank, fam in enumerate(ranked, 1):
            cat = CATEGORY[fam]
            tag = " ◆" if cat == "Lorentzian" else (" ●" if cat == "KR-family" else "")
            report.append(f"| {rank} | {fam}{tag} | {cat} | {means[fam]:.4f} |")

    # === Phase 4: LSD-W3 — N-scaled triple well ===
    report.append("\n\n## LSD-W3: N·[α·(C₁/C₀−c*)² + β·(w−w*)² + γ·(d_eff−d*)²]\n")
    report.append("**N-scaled well: O(N) penalty ensures wells compete at all scales.**\n")

    best_w3 = None
    best_w3_score = -1

    for cs, ws, ds, a, b, g in product(
        [c_star_ref], [w_star_ref],
        d_stars,
        [0.5, 1.0, 2.0, 5.0],
        [0.5, 1.0, 2.0],
        [0.2, 0.5, 1.0, 2.0],
    ):
        by_nf = defaultdict(list)
        for feat in all_feats:
            N = feat["N"]
            F = N * (a * (feat["c1_c0"] - cs)**2
                     + b * (feat["width_ratio"] - ws)**2
                     + g * (feat["d_eff"] - ds)**2)
            fc = dict(feat)
            fc["F"] = F
            by_nf[(feat["N"], feat["family"])].append(fc)
        ranks, frac = evaluate(by_nf, N_VALUES)
        if frac > best_w3_score:
            best_w3_score = frac
            best_w3 = (cs, ws, ds, a, b, g, ranks)

    cs, ws, ds, a, b, g, ranks = best_w3
    report.append(f"**Best**: c*={cs:.4f}, w*={ws:.4f}, d*={ds:.2f}, α={a}, β={b}, γ={g}")
    report.append(f"**Lor4D beats {best_w3_score*100:.1f}% of non-Lor across all N**\n")
    report.append("| N | Lor4D Rank |")
    report.append("|---|:----------:|")
    for N in N_VALUES:
        report.append(f"| {N} | #{ranks[N]}/17 |")

    by_nf_w3 = defaultdict(list)
    for feat in all_feats:
        N = feat["N"]
        F = N * (a * (feat["c1_c0"] - cs)**2
                 + b * (feat["width_ratio"] - ws)**2
                 + g * (feat["d_eff"] - ds)**2)
        fc = dict(feat)
        fc["F"] = F
        by_nf_w3[(feat["N"], feat["family"])].append(fc)

    for N in N_VALUES:
        means = {}
        for fam in FAMILIES:
            rows = by_nf_w3.get((N, fam), [])
            if rows:
                means[fam] = np.mean([r["F"] for r in rows])
        ranked = sorted(means, key=means.get)
        report.append(f"\n### LSD-W3 N={N}\n")
        report.append("| Rank | Family | Category | F |")
        report.append("|------|--------|----------|:-:|")
        for rank, fam in enumerate(ranked, 1):
            cat = CATEGORY[fam]
            tag = " ◆" if cat == "Lorentzian" else (" ●" if cat == "KR-family" else "")
            report.append(f"| {rank} | {fam}{tag} | {cat} | {means[fam]:.4f} |")

    # === Phase 5: Pairwise critical check for best overall ===
    # Find overall best
    results = [
        ("LSD-W1", best_w1_score, by_nf_w1),
        ("LSD-W2", best_w2_score, by_nf_w2),
        ("LSD-W3", best_w3_score, by_nf_w3),
    ]
    results.sort(key=lambda x: -x[1])
    winner_name, winner_score, winner_nf = results[0]

    report.append(f"\n\n## Overall Winner: {winner_name} ({winner_score*100:.1f}%)\n")

    report.append("### Lor4D vs KR_2layer across N\n")
    report.append("| N | Lor4D F | KR_2layer F | Lor4D wins? |")
    report.append("|---|:-------:|:-----------:|:-----------:|")
    for N in N_VALUES:
        l4 = winner_nf.get((N, "Lor4D"), [])
        k2 = winner_nf.get((N, "KR_2layer"), [])
        if l4 and k2:
            l4m = np.mean([r["F"] for r in l4])
            k2m = np.mean([r["F"] for r in k2])
            win = "✅" if l4m < k2m else "❌"
            report.append(f"| {N} | {l4m:.4f} | {k2m:.4f} | {win} |")

    report.append("\n### Lor4D vs ALL at N=64\n")
    report.append("| Competitor | Category | Lor4D F | Comp F | Win? |")
    report.append("|------------|----------|:-------:|:------:|:----:|")
    lor4d_vals = winner_nf.get((64, "Lor4D"), [])
    lor4d_64 = np.mean([r["F"] for r in lor4d_vals]) if lor4d_vals else 999
    for fam in sorted(FAMILIES.keys()):
        if fam == "Lor4D":
            continue
        comp = winner_nf.get((64, fam), [])
        if comp:
            cm = np.mean([r["F"] for r in comp])
            win = "✅" if lor4d_64 < cm else "❌"
            report.append(f"| {fam} | {CATEGORY[fam]} | {lor4d_64:.4f} | {cm:.4f} | {win} |")

    # === Phase 6: Comparison table ===
    report.append("\n\n## Comparison: All Generations\n")
    report.append("| Functional | Architecture | Lor4D best rank | N=64 rank | KR_2layer? | Score |")
    report.append("|------------|-------------|:---------------:|:---------:|:----------:|:-----:|")
    report.append("| F7 | logH + wall | #1 (N=16) | #8+ | ❌ N≥28 | ~60% |")
    report.append("| F10 | logH + d_eff-well | #1 (all N) | #1 | ✅ | ~95% |")
    report.append("| F_link | S_link + wall | #3 (N=16) | #8 | ❌ always | ~40% |")
    report.append("| LSD2 | C₁/C₀ + wall | #1 (N≤28) | #10 | ❌ N≥48 | 89.7% |")
    report.append(f"| LSD-W1 | shape+width well | ? | #{best_w1[4].get(64,'?')} | ? | {best_w1_score*100:.1f}% |")
    report.append(f"| LSD-W2 | triple well | ? | #{best_w2[6].get(64,'?')} | ? | {best_w2_score*100:.1f}% |")
    report.append(f"| LSD-W3 | N·triple well | ? | #{best_w3[6].get(64,'?')} | ? | {best_w3_score*100:.1f}% |")

    out = "\n".join(report)
    outdir = Path("outputs_carlip")
    outdir.mkdir(exist_ok=True)
    md_path = outdir / "lsd_well_17family_test.md"
    md_path.write_text(out, encoding="utf-8")
    print(f"\nSaved: {md_path}")
    print("\n" + out)


if __name__ == "__main__":
    main()
