"""
Lorentzian Structural Discriminator (LSD) — Post-Carlip Prototype
===================================================================
Uses INTERVAL SHAPE and TRANSVERSE STRUCTURE as primary features.
No logH. No d_eff-well. No link density as primary ordering.

Key insight from info-geometry analysis:
  - d_eff and f₂ CANNOT separate Lor4D from KR_2layer
  - C₁/C₀ and width_ratio CAN, with effect sizes growing with N

Design principles:
  - Lower F = more favored (convention matches F7)
  - Lorentzian structures should have LOW F
  - All features are pure causal-geometric (no logH)

Candidates:
  LSD1: -α·(C₁/C₀) - β·(1-width_ratio)
         → rewards non-trivial intervals, penalizes wide antichains
  LSD2: LSD1 + wall(R)
         → adds sparse-structure penalty  
  LSD3: -α·(C₁/C₀) - β·(1-width_ratio) - γ·chain_ratio + wall(R)
         → also rewards deep chains
  LSD4: -α·interval_richness - β·(1-width_ratio) + wall(R)
         → interval_richness = (C₁+C₂+C₃)/C₀ (total non-link fraction)

For each: scan over parameters, report Lor4D ranking and pairwise wins.
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
from unified_functional import compute_sigma_hist


def sigmoid(x: float) -> float:
    x = max(-30.0, min(30.0, x))
    return 1.0 / (1.0 + math.exp(-x))


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
    n_pairs = N * (N - 1) // 2

    R = 1.0 - float(C0) / float(total_rel) if total_rel > 0 else 0.0
    c1_c0 = C1 / max(1, C0)
    c2_c0 = C2 / max(1, C0)
    c3_c0 = C3 / max(1, C0)
    interval_richness = (C1 + C2 + C3) / max(1, C0)

    lc = longest_chain_length(poset)
    aw = max_antichain_width(poset)
    chain_ratio = lc / max(1, N)
    width_ratio = aw / max(1, N)

    return {
        "C0": C0, "C1": C1, "C2": C2, "C3": C3,
        "total_rel": total_rel,
        "R": R,
        "c1_c0": c1_c0,
        "c2_c0": c2_c0,
        "c3_c0": c3_c0,
        "interval_richness": interval_richness,
        "chain_ratio": chain_ratio,
        "width_ratio": width_ratio,
    }


def wall_term(R: float, N: int) -> float:
    """R-based wall: penalize sparse structures (high R)."""
    alpha_N = 2.0 * (20.0 / max(N, 1)) ** 0.5
    return alpha_N * sigmoid((R - 0.25) / 0.015)


def lsd1(feat: dict, alpha: float, beta: float) -> float:
    return -alpha * feat["c1_c0"] - beta * (1.0 - feat["width_ratio"])


def lsd2(feat: dict, alpha: float, beta: float, N: int) -> float:
    return lsd1(feat, alpha, beta) + wall_term(feat["R"], N)


def lsd3(feat: dict, alpha: float, beta: float, gamma: float, N: int) -> float:
    return (lsd1(feat, alpha, beta)
            - gamma * feat["chain_ratio"]
            + wall_term(feat["R"], N))


def lsd4(feat: dict, alpha: float, beta: float, N: int) -> float:
    return (-alpha * feat["interval_richness"]
            - beta * (1.0 - feat["width_ratio"])
            + wall_term(feat["R"], N))


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


def evaluate_functional(by_nf, N, func_fn, families, category):
    """Return (ranking_list, lor4d_rank, n_non_lor_beaten)."""
    means = {}
    for fam in families:
        rows = by_nf.get((N, fam), [])
        if rows:
            means[fam] = np.mean([r["F"] for r in rows])
    ranked = sorted(means, key=means.get)
    lor4d_rank = ranked.index("Lor4D") + 1 if "Lor4D" in ranked else -1

    lor4d_mean = means.get("Lor4D", float("inf"))
    non_lor = [f for f in families if category[f] != "Lorentzian"]
    beaten = sum(1 for f in non_lor if f in means and lor4d_mean < means[f])

    return ranked, lor4d_rank, beaten, len(non_lor), means


def main():
    N_VALUES = [16, 20, 28, 36, 48, 64]
    REPS = 12
    SEED_BASE = 99

    print("=" * 80)
    print("LORENTZIAN STRUCTURAL DISCRIMINATOR (LSD) — 17 Families")
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

    # === Phase 2: Parameter scan for each LSD variant ===
    report = []
    report.append("# Lorentzian Structural Discriminator (LSD) — 17 Family Test\n")
    report.append("**No logH. No d_eff-well. Pure interval shape + transverse structure.**\n")
    report.append(f"N = {N_VALUES}, reps = {REPS}\n")

    # LSD candidate definitions with parameter grids
    candidates = {
        "LSD1: -α·(C₁/C₀) - β·(1-w)": {
            "params": list(product(
                [0.5, 1.0, 2.0, 3.0, 5.0],     # alpha
                [0.5, 1.0, 2.0, 3.0],            # beta
            )),
            "fn": lambda feat, N, p: lsd1(feat, p[0], p[1]),
            "plabels": ["α", "β"],
        },
        "LSD2: -α·(C₁/C₀) - β·(1-w) + wall": {
            "params": list(product(
                [0.5, 1.0, 2.0, 3.0, 5.0],
                [0.5, 1.0, 2.0, 3.0],
            )),
            "fn": lambda feat, N, p: lsd2(feat, p[0], p[1], N),
            "plabels": ["α", "β"],
        },
        "LSD3: -α·(C₁/C₀) - β·(1-w) - γ·chain + wall": {
            "params": list(product(
                [1.0, 2.0, 3.0],
                [1.0, 2.0],
                [1.0, 2.0, 3.0],
            )),
            "fn": lambda feat, N, p: lsd3(feat, p[0], p[1], p[2], N),
            "plabels": ["α", "β", "γ"],
        },
        "LSD4: -α·IR - β·(1-w) + wall": {
            "params": list(product(
                [0.5, 1.0, 2.0, 3.0],
                [0.5, 1.0, 2.0],
            )),
            "fn": lambda feat, N, p: lsd4(feat, p[0], p[1], N),
            "plabels": ["α", "β"],
        },
    }

    best_overall = None
    best_overall_score = -1

    for cand_name, cand in candidates.items():
        report.append(f"\n## {cand_name}\n")
        best_params = None
        best_score = -1

        for params in cand["params"]:
            # Compute F for all samples
            by_nf = defaultdict(list)
            for feat in all_feats:
                feat_copy = dict(feat)
                feat_copy["F"] = cand["fn"](feat, feat["N"], params)
                by_nf[(feat["N"], feat["family"])].append(feat_copy)

            # Score: sum of (non-Lor families beaten by Lor4D) across all N
            total_beaten = 0
            total_possible = 0
            lor4d_ranks = {}
            for N in N_VALUES:
                ranked, r4d, beaten, n_non, means = evaluate_functional(
                    by_nf, N, None, FAMILIES, CATEGORY
                )
                lor4d_ranks[N] = r4d
                total_beaten += beaten
                total_possible += n_non

            score = total_beaten / max(1, total_possible)
            if score > best_score:
                best_score = score
                best_params = params
                best_ranks = dict(lor4d_ranks)

        # Report best params for this candidate
        param_str = ", ".join(f"{l}={v}" for l, v in zip(cand["plabels"], best_params))
        report.append(f"**Best params**: {param_str}")
        report.append(f"**Lor4D beats {best_score*100:.1f}% of non-Lor across all N**\n")
        report.append("| N | Lor4D Rank | Beats non-Lor |")
        report.append("|---|:----------:|:-------------:|")

        # Recompute with best params for detailed output
        by_nf_best = defaultdict(list)
        for feat in all_feats:
            feat_copy = dict(feat)
            feat_copy["F"] = cand["fn"](feat, feat["N"], best_params)
            by_nf_best[(feat["N"], feat["family"])].append(feat_copy)

        for N in N_VALUES:
            ranked, r4d, beaten, n_non, means = evaluate_functional(
                by_nf_best, N, None, FAMILIES, CATEGORY
            )
            report.append(f"| {N} | #{r4d}/17 | {beaten}/{n_non} |")

        if best_score > best_overall_score:
            best_overall_score = best_score
            best_overall = (cand_name, best_params, cand)

        # Full ranking table for best params at selected N values
        for N in [20, 36, 48, 64]:
            ranked, r4d, beaten, n_non, means = evaluate_functional(
                by_nf_best, N, None, FAMILIES, CATEGORY
            )
            report.append(f"\n### {cand_name} — N={N} Full Ranking\n")
            report.append("| Rank | Family | Category | F |")
            report.append("|------|--------|----------|:-:|")
            for rank, fam in enumerate(ranked, 1):
                cat = CATEGORY[fam]
                tag = " ◆" if cat == "Lorentzian" else (" ●" if cat == "KR-family" else "")
                report.append(f"| {rank} | {fam}{tag} | {cat} | {means[fam]:.4f} |")

    # === Phase 3: Winner analysis ===
    report.append("\n## Overall Winner\n")
    if best_overall:
        wname, wparams, wcand = best_overall
        wlabels = wcand["plabels"]
        wpstr = ", ".join(f"{l}={v}" for l, v in zip(wlabels, wparams))
        report.append(f"**{wname}** with {wpstr}")
        report.append(f"**Lor4D beats {best_overall_score*100:.1f}% of non-Lorentzian families**\n")

        # Detailed pairwise at N=48
        by_nf_win = defaultdict(list)
        for feat in all_feats:
            fc = dict(feat)
            fc["F"] = wcand["fn"](feat, feat["N"], wparams)
            by_nf_win[(feat["N"], feat["family"])].append(fc)

        report.append("### Lor4D Pairwise Wins at N=48\n")
        report.append("| Competitor | Lor4D F | Comp F | Lor4D wins? |")
        report.append("|------------|:-------:|:------:|:-----------:|")

        lor4d_vals = by_nf_win.get((48, "Lor4D"), [])
        lor4d_mean = np.mean([r["F"] for r in lor4d_vals]) if lor4d_vals else 0

        for fam in sorted(FAMILIES.keys()):
            if fam == "Lor4D":
                continue
            comp_vals = by_nf_win.get((48, fam), [])
            if comp_vals:
                comp_mean = np.mean([r["F"] for r in comp_vals])
                win = "✅" if lor4d_mean < comp_mean else "❌"
                report.append(f"| {fam} | {lor4d_mean:.4f} | {comp_mean:.4f} | {win} |")

        # KR_2layer specific check across all N
        report.append("\n### Critical Check: Lor4D vs KR_2layer across N\n")
        report.append("| N | Lor4D F | KR_2layer F | Lor4D wins? |")
        report.append("|---|:-------:|:-----------:|:-----------:|")
        for N in N_VALUES:
            l4 = by_nf_win.get((N, "Lor4D"), [])
            k2 = by_nf_win.get((N, "KR_2layer"), [])
            if l4 and k2:
                l4m = np.mean([r["F"] for r in l4])
                k2m = np.mean([r["F"] for r in k2])
                win = "✅" if l4m < k2m else "❌"
                report.append(f"| {N} | {l4m:.4f} | {k2m:.4f} | {win} |")

    out = "\n".join(report)
    outdir = Path("outputs_carlip")
    outdir.mkdir(exist_ok=True)
    md_path = outdir / "lsd_17family_test.md"
    md_path.write_text(out, encoding="utf-8")
    print(f"\nSaved: {md_path}")
    print("\n" + out)


if __name__ == "__main__":
    main()
