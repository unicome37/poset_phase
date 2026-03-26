"""
LSD-Well N-Adapted: well centers derived per-N from Lor4D centroids
====================================================================
Upgrade from LSD-W2 (constant well centers) to N-adapted version:

  F_LSD(N) = α·(d_eff − 4)² + β·(C₁/C₀ − c*(N))² + γ·(w − w*(N))²

where:
  d* = 4 (fixed, physics input)
  c*(N) = Lor4D centroid of C₁/C₀ at each N
  w*(N) = Lor4D centroid of width_ratio at each N

Only free parameters: α, β, γ (relative weights).

Two evaluation modes:
  MODE A: "Oracle" — c*(N), w*(N) from same-N Lor4D samples (upper bound)
  MODE B: "Extrapolated" — c*(N), w*(N) from power-law fit to smaller N,
          evaluated on held-out larger N (realistic performance)
"""
from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path
from itertools import product

import numpy as np
from scipy.optimize import curve_fit

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
    total_rel = counts.total_relations

    R = 1.0 - float(C0) / float(total_rel) if total_rel > 0 else 0.0
    c1_c0 = C1 / max(1, C0)

    xi_val, d_eff = compute_xi_dim(poset)

    lc = longest_chain_length(poset)
    aw = max_antichain_width(poset)
    width_ratio = aw / max(1, N)

    return {
        "R": R,
        "c1_c0": c1_c0,
        "d_eff": d_eff,
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
for f in ["AbsLayer", "MLR", "RLk4", "RLk6", "RLk8",
          "RLk6_tap", "RLk6_mid", "RLk6_lj"]:
    CATEGORY[f] = "Layered"
for f in ["TransPerc", "IntOrder"]:
    CATEGORY[f] = "Other"


def evaluate(by_nf, N_VALUES):
    """Return (lor4d_rank_dict, beaten_fraction)."""
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
        beaten = sum(1 for f in non_lor
                     if f in means and lor4d_mean < means[f])
        total_beaten += beaten
        total_possible += len(non_lor)
    frac = total_beaten / max(1, total_possible)
    return ranks, frac


def power_law(x, a, b):
    """f(N) = a + b/N for finite-size scaling."""
    return a + b / x


def main():
    N_VALUES = [16, 20, 28, 36, 48, 64]
    REPS = 15
    SEED_BASE = 42

    print("=" * 80)
    print("LSD-Well N-Adapted: N-dependent well centers from Lor4D statistics")
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

    # === Phase 2: Compute Lor4D centroids per N ===
    lor4d_centroids = {}
    for N in N_VALUES:
        lor_rows = [f for f in all_feats
                    if f["family"] == "Lor4D" and f["N"] == N]
        if lor_rows:
            lor4d_centroids[N] = {
                "c1_c0": np.mean([r["c1_c0"] for r in lor_rows]),
                "width_ratio": np.mean([r["width_ratio"] for r in lor_rows]),
                "d_eff": np.mean([r["d_eff"] for r in lor_rows]),
                "c1_c0_std": np.std([r["c1_c0"] for r in lor_rows]),
                "width_std": np.std([r["width_ratio"] for r in lor_rows]),
                "d_eff_std": np.std([r["d_eff"] for r in lor_rows]),
            }

    print("\nLor4D centroids per N:")
    for N in N_VALUES:
        c = lor4d_centroids[N]
        print(f"  N={N:3d}: c1/c0={c['c1_c0']:.4f}±{c['c1_c0_std']:.4f}, "
              f"width={c['width_ratio']:.4f}±{c['width_std']:.4f}, "
              f"d_eff={c['d_eff']:.3f}±{c['d_eff_std']:.3f}")

    # === Phase 3: Fit finite-size scaling for MODE B ===
    Ns = np.array(N_VALUES, dtype=float)
    c_vals = np.array([lor4d_centroids[N]["c1_c0"] for N in N_VALUES])
    w_vals = np.array([lor4d_centroids[N]["width_ratio"] for N in N_VALUES])

    try:
        popt_c, _ = curve_fit(power_law, Ns, c_vals, p0=[0.2, 1.0])
        c_inf, c_b = popt_c
        print(f"\nc*(N) fit: c_inf={c_inf:.4f}, b={c_b:.2f}  (c*(N) = {c_inf:.4f} + {c_b:.2f}/N)")
    except Exception:
        c_inf, c_b = np.mean(c_vals), 0.0
        print(f"\nc*(N) fit failed, using mean: {c_inf:.4f}")

    try:
        popt_w, _ = curve_fit(power_law, Ns, w_vals, p0=[0.4, 1.0])
        w_inf, w_b = popt_w
        print(f"w*(N) fit: w_inf={w_inf:.4f}, b={w_b:.2f}  (w*(N) = {w_inf:.4f} + {w_b:.2f}/N)")
    except Exception:
        w_inf, w_b = np.mean(w_vals), 0.0
        print(f"w*(N) fit failed, using mean: {w_inf:.4f}")

    report = []
    report.append("# LSD-Well N-Adapted: N-dependent Well Centers\n")
    report.append("**Upgrade**: well centers c*(N), w*(N) no longer fixed constants.\n")
    report.append(f"N = {N_VALUES}, reps = {REPS}\n")

    # --- Table: Lor4D centroid evolution ---
    report.append("\n## 1. Lor4D Centroid Evolution\n")
    report.append("| N | c*(N) = C₁/C₀ | σ(c) | w*(N) = width | σ(w) | d_eff | σ(d) |")
    report.append("|---|:-:|:-:|:-:|:-:|:-:|:-:|")
    for N in N_VALUES:
        c = lor4d_centroids[N]
        report.append(f"| {N} | {c['c1_c0']:.4f} | {c['c1_c0_std']:.4f} | "
                      f"{c['width_ratio']:.4f} | {c['width_std']:.4f} | "
                      f"{c['d_eff']:.3f} | {c['d_eff_std']:.3f} |")

    report.append(f"\n**Finite-size fit**: c*(N) = {c_inf:.4f} + {c_b:.2f}/N")
    report.append(f"**Finite-size fit**: w*(N) = {w_inf:.4f} + {w_b:.2f}/N")
    report.append("**d* = 4 (fixed, physics input)**\n")

    # =======================================================
    # MODE A: Oracle — use same-N Lor4D centroids
    # =======================================================
    report.append("\n## 2. MODE A: Oracle N-Adapted (same-N centroids)\n")
    report.append("F(N) = α·(d_eff − 4)² + β·(C₁/C₀ − c*(N))² + γ·(w − w*(N))²\n")
    report.append("where c*(N), w*(N) = Lor4D centroid at each N.\n")

    best_a = None
    best_a_score = -1

    for a, b, g in product(
        [0.5, 1.0, 2.0, 5.0, 10.0],   # α (d_eff weight)
        [1.0, 2.0, 5.0, 10.0, 20.0],   # β (c1_c0 weight)
        [1.0, 2.0, 5.0, 10.0, 20.0],   # γ (width weight)
    ):
        by_nf = defaultdict(list)
        for feat in all_feats:
            N = feat["N"]
            cN = lor4d_centroids[N]["c1_c0"]
            wN = lor4d_centroids[N]["width_ratio"]
            F = (a * (feat["d_eff"] - 4.0)**2
                 + b * (feat["c1_c0"] - cN)**2
                 + g * (feat["width_ratio"] - wN)**2)
            fc = dict(feat)
            fc["F"] = F
            by_nf[(feat["N"], feat["family"])].append(fc)
        ranks, frac = evaluate(by_nf, N_VALUES)
        if frac > best_a_score:
            best_a_score = frac
            best_a = (a, b, g, ranks)

    a_a, a_b, a_g, a_ranks = best_a
    report.append(f"**Best**: α={a_a}, β={a_b}, γ={a_g}")
    report.append(f"**Lor4D beats {best_a_score*100:.1f}% of non-Lor across all N**\n")
    report.append("| N | Lor4D Rank | c*(N) | w*(N) |")
    report.append("|---|:----------:|:-----:|:-----:|")
    for N in N_VALUES:
        cN = lor4d_centroids[N]["c1_c0"]
        wN = lor4d_centroids[N]["width_ratio"]
        report.append(f"| {N} | #{a_ranks[N]}/17 | {cN:.4f} | {wN:.4f} |")

    # Full rankings for MODE A
    by_nf_a = defaultdict(list)
    for feat in all_feats:
        N = feat["N"]
        cN = lor4d_centroids[N]["c1_c0"]
        wN = lor4d_centroids[N]["width_ratio"]
        F = (a_a * (feat["d_eff"] - 4.0)**2
             + a_b * (feat["c1_c0"] - cN)**2
             + a_g * (feat["width_ratio"] - wN)**2)
        fc = dict(feat)
        fc["F"] = F
        by_nf_a[(feat["N"], feat["family"])].append(fc)

    for N in N_VALUES:
        means = {}
        for fam in FAMILIES:
            rows = by_nf_a.get((N, fam), [])
            if rows:
                means[fam] = np.mean([r["F"] for r in rows])
        ranked = sorted(means, key=means.get)
        report.append(f"\n### MODE A: N={N}\n")
        report.append("| Rank | Family | Category | F |")
        report.append("|------|--------|----------|:-:|")
        for rank, fam in enumerate(ranked, 1):
            cat = CATEGORY[fam]
            tag = " ◆" if cat == "Lorentzian" else (" ●" if cat == "KR-family" else "")
            report.append(f"| {rank} | {fam}{tag} | {cat} | {means[fam]:.4f} |")

    # =======================================================
    # MODE B: Extrapolated — power-law fit from N<=36, test on N>=48
    # =======================================================
    report.append("\n\n## 3. MODE B: Extrapolated N-Adapted\n")
    report.append("Fit c*(N) = c_∞ + b_c/N and w*(N) = w_∞ + b_w/N from N ≤ 36,\n")
    report.append("then predict well centers for N = 48, 64.\n")

    train_Ns = [N for N in N_VALUES if N <= 36]
    test_Ns = [N for N in N_VALUES if N > 36]

    Ns_train = np.array(train_Ns, dtype=float)
    c_train = np.array([lor4d_centroids[N]["c1_c0"] for N in train_Ns])
    w_train = np.array([lor4d_centroids[N]["width_ratio"] for N in train_Ns])

    try:
        popt_c2, _ = curve_fit(power_law, Ns_train, c_train, p0=[0.2, 1.0])
        c_inf2, c_b2 = popt_c2
    except Exception:
        c_inf2, c_b2 = np.mean(c_train), 0.0

    try:
        popt_w2, _ = curve_fit(power_law, Ns_train, w_train, p0=[0.4, 1.0])
        w_inf2, w_b2 = popt_w2
    except Exception:
        w_inf2, w_b2 = np.mean(w_train), 0.0

    report.append(f"\n**Train fit (N ≤ 36)**: c*(N) = {c_inf2:.4f} + {c_b2:.2f}/N")
    report.append(f"**Train fit (N ≤ 36)**: w*(N) = {w_inf2:.4f} + {w_b2:.2f}/N\n")

    # Show predicted vs actual for test N
    report.append("| N | c*(N) pred | c*(N) actual | Δc | w*(N) pred | w*(N) actual | Δw |")
    report.append("|---|:-:|:-:|:-:|:-:|:-:|:-:|")
    for N in N_VALUES:
        c_pred = power_law(N, c_inf2, c_b2)
        w_pred = power_law(N, w_inf2, w_b2)
        c_actual = lor4d_centroids[N]["c1_c0"]
        w_actual = lor4d_centroids[N]["width_ratio"]
        tag = " (test)" if N in test_Ns else ""
        report.append(f"| {N}{tag} | {c_pred:.4f} | {c_actual:.4f} | "
                      f"{abs(c_pred-c_actual):.4f} | "
                      f"{w_pred:.4f} | {w_actual:.4f} | "
                      f"{abs(w_pred-w_actual):.4f} |")

    # Evaluate MODE B with extrapolated centers
    best_b = None
    best_b_score = -1

    for a, b, g in product(
        [0.5, 1.0, 2.0, 5.0, 10.0],
        [1.0, 2.0, 5.0, 10.0, 20.0],
        [1.0, 2.0, 5.0, 10.0, 20.0],
    ):
        by_nf = defaultdict(list)
        for feat in all_feats:
            N = feat["N"]
            cN = power_law(N, c_inf2, c_b2)
            wN = power_law(N, w_inf2, w_b2)
            F = (a * (feat["d_eff"] - 4.0)**2
                 + b * (feat["c1_c0"] - cN)**2
                 + g * (feat["width_ratio"] - wN)**2)
            fc = dict(feat)
            fc["F"] = F
            by_nf[(feat["N"], feat["family"])].append(fc)
        ranks, frac = evaluate(by_nf, N_VALUES)
        if frac > best_b_score:
            best_b_score = frac
            best_b = (a, b, g, ranks)

    b_a, b_b, b_g, b_ranks = best_b
    report.append(f"\n**Best**: α={b_a}, β={b_b}, γ={b_g}")
    report.append(f"**Lor4D beats {best_b_score*100:.1f}% of non-Lor across all N**\n")
    report.append("| N | Lor4D Rank |")
    report.append("|---|:----------:|")
    for N in N_VALUES:
        report.append(f"| {N} | #{b_ranks[N]}/17 |")

    # Full rankings for MODE B
    by_nf_b = defaultdict(list)
    for feat in all_feats:
        N = feat["N"]
        cN = power_law(N, c_inf2, c_b2)
        wN = power_law(N, w_inf2, w_b2)
        F = (b_a * (feat["d_eff"] - 4.0)**2
             + b_b * (feat["c1_c0"] - cN)**2
             + b_g * (feat["width_ratio"] - wN)**2)
        fc = dict(feat)
        fc["F"] = F
        by_nf_b[(feat["N"], feat["family"])].append(fc)

    for N in N_VALUES:
        means = {}
        for fam in FAMILIES:
            rows = by_nf_b.get((N, fam), [])
            if rows:
                means[fam] = np.mean([r["F"] for r in rows])
        ranked = sorted(means, key=means.get)
        report.append(f"\n### MODE B: N={N}\n")
        report.append("| Rank | Family | Category | F |")
        report.append("|------|--------|----------|:-:|")
        for rank, fam in enumerate(ranked, 1):
            cat = CATEGORY[fam]
            tag = " ◆" if cat == "Lorentzian" else (" ●" if cat == "KR-family" else "")
            report.append(f"| {rank} | {fam}{tag} | {cat} | {means[fam]:.4f} |")

    # =======================================================
    # MODE C: Constant-center baseline (original LSD-W2)
    # =======================================================
    report.append("\n\n## 4. Baseline: Constant-Center LSD-W2\n")
    report.append("For comparison: use N=48 centroid as fixed well center.\n")

    c48 = lor4d_centroids.get(48, lor4d_centroids.get(36, {}))
    c_const = c48.get("c1_c0", 0.20)
    w_const = c48.get("width_ratio", 0.41)

    report.append(f"Fixed centers: c* = {c_const:.4f}, w* = {w_const:.4f}, d* = 4\n")

    best_c = None
    best_c_score = -1

    for a, b, g in product(
        [0.5, 1.0, 2.0, 5.0, 10.0],
        [1.0, 2.0, 5.0, 10.0, 20.0],
        [1.0, 2.0, 5.0, 10.0, 20.0],
    ):
        by_nf = defaultdict(list)
        for feat in all_feats:
            F = (a * (feat["d_eff"] - 4.0)**2
                 + b * (feat["c1_c0"] - c_const)**2
                 + g * (feat["width_ratio"] - w_const)**2)
            fc = dict(feat)
            fc["F"] = F
            by_nf[(feat["N"], feat["family"])].append(fc)
        ranks, frac = evaluate(by_nf, N_VALUES)
        if frac > best_c_score:
            best_c_score = frac
            best_c = (a, b, g, ranks)

    c_a, c_b_w, c_g, c_ranks = best_c
    report.append(f"**Best**: α={c_a}, β={c_b_w}, γ={c_g}")
    report.append(f"**Lor4D beats {best_c_score*100:.1f}% of non-Lor across all N**\n")
    report.append("| N | Lor4D Rank |")
    report.append("|---|:----------:|")
    for N in N_VALUES:
        report.append(f"| {N} | #{c_ranks[N]}/17 |")

    # =======================================================
    # Final comparison
    # =======================================================
    report.append("\n\n## 5. Head-to-Head Comparison\n")
    report.append("| Mode | Description | Beat% | Worst Rank |")
    report.append("|------|------------|:-----:|:----------:|")

    worst_a = max(a_ranks.values())
    worst_b = max(b_ranks.values())
    worst_c = max(c_ranks.values())

    report.append(f"| A (Oracle) | c*(N), w*(N) from same-N | "
                  f"{best_a_score*100:.1f}% | #{worst_a}/17 |")
    report.append(f"| B (Extrap) | c*(N), w*(N) power-law fit N≤36 | "
                  f"{best_b_score*100:.1f}% | #{worst_b}/17 |")
    report.append(f"| C (Const)  | Fixed c*, w* from N=48 | "
                  f"{best_c_score*100:.1f}% | #{worst_c}/17 |")

    report.append("\n### Per-N Rank Comparison\n")
    report.append("| N | A (Oracle) | B (Extrap) | C (Const) |")
    report.append("|---|:----------:|:----------:|:---------:|")
    for N in N_VALUES:
        report.append(f"| {N} | #{a_ranks[N]} | #{b_ranks[N]} | #{c_ranks[N]} |")

    # === Critical pairwise checks ===
    report.append("\n### Lor4D vs KR_2layer (Oracle mode)\n")
    report.append("| N | Lor4D F | KR_2layer F | Gap | Win? |")
    report.append("|---|:------:|:----------:|:---:|:----:|")
    for N in N_VALUES:
        l4 = by_nf_a.get((N, "Lor4D"), [])
        k2 = by_nf_a.get((N, "KR_2layer"), [])
        if l4 and k2:
            l4m = np.mean([r["F"] for r in l4])
            k2m = np.mean([r["F"] for r in k2])
            win = "✅" if l4m < k2m else "❌"
            report.append(f"| {N} | {l4m:.4f} | {k2m:.4f} | "
                          f"{k2m-l4m:+.4f} | {win} |")

    # === Conclusion ===
    report.append("\n\n## 6. Conclusion\n")
    if best_a_score > best_c_score:
        delta = (best_a_score - best_c_score) * 100
        report.append(f"N-adapted wells improve discrimination by "
                      f"+{delta:.1f}pp over constant centers.\n")
    elif best_a_score == best_c_score:
        report.append("N-adapted wells match constant-center performance.\n")
    else:
        delta = (best_c_score - best_a_score) * 100
        report.append(f"Constant centers unexpectedly outperform by "
                      f"+{delta:.1f}pp — N-adaptation may not help here.\n")

    if best_b_score >= best_c_score:
        report.append("Extrapolated mode matches or exceeds constant baseline → "
                      "power-law scaling is viable for prediction.\n")
    else:
        delta = (best_c_score - best_b_score) * 100
        report.append(f"Extrapolated mode loses {delta:.1f}pp to constant baseline → "
                      "power-law fit may not capture true scaling.\n")

    report.append("**Final form:**\n")
    report.append("```")
    report.append("F_LSD(N) = α·(d_eff − 4)² + β·(C₁/C₀ − c*(N))² + γ·(w − w*(N))²")
    report.append(f"  d* = 4 (fixed)")
    report.append(f"  c*(N) = {c_inf:.4f} + {c_b:.2f}/N")
    report.append(f"  w*(N) = {w_inf:.4f} + {w_b:.2f}/N")
    report.append(f"  Best weights (oracle): α={a_a}, β={a_b}, γ={a_g}")
    report.append("```")

    out = "\n".join(report)
    outdir = Path("outputs_carlip")
    outdir.mkdir(exist_ok=True)
    md_path = outdir / "lsd_well_n_adapted.md"
    md_path.write_text(out, encoding="utf-8")
    print(f"\nSaved: {md_path}")
    print("\n" + out)


if __name__ == "__main__":
    main()
