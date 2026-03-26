"""
LSD-Well Large-N Scalability Test
===================================
Extend the N-adapted LSD-Well experiment from N=16–64 to N=16–256.

Key scientific questions:
  1. Does Lor4D remain #1/17 at N=96, 128, 192, 256?
  2. Does the finite-size scaling c*(N) = c_∞ + b/N extrapolate correctly?
  3. How does the Lor4D-to-runner-up margin evolve with N?
  4. At what N (if any) does discrimination break down?

Design:
  - N = [16, 20, 28, 36, 48, 64, 96, 128, 192, 256]
  - reps = 15 for N ≤ 64, 10 for N > 64 (runtime budget)
  - Uses best weights from prior experiment: α=0.5, β=1.0, γ=5.0
  - Also runs grid search to check if optimal weights shift at large N
  - MODE A (oracle) + MODE B (extrapolated from N ≤ 64)
"""
from __future__ import annotations

import math
import time
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


# ── helpers ──────────────────────────────────────────────────────

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


def power_law(x, a, b):
    return a + b / x


def evaluate(by_nf, N_VALUES):
    total_beaten = 0
    total_possible = 0
    ranks = {}
    margins = {}
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
        # margin to runner-up
        non_lor_means = [means[f] for f in non_lor if f in means]
        if non_lor_means:
            margins[N] = min(non_lor_means) - lor4d_mean
        else:
            margins[N] = 0.0
    frac = total_beaten / max(1, total_possible)
    return ranks, frac, margins


def main():
    N_VALUES = [16, 20, 28, 36, 48, 64, 96, 128, 192, 256]
    REPS_SMALL = 15   # for N <= 64
    REPS_LARGE = 10   # for N > 64
    SEED_BASE = 42

    print("=" * 80)
    print("LSD-Well Large-N Scalability Test")
    print(f"N = {N_VALUES}")
    print("=" * 80)

    # === Phase 1: Generate all raw features ===
    all_feats = []
    total = sum(len(FAMILIES) * (REPS_SMALL if N <= 64 else REPS_LARGE)
                for N in N_VALUES)
    done = 0
    t0 = time.time()

    for fam_name, gen_fn in FAMILIES.items():
        for N in N_VALUES:
            reps = REPS_SMALL if N <= 64 else REPS_LARGE
            for rep in range(reps):
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
                if done % 100 == 0:
                    elapsed = time.time() - t0
                    rate = done / max(elapsed, 0.01)
                    eta = (total - done) / max(rate, 0.01)
                    print(f"  [{done}/{total}] {elapsed:.0f}s elapsed, "
                          f"ETA {eta:.0f}s", flush=True)

    elapsed = time.time() - t0
    print(f"\n  Finished generation: {len(all_feats)} samples in {elapsed:.1f}s\n")

    # === Phase 2: Lor4D centroids per N ===
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

    print("Lor4D centroids per N:")
    for N in N_VALUES:
        c = lor4d_centroids.get(N)
        if c:
            print(f"  N={N:4d}: c1/c0={c['c1_c0']:.4f}±{c['c1_c0_std']:.4f}, "
                  f"width={c['width_ratio']:.4f}±{c['width_std']:.4f}, "
                  f"d_eff={c['d_eff']:.3f}±{c['d_eff_std']:.3f}")

    # === Phase 3: Finite-size scaling fits ===
    Ns = np.array(N_VALUES, dtype=float)
    c_vals = np.array([lor4d_centroids[N]["c1_c0"] for N in N_VALUES])
    w_vals = np.array([lor4d_centroids[N]["width_ratio"] for N in N_VALUES])

    # Full fit (all N)
    popt_c_all, _ = curve_fit(power_law, Ns, c_vals, p0=[0.2, 1.0])
    popt_w_all, _ = curve_fit(power_law, Ns, w_vals, p0=[0.4, 1.0])
    print(f"\nFull fit (all N): c*(N) = {popt_c_all[0]:.4f} + {popt_c_all[1]:.2f}/N")
    print(f"Full fit (all N): w*(N) = {popt_w_all[0]:.4f} + {popt_w_all[1]:.2f}/N")

    # Restricted fit (N <= 64, same as before)
    small_mask = Ns <= 64
    popt_c_small, _ = curve_fit(power_law, Ns[small_mask], c_vals[small_mask], p0=[0.2, 1.0])
    popt_w_small, _ = curve_fit(power_law, Ns[small_mask], w_vals[small_mask], p0=[0.4, 1.0])
    print(f"Small fit (N≤64): c*(N) = {popt_c_small[0]:.4f} + {popt_c_small[1]:.2f}/N")
    print(f"Small fit (N≤64): w*(N) = {popt_w_small[0]:.4f} + {popt_w_small[1]:.2f}/N")

    # Prediction accuracy: how well does N≤64 fit predict N>64?
    print("\n  Extrapolation accuracy (N≤64 fit → large N):")
    for N in N_VALUES:
        if N > 64:
            c_pred = power_law(N, popt_c_small[0], popt_c_small[1])
            w_pred = power_law(N, popt_w_small[0], popt_w_small[1])
            c_act = lor4d_centroids[N]["c1_c0"]
            w_act = lor4d_centroids[N]["width_ratio"]
            print(f"  N={N:4d}: c pred={c_pred:.4f} act={c_act:.4f} Δ={abs(c_pred-c_act):.4f}"
                  f"  |  w pred={w_pred:.4f} act={w_act:.4f} Δ={abs(w_pred-w_act):.4f}")

    # === Phase 4: Build report ===
    report = []
    report.append("# LSD-Well Large-N Scalability Test\n")
    report.append(f"N = {N_VALUES}")
    report.append(f"reps = {REPS_SMALL} (N≤64), {REPS_LARGE} (N>64)")
    report.append(f"Total samples: {len(all_feats)}")
    report.append(f"Generation time: {elapsed:.1f}s\n")

    # --- Table: Lor4D centroid evolution ---
    report.append("\n## 1. Lor4D Centroid Evolution (Extended Range)\n")
    report.append("| N | c*(N) | σ(c) | w*(N) | σ(w) | d_eff | σ(d) |")
    report.append("|---|:-:|:-:|:-:|:-:|:-:|:-:|")
    for N in N_VALUES:
        c = lor4d_centroids.get(N)
        if c:
            report.append(f"| {N} | {c['c1_c0']:.4f} | {c['c1_c0_std']:.4f} | "
                          f"{c['width_ratio']:.4f} | {c['width_std']:.4f} | "
                          f"{c['d_eff']:.3f} | {c['d_eff_std']:.3f} |")

    report.append(f"\n**Full fit (all N)**: c*(N) = {popt_c_all[0]:.4f} + {popt_c_all[1]:.2f}/N")
    report.append(f"**Full fit (all N)**: w*(N) = {popt_w_all[0]:.4f} + {popt_w_all[1]:.2f}/N")
    report.append(f"**Small fit (N≤64)**: c*(N) = {popt_c_small[0]:.4f} + {popt_c_small[1]:.2f}/N")
    report.append(f"**Small fit (N≤64)**: w*(N) = {popt_w_small[0]:.4f} + {popt_w_small[1]:.2f}/N\n")

    # Extrapolation table
    report.append("### Extrapolation Accuracy (N≤64 fit → large N)\n")
    report.append("| N | c pred | c actual | |Δc| | w pred | w actual | |Δw| |")
    report.append("|---|:-:|:-:|:-:|:-:|:-:|:-:|")
    for N in N_VALUES:
        if N > 64:
            c_pred = power_law(N, popt_c_small[0], popt_c_small[1])
            w_pred = power_law(N, popt_w_small[0], popt_w_small[1])
            c_act = lor4d_centroids[N]["c1_c0"]
            w_act = lor4d_centroids[N]["width_ratio"]
            report.append(f"| {N} | {c_pred:.4f} | {c_act:.4f} | {abs(c_pred-c_act):.4f} | "
                          f"{w_pred:.4f} | {w_act:.4f} | {abs(w_pred-w_act):.4f} |")

    # === Phase 5: Evaluate with fixed best weights α=0.5, β=1.0, γ=5.0 ===
    ALPHA, BETA, GAMMA = 0.5, 1.0, 5.0

    report.append(f"\n\n## 2. Fixed Weights (α={ALPHA}, β={BETA}, γ={GAMMA})\n")

    # MODE A: Oracle
    report.append("### MODE A: Oracle (same-N centroids)\n")
    by_nf_a = defaultdict(list)
    for feat in all_feats:
        N = feat["N"]
        cN = lor4d_centroids[N]["c1_c0"]
        wN = lor4d_centroids[N]["width_ratio"]
        F = (ALPHA * (feat["d_eff"] - 4.0)**2
             + BETA * (feat["c1_c0"] - cN)**2
             + GAMMA * (feat["width_ratio"] - wN)**2)
        fc = dict(feat)
        fc["F"] = F
        by_nf_a[(feat["N"], feat["family"])].append(fc)

    ranks_a, frac_a, margins_a = evaluate(by_nf_a, N_VALUES)
    report.append(f"**Beat rate**: {frac_a*100:.1f}%\n")
    report.append("| N | Lor4D Rank | Margin to runner-up |")
    report.append("|---|:----------:|:-------------------:|")
    for N in N_VALUES:
        report.append(f"| {N} | #{ranks_a[N]}/17 | {margins_a[N]:.4f} |")

    # MODE B: Extrapolated (N≤64 fit)
    report.append("\n### MODE B: Extrapolated (N≤64 fit applied to all N)\n")
    by_nf_b = defaultdict(list)
    for feat in all_feats:
        N = feat["N"]
        cN = power_law(N, popt_c_small[0], popt_c_small[1])
        wN = power_law(N, popt_w_small[0], popt_w_small[1])
        F = (ALPHA * (feat["d_eff"] - 4.0)**2
             + BETA * (feat["c1_c0"] - cN)**2
             + GAMMA * (feat["width_ratio"] - wN)**2)
        fc = dict(feat)
        fc["F"] = F
        by_nf_b[(feat["N"], feat["family"])].append(fc)

    ranks_b, frac_b, margins_b = evaluate(by_nf_b, N_VALUES)
    report.append(f"**Beat rate**: {frac_b*100:.1f}%\n")
    report.append("| N | Lor4D Rank | Margin to runner-up |")
    report.append("|---|:----------:|:-------------------:|")
    for N in N_VALUES:
        report.append(f"| {N} | #{ranks_b[N]}/17 | {margins_b[N]:.4f} |")

    # === Phase 6: Grid search over weights (full N range) ===
    report.append("\n\n## 3. Grid Search: Optimal Weights at Large N\n")
    report.append("Do optimal weights shift when large-N data is included?\n")

    best_grid = None
    best_grid_score = -1

    for a, b, g in product(
        [0.5, 1.0, 2.0, 5.0],
        [0.5, 1.0, 2.0, 5.0, 10.0],
        [1.0, 2.0, 5.0, 10.0, 20.0],
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
        ranks, frac, _ = evaluate(by_nf, N_VALUES)
        if frac > best_grid_score:
            best_grid_score = frac
            best_grid = (a, b, g, ranks)

    g_a, g_b, g_g, g_ranks = best_grid
    report.append(f"**Best weights (all N)**: α={g_a}, β={g_b}, γ={g_g}")
    report.append(f"**Beat rate**: {best_grid_score*100:.1f}%\n")
    report.append("| N | Lor4D Rank |")
    report.append("|---|:----------:|")
    for N in N_VALUES:
        report.append(f"| {N} | #{g_ranks[N]}/17 |")

    # === Phase 7: Full ranking at largest N for each mode ===
    for test_N in [128, 256]:
        report.append(f"\n\n## 4. Full Rankings at N={test_N}\n")
        report.append(f"### MODE A (Oracle, α={ALPHA} β={BETA} γ={GAMMA})\n")
        means = {}
        for fam in FAMILIES:
            rows = by_nf_a.get((test_N, fam), [])
            if rows:
                means[fam] = np.mean([r["F"] for r in rows])
        if means:
            ranked = sorted(means, key=means.get)
            report.append("| Rank | Family | Category | F |")
            report.append("|------|--------|----------|:-:|")
            for rank, fam in enumerate(ranked, 1):
                cat = CATEGORY[fam]
                tag = " ◆" if cat == "Lorentzian" else (" ●" if cat == "KR-family" else "")
                report.append(f"| {rank} | {fam}{tag} | {cat} | {means[fam]:.6f} |")

        report.append(f"\n### MODE B (Extrapolated, α={ALPHA} β={BETA} γ={GAMMA})\n")
        means = {}
        for fam in FAMILIES:
            rows = by_nf_b.get((test_N, fam), [])
            if rows:
                means[fam] = np.mean([r["F"] for r in rows])
        if means:
            ranked = sorted(means, key=means.get)
            report.append("| Rank | Family | Category | F |")
            report.append("|------|--------|----------|:-:|")
            for rank, fam in enumerate(ranked, 1):
                cat = CATEGORY[fam]
                tag = " ◆" if cat == "Lorentzian" else (" ●" if cat == "KR-family" else "")
                report.append(f"| {rank} | {fam}{tag} | {cat} | {means[fam]:.6f} |")

    # === Phase 8: Margin scaling analysis ===
    report.append("\n\n## 5. Margin Scaling: How Does Discrimination Improve with N?\n")
    report.append("| N | Margin A (Oracle) | Margin B (Extrap) | Margin ratio B/A |")
    report.append("|---|:-:|:-:|:-:|")
    for N in N_VALUES:
        ma = margins_a.get(N, 0)
        mb = margins_b.get(N, 0)
        ratio = mb / ma if ma > 0 else float("inf")
        report.append(f"| {N} | {ma:.6f} | {mb:.6f} | {ratio:.3f} |")

    # === Phase 9: Per-family mean feature at N=256 ===
    report.append("\n\n## 6. Feature Space at N=256\n")
    report.append("| Family | d_eff | c₁/c₀ | width | F (Oracle) |")
    report.append("|--------|:-:|:-:|:-:|:-:|")
    for fam in sorted(FAMILIES.keys()):
        rows_a = by_nf_a.get((256, fam), [])
        if rows_a:
            d = np.mean([r["d_eff"] for r in rows_a])
            c = np.mean([r["c1_c0"] for r in rows_a])
            w = np.mean([r["width_ratio"] for r in rows_a])
            f_val = np.mean([r["F"] for r in rows_a])
            report.append(f"| {fam} | {d:.3f} | {c:.4f} | {w:.4f} | {f_val:.6f} |")

    # === Summary ===
    report.append("\n\n## 7. Summary\n")
    all_a_one = all(r == 1 for r in ranks_a.values())
    all_b_one = all(r == 1 for r in ranks_b.values())
    report.append(f"- MODE A (Oracle): Lor4D #1 at ALL N? **{'YES' if all_a_one else 'NO'}**")
    report.append(f"- MODE B (Extrap): Lor4D #1 at ALL N? **{'YES' if all_b_one else 'NO'}**")
    report.append(f"- Grid search best: α={g_a}, β={g_b}, γ={g_g} (was α=0.5, β=1.0, γ=5.0)")
    report.append(f"- Finite-size scaling prediction error at N=256: "
                  f"Δc={abs(power_law(256, popt_c_small[0], popt_c_small[1]) - lor4d_centroids[256]['c1_c0']):.4f}, "
                  f"Δw={abs(power_law(256, popt_w_small[0], popt_w_small[1]) - lor4d_centroids[256]['width_ratio']):.4f}")

    # Write report
    out_dir = Path("outputs_carlip")
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "lsd_well_large_n_test.md"
    out_path.write_text("\n".join(report), encoding="utf-8")
    print(f"\nReport saved to {out_path}")

    # Console summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"MODE A (Oracle): Lor4D #1 at ALL N? {'YES' if all_a_one else 'NO'}")
    for N in N_VALUES:
        print(f"  N={N:4d}: rank #{ranks_a[N]}/17, margin {margins_a[N]:.6f}")
    print(f"\nMODE B (Extrap): Lor4D #1 at ALL N? {'YES' if all_b_one else 'NO'}")
    for N in N_VALUES:
        print(f"  N={N:4d}: rank #{ranks_b[N]}/17, margin {margins_b[N]:.6f}")
    print(f"\nGrid search best: α={g_a}, β={g_b}, γ={g_g}")
    print(f"Extrapolation (N≤64→256): Δc={abs(power_law(256, popt_c_small[0], popt_c_small[1]) - lor4d_centroids[256]['c1_c0']):.4f}, "
          f"Δw={abs(power_law(256, popt_w_small[0], popt_w_small[1]) - lor4d_centroids[256]['width_ratio']):.4f}")


if __name__ == "__main__":
    main()
