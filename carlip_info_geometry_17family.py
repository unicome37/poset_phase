"""
Post-Carlip Line 2: Information Geometry of 17 Families
=========================================================
Instead of constructing a single functional F[X] and ranking,
we map each poset into a multi-dimensional invariant space:

  Φ(X) = (ρ, d_eff, Σ_hist, R, link_frac, C₁/C₀, longest_chain/N, width/N, ...)

Then ask:
  1. Does Lor4D occupy a DISTINCT region in Φ-space?
  2. Is Lor4D separable from KR_2layer (the fatal counterexample)?
  3. Which features carry the most discriminative power?
  4. Is there a LOW-DIMENSIONAL manifold on which Lorentzian families lie?

No logH. No F7. Pure causal-geometric observables only.
"""
from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from bd_action import count_intervals_fast, bdg_action_d2_link, bdg_action_d4_standard
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
    """Longest chain (antichain-partition lower bound on height)."""
    n = poset.n
    if n <= 1:
        return n
    c = poset.closure
    dp = [1] * n
    order = list(range(n))
    for i in order:
        for j in range(i + 1, n):
            if c[i, j]:
                dp[j] = max(dp[j], dp[i] + 1)
    return max(dp)


def max_antichain_width(poset: Poset) -> int:
    """Size of largest antichain (greedy approximation via layer decomposition)."""
    n = poset.n
    if n <= 1:
        return n
    c = poset.closure
    assigned = [False] * n
    max_w = 0
    remaining = set(range(n))
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


def compute_features(poset: Poset, N: int) -> dict:
    counts = count_intervals_fast(poset, k_max=5)
    C0 = counts.get(0)
    C1 = counts.get(1)
    C2 = counts.get(2)
    C3 = counts.get(3)
    total_rel = counts.total_relations
    n_pairs = N * (N - 1) // 2

    # Basic density
    f2 = total_rel / max(1, n_pairs)  # comparability fraction (= our ρ)
    link_frac = C0 / max(1, n_pairs)  # link fraction
    R = 1.0 - float(C0) / float(total_rel) if total_rel > 0 else 0.0

    # Interval shape ratios
    c1_over_c0 = C1 / max(1, C0)
    c2_over_c0 = C2 / max(1, C0)
    c3_over_c0 = C3 / max(1, C0)

    # Effective dimension (Myrheim-Meyer)
    xi_dim_val, d_eff = compute_xi_dim(poset)

    # Sigma_hist
    sigma_hist = compute_sigma_hist(poset)

    # Chain and antichain
    lc = longest_chain_length(poset)
    aw = max_antichain_width(poset)
    chain_ratio = lc / max(1, N)
    width_ratio = aw / max(1, N)

    # BDG actions (normalized)
    s_link = bdg_action_d2_link(counts, N, normalized=True)
    s_d4 = bdg_action_d4_standard(counts, N, normalized=True)

    # Layer count proxy (from sigma_hist computation internals)
    # We approximate: n_layers ≈ longest chain
    n_layers = lc

    return {
        "f2": f2,
        "link_frac": link_frac,
        "R": R,
        "d_eff": d_eff,
        "sigma_hist": sigma_hist,
        "c1_c0": c1_over_c0,
        "c2_c0": c2_over_c0,
        "c3_c0": c3_over_c0,
        "chain_ratio": chain_ratio,
        "width_ratio": width_ratio,
        "s_link_n": s_link,
        "s_d4_n": s_d4,
        "n_layers": n_layers,
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


FEATURE_NAMES = [
    "f2", "link_frac", "R", "d_eff", "sigma_hist",
    "c1_c0", "c2_c0", "c3_c0",
    "chain_ratio", "width_ratio",
    "s_link_n", "s_d4_n",
]


def main():
    N_VALUES = [20, 36, 48]
    REPS = 15
    SEED_BASE = 77

    print("=" * 80)
    print("POST-CARLIP LINE 2: Information Geometry of 17 Families")
    print("=" * 80)

    all_rows = []
    total = len(FAMILIES) * len(N_VALUES) * REPS
    done = 0

    for fam_name, gen_fn in FAMILIES.items():
        for N in N_VALUES:
            for rep in range(REPS):
                seed = SEED_BASE + hash(fam_name) % 10000 + N * 100 + rep
                seed = seed % (2**31)
                try:
                    poset = gen_fn(N, seed=seed)
                    feat = compute_features(poset, N)
                    feat["family"] = fam_name
                    feat["category"] = CATEGORY[fam_name]
                    feat["N"] = N
                    feat["rep"] = rep
                    all_rows.append(feat)
                except Exception as e:
                    print(f"  ERROR: {fam_name} N={N} rep={rep}: {e}")
                done += 1
                if done % 100 == 0:
                    print(f"  [{done}/{total}]", flush=True)

    by_nf = defaultdict(list)
    for r in all_rows:
        by_nf[(r["N"], r["family"])].append(r)

    report = []
    report.append("# Post-Carlip Line 2: Information Geometry of 17 Families\n")
    report.append("**No logH. No F7. Pure causal-geometric observables.**\n")
    report.append(f"N = {N_VALUES}, reps = {REPS}\n")
    report.append("Features: f₂, link_frac, R, d_eff, Σ_hist, C₁/C₀, C₂/C₀, C₃/C₀,")
    report.append("chain_ratio, width_ratio, S_link/N, S_d4/N\n")

    # === Section 1: Feature profiles ===
    report.append("\n## 1. Family Feature Profiles\n")
    for N in N_VALUES:
        report.append(f"\n### N = {N}\n")
        header = "| Family | Cat | f₂ | link | R | d_eff | Σ_hist | C₁/C₀ | chain | width |"
        sep =    "|--------|-----|-----|------|---|-------|--------|-------|-------|-------|"
        report.append(header)
        report.append(sep)

        for fam in sorted(FAMILIES.keys()):
            rows = by_nf.get((N, fam), [])
            if not rows:
                continue
            cat = CATEGORY[fam][:3]
            f2 = np.mean([r["f2"] for r in rows])
            lf = np.mean([r["link_frac"] for r in rows])
            R = np.mean([r["R"] for r in rows])
            de = np.mean([r["d_eff"] for r in rows])
            sh = np.mean([r["sigma_hist"] for r in rows])
            c1 = np.mean([r["c1_c0"] for r in rows])
            cr = np.mean([r["chain_ratio"] for r in rows])
            wr = np.mean([r["width_ratio"] for r in rows])
            tag = "◆" if CATEGORY[fam] == "Lorentzian" else ""
            report.append(
                f"| {fam}{tag} | {cat} | {f2:.3f} | {lf:.3f} | {R:.3f} | "
                f"{de:.2f} | {sh:.3f} | {c1:.3f} | {cr:.3f} | {wr:.3f} |"
            )

    # === Section 2: Pairwise distances in feature space ===
    report.append("\n## 2. Lor4D vs Key Competitors (Euclidean Distance in Standardized Feature Space)\n")
    for N in N_VALUES:
        report.append(f"\n### N = {N}\n")

        # Compute standardized features
        all_vecs = {}
        for fam in FAMILIES:
            rows = by_nf.get((N, fam), [])
            if rows:
                vec = np.array([np.mean([r[f] for r in rows]) for f in FEATURE_NAMES])
                all_vecs[fam] = vec

        if not all_vecs:
            continue

        # Standardize
        all_mat = np.array(list(all_vecs.values()))
        mu = all_mat.mean(axis=0)
        sigma = all_mat.std(axis=0)
        sigma[sigma < 1e-10] = 1.0

        std_vecs = {f: (v - mu) / sigma for f, v in all_vecs.items()}

        lor4d = std_vecs.get("Lor4D")
        if lor4d is None:
            continue

        dists = {}
        for fam, vec in std_vecs.items():
            if fam != "Lor4D":
                dists[fam] = np.linalg.norm(lor4d - vec)

        report.append("| Family | Category | Distance to Lor4D | Nearest? |")
        report.append("|--------|----------|:-----------------:|:--------:|")
        for fam, d in sorted(dists.items(), key=lambda x: x[1]):
            cat = CATEGORY[fam]
            near = "★" if d < 1.5 else ""
            report.append(f"| {fam} | {cat} | {d:.3f} | {near} |")

    # === Section 3: Which features separate Lor4D from KR_2layer? ===
    report.append("\n## 3. Lor4D vs KR_2layer: Feature-by-Feature Comparison\n")
    report.append("**KR_2layer is the fatal counterexample for all single-functional approaches.**\n")

    for N in N_VALUES:
        report.append(f"\n### N = {N}\n")
        lor4d_rows = by_nf.get((N, "Lor4D"), [])
        kr2_rows = by_nf.get((N, "KR_2layer"), [])
        if not lor4d_rows or not kr2_rows:
            continue

        report.append("| Feature | Lor4D | KR_2layer | Δ | |Δ|/σ | Separable? |")
        report.append("|---------|:-----:|:---------:|:-:|:----:|:----------:|")

        for feat in FEATURE_NAMES:
            l4 = np.array([r[feat] for r in lor4d_rows])
            k2 = np.array([r[feat] for r in kr2_rows])
            mean_l4 = np.mean(l4)
            mean_k2 = np.mean(k2)
            delta = mean_l4 - mean_k2
            pooled_std = np.sqrt((np.var(l4) + np.var(k2)) / 2)
            effect_size = abs(delta) / pooled_std if pooled_std > 1e-10 else 0.0
            sep = "✅" if effect_size > 2.0 else ("⚠️" if effect_size > 1.0 else "❌")
            report.append(
                f"| {feat} | {mean_l4:.4f} | {mean_k2:.4f} | {delta:+.4f} | "
                f"{effect_size:.2f} | {sep} |"
            )

    # === Section 4: Lor4D vs ALL non-Lor: which features are unique? ===
    report.append("\n## 4. Lor4D Uniqueness: Features Where Lor4D is an Outlier\n")
    report.append("For each feature, compute Lor4D's z-score relative to all non-Lorentzian families.\n")

    for N in N_VALUES:
        report.append(f"\n### N = {N}\n")
        lor4d_rows = by_nf.get((N, "Lor4D"), [])
        if not lor4d_rows:
            continue

        non_lor_rows = []
        for fam in FAMILIES:
            if CATEGORY[fam] != "Lorentzian":
                non_lor_rows.extend(by_nf.get((N, fam), []))

        if not non_lor_rows:
            continue

        report.append("| Feature | Lor4D mean | non-Lor mean | non-Lor σ | z-score | Outlier? |")
        report.append("|---------|:----------:|:------------:|:---------:|:-------:|:--------:|")

        for feat in FEATURE_NAMES:
            l4 = np.mean([r[feat] for r in lor4d_rows])
            nl = np.array([r[feat] for r in non_lor_rows])
            nl_mean = np.mean(nl)
            nl_std = np.std(nl)
            z = (l4 - nl_mean) / nl_std if nl_std > 1e-10 else 0.0
            outlier = "✅" if abs(z) > 2.0 else ("⚠️" if abs(z) > 1.0 else "")
            report.append(
                f"| {feat} | {l4:.4f} | {nl_mean:.4f} | {nl_std:.4f} | {z:+.2f} | {outlier} |"
            )

    # === Section 5: Cluster Analysis (simple) ===
    report.append("\n## 5. Category Separation via d_eff + R Plane\n")
    report.append("The (d_eff, R) plane is the most interpretable 2D projection.\n")

    for N in N_VALUES:
        report.append(f"\n### N = {N}\n")
        report.append("| Family | Category | d_eff | R | Quadrant |")
        report.append("|--------|----------|:-----:|:-:|:--------:|")

        for fam in sorted(FAMILIES.keys()):
            rows = by_nf.get((N, fam), [])
            if not rows:
                continue
            de = np.mean([r["d_eff"] for r in rows])
            R = np.mean([r["R"] for r in rows])
            # Quadrant: high-d_eff + low-R = "Lorentzian corner"
            if de > 3.5 and R < 0.2:
                quad = "Lor-corner ★"
            elif de > 3.5:
                quad = "high-d"
            elif R < 0.2:
                quad = "low-R"
            else:
                quad = "bulk"
            report.append(f"| {fam} | {CATEGORY[fam]} | {de:.2f} | {R:.3f} | {quad} |")

    # === Conclusion ===
    report.append("\n## 6. Conclusions\n")
    report.append("### Key Question: Does Lor4D occupy a unique position in invariant space?\n")
    report.append("### Answer: Check the d_eff × R plane and feature separation tables above.\n")
    report.append("### If Lor4D is in the 'Lor-corner' (high d_eff, low R) and KR_2layer is NOT,")
    report.append("### then information geometry succeeds where single functionals failed.\n")

    out = "\n".join(report)
    outdir = Path("outputs_carlip")
    outdir.mkdir(exist_ok=True)
    md_path = outdir / "info_geometry_17family.md"
    md_path.write_text(out, encoding="utf-8")
    print(f"\nSaved: {md_path}")
    print("\n" + out)


if __name__ == "__main__":
    main()
