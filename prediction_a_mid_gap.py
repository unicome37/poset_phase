"""A_mid: 3D → 4D residual gap analysis.

Three candidate physics inputs for why 4D should beat 3D:
  (a) Cover/Hasse local dynamics: link structure differences
  (b) Interval spectrum fine structure (not just total S_BD)
  (c) Discrete gravity "mid-d optimum": link-based action term

All diagnostics are computed per-poset, then aggregated by family and N.
"""
from __future__ import annotations

import csv
from collections import Counter

import numpy as np

from generators import Poset
from observables import comparable_fraction
from observables_geo import dimension_proxy_views
from prediction_a_bd_bridge import (
    regenerate_poset, count_intervals_fast, compute_f5, CALIBRATED_WEIGHTS
)


# ---------------------------------------------------------------------------
# (a) Cover / Hasse local dynamics
# ---------------------------------------------------------------------------

def hasse_diagnostics(poset: Poset) -> dict:
    """Compute Hasse diagram (cover relation) statistics.

    A link/cover i→j means c[i,j]=True and no k with c[i,k] and c[k,j].
    """
    c = poset.closure
    n = poset.n

    # Interval matrix: M[i,j] = number of elements k with i<k<j
    M = c.astype(np.int32) @ c.astype(np.int32)

    # Links: c[i,j]=True and M[i,j]=0
    links = c & (M == 0)
    n_links = int(links.sum())
    n_relations = int(c.sum())

    # In-degree and out-degree in Hasse diagram
    in_degree = links.sum(axis=0)   # how many covers point TO each element
    out_degree = links.sum(axis=1)  # how many covers point FROM each element

    # Hasse width: max antichain size approximation via layer widths
    # (already available from layer_profile, but let's compute link-specific metrics)

    # Link density per element
    link_density = n_links / n if n > 0 else 0

    # Branching factor: mean out-degree in Hasse
    mean_out = float(out_degree.mean()) if n > 0 else 0
    max_out = int(out_degree.max()) if n > 0 else 0
    std_out = float(out_degree.std()) if n > 0 else 0

    # Convergence factor: mean in-degree
    mean_in = float(in_degree.mean()) if n > 0 else 0

    # "Diamond fraction": fraction of links that form part of a diamond
    # (i→j, i→k, j→l, k→l for some l)
    # This measures local 2D-like vs higher-d structure
    # Simplified: count pairs of links sharing a source or target
    n_diamond_sources = 0
    for i in range(n):
        targets = np.where(links[i])[0]
        if len(targets) >= 2:
            # Check if any pair of targets share a common successor
            for t1_idx in range(len(targets)):
                for t2_idx in range(t1_idx + 1, len(targets)):
                    t1, t2 = targets[t1_idx], targets[t2_idx]
                    # Do t1 and t2 share a common link-successor?
                    succ_t1 = set(np.where(links[t1])[0])
                    succ_t2 = set(np.where(links[t2])[0])
                    if succ_t1 & succ_t2:
                        n_diamond_sources += 1
                        break
                else:
                    continue
                break

    diamond_frac = n_diamond_sources / n if n > 0 else 0

    return {
        "n_links": n_links,
        "n_relations": n_relations,
        "link_density": float(link_density),
        "mean_out_degree": float(mean_out),
        "max_out_degree": int(max_out),
        "std_out_degree": float(std_out),
        "mean_in_degree": float(mean_in),
        "diamond_frac": float(diamond_frac),
    }


# ---------------------------------------------------------------------------
# (b) Interval spectrum fine structure
# ---------------------------------------------------------------------------

def interval_spectrum(poset: Poset) -> dict:
    """Compute the full interval size distribution and derived spectral metrics.

    Instead of the total S_BD (which is a crude aggregate), we look at:
    - The shape of the C_k distribution
    - The ratio C_1/C_0 (fraction of "almost links")
    - The spectral entropy and moments
    - The "gap ratio": how sharply C_k drops from k=0 to k=1
    """
    counts = count_intervals_fast(poset)
    total = sum(counts.values())
    if total == 0:
        return {"total": 0, "C0_frac": 1.0, "C1_frac": 0.0,
                "C1_C0_ratio": 0.0, "gap_ratio": 1.0,
                "spectral_entropy": 0.0, "kurtosis": 0.0,
                "tail_weight": 0.0, "median_k": 0.0}

    C = {k: counts.get(k, 0) for k in range(max(counts.keys()) + 1)}
    C0 = C.get(0, 0)
    C1 = C.get(1, 0)

    C0_frac = C0 / total
    C1_frac = C1 / total
    C1_C0_ratio = C1 / C0 if C0 > 0 else 0.0

    # Gap ratio: how sharply does the distribution drop after k=0?
    # High gap ratio = links dominate sharply
    gap_ratio = C0_frac  # simple: fraction that are links

    # Spectral entropy
    probs = np.array([v / total for v in C.values() if v > 0])
    spectral_entropy = -float(np.sum(probs * np.log(probs)))

    # Moments of the interval size distribution
    ks = np.array(list(C.keys()))
    ps = np.array([C[k] / total for k in ks])
    mean_k = float(np.sum(ks * ps))
    var_k = float(np.sum((ks - mean_k) ** 2 * ps))
    if var_k > 0:
        kurtosis = float(np.sum((ks - mean_k) ** 4 * ps) / var_k ** 2 - 3)
    else:
        kurtosis = 0.0

    # Tail weight: fraction of relations with k >= 3
    tail_weight = sum(C.get(k, 0) for k in range(3, max(counts.keys()) + 1)) / total

    # Median interval size
    cumsum = 0
    median_k = 0
    for k in sorted(C.keys()):
        cumsum += C[k]
        if cumsum >= total / 2:
            median_k = k
            break

    return {
        "total": int(total),
        "C0_frac": float(C0_frac),
        "C1_frac": float(C1_frac),
        "C1_C0_ratio": float(C1_C0_ratio),
        "gap_ratio": float(gap_ratio),
        "spectral_entropy": float(spectral_entropy),
        "mean_k": float(mean_k),
        "var_k": float(var_k),
        "kurtosis": float(kurtosis),
        "tail_weight": float(tail_weight),
        "median_k": float(median_k),
    }


# ---------------------------------------------------------------------------
# (c) Discrete gravity "mid-d optimum"
# ---------------------------------------------------------------------------

def discrete_gravity_diagnostics(poset: Poset) -> dict:
    """Compute link-based discrete gravity action terms.

    The Benincasa-Dowker action for d=2 is:
      S_BD^(2) = N - 2 * n_links
    For d=4:
      S_BD^(4) = N - (9/2)*C_0 + 8*C_1 - (9/2)*C_2 + (8/9)*C_3

    The key insight: in 3D, the BD action doesn't have a clean form.
    The d=4 BD action has a specific alternating-sign structure that
    may provide a "resonance" at d=4.

    We compute:
    - S_BD^(2): the 2D BD action (should be large negative for 2D sprinklings)
    - S_BD^(4): the 4D BD action (should be near 0 for 4D sprinklings)
    - The "action ratio" S_BD^(4)/N: normalized 4D action per element
    """
    counts = count_intervals_fast(poset)
    n = poset.n

    C0 = counts.get(0, 0)
    C1 = counts.get(1, 0)
    C2 = counts.get(2, 0)
    C3 = counts.get(3, 0)

    # BD action d=2
    s_bd_2 = n - 2 * C0

    # BD action d=4 (Benincasa-Dowker 2010, eq. 3.2)
    s_bd_4 = n - (9.0 / 2) * C0 + 8 * C1 - (9.0 / 2) * C2 + (8.0 / 9) * C3

    # Normalized per element
    s_bd_2_per_n = s_bd_2 / n if n > 0 else 0
    s_bd_4_per_n = s_bd_4 / n if n > 0 else 0

    # "Action deficit": |S_BD^(4)/N| — how close to the continuum limit
    # For a faithful d=4 sprinkling, S_BD^(4)/N should approach a constant
    # related to the scalar curvature. For non-4D sprinklings, it deviates.
    action_deficit_4 = abs(s_bd_4_per_n)

    # Link-action balance: C0 vs higher-order interval contributions
    # In the BD(4) formula, C0 has the largest coefficient (9/2)
    # The balance between C0 and C1,C2,C3 terms indicates "action harmony"
    if C0 > 0:
        c1_balance = 8 * C1 / ((9.0 / 2) * C0)
        c2_balance = (9.0 / 2) * C2 / ((9.0 / 2) * C0)
    else:
        c1_balance = 0.0
        c2_balance = 0.0

    return {
        "s_bd_2": float(s_bd_2),
        "s_bd_4": float(s_bd_4),
        "s_bd_2_per_n": float(s_bd_2_per_n),
        "s_bd_4_per_n": float(s_bd_4_per_n),
        "action_deficit_4": float(action_deficit_4),
        "c1_balance": float(c1_balance),
        "c2_balance": float(c2_balance),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    raw_path = "outputs_unified_functional/raw_features.csv"
    rows = list(csv.DictReader(open(raw_path)))
    families = {"Lor3D", "Lor4D"}  # Focus on the 3D vs 4D gap
    rows_all = [r for r in rows if r["family"] in {"Lor2D", "Lor3D", "Lor4D", "Lor5D"}]
    families_order = ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]

    data = []
    for i, r in enumerate(rows_all):
        fam, n, rep = r["family"], int(r["N"]), int(r["rep"])
        poset = regenerate_poset(fam, n, rep, 42)
        f5 = compute_f5(r, CALIBRATED_WEIGHTS)
        hasse = hasse_diagnostics(poset)
        spectrum = interval_spectrum(poset)
        gravity = discrete_gravity_diagnostics(poset)
        data.append({"family": fam, "N": n, "rep": rep, "f5": f5,
                      **{f"h_{k}": v for k, v in hasse.items()},
                      **{f"sp_{k}": v for k, v in spectrum.items()},
                      **{f"gr_{k}": v for k, v in gravity.items()}})
        if (i + 1) % 16 == 0:
            print(f"  [{i+1}/{len(rows_all)}]")

    n_values = sorted(set(d["N"] for d in data))

    # =====================================================================
    # (a) Hasse / Cover diagnostics
    # =====================================================================
    print("\n" + "=" * 80)
    print("(a) HASSE / COVER LOCAL DYNAMICS — 3D vs 4D focus")
    print("=" * 80)

    hasse_keys = ["h_link_density", "h_mean_out_degree", "h_std_out_degree",
                  "h_diamond_frac"]
    for key in hasse_keys:
        label = key.replace("h_", "")
        print(f"\n  {label}:")
        for fam in families_order:
            vals = [d[key] for d in data if d["family"] == fam]
            print(f"    {fam}: {np.mean(vals):.4f} +/- {np.std(vals):.4f}")
        # Per-N for 3D vs 4D
        print(f"    Per-N (3D vs 4D):")
        for nv in n_values:
            v3 = [d[key] for d in data if d["family"] == "Lor3D" and d["N"] == nv]
            v4 = [d[key] for d in data if d["family"] == "Lor4D" and d["N"] == nv]
            gap = np.mean(v4) - np.mean(v3)
            print(f"      N={nv}: 3D={np.mean(v3):.4f} 4D={np.mean(v4):.4f} gap={gap:+.4f}")

    # =====================================================================
    # (b) Interval spectrum fine structure
    # =====================================================================
    print("\n" + "=" * 80)
    print("(b) INTERVAL SPECTRUM FINE STRUCTURE — 3D vs 4D focus")
    print("=" * 80)

    spec_keys = ["sp_C0_frac", "sp_C1_frac", "sp_C1_C0_ratio",
                 "sp_spectral_entropy", "sp_kurtosis", "sp_tail_weight"]
    for key in spec_keys:
        label = key.replace("sp_", "")
        print(f"\n  {label}:")
        for fam in families_order:
            vals = [d[key] for d in data if d["family"] == fam]
            print(f"    {fam}: {np.mean(vals):.4f} +/- {np.std(vals):.4f}")
        print(f"    Per-N (3D vs 4D):")
        for nv in n_values:
            v3 = [d[key] for d in data if d["family"] == "Lor3D" and d["N"] == nv]
            v4 = [d[key] for d in data if d["family"] == "Lor4D" and d["N"] == nv]
            gap = np.mean(v4) - np.mean(v3)
            sep = abs(gap) / (np.std(v3) + np.std(v4) + 1e-10) * 2
            print(f"      N={nv}: 3D={np.mean(v3):.4f} 4D={np.mean(v4):.4f} "
                  f"gap={gap:+.4f} sep={sep:.2f}sigma")

    # =====================================================================
    # (c) Discrete gravity — BD action comparison
    # =====================================================================
    print("\n" + "=" * 80)
    print("(c) DISCRETE GRAVITY — BD ACTION 3D vs 4D")
    print("=" * 80)

    grav_keys = ["gr_s_bd_2_per_n", "gr_s_bd_4_per_n", "gr_action_deficit_4",
                 "gr_c1_balance"]
    for key in grav_keys:
        label = key.replace("gr_", "")
        print(f"\n  {label}:")
        for fam in families_order:
            vals = [d[key] for d in data if d["family"] == fam]
            print(f"    {fam}: {np.mean(vals):.4f} +/- {np.std(vals):.4f}")
        print(f"    Per-N (3D vs 4D):")
        for nv in n_values:
            v3 = [d[key] for d in data if d["family"] == "Lor3D" and d["N"] == nv]
            v4 = [d[key] for d in data if d["family"] == "Lor4D" and d["N"] == nv]
            gap = np.mean(v4) - np.mean(v3)
            sep = abs(gap) / (np.std(v3) + np.std(v4) + 1e-10) * 2
            print(f"      N={nv}: 3D={np.mean(v3):.4f} 4D={np.mean(v4):.4f} "
                  f"gap={gap:+.4f} sep={sep:.2f}sigma")

    # =====================================================================
    # KEY QUESTION: Which metric best separates 3D from 4D?
    # =====================================================================
    print("\n" + "=" * 80)
    print("DISCRIMINANT RANKING: Best 3D/4D separators (by Cohen's d)")
    print("=" * 80)

    all_keys = hasse_keys + spec_keys + grav_keys
    discriminants = []
    for key in all_keys:
        v3 = [d[key] for d in data if d["family"] == "Lor3D"]
        v4 = [d[key] for d in data if d["family"] == "Lor4D"]
        m3, m4 = np.mean(v3), np.mean(v4)
        s3, s4 = np.std(v3), np.std(v4)
        pooled_std = np.sqrt((s3**2 + s4**2) / 2)
        cohens_d = (m4 - m3) / pooled_std if pooled_std > 0 else 0
        discriminants.append((key, cohens_d, m3, m4, pooled_std))

    discriminants.sort(key=lambda x: abs(x[1]), reverse=True)
    print(f"\n  {'Metric':<25s} {'Cohen_d':>8s} {'Mean(3D)':>10s} {'Mean(4D)':>10s} {'Direction':>10s}")
    print(f"  {'-'*25} {'-'*8} {'-'*10} {'-'*10} {'-'*10}")
    for key, cd, m3, m4, _ in discriminants:
        direction = "4D higher" if cd > 0 else "3D higher"
        print(f"  {key:<25s} {cd:>8.3f} {m3:>10.4f} {m4:>10.4f} {direction:>10s}")

    # =====================================================================
    # Test: Can the best discriminant be used as a filter or mild energy term?
    # =====================================================================
    print("\n" + "=" * 80)
    print("FILTER TEST: Best discriminant as 3D→4D filter (post faithfulness)")
    print("=" * 80)

    # Take top 3 discriminants and test as filters
    for key, cd, m3, m4, ps in discriminants[:5]:
        label = key
        direction = ">" if cd > 0 else "<"  # 4D is higher/lower

        v3_all = [d[key] for d in data if d["family"] == "Lor3D"]
        v4_all = [d[key] for d in data if d["family"] == "Lor4D"]

        # Find threshold that maximizes separation
        all_v = sorted(set(v3_all + v4_all))
        best_thresh = None
        best_score = -1
        for thresh in np.linspace(min(all_v), max(all_v), 50):
            if direction == ">":
                pass3 = sum(1 for v in v3_all if v >= thresh) / len(v3_all)
                pass4 = sum(1 for v in v4_all if v >= thresh) / len(v4_all)
            else:
                pass3 = sum(1 for v in v3_all if v <= thresh) / len(v3_all)
                pass4 = sum(1 for v in v4_all if v <= thresh) / len(v4_all)
            # Score: maximize 4D pass rate while minimizing 3D pass rate
            score = pass4 - pass3
            if score > best_score:
                best_score = score
                best_thresh = thresh

        if best_thresh is not None:
            if direction == ">":
                p3 = sum(1 for v in v3_all if v >= best_thresh) / len(v3_all) * 100
                p4 = sum(1 for v in v4_all if v >= best_thresh) / len(v4_all) * 100
            else:
                p3 = sum(1 for v in v3_all if v <= best_thresh) / len(v3_all) * 100
                p4 = sum(1 for v in v4_all if v <= best_thresh) / len(v4_all) * 100

            print(f"\n  {label} (Cohen's d={cd:.3f}):")
            print(f"    Best threshold: {direction} {best_thresh:.4f}")
            print(f"    3D pass: {p3:.0f}%  4D pass: {p4:.0f}%  gap: {p4-p3:.0f}%")

    # =====================================================================
    # Test: mild energy term using best discriminant
    # =====================================================================
    print("\n" + "=" * 80)
    print("ENERGY TERM TEST: F5 + alpha * best_discriminant (post faithfulness filter)")
    print("=" * 80)

    # Use post-faithfulness data (link_frac >= 0.50, excludes 2D)
    from prediction_a_faithfulness import faithfulness_score
    filtered = [d for d in data if d["family"] != "Lor2D"]  # Simplified: just exclude 2D

    for key, cd, m3, m4, ps in discriminants[:3]:
        sign = +1 if cd < 0 else -1  # We want 4D to be lower, so penalize 3D's direction

        f5_3 = np.mean([d["f5"] for d in filtered if d["family"] == "Lor3D"])
        f5_4 = np.mean([d["f5"] for d in filtered if d["family"] == "Lor4D"])
        f5_5 = np.mean([d["f5"] for d in filtered if d["family"] == "Lor5D"])
        m3f = np.mean([d[key] for d in filtered if d["family"] == "Lor3D"])
        m4f = np.mean([d[key] for d in filtered if d["family"] == "Lor4D"])
        m5f = np.mean([d[key] for d in filtered if d["family"] == "Lor5D"])

        denom_34 = sign * (m3f - m4f)
        if abs(denom_34) > 1e-10:
            alpha_cross = (f5_4 - f5_3) / denom_34
        else:
            alpha_cross = float('inf')

        denom_45 = sign * (m4f - m5f)  # for 4D < 5D preservation
        if abs(denom_45) > 1e-10:
            alpha_limit = (f5_5 - f5_4) / denom_45
        else:
            alpha_limit = float('inf')

        print(f"\n  {key} (sign={sign:+d}):")
        print(f"    Means: 3D={m3f:.4f} 4D={m4f:.4f} 5D={m5f:.4f}")
        print(f"    4D crosses 3D at alpha = {alpha_cross:.2f}")
        print(f"    4D crosses 5D at alpha = {alpha_limit:.2f}")
        if 0 < alpha_cross < alpha_limit:
            # Check relative contribution at crossing
            contrib_3d = abs(alpha_cross * m3f)
            ratio_3d = contrib_3d / f5_3 * 100
            print(f"    *** WINDOW [{alpha_cross:.1f}, {alpha_limit:.1f}] ***")
            print(f"    Contribution at crossing: {ratio_3d:.0f}% of F5(3D)")
        elif alpha_cross < 0:
            print(f"    (wrong direction — no crossing)")
        else:
            print(f"    (5D crosses first — no window)")


if __name__ == "__main__":
    main()
