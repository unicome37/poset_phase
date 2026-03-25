#!/usr/bin/env python3
"""Prediction C — Reverse Intervention: Edge Removal → Layer Merge → Entropy Increase.

Symmetry test: if adding an edge that creates a layer split LOWERS entropy,
then REMOVING an edge that causes a layer merge should RAISE entropy.

Protocol:
1. Generate Lor2D posets at small N.
2. For each poset, find "critical cover edges" — edges in the Hasse diagram
   whose removal reduces layer_count (i.e., merges two layers).
3. Remove each such edge (and recompute transitive reduction), measure Δlog_H.
4. Compare with removal of non-critical edges (that don't affect layer_count).
5. If Δlog_H > 0 for critical removals and Δlog_H ≈ 0 or smaller for
   non-critical removals, the bidirectional causal link is established.

Additionally tests: does the MAGNITUDE of entropy increase under merge
match the magnitude of decrease under split?
"""
from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from generators import (
    Poset, transitive_closure,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
)
from entropy_exact import count_linear_extensions_exact
from matched_residual_freedom_check import layer_index_by_minima

OUT_DIR = Path("outputs_exploratory/prediction_c_reverse_intervention")
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_VALUES = [14, 16]
POSETS_PER_N = 50
MAX_REMOVALS_PER_TYPE = 12
SEED_BASE = 20260318


def layer_count(poset: Poset) -> int:
    idx = layer_index_by_minima(poset.closure)
    return int(idx.max() + 1) if len(idx) > 0 else 0


def log_H(poset: Poset) -> float:
    return math.log(count_linear_extensions_exact(poset, prefer_c=False))


def hasse_edges(closure: np.ndarray) -> list[tuple[int, int]]:
    """Return the Hasse diagram edges (transitive reduction)."""
    n = closure.shape[0]
    edges = []
    for i in range(n):
        succ = np.where(closure[i])[0]
        for j in succ:
            # Check if there's an intermediate: any k with i<k<j?
            intermediates = np.where(closure[i] & closure[:, j])[0]
            if len(intermediates) == 0:
                edges.append((i, j))
    return edges


def remove_edge_and_rebuild(poset: Poset, a: int, b: int) -> Poset:
    """Remove the direct relation a < b and recompute transitive closure.
    
    We rebuild from the Hasse diagram minus (a,b), then re-close.
    """
    # Get all Hasse edges except (a,b)
    edges = hasse_edges(poset.closure)
    n = poset.n
    adj = np.zeros((n, n), dtype=bool)
    for (u, v) in edges:
        if u == a and v == b:
            continue
        adj[u, v] = True
    return Poset(transitive_closure(adj))


def run():
    rng = np.random.default_rng(SEED_BASE)
    rows = []
    t0 = time.time()

    for N in N_VALUES:
        print(f"\n=== N = {N} ===")
        for idx in range(POSETS_PER_N):
            seed = SEED_BASE + N * 1000 + idx
            poset = generate_lorentzian_like_2d(N, seed=seed)
            lc0 = layer_count(poset)
            lh0 = log_H(poset)

            if lc0 <= 2:
                continue  # need at least 3 layers to have merge potential

            edges = hasse_edges(poset.closure)
            if len(edges) == 0:
                continue

            # Classify edges by whether their removal reduces layer count
            critical = []      # removal merges layers
            non_critical = []  # removal doesn't change layer count

            # Sample a subset of edges to classify (to bound compute)
            sample_size = min(len(edges), MAX_REMOVALS_PER_TYPE * 3)
            chosen_indices = rng.choice(len(edges), size=sample_size, replace=False)

            for ei in chosen_indices:
                a, b = edges[ei]
                new_poset = remove_edge_and_rebuild(poset, a, b)
                lc1 = layer_count(new_poset)
                lh1 = log_H(new_poset)

                edge_type = "critical" if lc1 < lc0 else "non_critical"
                if lc1 < lc0:
                    critical.append((a, b, lc1, lh1))
                else:
                    non_critical.append((a, b, lc1, lh1))

            # Record all classified removals
            for group_name, group_list in [("critical", critical),
                                            ("non_critical", non_critical)]:
                for (a, b, lc1, lh1) in group_list:
                    rows.append({
                        "N": N,
                        "poset_seed": seed,
                        "poset_idx": idx,
                        "original_layer_count": lc0,
                        "original_log_H": lh0,
                        "edge_removed": f"{a}->{b}",
                        "removal_type": group_name,
                        "new_layer_count": lc1,
                        "new_log_H": lh1,
                        "delta_layer_count": lc1 - lc0,
                        "delta_log_H": lh1 - lh0,
                        "abs_delta_log_H": abs(lh1 - lh0),
                        "merge_occurred": int(lc1 < lc0),
                    })

            if (idx + 1) % 10 == 0:
                elapsed = time.time() - t0
                print(f"  N={N}: {idx+1}/{POSETS_PER_N} ({elapsed:.0f}s)")

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "reverse_intervention_raw.csv", index=False)
    print(f"\nTotal edge removals: {len(df)}")

    # ══════════════════════════════════════════════════════════
    #  ANALYSIS
    # ══════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("ANALYSIS: Reverse Intervention (Edge Removal)")
    print("=" * 70)

    # A. Critical vs Non-critical comparison
    print("\n--- A. Critical (MERGE) vs Non-critical edge removal ---")
    for N in N_VALUES:
        sub = df[df["N"] == N]
        crit = sub[sub["removal_type"] == "critical"]
        noncrit = sub[sub["removal_type"] == "non_critical"]
        print(f"\n  N={N}:")
        print(f"    Critical (merge):     n={len(crit)}, mean Δlog_H={crit['delta_log_H'].mean():.4f}")
        print(f"    Non-critical:         n={len(noncrit)}, mean Δlog_H={noncrit['delta_log_H'].mean():.4f}")
        print(f"    Critical merge rate:  {crit['merge_occurred'].mean():.3f}")

        if len(crit) >= 3 and len(noncrit) >= 3:
            t_stat, p_val = stats.ttest_ind(
                crit["delta_log_H"], noncrit["delta_log_H"],
                equal_var=False
            )
            n1, n2 = len(crit), len(noncrit)
            s1, s2 = crit["delta_log_H"].std(), noncrit["delta_log_H"].std()
            pooled_sd = math.sqrt(((n1-1)*s1**2 + (n2-1)*s2**2) / (n1+n2-2)) if (n1+n2-2) > 0 else 1
            d = (crit["delta_log_H"].mean() - noncrit["delta_log_H"].mean()) / pooled_sd if pooled_sd > 0 else float('nan')
            print(f"    Welch t={t_stat:.3f}, p={p_val:.6f}, Cohen's d={d:.3f}")

    # B. Direction check: does removing critical edges INCREASE entropy?
    print("\n--- B. Direction check: Δlog_H sign for each group ---")
    crit_all = df[df["removal_type"] == "critical"]
    noncrit_all = df[df["removal_type"] == "non_critical"]
    print(f"  Critical:    mean Δlog_H = {crit_all['delta_log_H'].mean():.4f} "
          f"({'↑ INCREASES' if crit_all['delta_log_H'].mean() > 0 else '↓ decreases'})")
    print(f"  Non-critical: mean Δlog_H = {noncrit_all['delta_log_H'].mean():.4f} "
          f"({'↑ increases' if noncrit_all['delta_log_H'].mean() > 0 else '↓ decreases'})")
    frac_positive_crit = (crit_all["delta_log_H"] > 0).mean()
    frac_positive_noncrit = (noncrit_all["delta_log_H"] > 0).mean()
    print(f"  Fraction Δlog_H > 0: critical={frac_positive_crit:.3f}, non-critical={frac_positive_noncrit:.3f}")

    # C. Pooled test
    print("\n--- C. Pooled across N ---")
    if len(crit_all) >= 3 and len(noncrit_all) >= 3:
        t_stat, p_val = stats.ttest_ind(
            crit_all["delta_log_H"], noncrit_all["delta_log_H"],
            equal_var=False
        )
        n1, n2 = len(crit_all), len(noncrit_all)
        s1, s2 = crit_all["delta_log_H"].std(), noncrit_all["delta_log_H"].std()
        pooled_sd = math.sqrt(((n1-1)*s1**2 + (n2-1)*s2**2) / (n1+n2-2)) if (n1+n2-2) > 0 else 1
        d = (crit_all["delta_log_H"].mean() - noncrit_all["delta_log_H"].mean()) / pooled_sd if pooled_sd > 0 else float('nan')
        print(f"  Welch t={t_stat:.3f}, p={p_val:.6f}, Cohen's d={d:.3f}")

    # D. Symmetry test: compare |Δlog_H| of merge to |Δlog_H| of split
    print("\n--- D. Symmetry: merge vs split magnitude ---")
    split_file = Path("outputs_exploratory/prediction_c_intervention/intervention_raw.csv")
    if split_file.exists():
        split_df = pd.read_csv(split_file)
        split_data = split_df[split_df["split_occurred"] == 1]
        merge_data = crit_all

        for N in N_VALUES:
            sp = split_data[split_data["N"] == N]
            mg = merge_data[merge_data["N"] == N]
            if len(sp) > 0 and len(mg) > 0:
                print(f"\n  N={N}:")
                print(f"    Split  (add edge → +layer):  mean |Δlog_H| = {sp['abs_delta_log_H'].mean():.4f} (n={len(sp)})")
                print(f"    Merge  (remove edge → -layer): mean |Δlog_H| = {mg['abs_delta_log_H'].mean():.4f} (n={len(mg)})")
                ratio = mg["abs_delta_log_H"].mean() / sp["abs_delta_log_H"].mean() if sp["abs_delta_log_H"].mean() > 0 else float('nan')
                print(f"    Ratio merge/split: {ratio:.3f}")

        # Pooled
        if len(split_data) > 0 and len(merge_data) > 0:
            print(f"\n  Pooled: split |Δ|={split_data['abs_delta_log_H'].mean():.4f}, "
                  f"merge |Δ|={merge_data['abs_delta_log_H'].mean():.4f}, "
                  f"ratio={merge_data['abs_delta_log_H'].mean() / split_data['abs_delta_log_H'].mean():.3f}")

    # E. Save summary
    summary_rows = []
    for N in N_VALUES:
        for grp in ["critical", "non_critical"]:
            sub = df[(df["N"] == N) & (df["removal_type"] == grp)]
            if len(sub) == 0:
                continue
            summary_rows.append({
                "N": N, "type": grp,
                "n": len(sub),
                "mean_delta_layer_count": sub["delta_layer_count"].mean(),
                "mean_delta_log_H": sub["delta_log_H"].mean(),
                "mean_abs_delta_log_H": sub["abs_delta_log_H"].mean(),
                "frac_positive_delta_log_H": (sub["delta_log_H"] > 0).mean(),
            })
    pd.DataFrame(summary_rows).to_csv(OUT_DIR / "reverse_intervention_summary.csv", index=False)

    # F. Conclusion
    print("\n" + "=" * 70)
    if crit_all["delta_log_H"].mean() > 0 and crit_all["delta_log_H"].mean() > noncrit_all["delta_log_H"].mean():
        print("  ✓ REVERSE INTERVENTION SUPPORTED:")
        print("    Removing layer-merging edges INCREASES entropy")
        print("    Bidirectional causal link established:")
        print("    split → entropy ↓  AND  merge → entropy ↑")
    else:
        print("  ✗ Reverse intervention not supported")
    print("=" * 70)


if __name__ == "__main__":
    run()
