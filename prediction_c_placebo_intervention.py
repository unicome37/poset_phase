#!/usr/bin/env python3
"""Prediction C — Placebo-Controlled Intervention Experiment.

Extends the original intervention experiment with a crucial PLACEBO control:

  Treatment group : add edge between incomparable pair in SAME layer
                    (can increase layer_count)
  Placebo group   : add edge between incomparable pair in DIFFERENT layers
                    (cannot increase layer_count by construction)

If Prediction C's mechanism is correct, the treatment group should show
larger |Δlog_H| than the placebo group. This rules out the null hypothesis
that "any added comparability reduces entropy equally, regardless of whether
it creates new layers."

Additionally, we extend the original experiment to Lor3D and Lor4D to test
cross-dimensional generalizability.
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
    generate_lorentzian_like_4d,
)
from entropy_exact import count_linear_extensions_exact
from matched_residual_freedom_check import layer_index_by_minima

OUT_DIR = Path("outputs_exploratory/prediction_c_placebo_intervention")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Config ──────────────────────────────────────────────────
FAMILIES = {
    "Lor2D": (generate_lorentzian_like_2d, [14, 16]),
    "Lor3D": (generate_lorentzian_like_3d, [14, 16]),
    "Lor4D": (generate_lorentzian_like_4d, [12, 14]),  # 4D denser → smaller N
}
POSETS_PER_N = 40
MAX_INTERVENTIONS_PER_TYPE = 10  # per-type (same-layer / cross-layer)
SEED_BASE = 20260317


# ── Helpers ─────────────────────────────────────────────────
def layer_count(poset: Poset) -> int:
    idx = layer_index_by_minima(poset.closure)
    return int(idx.max() + 1) if len(idx) > 0 else 0


def log_H(poset: Poset) -> float:
    return math.log(count_linear_extensions_exact(poset, prefer_c=False))


def incomparable_pairs(poset: Poset):
    """Return (same_layer_pairs, cross_layer_pairs)."""
    c = poset.closure
    n = poset.n
    layer_idx = layer_index_by_minima(c)
    same, cross = [], []
    for i in range(n):
        for j in range(i + 1, n):
            if not c[i, j] and not c[j, i]:
                if layer_idx[i] == layer_idx[j]:
                    same.append((i, j))
                else:
                    # ensure i is in lower layer for consistent direction
                    if layer_idx[i] < layer_idx[j]:
                        cross.append((i, j))
                    else:
                        cross.append((j, i))
    return same, cross


def intervene(poset: Poset, a: int, b: int) -> Poset:
    new_adj = poset.closure.copy()
    new_adj[a, b] = True
    return Poset(transitive_closure(new_adj))


# ── Main ────────────────────────────────────────────────────
def run():
    rng = np.random.default_rng(SEED_BASE)
    rows = []
    t0 = time.time()

    for fam_name, (gen_fn, n_values) in FAMILIES.items():
        for N in n_values:
            print(f"\n=== {fam_name}, N={N} ===")
            for idx in range(POSETS_PER_N):
                seed = SEED_BASE + hash(fam_name) % 10000 + N * 1000 + idx
                # ensure positive seed
                seed = abs(seed)
                poset = gen_fn(N, seed=seed)
                lc0 = layer_count(poset)
                lh0 = log_H(poset)

                same_pairs, cross_pairs = incomparable_pairs(poset)

                # Sample from each type
                for group_name, pool in [("same_layer", same_pairs),
                                         ("cross_layer", cross_pairs)]:
                    if len(pool) == 0:
                        continue
                    sample_size = min(MAX_INTERVENTIONS_PER_TYPE, len(pool))
                    chosen = rng.choice(len(pool), size=sample_size, replace=False)
                    for ci in chosen:
                        a, b = pool[ci]
                        new_poset = intervene(poset, a, b)
                        lc1 = layer_count(new_poset)
                        lh1 = log_H(new_poset)
                        rows.append({
                            "family": fam_name,
                            "N": N,
                            "poset_seed": seed,
                            "original_layer_count": lc0,
                            "original_log_H": lh0,
                            "intervention_type": group_name,
                            "edge_added": f"{a}->{b}",
                            "new_layer_count": lc1,
                            "new_log_H": lh1,
                            "delta_layer_count": lc1 - lc0,
                            "delta_log_H": lh1 - lh0,
                            "abs_delta_log_H": abs(lh1 - lh0),
                            "split_occurred": int(lc1 > lc0),
                        })

                if (idx + 1) % 10 == 0:
                    elapsed = time.time() - t0
                    print(f"  {fam_name} N={N}: {idx+1}/{POSETS_PER_N} ({elapsed:.0f}s)")

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "placebo_intervention_raw.csv", index=False)
    print(f"\nTotal interventions: {len(df)}")

    # ══════════════════════════════════════════════════════════
    #  ANALYSIS
    # ══════════════════════════════════════════════════════════
    print("\n" + "=" * 70)
    print("ANALYSIS: Placebo-Controlled Intervention")
    print("=" * 70)

    # ── A. Same-layer vs Cross-layer comparison (main placebo test) ──
    print("\n--- A. Same-layer (TREATMENT) vs Cross-layer (PLACEBO) ---")
    results_rows = []

    for fam_name in FAMILIES:
        sub = df[df["family"] == fam_name]
        treat = sub[sub["intervention_type"] == "same_layer"]
        placebo = sub[sub["intervention_type"] == "cross_layer"]
        print(f"\n  [{fam_name}]")
        print(f"    Treatment (same-layer):   n={len(treat)}, "
              f"mean |Δlog_H|={treat['abs_delta_log_H'].mean():.4f}, "
              f"mean Δlog_H={treat['delta_log_H'].mean():.4f}")
        print(f"    Placebo   (cross-layer):  n={len(placebo)}, "
              f"mean |Δlog_H|={placebo['abs_delta_log_H'].mean():.4f}, "
              f"mean Δlog_H={placebo['delta_log_H'].mean():.4f}")
        print(f"    Treatment split rate:     {treat['split_occurred'].mean():.3f}")
        print(f"    Placebo split rate:       {placebo['split_occurred'].mean():.3f}")

        if len(treat) >= 3 and len(placebo) >= 3:
            t_stat, p_val = stats.ttest_ind(
                treat["abs_delta_log_H"], placebo["abs_delta_log_H"],
                equal_var=False
            )
            n1, n2 = len(treat), len(placebo)
            s1, s2 = treat["abs_delta_log_H"].std(), placebo["abs_delta_log_H"].std()
            pooled_sd = math.sqrt(((n1-1)*s1**2 + (n2-1)*s2**2) / (n1+n2-2)) if (n1+n2-2) > 0 else 1
            d = (treat["abs_delta_log_H"].mean() - placebo["abs_delta_log_H"].mean()) / pooled_sd if pooled_sd > 0 else float('nan')
            print(f"    Welch t={t_stat:.3f}, p={p_val:.6f}, Cohen's d={d:.3f}")

            results_rows.append({
                "family": fam_name,
                "n_treatment": n1, "n_placebo": n2,
                "mean_abs_dlogH_treatment": treat["abs_delta_log_H"].mean(),
                "mean_abs_dlogH_placebo": placebo["abs_delta_log_H"].mean(),
                "treatment_split_rate": treat["split_occurred"].mean(),
                "placebo_split_rate": placebo["split_occurred"].mean(),
                "welch_t": t_stat, "welch_p": p_val, "cohens_d": d,
            })

    results_df = pd.DataFrame(results_rows)
    results_df.to_csv(OUT_DIR / "placebo_comparison.csv", index=False)

    # ── B. Within treatment: split vs no-split (replication of original) ──
    print("\n\n--- B. Within same-layer interventions: split vs no-split ---")
    for fam_name in FAMILIES:
        treat = df[(df["family"] == fam_name) & (df["intervention_type"] == "same_layer")]
        split = treat[treat["split_occurred"] == 1]
        nosplit = treat[treat["split_occurred"] == 0]
        print(f"\n  [{fam_name}]")
        print(f"    Split:    n={len(split)}, mean |Δlog_H|={split['abs_delta_log_H'].mean():.4f}" if len(split)>0 else f"    Split:    n=0")
        print(f"    No-split: n={len(nosplit)}, mean |Δlog_H|={nosplit['abs_delta_log_H'].mean():.4f}" if len(nosplit)>0 else f"    No-split: n=0")
        if len(split) >= 3 and len(nosplit) >= 3:
            t_stat, p_val = stats.ttest_ind(
                split["abs_delta_log_H"], nosplit["abs_delta_log_H"],
                equal_var=False
            )
            print(f"    Welch t={t_stat:.3f}, p={p_val:.6f}")

    # ── C. Pooled placebo test ──────────────────────────────
    print("\n\n--- C. Pooled across all families ---")
    treat_all = df[df["intervention_type"] == "same_layer"]
    placebo_all = df[df["intervention_type"] == "cross_layer"]
    print(f"  Treatment: n={len(treat_all)}, mean |Δlog_H|={treat_all['abs_delta_log_H'].mean():.4f}")
    print(f"  Placebo:   n={len(placebo_all)}, mean |Δlog_H|={placebo_all['abs_delta_log_H'].mean():.4f}")
    if len(treat_all) >= 3 and len(placebo_all) >= 3:
        t_stat, p_val = stats.ttest_ind(
            treat_all["abs_delta_log_H"], placebo_all["abs_delta_log_H"],
            equal_var=False
        )
        n1, n2 = len(treat_all), len(placebo_all)
        s1, s2 = treat_all["abs_delta_log_H"].std(), placebo_all["abs_delta_log_H"].std()
        pooled_sd = math.sqrt(((n1-1)*s1**2 + (n2-1)*s2**2) / (n1+n2-2)) if (n1+n2-2) > 0 else 1
        d = (treat_all["abs_delta_log_H"].mean() - placebo_all["abs_delta_log_H"].mean()) / pooled_sd if pooled_sd > 0 else float('nan')
        print(f"  Welch t={t_stat:.3f}, p={p_val:.6f}, Cohen's d={d:.3f}")

    # ── D. Summary ──────────────────────────────────────────
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    if len(results_rows) > 0:
        all_positive_d = all(r["cohens_d"] > 0 for r in results_rows)
        all_sig = all(r["welch_p"] < 0.05 for r in results_rows)
        if all_positive_d and all_sig:
            print("  ✓ PLACEBO TEST PASSED: Same-layer interventions reduce entropy")
            print("    more than cross-layer interventions across ALL families.")
            print("    The layer-splitting mechanism is supported as causal pathway.")
        elif all_positive_d:
            print("  ~ Direction supported in all families, not all reach p<0.05")
        else:
            print("  ✗ Placebo test FAILED in at least one family")

    # Check cross-dimensional generality
    fam_names = list(FAMILIES.keys())
    non_2d = [r for r in results_rows if r["family"] != "Lor2D"]
    if non_2d and all(r["cohens_d"] > 0 for r in non_2d):
        print("  ✓ Cross-dimensional: Effect holds in Lor3D and/or Lor4D")
    print("=" * 70)


if __name__ == "__main__":
    run()
