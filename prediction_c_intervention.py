#!/usr/bin/env python3
"""Prediction C — Layer Split/Merge Intervention Experiment.

Quasi-causal test: does *adding* a comparability relation that increases
layer_count reduce log_H more than one that does not?

Protocol
--------
1. Generate Lor2D posets at small N where exact counting is feasible.
2. For each poset, locate ALL incomparable pairs in the same layer.
3. For a random sample of such pairs, add the relation a < b, recompute
   transitive closure, and measure new (layer_count, log_H).
4. Categorise each single-edge intervention by whether it increased
   layer_count (split) or not (no-split).
5. Compare groups:
   - Welch t-test on |Δlog_H| between split vs no-split.
   - Regression of Δlog_H on Δlayer_count.
"""
from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from generators import Poset, transitive_closure, generate_lorentzian_like_2d
from entropy_exact import count_linear_extensions_exact
from matched_residual_freedom_check import layer_index_by_minima

# ── Output dir ──────────────────────────────────────────────
OUT_DIR = Path("outputs_exploratory/prediction_c_intervention")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Config ──────────────────────────────────────────────────
N_VALUES = [14, 16]          # poset sizes (exact counting feasible)
POSETS_PER_N = 40            # how many Lor2D posets per N
MAX_INTERVENTIONS = 15       # max single-edge interventions per poset
SEED_BASE = 20260317


# ── Helpers ─────────────────────────────────────────────────
def layer_count(poset: Poset) -> int:
    idx = layer_index_by_minima(poset.closure)
    return int(idx.max() + 1) if len(idx) > 0 else 0


def log_H(poset: Poset) -> float:
    return math.log(count_linear_extensions_exact(poset, prefer_c=False))


def incomparable_pairs_same_layer(poset: Poset) -> list[tuple[int, int]]:
    """Return incomparable (i,j) pairs that share the same layer."""
    c = poset.closure
    n = poset.n
    layer_idx = layer_index_by_minima(c)
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            if not c[i, j] and not c[j, i] and layer_idx[i] == layer_idx[j]:
                pairs.append((i, j))
    return pairs


def intervene_add_edge(poset: Poset, a: int, b: int) -> Poset:
    """Add the relation a < b and recompute transitive closure."""
    new_adj = poset.closure.copy()
    new_adj[a, b] = True
    return Poset(transitive_closure(new_adj))


# ── Main experiment ─────────────────────────────────────────
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

            pairs = incomparable_pairs_same_layer(poset)
            if len(pairs) == 0:
                continue

            sample_size = min(MAX_INTERVENTIONS, len(pairs))
            chosen = rng.choice(len(pairs), size=sample_size, replace=False)

            for ci in chosen:
                a, b = pairs[ci]
                new_poset = intervene_add_edge(poset, a, b)
                lc1 = layer_count(new_poset)
                lh1 = log_H(new_poset)
                rows.append({
                    "N": N,
                    "poset_seed": seed,
                    "poset_idx": idx,
                    "original_layer_count": lc0,
                    "original_log_H": lh0,
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
                print(f"  N={N}: {idx+1}/{POSETS_PER_N} posets done  ({elapsed:.0f}s)")

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "intervention_raw.csv", index=False)
    print(f"\nTotal interventions: {len(df)}")

    # ── Analysis ────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("ANALYSIS: Layer Split Intervention Experiment")
    print("=" * 60)

    for N in N_VALUES:
        sub = df[df["N"] == N]
        split = sub[sub["split_occurred"] == 1]
        nosplit = sub[sub["split_occurred"] == 0]
        print(f"\n--- N = {N} ---")
        print(f"  Total interventions: {len(sub)}")
        print(f"  Split (Δlayer > 0): {len(split)}  |  No-split: {len(nosplit)}")

        if len(split) < 3 or len(nosplit) < 3:
            print("  ⚠ Too few observations in one group; skipping t-test")
            continue

        # Welch t-test on |Δlog_H|
        t_stat, p_val = stats.ttest_ind(
            split["abs_delta_log_H"], nosplit["abs_delta_log_H"],
            equal_var=False
        )
        print(f"  Mean |Δlog_H|  (split):    {split['abs_delta_log_H'].mean():.4f}")
        print(f"  Mean |Δlog_H|  (no-split): {nosplit['abs_delta_log_H'].mean():.4f}")
        print(f"  Welch t = {t_stat:.3f}, p = {p_val:.4f}")

        # Also: Δlog_H mean (should be negative for both, more negative for split)
        print(f"  Mean Δlog_H    (split):    {split['delta_log_H'].mean():.4f}")
        print(f"  Mean Δlog_H    (no-split): {nosplit['delta_log_H'].mean():.4f}")

        # Regression: Δlog_H ~ Δlayer_count
        if sub["delta_layer_count"].nunique() > 1:
            slope, intercept, r, p_reg, se = stats.linregress(
                sub["delta_layer_count"], sub["delta_log_H"]
            )
            print(f"  Regression Δlog_H ~ Δlayer_count: slope={slope:.4f}, r={r:.3f}, p={p_reg:.4f}")

    # ── Pooled analysis ────────────────────────────────────
    print("\n--- Pooled across N ---")
    split_all = df[df["split_occurred"] == 1]
    nosplit_all = df[df["split_occurred"] == 0]
    print(f"  Split: n={len(split_all)}, No-split: n={len(nosplit_all)}")
    if len(split_all) >= 3 and len(nosplit_all) >= 3:
        t_stat, p_val = stats.ttest_ind(
            split_all["abs_delta_log_H"], nosplit_all["abs_delta_log_H"],
            equal_var=False
        )
        print(f"  Mean |Δlog_H|  (split):    {split_all['abs_delta_log_H'].mean():.4f}")
        print(f"  Mean |Δlog_H|  (no-split): {nosplit_all['abs_delta_log_H'].mean():.4f}")
        print(f"  Welch t = {t_stat:.3f}, p = {p_val:.4f}")

        # Effect size (Cohen's d)
        n1, n2 = len(split_all), len(nosplit_all)
        s1, s2 = split_all["abs_delta_log_H"].std(), nosplit_all["abs_delta_log_H"].std()
        pooled_sd = math.sqrt(((n1 - 1) * s1**2 + (n2 - 1) * s2**2) / (n1 + n2 - 2))
        cohen_d = (split_all["abs_delta_log_H"].mean() - nosplit_all["abs_delta_log_H"].mean()) / pooled_sd if pooled_sd > 0 else float('nan')
        print(f"  Cohen's d = {cohen_d:.3f}")

    # Regression pooled
    if df["delta_layer_count"].nunique() > 1:
        slope, intercept, r, p_reg, se = stats.linregress(
            df["delta_layer_count"], df["delta_log_H"]
        )
        print(f"  Regression pooled: slope={slope:.4f}, r={r:.3f}, p={p_reg:.4f}")

    # ── Summary table ───────────────────────────────────────
    summary_rows = []
    for N in N_VALUES:
        sub = df[df["N"] == N]
        for grp_name, grp in [("split", sub[sub["split_occurred"] == 1]),
                               ("no-split", sub[sub["split_occurred"] == 0])]:
            if len(grp) == 0:
                continue
            summary_rows.append({
                "N": N,
                "group": grp_name,
                "n_interventions": len(grp),
                "mean_delta_layer_count": grp["delta_layer_count"].mean(),
                "mean_delta_log_H": grp["delta_log_H"].mean(),
                "mean_abs_delta_log_H": grp["abs_delta_log_H"].mean(),
                "std_abs_delta_log_H": grp["abs_delta_log_H"].std(),
            })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT_DIR / "intervention_summary.csv", index=False)
    print(f"\n  Saved: intervention_summary.csv")

    # ── Conclusion ──────────────────────────────────────────
    print("\n" + "=" * 60)
    if len(split_all) >= 3 and len(nosplit_all) >= 3:
        direction = "SUPPORTED" if split_all["abs_delta_log_H"].mean() > nosplit_all["abs_delta_log_H"].mean() else "NOT supported"
        print(f"  Prediction C direction: {direction}")
        print(f"  (Adding an edge that increases layer count should")
        print(f"   reduce entropy MORE than one that does not.)")
    print("=" * 60)


if __name__ == "__main__":
    run()
