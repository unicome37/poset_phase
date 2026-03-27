"""
Triple Screening Quantitative Test
====================================
Tests whether Lor4D is uniquely selected by the three-fold screening:
  E₁  Existability  — non-trivial causal depth
  E₂  Sustainability — feature stability under coarse-graining
  E₃  Manifestability — low variance of features across realizations

Composite score  S_triple = E₁ × E₂ × E₃   (must pass ALL three).
"""
from __future__ import annotations

import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

from generators import (
    Poset,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
    generate_kr_like,
    generate_kr_2layer,
    generate_kr_4layer,
    generate_absolute_layered,
    generate_multi_layer_random,
    generate_random_layered_k4_uniform,
    generate_random_layered_k6_uniform,
    generate_random_layered_k8_uniform,
    generate_random_layered_k6_tapered,
    generate_random_layered_k6_middle_heavy,
    generate_random_layered_k6_longjump,
    generate_transitive_percolation,
    generate_interval_order,
)
from bd_action import count_intervals_fast
from unified_functional import compute_xi_dim

# ---------------------------------------------------------------------------
# Family registry
# ---------------------------------------------------------------------------
FAMILIES: dict[str, callable] = {
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

N_VALUES = [16, 28, 48, 64, 96, 128]
N_REPS = 20
N_CG_SAMPLES = 5
CG_KEEP_RATIO = 0.7

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

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


def compute_features(poset: Poset) -> tuple[float, float, float]:
    """Return (d_eff, c1_c0, width_ratio)."""
    counts = count_intervals_fast(poset, k_max=5)
    C0 = counts.get(0)
    C1 = counts.get(1)
    c1_c0 = C1 / max(1, C0)
    _, d_eff = compute_xi_dim(poset)
    aw = max_antichain_width(poset)
    width_ratio = aw / max(1, poset.n)
    return d_eff, c1_c0, width_ratio


def coarse_grain(poset: Poset, keep_ratio: float, rng: np.random.Generator) -> Poset:
    """Random node deletion coarse-graining."""
    n = poset.n
    n_keep = max(3, int(n * keep_ratio))
    keep_idx = np.sort(rng.choice(n, size=n_keep, replace=False))
    sub = poset.closure[np.ix_(keep_idx, keep_idx)]
    return Poset(closure=sub)


# ---------------------------------------------------------------------------
# Screening Layer 1: Existability
# ---------------------------------------------------------------------------

def existability(poset: Poset) -> float:
    """Continuous existability score based on causal structure richness.
    
    A physically meaningful "existability" should measure whether the structure
    has non-trivial causal depth AND non-trivial lateral extent — i.e., it is
    neither a trivial antichain nor a trivial chain.
    
    Score = geometric mean of:
      - depth_score: normalized layer count (n_layers / N^alpha, alpha chosen so Lor4D ~ 1)
      - richness_score: interval richness C1/C0 (higher = richer causal structure)
      - non_degeneracy: 1 - max(layer_size)/N (not dominated by a single layer)
    """
    from observables import layer_profile as _layer_profile
    n = poset.n
    if n < 3:
        return 0.0
    
    profile = _layer_profile(poset)
    n_layers = len(profile)
    
    # Depth score: n_layers normalized by sqrt(N) 
    # For Lor4D d=4: n_layers ~ N^{3/4}, so n_layers/sqrt(N) ~ N^{1/4} → grows
    # For KR_2layer: n_layers ~ 2-3, so this → 0
    # Use sigmoid: score = 1 / (1 + exp(-k*(x - x0)))
    depth_raw = n_layers / max(1, n ** 0.5)
    depth_score = 1.0 / (1.0 + np.exp(-5.0 * (depth_raw - 0.5)))
    
    # Non-degeneracy: 1 - max_fraction (fraction of elements in largest layer)
    max_frac = float(np.max(profile)) / n
    non_degen = 1.0 - max_frac
    
    # Interval richness: C1/C0 > 0 means non-trivial intervals exist
    counts = count_intervals_fast(poset, k_max=3)
    C0 = max(1, counts.get(0))
    C1 = counts.get(1)
    richness_raw = C1 / C0
    richness_score = 1.0 / (1.0 + np.exp(-10.0 * (richness_raw - 0.05)))
    
    return float((depth_score * non_degen * richness_score) ** (1.0 / 3.0))


# ---------------------------------------------------------------------------
# Screening Layer 2: Sustainability (coarse-graining drift)
# ---------------------------------------------------------------------------

def sustainability(poset: Poset, rng: np.random.Generator) -> float:
    """1 / (1 + mean_drift²)  over N_CG_SAMPLES coarse-grainings."""
    try:
        feat_orig = np.array(compute_features(poset))
    except Exception:
        return 0.0
    drifts = []
    for _ in range(N_CG_SAMPLES):
        try:
            cg = coarse_grain(poset, CG_KEEP_RATIO, rng)
            feat_cg = np.array(compute_features(cg))
            drift = np.linalg.norm(feat_orig - feat_cg)
            drifts.append(drift)
        except Exception:
            drifts.append(10.0)  # large penalty
    mean_drift = float(np.mean(drifts))
    return 1.0 / (1.0 + mean_drift ** 2)


# ---------------------------------------------------------------------------
# Screening Layer 3: Manifestability (feature CV across realizations)
# ---------------------------------------------------------------------------

def manifestability_from_features(feat_array: np.ndarray) -> float:
    """Given (n_reps, 3) array of features, return 1 / (1 + mean_CV)."""
    means = np.mean(feat_array, axis=0)
    stds = np.std(feat_array, axis=0)
    # CV = std / |mean|, guard against zero mean
    cvs = stds / np.maximum(np.abs(means), 1e-12)
    mean_cv = float(np.mean(cvs))
    return 1.0 / (1.0 + mean_cv)


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

def run_experiment():
    rng = np.random.default_rng(42)
    out_dir = Path("outputs_carlip")
    out_dir.mkdir(exist_ok=True)

    # results[family][N] = {"E1": ..., "E2": ..., "E3": ..., "S": ...}
    results: dict[str, dict[int, dict]] = defaultdict(dict)

    total_configs = len(FAMILIES) * len(N_VALUES)
    config_idx = 0

    for fam_name, gen_fn in FAMILIES.items():
        for N in N_VALUES:
            config_idx += 1
            t0 = time.time()

            e1_vals = []
            e2_vals = []
            feat_list = []  # for E3

            for rep in range(N_REPS):
                try:
                    poset = gen_fn(N, rng=rng)
                except TypeError:
                    poset = gen_fn(N)

                # E1: existability
                e1_vals.append(existability(poset))

                # E2: sustainability
                e2_vals.append(sustainability(poset, rng))

                # Features for E3
                try:
                    feat = compute_features(poset)
                    feat_list.append(feat)
                except Exception:
                    feat_list.append((0.0, 0.0, 0.0))

            # Aggregate
            E1 = float(np.mean(e1_vals))
            E2 = float(np.mean(e2_vals))

            feat_array = np.array(feat_list)
            E3 = manifestability_from_features(feat_array)

            S_triple = E1 * E2 * E3

            results[fam_name][N] = {
                "E1": E1, "E2": E2, "E3": E3, "S": S_triple,
            }

            elapsed = time.time() - t0
            print(f"  [{config_idx:3d}/{total_configs}] {fam_name:12s} N={N:4d}  "
                  f"E1={E1:.3f} E2={E2:.3f} E3={E3:.3f} S={S_triple:.4f}  "
                  f"({elapsed:.1f}s)")

    # ------------------------------------------------------------------
    # Report
    # ------------------------------------------------------------------
    lines: list[str] = []
    lines.append("=" * 90)
    lines.append("Triple Screening Quantitative Results")
    lines.append("=" * 90)
    lines.append("")

    # Per-N tables
    for N in N_VALUES:
        lines.append(f"--- N = {N} ---")
        header = f"{'Family':14s} {'E1(exist)':>10s} {'E2(sustain)':>12s} {'E3(manifest)':>13s} {'S_triple':>10s}"
        lines.append(header)
        lines.append("-" * len(header))

        # Sort by S_triple descending
        ranking = sorted(FAMILIES.keys(),
                         key=lambda f: results[f][N]["S"], reverse=True)
        for fam in ranking:
            r = results[fam][N]
            lines.append(f"{fam:14s} {r['E1']:10.3f} {r['E2']:12.4f} {r['E3']:13.4f} {r['S']:10.5f}")
        lines.append("")

    # Summary: how many families pass each layer at each N
    lines.append("=" * 90)
    lines.append("Screening Funnel: families with score > 0.5 at each layer")
    lines.append("=" * 90)
    for N in N_VALUES:
        pass_e1 = [f for f in FAMILIES if results[f][N]["E1"] > 0.5]
        pass_e2 = [f for f in FAMILIES if results[f][N]["E2"] > 0.5]
        pass_e3 = [f for f in FAMILIES if results[f][N]["E3"] > 0.5]
        pass_all = [f for f in FAMILIES
                    if results[f][N]["E1"] > 0.5
                    and results[f][N]["E2"] > 0.5
                    and results[f][N]["E3"] > 0.5]
        lines.append(f"N={N:4d}: E1 pass {len(pass_e1):2d}/17  "
                      f"E2 pass {len(pass_e2):2d}/17  "
                      f"E3 pass {len(pass_e3):2d}/17  "
                      f"ALL pass {len(pass_all):2d}/17  {pass_all}")

    # Lor4D uniqueness check
    lines.append("")
    lines.append("=" * 90)
    lines.append("Lor4D Uniqueness Check: does Lor4D rank #1 in S_triple?")
    lines.append("=" * 90)
    for N in N_VALUES:
        ranking = sorted(FAMILIES.keys(),
                         key=lambda f: results[f][N]["S"], reverse=True)
        top = ranking[0]
        top_s = results[top][N]["S"]
        lor4d_s = results["Lor4D"][N]["S"]
        lor4d_rank = ranking.index("Lor4D") + 1
        is_unique = (lor4d_rank == 1)
        gap = lor4d_s - results[ranking[1]][N]["S"] if is_unique else top_s - lor4d_s
        lines.append(f"N={N:4d}: Lor4D rank={lor4d_rank:2d}  "
                      f"S={lor4d_s:.5f}  top={top}({top_s:.5f})  "
                      f"gap={'+'if is_unique else '-'}{abs(gap):.5f}  "
                      f"{'*** UNIQUE #1 ***' if is_unique else ''}")

    # Composite ranking across all N (geometric mean of S)
    lines.append("")
    lines.append("=" * 90)
    lines.append("Composite Ranking: geometric mean of S_triple across all N")
    lines.append("=" * 90)
    geo_means = {}
    for fam in FAMILIES:
        scores = [results[fam][N]["S"] for N in N_VALUES]
        # geometric mean, treating 0 as epsilon to allow log
        scores_safe = [max(s, 1e-15) for s in scores]
        geo_means[fam] = float(np.exp(np.mean(np.log(scores_safe))))
    geo_ranking = sorted(geo_means.items(), key=lambda x: x[1], reverse=True)
    for rank, (fam, gm) in enumerate(geo_ranking, 1):
        marker = " <<<" if fam == "Lor4D" else ""
        lines.append(f"  #{rank:2d}  {fam:14s}  geo_mean(S) = {gm:.6f}{marker}")

    report = "\n".join(lines)
    print("\n" + report)

    out_path = out_dir / "triple_screening_results.txt"
    out_path.write_text(report, encoding="utf-8")
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    run_experiment()
