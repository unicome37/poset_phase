"""
Experiment: 4D de Sitter sprinkling — A2 gap closure.

Tests whether the layered screening architecture (S_BD + S_MD centred
on flat Lor4D) remains valid when curved-spacetime (de Sitter) families
are added to the 25-family poset library.

Key questions:
  Q1. Where does dS4D(H) rank under flat-Lor4D-centred S_MD?
  Q2. Does flat Lor4D remain #1 when dS4D families are present?
  Q3. How do dS4D features (d_eff, C1/C0, w/N) differ from flat Lor4D?
  Q4. At what Hubble H does dS4D become indistinguishable from flat Lor4D?
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

# ── project imports ──────────────────────────────────────────────────────
from generators import Poset, transitive_closure
from curvature_layer2_de_sitter_scan import (
    sprinkle_de_sitter_like_diamond,
    poset_from_de_sitter_points,
)
from expanded_family_robustness import (
    ALL_FAMILIES,
    compute_features,
    mahalanobis_score,
)
from bd_action import count_intervals_fast, bdg_action_d4_standard

# ── constants ────────────────────────────────────────────────────────────
SEED_BASE = 42
REPS = 40                          # realizations per family per N
N_GRID = [16, 28, 48, 64, 96, 128]
HUBBLE_VALUES = [0.0, 0.1, 0.3, 0.5, 1.0, 2.0]

OUT_DIR = Path("outputs_carlip")
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ── helpers ──────────────────────────────────────────────────────────────
def generate_ds4d(n: int, hubble: float, seed: int) -> Poset:
    """Generate a 4D de Sitter causal set (1 time + 3 spatial)."""
    pts = sprinkle_de_sitter_like_diamond(n, d_spatial=3, hubble=hubble, seed=seed)
    return poset_from_de_sitter_points(pts, hubble)


def feature_ensemble(gen_fn, N: int, reps: int, seed_base: int):
    """Compute feature matrix (reps × 3) for a generator function."""
    feats = []
    for r in range(reps):
        poset = gen_fn(N, seed=seed_base + r)
        f = compute_features(poset, N)
        feats.append(f)
    return np.array(feats)


def ds4d_feature_ensemble(N: int, hubble: float, reps: int, seed_base: int):
    """Compute feature matrix for dS4D at given Hubble."""
    feats = []
    for r in range(reps):
        poset = generate_ds4d(N, hubble, seed=seed_base + r)
        f = compute_features(poset, N)
        feats.append(f)
    return np.array(feats)


# ══════════════════════════════════════════════════════════════════════════
#  Experiment A: dS4D feature profiles across Hubble and N
# ══════════════════════════════════════════════════════════════════════════
def exp_a_feature_profiles(report):
    """Compare dS4D(H) features vs flat Lor4D across N and H."""
    report.append("=" * 72)
    report.append("# Exp A: dS4D Feature Profiles vs Flat Lor4D")
    report.append("=" * 72)
    report.append("")
    report.append("| N | H | d_eff | C1/C0 | w/N | "
                  "flat_d_eff | flat_C1C0 | flat_wN |")
    report.append("|---|---|-------|-------|-----|"
                  "-----------|-----------|---------|")

    for N in N_GRID:
        # flat Lor4D reference
        flat_feats = feature_ensemble(
            ALL_FAMILIES["Lor4D"], N, REPS, SEED_BASE + 7000
        )
        flat_mu = flat_feats.mean(axis=0)

        for H in HUBBLE_VALUES:
            ds_feats = ds4d_feature_ensemble(N, H, REPS, SEED_BASE + 8000)
            ds_mu = ds_feats.mean(axis=0)
            report.append(
                f"| {N} | {H:.1f} | {ds_mu[0]:.4f} | {ds_mu[1]:.4f} | "
                f"{ds_mu[2]:.4f} | {flat_mu[0]:.4f} | {flat_mu[1]:.4f} | "
                f"{flat_mu[2]:.4f} |"
            )
    report.append("")


# ══════════════════════════════════════════════════════════════════════════
#  Experiment B: S_MD rankings when dS4D families are added
# ══════════════════════════════════════════════════════════════════════════
def exp_b_smd_rankings(report):
    """Compute S_MD (flat-Lor4D centre) for 25 + dS4D families."""
    report.append("=" * 72)
    report.append("# Exp B: S_MD Rankings with dS4D Families Added")
    report.append("=" * 72)
    report.append("")

    for N in N_GRID:
        # Build flat Lor4D reference (μ, Σ)
        ref_feats = feature_ensemble(
            ALL_FAMILIES["Lor4D"], N, REPS, SEED_BASE + 9000
        )
        mu = ref_feats.mean(axis=0)
        cov = np.cov(ref_feats.T)
        cov += 1e-8 * np.eye(3)
        cov_inv = np.linalg.inv(cov)

        # Compute scores for all 25 standard families
        family_scores = {}
        for name, gen_fn in ALL_FAMILIES.items():
            feats = feature_ensemble(gen_fn, N, REPS, SEED_BASE + 5000)
            feat_mu = feats.mean(axis=0)
            family_scores[name] = mahalanobis_score(feat_mu, mu, cov_inv)

        # Compute scores for dS4D at each Hubble
        for H in HUBBLE_VALUES:
            if H == 0.0:
                continue  # skip H=0 (= flat Lor4D)
            ds_feats = ds4d_feature_ensemble(N, H, REPS, SEED_BASE + 8000)
            ds_mu = ds_feats.mean(axis=0)
            family_scores[f"dS4D_H{H:.1f}"] = mahalanobis_score(
                ds_mu, mu, cov_inv
            )

        # Rank all families
        sorted_fams = sorted(family_scores.items(), key=lambda x: x[1])
        lor4d_score = family_scores["Lor4D"]
        lor4d_rank = [name for name, _ in sorted_fams].index("Lor4D") + 1

        report.append(f"## N = {N}, total families = {len(sorted_fams)}")
        report.append(f"Lor4D rank: #{lor4d_rank}/{len(sorted_fams)}, "
                      f"score = {lor4d_score:.4f}")
        report.append("")
        report.append("| Rank | Family | S_MD | Gap from Lor4D |")
        report.append("|------|--------|------|----------------|")
        for rank, (name, score) in enumerate(sorted_fams[:10], 1):
            gap = score - lor4d_score
            marker = " ← *" if name == "Lor4D" else ""
            report.append(
                f"| {rank} | {name} | {score:.4f} | "
                f"{gap:+.4f}{marker} |"
            )
        # Also show all dS4D entries
        report.append("")
        report.append("**dS4D entries:**")
        report.append("| Family | S_MD | Rank |")
        report.append("|--------|------|------|")
        for rank, (name, score) in enumerate(sorted_fams, 1):
            if name.startswith("dS4D"):
                report.append(f"| {name} | {score:.4f} | #{rank} |")
        report.append("")


# ══════════════════════════════════════════════════════════════════════════
#  Experiment C: Feature-space distance from dS4D to flat Lor4D
# ══════════════════════════════════════════════════════════════════════════
def exp_c_feature_distances(report):
    """Euclidean and Mahalanobis distance from dS4D(H) to flat Lor4D."""
    report.append("=" * 72)
    report.append("# Exp C: Feature-Space Distance dS4D → Flat Lor4D")
    report.append("=" * 72)
    report.append("")
    report.append("| N | H | Euclid_dist | Mahal_dist | "
                  "Δ(d_eff) | Δ(C1/C0) | Δ(w/N) |")
    report.append("|---|---|-------------|------------|"
                  "---------|----------|--------|")

    for N in N_GRID:
        ref_feats = feature_ensemble(
            ALL_FAMILIES["Lor4D"], N, REPS, SEED_BASE + 9000
        )
        mu = ref_feats.mean(axis=0)
        cov = np.cov(ref_feats.T)
        cov += 1e-8 * np.eye(3)
        cov_inv = np.linalg.inv(cov)

        for H in HUBBLE_VALUES:
            if H == 0.0:
                continue
            ds_feats = ds4d_feature_ensemble(N, H, REPS, SEED_BASE + 8000)
            ds_mu = ds_feats.mean(axis=0)
            delta = ds_mu - mu
            euclid = np.linalg.norm(delta)
            mahal = mahalanobis_score(ds_mu, mu, cov_inv)

            report.append(
                f"| {N} | {H:.1f} | {euclid:.6f} | {mahal:.4f} | "
                f"{delta[0]:+.4f} | {delta[1]:+.4f} | {delta[2]:+.4f} |"
            )
    report.append("")


# ══════════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════════
def main():
    t0 = time.time()
    report = [
        "# 4D de Sitter Sprinkling Experiment",
        f"Date: {time.strftime('%Y-%m-%d %H:%M')}",
        f"N grid: {N_GRID}",
        f"Hubble values: {HUBBLE_VALUES}",
        f"REPS per family: {REPS}",
        f"Flat families: {len(ALL_FAMILIES)}",
        "",
    ]

    print("=== Exp A: Feature Profiles ===")
    exp_a_feature_profiles(report)

    print("=== Exp B: S_MD Rankings ===")
    exp_b_smd_rankings(report)

    print("=== Exp C: Feature Distances ===")
    exp_c_feature_distances(report)

    elapsed = time.time() - t0
    report.append(f"---\nTotal runtime: {elapsed:.0f}s")

    out_file = OUT_DIR / "experiment_de_sitter.md"
    out_file.write_text("\n".join(report), encoding="utf-8")
    print(f"Report written to: {out_file}")
    print(f"Total runtime: {elapsed:.0f}s")


if __name__ == "__main__":
    main()
