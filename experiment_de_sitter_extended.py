"""
Extended de Sitter experiment — finer H grid + larger N.
Builds on experiment_de_sitter.py with:
  - H grid: 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0
  - N grid: 16, 28, 48, 64, 96, 128, 192, 256
  - Focus on the H_c boundary where dS4D goes from "inside basin" to "outside"
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

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

SEED_BASE = 42
REPS = 40
N_GRID = [16, 28, 48, 64, 96, 128, 192, 256]
HUBBLE_FINE = [0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0]

OUT_DIR = Path("outputs_carlip")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def generate_ds4d(n, hubble, seed):
    pts = sprinkle_de_sitter_like_diamond(n, d_spatial=3, hubble=hubble, seed=seed)
    return poset_from_de_sitter_points(pts, hubble)


def feature_ensemble_std(gen_fn, N, reps, seed_base):
    feats = []
    for r in range(reps):
        poset = gen_fn(N, seed=seed_base + r)
        f = compute_features(poset, N)
        feats.append(f)
    return np.array(feats)


def ds4d_feature_ensemble(N, hubble, reps, seed_base):
    feats = []
    for r in range(reps):
        poset = generate_ds4d(N, hubble, seed=seed_base + r)
        f = compute_features(poset, N)
        feats.append(f)
    return np.array(feats)


def main():
    t0 = time.time()
    report = ["# Extended de Sitter Experiment\n"]
    report.append(f"N grid: {N_GRID}")
    report.append(f"H grid: {HUBBLE_FINE}")
    report.append(f"REPS per (family, N): {REPS}\n")

    # ================================================================
    # Part 1: S_MD score of dS4D(H) — fine H scan at all N
    # ================================================================
    report.append("=" * 72)
    report.append("# Part 1: S_MD(dS4D_H) vs flat Lor4D — Fine Hubble Scan")
    report.append("=" * 72)
    report.append("")
    report.append("| N | H | S_MD(dS4D) | S_MD(Lor4D) | Δ | Lor4D_rank | dS4D_rank | total_fams |")
    report.append("|---|---|------------|-------------|---|------------|-----------|------------|")

    for N in N_GRID:
        print(f"  N={N} ...", flush=True)
        ref_feats = feature_ensemble_std(ALL_FAMILIES["Lor4D"], N, REPS, SEED_BASE + 9000)
        mu = ref_feats.mean(axis=0)
        cov = np.cov(ref_feats.T) + 1e-8 * np.eye(3)
        cov_inv = np.linalg.inv(cov)

        # Standard families scores
        std_scores = {}
        for name, gen_fn in ALL_FAMILIES.items():
            feats = feature_ensemble_std(gen_fn, N, REPS, SEED_BASE + 5000)
            std_scores[name] = mahalanobis_score(feats.mean(axis=0), mu, cov_inv)

        lor4d_score = std_scores["Lor4D"]

        for H in HUBBLE_FINE:
            ds_feats = ds4d_feature_ensemble(N, H, REPS, SEED_BASE + 8000)
            ds_mu = ds_feats.mean(axis=0)
            ds_score = mahalanobis_score(ds_mu, mu, cov_inv)

            all_scores = {**std_scores, f"dS4D_H{H}": ds_score}
            sorted_all = sorted(all_scores.items(), key=lambda x: x[1])
            names = [n for n, _ in sorted_all]
            lor4d_rank = names.index("Lor4D") + 1
            ds_rank = names.index(f"dS4D_H{H}") + 1

            delta = ds_score - lor4d_score
            report.append(
                f"| {N} | {H} | {ds_score:.4f} | {lor4d_score:.4f} | "
                f"{delta:+.4f} | #{lor4d_rank} | #{ds_rank} | {len(all_scores)} |"
            )

    report.append("")

    # ================================================================
    # Part 2: Critical Hubble H_c — where dS4D leaves top-2
    # ================================================================
    report.append("=" * 72)
    report.append("# Part 2: Critical Hubble H_c Summary")
    report.append("=" * 72)
    report.append("")
    report.append("H_c is defined as the largest H where dS4D still ranks #2 (or better).")
    report.append("")
    report.append("| N | H_c (dS4D ≤ #2) | S_MD at H_c | Nearest non-dS4D competitor | Competitor S_MD |")
    report.append("|---|------------------|-------------|----------------------------|-----------------|")

    for N in N_GRID:
        ref_feats = feature_ensemble_std(ALL_FAMILIES["Lor4D"], N, REPS, SEED_BASE + 9000)
        mu = ref_feats.mean(axis=0)
        cov = np.cov(ref_feats.T) + 1e-8 * np.eye(3)
        cov_inv = np.linalg.inv(cov)

        std_scores = {}
        for name, gen_fn in ALL_FAMILIES.items():
            feats = feature_ensemble_std(gen_fn, N, REPS, SEED_BASE + 5000)
            std_scores[name] = mahalanobis_score(feats.mean(axis=0), mu, cov_inv)

        # Find nearest non-Lor4D standard family
        non_lor4d = {k: v for k, v in std_scores.items() if k != "Lor4D"}
        nearest_name = min(non_lor4d, key=non_lor4d.get)
        nearest_score = non_lor4d[nearest_name]

        h_c = None
        h_c_score = None
        for H in HUBBLE_FINE:
            ds_feats = ds4d_feature_ensemble(N, H, REPS, SEED_BASE + 8000)
            ds_score = mahalanobis_score(ds_feats.mean(axis=0), mu, cov_inv)

            all_scores = {**std_scores, f"dS4D_H{H}": ds_score}
            sorted_all = sorted(all_scores.items(), key=lambda x: x[1])
            names = [n for n, _ in sorted_all]
            ds_rank = names.index(f"dS4D_H{H}") + 1

            if ds_rank <= 2:
                h_c = H
                h_c_score = ds_score

        if h_c is not None:
            report.append(
                f"| {N} | {h_c} | {h_c_score:.4f} | {nearest_name} | {nearest_score:.4f} |"
            )
        else:
            report.append(
                f"| {N} | < {HUBBLE_FINE[0]} | — | {nearest_name} | {nearest_score:.4f} |"
            )

    report.append("")

    # ================================================================
    # Part 3: Scaling of S_MD(dS4D_H) with N at fixed H
    # ================================================================
    report.append("=" * 72)
    report.append("# Part 3: S_MD Scaling with N at Fixed H")
    report.append("=" * 72)
    report.append("")
    report.append("Fit: S_MD(N) ~ A · N^β for each H")
    report.append("")

    from scipy.stats import linregress

    for H in [0.01, 0.05, 0.1, 0.3, 0.5, 1.0]:
        scores = []
        Ns_used = []
        for N in N_GRID:
            if N < 28:
                continue  # skip noisy small N
            ref_feats = feature_ensemble_std(ALL_FAMILIES["Lor4D"], N, REPS, SEED_BASE + 9000)
            mu = ref_feats.mean(axis=0)
            cov = np.cov(ref_feats.T) + 1e-8 * np.eye(3)
            cov_inv = np.linalg.inv(cov)
            ds_feats = ds4d_feature_ensemble(N, H, REPS, SEED_BASE + 8000)
            sc = mahalanobis_score(ds_feats.mean(axis=0), mu, cov_inv)
            scores.append(sc)
            Ns_used.append(N)

        if len(scores) >= 3:
            log_n = np.log(Ns_used)
            log_s = np.log(np.maximum(scores, 1e-12))
            res = linregress(log_n, log_s)
            report.append(f"H = {H}: β = {res.slope:.3f} ± {res.stderr:.3f}, "
                          f"R² = {res.rvalue**2:.3f}, "
                          f"range S_MD = [{min(scores):.4f}, {max(scores):.4f}]")

    report.append("")

    elapsed = time.time() - t0
    report.append(f"---\nTotal runtime: {elapsed:.0f}s")

    out_path = OUT_DIR / "experiment_de_sitter_extended.md"
    out_path.write_text("\n".join(report), encoding="utf-8")
    print(f"\nReport written to: {out_path}")
    print(f"Total runtime: {elapsed:.0f}s")


if __name__ == "__main__":
    main()
