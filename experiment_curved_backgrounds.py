"""
Curved background experiment: FLRW matter-dominated + Schwarzschild weak-field.

Tests whether the layered screening architecture (S_MD ranking)
survives non-de-Sitter curvature backgrounds.

FLRW_matter: kappa ∈ {0.0, 0.3, 1.0, 3.0}
  kappa→0: flat Minkowski;  H₀ = 2κ/3 ≈ 0.2 – 2.0
Schwarzschild weak-field: phi0 ∈ {0.0, 0.01, 0.05, 0.1}
  phi0 = M/(2R);  phi0→0: flat;  phi0=0.1: moderate compactness

N ∈ {64, 128, 256},  15 reps per condition.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import linregress

from generators import Poset, transitive_closure
from curvature_backgrounds import (
    sprinkle_flrw_matter_diamond,
    poset_from_flrw_matter_points,
    sprinkle_schwarzschild_diamond,
    poset_from_schwarzschild_points,
)
from expanded_family_robustness import (
    ALL_FAMILIES,
    compute_features,
    mahalanobis_score,
)

SEED_BASE = 42
REPS = 15
N_GRID = [64, 128, 256]
KAPPA_GRID = [0.0, 0.3, 1.0, 3.0]
PHI0_GRID = [0.0, 0.01, 0.05, 0.1]

OUT_DIR = Path("outputs_carlip")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def feature_ensemble(gen_fn, N, reps, seed_base):
    feats = []
    for r in range(reps):
        poset = gen_fn(N, seed=seed_base + r)
        f = compute_features(poset, N)
        feats.append(f)
    return np.array(feats)


def curved_feature_ensemble(sprinkle_fn, poset_fn, N, param, reps, seed_base):
    feats = []
    for r in range(reps):
        pts = sprinkle_fn(N, d_spatial=3, **{param_name(sprinkle_fn): param}, seed=seed_base + r)
        poset = poset_fn(pts, param)
        f = compute_features(poset, N)
        feats.append(f)
    return np.array(feats)


def param_name(fn):
    if "flrw" in fn.__name__:
        return "kappa"
    return "phi0"


def main():
    if hasattr(sys.stdout, "reconfigure"):
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    t0 = time.time()
    report = ["# Curved Background Experiment: FLRW + Schwarzschild\n"]
    report.append(f"N grid: {N_GRID}")
    report.append(f"FLRW kappa grid: {KAPPA_GRID}")
    report.append(f"Schwarz phi0 grid: {PHI0_GRID}")
    report.append(f"REPS: {REPS}\n")

    # ── build reference per N (same as all other experiments) ─────
    cached_ref = {}
    for N in N_GRID:
        tN = time.time()
        print(f"[Ref] N={N} — building 25-family reference ({REPS} reps) ...", flush=True)
        ref_feats = feature_ensemble(ALL_FAMILIES["Lor4D"], N, REPS, SEED_BASE + 9000)
        mu = ref_feats.mean(axis=0)
        cov = np.cov(ref_feats.T) + 1e-8 * np.eye(3)
        cov_inv = np.linalg.inv(cov)

        std_scores = {}
        for name, gen_fn in ALL_FAMILIES.items():
            feats = feature_ensemble(gen_fn, N, REPS, SEED_BASE + 5000)
            std_scores[name] = mahalanobis_score(feats.mean(axis=0), mu, cov_inv)

        cached_ref[N] = (mu, cov_inv, std_scores)
        print(f"  N={N} reference done in {time.time() - tN:.0f}s", flush=True)

    # ================================================================
    #  FLRW matter-dominated
    # ================================================================
    report.append("=" * 72)
    report.append("# FLRW Matter-Dominated (kappa scan)")
    report.append("=" * 72)
    report.append("")
    report.append("| N | kappa | H_0=2k/3 | S_MD(FLRW) | S_MD(Lor4D) | Delta | Lor4D_rank | FLRW_rank | total |")
    report.append("|---|-------|----------|------------|-------------|-------|------------|-----------|-------|")

    for N in N_GRID:
        mu, cov_inv, std_scores = cached_ref[N]
        lor4d_score = std_scores["Lor4D"]

        for kappa in KAPPA_GRID:
            feats = []
            for r in range(REPS):
                pts = sprinkle_flrw_matter_diamond(N, d_spatial=3, kappa=kappa, seed=SEED_BASE + 7000 + r)
                poset = poset_from_flrw_matter_points(pts, kappa)
                f = compute_features(poset, N)
                feats.append(f)
            feats = np.array(feats)
            score = mahalanobis_score(feats.mean(axis=0), mu, cov_inv)

            label = f"FLRW_k{kappa}"
            all_scores = {**std_scores, label: score}
            sorted_all = sorted(all_scores.items(), key=lambda x: x[1])
            names = [n for n, _ in sorted_all]
            lor4d_rank = names.index("Lor4D") + 1
            flrw_rank = names.index(label) + 1
            delta = score - lor4d_score
            h0 = 2.0 * kappa / 3.0

            report.append(
                f"| {N} | {kappa} | {h0:.2f} | {score:.4f} | {lor4d_score:.4f} | "
                f"{delta:+.4f} | #{lor4d_rank} | #{flrw_rank} | {len(all_scores)} |"
            )
        print(f"  FLRW N={N} done", flush=True)

    report.append("")

    # ================================================================
    #  Schwarzschild weak-field
    # ================================================================
    report.append("=" * 72)
    report.append("# Schwarzschild Weak-Field (phi0 scan)")
    report.append("=" * 72)
    report.append("")
    report.append("| N | phi0 | S_MD(Schw) | S_MD(Lor4D) | Delta | Lor4D_rank | Schw_rank | total |")
    report.append("|---|------|------------|-------------|-------|------------|-----------|-------|")

    for N in N_GRID:
        mu, cov_inv, std_scores = cached_ref[N]
        lor4d_score = std_scores["Lor4D"]

        for phi0 in PHI0_GRID:
            feats = []
            for r in range(REPS):
                pts = sprinkle_schwarzschild_diamond(N, d_spatial=3, phi0=phi0, seed=SEED_BASE + 6000 + r)
                poset = poset_from_schwarzschild_points(pts, phi0)
                f = compute_features(poset, N)
                feats.append(f)
            feats = np.array(feats)
            score = mahalanobis_score(feats.mean(axis=0), mu, cov_inv)

            label = f"Schw_p{phi0}"
            all_scores = {**std_scores, label: score}
            sorted_all = sorted(all_scores.items(), key=lambda x: x[1])
            names = [n for n, _ in sorted_all]
            lor4d_rank = names.index("Lor4D") + 1
            schw_rank = names.index(label) + 1
            delta = score - lor4d_score

            report.append(
                f"| {N} | {phi0} | {score:.4f} | {lor4d_score:.4f} | "
                f"{delta:+.4f} | #{lor4d_rank} | #{schw_rank} | {len(all_scores)} |"
            )
        print(f"  Schwarzschild N={N} done", flush=True)

    report.append("")

    # ================================================================
    #  Physical interpretation summary
    # ================================================================
    report.append("=" * 72)
    report.append("# Summary")
    report.append("=" * 72)
    report.append("")
    report.append("FLRW curvature scale: H_0 = 2*kappa/3")
    report.append("  kappa=0.3 → H_0≈0.2 (mild, like slow expansion)")
    report.append("  kappa=1.0 → H_0≈0.67 (moderate)")
    report.append("  kappa=3.0 → H_0≈2.0 (strong, Hubble radius ~ diamond size)")
    report.append("")
    report.append("Schwarzschild compactness: phi0 = M/(2R)")
    report.append("  phi0=0.01 → very weak field")
    report.append("  phi0=0.05 → moderate (c_eff deviation ~10%)")
    report.append("  phi0=0.10 → strong weak-field limit (c_eff deviation ~20%)")
    report.append("")

    elapsed = time.time() - t0
    report.append(f"\nTotal elapsed: {elapsed:.0f}s")
    print(f"\nTotal elapsed: {elapsed:.0f}s", flush=True)

    out_path = OUT_DIR / "experiment_curved_backgrounds.md"
    out_path.write_text("\n".join(report), encoding="utf-8")
    print(f"Report saved to {out_path}", flush=True)


if __name__ == "__main__":
    main()
