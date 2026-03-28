"""
Curved background large-N extension (includes N=512).

Goal: strengthen FLRW/Schwarzschild robustness evidence and provide
Table-9-ready numbers for manuscript.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

from curvature_backgrounds import (
    sprinkle_flrw_matter_diamond,
    poset_from_flrw_matter_points,
    sprinkle_schwarzschild_diamond,
    poset_from_schwarzschild_points,
)
from expanded_family_robustness import ALL_FAMILIES, compute_features, mahalanobis_score

SEED_BASE = 42
REPS = 15
N_GRID = [256, 512]
KAPPA_GRID = [0.0, 0.3, 1.0, 3.0]
PHI0_GRID = [0.0, 0.01, 0.05, 0.1]

OUT_DIR = Path("outputs_carlip")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def feature_ensemble(gen_fn, N, reps, seed_base):
    feats = []
    for r in range(reps):
        poset = gen_fn(N, seed=seed_base + r)
        feats.append(compute_features(poset, N))
    return np.array(feats)


def main():
    reconfigure = getattr(sys.stdout, "reconfigure", None)
    if callable(reconfigure):
        try:
            reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    t0 = time.time()
    report = ["# Curved Background Large-N Extension (N=256,512)", ""]
    report.append(f"REPS={REPS}, N={N_GRID}")
    report.append(f"FLRW kappa={KAPPA_GRID}, Schwarz phi0={PHI0_GRID}")
    report.append("")

    # build Lor4D reference + 25-family baseline per N
    cached = {}
    for N in N_GRID:
        tN = time.time()
        print(f"[Ref] N={N} building 25-family baseline ...", flush=True)
        ref_feats = feature_ensemble(ALL_FAMILIES["Lor4D"], N, REPS, SEED_BASE + 9000)
        mu = ref_feats.mean(axis=0)
        cov = np.cov(ref_feats.T) + 1e-8 * np.eye(3)
        cov_inv = np.linalg.inv(cov)

        std_scores = {}
        for name, gen_fn in ALL_FAMILIES.items():
            feats = feature_ensemble(gen_fn, N, REPS, SEED_BASE + 5000)
            std_scores[name] = mahalanobis_score(feats.mean(axis=0), mu, cov_inv)

        cached[N] = (mu, cov_inv, std_scores)
        print(f"  N={N} baseline done in {time.time() - tN:.0f}s", flush=True)

    report.append("## FLRW (matter-dominated)")
    report.append("")
    report.append("| N | kappa | H0=2k/3 | S_MD(FLRW) | S_MD(Lor4D) | Delta | FLRW rank | Lor4D rank |")
    report.append("|---|---|---|---|---|---|---|---|")

    for N in N_GRID:
        mu, cov_inv, std_scores = cached[N]
        lor4d = std_scores["Lor4D"]
        for kappa in KAPPA_GRID:
            feats = []
            for r in range(REPS):
                pts = sprinkle_flrw_matter_diamond(N, d_spatial=3, kappa=kappa, seed=SEED_BASE + 7000 + r)
                p = poset_from_flrw_matter_points(pts, kappa)
                feats.append(compute_features(p, N))
            score = mahalanobis_score(np.array(feats).mean(axis=0), mu, cov_inv)
            label = f"FLRW_k{kappa}"
            all_scores = {**std_scores, label: score}
            names = [n for n, _ in sorted(all_scores.items(), key=lambda x: x[1])]
            flrw_rank = names.index(label) + 1
            lor_rank = names.index("Lor4D") + 1
            delta = score - lor4d
            h0 = 2.0 * kappa / 3.0
            report.append(f"| {N} | {kappa} | {h0:.2f} | {score:.4f} | {lor4d:.4f} | {delta:+.4f} | #{flrw_rank} | #{lor_rank} |")
        print(f"  FLRW N={N} done", flush=True)

    report.append("")
    report.append("## Schwarzschild (weak-field proxy)")
    report.append("")
    report.append("| N | phi0 | S_MD(Schw) | S_MD(Lor4D) | Delta | Schw rank | Lor4D rank |")
    report.append("|---|---|---|---|---|---|---|")

    for N in N_GRID:
        mu, cov_inv, std_scores = cached[N]
        lor4d = std_scores["Lor4D"]
        for phi0 in PHI0_GRID:
            feats = []
            for r in range(REPS):
                pts = sprinkle_schwarzschild_diamond(N, d_spatial=3, phi0=phi0, seed=SEED_BASE + 6000 + r)
                p = poset_from_schwarzschild_points(pts, phi0)
                feats.append(compute_features(p, N))
            score = mahalanobis_score(np.array(feats).mean(axis=0), mu, cov_inv)
            label = f"Schw_p{phi0}"
            all_scores = {**std_scores, label: score}
            names = [n for n, _ in sorted(all_scores.items(), key=lambda x: x[1])]
            schw_rank = names.index(label) + 1
            lor_rank = names.index("Lor4D") + 1
            delta = score - lor4d
            report.append(f"| {N} | {phi0} | {score:.4f} | {lor4d:.4f} | {delta:+.4f} | #{schw_rank} | #{lor_rank} |")
        print(f"  Schwarzschild N={N} done", flush=True)

    elapsed = time.time() - t0
    report.append("")
    report.append(f"Total elapsed: {elapsed:.0f}s")

    out = OUT_DIR / "experiment_curved_backgrounds_largeN.md"
    out.write_text("\n".join(report), encoding="utf-8")
    print(f"\nTotal elapsed: {elapsed:.0f}s", flush=True)
    print(f"Report saved to {out}", flush=True)


if __name__ == "__main__":
    main()
