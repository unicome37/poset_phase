"""
N=16 Mahalanobis CV Failure — Root Cause Diagnosis
====================================================

Goal: Understand WHY Mahalanobis fails at N=16 under cross-validation.

Hypotheses:
  H1: Small-sample covariance estimation instability (8 train → 3×3 Σ)
  H2: Condition number explosion → Σ⁻¹ noise amplification
  H3: μ(N=16) itself not well-separated from competitors
  H4: Specific competitor(s) intrude into Lor4D's covariance ellipsoid

Diagnostics:
  1. Eigenvalue spectrum of Σ at small vs large N
  2. Condition number κ(Σ) vs N
  3. μ shift between train/test splits
  4. Which competitor wins when Lor4D fails
  5. Regularized Mahalanobis: Σ + λI (shrinkage estimator)
"""
from __future__ import annotations

import time
from collections import defaultdict
from pathlib import Path

import numpy as np

from bd_action import count_intervals_fast
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
from unified_functional import compute_xi_dim


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


def compute_features(poset: Poset, N: int) -> np.ndarray:
    counts = count_intervals_fast(poset, k_max=5)
    C0 = counts.get(0)
    C1 = counts.get(1)
    c1_c0 = C1 / max(1, C0)
    _, d_eff = compute_xi_dim(poset)
    aw = max_antichain_width(poset)
    width_ratio = aw / max(1, N)
    return np.array([d_eff, c1_c0, width_ratio])


def main():
    SEED_BASES = [42, 137, 271, 500, 777, 1001, 2023, 3141, 5000, 8888]
    N_ALL = [16, 20, 28, 64, 128]  # focused set: problem zone + reference
    N_VALUES = N_ALL
    REPS = 30
    N_FOLDS = 5

    report = []
    report.append("# N=16 Mahalanobis CV Failure — Root Cause Diagnosis\n")

    t0 = time.time()

    # Collect data for all seeds
    all_data = {}
    for seed_base in SEED_BASES:
        data = defaultdict(lambda: defaultdict(list))
        for fam_name, gen_fn in FAMILIES.items():
            for N in N_VALUES:
                for rep in range(REPS):
                    seed = seed_base + hash(fam_name) % 10000 + N * 100 + rep
                    seed = seed % (2**31)
                    try:
                        poset = gen_fn(N, seed=seed)
                        feat = compute_features(poset, N)
                        data[N][fam_name].append(feat)
                    except Exception:
                        pass
        all_data[seed_base] = data

    elapsed_gen = time.time() - t0
    print(f"Data generation: {elapsed_gen:.1f}s")

    # ══════════════════════════════════════════════════════════
    # Diagnostic 1: Eigenvalue spectrum & condition number
    # ══════════════════════════════════════════════════════════
    report.append("## 1. Covariance Eigenvalue Spectrum\n")
    report.append("| N | λ₁ | λ₂ | λ₃ | κ(Σ) = λ_max/λ_min | det(Σ) |")
    report.append("|---|:-:|:-:|:-:|:-:|:-:|")

    # Average over all seeds
    for N in N_VALUES:
        eig_lists = [[], [], []]
        conds = []
        dets = []
        for seed_base in SEED_BASES:
            lor4d = np.array(all_data[seed_base][N]["Lor4D"])
            if len(lor4d) < 4:
                continue
            cov = np.cov(lor4d.T)
            eigvals = np.sort(np.linalg.eigvalsh(cov))[::-1]
            for k in range(3):
                eig_lists[k].append(eigvals[k])
            conds.append(eigvals[0] / max(eigvals[-1], 1e-15))
            dets.append(np.linalg.det(cov))

        report.append(f"| {N} | {np.mean(eig_lists[0]):.4f} | "
                      f"{np.mean(eig_lists[1]):.4f} | "
                      f"{np.mean(eig_lists[2]):.6f} | "
                      f"{np.mean(conds):.1f} | {np.mean(dets):.2e} |")
    report.append("")

    # ══════════════════════════════════════════════════════════
    # Diagnostic 2: Train/test μ shift at N=16
    # ══════════════════════════════════════════════════════════
    report.append("## 2. Train/Test μ Shift at N=16\n")
    report.append("How much does the centroid shift when we remove 20% of data?\n")
    report.append("| Seed | Full μ [d, c, w] | Max fold Δμ | Max fold ‖Δμ‖/‖σ‖ |")
    report.append("|------|:-:|:-:|:-:|")

    for seed_base in SEED_BASES:
        lor4d = np.array(all_data[seed_base][16]["Lor4D"])
        n = len(lor4d)
        mu_full = np.mean(lor4d, axis=0)
        sigma_full = np.std(lor4d, axis=0)
        rng = np.random.RandomState(seed_base + 16)
        indices = rng.permutation(n)
        fold_size = n // N_FOLDS

        max_shift = 0
        max_rel = 0
        for fold in range(N_FOLDS):
            test_start = fold * fold_size
            test_end = test_start + fold_size if fold < N_FOLDS - 1 else n
            test_idx = set(indices[test_start:test_end])
            train_idx = [i for i in range(n) if i not in test_idx]
            mu_train = np.mean(lor4d[train_idx], axis=0)
            shift = np.abs(mu_train - mu_full)
            max_shift = max(max_shift, np.max(shift))
            rel = np.max(shift / np.maximum(sigma_full, 1e-6))
            max_rel = max(max_rel, rel)

        report.append(f"| {seed_base} | [{mu_full[0]:.3f}, {mu_full[1]:.3f}, {mu_full[2]:.3f}] | "
                      f"{max_shift:.4f} | {max_rel:.3f} |")
    report.append("")

    # ══════════════════════════════════════════════════════════
    # Diagnostic 3: Who wins when Lor4D fails?
    # ══════════════════════════════════════════════════════════
    report.append("## 3. Competitor Identity at Failure Points\n")
    report.append("5-fold CV at N=16 and N=20: who beats Lor4D?\n")
    report.append("| Seed | N | Fold | Winner | Winner S_M | Lor4D S_M | Gap |")
    report.append("|------|---|------|--------|:-:|:-:|:-:|")

    failure_count = 0
    winner_census = defaultdict(int)

    for seed_base in SEED_BASES:
        for N in [16, 20]:
            lor4d_all = np.array(all_data[seed_base][N]["Lor4D"])
            n_lor = len(lor4d_all)
            rng = np.random.RandomState(seed_base + N)
            indices = rng.permutation(n_lor)
            fold_size = n_lor // N_FOLDS

            for fold in range(N_FOLDS):
                test_start = fold * fold_size
                test_end = test_start + fold_size if fold < N_FOLDS - 1 else n_lor
                test_idx = set(indices[test_start:test_end])
                train_idx = [i for i in range(n_lor) if i not in test_idx]

                lor4d_train = lor4d_all[train_idx]
                mu = np.mean(lor4d_train, axis=0)
                cov = np.cov(lor4d_train.T)
                cov_inv = np.linalg.inv(cov + 1e-12 * np.eye(3))

                fam_scores = {}
                for fam in FAMILIES:
                    if fam == "Lor4D":
                        feats = lor4d_all[sorted(test_idx)]
                    else:
                        feats = all_data[seed_base][N][fam]
                        if not feats:
                            continue
                        feats = np.array(feats)
                    scores = [float((f - mu) @ cov_inv @ (f - mu)) for f in feats]
                    fam_scores[fam] = np.mean(scores)

                ranked = sorted(fam_scores.items(), key=lambda x: x[1])
                lor_rank = next(i+1 for i, (f, _) in enumerate(ranked) if f == "Lor4D")

                if lor_rank > 1:
                    failure_count += 1
                    winner = ranked[0][0]
                    winner_census[winner] += 1
                    report.append(f"| {seed_base} | {N} | {fold+1} | {winner} | "
                                  f"{fam_scores[winner]:.2f} | {fam_scores['Lor4D']:.2f} | "
                                  f"{fam_scores['Lor4D'] - fam_scores[winner]:.2f} |")
    report.append("")
    report.append(f"Total failures: {failure_count}")
    report.append(f"Winner census: {dict(winner_census)}\n")

    # ══════════════════════════════════════════════════════════
    # Diagnostic 4: Shrinkage estimator (Ledoit-Wolf style)
    # ══════════════════════════════════════════════════════════
    report.append("## 4. Regularized Mahalanobis (Shrinkage)\n")
    report.append("Replace Σ with (1-α)Σ + α·tr(Σ)/3·I (shrinkage toward spherical)\n")
    report.append("| α | N=16 CV #1 rate | N=20 CV #1 rate | N=28+ CV #1 rate |")
    report.append("|---|:-:|:-:|:-:|")

    for alpha in [0.0, 0.1, 0.2, 0.5, 1.0]:
        counts = {16: [0, 0], 20: [0, 0], 28: [0, 0]}
        for seed_base in SEED_BASES:
            for N in [16, 20, 28]:
                lor4d_all = np.array(all_data[seed_base][N]["Lor4D"])
                n_lor = len(lor4d_all)
                rng = np.random.RandomState(seed_base + N)
                indices = rng.permutation(n_lor)
                fold_size = n_lor // N_FOLDS

                for fold in range(N_FOLDS):
                    test_start = fold * fold_size
                    test_end = test_start + fold_size if fold < N_FOLDS - 1 else n_lor
                    test_idx = set(indices[test_start:test_end])
                    train_idx = [i for i in range(n_lor) if i not in test_idx]

                    lor4d_train = lor4d_all[train_idx]
                    mu = np.mean(lor4d_train, axis=0)
                    cov = np.cov(lor4d_train.T)
                    # Shrinkage
                    target = np.trace(cov) / 3.0 * np.eye(3)
                    cov_reg = (1 - alpha) * cov + alpha * target
                    cov_inv = np.linalg.inv(cov_reg + 1e-12 * np.eye(3))

                    fam_scores = {}
                    for fam in FAMILIES:
                        if fam == "Lor4D":
                            feats = lor4d_all[sorted(test_idx)]
                        else:
                            feats = all_data[seed_base][N][fam]
                            if not feats:
                                continue
                            feats = np.array(feats)
                        scores = [float((f - mu) @ cov_inv @ (f - mu)) for f in feats]
                        fam_scores[fam] = np.mean(scores)

                    ranked = sorted(fam_scores.items(), key=lambda x: x[1])
                    lor_rank = next(i+1 for i, (f, _) in enumerate(ranked) if f == "Lor4D")
                    bucket = 28 if N >= 28 else N
                    counts[bucket][1] += 1
                    if lor_rank == 1:
                        counts[bucket][0] += 1

        r16 = f"{counts[16][0]}/{counts[16][1]}" if counts[16][1] > 0 else "N/A"
        r20 = f"{counts[20][0]}/{counts[20][1]}" if counts[20][1] > 0 else "N/A"
        r28 = f"{counts[28][0]}/{counts[28][1]}" if counts[28][1] > 0 else "N/A"
        report.append(f"| {alpha} | {r16} | {r20} | {r28} |")
    report.append("")

    # ══════════════════════════════════════════════════════════
    # Diagnostic 5: Minimum train size needed
    # ══════════════════════════════════════════════════════════
    report.append("## 5. Minimum Training Size for Stable Mahalanobis at N=16\n")
    report.append("Vary train fraction from 50% to 95%.\n")

    train_fracs = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    report.append("| Train fraction | n_train | CV #1 rate (N=16) |")
    report.append("|:-:|:-:|:-:|")

    for frac in train_fracs:
        total_1 = 0
        total_tests = 0
        for seed_base in SEED_BASES:
            lor4d_all = np.array(all_data[seed_base][16]["Lor4D"])
            n_lor = len(lor4d_all)
            n_train = int(frac * n_lor)
            n_test = n_lor - n_train
            if n_test < 2 or n_train < 4:
                continue

            # 10 random splits
            for split_rep in range(10):
                rng = np.random.RandomState(seed_base + split_rep * 1000)
                perm = rng.permutation(n_lor)
                train_idx = perm[:n_train]
                test_idx = perm[n_train:]

                mu = np.mean(lor4d_all[train_idx], axis=0)
                cov = np.cov(lor4d_all[train_idx].T)
                cov_inv = np.linalg.inv(cov + 1e-12 * np.eye(3))

                fam_scores = {}
                for fam in FAMILIES:
                    if fam == "Lor4D":
                        feats = lor4d_all[test_idx]
                    else:
                        feats = all_data[seed_base][16][fam]
                        if not feats:
                            continue
                        feats = np.array(feats)
                    scores = [float((f - mu) @ cov_inv @ (f - mu)) for f in feats]
                    fam_scores[fam] = np.mean(scores)

                ranked = sorted(fam_scores.items(), key=lambda x: x[1])
                lor_rank = next(i+1 for i, (f, _) in enumerate(ranked) if f == "Lor4D")
                total_tests += 1
                if lor_rank == 1:
                    total_1 += 1

        rate = f"{total_1}/{total_tests}" if total_tests > 0 else "N/A"
        pct = f"({100*total_1/total_tests:.0f}%)" if total_tests > 0 else ""
        report.append(f"| {frac:.0%} | ~{int(frac*40)} | {rate} {pct} |")
    report.append("")

    # ══════════════════════════════════════════════════════════
    # Synthesis
    # ══════════════════════════════════════════════════════════
    report.append("## 6. Root Cause Synthesis\n")
    report.append("The analysis above identifies the mechanism behind N=16 CV failures.\n")

    elapsed = time.time() - t0
    print(f"\nTotal: {elapsed:.1f}s")

    out = "\n".join(report)
    outdir = Path("outputs_carlip")
    outdir.mkdir(exist_ok=True)
    md_path = outdir / "n16_mahalanobis_cv_rootcause.md"
    md_path.write_text(out, encoding="utf-8")
    print(f"Saved: {md_path}")
    print("\n" + out)


if __name__ == "__main__":
    main()
