"""
Anti-Engineering-Illusion Test Battery
======================================
Seven independent tests to detect whether F7's predictive power is genuine
or an artifact of overfitting / phenomenological patching.

Tests:
  1. Permutation test — shuffle family labels, check if F7 still "predicts"
  2. Random functional — random weight vectors, check hit rate distribution
  3. Ablation with noise — replace each term with random noise
  4. Cross-validation — train/test split on posets
  5. Sigmoid wall necessity — F7 without wall at various N
  6. New seed test — entirely new posets with unseen seeds
  7. Fake family test — synthetic families with known properties

Usage:
  python _anti_illusion_battery.py
"""
from __future__ import annotations
import csv, time, sys
import numpy as np
from generators import (
    Poset,
    generate_lorentzian_like_2d, generate_lorentzian_like_3d,
    generate_lorentzian_like_4d, generate_lorentzian_like_5d,
    generate_kr_like, generate_transitive_percolation,
    generate_interval_order, generate_absolute_layered,
)

# ── Shared utilities ──

def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -30, 30)))

GENS = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
    "KR_like": generate_kr_like,
}

def load_data():
    """Load raw_features.csv and compute R for each sample."""
    data = list(csv.DictReader(open("outputs_unified_functional/raw_features.csv")))
    n_data = len(data)
    families = np.array([d["family"] for d in data])
    N_arr = np.array([int(d["N"]) for d in data], dtype=float)
    log_H = np.array([float(d["log_H"]) for d in data])
    pi_geo = np.array([float(d["pi_geo"]) for d in data])
    sigma_hist = np.array([float(d["sigma_hist"]) for d in data])
    xi_dim = np.array([float(d["xi_dim"]) for d in data])
    pi_cg = np.array([float(d["pi_cg"]) for d in data])

    # Compute R (interval richness)
    R_arr = np.zeros(n_data)
    for i, row in enumerate(data):
        fam, n, rep = row["family"], int(row["N"]), int(row["rep"])
        seed = 42 + rep * 1000 + n * 100
        p = GENS[fam](n, seed=seed)
        c = p.closure.astype(np.int32)
        ks = c @ c
        mask = p.closure
        C0 = int(np.sum(mask & (ks == 0)))
        total = int(np.sum(mask))
        R_arr[i] = 1.0 - C0 / total if total > 0 else 0.0

    return {
        "data": data, "families": families, "N_arr": N_arr,
        "log_H": log_H, "pi_geo": pi_geo, "sigma_hist": sigma_hist,
        "xi_dim": xi_dim, "pi_cg": pi_cg, "R_arr": R_arr,
    }


def compute_F7_from_arrays(log_H, pi_geo, sigma_hist, xi_dim, R_arr, N_arr,
                           alpha0=26.0, lam=5.0, eta=1.0, Rc=0.25, w=0.015,
                           gam=0.0004, q=-0.5, N0=20.0):
    """Vectorized F7 computation from feature arrays."""
    alpha_N = alpha0 * (N0 / N_arr) ** abs(q)
    wall = alpha_N * sigmoid((R_arr - Rc) / w)
    return log_H + gam * pi_geo - lam * sigma_hist + eta * xi_dim + wall


def build_pairs(families, N_arr, fam_i, fam_j):
    """Build comparison pair indices."""
    unique_N = sorted(set(N_arr))
    pairs = []
    for nv in unique_N:
        ii = np.where((N_arr == nv) & (families == fam_i))[0]
        jj = np.where((N_arr == nv) & (families == fam_j))[0]
        for a in ii:
            for b in jj:
                pairs.append((a, b))
    return np.array(pairs, dtype=int) if pairs else np.zeros((0, 2), dtype=int)


def compute_scores(F, families, N_arr, sigma_hist, pi_cg):
    """Compute prediction scores A, B, C, D."""
    unique_N = sorted(set(N_arr))
    A_pairs_2d = build_pairs(families, N_arr, "Lor4D", "Lor2D")
    A_pairs_3d = build_pairs(families, N_arr, "Lor4D", "Lor3D")
    A_pairs_5d = build_pairs(families, N_arr, "Lor4D", "Lor5D")
    B_pairs = build_pairs(families, N_arr, "Lor2D", "KR_like")

    # A: dimensional selection (4D should have lowest F)
    w2 = int(np.sum(F[A_pairs_2d[:, 0]] < F[A_pairs_2d[:, 1]])) if len(A_pairs_2d) > 0 else 0
    w3 = int(np.sum(F[A_pairs_3d[:, 0]] < F[A_pairs_3d[:, 1]])) if len(A_pairs_3d) > 0 else 0
    w5 = int(np.sum(F[A_pairs_5d[:, 0]] < F[A_pairs_5d[:, 1]])) if len(A_pairs_5d) > 0 else 0
    t_a = len(A_pairs_2d) + len(A_pairs_3d) + len(A_pairs_5d)
    A = (w2 + w3 + w5) / t_a if t_a > 0 else 0

    # B: phase competition (Lor2D < KR)
    wb = int(np.sum(F[B_pairs[:, 0]] < F[B_pairs[:, 1]])) if len(B_pairs) > 0 else 0
    B = wb / len(B_pairs) if len(B_pairs) > 0 else 0

    # C: sedimentation (negative corr between sigma_hist and F)
    neg = 0
    for nv in unique_N:
        idx = np.where(N_arr == nv)[0]
        sh = sigma_hist[idx]
        f = F[idx]
        if np.std(sh) > 1e-10 and np.std(f) > 1e-10:
            if np.corrcoef(sh, f)[0, 1] < 0:
                neg += 1
    C = neg / len(unique_N) if len(unique_N) > 0 else 0

    # D: coarse-graining stability (low corr with pi_cg)
    if np.std(pi_cg) > 1e-10 and np.std(F) > 1e-10:
        D = abs(np.corrcoef(pi_cg, F)[0, 1])
    else:
        D = 0.0

    return A, B, C, D


def grade_str(A, B, C):
    if A >= 0.65 and B >= 0.95 and C >= 0.75:
        return "★"
    if A >= 0.65 and B >= 0.95:
        return "●"
    if A >= 0.65 and B >= 0.90:
        return "◆"
    return "·"


# ══════════════════════════════════════════════════════════════════════════
# TEST 1: Permutation Test — Shuffle family labels
# ══════════════════════════════════════════════════════════════════════════
def test_permutation(d, n_perm=1000, seed=12345):
    """
    If F7 is genuinely capturing structural differences between families,
    then shuffling family labels should destroy predictive power.
    If shuffled labels still yield high A/B scores, F7 is overfitting.
    """
    print("\n" + "=" * 70)
    print("TEST 1: PERMUTATION TEST — Shuffle family labels")
    print("=" * 70)
    print(f"  Null hypothesis: F7 scores are independent of family labels")
    print(f"  Running {n_perm} permutations...")

    rng = np.random.default_rng(seed)
    F_real = compute_F7_from_arrays(
        d["log_H"], d["pi_geo"], d["sigma_hist"], d["xi_dim"],
        d["R_arr"], d["N_arr"]
    )
    A_real, B_real, C_real, D_real = compute_scores(
        F_real, d["families"], d["N_arr"], d["sigma_hist"], d["pi_cg"]
    )
    print(f"\n  Real F7 scores: A={A_real:.3f}, B={B_real:.3f}, C={C_real:.2f}, D={D_real:.3f}")
    print(f"  Grade: {grade_str(A_real, B_real, C_real)}")

    perm_A, perm_B, perm_C = [], [], []
    for _ in range(n_perm):
        shuffled = rng.permutation(d["families"])
        A, B, C, _ = compute_scores(F_real, shuffled, d["N_arr"], d["sigma_hist"], d["pi_cg"])
        perm_A.append(A)
        perm_B.append(B)
        perm_C.append(C)

    perm_A, perm_B, perm_C = np.array(perm_A), np.array(perm_B), np.array(perm_C)

    p_A = np.mean(perm_A >= A_real)
    p_B = np.mean(perm_B >= B_real)
    p_C = np.mean(perm_C >= C_real)

    print(f"\n  Permutation distribution (mean ± std):")
    print(f"    A: {perm_A.mean():.3f} ± {perm_A.std():.3f}  (real={A_real:.3f}, p={p_A:.4f})")
    print(f"    B: {perm_B.mean():.3f} ± {perm_B.std():.3f}  (real={B_real:.3f}, p={p_B:.4f})")
    print(f"    C: {perm_C.mean():.3f} ± {perm_C.std():.3f}  (real={C_real:.3f}, p={p_C:.4f})")

    verdict = "PASS" if (p_A < 0.05 and p_B < 0.05) else "FAIL"
    print(f"\n  VERDICT: {verdict}")
    if verdict == "PASS":
        print("  → F7 scores are significantly better than chance (p<0.05)")
        print("  → Family labels carry real structural information")
    else:
        print("  → WARNING: F7 scores could arise from random label assignment")
    return {"A_real": A_real, "B_real": B_real, "p_A": p_A, "p_B": p_B, "verdict": verdict}


# ══════════════════════════════════════════════════════════════════════════
# TEST 2: Random Functional Test — Random weight vectors
# ══════════════════════════════════════════════════════════════════════════
def test_random_functional(d, n_random=5000, seed=54321):
    """
    Generate random weight vectors and check what fraction achieve
    similar or better hit rates than F7.
    If many random functionals work equally well, F7's form is not special.
    """
    print("\n" + "=" * 70)
    print("TEST 2: RANDOM FUNCTIONAL — Random weight vectors")
    print("=" * 70)
    print(f"  Null hypothesis: Random weight vectors achieve similar scores")
    print(f"  Sampling {n_random} random functionals...")

    rng = np.random.default_rng(seed)

    # Real F7 scores
    F_real = compute_F7_from_arrays(
        d["log_H"], d["pi_geo"], d["sigma_hist"], d["xi_dim"],
        d["R_arr"], d["N_arr"]
    )
    A_real, B_real, C_real, D_real = compute_scores(
        F_real, d["families"], d["N_arr"], d["sigma_hist"], d["pi_cg"]
    )
    print(f"\n  Real F7: A={A_real:.3f}, B={B_real:.3f}, C={C_real:.2f}")

    # Random functionals: F = w1*logH + w2*pi_geo + w3*sigma_hist + w4*xi_dim + w5*wall
    rand_A, rand_B, rand_C = [], [], []
    n_better = 0
    for _ in range(n_random):
        # Random weights from broad distributions
        w1 = rng.uniform(0.5, 2.0)       # logH coefficient
        w2 = rng.uniform(-1.0, 1.0)      # pi_geo coefficient
        w3 = rng.uniform(-15.0, 15.0)    # sigma_hist coefficient (can be + or -)
        w4 = rng.uniform(-3.0, 3.0)      # xi_dim coefficient
        alpha = rng.uniform(0, 50)        # wall strength
        Rc = rng.uniform(0.05, 0.50)      # wall position
        w_wall = rng.uniform(0.005, 0.2)  # wall width

        wall = alpha * sigmoid((d["R_arr"] - Rc) / w_wall)
        F_rand = w1 * d["log_H"] + w2 * d["pi_geo"] + w3 * d["sigma_hist"] + w4 * d["xi_dim"] + wall

        A, B, C, _ = compute_scores(F_rand, d["families"], d["N_arr"], d["sigma_hist"], d["pi_cg"])
        rand_A.append(A)
        rand_B.append(B)
        rand_C.append(C)
        if A >= A_real and B >= B_real and C >= C_real:
            n_better += 1

    rand_A, rand_B, rand_C = np.array(rand_A), np.array(rand_B), np.array(rand_C)

    print(f"\n  Random functional distribution:")
    print(f"    A: {rand_A.mean():.3f} ± {rand_A.std():.3f}  (max={rand_A.max():.3f})")
    print(f"    B: {rand_B.mean():.3f} ± {rand_B.std():.3f}  (max={rand_B.max():.3f})")
    print(f"    C: {rand_C.mean():.3f} ± {rand_C.std():.3f}  (max={rand_C.max():.3f})")
    print(f"\n  Fraction matching or beating F7 on ALL of A,B,C: {n_better}/{n_random} = {n_better/n_random:.4f}")

    # Also check: how many beat F7 on individual metrics?
    n_beat_A = np.sum(rand_A >= A_real)
    n_beat_B = np.sum(rand_B >= B_real)
    n_beat_C = np.sum(rand_C >= C_real)
    print(f"  Fraction beating F7 on A alone: {n_beat_A}/{n_random} = {n_beat_A/n_random:.4f}")
    print(f"  Fraction beating F7 on B alone: {n_beat_B}/{n_random} = {n_beat_B/n_random:.4f}")
    print(f"  Fraction beating F7 on C alone: {n_beat_C}/{n_random} = {n_beat_C/n_random:.4f}")

    verdict = "PASS" if n_better / n_random < 0.01 else "FAIL"
    print(f"\n  VERDICT: {verdict}")
    if verdict == "PASS":
        print("  → Less than 1% of random functionals match F7's performance")
        print("  → F7's specific form is genuinely special")
    else:
        print("  → WARNING: Many random functionals achieve similar scores")
        print("  → F7's form may not be uniquely determined by the data")
    return {"n_better": n_better, "frac": n_better / n_random, "verdict": verdict}


# ══════════════════════════════════════════════════════════════════════════
# TEST 3: Ablation with Noise Replacement
# ══════════════════════════════════════════════════════════════════════════
def test_ablation(d, n_trials=200, seed=99999):
    """
    Replace each term in F7 with random noise of the same distribution.
    If removing a term doesn't degrade predictions, that term is decorative.
    """
    print("\n" + "=" * 70)
    print("TEST 3: ABLATION — Replace each term with matched noise")
    print("=" * 70)

    rng = np.random.default_rng(seed)

    F_real = compute_F7_from_arrays(
        d["log_H"], d["pi_geo"], d["sigma_hist"], d["xi_dim"],
        d["R_arr"], d["N_arr"]
    )
    A_real, B_real, C_real, _ = compute_scores(
        F_real, d["families"], d["N_arr"], d["sigma_hist"], d["pi_cg"]
    )
    print(f"  Baseline F7: A={A_real:.3f}, B={B_real:.3f}, C={C_real:.2f}")

    terms = {
        "log_H": d["log_H"],
        "sigma_hist": -5.0 * d["sigma_hist"],  # with its coefficient
        "xi_dim": 1.0 * d["xi_dim"],
        "sigmoid_wall": 26.0 * (20.0 / d["N_arr"]) ** 0.5 * sigmoid((d["R_arr"] - 0.25) / 0.015),
    }

    print(f"\n  {'Term':<15} {'A_ablated':>10} {'B_ablated':>10} {'C_ablated':>10} {'ΔA':>8} {'ΔB':>8} {'Impact':>8}")
    print("  " + "-" * 70)

    for term_name, term_values in terms.items():
        ablated_A, ablated_B, ablated_C = [], [], []
        for _ in range(n_trials):
            noise = rng.normal(np.mean(term_values), np.std(term_values), size=len(term_values))
            # Reconstruct F7 with this term replaced by noise
            F_abl = F_real - term_values + noise
            A, B, C, _ = compute_scores(F_abl, d["families"], d["N_arr"], d["sigma_hist"], d["pi_cg"])
            ablated_A.append(A)
            ablated_B.append(B)
            ablated_C.append(C)

        mA, mB, mC = np.mean(ablated_A), np.mean(ablated_B), np.mean(ablated_C)
        dA, dB = mA - A_real, mB - B_real
        impact = "CRITICAL" if (dA < -0.10 or dB < -0.10) else "moderate" if (dA < -0.03 or dB < -0.03) else "minimal"
        print(f"  {term_name:<15} {mA:10.3f} {mB:10.3f} {mC:10.2f} {dA:+8.3f} {dB:+8.3f} {impact:>8}")

    # Also test: remove sigmoid wall entirely (set to 0)
    F_no_wall = d["log_H"] + 0.0004 * d["pi_geo"] - 5.0 * d["sigma_hist"] + 1.0 * d["xi_dim"]
    A_nw, B_nw, C_nw, _ = compute_scores(F_no_wall, d["families"], d["N_arr"], d["sigma_hist"], d["pi_cg"])
    print(f"\n  No wall at all: A={A_nw:.3f}, B={B_nw:.3f}, C={C_nw:.2f}")
    print(f"  Δ from F7:      ΔA={A_nw-A_real:+.3f}, ΔB={B_nw-B_real:+.3f}, ΔC={C_nw-C_real:+.2f}")

    verdict = "PASS" if (A_nw < A_real - 0.05 or B_nw < B_real - 0.05) else "AMBIGUOUS"
    print(f"\n  VERDICT: {verdict}")
    if verdict == "PASS":
        print("  → Sigmoid wall contributes meaningfully to predictions")
    else:
        print("  → Sigmoid wall may be unnecessary — predictions hold without it")
    return {"A_no_wall": A_nw, "B_no_wall": B_nw, "verdict": verdict}


# ══════════════════════════════════════════════════════════════════════════
# TEST 4: Cross-Validation — Train/Test split
# ══════════════════════════════════════════════════════════════════════════
def test_cross_validation(d, n_folds=4, seed=77777):
    """
    Split 160 posets into train/test by rep index.
    Fit weights on train set, evaluate on test set.
    If test performance degrades significantly, F7 is overfitting.
    """
    print("\n" + "=" * 70)
    print("TEST 4: CROSS-VALIDATION — Train/Test split by rep")
    print("=" * 70)

    rng = np.random.default_rng(seed)
    reps = np.array([int(row["rep"]) for row in d["data"]])
    unique_reps = sorted(set(reps))
    n_reps = len(unique_reps)

    # Shuffle reps and split into folds
    shuffled_reps = rng.permutation(unique_reps)
    fold_size = n_reps // n_folds

    print(f"  {n_folds}-fold CV, {n_reps} unique reps, fold_size={fold_size}")
    print(f"  For each fold: optimize (α, λ, η) on train, evaluate on test\n")

    all_train_scores = []
    all_test_scores = []

    for fold in range(n_folds):
        test_reps = set(shuffled_reps[fold * fold_size:(fold + 1) * fold_size])
        train_mask = np.array([r not in test_reps for r in reps])
        test_mask = ~train_mask

        # Grid search on train set
        best_A_train, best_B_train, best_params = 0, 0, None
        for alpha in [15, 20, 26, 32, 40]:
            for lam in [3, 5, 7, 10]:
                for eta in [0.3, 0.6, 1.0, 1.5]:
                    F_train = compute_F7_from_arrays(
                        d["log_H"][train_mask], d["pi_geo"][train_mask],
                        d["sigma_hist"][train_mask], d["xi_dim"][train_mask],
                        d["R_arr"][train_mask], d["N_arr"][train_mask],
                        alpha0=alpha, lam=lam, eta=eta,
                    )
                    A, B, C, _ = compute_scores(
                        F_train, d["families"][train_mask], d["N_arr"][train_mask],
                        d["sigma_hist"][train_mask], d["pi_cg"][train_mask],
                    )
                    score = A + B + C
                    if score > best_A_train + best_B_train + (best_params[2] if best_params else 0):
                        best_A_train, best_B_train = A, B
                        best_params = (alpha, lam, eta, A, B, C)

        # Evaluate on test set with best params
        alpha, lam, eta = best_params[0], best_params[1], best_params[2]
        F_test = compute_F7_from_arrays(
            d["log_H"][test_mask], d["pi_geo"][test_mask],
            d["sigma_hist"][test_mask], d["xi_dim"][test_mask],
            d["R_arr"][test_mask], d["N_arr"][test_mask],
            alpha0=alpha, lam=lam, eta=eta,
        )
        A_test, B_test, C_test, _ = compute_scores(
            F_test, d["families"][test_mask], d["N_arr"][test_mask],
            d["sigma_hist"][test_mask], d["pi_cg"][test_mask],
        )

        train_A, train_B, train_C = best_params[3], best_params[4], best_params[5]
        all_train_scores.append((train_A, train_B, train_C))
        all_test_scores.append((A_test, B_test, C_test))

        print(f"  Fold {fold}: params=(α={alpha}, λ={lam}, η={eta:.1f})")
        print(f"    Train: A={train_A:.3f}, B={train_B:.3f}, C={train_C:.2f}")
        print(f"    Test:  A={A_test:.3f}, B={B_test:.3f}, C={C_test:.2f}")
        print(f"    Gap:   ΔA={A_test-train_A:+.3f}, ΔB={B_test-train_B:+.3f}")

    # Summary
    mean_train_A = np.mean([s[0] for s in all_train_scores])
    mean_test_A = np.mean([s[0] for s in all_test_scores])
    mean_train_B = np.mean([s[1] for s in all_train_scores])
    mean_test_B = np.mean([s[1] for s in all_test_scores])

    gap_A = mean_test_A - mean_train_A
    gap_B = mean_test_B - mean_train_B

    print(f"\n  Summary:")
    print(f"    Mean train: A={mean_train_A:.3f}, B={mean_train_B:.3f}")
    print(f"    Mean test:  A={mean_test_A:.3f}, B={mean_test_B:.3f}")
    print(f"    Mean gap:   ΔA={gap_A:+.3f}, ΔB={gap_B:+.3f}")

    verdict = "PASS" if (gap_A > -0.15 and gap_B > -0.15) else "FAIL"
    print(f"\n  VERDICT: {verdict}")
    if verdict == "PASS":
        print("  → Test performance close to train — no significant overfitting")
    else:
        print("  → WARNING: Large train-test gap suggests overfitting")
    return {"gap_A": gap_A, "gap_B": gap_B, "verdict": verdict}


# ══════════════════════════════════════════════════════════════════════════
# TEST 5: Sigmoid Wall Necessity — F7 without wall at various N
# ══════════════════════════════════════════════════════════════════════════
def test_sigmoid_necessity(d):
    """
    Run F7 with and without sigmoid wall, broken down by N.
    Check if predictions degrade gracefully or catastrophically.
    Also test: does the wall's importance decrease with N (finite-size effect)?
    """
    print("\n" + "=" * 70)
    print("TEST 5: SIGMOID WALL NECESSITY — With vs without, by N")
    print("=" * 70)

    unique_N = sorted(set(d["N_arr"]))

    print(f"\n  {'N':>4} | {'A_with':>7} {'B_with':>7} {'C_with':>7} | {'A_without':>9} {'B_without':>9} {'C_without':>9} | {'ΔA':>6} {'ΔB':>6}")
    print("  " + "-" * 80)

    deltas_A = []
    deltas_B = []

    for nv in unique_N:
        mask = d["N_arr"] == nv
        # With wall
        F_with = compute_F7_from_arrays(
            d["log_H"][mask], d["pi_geo"][mask], d["sigma_hist"][mask],
            d["xi_dim"][mask], d["R_arr"][mask], d["N_arr"][mask],
        )
        A_w, B_w, C_w, _ = compute_scores(
            F_with, d["families"][mask], d["N_arr"][mask],
            d["sigma_hist"][mask], d["pi_cg"][mask],
        )
        # Without wall
        F_without = d["log_H"][mask] + 0.0004 * d["pi_geo"][mask] - 5.0 * d["sigma_hist"][mask] + 1.0 * d["xi_dim"][mask]
        A_wo, B_wo, C_wo, _ = compute_scores(
            F_without, d["families"][mask], d["N_arr"][mask],
            d["sigma_hist"][mask], d["pi_cg"][mask],
        )
        dA = A_wo - A_w
        dB = B_wo - B_w
        deltas_A.append(dA)
        deltas_B.append(dB)
        print(f"  {int(nv):4d} | {A_w:7.3f} {B_w:7.3f} {C_w:7.2f} | {A_wo:9.3f} {B_wo:9.3f} {C_wo:9.2f} | {dA:+6.3f} {dB:+6.3f}")

    # Check if wall importance decreases with N (finite-size effect hypothesis)
    if len(unique_N) >= 3:
        corr_A = np.corrcoef(list(unique_N), deltas_A)[0, 1]
        corr_B = np.corrcoef(list(unique_N), deltas_B)[0, 1]
        print(f"\n  Correlation of ΔA with N: {corr_A:+.3f}")
        print(f"  Correlation of ΔB with N: {corr_B:+.3f}")
        if corr_A > 0.3 or corr_B > 0.3:
            print("  → Wall becomes less important at larger N (finite-size effect)")
        else:
            print("  → Wall importance does NOT decrease with N")

    verdict = "INFORMATIVE"
    print(f"\n  VERDICT: {verdict} (this test characterizes, not pass/fail)")
    return {"deltas_A": deltas_A, "deltas_B": deltas_B}


# ══════════════════════════════════════════════════════════════════════════
# TEST 6: New Seed Test — Entirely new posets
# ══════════════════════════════════════════════════════════════════════════
def test_new_seeds(n_reps=8, seed_offset=500000):
    """
    Generate entirely new posets with seeds never used in the original dataset.
    Compute features from scratch and evaluate F7.
    If predictions hold on unseen data, F7 generalizes.
    """
    print("\n" + "=" * 70)
    print("TEST 6: NEW SEED TEST — Entirely new posets")
    print("=" * 70)
    print(f"  Generating {n_reps} new posets per (family, N) with seed_offset={seed_offset}")

    from unified_functional import compute_log_H, compute_pi_geo, compute_sigma_hist, compute_xi_dim, compute_pi_cg

    N_values = [16, 20, 28, 36]
    all_families = []
    all_N = []
    all_logH = []
    all_pi_geo = []
    all_sigma = []
    all_xi = []
    all_pi_cg_arr = []
    all_R = []

    for n_val in N_values:
        for fam_name, gen_func in GENS.items():
            for rep in range(n_reps):
                s = seed_offset + rep * 1000 + n_val * 100
                p = gen_func(n_val, seed=s)

                all_families.append(fam_name)
                all_N.append(float(n_val))
                all_logH.append(compute_log_H(p))
                all_pi_geo.append(compute_pi_geo(p))
                all_sigma.append(compute_sigma_hist(p))
                xi, _ = compute_xi_dim(p)
                all_xi.append(xi)
                all_pi_cg_arr.append(compute_pi_cg(p, keep_ratio=0.7, n_cg_samples=3))

                # Compute R
                c = p.closure.astype(np.int32)
                ks = c @ c
                mask = p.closure
                C0 = int(np.sum(mask & (ks == 0)))
                total = int(np.sum(mask))
                all_R.append(1.0 - C0 / total if total > 0 else 0.0)

    families = np.array(all_families)
    N_arr = np.array(all_N)
    log_H = np.array(all_logH)
    pi_geo = np.array(all_pi_geo)
    sigma_hist = np.array(all_sigma)
    xi_dim = np.array(all_xi)
    pi_cg = np.array(all_pi_cg_arr)
    R_arr = np.array(all_R)

    print(f"  Generated {len(families)} new posets")

    F_new = compute_F7_from_arrays(log_H, pi_geo, sigma_hist, xi_dim, R_arr, N_arr)
    A, B, C, D = compute_scores(F_new, families, N_arr, sigma_hist, pi_cg)

    print(f"\n  New-seed F7 scores: A={A:.3f}, B={B:.3f}, C={C:.2f}, D={D:.3f}")
    print(f"  Grade: {grade_str(A, B, C)}")

    # Compare with original dataset scores
    print(f"\n  (Compare with original dataset to check generalization)")

    verdict = "PASS" if (A >= 0.55 and B >= 0.85) else "FAIL"
    print(f"\n  VERDICT: {verdict}")
    if verdict == "PASS":
        print("  → F7 generalizes to entirely unseen posets")
    else:
        print("  → WARNING: F7 fails on new data — possible overfitting to original seeds")
    return {"A": A, "B": B, "C": C, "verdict": verdict}


# ══════════════════════════════════════════════════════════════════════════
# TEST 7: Fake Family Test — Synthetic families
# ══════════════════════════════════════════════════════════════════════════
def test_fake_families(n_reps=8, seed_base=888888):
    """
    Create synthetic "families" with known properties:
    - TransPerc: transitive percolation (should NOT be Lorentzian)
    - IntOrder: interval orders (should NOT be Lorentzian)
    - AbsLayer: absolute layered (should NOT be Lorentzian)
    
    F7 should correctly rank these as higher (worse) than Lor2D.
    If F7 ranks them as better than Lor2D, it's not capturing real physics.
    """
    print("\n" + "=" * 70)
    print("TEST 7: FAKE FAMILY TEST — Synthetic non-Lorentzian families")
    print("=" * 70)

    from unified_functional import compute_log_H, compute_pi_geo, compute_sigma_hist, compute_xi_dim, compute_pi_cg

    fake_gens = {
        "TransPerc": lambda n, seed: generate_transitive_percolation(n, p=0.08, seed=seed),
        "IntOrder": lambda n, seed: generate_interval_order(n, mean_length=0.3, seed=seed),
        "AbsLayer": lambda n, seed: generate_absolute_layered(n, n_layers=4, seed=seed),
    }

    N_values = [16, 20, 28]

    # First compute Lor2D baseline
    lor2d_F = {}
    for n_val in N_values:
        fs = []
        for rep in range(n_reps):
            s = seed_base + rep * 1000 + n_val * 100
            p = generate_lorentzian_like_2d(n_val, seed=s)
            log_H = compute_log_H(p)
            pi_geo = compute_pi_geo(p)
            sigma = compute_sigma_hist(p)
            xi, _ = compute_xi_dim(p)
            c = p.closure.astype(np.int32)
            ks = c @ c; mask = p.closure
            C0 = int(np.sum(mask & (ks == 0))); total = int(np.sum(mask))
            R = 1.0 - C0 / total if total > 0 else 0.0
            alpha_N = 26.0 * (20.0 / n_val) ** 0.5
            wall = alpha_N * sigmoid((R - 0.25) / 0.015)
            f7 = log_H + 0.0004 * pi_geo - 5.0 * sigma + 1.0 * xi + wall
            fs.append(f7)
        lor2d_F[n_val] = np.mean(fs)

    print(f"\n  Lor2D baseline F7: {lor2d_F}")

    print(f"\n  {'Family':<12} {'N':>4} {'mean_F7':>10} {'Lor2D_F7':>10} {'F7>Lor2D?':>10}")
    print("  " + "-" * 50)

    all_correct = 0
    all_total = 0

    for fake_name, gen_func in fake_gens.items():
        for n_val in N_values:
            fs = []
            for rep in range(n_reps):
                s = seed_base + 100000 + rep * 1000 + n_val * 100
                try:
                    p = gen_func(n_val, seed=s)
                    log_H = compute_log_H(p)
                    pi_geo = compute_pi_geo(p)
                    sigma = compute_sigma_hist(p)
                    xi, _ = compute_xi_dim(p)
                    c = p.closure.astype(np.int32)
                    ks = c @ c; mask = p.closure
                    C0 = int(np.sum(mask & (ks == 0))); total = int(np.sum(mask))
                    R = 1.0 - C0 / total if total > 0 else 0.0
                    alpha_N = 26.0 * (20.0 / n_val) ** 0.5
                    wall = alpha_N * sigmoid((R - 0.25) / 0.015)
                    f7 = log_H + 0.0004 * pi_geo - 5.0 * sigma + 1.0 * xi + wall
                    fs.append(f7)
                except Exception as e:
                    print(f"  WARNING: {fake_name} N={n_val} rep={rep} failed: {e}")

            if fs:
                mean_f7 = np.mean(fs)
                correct = mean_f7 > lor2d_F[n_val]
                all_total += 1
                if correct:
                    all_correct += 1
                print(f"  {fake_name:<12} {n_val:4d} {mean_f7:10.2f} {lor2d_F[n_val]:10.2f} {'YES ✓' if correct else 'NO ✗':>10}")

    accuracy = all_correct / all_total if all_total > 0 else 0
    print(f"\n  Fake family rejection rate: {all_correct}/{all_total} = {accuracy:.1%}")

    verdict = "PASS" if accuracy >= 0.75 else "FAIL"
    print(f"\n  VERDICT: {verdict}")
    if verdict == "PASS":
        print("  → F7 correctly identifies non-Lorentzian structures as higher-action")
    else:
        print("  → WARNING: F7 fails to distinguish fake families from Lorentzian")
    return {"accuracy": accuracy, "verdict": verdict}


# ══════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════
def main():
    print("=" * 70)
    print("ANTI-ENGINEERING-ILLUSION TEST BATTERY")
    print("Detecting overfitting / phenomenological patching in F7")
    print("=" * 70)

    t0 = time.time()

    print("\nLoading data and computing R...")
    d = load_data()
    print(f"  Loaded {len(d['data'])} samples in {time.time()-t0:.1f}s")

    results = {}

    # Test 1: Permutation
    results["permutation"] = test_permutation(d)

    # Test 2: Random functional
    results["random_functional"] = test_random_functional(d)

    # Test 3: Ablation
    results["ablation"] = test_ablation(d)

    # Test 4: Cross-validation
    results["cross_validation"] = test_cross_validation(d)

    # Test 5: Sigmoid necessity
    results["sigmoid_necessity"] = test_sigmoid_necessity(d)

    # Test 6: New seeds
    print("\n  (Test 6 generates new posets — may take a minute...)")
    results["new_seeds"] = test_new_seeds()

    # Test 7: Fake families
    print("\n  (Test 7 generates fake families — may take a minute...)")
    results["fake_families"] = test_fake_families()

    # ── Final Summary ──
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    verdicts = {k: v["verdict"] for k, v in results.items() if "verdict" in v}
    for test_name, verdict in verdicts.items():
        symbol = "✓" if verdict == "PASS" else "✗" if verdict == "FAIL" else "?"
        print(f"  {symbol} {test_name}: {verdict}")

    n_pass = sum(1 for v in verdicts.values() if v == "PASS")
    n_fail = sum(1 for v in verdicts.values() if v == "FAIL")
    n_other = sum(1 for v in verdicts.values() if v not in ("PASS", "FAIL"))

    print(f"\n  Total: {n_pass} PASS, {n_fail} FAIL, {n_other} INFORMATIVE")
    print(f"  Elapsed: {time.time()-t0:.1f}s")

    if n_fail == 0:
        print("\n  CONCLUSION: No evidence of engineering illusion detected.")
        print("  F7's predictive power appears genuine and robust.")
    elif n_fail <= 2:
        print("\n  CONCLUSION: Some concerns detected — review failing tests.")
        print("  F7 may have partial overfitting in specific aspects.")
    else:
        print("\n  CONCLUSION: Significant evidence of engineering illusion.")
        print("  F7's predictions may be artifacts of data fitting.")


if __name__ == "__main__":
    main()
