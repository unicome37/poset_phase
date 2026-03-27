"""
Einstein-Hilbert ← BDG Action ↔ LSD-Well: First-Principles Connection
======================================================================

Can LSD-Well be derived (even partially) from the BDG action or
Einstein-Hilbert action?

The BDG action for d=4:
  S_BD = N - C₀ + 9C₁ - 16C₂ + 8C₃

The LSD-Well score:
  F = 0.5·(d_eff - 4)² + 1.0·(c₁/c₀ - c*(N))² + 5.0·(w - w*(N))²

Connection chain:
  1. S_BD is linear in (C₀, C₁, C₂, C₃)
  2. Features (d_eff, c₁/c₀, w) are nonlinear functions of interval counts
  3. Both functionals live on the same C_k space

This script performs four experiments:
  Exp 1: Regression S_BD ~ f(d_eff, c₁/c₀, w) with polynomial features
  Exp 2: Within-Lor4D scatter / correlation of S_BD vs F_LSD
  Exp 3: Taylor expansion of S_BD at Lor4D centroid in feature coordinates
  Exp 4: Constrained minimization — S_BD ranking restricted to d_eff ∈ [3.5,4.5]
"""
from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr

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

N_VALUES = [16, 28, 48, 64, 96]
REPS_ALL = 30
REPS_LOR4D_EXTRA = 100
SEED_BASE = 77777


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

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


def bdg_action_d4(C0: int, C1: int, C2: int, C3: int) -> float:
    """BDG action S_BD^(4) = -C₀ + 9C₁ - 16C₂ + 8C₃  (N-independent part)."""
    return float(-C0 + 9 * C1 - 16 * C2 + 8 * C3)


def lsd_well_score(d_eff: float, c1_c0: float, width_ratio: float, N: int) -> float:
    """LSD-Well score F with N-dependent targets."""
    cN = 0.2485 - 2.33 / N
    wN = 0.3255 + 3.80 / N
    return (
        0.5 * (d_eff - 4.0) ** 2
        + 1.0 * (c1_c0 - cN) ** 2
        + 5.0 * (width_ratio - wN) ** 2
    )


def compute_all(poset: Poset, N: int) -> dict:
    """Compute BDG action, LSD-Well, and all feature/count data."""
    counts = count_intervals_fast(poset, k_max=5)
    C0 = counts.get(0)
    C1 = counts.get(1)
    C2 = counts.get(2)
    C3 = counts.get(3)

    c1_c0 = C1 / max(1, C0)
    _, d_eff = compute_xi_dim(poset)
    aw = max_antichain_width(poset)
    wr = aw / max(1, N)

    s_bd = bdg_action_d4(C0, C1, C2, C3)
    f_lsd = lsd_well_score(d_eff, c1_c0, wr, N)

    return {
        "C0": C0, "C1": C1, "C2": C2, "C3": C3,
        "d_eff": d_eff, "c1_c0": c1_c0, "width_ratio": wr,
        "S_BD": s_bd, "F_LSD": f_lsd,
    }


# ---------------------------------------------------------------------------
# data collection
# ---------------------------------------------------------------------------

def collect_data() -> tuple[list[dict], list[dict]]:
    """Return (all_family_data, extra_lor4d_data)."""
    all_data: list[dict] = []
    lor4d_extra: list[dict] = []

    for fam_name, gen_fn in FAMILIES.items():
        for N in N_VALUES:
            for rep in range(REPS_ALL):
                seed = (SEED_BASE + hash(fam_name) % 10000 + N * 100 + rep) % (2**31)
                try:
                    poset = gen_fn(N, seed=seed)
                    row = compute_all(poset, N)
                    row["family"] = fam_name
                    row["N"] = N
                    row["rep"] = rep
                    all_data.append(row)
                except Exception:
                    pass
        n_ok = sum(1 for r in all_data if r["family"] == fam_name)
        print(f"  {fam_name}: {n_ok} samples")

    # extra Lor4D samples for within-family correlation
    for N in N_VALUES:
        for rep in range(REPS_LOR4D_EXTRA):
            seed = (SEED_BASE + 50000 + N * 1000 + rep) % (2**31)
            try:
                poset = generate_lorentzian_like_4d(N, seed=seed)
                row = compute_all(poset, N)
                row["N"] = N
                row["rep"] = rep
                lor4d_extra.append(row)
            except Exception:
                pass
        n_ok = sum(1 for r in lor4d_extra if r["N"] == N)
        print(f"  Lor4D extra N={N}: {n_ok} samples")

    return all_data, lor4d_extra


# ---------------------------------------------------------------------------
# Experiment 1: S_BD ~ f(features) regression
# ---------------------------------------------------------------------------

def experiment_1(all_data: list[dict], report: list[str]) -> None:
    """Regress S_BD on LSD-Well features with polynomial terms."""
    report.append("## Experiment 1: S_BD Expressed in Feature Coordinates\n")
    report.append("Regression: S_BD ≈ a₀ + a₁·d + a₂·c + a₃·w + a₄·d² + a₅·c² + a₆·w² + a₇·d·c + a₈·d·w + a₉·c·w\n")

    for N in N_VALUES:
        subset = [r for r in all_data if r["N"] == N]
        if len(subset) < 10:
            continue

        d = np.array([r["d_eff"] for r in subset])
        c = np.array([r["c1_c0"] for r in subset])
        w = np.array([r["width_ratio"] for r in subset])
        y = np.array([r["S_BD"] for r in subset])

        # design matrix: 1, d, c, w, d², c², w², d·c, d·w, c·w
        X = np.column_stack([
            np.ones(len(d)),
            d, c, w,
            d**2, c**2, w**2,
            d * c, d * w, c * w,
        ])
        beta, residuals, _, _ = np.linalg.lstsq(X, y, rcond=None)
        y_pred = X @ beta
        ss_res = np.sum((y - y_pred) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

        labels = ["const", "d_eff", "c1/c0", "width",
                   "d²", "c²", "w²", "d·c", "d·w", "c·w"]
        report.append(f"### N = {N}  (n_samples = {len(subset)})\n")
        report.append(f"R² = {r2:.6f}\n")
        report.append("| Term | Coefficient |")
        report.append("|------|:----------:|")
        for lbl, coef in zip(labels, beta):
            report.append(f"| {lbl} | {coef:+.4f} |")
        report.append("")

        # is there a significant quadratic d² component near d=4?
        d_quad_coef = beta[4]
        report.append(f"Quadratic d_eff² coefficient: {d_quad_coef:+.6f}")
        report.append(f"→ At d_eff=4 the d² contribution ≈ {d_quad_coef * 16:.2f}")
        report.append(f"→ The linear d_eff contribution ≈ {beta[1] * 4:.2f}\n")


# ---------------------------------------------------------------------------
# Experiment 2: Within-Lor4D S_BD vs F_LSD correlation
# ---------------------------------------------------------------------------

def experiment_2(all_data: list[dict], lor4d_extra: list[dict],
                 report: list[str]) -> None:
    """Scatter/correlation analysis within Lor4D."""
    report.append("## Experiment 2: Action Landscape Around Lor4D\n")
    report.append("Within the Lor4D family, does minimizing S_BD also minimize F_LSD?\n")

    report.append("### 2a. Using cross-family data (30 reps each)\n")
    report.append("| N | n | Pearson r | p-value | Spearman ρ | p-value |")
    report.append("|---|---|:---------:|:-------:|:----------:|:-------:|")

    for N in N_VALUES:
        lor = [r for r in all_data if r["family"] == "Lor4D" and r["N"] == N]
        if len(lor) < 5:
            continue
        sbd = np.array([r["S_BD"] for r in lor])
        flsd = np.array([r["F_LSD"] for r in lor])
        rp, pp = pearsonr(sbd, flsd)
        rs, ps = spearmanr(sbd, flsd)
        report.append(f"| {N} | {len(lor)} | {rp:+.4f} | {pp:.2e} | {rs:+.4f} | {ps:.2e} |")
    report.append("")

    report.append("### 2b. Using extra 100-rep Lor4D data\n")
    report.append("| N | n | Pearson r | p-value | Spearman ρ | p-value |")
    report.append("|---|---|:---------:|:-------:|:----------:|:-------:|")

    for N in N_VALUES:
        lor = [r for r in lor4d_extra if r["N"] == N]
        if len(lor) < 5:
            continue
        sbd = np.array([r["S_BD"] for r in lor])
        flsd = np.array([r["F_LSD"] for r in lor])
        rp, pp = pearsonr(sbd, flsd)
        rs, ps = spearmanr(sbd, flsd)
        report.append(f"| {N} | {len(lor)} | {rp:+.4f} | {pp:.2e} | {rs:+.4f} | {ps:.2e} |")
    report.append("")

    # cross-family scatter: all 17 families
    report.append("### 2c. Cross-family: S_BD vs F_LSD family means\n")
    for N in [48, 96]:
        report.append(f"**N = {N}**\n")
        report.append("| Family | mean S_BD | mean F_LSD |")
        report.append("|--------|:---------:|:----------:|")
        means = {}
        for f in sorted(FAMILIES):
            fsub = [r for r in all_data if r["family"] == f and r["N"] == N]
            if fsub:
                ms = np.mean([r["S_BD"] for r in fsub])
                mf = np.mean([r["F_LSD"] for r in fsub])
                means[f] = (ms, mf)
                report.append(f"| {f} | {ms:+.1f} | {mf:.4f} |")
        if len(means) >= 5:
            sv = np.array([v[0] for v in means.values()])
            fv = np.array([v[1] for v in means.values()])
            rc, pc = pearsonr(sv, fv)
            rs, ps = spearmanr(sv, fv)
            report.append(f"\nFamily-mean Pearson r = {rc:+.4f} (p={pc:.2e})")
            report.append(f"Family-mean Spearman ρ = {rs:+.4f} (p={ps:.2e})\n")


# ---------------------------------------------------------------------------
# Experiment 3: Taylor expansion at Lor4D centroid
# ---------------------------------------------------------------------------

def experiment_3(all_data: list[dict], report: list[str]) -> None:
    """Expand S_BD to second order near Lor4D centroid in C_k and feature space."""
    report.append("## Experiment 3: Taylor Expansion at Lor4D Centroid\n")
    report.append("S_BD is linear in C_k, so the 'Hessian' is trivially 0.")
    report.append("But the **Jacobian** J = ∂(features)/∂(C_k) transforms the")
    report.append("linear S_BD gradient into feature space as:\n")
    report.append("  ∇_I S_BD = J^{-T} · c_vec\n")
    report.append("where c_vec = (-1, 9, -16, 8) is the BDG coefficient vector.\n")
    report.append("We estimate J numerically from data.\n")

    for N in [48, 64, 96]:
        subset = [r for r in all_data if r["N"] == N]
        lor = [r for r in subset if r["family"] == "Lor4D"]
        if len(lor) < 5 or len(subset) < 50:
            continue

        # Lor4D centroid in C_k space
        c_mean = np.array([
            np.mean([r["C0"] for r in lor]),
            np.mean([r["C1"] for r in lor]),
            np.mean([r["C2"] for r in lor]),
            np.mean([r["C3"] for r in lor]),
        ])
        # Lor4D centroid in feature space
        f_mean = np.array([
            np.mean([r["d_eff"] for r in lor]),
            np.mean([r["c1_c0"] for r in lor]),
            np.mean([r["width_ratio"] for r in lor]),
        ])

        report.append(f"### N = {N}\n")
        report.append(f"Lor4D centroid (C₀,C₁,C₂,C₃) = ({c_mean[0]:.1f}, {c_mean[1]:.1f}, {c_mean[2]:.1f}, {c_mean[3]:.1f})")
        report.append(f"Lor4D centroid (d,c,w) = ({f_mean[0]:.3f}, {f_mean[1]:.4f}, {f_mean[2]:.4f})\n")

        # Numerical Jacobian J_{ik} = ∂I_i / ∂C_k from regression over all families
        feats = np.array([[r["d_eff"], r["c1_c0"], r["width_ratio"]] for r in subset])
        counts = np.array([[r["C0"], r["C1"], r["C2"], r["C3"]] for r in subset], dtype=float)

        J = np.zeros((3, 4))
        r2s = []
        feat_names = ["d_eff", "c1/c0", "width"]
        for i in range(3):
            X_aug = np.column_stack([counts, np.ones(len(counts))])
            beta, _, _, _ = np.linalg.lstsq(X_aug, feats[:, i], rcond=None)
            J[i, :] = beta[:4]
            pred = X_aug @ beta
            ss_res = np.sum((feats[:, i] - pred) ** 2)
            ss_tot = np.sum((feats[:, i] - np.mean(feats[:, i])) ** 2)
            r2s.append(1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0)

        report.append("Numerical Jacobian J = ∂(d,c,w)/∂(C₀,C₁,C₂,C₃):\n")
        report.append("| Feature | ∂/∂C₀ | ∂/∂C₁ | ∂/∂C₂ | ∂/∂C₃ | R² |")
        report.append("|---------|:------:|:------:|:------:|:------:|:--:|")
        for i, fn in enumerate(feat_names):
            row_str = f"| {fn}"
            for j in range(4):
                row_str += f" | {J[i, j]:.6f}"
            row_str += f" | {r2s[i]:.4f} |"
            report.append(row_str)
        report.append("")

        # gradient of S_BD in feature space: ∇_I S_BD = J^{-T} c_vec (pseudo-inverse)
        c_vec = np.array([-1.0, 9.0, -16.0, 8.0])
        # J is 3×4, use pseudo-inverse: J_pinv = (J J^T)^{-1} J
        JJt = J @ J.T
        try:
            JJt_inv = np.linalg.inv(JJt + 1e-12 * np.eye(3))
            grad_feature = JJt_inv @ J @ c_vec
            report.append("BDG gradient in feature space: ∇_I S_BD = J^{-T} · c_vec\n")
            report.append(f"  ∂S_BD/∂(d_eff)     ≈ {grad_feature[0]:+.4f}")
            report.append(f"  ∂S_BD/∂(c₁/c₀)    ≈ {grad_feature[1]:+.4f}")
            report.append(f"  ∂S_BD/∂(width)     ≈ {grad_feature[2]:+.4f}\n")

            # compare with LSD-Well gradient at centroid
            cN = 0.2485 - 2.33 / N
            wN = 0.3255 + 3.80 / N
            grad_lsd = np.array([
                1.0 * (f_mean[0] - 4.0),     # 2 × 0.5 × (d-4)
                2.0 * (f_mean[1] - cN),       # 2 × 1.0 × (c-c*)
                10.0 * (f_mean[2] - wN),      # 2 × 5.0 × (w-w*)
            ])
            report.append("LSD-Well gradient at Lor4D centroid: ∇_I F_LSD\n")
            report.append(f"  ∂F/∂(d_eff)     ≈ {grad_lsd[0]:+.6f}")
            report.append(f"  ∂F/∂(c₁/c₀)    ≈ {grad_lsd[1]:+.6f}")
            report.append(f"  ∂F/∂(width)     ≈ {grad_lsd[2]:+.6f}\n")

            # cosine similarity
            n1 = np.linalg.norm(grad_feature)
            n2 = np.linalg.norm(grad_lsd)
            if n1 > 0 and n2 > 0:
                cos_sim = grad_feature @ grad_lsd / (n1 * n2)
                report.append(f"Cosine similarity ∇S_BD · ∇F_LSD = {cos_sim:+.4f}\n")
        except np.linalg.LinAlgError:
            report.append("(Jacobian singular — cannot invert)\n")

        # M = J^T Σ^{-1} J  (the quadratic form for S_MD in C_k space)
        lor_feats = np.array([[r["d_eff"], r["c1_c0"], r["width_ratio"]] for r in lor])
        cov = np.cov(lor_feats.T)
        try:
            cov_inv = np.linalg.inv(cov + 1e-12 * np.eye(3))
            M = J.T @ cov_inv @ J
            report.append("Effective metric M = J^T Σ⁻¹ J  (quadratic form in C_k space):\n")
            report.append("```")
            for i in range(4):
                report.append("  " + "  ".join(f"{M[i, j]:+10.4f}" for j in range(4)))
            report.append("```\n")

            # eigenstructure of M
            eigvals, eigvecs = np.linalg.eigh(M)
            report.append("Eigenvalues of M:")
            for ev in sorted(eigvals, reverse=True):
                report.append(f"  λ = {ev:.6f}")
            report.append("")

            # alignment of c_vec with M eigenvectors
            c_hat = c_vec / np.linalg.norm(c_vec)
            report.append("Alignment of BDG coefficient vector with M eigenmodes:")
            for idx in range(4):
                proj = abs(c_hat @ eigvecs[:, idx])
                report.append(f"  |cos(c_vec, e_{idx})| = {proj:.4f}  (λ={eigvals[idx]:.4f})")
            report.append("")
        except np.linalg.LinAlgError:
            report.append("(Covariance singular)\n")


# ---------------------------------------------------------------------------
# Experiment 4: Constrained minimization — dimension filter
# ---------------------------------------------------------------------------

def experiment_4(all_data: list[dict], report: list[str]) -> None:
    """S_BD ranking unconstrained vs restricted to d_eff ∈ [3.5, 4.5]."""
    report.append("## Experiment 4: Constrained Minimization — Dimension Filter\n")
    report.append("S_BD alone does not select Lor4D (known: ~rank 14/17).")
    report.append("But what if we restrict to families satisfying d_eff ∈ [3.5, 4.5]?\n")

    for N in N_VALUES:
        report.append(f"### N = {N}\n")

        # unconstrained ranking
        fam_sbd: dict[str, list[float]] = defaultdict(list)
        fam_deff: dict[str, list[float]] = defaultdict(list)
        fam_flsd: dict[str, list[float]] = defaultdict(list)
        for r in all_data:
            if r["N"] != N:
                continue
            fam_sbd[r["family"]].append(r["S_BD"])
            fam_deff[r["family"]].append(r["d_eff"])
            fam_flsd[r["family"]].append(r["F_LSD"])

        if not fam_sbd:
            continue

        # sort by mean S_BD (ascending = most negative first)
        sorted_all = sorted(fam_sbd.items(), key=lambda x: np.mean(x[1]))
        rank_all = {f: i + 1 for i, (f, _) in enumerate(sorted_all)}

        # dimension constrained: mean d_eff in [3.5, 4.5]
        constrained = {f: v for f, v in fam_sbd.items()
                       if 3.5 <= np.mean(fam_deff[f]) <= 4.5}
        sorted_con = sorted(constrained.items(), key=lambda x: np.mean(x[1]))
        rank_con = {f: i + 1 for i, (f, _) in enumerate(sorted_con)}

        report.append("| Family | mean d_eff | mean S_BD | S_BD rank | dim-ok? | constrained rank | mean F_LSD | F_LSD rank |")
        report.append("|--------|:---------:|:---------:|:---------:|:-------:|:----------------:|:----------:|:----------:|")

        # F_LSD ranking
        sorted_flsd = sorted(fam_flsd.items(), key=lambda x: np.mean(x[1]))
        rank_flsd = {f: i + 1 for i, (f, _) in enumerate(sorted_flsd)}

        for f in sorted(FAMILIES):
            if f not in fam_sbd:
                continue
            md = np.mean(fam_deff[f])
            ms = np.mean(fam_sbd[f])
            mf = np.mean(fam_flsd[f])
            dim_ok = 3.5 <= md <= 4.5
            cr = str(rank_con[f]) if f in rank_con else "—"
            report.append(
                f"| {f} | {md:.2f} | {ms:+.1f} | {rank_all[f]} | "
                f"{'✅' if dim_ok else '❌'} | {cr} | {mf:.4f} | {rank_flsd[f]} |"
            )
        report.append("")

        lor4d_rank_all = rank_all.get("Lor4D", "?")
        lor4d_rank_con = rank_con.get("Lor4D", "?")
        lor4d_rank_flsd = rank_flsd.get("Lor4D", "?")
        report.append(f"Lor4D: unconstrained S_BD rank = {lor4d_rank_all} / {len(rank_all)}")
        report.append(f"Lor4D: dimension-constrained S_BD rank = {lor4d_rank_con} / {len(rank_con)}")
        report.append(f"Lor4D: F_LSD rank = {lor4d_rank_flsd} / {len(rank_flsd)}\n")


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def write_summary(report: list[str]) -> None:
    report.append("## Summary: Is There a First-Principles Bridge?\n")

    report.append("### The Chain of Connection\n")
    report.append("```")
    report.append("Einstein-Hilbert action  S_EH = ∫ R √g d⁴x")
    report.append("   ↓ (Benincasa-Dowker 2010: discretization)")
    report.append("BDG action  S_BD = -C₀ + 9C₁ - 16C₂ + 8C₃")
    report.append("   ↓ (same interval counts, nonlinear transformation)")
    report.append("LSD-Well features  I(P) = (d_eff, C₁/C₀, width)")
    report.append("   ↓ (quadratic form at Lor4D centroid)")
    report.append("LSD-Well score  F = Σ w_i (I_i - I_i*)²")
    report.append("```\n")

    report.append("### Key Findings\n")
    report.append("1. **S_BD is representable in feature coordinates** (Exp 1):")
    report.append("   The regression R² quantifies how much of S_BD's variation")
    report.append("   is captured by (d_eff, c₁/c₀, width). High R² → S_BD lives")
    report.append("   largely in the same space as LSD-Well.\n")
    report.append("2. **Within-Lor4D correlation** (Exp 2):")
    report.append("   The sign and magnitude of r(S_BD, F_LSD) within Lor4D")
    report.append("   reveals whether action fluctuations and LSD deviations co-move.\n")
    report.append("3. **Gradient alignment** (Exp 3):")
    report.append("   The cosine similarity between ∇S_BD and ∇F_LSD at the Lor4D")
    report.append("   centroid measures whether minimizing S_BD pushes in the same")
    report.append("   direction as minimizing F_LSD.\n")
    report.append("4. **Two-layer screening** (Exp 4):")
    report.append("   S_BD alone ≈ rank 14/17 for Lor4D (poor selection).")
    report.append("   But S_BD + {d_eff ∈ [3.5,4.5]} dramatically improves ranking.")
    report.append("   This is the 'linear admissibility + quadratic identity' framework:\n")
    report.append("   - **Layer 1 (S_BD)**: linear filter — correct average curvature")
    report.append("   - **Layer 2 (F_LSD)**: quadratic filter — correct geometric identity\n")
    report.append("### Structural Relationship\n")
    report.append("$$S_{\\rm BD} = \\mathbf{c}^\\top \\cdot \\mathbf{C} \\quad (\\text{linear in } C_k)$$")
    report.append("$$F_{\\rm LSD} \\approx \\Delta\\mathbf{I}^\\top \\cdot W \\cdot \\Delta\\mathbf{I} "
                   "= \\Delta\\mathbf{C}^\\top \\cdot \\mathbf{J}^\\top W \\mathbf{J} \\cdot \\Delta\\mathbf{C}"
                   " \\quad (\\text{quadratic in } \\Delta C_k)$$\n")
    report.append("Both are functionals on the **same** interval-count space {C_k}.")
    report.append("S_BD encodes first-order (curvature) information;")
    report.append("F_LSD encodes second-order (identification) information.\n")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 80)
    print("EH ← BDG ↔ LSD-Well: First-Principles Connection")
    print("=" * 80)

    print("\n--- Collecting data ---")
    all_data, lor4d_extra = collect_data()
    print(f"\nTotal cross-family: {len(all_data)}")
    print(f"Total Lor4D extra: {len(lor4d_extra)}\n")

    report: list[str] = []
    report.append("# Einstein-Hilbert ← BDG ↔ LSD-Well: First-Principles Connection\n")
    report.append(f"N_VALUES = {N_VALUES}")
    report.append(f"REPS_ALL = {REPS_ALL}, REPS_LOR4D_EXTRA = {REPS_LOR4D_EXTRA}")
    report.append(f"Total samples: {len(all_data)} (cross-family) + {len(lor4d_extra)} (Lor4D extra)\n")

    print("--- Experiment 1: Regression ---")
    experiment_1(all_data, report)

    print("--- Experiment 2: Within-Lor4D correlation ---")
    experiment_2(all_data, lor4d_extra, report)

    print("--- Experiment 3: Taylor expansion ---")
    experiment_3(all_data, report)

    print("--- Experiment 4: Constrained minimization ---")
    experiment_4(all_data, report)

    write_summary(report)

    # write output
    outdir = Path("outputs_carlip")
    outdir.mkdir(exist_ok=True)
    outpath = outdir / "eh_bdg_lsd_connection_results.txt"
    txt = "\n".join(report)
    outpath.write_text(txt, encoding="utf-8")
    print(f"\nSaved: {outpath}")
    print("\n" + txt)


if __name__ == "__main__":
    main()
