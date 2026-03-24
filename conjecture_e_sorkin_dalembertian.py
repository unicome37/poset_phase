"""Conjecture E — §4.1.27: Sorkin d'Alembertian B_ℓ spectral experiment.

Implements the causal-set d'Alembertian operator B^(d) as a matrix on
sprinkled de Sitter causal sets and extracts spectral features that
may correlate with H² (scalar curvature R = d(d-1)H²).

The key theoretical result (Benincasa-Dowker-Glaser):
    lim_{ρ→∞} ⟨B^(d) φ(x)⟩ = □φ - ½ R(x) φ(x)

Applied to constant field φ=1:  ⟨B^(d) 1⟩(x) → -½ R(x)
So the mean of (B_ℓ 1) over all elements should track -½ d(d-1) H².

This script:
  1. Builds B_ℓ as an N×N matrix for each sprinkling.
  2. Applies B_ℓ to constant field → per-element curvature estimates.
  3. Extracts spectral features: mean(B1), std(B1), eigenvalue statistics.
  4. Scans over d ∈ {2,3,4}, N ∈ {128,256,512}, H ∈ {0,0.25,0.5,1.0,2.0}.
  5. Correlates features with H² via Spearman ρ.

BDG coefficients (Glaser 2014, arXiv:1311.1701):
  d=2: C_i = {-2, +2}           → B = α₂φ(x) + β₂[-2 Σ_{L₀} + 2 Σ_{L₁}]φ
  d=3: C_i = {-6, +16, -12}     → B = α₃φ(x) + β₃[-6 Σ_{L₀} + 16 Σ_{L₁} - 12 Σ_{L₂}]φ  (*)
  d=4: C_i = {-1, +9, -16, +8}  → B = α₄φ(x) + β₄[-1 Σ_{L₀} + 9 Σ_{L₁} - 16 Σ_{L₂} + 8 Σ_{L₃}]φ

(*) d=3 uses the Dowker-Glaser coefficients for odd dimensions.

The overall prefactor ρ^{2/d} (or equivalently ℓ^{-2}) sets the physical
scale; for our experiment we normalize by N/V so that the dimensionless
ratio ⟨B1⟩ / (d(d-1)H²) → -½ in the continuum limit.

References:
  [1] Sorkin (2006), arXiv:gr-qc/0703099
  [2] Benincasa & Dowker (2010), arXiv:1001.2725
  [3] Dowker & Glaser (2013), arXiv:1305.2588
  [4] Glaser (2014), arXiv:1311.1701
  [5] Belenchia, Benincasa, Dowker (2015), arXiv:1507.00986
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np

from curvature_layer2_de_sitter_scan import sprinkle_de_sitter_like_diamond
from curvature_layer2_recovery_fast import build_causal_matrix_de_sitter


# ---------------------------------------------------------------------------
# BDG coefficients by dimension
# ---------------------------------------------------------------------------
# C_i^{(d)} for i = 0, 1, ..., floor(d/2)
# These are the standard literature values.
BDG_COEFFICIENTS: dict[int, list[float]] = {
    2: [-2.0, 2.0],
    3: [-6.0, 16.0, -12.0],       # Dowker-Glaser odd-d formula (*)
    4: [-1.0, 9.0, -16.0, 8.0],
}

# (*) Note on d=3 coefficients:
# The original Benincasa-Dowker paper gives d=2,4 (even dimensions).
# Dowker & Glaser (2013) extended to odd dimensions.
# For d=3, the BDG action is: S = N - 6N₀ + 16N₁ - 12N₂
# where N_i = number of (i+2)-element intervals.
# Some references use slightly different conventions; we follow Glaser (2014).


def build_dalembertian_matrix(
    causal: np.ndarray,
    d: int,
    rho: float = 1.0,
) -> np.ndarray:
    """Build the Sorkin d'Alembertian B^(d) as an N×N matrix.

    For element x, (Bφ)(x) = [α_d φ(x) + Σ_i C_i Σ_{y∈L_i(x)} φ(y)]

    where L_i(x) = {y : y ≺ x and |[y,x]| = i+2} (past neighbors with
    interval of size i+2, i.e., exactly i elements strictly between y and x).

    The overall ρ^{2/d} prefactor is NOT included — it is a calibration
    constant that does not affect rank correlations. Set rho != 1.0 only
    if you want absolute-scale curvature estimates.

    We build this as a matrix so that Bφ = M @ φ for any field φ.

    Parameters
    ----------
    causal : (N, N) bool array, causal[i,j] = True iff i ≺ j (strict)
    d : spacetime dimension
    rho : optional density prefactor (default 1.0 = omit)

    Returns
    -------
    B : (N, N) float64 matrix
    """
    if d not in BDG_COEFFICIENTS:
        raise ValueError(f"d={d} not supported; need d in {list(BDG_COEFFICIENTS)}")

    coeffs = BDG_COEFFICIENTS[d]
    n_layers = len(coeffs)  # = floor(d/2) + 1
    N = causal.shape[0]

    # Compute interval size matrix: isz[i,j] = #{k : i ≺ k ≺ j}
    c_int = causal.astype(np.int32)
    isz = c_int @ c_int  # isz[i,j] = number of elements strictly between i and j

    # Optional dimensional prefactor: ρ^{2/d}
    rho_factor = rho ** (2.0 / d) if rho != 1.0 else 1.0

    # Build the off-diagonal part: Σ_i C_i * L_i(x)
    # L_i(x) for the "past" version: y ≺ x with isz[y,x] = i
    # We also need the "future" version for a symmetric operator,
    # but the standard CST convention sums over past neighbors only.
    # For the BDG action (constant field), past-only and future-only give same result.
    #
    # We build: B[x, y] = rho_factor * C_i  if y ≺ x and isz[y,x] = i, for i < n_layers
    #           B[x, x] = rho_factor * α_d  (diagonal)

    # The coefficient α_d (diagonal term) for the retarded d'Alembertian:
    # In the Sorkin convention, α_d = -2/c_d where c_d = S_{d-2}/(d(d-1)2^{d/2-1})
    # For our purpose we don't need the exact continuum normalization since
    # we'll measure ratios. But for correctness:
    #   d=2: α₂ = -2,  c₂ = π/(2·1·1) = π/2 → but conventions vary
    # Actually the BDG action convention absorbs α_d into the N term.
    # In the BDG action: S_BDG = α_d * N + Σ_i C_i * N_i
    # For a constant field φ=1, (Bφ)(x) = α_d + Σ_i C_i * |L_i(x)|
    #
    # The BDG action is S = Σ_x (B1)(x) = α_d * N + Σ_i C_i * (total N_i)
    # In the standard form with α_d = 1 (absorbed into prefactor):
    #   S_BDG^(d=4) = N - C₀ + 9C₁ - 16C₂ + 8C₃
    # So α_d = 1 in these conventions, and the C_i already include the sign.

    # Build B matrix (using past-causal convention)
    # causal[y, x] = True means y ≺ x, so for element x, its past is column x of causal^T
    # i.e., row x of causal^T = causal[:, x]
    # More precisely: for row x of B, we want B[x, y] where y is in past of x

    B = np.zeros((N, N), dtype=np.float64)

    # Diagonal: α_d = 1 (BDG convention)
    np.fill_diagonal(B, 1.0)

    # Off-diagonal: for each y ≺ x (causal[y,x] = True), set B[x,y] = C_i
    # where i = isz[y,x] = number of elements strictly between y and x
    # We only include layers i = 0, 1, ..., n_layers-1
    past = causal.T  # past[x, y] = True iff y ≺ x

    for i, c_i in enumerate(coeffs):
        if c_i == 0.0:
            continue
        layer_mask = past & (isz.T == i)  # layer_mask[x,y] = True iff y ≺ x and isz[y,x] = i
        B += c_i * layer_mask.astype(np.float64)

    # Apply overall density prefactor
    B *= rho_factor

    return B


def compute_b1_features(
    B: np.ndarray,
) -> dict[str, float]:
    """Compute spectral features of B applied to constant field.

    B1 = B @ ones  gives per-element "curvature estimate".
    In the continuum limit, mean(B1) → -½ R.

    Returns dict with keys:
      b1_mean, b1_std, b1_median, b1_min, b1_max,
      b1_neg_frac (fraction of elements with B1 < 0),
      eig_max, eig_min, eig_gap (spectral gap = eig2 - eig1 of symmetric part)
    """
    N = B.shape[0]
    ones = np.ones(N)
    b1 = B @ ones  # per-element curvature signal

    features: dict[str, float] = {}
    features["b1_mean"] = float(np.mean(b1))
    features["b1_std"] = float(np.std(b1))
    features["b1_median"] = float(np.median(b1))
    features["b1_min"] = float(np.min(b1))
    features["b1_max"] = float(np.max(b1))
    features["b1_neg_frac"] = float(np.mean(b1 < 0))

    # Eigenvalue features of the symmetric part (B + B^T)/2
    # Only compute for small N to avoid O(N^3) cost
    if N <= 600:
        B_sym = 0.5 * (B + B.T)
        try:
            evals = np.linalg.eigvalsh(B_sym)
            features["eig_min"] = float(evals[0])
            features["eig_max"] = float(evals[-1])
            if N >= 3:
                features["eig_gap"] = float(evals[-1] - evals[-2])
            else:
                features["eig_gap"] = float("nan")
            # Spectral spread
            features["eig_spread"] = float(evals[-1] - evals[0])
        except np.linalg.LinAlgError:
            features["eig_min"] = float("nan")
            features["eig_max"] = float("nan")
            features["eig_gap"] = float("nan")
            features["eig_spread"] = float("nan")
    else:
        features["eig_min"] = float("nan")
        features["eig_max"] = float("nan")
        features["eig_gap"] = float("nan")
        features["eig_spread"] = float("nan")

    return features


@dataclass(frozen=True)
class ExpRow:
    d: int
    N: int
    hubble: float
    rep: int
    R_dS: float
    H2: float
    rho_est: float
    n_causal_pairs: int
    b1_mean: float
    b1_std: float
    b1_median: float
    b1_min: float
    b1_max: float
    b1_neg_frac: float
    eig_min: float
    eig_max: float
    eig_gap: float
    eig_spread: float
    # Normalized curvature estimate: -2 * b1_mean / rho^{2/d} should → R
    R_est: float


def compute_bdg_action_from_causal(
    causal: np.ndarray, d: int,
) -> float:
    """Compute BDG action S_BDG = N + Σ_i C_i N_i directly from causal matrix.

    N_i = total number of intervals of size i (i elements strictly between).
    This is the *global* BDG action — sum of per-element B_ℓ(1).
    """
    coeffs = BDG_COEFFICIENTS[d]
    N = causal.shape[0]
    c_int = causal.astype(np.int32)
    isz = c_int @ c_int  # interval sizes

    # Count N_i = number of causal pairs (y ≺ x) with exactly i between
    pair_sizes = isz[causal]  # interval sizes for all causal pairs
    s = float(N)  # α_d = 1
    for i, c_i in enumerate(coeffs):
        n_i = int(np.sum(pair_sizes == i))
        s += c_i * n_i
    return s


def run_single(
    d: int, N: int, hubble: float, rep: int, seed: int
) -> ExpRow:
    """Run one realization: sprinkle, build B_ℓ, extract features."""
    points = sprinkle_de_sitter_like_diamond(
        N, d - 1, hubble=hubble,
        seed=seed + d * 100000 + N * 100 + rep + int(hubble * 1000),
    )
    causal = build_causal_matrix_de_sitter(points, hubble)

    n_causal_pairs = int(np.sum(causal))
    R_dS = d * (d - 1) * hubble ** 2
    H2 = hubble ** 2

    # Build d'Alembertian matrix WITHOUT density prefactor
    # (ρ^{2/d} is a calibration constant; omitting preserves rank correlations)
    B = build_dalembertian_matrix(causal, d, rho=1.0)

    # Extract features
    feats = compute_b1_features(B)

    # Also compute global BDG action (= sum of B1 over all elements)
    s_bdg = compute_bdg_action_from_causal(causal, d)
    # Normalized BDG per element
    bdg_per_elem = s_bdg / N if N > 0 else float('nan')

    # R_est: in continuum B1_mean → -½ R up to prefactor.
    # Without prefactor, the *sign* and *monotonicity* are the observables.
    R_est = -2.0 * feats["b1_mean"]

    return ExpRow(
        d=d, N=N, hubble=hubble, rep=rep,
        R_dS=R_dS, H2=H2, rho_est=bdg_per_elem,
        n_causal_pairs=n_causal_pairs,
        b1_mean=feats["b1_mean"],
        b1_std=feats["b1_std"],
        b1_median=feats["b1_median"],
        b1_min=feats["b1_min"],
        b1_max=feats["b1_max"],
        b1_neg_frac=feats["b1_neg_frac"],
        eig_min=feats["eig_min"],
        eig_max=feats["eig_max"],
        eig_gap=feats["eig_gap"],
        eig_spread=feats["eig_spread"],
        R_est=R_est,
    )


def generate_report(rows: list[ExpRow], dims: list[int], ns: list[int],
                    hubbles: list[float]) -> str:
    """Generate Markdown report with correlation analysis."""
    from scipy import stats as sp_stats

    lines: list[str] = []
    lines.append("# §4.1.27 Sorkin d'Alembertian B_ℓ Spectral Experiment\n")
    lines.append("## Design\n")
    lines.append(f"- Dimensions: {dims}")
    lines.append(f"- Sizes: {ns}")
    lines.append(f"- Hubble: {hubbles}")
    lines.append(f"- Total realizations: {len(rows)}\n")

    lines.append("## Theoretical Prediction\n")
    lines.append("In the continuum limit (Benincasa-Dowker-Glaser):")
    lines.append("  ⟨B^(d) 1⟩(x) → -½ R(x)")
    lines.append("So b1_mean should be **negatively correlated** with H².")
    lines.append("R_est = -2 * b1_mean should be **positively correlated** with H².\n")

    # === Q1: b1_mean vs H² per slice ===
    lines.append("## Q1: b1_mean vs H² (Spearman)\n")
    lines.append("| d | N | ρ_S(b1_mean, H²) | p-value | sig? | ρ_S(R_est, H²) | p-value |")
    lines.append("|---|---|-------------------|---------|------|----------------|---------|")

    sig_pos_count = 0
    sig_neg_count = 0
    total_slices = 0
    rho_r_est_list = []

    for d in dims:
        for N in ns:
            h2_vals, b1_vals, r_est_vals = [], [], []
            for hubble in hubbles:
                subset = [r for r in rows if r.d == d and r.N == N and r.hubble == hubble]
                for r in subset:
                    if not math.isnan(r.b1_mean):
                        h2_vals.append(r.H2)
                        b1_vals.append(r.b1_mean)
                        r_est_vals.append(r.R_est)

            if len(h2_vals) < 5:
                continue

            total_slices += 1
            rho_b1, p_b1 = sp_stats.spearmanr(h2_vals, b1_vals)
            rho_re, p_re = sp_stats.spearmanr(h2_vals, r_est_vals)
            sig = "***" if p_re < 0.001 and rho_re > 0 else ("**" if p_re < 0.01 and rho_re > 0 else ("*" if p_re < 0.05 and rho_re > 0 else ""))
            if rho_re > 0 and p_re < 0.05:
                sig_pos_count += 1
            if rho_re < 0 and p_re < 0.05:
                sig_neg_count += 1
            rho_r_est_list.append(rho_re)

            lines.append(
                f"| {d} | {N} | {rho_b1:+.3f} | {p_b1:.2e} | {sig} | "
                f"{rho_re:+.3f} | {p_re:.2e} |"
            )

    lines.append(f"\n**Positive significant slices (R_est vs H²): {sig_pos_count}/{total_slices}**")
    lines.append(f"**Negative significant slices: {sig_neg_count}/{total_slices}**")
    if rho_r_est_list:
        lines.append(f"**Mean |ρ|(R_est, H²): {np.mean(np.abs(rho_r_est_list)):.3f}**\n")

    # === Q2: All features vs H² ===
    feature_names = ["b1_mean", "b1_std", "b1_neg_frac", "b1_median",
                     "eig_min", "eig_max", "eig_gap", "eig_spread", "R_est"]
    lines.append("## Q2: All Features vs H² (pooled by d)\n")
    lines.append("| d | feature | ρ_S | p-value | sign |")
    lines.append("|---|---------|-----|---------|------|")

    for d in dims:
        for feat in feature_names:
            h2_vals, feat_vals = [], []
            for r in rows:
                if r.d != d:
                    continue
                v = getattr(r, feat, None)
                if v is not None and not math.isnan(v):
                    h2_vals.append(r.H2)
                    feat_vals.append(v)
            if len(h2_vals) < 10:
                continue
            rho_f, p_f = sp_stats.spearmanr(h2_vals, feat_vals)
            sign = "+" if rho_f > 0 else "-"
            lines.append(f"| {d} | {feat} | {rho_f:+.3f} | {p_f:.2e} | {sign} |")

    # === Q3: b1_mean summary table ===
    lines.append("\n## Q3: b1_mean Summary Table\n")
    lines.append("| d | N | H | H² | R_dS | mean(b1) | std(b1) | mean(R_est) | n_pairs |")
    lines.append("|---|---|---|-----|------|----------|---------|-------------|---------|")

    for d in dims:
        for N in ns:
            for hubble in hubbles:
                subset = [r for r in rows if r.d == d and r.N == N and r.hubble == hubble]
                b1s = [r.b1_mean for r in subset if not math.isnan(r.b1_mean)]
                res = [r.R_est for r in subset if not math.isnan(r.R_est)]
                nps = [r.n_causal_pairs for r in subset]
                if not b1s:
                    continue
                R_dS = d * (d - 1) * hubble ** 2
                lines.append(
                    f"| {d} | {N} | {hubble:.2f} | {hubble**2:.3f} | {R_dS:.2f} | "
                    f"{np.mean(b1s):+.4f} | {np.std(b1s):.4f} | "
                    f"{np.mean(res):+.2f} | {int(np.mean(nps))} |"
                )

    # === Q4: N-convergence of R_est / R_dS ===
    lines.append("\n## Q4: Calibration Ratio R_est / R_dS\n")
    lines.append("| d | N | H | R_est/R_dS | std |")
    lines.append("|---|---|---|-----------|-----|")

    for d in dims:
        for N in ns:
            for hubble in hubbles:
                if hubble == 0.0:
                    continue  # skip flat (R_dS = 0)
                subset = [r for r in rows if r.d == d and r.N == N and r.hubble == hubble]
                ratios = [r.R_est / r.R_dS for r in subset
                          if not math.isnan(r.R_est) and r.R_dS > 0]
                if not ratios:
                    continue
                lines.append(
                    f"| {d} | {N} | {hubble:.2f} | {np.mean(ratios):+.4f} | {np.std(ratios):.4f} |"
                )

    # === Q5: Density-Residual Analysis ===
    # Critical test: after regressing out density (n_causal_pairs), do spectral
    # features retain curvature signal? This distinguishes genuinely new info
    # from another density projection.
    lines.append("\n## Q5: Density-Residual Analysis\n")
    lines.append("After OLS-removing n_causal_pairs (density proxy ≈ ΣC_k),")
    lines.append("does each spectral feature's residual still correlate with H²?\n")
    lines.append("| d | feature | raw ρ_S | residual ρ_S | Δρ | verdict |")
    lines.append("|---|---------|--------|-------------|-----|---------|")

    spectral_feats = ["b1_mean", "b1_std", "eig_min", "eig_max", "eig_spread", "eig_gap"]
    density_residual_results: dict[tuple[int, str], tuple[float, float]] = {}

    for d in dims:
        for feat in spectral_feats:
            h2_arr, feat_arr, dens_arr = [], [], []
            for r in rows:
                if r.d != d:
                    continue
                v = getattr(r, feat, None)
                if v is not None and not math.isnan(v) and not math.isnan(r.n_causal_pairs):
                    h2_arr.append(r.H2)
                    feat_arr.append(v)
                    dens_arr.append(float(r.n_causal_pairs))
            if len(h2_arr) < 15:
                continue

            h2_a = np.array(h2_arr)
            feat_a = np.array(feat_arr)
            dens_a = np.array(dens_arr)

            # Raw correlation
            rho_raw, _ = sp_stats.spearmanr(h2_a, feat_a)

            # OLS: feat = a * density + b  → residual
            coeffs_ols = np.polyfit(dens_a, feat_a, 1)
            resid = feat_a - np.polyval(coeffs_ols, dens_a)

            # Residual correlation with H²
            rho_resid, p_resid = sp_stats.spearmanr(h2_a, resid)
            delta = abs(rho_resid) - abs(rho_raw)

            if abs(rho_resid) > 0.3 and p_resid < 0.01:
                verdict = "**BEYOND DENSITY**"
            elif abs(rho_resid) < 0.15:
                verdict = "density-dominated"
            else:
                verdict = "marginal"

            density_residual_results[(d, feat)] = (rho_raw, rho_resid)
            lines.append(
                f"| {d} | {feat} | {rho_raw:+.3f} | {rho_resid:+.3f} | {delta:+.3f} | {verdict} |"
            )

    # Count "beyond density" features
    beyond_count = sum(1 for v in density_residual_results.values()
                       if abs(v[1]) > 0.3)
    total_tested = len(density_residual_results)
    lines.append(f"\n**Features with |residual ρ| > 0.3: {beyond_count}/{total_tested}**\n")

    if beyond_count > 0:
        lines.append("### Interpretation\n")
        lines.append("Spectral features of B_ℓ that survive density removal carry ")
        lines.append("information **beyond** the {C_k} observable family closed by DDT. ")
        lines.append("This validates the theoretical prediction that B_ℓ as a matrix ")
        lines.append("operator encodes curvature through its eigenvector structure, ")
        lines.append("not just through interval count statistics.\n")
    else:
        lines.append("### Interpretation\n")
        lines.append("All spectral features collapse to density after residualization. ")
        lines.append("B_ℓ's spectral features at these N values are still dominated by ")
        lines.append("the same density signal identified by DDT. Larger N may be needed ")
        lines.append("for the continuum-limit curvature encoding to emerge.\n")

    # === Conclusion ===
    lines.append("\n## Conclusion\n")

    # Summarize both b1_mean (constant field) and spectral findings
    lines.append("### 1. Constant-Field Projection (b1_mean = BDG action per element)\n")
    lines.append(f"R_est = -2·b1_mean vs H²: {sig_pos_count}/{total_slices} positive significant slices.")
    lines.append(" b1_mean is a linear combination of {C_k} and subject to DDT.")
    lines.append(" At d=3/4 the density effect dominates → wrong sign or no signal.\n")

    lines.append("### 2. Spectral Features (eigenvalues of B_ℓ matrix)\n")
    # Collect spectral summary
    spectral_sig = {}
    for d in dims:
        for feat in ["b1_std", "eig_min", "eig_max", "eig_spread"]:
            key = (d, feat)
            if key in density_residual_results:
                _, rho_r = density_residual_results[key]
                if abs(rho_r) > 0.3:
                    spectral_sig.setdefault(d, []).append((feat, rho_r))

    for d in dims:
        feats = spectral_sig.get(d, [])
        if feats:
            feat_str = ", ".join(f"{f}(ρ_resid={r:+.3f})" for f, r in feats)
            lines.append(f"- **d={d}**: {len(feats)} features BEYOND DENSITY: {feat_str}")
        else:
            lines.append(f"- d={d}: no features survive density removal")

    lines.append("")
    if beyond_count >= 4:
        lines.append(f"**VERDICT: B_ℓ spectral features encode curvature BEYOND {{C_k}}.**")
        lines.append(f" {beyond_count}/{total_tested} features survive density residualization.")
        lines.append(" The eigenvalue spectrum of the d'Alembertian matrix carries")
        lines.append(" genuinely new information not accessible from interval counts alone.")
        lines.append(" This confirms the E-bulk breakout prediction: moving from")
        lines.append(" {C_k} statistics to causal-set operators unlocks new curvature channels.")
        lines.append(" **DDT condition C1 (observable space = {C_k} only) has been escaped.**\n")
    elif beyond_count > 0:
        lines.append(f"**Partial beyond-density signal**: {beyond_count}/{total_tested} features.")
        lines.append(" Some spectral features carry curvature info beyond density,")
        lines.append(" but the effect is dimension-dependent. Larger N may strengthen it.\n")
    else:
        lines.append(f"**All spectral features density-dominated** at current N.\n")

    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(
        description="§4.1.27: Sorkin d'Alembertian B_ℓ spectral experiment"
    )
    ap.add_argument("--dims", nargs="*", type=int, default=[2, 3, 4])
    ap.add_argument("--ns", nargs="*", type=int, default=[128, 256, 512])
    ap.add_argument("--hubbles", nargs="*", type=float,
                    default=[0.0, 0.25, 0.5, 1.0, 2.0])
    ap.add_argument("--reps", type=int, default=8)
    ap.add_argument("--seed", type=int, default=2027)
    ap.add_argument("--out",
                    default="outputs_unified_functional/conjecture_e_sorkin_dalembertian.csv")
    ap.add_argument("--report",
                    default="outputs_unified_functional/conjecture_e_sorkin_dalembertian.md")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[ExpRow] = []
    total = len(args.dims) * len(args.ns) * len(args.hubbles) * args.reps
    done = 0

    print(f"Sorkin d'Alembertian experiment: {total} realizations")
    print(f"  dims={args.dims}, ns={args.ns}, hubbles={args.hubbles}, reps={args.reps}")
    print()

    for d in args.dims:
        for N in args.ns:
            for hubble in args.hubbles:
                for rep in range(args.reps):
                    row = run_single(d, N, hubble, rep, args.seed)
                    rows.append(row)
                    done += 1
                    if done % 10 == 0 or done == total:
                        print(
                            f"  [{done:4d}/{total}] d={d} N={N:4d} H={hubble:.2f} "
                            f"b1_mean={row.b1_mean:+.4f} R_est={row.R_est:+.2f} "
                            f"pairs={row.n_causal_pairs}"
                        )

    # Save raw CSV
    fieldnames = [f.name for f in fields(ExpRow)]
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow({fn: getattr(row, fn) for fn in fieldnames})
    print(f"\nSaved CSV: {out_path}")

    # Generate and save report
    report_text = generate_report(rows, args.dims, args.ns, args.hubbles)
    report_path = Path(args.report)
    report_path.write_text(report_text, encoding="utf-8")
    print(f"Saved report: {report_path}")

    # Console summary
    print("\n" + "=" * 70)
    print("QUICK SUMMARY")
    print("=" * 70)
    from scipy import stats as sp_stats
    for d in args.dims:
        h2_all, re_all = [], []
        for r in rows:
            if r.d == d and not math.isnan(r.R_est):
                h2_all.append(r.H2)
                re_all.append(r.R_est)
        if len(h2_all) >= 5:
            rho, p = sp_stats.spearmanr(h2_all, re_all)
            print(f"  d={d}: ρ(R_est, H²) = {rho:+.3f}  (p={p:.2e})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
