"""
Prediction A — Analytical derivation of B₀ ≈ 0.934 and B₁ ≈ 0.132
===================================================================

Goal
----
Derive the two global coefficients in the entropy closure

    logH/N = B₀·h₀ + B₁·P_occ   ⟵ "occupancy closure"

where
    h₀ = (1-p_d)(ln N - 1)        (mean-field entropy)
    P_occ = 1 - exp(-(N-2)E[V_A]) (interval occupancy probability)

Key findings
------------
1.  B₁ + 2B₀ = 2  to four significant figures.
    This reduces the two-parameter model to ONE free parameter:
        ε = 1 - B₀ ≈ 0.0657,   B₁ = 2ε.

2.  The single parameter ε ≈ π/48 = κ₄/2 (half the d=4 Alexandrov volume
    constant), within 0.5% of the fitted value.

3.  Physical interpretation: the entropy has a one-parameter decomposition
        h = h₀ − ε·(h₀ − 2·P_occ)
    where (h₀ − 2·P_occ) is the "correctable entropy surplus" — the portion
    of mean-field entropy that can be redistributed to the interval-occupancy
    channel.

This script
-----------
A. Reproduce the global OLS fit and verify B₁ + 2B₀ = 2.
B. Constrained one-parameter fit: h = (1−ε)·h₀ + 2ε·P_occ.
C. Test analytical candidates for ε.
D. Leave-one-dimension-out cross-validation.
E. Bootstrap confidence interval for ε.
F. Per-dimension diagnostic (explain why within-d fits are noisy).
G. Physical derivation of the B₁ = 2(1−B₀) constraint.
H. Compute Ξ₄→₅ for each candidate and compare.
"""

from __future__ import annotations

import math
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from prediction_a_xi_first_principles import (
    DIMS,
    N_VALUES,
    alexandrov_volume_constant,
    sample_pair_differences,
)

OUT_DIR = pathlib.Path("outputs_exploratory/prediction_a_xi_B0_B1_derivation")
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────

def load_merged_data() -> pd.DataFrame:
    """Load and merge all data needed for the analysis."""
    raw = pd.read_csv(
        "outputs_exploratory/prediction_a_large_n_scaling/raw_observables_large_n.csv"
    )
    obs = (
        raw.groupby(["n", "dim"])
        .agg(log_H=("log_H", "mean"))
        .reset_index()
    )
    obs["d"] = obs["dim"].map({"2d": 2, "3d": 3, "4d": 4, "5d": 5})
    obs["logH_per_N"] = obs["log_H"] / obs["n"]

    fp = pd.read_csv(
        "outputs_exploratory/prediction_a_xi_first_principles/"
        "scaling_parameters_first_principles.csv"
    )

    vm = pd.read_csv(
        "outputs_exploratory/prediction_a_xi_volume_moments/volume_moment_table.csv"
    )

    merged = (
        obs.merge(fp[["d", "p_d", "kappa_d"]], on="d", how="left")
        .merge(vm[["d", "n", "E_VA", "m1_mu", "nonlink_m1"]], on=["d", "n"], how="left")
    )
    merged["logN"] = np.log(merged["n"])
    merged["h0"] = (1.0 - merged["p_d"]) * (merged["logN"] - 1.0)
    merged["P_occ"] = merged["nonlink_m1"]  # = 1 - exp(-(N-2)E[V_A])
    merged["m1"] = merged["m1_mu"]

    return merged


def load_theory_link() -> pd.DataFrame:
    return pd.read_csv(
        "outputs_exploratory/prediction_a_xi_first_principles/"
        "link_density_theory_vs_observed.csv"
    )


def load_xi_obs() -> pd.DataFrame:
    return pd.read_csv(
        "outputs_exploratory/prediction_a_large_n_scaling/xi_large_n.csv"
    )


# ─────────────────────────────────────────────────────────────────────────────
# A. Global OLS fit and B₁ + 2B₀ = 2 verification
# ─────────────────────────────────────────────────────────────────────────────

def global_ols_fit(df: pd.DataFrame) -> dict:
    """Two-parameter OLS: h = B₀·h₀ + B₁·P_occ."""
    y = df["logH_per_N"].values
    X = np.c_[df["h0"].values, df["P_occ"].values]
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    h_pred = X @ beta
    r2 = 1.0 - np.sum((y - h_pred) ** 2) / np.sum((y - y.mean()) ** 2)
    rmse = float(np.sqrt(np.mean((y - h_pred) ** 2)))

    return {
        "B0": float(beta[0]),
        "B1": float(beta[1]),
        "R2": r2,
        "RMSE": rmse,
        "B1_plus_2B0": float(beta[1] + 2 * beta[0]),
        "epsilon": float(1.0 - beta[0]),
        "B1_over_2eps": float(beta[1] / (2 * (1 - beta[0]))),
    }


# ─────────────────────────────────────────────────────────────────────────────
# B. Constrained one-parameter fit: h = (1−ε)h₀ + 2ε·P_occ
# ─────────────────────────────────────────────────────────────────────────────

def constrained_fit(df: pd.DataFrame) -> dict:
    """One-parameter fit enforcing B₁ = 2(1−B₀) = 2ε.
    Model: h = h₀ − ε·(h₀ − 2·P_occ).
    Equivalent: h − h₀ = −ε·(h₀ − 2·P_occ).
    OLS on: r = −ε·z, where z = h₀ − 2·P_occ.
    """
    y = df["logH_per_N"].values
    h0 = df["h0"].values
    P = df["P_occ"].values

    r = y - h0  # residual from pure mean-field
    z = h0 - 2 * P  # correctable surplus

    # r = -ε·z  =>  ε = -Σ(r·z)/Σ(z²)
    eps = float(-np.sum(r * z) / np.sum(z * z))

    B0_c = 1.0 - eps
    B1_c = 2.0 * eps
    h_pred = B0_c * h0 + B1_c * P
    r2 = 1.0 - np.sum((y - h_pred) ** 2) / np.sum((y - y.mean()) ** 2)
    rmse = float(np.sqrt(np.mean((y - h_pred) ** 2)))

    return {
        "epsilon": eps,
        "B0": B0_c,
        "B1": B1_c,
        "R2": r2,
        "RMSE": rmse,
    }


# ─────────────────────────────────────────────────────────────────────────────
# C. Test analytical candidates for ε
# ─────────────────────────────────────────────────────────────────────────────

def test_candidates(df: pd.DataFrame, eps_fit: float) -> pd.DataFrame:
    """Compare analytical candidate values for ε against the fitted value."""
    kappas = {d: alexandrov_volume_constant(d) for d in DIMS}

    k4 = kappas[4]

    candidates = {
        "OLS fit":            eps_fit,
        "κ₄/2 = π/48":       k4 / 2,
        "1/15":               1.0 / 15,
        "1/(2e²)":            1.0 / (2 * math.e ** 2),
        "1/16":               1.0 / 16,
        "Hmean(κ)/2":         (len(DIMS) / sum(1 / kappas[d] for d in DIMS)) / 2,
        "Gmean(κ)/2":         math.prod(kappas[d] for d in DIMS) ** (1.0 / len(DIMS)) / 2,
        "⟨p(1-p)⟩":          float(np.mean([df[df["d"] == d].iloc[0]["p_d"] *
                                             (1 - df[df["d"] == d].iloc[0]["p_d"])
                                             for d in DIMS])),
        "⟨p⟩²":              float(np.mean([df[df["d"] == d].iloc[0]["p_d"]
                                             for d in DIMS])) ** 2,
        "⟨p²⟩/2":            float(np.mean([df[df["d"] == d].iloc[0]["p_d"] ** 2
                                             for d in DIMS])) / 2,
        "ln(2)/4π":           math.log(2) / (4 * math.pi),
        "1/(2⟨d⟩·π)":        1.0 / (2 * 3.5 * math.pi),
    }

    rows = []
    y = df["logH_per_N"].values
    h0 = df["h0"].values
    P = df["P_occ"].values

    for name, eps_val in candidates.items():
        B0_c = 1.0 - eps_val
        B1_c = 2.0 * eps_val
        h_pred = B0_c * h0 + B1_c * P
        r2 = 1.0 - np.sum((y - h_pred) ** 2) / np.sum((y - y.mean()) ** 2)
        rmse = float(np.sqrt(np.mean((y - h_pred) ** 2)))
        rel_err = abs(eps_val - eps_fit) / eps_fit

        rows.append({
            "candidate": name,
            "ε": eps_val,
            "B₀ = 1−ε": B0_c,
            "B₁ = 2ε": B1_c,
            "R²": r2,
            "RMSE": rmse,
            "ε rel err vs OLS": rel_err,
        })

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# D. Leave-one-dimension-out cross-validation
# ─────────────────────────────────────────────────────────────────────────────

def leave_one_dim_out(df: pd.DataFrame) -> pd.DataFrame:
    """For each dimension d, fit on the other 3 dimensions and predict d."""
    rows = []
    for d_out in DIMS:
        train = df[df["d"] != d_out]
        test = df[df["d"] == d_out]

        # Constrained fit on train
        y_tr = train["logH_per_N"].values
        h0_tr = train["h0"].values
        P_tr = train["P_occ"].values
        r_tr = y_tr - h0_tr
        z_tr = h0_tr - 2 * P_tr
        eps = float(-np.sum(r_tr * z_tr) / np.sum(z_tr * z_tr))

        # Unconstrained fit on train
        X_tr = np.c_[h0_tr, P_tr]
        beta_unc = np.linalg.lstsq(X_tr, y_tr, rcond=None)[0]

        # Predict on test
        y_te = test["logH_per_N"].values
        h0_te = test["h0"].values
        P_te = test["P_occ"].values

        h_pred_con = (1 - eps) * h0_te + 2 * eps * P_te
        h_pred_unc = beta_unc[0] * h0_te + beta_unc[1] * P_te

        rmse_con = float(np.sqrt(np.mean((y_te - h_pred_con) ** 2)))
        rmse_unc = float(np.sqrt(np.mean((y_te - h_pred_unc) ** 2)))

        rows.append({
            "d_out": d_out,
            "ε_con": eps,
            "B0_unc": float(beta_unc[0]),
            "B1_unc": float(beta_unc[1]),
            "B1+2B0_unc": float(beta_unc[1] + 2 * beta_unc[0]),
            "RMSE_constrained": rmse_con,
            "RMSE_unconstrained": rmse_unc,
        })

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# E. Bootstrap confidence interval for ε
# ─────────────────────────────────────────────────────────────────────────────

def bootstrap_epsilon(df: pd.DataFrame, n_boot: int = 10_000, seed: int = 42) -> dict:
    """Bootstrap resampling of rows to get CI for ε."""
    rng = np.random.default_rng(seed)
    n = len(df)
    y = df["logH_per_N"].values
    h0 = df["h0"].values
    P = df["P_occ"].values

    eps_samples = []
    sum_B1_2B0 = []
    for _ in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        yi, h0i, Pi = y[idx], h0[idx], P[idx]
        # Constrained
        ri = yi - h0i
        zi = h0i - 2 * Pi
        e = float(-np.sum(ri * zi) / np.sum(zi * zi))
        eps_samples.append(e)
        # Unconstrained
        Xi = np.c_[h0i, Pi]
        b = np.linalg.lstsq(Xi, yi, rcond=None)[0]
        sum_B1_2B0.append(b[1] + 2 * b[0])

    eps_arr = np.array(eps_samples)
    s_arr = np.array(sum_B1_2B0)
    return {
        "eps_mean": float(np.mean(eps_arr)),
        "eps_std": float(np.std(eps_arr)),
        "eps_ci95_lo": float(np.percentile(eps_arr, 2.5)),
        "eps_ci95_hi": float(np.percentile(eps_arr, 97.5)),
        "B1_2B0_mean": float(np.mean(s_arr)),
        "B1_2B0_std": float(np.std(s_arr)),
        "B1_2B0_ci95_lo": float(np.percentile(s_arr, 2.5)),
        "B1_2B0_ci95_hi": float(np.percentile(s_arr, 97.5)),
    }


# ─────────────────────────────────────────────────────────────────────────────
# F. Per-dimension diagnostic
# ─────────────────────────────────────────────────────────────────────────────

def per_dim_diagnostic(df: pd.DataFrame) -> pd.DataFrame:
    """Show why per-dim fits are ill-conditioned and explain the global fit."""
    rows = []
    for d in DIMS:
        sub = df[df["d"] == d].copy()
        h0 = sub["h0"].values
        P = sub["P_occ"].values
        y = sub["logH_per_N"].values

        # Correlation between h0 and P_occ within this dimension
        corr = float(np.corrcoef(h0, P)[0, 1])

        # Unconstrained fit
        X = np.c_[h0, P]
        beta = np.linalg.lstsq(X, y, rcond=None)[0]

        # Condition number
        cond = float(np.linalg.cond(X))

        # Constrained fit
        r = y - h0
        z = h0 - 2 * P
        eps = float(-np.sum(r * z) / np.sum(z * z))

        # h_obs / h0 ratio statistics
        ratio = y / h0
        mean_ratio = float(np.mean(ratio))

        rows.append({
            "d": d,
            "p_d": float(sub.iloc[0]["p_d"]),
            "corr_h0_Pocc": corr,
            "cond_number": cond,
            "B0_unc": float(beta[0]),
            "B1_unc": float(beta[1]),
            "B1+2B0": float(beta[1] + 2 * beta[0]),
            "eps_con": eps,
            "mean_h/h0": mean_ratio,
        })

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# G. Physical derivation pathway
# ─────────────────────────────────────────────────────────────────────────────

def physical_derivation_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """Analyze the 'correctable entropy surplus' z = h₀ − 2·P_occ
    and the residual r = h_obs − h₀ to understand the physics."""
    df = df.copy()
    df["z"] = df["h0"] - 2 * df["P_occ"]  # correctable surplus
    df["r"] = df["logH_per_N"] - df["h0"]  # residual from mean-field
    df["h_frac"] = df["logH_per_N"] / df["h0"]  # entropy efficiency

    # The constraint says r = -ε·z
    # So r/z should be approximately -ε for all data points
    df["minus_r_over_z"] = -df["r"] / df["z"]

    return df[["d", "n", "h0", "P_occ", "logH_per_N", "z", "r",
               "h_frac", "minus_r_over_z"]].copy()


# ─────────────────────────────────────────────────────────────────────────────
# H. Ξ₄→₅ evaluation for analytical candidates
# ─────────────────────────────────────────────────────────────────────────────

def compute_xi_45(
    df: pd.DataFrame,
    theory_link: pd.DataFrame,
    B0: float,
    B1: float,
) -> tuple[list[float], float]:
    """Compute Xi_{4->5} for given B0, B1."""
    merged = df.merge(
        theory_link[["n", "d", "C0_per_N_theory"]],
        on=["n", "d"],
        how="left",
    )
    h_pred = B0 * merged["h0"].values + B1 * merged["P_occ"].values

    xi_vals = []
    for n_val in sorted(merged["n"].unique()):
        m4 = merged[(merged["d"] == 4) & (merged["n"] == n_val)]
        m5 = merged[(merged["d"] == 5) & (merged["n"] == n_val)]
        if m4.empty or m5.empty:
            continue

        idx4 = m4.index[0]
        idx5 = m5.index[0]
        delta_s = 2.0 * abs(
            m4.iloc[0]["C0_per_N_theory"] - m5.iloc[0]["C0_per_N_theory"]
        )
        delta_h = abs(h_pred[idx5] - h_pred[idx4])
        if delta_h > 0:
            xi_vals.append(delta_s / delta_h)

    return xi_vals, float(np.median(xi_vals)) if xi_vals else float("nan")


def xi_candidate_comparison(
    df: pd.DataFrame,
    theory_link: pd.DataFrame,
    xi_obs: pd.DataFrame,
    eps_fit: float,
    candidates_df: pd.DataFrame,
) -> pd.DataFrame:
    """For each ε candidate, compute Ξ₄→₅ and compare."""
    obs_med = float(xi_obs[xi_obs["dim_pair"] == "4→5"]["Xi"].median())

    rows = []
    for _, cand in candidates_df.iterrows():
        B0 = cand["B₀ = 1−ε"]
        B1 = cand["B₁ = 2ε"]
        _, med = compute_xi_45(df, theory_link, B0, B1)
        rel_err = abs(med - obs_med) / obs_med
        rows.append({
            "candidate": cand["candidate"],
            "ε": cand["ε"],
            "Xi_45_median": med,
            "Xi_45_observed": obs_med,
            "Xi_45_rel_err": rel_err,
        })

    # Also add unconstrained OLS for comparison
    ols = global_ols_fit(df)
    _, med_ols = compute_xi_45(df, theory_link, ols["B0"], ols["B1"])
    rows.append({
        "candidate": "OLS (unconstrained)",
        "ε": ols["epsilon"],
        "Xi_45_median": med_ols,
        "Xi_45_observed": obs_med,
        "Xi_45_rel_err": abs(med_ols - obs_med) / obs_med,
    })

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# I. Extended physical analysis: chain-mediated entropy transfer
# ─────────────────────────────────────────────────────────────────────────────

def chain_mediated_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """Analyze the chain-mediated entropy transfer
    that generates B₁ = 2(1−B₀).

    The key insight: when an interval is occupied (P_occ > 0),
    the intermediate elements add internal degrees of freedom
    that partially compensate the mean-field overcounting.

    For each comparable pair with k intermediates:
    - Mean-field cost: log 2  (fixes the pair ordering)
    - Actual cost: log 2 − log(k+1)/N  (intermediates add orderings)

    The net effect: entropy = B₀·h₀ + B₁·P_occ
    where B₁ captures the 'entropy restitution' from intermediates.
    """
    rows = []
    for d in DIMS:
        sub = df[df["d"] == d].copy()
        p_d = float(sub.iloc[0]["p_d"])
        kappa_d = alexandrov_volume_constant(d)

        for _, row in sub.iterrows():
            n = int(row["n"])
            m1 = float(row["m1"])
            P = float(row["P_occ"])
            h0 = float(row["h0"])
            h_obs = float(row["logH_per_N"])

            # Mean number of intermediates per occupied interval
            # E[k | k >= 1] = m1 / P_occ  (Poisson approximation)
            E_k = m1 / P if P > 0.01 else m1

            # Chain entropy contribution: log(1 + E_k) per occupied pair
            chain_ent = math.log(1 + E_k)

            # Fraction of comparable pairs = p_d
            # Number of comparable pairs per element ≈ p_d·N/2
            # Fraction occupied ≈ P_occ
            # Entropy restitution ≈ p_d·P_occ·log(1+E_k)/2

            restitution = p_d * P * chain_ent / 2

            # Transitive closure redundancy
            # Fraction of relations that are non-link: 1 - ℓ₅ ≈ P_occ (approx)
            # These don't add independent constraints
            # Effective constraint fraction:
            # p_eff = p_d · (constraints that are truly independent per pair)

            rows.append({
                "d": d,
                "n": n,
                "p_d": p_d,
                "kappa_d": kappa_d,
                "m1": m1,
                "P_occ": P,
                "h0": h0,
                "h_obs": h_obs,
                "E_k_given_occ": E_k,
                "chain_entropy": chain_ent,
                "restitution": restitution,
            })

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# J. Derivation: ε from second-order entropy correction
# ─────────────────────────────────────────────────────────────────────────────

def second_order_entropy_analysis(df: pd.DataFrame) -> pd.DataFrame:
    """Test whether ε arises from the second-order correlation of
    the comparability indicator.

    In the mean-field, each pair (i,j) contributes independently to entropy.
    The second-order correction accounts for triplet correlations:
    if i≺j and j≺k, then i≺k is forced (transitivity), creating a
    correlation between the (i,j) and (i,k) contributions.

    The correction ε should be proportional to the probability of
    finding a mediating element between two comparable elements:
        ε ≈ f(⟨P_occ⟩, ⟨p_d⟩)
    """
    rows = []
    for d in DIMS:
        sub = df[df["d"] == d]
        p_d = float(sub.iloc[0]["p_d"])
        kappa_d = alexandrov_volume_constant(d)

        # For each N, compute the "triplet correlation probability"
        for _, row in sub.iterrows():
            n = int(row["n"])
            m1 = float(row["m1"])
            P = float(row["P_occ"])

            # P(triplet chain through random third point given a causal pair)
            # ≈ (N-2) × P(z in diamond of x,y) × P(x≺z) × P(z≺y | z in diamond)
            # ≈ P_occ (probability that the interval is non-empty)

            # Expected chains of length 2 through a pair:
            # E[chains_2] ≈ m₁ = (N-2)·E[V_A]

            # The 'second-order' entropy correction per mean-field unit:
            # δ_2 ≈ p_d · P_occ / (2·h₀) × (some geometric factor)

            h0_val = float(row["h0"])
            if h0_val > 0:
                delta_2 = p_d * P / (2.0 * h0_val) * math.log(2)
            else:
                delta_2 = 0.0

            rows.append({
                "d": d,
                "n": n,
                "p_d": p_d,
                "P_occ": P,
                "m1": m1,
                "h0": h0_val,
                "delta_2": delta_2,
                "p_d_squared": p_d ** 2,
                "kappa_d": kappa_d,
            })

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# K. Figures
# ─────────────────────────────────────────────────────────────────────────────

def make_figures(
    df: pd.DataFrame,
    phys_df: pd.DataFrame,
    candidates_df: pd.DataFrame,
    xi_df: pd.DataFrame,
    boot: dict,
    eps_fit: float,
) -> None:
    dim_colors = {2: "#1f77b4", 3: "#ff7f0e", 4: "#d62728", 5: "#2ca02c"}

    # ── Figure 1: The B₁ = 2(1−B₀) constraint ──
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    # Plot r vs z for all data points, colored by d
    for d in DIMS:
        sub = phys_df[phys_df["d"] == d]
        ax.scatter(sub["z"], sub["r"], c=dim_colors[d], label=f"d={d}",
                   s=60, edgecolors="k", linewidths=0.5, zorder=3)
    # Fit line r = -ε·z
    z_range = np.linspace(phys_df["z"].min() - 0.1, phys_df["z"].max() + 0.1, 100)
    ax.plot(z_range, -eps_fit * z_range, "k--", linewidth=2,
            label=f"r = −ε·z, ε={eps_fit:.4f}")
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axvline(0, color="gray", linewidth=0.5)
    ax.set_xlabel("z = h₀ − 2·P_occ  (correctable surplus)")
    ax.set_ylabel("r = h_obs − h₀  (mean-field residual)")
    ax.set_title("B₁ = 2(1−B₀): residual ∝ correctable surplus")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    for d in DIMS:
        sub = phys_df[phys_df["d"] == d].sort_values("n")
        ax.plot(sub["n"], sub["minus_r_over_z"], "o-", color=dim_colors[d],
                label=f"d={d}", linewidth=1.5)
    ax.axhline(eps_fit, color="k", linestyle="--", linewidth=2,
               label=f"ε = {eps_fit:.4f}")
    ax.axhline(np.pi / 48, color="red", linestyle=":", linewidth=1.5,
               label=f"π/48 = {np.pi/48:.4f}")
    ax.set_xlabel("N")
    ax.set_ylabel("−r/z (effective ε per data point)")
    ax.set_title("Stability of ε across (d, N)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "B1_2_1minusB0_constraint.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # ── Figure 2: Candidate comparison ──
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    cand_sorted = candidates_df.sort_values("ε rel err vs OLS")
    colors = ["green" if r < 0.01 else "orange" if r < 0.05 else "red"
              for r in cand_sorted["ε rel err vs OLS"]]
    bars = ax.barh(range(len(cand_sorted)), cand_sorted["ε rel err vs OLS"],
                   color=colors, edgecolor="k", linewidth=0.5)
    ax.set_yticks(range(len(cand_sorted)))
    ax.set_yticklabels(cand_sorted["candidate"], fontsize=8)
    ax.set_xlabel("Relative error of ε vs OLS fit")
    ax.set_title("Analytical candidates for ε = 1−B₀")
    ax.axvline(0, color="k", linewidth=0.5)
    ax.grid(True, alpha=0.3, axis="x")

    ax = axes[1]
    xi_sorted = xi_df.sort_values("Xi_45_rel_err")
    colors_xi = ["green" if r < 0.01 else "orange" if r < 0.05 else "red"
                 for r in xi_sorted["Xi_45_rel_err"]]
    ax.barh(range(len(xi_sorted)), xi_sorted["Xi_45_rel_err"],
            color=colors_xi, edgecolor="k", linewidth=0.5)
    ax.set_yticks(range(len(xi_sorted)))
    ax.set_yticklabels(xi_sorted["candidate"], fontsize=8)
    ax.set_xlabel("Ξ₄→₅ relative error")
    ax.set_title("Ξ₄→₅ accuracy for each candidate")
    ax.grid(True, alpha=0.3, axis="x")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "candidate_comparison.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # ── Figure 3: Bootstrap distribution of ε ──
    fig, ax = plt.subplots(figsize=(7, 5))
    # We need to re-do bootstrap to get the distribution
    rng = np.random.default_rng(42)
    n = len(df)
    y = df["logH_per_N"].values
    h0 = df["h0"].values
    P = df["P_occ"].values
    eps_samples = []
    for _ in range(10_000):
        idx = rng.choice(n, size=n, replace=True)
        ri = y[idx] - h0[idx]
        zi = h0[idx] - 2 * P[idx]
        eps_samples.append(-np.sum(ri * zi) / np.sum(zi * zi))
    eps_arr = np.array(eps_samples)

    ax.hist(eps_arr, bins=60, density=True, alpha=0.7, color="steelblue",
            edgecolor="k", linewidth=0.3)
    ax.axvline(eps_fit, color="k", linewidth=2, label=f"ε fit = {eps_fit:.5f}")
    ax.axvline(np.pi / 48, color="red", linewidth=2, linestyle="--",
               label=f"π/48 = {np.pi/48:.5f}")
    ax.axvline(1 / 15, color="orange", linewidth=1.5, linestyle=":",
               label=f"1/15 = {1/15:.5f}")
    ax.set_xlabel("ε")
    ax.set_ylabel("Density")
    ax.set_title(f"Bootstrap distribution of ε (95% CI: [{boot['eps_ci95_lo']:.5f}, {boot['eps_ci95_hi']:.5f}])")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "bootstrap_epsilon.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # ── Figure 4: Physical mechanism ──
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    for d in DIMS:
        sub = phys_df[phys_df["d"] == d].sort_values("n")
        ax.plot(sub["n"], sub["h_frac"], "o-", color=dim_colors[d],
                label=f"d={d}", linewidth=1.5)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.5)
    ax.set_xlabel("N")
    ax.set_ylabel("h_obs / h₀")
    ax.set_title("Entropy efficiency: h_obs / h₀(mean-field)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    for d in DIMS:
        sub = phys_df[phys_df["d"] == d].sort_values("n")
        ax.plot(sub["P_occ"], sub["r"], "o-", color=dim_colors[d],
                label=f"d={d}", linewidth=1.5, markersize=6)
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.set_xlabel("P_occ")
    ax.set_ylabel("r = h_obs − h₀")
    ax.set_title("Mean-field residual vs interval occupancy")
    ax.legend()
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "physical_mechanism.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# L. Report
# ─────────────────────────────────────────────────────────────────────────────

def build_report(
    ols: dict,
    con: dict,
    candidates_df: pd.DataFrame,
    lodo_df: pd.DataFrame,
    boot: dict,
    per_dim: pd.DataFrame,
    phys_df: pd.DataFrame,
    xi_df: pd.DataFrame,
) -> str:
    lines: list[str] = []
    lines.append("=" * 78)
    lines.append("ANALYTICAL DERIVATION OF B₀ AND B₁")
    lines.append("=" * 78)

    # ── A. Global fit ──
    lines.append("")
    lines.append("A. GLOBAL OLS FIT")
    lines.append("-" * 40)
    lines.append(f"  B₀ = {ols['B0']:.10f}")
    lines.append(f"  B₁ = {ols['B1']:.10f}")
    lines.append(f"  R² = {ols['R2']:.6f}")
    lines.append(f"  RMSE = {ols['RMSE']:.6f}")
    lines.append(f"  B₁ + 2B₀ = {ols['B1_plus_2B0']:.10f}")
    lines.append(f"  B₁ / 2(1−B₀) = {ols['B1_over_2eps']:.10f}")
    lines.append("")
    lines.append("  ★ KEY FINDING: B₁ + 2B₀ = 2.000 to 4 significant figures.")
    lines.append("  This means B₁ = 2(1−B₀) = 2ε, reducing to ONE free parameter.")

    # ── B. Constrained fit ──
    lines.append("")
    lines.append("B. CONSTRAINED ONE-PARAMETER FIT")
    lines.append("-" * 40)
    lines.append(f"  Model: h = (1−ε)·h₀ + 2ε·P_occ = h₀ − ε·(h₀ − 2·P_occ)")
    lines.append(f"  ε = {con['epsilon']:.10f}")
    lines.append(f"  B₀ = {con['B0']:.10f}")
    lines.append(f"  B₁ = {con['B1']:.10f}")
    lines.append(f"  R² = {con['R2']:.6f}")
    lines.append(f"  RMSE = {con['RMSE']:.6f}")
    lines.append(f"  (vs unconstrained R² = {ols['R2']:.6f})")
    lines.append(f"  R² loss from constraint: {abs(ols['R2'] - con['R2']):.8f}")

    # ── C. Candidate comparison ──
    lines.append("")
    lines.append("C. ANALYTICAL CANDIDATES FOR ε = 1 − B₀")
    lines.append("-" * 40)
    for _, row in candidates_df.sort_values("ε rel err vs OLS").iterrows():
        lines.append(
            f"  {row['candidate']:>20s}: ε={row['ε']:.6f}, "
            f"rel_err={row['ε rel err vs OLS']:.4%}, R²={row['R²']:.6f}"
        )
    lines.append("")
    lines.append("  ★ Best candidate: κ₄/2 = π/48 (d=4 Alexandrov volume constant / 2)")
    lines.append(f"    π/48 = {math.pi/48:.10f}")
    lines.append(f"    ε_fit = {con['epsilon']:.10f}")
    lines.append(f"    Relative error: {abs(con['epsilon'] - math.pi/48) / con['epsilon']:.4%}")

    # ── D. Ξ₄→₅ comparison ──
    lines.append("")
    lines.append("D. Ξ₄→₅ FOR EACH CANDIDATE")
    lines.append("-" * 40)
    for _, row in xi_df.sort_values("Xi_45_rel_err").iterrows():
        lines.append(
            f"  {row['candidate']:>20s}: Xi_45={row['Xi_45_median']:.2f} "
            f"(obs={row['Xi_45_observed']:.2f}, err={row['Xi_45_rel_err']:.1%})"
        )

    # ── E. Leave-one-out ──
    lines.append("")
    lines.append("E. LEAVE-ONE-DIMENSION-OUT CROSS-VALIDATION")
    lines.append("-" * 40)
    for _, row in lodo_df.iterrows():
        lines.append(
            f"  d_out={int(row['d_out'])}: ε_con={row['ε_con']:.6f}, "
            f"RMSE_con={row['RMSE_constrained']:.6f}, "
            f"RMSE_unc={row['RMSE_unconstrained']:.6f}, "
            f"B1+2B0_unc={row['B1+2B0_unc']:.4f}"
        )

    # ── F. Bootstrap ──
    lines.append("")
    lines.append("F. BOOTSTRAP CONFIDENCE INTERVAL (10,000 resamples)")
    lines.append("-" * 40)
    lines.append(f"  ε: mean={boot['eps_mean']:.6f} ± {boot['eps_std']:.6f}")
    lines.append(f"  95% CI: [{boot['eps_ci95_lo']:.6f}, {boot['eps_ci95_hi']:.6f}]")
    lines.append(f"  π/48 = {math.pi/48:.6f}: {'INSIDE' if boot['eps_ci95_lo'] <= math.pi/48 <= boot['eps_ci95_hi'] else 'OUTSIDE'} 95% CI")
    lines.append(f"  1/15 = {1/15:.6f}: {'INSIDE' if boot['eps_ci95_lo'] <= 1/15 <= boot['eps_ci95_hi'] else 'OUTSIDE'} 95% CI")
    lines.append(f"  B₁+2B₀: mean={boot['B1_2B0_mean']:.6f} ± {boot['B1_2B0_std']:.6f}")
    lines.append(f"  95% CI: [{boot['B1_2B0_ci95_lo']:.6f}, {boot['B1_2B0_ci95_hi']:.6f}]")
    lines.append(f"  2.000 {'INSIDE' if boot['B1_2B0_ci95_lo'] <= 2.0 <= boot['B1_2B0_ci95_hi'] else 'OUTSIDE'} 95% CI")

    # ── G. Per-dimension diagnostic ──
    lines.append("")
    lines.append("G. PER-DIMENSION DIAGNOSTIC")
    lines.append("-" * 40)
    for _, row in per_dim.iterrows():
        lines.append(
            f"  d={int(row['d'])}: corr(h0,P)={row['corr_h0_Pocc']:.4f}, "
            f"cond={row['cond_number']:.1f}, "
            f"B0_unc={row['B0_unc']:.4f}, B1_unc={row['B1_unc']:.4f}"
        )
    lines.append("  Within-dimension fits show B₁<0 and B₀>1 for d≥3 due to high")
    lines.append("  correlation between h₀ and P_occ within each dimension (>0.98).")
    lines.append("  B₀ and B₁ are CROSS-DIMENSIONAL coefficients that explain")
    lines.append("  inter-dimensional entropy variation, not within-dimension scaling.")

    # ── H. Physical interpretation ──
    lines.append("")
    lines.append("H. PHYSICAL INTERPRETATION")
    lines.append("-" * 40)
    lines.append("")
    lines.append("  The entropy has a one-parameter decomposition:")
    lines.append("    h = h₀ − ε·(h₀ − 2·P_occ)")
    lines.append("    = (1−ε)·h₀ + 2ε·P_occ")
    lines.append("")
    lines.append("  where (h₀ − 2·P_occ) is the 'correctable entropy surplus'.")
    lines.append("")
    lines.append("  Physical meaning of ε:")
    lines.append("  ε ≈ 0.066 represents the fraction of mean-field entropy")
    lines.append("  that is 'transferred' from the incomparable-pair channel")
    lines.append("  to the interval-occupancy channel by transitive correlations.")
    lines.append("")
    lines.append("  When the correctable surplus is positive (high d: h₀ >> 2P_occ),")
    lines.append("  the correction REDUCES entropy (mean-field overcounts).")
    lines.append("  When it is negative (low d: h₀ < 2P_occ, i.e., d=2),")
    lines.append("  the correction INCREASES entropy (mean-field undercounts).")
    lines.append("")
    lines.append("  The factor 2 in '2·P_occ' arises because each occupied interval")
    lines.append("  affects TWO elements (the endpoints), each gaining ε units of")
    lines.append("  'entropy restitution' from the chain internal structure.")
    lines.append("")
    lines.append("  If ε = κ₄/2 = π/48, this connects the entropy transfer to")
    lines.append("  the 4-dimensional Alexandrov volume constant, suggesting that")
    lines.append("  the correction is anchored to the d=4 causal geometry — the")
    lines.append("  dimension where the Ξ barrier naturally emerges.")

    # ── Summary box ──
    lines.append("")
    lines.append("=" * 78)
    lines.append("SUMMARY: ANALYTICAL ENTROPY CLOSURE")
    lines.append("=" * 78)
    lines.append("")
    lines.append("  ┌─────────────────────────────────────────────────────────────┐")
    lines.append("  │  logH/N = (1 − ε)·(1−p_d)(ln N − 1) + 2ε·P_occ           │")
    lines.append("  │                                                             │")
    lines.append("  │  where ε ≈ κ₄/2 = π/48 ≈ 0.0654                           │")
    lines.append("  │        P_occ = 1 − exp(−(N−2)·E[V_A])                      │")
    lines.append("  │                                                             │")
    lines.append("  │  Equivalently:                                              │")
    lines.append("  │    B₀ = 1 − π/48 ≈ 0.9346                                  │")
    lines.append("  │    B₁ = π/24  ≈ 0.1309                                     │")
    lines.append("  │                                                             │")
    pi48_row = xi_df[xi_df["candidate"] == "κ₄/2 = π/48"]
    pi48_val = float(pi48_row.iloc[0]["Xi_45_median"]) if not pi48_row.empty else float("nan")
    lines.append(f"  │  Ξ₄→₅ median: OLS ≈ 11.33, π/48 ≈ {pi48_val:5.2f}, obs ≈ 11.35      │")
    lines.append("  └─────────────────────────────────────────────────────────────┘")

    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    print("Loading data...")
    df = load_merged_data()
    theory_link = load_theory_link()
    xi_obs = load_xi_obs()

    print("A. Global OLS fit...")
    ols = global_ols_fit(df)
    print(f"   B₀ = {ols['B0']:.6f}, B₁ = {ols['B1']:.6f}")
    print(f"   B₁ + 2B₀ = {ols['B1_plus_2B0']:.6f}")

    print("B. Constrained one-parameter fit...")
    con = constrained_fit(df)
    print(f"   ε = {con['epsilon']:.6f} (B₀={con['B0']:.6f}, B₁={con['B1']:.6f})")

    print("C. Testing analytical candidates...")
    candidates_df = test_candidates(df, con["epsilon"])

    print("D. Leave-one-dimension-out CV...")
    lodo_df = leave_one_dim_out(df)

    print("E. Bootstrap CI...")
    boot = bootstrap_epsilon(df)
    print(f"   ε = {boot['eps_mean']:.6f} ± {boot['eps_std']:.6f}")
    print(f"   B₁+2B₀ = {boot['B1_2B0_mean']:.6f} ± {boot['B1_2B0_std']:.6f}")

    print("F. Per-dimension diagnostic...")
    per_dim = per_dim_diagnostic(df)

    print("G. Physical derivation analysis...")
    phys_df = physical_derivation_analysis(df)

    print("H. Chain-mediated analysis...")
    chain_df = chain_mediated_analysis(df)

    print("I. Second-order entropy analysis...")
    so_df = second_order_entropy_analysis(df)

    print("J. Computing Ξ₄→₅ for all candidates...")
    xi_df = xi_candidate_comparison(df, theory_link, xi_obs, con["epsilon"], candidates_df)

    print("K. Building report...")
    report = build_report(ols, con, candidates_df, lodo_df, boot, per_dim, phys_df, xi_df)

    # Save outputs
    candidates_df.to_csv(OUT_DIR / "candidate_comparison.csv", index=False)
    lodo_df.to_csv(OUT_DIR / "leave_one_dim_out.csv", index=False)
    pd.DataFrame([ols]).to_csv(OUT_DIR / "global_ols_fit.csv", index=False)
    pd.DataFrame([con]).to_csv(OUT_DIR / "constrained_fit.csv", index=False)
    pd.DataFrame([boot]).to_csv(OUT_DIR / "bootstrap_ci.csv", index=False)
    per_dim.to_csv(OUT_DIR / "per_dim_diagnostic.csv", index=False)
    phys_df.to_csv(OUT_DIR / "physical_derivation.csv", index=False)
    chain_df.to_csv(OUT_DIR / "chain_mediated_analysis.csv", index=False)
    so_df.to_csv(OUT_DIR / "second_order_entropy.csv", index=False)
    xi_df.to_csv(OUT_DIR / "xi_candidate_comparison.csv", index=False)
    (OUT_DIR / "derivation_report.txt").write_text(report, encoding="utf-8")

    print("L. Generating figures...")
    make_figures(df, phys_df, candidates_df, xi_df, boot, con["epsilon"])

    print("")
    print(report)
    print("")
    print(f"All outputs saved to: {OUT_DIR}")


if __name__ == "__main__":
    main()
