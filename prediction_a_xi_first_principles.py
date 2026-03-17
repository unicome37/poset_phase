"""
Prediction A — First-principles Xi derivation from Minkowski light-cone geometry
===============================================================================

Goal
----
Replace the purely empirical scaling-law fit with a geometry-first derivation:

1. Compute the ordering fraction p_d directly from the pair-difference
   distribution induced by the cube-sprinkle generator.
2. Compute the expected link density C0/N from the Alexandrov-volume emptiness
   probability:

       P(link | separation) ~= (1 - V_diamond)^ (N-2)

   with the flat-space interval volume V_diamond = kappa_d * tau^d.
3. Derive theory-side scaling parameters (a_d, alpha_d) from this geometric
   prediction, instead of fitting them to the observed C0/N data.
4. Use a minimal mean-field entropy closure,

       logH / N ~= (1 - p_d) * (log N - 1),

   to obtain theory-side (b_d, c_d) with no data fit:

       b_d = 1 - p_d,   c_d = -(1 - p_d).

This is not an exact theorem for the finite cube generator. It is a bulk
Minkowski approximation meant to test whether the 4->5 barrier already follows
from light-cone thinning plus empty-interval statistics.

Update
------
The minimal entropy closure is useful but systematically underestimates the
4->5 barrier. We therefore also include a two-constant renormalized closure,

    logH / N ~= A + B * (1 - p_d) * (log N - 1),

where A and B are global constants shared across all tested dimensions.
This still avoids per-dimension fitting of (b_d, c_d), while capturing the
finite-generator normalization missed by the raw mean-field expression.

An even more geometric alternative is to replace the free intercept by the
link-sparsity correction:

    logH / N ~= B0 * (1 - p_d) * (log N - 1) + B1 * (1 - ell_d),

where ell_d is the theory-side link fraction. This keeps two global
coefficients, but both multiply explicit geometric observables rather than
an abstract additive offset.
"""

from __future__ import annotations

import math
import pathlib
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

OUT_DIR = pathlib.Path("outputs_exploratory/prediction_a_xi_first_principles")
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_VALUES = np.array([20, 36, 52, 68, 80, 96, 112], dtype=int)
DIMS = [2, 3, 4, 5]


def unit_ball_volume(m: int) -> float:
    return math.pi ** (m / 2.0) / math.gamma(m / 2.0 + 1.0)


def alexandrov_volume_constant(d: int) -> float:
    """Volume of a flat d-dimensional Alexandrov interval of proper time tau:
    V_diamond(tau) = kappa_d * tau^d.
    """
    return unit_ball_volume(d - 1) / (2 ** (d - 1) * d)


def sample_pair_differences(
    d: int,
    n_samples: int = 400_000,
    seed: int = 12345,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Sample |dt| and spatial differences from the exact cube-difference law."""
    rng = np.random.default_rng(seed + d * 1000)
    a = rng.random((n_samples, d))
    b = rng.random((n_samples, d))
    diff = b - a
    dt = np.abs(diff[:, 0])
    spatial = diff[:, 1:]
    if spatial.size == 0:
        r = np.zeros(n_samples)
    else:
        r = np.sqrt(np.sum(spatial ** 2, axis=1))
    tau2 = np.clip(dt * dt - r * r, 0.0, None)
    tau = np.sqrt(tau2)
    related = dt >= r
    return dt, r, tau, related


def geometric_summary(
    d: int,
    n_samples: int = 400_000,
    seed: int = 12345,
) -> dict:
    dt, r, tau, related = sample_pair_differences(d, n_samples=n_samples, seed=seed)
    p_d = float(np.mean(related))
    tau_related = tau[related]
    kappa_d = alexandrov_volume_constant(d)

    return {
        "d": d,
        "kappa_d": kappa_d,
        "p_d": p_d,
        "tau_mean_related": float(np.mean(tau_related)) if tau_related.size else 0.0,
        "tau_d_mean_related": float(np.mean(tau_related ** d)) if tau_related.size else 0.0,
        "samples": n_samples,
    }


def theory_link_density_curve(
    d: int,
    n_values: np.ndarray,
    n_samples: int = 400_000,
    seed: int = 12345,
) -> pd.DataFrame:
    """Bulk Minkowski prediction for C0/N from empty Alexandrov intervals."""
    _, _, tau, related = sample_pair_differences(d, n_samples=n_samples, seed=seed)
    kappa_d = alexandrov_volume_constant(d)
    vol = kappa_d * tau[related] ** d

    rows = []
    for n in n_values:
        emptiness = np.power(np.clip(1.0 - vol, 0.0, 1.0), n - 2)
        pair_link_prob = float(np.mean(emptiness) * np.mean(related))
        c0_per_n = 0.5 * (n - 1) * pair_link_prob
        relation_per_n = 0.5 * (n - 1) * float(np.mean(related))
        link_fraction = c0_per_n / relation_per_n if relation_per_n > 0 else 0.0
        rows.append(
            {
                "d": d,
                "n": int(n),
                "C0_per_N_theory": c0_per_n,
                "R_per_N_theory": relation_per_n,
                "link_fraction_theory": link_fraction,
            }
        )
    return pd.DataFrame(rows)


def fit_power_law_from_theory(curve: pd.DataFrame) -> tuple[float, float]:
    ns = curve["n"].to_numpy(dtype=float)
    ys = curve["C0_per_N_theory"].to_numpy(dtype=float)
    logn = np.log(ns)
    logy = np.log(ys)
    slope, intercept = np.polyfit(logn, logy, 1)
    return float(math.exp(intercept)), float(slope)


def entropy_params_from_p(p_d: float) -> tuple[float, float]:
    b_d = 1.0 - p_d
    c_d = -(1.0 - p_d)
    return b_d, c_d


def fit_entropy_renormalization(obs_df: pd.DataFrame, params_df: pd.DataFrame) -> tuple[float, float]:
    merged = obs_df.merge(params_df[["d", "p_d"]], on="d", how="left").copy()
    merged["logN"] = np.log(merged["n"])
    merged["h0"] = (1.0 - merged["p_d"]) * (merged["logN"] - 1.0)
    X = np.c_[np.ones(len(merged)), merged["h0"].to_numpy()]
    y = merged["logH_per_N_obs"].to_numpy()
    A, B = np.linalg.lstsq(X, y, rcond=None)[0]
    return float(A), float(B)


def fit_entropy_geometric_mixture(obs_df: pd.DataFrame, curve_df: pd.DataFrame, params_df: pd.DataFrame) -> tuple[float, float]:
    merged = (
        obs_df.merge(curve_df[["n", "d", "link_fraction_theory"]], on=["n", "d"], how="left")
        .merge(params_df[["d", "p_d"]], on="d", how="left")
        .copy()
    )
    merged["logN"] = np.log(merged["n"])
    merged["h0"] = (1.0 - merged["p_d"]) * (merged["logN"] - 1.0)
    merged["sparsity"] = 1.0 - merged["link_fraction_theory"]
    X = np.c_[merged["h0"].to_numpy(), merged["sparsity"].to_numpy()]
    y = merged["logH_per_N_obs"].to_numpy()
    B0, B1 = np.linalg.lstsq(X, y, rcond=None)[0]
    return float(B0), float(B1)


def load_observed_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(
        "outputs_exploratory/prediction_a_large_n_scaling/raw_observables_large_n.csv"
    )
    agg = (
        raw.groupby(["n", "dim"])
        .agg(C0=("C0_links", "mean"), logH=("log_H", "mean"))
        .reset_index()
    )
    agg["d"] = agg["dim"].map({"2d": 2, "3d": 3, "4d": 4, "5d": 5})
    agg["C0_per_N_obs"] = agg["C0"] / agg["n"]
    agg["logH_per_N_obs"] = agg["logH"] / agg["n"]

    xi = pd.read_csv("outputs_exploratory/prediction_a_large_n_scaling/xi_large_n.csv")
    return agg, xi


def build_first_principles_tables(
    n_samples: int = 400_000,
    seed: int = 12345,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    geom_rows = []
    curve_parts = []
    param_rows = []

    for d in DIMS:
        geom = geometric_summary(d, n_samples=n_samples, seed=seed)
        geom_rows.append(geom)
        curve = theory_link_density_curve(d, N_VALUES, n_samples=n_samples, seed=seed)
        curve_parts.append(curve)
        a_d, alpha_d = fit_power_law_from_theory(curve)
        b_d, c_d = entropy_params_from_p(geom["p_d"])
        param_rows.append(
            {
                "d": d,
                "a_d_theory": a_d,
                "alpha_d_theory": alpha_d,
                "b_d_theory": b_d,
                "c_d_theory": c_d,
                "p_d": geom["p_d"],
                "kappa_d": geom["kappa_d"],
            }
        )

    geom_df = pd.DataFrame(geom_rows)
    curve_df = pd.concat(curve_parts, ignore_index=True)
    params_df = pd.DataFrame(param_rows)
    return geom_df, curve_df, params_df


def build_xi_table(
    curve_df: pd.DataFrame,
    params_df: pd.DataFrame,
    entropy_mode: str = "minimal",
    entropy_renorm: tuple[float, float] | None = None,
    entropy_geom_mix: tuple[float, float] | None = None,
) -> pd.DataFrame:
    curve_map = {
        (int(row["d"]), int(row["n"])): row["C0_per_N_theory"]
        for _, row in curve_df.iterrows()
    }
    p_map = {int(row["d"]): row["p_d"] for _, row in params_df.iterrows()}
    ell_map = {
        (int(row["d"]), int(row["n"])): row["link_fraction_theory"]
        for _, row in curve_df.iterrows()
    }

    rows = []
    for d_lo, d_hi in [(2, 3), (3, 4), (4, 5)]:
        for n in N_VALUES:
            c_lo = curve_map[(d_lo, int(n))]
            c_hi = curve_map[(d_hi, int(n))]
            delta_s = 2.0 * abs(c_lo - c_hi)
            h0_lo = (1.0 - p_map[d_lo]) * (math.log(n) - 1.0)
            h0_hi = (1.0 - p_map[d_hi]) * (math.log(n) - 1.0)
            if entropy_mode == "minimal":
                h_lo = h0_lo
                h_hi = h0_hi
            elif entropy_mode == "renormalized":
                if entropy_renorm is None:
                    raise ValueError("entropy_renorm must be provided for renormalized mode")
                A, B = entropy_renorm
                h_lo = A + B * h0_lo
                h_hi = A + B * h0_hi
            elif entropy_mode == "geom_mix":
                if entropy_geom_mix is None:
                    raise ValueError("entropy_geom_mix must be provided for geom_mix mode")
                B0, B1 = entropy_geom_mix
                h_lo = B0 * h0_lo + B1 * (1.0 - ell_map[(d_lo, int(n))])
                h_hi = B0 * h0_hi + B1 * (1.0 - ell_map[(d_hi, int(n))])
            else:
                raise ValueError(f"Unknown entropy_mode: {entropy_mode}")
            delta_h = abs(h_hi - h_lo)
            xi = delta_s / delta_h if delta_h > 0 else float("inf")
            rows.append(
                {
                    "entropy_mode": entropy_mode,
                    "dim_pair": f"{d_lo}→{d_hi}",
                    "n": int(n),
                    "delta_S_theory": delta_s,
                    "delta_H_theory": delta_h,
                    "Xi_theory": xi,
                }
            )
    return pd.DataFrame(rows)


def build_report(
    geom_df: pd.DataFrame,
    params_df: pd.DataFrame,
    xi_min_df: pd.DataFrame,
    xi_renorm_df: pd.DataFrame,
    xi_geom_df: pd.DataFrame,
    xi_obs: pd.DataFrame,
    entropy_renorm: tuple[float, float],
    entropy_geom_mix: tuple[float, float],
) -> str:
    A, B = entropy_renorm
    B0, B1 = entropy_geom_mix
    lines: list[str] = []
    lines.append("=" * 78)
    lines.append("FIRST-PRINCIPLES XI DERIVATION FROM MINKOWSKI LIGHT-CONE GEOMETRY")
    lines.append("=" * 78)
    lines.append("")
    lines.append("Model assumptions")
    lines.append("  1. Pair differences are sampled exactly from the cube-sprinkle generator.")
    lines.append("  2. Link probability uses the flat Alexandrov volume emptiness factor:")
    lines.append("       P_link(Δ) ~= (1 - kappa_d * tau^d)^(N-2).")
    lines.append("  3. Entropy uses the mean-field closure:")
    lines.append("       logH/N ~= (1 - p_d) * (log N - 1).")
    lines.append("")
    lines.append("Geometric constants")
    for _, row in geom_df.iterrows():
        lines.append(
            f"  d={int(row['d'])}: p_d={row['p_d']:.5f}, "
            f"kappa_d={row['kappa_d']:.6f}, "
            f"<tau^d>_related={row['tau_d_mean_related']:.6f}"
        )
    lines.append("")
    lines.append("Theory-side scaling parameters")
    for _, row in params_df.iterrows():
        lines.append(
            f"  d={int(row['d'])}: a_d={row['a_d_theory']:.4f}, "
            f"alpha_d={row['alpha_d_theory']:.4f}, "
            f"b_d={row['b_d_theory']:.4f}, c_d={row['c_d_theory']:.4f}"
        )
    lines.append("")
    lines.append("Entropy renormalization")
    lines.append(
        f"  logH/N ~= A + B * (1 - p_d) * (log N - 1), "
        f"with A={A:.6f}, B={B:.6f}"
    )
    lines.append(
        f"  logH/N ~= B0 * (1 - p_d) * (log N - 1) + B1 * (1 - ell_d), "
        f"with B0={B0:.6f}, B1={B1:.6f}"
    )
    lines.append("")
    lines.append("Xi comparison (theory vs observed)")
    xi_obs_45 = xi_obs[xi_obs["dim_pair"] == "4→5"].sort_values("n")
    xi_min_45 = xi_min_df[xi_min_df["dim_pair"] == "4→5"].sort_values("n")
    xi_renorm_45 = xi_renorm_df[xi_renorm_df["dim_pair"] == "4→5"].sort_values("n")
    xi_geom_45 = xi_geom_df[xi_geom_df["dim_pair"] == "4→5"].sort_values("n")
    for ((_, th_min), (_, th_renorm), (_, th_geom), (_, obs)) in zip(
        xi_min_45.iterrows(),
        xi_renorm_45.iterrows(),
        xi_geom_45.iterrows(),
        xi_obs_45.iterrows(),
    ):
        lines.append(
            f"  N={int(obs['n']):>3}: Xi_min={th_min['Xi_theory']:.2f}, "
            f"Xi_renorm={th_renorm['Xi_theory']:.2f}, "
            f"Xi_geom={th_geom['Xi_theory']:.2f}, Xi_obs={obs['Xi']:.2f}"
        )
    med_min = float(xi_min_45["Xi_theory"].median())
    med_renorm = float(xi_renorm_45["Xi_theory"].median())
    med_geom = float(xi_geom_45["Xi_theory"].median())
    med_obs = float(xi_obs_45["Xi"].median())
    rel_err_min = abs(med_min - med_obs) / med_obs if med_obs else float("nan")
    rel_err_renorm = abs(med_renorm - med_obs) / med_obs if med_obs else float("nan")
    rel_err_geom = abs(med_geom - med_obs) / med_obs if med_obs else float("nan")
    lines.append("")
    lines.append(
        f"  Median Xi_45 (minimal)={med_min:.2f}, observed={med_obs:.2f}, "
        f"relative_error={rel_err_min:.1%}"
    )
    lines.append(
        f"  Median Xi_45 (renormalized)={med_renorm:.2f}, observed={med_obs:.2f}, "
        f"relative_error={rel_err_renorm:.1%}"
    )
    lines.append(
        f"  Median Xi_45 (geom_mix)={med_geom:.2f}, observed={med_obs:.2f}, "
        f"relative_error={rel_err_geom:.1%}"
    )
    lines.append("")
    lines.append("Interpretation")
    lines.append("  - The numerator is now geometry-driven: no empirical C0/N fit is used.")
    lines.append("  - The denominator is controlled by p_d alone; the 4->5 entropy gap shrinks")
    lines.append("    because p_d already becomes small by d=4.")
    lines.append("  - The raw mean-field closure is too stiff; a single global renormalization")
    lines.append("    already absorbs most of the finite-generator mismatch.")
    lines.append("  - Replacing the additive intercept by link sparsity produces an even tighter")
    lines.append("    geometry-side closure, suggesting that the missing entropy correction tracks")
    lines.append("    the fraction of causal relations that are still non-link mediated.")
    return "\n".join(lines)


def make_figures(
    curve_df: pd.DataFrame,
    params_df: pd.DataFrame,
    obs_df: pd.DataFrame,
    xi_min_df: pd.DataFrame,
    xi_renorm_df: pd.DataFrame,
    xi_geom_df: pd.DataFrame,
    xi_obs: pd.DataFrame,
) -> None:
    dim_colors = {2: "#1f77b4", 3: "#ff7f0e", 4: "#d62728", 5: "#2ca02c"}
    pair_colors = {"2→3": "#1f77b4", "3→4": "#ff7f0e", "4→5": "#d62728"}

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    ax = axes[0]
    for d in DIMS:
        sub_obs = obs_df[obs_df["d"] == d].sort_values("n")
        sub_th = curve_df[curve_df["d"] == d].sort_values("n")
        ax.plot(
            sub_obs["n"],
            sub_obs["C0_per_N_obs"],
            "o",
            color=dim_colors[d],
            label=f"{d}D observed",
        )
        ax.plot(
            sub_th["n"],
            sub_th["C0_per_N_theory"],
            "-",
            color=dim_colors[d],
            alpha=0.8,
            label=f"{d}D theory",
        )
    ax.set_xlabel("N")
    ax.set_ylabel("C0 / N")
    ax.set_title("Link Density: theory vs observed")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2)

    ax = axes[1]
    for pair in ["2→3", "3→4", "4→5"]:
        sub_min = xi_min_df[xi_min_df["dim_pair"] == pair].sort_values("n")
        sub_renorm = xi_renorm_df[xi_renorm_df["dim_pair"] == pair].sort_values("n")
        sub_geom = xi_geom_df[xi_geom_df["dim_pair"] == pair].sort_values("n")
        sub_obs = xi_obs[xi_obs["dim_pair"] == pair].sort_values("n")
        ax.plot(sub_min["n"], sub_min["Xi_theory"], "--", color=pair_colors[pair], alpha=0.55, label=f"{pair} minimal")
        ax.plot(sub_renorm["n"], sub_renorm["Xi_theory"], "-", color=pair_colors[pair], label=f"{pair} renorm")
        ax.plot(sub_geom["n"], sub_geom["Xi_theory"], ":", color=pair_colors[pair], linewidth=2.4, label=f"{pair} geom")
        ax.plot(sub_obs["n"], sub_obs["Xi"], "o", color=pair_colors[pair], alpha=0.75, label=f"{pair} observed")
    ax.set_xlabel("N")
    ax.set_ylabel("Xi")
    ax.set_title("First-principles Xi vs observed")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "first_principles_overview.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    pd_vals = params_df.sort_values("d")["p_d"].to_numpy()
    dims = params_df.sort_values("d")["d"].to_numpy()
    ax.plot(dims, pd_vals, "o-", color="#d62728", linewidth=2)
    for d, p_d in zip(dims, pd_vals):
        ax.annotate(f"{p_d:.3f}", (d, p_d), textcoords="offset points", xytext=(8, 6))
    ax.set_xlabel("Dimension d")
    ax.set_ylabel("p_d")
    ax.set_title("Ordering fraction from geometry")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ordering_fraction_first_principles.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    geom_df, curve_df, params_df = build_first_principles_tables()
    obs_df, xi_obs = load_observed_data()
    entropy_renorm = fit_entropy_renormalization(obs_df, params_df)
    entropy_geom_mix = fit_entropy_geometric_mixture(obs_df, curve_df, params_df)
    xi_min_df = build_xi_table(curve_df, params_df, entropy_mode="minimal")
    xi_renorm_df = build_xi_table(
        curve_df,
        params_df,
        entropy_mode="renormalized",
        entropy_renorm=entropy_renorm,
    )
    xi_geom_df = build_xi_table(
        curve_df,
        params_df,
        entropy_mode="geom_mix",
        entropy_geom_mix=entropy_geom_mix,
    )

    merged_curve = curve_df.merge(
        obs_df[["n", "d", "C0_per_N_obs", "logH_per_N_obs"]],
        on=["n", "d"],
        how="left",
    )
    merged_xi_min = xi_min_df.merge(
        xi_obs[["n", "dim_pair", "Xi"]],
        on=["n", "dim_pair"],
        how="left",
    )
    merged_xi_renorm = xi_renorm_df.merge(
        xi_obs[["n", "dim_pair", "Xi"]],
        on=["n", "dim_pair"],
        how="left",
    )
    merged_xi_geom = xi_geom_df.merge(
        xi_obs[["n", "dim_pair", "Xi"]],
        on=["n", "dim_pair"],
        how="left",
    )

    report = build_report(
        geom_df,
        params_df,
        xi_min_df,
        xi_renorm_df,
        xi_geom_df,
        xi_obs,
        entropy_renorm,
        entropy_geom_mix,
    )

    geom_df.to_csv(OUT_DIR / "geometric_constants.csv", index=False)
    merged_curve.to_csv(OUT_DIR / "link_density_theory_vs_observed.csv", index=False)
    params_df.to_csv(OUT_DIR / "scaling_parameters_first_principles.csv", index=False)
    merged_xi_min.to_csv(OUT_DIR / "xi_first_principles_vs_observed.csv", index=False)
    merged_xi_renorm.to_csv(OUT_DIR / "xi_first_principles_renorm_vs_observed.csv", index=False)
    merged_xi_geom.to_csv(OUT_DIR / "xi_first_principles_geom_mix_vs_observed.csv", index=False)
    pd.DataFrame(
        [{"A_entropy": entropy_renorm[0], "B_entropy": entropy_renorm[1]}]
    ).to_csv(OUT_DIR / "entropy_renormalization.csv", index=False)
    pd.DataFrame(
        [{"B0_geom_mix": entropy_geom_mix[0], "B1_geom_mix": entropy_geom_mix[1]}]
    ).to_csv(OUT_DIR / "entropy_geometric_mixture.csv", index=False)
    (OUT_DIR / "first_principles_report.txt").write_text(report, encoding="utf-8")

    make_figures(curve_df, params_df, obs_df, xi_min_df, xi_renorm_df, xi_geom_df, xi_obs)

    print(report)
    print("")
    print(f"Outputs saved to: {OUT_DIR}")


if __name__ == "__main__":
    main()
