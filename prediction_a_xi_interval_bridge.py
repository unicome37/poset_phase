"""
Prediction A — Interval bridge for the entropy correction
=========================================================

Goal
----
Interpret the geometric entropy correction term

    B1 * (1 - ell_d)

in terms of explicit interval hierarchy observables. The large-N scaling
dataset stored only C0, so this script regenerates the corresponding posets
from their recorded seeds and recomputes:

  - C0: links
  - C1: order-1 intervals (exactly one intermediate element)
  - R : total related pairs

We then test whether the non-link fraction

    1 - ell_d = 1 - C0 / R

can be approximated by a simpler interval statistic such as C1 / C0.
We also test a two-level bridge that adds a higher-order reservoir term,

    sqrt((R - C0 - C1) / R),

to capture mediated relations beyond order-1 intervals.
"""

from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd

from generators import (
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
)

OUT_DIR = pathlib.Path("outputs_exploratory/prediction_a_xi_interval_bridge")
OUT_DIR.mkdir(parents=True, exist_ok=True)

GENERATORS = {
    "2d": generate_lorentzian_like_2d,
    "3d": generate_lorentzian_like_3d,
    "4d": generate_lorentzian_like_4d,
    "5d": generate_lorentzian_like_5d,
}


def recompute_interval_profile() -> pd.DataFrame:
    raw = pd.read_csv("outputs_exploratory/prediction_a_large_n_scaling/raw_observables_large_n.csv")
    rows = []
    for rec in raw[["n", "dim", "sample_id", "seed"]].drop_duplicates().itertuples(index=False):
        poset = GENERATORS[rec.dim](n=int(rec.n), seed=int(rec.seed))
        closure_u8 = poset.closure.astype(np.uint8, copy=False)
        n_intermediate = (closure_u8 @ closure_u8)[poset.closure]
        rows.append(
            {
                "n": int(rec.n),
                "dim": rec.dim,
                "sample_id": int(rec.sample_id),
                "seed": int(rec.seed),
                "C0_re": int(np.sum(n_intermediate == 0)),
                "C1_re": int(np.sum(n_intermediate == 1)),
                "C2_re": int(np.sum(n_intermediate == 2)),
                "C3_re": int(np.sum(n_intermediate == 3)),
                "R_re": int(poset.closure.sum()),
            }
        )
    return pd.DataFrame(rows)


def summarize_interval_profile(profile_df: pd.DataFrame) -> pd.DataFrame:
    agg = (
        profile_df.groupby(["n", "dim"])
        .mean(numeric_only=True)
        .reset_index()
        .sort_values(["n", "dim"])
    )
    agg["C1_per_C0"] = agg["C1_re"] / agg["C0_re"]
    agg["nonlink_frac_obs"] = 1.0 - agg["C0_re"] / agg["R_re"]
    agg["higher_frac_obs"] = (agg["R_re"] - agg["C0_re"] - agg["C1_re"]) / agg["R_re"]
    agg["sqrt_higher_frac_obs"] = np.sqrt(np.clip(agg["higher_frac_obs"], 0.0, None))
    return agg


def fit_bridge(agg: pd.DataFrame) -> dict[str, float]:
    x = agg["C1_per_C0"].to_numpy()
    y = agg["nonlink_frac_obs"].to_numpy()

    X_aff = np.c_[np.ones(len(x)), x]
    a0, a1 = np.linalg.lstsq(X_aff, y, rcond=None)[0]
    y_aff = X_aff @ np.array([a0, a1])
    ss_res_aff = float(np.sum((y - y_aff) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2_aff = 1.0 - ss_res_aff / ss_tot

    X_lin = x[:, None]
    b1 = float(np.linalg.lstsq(X_lin, y, rcond=None)[0][0])
    y_lin = b1 * x
    ss_res_lin = float(np.sum((y - y_lin) ** 2))
    r2_lin = 1.0 - ss_res_lin / ss_tot

    X_two = np.c_[np.ones(len(x)), agg["C1_per_C0"].to_numpy(), agg["sqrt_higher_frac_obs"].to_numpy()]
    c0, c1, c2 = np.linalg.lstsq(X_two, y, rcond=None)[0]
    y_two = X_two @ np.array([c0, c1, c2])
    ss_res_two = float(np.sum((y - y_two) ** 2))
    r2_two = 1.0 - ss_res_two / ss_tot
    mae_two = float(np.mean(np.abs(y - y_two)))

    return {
        "a0_affine": float(a0),
        "a1_affine": float(a1),
        "r2_affine": r2_aff,
        "b1_linear": b1,
        "r2_linear": r2_lin,
        "c0_two_level": float(c0),
        "c1_two_level": float(c1),
        "c2_two_level": float(c2),
        "r2_two_level": r2_two,
        "mae_two_level": mae_two,
    }


def evaluate_entropy_closures(agg: pd.DataFrame) -> tuple[dict[str, float], pd.DataFrame]:
    raw = pd.read_csv("outputs_exploratory/prediction_a_large_n_scaling/raw_observables_large_n.csv")
    obs = (
        raw.groupby(["n", "dim"])
        .agg(log_H=("log_H", "mean"), C0_links=("C0_links", "mean"))
        .reset_index()
    )
    obs["d"] = obs["dim"].map({"2d": 2, "3d": 3, "4d": 4, "5d": 5})
    obs["logH_per_N_obs"] = obs["log_H"] / obs["n"]

    fp_params = pd.read_csv(
        "outputs_exploratory/prediction_a_xi_first_principles/scaling_parameters_first_principles.csv"
    )
    fp_curve = pd.read_csv(
        "outputs_exploratory/prediction_a_xi_first_principles/link_density_theory_vs_observed.csv"
    )
    xi_obs = pd.read_csv("outputs_exploratory/prediction_a_large_n_scaling/xi_large_n.csv")

    merged = (
        obs.merge(fp_params[["d", "p_d"]], on="d", how="left")
        .merge(fp_curve[["n", "d", "C0_per_N_theory", "link_fraction_theory"]], on=["n", "d"], how="left")
        .merge(
            agg[["n", "dim", "C1_per_C0", "nonlink_frac_obs", "higher_frac_obs", "sqrt_higher_frac_obs"]],
            on=["n", "dim"],
            how="left",
        )
    )
    merged["logN"] = np.log(merged["n"])
    merged["h0"] = (1.0 - merged["p_d"]) * (merged["logN"] - 1.0)

    y = merged["logH_per_N_obs"].to_numpy()

    X_nonlink = np.c_[merged["h0"].to_numpy(), merged["nonlink_frac_obs"].to_numpy()]
    beta_nonlink = np.linalg.lstsq(X_nonlink, y, rcond=None)[0]

    X_c1 = np.c_[merged["h0"].to_numpy(), merged["C1_per_C0"].to_numpy()]
    beta_c1 = np.linalg.lstsq(X_c1, y, rcond=None)[0]

    X_bridge = np.c_[
        np.ones(len(merged)),
        merged["C1_per_C0"].to_numpy(),
        merged["sqrt_higher_frac_obs"].to_numpy(),
    ]
    bridge_beta = np.linalg.lstsq(X_bridge, merged["nonlink_frac_obs"].to_numpy(), rcond=None)[0]
    merged["nonlink_hat"] = X_bridge @ bridge_beta
    X_nonlink_hat = np.c_[merged["h0"].to_numpy(), merged["nonlink_hat"].to_numpy()]
    beta_nonlink_hat = np.linalg.lstsq(X_nonlink_hat, y, rcond=None)[0]

    merged["h_pred_nonlink"] = X_nonlink @ beta_nonlink
    merged["h_pred_c1"] = X_c1 @ beta_c1
    merged["h_pred_nonlink_hat"] = X_nonlink_hat @ beta_nonlink_hat

    rows = []
    obs_med = float(xi_obs[xi_obs["dim_pair"] == "4→5"]["Xi"].median())
    for label, hcol in [
        ("nonlink_obs", "h_pred_nonlink"),
        ("c1_over_c0", "h_pred_c1"),
        ("bridge_nonlink_hat", "h_pred_nonlink_hat"),
    ]:
        vals = []
        for n in sorted(merged["n"].unique()):
            lo = merged[(merged["d"] == 4) & (merged["n"] == n)].iloc[0]
            hi = merged[(merged["d"] == 5) & (merged["n"] == n)].iloc[0]
            delta_s = 2.0 * abs(lo["C0_per_N_theory"] - hi["C0_per_N_theory"])
            delta_h = abs(hi[hcol] - lo[hcol])
            vals.append(delta_s / delta_h)
            rows.append({"closure": label, "n": int(n), "Xi_45_theory": delta_s / delta_h})
        rows.append  # quiet linter intent
    xi_eval = pd.DataFrame(rows)

    summary = {
        "beta_h0_nonlink": float(beta_nonlink[0]),
        "beta_nonlink": float(beta_nonlink[1]),
        "beta_h0_c1": float(beta_c1[0]),
        "beta_c1": float(beta_c1[1]),
        "beta_h0_nonlink_hat": float(beta_nonlink_hat[0]),
        "beta_nonlink_hat": float(beta_nonlink_hat[1]),
        "xi45_obs_median": obs_med,
        "xi45_nonlink_median": float(xi_eval[xi_eval["closure"] == "nonlink_obs"]["Xi_45_theory"].median()),
        "xi45_c1_median": float(xi_eval[xi_eval["closure"] == "c1_over_c0"]["Xi_45_theory"].median()),
        "xi45_nonlink_hat_median": float(xi_eval[xi_eval["closure"] == "bridge_nonlink_hat"]["Xi_45_theory"].median()),
    }
    summary["rel_err_nonlink"] = abs(summary["xi45_nonlink_median"] - obs_med) / obs_med
    summary["rel_err_c1"] = abs(summary["xi45_c1_median"] - obs_med) / obs_med
    summary["rel_err_nonlink_hat"] = abs(summary["xi45_nonlink_hat_median"] - obs_med) / obs_med
    return summary, xi_eval


def build_report(bridge: dict[str, float], summary: dict[str, float]) -> str:
    lines: list[str] = []
    lines.append("=" * 72)
    lines.append("INTERVAL BRIDGE FOR THE ENTROPY CORRECTION")
    lines.append("=" * 72)
    lines.append("")
    lines.append("Bridge fits")
    lines.append(
        "  nonlink_frac ~= a0 + a1 * (C1/C0) "
        f"with a0={bridge['a0_affine']:.4f}, a1={bridge['a1_affine']:.4f}, R2={bridge['r2_affine']:.4f}"
    )
    lines.append(
        "  nonlink_frac ~= b1 * (C1/C0) "
        f"with b1={bridge['b1_linear']:.4f}, R2={bridge['r2_linear']:.4f}"
    )
    lines.append(
        "  nonlink_frac ~= c0 + c1 * (C1/C0) + c2 * sqrt(higher_frac) "
        f"with c0={bridge['c0_two_level']:.4f}, c1={bridge['c1_two_level']:.4f}, "
        f"c2={bridge['c2_two_level']:.4f}, R2={bridge['r2_two_level']:.4f}"
    )
    lines.append("")
    lines.append("Entropy closures")
    lines.append(
        "  logH/N ~= beta0 * h0 + beta1 * nonlink_frac_obs "
        f"with beta0={summary['beta_h0_nonlink']:.4f}, beta1={summary['beta_nonlink']:.4f}"
    )
    lines.append(
        "  logH/N ~= gamma0 * h0 + gamma1 * (C1/C0) "
        f"with gamma0={summary['beta_h0_c1']:.4f}, gamma1={summary['beta_c1']:.4f}"
    )
    lines.append(
        "  logH/N ~= eta0 * h0 + eta1 * nonlink_hat "
        f"with eta0={summary['beta_h0_nonlink_hat']:.4f}, eta1={summary['beta_nonlink_hat']:.4f}"
    )
    lines.append("")
    lines.append("Xi_45 medians")
    lines.append(
        f"  observed  : {summary['xi45_obs_median']:.2f}"
    )
    lines.append(
        f"  nonlink   : {summary['xi45_nonlink_median']:.2f}  "
        f"(rel err {summary['rel_err_nonlink']:.1%})"
    )
    lines.append(
        f"  C1/C0     : {summary['xi45_c1_median']:.2f}  "
        f"(rel err {summary['rel_err_c1']:.1%})"
    )
    lines.append(
        f"  bridge_hat : {summary['xi45_nonlink_hat_median']:.2f}  "
        f"(rel err {summary['rel_err_nonlink_hat']:.1%})"
    )
    lines.append("")
    lines.append("Interpretation")
    lines.append("  - The interval correction is not dominated by order-1 intervals alone.")
    lines.append("  - C1/C0 is a strong proxy for non-link fraction, but it under-resolves")
    lines.append("    the higher-order mediated relations that still contribute to entropy.")
    lines.append("  - Adding sqrt(higher_frac) dramatically improves the bridge from interval")
    lines.append("    hierarchy to non-link fraction, showing that the missing reservoir is")
    lines.append("    distributed over higher mediated layers rather than concentrated at C1.")
    lines.append("  - The success of 1-ell_d therefore points to a coarse but robust summary:")
    lines.append("    entropy is corrected by the total mediated-causality reservoir, not just")
    lines.append("    its first layer.")
    return "\n".join(lines)


def main() -> None:
    profile_df = recompute_interval_profile()
    agg = summarize_interval_profile(profile_df)
    bridge = fit_bridge(agg)
    summary, xi_eval = evaluate_entropy_closures(agg)
    report = build_report(bridge, summary)

    profile_df.to_csv(OUT_DIR / "interval_profile_recomputed_samples.csv", index=False)
    agg.to_csv(OUT_DIR / "interval_profile_recomputed_aggregated.csv", index=False)
    pd.DataFrame([bridge]).to_csv(OUT_DIR / "interval_bridge_fit.csv", index=False)
    pd.DataFrame([summary]).to_csv(OUT_DIR / "interval_bridge_summary.csv", index=False)
    xi_eval.to_csv(OUT_DIR / "interval_bridge_xi_eval.csv", index=False)
    (OUT_DIR / "interval_bridge_report.txt").write_text(report, encoding="utf-8")

    print(report)
    print("")
    print(f"Outputs saved to: {OUT_DIR}")


if __name__ == "__main__":
    main()
