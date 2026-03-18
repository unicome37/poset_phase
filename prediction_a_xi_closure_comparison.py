"""
Prediction A — Unified closure comparison
=========================================

Collect in one place:
  1. alpha scan for P_occ(alpha) = 1 - exp(-alpha * m1)
  2. Xi_45 comparison table across the main closure candidates
"""

from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd

from prediction_a_xi_first_principles import (
    N_VALUES,
    alexandrov_volume_constant,
    sample_pair_differences,
)

OUT_DIR = pathlib.Path("outputs_exploratory/prediction_a_xi_closure_comparison")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def build_base_table() -> pd.DataFrame:
    rows = []
    for d in [2, 3, 4, 5]:
        _, _, tau, related = sample_pair_differences(d)
        vol = alexandrov_volume_constant(d) * tau[related] ** d
        for n in N_VALUES:
            x = (n - 2) * vol
            m1 = float(np.mean(x))
            var = float(np.var(x))
            ell_exact = float(np.mean((1.0 - vol) ** (n - 2)))
            ell_gamma = float((1.0 + var / m1) ** (-(m1 * m1 / var))) if (m1 > 0 and var > 0) else float("nan")
            ell2 = float(np.exp(-m1 + 0.5 * var))
            ell2_clip = min(1.0, max(0.0, ell2))
            rows.append(
                {
                    "d": d,
                    "n": int(n),
                    "m1": m1,
                    "var": var,
                    "P_occ": 1.0 - np.exp(-m1),
                    "one_minus_ell_mc": 1.0 - ell_exact,
                    "P_gamma": 1.0 - ell_gamma,
                    "one_minus_ell_2_clip": 1.0 - ell2_clip,
                }
            )
    return pd.DataFrame(rows)


def load_obs_frame() -> tuple[pd.DataFrame, float]:
    raw = pd.read_csv("outputs_exploratory/prediction_a_large_n_scaling/raw_observables_large_n.csv")
    obs = raw.groupby(["n", "dim"]).agg(log_H=("log_H", "mean")).reset_index()
    obs["d"] = obs["dim"].map({"2d": 2, "3d": 3, "4d": 4, "5d": 5})
    obs["logH_per_N_obs"] = obs["log_H"] / obs["n"]

    params = pd.read_csv("outputs_exploratory/prediction_a_xi_first_principles/scaling_parameters_first_principles.csv")
    fp_curve = pd.read_csv("outputs_exploratory/prediction_a_xi_first_principles/link_density_theory_vs_observed.csv")
    xi_obs = pd.read_csv("outputs_exploratory/prediction_a_large_n_scaling/xi_large_n.csv")

    merged = obs.merge(params[["d", "p_d"]], on="d", how="left").merge(
        fp_curve[["n", "d", "C0_per_N_theory"]],
        on=["n", "d"],
        how="left",
    )
    merged["logN"] = np.log(merged["n"])
    merged["h0"] = (1.0 - merged["p_d"]) * (merged["logN"] - 1.0)
    obs_med = float(xi_obs[xi_obs["dim_pair"] == "4→5"]["Xi"].median())
    return merged, obs_med


def evaluate_xi45(merged: pd.DataFrame, correction: np.ndarray, obs_med: float) -> tuple[np.ndarray, float, np.ndarray]:
    X = np.c_[merged["h0"].to_numpy(), correction]
    beta = np.linalg.lstsq(X, merged["logH_per_N_obs"].to_numpy(), rcond=None)[0]
    pred = X @ beta
    tmp = merged.copy()
    tmp["h_pred"] = pred
    vals = []
    for n in sorted(tmp["n"].unique()):
        lo = tmp[(tmp["d"] == 4) & (tmp["n"] == n)].iloc[0]
        hi = tmp[(tmp["d"] == 5) & (tmp["n"] == n)].iloc[0]
        delta_s = 2.0 * abs(lo["C0_per_N_theory"] - hi["C0_per_N_theory"])
        delta_h = abs(hi["h_pred"] - lo["h_pred"])
        vals.append(delta_s / delta_h)
    vals = np.array(vals, dtype=float)
    med = float(np.median(vals))
    return beta, med, vals


def main() -> None:
    base = build_base_table()
    merged, obs_med = load_obs_frame()
    merged = merged.merge(base, on=["n", "d"], how="left")

    alpha_rows = []
    for alpha in np.round(np.linspace(0.3, 2.5, 221), 3):
        corr = 1.0 - np.exp(-alpha * merged["m1"].to_numpy())
        beta, med, _ = evaluate_xi45(merged, corr, obs_med)
        alpha_rows.append(
            {
                "alpha": float(alpha),
                "xi45_median": med,
                "rel_err": abs(med - obs_med) / obs_med,
                "beta_h0": float(beta[0]),
                "beta_occ": float(beta[1]),
            }
        )
    alpha_df = pd.DataFrame(alpha_rows).sort_values(["rel_err", "alpha"])

    closures = []
    for label, col in [
        ("P_occ", "P_occ"),
        ("1-ell_MC", "one_minus_ell_mc"),
        ("P_gamma", "P_gamma"),
        ("1-ell_(2)", "one_minus_ell_2_clip"),
    ]:
        beta, med, vals = evaluate_xi45(merged, merged[col].to_numpy(), obs_med)
        closures.append(
            {
                "closure": label,
                "xi45_median": med,
                "rel_err": abs(med - obs_med) / obs_med,
                "beta_h0": float(beta[0]),
                "beta_corr": float(beta[1]),
                "xi45_values": ", ".join(f"{v:.2f}" for v in vals),
            }
        )
    closure_df = pd.DataFrame(closures).sort_values("rel_err")

    alpha_df.to_csv(OUT_DIR / "alpha_scan.csv", index=False)
    closure_df.to_csv(OUT_DIR / "closure_comparison.csv", index=False)

    print("Top alpha values:")
    print(alpha_df.head(12).to_string(index=False))
    print("\nClosure comparison:")
    print(closure_df.to_string(index=False))
    print(f"\nOutputs saved to: {OUT_DIR}")


if __name__ == "__main__":
    main()
