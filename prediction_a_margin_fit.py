"""
Prediction A — Margin of Victory Statistical Fit
OLS regression: margin(N) = a + b*N for consistency action, 4D vs {2D, 3D}.
Reports: slope ± SE, R², p-value.
"""
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

OUT_DIR = Path("outputs_exploratory/prediction_a_margin_fit")
OUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv("outputs_exploratory/prediction_a_margin_summary/prediction_a_margin_summary.csv")

# Filter consistency action only
df_con = df[df["variant"] == "A2_replace_dim_with_consistency"].copy()

results = []

for opponent in ["lorentzian_like_2d", "lorentzian_like_3d"]:
    sub = df_con[df_con["right_family"] == opponent].sort_values("n")
    N = sub["n"].values
    margin = sub["mean_margin"].values

    # OLS: margin = a + b * N
    slope, intercept, r_value, p_value, std_err = stats.linregress(N, margin)

    print(f"\n=== 4D vs {opponent} (consistency action) ===")
    print(f"  N range: {N.min()}–{N.max()} ({len(N)} points)")
    print(f"  Slope:     {slope:.4f} ± {std_err:.4f} per unit N")
    print(f"  Intercept: {intercept:.4f}")
    print(f"  R²:        {r_value**2:.6f}")
    print(f"  p-value:   {p_value:.2e}")

    # Also try power-law fit: margin = c * N^d
    log_N = np.log(N)
    log_m = np.log(margin)
    valid = np.isfinite(log_m)
    if valid.sum() >= 3:
        sl_pw, int_pw, r_pw, p_pw, se_pw = stats.linregress(log_N[valid], log_m[valid])
        print(f"\n  Power-law fit: margin ∝ N^{sl_pw:.3f}")
        print(f"  Power-law R²: {r_pw**2:.6f}")
    else:
        sl_pw, r_pw, p_pw = np.nan, np.nan, np.nan

    results.append({
        "opponent": opponent,
        "n_points": len(N),
        "n_min": int(N.min()),
        "n_max": int(N.max()),
        "linear_slope": slope,
        "linear_slope_se": std_err,
        "linear_intercept": intercept,
        "linear_R2": r_value**2,
        "linear_p_value": p_value,
        "powerlaw_exponent": sl_pw,
        "powerlaw_R2": r_pw**2 if not np.isnan(r_pw) else np.nan,
    })

rdf = pd.DataFrame(results)
rdf.to_csv(OUT_DIR / "margin_ols_fit.csv", index=False)
print(f"\nSaved to {OUT_DIR / 'margin_ols_fit.csv'}")

# --- Visualization ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, opponent in zip(axes, ["lorentzian_like_2d", "lorentzian_like_3d"]):
    sub = df_con[df_con["right_family"] == opponent].sort_values("n")
    N = sub["n"].values
    margin = sub["mean_margin"].values
    margin_min = sub["min_margin"].values
    margin_max = sub["max_margin"].values

    slope, intercept, r_value, p_value, std_err = stats.linregress(N, margin)

    ax.errorbar(N, margin,
                yerr=[margin - margin_min, margin_max - margin],
                fmt='o-', capsize=3, label="Mean margin ± range")
    N_fit = np.linspace(N.min(), N.max(), 100)
    ax.plot(N_fit, intercept + slope * N_fit, 'r--',
            label=f"OLS: {slope:.3f}·N + {intercept:.2f}\nR²={r_value**2:.4f}, p={p_value:.1e}")
    ax.set_xlabel("N")
    ax.set_ylabel("Margin (4D − opponent)")
    short = "2D" if "2d" in opponent else "3D"
    ax.set_title(f"4D vs {short} (consistency)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

plt.tight_layout()
fig.savefig(OUT_DIR / "margin_ols_fit.png", dpi=150)
print(f"Saved figure to {OUT_DIR / 'margin_ols_fit.png'}")
