"""
Prediction B — γ_c(N) Finite-Size Scaling Fit

Fits:
  1) Constant: γ_c = a              (null hypothesis: truly bounded)
  2) Power-law decay: γ_c = a + b / N^c
  3) Linear: γ_c = a + b * N       (alternative: drift)

Reports asymptotic limit, AIC comparison, 95% CI.
Uses A2-action lor2d vs KR_like data from both frozen_exact and medium_exact_scan.
"""
import numpy as np
import pandas as pd
from scipy import stats, optimize
from pathlib import Path

OUT_DIR = Path("outputs_exploratory/prediction_b_gamma_c_scaling")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Load γ_c data from both confirmatory sources
df1 = pd.read_csv("outputs_confirmatory/frozen_exact/gamma_c_report.csv")
df2 = pd.read_csv("outputs_confirmatory/medium_exact_scan/gamma_c_report.csv")
df = pd.concat([df1, df2], ignore_index=True)

# Filter: A2, lor2d vs KR_like
mask = (df["action_mode"] == "A2") & \
       (df["family_a"] == "lorentzian_like_2d") & \
       (df["family_b"] == "KR_like")
sub = df[mask].drop_duplicates(subset="n").sort_values("n")

N = sub["n"].values
gc = sub["gamma_c_est"].values

print(f"γ_c data points: {len(N)}")
print(f"N range: {N.min()} – {N.max()}")
print(f"γ_c range: {gc.min():.4f} – {gc.max():.4f}")
print()

# --- Exclude N=10 outlier (known small-N artifact, γ_c=9.19) ---
mask_valid = N >= 12
N_fit = N[mask_valid]
gc_fit = gc[mask_valid]

print(f"Fitting on N >= 12: {len(N_fit)} points")
print()

# Model 1: Constant
mean_gc = np.mean(gc_fit)
std_gc = np.std(gc_fit, ddof=1)
se_gc = std_gc / np.sqrt(len(gc_fit))
ci_95 = stats.t.interval(0.95, df=len(gc_fit)-1, loc=mean_gc, scale=se_gc)
resid1 = gc_fit - mean_gc
ss1 = np.sum(resid1**2)
k1 = 1
aic1 = len(gc_fit) * np.log(ss1 / len(gc_fit)) + 2 * k1

print(f"Model 1: Constant")
print(f"  γ_c = {mean_gc:.4f} ± {se_gc:.4f}")
print(f"  95% CI: [{ci_95[0]:.4f}, {ci_95[1]:.4f}]")
print(f"  AIC = {aic1:.2f}")
print()

# Model 2: γ_c = a + b / N^c
def power_model(N, a, b, c):
    return a + b / np.power(N, c)

try:
    popt, pcov = optimize.curve_fit(power_model, N_fit, gc_fit,
                                     p0=[0.3, 5.0, 1.0],
                                     maxfev=10000)
    perr = np.sqrt(np.diag(pcov))
    resid2 = gc_fit - power_model(N_fit, *popt)
    ss2 = np.sum(resid2**2)
    k2 = 3
    aic2 = len(gc_fit) * np.log(ss2 / len(gc_fit)) + 2 * k2

    print(f"Model 2: γ_c = a + b / N^c")
    print(f"  a (asymptotic) = {popt[0]:.4f} ± {perr[0]:.4f}")
    print(f"  b = {popt[1]:.4f} ± {perr[1]:.4f}")
    print(f"  c = {popt[2]:.4f} ± {perr[2]:.4f}")
    print(f"  Asymptotic γ_c(∞) = {popt[0]:.4f}")
    print(f"  AIC = {aic2:.2f}")
    power_fit_ok = True
except Exception as e:
    print(f"Model 2 fit failed: {e}")
    power_fit_ok = False
print()

# Model 3: Linear γ_c = a + b * N
slope, intercept, r_value, p_value, std_err = stats.linregress(N_fit, gc_fit)
resid3 = gc_fit - (intercept + slope * N_fit)
ss3 = np.sum(resid3**2)
k3 = 2
aic3 = len(gc_fit) * np.log(ss3 / len(gc_fit)) + 2 * k3

print(f"Model 3: Linear γ_c = a + b·N")
print(f"  slope = {slope:.6f} ± {std_err:.6f}")
print(f"  intercept = {intercept:.4f}")
print(f"  R² = {r_value**2:.6f}")
print(f"  p(slope≠0) = {p_value:.4f}")
print(f"  AIC = {aic3:.2f}")
print()

# --- Summary table ---
rows = [
    {"model": "constant", "params": 1, "AIC": aic1,
     "asymptotic_gamma_c": mean_gc, "asymptotic_95ci_lo": ci_95[0],
     "asymptotic_95ci_hi": ci_95[1]},
    {"model": "linear", "params": 2, "AIC": aic3,
     "slope": slope, "slope_se": std_err, "slope_p": p_value},
]
if power_fit_ok:
    rows.insert(1, {"model": "power_decay", "params": 3, "AIC": aic2,
                     "asymptotic_gamma_c": popt[0], "decay_exponent": popt[2]})

rdf = pd.DataFrame(rows)
rdf.to_csv(OUT_DIR / "gamma_c_scaling_fit.csv", index=False)

# Best model
aics = {"constant": aic1, "linear": aic3}
if power_fit_ok:
    aics["power_decay"] = aic2
best = min(aics, key=aics.get)
print(f"** Best model by AIC: {best} (AIC={aics[best]:.2f}) **")
print()

# --- Visualization ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(N, gc, 'ko', markersize=8, label="Data (all N)")
ax.plot(N_fit, gc_fit, 'bs', markersize=6, label="Fit data (N≥12)")

# Constant
ax.axhline(mean_gc, color='green', linestyle='--',
           label=f"Constant: γ_c={mean_gc:.3f} [{ci_95[0]:.3f}, {ci_95[1]:.3f}]")
ax.axhspan(ci_95[0], ci_95[1], alpha=0.1, color='green')

# Power decay
if power_fit_ok:
    N_dense = np.linspace(12, 100, 200)
    ax.plot(N_dense, power_model(N_dense, *popt), 'r-',
            label=f"Power: γ_c(∞)={popt[0]:.3f}, exponent={popt[2]:.2f}")

# Linear
N_dense = np.linspace(10, 100, 200)
ax.plot(N_dense, intercept + slope * N_dense, 'purple', linestyle=':',
        label=f"Linear: slope={slope:.4f}, p={p_value:.3f}")

ax.set_xlabel("N", fontsize=12)
ax.set_ylabel("γ_c", fontsize=12)
ax.set_title("Prediction B: γ_c(N) Finite-Size Scaling (A2, lor2d vs KR)")
ax.legend(fontsize=8)
ax.set_xlim(8, 52)
ax.set_ylim(-0.5, 2.0)
ax.grid(alpha=0.3)

plt.tight_layout()
fig.savefig(OUT_DIR / "gamma_c_scaling_fit.png", dpi=150)
print(f"Saved to {OUT_DIR}")
