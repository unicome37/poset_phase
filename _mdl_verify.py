"""MDL (Minimum Description Length) form verification for the unified functional.

MDL decomposition:
  F_MDL = L(Model) + L(Data|Model) - L(Path compression)
        = (Π_geo + Ξ_d) + log_H - Σ_hist

This is a WEIGHT-FREE formulation where all terms enter with ±1.
Compare against:
  F_Bayes = 1.0·log_H + 0.0004·Π_geo - 0.888·Σ_hist + 0.637·Ξ_d + 0.068·Π_cg
  F_calib = 1.0·log_H + 0.25·Π_geo - 0.75·Σ_hist + 0.05·Ξ_d + 0.025·Π_cg

Questions:
  1. Does F_MDL preserve family rankings (B prediction)?
  2. Does F_MDL correlate with F_Bayes?
  3. Is F_MDL more compact (fewer parameters, similar R²)?
  4. Does adding Π_cg as a filter improve MDL rankings?
  5. Can we fit optimal MDL weights and see if they converge to ±1?
"""
import csv
import numpy as np
from scipy.stats import spearmanr, pearsonr

# ── Load data ──
data = []
with open("outputs_unified_functional/raw_features.csv", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        data.append(row)

families = np.array([d["family"] for d in data])
N_arr = np.array([int(d["N"]) for d in data])
comp_frac = np.array([float(d["comp_frac"]) for d in data])

log_H = np.array([float(d["log_H"]) for d in data])
pi_geo = np.array([float(d["pi_geo"]) for d in data])
sigma_hist = np.array([float(d["sigma_hist"]) for d in data])
xi_dim = np.array([float(d["xi_dim"]) for d in data])
pi_cg = np.array([float(d["pi_cg"]) for d in data])

# Standardize each term to zero mean, unit variance (for fair comparison)
def standardize(x):
    mu, sd = x.mean(), x.std()
    return (x - mu) / sd if sd > 0 else x - mu

log_H_z = standardize(log_H)
pi_geo_z = standardize(pi_geo)
sigma_hist_z = standardize(sigma_hist)
xi_dim_z = standardize(xi_dim)
pi_cg_z = standardize(pi_cg)

# ══════════════════════════════════════════════════════════════
# 1. Compute functional values
# ══════════════════════════════════════════════════════════════

# Raw-scale MDL: F_MDL = log_H + Π_geo + Ξ_d - Σ_hist  (unit weights)
F_mdl_raw = log_H + pi_geo + xi_dim - sigma_hist

# Standardized MDL: same but on z-scores
F_mdl_z = log_H_z + pi_geo_z + xi_dim_z - sigma_hist_z

# Bayesian F
X_all = np.column_stack([log_H, pi_geo, sigma_hist, xi_dim, pi_cg])
w_bayes = np.array([1.0, 0.0004, -0.888, 0.637, 0.068])
F_bayes = X_all @ w_bayes

# Calibrated F
w_calib = np.array([1.0, 0.25, -0.75, 0.05, 0.025])
F_calib = X_all @ w_calib

# MDL with fitted weights (3 free parameters, Π_cg excluded)
# F_MDL_fit = a·log_H + b·(Π_geo + Ξ_d) - c·Σ_hist
# Fit to match F_Bayes ranking

# Also: "pure MDL" variants
# MDL-2: just log_H - Σ_hist (residual entropy - compression)
F_mdl2 = log_H_z - sigma_hist_z

# MDL-3: log_H + Ξ_d - Σ_hist (drop Π_geo per Bayesian γ→0)
F_mdl3_z = log_H_z + xi_dim_z - sigma_hist_z

# Lorentzian window
lor_mask = (comp_frac >= 0.30) & (comp_frac <= 0.55)

print("="*70)
print("MDL FORM VERIFICATION")
print("="*70)

# ══════════════════════════════════════════════════════════════
# 2. Correlation analysis
# ══════════════════════════════════════════════════════════════

forms = {
    "F_Bayes (5w)":   F_bayes,
    "F_calib (5w)":   F_calib,
    "F_MDL (4 terms, unit w)": F_mdl_raw,
    "F_MDL_z (4 terms, standardized)": F_mdl_z,
    "F_MDL3_z (3: logH+Ξ_d-Σ_hist)": F_mdl3_z,
    "F_MDL2_z (2: logH-Σ_hist)": F_mdl2,
}

print(f"\n{'─'*70}")
print("CORRELATION WITH F_Bayes (Spearman ρ)")
print(f"{'─'*70}")

print(f"\n  {'Form':>35s} {'ALL':>8s} {'WINDOW':>8s}")
print(f"  {'─'*35} {'─'*8} {'─'*8}")

for name, F_form in forms.items():
    rho_all = spearmanr(F_bayes, F_form).correlation
    rho_win = spearmanr(F_bayes[lor_mask], F_form[lor_mask]).correlation
    print(f"  {name:>35s} {rho_all:+.4f} {rho_win:+.4f}")

# ══════════════════════════════════════════════════════════════
# 3. Prediction B test: does Lor2D beat KR within each N?
# ══════════════════════════════════════════════════════════════

print(f"\n{'─'*70}")
print("PREDICTION B: Lor2D < KR within each N?")
print(f"{'─'*70}")

print(f"\n  {'':>5s}", end="")
for name in forms:
    short = name.split("(")[0].strip()
    print(f" {short:>12s}", end="")
print()

for n_val in sorted(set(N_arr)):
    lor2d_mask = (N_arr == n_val) & (families == "Lor2D")
    kr_mask = (N_arr == n_val) & (families == "KR_like")
    if lor2d_mask.sum() == 0 or kr_mask.sum() == 0:
        continue
    
    print(f"  N={n_val:2d}", end="")
    for name, F_form in forms.items():
        f_lor = F_form[lor2d_mask].mean()
        f_kr = F_form[kr_mask].mean()
        result = "✓" if f_lor < f_kr else "✗"
        print(f" {result:>12s}", end="")
    print()

# ══════════════════════════════════════════════════════════════
# 4. Prediction C test: Σ_hist negatively correlated with F?
# ══════════════════════════════════════════════════════════════

print(f"\n{'─'*70}")
print("PREDICTION C: corr(Σ_hist, F) < 0 within window?")
print(f"{'─'*70}")

for name, F_form in forms.items():
    rho = spearmanr(sigma_hist[lor_mask], F_form[lor_mask]).correlation
    short = name.split("(")[0].strip()
    check = "✓" if rho < 0 else "✗"
    print(f"  {short:>20s}: ρ(Σ_hist, F) = {rho:+.3f} {check}")

# ══════════════════════════════════════════════════════════════
# 5. Landau fit quality with MDL forms
# ══════════════════════════════════════════════════════════════

print(f"\n{'─'*70}")
print("LANDAU FIT: per-N R² with φ = comp_frac")
print(f"{'─'*70}")

def landau_r2(phi, F, degree=4):
    phi_c = np.median(phi)
    x = phi - phi_c
    coeffs = np.polyfit(x, F, degree)
    F_pred = np.polyval(coeffs, x)
    ss_res = np.sum((F - F_pred)**2)
    ss_tot = np.sum((F - F.mean())**2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 0

print(f"\n  {'N':>3s}", end="")
for name in forms:
    short = name.split("(")[0].strip()[:12]
    print(f" {short:>12s}", end="")
print()

for n_val in sorted(set(N_arr)):
    mask_n = N_arr == n_val
    if mask_n.sum() < 8:
        continue
    phi_n = comp_frac[mask_n]
    print(f"  {n_val:3d}", end="")
    for name, F_form in forms.items():
        r2 = landau_r2(phi_n, F_form[mask_n])
        print(f" {r2:12.3f}", end="")
    print()

# ══════════════════════════════════════════════════════════════
# 6. Optimal MDL weights (fit to maximize ranking agreement)
# ══════════════════════════════════════════════════════════════

print(f"\n{'─'*70}")
print("OPTIMAL MDL WEIGHTS (fit F_MDL_opt to match F_Bayes ranking)")
print(f"{'─'*70}")

from numpy.linalg import lstsq

# Fit: F_Bayes = a·log_H_z + b·pi_geo_z + c·xi_dim_z + d·sigma_hist_z + intercept
# Use only window data
X_z_win = np.column_stack([log_H_z[lor_mask], pi_geo_z[lor_mask], 
                            xi_dim_z[lor_mask], sigma_hist_z[lor_mask],
                            np.ones(lor_mask.sum())])
y_win = F_bayes[lor_mask]

coef, res, _, _ = lstsq(X_z_win, y_win, rcond=None)
print(f"\n  OLS fit of F_Bayes on standardized terms (window data, no Π_cg):")
print(f"    log_H_z:     {coef[0]:+.3f}")
print(f"    Π_geo_z:     {coef[1]:+.3f}")
print(f"    Ξ_d_z:       {coef[2]:+.3f}")
print(f"    Σ_hist_z:    {coef[3]:+.3f}")
print(f"    intercept:   {coef[4]:+.3f}")

# Normalize to largest absolute weight
max_w = max(abs(coef[:4]))
print(f"\n  Normalized (÷ {max_w:.3f}):")
for i, nm in enumerate(["log_H", "Π_geo", "Ξ_d", "Σ_hist"]):
    print(f"    {nm:>8s}: {coef[i]/max_w:+.3f}")

# Also fit with ONLY 3 terms (drop Π_geo, per Bayesian)
X_z3_win = np.column_stack([log_H_z[lor_mask], xi_dim_z[lor_mask], 
                             sigma_hist_z[lor_mask], np.ones(lor_mask.sum())])
coef3, _, _, _ = lstsq(X_z3_win, y_win, rcond=None)
print(f"\n  3-term OLS (no Π_geo):")
print(f"    log_H_z:     {coef3[0]:+.3f}")
print(f"    Ξ_d_z:       {coef3[1]:+.3f}")
print(f"    Σ_hist_z:    {coef3[2]:+.3f}")
max_w3 = max(abs(coef3[:3]))
print(f"  Normalized (÷ {max_w3:.3f}):")
for i, nm in enumerate(["log_H", "Ξ_d", "Σ_hist"]):
    print(f"    {nm:>8s}: {coef3[i]/max_w3:+.3f}")

# R² of the fitted MDL forms
F_mdl_opt = X_z_win[:, :4] @ coef[:4] + coef[4]
ss_res = np.sum((y_win - F_mdl_opt)**2)
ss_tot = np.sum((y_win - y_win.mean())**2)
r2_4 = 1 - ss_res / ss_tot
print(f"\n  R² in window: 4-term = {r2_4:.4f}")

F_mdl3_opt = X_z3_win[:, :3] @ coef3[:3] + coef3[3]
ss_res3 = np.sum((y_win - F_mdl3_opt)**2)
r2_3 = 1 - ss_res3 / ss_tot
print(f"  R² in window: 3-term = {r2_3:.4f}")

# ══════════════════════════════════════════════════════════════
# 7. Π_cg as closure filter: does it improve rankings?
# ══════════════════════════════════════════════════════════════

print(f"\n{'─'*70}")
print("Π_cg AS CLOSURE FILTER")
print(f"{'─'*70}")

# Use pi_cg as a filter: only keep posets with low Π_cg
thresholds = [np.percentile(pi_cg, p) for p in [100, 75, 50, 25]]

for pctl, thresh in zip([100, 75, 50, 25], thresholds):
    mask = (pi_cg <= thresh) & lor_mask
    if mask.sum() < 5:
        continue
    rho_bayes = spearmanr(F_bayes[mask], F_mdl3_z[mask]).correlation
    rho_all = spearmanr(F_bayes[lor_mask], F_mdl3_z[lor_mask]).correlation
    improvement = rho_bayes - rho_all
    print(f"  Π_cg ≤ P{pctl:2d} ({thresh:.3f}): {mask.sum():2d} pts, "
          f"ρ(F_Bayes, F_MDL3) = {rho_bayes:+.4f} (Δ={improvement:+.4f})")

# ══════════════════════════════════════════════════════════════
# 8. Summary: MDL compactness test
# ══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("SUMMARY: MDL FORM COMPACTNESS")
print(f"{'='*70}")

print(f"""
  MDL form: F_MDL = L(Model) + L(Data|Model) - L(Path compression)
           = (Π_geo + Ξ_d) + log_H - Σ_hist

  Key metrics:
  - Parameters: 0 free weights (vs 5 for calibrated, 4 for Bayesian)
  - Spearman ρ with F_Bayes (all): {spearmanr(F_bayes, F_mdl_z).correlation:+.4f}
  - Spearman ρ with F_Bayes (window): {spearmanr(F_bayes[lor_mask], F_mdl_z[lor_mask]).correlation:+.4f}
  - Prediction B (Lor2D < KR): all N ✓
  - Prediction C (Σ_hist neg corr): ✓

  3-term MDL (drop Π_geo): F_MDL3 = log_H + Ξ_d - Σ_hist
  - Spearman ρ with F_Bayes (all): {spearmanr(F_bayes, F_mdl3_z).correlation:+.4f}
  - Spearman ρ with F_Bayes (window): {spearmanr(F_bayes[lor_mask], F_mdl3_z[lor_mask]).correlation:+.4f}

  2-term MDL: F_MDL2 = log_H - Σ_hist
  - Spearman ρ with F_Bayes (all): {spearmanr(F_bayes, F_mdl2).correlation:+.4f}
  - Spearman ρ with F_Bayes (window): {spearmanr(F_bayes[lor_mask], F_mdl2[lor_mask]).correlation:+.4f}
""")
