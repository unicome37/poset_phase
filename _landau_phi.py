"""Landau order parameter analysis for the unified structural functional.

Goal: Identify a macroscopic order parameter φ such that F can be approximately
written as F_macro[φ] = a·φ² + b·φ⁴ (+ higher order), where the Lorentzian
window corresponds to a specific φ range (potential minimum).

Candidates for φ:
  1. comp_frac (comparable fraction) — directly measures "Lorentzian-ness"
  2. d_eff (effective dimension) — measures dimensional consistency
  3. n_layers / N (normalized layer depth) — measures causal depth

For each candidate:
  - Plot F vs φ scatter (colored by family)
  - Fit polynomial F(φ) = a₀ + a₂(φ-φ₀)² + a₄(φ-φ₀)⁴
  - Check if minimum falls in Lorentzian window
  - Assess R² and residual structure
"""
import csv
import numpy as np
from numpy.polynomial import polynomial as P

# ── Load data ──
data = []
with open("outputs_unified_functional/raw_features.csv", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        data.append(row)

families = np.array([d["family"] for d in data])
N_arr = np.array([int(d["N"]) for d in data])

# Five functional terms
log_H = np.array([float(d["log_H"]) for d in data])
pi_geo = np.array([float(d["pi_geo"]) for d in data])
sigma_hist = np.array([float(d["sigma_hist"]) for d in data])
xi_dim = np.array([float(d["xi_dim"]) for d in data])
pi_cg = np.array([float(d["pi_cg"]) for d in data])

# Order parameter candidates
comp_frac = np.array([float(d["comp_frac"]) for d in data])
d_eff = np.array([float(d["d_eff"]) for d in data])
n_layers = np.array([int(d["n_layers"]) for d in data])
layer_density = n_layers / N_arr  # normalized

# Compute F with Bayesian posterior weights (β≡1)
X_all = np.column_stack([log_H, pi_geo, sigma_hist, xi_dim, pi_cg])
w_bayes = np.array([1.0, 0.0004, -0.888, 0.637, 0.068])
F_bayes = X_all @ w_bayes

# Also calibrated weights for comparison
w_calib = np.array([1.0, 0.25, -0.75, 0.05, 0.025])
F_calib = X_all @ w_calib

# ══════════════════════════════════════════════════════════════
# Analysis for each candidate φ
# ══════════════════════════════════════════════════════════════

candidates = {
    "comp_frac": comp_frac,
    "d_eff": d_eff, 
    "layer_density": layer_density,
}

# Lorentzian window
lor_mask = (comp_frac >= 0.30) & (comp_frac <= 0.55)

family_list = sorted(set(families))

def fit_landau(phi, F, label, degree=4):
    """Fit F(φ) = Σ aₖ φᵏ and report."""
    # Center φ
    phi_c = np.median(phi)
    x = phi - phi_c
    
    # Fit polynomial up to degree
    coeffs = np.polyfit(x, F, degree)  # highest degree first
    F_pred = np.polyval(coeffs, x)
    
    ss_res = np.sum((F - F_pred)**2)
    ss_tot = np.sum((F - F.mean())**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    
    # Find minimum
    # For degree-4 polynomial, find critical points
    deriv_coeffs = np.polyder(coeffs)
    roots = np.roots(deriv_coeffs)
    real_roots = roots[np.isreal(roots)].real
    
    # Evaluate at real roots and endpoints
    phi_range = np.linspace(x.min(), x.max(), 1000)
    F_curve = np.polyval(coeffs, phi_range)
    min_idx = np.argmin(F_curve)
    phi_min = phi_range[min_idx] + phi_c
    F_min = F_curve[min_idx]
    
    return coeffs, phi_c, r2, phi_min, F_min, F_pred


print("="*70)
print("LANDAU ORDER PARAMETER ANALYSIS")
print("="*70)

for name, phi in candidates.items():
    print(f"\n{'─'*70}")
    print(f"Candidate: φ = {name}")
    print(f"{'─'*70}")
    
    # Basic stats
    print(f"  Range: [{phi.min():.3f}, {phi.max():.3f}]")
    print(f"  In Lorentzian window: mean={phi[lor_mask].mean():.3f}, std={phi[lor_mask].std():.3f}")
    
    # Per-family distribution
    print(f"\n  Per-family φ means:")
    for fam in family_list:
        mask = families == fam
        print(f"    {fam:>10s}: φ={phi[mask].mean():.3f} ± {phi[mask].std():.3f}  F={F_bayes[mask].mean():.2f}")
    
    # Correlation with F
    from scipy.stats import spearmanr, pearsonr
    rho_s, p_s = spearmanr(phi, F_bayes)
    rho_p, p_p = pearsonr(phi, F_bayes)
    print(f"\n  Correlation with F (Bayes):")
    print(f"    Pearson r = {rho_p:.3f} (p={p_p:.2e})")
    print(f"    Spearman ρ = {rho_s:.3f} (p={p_s:.2e})")
    
    # Fit Landau polynomial (degree 2 and 4)
    for deg in [2, 4]:
        coeffs, phi_c, r2, phi_min, F_min, F_pred = fit_landau(phi, F_bayes, name, degree=deg)
        
        in_window = 0.30 <= phi_min <= 0.55 if name == "comp_frac" else True
        
        print(f"\n  Degree-{deg} fit (centered at φ₀={phi_c:.3f}):")
        print(f"    R² = {r2:.4f}")
        print(f"    Minimum at φ* = {phi_min:.3f}, F* = {F_min:.2f}")
        if name == "comp_frac":
            print(f"    φ* in Lorentzian window [0.30, 0.55]? {'YES' if in_window else 'NO'}")
        
        # Print coefficients
        print(f"    Coefficients (highest to lowest):")
        for i, c in enumerate(coeffs):
            power = deg - i
            print(f"      a_{power} = {c:+.4f}")
    
    # Within-N analysis: does φ separate families?
    print(f"\n  Within-N family separation by φ:")
    for n_val in sorted(set(N_arr)):
        mask_n = N_arr == n_val
        fams_at_n = sorted(set(families[mask_n]))
        
        # Check if Lor2D has lower F than KR within this N
        lor2d_mask = mask_n & (families == "Lor2D")
        kr_mask = mask_n & (families == "KR_like")
        
        if lor2d_mask.sum() > 0 and kr_mask.sum() > 0:
            phi_lor = phi[lor2d_mask].mean()
            phi_kr = phi[kr_mask].mean()
            F_lor = F_bayes[lor2d_mask].mean()
            F_kr = F_bayes[kr_mask].mean()
            sep = "✓" if F_lor < F_kr else "✗"
            print(f"    N={n_val:2d}: Lor2D φ={phi_lor:.3f} F={F_lor:.1f} | KR φ={phi_kr:.3f} F={F_kr:.1f} {sep}")


# ══════════════════════════════════════════════════════════════
# Best candidate analysis: deeper Landau fit
# ══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("BEST CANDIDATE DEEP ANALYSIS")
print(f"{'='*70}")

# comp_frac is the most natural Lorentzian-ness measure
# But let's also try PC1 score as a composite order parameter
from numpy.linalg import eigh

X_std = (X_all - X_all.mean(axis=0)) / X_all.std(axis=0)
C = np.cov(X_std.T)
eigvals, eigvecs = eigh(C)
idx = np.argsort(eigvals)[::-1]
eigvals = eigvals[idx]
eigvecs = eigvecs[:, idx]

pc1_score = X_std @ eigvecs[:, 0]

# Also try: composite "Lorentzian-ness" = weighted combination
# Normalize comp_frac to [0,1] range and combine with d_eff proximity to 4
lor_ness = comp_frac  # simplest: just comp_frac itself

candidates_extended = {
    "comp_frac": comp_frac,
    "PC1_score": pc1_score,
}

for name, phi in candidates_extended.items():
    print(f"\n  === {name} as order parameter ===")
    
    # Fit degree-4 Landau for ALL data
    coeffs, phi_c, r2, phi_min, F_min, F_pred = fit_landau(phi, F_bayes, name, degree=4)
    print(f"  All data: R²={r2:.4f}, φ*={phi_min:.3f}")
    
    # Within-N fits (N-dependent Landau coefficients)
    print(f"\n  N-dependent Landau fit F(φ) = a₂(φ-φ₀)² + a₄(φ-φ₀)⁴:")
    for n_val in sorted(set(N_arr)):
        mask_n = N_arr == n_val
        if mask_n.sum() < 6:
            continue
        phi_n = phi[mask_n]
        F_n = F_bayes[mask_n]
        
        coeffs_n, phi_c_n, r2_n, phi_min_n, F_min_n, _ = fit_landau(phi_n, F_n, name, degree=4)
        
        # Extract effective a₂ coefficient (controls curvature at minimum)
        # For centered polynomial a₄x⁴ + a₃x³ + a₂x² + a₁x + a₀
        a4 = coeffs_n[0]
        a2 = coeffs_n[2] if len(coeffs_n) > 2 else 0
        
        print(f"    N={n_val:2d}: R²={r2_n:.3f}, φ*={phi_min_n:.3f}, a₂={a2:+.2f}, a₄={a4:+.4f}")

# ══════════════════════════════════════════════════════════════
# Phase diagram: F vs comp_frac with Lorentzian window marked
# ══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("PHASE STRUCTURE IN φ = comp_frac SPACE")
print(f"{'='*70}")

# Bin comp_frac and compute mean F
bins = np.linspace(0, 1, 21)
bin_centers = (bins[:-1] + bins[1:]) / 2

print(f"\n  {'bin_center':>10s} {'count':>6s} {'mean_F':>8s} {'std_F':>8s} {'phase':>12s}")
print(f"  {'-'*10} {'-'*6} {'-'*8} {'-'*8} {'-'*12}")

for i in range(len(bins)-1):
    mask = (comp_frac >= bins[i]) & (comp_frac < bins[i+1])
    if mask.sum() == 0:
        continue
    mf = F_bayes[mask].mean()
    sf = F_bayes[mask].std()
    
    # Phase classification
    if bin_centers[i] < 0.15:
        phase = "random"
    elif bin_centers[i] < 0.30:
        phase = "transition"
    elif bin_centers[i] <= 0.55:
        phase = "Lorentzian"
    elif bin_centers[i] <= 0.70:
        phase = "transition"
    else:
        phase = "over-ordered"
    
    print(f"  {bin_centers[i]:10.2f} {mask.sum():6d} {mf:8.2f} {sf:8.2f} {phase:>12s}")


# ══════════════════════════════════════════════════════════════
# Key test: does the Landau minimum actually fall in the window?
# ══════════════════════════════════════════════════════════════

print(f"\n{'='*70}")
print("LANDAU MINIMUM LOCATION TEST")
print(f"{'='*70}")

phi = comp_frac
coeffs, phi_c, r2, phi_min, F_min, F_pred = fit_landau(phi, F_bayes, "comp_frac", degree=4)

print(f"\n  Degree-4 Landau fit of F(comp_frac):")
print(f"    R² = {r2:.4f}")
print(f"    Minimum at comp_frac* = {phi_min:.3f}")
print(f"    Lorentzian window = [0.30, 0.55]")
print(f"    Minimum IN window? {'YES ✓' if 0.30 <= phi_min <= 0.55 else 'NO ✗'}")

# Also check: does the Lor4D family sit near the minimum?
lor4d_mask = families == "Lor4D"
if lor4d_mask.sum() > 0:
    print(f"\n  Lor4D:  comp_frac = {phi[lor4d_mask].mean():.3f}, F = {F_bayes[lor4d_mask].mean():.2f}")
lor2d_mask = families == "Lor2D"
if lor2d_mask.sum() > 0:
    print(f"  Lor2D:  comp_frac = {phi[lor2d_mask].mean():.3f}, F = {F_bayes[lor2d_mask].mean():.2f}")
kr_mask = families == "KR_like"
if kr_mask.sum() > 0:
    print(f"  KR:     comp_frac = {phi[kr_mask].mean():.3f}, F = {F_bayes[kr_mask].mean():.2f}")

# Degree-6 for robustness
coeffs6, phi_c6, r26, phi_min6, F_min6, _ = fit_landau(phi, F_bayes, "comp_frac", degree=6)
print(f"\n  Degree-6 fit: R²={r26:.4f}, minimum at {phi_min6:.3f}")

# Per-N Landau minima
print(f"\n  Per-N Landau minima (degree-4):")
for n_val in sorted(set(N_arr)):
    mask_n = N_arr == n_val
    if mask_n.sum() < 8:
        continue
    phi_n = phi[mask_n]
    F_n = F_bayes[mask_n]
    _, _, r2_n, phi_min_n, _, _ = fit_landau(phi_n, F_n, "comp_frac", degree=4)
    in_win = "✓" if 0.30 <= phi_min_n <= 0.55 else "✗"
    print(f"    N={n_val:2d} ({mask_n.sum():2d} pts): φ*={phi_min_n:.3f} R²={r2_n:.3f} {in_win}")
