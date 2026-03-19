"""PCA / rank-degeneracy analysis of five functional terms within the Lorentzian window.

Goal: Test the "renormalization conjecture" — do the 5 micro terms collapse to
2-3 effective directions inside the competitive window (comp_frac ∈ [0.3, 0.55])?

Outputs:
1. Correlation matrix of 5 terms (all data + window-only)
2. PCA eigenvalues and cumulative variance
3. Loadings: which terms cluster into "C_maint" vs "S_res" vs "R_hist"
4. Ranking preservation: can top-2 or top-3 PCs reproduce F ranking?
"""
import csv
import numpy as np
from itertools import combinations

# ── Load data ──
data = []
with open("outputs_unified_functional/raw_features.csv", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        data.append(row)

families = np.array([d["family"] for d in data])
N_arr = np.array([int(d["N"]) for d in data])
comp_frac = np.array([float(d["comp_frac"]) for d in data])

# Five functional terms
names = ["log_H", "Π_geo", "Σ_hist", "Ξ_d", "Π_cg"]
log_H = np.array([float(d["log_H"]) for d in data])
pi_geo = np.array([float(d["pi_geo"]) for d in data])
sigma_hist = np.array([float(d["sigma_hist"]) for d in data])
xi_dim = np.array([float(d["xi_dim"]) for d in data])
pi_cg = np.array([float(d["pi_cg"]) for d in data])

X_all = np.column_stack([log_H, pi_geo, sigma_hist, xi_dim, pi_cg])

# Lorentzian window mask
window_mask = (comp_frac >= 0.30) & (comp_frac <= 0.55)
X_win = X_all[window_mask]

print(f"Total posets: {len(data)}")
print(f"In Lorentzian window (comp_frac ∈ [0.30, 0.55]): {window_mask.sum()}")
print(f"  Families in window: ", end="")
for fam in sorted(set(families[window_mask])):
    count = np.sum((families == fam) & window_mask)
    print(f"{fam}({count}) ", end="")
print()

# ══════════════════════════════════════════════════════════════
# 1. Correlation matrices
# ══════════════════════════════════════════════════════════════
def print_corr(X, label, col_names):
    print(f"\n{'='*70}")
    print(f"CORRELATION MATRIX — {label}")
    print(f"{'='*70}")
    corr = np.corrcoef(X.T)
    print(f"  {'':>8s}", end="")
    for n in col_names:
        print(f" {n:>8s}", end="")
    print()
    for i, n in enumerate(col_names):
        print(f"  {n:>8s}", end="")
        for j in range(len(col_names)):
            print(f" {corr[i,j]:+8.3f}", end="")
        print()

print_corr(X_all, "ALL DATA (N=160)", names)
print_corr(X_win, f"LORENTZIAN WINDOW (N={window_mask.sum()})", names)


# ══════════════════════════════════════════════════════════════
# 2. PCA
# ══════════════════════════════════════════════════════════════
def do_pca(X, label, col_names):
    # Standardize
    mu = X.mean(axis=0)
    std = X.std(axis=0)
    std[std < 1e-12] = 1.0
    Z = (X - mu) / std
    
    # Covariance
    C = np.cov(Z.T)
    eigvals, eigvecs = np.linalg.eigh(C)
    
    # Sort descending
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    
    total = eigvals.sum()
    cumvar = np.cumsum(eigvals) / total
    
    print(f"\n{'='*70}")
    print(f"PCA — {label}")
    print(f"{'='*70}")
    print(f"\n  {'PC':>4s} {'Eigenvalue':>12s} {'Var%':>8s} {'CumVar%':>8s}")
    print(f"  {'-'*4} {'-'*12} {'-'*8} {'-'*8}")
    for k in range(len(eigvals)):
        pct = 100 * eigvals[k] / total
        cum = 100 * cumvar[k]
        print(f"  PC{k+1:d} {eigvals[k]:12.4f} {pct:7.1f}% {cum:7.1f}%")
    
    # Effective rank
    eff_rank = np.exp(-np.sum((eigvals/total) * np.log(eigvals/total + 1e-20)))
    print(f"\n  Effective rank (Shannon entropy): {eff_rank:.2f}")
    print(f"  Variance explained by PC1-2: {100*cumvar[1]:.1f}%")
    print(f"  Variance explained by PC1-3: {100*cumvar[2]:.1f}%")
    
    # Loadings
    print(f"\n  Loadings (eigenvectors):")
    print(f"  {'':>8s}", end="")
    for k in range(min(4, len(eigvals))):
        print(f"   {'PC'+str(k+1):>6s}", end="")
    print()
    for i, n in enumerate(col_names):
        print(f"  {n:>8s}", end="")
        for k in range(min(4, len(eigvals))):
            print(f"   {eigvecs[i,k]:+6.3f}", end="")
        print()
    
    return eigvals, eigvecs, mu, std, cumvar

eigvals_all, eigvecs_all, mu_all, std_all, _ = do_pca(X_all, "ALL DATA", names)
eigvals_win, eigvecs_win, mu_win, std_win, cumvar_win = do_pca(X_win, "LORENTZIAN WINDOW", names)


# ══════════════════════════════════════════════════════════════
# 3. Collinearity check: Π_geo, Ξ_d, Π_cg
# ══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("COLLINEARITY CHECK: Π_geo, Ξ_d, Π_cg (structure maintenance candidates)")
print(f"{'='*70}")

# Pairwise correlations in window
struct_idx = [1, 3, 4]  # pi_geo, xi_dim, pi_cg
struct_names = ["Π_geo", "Ξ_d", "Π_cg"]

for subset_label, X_sub in [("ALL", X_all), ("WINDOW", X_win)]:
    print(f"\n  {subset_label}:")
    corr = np.corrcoef(X_sub[:, struct_idx].T)
    for i in range(3):
        for j in range(i+1, 3):
            print(f"    corr({struct_names[i]}, {struct_names[j]}) = {corr[i,j]:+.3f}")
    
    # VIF-like: regress each on the other two
    from numpy.linalg import lstsq
    for i in range(3):
        y = X_sub[:, struct_idx[i]]
        others = np.column_stack([X_sub[:, struct_idx[j]] for j in range(3) if j != i])
        others = np.column_stack([others, np.ones(len(y))])
        coef, res, _, _ = lstsq(others, y, rcond=None)
        ss_res = np.sum((y - others @ coef)**2)
        ss_tot = np.sum((y - y.mean())**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        vif = 1 / (1 - r2) if r2 < 1 else float('inf')
        print(f"    {struct_names[i]} → R² from others = {r2:.3f}, VIF = {vif:.1f}")


# ══════════════════════════════════════════════════════════════
# 4. Σ_hist independence check
# ══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("Σ_hist INDEPENDENCE CHECK")
print(f"{'='*70}")

for subset_label, X_sub in [("ALL", X_all), ("WINDOW", X_win)]:
    # Regress Σ_hist on all other 4
    y = X_sub[:, 2]  # sigma_hist
    others = np.column_stack([X_sub[:, j] for j in [0, 1, 3, 4]])
    others = np.column_stack([others, np.ones(len(y))])
    coef, res, _, _ = np.linalg.lstsq(others, y, rcond=None)
    ss_res = np.sum((y - others @ coef)**2)
    ss_tot = np.sum((y - y.mean())**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    print(f"  {subset_label}: Σ_hist ~ (log_H + Π_geo + Ξ_d + Π_cg) → R² = {r2:.3f}")
    print(f"    Σ_hist retains {100*(1-r2):.1f}% unique variance")


# ══════════════════════════════════════════════════════════════
# 5. Ranking preservation with reduced PCs
# ══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("RANKING PRESERVATION TEST")
print(f"{'='*70}")

# Use Bayesian posterior weights for F
# w = [β=1, γ~0, -λ=-0.89, η=0.64, κ=0.07]
w_bayes = np.array([1.0, 0.0004, -0.888, 0.637, 0.068])
F_full = X_all @ w_bayes

# Also test with calibrated weights (β=2 scale, normalized to β=1)
w_calib = np.array([1.0, 0.25, -0.75, 0.05, 0.025])
F_calib = X_all @ w_calib

def rank_corr(a, b):
    """Spearman rank correlation."""
    from scipy.stats import spearmanr
    return spearmanr(a, b).correlation

# For window data
F_full_win = X_win @ w_bayes
F_calib_win = X_win @ w_calib

# Project onto top-k PCs and reconstruct F
Z_win = (X_win - mu_win) / std_win

for k in [1, 2, 3, 4]:
    # Project onto top-k PCs
    scores = Z_win @ eigvecs_win[:, :k]
    # Reconstruct in standardized space
    Z_recon = scores @ eigvecs_win[:, :k].T
    # Back to original scale
    X_recon = Z_recon * std_win + mu_win
    
    F_recon_bayes = X_recon @ w_bayes
    F_recon_calib = X_recon @ w_calib
    
    rho_bayes = rank_corr(F_full_win, F_recon_bayes)
    rho_calib = rank_corr(F_calib_win, F_recon_calib)
    
    print(f"  PC1-{k}: Spearman ρ(F_full, F_recon) = {rho_bayes:.4f} (Bayes w), {rho_calib:.4f} (calib w)")

# Also check: within each N, does the family ranking preserve?
print(f"\n  Within-N family ranking preservation (Bayes weights, window data):")
for n_val in sorted(set(N_arr)):
    mask = window_mask & (N_arr == n_val)
    if mask.sum() < 3:
        continue
    X_n = X_all[mask]
    F_n_full = X_n @ w_bayes
    
    Z_n = (X_n - mu_win) / std_win
    for k in [2, 3]:
        scores_n = Z_n @ eigvecs_win[:, :k]
        Z_n_recon = scores_n @ eigvecs_win[:, :k].T
        X_n_recon = Z_n_recon * std_win + mu_win
        F_n_recon = X_n_recon @ w_bayes
        
        rho = rank_corr(F_n_full, F_n_recon)
        print(f"    N={int(n_val):2d} ({mask.sum():2d} posets): PC1-{k} ρ = {rho:.3f}")


# ══════════════════════════════════════════════════════════════
# 6. Effective functional interpretation
# ══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print("EFFECTIVE FUNCTIONAL INTERPRETATION")
print(f"{'='*70}")

# In window PCA, interpret top 3 PCs
pc_labels_guess = []
for k in range(min(3, len(names))):
    loadings = eigvecs_win[:, k]
    # Find dominant terms
    dominant = [(names[i], loadings[i]) for i in np.argsort(np.abs(loadings))[::-1]]
    
    print(f"\n  PC{k+1} (var = {100*eigvals_win[k]/eigvals_win.sum():.1f}%):")
    for name, load in dominant:
        bar = "█" * int(abs(load) * 20)
        sign = "+" if load > 0 else "-"
        print(f"    {name:>8s}: {sign}{abs(load):.3f} {bar}")
    
    # Interpret
    top2_names = [dominant[0][0], dominant[1][0]]
    if "Σ_hist" in top2_names and "log_H" in top2_names:
        label = "S_res - R_hist (residual entropy vs history)"
    elif "Ξ_d" in top2_names or "Π_geo" in top2_names:
        label = "C_maint (structure maintenance cost)"
    elif "Π_cg" in top2_names:
        label = "RG closure / scale consistency"
    else:
        label = "(mixed)"
    pc_labels_guess.append(label)
    print(f"    → Interpretation: {label}")

print(f"\n  Summary: 5-term micro functional → {min(3, len(pc_labels_guess))}-PC effective functional")
print(f"    This {'supports' if cumvar_win[2] > 0.90 else 'partially supports'} the renormalization conjecture")
print(f"    (PC1-3 explain {100*cumvar_win[2]:.1f}% of variance)")
