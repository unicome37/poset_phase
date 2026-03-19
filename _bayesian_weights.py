"""Bayesian inference for unified functional weights (β, γ, λ, η, κ).

Model (scale-fixed, β ≡ 1):
  F[X] = 1·log_H + γ·Π_geo - λ·Σ_hist + η·Ξ_d + κ·Π_cg

The overall scale of F is unidentifiable from ranking data alone
(only relative F matters), so we fix β=1 and infer (γ, λ, η, κ).

Likelihood from prediction constraints:
  B: F(Lor2D) < F(KR_like) within each N        →  Bradley-Terry pairwise
  C: corr(Σ_hist, F) < 0 within each N           →  sign constraint
  A: F(Lor4D) < F(Lor2D/3D/5D) at large N        →  ranking (weak)
  R: within-N residual variance minimization      →  Gaussian regression

Implementation: manual Metropolis-Hastings.
"""
import csv
import numpy as np

# ── Load data ──
data = []
with open("outputs_unified_functional/raw_features.csv", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        data.append(row)

families = np.array([d["family"] for d in data])
N_arr = np.array([int(d["N"]) for d in data], dtype=float)
log_H = np.array([float(d["log_H"]) for d in data])
pi_geo = np.array([float(d["pi_geo"]) for d in data])
sigma_hist = np.array([float(d["sigma_hist"]) for d in data])
xi_dim = np.array([float(d["xi_dim"]) for d in data])
pi_cg = np.array([float(d["pi_cg"]) for d in data])

n_data = len(data)
unique_N = sorted(set(N_arr))

# Features matrix (β=1 fixed, so log_H is always included with coeff 1)
# w = [γ, λ, η, κ] → F = log_H + γ·pi_geo - λ·sigma_hist + η·xi_dim + κ·pi_cg
other_features = np.column_stack([pi_geo, sigma_hist, xi_dim, pi_cg])

def compute_F(w):
    """F = log_H + γ·pi_geo - λ·sigma_hist + η·xi_dim + κ·pi_cg.
    w = [γ, -λ, η, κ] (note: -λ so all enter as +)."""
    return log_H + other_features @ w


# ── Precompute pair indices ──

B_pairs = []
for n_val in unique_N:
    lor_idx = np.where((N_arr == n_val) & (families == "Lor2D"))[0]
    kr_idx = np.where((N_arr == n_val) & (families == "KR_like"))[0]
    for i in lor_idx:
        for j in kr_idx:
            B_pairs.append((i, j))

A_pairs = []
for n_val in unique_N:
    if n_val < 28:
        continue
    lor4_idx = np.where((N_arr == n_val) & (families == "Lor4D"))[0]
    for other in ["Lor2D", "Lor3D", "Lor5D"]:
        other_idx = np.where((N_arr == n_val) & (families == other))[0]
        for i in lor4_idx:
            for j in other_idx:
                A_pairs.append((i, j))

C_groups = {}
for n_val in unique_N:
    C_groups[n_val] = np.where(N_arr == n_val)[0]

print(f"Data: {n_data} posets, {len(B_pairs)} B-pairs, {len(A_pairs)} A-pairs, {len(C_groups)} C-groups")
print(f"Scale fixed: β ≡ 1. Inferring (γ, λ, η, κ).")


# ── Log-likelihood ──

def log_likelihood(w, scale_B=2.0, scale_A=0.5, scale_C=2.0):
    F = compute_F(w)
    ll = 0.0
    
    # B: Lor2D should have lower F than KR → sigmoid likelihood
    for i, j in B_pairs:
        diff = (F[j] - F[i]) / scale_B
        diff = np.clip(diff, -30, 30)
        ll += -np.log1p(np.exp(-diff))
    
    # A: Lor4D should have lower F than others (weaker signal)
    for i, j in A_pairs:
        diff = (F[j] - F[i]) / scale_B
        diff = np.clip(diff, -30, 30)
        ll += scale_A * (-np.log1p(np.exp(-diff)))
    
    # C: within-N correlation of Σ_hist with F should be negative
    for n_val, idx in C_groups.items():
        sh = sigma_hist[idx]
        f = F[idx]
        if np.std(sh) > 1e-10 and np.std(f) > 1e-10:
            r = np.corrcoef(sh, f)[0, 1]
            ll += scale_C * (-r)
    
    # Regularization: within-N residual variance of F should be small
    # (tighter F distribution → more deterministic predictions)
    for n_val, idx in C_groups.items():
        f = F[idx]
        # Reward smaller within-family variance
        for fam in ["Lor2D", "Lor3D", "Lor4D", "Lor5D", "KR_like"]:
            fam_idx = idx[families[idx] == fam]
            if len(fam_idx) > 2:
                v = np.var(f[fam_idx - idx[0]] if False else F[fam_idx])
                ll += -0.01 * v  # weak regularizer
    
    return ll


# ── Log-prior (on [γ, -λ, η, κ]) ──

def log_prior(w):
    gamma, neg_lam, eta, kappa = w
    lam = -neg_lam
    
    if gamma < 0 or lam < 0 or eta < 0 or kappa < 0:
        return -np.inf
    if gamma > 10 or lam > 10 or eta > 5 or kappa > 5:
        return -np.inf
    
    lp = 0.0
    lp += -0.5 * ((gamma - 0.25) / 0.5)**2
    lp += -0.5 * ((lam - 0.75) / 1.0)**2
    lp += -0.5 * ((eta - 0.05) / 0.3)**2
    lp += -0.5 * ((kappa - 0.025) / 0.2)**2
    return lp


def log_posterior(w):
    lp = log_prior(w)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(w)
    if not np.isfinite(ll):
        return -np.inf
    return lp + ll


# ── MCMC ──

def run_mcmc(n_steps=80000, burn_in=20000, seed=42):
    rng = np.random.default_rng(seed)
    
    # Start: [γ, -λ, η, κ] = [0.25, -0.75, 0.05, 0.025]
    w0 = np.array([0.25, -0.75, 0.05, 0.025])
    proposal_scale = np.array([0.04, 0.06, 0.015, 0.01])
    
    chain = np.zeros((n_steps, 4))
    log_posts = np.zeros(n_steps)
    w_current = w0.copy()
    lp_current = log_posterior(w_current)
    n_accept = 0
    
    for step in range(n_steps):
        w_prop = w_current + rng.normal(0, proposal_scale)
        lp_prop = log_posterior(w_prop)
        
        if np.log(rng.random()) < (lp_prop - lp_current):
            w_current = w_prop
            lp_current = lp_prop
            n_accept += 1
        
        chain[step] = w_current
        log_posts[step] = lp_current
        
        if step < burn_in and step > 0 and step % 1000 == 0:
            rate = n_accept / (step + 1)
            if rate < 0.15:
                proposal_scale *= 0.8
            elif rate > 0.40:
                proposal_scale *= 1.2
    
    return chain[burn_in:], n_accept / n_steps, chain, log_posts


print("\nRunning MCMC (200000 steps, 40000 burn-in)...")
samples, acc_rate, full_chain, log_posts = run_mcmc(n_steps=200000, burn_in=40000)
print(f"Acceptance rate: {acc_rate:.3f}")
print(f"Posterior samples (before thinning): {len(samples)}")

# Thin by 5 to reduce autocorrelation
samples = samples[::5]
print(f"Posterior samples (after thin-5): {len(samples)}")

# ── Convert to display form ──
param_names = ["γ", "λ", "η", "κ"]
calibrated = [0.25, 0.75, 0.05, 0.025]

display = samples.copy()
display[:, 1] = -display[:, 1]  # -λ → λ

print("\n" + "="*70)
print("POSTERIOR SUMMARY (β ≡ 1 fixed)")
print("="*70)
print(f"\n  {'Param':>6s} {'Calib':>8s} {'Mean':>8s} {'Median':>8s} {'Std':>8s} {'2.5%':>8s} {'97.5%':>8s}")
print(f"  {'-'*6} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")

for i, (name, cal) in enumerate(zip(param_names, calibrated)):
    s = display[:, i]
    print(f"  {name:>6s} {cal:8.4f} {np.mean(s):8.4f} {np.median(s):8.4f} "
          f"{np.std(s):8.4f} {np.percentile(s, 2.5):8.4f} {np.percentile(s, 97.5):8.4f}")

# Full weight vector with β=1
print(f"\n  Full posterior weights (β ≡ 1):")
print(f"    β = 1.000 (fixed)")
for i, name in enumerate(param_names):
    s = display[:, i]
    print(f"    {name} = {np.median(s):.4f}  [{np.percentile(s, 2.5):.4f}, {np.percentile(s, 97.5):.4f}]")

# Original scale (multiply by 2 to recover β=2 calibration)
print(f"\n  Rescaled to β = 2.0 for comparison with calibrated values:")
print(f"    β = 2.000")
for i, (name, cal) in enumerate(zip(param_names, [0.5, 1.5, 0.1, 0.05])):
    s = display[:, i] * 2.0
    print(f"    {name} = {np.median(s):.3f}  [{np.percentile(s, 2.5):.3f}, {np.percentile(s, 97.5):.3f}]  (calib: {cal:.3f})")


# ── Correlation matrix ──
print("\n  Posterior correlation matrix:")
corr = np.corrcoef(display.T)
print(f"  {'':>6s}", end="")
for name in param_names:
    print(f" {name:>6s}", end="")
print()
for i, name in enumerate(param_names):
    print(f"  {name:>6s}", end="")
    for j in range(4):
        print(f" {corr[i,j]:+6.2f}", end="")
    print()


# ── Identifiability ──
print("\n" + "="*70)
print("IDENTIFIABILITY DIAGNOSTICS")
print("="*70)

prior_stds = [0.5, 1.0, 0.3, 0.2]
for i, (name, ps) in enumerate(zip(param_names, prior_stds)):
    post_std = np.std(display[:, i])
    ratio = post_std / ps
    status = "WELL-CONSTRAINED" if ratio < 0.5 else "partially constrained" if ratio < 0.9 else "poorly constrained"
    print(f"  {name}: post_σ={post_std:.4f}, prior_σ={ps:.3f}, ratio={ratio:.2f} → {status}")

# ESS
def ess(x):
    n = len(x)
    x = x - np.mean(x)
    var = np.var(x)
    if var < 1e-20:
        return n
    acf = np.correlate(x, x, mode='full')[n-1:]
    acf = acf / (var * n)
    tau = 1.0
    for k in range(1, min(n//2, 500)):
        if acf[k] < 0:
            break
        tau += 2 * acf[k]
    return n / tau

print(f"\n  Effective sample sizes:")
for i, name in enumerate(param_names):
    print(f"    {name}: ESS = {ess(display[:, i]):.0f}")


# ── Posterior predictive check ──
print("\n" + "="*70)
print("POSTERIOR PREDICTIVE CHECK")
print("="*70)

rng = np.random.default_rng(123)
check_idx = rng.choice(len(samples), size=min(500, len(samples)), replace=False)

B_rates, C_rates, A_rates = [], [], []
for idx in check_idx:
    w = samples[idx]
    F = compute_F(w)
    
    B_rates.append(sum(1 for i, j in B_pairs if F[i] < F[j]) / max(len(B_pairs), 1))
    
    c_neg, c_tot = 0, 0
    for n_val, grp_idx in C_groups.items():
        sh = sigma_hist[grp_idx]
        f = F[grp_idx]
        if np.std(sh) > 1e-10 and np.std(f) > 1e-10:
            c_tot += 1
            if np.corrcoef(sh, f)[0, 1] < 0:
                c_neg += 1
    C_rates.append(c_neg / max(c_tot, 1))
    
    a_wins = sum(1 for i, j in A_pairs if F[i] < F[j])
    A_rates.append(a_wins / max(len(A_pairs), 1))

print(f"  B (Lor2D < KR):   {np.mean(B_rates):.3f} ± {np.std(B_rates):.3f}")
print(f"  C (corr < 0):     {np.mean(C_rates):.3f} ± {np.std(C_rates):.3f}")
print(f"  A (Lor4D best):   {np.mean(A_rates):.3f} ± {np.std(A_rates):.3f}")


# ── Sign constraints check ──
print("\n" + "="*70)
print("SIGN CONSTRAINTS (percentage of posterior satisfying each)")
print("="*70)

pct_gamma_pos = np.mean(display[:, 0] > 0) * 100
pct_lambda_pos = np.mean(display[:, 1] > 0) * 100
pct_eta_pos = np.mean(display[:, 2] > 0) * 100
pct_kappa_pos = np.mean(display[:, 3] > 0) * 100

print(f"  γ > 0: {pct_gamma_pos:.1f}%")
print(f"  λ > 0: {pct_lambda_pos:.1f}%")
print(f"  η ≥ 0: {pct_eta_pos:.1f}%")
print(f"  κ ≥ 0: {pct_kappa_pos:.1f}%")
