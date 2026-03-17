#!/usr/bin/env python3
"""
Prediction A — d>=6 extrapolation & predictive validation
==========================================================

Workflow:
  Phase 1 (BLIND):  Extrapolate scaling-law parameters to d=6,7
                     using d=2..5 fits. Predict Xi_5->6, Xi_6->7.
  Phase 2 (VERIFY): Run numerical experiments at d=6 (and 7 if feasible)
                     and compute *measured* Xi.
  Phase 3 (REPORT): Compare blind prediction vs measurement.

Depends on: generators.py (6d/7d), runtime_utils.py, observables_geo.py
"""

from __future__ import annotations
import os, sys, time, csv, pathlib
import numpy as np
from scipy.optimize import curve_fit

# ── output ──────────────────────────────────────────────────────────
OUT = pathlib.Path("outputs_exploratory/prediction_a_d6_extrapolation")
OUT.mkdir(parents=True, exist_ok=True)

# ── imports from project ────────────────────────────────────────────
from generators import (
    Poset, transitive_closure,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
    generate_lorentzian_like_6d,
    generate_lorentzian_like_7d,
)
from runtime_utils import estimate_entropy
from observables_geo import geometric_components

# ── constants ───────────────────────────────────────────────────────
SAMPLES = 4            # posets per (d, N, lambda)
SIS_RUNS = 4096        # entropy samples
LAMBDAS = [5, 6, 7, 8, 10]
N_LIST = [20, 36, 52, 68]   # same as original experiments

# d=2..5 fitted parameters from xi_derivation (commit 8c36896)
FIT_D = np.array([2, 3, 4, 5])
FIT_a = np.array([0.5729, 0.3014, 0.1138, 0.0336])
FIT_alpha = np.array([0.3715, 0.6175, 0.8392, 1.0486])
FIT_b = np.array([0.4873, 0.6210, 0.7308, 0.7727])
FIT_c = np.array([-0.3990, -0.3933, -0.5551, -0.5410])
# MC fractions
FIT_p = np.array([0.503, 0.289, 0.175, 0.109])
FIT_l = np.array([0.288, 0.593, 0.840, 0.947])


# ====================================================================
#  PHASE 1: BLIND EXTRAPOLATION
# ====================================================================
def extrapolate_parameters():
    """Fit d-dependence of (a,alpha,b,c,p,l) and extrapolate to d=6,7."""
    results = {}

    # --- a_d: exponential decay  a = A * exp(-kappa * d) ---
    def exp_decay(d, A, kappa):
        return A * np.exp(-kappa * d)
    popt_a, _ = curve_fit(exp_decay, FIT_D, FIT_a, p0=[3.0, 0.9])
    results['a_model'] = f"a_d = {popt_a[0]:.4f} * exp(-{popt_a[1]:.4f} * d)"
    a6 = exp_decay(6, *popt_a)
    a7 = exp_decay(7, *popt_a)
    results['a6'], results['a7'] = a6, a7

    # --- alpha_d: linear  alpha = m*d + q ---
    popt_al = np.polyfit(FIT_D, FIT_alpha, 1)   # [m, q]
    results['alpha_model'] = f"alpha_d = {popt_al[0]:.4f}*d + ({popt_al[1]:.4f})"
    alpha6 = np.polyval(popt_al, 6)
    alpha7 = np.polyval(popt_al, 7)
    results['alpha6'], results['alpha7'] = alpha6, alpha7

    # --- b_d: saturating  b = b_inf - B*exp(-mu*d) ---
    def saturating(d, b_inf, B, mu):
        return b_inf - B * np.exp(-mu * d)
    popt_b, _ = curve_fit(saturating, FIT_D, FIT_b, p0=[0.85, 1.0, 0.5])
    results['b_model'] = (f"b_d = {popt_b[0]:.4f} - {popt_b[1]:.4f}*exp(-{popt_b[2]:.4f}*d)")
    b6 = saturating(6, *popt_b)
    b7 = saturating(7, *popt_b)
    results['b6'], results['b7'] = b6, b7

    # --- c_d: second-order polynomial ---
    popt_c = np.polyfit(FIT_D, FIT_c, 2)
    results['c_model'] = (f"c_d = {popt_c[0]:.4f}*d^2 + ({popt_c[1]:.4f})*d + ({popt_c[2]:.4f})")
    c6 = np.polyval(popt_c, 6)
    c7 = np.polyval(popt_c, 7)
    results['c6'], results['c7'] = c6, c7

    # --- p_d: geometric ratio p_{d+1}/p_d ---
    ratios_p = FIT_p[1:] / FIT_p[:-1]
    r_mean = np.mean(ratios_p)
    p6 = FIT_p[-1] * r_mean
    p7 = p6 * r_mean
    results['p_model'] = f"p_d ~ geometric (ratio={r_mean:.4f})"
    results['p6'], results['p7'] = p6, p7

    # --- l_d: approach 1 from above  l = 1 - L*exp(-nu*d) ---
    def link_sat(d, L, nu):
        return 1.0 - L * np.exp(-nu * d)
    popt_l, _ = curve_fit(link_sat, FIT_D, FIT_l, p0=[2.0, 0.5])
    results['l_model'] = f"l_d = 1 - {popt_l[0]:.4f}*exp(-{popt_l[1]:.4f}*d)"
    l6 = link_sat(6, *popt_l)
    l7 = link_sat(7, *popt_l)
    results['l6'], results['l7'] = l6, l7

    # Store all extrapolated values
    results['params'] = {
        6: {'a': a6, 'alpha': alpha6, 'b': b6, 'c': c6, 'p': p6, 'l': l6},
        7: {'a': a7, 'alpha': alpha7, 'b': b7, 'c': c7, 'p': p7, 'l': l7},
    }

    return results, popt_a, popt_al, popt_b, popt_c, popt_l


def predict_xi(a_lo, alpha_lo, b_lo, c_lo, a_hi, alpha_hi, b_hi, c_hi, N):
    """Closed-form Xi from scaling-law parameters."""
    C0_lo = a_lo * N ** alpha_lo
    C0_hi = a_hi * N ** alpha_hi
    logH_lo = b_lo * np.log(N) + c_lo
    logH_hi = b_hi * np.log(N) + c_hi
    numerator = 2.0 * abs(C0_lo - C0_hi)
    denominator = abs(logH_hi - logH_lo)
    if denominator < 1e-12:
        return float('inf')
    return numerator / denominator


def blind_predictions(ext):
    """Compute blind Xi predictions for d->d+1 transitions."""
    N_ref = 68  # reference N for Xi
    preds = {}

    # Xi_4->5 (sanity check against known 11.8)
    xi45 = predict_xi(
        FIT_a[2], FIT_alpha[2], FIT_b[2], FIT_c[2],
        FIT_a[3], FIT_alpha[3], FIT_b[3], FIT_c[3], N_ref)
    preds['Xi_4->5'] = xi45

    # Xi_5->6
    p6 = ext['params'][6]
    xi56 = predict_xi(
        FIT_a[3], FIT_alpha[3], FIT_b[3], FIT_c[3],
        p6['a'], p6['alpha'], p6['b'], p6['c'], N_ref)
    preds['Xi_5->6'] = xi56

    # Xi_6->7
    p7 = ext['params'][7]
    xi67 = predict_xi(
        p6['a'], p6['alpha'], p6['b'], p6['c'],
        p7['a'], p7['alpha'], p7['b'], p7['c'], N_ref)
    preds['Xi_6->7'] = xi67

    return preds


# ====================================================================
#  PHASE 2: NUMERICAL VERIFICATION
# ====================================================================
def interval_counts_fast(poset: Poset, n: int):
    """Return C0 only (links per element)."""
    closure = poset.closure
    link_mat = closure.copy()
    for k in range(n):
        for i in range(n):
            if not closure[i, k]:
                continue
            for j in range(n):
                if closure[k, j] and i != k and k != j:
                    link_mat[i, j] = False
    C0 = link_mat.sum() / n
    return C0


def run_experiment_for_dim(dim, gen_func, N_list, lambdas, samples, sis_runs):
    """Run link-action + entropy experiment for a given dimension."""
    print(f"\n{'='*60}")
    print(f"  NUMERICAL EXPERIMENT: d={dim}")
    print(f"{'='*60}")

    rows = []
    for N in N_list:
        lam_N = int(round(N * 0.5))  # natural lambda ~ N/2
        for lam in lambdas:
            c0_list = []
            logH_list = []
            for s in range(samples):
                seed = dim * 10000 + N * 100 + lam * 10 + s
                poset = gen_func(N, seed=seed)
                n_eff = poset.closure.shape[0]
                c0 = interval_counts_fast(poset, n_eff)
                c0_per_n = c0 / n_eff

                # entropy (exact if small, SIS otherwise)
                exact_thr = 24  # same as derivation script (d>=3)
                H, _ = estimate_entropy(
                    poset, sis_runs=sis_runs, seed=seed + 7777,
                    exact_threshold=exact_thr)
                logH_per_n = np.log(H) / n_eff if H > 1 else 0.0

                c0_list.append(c0_per_n)
                logH_list.append(logH_per_n)

            row = {
                'd': dim, 'N': N, 'lambda': lam,
                'C0_N_mean': np.mean(c0_list),
                'C0_N_std': np.std(c0_list),
                'logH_N_mean': np.mean(logH_list),
                'logH_N_std': np.std(logH_list),
            }
            rows.append(row)
            print(f"  d={dim} N={N} lam={lam}: "
                  f"C0/N={row['C0_N_mean']:.4f}+/-{row['C0_N_std']:.4f}  "
                  f"logH/N={row['logH_N_mean']:.4f}+/-{row['logH_N_std']:.4f}")
    return rows


def compute_xi_numerical(rows_lo, rows_hi, N_ref=68, lam_ref=8):
    """Compute numerical Xi from experimental data at reference N, lambda."""
    def get_mean(rows, field, N, lam):
        for r in rows:
            if r['N'] == N and r['lambda'] == lam:
                return r[field]
        # fallback: closest N
        dists = [(abs(r['N'] - N), r) for r in rows if r['lambda'] == lam]
        if dists:
            _, best = min(dists, key=lambda x: x[0])
            return best[field]
        return None

    c0_lo = get_mean(rows_lo, 'C0_N_mean', N_ref, lam_ref)
    c0_hi = get_mean(rows_hi, 'C0_N_mean', N_ref, lam_ref)
    h_lo = get_mean(rows_lo, 'logH_N_mean', N_ref, lam_ref)
    h_hi = get_mean(rows_hi, 'logH_N_mean', N_ref, lam_ref)

    if any(v is None for v in [c0_lo, c0_hi, h_lo, h_hi]):
        return None, {}
    num = 2.0 * abs(c0_lo - c0_hi)
    den = abs(h_hi - h_lo)
    xi = num / den if den > 1e-12 else float('inf')
    details = {
        'C0_N_lo': c0_lo, 'C0_N_hi': c0_hi,
        'logH_N_lo': h_lo, 'logH_N_hi': h_hi,
        'numerator': num, 'denominator': den,
    }
    return xi, details


# ====================================================================
#  PHASE 3: FIGURES
# ====================================================================
def make_figures(ext, preds, results_5d, results_6d, results_7d, xi_meas):
    """Create publication figures."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # --- Fig 1: Parameter extrapolation ---
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    d_range = np.linspace(2, 7, 100)

    # a_d
    ax = axes[0, 0]
    ax.semilogy(FIT_D, FIT_a, 'ko', ms=8, label='d=2..5 (fitted)')
    ax.semilogy([6, 7], [ext['params'][6]['a'], ext['params'][7]['a']],
                'rs', ms=10, label='d=6,7 (predicted)')
    ax.set_xlabel('d'); ax.set_ylabel('a_d')
    ax.set_title('Link density amplitude')
    ax.legend()

    # alpha_d
    ax = axes[0, 1]
    ax.plot(FIT_D, FIT_alpha, 'ko', ms=8, label='fitted')
    ax.plot([6, 7], [ext['params'][6]['alpha'], ext['params'][7]['alpha']],
            'rs', ms=10, label='predicted')
    ax.set_xlabel('d'); ax.set_ylabel(r'$\alpha_d$')
    ax.set_title('Link density exponent')
    ax.legend()

    # b_d
    ax = axes[0, 2]
    ax.plot(FIT_D, FIT_b, 'ko', ms=8, label='fitted')
    ax.plot([6, 7], [ext['params'][6]['b'], ext['params'][7]['b']],
            'rs', ms=10, label='predicted')
    ax.set_xlabel('d'); ax.set_ylabel(r'$b_d$')
    ax.set_title('Entropy saturation')
    ax.legend()

    # c_d
    ax = axes[1, 0]
    ax.plot(FIT_D, FIT_c, 'ko', ms=8, label='fitted')
    ax.plot([6, 7], [ext['params'][6]['c'], ext['params'][7]['c']],
            'rs', ms=10, label='predicted')
    ax.set_xlabel('d'); ax.set_ylabel(r'$c_d$')
    ax.set_title('Entropy offset')
    ax.legend()

    # p_d
    ax = axes[1, 1]
    ax.semilogy(FIT_D, FIT_p, 'ko', ms=8, label='fitted')
    ax.semilogy([6, 7], [ext['params'][6]['p'], ext['params'][7]['p']],
                'rs', ms=10, label='predicted')
    ax.set_xlabel('d'); ax.set_ylabel(r'$p_d$')
    ax.set_title('MC ordering fraction')
    ax.legend()

    # l_d
    ax = axes[1, 2]
    ax.plot(FIT_D, FIT_l, 'ko', ms=8, label='fitted')
    ax.plot([6, 7], [ext['params'][6]['l'], ext['params'][7]['l']],
            'rs', ms=10, label='predicted')
    ax.axhline(1.0, color='gray', ls='--', alpha=0.5)
    ax.set_xlabel('d'); ax.set_ylabel(r'$\ell_d$')
    ax.set_title('Link saturation')
    ax.legend()

    fig.suptitle('Scaling-law parameter extrapolation (d=2..5 -> d=6,7)',
                 fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUT / 'fig_parameter_extrapolation.png', dpi=200)
    fig.savefig(OUT / 'fig_parameter_extrapolation.pdf')
    plt.close(fig)
    print(f"  [Fig 1] Parameter extrapolation saved")

    # --- Fig 2: Xi staircase (blind + measured) ---
    fig, ax = plt.subplots(figsize=(8, 5))
    transitions = ['3->4', '4->5', '5->6', '6->7']
    xi_blind = [
        predict_xi(FIT_a[1], FIT_alpha[1], FIT_b[1], FIT_c[1],
                   FIT_a[2], FIT_alpha[2], FIT_b[2], FIT_c[2], 68),
        preds['Xi_4->5'],
        preds['Xi_5->6'],
        preds['Xi_6->7'],
    ]
    x_pos = np.arange(len(transitions))
    bars = ax.bar(x_pos - 0.15, xi_blind, 0.3, color='steelblue',
                  label='Blind prediction (formula)', alpha=0.8)

    # measured
    if xi_meas:
        xi_m = [xi_meas.get(t, 0) for t in transitions]
        ax.bar(x_pos + 0.15, xi_m, 0.3, color='coral',
               label='Numerical measurement', alpha=0.8)

    ax.set_xticks(x_pos)
    ax.set_xticklabels([f'd={t}' for t in transitions])
    ax.set_ylabel(r'$\Xi_{d \to d+1}$')
    ax.set_title('Dimensional barrier staircase: blind prediction vs measurement')
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / 'fig_xi_staircase.png', dpi=200)
    fig.savefig(OUT / 'fig_xi_staircase.pdf')
    plt.close(fig)
    print(f"  [Fig 2] Xi staircase saved")

    # --- Fig 3: Raw data comparison d=5,6,7 ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    colors = {5: 'green', 6: 'blue', 7: 'purple'}
    for d_val, rows in [(5, results_5d), (6, results_6d), (7, results_7d)]:
        if not rows:
            continue
        ns = sorted(set(r['N'] for r in rows))
        # plot at lambda=8
        c0_vals = [np.mean([r['C0_N_mean'] for r in rows
                            if r['N'] == n and r['lambda'] == 8]) for n in ns]
        h_vals = [np.mean([r['logH_N_mean'] for r in rows
                           if r['N'] == n and r['lambda'] == 8]) for n in ns]
        ax1.plot(ns, c0_vals, 'o-', color=colors[d_val], ms=7, label=f'd={d_val}')
        ax2.plot(ns, h_vals, 's-', color=colors[d_val], ms=7, label=f'd={d_val}')

    ax1.set_xlabel('N'); ax1.set_ylabel(r'$C_0 / N$')
    ax1.set_title(r'Link density at $\lambda=8$')
    ax1.legend()
    ax2.set_xlabel('N'); ax2.set_ylabel(r'$\log H / N$')
    ax2.set_title(r'Entropy density at $\lambda=8$')
    ax2.legend()
    fig.suptitle('Observable comparison: d=5,6,7', fontsize=13, fontweight='bold')
    fig.tight_layout()
    fig.savefig(OUT / 'fig_d567_comparison.png', dpi=200)
    fig.savefig(OUT / 'fig_d567_comparison.pdf')
    plt.close(fig)
    print(f"  [Fig 3] d=5,6,7 comparison saved")


# ====================================================================
#  MAIN
# ====================================================================
def main():
    t0 = time.time()

    # ── Phase 1: Blind extrapolation ──
    print("=" * 60)
    print("  PHASE 1: BLIND EXTRAPOLATION (d=6, d=7)")
    print("=" * 60)

    ext, *_ = extrapolate_parameters()
    preds = blind_predictions(ext)

    print("\n--- Extrapolated parameters ---")
    for key in ['a_model', 'alpha_model', 'b_model', 'c_model', 'p_model', 'l_model']:
        print(f"  {key}: {ext[key]}")
    print()
    for d in [6, 7]:
        p = ext['params'][d]
        print(f"  d={d}: a={p['a']:.5f}  alpha={p['alpha']:.4f}  "
              f"b={p['b']:.4f}  c={p['c']:.4f}  p={p['p']:.4f}  l={p['l']:.4f}")

    print("\n--- BLIND PREDICTIONS ---")
    for k, v in preds.items():
        print(f"  {k} = {v:.2f}")
    print(f"\n  (Sanity: Xi_4->5 should be ~11.8, got {preds['Xi_4->5']:.2f})")

    # Save blind predictions to file (timestamp for audit trail)
    blind_file = OUT / 'blind_predictions.txt'
    with open(blind_file, 'w', encoding='utf-8') as f:
        f.write(f"Blind predictions generated at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Based on d=2..5 scaling-law fits\n\n")
        for key in ['a_model', 'alpha_model', 'b_model', 'c_model', 'p_model', 'l_model']:
            f.write(f"{key}: {ext[key]}\n")
        f.write("\nExtrapolated parameters:\n")
        for d in [6, 7]:
            p = ext['params'][d]
            f.write(f"  d={d}: a={p['a']:.6f} alpha={p['alpha']:.4f} "
                    f"b={p['b']:.4f} c={p['c']:.4f} p={p['p']:.4f} l={p['l']:.4f}\n")
        f.write(f"\nBLIND PREDICTIONS:\n")
        for k, v in preds.items():
            f.write(f"  {k} = {v:.4f}\n")

    # ── Phase 2: Numerical experiments ──
    print("\n" + "=" * 60)
    print("  PHASE 2: NUMERICAL VERIFICATION")
    print("=" * 60)

    gen_map = {
        5: generate_lorentzian_like_5d,
        6: generate_lorentzian_like_6d,
        7: generate_lorentzian_like_7d,
    }

    # --- 5D reference (recompute for consistency) ---
    print("\n--- Running d=5 (reference) ---")
    results_5d = run_experiment_for_dim(
        5, gen_map[5], N_LIST, LAMBDAS, SAMPLES, SIS_RUNS)

    # --- 6D ---
    print("\n--- Running d=6 ---")
    results_6d = run_experiment_for_dim(
        6, gen_map[6], N_LIST, LAMBDAS, SAMPLES, SIS_RUNS)

    # --- 7D (smaller N due to sparsity) ---
    print("\n--- Running d=7 ---")
    N_LIST_7 = [20, 36, 52, 68]
    results_7d = run_experiment_for_dim(
        7, gen_map[7], N_LIST_7, LAMBDAS, SAMPLES, SIS_RUNS)

    # ── Compute measured Xi ──
    print("\n" + "=" * 60)
    print("  PHASE 3: COMPARISON")
    print("=" * 60)

    xi_meas = {}

    # Xi_4->5: use stored d=4 values + new d=5
    # (We don't re-run d=4; use known fit from derivation)
    # For d=4, C0/N ~ a4*N^alpha4/N at N=68 => 0.1138*68^0.8392/68
    c0_4_pred = FIT_a[2] * 68 ** FIT_alpha[2] / 68   # per-element
    # Wait, C0/N = a_d * N^alpha_d is already per-element? Let me check.
    # From derivation: C0/N = a_d * N^alpha_d, where C0 is total links
    # Actually it was C0 = a_d * N^(1+alpha_d) => C0/N = a_d * N^alpha_d
    # No, let me re-check. The fit was C0/N = a_d * N^alpha_d.
    # So at N=68, C0_per_N_d4 = 0.1138 * 68^0.8392
    # OK but for measured Xi, we compare *the same* data. Let me just use
    # measured 5d vs 6d vs 7d directly.

    # use lambda=8 as reference
    lam_ref = 8

    xi_56, det_56 = compute_xi_numerical(results_5d, results_6d, N_ref=68, lam_ref=lam_ref)
    xi_67, det_67 = compute_xi_numerical(results_6d, results_7d, N_ref=68, lam_ref=lam_ref)

    if xi_56 is not None:
        xi_meas['5->6'] = xi_56
    if xi_67 is not None:
        xi_meas['6->7'] = xi_67

    # Also compute Xi_4->5 numerically from 5d data vs known 4d curve
    # We'll use the formula values for 4d as "measured" baseline
    # (since d=4 data is known and stable from previous experiments)

    print("\n--- MEASURED Xi VALUES ---")
    for k, v in xi_meas.items():
        print(f"  Xi_{k} (measured) = {v:.2f}")

    print("\n--- COMPARISON: BLIND vs MEASURED ---")
    print(f"  Xi_5->6:  blind = {preds['Xi_5->6']:.2f},  measured = {xi_meas.get('5->6', 'N/A')}")
    print(f"  Xi_6->7:  blind = {preds['Xi_6->7']:.2f},  measured = {xi_meas.get('6->7', 'N/A')}")

    if '5->6' in xi_meas and preds['Xi_5->6'] > 0:
        err_56 = abs(preds['Xi_5->6'] - xi_meas['5->6']) / preds['Xi_5->6'] * 100
        print(f"  Xi_5->6 error: {err_56:.1f}%")
    if '6->7' in xi_meas and preds['Xi_6->7'] > 0:
        err_67 = abs(preds['Xi_6->7'] - xi_meas['6->7']) / preds['Xi_6->7'] * 100
        print(f"  Xi_6->7 error: {err_67:.1f}%")

    # ── Make figures ──
    print("\n--- Generating figures ---")
    make_figures(ext, preds, results_5d, results_6d, results_7d, xi_meas)

    # ── Save full report ──
    report_file = OUT / 'extrapolation_report.txt'
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("="*60 + "\n")
        f.write("  d>=6 EXTRAPOLATION REPORT\n")
        f.write("="*60 + "\n\n")
        f.write("PHASE 1: BLIND PREDICTIONS\n")
        for key in ['a_model', 'alpha_model', 'b_model', 'c_model', 'p_model', 'l_model']:
            f.write(f"  {key}: {ext[key]}\n")
        f.write("\nExtrapolated parameters:\n")
        for d in [6, 7]:
            p = ext['params'][d]
            f.write(f"  d={d}: a={p['a']:.6f} alpha={p['alpha']:.4f} "
                    f"b={p['b']:.4f} c={p['c']:.4f}\n")
        f.write(f"\nBlind predictions:\n")
        for k, v in preds.items():
            f.write(f"  {k} = {v:.4f}\n")

        f.write(f"\nPHASE 2: NUMERICAL RESULTS\n")
        for label, rows in [('d=5', results_5d), ('d=6', results_6d), ('d=7', results_7d)]:
            f.write(f"\n  {label}:\n")
            for r in rows:
                f.write(f"    N={r['N']:3d} lam={r['lambda']:2d}: "
                        f"C0/N={r['C0_N_mean']:.4f}+/-{r['C0_N_std']:.4f} "
                        f"logH/N={r['logH_N_mean']:.4f}+/-{r['logH_N_std']:.4f}\n")

        f.write(f"\nPHASE 3: COMPARISON\n")
        for k, v in xi_meas.items():
            f.write(f"  Xi_{k} (measured) = {v:.4f}\n")
        f.write(f"\n  Xi_5->6:  blind = {preds['Xi_5->6']:.4f}")
        if '5->6' in xi_meas:
            f.write(f",  measured = {xi_meas['5->6']:.4f}")
            err = abs(preds['Xi_5->6'] - xi_meas['5->6']) / preds['Xi_5->6'] * 100
            f.write(f",  error = {err:.1f}%")
        f.write(f"\n  Xi_6->7:  blind = {preds['Xi_6->7']:.4f}")
        if '6->7' in xi_meas:
            f.write(f",  measured = {xi_meas['6->7']:.4f}")
            err = abs(preds['Xi_6->7'] - xi_meas['6->7']) / preds['Xi_6->7'] * 100
            f.write(f",  error = {err:.1f}%")
        f.write("\n")

    # ── Save CSV ──
    csv_file = OUT / 'all_results.csv'
    all_rows = results_5d + results_6d + results_7d
    if all_rows:
        keys = all_rows[0].keys()
        with open(csv_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(all_rows)

    elapsed = time.time() - t0
    print(f"\n  Total elapsed: {elapsed:.0f}s ({elapsed/60:.1f} min)")
    print(f"  Output: {OUT}")
    print("  DONE.")


if __name__ == '__main__':
    main()
