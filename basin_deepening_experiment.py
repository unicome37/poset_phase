"""
Basin Deepening Experiment
===========================

Tracks the "basin depth" of Lor4D's well-bottom position in feature space
as a function of N. The theory predicts that as N grows, the well deepens —
separation between Lor4D and competitors increases (historical sedimentation).

Five complementary quantities:
  1. LSD-Well Margin Ratio: mean_score(runner-up) / mean_score(Lor4D)
  2. Mahalanobis Separation Gap: min(Mahal_dist of non-Lor) - Lor4D self-dist
  3. Effective Volume: det(Σ)^{1/2} of Lor4D uncertainty ellipsoid
  4. Fisher Information: tr(Σ⁻¹) growth
  5. Power-law fit: margin_ratio(N) = α·N^β + c₀

Previously established:
  μ̂(N) = μ(∞) + a/N + b/N²
  Σ̂(N) = diag(A_i · N^{-p_i})
  det(Σ) ∝ N^{-3.38}
"""
from __future__ import annotations

import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

from bd_action import count_intervals_fast
from generators import (
    Poset,
    generate_lorentzian_like_2d,
    generate_lorentzian_like_3d,
    generate_lorentzian_like_4d,
    generate_lorentzian_like_5d,
    generate_kr_like,
    generate_kr_2layer,
    generate_kr_4layer,
    generate_absolute_layered,
    generate_multi_layer_random,
    generate_random_layered_k4_uniform,
    generate_random_layered_k6_uniform,
    generate_random_layered_k8_uniform,
    generate_random_layered_k6_tapered,
    generate_random_layered_k6_middle_heavy,
    generate_random_layered_k6_longjump,
    generate_transitive_percolation,
    generate_interval_order,
)
from unified_functional import compute_xi_dim


# ══════════════════════════════════════════════════════════════════════════
#  Families
# ══════════════════════════════════════════════════════════════════════════

FAMILIES = {
    "Lor2D": generate_lorentzian_like_2d,
    "Lor3D": generate_lorentzian_like_3d,
    "Lor4D": generate_lorentzian_like_4d,
    "Lor5D": generate_lorentzian_like_5d,
    "KR_like": generate_kr_like,
    "KR_2layer": generate_kr_2layer,
    "KR_4layer": generate_kr_4layer,
    "AbsLayer": generate_absolute_layered,
    "MLR": generate_multi_layer_random,
    "RLk4": generate_random_layered_k4_uniform,
    "RLk6": generate_random_layered_k6_uniform,
    "RLk8": generate_random_layered_k8_uniform,
    "RLk6_tap": generate_random_layered_k6_tapered,
    "RLk6_mid": generate_random_layered_k6_middle_heavy,
    "RLk6_lj": generate_random_layered_k6_longjump,
    "TransPerc": generate_transitive_percolation,
    "IntOrder": generate_interval_order,
}


# ══════════════════════════════════════════════════════════════════════════
#  Feature computation (project standard)
# ══════════════════════════════════════════════════════════════════════════

def max_antichain_width(poset: Poset) -> int:
    n = poset.n
    if n <= 1:
        return n
    c = poset.closure
    remaining = set(range(n))
    max_w = 0
    while remaining:
        minimals = []
        for i in remaining:
            is_min = True
            for j in remaining:
                if j != i and c[j, i] and not c[i, j]:
                    is_min = False
                    break
            if is_min:
                minimals.append(i)
        if not minimals:
            break
        max_w = max(max_w, len(minimals))
        for m in minimals:
            remaining.discard(m)
    return max_w


def compute_features(poset: Poset, N: int) -> np.ndarray:
    counts = count_intervals_fast(poset, k_max=5)
    C0 = counts.get(0)
    C1 = counts.get(1)
    c1_c0 = C1 / max(1, C0)
    _, d_eff = compute_xi_dim(poset)
    aw = max_antichain_width(poset)
    width_ratio = aw / max(1, N)
    return np.array([d_eff, c1_c0, width_ratio])


# ══════════════════════════════════════════════════════════════════════════
#  Scoring functions
# ══════════════════════════════════════════════════════════════════════════

def cstar(N: int) -> float:
    return 0.2485 - 2.33 / N


def wstar(N: int) -> float:
    return 0.3255 + 3.80 / N


def lsd_well_score(feats: np.ndarray, N: int) -> float:
    d, c, w = feats
    return 0.5 * (d - 4) ** 2 + 1.0 * (c - cstar(N)) ** 2 + 5.0 * (w - wstar(N)) ** 2


def mahalanobis_score(feats: np.ndarray, mu: np.ndarray, Sigma_inv: np.ndarray) -> float:
    delta = feats - mu
    return float(delta @ Sigma_inv @ delta)


# ══════════════════════════════════════════════════════════════════════════
#  Data generation with adaptive reps
# ══════════════════════════════════════════════════════════════════════════

N_VALUES = [12, 16, 20, 28, 36, 48, 64, 96, 128, 192, 256]
REPS_DEFAULT = 30
REPS_REDUCED = 10
SEED_BASE = 42
TIME_THRESHOLD = 10.0  # seconds per sample → reduce reps


def generate_all_data():
    """Generate feature data for all families × N values with adaptive reps."""
    data = defaultdict(lambda: defaultdict(list))  # data[N][fam] = [feat, ...]
    reps_used = {}  # (N, fam) → actual reps

    for N in N_VALUES:
        print(f"\n--- N = {N} ---")
        # Determine reps: probe timing with one sample from Lor4D
        t_probe = time.time()
        seed_probe = SEED_BASE + N * 100
        try:
            p = generate_lorentzian_like_4d(N, seed=seed_probe)
            _ = compute_features(p, N)
        except Exception:
            pass
        dt_probe = time.time() - t_probe
        reps = REPS_REDUCED if dt_probe > TIME_THRESHOLD else REPS_DEFAULT
        print(f"  Probe time: {dt_probe:.2f}s → reps = {reps}")

        for fam_name, gen_fn in FAMILIES.items():
            t_fam = time.time()
            success = 0
            for rep in range(reps):
                seed = SEED_BASE + hash(fam_name) % 10000 + N * 100 + rep
                seed = seed % (2**31)
                try:
                    poset = gen_fn(N, seed=seed)
                    feat = compute_features(poset, N)
                    if np.all(np.isfinite(feat)):
                        data[N][fam_name].append(feat)
                        success += 1
                except Exception:
                    pass
            dt_fam = time.time() - t_fam
            reps_used[(N, fam_name)] = success
            print(f"  {fam_name:12s}: {success}/{reps} ok  ({dt_fam:.1f}s)")

    return data, reps_used


# ══════════════════════════════════════════════════════════════════════════
#  Analysis
# ══════════════════════════════════════════════════════════════════════════

def analyze(data, reps_used):
    results = []

    for N in N_VALUES:
        fam_data = data[N]
        if "Lor4D" not in fam_data or len(fam_data["Lor4D"]) < 3:
            print(f"  [SKIP] N={N}: insufficient Lor4D data")
            continue

        lor4d_feats = np.array(fam_data["Lor4D"])  # (n_samples, 3)
        mu_lor4d = lor4d_feats.mean(axis=0)
        cov_lor4d = np.cov(lor4d_feats.T)

        # Regularize covariance
        cov_lor4d += np.eye(3) * 1e-10
        try:
            Sigma_inv = np.linalg.inv(cov_lor4d)
        except np.linalg.LinAlgError:
            Sigma_inv = np.linalg.pinv(cov_lor4d)

        # --- Quantity 1: LSD-Well Margin Ratio ---
        cstar_N = cstar(N)
        wstar_N = wstar(N)
        fam_mean_scores = {}
        for fam_name, feats_list in fam_data.items():
            if len(feats_list) == 0:
                continue
            scores = [lsd_well_score(f, N) for f in feats_list]
            fam_mean_scores[fam_name] = np.mean(scores)

        lor4d_score = fam_mean_scores.get("Lor4D", np.nan)
        non_lor4d_scores = {k: v for k, v in fam_mean_scores.items() if k != "Lor4D"}
        if non_lor4d_scores and lor4d_score > 0:
            runner_up_score = min(non_lor4d_scores.values())
            margin_ratio = runner_up_score / lor4d_score
        else:
            margin_ratio = np.nan

        # --- Quantity 2: Mahalanobis Separation Gap ---
        lor4d_self_dists = [mahalanobis_score(f, mu_lor4d, Sigma_inv) for f in lor4d_feats]
        lor4d_mean_self = np.mean(lor4d_self_dists)

        min_non_lor_mahal = np.inf
        for fam_name, feats_list in fam_data.items():
            if fam_name == "Lor4D" or len(feats_list) == 0:
                continue
            dists = [mahalanobis_score(f, mu_lor4d, Sigma_inv) for f in feats_list]
            fam_mean_mahal = np.mean(dists)
            if fam_mean_mahal < min_non_lor_mahal:
                min_non_lor_mahal = fam_mean_mahal

        mahal_gap = min_non_lor_mahal - lor4d_mean_self if np.isfinite(min_non_lor_mahal) else np.nan

        # --- Quantity 3: Effective Volume ---
        det_cov = np.linalg.det(cov_lor4d)
        V_eff = np.sqrt(max(det_cov, 0.0))

        # --- Quantity 4: Fisher Information ---
        tr_Fisher = np.trace(Sigma_inv)

        results.append({
            "N": N,
            "margin_ratio": margin_ratio,
            "mahal_gap": mahal_gap,
            "V_eff": V_eff,
            "tr_Fisher": tr_Fisher,
            "lor4d_score": lor4d_score,
            "runner_up_score": runner_up_score if non_lor4d_scores else np.nan,
            "lor4d_self_mahal": lor4d_mean_self,
            "n_lor4d": len(lor4d_feats),
        })

    return results


# ══════════════════════════════════════════════════════════════════════════
#  Power-law fits
# ══════════════════════════════════════════════════════════════════════════

def fit_power_laws(results):
    fits = {}
    Ns = np.array([r["N"] for r in results], dtype=float)

    # --- Quantity 5: margin_ratio(N) = alpha * N^beta + c0 ---
    mrs = np.array([r["margin_ratio"] for r in results])
    mask = np.isfinite(mrs)
    if mask.sum() >= 3:
        try:
            def margin_model(N, alpha, beta, c0):
                return alpha * N**beta + c0

            popt, _ = curve_fit(margin_model, Ns[mask], mrs[mask],
                                p0=[1.0, 0.5, 1.0], maxfev=10000)
            fits["margin"] = {"alpha": popt[0], "beta": popt[1], "c0": popt[2]}
        except Exception as e:
            fits["margin"] = {"error": str(e)}
    else:
        fits["margin"] = {"error": "insufficient data"}

    # --- V_eff(N) = A * N^{-p} ---
    veffs = np.array([r["V_eff"] for r in results])
    mask_v = np.isfinite(veffs) & (veffs > 0)
    if mask_v.sum() >= 3:
        try:
            log_N = np.log(Ns[mask_v])
            log_V = np.log(veffs[mask_v])

            def linear(x, a, b):
                return a * x + b

            popt_v, _ = curve_fit(linear, log_N, log_V)
            fits["V_eff"] = {"A": np.exp(popt_v[1]), "p": -popt_v[0]}
        except Exception as e:
            fits["V_eff"] = {"error": str(e)}
    else:
        fits["V_eff"] = {"error": "insufficient data"}

    # --- tr_Fisher(N) = A_F * N^{p_F} ---
    tfs = np.array([r["tr_Fisher"] for r in results])
    mask_f = np.isfinite(tfs) & (tfs > 0)
    if mask_f.sum() >= 3:
        try:
            log_N = np.log(Ns[mask_f])
            log_F = np.log(tfs[mask_f])

            def linear(x, a, b):
                return a * x + b

            popt_f, _ = curve_fit(linear, log_N, log_F)
            fits["Fisher"] = {"A_F": np.exp(popt_f[1]), "p_F": popt_f[0]}
        except Exception as e:
            fits["Fisher"] = {"error": str(e)}
    else:
        fits["Fisher"] = {"error": "insufficient data"}

    return fits


# ══════════════════════════════════════════════════════════════════════════
#  Report
# ══════════════════════════════════════════════════════════════════════════

def generate_report(results, fits):
    lines = []
    lines.append("# Basin Deepening Experiment Results")
    lines.append("=" * 70)
    lines.append("")
    lines.append("Theory: As N grows, the Lor4D well deepens — separation from")
    lines.append("competitors increases (historical sedimentation of structure).")
    lines.append("")

    # Table header
    lines.append("## Quantity Table")
    lines.append("")
    hdr = f"{'N':>5s}  {'margin_ratio':>13s}  {'mahal_gap':>10s}  {'V_eff':>12s}  {'tr_Fisher':>12s}  {'n_lor4d':>8s}"
    lines.append(hdr)
    lines.append("-" * len(hdr))

    for r in results:
        line = (
            f"{r['N']:5d}  "
            f"{r['margin_ratio']:13.4f}  "
            f"{r['mahal_gap']:10.4f}  "
            f"{r['V_eff']:12.6e}  "
            f"{r['tr_Fisher']:12.2f}  "
            f"{r['n_lor4d']:8d}"
        )
        lines.append(line)

    lines.append("")

    # Detail: LSD-Well scores
    lines.append("## LSD-Well Score Details")
    lines.append("")
    for r in results:
        lines.append(f"  N={r['N']:4d}: Lor4D={r['lor4d_score']:.6f}  "
                      f"runner-up={r['runner_up_score']:.6f}  "
                      f"ratio={r['margin_ratio']:.4f}")
    lines.append("")

    # Power-law fits
    lines.append("## Power-Law Fits")
    lines.append("")

    mf = fits.get("margin", {})
    if "error" in mf:
        lines.append(f"  margin_ratio(N) = alpha * N^beta + c0 : FIT FAILED ({mf['error']})")
    else:
        lines.append(f"  margin_ratio(N) = {mf['alpha']:.4f} * N^{mf['beta']:.4f} + {mf['c0']:.4f}")
        if mf['beta'] > 0:
            lines.append("    → beta > 0: BASIN IS DEEPENING (margin grows with N)")
        else:
            lines.append("    → beta <= 0: margin not growing as power law")
    lines.append("")

    vf = fits.get("V_eff", {})
    if "error" in vf:
        lines.append(f"  V_eff(N) = A * N^(-p) : FIT FAILED ({vf['error']})")
    else:
        lines.append(f"  V_eff(N) = {vf['A']:.4e} * N^(-{vf['p']:.4f})")
        lines.append(f"    → well narrows as N^(-{vf['p']:.2f})")
    lines.append("")

    ff = fits.get("Fisher", {})
    if "error" in ff:
        lines.append(f"  tr(Fisher)(N) = A * N^p : FIT FAILED ({ff['error']})")
    else:
        lines.append(f"  tr(Fisher)(N) = {ff['A_F']:.4e} * N^{ff['p_F']:.4f}")
        lines.append(f"    → Fisher information grows as N^{ff['p_F']:.2f}")
    lines.append("")

    # Monotonicity check
    lines.append("## Monotonicity Check")
    lines.append("")
    mrs = [r["margin_ratio"] for r in results]
    mono_margin = all(mrs[i] <= mrs[i + 1] for i in range(len(mrs) - 1) if
                      np.isfinite(mrs[i]) and np.isfinite(mrs[i + 1]))
    lines.append(f"  margin_ratio monotonically increasing? {'YES' if mono_margin else 'NO'}")

    gaps = [r["mahal_gap"] for r in results]
    mono_gap = all(gaps[i] <= gaps[i + 1] for i in range(len(gaps) - 1) if
                   np.isfinite(gaps[i]) and np.isfinite(gaps[i + 1]))
    lines.append(f"  mahal_gap monotonically increasing?    {'YES' if mono_gap else 'NO'}")

    veffs = [r["V_eff"] for r in results]
    mono_veff = all(veffs[i] >= veffs[i + 1] for i in range(len(veffs) - 1) if
                    np.isfinite(veffs[i]) and np.isfinite(veffs[i + 1]))
    lines.append(f"  V_eff monotonically decreasing?        {'YES' if mono_veff else 'NO'}")

    tfs = [r["tr_Fisher"] for r in results]
    mono_fisher = all(tfs[i] <= tfs[i + 1] for i in range(len(tfs) - 1) if
                      np.isfinite(tfs[i]) and np.isfinite(tfs[i + 1]))
    lines.append(f"  tr(Fisher) monotonically increasing?   {'YES' if mono_fisher else 'NO'}")
    lines.append("")

    # Conclusion
    lines.append("## Conclusion")
    lines.append("")
    n_mono = sum([mono_margin, mono_gap, mono_veff, mono_fisher])
    if n_mono == 4:
        lines.append("  ALL FOUR indicators confirm monotonic basin deepening.")
    elif n_mono >= 3:
        lines.append(f"  {n_mono}/4 indicators confirm monotonic deepening (strong evidence).")
    elif n_mono >= 2:
        lines.append(f"  {n_mono}/4 indicators confirm monotonic deepening (moderate evidence).")
    else:
        lines.append(f"  Only {n_mono}/4 indicators monotonic — deepening not clearly established.")

    beta_val = fits.get("margin", {}).get("beta", None)
    if beta_val is not None:
        lines.append(f"  Deepening rate exponent: beta = {beta_val:.4f}")
    p_val = fits.get("V_eff", {}).get("p", None)
    if p_val is not None:
        lines.append(f"  Volume shrinkage exponent: p = {p_val:.4f}")
    lines.append("")

    return "\n".join(lines)


# ══════════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("BASIN DEEPENING EXPERIMENT")
    print(f"  N:    {N_VALUES}")
    print(f"  Fams: {len(FAMILIES)}")
    print(f"  Reps: {REPS_DEFAULT} (reduced to {REPS_REDUCED} if >10s/sample)")
    print("=" * 70)

    t0 = time.time()

    # Step 1: Generate data
    data, reps_used = generate_all_data()
    t_gen = time.time() - t0
    print(f"\nData generation complete: {t_gen:.1f}s")

    # Step 2: Analyze
    print("\n--- Analysis ---")
    results = analyze(data, reps_used)

    # Step 3: Power-law fits
    print("\n--- Power-law fits ---")
    fits = fit_power_laws(results)

    # Step 4: Report
    report = generate_report(results, fits)
    print("\n" + report)

    # Step 5: Save
    out_dir = Path("outputs_carlip")
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "basin_deepening_results.txt"
    out_path.write_text(report, encoding="utf-8")
    print(f"\nSaved to {out_path}")
    print(f"Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
