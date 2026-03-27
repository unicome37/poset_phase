"""
Manuscript Supplement Experiments — Addressing Reviewer Concerns
================================================================

Runs 6 focused experiments in one pass to fill argumentation gaps:

  Exp 1 (A1): Counter-factual S_MD — center on Lor2D/Lor3D/Lor5D,
              check if each can uniquely self-select.  The claim:
              only Lor4D's own covariance yields self-selection.

  Exp 2 (A3): Re-fit all scaling laws with bootstrap error bars and R².

  Exp 3 (B1): Functional separation — more N values, bootstrap CI on r.

  Exp 4 (B4): Normalized gap — separate statistical precision from
              physical separation at large N.

  Exp 5 (B5): Feature ablation — all 3 choose-2 subsets + single features.

  Exp 6 (B2/B3): Turn-on CI and gap error bars from existing data.

Output: outputs_carlip/manuscript_supplement_experiments.md
"""
from __future__ import annotations

import sys
import time
from collections import defaultdict
from pathlib import Path
from itertools import combinations

import numpy as np
from scipy import stats as sp_stats

from expanded_family_robustness import (
    ALL_FAMILIES,
    compute_features,
    mahalanobis_score,
)


# ═══════════════════════════════════════════════════════════════════
#  Shared infrastructure
# ═══════════════════════════════════════════════════════════════════

SEED_BASE = 42
REPS = 80  # Reference ensemble size (matches turn-on protocol)
EVAL_REPS = 30  # Evaluation samples per family

# N grid covering all relevant scales
N_GRID_DENSE = [12, 14, 16, 18, 20, 24, 28, 32, 36, 48, 64, 96, 128]
N_GRID_CORE = [16, 28, 48, 64, 96, 128]
N_GRID_BASIN = [12, 16, 20, 28, 36, 48, 64, 96, 128, 192, 256]

LORENTZIAN_CENTERS = ["Lor2D", "Lor3D", "Lor4D", "Lor5D"]

OUT_DIR = Path("outputs_carlip")
OUT_DIR.mkdir(exist_ok=True)

report: list[str] = []


def _header(title: str) -> None:
    report.append(f"\n\n{'='*72}")
    report.append(f"# {title}")
    report.append(f"{'='*72}\n")


def _collect_features(N: int, reps: int, seed_base: int = SEED_BASE) -> dict[str, np.ndarray]:
    """Collect feature vectors for all families at given N."""
    data: dict[str, list[np.ndarray]] = defaultdict(list)
    for fam_name, gen_fn in ALL_FAMILIES.items():
        for rep in range(reps):
            seed = (seed_base + hash(fam_name) % 10000 + N * 100 + rep) % (2**31)
            try:
                poset = gen_fn(N, seed=seed)
                feat = compute_features(poset, N)
                data[fam_name].append(feat)
            except Exception:
                pass
    return {k: np.array(v) for k, v in data.items() if len(v) > 0}


def _make_mahalanobis_scorer(ref_feats: np.ndarray):
    """Build (mu, cov_inv) from reference features."""
    mu = np.mean(ref_feats, axis=0)
    cov = np.cov(ref_feats.T)
    reg = 1e-12 * np.eye(cov.shape[0])
    cov_inv = np.linalg.inv(cov + reg)
    det_cov = float(np.linalg.det(cov + reg))
    return mu, cov_inv, cov, det_cov


def _rank_families(data: dict[str, np.ndarray], mu: np.ndarray, cov_inv: np.ndarray) -> list[tuple[str, float]]:
    """Score and rank families. Returns sorted list of (name, mean_score)."""
    scores = {}
    for fam_name, feats in data.items():
        s = [mahalanobis_score(f, mu, cov_inv) for f in feats]
        scores[fam_name] = float(np.mean(s))
    return sorted(scores.items(), key=lambda x: x[1])


# ═══════════════════════════════════════════════════════════════════
#  Experiment 1: Counter-factual S_MD
# ═══════════════════════════════════════════════════════════════════

def exp1_counterfactual_smd() -> None:
    _header("Exp 1: Counter-Factual S_MD — Self-Selection Test")
    report.append("For each Lorentzian center family X, compute S_MD centered on X")
    report.append("and check whether X uniquely ranks #1 across all 25 families.\n")

    report.append("| N | Center | Center-Rank | #1 Family | Margin | Self-Select? |")
    report.append("|---|--------|-------------|-----------|--------|--------------|")

    for N in N_GRID_CORE:
        print(f"  Exp1: N={N}", flush=True)
        data = _collect_features(N, EVAL_REPS)

        for center_name in LORENTZIAN_CENTERS:
            if center_name not in data or len(data[center_name]) < 5:
                report.append(f"| {N} | {center_name} | — | — | — | SKIP |")
                continue

            # Build S_MD centered on this family
            ref_feats = data[center_name]
            mu, cov_inv, _, _ = _make_mahalanobis_scorer(ref_feats)

            ranked = _rank_families(data, mu, cov_inv)
            center_rank = next(
                (i + 1 for i, (name, _) in enumerate(ranked) if name == center_name), 99
            )
            winner_name, winner_score = ranked[0]
            center_score = next(s for n, s in ranked if n == center_name)

            if winner_name == center_name:
                margin = ranked[1][1] - ranked[0][1] if len(ranked) > 1 else float("inf")
                self_sel = "YES"
            else:
                margin = center_score - winner_score
                self_sel = "NO"

            report.append(
                f"| {N} | {center_name} | {center_rank} | {winner_name} | "
                f"{margin:+.4f} | {self_sel} |"
            )


# ═══════════════════════════════════════════════════════════════════
#  Experiment 2: Scaling Laws with Error Bars
# ═══════════════════════════════════════════════════════════════════

def exp2_scaling_laws() -> None:
    _header("Exp 2: Scaling Laws with Bootstrap Uncertainties")

    N_FIT = [N for N in N_GRID_BASIN if N >= 16]
    n_boot = 500

    # Collect data at all N values
    det_covs = []
    fisher_infos = []
    v_effs = []
    gaps = []
    ns_used = []

    for N in N_FIT:
        print(f"  Exp2: N={N}", flush=True)
        data = _collect_features(N, EVAL_REPS)
        if "Lor4D" not in data or len(data["Lor4D"]) < 5:
            continue

        mu, cov_inv, cov, det_cov = _make_mahalanobis_scorer(data["Lor4D"])
        fisher_info = float(np.trace(cov_inv))
        v_eff = float(np.sqrt(max(det_cov, 1e-100)))

        # Gap to nearest competitor
        ranked = _rank_families(data, mu, cov_inv)
        lor_score = next(s for n, s in ranked if n == "Lor4D")
        runner_scores = [(n, s) for n, s in ranked if n != "Lor4D"]
        gap = runner_scores[0][1] - lor_score if runner_scores else 0.0

        ns_used.append(N)
        det_covs.append(det_cov)
        fisher_infos.append(fisher_info)
        v_effs.append(v_eff)
        gaps.append(gap)

    ns_arr = np.array(ns_used, dtype=float)
    log_n = np.log(ns_arr)

    def _power_law_fit(x_log, y_log):
        """Fit log(y) = a + b*log(x), return (exponent, intercept, R², std_err)."""
        mask = np.isfinite(y_log) & np.isfinite(x_log)
        if mask.sum() < 3:
            return np.nan, np.nan, np.nan, np.nan
        result = sp_stats.linregress(x_log[mask], y_log[mask])
        return result.slope, result.intercept, result.rvalue**2, result.stderr

    def _bootstrap_exponent(x_log, y_log, n_boot=500):
        """Bootstrap the power-law exponent."""
        mask = np.isfinite(y_log)
        xx, yy = x_log[mask], y_log[mask]
        n = len(xx)
        if n < 3:
            return np.nan, np.nan
        exponents = []
        rng = np.random.default_rng(42)
        for _ in range(n_boot):
            idx = rng.integers(0, n, n)
            slope, _, r, _, _ = sp_stats.linregress(xx[idx], yy[idx])
            exponents.append(slope)
        return float(np.mean(exponents)), float(np.std(exponents))

    report.append("## Power-law fits: quantity ~ N^exponent\n")
    report.append("| Quantity | Exponent | Bootstrap SE | R² | Data Points |")
    report.append("|----------|----------|--------------|----|-------------|")

    quantities = [
        ("det(Sigma)", np.log(np.array(det_covs))),
        ("V_eff = sqrt(det)", np.log(np.array(v_effs))),
        ("Fisher I_F = tr(Sigma^-1)", np.log(np.array(fisher_infos))),
        ("Delta_hist (gap)", np.log(np.maximum(np.array(gaps), 1e-30))),
    ]

    for name, log_y in quantities:
        slope, intercept, r2, se = _power_law_fit(log_n, log_y)
        boot_mean, boot_se = _bootstrap_exponent(log_n, log_y, n_boot)
        report.append(
            f"| {name} | {slope:.3f} ± {boot_se:.3f} | {boot_se:.4f} | "
            f"{r2:.4f} | {len(ns_used)} |"
        )

    report.append(f"\nN values used: {ns_used}")

    # Consistency check: V_eff vs sqrt(det)
    if len(det_covs) > 0:
        det_exp, _, _, _ = _power_law_fit(log_n, np.log(np.array(det_covs)))
        veff_exp, _, _, _ = _power_law_fit(log_n, np.log(np.array(v_effs)))
        report.append(f"\nConsistency: det exponent = {det_exp:.3f}, "
                      f"V_eff exponent = {veff_exp:.3f}, "
                      f"½·det_exp = {det_exp/2:.3f} (should match V_eff)")


# ═══════════════════════════════════════════════════════════════════
#  Experiment 3: Functional Separation with more N and CI
# ═══════════════════════════════════════════════════════════════════

def exp3_functional_separation() -> None:
    _header("Exp 3: S_BD vs S_MD Functional Separation — Extended")
    report.append("Pearson r and Spearman rho between S_BD and S_MD across families,")
    report.append("with bootstrap 95% CI.\n")

    from bd_action import count_intervals_fast as _cif, bdg_action_d4_standard

    report.append("| N | Pearson r | 95% CI | Spearman rho | 95% CI | n_families |")
    report.append("|---|-----------|--------|--------------|--------|------------|")

    for N in N_GRID_DENSE:
        if N < 16:
            continue
        print(f"  Exp3: N={N}", flush=True)
        data = _collect_features(N, EVAL_REPS)
        if "Lor4D" not in data or len(data["Lor4D"]) < 5:
            continue

        mu, cov_inv, _, _ = _make_mahalanobis_scorer(data["Lor4D"])

        # Compute S_BD and S_MD for each family
        sbd_means = []
        smd_means = []
        for fam_name, gen_fn in ALL_FAMILIES.items():
            if fam_name not in data:
                continue
            # S_MD
            feats = data[fam_name]
            smd_vals = [mahalanobis_score(f, mu, cov_inv) for f in feats]
            smd_means.append(float(np.mean(smd_vals)))
            # S_BD (compute from posets)
            bd_vals = []
            for rep in range(min(EVAL_REPS, len(feats))):
                seed = (SEED_BASE + hash(fam_name) % 10000 + N * 100 + rep) % (2**31)
                try:
                    poset = gen_fn(N, seed=seed)
                    counts = _cif(poset, k_max=5)
                    bd_val = bdg_action_d4_standard(counts, N)
                    bd_vals.append(bd_val)
                except Exception:
                    pass
            sbd_means.append(float(np.mean(bd_vals)) if bd_vals else 0.0)

        sbd_arr = np.array(sbd_means)
        smd_arr = np.array(smd_means)

        if len(sbd_arr) < 5:
            continue

        r_pearson, p_pearson = sp_stats.pearsonr(sbd_arr, smd_arr)
        r_spearman, p_spearman = sp_stats.spearmanr(sbd_arr, smd_arr)

        # Bootstrap CI
        rng = np.random.default_rng(42)
        n_fam = len(sbd_arr)
        boot_r = []
        boot_rho = []
        for _ in range(1000):
            idx = rng.integers(0, n_fam, n_fam)
            if np.std(sbd_arr[idx]) < 1e-15 or np.std(smd_arr[idx]) < 1e-15:
                continue
            br, _ = sp_stats.pearsonr(sbd_arr[idx], smd_arr[idx])
            bh, _ = sp_stats.spearmanr(sbd_arr[idx], smd_arr[idx])
            boot_r.append(br)
            boot_rho.append(bh)

        r_ci = np.percentile(boot_r, [2.5, 97.5]) if boot_r else [np.nan, np.nan]
        rho_ci = np.percentile(boot_rho, [2.5, 97.5]) if boot_rho else [np.nan, np.nan]

        report.append(
            f"| {N} | {r_pearson:+.3f} | [{r_ci[0]:+.3f}, {r_ci[1]:+.3f}] | "
            f"{r_spearman:+.3f} | [{rho_ci[0]:+.3f}, {rho_ci[1]:+.3f}] | {n_fam} |"
        )


# ═══════════════════════════════════════════════════════════════════
#  Experiment 4: Normalized Gap (Physical vs Statistical)
# ═══════════════════════════════════════════════════════════════════

def exp4_normalized_gap() -> None:
    _header("Exp 4: Physical vs Statistical Gap Decomposition")
    report.append("Separating the Mahalanobis gap into:")
    report.append("  - Euclidean distance (absolute feature deviation)")
    report.append("  - Precision amplification (Sigma^-1 scaling)")
    report.append("If Euclidean distance also grows → genuine physical separation.")
    report.append("If only Mahalanobis grows → statistical precision artifact.\n")

    report.append("| N | Runner-Up | Euclid_dist | Mahal_gap | Sigma_amplif | "
                  "Euclid_trend |")
    report.append("|---|-----------|-------------|-----------|--------------|"
                  "-------------|")

    euclid_dists = []
    mahal_gaps = []
    ns_used = []

    for N in N_GRID_BASIN:
        if N < 14:
            continue
        print(f"  Exp4: N={N}", flush=True)
        data = _collect_features(N, EVAL_REPS)
        if "Lor4D" not in data or len(data["Lor4D"]) < 5:
            continue

        mu, cov_inv, cov, det_cov = _make_mahalanobis_scorer(data["Lor4D"])

        # Find nearest competitor
        ranked = _rank_families(data, mu, cov_inv)
        lor_mahal = next(s for n, s in ranked if n == "Lor4D")
        runner = next(((n, s) for n, s in ranked if n != "Lor4D"), None)
        if runner is None:
            continue

        runner_name, runner_mahal = runner
        mahal_gap = runner_mahal - lor_mahal

        # Euclidean distance of runner-up from Lor4D center
        runner_feats = data[runner_name]
        runner_mean = np.mean(runner_feats, axis=0)
        euclid_dist = float(np.linalg.norm(runner_mean - mu))

        # Amplification factor: how much Sigma^-1 amplifies the gap
        amplification = mahal_gap / max(euclid_dist**2, 1e-30)

        ns_used.append(N)
        euclid_dists.append(euclid_dist)
        mahal_gaps.append(mahal_gap)

        report.append(
            f"| {N} | {runner_name} | {euclid_dist:.6f} | {mahal_gap:.4f} | "
            f"{amplification:.2f} | — |"
        )

    # Fit Euclidean distance trend
    if len(ns_used) >= 3:
        log_n = np.log(np.array(ns_used, dtype=float))
        log_euclid = np.log(np.maximum(np.array(euclid_dists), 1e-30))
        res_e = sp_stats.linregress(log_n, log_euclid)
        report.append(f"\nEuclidean distance scaling: exponent = {res_e.slope:.3f}, "
                      f"R² = {res_e.rvalue**2:.4f}")
        if res_e.slope > -0.1:
            report.append("→ Euclidean distance does NOT vanish — physical separation is genuine")
        else:
            report.append("→ Euclidean distance shrinks — gap growth is partly statistical artifact")

        log_mahal = np.log(np.maximum(np.array(mahal_gaps), 1e-30))
        res_m = sp_stats.linregress(log_n, log_mahal)
        report.append(f"Mahalanobis gap scaling: exponent = {res_m.slope:.3f}, "
                      f"R² = {res_m.rvalue**2:.4f}")


# ═══════════════════════════════════════════════════════════════════
#  Experiment 5: Feature Ablation
# ═══════════════════════════════════════════════════════════════════

def exp5_feature_ablation() -> None:
    _header("Exp 5: Feature Space Ablation")
    report.append("Test S_MD with all subsets of the 3 features:")
    report.append("  Full: (d_eff, C1/C0, w/N)")
    report.append("  Drop-1: each pair of 2 features")
    report.append("  Single: each feature alone\n")

    feature_names = ["d_eff", "C1/C0", "w/N"]

    report.append("| N | Feature Set | Lor4D Rank | Runner-Up | Margin |")
    report.append("|---|-------------|------------|-----------|--------|")

    for N in [28, 64, 128]:
        print(f"  Exp5: N={N}", flush=True)
        data = _collect_features(N, EVAL_REPS)
        if "Lor4D" not in data or len(data["Lor4D"]) < 5:
            continue

        # All subsets: singles (3), pairs (3), full (1)
        subsets = []
        # Singles
        for i in range(3):
            subsets.append(([i], feature_names[i]))
        # Pairs
        for i, j in combinations(range(3), 2):
            subsets.append(([i, j], f"{feature_names[i]}+{feature_names[j]}"))
        # Full
        subsets.append(([0, 1, 2], "ALL 3"))

        for indices, label in subsets:
            # Project features
            data_proj = {k: v[:, indices] for k, v in data.items()}
            lor_proj = data_proj["Lor4D"]

            mu = np.mean(lor_proj, axis=0)
            cov = np.cov(lor_proj.T) if lor_proj.shape[1] > 1 else np.array([[np.var(lor_proj)]])
            if cov.ndim == 0:
                cov = np.array([[float(cov)]])
            reg = 1e-12 * np.eye(cov.shape[0])
            cov_inv = np.linalg.inv(cov + reg)

            ranked = _rank_families(data_proj, mu, cov_inv)
            lor_rank = next(
                (i + 1 for i, (n, _) in enumerate(ranked) if n == "Lor4D"), 99
            )
            lor_score = next(s for n, s in ranked if n == "Lor4D")
            winner_name = ranked[0][0]
            if winner_name == "Lor4D" and len(ranked) > 1:
                margin = ranked[1][1] - ranked[0][1]
            else:
                margin = lor_score - ranked[0][1]

            report.append(
                f"| {N} | {label:20s} | {lor_rank} | {winner_name:12s} | "
                f"{margin:+.4f} |"
            )


# ═══════════════════════════════════════════════════════════════════
#  Experiment 6: Turn-On CI and Gap Error Bars
# ═══════════════════════════════════════════════════════════════════

def exp6_turnon_error_bars() -> None:
    _header("Exp 6: Turn-On CI and Gap Error Bars (Multi-Seed)")

    N_GRID = [12, 14, 16, 18, 20, 24, 28, 32]
    SEED_BASES = [42, 137, 271, 500, 777, 1001, 2023, 3141, 5000, 8888]
    REPS_PER = 80

    report.append("10 independent seeds × REPS=80, full 25-family library\n")
    report.append("| N | #1_rate | Mean_margin | Margin_SE | Min_margin | "
                  "95% CI_margin | Runner_census |")
    report.append("|---|---------|-------------|-----------|------------|"
                  "---------------|---------------|")

    for N in N_GRID:
        print(f"  Exp6: N={N}", flush=True)
        ranks = []
        margins = []
        runners = []

        for seed_base in SEED_BASES:
            data: dict[str, list[np.ndarray]] = defaultdict(list)
            for fam_name, gen_fn in ALL_FAMILIES.items():
                for rep in range(REPS_PER):
                    seed = (seed_base + hash(fam_name) % 10000 + N * 100 + rep) % (2**31)
                    try:
                        poset = gen_fn(N, seed=seed)
                        feat = compute_features(poset, N)
                        data[fam_name].append(feat)
                    except Exception:
                        pass

            lor_feats = data.get("Lor4D", [])
            if len(lor_feats) < 5:
                continue
            lor_arr = np.array(lor_feats)
            mu = np.mean(lor_arr, axis=0)
            cov = np.cov(lor_arr.T)
            cov_inv = np.linalg.inv(cov + 1e-12 * np.eye(3))

            mean_scores = {}
            for fam_name, feats in data.items():
                if not feats:
                    continue
                farr = np.array(feats)
                ss = [mahalanobis_score(f, mu, cov_inv) for f in farr]
                mean_scores[fam_name] = float(np.mean(ss))

            ranked_names = sorted(mean_scores, key=mean_scores.get)
            lor_rank = ranked_names.index("Lor4D") + 1 if "Lor4D" in ranked_names else 99
            ranks.append(lor_rank)

            runner = next((f for f in ranked_names if f != "Lor4D"), "N/A")
            margin = mean_scores.get(runner, 0) - mean_scores.get("Lor4D", 0)
            margins.append(margin)
            if lor_rank != 1:
                runners.append(ranked_names[0])

        if not ranks:
            continue

        n_rank1 = sum(1 for r in ranks if r == 1)
        rate_str = f"{n_rank1}/{len(ranks)}"
        mean_m = float(np.mean(margins))
        se_m = float(np.std(margins) / np.sqrt(len(margins)))
        min_m = float(np.min(margins))
        ci_lo = mean_m - 1.96 * se_m
        ci_hi = mean_m + 1.96 * se_m
        runner_str = ", ".join(runners) if runners else "—"

        report.append(
            f"| {N} | {rate_str} | {mean_m:+.4f} | {se_m:.4f} | "
            f"{min_m:+.4f} | [{ci_lo:+.3f}, {ci_hi:+.3f}] | {runner_str} |"
        )


# ═══════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════

def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        try:
            sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    report.append("# Manuscript Supplement Experiments")
    report.append(f"Date: {time.strftime('%Y-%m-%d %H:%M')}")
    report.append(f"Families: {len(ALL_FAMILIES)} (17 standard + 8 adversarial)")
    report.append(f"Reference REPS: {REPS}\n")

    t0 = time.time()

    exp1_counterfactual_smd()
    exp2_scaling_laws()
    exp3_functional_separation()
    exp4_normalized_gap()
    exp5_feature_ablation()
    exp6_turnon_error_bars()

    elapsed = time.time() - t0
    report.append(f"\n\n---\nTotal runtime: {elapsed:.0f}s")

    outpath = OUT_DIR / "manuscript_supplement_experiments.md"
    outpath.write_text("\n".join(report), encoding="utf-8")
    print(f"\nReport written to: {outpath}")
    print(f"Total runtime: {elapsed:.0f}s")


if __name__ == "__main__":
    main()
