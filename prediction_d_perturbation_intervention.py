"""Prediction D: CG Perturbation Quasi-Intervention Experiment.

Analog of Prediction C's Layer Split — a quasi-causal test asking:
    "If we *intervene* on poset structure, does the resulting ΔI_cg
     predict the resulting Δimprove_rank?"

Design:
  1. Reconstruct original posets from seeds (matching confirmatory rep3).
  2. Perturb structure: randomly remove cover relations from the Hasse diagram,
     making each poset "less ordered" in a controlled dose.
  3. Run the full CG pipeline on both original and perturbed posets.
  4. Compute ΔI_cg and Δimprove_rank per family.
  5. Three tests:
       A. Rigidity: higher I_cg_orig → smaller |Δpenalty_cg|
       B. Dose-intervention: ΔI_cg predicts Δimprove_rank
       C. Dose-strength: higher perturbation → larger |ΔI_cg|

Perturbation levels: remove 5%, 10%, 20% of cover relations.
Frozen window: N=30, keep_ratio=0.6, gamma=0.2 (primary).
Also: N=40 and N=52 for scaling check.
"""

from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Imports from existing codebase
# ---------------------------------------------------------------------------
from experiment import FAMILIES
from coarse_grain import coarse_grain_delete_nodes
from stability import (
    signature_dict,
    self_drift,
    coarse_grain_penalty,
    family_centroids,
    nearest_family,
    SIGNATURE_COLUMNS,
)
from runtime_utils import estimate_entropy
from action import action_value, get_action_penalty
from generators import Poset

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SEED_OFFSET = 172000  # rep3
N_VALUES = (30, 40, 52)
FAMILY_NAMES = (
    "lorentzian_like_2d",
    "multi_layer_random",
    "KR_like",
    "interval_order",
    "random_layered_k6_uniform",
    "transitive_percolation",
)
SAMPLES_PER_FAMILY = 8
SIS_RUNS = 64  # sufficient for N≤52
KEEP_RATIO = 0.6
GAMMA = 0.2
BETA = 1.0
ACTION_MODE = "A2"
CG_REPEATS = 5  # reduced from 10 for speed; still enough for mean
PERTURB_FRACS = (0.0, 0.05, 0.10, 0.20)  # 0.0 = original baseline
OUT_DIR = Path("outputs_exploratory/prediction_d_perturbation")

# Zeta values for ranking (subset of full scan, matches frozen window)
ZETA_VALUES = [0.5, 1.0, 1.5, 2.0]


# ===================================================================
# Hasse diagram and perturbation
# ===================================================================

def hasse_covers(closure: np.ndarray) -> np.ndarray:
    """Compute the Hasse diagram (cover matrix) from a transitive closure.

    covers[i, j] = True  iff  closure[i, j] = True  AND  there is no
    intermediate k with closure[i, k] and closure[k, j].
    """
    n = closure.shape[0]
    covers = closure.copy()
    np.fill_diagonal(covers, False)
    # Remove non-cover relations: if i<k<j then i→j is not a cover
    for k in range(n):
        for i in range(n):
            if not covers[i, k]:
                continue
            for j in range(n):
                if covers[k, j]:
                    covers[i, j] = False
    return covers


def transitive_closure(adj: np.ndarray) -> np.ndarray:
    """Compute transitive closure from an adjacency/cover matrix via Warshall."""
    n = adj.shape[0]
    reach = adj.copy()
    for k in range(n):
        reach |= np.outer(reach[:, k], reach[k, :])
    return reach


def perturb_remove_covers(
    closure: np.ndarray,
    frac: float,
    seed: int,
) -> np.ndarray:
    """Remove a fraction of cover relations, then recompute transitive closure.

    This makes the poset "less ordered" — some comparable pairs become
    incomparable.  The coverage removal is the intervention; the resulting
    ΔI_cg is the treatment effect.
    """
    if frac <= 0.0:
        return closure.copy()

    rng = np.random.default_rng(seed)
    covers = hasse_covers(closure)
    # Collect cover edges as (i, j) pairs
    edges = list(zip(*np.where(covers)))
    n_edges = len(edges)
    if n_edges == 0:
        return closure.copy()

    n_remove = max(1, int(round(frac * n_edges)))
    n_remove = min(n_remove, n_edges - 1)  # keep at least one edge

    remove_idx = rng.choice(n_edges, size=n_remove, replace=False)
    for idx in remove_idx:
        i, j = edges[idx]
        covers[i, j] = False

    # Recompute transitive closure from the remaining covers
    new_closure = transitive_closure(covers)
    np.fill_diagonal(new_closure, False)
    return new_closure


# ===================================================================
# CG pipeline (lightweight, self-contained)
# ===================================================================

def run_cg_pipeline_for_posets(
    posets: dict[tuple[int, str, int], Poset],
    centroids: dict[tuple[int, str], np.ndarray],
    base_family_ranks: dict[tuple[int, str], int],
    *,
    keep_ratio: float,
    cg_repeats: int,
    sis_runs: int,
    gamma: float,
    beta: float,
    action_mode: str,
    zeta_values: list[float],
    perturb_label: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run CG → compute I_cg, score, ranking, improve_rank for a set of posets.

    Returns:
        (family_rankings, sample_cg): family-level rankings + sample-level CG data.
    """
    # Stage 1: compute CG penalty for each (family, n, sample, cg_repeat)
    cg_rows = []
    total = sum(1 for _ in posets) * cg_repeats
    done = 0
    for (n, family, sample_id), poset in sorted(posets.items()):
        sig_before = signature_dict(poset)
        # Original entropy and score (for local ranking)
        log_h_orig, _ = estimate_entropy(poset, sis_runs=sis_runs, seed=n * 1000 + sample_id, exact_threshold=0)
        penalty_local_orig = get_action_penalty(poset, action_mode)
        score_local_orig = action_value(log_h_orig, penalty_local_orig, beta, gamma)

        for cg_rep in range(cg_repeats):
            cg_seed = 500_000 + n * 1000 + sample_id * 100 + int(round(keep_ratio * 100)) * 10 + cg_rep
            cg_poset = coarse_grain_delete_nodes(poset, keep_ratio=keep_ratio, seed=cg_seed)
            log_h_cg, _ = estimate_entropy(cg_poset, sis_runs=sis_runs, seed=cg_seed, exact_threshold=0)
            sig_after = signature_dict(cg_poset)
            drift_val = self_drift(sig_before, sig_after)
            nearest_fam, nearest_dist = nearest_family(sig_after, centroids, n_value=n)
            fam_switch = 0.0 if nearest_fam == family else 1.0
            rank_before = base_family_ranks.get((n, family), 0)
            rank_after = base_family_ranks.get((n, nearest_fam), 0)
            rank_shift = abs(rank_before - rank_after)
            pcg = coarse_grain_penalty(drift_val, fam_switch, float(rank_shift))
            penalty_local_cg = get_action_penalty(cg_poset, action_mode)
            score_local_cg = action_value(log_h_cg, penalty_local_cg, beta, gamma)

            cg_rows.append({
                "n": n, "family": family, "sample_id": sample_id,
                "cg_repeat": cg_rep, "perturb": perturb_label,
                "penalty_cg": pcg, "self_drift": drift_val,
                "family_switch_penalty": fam_switch,
                "rank_shift_penalty": float(rank_shift),
                "score_local_orig": score_local_orig,
                "score_local_cg": score_local_cg,
            })
            done += 1
            if done % 50 == 0:
                print(f"    [{perturb_label}] CG progress: {done}/{total}", flush=True)

    cg_df = pd.DataFrame(cg_rows)

    # Stage 2: average over CG repeats to get mean penalty per (n, family, sample)
    mean_df = (
        cg_df
        .groupby(["n", "family", "sample_id", "perturb"])
        .agg(
            mean_penalty_cg=("penalty_cg", "mean"),
            mean_self_drift=("self_drift", "mean"),
            mean_fam_switch=("family_switch_penalty", "mean"),
            mean_score_local_orig=("score_local_orig", "first"),
        )
        .reset_index()
    )

    # Stage 3: for each zeta, compute I_cg and rankings within blocks
    frames = []
    for zeta in zeta_values:
        block = mean_df.copy()
        block["zeta"] = zeta
        # I_cg = -zscore(mean_penalty_cg) within each (n,)
        for n_val in block["n"].unique():
            mask = block["n"] == n_val
            pcg_arr = block.loc[mask, "mean_penalty_cg"].to_numpy(float)
            std = pcg_arr.std(ddof=0)
            std = std if std > 1e-12 else 1.0
            block.loc[mask, "I_cg"] = -(pcg_arr - pcg_arr.mean()) / std

        # Score_eval = score_local_orig + zeta * mean_penalty_cg (per family-sample)
        block["score_eval"] = block["mean_score_local_orig"] + zeta * block["mean_penalty_cg"]

        # Ranking within each (n, zeta) block across families
        def _rank_families(sub):
            # Average score across samples per family
            fam_scores = sub.groupby("family")["score_eval"].mean()
            fam_scores_local = sub.groupby("family")["mean_score_local_orig"].mean()
            fam_icg = sub.groupby("family")["I_cg"].mean()
            out = pd.DataFrame({
                "mean_score_eval": fam_scores,
                "mean_score_local": fam_scores_local,
                "mean_I_cg": fam_icg,
            })
            out["rank_eval"] = out["mean_score_eval"].rank(method="dense").astype(int)
            out["rank_local_eval"] = out["mean_score_local"].rank(method="dense").astype(int)
            out["improve_rank"] = -(out["rank_eval"] - out["rank_local_eval"])
            return out.reset_index()

        for n_val in block["n"].unique():
            sub = block[block["n"] == n_val]
            ranked = _rank_families(sub)
            ranked["n"] = n_val
            ranked["zeta"] = zeta
            ranked["perturb"] = perturb_label
            frames.append(ranked)

    # Also save sample-level mean CG data (for continuous analysis)
    sample_cg = mean_df.copy()
    # Add I_cg at sample level (within-N z-score of mean_penalty_cg)
    for n_val in sample_cg["n"].unique():
        mask = sample_cg["n"] == n_val
        pcg_arr = sample_cg.loc[mask, "mean_penalty_cg"].to_numpy(float)
        std = pcg_arr.std(ddof=0)
        std = std if std > 1e-12 else 1.0
        sample_cg.loc[mask, "I_cg_sample"] = -(pcg_arr - pcg_arr.mean()) / std

    return pd.concat(frames, ignore_index=True), sample_cg


# ===================================================================
# Main experiment
# ===================================================================

def main():
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║   Prediction D: CG Perturbation Quasi-Intervention                ║")
    print("║   Does ΔI_cg (from structure perturbation) predict Δimprove_rank?  ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")
    print()

    seed_base = SEED_OFFSET * 10_000_000

    # ------------------------------------------------------------------
    # Step 1: Reconstruct original posets and compute centroids
    # ------------------------------------------------------------------
    print("Step 1: Reconstructing original posets...", flush=True)
    posets_orig: dict[tuple[int, str, int], Poset] = {}
    sig_rows = []
    for n in N_VALUES:
        for family in FAMILY_NAMES:
            generator = FAMILIES[family]
            for sid in range(SAMPLES_PER_FAMILY):
                seed = seed_base + 1000 * n + sid
                poset = generator(n=n, seed=seed)
                posets_orig[(n, family, sid)] = poset
                sig = signature_dict(poset)
                sig["n"] = n
                sig["family"] = family
                sig["sample_id"] = sid
                sig_rows.append(sig)

    sig_df = pd.DataFrame(sig_rows)
    centroids = family_centroids(sig_df[["n", "family"] + SIGNATURE_COLUMNS])
    print(f"  Reconstructed {len(posets_orig)} posets across N={N_VALUES}")

    # Compute base family ranks from original scores
    base_ranks: dict[tuple[int, str], int] = {}
    for n in N_VALUES:
        fam_scores = {}
        for family in FAMILY_NAMES:
            scores = []
            for sid in range(SAMPLES_PER_FAMILY):
                poset = posets_orig[(n, family, sid)]
                log_h, _ = estimate_entropy(poset, sis_runs=SIS_RUNS, seed=n * 1000 + sid, exact_threshold=0)
                pen = get_action_penalty(poset, ACTION_MODE)
                sc = action_value(log_h, pen, BETA, GAMMA)
                scores.append(sc)
            fam_scores[family] = np.mean(scores)
        sorted_fams = sorted(fam_scores.items(), key=lambda x: x[1])
        for rank_i, (fam, _) in enumerate(sorted_fams, 1):
            base_ranks[(n, fam)] = rank_i

    print(f"  Base family ranks computed. Time: {time.time() - t0:.1f}s")

    # ------------------------------------------------------------------
    # Step 2: For each perturbation level, create perturbed posets
    # ------------------------------------------------------------------
    all_ranking_frames = []
    all_sample_frames = []

    for perturb_frac in PERTURB_FRACS:
        label = f"p{int(perturb_frac * 100):02d}"
        print(f"\nStep 2-{label}: Perturbation = {perturb_frac:.0%} cover removal...", flush=True)

        posets_pert: dict[tuple[int, str, int], Poset] = {}
        n_covers_removed = []

        for (n, family, sid), poset_orig in sorted(posets_orig.items()):
            if perturb_frac == 0.0:
                posets_pert[(n, family, sid)] = poset_orig
            else:
                pert_seed = 7_777_000 + n * 1000 + sid * 100 + int(perturb_frac * 1000)
                new_closure = perturb_remove_covers(
                    poset_orig.closure, frac=perturb_frac, seed=pert_seed,
                )
                posets_pert[(n, family, sid)] = Poset(new_closure)

                # Track how many covers were removed
                orig_covers = hasse_covers(poset_orig.closure).sum()
                new_covers = hasse_covers(new_closure).sum()
                n_covers_removed.append(orig_covers - new_covers)

        if n_covers_removed:
            print(f"  Covers removed: mean={np.mean(n_covers_removed):.1f}, "
                  f"median={np.median(n_covers_removed):.0f}, "
                  f"range=[{min(n_covers_removed)}, {max(n_covers_removed)}]")

        # Run CG pipeline on the perturbed posets
        ranked, sample_cg = run_cg_pipeline_for_posets(
            posets_pert, centroids, base_ranks,
            keep_ratio=KEEP_RATIO,
            cg_repeats=CG_REPEATS,
            sis_runs=SIS_RUNS,
            gamma=GAMMA,
            beta=BETA,
            action_mode=ACTION_MODE,
            zeta_values=ZETA_VALUES,
            perturb_label=label,
        )
        all_ranking_frames.append(ranked)
        all_sample_frames.append(sample_cg)
        print(f"  Done. Time: {time.time() - t0:.1f}s")

    all_df = pd.concat(all_ranking_frames, ignore_index=True)
    all_df.to_csv(OUT_DIR / "perturbation_rankings.csv", index=False, encoding="utf-8-sig")
    all_sample_df = pd.concat(all_sample_frames, ignore_index=True)
    all_sample_df.to_csv(OUT_DIR / "perturbation_sample_cg.csv", index=False, encoding="utf-8-sig")

    # ------------------------------------------------------------------
    # Step 3: Analysis
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("ANALYSIS")
    print("=" * 70)

    # Get baseline (p00) and merge with perturbed
    baseline = all_df[all_df["perturb"] == "p00"].copy()
    baseline = baseline.rename(columns={
        "mean_I_cg": "I_cg_orig",
        "improve_rank": "improve_rank_orig",
        "mean_score_eval": "score_eval_orig",
    })
    merge_cols = ["family", "n", "zeta"]

    results_rows = []

    for perturb_label in ["p05", "p10", "p20"]:
        frac_val = {"p05": 0.05, "p10": 0.10, "p20": 0.20}[perturb_label]
        pert_data = all_df[all_df["perturb"] == perturb_label].copy()
        merged = pert_data.merge(
            baseline[merge_cols + ["I_cg_orig", "improve_rank_orig"]],
            on=merge_cols,
            suffixes=("_pert", "_orig"),
        )
        merged["delta_I_cg"] = merged["mean_I_cg"] - merged["I_cg_orig"]
        merged["delta_improve_rank"] = merged["improve_rank"] - merged["improve_rank_orig"]

        print(f"\n--- Perturbation: {perturb_label} ({frac_val:.0%} covers removed) ---")

        # Test A: Rigidity — higher I_cg_orig → smaller |delta_I_cg|?
        for n_val in N_VALUES:
            sub = merged[merged["n"] == n_val]
            if len(sub) < 4:
                continue
            # Average across zeta values per family
            fam_agg = sub.groupby("family").agg(
                I_cg_orig=("I_cg_orig", "mean"),
                abs_delta_I_cg=("delta_I_cg", lambda x: np.mean(np.abs(x))),
                delta_I_cg=("delta_I_cg", "mean"),
                delta_improve_rank=("delta_improve_rank", "mean"),
            ).reset_index()

            # Rigidity: corr(I_cg_orig, |delta_I_cg|) should be negative
            from prediction_d_dynamic_validation import spearman
            if len(fam_agg) >= 4:
                rho_rigid = spearman(
                    fam_agg["I_cg_orig"].values,
                    fam_agg["abs_delta_I_cg"].values,
                )
                # Dose-intervention: corr(delta_I_cg, delta_improve_rank) should be positive
                rho_dose = spearman(
                    fam_agg["delta_I_cg"].values,
                    fam_agg["delta_improve_rank"].values,
                )

                print(f"  N={n_val}:")
                print(f"    Test A (rigidity):      rho(I_cg_orig, |ΔI_cg|) = {rho_rigid:+.4f}")
                print(f"    Test B (dose→outcome):  rho(ΔI_cg, Δimprove)    = {rho_dose:+.4f}")

                for _, row in fam_agg.iterrows():
                    print(f"      {row['family']:30s}  I_cg_o={row['I_cg_orig']:+.3f}  "
                          f"|ΔI_cg|={row['abs_delta_I_cg']:.3f}  "
                          f"ΔI_cg={row['delta_I_cg']:+.3f}  "
                          f"Δimpr={row['delta_improve_rank']:+.3f}")

                results_rows.append({
                    "perturb": perturb_label,
                    "perturb_frac": frac_val,
                    "n": n_val,
                    "n_families": len(fam_agg),
                    "rho_rigidity": rho_rigid,
                    "rho_dose_intervention": rho_dose,
                })

    results_df = pd.DataFrame(results_rows)
    results_df.to_csv(OUT_DIR / "perturbation_results.csv", index=False, encoding="utf-8-sig")

    # ------------------------------------------------------------------
    # Summary table
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)
    print(f"{'perturb':>8s}  {'N':>4s}  {'rho_rigidity':>14s}  {'rho_dose_inter':>14s}")
    for _, r in results_df.iterrows():
        print(f"  {r['perturb']:>6s}  {int(r['n']):>4d}  {r['rho_rigidity']:>+14.4f}  {r['rho_dose_intervention']:>+14.4f}")

    # Test C: dose-strength scaling — larger perturbation → larger mean |ΔI_cg|?
    print("\n--- Test C: Dose-strength scaling ---")
    for n_val in N_VALUES:
        sub = results_df[results_df["n"] == n_val]
        if len(sub) < 2:
            continue
        baseline_data = all_df[(all_df["perturb"] == "p00") & (all_df["n"] == n_val)]
        print(f"  N={n_val}:")
        for _, r in sub.iterrows():
            pert_data = all_df[(all_df["perturb"] == r["perturb"]) & (all_df["n"] == n_val)]
            merged_temp = pert_data.merge(
                baseline_data[["family", "zeta", "mean_I_cg"]].rename(columns={"mean_I_cg": "I_cg_orig"}),
                on=["family", "zeta"],
            )
            mean_abs_delta = (merged_temp["mean_I_cg"] - merged_temp["I_cg_orig"]).abs().mean()
            print(f"    {r['perturb']}: mean|ΔI_cg| = {mean_abs_delta:.4f}")

    # ------------------------------------------------------------------
    # Pooled cross-N permutation test for Test B
    # ------------------------------------------------------------------
    print("\n--- Pooled Permutation Test (Test B) ---")
    rng = np.random.default_rng(42)
    n_perm = 50000

    for perturb_label in ["p05", "p10", "p20"]:
        pert_data = all_df[all_df["perturb"] == perturb_label]
        merged_all = pert_data.merge(
            baseline[merge_cols + ["I_cg_orig", "improve_rank_orig"]],
            on=merge_cols,
        )
        merged_all["delta_I_cg"] = merged_all["mean_I_cg"] - merged_all["I_cg_orig"]
        merged_all["delta_improve_rank"] = merged_all["improve_rank"] - merged_all["improve_rank_orig"]

        # Aggregate per family (across N and zeta) for the test
        fam_pool = merged_all.groupby("family").agg(
            delta_I_cg=("delta_I_cg", "mean"),
            delta_improve_rank=("delta_improve_rank", "mean"),
        ).reset_index()

        x = fam_pool["delta_I_cg"].values
        y = fam_pool["delta_improve_rank"].values
        obs_rho = spearman(x, y)
        count = 1
        for _ in range(n_perm):
            if abs(spearman(rng.permutation(x), y)) >= abs(obs_rho):
                count += 1
        p = count / (n_perm + 1)
        print(f"  {perturb_label}: rho(ΔI_cg, Δimprove) = {obs_rho:+.4f}, p = {p:.5f} (n_perm={n_perm})")

    # ------------------------------------------------------------------
    # SAMPLE-LEVEL CONTINUOUS ANALYSIS (the powerful version)
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("SAMPLE-LEVEL CONTINUOUS QUASI-INTERVENTION ANALYSIS")
    print("=" * 70)
    print("(Addresses family-level underpowering: 48 samples per N instead of 6)")

    # Baseline sample CG data
    sample_baseline = all_sample_df[all_sample_df["perturb"] == "p00"].copy()
    sample_baseline = sample_baseline.rename(columns={
        "mean_penalty_cg": "penalty_cg_orig",
        "mean_self_drift": "drift_orig",
        "mean_fam_switch": "switch_orig",
        "I_cg_sample": "I_cg_orig_sample",
        "mean_score_local_orig": "score_local_orig",
    })
    sample_merge_cols = ["n", "family", "sample_id"]

    sample_results = []

    for perturb_label in ["p05", "p10", "p20"]:
        frac_val = {"p05": 0.05, "p10": 0.10, "p20": 0.20}[perturb_label]
        sample_pert = all_sample_df[all_sample_df["perturb"] == perturb_label].copy()
        sm = sample_pert.merge(
            sample_baseline[sample_merge_cols + [
                "penalty_cg_orig", "drift_orig", "switch_orig",
                "I_cg_orig_sample", "score_local_orig",
            ]],
            on=sample_merge_cols,
        )
        sm["delta_penalty"] = sm["mean_penalty_cg"] - sm["penalty_cg_orig"]
        sm["delta_drift"] = sm["mean_self_drift"] - sm["drift_orig"]
        sm["delta_switch"] = sm["mean_fam_switch"] - sm["switch_orig"]
        sm["abs_delta_penalty"] = sm["delta_penalty"].abs()

        print(f"\n--- {perturb_label} ({frac_val:.0%}) ---")

        for n_val in N_VALUES:
            sub = sm[sm["n"] == n_val]
            n_samples = len(sub)
            if n_samples < 6:
                continue

            # Test A_s: Rigidity (sample level)
            # Does higher I_cg → smaller |Δpenalty_cg|?
            rho_rigid = spearman(sub["I_cg_orig_sample"].values, sub["abs_delta_penalty"].values)

            # Test B_s: Dose-intervention (continuous)
            # Does Δpenalty predict how much the CG score changes?
            # score_eval = score_local + zeta * penalty_cg
            # For zeta=1: Δscore_eval = Δscore_local + 1 * Δpenalty_cg
            # So the rank-relevant quantity is the RELATIVE change in score_eval
            # We compute: ΔI_cg_sample = I_cg_pert - I_cg_orig (continuous)
            rho_delta_icg_penalty = spearman(
                sub["I_cg_sample"].values - sub["I_cg_orig_sample"].values,
                sub["delta_penalty"].values,
            )

            # Test B_s2: Within-family — do samples that lose more I_cg also
            # lose more CG benefit? (controls for family identity)
            rho_within_family = []
            for fam in FAMILY_NAMES:
                fsub = sub[sub["family"] == fam]
                if len(fsub) >= 4:
                    delta_icg = fsub["I_cg_sample"].values - fsub["I_cg_orig_sample"].values
                    delta_pen = fsub["delta_penalty"].values
                    rho_within_family.append(spearman(delta_icg, delta_pen))
            mean_within_rho = np.mean(rho_within_family) if rho_within_family else float("nan")

            # Permutation test for rigidity at sample level
            rng_s = np.random.default_rng(n_val * 100 + int(frac_val * 100))
            x_rig = sub["I_cg_orig_sample"].values
            y_rig = sub["abs_delta_penalty"].values
            obs_abs = abs(rho_rigid)
            count_rig = 1
            n_perm_s = 20000
            for _ in range(n_perm_s):
                if abs(spearman(rng_s.permutation(x_rig), y_rig)) >= obs_abs:
                    count_rig += 1
            p_rigid = count_rig / (n_perm_s + 1)

            print(f"  N={n_val} (n={n_samples} samples):")
            print(f"    A_s (rigidity): rho(I_cg_orig, |Δpenalty|) = {rho_rigid:+.4f}, p = {p_rigid:.5f}")
            print(f"    B_s (ΔI_cg~Δpenalty): rho = {rho_delta_icg_penalty:+.4f}")
            print(f"    B_s2 (within-family): mean rho = {mean_within_rho:+.4f}")

            # Per-family detail
            for fam in FAMILY_NAMES:
                fsub = sub[sub["family"] == fam]
                if len(fsub) > 0:
                    print(f"      {fam:30s}  I_cg_o={fsub['I_cg_orig_sample'].mean():+.3f}  "
                          f"Δpen={fsub['delta_penalty'].mean():+.4f}  "
                          f"|Δpen|={fsub['abs_delta_penalty'].mean():.4f}  "
                          f"Δdrift={fsub['delta_drift'].mean():+.4f}  "
                          f"Δswitch={fsub['delta_switch'].mean():+.3f}")

            sample_results.append({
                "perturb": perturb_label,
                "perturb_frac": frac_val,
                "n": n_val,
                "n_samples": n_samples,
                "rho_rigidity_sample": rho_rigid,
                "p_rigidity_sample": p_rigid,
                "rho_delta_icg_penalty": rho_delta_icg_penalty,
                "mean_within_family_rho": mean_within_rho,
            })

    sample_results_df = pd.DataFrame(sample_results)
    sample_results_df.to_csv(OUT_DIR / "perturbation_sample_results.csv", index=False, encoding="utf-8-sig")

    # Summary of sample-level results
    print("\n" + "=" * 70)
    print("SAMPLE-LEVEL SUMMARY")
    print("=" * 70)
    print(f"{'pert':>6s}  {'N':>4s}  {'n':>4s}  {'rho_rigid':>10s}  {'p':>8s}  {'rho_ΔIcg~Δpen':>14s}  {'within_fam':>10s}")
    for _, r in sample_results_df.iterrows():
        sig = "***" if r["p_rigidity_sample"] < 0.001 else ("**" if r["p_rigidity_sample"] < 0.01 else ("*" if r["p_rigidity_sample"] < 0.05 else "ns"))
        print(f"  {r['perturb']:>4s}  {int(r['n']):>4d}  {int(r['n_samples']):>4d}  "
              f"{r['rho_rigidity_sample']:>+10.4f}  {r['p_rigidity_sample']:>7.5f}{sig:>3s}  "
              f"{r['rho_delta_icg_penalty']:>+14.4f}  {r['mean_within_family_rho']:>+10.4f}")

    elapsed = time.time() - t0
    print(f"\nTotal runtime: {elapsed:.1f}s ({elapsed / 60:.1f} min)")
    print(f"Outputs: {OUT_DIR}/")


if __name__ == "__main__":
    main()
