"""Prediction D — Continuous-Y & Rich-Residualized Within-Stratum Test.

Goal
----
Previous n=32 high-power test confirmed that the pooled ρ(Δpenalty_cg,
Δscore_local) is driven by between-stratum composition, not within-stratum
mechanism. This script attempts to recover a within-stratum signal by:

  1. **Multiple continuous Y targets** — Δlog_H, Δpenalty_local, Δgeo_total,
     Δneutral_penalty, Δscore_local  (all continuous, no rank/tie issues).
  2. **Rich residualization** — within each (perturb, N, family) stratum,
     regress out baseline structural features (sig_comp, sig_d_eff,
     sig_height_ratio, sig_width_ratio, sig_degree_var, layer_count,
     mean_layer_gap) and baseline score components (log_H₀, penalty_local₀,
     geo_total₀, penalty_cg₀).
  3. **Boundary-sample sub-analysis** — restrict to samples near the median
     |Δpenalty_cg| within each stratum (the "decision boundary"), where a true
     within-stratum mechanism has the best chance of being visible.

If the residualized within-stratum signal is non-null on continuous targets
(especially on boundary samples), D can be upgraded toward quasi-causal.
If it remains null, the current perturbation design is exhausted for this
purpose.

Uses the same SEED_OFFSET / N_VALUES / FAMILY_NAMES / SAMPLES_PER_FAMILY=32
as the high-power test, so the generated posets are identical and results
are directly comparable.
"""

from __future__ import annotations

import math
import time
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

# ── codebase imports ────────────────────────────────────────────────
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
from observables import (
    comparable_fraction,
    layer_profile,
    neutral_penalty as compute_neutral_penalty,
)
from observables_geo import geometric_components, geometric_penalty

from prediction_d_perturbation_intervention import (
    hasse_covers,
    transitive_closure,
    perturb_remove_covers,
    FAMILY_NAMES,
    N_VALUES,
    KEEP_RATIO,
    GAMMA,
    BETA,
    ACTION_MODE,
    CG_REPEATS,
    SIS_RUNS,
    ZETA_VALUES,
    SEED_OFFSET,
)

# ── configuration ──────────────────────────────────────────────────
SAMPLES_PER_FAMILY = 32
PERTURB_FRACS = (0.0, 0.05, 0.10, 0.20)
OUT_DIR = Path("outputs_exploratory/prediction_d_continuous_residualized")
N_PERM = 100_000
SEED = 42
BOUNDARY_QUANTILE = (0.25, 0.75)  # middle 50% of |Δpenalty_cg| per stratum


# ═══════════════════════════════════════════════════════════════════
# Statistics
# ═══════════════════════════════════════════════════════════════════

def _tieaware_rank(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    n = len(x)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1, dtype=float)
    xs = x[order]
    i = 0
    while i < n:
        j = i + 1
        while j < n and xs[j] == xs[i]:
            j += 1
        if j - i > 1:
            avg = 0.5 * ((i + 1) + j)
            ranks[order[i:j]] = avg
        i = j
    return ranks


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx = _tieaware_rank(x)
    ry = _tieaware_rank(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = math.sqrt(float((rx * rx).sum() * (ry * ry).sum()))
    return float((rx * ry).sum() / denom) if denom > 1e-12 else 0.0


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    xc = x - x.mean()
    yc = y - y.mean()
    denom = math.sqrt(float((xc * xc).sum() * (yc * yc).sum()))
    return float((xc * yc).sum() / denom) if denom > 1e-12 else 0.0


def residualize(y: np.ndarray, design: np.ndarray) -> np.ndarray:
    """OLS residualize y against design matrix (intercept added automatically)."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(design, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    x = np.column_stack([np.ones(len(y)), x])
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ beta


def perm_p_abs_spearman(x: np.ndarray, y: np.ndarray, n_perm: int, seed: int) -> float:
    """Two-sided permutation p-value for |Spearman ρ|."""
    rng = np.random.default_rng(seed)
    obs = abs(spearman(x, y))
    count = 1
    for _ in range(n_perm):
        if abs(spearman(rng.permutation(x), y)) >= obs:
            count += 1
    return count / (n_perm + 1)


def stratified_perm_test(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    strata_cols: list[str],
    n_perm: int,
    seed: int,
) -> tuple[float, float, int]:
    """Stratified weighted-mean Spearman ρ with permutation p-value."""
    rng = np.random.default_rng(seed)
    blocks: list[tuple[np.ndarray, np.ndarray, int]] = []
    for _, sub in df.groupby(strata_cols, sort=False):
        x = sub[x_col].to_numpy(float)
        y = sub[y_col].to_numpy(float)
        if len(x) < 4:
            continue
        blocks.append((x, y, len(x)))
    if not blocks:
        return float("nan"), float("nan"), 0

    total_n = sum(b[2] for b in blocks)
    obs = sum(spearman(x, y) * nn for x, y, nn in blocks) / total_n

    count = 1
    for _ in range(n_perm):
        perm_val = sum(
            spearman(rng.permutation(x), y) * nn for x, y, nn in blocks
        ) / total_n
        if abs(perm_val) >= abs(obs):
            count += 1
    return obs, count / (n_perm + 1), len(blocks)


# ═══════════════════════════════════════════════════════════════════
# Extended feature extraction
# ═══════════════════════════════════════════════════════════════════

def extended_features(poset: Poset, sis_runs: int, seed: int) -> dict[str, float]:
    """Compute all structural + score features for a single poset."""
    sig = signature_dict(poset)
    profile = layer_profile(poset)
    layer_count = len(profile)
    gaps = np.diff(profile)
    mean_layer_gap = float(np.mean(np.abs(gaps))) if len(gaps) > 0 else 0.0

    log_h, _ = estimate_entropy(poset, sis_runs=sis_runs, seed=seed, exact_threshold=0)
    pen_neutral = compute_neutral_penalty(poset)
    geo_comp = geometric_components(poset)
    geo_total = geo_comp["geo_total"]
    pen_local = get_action_penalty(poset, ACTION_MODE)
    score_local = action_value(log_h, pen_local, BETA, GAMMA)

    return {
        # structural signature
        "sig_comp": sig["sig_comp"],
        "sig_d_eff": sig["sig_d_eff"],
        "sig_height_ratio": sig["sig_height_ratio"],
        "sig_width_ratio": sig["sig_width_ratio"],
        "sig_degree_var": sig["sig_degree_var"],
        "layer_count": layer_count,
        "mean_layer_gap": mean_layer_gap,
        # score components
        "log_H": log_h,
        "penalty_neutral": pen_neutral,
        "geo_total": geo_total,
        "penalty_local": pen_local,
        "score_local": score_local,
    }


# ═══════════════════════════════════════════════════════════════════
# CG pipeline (extended — records per-sample features)
# ═══════════════════════════════════════════════════════════════════

def run_extended_cg_pipeline(
    posets: dict[tuple[int, str, int], Poset],
    centroids: dict[tuple[int, str], np.ndarray],
    base_family_ranks: dict[tuple[int, str], int],
    perturb_label: str,
) -> pd.DataFrame:
    """Like the original CG pipeline but records ALL extended features."""
    rows: list[dict] = []
    total = len(posets) * CG_REPEATS
    done = 0

    for (n, family, sample_id), poset in sorted(posets.items()):
        # Extended features of the poset itself
        feat = extended_features(poset, sis_runs=SIS_RUNS, seed=n * 1000 + sample_id)
        sig_before = signature_dict(poset)

        # CG repeats → mean penalty_cg
        cg_penalties = []
        for cg_rep in range(CG_REPEATS):
            cg_seed = 500_000 + n * 1000 + sample_id * 100 + int(round(KEEP_RATIO * 100)) * 10 + cg_rep
            cg_poset = coarse_grain_delete_nodes(poset, keep_ratio=KEEP_RATIO, seed=cg_seed)
            sig_after = signature_dict(cg_poset)
            drift_val = self_drift(sig_before, sig_after)
            nearest_fam, _ = nearest_family(sig_after, centroids, n_value=n)
            fam_switch = 0.0 if nearest_fam == family else 1.0
            rank_before = base_family_ranks.get((n, family), 0)
            rank_after = base_family_ranks.get((n, nearest_fam), 0)
            rank_shift = abs(rank_before - rank_after)
            pcg = coarse_grain_penalty(drift_val, fam_switch, float(rank_shift))
            cg_penalties.append(pcg)
            done += 1
            if done % 200 == 0:
                print(f"    [{perturb_label}] CG progress: {done}/{total}", flush=True)

        row = {
            "n": n,
            "family": family,
            "sample_id": sample_id,
            "perturb": perturb_label,
            "mean_penalty_cg": float(np.mean(cg_penalties)),
        }
        row.update(feat)
        rows.append(row)

    return pd.DataFrame(rows)


# ═══════════════════════════════════════════════════════════════════
# Delta-frame construction
# ═══════════════════════════════════════════════════════════════════

# Columns that represent score/feature targets
SCORE_COLS = ["log_H", "penalty_neutral", "geo_total", "penalty_local", "score_local"]
# Columns used as confounders for residualization (baseline structural features)
CONFOUNDER_COLS = [
    "sig_comp", "sig_d_eff", "sig_height_ratio", "sig_width_ratio",
    "sig_degree_var", "layer_count", "mean_layer_gap",
    "log_H", "penalty_local", "geo_total", "mean_penalty_cg",
]

# Y targets: delta of each score column
Y_TARGETS = {
    "delta_score_local": "score_local",
    "delta_log_H": "log_H",
    "delta_penalty_local": "penalty_local",
    "delta_geo_total": "geo_total",
    "delta_penalty_neutral": "penalty_neutral",
}

X_COL = "delta_penalty_cg"


def build_delta_frame(all_samples: pd.DataFrame) -> pd.DataFrame:
    """Pivot baseline vs perturbed. Return one row per (n, family, sample_id,
    perturb) with X = Δpenalty_cg, Y = each Δscore component, plus baseline
    confounders."""
    base = all_samples[all_samples["perturb"] == "p00"].copy()
    base_cols = ["n", "family", "sample_id"]

    # Rename baseline columns → _0 suffix
    rename_map = {}
    for c in SCORE_COLS + CONFOUNDER_COLS + ["mean_penalty_cg"]:
        if c not in rename_map:
            rename_map[c] = f"{c}_0"
    base_ren = base[base_cols + list(rename_map.keys())].rename(columns=rename_map)

    parts: list[pd.DataFrame] = []
    for pert in ["p05", "p10", "p20"]:
        sub = all_samples[all_samples["perturb"] == pert].copy()
        merged = sub.merge(base_ren, on=base_cols, how="inner")
        merged["perturb"] = pert

        merged[X_COL] = merged["mean_penalty_cg"] - merged["mean_penalty_cg_0"]
        for delta_name, src_col in Y_TARGETS.items():
            merged[delta_name] = merged[src_col] - merged[f"{src_col}_0"]

        parts.append(merged)

    return pd.concat(parts, ignore_index=True)


# ═══════════════════════════════════════════════════════════════════
# Residualization + boundary selection
# ═══════════════════════════════════════════════════════════════════

def add_residualized_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Within each (perturb, n, family) stratum, residualize X and all Y targets
    against the baseline confounder set."""
    confounder_0_cols = [f"{c}_0" for c in CONFOUNDER_COLS if f"{c}_0" in df.columns]
    out_frames: list[pd.DataFrame] = []

    for _, sub in df.groupby(["perturb", "n", "family"], sort=False):
        sub = sub.copy()
        if len(sub) < 5:
            # too few to residualize; keep raw
            sub[f"{X_COL}_resid"] = sub[X_COL]
            for y_name in Y_TARGETS:
                sub[f"{y_name}_resid"] = sub[y_name]
            out_frames.append(sub)
            continue

        design = sub[confounder_0_cols].to_numpy(float)
        # Remove constant columns from design (avoid rank-deficient X)
        col_std = design.std(axis=0)
        keep = col_std > 1e-12
        design = design[:, keep]

        sub[f"{X_COL}_resid"] = residualize(sub[X_COL].to_numpy(float), design)
        for y_name in Y_TARGETS:
            sub[f"{y_name}_resid"] = residualize(sub[y_name].to_numpy(float), design)
        out_frames.append(sub)

    return pd.concat(out_frames, ignore_index=True)


def flag_boundary_samples(df: pd.DataFrame, q_lo: float, q_hi: float) -> pd.DataFrame:
    """Flag samples in the middle quantile range of |Δpenalty_cg| per stratum."""
    flags: list[pd.DataFrame] = []
    for _, sub in df.groupby(["perturb", "n", "family"], sort=False):
        sub = sub.copy()
        abs_x = sub[X_COL].abs()
        lo = abs_x.quantile(q_lo)
        hi = abs_x.quantile(q_hi)
        sub["is_boundary"] = (abs_x >= lo) & (abs_x <= hi)
        flags.append(sub)
    return pd.concat(flags, ignore_index=True)


# ═══════════════════════════════════════════════════════════════════
# Analysis engine
# ═══════════════════════════════════════════════════════════════════

def analyse_one_y(
    df: pd.DataFrame,
    perturb: str,
    y_name: str,
    label: str,
    n_perm: int,
    seed: int,
) -> dict:
    """Run pooled + stratified tests for one Y target on one subsample."""
    sub = df[df["perturb"] == perturb].copy()
    if len(sub) < 10:
        return {"perturb": perturb, "y_target": y_name, "subset": label,
                "n_obs": len(sub), "rho_pooled": float("nan"),
                "p_pooled": float("nan"), "rho_strat": float("nan"),
                "p_strat": float("nan"), "n_blocks": 0,
                "mean_within_rho": float("nan"), "neg_frac": float("nan")}

    x = sub[f"{X_COL}_resid"].to_numpy(float)
    y = sub[f"{y_name}_resid"].to_numpy(float)

    rho_pool = spearman(x, y)
    p_pool = perm_p_abs_spearman(x, y, n_perm, seed)

    rho_strat, p_strat, n_blocks = stratified_perm_test(
        sub, f"{X_COL}_resid", f"{y_name}_resid",
        ["n", "family"], n_perm, seed + 1,
    )

    within_rhos = []
    for _, blk in sub.groupby(["n", "family"], sort=False):
        if len(blk) >= 4:
            within_rhos.append(spearman(
                blk[f"{X_COL}_resid"].to_numpy(float),
                blk[f"{y_name}_resid"].to_numpy(float),
            ))
    mean_wr = float(np.mean(within_rhos)) if within_rhos else float("nan")
    neg_fr = float(np.mean([1.0 if r < 0 else 0.0 for r in within_rhos])) if within_rhos else float("nan")

    return {
        "perturb": perturb,
        "y_target": y_name,
        "subset": label,
        "n_obs": len(sub),
        "rho_pooled": rho_pool,
        "p_pooled": p_pool,
        "rho_strat": rho_strat,
        "p_strat": p_strat,
        "n_blocks": n_blocks,
        "mean_within_rho": mean_wr,
        "neg_frac": neg_fr,
    }


def run_full_analysis(df: pd.DataFrame, n_perm: int, seed: int) -> pd.DataFrame:
    """Run all combinations of perturb × Y target × {all, boundary}."""
    results: list[dict] = []
    total_combos = 3 * len(Y_TARGETS) * 2  # 3 perturbs × 5 targets × 2 subsets
    done = 0

    for pert in ["p05", "p10", "p20"]:
        for y_name in Y_TARGETS:
            for label, mask_col in [("all", None), ("boundary", "is_boundary")]:
                sub = df if mask_col is None else df[df[mask_col]]
                s = seed + hash((pert, y_name, label)) % 100_000
                row = analyse_one_y(sub, pert, y_name, label, n_perm, s)
                results.append(row)
                done += 1
                sig_p = "***" if row["p_pooled"] < 0.001 else (
                    "**" if row["p_pooled"] < 0.01 else (
                    "*" if row["p_pooled"] < 0.05 else "ns"))
                sig_s = "***" if row["p_strat"] < 0.001 else (
                    "**" if row["p_strat"] < 0.01 else (
                    "*" if row["p_strat"] < 0.05 else "ns"))
                print(f"  [{done}/{total_combos}] {pert} {y_name:25s} {label:10s}  "
                      f"ρ_pool={row['rho_pooled']:+.4f} p={row['p_pooled']:.5f}{sig_p:>4s}  "
                      f"ρ_strat={row['rho_strat']:+.4f} p={row['p_strat']:.5f}{sig_s:>4s}  "
                      f"ρ̄_within={row['mean_within_rho']:+.4f} %neg={row['neg_frac']:.0%}")

    return pd.DataFrame(results)


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

def main():
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("╔═══════════════════════════════════════════════════════════════════════╗")
    print("║  Prediction D: Continuous-Y + Rich-Residualized Within-Stratum Test ║")
    print("║  X = Δpenalty_cg (residualized)                                     ║")
    print("║  Y = 5 continuous targets (residualized)                            ║")
    print("║  + boundary-sample sub-analysis                                     ║")
    print("╚═══════════════════════════════════════════════════════════════════════╝")

    seed_base = SEED_OFFSET * 10_000_000

    # ─── Step 1: Generate posets (same seeds as high-power test) ────
    print(f"\nStep 1: Generating {SAMPLES_PER_FAMILY} posets per family "
          f"({len(N_VALUES)} N × {len(FAMILY_NAMES)} families)...", flush=True)

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
    print(f"  Generated {len(posets_orig)} posets. ({time.time() - t0:.1f}s)")

    # Base family ranks (for CG penalty)
    base_ranks: dict[tuple[int, str], int] = {}
    for n in N_VALUES:
        fam_scores = {}
        for family in FAMILY_NAMES:
            scores = []
            for sid in range(SAMPLES_PER_FAMILY):
                poset = posets_orig[(n, family, sid)]
                log_h, _ = estimate_entropy(poset, sis_runs=SIS_RUNS,
                                            seed=n * 1000 + sid, exact_threshold=0)
                pen = get_action_penalty(poset, ACTION_MODE)
                scores.append(action_value(log_h, pen, BETA, GAMMA))
            fam_scores[family] = float(np.mean(scores))
        sorted_fams = sorted(fam_scores.items(), key=lambda x: x[1])
        for rank_i, (fam, _) in enumerate(sorted_fams, 1):
            base_ranks[(n, fam)] = rank_i
    print(f"  Base ranks computed. ({time.time() - t0:.1f}s)")

    # ─── Step 2: Perturbation + extended CG pipeline ───────────────
    all_sample_frames: list[pd.DataFrame] = []
    for perturb_frac in PERTURB_FRACS:
        label = f"p{int(perturb_frac * 100):02d}"
        print(f"\nStep 2-{label}: {perturb_frac:.0%} cover removal + extended features...",
              flush=True)

        posets_pert: dict[tuple[int, str, int], Poset] = {}
        for (n, family, sid), poset_orig in sorted(posets_orig.items()):
            if perturb_frac == 0.0:
                posets_pert[(n, family, sid)] = poset_orig
            else:
                closure = poset_orig.closure.copy()
                new_closure = perturb_remove_covers(
                    closure, perturb_frac,
                    seed=SEED_OFFSET + n * 10000 + sid * 100 + int(perturb_frac * 100),
                )
                posets_pert[(n, family, sid)] = Poset(new_closure)

        sample_df = run_extended_cg_pipeline(
            posets_pert, centroids, base_ranks, perturb_label=label,
        )
        all_sample_frames.append(sample_df)
        print(f"  Done {label}: {len(sample_df)} rows. ({time.time() - t0:.1f}s)")

    all_samples = pd.concat(all_sample_frames, ignore_index=True)
    all_samples.to_csv(OUT_DIR / "extended_sample_features.csv",
                       index=False, encoding="utf-8-sig")
    print(f"\nSaved {len(all_samples)} rows to extended_sample_features.csv")

    # ─── Step 3: Build delta frame ─────────────────────────────────
    print("\nStep 3: Building delta frame (baseline-vs-perturbed)...", flush=True)
    delta = build_delta_frame(all_samples)
    print(f"  Delta frame: {len(delta)} rows, {len(delta.columns)} columns")

    # ─── Step 4: Residualize + boundary flag ───────────────────────
    print("Step 4: Within-stratum residualization + boundary flagging...", flush=True)
    delta = add_residualized_columns(delta)
    delta = flag_boundary_samples(delta, *BOUNDARY_QUANTILE)
    n_boundary = delta["is_boundary"].sum()
    print(f"  Boundary samples: {n_boundary}/{len(delta)} ({n_boundary / len(delta):.0%})")

    delta.to_csv(OUT_DIR / "delta_residualized.csv", index=False, encoding="utf-8-sig")

    # ─── Step 5: Full analysis ─────────────────────────────────────
    print(f"\nStep 5: Running {3 * len(Y_TARGETS) * 2} test combinations "
          f"({N_PERM:,} permutations each)...\n", flush=True)
    results = run_full_analysis(delta, N_PERM, SEED)
    results.to_csv(OUT_DIR / "continuous_residualized_results.csv",
                   index=False, encoding="utf-8-sig")

    # ─── Step 6: Summary ───────────────────────────────────────────
    print(f"\n{'═' * 75}")
    print("SUMMARY — Continuous-Y / Rich-Residualized / Boundary-Subsample")
    print(f"{'═' * 75}")

    for subset_label in ["all", "boundary"]:
        sub = results[results["subset"] == subset_label]
        print(f"\n  ── {subset_label.upper()} samples ──")
        print(f"  {'perturb':>7s}  {'y_target':>25s}  {'n':>5s}  "
              f"{'ρ_pool':>8s}  {'p_pool':>10s}  {'ρ_strat':>8s}  {'p_strat':>10s}  "
              f"{'ρ̄_within':>8s}  {'%neg':>5s}")
        for _, r in sub.iterrows():
            sig_p = ("***" if r["p_pooled"] < 0.001 else
                     "**" if r["p_pooled"] < 0.01 else
                     "*" if r["p_pooled"] < 0.05 else "ns")
            sig_s = ("***" if r["p_strat"] < 0.001 else
                     "**" if r["p_strat"] < 0.01 else
                     "*" if r["p_strat"] < 0.05 else "ns")
            print(f"  {r['perturb']:>7s}  {r['y_target']:>25s}  {int(r['n_obs']):>5d}  "
                  f"{r['rho_pooled']:>+8.4f}  {r['p_pooled']:>9.5f}{sig_p:>4s}  "
                  f"{r['rho_strat']:>+8.4f}  {r['p_strat']:>9.5f}{sig_s:>4s}  "
                  f"{r['mean_within_rho']:>+8.4f}  {r['neg_frac']:>5.0%}")

    # Any significant within-stratum results?
    sig_strat = results[(results["p_strat"] < 0.05)]
    if len(sig_strat) > 0:
        print(f"\n  ★ {len(sig_strat)} test(s) reached within-stratum significance (p<0.05):")
        for _, r in sig_strat.iterrows():
            print(f"    {r['perturb']} / {r['y_target']} / {r['subset']}: "
                  f"ρ_strat={r['rho_strat']:+.4f}, p={r['p_strat']:.5f}")
    else:
        print("\n  ☐ No test reached within-stratum significance at α=0.05.")
        print("    → Current perturbation design exhausted for quasi-causal upgrade.")

    elapsed = time.time() - t0
    print(f"\nTotal runtime: {elapsed:.0f}s ({elapsed / 60:.1f} min)")

    # ─── Step 7: Write report ──────────────────────────────────────
    lines = [
        "# Prediction D: Continuous-Y + Rich-Residualized Within-Stratum Test",
        "",
        "## Design",
        f"- **Samples**: {SAMPLES_PER_FAMILY}/family × {len(FAMILY_NAMES)} families × {len(N_VALUES)} N = "
        f"{SAMPLES_PER_FAMILY * len(FAMILY_NAMES) * len(N_VALUES)} per perturbation level",
        f"- **X**: Δpenalty_cg (cover-removal perturbation)",
        f"- **Y targets**: {', '.join(Y_TARGETS.keys())}",
        f"- **Residualization confounders**: baseline {', '.join(CONFOUNDER_COLS)}",
        f"- **Boundary samples**: |Δpenalty_cg| in [{BOUNDARY_QUANTILE[0]:.0%}, {BOUNDARY_QUANTILE[1]:.0%}] quantile per stratum",
        f"- **Permutations**: {N_PERM:,}",
        "",
        "## Results",
        "",
        results.to_markdown(index=False),
        "",
        "## Interpretation",
        "",
    ]
    if len(sig_strat) > 0:
        lines.append(f"**{len(sig_strat)} test(s) reached within-stratum significance.**")
        lines.append("This suggests a genuine within-stratum mechanism channel worth pursuing.")
    else:
        lines.append("**No test reached within-stratum significance.**")
        lines.append("The current cover-removal perturbation design does not produce detectable")
        lines.append("within-stratum quasi-causal signal, even with continuous targets, rich")
        lines.append("residualization, and boundary-sample focus. D remains 'structural covariation'.")
        lines.append("")
        lines.append("Next steps if quasi-causal upgrade is still desired:")
        lines.append("1. Design CG-specific perturbation (matched rewiring that only affects CG channel)")
        lines.append("2. Expand family count to increase within-stratum heterogeneity")
        lines.append("3. Targeted boundary-sample generation (near family-switching threshold)")

    report_path = OUT_DIR / "continuous_residualized_report.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"\nReport written to {report_path.as_posix()}")


if __name__ == "__main__":
    main()
