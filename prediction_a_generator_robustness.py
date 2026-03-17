"""
Prediction A — Generator Robustness Verification

Tests whether the 4D selection window (λ=6-8 under link action)
survives under:
  1. Causal diamond sprinkle (Alexandrov set) — more physical than cube sprinkle
  2. Independent seed family (base 1234567 instead of 980000)
  3. Both combined

If 4D still wins unanimously at λ=6-8, the result is generator-independent.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd

from generators import Poset, transitive_closure
from runtime_utils import estimate_entropy_by_family

OUT_DIR = Path("outputs_exploratory/prediction_a_generator_robustness")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ===== Causal Diamond Generators =====
# Alexandrov set J+(p) ∩ J-(q) in d-dim Minkowski space
# p = (0, 0, ..., 0), q = (1, 0, ..., 0)
# Interior: t ∈ (0,1), t² ≥ Σx_i² AND (1-t)² ≥ Σx_i²
# → spatial radius r ≤ min(t, 1-t)

def _sprinkle_causal_diamond(n: int, d_spatial: int, rng: np.random.Generator) -> np.ndarray:
    """Rejection-sample n points uniformly inside d-dim causal diamond.
    Returns array of shape (n, 1+d_spatial): [t, x1, x2, ...]
    """
    points = np.empty((0, 1 + d_spatial))
    while points.shape[0] < n:
        batch_size = max(n * 20, 500)  # oversample for rejection
        t = rng.random(batch_size)
        # spatial: uniform in hypercube [-0.5, 0.5]^d_spatial
        x = rng.random((batch_size, d_spatial)) - 0.5
        r2 = (x ** 2).sum(axis=1)
        # diamond condition: r² ≤ min(t, 1-t)²
        r_max = np.minimum(t, 1.0 - t)
        accept = r2 <= r_max ** 2
        accepted = np.column_stack([t[accept], x[accept]])
        points = np.vstack([points, accepted])
    return points[:n]


def generate_diamond_2d(n: int, seed: int | None = None) -> Poset:
    """2D causal diamond sprinkle (1+1 Minkowski)."""
    rng = np.random.default_rng(seed)
    pts = _sprinkle_causal_diamond(n, 1, rng)
    t = pts[:, 0]
    x = pts[:, 1]
    dt = t[None, :] - t[:, None]
    dx = np.abs(x[None, :] - x[:, None])
    adj = (dt > 0.0) & (dt >= dx)
    return Poset(transitive_closure(adj))


def generate_diamond_3d(n: int, seed: int | None = None) -> Poset:
    """3D causal diamond sprinkle (1+2 Minkowski)."""
    rng = np.random.default_rng(seed)
    pts = _sprinkle_causal_diamond(n, 2, rng)
    t = pts[:, 0]
    dt = t[None, :] - t[:, None]
    dx2 = (pts[None, :, 1] - pts[:, None, 1]) ** 2
    dy2 = (pts[None, :, 2] - pts[:, None, 2]) ** 2
    adj = (dt > 0.0) & (dt ** 2 >= dx2 + dy2)
    return Poset(transitive_closure(adj))


def generate_diamond_4d(n: int, seed: int | None = None) -> Poset:
    """4D causal diamond sprinkle (1+3 Minkowski)."""
    rng = np.random.default_rng(seed)
    pts = _sprinkle_causal_diamond(n, 3, rng)
    t = pts[:, 0]
    dt = t[None, :] - t[:, None]
    spatial_d2 = sum((pts[None, :, i] - pts[:, None, i]) ** 2 for i in range(1, 4))
    adj = (dt > 0.0) & (dt ** 2 >= spatial_d2)
    return Poset(transitive_closure(adj))


def generate_diamond_5d(n: int, seed: int | None = None) -> Poset:
    """5D causal diamond sprinkle (1+4 Minkowski)."""
    rng = np.random.default_rng(seed)
    pts = _sprinkle_causal_diamond(n, 4, rng)
    t = pts[:, 0]
    dt = t[None, :] - t[:, None]
    spatial_d2 = sum((pts[None, :, i] - pts[:, None, i]) ** 2 for i in range(1, 5))
    adj = (dt > 0.0) & (dt ** 2 >= spatial_d2)
    return Poset(transitive_closure(adj))


# ===== Shared utilities =====

def hasse_links(poset: Poset) -> int:
    c = poset.closure.astype(np.uint8, copy=False)
    has_intermediate = (c @ c).astype(bool, copy=False)
    cover = poset.closure & ~has_intermediate
    np.fill_diagonal(cover, False)
    return int(cover.sum())


def interval_counts(poset: Poset, max_k: int = 3) -> dict[str, int]:
    c = poset.closure
    n = poset.n
    counts = {f"C{k}": 0 for k in range(max_k + 1)}
    for i in range(n):
        for j in range(n):
            if not c[i, j] or i == j:
                continue
            between = int(c[i, :].astype(np.uint8) @ c[:, j].astype(np.uint8))
            k = between
            if k <= max_k:
                counts[f"C{k}"] += 1
    return counts


# ===== Experiment configuration =====

N_VALUES = [20, 28, 36, 44, 52, 60, 68]
SAMPLES = 4
SIS_RUNS = 4096
BETA = 1.0
LAMBDA_VALUES = [0.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0]

FAMILY_EXACT_THRESHOLDS = {
    "2d": 104, "3d": 24, "4d": 24, "5d": 24,
}

# --- Generator suites ---
# Suite A: Original cube sprinkle, INDEPENDENT seeds (base 1234567)
from generators import (generate_lorentzian_like_2d, generate_lorentzian_like_3d,
                         generate_lorentzian_like_4d, generate_lorentzian_like_5d)

GENERATOR_SUITES = {
    "cube_indep_seed": {
        "2d": generate_lorentzian_like_2d,
        "3d": generate_lorentzian_like_3d,
        "4d": generate_lorentzian_like_4d,
        "5d": generate_lorentzian_like_5d,
        "seed_base": 1234567,
        "label": "Cube sprinkle, independent seeds",
    },
    "diamond": {
        "2d": generate_diamond_2d,
        "3d": generate_diamond_3d,
        "4d": generate_diamond_4d,
        "5d": generate_diamond_5d,
        "seed_base": 7770000,
        "label": "Causal diamond sprinkle",
    },
}


def run_suite(suite_name: str, suite: dict) -> pd.DataFrame:
    """Run one generator suite and return base observables."""
    print(f"\n{'='*70}")
    print(f"SUITE: {suite['label']} ({suite_name})")
    print(f"{'='*70}")

    rows = []
    seed_base = suite["seed_base"]

    for n in N_VALUES:
        print(f"  N={n}", end="", flush=True)
        for dim_label in ["2d", "3d", "4d", "5d"]:
            gen = suite[dim_label]
            fam_name = f"lorentzian_like_{dim_label}"
            exact_thr = FAMILY_EXACT_THRESHOLDS[dim_label]

            for sid in range(SAMPLES):
                seed = seed_base + 1000 * n + sid
                poset = gen(n=n, seed=seed)

                log_h, method = estimate_entropy_by_family(
                    poset, family=fam_name, sis_runs=SIS_RUNS, seed=seed,
                    default_exact_threshold=24,
                    family_exact_thresholds={fam_name: exact_thr},
                )

                n_links = hasse_links(poset)
                ic = interval_counts(poset, max_k=3)
                n_rel = int(poset.closure.sum())

                rows.append({
                    "suite": suite_name, "n": n, "dim": dim_label,
                    "sample_id": sid, "seed": seed,
                    "log_H": float(log_h), "method": method,
                    "C0_links": n_links,
                    "C1": ic["C1"], "C2": ic["C2"], "C3": ic["C3"],
                    "n_relations": n_rel,
                    "rel_density": n_rel / (n*(n-1)/2) if n > 1 else 0,
                    "S_link_d2": n - 2 * n_links,
                    "S_link_d2_norm": (n - 2 * n_links) / n,
                })
            print(f" [{dim_label}]", end="", flush=True)
        print()

    return pd.DataFrame(rows)


def analyze_winners(df: pd.DataFrame, suite_name: str) -> pd.DataFrame:
    """Compute winners at each (N, λ) for a given suite."""
    dims = ["2d", "3d", "4d", "5d"]
    winner_rows = []

    for n in N_VALUES:
        for lam in LAMBDA_VALUES:
            dim_scores = {}
            for dim in dims:
                sub = df[(df["n"] == n) & (df["dim"] == dim)]
                if len(sub) == 0:
                    continue
                mean_score = (-BETA * sub["log_H"] + lam * sub["S_link_d2_norm"]).mean()
                dim_scores[dim] = mean_score

            if not dim_scores:
                continue
            ranked = sorted(dim_scores, key=dim_scores.get)
            winner = ranked[0]
            runner = ranked[1] if len(ranked) > 1 else winner
            margin = dim_scores[runner] - dim_scores[winner]

            winner_rows.append({
                "suite": suite_name, "n": n, "lambda": lam,
                "winner": winner, "margin": margin, "runner_up": runner,
                "ranking": ">".join(ranked),
            })

    return pd.DataFrame(winner_rows)


# ===== Main =====
if __name__ == "__main__":
    all_obs = []
    all_winners = []

    for sname, suite in GENERATOR_SUITES.items():
        obs_df = run_suite(sname, suite)
        all_obs.append(obs_df)

        win_df = analyze_winners(obs_df, sname)
        all_winners.append(win_df)

    obs_all = pd.concat(all_obs, ignore_index=True)
    win_all = pd.concat(all_winners, ignore_index=True)

    obs_all.to_csv(OUT_DIR / "base_observables_robustness.csv", index=False)
    win_all.to_csv(OUT_DIR / "winners_robustness.csv", index=False)

    # ===== Print Results =====
    print("\n" + "=" * 76)
    print("ROBUSTNESS VERIFICATION: WINNER TABLES")
    print("=" * 76)

    for sname in GENERATOR_SUITES:
        suite = GENERATOR_SUITES[sname]
        wdf = win_all[win_all["suite"] == sname]
        print(f"\n--- {suite['label']} ---")

        # Pivot table
        pivot = wdf.pivot(index="lambda", columns="n", values="winner")
        print(f"\n{'λ':>6}", end="")
        for n in N_VALUES:
            print(f"  {n:>4}", end="")
        print()

        for lam in LAMBDA_VALUES:
            row = wdf[wdf["lambda"] == lam]
            line = f"{lam:>6.1f}"
            for n in N_VALUES:
                cell = row[row["n"] == n]
                if len(cell) > 0:
                    w = cell.iloc[0]["winner"].upper()
                    line += f"  {w:>4}"
                else:
                    line += "     -"
            # Mark unanimous 4D
            vals = [row[row["n"] == n].iloc[0]["winner"] for n in N_VALUES if len(row[row["n"] == n]) > 0]
            if all(v == "4d" for v in vals):
                line += "  ★ 4D unanimous"
            print(line)

    # ===== Critical Comparison =====
    print("\n" + "=" * 76)
    print("4D WIN COUNT AT KEY λ VALUES")
    print("=" * 76)

    print(f"\n{'λ':>6}", end="")
    for sname in GENERATOR_SUITES:
        print(f"  {sname:>20}", end="")
    print("  {'original (reference)':>20}")

    # Reference data from previous experiment
    ref_4d_counts = {
        0: 0, 2: 0, 3: 2, 4: 4, 5: 4, 6: 7, 7: 7, 8: 7, 10: 6, 12: 5, 15: 3, 20: 2
    }

    for lam in LAMBDA_VALUES:
        line = f"{lam:>6.1f}"
        for sname in GENERATOR_SUITES:
            wdf = win_all[(win_all["suite"] == sname) & (win_all["lambda"] == lam)]
            n4d = sum(1 for _, r in wdf.iterrows() if r["winner"] == "4d")
            total = len(wdf)
            marker = "★" if n4d == total and total > 0 else " "
            line += f"  {n4d:>2}/{total}{marker:>16}"
            # fix formatting
        ref = ref_4d_counts.get(lam, "?")
        line_fixed = f"{lam:>6.1f}"
        for sname in GENERATOR_SUITES:
            wdf = win_all[(win_all["suite"] == sname) & (win_all["lambda"] == lam)]
            n4d = sum(1 for _, r in wdf.iterrows() if r["winner"] == "4d")
            total = len(wdf)
            marker = "★" if n4d == total and total > 0 else ""
            line_fixed += f"    {n4d}/{total}{marker}"
        line_fixed += f"    {ref}/7" if isinstance(ref, int) else f"    {ref}"
        print(line_fixed)

    # ===== Detailed margins at λ=6 =====
    print("\n" + "=" * 76)
    print("DETAILED MARGINS AT λ=6")
    print("=" * 76)

    for sname in GENERATOR_SUITES:
        suite = GENERATOR_SUITES[sname]
        print(f"\n  {suite['label']}:")
        wdf = win_all[(win_all["suite"] == sname) & (win_all["lambda"] == 6.0)]
        for _, row in wdf.iterrows():
            marker = "★" if row["winner"] == "4d" else " "
            print(f"    N={row['n']:>2}: {row['winner'].upper()} beats "
                  f"{row['runner_up'].upper()} by {row['margin']:.3f} {marker}")

    # ===== Interval profile comparison =====
    print("\n" + "=" * 76)
    print("INTERVAL PROFILE COMPARISON (N=52)")
    print("=" * 76)

    for sname in GENERATOR_SUITES:
        suite = GENERATOR_SUITES[sname]
        sub = obs_all[(obs_all["suite"] == sname) & (obs_all["n"] == 52)]
        print(f"\n  {suite['label']}:")
        print(f"    {'Dim':>4} {'C0':>8} {'C1':>8} {'C2':>8} {'C3':>8} | {'logH':>8} {'S_link/N':>10}")
        for dim in ["2d", "3d", "4d", "5d"]:
            ds = sub[sub["dim"] == dim]
            print(f"    {dim:>4} {ds['C0_links'].mean():>8.1f} {ds['C1'].mean():>8.1f} "
                  f"{ds['C2'].mean():>8.1f} {ds['C3'].mean():>8.1f} | "
                  f"{ds['log_H'].mean():>8.1f} {ds['S_link_d2_norm'].mean():>+10.3f}")

    print(f"\nAll outputs saved to {OUT_DIR}")
