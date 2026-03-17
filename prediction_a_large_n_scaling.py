"""
Prediction A — Large-N Finite-Size Scaling + Ξ Closure

Push link-action 4D selection to N = 80, 96, 112 and verify:
  1. Does 4D still win at λ = 6-8?
  2. Does margin grow, plateau, or shrink?
  3. Does Ξ₄→₅ ≈ 10 hold at large N?

This is the confirmatory finite-size scaling needed to close the loop.

Also includes a unification test: correlation between link density
(C₀/N) and geometric consistency penalty across dimensions.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd

from generators import (
    Poset, transitive_closure,
    generate_lorentzian_like_2d, generate_lorentzian_like_3d,
    generate_lorentzian_like_4d, generate_lorentzian_like_5d,
)
from runtime_utils import estimate_entropy_by_family
from observables_geo import geometric_components

OUT_DIR = Path("outputs_exploratory/prediction_a_large_n_scaling")
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ===== Core functions =====

def hasse_links(poset: Poset) -> int:
    c = poset.closure.astype(np.uint8, copy=False)
    has_intermediate = (c @ c).astype(bool, copy=False)
    cover = poset.closure & ~has_intermediate
    np.fill_diagonal(cover, False)
    return int(cover.sum())


def interval_counts_fast(poset: Poset) -> dict[str, int]:
    """Only compute C0 (links) — skip expensive C1/C2/C3 for large N.
    For Ξ analysis we only need C0 and log_H."""
    return {"C0": hasse_links(poset), "C1": 0, "C2": 0, "C3": 0}


# ===== Configuration =====

# Extend beyond N=68 to N=80, 96, 112
N_VALUES = [20, 36, 52, 68, 80, 96, 112]
SAMPLES = 4
SIS_RUNS = 4096
BETA = 1.0
LAMBDA_VALUES = [5.0, 6.0, 7.0, 8.0, 10.0]
SEED_BASE = 980000  # Same as original for consistency

GENERATORS = {
    "2d": generate_lorentzian_like_2d,
    "3d": generate_lorentzian_like_3d,
    "4d": generate_lorentzian_like_4d,
    "5d": generate_lorentzian_like_5d,
}

FAMILY_EXACT_THRESHOLDS = {
    "lorentzian_like_2d": 68,   # Use SIS for N>68 (exact is too slow)
    "lorentzian_like_3d": 24,
    "lorentzian_like_4d": 24,
    "lorentzian_like_5d": 24,
}


def collect_observables() -> pd.DataFrame:
    """Generate posets and collect all observables."""
    rows = []
    for n in N_VALUES:
        print(f"\n  N={n}", end="", flush=True)
        for dim_label, gen in GENERATORS.items():
            fam_name = f"lorentzian_like_{dim_label}"
            exact_thr = FAMILY_EXACT_THRESHOLDS[fam_name]

            for sid in range(SAMPLES):
                seed = SEED_BASE + 1000 * n + sid
                poset = gen(n=n, seed=seed)

                log_h, method = estimate_entropy_by_family(
                    poset, family=fam_name, sis_runs=SIS_RUNS, seed=seed,
                    default_exact_threshold=24,
                    family_exact_thresholds={fam_name: exact_thr},
                )

                n_links = hasse_links(poset)
                ic = interval_counts_fast(poset)

                # Geometric consistency (for unification test — skip at large N for speed)
                if n <= 80:
                    try:
                        geo = geometric_components(poset)
                        p_consistency = geo.get("consistency", 0.0)
                    except Exception:
                        p_consistency = np.nan
                else:
                    p_consistency = np.nan

                rows.append({
                    "n": n, "dim": dim_label,
                    "sample_id": sid, "seed": seed,
                    "log_H": float(log_h), "method": method,
                    "C0_links": n_links,
                    "C1": ic["C1"], "C2": ic["C2"], "C3": ic["C3"],
                    "S_link_d2": n - 2 * n_links,
                    "S_link_d2_norm": (n - 2 * n_links) / n,
                    "link_density": n_links / (n * (n - 1) / 2),
                    "p_consistency": p_consistency,
                })
            print(f" [{dim_label}]", end="", flush=True)
        print()
    return pd.DataFrame(rows)


def compute_winners(df: pd.DataFrame) -> pd.DataFrame:
    """Compute winners at each (N, λ)."""
    agg = df.groupby(["n", "dim"]).agg(
        log_H_mean=("log_H", "mean"),
        S_link_norm_mean=("S_link_d2_norm", "mean"),
    ).reset_index()

    rows = []
    for n in sorted(agg["n"].unique()):
        for lam in LAMBDA_VALUES:
            sub = agg[agg["n"] == n].copy()
            sub["score"] = -BETA * sub["log_H_mean"] + lam * sub["S_link_norm_mean"]
            sub = sub.sort_values("score")
            winner = sub.iloc[0]
            runner_up = sub.iloc[1]
            margin = runner_up["score"] - winner["score"]
            rows.append({
                "n": n, "lambda": lam,
                "winner": winner["dim"],
                "runner_up": runner_up["dim"],
                "margin": margin,
            })
    return pd.DataFrame(rows)


def compute_xi(df: pd.DataFrame) -> pd.DataFrame:
    """Compute Ξ between adjacent dimension pairs."""
    agg = df.groupby(["n", "dim"]).agg(
        log_H_mean=("log_H", "mean"),
        S_link_norm_mean=("S_link_d2_norm", "mean"),
        C0_links_mean=("C0_links", "mean"),
    ).reset_index()

    dim_order = ["2d", "3d", "4d", "5d"]
    dim_num = {"2d": 2, "3d": 3, "4d": 4, "5d": 5}

    rows = []
    for n in sorted(agg["n"].unique()):
        sub = agg[agg["n"] == n].set_index("dim")
        for i in range(len(dim_order) - 1):
            d_lo, d_hi = dim_order[i], dim_order[i + 1]
            if d_lo not in sub.index or d_hi not in sub.index:
                continue
            s_lo = sub.loc[d_lo, "S_link_norm_mean"]
            s_hi = sub.loc[d_hi, "S_link_norm_mean"]
            h_lo = sub.loc[d_lo, "log_H_mean"] / n
            h_hi = sub.loc[d_hi, "log_H_mean"] / n
            delta_S = s_hi - s_lo
            delta_H = h_hi - h_lo
            xi = abs(delta_S) / abs(delta_H) if abs(delta_H) > 1e-12 else np.nan
            rows.append({
                "n": n,
                "dim_pair": f"{dim_num[d_lo]}→{dim_num[d_hi]}",
                "Xi": xi,
                "delta_S_link_norm": delta_S,
                "delta_logH_per_N": delta_H,
            })
    return pd.DataFrame(rows)


def unification_test(df: pd.DataFrame) -> str:
    """Test correlation between link density and consistency penalty."""
    lines = []
    lines.append("\n" + "=" * 60)
    lines.append("UNIFICATION TEST: Link Density vs Consistency Penalty")
    lines.append("=" * 60)

    # For each N, compute correlation across dimensions
    agg = df.groupby(["n", "dim"]).agg(
        link_density_mean=("link_density", "mean"),
        C0_links_mean=("C0_links", "mean"),
        p_consistency_mean=("p_consistency", "mean"),
        log_H_mean=("log_H", "mean"),
    ).reset_index()

    lines.append(f"\n  {'N':>4}  {'r(link_dens, consistency)':>28}  {'p-value':>10}  {'interpretation':>20}")
    lines.append("  " + "-" * 70)

    for n in sorted(agg["n"].unique()):
        sub = agg[agg["n"] == n].dropna(subset=["p_consistency_mean"])
        if len(sub) < 3:
            continue
        x = sub["link_density_mean"].values
        y = sub["p_consistency_mean"].values
        if np.std(x) < 1e-12 or np.std(y) < 1e-12:
            continue

        from scipy.stats import pearsonr
        r, p = pearsonr(x, y)
        interp = "STRONG" if abs(r) > 0.9 else "moderate" if abs(r) > 0.7 else "weak"
        lines.append(f"  {n:>4}  {r:>28.4f}  {p:>10.4f}  {interp:>20}")

    # Also show the raw values
    lines.append("\n  Detailed data (mean over samples):")
    lines.append(f"  {'N':>4}  {'dim':>4}  {'link_dens':>10}  {'C0/N':>8}  {'consistency':>12}  {'logH/N':>8}")
    for _, row in agg.sort_values(["n", "dim"]).iterrows():
        lines.append(f"  {int(row['n']):>4}  {row['dim']:>4}  {row['link_density_mean']:>10.4f}  {row['C0_links_mean']/row['n']:>8.4f}  {row['p_consistency_mean']:>12.4f}  {row['log_H_mean']/row['n']:>8.4f}")

    return "\n".join(lines)


def main():
    print("=" * 72)
    print("Prediction A — Large-N Finite-Size Scaling + Ξ Closure")
    print("=" * 72)

    # Step 1: Collect observables
    print("\n[1/4] Collecting observables...")
    df = collect_observables()
    df.to_csv(OUT_DIR / "raw_observables_large_n.csv", index=False)
    print(f"\n  Saved {len(df)} observations")

    # Step 2: Winners
    print("\n[2/4] Computing winners...")
    winners = compute_winners(df)
    winners.to_csv(OUT_DIR / "winners_large_n.csv", index=False)

    print("\n  Winner table (λ=7):")
    w7 = winners[winners["lambda"] == 7.0].sort_values("n")
    print(f"  {'N':>4}  {'Winner':>8}  {'Runner-up':>10}  {'Margin':>8}")
    count_4d = 0
    for _, r in w7.iterrows():
        is4d = "★" if r["winner"] == "4d" else " "
        if r["winner"] == "4d":
            count_4d += 1
        print(f"  {int(r['n']):>4}  {r['winner']:>8}  {r['runner_up']:>10}  {r['margin']:>8.3f}  {is4d}")
    print(f"\n  4D wins: {count_4d}/{len(w7)}")

    # Print all λ values
    print("\n  Winner summary across λ:")
    print(f"  {'λ':>4}", end="")
    for n in sorted(winners["n"].unique()):
        print(f"  {f'N={int(n)}':>8}", end="")
    print(f"  {'4D count':>10}")
    for lam in LAMBDA_VALUES:
        wl = winners[winners["lambda"] == lam].sort_values("n")
        print(f"  {lam:>4.0f}", end="")
        c4d = 0
        for _, r in wl.iterrows():
            w = r["winner"]
            if w == "4d":
                c4d += 1
            print(f"  {w:>8}", end="")
        star = " ★" if c4d == len(wl) else ""
        print(f"  {c4d}/{len(wl):>8}{star}")

    # Step 3: Ξ analysis
    print("\n[3/4] Computing Ξ...")
    xi_df = compute_xi(df)
    xi_df.to_csv(OUT_DIR / "xi_large_n.csv", index=False)

    print("\n  Ξ₄→₅ across N:")
    xi45 = xi_df[xi_df["dim_pair"] == "4→5"].sort_values("n")
    print(f"  {'N':>4}  {'Ξ₄→₅':>8}  {'ΔS_link/N':>10}  {'Δ(logH/N)':>10}")
    for _, r in xi45.iterrows():
        print(f"  {int(r['n']):>4}  {r['Xi']:>8.4f}  {r['delta_S_link_norm']:>10.4f}  {r['delta_logH_per_N']:>10.4f}")
    print(f"\n  Median Ξ₄→₅: {xi45['Xi'].median():.2f}")
    print(f"  Std Ξ₄→₅: {xi45['Xi'].std():.2f}")
    print(f"  CV: {xi45['Xi'].std()/xi45['Xi'].mean()*100:.1f}%")

    # Check large-N stability
    small_n_xi = xi45[xi45["n"] <= 68]["Xi"]
    large_n_xi = xi45[xi45["n"] > 68]["Xi"]
    if len(large_n_xi) > 0:
        print(f"\n  Small N (≤68) median: {small_n_xi.median():.2f}")
        print(f"  Large N (>68) median: {large_n_xi.median():.2f}")
        drift = abs(large_n_xi.median() - small_n_xi.median()) / small_n_xi.median() * 100
        print(f"  Drift: {drift:.1f}%")
        if drift < 20:
            print("  → Ξ₄→₅ is STABLE at large N ✓")
        else:
            print("  → Ξ₄→₅ shows drift at large N — needs investigation")

    # Margin scaling
    print("\n  4D margin scaling (λ=7):")
    margins = winners[(winners["lambda"] == 7.0) & (winners["winner"] == "4d")].sort_values("n")
    if len(margins) > 2:
        ns = margins["n"].values
        ms = margins["margin"].values
        print(f"  {'N':>4}  {'Margin':>8}")
        for n, m in zip(ns, ms):
            print(f"  {int(n):>4}  {m:>8.3f}")

        # Linear fit of margin vs N
        if len(ns) >= 3:
            from numpy.polynomial import polynomial as P
            coeffs = np.polyfit(ns, ms, 1)
            print(f"\n  Linear fit: margin ≈ {coeffs[0]:.4f} * N + {coeffs[1]:.3f}")
            if coeffs[0] > 0:
                print("  → Margin GROWS with N ✓ (4D selection strengthens)")
            else:
                print("  → Margin shrinks with N ⚠ (may destabilize)")

    # Step 4: Unification test
    print("\n[4/4] Unification test...")
    unif_report = unification_test(df)
    print(unif_report)

    with open(OUT_DIR / "unification_test.txt", "w", encoding="utf-8") as f:
        f.write(unif_report)

    # Summary
    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)

    # Count 4D wins at λ=7 across all N
    w7_all = winners[winners["lambda"] == 7.0]
    n4d_all = (w7_all["winner"] == "4d").sum()
    print(f"\n  4D wins at λ=7: {n4d_all}/{len(w7_all)}")

    # New N specifically
    for n in [80, 96, 112]:
        w = w7_all[w7_all["n"] == n]
        if not w.empty:
            r = w.iloc[0]
            star = "✓" if r["winner"] == "4d" else "✗"
            print(f"  N={n}: winner={r['winner']}, margin={r['margin']:.3f} {star}")

    print(f"\n  Ξ₄→₅ overall median: {xi45['Xi'].median():.2f}")
    print(f"  Ξ₄→₅ overall CV: {xi45['Xi'].std()/xi45['Xi'].mean()*100:.1f}%")

    # Closure verdict
    all_new_n_4d = all(
        w7_all[w7_all["n"] == n]["winner"].iloc[0] == "4d"
        for n in [80, 96, 112]
        if not w7_all[w7_all["n"] == n].empty
    )
    xi_stable = True
    if len(large_n_xi) > 0:
        xi_stable = abs(large_n_xi.median() - small_n_xi.median()) / small_n_xi.median() < 0.2

    print(f"\n  === CLOSURE ASSESSMENT ===")
    print(f"  4D wins all new N at λ=7: {'YES ✓' if all_new_n_4d else 'NO ✗'}")
    print(f"  Ξ₄→₅ stable at large N: {'YES ✓' if xi_stable else 'NO ✗'}")
    if all_new_n_4d and xi_stable:
        print(f"\n  ★★★ PREDICTION A CLOSURE CRITERIA MET ★★★")
    else:
        print(f"\n  Closure criteria NOT fully met — see details above")


if __name__ == "__main__":
    main()
