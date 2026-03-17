"""
Prediction A — Dimensionless Control Parameter Ξ Analysis

Computes Ξ_d(N, geom) ≈ Δ(link penalty per element) / Δ(entropy per element)
across all three generator types (original cube, independent-seed cube,
causal diamond) to test whether 4D wins when Ξ crosses a
generator-independent threshold.

The hypothesis (from GPT review):
  "If the 4D selection window exists wherever a normalized
   connectivity-vs-entropy ratio crosses a universal threshold,
   then the mechanism is truly geometry-independent."

Definition:
  Ξ_{d→d+1}(N, geom) = |S_link(d)/N - S_link(d+1)/N| / |log H(d)/N - log H(d+1)/N|
                       = |ΔS_link/N| / |Δlog_H/N|

  This measures how much link-penalty separation you get per unit of
  entropy separation between adjacent dimensions.
  Large Ξ → link action strongly differentiates (favors lower dim)
  Small Ξ → entropy dominates (favors higher dim)
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd

OUT_DIR = Path("outputs_exploratory/prediction_a_xi_parameter")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def load_all_data() -> pd.DataFrame:
    """Load the raw CSV data from all three generator runs."""
    frames = []

    # 1. Original cube (from bd_extended output — seed 980000 family)
    p1 = Path("outputs_exploratory/prediction_a_bd_extended/base_observables.csv")
    if p1.exists():
        df1 = pd.read_csv(p1)
        # Filter to Lorentzian families only, extract dim label
        lor_mask = df1["family"].str.startswith("lorentzian_like_")
        df1 = df1[lor_mask].copy()
        df1["suite"] = "cube_original"
        df1["dim"] = df1["family"].str.replace("lorentzian_like_", "")
        # Harmonize column names
        df1.rename(columns={"S_BD_d2": "S_link_d2", "S_BD_d2_norm": "S_link_d2_norm"}, inplace=True)
        frames.append(df1)
        print(f"  Loaded {len(df1)} rows from cube_original (bd_extended)")

    # 2. Independent-seed cube + diamond (from generator_robustness output)
    p2 = Path("outputs_exploratory/prediction_a_generator_robustness/base_observables_robustness.csv")
    if p2.exists():
        df2 = pd.read_csv(p2)
        frames.append(df2)
        for s in df2["suite"].unique():
            cnt = (df2["suite"] == s).sum()
            print(f"  Loaded {cnt} rows from {s}")

    if not frames:
        raise FileNotFoundError("No raw observable CSVs found!")

    combined = pd.concat(frames, ignore_index=True)
    print(f"\n  Total: {len(combined)} observations across {combined['suite'].nunique()} suites")
    return combined


def compute_xi(df: pd.DataFrame) -> pd.DataFrame:
    """Compute Ξ between adjacent dimension pairs for each (suite, N)."""
    # Average over samples
    agg = df.groupby(["suite", "n", "dim"]).agg(
        log_H_mean=("log_H", "mean"),
        S_link_norm_mean=("S_link_d2_norm", "mean"),
        C0_links_mean=("C0_links", "mean"),
    ).reset_index()

    dim_order = ["2d", "3d", "4d", "5d"]
    dim_num = {"2d": 2, "3d": 3, "4d": 4, "5d": 5}

    rows = []
    for suite in agg["suite"].unique():
        for n in sorted(agg["n"].unique()):
            sub = agg[(agg["suite"] == suite) & (agg["n"] == n)].set_index("dim")
            for i in range(len(dim_order) - 1):
                d_lo, d_hi = dim_order[i], dim_order[i + 1]
                if d_lo not in sub.index or d_hi not in sub.index:
                    continue

                s_lo = sub.loc[d_lo, "S_link_norm_mean"]
                s_hi = sub.loc[d_hi, "S_link_norm_mean"]
                h_lo = sub.loc[d_lo, "log_H_mean"] / n
                h_hi = sub.loc[d_hi, "log_H_mean"] / n

                # Link penalty per element: S_link/N is negative; higher dim = less negative
                # ΔS = S_hi - S_lo > 0 (higher dim has less penalty = higher S)
                delta_S = s_hi - s_lo  # positive: higher dim less penalized
                delta_H = h_hi - h_lo  # positive: higher dim has more entropy

                # Ξ = |ΔS_link/N| / |Δ(logH/N)|
                if abs(delta_H) < 1e-12:
                    xi = np.nan
                else:
                    xi = abs(delta_S) / abs(delta_H)

                rows.append({
                    "suite": suite,
                    "n": n,
                    "dim_pair": f"{dim_num[d_lo]}→{dim_num[d_hi]}",
                    "d_lo": dim_num[d_lo],
                    "d_hi": dim_num[d_hi],
                    "S_link_lo": s_lo,
                    "S_link_hi": s_hi,
                    "delta_S_link_norm": delta_S,
                    "logH_per_N_lo": h_lo,
                    "logH_per_N_hi": h_hi,
                    "delta_logH_per_N": delta_H,
                    "Xi": xi,
                    "C0_lo": sub.loc[d_lo, "C0_links_mean"],
                    "C0_hi": sub.loc[d_hi, "C0_links_mean"],
                })

    return pd.DataFrame(rows)


def analyze_and_report(xi_df: pd.DataFrame) -> str:
    """Analyze Ξ patterns and generate a text report."""
    lines = []
    lines.append("=" * 72)
    lines.append("Ξ (Xi) Dimensionless Control Parameter Analysis")
    lines.append("=" * 72)

    # Focus on 4→5 transition (the critical boundary)
    lines.append("\n\n### 4→5 Transition (Key Boundary)")
    lines.append("-" * 60)
    t45 = xi_df[xi_df["dim_pair"] == "4→5"].copy()
    for suite in sorted(t45["suite"].unique()):
        sub = t45[t45["suite"] == suite].sort_values("n")
        lines.append(f"\n  Suite: {suite}")
        lines.append(f"  {'N':>4}  {'Ξ_4→5':>8}  {'ΔS_link/N':>10}  {'Δ(logH/N)':>10}  {'C0_4D':>6}  {'C0_5D':>6}")
        for _, r in sub.iterrows():
            lines.append(f"  {int(r['n']):>4}  {r['Xi']:>8.4f}  {r['delta_S_link_norm']:>10.4f}  {r['delta_logH_per_N']:>10.4f}  {r['C0_lo']:>6.1f}  {r['C0_hi']:>6.1f}")

    # Focus on 3→4 transition
    lines.append("\n\n### 3→4 Transition")
    lines.append("-" * 60)
    t34 = xi_df[xi_df["dim_pair"] == "3→4"].copy()
    for suite in sorted(t34["suite"].unique()):
        sub = t34[t34["suite"] == suite].sort_values("n")
        lines.append(f"\n  Suite: {suite}")
        lines.append(f"  {'N':>4}  {'Ξ_3→4':>8}  {'ΔS_link/N':>10}  {'Δ(logH/N)':>10}")
        for _, r in sub.iterrows():
            lines.append(f"  {int(r['n']):>4}  {r['Xi']:>8.4f}  {r['delta_S_link_norm']:>10.4f}  {r['delta_logH_per_N']:>10.4f}")

    # All transitions side by side for each suite
    lines.append("\n\n### All Transitions — Ξ Summary Table")
    lines.append("-" * 60)
    for suite in sorted(xi_df["suite"].unique()):
        sub = xi_df[xi_df["suite"] == suite]
        lines.append(f"\n  Suite: {suite}")
        pivot = sub.pivot_table(index="n", columns="dim_pair", values="Xi", aggfunc="first")
        lines.append(f"  {'N':>4}" + "".join(f"  {col:>8}" for col in pivot.columns))
        for n in sorted(pivot.index):
            vals = "".join(f"  {pivot.loc[n, c]:>8.4f}" if not np.isnan(pivot.loc[n, c]) else f"  {'NaN':>8}" for c in pivot.columns)
            lines.append(f"  {int(n):>4}{vals}")

    # Cross-generator comparison: Ξ_4→5 at same N
    lines.append("\n\n### Cross-Generator Ξ₄→₅ Comparison")
    lines.append("-" * 60)
    t45_pivot = t45.pivot_table(index="n", columns="suite", values="Xi", aggfunc="first")
    if len(t45_pivot.columns) > 1:
        lines.append(f"  {'N':>4}" + "".join(f"  {s:>16}" for s in t45_pivot.columns))
        for n in sorted(t45_pivot.index):
            vals = "".join(f"  {t45_pivot.loc[n, s]:>16.4f}" if not np.isnan(t45_pivot.loc[n, s]) else f"  {'NaN':>16}" for s in t45_pivot.columns)
            lines.append(f"  {int(n):>4}{vals}")

    # Key insight: the relevant comparison for 4D selection
    # At λ where 4D wins, what is Ξ_4→5 relative to Ξ_3→4?
    lines.append("\n\n### Ξ Ratio: Ξ₄→₅ / Ξ₃→₄ (Selection Asymmetry)")
    lines.append("-" * 60)
    lines.append("  If Ξ₄→₅ / Ξ₃→₄ > 1, then 4→5 link penalty gap is relatively larger")
    lines.append("  per unit entropy than 3→4, which helps 4D win.")
    for suite in sorted(xi_df["suite"].unique()):
        sub34 = xi_df[(xi_df["suite"] == suite) & (xi_df["dim_pair"] == "3→4")].set_index("n")
        sub45 = xi_df[(xi_df["suite"] == suite) & (xi_df["dim_pair"] == "4→5")].set_index("n")
        common_n = sorted(set(sub34.index) & set(sub45.index))
        lines.append(f"\n  Suite: {suite}")
        lines.append(f"  {'N':>4}  {'Ξ₃→₄':>8}  {'Ξ₄→₅':>8}  {'Ratio':>8}")
        for n in common_n:
            xi34 = sub34.loc[n, "Xi"]
            xi45 = sub45.loc[n, "Xi"]
            ratio = xi45 / xi34 if xi34 > 1e-12 else np.nan
            lines.append(f"  {int(n):>4}  {xi34:>8.4f}  {xi45:>8.4f}  {ratio:>8.4f}")

    # Summary hypothesis test
    lines.append("\n\n### Summary: Is there a universal Ξ threshold?")
    lines.append("-" * 60)
    medians = t45.groupby("suite")["Xi"].median()
    lines.append("  Median Ξ₄→₅ across generators:")
    for suite, med in medians.items():
        lines.append(f"    {suite}: {med:.4f}")
    if len(medians) > 1:
        spread = medians.max() - medians.min()
        mean_val = medians.mean()
        lines.append(f"\n  Spread: {spread:.4f} (CV = {spread/mean_val*100:.1f}%)")
        if spread / mean_val < 0.3:
            lines.append("  → Ξ₄→₅ is CONSISTENT across generators (CV < 30%)")
            lines.append("  → Supports universal mechanism hypothesis!")
        else:
            lines.append("  → Ξ₄→₅ shows significant variation across generators")
            lines.append("  → Mechanism is the same but quantitative threshold is geometry-dependent")

    return "\n".join(lines)


def main():
    print("Loading data from all generator suites...")
    df = load_all_data()

    # Standardize column names
    if "C0_links" not in df.columns and "C0" in df.columns:
        df.rename(columns={"C0": "C0_links"}, inplace=True)

    # Need S_link_d2_norm
    if "S_link_d2_norm" not in df.columns:
        if "S_link_d2" in df.columns:
            df["S_link_d2_norm"] = df["S_link_d2"] / df["n"]
        elif "C0_links" in df.columns:
            df["S_link_d2"] = df["n"] - 2 * df["C0_links"]
            df["S_link_d2_norm"] = df["S_link_d2"] / df["n"]

    print("\nComputing Ξ for all dimension pairs...")
    xi_df = compute_xi(df)
    xi_df.to_csv(OUT_DIR / "xi_all_transitions.csv", index=False)
    print(f"  Saved {len(xi_df)} Ξ values to xi_all_transitions.csv")

    report = analyze_and_report(xi_df)
    print("\n" + report)

    with open(OUT_DIR / "xi_analysis_report.txt", "w", encoding="utf-8") as f:
        f.write(report)
    print(f"\n  Report saved to {OUT_DIR / 'xi_analysis_report.txt'}")


if __name__ == "__main__":
    main()
