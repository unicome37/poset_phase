"""
Prediction A — BDG Action Component Contribution Figure

Creates a publication-quality stacked bar chart showing how each BDG d=4 term
contributes to the total action score across dimensions.

Key diagnostic: WHY does BDG_d4 select 5D instead of 4D?
Answer: The +9*C1 term massively boosts 2D/3D, while 5D has so few intervals
that it's barely penalized. Combined with 5D's maximum entropy → 5D always wins.

Also creates a link-action comparison panel showing the clean 4D selection.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

DATA_DIR = Path("outputs_exploratory/prediction_a_bd_extended")
OUT_DIR = Path("outputs_exploratory/prediction_a_bdg_component_figure")
OUT_DIR.mkdir(parents=True, exist_ok=True)

base = pd.read_csv(DATA_DIR / "base_observables.csv")

# Select representative N values
N_REPR = [20, 36, 52, 68]
DIM_ORDER = ["2d", "3d", "4d", "5d"]
DIM_LABELS = ["2D", "3D", "4D", "5D"]
DIM_COLORS = {"2d": "#e74c3c", "3d": "#f39c12", "4d": "#2ecc71", "5d": "#3498db"}

lor = base[base["family"].str.startswith("lorentzian_like_")].copy()
lor["dim"] = lor["family"].str.replace("lorentzian_like_", "")

# ========== Table 1: Component breakdown ==========
print("=" * 80)
print("BDG d=4 ACTION COMPONENT BREAKDOWN BY DIMENSION")
print("=" * 80)

print(f"\n{'N':>4} {'Dim':>4} | {'C0':>8} {'C1':>8} {'C2':>8} {'C3':>8} | "
      f"{'-C0':>8} {'+9C1':>8} {'-16C2':>8} {'+8C3':>8} | {'S_bdg':>8} {'S_link':>8} {'logH':>8}")
print("-" * 120)

rows_for_csv = []
for n in N_REPR:
    for dim in DIM_ORDER:
        sub = lor[(lor["n"] == n) & (lor["dim"] == dim)]
        if len(sub) == 0:
            continue
        c0 = sub["C0_links"].mean()
        c1 = sub["C1"].mean()
        c2 = sub["C2"].mean()
        c3 = sub["C3"].mean()
        log_h = sub["log_H"].mean()
        
        # BDG d=4 components
        term_c0 = -c0
        term_c1 = 9 * c1
        term_c2 = -16 * c2
        term_c3 = 8 * c3
        s_bdg = n + term_c0 + term_c1 + term_c2 + term_c3
        
        # Link action
        s_link = n - 2 * c0
        
        print(f"{n:>4} {dim:>4} | {c0:>8.1f} {c1:>8.1f} {c2:>8.1f} {c3:>8.1f} | "
              f"{term_c0:>8.1f} {term_c1:>8.1f} {term_c2:>8.1f} {term_c3:>8.1f} | "
              f"{s_bdg:>8.1f} {s_link:>8.1f} {log_h:>8.1f}")
        
        rows_for_csv.append({
            "N": n, "dim": dim,
            "C0": c0, "C1": c1, "C2": c2, "C3": c3,
            "term_minus_C0": term_c0, "term_plus_9C1": term_c1,
            "term_minus_16C2": term_c2, "term_plus_8C3": term_c3,
            "S_BDG_d4": s_bdg, "S_link_d2": s_link, "log_H": log_h,
            "S_BDG_d4_norm": s_bdg / n, "S_link_d2_norm": s_link / n,
        })
    print()

df_comp = pd.DataFrame(rows_for_csv)
df_comp.to_csv(OUT_DIR / "component_breakdown.csv", index=False)

# ========== Table 2: Score analysis at λ=6 ==========
print("\n" + "=" * 80)
print("SCORE COMPARISON AT λ=6 (the 4D selection window)")
print("=" * 80)

LAM = 6.0
BETA = 1.0

print(f"\n{'N':>4} {'Dim':>4} | {'logH':>8} {'S_link/N':>10} {'S_bdg/N':>10} | "
      f"{'score_link':>12} {'score_bdg':>12} | {'Δ_link':>10} {'Δ_bdg':>10}")
print("-" * 110)

for n in N_REPR:
    scores_link = {}
    scores_bdg = {}
    for dim in DIM_ORDER:
        row = df_comp[(df_comp["N"] == n) & (df_comp["dim"] == dim)].iloc[0]
        s_link = -BETA * row["log_H"] + LAM * row["S_link_d2_norm"]
        s_bdg = -BETA * row["log_H"] + LAM * row["S_BDG_d4_norm"]
        scores_link[dim] = s_link
        scores_bdg[dim] = s_bdg
    
    # Find 4D rank
    for dim in DIM_ORDER:
        row = df_comp[(df_comp["N"] == n) & (df_comp["dim"] == dim)].iloc[0]
        delta_link = scores_link[dim] - scores_link["4d"]
        delta_bdg = scores_bdg[dim] - scores_bdg["4d"]
        marker = " ★" if dim == "4d" else ""
        print(f"{n:>4} {dim:>4} | {row['log_H']:>8.1f} {row['S_link_d2_norm']:>10.3f} "
              f"{row['S_BDG_d4_norm']:>10.3f} | {scores_link[dim]:>12.2f} "
              f"{scores_bdg[dim]:>12.2f} | {delta_link:>10.2f} {delta_bdg:>10.2f}{marker}")
    print()

# ========== Figure generation ==========
if not HAS_MPL:
    print("\n[WARN] matplotlib not available — skipping figure generation")
    print("       Install with: pip install matplotlib")
    print(f"\nAll CSV data saved to {OUT_DIR}")
else:
    fig, axes = plt.subplots(2, len(N_REPR), figsize=(4.5 * len(N_REPR), 10),
                              gridspec_kw={"hspace": 0.35, "wspace": 0.3})
    
    TERM_COLORS = {
        "-C₀": "#e74c3c",
        "+9C₁": "#2ecc71",
        "-16C₂": "#3498db",
        "+8C₃": "#f39c12",
        "+N": "#95a5a6",
    }
    
    for col_idx, n in enumerate(N_REPR):
        sub = df_comp[df_comp["N"] == n]
        
        # --- Top row: BDG d=4 component stacked bar ---
        ax = axes[0, col_idx]
        x = np.arange(len(DIM_ORDER))
        bar_width = 0.6
        
        # Separate positive and negative components for stacked bar
        for i, dim in enumerate(DIM_ORDER):
            row = sub[sub["dim"] == dim].iloc[0]
            components = {
                "+N": float(n),
                "-C₀": row["term_minus_C0"],
                "+9C₁": row["term_plus_9C1"],
                "-16C₂": row["term_minus_16C2"],
                "+8C₃": row["term_plus_8C3"],
            }
            
            # Stack positive above 0
            bottom_pos = 0
            for label in ["+N", "+9C₁", "+8C₃"]:
                val = components[label]
                if val > 0:
                    ax.bar(i, val, bar_width, bottom=bottom_pos,
                           color=TERM_COLORS[label], edgecolor="white",
                           linewidth=0.5, label=label if i == 0 else "")
                    bottom_pos += val
            
            # Stack negative below 0
            bottom_neg = 0
            for label in ["-C₀", "-16C₂"]:
                val = components[label]
                if val < 0:
                    ax.bar(i, val, bar_width, bottom=bottom_neg,
                           color=TERM_COLORS[label], edgecolor="white",
                           linewidth=0.5, label=label if i == 0 else "")
                    bottom_neg += val
            
            # Net S_BDG marker
            net = row["S_BDG_d4"]
            ax.plot(i, net, "kD", markersize=6, zorder=5)
        
        ax.set_xticks(x)
        ax.set_xticklabels(DIM_LABELS)
        ax.set_title(f"N = {n}", fontsize=12, fontweight="bold")
        ax.axhline(0, color="black", linewidth=0.8, linestyle="-")
        ax.set_ylabel("Action component value" if col_idx == 0 else "")
        if col_idx == 0:
            ax.legend(fontsize=7, loc="upper left", ncol=2)
        
        # --- Bottom row: Score comparison at λ=6 ---
        ax2 = axes[1, col_idx]
        scores_l = []
        scores_b = []
        for dim in DIM_ORDER:
            row = sub[sub["dim"] == dim].iloc[0]
            scores_l.append(-BETA * row["log_H"] + LAM * row["S_link_d2_norm"])
            scores_b.append(-BETA * row["log_H"] + LAM * row["S_BDG_d4_norm"])
        
        x_off = np.arange(len(DIM_ORDER))
        w = 0.35
        bars1 = ax2.bar(x_off - w/2, scores_l, w, color="#2ecc71", alpha=0.8,
                        label="Link action" if col_idx == 0 else "")
        bars2 = ax2.bar(x_off + w/2, scores_b, w, color="#e74c3c", alpha=0.8,
                        label="BDG d=4" if col_idx == 0 else "")
        
        # Highlight winner
        winner_l = DIM_ORDER[np.argmin(scores_l)]
        winner_b = DIM_ORDER[np.argmin(scores_b)]
        for i, dim in enumerate(DIM_ORDER):
            if dim == winner_l:
                ax2.bar(x_off[i] - w/2, scores_l[i], w, color="#2ecc71",
                        edgecolor="black", linewidth=2)
            if dim == winner_b:
                ax2.bar(x_off[i] + w/2, scores_b[i], w, color="#e74c3c",
                        edgecolor="black", linewidth=2)
        
        ax2.set_xticks(x_off)
        ax2.set_xticklabels(DIM_LABELS)
        ax2.set_title(f"Score at λ={LAM:.0f}, N={n}", fontsize=11)
        ax2.set_ylabel("Score (lower = selected)" if col_idx == 0 else "")
        if col_idx == 0:
            ax2.legend(fontsize=8, loc="upper right")
    
    fig.suptitle("BDG d=4 Component Breakdown & Link vs BDG Score Comparison",
                 fontsize=14, fontweight="bold", y=0.98)
    
    fig.savefig(OUT_DIR / "bdg_component_breakdown.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT_DIR / "bdg_component_breakdown.pdf", bbox_inches="tight")
    print(f"\n[OK] Figure saved to {OUT_DIR}/bdg_component_breakdown.png/pdf")
    plt.close()

    # ===== Figure 2: Link action λ scan heatmap =====
    fig2, (ax_l, ax_b) = plt.subplots(1, 2, figsize=(14, 5))
    
    LAMBDA_GRID = [0, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20]
    ALL_N = sorted(lor["n"].unique())
    
    # Compute mean observables per (n, dim)
    agg = lor.groupby(["n", "dim"]).agg({
        "log_H": "mean", "C0_links": "mean", "C1": "mean", "C2": "mean", "C3": "mean"
    }).reset_index()
    
    dim_to_num = {"2d": 0, "3d": 1, "4d": 2, "5d": 3}
    
    for ax, action_name, action_func, cmap_name in [
        (ax_l, "Link action (N−2C₀)", 
         lambda r: (r["n"] - 2*r["C0_links"]) / r["n"], "Greens"),
        (ax_b, "BDG d=4 standard", 
         lambda r: (r["n"] - r["C0_links"] + 9*r["C1"] - 16*r["C2"] + 8*r["C3"]) / r["n"], "Reds"),
    ]:
        winner_matrix = np.full((len(LAMBDA_GRID), len(ALL_N)), np.nan)
        
        for li, lam in enumerate(LAMBDA_GRID):
            for ni, n in enumerate(ALL_N):
                best_score = np.inf
                best_dim = -1
                for dim in DIM_ORDER:
                    row = agg[(agg["n"] == n) & (agg["dim"] == dim)]
                    if len(row) == 0:
                        continue
                    r = row.iloc[0]
                    s_norm = action_func(r)
                    score = -BETA * r["log_H"] + lam * s_norm
                    if score < best_score:
                        best_score = score
                        best_dim = dim_to_num[dim]
                winner_matrix[li, ni] = best_dim
        
        im = ax.imshow(winner_matrix, aspect="auto", origin="lower",
                       cmap=plt.cm.get_cmap("RdYlGn", 4), vmin=-0.5, vmax=3.5)
        ax.set_xticks(range(len(ALL_N)))
        ax.set_xticklabels(ALL_N, fontsize=8)
        ax.set_yticks(range(len(LAMBDA_GRID)))
        ax.set_yticklabels(LAMBDA_GRID)
        ax.set_xlabel("N")
        ax.set_ylabel("λ")
        ax.set_title(action_name, fontsize=12, fontweight="bold")
        
        # Annotate winners
        for li in range(len(LAMBDA_GRID)):
            for ni in range(len(ALL_N)):
                w = int(winner_matrix[li, ni])
                txt = DIM_LABELS[w]
                ax.text(ni, li, txt, ha="center", va="center", fontsize=7,
                        fontweight="bold" if w == 2 else "normal",
                        color="white" if w in [0, 3] else "black")
    
    fig2.suptitle("Winner Heatmap: Link Action vs BDG d=4 Standard",
                  fontsize=14, fontweight="bold")
    fig2.tight_layout()
    fig2.savefig(OUT_DIR / "winner_heatmap_comparison.png", dpi=200, bbox_inches="tight")
    fig2.savefig(OUT_DIR / "winner_heatmap_comparison.pdf", bbox_inches="tight")
    print(f"[OK] Heatmap saved to {OUT_DIR}/winner_heatmap_comparison.png/pdf")
    plt.close()

print(f"\nAll outputs saved to {OUT_DIR}")
