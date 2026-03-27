"""
Margin vs N Visualization for Letter Figure
============================================

Reads the turn-on experiment results and produces publication-quality plots:
  1. Mean margin ± std vs N  (bar chart + error bars)
  2. #1 rate vs N  (step diagram)
  3. Reference manifold μ(N) components

Output: outputs_carlip/fig_margin_vs_n.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# ── Data from turn-on experiment (2026-03-27) ────────────────────────────
N_GRID    = np.array([12,   14,   16,   18,   20,   24,   28,   32])
RATE_1    = np.array([6,    10,   10,   10,   10,   10,   10,   10])  # out of 10
MEAN_MRG  = np.array([0.1041, 1.2775, 1.7242, 2.2009, 2.0535, 3.1672, 4.5874, 3.9896])
STD_MRG   = np.array([0.3567, 0.4335, 0.7373, 0.6575, 0.9341, 0.9963, 1.0272, 1.6748])
MIN_MRG   = np.array([-0.4668, 0.3349, 0.4207, 1.4299, 1.1691, 1.5838, 3.0851, 1.3707])
MAX_MRG   = np.array([0.6508, 1.8422, 2.4245, 3.6731, 4.4253, 4.4035, 6.0192, 7.5554])

# Reference manifold means
MU_D = np.array([3.9556, 3.9515, 3.9582, 3.9501, 3.9381, 3.9592, 3.9468, 3.9581])
MU_C = np.array([0.0691, 0.0820, 0.0885, 0.0985, 0.1136, 0.1256, 0.1396, 0.1522])
MU_W = np.array([0.6049, 0.5736, 0.5558, 0.5367, 0.5196, 0.4959, 0.4723, 0.4618])

N_SEEDS = 10

# ── Figure Setup ─────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "legend.fontsize": 10,
    "figure.dpi": 150,
})

fig, axes = plt.subplots(2, 2, figsize=(10, 8))
fig.suptitle(r"$S_{\mathrm{MD}}$ Turn-On Diagnostics  (REPS=80, 25 families, 10 seeds)",
             fontsize=14, fontweight="bold")

# ── Panel (a): Margin vs N ──────────────────────────────────────────────
ax = axes[0, 0]
ax.errorbar(N_GRID, MEAN_MRG, yerr=STD_MRG, fmt="o-", color="#2563EB",
            capsize=4, capthick=1.5, linewidth=1.5, markersize=6, label="Mean ± std")
ax.fill_between(N_GRID, MIN_MRG, MAX_MRG, alpha=0.12, color="#2563EB", label="Min–Max range")
ax.axhline(0, color="r", linestyle="--", linewidth=0.8, alpha=0.6)
ax.axvline(14, color="green", linestyle=":", linewidth=1.2, alpha=0.7, label="Turn-on: N=14")
ax.set_xlabel("N (sprinkled points)")
ax.set_ylabel("Margin (runner-up − Lor4D)")
ax.set_title("(a) Separation margin vs N")
ax.legend(loc="upper left", fontsize=9)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

# ── Panel (b): #1 rate vs N ────────────────────────────────────────────
ax = axes[0, 1]
pct = 100.0 * RATE_1 / N_SEEDS
colors = ["#EF4444" if r < N_SEEDS else "#22C55E" for r in RATE_1]
bars = ax.bar(N_GRID, pct, width=1.5, color=colors, edgecolor="white", linewidth=0.5)
ax.axhline(100, color="green", linestyle="--", linewidth=0.8, alpha=0.5)
ax.set_xlabel("N")
ax.set_ylabel("Lor4D rank #1 rate (%)")
ax.set_title("(b) Identification success rate")
ax.set_ylim(0, 115)
for bar, rate in zip(bars, RATE_1):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
            f"{rate}/{N_SEEDS}", ha="center", va="bottom", fontsize=9)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

# ── Panel (c): Min margin (worst case) ─────────────────────────────────
ax = axes[1, 0]
colors_min = ["#EF4444" if m < 0 else "#2563EB" for m in MIN_MRG]
ax.bar(N_GRID, MIN_MRG, width=1.5, color=colors_min, edgecolor="white", linewidth=0.5)
ax.axhline(0, color="r", linestyle="--", linewidth=0.8, alpha=0.6)
ax.axvline(14, color="green", linestyle=":", linewidth=1.2, alpha=0.7)
ax.set_xlabel("N")
ax.set_ylabel("Minimum margin across seeds")
ax.set_title("(c) Worst-case margin (conservative bound)")
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

# ── Panel (d): Reference manifold μ(N) ─────────────────────────────────
ax = axes[1, 1]
ax.plot(N_GRID, MU_D, "s-", color="#8B5CF6", label=r"$\mu_{d_\mathrm{eff}}$", markersize=5)
ax.set_ylabel(r"$\mu_{d_\mathrm{eff}}$", color="#8B5CF6")
ax.tick_params(axis="y", labelcolor="#8B5CF6")
ax.axhline(4.0, color="#8B5CF6", linestyle=":", linewidth=0.8, alpha=0.4)
ax.set_ylim(3.90, 4.05)

ax2 = ax.twinx()
ax2.plot(N_GRID, MU_C, "^-", color="#F59E0B", label=r"$\mu_{C_1/C_0}$", markersize=5)
ax2.plot(N_GRID, MU_W, "D-", color="#10B981", label=r"$\mu_{w/N}$", markersize=5)
ax2.set_ylabel(r"$\mu_{C_1/C_0}$, $\mu_{w/N}$")

ax.set_xlabel("N")
ax.set_title("(d) Reference manifold μ(N)")

# Combined legend
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, loc="center right", fontsize=9)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

plt.tight_layout()

outdir = Path(__file__).parent / "outputs_carlip"
outdir.mkdir(exist_ok=True)
outpath = outdir / "fig_margin_vs_n.png"
fig.savefig(outpath, dpi=200, bbox_inches="tight")
print(f"Saved: {outpath}")
plt.close()
