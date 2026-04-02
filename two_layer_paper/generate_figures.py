"""
Generate publication-quality figures for the two-layer screening paper.

Figures:
  Fig 1: S_MD margin vs N (turn-on + basin deepening)
  Fig 2: Basin deepening diagnostics (V_eff, Fisher, det(Sigma))
  Fig 3: Gap decomposition (Euclidean vs Mahalanobis)
  Fig 4: S_BD vs S_MD rank discordance at N=128
  Fig 5: Feature ablation bar chart

Run from: d:\Kiro\理论体系\poset_phase\two_layer_paper\
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'font.family': 'serif',
})

outdir = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(outdir, exist_ok=True)

# ============================================================
# Fig 1: S_MD margin vs N (turn-on + basin deepening)
# ============================================================
def fig1_margin_vs_n():
    # F2 onset data (20 seeds x 120 reps)
    N_f2 = [10, 12, 14, 16, 18, 20, 22, 24]
    margin_f2 = [0.308, 1.161, 1.487, 1.771, 2.009, 2.615, 2.988, 3.140]
    ci_lo = [0.268, 1.034, 1.280, 1.628, 1.735, 2.342, 2.731, 2.787]
    ci_hi = [0.348, 1.287, 1.694, 1.914, 2.282, 2.888, 3.245, 3.494]

    # Basin deepening + extreme N
    N_basin = [12, 16, 20, 28, 48, 64, 96, 128, 256]
    gap_basin = [-0.816, 1.888, 2.462, 4.604, 14.833, 17.467, 32.707, 37.227, 94.052]

    # Extreme N
    N_extreme = [128, 256, 512, 1024]
    gap_extreme = [93.0, 243.3, 5689.6, 192998805.0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    # Left: linear scale, onset region
    err_lo = [m - l for m, l in zip(margin_f2, ci_lo)]
    err_hi = [h - m for m, h in zip(margin_f2, ci_hi)]
    ax1.errorbar(N_f2, margin_f2, yerr=[err_lo, err_hi],
                 fmt='o-', color='#2166ac', capsize=3, markersize=5, linewidth=1.5,
                 label='F2 protocol (20 seeds)')
    ax1.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
    ax1.axvline(x=10, color='#b2182b', linestyle=':', linewidth=1.2, label='$N_{\\mathrm{id}}=10$')
    ax1.set_xlabel('$N$')
    ax1.set_ylabel('Mean Mahalanobis margin')
    ax1.set_title('(a) Identity turn-on')
    ax1.legend(loc='upper left')
    ax1.set_xlim(8, 26)

    # Right: log-log scale, full range
    N_all = N_basin + [512, 1024]
    gap_all = gap_basin + [5689.6, 192998805.0]
    # filter positive only
    mask = [g > 0 for g in gap_all]
    N_pos = [n for n, m in zip(N_all, mask) if m]
    gap_pos = [g for g, m in zip(gap_all, mask) if m]

    ax2.loglog(N_pos, gap_pos, 's-', color='#d6604d', markersize=5, linewidth=1.5)
    ax2.set_xlabel('$N$')
    ax2.set_ylabel('Mahalanobis gap $\\Delta_{\\mathrm{hist}}$')
    ax2.set_title('(b) Basin deepening (log-log)')
    ax2.annotate('$1.93\\times 10^8$', xy=(1024, 1.93e8), fontsize=9,
                 ha='right', va='bottom', color='#d6604d')

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'fig01_margin_and_deepening.png'))
    plt.savefig(os.path.join(outdir, 'fig01_margin_and_deepening.pdf'))
    plt.close()
    print('Fig 1 done.')


# ============================================================
# Fig 2: Basin deepening diagnostics (3-panel)
# ============================================================
def fig2_basin_diagnostics():
    N = np.array([12, 16, 20, 28, 36, 48, 64, 96, 128, 192, 256])
    V_eff = np.array([3.321e-3, 1.582e-3, 2.393e-3, 1.075e-3, 8.521e-4,
                      2.839e-4, 2.366e-4, 1.312e-4, 9.719e-5, 6.290e-5, 2.126e-5])
    Fisher = np.array([197, 469, 368, 540, 659, 1430, 1589, 1673, 2355, 3315, 6533])
    gap = np.array([-0.816, 1.888, 2.462, 4.604, 6.171, 14.833, 17.467,
                    32.707, 37.227, 35.067, 94.052])

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))

    # V_eff
    axes[0].loglog(N, V_eff, 'o-', color='#4393c3', markersize=4)
    # fit line
    mask = N >= 16
    p = np.polyfit(np.log(N[mask]), np.log(V_eff[mask]), 1)
    Nfit = np.linspace(16, 256, 100)
    axes[0].loglog(Nfit, np.exp(p[1]) * Nfit**p[0], '--', color='gray', linewidth=0.8)
    axes[0].set_xlabel('$N$')
    axes[0].set_ylabel('$V_{\\mathrm{eff}}$')
    axes[0].set_title(f'(a) $V_{{\\mathrm{{eff}}}} \\propto N^{{{p[0]:.2f}}}$')

    # Fisher
    axes[1].loglog(N, Fisher, 's-', color='#d6604d', markersize=4)
    p2 = np.polyfit(np.log(N[mask]), np.log(Fisher[mask]), 1)
    axes[1].loglog(Nfit, np.exp(p2[1]) * Nfit**p2[0], '--', color='gray', linewidth=0.8)
    axes[1].set_xlabel('$N$')
    axes[1].set_ylabel('$I_F = \\mathrm{tr}(\\Sigma^{-1})$')
    axes[1].set_title(f'(b) $I_F \\propto N^{{+{p2[0]:.2f}}}$')

    # Gap
    mask_pos = gap > 0
    axes[2].loglog(N[mask_pos], gap[mask_pos], 'D-', color='#762a83', markersize=4)
    mask3 = mask_pos & mask
    p3 = np.polyfit(np.log(N[mask3]), np.log(gap[mask3]), 1)
    axes[2].loglog(Nfit, np.exp(p3[1]) * Nfit**p3[0], '--', color='gray', linewidth=0.8)
    axes[2].set_xlabel('$N$')
    axes[2].set_ylabel('$\\Delta_{\\mathrm{hist}}$')
    axes[2].set_title(f'(c) Gap $\\propto N^{{+{p3[0]:.2f}}}$')

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'fig02_basin_diagnostics.png'))
    plt.savefig(os.path.join(outdir, 'fig02_basin_diagnostics.pdf'))
    plt.close()
    print('Fig 2 done.')


# ============================================================
# Fig 3: Gap decomposition (Euclidean vs Mahalanobis)
# ============================================================
def fig3_gap_decomposition():
    N = [16, 28, 64, 128, 256]
    d_E = [0.511, 0.369, 0.488, 0.459, 0.452]
    gap = [3.65, 4.47, 20.1, 40.5, 92.6]
    sigma_A = [14.0, 32.8, 84.3, 191.9, 454.0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.5))

    # Left: Euclidean distance (flat)
    ax1.plot(N, d_E, 'o-', color='#2166ac', markersize=6, linewidth=1.5)
    ax1.axhline(y=np.mean(d_E), color='gray', linestyle='--', linewidth=0.8,
                label=f'mean = {np.mean(d_E):.3f}')
    ax1.set_xlabel('$N$')
    ax1.set_ylabel('Euclidean distance $d_E$')
    ax1.set_title('(a) Physical separation (flat)')
    ax1.legend()
    ax1.set_ylim(0, 0.7)

    # Right: Precision amplification
    ax2.semilogy(N, sigma_A, 's-', color='#d6604d', markersize=6, linewidth=1.5,
                 label='$\\sigma_A = \\Delta_{\\mathrm{hist}} / d_E$')
    ax2.semilogy(N, gap, 'D-', color='#762a83', markersize=5, linewidth=1.2,
                 label='$\\Delta_{\\mathrm{hist}}$ (Mahalanobis)')
    ax2.set_xlabel('$N$')
    ax2.set_ylabel('Value')
    ax2.set_title('(b) Precision amplification')
    ax2.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'fig03_gap_decomposition.png'))
    plt.savefig(os.path.join(outdir, 'fig03_gap_decomposition.pdf'))
    plt.close()
    print('Fig 3 done.')


# ============================================================
# Fig 4: S_BD vs S_MD rank discordance at N=128
# ============================================================
def fig4_rank_discordance():
    families = ['KR_like', 'RLk6', 'KR_4layer', 'RLk6_lj', 'RLk6_tap',
                'KR_2layer', 'RLk6_mid', 'AbsLayer', 'RLk4', 'RLk8',
                'MLR', 'TransPerc', 'Lor5D', 'Lor4D', 'IntOrder', 'Lor3D', 'Lor2D']
    sbd_rank = list(range(1, 18))
    smd_rank = [5, 14, 9, 12, 11, 6, 10, 17, 13, 16, 7, 15, 2, 1, 8, 3, 4]

    colors = []
    for f in families:
        if f == 'Lor4D':
            colors.append('#d6604d')
        elif f.startswith('Lor'):
            colors.append('#4393c3')
        elif f.startswith('KR'):
            colors.append('#762a83')
        else:
            colors.append('#999999')

    fig, ax = plt.subplots(figsize=(7, 5))
    for i in range(len(families)):
        ax.annotate('', xy=(smd_rank[i], sbd_rank[i]),
                    xytext=(smd_rank[i], sbd_rank[i]))
        ax.scatter(smd_rank[i], sbd_rank[i], c=colors[i], s=60, zorder=3,
                   edgecolors='black', linewidth=0.5)
        offset = (5, 3) if families[i] != 'Lor4D' else (5, -12)
        ax.annotate(families[i], (smd_rank[i], sbd_rank[i]),
                    fontsize=7, textcoords='offset points', xytext=offset)

    ax.set_xlabel('$S_{\\mathrm{MD}}$ rank (identity)')
    ax.set_ylabel('$S_{\\mathrm{BD}}$ rank (admissibility)')
    ax.set_title('Rank discordance at $N=128$')
    ax.set_xlim(0, 18)
    ax.set_ylim(0, 18)
    ax.plot([0, 18], [0, 18], '--', color='lightgray', linewidth=0.8)
    ax.invert_yaxis()

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#d6604d',
               markersize=8, label='Lor4D'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#4393c3',
               markersize=8, label='Other Lorentzian'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#762a83',
               markersize=8, label='KR-type'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#999999',
               markersize=8, label='Layered/Other'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'fig04_rank_discordance.png'))
    plt.savefig(os.path.join(outdir, 'fig04_rank_discordance.pdf'))
    plt.close()
    print('Fig 4 done.')


# ============================================================
# Fig 5: Feature ablation bar chart
# ============================================================
def fig5_feature_ablation():
    labels = ['$d_{\\mathrm{eff}}$', '$C_1/C_0$', '$w/N$',
              '$d+c$', '$d+w$', '$c+w$', 'All 3']
    margins = [0.17, 7.17, 0.04, 36.2, 28.8, 8.31, 40.5]
    ranks = [2, 1, 1, 1, 1, 1, 1]
    colors_bar = ['#fee0d2' if r > 1 else '#deebf7' for r in ranks]
    colors_bar[-1] = '#2166ac'

    fig, ax = plt.subplots(figsize=(7, 3.5))
    bars = ax.bar(labels, margins, color=colors_bar, edgecolor='black', linewidth=0.5)

    for i, (m, r) in enumerate(zip(margins, ranks)):
        txt = f'{m:.1f}' if m >= 1 else f'{m:.2f}'
        if r > 1:
            txt += ' (rank 2)'
        ax.text(i, m + 1.2, txt, ha='center', fontsize=8)

    ax.set_ylabel('Mahalanobis margin at $N=128$')
    ax.set_title('Feature ablation: minimal complete basis')
    ax.set_ylim(0, 48)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'fig05_feature_ablation.png'))
    plt.savefig(os.path.join(outdir, 'fig05_feature_ablation.pdf'))
    plt.close()
    print('Fig 5 done.')


# ============================================================
if __name__ == '__main__':
    fig1_margin_vs_n()
    fig2_basin_diagnostics()
    fig3_gap_decomposition()
    fig4_rank_discordance()
    fig5_feature_ablation()
    print('All figures generated in figures/')
