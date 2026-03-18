CQG submission package (Unified B + A)
Date: 2026-03-18

Entry point:
- manuscript_cqg.tex

Figures:
- manuscript_figures/fig1_gamma_c_curve.png
- manuscript_figures/fig2_ablation_summary.png
- manuscript_figures/fig3_noncircular_replacement.png
- manuscript_figures/fig4_exact_timing_frontier.png
- manuscript_figures/fig5_mixed_lor2d_vs_kr.png
- preA/manuscript_figures/fig1_margin_of_victory.png
- preA/manuscript_figures/fig2_winner_phase_comparison.png

Build (MiKTeX):
- latexmk -pdf -interaction=nonstopmode -halt-on-error manuscript_cqg.tex

Notes:
- manuscript uses a generic LaTeX article class. If CQG requests the IOP class (iopart),
  we can switch after acceptance or on request (content should stay stable).
