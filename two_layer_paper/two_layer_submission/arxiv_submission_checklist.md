# arXiv Submission Checklist (Two-Layer Screening)

## A. Files
- [ ] `two_layer_screening.tex` (main source)
- [ ] `two_layer_screening.pdf` (reference PDF)
- [ ] `figures/fig01_margin_and_deepening.pdf`
- [ ] `figures/fig02_basin_diagnostics.pdf`
- [ ] `figures/fig03_gap_decomposition.pdf`
- [ ] `figures/fig04_rank_discordance.pdf`
- [ ] `figures/fig05_feature_ablation.pdf`
- [ ] bibliography contained in source / no missing refs

## B. Build sanity
- [ ] clean compile in submission folder (2-pass xelatex)
- [ ] no missing-file errors
- [ ] cross-references resolved
- [ ] output page count confirmed

## C. Metadata
- [ ] title matches manuscript
- [ ] author list verified
- [ ] abstract pasted as plain text for arXiv form
- [ ] primary category selected (proposed: gr-qc)
- [ ] comments field updated (pages/figures/version)

## D. Scientific consistency
- [ ] conclusion wording consistent with Report final state
- [ ] κ=1.0 statement uses final finite-size interpretation
- [ ] §6.11 μ(N,κ) prototype language marked as prototype (not theorem)

## E. Final release hygiene
- [ ] `README.md` updated with latest compile status
- [ ] git tag for submission snapshot (optional)
- [ ] archive package test-upload locally (optional)
