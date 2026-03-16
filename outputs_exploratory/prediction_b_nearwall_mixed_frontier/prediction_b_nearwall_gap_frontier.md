# Prediction B Near-Wall Mixed Frontier

> Scope: exploratory mixed continuation beyond the confirmatory exact window.
> Families: `lorentzian_like_2d` vs `KR_like`.
> Action: `A2_full`.

No crossing is observed up to the scanned `gamma_max`, but the residual KR-vs-Lor2D gap stays finite and contracts substantially toward the right boundary.

| N | status | ΔA at γ_max | |ΔA| min | γ at min gap |
|---|---|---:|---:|---:|
| 52 | no_cross_kr_favored | -0.201 | 0.201 | 2.0 |
| 56 | no_cross_kr_favored | -0.839 | 0.839 | 2.0 |
| 60 | no_cross_kr_favored | -6.442 | 6.442 | 2.0 |
| 64 | no_cross_kr_favored | -6.064 | 6.064 | 2.0 |
| 68 | no_cross_kr_favored | -6.892 | 6.892 | 2.0 |
| 72 | no_cross_kr_favored | -12.772 | 12.772 | 2.0 |

Interpretation:
- Negative `ΔA_KR-L2D` means KR still wins; values near zero indicate near-degeneracy.
- These mixed results extend the competition frontier, but they do not by themselves establish persistence of a bounded `gamma_c` beyond the exact `N<=44` confirmatory window.
- The practical value of this report is to separate a stronger claim from a weaker one: the exact bounded-window result remains confirmatory, while the large-`N` mixed frontier shows that the competition remains active rather than collapsing.

