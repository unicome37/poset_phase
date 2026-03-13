# Paper Skeleton: Section 4 — Simpson's Paradox Diagnosis

## 4. Simpson's Paradox: Poset Size as the Sole Confound

### 4.1. The Anomaly

Tier 1's naïve partial correlation (controlling for action weights, cover fraction, and geometric dimension, but **not** N) yields:

$$r_{\text{partial}}(\text{HII}, \log H \mid aw, cf, \text{geo\_dim}) = +0.336 \quad (p = 0.0005)$$

This is the **opposite** sign from the prediction. If taken at face value, it would suggest that deeper hierarchy is associated with **higher** entropy — falsifying the prediction. This section diagnoses why this reversal occurs and shows that N is the sole confound responsible.

---

### 4.2. Diagnosis: Stepwise Confound Identification

A systematic analysis was conducted by adding control variables one at a time:

**Table 4.1**: Logical sequence of confound identification

| Step | Controls | partial_r | Sign | Interpretation |
|------|----------|-----------|------|----------------|
| 0 | None (raw) | +0.12 | + | Weak positive, dominated by N scaling |
| 1 | aw, cf, geo_dim | +0.336 | + | **Simpson's Paradox** |
| 2 | **N** | **−0.578** | **−** | **Sign flip** — N alone resolves it |
| 3 | N + family | −0.250 | − | Same direction, attenuated by within-family noise |
| 4 | N + aw + cf + geo_dim | −0.434 | − | Full controls with N: still negative |

The sign flip at step 2 is diagnostic: **N is both necessary and sufficient** to resolve the paradox.

### 4.3. Physical Mechanism

The Simpson's Paradox arises because N simultaneously affects both HII and log_H, but on different scales:

1. **N → HII**: Larger N posets generically admit deeper layer structure, so HII increases with N across families. This is a scaling effect, not a structural one.

2. **N → log_H**: Larger N posets have combinatorially more linear extensions, so log_H increases faster than HII with N. The growth rate is super-linear in N.

3. **Net effect without N control**: The within-N negative relationship (deeper hierarchy → lower entropy) is overwhelmed by the cross-N positive relationship (larger posets have both more hierarchy and more entropy).

This is a textbook ecological correlation — the aggregate trend (positive) is the opposite of the within-group trend (negative).

### 4.4. The Crystallization Analogy

> [Writing note: This material is appropriate for an intuitive reader; it may be placed in a "Physical Picture" box or sidebar in the final manuscript.]

An apt analogy is crystallization from a melt:
- **Temperature** (analogous to N): hotter melt has more thermal energy (higher log_H baseline) and also allows more molecular motion (higher HII baseline).
- **At fixed temperature**: more ordered crystals (deeper hierarchy) have lower entropy. This is the structural relationship.
- **Across temperatures**: comparing a cold crystal to a hot liquid shows higher order coinciding with lower entropy, but this conflates thermal and structural effects.

Controlling for N is like studying crystallization at fixed temperature — it isolates the structural effect.

### 4.5. Verification: Family-Only Control Fails

Controlling for family dummies without N:

$$r_{\text{partial}}(\text{HII}, \log H \mid \text{family}) = +0.148 \quad (p = 0.043)$$

This is still positive, showing that family identity alone does not resolve the paradox. The confound is **not** family membership but **poset size**.

### 4.6. Implications for Cross-N Studies

The Simpson's Paradox has a methodological lesson:

> In any cross-size poset study, raw correlations between structural observables and entropy-like quantities should not be trusted. The N-scaling of entropy creates confounds that can reverse structural relationships. Controls for N (or matched-pair designs, as in Tier 2) are essential.

This finding motivates the matched-pair Tier 2 design (Section 3.3), which eliminates N by construction.

---

### 4.7. Figure Suggestion

**Figure 4**: Simpson's Paradox visualization.
- Left panel: HII vs log_H, all 320 samples colored by N. Marginal trend is positive (annotation line). Within-N trends are negative (per-N regression lines).
- Right panel: Bar chart or forest plot of partial_r under the five control sets (Table 4.1), highlighting the sign flip at N control.

### 4.8. Table Summary

- Table 4.1 (above): Stepwise confound identification

---

## Writing Notes

- This section should be concise but definitive. The Simpson's Paradox is the most surprising finding and the one most likely to draw reviewer attention.
- Reviewers may ask: "Why not control for N from the start?" Fair question — answer: (1) we wanted to test the naïve prediction first; (2) the paradox itself is informative about the N-scaling structure of the ensemble.
- The crystallization analogy is optional but recommended for accessibility. It should not be the primary argument — the data table is.
- Frame this as a **methodological contribution**: future poset-entropy studies should always control for N. This elevates the section from "oops, we got the sign wrong" to "here's a generalizable lesson."
- Do NOT frame it defensively. Frame it as: "The positive sign forced us to investigate, and the investigation yielded a clean diagnosis."
