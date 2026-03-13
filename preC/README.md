# Prediction C — Paper Skeleton

## Working Title

**Hierarchy Integration Lowers Combinatorial Entropy: A Three-Tier Correlational Study in Finite Causal Posets**

## Status

Skeleton stage. Three-tier validation completed; Simpson's Paradox analyzed; GPT review incorporated.

## File Structure

| File | Section | Status |
|------|---------|--------|
| `PredictionC_Section1_Abstract_Introduction.md` | Abstract + §1 Introduction | skeleton |
| `PredictionC_Section2_Methods.md` | §2 Methods | skeleton |
| `PredictionC_Section3_Results.md` | §3 Three-Tier Results | skeleton |
| `PredictionC_Section4_Simpson.md` | §4 Simpson's Paradox | skeleton |
| `PredictionC_Section5_Components.md` | §5 Component Decomposition & HII Narrowing | skeleton |
| `PredictionC_Section6_Discussion.md` | §6 Discussion, Limitations, Outlook | skeleton |

## Relation to Other Papers

| Paper | Core claim | Relation to C |
|-------|-----------|---------------|
| **Prediction A** | gcon selects higher-d Lorentzian; 4D wins under consistency | C provides the mechanism: deeper hierarchy → lower entropy |
| **Prediction B** | Bounded γ_c for Lor2D vs KR; non-circular dim_consistency replacement | C explains *why* the structural advantage exists |

### B→C Logical Dependency

- **Logically independent**: C's HII–log_H correlation can be verified without B being true.
- **Semantically layered**: If B+C both hold, C elevates B from "numerical coincidence" to "structurally supported regularity."
- **Three-prediction tower**: A (dimension selection) → B (self-consistency) → C (structural mechanism). Each layer adds explanatory depth.

## Key Quantitative Results

- Tier 1 (fixed N): partial_r(HII, log_H | N) = −0.578
- Tier 2 (46 pairs, N=30–56): HII_delta vs log_H_delta r = −0.834, p < 0.001
- Tier 3 (92 samples, N=30–56): layer_count → cg_switch_rate r = −0.874, p < 0.001
- Simpson's Paradox: N is the sole confound; original controls yield +0.336 (spurious)
- Near-wall: P5~P95 moderate filter rescued 6 MLR survivors each at N=52 and N=56; correlation unchanged (−0.836 → −0.834)

## Data Paths

- Tier 1: `outputs_exploratory/prediction_c_comprehensive/tier1_*.csv`
- Tier 2: `outputs_exploratory/prediction_c_comprehensive/tier2_*.csv`
- Tier 3: `outputs_exploratory/prediction_c_comprehensive/tier3_*.csv`
- Full validation report: `PredictionC_精确表述与验证报告.md`
