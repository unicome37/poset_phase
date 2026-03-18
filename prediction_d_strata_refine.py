"""Refine Prediction D stratified block-permutation p-values on selected strata.

Goal: take v6 stratified table, select "candidate freeze windows" where multiple
I_cg variants are simultaneously positive & significant, then re-run the same
stratified block permutation with larger n_perm only for those strata.

This is intentionally a narrow, fast step:
- It reuses the already materialized zeta-scan ranking table (variants) from v7.
- It avoids recomputing the full pipeline (independent gains, etc.).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import prediction_d_dynamic_validation as pdv


def _parse_csv_path(s: str) -> Path:
    p = Path(s)
    if not p.exists():
        raise FileNotFoundError(str(p))
    return p


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Refine Prediction D stratified blockperm on selected strata.")
    p.add_argument(
        "--strata-csv",
        type=_parse_csv_path,
        default=Path("outputs_exploratory/prediction_d_dynamic_v6/cg_blockperm_delta_rank_stratified.csv"),
        help="v6 stratified p-table CSV path.",
    )
    p.add_argument(
        "--zeta-rank-csv",
        type=_parse_csv_path,
        default=Path("outputs_exploratory/prediction_d_dynamic_v7_blockperm/cg_zeta_scan_rankings_variants.csv"),
        help="zeta ranking (variants) CSV path to re-test.",
    )
    p.add_argument(
        "--out",
        default="outputs_exploratory/prediction_d_dynamic_v8_refine",
        help="Output directory.",
    )
    p.add_argument("--n-perm", type=int, default=20000, help="Permutation count for refined stratified p-values.")
    p.add_argument("--seed", type=int, default=11, help="RNG seed for permutations.")
    p.add_argument(
        "--require-variants",
        default="full,switch,no_switch",
        help="Comma-separated icg_variant set required for a stratum to be selected.",
    )
    p.add_argument("--alpha", type=float, default=0.05, help="Selection threshold on v6 p-values.")
    p.add_argument(
        "--require-positive",
        action="store_true",
        help="If set, require obs_mean_spearman > 0 in v6 selection (recommended).",
    )
    p.add_argument(
        "--keys",
        default="",
        help=(
            "Optional explicit strata keys to refine, as comma-separated triples 'n:keep:gamma'. "
            "If provided, v6 selection is skipped and only these (n,keep_ratio,gamma) are refined."
        ),
    )
    return p


def main() -> None:
    args = build_arg_parser().parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    require_variants = [x.strip() for x in str(args.require_variants).split(",") if x.strip()]
    alpha = float(args.alpha)

    key_cols = ["n", "keep_ratio", "gamma"]

    explicit_keys = [x.strip() for x in str(args.keys).split(",") if x.strip()]
    if explicit_keys:
        triples = []
        for item in explicit_keys:
            parts = item.split(":")
            if len(parts) != 3:
                raise ValueError(f"Invalid --keys item (expected n:keep:gamma): {item!r}")
            n, keep, gamma = int(parts[0]), float(parts[1]), float(parts[2])
            triples.append({"n": n, "keep_ratio": keep, "gamma": gamma})
        chosen_keys = pd.DataFrame(triples, columns=key_cols).drop_duplicates().reset_index(drop=True)
        chosen = pd.DataFrame()
    else:
        v6 = pd.read_csv(args.strata_csv)
        # Normalize types
        for c in ["n_perm", "n_blocks", "obs_mean_spearman", "p_abs_mean_spearman", "n", "keep_ratio", "gamma"]:
            if c in v6.columns:
                v6[c] = pd.to_numeric(v6[c], errors="coerce")
        v6["icg_variant"] = v6["icg_variant"].astype(str)

        sel = v6[v6["icg_variant"].isin(require_variants)].copy()
        sel = sel[sel["p_abs_mean_spearman"] < alpha].copy()
        if bool(args.require_positive):
            sel = sel[sel["obs_mean_spearman"] > 0.0].copy()

        # Find (n, keep_ratio, gamma) strata where all required variants pass.
        have = (
            sel.groupby(key_cols)["icg_variant"]
            .apply(lambda s: set(s.tolist()))
            .reset_index()
            .rename(columns={"icg_variant": "variants_present"})
        )
        required = set(require_variants)
        have["is_selected"] = have["variants_present"].apply(lambda s: required.issubset(s))
        chosen_keys = have[have["is_selected"]][key_cols].copy()

    if chosen_keys.empty:
        (out_dir / "report.md").write_text(
            "\n".join(
                [
                    "# Prediction D v8 Refine",
                    "",
                    "No strata selected from v6 under current criteria.",
                    f"- require_variants: `{require_variants}`",
                    f"- alpha: `{alpha}`",
                    f"- require_positive: `{bool(args.require_positive)}`",
                    f"- keys: `{str(args.keys)}`",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        return

    if not explicit_keys:
        chosen = sel.merge(chosen_keys, on=key_cols, how="inner").copy()
        chosen = chosen.sort_values(["n", "keep_ratio", "gamma", "icg_variant"]).reset_index(drop=True)
        chosen.to_csv(out_dir / "selected_strata_v6.csv", index=False, encoding="utf-8-sig")

    zeta_rank = pd.read_csv(args.zeta_rank_csv)
    # Filter to the chosen strata + required variants only, to keep the permutation small/targeted.
    zeta_rank["icg_variant"] = zeta_rank["icg_variant"].astype(str)
    zeta_rank = zeta_rank[zeta_rank["icg_variant"].isin(require_variants)].copy()
    zeta_rank = zeta_rank.merge(chosen_keys, on=key_cols, how="inner")

    refined = pdv.stratified_blockperm_delta_rank_test(
        zeta_rank,
        strata_cols=("icg_variant", "n", "keep_ratio", "gamma"),
        n_perm=int(args.n_perm),
        seed=int(args.seed),
    )
    refined.to_csv(out_dir / "cg_blockperm_delta_rank_stratified_refined.csv", index=False, encoding="utf-8-sig")

    # Summaries for the report.
    n_strata = int(chosen_keys.shape[0])
    worst = refined.sort_values("p_abs_mean_spearman", ascending=False).head(10)
    best = refined.sort_values("p_abs_mean_spearman", ascending=True).head(10)

    lines = ["# Prediction D v8 Refine", ""]
    lines.append(f"- v6 strata CSV: `{Path(args.strata_csv).as_posix()}`")
    lines.append(f"- zeta-rank CSV: `{Path(args.zeta_rank_csv).as_posix()}`")
    lines.append(f"- require_variants: `{require_variants}`")
    lines.append(f"- alpha (v6 selection): `{alpha}`")
    lines.append(f"- require_positive: `{bool(args.require_positive)}`")
    if explicit_keys:
        lines.append(f"- explicit keys: `{explicit_keys}`")
    lines.append(f"- selected (n,keep_ratio,gamma) strata: `{n_strata}`")
    lines.append(f"- refined n_perm: `{int(args.n_perm)}`")
    lines.append("")
    lines.append("Outputs:")
    if not explicit_keys:
        lines.append("- `selected_strata_v6.csv`")
    lines.append("- `cg_blockperm_delta_rank_stratified_refined.csv`")
    lines.append("")
    lines.append("Best (lowest p) strata (refined):")
    for _, r in best.iterrows():
        lines.append(
            f"- `{r['icg_variant']}` n={int(r['n'])} keep={float(r['keep_ratio']):.2f} gamma={float(r['gamma']):.1f}: "
            f"rho={float(r['obs_mean_spearman']):.3f}, p={float(r['p_abs_mean_spearman']):.4g}"
        )
    lines.append("")
    lines.append("Worst (highest p) strata (refined):")
    for _, r in worst.iterrows():
        lines.append(
            f"- `{r['icg_variant']}` n={int(r['n'])} keep={float(r['keep_ratio']):.2f} gamma={float(r['gamma']):.1f}: "
            f"rho={float(r['obs_mean_spearman']):.3f}, p={float(r['p_abs_mean_spearman']):.4g}"
        )
    (out_dir / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
