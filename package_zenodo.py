"""Package for Zenodo release.

Creates a zenodo_release/ directory with essential files:
  - Source code (.py)
  - Configuration files (.yaml)
  - Output CSVs (key results only)
  - Manuscript figures
  - Documentation (README, LICENSE, MANUSCRIPT_DRAFT, etc.)
"""
import shutil
from pathlib import Path

SRC = Path(".")
DST = Path("zenodo_release/poset_phase")

# Clean
if DST.exists():
    shutil.rmtree(DST)
DST.mkdir(parents=True)

# ── Core source code ──
for f in SRC.glob("*.py"):
    shutil.copy2(f, DST / f.name)

# ── Config files ──
for f in SRC.glob("*.yaml"):
    shutil.copy2(f, DST / f.name)

# ── Documents ──
for name in [
    "README.md",
    "README_old.md",
    "LICENSE",
    "requirements.txt",
    ".zenodo.json",
    ".gitignore",
    "MANUSCRIPT_DRAFT.md",
    "MANUSCRIPT_OUTLINE.md",
    "AUDIT_REPORT.md",
]:
    src = SRC / name
    if src.exists():
        shutil.copy2(src, DST / name)

# ── Manuscript figures ──
fig_src = SRC / "manuscript_figures"
if fig_src.exists():
    shutil.copytree(fig_src, DST / "manuscript_figures")

# ── Key output CSVs (confirmatory) ──
csv_dirs = [
    "outputs_confirmatory/frozen_exact",
    "outputs_confirmatory/medium_exact_scan",
]
for d in csv_dirs:
    src_dir = SRC / d
    if src_dir.exists():
        dst_dir = DST / d
        dst_dir.mkdir(parents=True, exist_ok=True)
        for csv in src_dir.glob("*.csv"):
            shutil.copy2(csv, dst_dir / csv.name)

# ── Ablation + noncyclic CSVs (exploratory — at repo root level) ──
ext_dirs = [
    (Path("d:/Kiro/outputs_exploratory/geometric_ablation_gamma_c"),
     DST / "outputs_exploratory/geometric_ablation_gamma_c"),
    (Path("d:/Kiro/outputs_exploratory/noncyclic_dim_replacement_gamma_c"),
     DST / "outputs_exploratory/noncyclic_dim_replacement_gamma_c"),
]
for src_dir, dst_dir in ext_dirs:
    if src_dir.exists():
        dst_dir.mkdir(parents=True, exist_ok=True)
        for csv in src_dir.glob("*.csv"):
            shutil.copy2(csv, dst_dir / csv.name)

# ── Create zip archive ──
zip_path = Path("zenodo_release/poset_phase")
shutil.make_archive(str(Path("zenodo_release/poset_phase_v1.0")), "zip", "zenodo_release", "poset_phase")

# ── Summary ──
total = sum(1 for _ in DST.rglob("*") if _.is_file())
print(f"Packaged {total} files into {DST}")
print(f"Archive: zenodo_release/poset_phase_v1.0.zip")
