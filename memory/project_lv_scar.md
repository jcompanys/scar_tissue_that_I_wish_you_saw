---
name: LV Scar Segmentation Project
description: Current context for the TFG project: repo organization, live project structure, analysis notebooks, generated outputs, and working conventions
type: project
---

## Project goal
Characterize LV scar geometry from ADAS3D cardiac MRI exports and correlate scar shape and distribution with patient demographics and clinical outcomes.

## Repository layout
- Repo root: `C:\Users\joan\Desktop\FEINA\UPF\TFG\repo\scar-char-v0.1`
- Main Git branch currently in use: `claude/charming-buck`
- Main project code lives in: `lv-scar-segmentation/`
- Project memory lives in: `memory/`
- Auxiliary repo folders also present at root: `.claude/`, `backup/`, and `tools/`

## Branch and worktree model
- `scar-char-v0.1` is now the one real Git repository folder to use for normal work.
- Commits are made from the repo root `scar-char-v0.1`, even when the edited files are inside `lv-scar-segmentation/`.
- The branch is attached to the whole repo/worktree, not to the `lv-scar-segmentation/` subfolder.
- `.claude/worktrees/` may contain extra Claude-created worktrees for other sessions or branches, but they are secondary. The canonical working copy is the repo root above.
- `memory/` is the shared project memory folder that should be updated to reflect the current main organization.

## Live project structure
Current `lv-scar-segmentation/` layout:

```text
lv-scar-segmentation/
|-- .gitignore
|-- brief_report.tex
|-- data_loading.py
|-- clinical_data.py
|-- requirements.txt
|-- __pycache__/
|-- notebooks/
|   |-- .gitkeep
|   |-- 01_scar_exploration.ipynb
|   |-- 02_clinical_eda.ipynb
|   |-- 03_xyz_shape_analysis.ipynb
|   `-- 04_bspline_polar_test.ipynb
`-- results/
    |-- split.json
    |-- 01_scar_exploration/
    |-- 02_clinical_eda/
    |-- 03_xyz_analysis/
    `-- 04_bspline_polar_test/
```

## File and folder roles
- `data_loading.py`: dataset scanner and loader utilities for ADAS3D-derived LV anatomy and scar surfaces.
- `clinical_data.py`: clinical registry ingestion, cleaning, column normalization, and grouped EDA helpers.
- `brief_report.tex`: LaTeX report summarizing the scar exploration and clinical EDA outputs.
- `notebooks/`: main exploratory and method-development workflows.
- `results/`: generated figures, caches, and split metadata exported by the notebooks.
- `__pycache__/`: generated Python bytecode cache; not project source.

## Dataset
- External HD: `F:/RM_TEKNON_DEVELOP`
- Expected layout: `<root>/<year>/<patient_id>/Basal/Data/DE-MRI/LV/`
- 168 total patients scanned
- 105 have valid Core Zone scar data
- 80/20 split frozen in `results/split.json`: 84 known and 21 held-out
- `results/01_scar_exploration/split.json` also exists as a notebook-level copy of the split metadata

## data_loading.py
Scans the HD dataset and returns a `PatientCase` dataclass per patient with paths to:
- VTK anatomy surfaces: Endo, Epi, LV, Myocardium
- Tissue surfaces: Core, BZ, Healthy, Scar, Power Paths
- Layer files: `Layer_10` through `Layer_90`
- Stats CSV files

## clinical_data.py
Reads the study clinical registry CSV with Spanish headers, comma decimals, and `DD/MM/YYYY` dates.
- `load(csv_path)` returns a cleaned DataFrame and auto-detects encoding and separator
- Around 90 columns are mapped from Spanish names to `snake_case`
- Handles dates, nullable integers, decimal conversion, channels yes/no, and NYHA Roman numerals
- Splits duplicate `FECHA IAM` information into `date_mi` and `date_mi_history`
- Exposes `COLUMN_GROUPS` for grouped EDA summaries

## Notebook inventory
- `01_scar_exploration.ipynb`: dataset scan, valid scar-case selection, split creation, LV/scar visualization (raw + anatomically aligned side/top views), pre-flattening polar-frame diagnostic (blue dot = septal ref), polar mapping (individual + superposition + sex-stratified + ring density), AHA summaries, fragmentation plots, and ISOMAP diagnostics.
- `02_clinical_eda.ipynb`: clinical registry cleaning checks plus demographics, scar burden summaries, risk factors, echo variables, subgroup comparisons, correlations, arrhythmias, medications, and outcomes plots.
- `03_xyz_shape_analysis.ipynb`: XYZ-coordinate shape analysis and radial-envelope comparisons, with cached intermediate analyses.
- `04_bspline_polar_test.ipynb`: newer B-spline / polar-method experimentation notebook currently present in the working tree.

## Results inventory
- `results/01_scar_exploration/`: scar exploration figures including:
  - `3d_grid_known.png` — raw side-view grid (pre-alignment)
  - `3d_grid_known_anatomic.png` — anatomically aligned side view (long axis +Z, septal ref +X)
  - `3d_grid_known_anatomic_top.png` — apex-to-base top view (septal ref upper-left, AHA convention)
  - `polar_grid_known.png` — per-patient individual bull's-eye polar maps
  - `polar_superposition_known.png` — any-patient footprint + scar density (with ANT/LAT/INF/SEPT labels)
  - `polar_sex_stratified.png` — scar density split by sex (M / F / All, shared colorscale)
  - `aha_scar_density_rings.png` — vertex-normalised scar coverage fraction per LV ring zone
  - `aha_scar_density*.png`, `fragmentation_summary.png`, `top_fragmented_3d.png`, `xyz_candidates.png`, ISOMAP/method-comparison diagnostics
- `results/02_clinical_eda/`: large set of clinical EDA outputs including missingness, demographics, scar summaries, risk factors, echo plots, subgroup analyses, device/arrhythmia/medication/outcome figures, and correlation CSV/PNG exports.
- `results/03_xyz_analysis/`: `xyz_cases_overview.png`, `xyz_radial_envelopes.png`, and `analyses_cache.pkl`.
- `results/04_bspline_polar_test/`: folder exists but was empty at the time this memory was updated.
- Root-level `results/split.json`: canonical frozen known/held-out partition used across analyses.

## Working tree status noted during memory refresh
- `notebooks/04_bspline_polar_test.ipynb` exists locally and was untracked in Git when this memory was refreshed.

## Working rule
Always treat `C:\Users\joan\Desktop\FEINA\UPF\TFG\repo\scar-char-v0.1` as the source of truth:
- start Claude sessions there
- use its `memory/` folder as the project memory
- treat `.claude/worktrees/` as auxiliary branch or session folders, not as the main home of the project

## Next steps
- Link `patient_id` from the clinical CSV to `PatientCase.patient_id`
- Consolidate the newer geometry experiments (`03_xyz_shape_analysis.ipynb` and `04_bspline_polar_test.ipynb`) into the documented workflow as methods stabilize
- Continue correlating imaging-derived scar characteristics with demographic and clinical outcomes
- Keep `claude/charming-buck` as the current reference branch unless a new task needs an isolated worktree
