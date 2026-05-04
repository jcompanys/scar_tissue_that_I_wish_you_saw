---
name: LV Scar Segmentation Project
description: Current context for the TFG project: repo organization, live project structure, analysis notebooks, generated outputs, geometry utilities, and working conventions
type: project
---

## Project goal
Characterize LV scar geometry from ADAS3D cardiac MRI exports and correlate scar shape and distribution with patient demographics and clinical outcomes.

## Repository layout
- Repo root: `C:\Users\joan\Desktop\FEINA\UPF\TFG\repo\scar-char-v0.1`
- Main Git branch currently in use: `claude/charming-buck`
- Main project code lives in: `lv-scar-segmentation/`
- Project memory lives in: `memory/`
- Auxiliary repo folders also present at root: `.claude/`, `backup/`, `slides/`, and `tools/`

## Branch and worktree model
- `scar-char-v0.1` is the canonical Git repository folder for normal work.
- Commits are made from the repo root `scar-char-v0.1`, even when edited files are inside `lv-scar-segmentation/`.
- The branch is attached to the whole repo/worktree, not to the `lv-scar-segmentation/` subfolder.
- `.claude/worktrees/` may contain extra assistant-created worktrees for other sessions or branches, but they are secondary.
- `memory/` is the shared project memory folder that should be updated when the project organization or pipeline changes.

## Live project structure
Current `lv-scar-segmentation/` layout:

```text
lv-scar-segmentation/
|-- .gitignore
|-- brief_report.tex
|-- requirements.txt
|-- __pycache__/
|-- notebooks/
|   |-- .gitkeep
|   |-- 01_ISOMAPS.ipynb
|   |-- 01_scar_exploration.ipynb
|   |-- 02_clinical_eda.ipynb
|   |-- 03_xyz_shape_analysis.ipynb
|   |-- 04_bspline_polar_test.ipynb
|   |-- 05_inspect_patient_50122488.ipynb
|   |-- 06_polar_diagnostics.ipynb
|   `-- 07_area_surface_polar_maps.ipynb
|-- src/
|   |-- __init__.py
|   |-- clinical_data.py
|   |-- cone_bspline_simple.py
|   |-- data_loading.py
|   |-- lv_geometry.py
|   `-- mesh_utils.py
`-- results/
    |-- split.json
    |-- 01_scar_exploration/
    |-- 02_clinical_eda/
    |-- 03_xyz_analysis/
    |-- 04_bspline_polar_test/
    |-- 06_polar_diagnostics/
    `-- 07_area_surface_polar_maps/
```

## File and folder roles
- `src/data_loading.py`: dataset scanner and loader utilities for ADAS3D-derived LV anatomy, optional RV anatomy, scar surfaces, layer files, and stats CSVs.
- `src/clinical_data.py`: clinical registry ingestion, cleaning, column normalization, and grouped EDA helpers.
- `src/mesh_utils.py`: shared mesh loaders and connected-fragment filtering utilities. It includes `best_shell`, `best_rv_mesh`, `get_cz_mesh`, `connected_fragments`, and `filter_fragments`.
- `src/lv_geometry.py`: shared LV geometry source of truth. It includes long-axis estimation, cross-section basis construction, PCA-based septal reference, RV-mesh-based septal/RV reference, anatomical transforms, and polar-coordinate projection.
- `src/cone_bspline_simple.py`: B-spline/cone projection logic; now imports shared mesh and LV geometry helpers instead of redefining them.
- `brief_report.tex`: LaTeX report summarizing the scar exploration and clinical EDA outputs.
- `notebooks/`: main exploratory and method-development workflows.
- `results/`: generated figures, caches, and split metadata exported by the notebooks.
- `__pycache__/`: generated Python bytecode cache; not project source.

## Dataset
- External HD: `F:/RM_TEKNON_DEVELOP`
- Expected LV layout: `<root>/<year>/<patient_id>/Basal/Data/<DE-MRI variant>/LV/`
- Optional RV layout: `<root>/<year>/<patient_id>/Basal/Data/<DE-MRI variant>/RV/` or `Right Ventricle/`
- 168 total patients scanned in the current working dataset
- 105 have valid Core Zone scar data
- 80/20 split frozen in `results/split.json`: 84 known and 21 held-out
- `results/01_scar_exploration/split.json` also exists as a notebook-level copy of the split metadata

## data_loading.py
Scans the HD dataset and returns a `PatientCase` dataclass per patient with paths to:
- LV VTK anatomy surfaces: Endo, Epi, LV, Myocardium
- Optional RV VTK anatomy surfaces: RV, RV Endo, RV Epi
- Tissue surfaces: Core, BZ, Healthy, Scar, Power Paths
- Layer files: `Layer_10` through `Layer_90`
- Stats CSV files

Current `PatientCase` fields include `rv_dir` and `rv_anatomy` so `best_rv_mesh(case)` can work with newly scanned cases. `mesh_utils.best_rv_mesh` also remains tolerant of older `PatientCase` objects that do not yet have `rv_anatomy`, avoiding the previous `AttributeError`.

## clinical_data.py
Reads the study clinical registry CSV with Spanish headers, comma decimals, and `DD/MM/YYYY` dates.
- `load(csv_path)` returns a cleaned DataFrame and auto-detects encoding and separator
- Around 90 columns are mapped from Spanish names to `snake_case`
- Handles dates, nullable integers, decimal conversion, channels yes/no, and NYHA Roman numerals
- Splits duplicate `FECHA IAM` information into `date_mi` and `date_mi_history`
- Exposes `COLUMN_GROUPS` for grouped EDA summaries

## Geometry and orientation conventions
- Shared geometry code should live in `src/lv_geometry.py`; notebooks should import it rather than copy helper definitions.
- Shared mesh-loading and fragment utilities should live in `src/mesh_utils.py`.
- For 3-D diagnostic grids, the preferred convention is `+Z = apex -> base`, with the camera looking from apex toward base and the septal/RV direction on screen-right.
- Putting every apex at the same origin is useful for coordinate computation but is not required for thumbnail visualization; view centering can be based on mesh bounds or centroid.
- For polar maps, the critical invariant is a consistent `theta=0` anatomical reference. If `theta=0` changes across patients, population density maps become unreliable.
- ADAS/scanner coordinates may look visually organized because the exported coordinate frame is already consistent. Mesh-derived anatomical alignment can look random if the septal/RV reference estimate is unstable.
- The project now supports both PCA basal-ring reference estimation and RV-mesh-based reference estimation when RV anatomy is available.

## Notebook inventory
- `01_ISOMAPS.ipynb`: split-out ISOMAP/fragmentation notebook for intrinsic scar-fragment shape analysis.
- `01_scar_exploration.ipynb`: dataset scan, valid scar-case selection, split creation, LV/scar visualization, anatomically aligned side/top views, pre-flattening polar-frame diagnostic, polar mapping, superposition, sex-stratified maps, ring density, and AHA summaries. Recent updates import shared helpers from `src.mesh_utils` and `src.lv_geometry`.
- `02_clinical_eda.ipynb`: clinical registry cleaning checks plus demographics, scar burden summaries, risk factors, echo variables, subgroup comparisons, correlations, arrhythmias, medications, and outcomes plots.
- `03_xyz_shape_analysis.ipynb`: XYZ-coordinate shape analysis and radial-envelope comparisons, with cached intermediate analyses.
- `04_bspline_polar_test.ipynb`: B-spline / polar-method experimentation notebook. Recent updates import shared mesh and geometry helpers while preserving notebook-facing wrapper outputs.
- `05_inspect_patient_50122488.ipynb`: patient-specific inspection notebook for LV/scar orientation and PCA long-axis debugging.
- `06_polar_diagnostics.ipynb`: polar-map reference diagnostics. It includes D1 basal-ring arrows, D2 angular disagreement across reference methods, D3 ADAS direction clustering, D4 three-panel 3-D grids, and D5 split-half polar density checks.
- `07_area_surface_polar_maps.ipynb`: area-based and surface-normalized polar maps. It accumulates CZ and LV shell cell areas into polar bins, computes CZ/LV coverage maps, and exports AHA/ring summaries.

## Results inventory
- `results/01_scar_exploration/`: scar exploration figures including raw and anatomically aligned 3-D grids, polar grids, population superposition, sex-stratified polar maps, ring-density maps, AHA summaries, fragmentation figures, and method diagnostics.
- `results/02_clinical_eda/`: clinical EDA outputs including missingness, demographics, scar summaries, risk factors, echo plots, subgroup analyses, device/arrhythmia/medication/outcome figures, and correlation CSV/PNG exports.
- `results/03_xyz_analysis/`: `xyz_cases_overview.png`, `xyz_radial_envelopes.png`, and `analyses_cache.pkl`.
- `results/04_bspline_polar_test/`: B-spline polar test image outputs (`embedded_output_*.png`).
- `results/06_polar_diagnostics/`: orientation/QC figures including `d1_basal_ring_arrows.png`, `d2_angular_disagreement.png`, `d3_adas_direction_clustering.png`, `d4_three_panel_3d_grid.png`, `d5_split_half_density.png`, and polar-map QC images.
- `results/07_area_surface_polar_maps/`: generated by `07_area_surface_polar_maps.ipynb`; stores area-based population bullseyes, patient coverage grids, AHA/ring area coverage CSVs, and summary plots.
- Root-level `results/split.json`: canonical frozen known/held-out partition used across analyses.

## Working tree status noted during memory refresh
- `README.md` was refreshed to mirror the current source/notebook/results layout and to document the RV-aware polar-map pipeline.
- `notebooks/01_ISOMAPS.ipynb` and `notebooks/06_polar_diagnostics.ipynb` existed as untracked notebooks.
- `src/lv_geometry.py` and `src/mesh_utils.py` existed as untracked shared modules.
- `notebooks/01_scar_exploration.ipynb`, `notebooks/04_bspline_polar_test.ipynb`, `src/cone_bspline_simple.py`, and `src/data_loading.py` had local edits.
- `src/data_loading.py` had been updated with RV-aware scanning, and `src/mesh_utils.best_rv_mesh` had been made tolerant of older case objects without `rv_anatomy`.
- `notebooks/_patch_qc2.py` was a temporary one-off notebook patch helper and should not be committed as project source; `lv-scar-segmentation/.gitignore` now ignores `notebooks/_patch_*.py`.
- Do not revert these local edits unless explicitly requested.

## Working rule
Always treat `C:\Users\joan\Desktop\FEINA\UPF\TFG\repo\scar-char-v0.1` as the source of truth:
- start assistant sessions there
- use its `memory/` folder as the project memory
- treat `.claude/worktrees/` as auxiliary branch or session folders, not as the main home of the project

## Next steps
- Reload notebook modules after code edits and rescan cases so old in-memory `PatientCase` objects gain `rv_anatomy`.
- Validate whether the final polar `theta=0` reference should use ADAS/scanner coordinates, PCA basal-ring geometry, or RV-mesh direction.
- Fix or explicitly document the bull's-eye display direction if the plotted AHA orientation still uses `set_theta_direction(1)` while expecting anterior at 12 o'clock.
- Consolidate any remaining duplicated geometry, mesh-loading, or plotting helpers from notebooks into `src/`.
- Continue correlating imaging-derived scar characteristics with demographic and clinical outcomes.
