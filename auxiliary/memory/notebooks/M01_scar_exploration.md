# M01 — LV Scar Exploration

**Focus:** Core Zone (CZ) meshes only.  
**Output directory:** `results/01_scar_exploration/`

---

## Overview

This notebook is the entry point of the LV scar analysis pipeline. It scans the raw dataset stored on an external hard drive, performs a basic Exploratory Data Analysis (EDA), establishes a reproducible train/held-out split, and produces visual quality-control outputs for the known cohort.

The declared pipeline is:

1. Scan + EDA  
2. 80/20 train/held-out split  
3. 3-D visualisation (aligned overlay + per-patient grid)  

---

## Section 0 — Imports & Configuration

All third-party libraries and project-local modules are imported here. The notebook locates the project root dynamically by walking up the directory tree until it finds `src/data_loading.py`, then appends that root to `sys.path`.

Key configuration constants:

| Constant | Value | Purpose |
|---|---|---|
| `HD_ROOT` | `F:/RM_TEKNON_DEVELOP` | Root path of the raw dataset on the external drive |
| `SEED` | `42` | Random seed for reproducible splits |
| `TRAIN_RATIO` | `0.80` | Fraction of valid cases assigned to the known set |
| `N_R` / `N_THETA` | `64` / `128` | Polar grid resolution (radial / angular) for future polar maps |
| `GRID3D_MAX` | `24` | Maximum cases shown in the 3-D grid |
| `ISOMAP_PTS` | `800` | Number of points sampled per mesh for ISOMAP |
| `ISOMAP_K` | `10` | Number of neighbours for the ISOMAP manifold |
| `MIN_VERTS_FRAGMENT` | `30` | Minimum vertex count to keep a disconnected mesh fragment |
| `FLIP_SEPTAL_DIRECTION` | `False` | Flips the septal rotation if the RV appears on the left side |
| `CLINICAL_CSV` | path to `develop-vt.csv` | External clinical registry used for sex stratification |

Output folders `results/` and `results/01_scar_exploration/` are created if they do not exist. PyVista is set to static Jupyter backend and matplotlib DPI is fixed at 120.

---

## Section 1 — Scan & EDA

### What it does

`scan(HD_ROOT)` (from `src.data_loading`) recursively discovers all patient cases on the external drive. Each `PatientCase` object exposes:

- `lv_dir` — path to the LV subfolder under `DE-MRI`
- `tissue_surfaces` — dictionary mapping tissue labels (`"core"`, `"border"`, `"scar"`) to `.vtk` file paths

Two filter predicates are defined:

- `has_lv(c)` — returns `True` if the case has a resolved LV directory
- `has_cz(c)` — returns `True` if a non-empty Core Zone `.vtk` file exists

`valid` is the list of cases that pass `has_cz`.

### Dataset overview table

The notebook prints and saves a breakdown of cases per year with three columns: total cases, cases with an LV directory, and cases with a Core Zone. The result is written to `eda_dataset_overview.csv`.

### Missing-case diagnostics

Cases that fail either filter are passed through `_diagnose_case(c)`, which walks the expected folder hierarchy step by step and returns the first failure reason found:

1. No `Basal/Data` path
2. No `DE-MRI` folder (tries three naming variants)
3. `DE-MRI` found but no LV subfolder
4. No `TISSUE` folder under the LV dir
5. `TISSUE` exists but contains no `.vtk` files
6. No `.vtk` file matching the Core Zone naming convention
7. CZ file exists but is empty (0 bytes)

Counts are grouped by failure reason and printed alongside the affected patient IDs. The full list is saved to `eda_missing_cases.csv`.

### Outputs

| File | Description |
|---|---|
| `eda_dataset_overview.csv` | Year-by-year case counts (total / has LV / has CZ) |
| `eda_missing_cases.csv` | Per-patient failure reasons for cases without LV or CZ |

---

## Section 2 — 80/20 Split

### What it does

The `valid` list is shuffled with a fixed seed (`SEED = 42`) and split into:

- **known** (~80 %) — used for all subsequent EDA, visualisation and analysis
- **held_out** (~20 %) — reserved, untouched until final evaluation

The split is immediately serialised to `results/split.json` as a JSON file containing the seed, the ratio, and the `year/patient_id` identifiers for each partition. This guarantees full reproducibility across runs and notebooks.

### Clinical sex loading

After the split, the notebook attempts to load patient sex from the external clinical CSV via `src.clinical_data.load()`. The resulting `sex_map` dictionary (`patient_id → "M"/"F"`) is used downstream for stratified plots. If the CSV path is `None` or the load fails, sex stratification is silently skipped.

### Outputs

| File | Description |
|---|---|
| `results/split.json` | Serialised train/held-out partition with seed and ratio |

---

## Section 3 — 3-D Visualisation

### Geometry helpers (Section 3 setup cell)

Before rendering, the notebook imports geometry utilities from `src.mesh_utils` and `src.lv_geometry` and wraps two of them in local functions that automatically inject the global `FLIP_SEPTAL_DIRECTION` flag:

- `_oriented_septal_reference(ref)` — computes a septal reference frame, optionally flipped
- `_compute_transform(lv_mesh, align_septum)` — computes the rigid alignment transform for a given LV mesh

### 3b — Per-patient 3-D grid

For every case in the `known` set, an off-screen PyVista plotter renders:

- The LV shell (endo/epi/myocardium) in light grey at 15 % opacity — obtained via `best_shell(case)` from `src.mesh_utils`
- The Core Zone mesh in crimson at 95 % opacity — obtained via `get_cz_mesh(case)`

Both meshes use smooth shading. Cases that raise any exception during loading or rendering are skipped with a warning. All thumbnails (300 × 260 px) are assembled into a matplotlib grid (6 columns, as many rows as needed) titled with year and patient ID, saved to `3d_grid_known.png`, and displayed inline.

### Outputs

| File | Description |
|---|---|
| `3d_grid_known.png` | Isometric 3-D thumbnail grid for all known-set cases |


---

## Dependencies

### External libraries

`numpy`, `matplotlib`, `pyvista`, `scipy`, `scikit-learn`, `tqdm`, `pandas`

### Project modules (`src/`)

| Module | Used for |
|---|---|
| `src.data_loading` | `scan()`, `PatientCase` — dataset discovery |
| `src.mesh_utils` | `best_shell()`, `get_cz_mesh()` — mesh loading helpers |
| `src.lv_geometry` | Rigid alignment transforms, septal reference, plane basis |
| `src.clinical_data` | Loading and normalising the clinical registry CSV |
