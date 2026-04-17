---
name: LV Scar Segmentation Project
description: Main context for the TFG project: repo organization, branch/worktree setup, data pipeline, and next steps
type: project
---

## Project goal
Characterize LV scar geometry from ADAS3D cardiac MRI exports and correlate scar shape and distribution with patient demographics and clinical outcomes.

## Repository layout
- Repo root: `C:\Users\joan\Desktop\FEINA\UPF\TFG\repo\scar-char-v0.1`
- Main Git branch currently in use: `claude/charming-buck`
- Main project code lives in: `lv-scar-segmentation/`
- Project memory lives in: `memory/`
- Claude metadata and auxiliary worktrees live in: `.claude/`

## Branch and worktree model
- `scar-char-v0.1` is now the one real Git repository folder to use for normal work.
- Commits are made from the repo root `scar-char-v0.1`, even when the edited files are inside `lv-scar-segmentation/`.
- The branch is attached to the whole repo/worktree, not to the `lv-scar-segmentation/` subfolder.
- `.claude/worktrees/` may contain extra Claude-created worktrees for other sessions or branches, but they are secondary. The canonical working copy is the repo root above.
- `memory/` is the shared project memory folder that should be updated to reflect the current main organization.

## Main code layout
Current main layout in this branch is flat:

```text
lv-scar-segmentation/
|-- data_loading.py
|-- clinical_data.py
|-- requirements.txt
|-- notebooks/
|   |-- 01_scar_exploration.ipynb
|   `-- 02_clinical_eda.ipynb
`-- results/
```

## Dataset
- External HD: `F:/RM_TEKNON_DEVELOP`
- Expected layout: `<root>/<year>/<patient_id>/Basal/Data/DE-MRI/LV/`
- 168 total patients scanned
- 105 have valid Core Zone scar data
- 80/20 split frozen in `results/split.json`: 84 known and 21 held-out

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

## Notebook status
`01_scar_exploration.ipynb` includes:
1. Dataset scan and missing-case diagnostics
2. Frozen 80/20 split
3. 3D visualization of LV and scar
4. Polar map generation and overlays
5. ISOMAP-based scar shape analysis

`02_clinical_eda.ipynb` includes:
1. Registry completeness review
2. Demographics and risk factors
3. Scar geometry and echo variables
4. Correlations, subgroup comparisons, and outcomes preview

## Working rule
Always treat `C:\Users\joan\Desktop\FEINA\UPF\TFG\repo\scar-char-v0.1` as the source of truth:
- start Claude sessions there
- use its `memory/` folder as the project memory
- treat `.claude/worktrees/` as auxiliary branch or session folders, not as the main home of the project

## Next steps
- Link `patient_id` from the clinical CSV to `PatientCase.patient_id`
- Continue correlating imaging-derived scar characteristics with demographic and clinical outcomes
- Keep `claude/charming-buck` as the main reference branch unless a new task needs an isolated worktree
