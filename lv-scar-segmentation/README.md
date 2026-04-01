# LV Scar Segmentation

Minimal, scalable repository for loading and processing LV scar study exports.

This first implementation step focuses on dataset discovery and structure validation
for a root layout like:

F:/RM_TEKNON_DEVELOP/<YEAR>/<PATIENT_ID>/Basal/Data/DE-MRI/LV/...

## Quick Start

1) Install package in editable mode

pip install -e .

2) Run dataset scan

lv-scar-scan --dataset-root F:/RM_TEKNON_DEVELOP --out results/dataset_manifest.json

3) Review report

results/dataset_manifest.json contains all discovered year/patient cases and
whether each case passes minimum structure checks.

## What The Loader Checks

For each patient case, the loader validates:

- Required directories under Basal/Data/DE-MRI/LV
- Required key VTK files
- Presence of TISSUE statistics CSV files
- TRANSMURALITY or TRANSMURABILITY directory naming mismatch

Primary module:

src/lv_scar_segmentation/data_loading.py

## Current Source Files

- src/lv_scar_segmentation/data_loading.py
- src/lv_scar_segmentation/cli.py
- configs/default.yaml
- tests/test_data_loading.py

## Data Protection

### ⚠️ Critical Rules

1. **NEVER** commit to `data/` folders
2. **NEVER** push sensitive formats (`.nii.gz`, `.dcm`, `.h5`, etc.)
3. **ALWAYS** verify with `.gitignore` before `git add`
4. **ALWAYS** use `.gitkeep` placeholder files

### Check Before Committing

```bash
# See what's staged
git diff --cached

# See all untracked files in data/
git clean -n data/

# Remove accidentally added data file
git rm --cached data/raw/patient_01.nii.gz
echo "data/raw/*.nii.gz" >> .gitignore
git commit "Remove sensitive data file"
```

## Next Data Loading Steps

Planned next increments after structure scan:

- Parse selected metadata from .lmi/.lses files
- Index TISSUE surfaces and layer files
- Add a normalized internal manifest schema
- Add optional copy/link step into local data/raw staging

## Sensitive Data Reminder

Keep all source patient files outside git tracking. Current .gitignore already
blocks medical and derived formats such as .vtk, .nii.gz, .dcm, .csv.
