# LV Scar Characterization

Research notebooks and helper code for left-ventricular scar characterization
from ADAS3D cardiac MRI exports.

The repository stores code only. Patient data, generated figures, CSV outputs,
and old exploratory notebooks are ignored by git.

## Project Layout

```text
scar-char-v0.1/
|-- lv-scar-segmentation/
|   |-- notebooks/
|   |   |-- M01_scar_exploration.ipynb
|   |   |-- M02_clinical_eda.ipynb
|   |   |-- M02_isomaps.ipynb
|   |   |-- M03_polar_diagnostics.ipynb
|   |   |-- M04_scar_analysis.ipynb
|   |   `-- M05_retrieval_based_generation.ipynb
|   |-- src/
|   |   |-- data_loading.py
|   |   |-- mesh_utils.py
|   |   |-- lv_geometry.py
|   |   |-- clinical_data.py
|   |   |-- scar_characterization.py
|   |   `-- scar_analysis.py
|   |-- results/
|   |   `-- .gitkeep
|   |-- old_results/        # ignored
|   |-- notebooks/old/      # ignored
|   `-- requirements.txt
|-- auxiliary/
|   |-- memory/
|   |-- slides/
|   |-- tools/
|   `-- writting/
|-- backup/
|   `-- lv-scar-segmentation/
`-- README.md
```

## Notebook Pipeline

Run notebooks from `lv-scar-segmentation/` so imports like `from src...` work.

1. `M01_scar_exploration.ipynb` scans the external ADAS3D dataset and builds the
   full scar cohort. No train/test split is created.
2. `M02_clinical_eda.ipynb` summarizes the clinical registry and marks patients
   with available MRI segmentation.
3. `M02_isomaps.ipynb` explores scar fragment shape with Isomap diagnostics.
4. `M03_polar_diagnostics.ipynb` validates LV axis, RV/septal reference, and
   polar-map orientation.
5. `M04_scar_analysis.ipynb` produces cohort-level scar descriptors, regional
   summaries, and statistical outputs.
6. `M05_retrieval_based_generation.ipynb` runs retrieval-based scar generation
   and leave-one-out validation on the full cohort.

All current analysis notebooks use the full available cohort. The old numbered
notebooks (`02`, `03`, `04`, `05`, `07`) live in `notebooks/old/` and are not
tracked.

## Data And Outputs

Default local paths used by the notebooks:

```text
MRI data:       F:/RM_TEKNON_DEVELOP
Clinical CSV:   C:/Users/joan/Desktop/FEINA/UPF/TFG/develop-vt.csv
Outputs:        lv-scar-segmentation/results/
```

Generated outputs are ignored. The previous output folder was renamed to
`old_results/`; new runs should write fresh files into `results/`.

## Auxiliary Files And Backup

Auxiliary project materials live under `auxiliary/`:

- `auxiliary/memory/` stores project memory and reference notes.
- `auxiliary/slides/` stores presentation materials.
- `auxiliary/tools/` stores helper tooling outside the main package.
- `auxiliary/writting/` stores writing materials.

The `backup/lv-scar-segmentation/` folder is a working copy snapshot of the main
`lv-scar-segmentation/` project folder.

## Setup

```bash
cd lv-scar-segmentation
python -m pip install -r requirements.txt
```

Then open JupyterLab and run the notebooks in order.
