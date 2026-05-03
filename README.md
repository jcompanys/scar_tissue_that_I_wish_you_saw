# LV Scar Segmentation

Research code for LV scar characterisation from ADAS3D cardiac MRI exports and
the local clinical registry.

Data stays local and is never committed to this repo.

## Organization

```text
scar-char-v0.1/
|-- lv-scar-segmentation/
|   |-- notebooks/
|   |   |-- 01_scar_exploration.ipynb
|   |   |-- 01_ISOMAPS.ipynb
|   |   |-- 02_clinical_eda.ipynb
|   |   |-- 03_xyz_shape_analysis.ipynb
|   |   |-- 04_bspline_polar_test.ipynb
|   |   |-- 05_inspect_patient_50122488.ipynb
|   |   `-- 06_polar_diagnostics.ipynb
|   |-- results/
|   |   |-- 01_scar_exploration/
|   |   |-- 02_clinical_eda/
|   |   |-- 03_xyz_analysis/
|   |   |-- 04_bspline_polar_test/
|   |   |-- 06_polar_diagnostics/
|   |   `-- split.json
|   |-- src/
|   |   |-- __init__.py
|   |   |-- clinical_data.py
|   |   |-- cone_bspline_simple.py
|   |   |-- data_loading.py
|   |   |-- lv_geometry.py
|   |   `-- mesh_utils.py
|   `-- requirements.txt
|-- memory/
`-- slides/
```

Notebook imports use the package under `lv-scar-segmentation/src`, for example:

```python
from src.data_loading import scan
from src.mesh_utils import best_shell, best_rv_mesh, get_cz_mesh
from src.lv_geometry import prepare_polar_geometry, project_points_to_polar
```

## Current Pipeline

1. `src/data_loading.py` scans the external ADAS3D dataset and builds one
   `PatientCase` per patient. It resolves LV anatomy, optional RV anatomy,
   tissue surfaces, layer files, and stats CSVs.
2. `src/mesh_utils.py` centralises mesh selection and scar-fragment helpers:
   best LV shell, optional RV mesh, Core Zone mesh, connected fragments, and
   fragment filtering.
3. `src/lv_geometry.py` is the shared geometry source of truth: LV long axis,
   cross-section basis, septal/RV reference estimation, anatomical transform,
   and polar-coordinate projection.
4. `src/cone_bspline_simple.py` contains the B-spline/cone polar workflow and
   imports the shared helpers instead of redefining them.
5. `notebooks/06_polar_diagnostics.ipynb` validates the angular reference used
   by polar maps before population-level density maps are trusted.

## Data Conventions

- Main external dataset: `F:/RM_TEKNON_DEVELOP`
- Expected LV layout: `<root>/<year>/<patient_id>/Basal/Data/<DE-MRI variant>/LV/`
- Optional RV layouts: `<root>/<year>/<patient_id>/Basal/Data/<DE-MRI variant>/RV/`
  or `<root>/<year>/<patient_id>/Basal/Data/<DE-MRI variant>/Right Ventricle/`
- Frozen known/held-out split: `lv-scar-segmentation/results/split.json`

For 3-D diagnostic views, the current preferred convention is:

- `+Z` means apex to base.
- The camera looks from apex toward base when checking polar orientation.
- Septal/RV direction should appear on screen-right for the diagnostic view.
- The extra AHA display roll is only for matching the bull's-eye layout.

## Setup

```bash
pip install -r lv-scar-segmentation/requirements.txt
```

Run notebooks from `lv-scar-segmentation/` or add that folder to `sys.path` so
`src.*` imports resolve correctly.
