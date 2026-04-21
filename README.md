# LV Scar Segmentation

Research code for LV scar characterisation from ADAS3D exports and the
local clinical registry.

Data stays local and is never committed to this repo.

## Organization

```text
scar-char-v0.1/
|-- lv-scar-segmentation/
|   |-- notebooks/
|   |-- results/
|   |-- src/
|   |   |-- __init__.py
|   |   |-- clinical_data.py
|   |   |-- cone_bspline_simple.py
|   |   `-- data_loading.py
|   `-- requirements.txt
|-- memory/          (gitignored)
`-- slides/          (gitignored)
```

Notebook imports now use the package under `lv-scar-segmentation/src`, for
example `from src.data_loading import scan`.

## Setup

```bash
pip install -r lv-scar-segmentation/requirements.txt
```
