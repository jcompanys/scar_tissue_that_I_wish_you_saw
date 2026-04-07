# LV Scar Segmentation

Research code for LV scar characterisation from ADAS3D exports.

Data lives on an external hard drive and is **never** committed to this repo.

## Setup

```bash
pip install -r requirements.txt
```

## Usage

```python
import data_loading

cases = data_loading.scan(r"F:/RM_TEKNON_DEVELOP")
for c in cases:
    print(c.patient_id, len(c.vtk_files), "VTK |", len(c.csv_files), "CSV")
```

## Structure

```
lv-scar-segmentation/
├── data_loading.py   # scan ADAS3D export root → list of PatientCase
├── notebooks/        # exploration notebooks (outputs gitignored)
├── results/          # outputs: manifests, figures (gitignored)
└── requirements.txt
```
