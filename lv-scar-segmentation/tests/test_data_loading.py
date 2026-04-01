from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from lv_scar_segmentation.data_loading import scan_dataset


def test_scan_dataset_minimum_structure(tmp_path: Path) -> None:
    case_root = tmp_path / "2021" / "02026677"

    # Create expected folders
    (case_root / "Basal" / "Data" / "DE-MRI" / "LV" / "EAM").mkdir(parents=True)
    (case_root / "Basal" / "Data" / "DE-MRI" / "LV" / "THICKNESS").mkdir(parents=True)
    (case_root / "Basal" / "Data" / "DE-MRI" / "LV" / "TISSUE").mkdir(parents=True)
    (case_root / "Basal" / "Data" / "DE-MRI" / "LV" / "TISSUE_CE").mkdir(parents=True)
    (case_root / "Basal" / "Data" / "DE-MRI" / "LV" / "TRANSMURALITY").mkdir(parents=True)

    # Create key files
    for rel in (
        "Basal/Data/DE-MRI/CardiacModelLVRV.vtk",
        "Basal/Data/DE-MRI/LV/Left Ventricle.vtk",
        "Basal/Data/DE-MRI/LV/Myocardium.vtk",
        "Basal/Data/DE-MRI/LV/Endo Layer.vtk",
        "Basal/Data/DE-MRI/LV/Epi Layer.vtk",
    ):
        full = case_root / rel
        full.parent.mkdir(parents=True, exist_ok=True)
        full.write_text("", encoding="utf-8")

    # Tissue CSV example
    tissue_csv = case_root / "Basal" / "Data" / "DE-MRI" / "LV" / "TISSUE" / "Histogram statistics.csv"
    tissue_csv.write_text("metric,value\nexample,1\n", encoding="utf-8")

    report = scan_dataset(tmp_path)

    assert report.total_years == 1
    assert report.total_cases == 1
    assert report.valid_cases == 1
    assert report.invalid_cases == 0
    assert len(report.cases[0].tissue_csv_files) == 1
