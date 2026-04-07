"""Scan an ADAS3D export root and collect per-patient file paths.

Expected layout on the external HD:
    <root>/<patient_id>/   (one folder per patient)
        *.vtk              (cardiac meshes)
        *.csv              (scar statistics)
        ...                (other ADAS3D outputs)

Usage:
    import data_loading
    cases = data_loading.scan(r"F:/RM_TEKNON_DEVELOP")
    for c in cases:
        print(c.patient_id, c.vtk_files, c.csv_files)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List


@dataclass
class PatientCase:
    patient_id: str
    root: Path
    vtk_files: List[Path] = field(default_factory=list)
    csv_files: List[Path] = field(default_factory=list)


def scan(dataset_root: str | Path) -> List[PatientCase]:
    """Return one PatientCase per patient folder found under dataset_root."""
    root = Path(dataset_root)
    if not root.exists():
        raise FileNotFoundError(f"Dataset root not found: {root}")

    cases = []
    for patient_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        cases.append(
            PatientCase(
                patient_id=patient_dir.name,
                root=patient_dir,
                vtk_files=sorted(patient_dir.rglob("*.vtk")),
                csv_files=sorted(patient_dir.rglob("*.csv")),
            )
        )
    return cases
