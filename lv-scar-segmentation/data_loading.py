"""ADAS3D dataset loader for RM_TEKNON_DEVELOP.

Dataset layout on the external HD:
    <root>/<year>/<patient_id>/Basal/Data/DE-MRI/LV/
        Endo Layer.vtk
        Epi Layer.vtk
        Left Ventricle.vtk
        Myocardium.vtk
        TISSUE/
            Core Surface.vtk
            Border Zone Surface.vtk
            Healthy Surface.vtk
            Scar Surface.vtk
            Layer_10.vtk ... Layer_90.vtk
            *.csv  (quantitative statistics)
        EAM/
        THICKNESS/
        TISSUE_CE/
        TRANSMURALITY/   (or TRANSMURABILITY/ in some exports)

Usage:
    from data_loading import scan

    cases = scan(r"F:/RM_TEKNON_DEVELOP")
    for c in cases:
        print(c)
        print(c.stats_csv)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional


# ── Named surface files inside LV/ ────────────────────────────────────────────

_ANATOMY_FILES: Dict[str, str] = {
    "endo":        "Endo Layer.vtk",
    "epi":         "Epi Layer.vtk",
    "lv":          "Left Ventricle.vtk",
    "myocardium":  "Myocardium.vtk",
}

_TISSUE_SURFACES: Dict[str, str] = {
    "core":        "Core Surface.vtk",
    "border_zone": "Border Zone Surface.vtk",
    "healthy":     "Healthy Surface.vtk",
    "scar":        "Scar Surface.vtk",
}


# ── Data model ─────────────────────────────────────────────────────────────────

@dataclass
class PatientCase:
    year: str
    patient_id: str
    root: Path                               # <root>/<year>/<patient_id>/
    lv_dir: Optional[Path]                   # .../Basal/Data/DE-MRI/LV/

    # Anatomical surfaces (key → Path or None if missing)
    anatomy: Dict[str, Optional[Path]] = field(default_factory=dict)

    # Scar/tissue surfaces (key → Path or None if missing)
    tissue_surfaces: Dict[str, Optional[Path]] = field(default_factory=dict)

    # Transmural layers: Layer_10.vtk … Layer_90.vtk
    tissue_layers: List[Path] = field(default_factory=list)

    # Quantitative CSV statistics from TISSUE/
    stats_csv: List[Path] = field(default_factory=list)

    def __str__(self) -> str:
        n_anatomy = sum(1 for v in self.anatomy.values() if v is not None)
        n_tissue  = sum(1 for v in self.tissue_surfaces.values() if v is not None)
        return (
            f"[{self.year}] {self.patient_id} | "
            f"anatomy {n_anatomy}/{len(_ANATOMY_FILES)}  "
            f"tissue surfaces {n_tissue}/{len(_TISSUE_SURFACES)}  "
            f"layers {len(self.tissue_layers)}  "
            f"CSVs {len(self.stats_csv)}"
        )


# ── Public API ─────────────────────────────────────────────────────────────────

def scan(dataset_root: str | Path) -> List[PatientCase]:
    """Scan the full dataset and return one PatientCase per patient.

    Parameters
    ----------
    dataset_root:
        Root of the external HD, e.g. ``r"F:/RM_TEKNON_DEVELOP"``.

    Returns
    -------
    List[PatientCase]
        Sorted by year then patient ID.
    """
    root = Path(dataset_root)
    if not root.exists():
        raise FileNotFoundError(f"Dataset root not found: {root}")

    cases: List[PatientCase] = []
    for year_dir in _year_dirs(root):
        for patient_dir in sorted(p for p in year_dir.iterdir() if p.is_dir()):
            cases.append(_load_case(year_dir.name, patient_dir))

    return cases


# ── Internal helpers ───────────────────────────────────────────────────────────

def _year_dirs(root: Path) -> List[Path]:
    return sorted(
        p for p in root.iterdir()
        if p.is_dir() and p.name.isdigit() and len(p.name) == 4
    )


def _load_case(year: str, patient_dir: Path) -> PatientCase:
    lv_dir = _find_lv_dir(patient_dir)

    anatomy: Dict[str, Optional[Path]] = {}
    tissue_surfaces: Dict[str, Optional[Path]] = {}
    tissue_layers: List[Path] = []
    stats_csv: List[Path] = []

    if lv_dir is not None:
        for key, filename in _ANATOMY_FILES.items():
            p = lv_dir / filename
            anatomy[key] = p if p.exists() else None

        tissue_dir = lv_dir / "TISSUE"
        if tissue_dir.is_dir():
            for key, filename in _TISSUE_SURFACES.items():
                p = tissue_dir / filename
                tissue_surfaces[key] = p if p.exists() else None

            tissue_layers = sorted(tissue_dir.glob("Layer_*.vtk"))
            stats_csv = sorted(tissue_dir.glob("*.csv"))

    return PatientCase(
        year=year,
        patient_id=patient_dir.name,
        root=patient_dir,
        lv_dir=lv_dir,
        anatomy=anatomy,
        tissue_surfaces=tissue_surfaces,
        tissue_layers=tissue_layers,
        stats_csv=stats_csv,
    )


def _find_lv_dir(patient_dir: Path) -> Optional[Path]:
    """Resolve .../Basal/Data/DE-MRI/LV/ for a patient folder."""
    lv = patient_dir / "Basal" / "Data" / "DE-MRI" / "LV"
    return lv if lv.is_dir() else None
