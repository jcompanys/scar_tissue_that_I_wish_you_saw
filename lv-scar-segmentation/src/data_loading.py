"""ADAS3D dataset loader for RM_TEKNON_DEVELOP.

Dataset layout on the external HD:
    <root>/<year>/<patient_id>/Basal/Data/<de_mri_variant>/LV/
        Endo Layer.vtk
        Epi Layer.vtk
        Left Ventricle.vtk
        Myocardium.vtk
        TISSUE/
            Core Surface.vtk          (exact name varies across ADAS3D versions)
            Border Zone Surface.vtk   (matched via rglob + ordered patterns)
            Healthy Surface.vtk
            Scar Surface.vtk
            Power Paths.vtk           (newer exports: critical corridor isthmuses)
            Layer_10.vtk ... Layer_90.vtk
            *.csv
        EAM/
        THICKNESS/
        TISSUE_CE/
        TRANSMURALITY/   (or TRANSMURABILITY/ in some exports)

    <de_mri_variant> is one of: "DE-MRI", "DE-MRI 3D", "DE-MRI 2D"
    (tried in that order; first match wins)

Usage:
    from src.data_loading import scan

    cases = scan(r"F:/RM_TEKNON_DEVELOP")
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# -- Anatomy: exact filenames (stable across ADAS3D versions) -----------------

_ANATOMY_FILES: Dict[str, str] = {
    "endo":       "Endo Layer.vtk",
    "epi":        "Epi Layer.vtk",
    "lv":         "Left Ventricle.vtk",
    "myocardium": "Myocardium.vtk",
}

_RV_ANATOMY_GLOBS: Dict[str, Tuple[str, ...]] = {
    "rv": (
        "Right Ventricle.vtk",
        "Right Ventricle Surface.vtk",
        "*Right*Ventricle*.vtk",
        "RV.vtk",
        "RV Surface.vtk",
    ),
    "rv_endo": (
        "RV Endo Layer.vtk",
        "*RV*Endo*.vtk",
        "*Endo*Layer*.vtk",
    ),
    "rv_epi": (
        "RV Epi Layer.vtk",
        "*RV*Epi*.vtk",
        "*Epi*Layer*.vtk",
    ),
}

# -- Tissue surfaces: glob patterns tried in order, searched recursively -------
#
# rglob is used so files inside sub-folders (e.g. TISSUE/Core Surface/*.vtk)
# are found as well as files directly in TISSUE/.

_TISSUE_SURFACE_GLOBS: Dict[str, Tuple[str, ...]] = {
    "core": (
        "Core Surface.vtk",
        "*Core*Surface*.vtk",
        "*Core*Zone*.vtk",
        "*CZ*.vtk",
    ),
    "border_zone": (
        "Border Zone Surface.vtk",
        "*Border*Zone*Surface*.vtk",
        "*Border*Zone*.vtk",
        "*BZ*.vtk",
    ),
    "scar": (
        "Scar Surface.vtk",
        "*Scar*Surface*.vtk",
        "*Scar*.vtk",
    ),
    "corridors": (
        "Automatic 3D Corridors.vtk",
        "*Automatic*Corridor*.vtk",
        "*Automatic*3D*Corridor*.vtk",
    ),
    "corridor_labels": (
        "3D Corridor Labels.vtk",
        "*Corridor*Labels*.vtk",
    ),
    "corridor_centerlines_unified": (
        "3D Corridor Unified Centerlines.vtk",
        "*Corridor*Unified*Centerlines*.vtk",
    ),
    "corridor_bz": (
        "3D Corridor 3D Border Zone.vtk",
        "*Corridor*3D*Border*Zone*.vtk",
        "*Corridor*Border*Zone*.vtk",
    ),
}

# -- DE-MRI folder name variants (tried in order, first match wins) -----------

DE_MRI_VARIANTS = ("DE-MRI", "DE-MRI 3D", "DE-MRI 2D")


# -- Data model ----------------------------------------------------------------

@dataclass
class PatientCase:
    year: str
    patient_id: str
    root: Path
    lv_dir: Optional[Path]

    anatomy: Dict[str, Optional[Path]]         = field(default_factory=dict)
    rv_anatomy: Dict[str, Optional[Path]]      = field(default_factory=dict)
    tissue_surfaces: Dict[str, Optional[Path]] = field(default_factory=dict)


# -- Public API ----------------------------------------------------------------

def scan(dataset_root: str | Path) -> List[PatientCase]:
    """Scan the full dataset and return one PatientCase per patient."""
    root = Path(dataset_root)
    if not root.exists():
        raise FileNotFoundError(f"Dataset root not found: {root}")

    cases: List[PatientCase] = []
    for year_dir in _year_dirs(root):
        for patient_dir in sorted(p for p in year_dir.iterdir() if p.is_dir()):
            cases.append(_load_case(year_dir.name, patient_dir))
    return cases


# -- Internal helpers ----------------------------------------------------------

def _year_dirs(root: Path) -> List[Path]:
    return sorted(
        p for p in root.iterdir()
        if p.is_dir() and p.name.isdigit() and len(p.name) == 4
    )


def _find_surface(root: Path, globs: Tuple[str, ...]) -> Optional[Path]:
    """Search root recursively for the first file matching any pattern."""
    for pattern in globs:
        matches = sorted(root.rglob(pattern))
        if matches:
            return matches[0]
    return None


def _load_case(year: str, patient_dir: Path) -> PatientCase:
    lv_dir = _find_lv_dir(patient_dir)
    rv_dir = _find_rv_dir(patient_dir)

    anatomy: Dict[str, Optional[Path]] = {}
    rv_anatomy: Dict[str, Optional[Path]] = {}
    tissue_surfaces: Dict[str, Optional[Path]] = {}

    if lv_dir is not None:
        for key, filename in _ANATOMY_FILES.items():
            p = lv_dir / filename
            anatomy[key] = p if p.exists() else None

        tissue_dir = lv_dir / "TISSUE"
        if tissue_dir.is_dir():
            for key, globs in _TISSUE_SURFACE_GLOBS.items():
                tissue_surfaces[key] = _find_surface(tissue_dir, globs)

    if rv_dir is not None:
        for key, globs in _RV_ANATOMY_GLOBS.items():
            rv_anatomy[key] = _find_surface(rv_dir, globs)

    return PatientCase(
        year=year,
        patient_id=patient_dir.name,
        root=patient_dir,
        lv_dir=lv_dir,
        anatomy=anatomy,
        rv_anatomy=rv_anatomy,
        tissue_surfaces=tissue_surfaces,
    )


def _find_lv_dir(patient_dir: Path) -> Optional[Path]:
    base = patient_dir / "Basal" / "Data"
    for variant in DE_MRI_VARIANTS:
        lv = base / variant / "LV"
        if lv.is_dir():
            return lv
    return None


def _find_rv_dir(patient_dir: Path) -> Optional[Path]:
    base = patient_dir / "Basal" / "Data"
    for variant in DE_MRI_VARIANTS:
        variant_dir = base / variant
        for name in ("RV", "Right Ventricle"):
            rv = variant_dir / name
            if rv.is_dir():
                return rv
    return None
