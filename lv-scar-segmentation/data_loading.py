"""ADAS3D dataset loader for RM_TEKNON_DEVELOP.

Dataset layout on the external HD:
    <root>/<year>/<patient_id>/Basal/Data/DE-MRI/LV/
        Endo Layer.vtk
        Epi Layer.vtk
        Left Ventricle.vtk
        Myocardium.vtk
        TISSUE/
            Core Surface.vtk          (or variant names — matched by glob)
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
    from data_loading import scan, diagnose

    cases = scan(r"F:/RM_TEKNON_DEVELOP")
    for c in cases:
        print(c)

    # inspect what files are actually present for one case:
    diagnose(cases[0])
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ── Anatomy: exact filenames (stable across ADAS3D versions) ──────────────────

_ANATOMY_FILES: Dict[str, str] = {
    "endo":       "Endo Layer.vtk",
    "epi":        "Epi Layer.vtk",
    "lv":         "Left Ventricle.vtk",
    "myocardium": "Myocardium.vtk",
}

# ── Tissue surfaces: glob patterns (filenames vary across ADAS3D exports) ──────
#
# Each key maps to a tuple of glob patterns tried in order.
# The first file that matches wins.  Patterns are case-insensitive on Windows
# (the default NTFS filesystem ignores case).

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
    "healthy": (
        "Healthy Surface.vtk",
        "*Healthy*Surface*.vtk",
        "*Healthy*.vtk",
    ),
    "scar": (
        "Scar Surface.vtk",
        "*Scar*Surface*.vtk",
        "*Scar*.vtk",
    ),
}


# ── Data model ─────────────────────────────────────────────────────────────────

@dataclass
class PatientCase:
    year: str
    patient_id: str
    root: Path                             # <root>/<year>/<patient_id>/
    lv_dir: Optional[Path]                 # .../Basal/Data/DE-MRI/LV/

    anatomy: Dict[str, Optional[Path]]        = field(default_factory=dict)
    tissue_surfaces: Dict[str, Optional[Path]] = field(default_factory=dict)
    tissue_layers: List[Path]                  = field(default_factory=list)
    stats_csv: List[Path]                      = field(default_factory=list)

    def __str__(self) -> str:
        n_anatomy = sum(1 for v in self.anatomy.values() if v is not None)
        n_tissue  = sum(1 for v in self.tissue_surfaces.values() if v is not None)
        return (
            f"[{self.year}] {self.patient_id} | "
            f"anatomy {n_anatomy}/{len(_ANATOMY_FILES)}  "
            f"tissue {n_tissue}/{len(_TISSUE_SURFACE_GLOBS)}  "
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


def diagnose(case: PatientCase) -> None:
    """Print all files found inside this patient's TISSUE folder.

    Call this when the scan finds fewer cases than expected — it reveals
    the actual filenames on disk so you can adjust the glob patterns.
    """
    print(f"\n{'─'*60}")
    print(f"Patient : {case.year}/{case.patient_id}")
    print(f"LV dir  : {case.lv_dir}")

    if case.lv_dir is None:
        print("  ⚠  LV directory not found — check Basal/Data/DE-MRI/LV path")
        return

    tissue_dir = case.lv_dir / "TISSUE"
    if not tissue_dir.is_dir():
        print(f"  ⚠  TISSUE folder not found at {tissue_dir}")
        return

    vtk_files = sorted(tissue_dir.glob("*.vtk"))
    print(f"\n  TISSUE/ contains {len(vtk_files)} VTK files:")
    for f in vtk_files:
        matched_as = _which_key(f.name)
        tag = f"  ← matched as '{matched_as}'" if matched_as else ""
        print(f"    {f.name}{tag}")

    print(f"\n  Resolved tissue_surfaces:")
    for key, path in case.tissue_surfaces.items():
        status = str(path.name) if path else "NOT FOUND"
        print(f"    {key:12s} → {status}")
    print(f"{'─'*60}\n")


# ── Internal helpers ───────────────────────────────────────────────────────────

def _year_dirs(root: Path) -> List[Path]:
    return sorted(
        p for p in root.iterdir()
        if p.is_dir() and p.name.isdigit() and len(p.name) == 4
    )


def _find_tissue_surface(tissue_dir: Path,
                          globs: Tuple[str, ...]) -> Optional[Path]:
    """Return the first file in tissue_dir that matches any of the glob patterns."""
    for pattern in globs:
        matches = sorted(tissue_dir.glob(pattern))
        if matches:
            return matches[0]
    return None


def _which_key(filename: str) -> Optional[str]:
    """Reverse-lookup: which tissue key would this filename match?"""
    for key, globs in _TISSUE_SURFACE_GLOBS.items():
        for pattern in globs:
            # simple substring check for the diagnosis display
            stem = pattern.replace("*", "").replace(".vtk", "").lower()
            if stem and stem in filename.lower():
                return key
    return None


def _load_case(year: str, patient_dir: Path) -> PatientCase:
    lv_dir = _find_lv_dir(patient_dir)

    anatomy: Dict[str, Optional[Path]] = {}
    tissue_surfaces: Dict[str, Optional[Path]] = {}
    tissue_layers: List[Path] = []
    stats_csv: List[Path] = []

    if lv_dir is not None:
        # anatomy — exact names (stable)
        for key, filename in _ANATOMY_FILES.items():
            p = lv_dir / filename
            anatomy[key] = p if p.exists() else None

        # tissue surfaces — flexible glob matching
        tissue_dir = lv_dir / "TISSUE"
        if tissue_dir.is_dir():
            for key, globs in _TISSUE_SURFACE_GLOBS.items():
                tissue_surfaces[key] = _find_tissue_surface(tissue_dir, globs)

            tissue_layers = sorted(tissue_dir.glob("Layer_*.vtk"))
            stats_csv     = sorted(tissue_dir.glob("*.csv"))

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
