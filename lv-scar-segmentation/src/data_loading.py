"""ADAS3D dataset loader for RM_TEKNON_DEVELOP.

Finds, for every patient on the external drive, the paths of the ADAS3D
meshes the notebooks need (LV anatomy, RV anatomy, scar tissue surfaces and
corridor outputs). Only paths are resolved here; meshes are read later,
mostly through ``src.mesh_utils``.

Main API
--------
scan(dataset_root) -> list[PatientCase]
    One ``PatientCase`` per patient folder, in sorted (year, patient) order.
PatientCase
    Resolved paths for one patient. Missing files are stored as ``None``.
DE_MRI_VARIANTS
    DE-MRI folder names tried for each patient, in priority order.

Used by
-------
``scan``: F01, M01, M02_clinical_eda, exploratory/M02_isomaps, M03,
M04_scar_analysis, M04_BZ_SCC, M05. ``DE_MRI_VARIANTS``: M01, M06.
``PatientCase`` fields are read by the notebooks and by ``src.mesh_utils``.

Dataset layout
--------------
::

    <root>/<year>/<patient_id>/Basal/Data/<de_mri_variant>/
        LV/
            Endo Layer.vtk
            Epi Layer.vtk
            Left Ventricle.vtk
            Myocardium.vtk
            TISSUE/                   (searched recursively)
                Core Surface.vtk
                Border Zone Surface.vtk
                Scar Surface.vtk
                Automatic 3D Corridors.vtk
                3D Corridor Labels.vtk
                3D Corridor Unified Centerlines.vtk
                3D Corridor 3D Border Zone.vtk
        RV/  (or "Right Ventricle/")
            Right Ventricle.vtk, RV Endo Layer.vtk, RV Epi Layer.vtk, ...

``<de_mri_variant>`` is one of ``DE_MRI_VARIANTS``. Tissue and RV file names
vary across ADAS3D versions, so they are matched with ordered glob patterns
(see ``_TISSUE_SURFACE_GLOBS`` and ``_RV_ANATOMY_GLOBS``).

Usage
-----
::

    from src.data_loading import scan

    cases = scan(r"F:/RM_TEKNON_DEVELOP")
    core_path = cases[0].tissue_surfaces.get("core")
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# -- LV anatomy: exact filenames (stable across ADAS3D versions) --------------
#
# Keys are the ones ``mesh_utils.SHELL_PRIORITY`` asks for.

_ANATOMY_FILES: Dict[str, str] = {
    "endo":       "Endo Layer.vtk",
    "epi":        "Epi Layer.vtk",
    "lv":         "Left Ventricle.vtk",
    "myocardium": "Myocardium.vtk",
}

# -- RV anatomy: glob patterns tried in order, searched recursively -----------
#
# RV exports are less consistent than LV ones, so each key lists the exact
# name first and broader patterns after it. Keys match
# ``mesh_utils.RV_PRIORITY``.
#
# Flagged (not fixed): the last "rv_endo" / "rv_epi" patterns
# ("*Endo*Layer*.vtk", "*Epi*Layer*.vtk") would also match LV layer files.
# Safe only because the search is restricted to the RV folder.

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
# are found as well as files directly in TISSUE/. Patterns go from the exact
# ADAS3D name to broader ones; the first pattern with any match wins.
#
# Flagged (not fixed): broad fallbacks can pick the wrong file when the exact
# name is missing.
#   - "border_zone": "*Border*Zone*.vtk" also matches
#     "3D Corridor 3D Border Zone.vtk" (the corridor BZ, not the tissue BZ).
#   - "scar": "*Scar*.vtk" matches any file with "Scar" in its name.
# No record is kept of which pattern matched, so it is not known whether
# these fallbacks ever fired on the dataset.

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
#
# Some patients have more than one DE-MRI acquisition. The plain "DE-MRI"
# export is preferred, then the 3D one, then the 2D one.

DE_MRI_VARIANTS = ("DE-MRI", "DE-MRI 3D", "DE-MRI 2D")


# -- Data model ----------------------------------------------------------------

@dataclass
class PatientCase:
    """Resolved file paths for one patient.

    Attributes
    ----------
    year : str
        Four-digit year folder name, e.g. ``"2023"``.
    patient_id : str
        Patient folder name, kept as a string (IDs have leading zeros).
    root : Path
        Patient folder ``<dataset_root>/<year>/<patient_id>``.
    lv_dir : Path or None
        ``Basal/Data/<variant>/LV`` folder, ``None`` if no variant has one.
    anatomy : dict of str to Path or None
        LV anatomy meshes, keys ``"endo"``, ``"epi"``, ``"lv"``,
        ``"myocardium"``. Empty if ``lv_dir`` is ``None``.
    rv_anatomy : dict of str to Path or None
        RV anatomy meshes, keys ``"rv"``, ``"rv_endo"``, ``"rv_epi"``.
        Empty if no RV folder is found.
    tissue_surfaces : dict of str to Path or None
        Meshes from ``LV/TISSUE``, keys of ``_TISSUE_SURFACE_GLOBS``
        (``"core"``, ``"border_zone"``, ``"scar"``, corridor outputs).
        Empty if there is no TISSUE folder.

    Notes
    -----
    The LV has no ``lv_anatomy`` field: ``anatomy`` is the LV anatomy (the
    name predates RV support). A value of ``None`` means the file was not
    found; callers check this before reading.
    """

    year: str
    patient_id: str
    root: Path
    lv_dir: Optional[Path]

    anatomy: Dict[str, Optional[Path]]         = field(default_factory=dict)
    rv_anatomy: Dict[str, Optional[Path]]      = field(default_factory=dict)
    tissue_surfaces: Dict[str, Optional[Path]] = field(default_factory=dict)


# -- Public API ----------------------------------------------------------------

def scan(dataset_root: str | Path) -> List[PatientCase]:
    """Scan the dataset and return one ``PatientCase`` per patient.

    Parameters
    ----------
    dataset_root : str or Path
        Dataset root containing the year folders, e.g.
        ``"F:/RM_TEKNON_DEVELOP"``.

    Returns
    -------
    list of PatientCase
        Every patient folder, including those with missing files. Sorted by
        year, then by patient folder name.

    Raises
    ------
    FileNotFoundError
        If ``dataset_root`` does not exist.

    Notes
    -----
    The order is deterministic, so case indices and random subsets drawn
    with a fixed seed in the notebooks are reproducible.
    """
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
    """Return the sorted four-digit year folders under ``root``.

    Other folders at the root (exports, notes) are ignored.
    """
    return sorted(
        p for p in root.iterdir()
        if p.is_dir() and p.name.isdigit() and len(p.name) == 4
    )


def _find_surface(root: Path, globs: Tuple[str, ...]) -> Optional[Path]:
    """Return the first file under ``root`` matching one of ``globs``.

    Parameters
    ----------
    root : Path
        Folder searched recursively.
    globs : tuple of str
        Glob patterns, most specific first.

    Returns
    -------
    Path or None
        Alphabetically first match of the first pattern with any match,
        ``None`` if no pattern matches. Sorting makes the pick independent
        of filesystem listing order.
    """
    for pattern in globs:
        matches = sorted(root.rglob(pattern))
        if matches:
            return matches[0]
    return None


def _load_case(year: str, patient_dir: Path) -> PatientCase:
    """Resolve all paths for one patient folder.

    Parameters
    ----------
    year : str
        Year folder name.
    patient_dir : Path
        Patient folder.

    Returns
    -------
    PatientCase
        Paths found; dicts stay empty when their parent folder is missing.
    """
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
    """Return ``Basal/Data/<variant>/LV`` for the first variant that has it."""
    base = patient_dir / "Basal" / "Data"
    for variant in DE_MRI_VARIANTS:
        lv = base / variant / "LV"
        if lv.is_dir():
            return lv
    return None


def _find_rv_dir(patient_dir: Path) -> Optional[Path]:
    """Return the RV folder (``RV`` or ``Right Ventricle``) of the first variant.

    Notes
    -----
    Flagged (not fixed): the RV folder is searched independently of the LV
    one, so LV and RV may come from different DE-MRI variants when a
    patient has several acquisitions.
    """
    base = patient_dir / "Basal" / "Data"
    for variant in DE_MRI_VARIANTS:
        variant_dir = base / variant
        for name in ("RV", "Right Ventricle"):
            rv = variant_dir / name
            if rv.is_dir():
                return rv
    return None
