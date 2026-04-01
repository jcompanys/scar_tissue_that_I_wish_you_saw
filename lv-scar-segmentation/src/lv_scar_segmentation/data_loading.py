"""Dataset discovery and structure validation utilities.

This module is intentionally lightweight and only handles file discovery:
- years at dataset root
- patient IDs under each year
- expected folder and key file checks for Basal/Data/DE-MRI/LV

No clinical content is parsed here yet.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


EXPECTED_RELATIVE_DIRS: Tuple[str, ...] = (
    "Basal/Data",
    "Basal/Data/DE-MRI",
    "Basal/Data/DE-MRI/LV",
    "Basal/Data/DE-MRI/LV/EAM",
    "Basal/Data/DE-MRI/LV/THICKNESS",
    "Basal/Data/DE-MRI/LV/TISSUE",
    "Basal/Data/DE-MRI/LV/TISSUE_CE",
)

# Accommodates the observed naming mismatch: TRANSMURALITY folder with
# TRANSMURABILITY files in many exports.
TRANSMURAL_DIR_OPTIONS: Tuple[str, ...] = (
    "Basal/Data/DE-MRI/LV/TRANSMURALITY",
    "Basal/Data/DE-MRI/LV/TRANSMURABILITY",
)

EXPECTED_KEY_FILES: Tuple[str, ...] = (
    "Basal/Data/DE-MRI/CardiacModelLVRV.vtk",
    "Basal/Data/DE-MRI/LV/Left Ventricle.vtk",
    "Basal/Data/DE-MRI/LV/Myocardium.vtk",
    "Basal/Data/DE-MRI/LV/Endo Layer.vtk",
    "Basal/Data/DE-MRI/LV/Epi Layer.vtk",
)

TISSUE_CSV_PATTERNS: Tuple[str, ...] = (
    "Basal/Data/DE-MRI/LV/TISSUE/*statistics*.csv",
    "Basal/Data/DE-MRI/LV/TISSUE/*Statistics*.csv",
)


@dataclass(slots=True)
class DataCase:
    """One patient case discovered in the dataset tree."""

    year: str
    patient_id: str
    root: Path
    available_dirs: List[str]
    missing_dirs: List[str]
    missing_key_files: List[str]
    tissue_csv_files: List[Path]

    @property
    def is_valid_minimum(self) -> bool:
        """Return True when minimum required structure is present."""
        return not self.missing_dirs and not self.missing_key_files


@dataclass(slots=True)
class DataStructureReport:
    """Aggregate report for all discovered cases."""

    dataset_root: Path
    total_years: int
    total_cases: int
    valid_cases: int
    invalid_cases: int
    cases: List[DataCase] = field(default_factory=list)

    def to_dict(self) -> Dict[str, object]:
        """Serialize report into plain Python objects."""
        return {
            "dataset_root": str(self.dataset_root),
            "total_years": self.total_years,
            "total_cases": self.total_cases,
            "valid_cases": self.valid_cases,
            "invalid_cases": self.invalid_cases,
            "cases": [
                {
                    "year": c.year,
                    "patient_id": c.patient_id,
                    "root": str(c.root),
                    "is_valid_minimum": c.is_valid_minimum,
                    "available_dirs": c.available_dirs,
                    "missing_dirs": c.missing_dirs,
                    "missing_key_files": c.missing_key_files,
                    "tissue_csv_files": [str(p) for p in c.tissue_csv_files],
                }
                for c in self.cases
            ],
        }


def scan_dataset(dataset_root: Path | str) -> DataStructureReport:
    """Scan dataset root with structure <root>/<year>/<patient_id>/... .

    Parameters
    ----------
    dataset_root:
        Dataset root path, for example:
        F:/RM_TEKNON_DEVELOP

    Returns
    -------
    DataStructureReport
        Full discovery report with per-case validity checks.
    """
    root = Path(dataset_root)
    if not root.exists():
        raise FileNotFoundError(f"Dataset root not found: {root}")

    years = _discover_year_dirs(root)
    cases: List[DataCase] = []

    for year_dir in years:
        for patient_dir in _discover_patient_dirs(year_dir):
            cases.append(_inspect_case(year_dir.name, patient_dir))

    valid_cases = sum(1 for c in cases if c.is_valid_minimum)
    invalid_cases = len(cases) - valid_cases

    return DataStructureReport(
        dataset_root=root,
        total_years=len(years),
        total_cases=len(cases),
        valid_cases=valid_cases,
        invalid_cases=invalid_cases,
        cases=cases,
    )


def _discover_year_dirs(root: Path) -> List[Path]:
    """Return directories whose names look like 4-digit years."""
    return sorted(
        [p for p in root.iterdir() if p.is_dir() and p.name.isdigit() and len(p.name) == 4],
        key=lambda p: p.name,
    )


def _discover_patient_dirs(year_dir: Path) -> List[Path]:
    """Return patient directories under one year directory."""
    return sorted([p for p in year_dir.iterdir() if p.is_dir()], key=lambda p: p.name)


def _inspect_case(year: str, patient_dir: Path) -> DataCase:
    """Inspect one patient folder and validate expected layout."""
    available_dirs: List[str] = []
    missing_dirs: List[str] = []

    for rel in EXPECTED_RELATIVE_DIRS:
        if (patient_dir / rel).is_dir():
            available_dirs.append(rel)
        else:
            missing_dirs.append(rel)

    # TRANSMURALITY / TRANSMURABILITY accepted as either directory name.
    if any((patient_dir / candidate).is_dir() for candidate in TRANSMURAL_DIR_OPTIONS):
        available_dirs.append("Basal/Data/DE-MRI/LV/TRANSMURALITY|TRANSMURABILITY")
    else:
        missing_dirs.append("Basal/Data/DE-MRI/LV/TRANSMURALITY|TRANSMURABILITY")

    missing_key_files: List[str] = []
    for rel_file in EXPECTED_KEY_FILES:
        if not (patient_dir / rel_file).exists():
            missing_key_files.append(rel_file)

    tissue_csv_files = _find_files(patient_dir, TISSUE_CSV_PATTERNS)

    return DataCase(
        year=year,
        patient_id=patient_dir.name,
        root=patient_dir,
        available_dirs=available_dirs,
        missing_dirs=missing_dirs,
        missing_key_files=missing_key_files,
        tissue_csv_files=tissue_csv_files,
    )


def _find_files(root: Path, patterns: Iterable[str]) -> List[Path]:
    """Find files matching a list of glob patterns relative to case root."""
    matches: List[Path] = []
    for pattern in patterns:
        matches.extend(root.glob(pattern))

    # Deduplicate while preserving deterministic ordering.
    unique = sorted(set(matches), key=lambda p: str(p).lower())
    return unique
