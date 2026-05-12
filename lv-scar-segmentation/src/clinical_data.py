"""Clinical registry loader for the LV scar characterization study.

Reads the study's main CSV (one row per patient) and returns a clean
pandas DataFrame with:

  - Python-friendly snake_case column names
  - Correct dtypes: datetime64, float64, Int64 (nullable integer), object
  - Grouped column sets via COLUMN_GROUPS for structured EDA / correlation

The source file is a patient-level clinical registry exported from Excel,
with Spanish column headers, comma decimal separators, and DD/MM/YYYY
dates. The columns are organized around:

  - identification and inclusion
  - MRI scar geometry
  - demographics and risk factors
  - infarction and revascularization history
  - echocardiography
  - devices
  - ventricular arrhythmias
  - ICD therapies / ablation
  - 6-month follow-up outcomes
  - medication

Most yes/no registry fields appear to use binary 0/1 coding. The main
known exception is ``sex`` where the local study notes indicate
``1 = male`` and ``2 = female``.

The loader handles encoding detection, separator inference, duplicate
headers, and the naming quirks observed in the local registry export.
Some fields remain provisional until checked directly against the local
CSV, especially ``revasc_type``, ``enhancement_distribution``,
``enhancement_grade``, ``soo_clinical_va``, ``soo_inducible_va``, and
``notes_misc`` (raw ``@`` column).

Typical usage::

    from src.clinical_data import load, COLUMN_GROUPS, summary

    df = load("path/to/registry.csv")
    summary(df)

    # Correlation: scar geometry vs demographics
    cols = COLUMN_GROUPS["scar_geometry"] + COLUMN_GROUPS["demographics"]
    print(df[cols].corr(numeric_only=True))
"""

from __future__ import annotations

import io
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


# ---------------------------------------------------------------------------
# Column name mapping: raw header (Spanish) → clean Python name
# ---------------------------------------------------------------------------
# Empty / separator columns (unnamed) are dropped automatically.
# Duplicate "FECHA IAM" columns are disambiguated during loading:
#   first occurrence  → date_mi          (date of the index MI)
#   second occurrence → date_mi_history  (date in the clinical-history block)

# The mapping below reflects the current local understanding of the CSV.
# A few semantic interpretations remain provisional and are documented above.
_COLUMN_MAP: Dict[str, str] = {
    # --- Admin / identification ------------------------------------------------
    "INCLUSION":                        "inclusion",
    "ID":                               "patient_id",
    "@":                                "notes_misc",

    # --- MRI session ----------------------------------------------------------
    "FECHA RM":                         "date_mri",
    "TIPO RM":                          "mri_type",
    "DIAGNÓSTICO":                      "diagnosis",

    # --- Scar geometry (MRI-derived) -----------------------------------------
    "LV MASS (g)":                      "lv_mass_g",
    "BZ + CORE (g)":                    "bz_core_g",
    "BZ + CORE (%)":                    "bz_core_pct",
    "BZ (g)":                           "bz_g",
    "BZ (%)":                           "bz_pct",
    "CORE (g)":                         "core_g",
    "CORE (%)":                         "core_pct",
    "CHANNELS":                         "channels",
    "CHANNEL MASS (g)":                 "channel_mass_g",
    "DISTRIBUCIÓN REALCE":              "enhancement_distribution",
    "GRADO REALCE":                     "enhancement_grade",
    "TERRITORIO REALCE":                "enhancement_territory",

    # --- Demographics ---------------------------------------------------------
    "FECHA NAC":                        "date_birth",
    "EDAD":                             "age",
    "SEXO":                             "sex",           # local study notes: 1 = male, 2 = female

    # --- Cardiovascular risk factors ------------------------------------------
    "HTA":                              "hta",
    "DLP":                              "dlp",
    "DM":                               "dm",
    "TABACO":                           "smoking",
    "AFIB":                             "afib",
    "NYHA":                             "nyha",
    "IAM SÍ/NO":                        "mi_yn",

    # --- MI & revascularization history ----------------------------------------
    # First "FECHA IAM" → date_mi;  second → date_mi_history  (pandas mangles to .1)
    "FECHA IAM":                        "date_mi",
    "FECHA IAM.1":                      "date_mi_history",
    "FECHA REVASC":                     "date_revasc",
    "TIPO REVASC":                      "revasc_type",   # provisional meaning until rechecked in CSV
    "NOTAS REVASC":                     "revasc_notes",

    # --- Echocardiography -----------------------------------------------------
    "FEVI":                             "lvef_pct",
    "DTDVI":                            "lvedd_mm",
    "DTSVI":                            "lvesd_mm",
    "VTDVI":                            "lvedv_ml",
    "VTDVIindex":                       "lvedv_index",
    "VTSVI":                            "lvesv_ml",
    "VTSVIindex":                       "lvesv_index",
    "AI":                               "la_mm",
    "SEPTO":                            "septum_mm",
    "PARED POST":                       "post_wall_mm",

    # --- Devices --------------------------------------------------------------
    "PORTADOR DISPOSITIVO":             "device_carrier",
    "FECHA IMPLANTE DISPOSITIVO":       "date_device_implant",
    "TIPO DISPOSITIVO":                 "device_type",
    "TIPO PREVENCIÓN DAI":              "icd_prevention_type",
    "BIOMONITOR":                       "biomonitor",

    # --- Clinical arrhythmias -------------------------------------------------
    "AV CLÍNICAS":                      "clinical_va",
    "FECHA AV CLÍNICA":                 "date_clinical_va",
    "TIPO AV CLÍNICA":                  "type_clinical_va",
    "SOO AV CLÍNICA":                   "soo_clinical_va",
    "CARGA EV":                         "pvc_burden",
    "NOTAS AV":                         "va_notes",
    "AV INDUCIBLES":                    "inducible_va",
    "FECHA EEF INDUCCIÓN":              "date_eef",
    "TIPO AV INDUCIBLES":               "type_inducible_va",
    "SOO AV INDUCIBLE":                 "soo_inducible_va",
    "TIPO CARDIOVERSIÓN AV INDUCIBLE":  "cardioversion_type_inducible_va",
    "TERAPIAS DAI":                     "icd_therapies",
    "FECHA TERAPIA DAI":                "date_icd_therapy",
    "TIPO TERAPIA DAI":                 "type_icd_therapy",
    "ABLACIÓN AV":                      "va_ablation",
    "FECHA ABLACIÓN":                   "date_ablation",
    "NOTAS ABLACIÓN":                   "ablation_notes",
    "OTRAS NOTAS":                      "other_notes",

    # --- 6-month follow-up ----------------------------------------------------
    "6 MONTHS FOLLOW-UP DATE":          "date_6m_followup",
    "ADMISSION FOR HEART FAILURE":      "admission_hf",
    "ADMISSION FOR ANGINA":             "admission_angina",
    "REVASC. (Y/N)":                    "revasc_followup_yn",
    "CLINICAL VA (Y/N)":                "followup_va_yn",
    "CLINICAL VA DATE":                 "date_followup_va",
    "CLINICAL VA TYPE":                 "type_followup_va",
    "CLINICAL VA NOTES":                "notes_followup_va",
    "DEATH (Y/N)":                      "death_yn",
    "CARDIAC DEATH (Y/N)":              "cardiac_death_yn",
    "ARRHYTHMIC DEATH (Y/N)":           "arrhythmic_death_yn",
    "DEATH DATE":                       "date_death",

    # --- Medications (text: name + dose, or "0" when not prescribed) ----------
    "IECAS (Y/N) + DOSE":               "med_acei",
    "ARA 2 (Y/N) + DOSE":              "med_arb",
    "BETA-BLOCKERS + DOSE":             "med_bb",
    "SPIRONOLACTONE (Y/N) + DOSE":      "med_spiro",
    "FUROSEMIDE (Y/N)":                 "med_furosemide",
    "DIGOXIN (Y/N)":                    "med_digoxin",
    "AMIODARONE (Y/N)":                 "med_amiodarone",
    "NYHA CLASS":                       "nyha_class",
}

# ---------------------------------------------------------------------------
# Column type catalogues
# ---------------------------------------------------------------------------

#: Columns that should be parsed as DD/MM/YYYY dates.
DATE_COLS: List[str] = [
    "date_mri",
    "date_birth",
    "date_mi",
    "date_mi_history",
    "date_revasc",
    "date_device_implant",
    "date_clinical_va",
    "date_eef",
    "date_icd_therapy",
    "date_ablation",
    "date_6m_followup",
    "date_followup_va",
    "date_death",
]

#: Columns that should be float64 (comma-separated decimal in source).
FLOAT_COLS: List[str] = [
    "lv_mass_g",
    "bz_core_g",
    "bz_core_pct",
    "bz_g",
    "bz_pct",
    "core_g",
    "core_pct",
    "channel_mass_g",
    "age",
    "lvef_pct",
    "lvedd_mm",
    "lvesd_mm",
    "lvedv_ml",
    "lvedv_index",
    "lvesv_ml",
    "lvesv_index",
    "la_mm",
    "septum_mm",
    "post_wall_mm",
    "pvc_burden",
    "enhancement_grade",
]

#: Coded categorical columns stored as nullable Int64.
CATEGORICAL_CODE_COLS: List[str] = [
    "sex",
]

#: Human-readable labels for coded categorical columns.
CATEGORICAL_LABELS: Dict[str, Dict[int, str]] = {
    "sex": {1: "Male", 2: "Female"},
}

#: Fields whose clinical meaning or coding should be confirmed with doctors.
CLINICAL_REVIEW_NOTES: Dict[str, str] = {
    "risk_factors": (
        "Risk-factor definitions and 0/1 coding should be confirmed with "
        "the clinical team before interpreting prevalence as final."
    ),
    "enhancement_distribution": (
        "Numeric distribution codes need the registry coding legend."
    ),
    "enhancement_grade": (
        "The direction and clinical meaning of the grade should be confirmed."
    ),
}

#: Binary / flag columns stored as nullable Int64 (0/1, possibly NaN).
BINARY_COLS: List[str] = [
    "hta",
    "dlp",
    "dm",
    "smoking",
    "afib",
    "mi_yn",
    "device_carrier",
    "biomonitor",
    "clinical_va",
    "inducible_va",
    "icd_therapies",
    "va_ablation",
    "admission_hf",
    "admission_angina",
    "revasc_followup_yn",
    "followup_va_yn",
    "death_yn",
    "cardiac_death_yn",
    "arrhythmic_death_yn",
]

# ---------------------------------------------------------------------------
# Column groups for EDA
# ---------------------------------------------------------------------------

#: Named column groups for structured access and correlation analysis.
COLUMN_GROUPS: Dict[str, List[str]] = {
    "demographics": [
        "patient_id", "date_birth", "age", "sex",
    ],
    "risk_factors": [
        "hta", "dlp", "dm", "smoking", "afib", "nyha", "mi_yn",
    ],
    "scar_geometry": [
        "lv_mass_g",
        "bz_core_g", "bz_core_pct",
        "bz_g", "bz_pct",
        "core_g", "core_pct",
        "channels", "channel_mass_g",
        "enhancement_distribution", "enhancement_grade", "enhancement_territory",
    ],
    "echo": [
        "lvef_pct",
        "lvedd_mm", "lvesd_mm",
        "lvedv_ml", "lvedv_index",
        "lvesv_ml", "lvesv_index",
        "la_mm", "septum_mm", "post_wall_mm",
    ],
    "devices": [
        "device_carrier", "date_device_implant", "device_type",
        "icd_prevention_type", "biomonitor",
    ],
    "arrhythmias": [
        "clinical_va", "date_clinical_va", "type_clinical_va", "soo_clinical_va",
        "pvc_burden",
        "inducible_va", "date_eef", "type_inducible_va", "soo_inducible_va",
        "icd_therapies", "va_ablation",
    ],
    "followup": [
        "date_6m_followup",
        "admission_hf", "admission_angina",
        "revasc_followup_yn",
        "followup_va_yn", "date_followup_va", "type_followup_va",
        "death_yn", "cardiac_death_yn", "arrhythmic_death_yn", "date_death",
        "nyha_class",
    ],
    "medications": [
        "med_acei", "med_arb", "med_bb", "med_spiro",
        "med_furosemide", "med_digoxin", "med_amiodarone",
    ],
    "mri_session": [
        "date_mri", "mri_type", "diagnosis",
    ],
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def load(
    csv_path: str | Path,
    encoding: Optional[str] = None,
    separator: Optional[str] = None,
) -> pd.DataFrame:
    """Load the clinical registry CSV and return a cleaned DataFrame.

    Parameters
    ----------
    csv_path:
        Path to the registry CSV file.
    encoding:
        Character encoding.  ``None`` triggers auto-detection (tries
        utf-8-sig, utf-8, latin-1, cp1252 in order).
    separator:
        Column delimiter.  ``None`` triggers auto-detection (tries
        ``\\t``, ``;``, ``|``, ``,``).

    Returns
    -------
    pd.DataFrame
        One row per patient with clean column names and correct dtypes.
        Index is reset to 0-based integers; ``patient_id`` is preserved
        as a regular column.
    """
    path = Path(csv_path)
    if not path.exists():
        raise FileNotFoundError(f"Registry CSV not found: {path}")

    raw_text = _read_text(path, encoding)
    sep = separator or _infer_separator(raw_text)
    # pandas ≥2.0 removed mangle_dupe_cols; duplicate headers are renamed
    # automatically with .1/.2 suffixes, so no extra flag is needed.
    df = pd.read_csv(
        io.StringIO(raw_text),
        sep=sep,
        dtype=str,          # read everything as text first; type-cast below
        keep_default_na=False,
        na_values=["", " ", "N/A", "NA", "n/a", "-"],
    )

    df = _drop_unnamed(df)
    df = _rename_columns(df)
    df = _cast_types(df)
    df = df.reset_index(drop=True)
    return df


def summary(df: pd.DataFrame) -> None:
    """Print a structured overview of the loaded DataFrame.

    Covers completeness, scar geometry statistics, demographic breakdown,
    and echo parameters.
    """
    n = len(df)
    print(f"=== Clinical registry: {n} patients ===\n")

    # -- Completeness per group -----------------------------------------------
    print("Completeness by column group:")
    for group, cols in COLUMN_GROUPS.items():
        present = [c for c in cols if c in df.columns]
        if not present:
            continue
        non_null = df[present].notna().mean().mean() * 100
        print(f"  {group:<20} {len(present):>2} cols  {non_null:5.1f}% non-null")

    # -- Demographics ---------------------------------------------------------
    print("\nDemographics:")
    if "age" in df.columns:
        print(f"  Age          mean={df['age'].mean():.1f}  "
              f"sd={df['age'].std():.1f}  "
              f"[{df['age'].min():.0f}–{df['age'].max():.0f}]")
    if "sex" in df.columns:
        counts = df["sex"].value_counts(dropna=False)
        print(f"  Sex (raw)    {dict(counts)}")

    # -- Scar geometry --------------------------------------------------------
    print("\nScar geometry (non-null patients):")
    scar_num = [c for c in ["lv_mass_g", "bz_core_g", "bz_core_pct",
                             "core_g", "core_pct", "bz_g", "bz_pct"]
                if c in df.columns]
    if scar_num:
        print(df[scar_num].describe().to_string())

    # -- Echo -----------------------------------------------------------------
    print("\nEchocardiography (non-null patients):")
    echo_num = [c for c in COLUMN_GROUPS["echo"] if c in df.columns]
    if echo_num:
        print(df[echo_num].describe().to_string())

    print()


def available_columns(df: pd.DataFrame) -> Dict[str, List[str]]:
    """Return COLUMN_GROUPS filtered to columns actually present in *df*."""
    return {
        group: [c for c in cols if c in df.columns]
        for group, cols in COLUMN_GROUPS.items()
    }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_ENCODINGS = ("utf-8-sig", "utf-8", "latin-1", "cp1252")
_SEPARATORS = ("\t", ";", "|", ",")


def _read_text(path: Path, encoding: Optional[str]) -> str:
    """Read file text, trying encodings in order if encoding is None."""
    candidates = [encoding] if encoding else list(_ENCODINGS)
    for enc in candidates:
        try:
            return path.read_text(encoding=enc)
        except (UnicodeDecodeError, LookupError):
            continue
    raise ValueError(
        f"Could not decode {path} with any of: {candidates}"
    )


def _infer_separator(text: str) -> str:
    """Pick the delimiter that produces the most columns on the first line."""
    first_line = text.split("\n", 1)[0]
    best_sep = "\t"
    best_count = 0
    for sep in _SEPARATORS:
        count = first_line.count(sep)
        if count > best_count:
            best_count = count
            best_sep = sep
    return best_sep


def _drop_unnamed(df: pd.DataFrame) -> pd.DataFrame:
    """Drop columns with empty or pandas-generated 'Unnamed: N' headers."""
    keep = [
        c for c in df.columns
        if not (str(c).startswith("Unnamed:") or str(c).strip() == "")
    ]
    return df[keep]


def _rename_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Apply _COLUMN_MAP; strip whitespace from all headers first."""
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    return df.rename(columns=_COLUMN_MAP)


def _cast_types(df: pd.DataFrame) -> pd.DataFrame:
    """Apply dtype conversions: dates, floats, nullable ints, channels."""
    df = df.copy()

    # Dates ----------------------------------------------------------------
    for col in DATE_COLS:
        if col in df.columns:
            df[col] = pd.to_datetime(df[col], format="%d/%m/%Y", errors="coerce")

    # Floats ---------------------------------------------------------------
    for col in FLOAT_COLS:
        if col in df.columns:
            # Source uses comma as decimal; already handled by read_csv(decimal=",")
            # but values read as str may still need coercion.
            df[col] = pd.to_numeric(
                df[col].astype(str).str.replace(",", ".", regex=False),
                errors="coerce",
            )

    # Binary / nullable int ------------------------------------------------
    for col in BINARY_COLS:
        if col in df.columns:
            df[col] = (
                pd.to_numeric(df[col], errors="coerce")
                .astype("Int64")
            )

    # Coded categorical / nullable int ------------------------------------
    for col in CATEGORICAL_CODE_COLS:
        if col in df.columns:
            df[col] = (
                pd.to_numeric(df[col], errors="coerce")
                .astype("Int64")
            )

    # CHANNELS: "NO" / "SÍ" / "SI" / "YES" → 0/1 Int64 -------------------
    if "channels" in df.columns:
        ch = df["channels"].astype(str).str.strip().str.upper()
        mapping = {"NO": 0, "SI": 1, "SÍ": 1, "YES": 1}
        numeric = pd.to_numeric(df["channels"], errors="coerce")
        text_mapped = ch.map(mapping)
        df["channels"] = numeric.fillna(text_mapped).astype("Int64")

    # NYHA: may be written "I", "II", "III", "IV" or 1–4 -----------------
    for col in ("nyha", "nyha_class"):
        if col in df.columns:
            df[col] = _parse_nyha(df[col])

    return df


_NYHA_MAP = {"I": 1, "II": 2, "III": 3, "IV": 4}


def _parse_nyha(series: pd.Series) -> pd.Series:
    """Convert Roman or Arabic NYHA values to nullable Int64."""
    roman = series.astype(str).str.strip().str.upper().map(_NYHA_MAP)
    numeric = pd.to_numeric(series, errors="coerce")
    return numeric.fillna(roman).astype("Int64")
