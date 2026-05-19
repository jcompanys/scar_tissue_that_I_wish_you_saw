# Data and Cohort Description

This chapter describes the clinical and imaging-derived dataset used throughout this project, the selection criteria applied to obtain the analysis cohort, and the data structures available for each patient case. All cases originate from the DEVELOP registry, a prospective collection of post-myocardial infarction patients studied at Centro Medico Teknon (Barcelona, Spain) in the context of ventricular tachycardia (VT) risk assessment.

The chapter is organised as follows. Section 2.1 introduces the source registry and the available data structure. Section 2.2 details the cohort selection process and the resulting exclusion cascade. Section 2.3 describes the data partitioning strategy adopted for modelling. Section 2.4 summarises the anatomical, tissue and clinical representations available for each patient.

## 2.1 Source Data and Available File Structure

The source data consist of anonymised cardiac magnetic resonance (CMR) studies from the DEVELOP registry. The original DICOM acquisitions are stored in Teknon's internal storage system. The working dataset used in this project does not operate directly on the DICOM files; instead, it uses the three-dimensional anatomical and tissue meshes exported after post-processing.

For each patient, image post-processing was carried out using ADAS 3D (Galgo Medical, Barcelona), a software platform for three-dimensional cardiac image analysis. ADAS 3D reconstructs the left ventricular (LV) shell from the late gadolinium enhancement CMR study and performs semi-automatic tissue characterisation based on signal intensity. The software distinguishes dense scar, also referred to as core zone (CZ), heterogeneous peri-infarct tissue, also referred to as border zone (BZ), and healthy myocardium. The resulting surfaces are exported as VTK meshes, which constitute the primary geometric input to the analyses in this work.

The local dataset is organised by acquisition year and patient identifier. The expected structure used by the loading code is:

```text
<dataset_root>/<year>/<patient_id>/Basal/Data/<de_mri_variant>/LV/
```

where `<de_mri_variant>` is searched in the following order: `DE-MRI`, `DE-MRI 3D`, and `DE-MRI 2D`. Within the `LV` directory, the loader searches for the following anatomical meshes:

| File | Meaning |
|---|---|
| `Endo Layer.vtk` | LV endocardial surface |
| `Epi Layer.vtk` | LV epicardial surface |
| `Left Ventricle.vtk` | LV surface export |
| `Myocardium.vtk` | LV myocardial wall surface |

The `LV/TISSUE/` subdirectory contains the tissue and corridor exports. File names vary slightly between ADAS 3D versions, so the pipeline searches recursively using ordered filename patterns rather than relying on one exact name. The main expected tissue files are:

| Tissue export | Main filename or pattern |
|---|---|
| Core zone / dense scar | `Core Surface.vtk`, `*Core*Surface*.vtk`, `*Core*Zone*.vtk`, `*CZ*.vtk` |
| Border zone | `Border Zone Surface.vtk`, `*Border*Zone*Surface*.vtk`, `*Border*Zone*.vtk`, `*BZ*.vtk` |
| Healthy myocardium | `Healthy Surface.vtk`, `*Healthy*Surface*.vtk`, `*Healthy*.vtk` |
| Total scar | `Scar Surface.vtk`, `*Scar*Surface*.vtk`, `*Scar*.vtk` |
| Corridor paths | `Power Paths.vtk`, `*Power*Path*.vtk`, `*PowerPath*.vtk` |
| Transmural layers | `Layer_*.vtk` |
| Tissue statistics | `*.csv` files inside `TISSUE/` |

Some ADAS exports may also contain additional folders under the same CMR variant, such as `EAM/`, `THICKNESS/`, `TISSUE_CE/`, and `TRANSMURALITY/`. In some cases, the transmurality folder appears with the spelling `TRANSMURABILITY/`. These auxiliary folders were not the central input for the present project, but they indicate that the export may contain richer ADAS 3D-derived information than the core CZ mesh alone.

The loader can also identify right ventricular exports when present. It searches for an `RV/` or `Right Ventricle/` directory under the same CMR variant and then looks for files matching names such as `Right Ventricle.vtk`, `Right Ventricle Surface.vtk`, `RV.vtk`, `RV Endo Layer.vtk`, and `RV Epi Layer.vtk`. However, the current project is restricted to LV scar geometry.

## 2.2 Cohort Selection and Exclusion Criteria

The scope of this project is the geometric characterisation of dense scar within the left ventricle. A case was therefore considered valid only if it contained both:

| Requirement | Operational definition |
|---|---|
| LV reconstruction | A valid `Basal/Data/<de_mri_variant>/LV/` directory |
| Core zone mesh | A CZ surface file found inside `LV/TISSUE/` |

An automated scan of the dataset was performed across all patient folders. The scan first verified whether an LV reconstruction directory existed, then checked whether a core zone surface could be resolved inside the corresponding `TISSUE` directory.

Table 2.1 summarises the scan results by acquisition year. Of the 168 total patient folders, 118 (70.2%) contained a valid LV reconstruction directory, and 105 (62.5%) additionally contained a CZ mesh. These 105 cases constitute the valid cohort used in the subsequent geometric analysis. The patient-level list of excluded cases is provided in Annex A.

**Table 2.1. Dataset overview by acquisition year.**

| Year | Total | Has LV | Has CZ | Yield |
|---|---:|---:|---:|---:|
| 2021 | 48 | 12 | 12 | 25.0% |
| 2022 | 27 | 23 | 23 | 85.2% |
| 2023 | 47 | 45 | 32 | 68.1% |
| 2024 | 31 | 31 | 31 | 100.0% |
| 2025 | 15 | 7 | 7 | 46.7% |
| Total | 168 | 118 | 105 | 62.5% |

The 63 excluded cases fall into two operational categories. Fifty cases lacked an LV reconstruction directory. Among these, 35 had a `Basal/Data` path but no recognised DE-MRI folder among `DE-MRI`, `DE-MRI 3D`, or `DE-MRI 2D`; the remaining 15 had no `Basal/Data` path at all in the expected local folder structure. This does not imply that the original DICOM studies do not exist, since the DICOM acquisitions are kept in Teknon's internal storage. It only means that the corresponding processed ADAS 3D export was not available in the working dataset used here.

The remaining 13 excluded cases had an LV reconstruction but no usable CZ mesh. In 12 cases, the `TISSUE` directory existed but contained no VTK files. In one case, patient `96024662`, the `TISSUE` directory contained corridor-related meshes, including `3D Corridor Centerlines.vtk`, `3D Corridor Filtered Centerlines.vtk`, and `Automatic 3D Corridors.vtk`, but no file matching the CZ surface patterns. This case was therefore excluded because the dense scar surface required by this project could not be identified.

## 2.3 Data Partitioning

The 105 valid cases were split into two disjoint subsets using a random permutation with a fixed seed (`seed = 42`) to ensure reproducibility. Eighty percent of the cases (`n = 84`) were assigned to the known set, used for exploratory analysis, descriptor extraction, and model fitting. The remaining 20% (`n = 21`) were reserved as a held-out set for future validation of generative or predictive models. The split was serialised to `results/split.json`, ensuring that subsequent pipeline stages operate on identical partitions.

No stratification by clinical variables was applied at this stage. The main objective of this work is methodological: to build and evaluate a geometric representation of LV scar morphology. Here, "prognostic" would mean using the data to predict future clinical outcomes, such as arrhythmia recurrence, ICD therapy, hospital admission, or death. Since this project does not train an outcome-prediction model, the split was defined only at the case level and not balanced by outcome labels.

**Table 2.2. Data partition summary.**

| Subset | n | Proportion |
|---|---:|---:|
| Known set | 84 | 80% |
| Held-out set | 21 | 20% |
| Total valid | 105 | 100% |

Clinical metadata, including patient sex, were loaded from the DEVELOP registry CSV and linked to the imaging cases by patient identifier. Sex information was available for 94 of the 105 valid patients. Among these, 79 patients (84.0%) were coded as male and 15 patients (16.0%) as female, using the registry encoding `1 = male` and `2 = female`.

This male-skewed distribution is plausible for a post-myocardial infarction cohort. Epidemiological studies report a higher incidence of myocardial infarction in men than in women, although women may have different symptom profiles and, in some contexts, worse outcomes after infarction. For example, the Tromso Study reported a lifelong higher risk of incident myocardial infarction in men compared with women, and clinical summaries such as StatPearls note male sex as a non-modifiable MI risk factor, with men tending to experience MI earlier in life. This project therefore treats sex as contextual cohort information, not as a modelling input.

## 2.4 Available Data per Patient

For each valid case, the following data are available and constitute the inputs or contextual information for the geometric analysis pipeline described in later chapters.

### Left Ventricular Shell

The LV shell consists of triangular surface meshes representing the endocardial and epicardial boundaries of the left ventricle, together with additional LV or myocardial wall exports when available. These meshes define the anatomical reference frame within which scar geometry is described. In the visualisation pipeline, the best available shell mesh is rendered as a semi-transparent surface to provide spatial context for the scar geometry.

### Tissue Classification Surfaces

The `TISSUE` subdirectory contains VTK surface meshes for the tissue classes identified by ADAS 3D. The primary surface used in this work is the core zone, which represents dense scar. Additional surfaces, when present, include the border zone, healthy myocardium, total scar, and corridor-related structures such as centreline and filtered centreline meshes. While this project focuses on CZ geometry, the presence of BZ and corridor surfaces may support future multi-class geometric analyses.

### Clinical Metadata

A patient-level clinical registry CSV accompanies the imaging-derived data. The source file was exported from Excel and uses Spanish column headers, comma decimal separators, and dates in `DD/MM/YYYY` format. The loader in `src/clinical_data.py` standardises this file by detecting the encoding and delimiter, removing empty or unnamed separator columns, renaming columns to Python-friendly `snake_case`, resolving duplicated headers such as the repeated `FECHA IAM`, and casting dates, continuous variables, binary indicators, and coded categorical fields to appropriate data types.

The CSV contains one row per patient and includes the following groups of variables:

| Group | Examples |
|---|---|
| Identification and inclusion | `inclusion`, `patient_id`, `notes_misc` |
| MRI session | `date_mri`, `mri_type`, `diagnosis` |
| MRI-derived scar geometry | `lv_mass_g`, `bz_core_g`, `bz_core_pct`, `bz_g`, `bz_pct`, `core_g`, `core_pct`, `channels`, `channel_mass_g` |
| Demographics | `date_birth`, `age`, `sex` |
| Cardiovascular risk factors | `hta`, `dlp`, `dm`, `smoking`, `afib`, `nyha`, `mi_yn` |
| Infarction and revascularisation history | `date_mi`, `date_revasc`, `revasc_type`, `revasc_notes` |
| Echocardiography | `lvef_pct`, `lvedd_mm`, `lvesd_mm`, `lvedv_ml`, `lvesv_ml`, `la_mm`, `septum_mm`, `post_wall_mm` |
| Devices | `device_carrier`, `date_device_implant`, `device_type`, `icd_prevention_type`, `biomonitor` |
| Ventricular arrhythmias | `clinical_va`, `inducible_va`, `pvc_burden`, `icd_therapies`, `va_ablation` |
| Six-month follow-up | `date_6m_followup`, `admission_hf`, `admission_angina`, `followup_va_yn`, `death_yn` |
| Medication | `med_acei`, `med_arb`, `med_bb`, `med_spiro`, `med_furosemide`, `med_digoxin`, `med_amiodarone` |

Several clinical fields are incomplete. In particular, missingness is concentrated in follow-up dates and event-detail variables, such as dates of ventricular arrhythmia follow-up, ICD therapy, death, ablation, device implantation, and electrophysiological study. This is expected in a registry where many event-specific fields are only filled when the event occurs. The clinical registry is therefore used in this chapter mainly for cohort characterisation. It is not used as a modelling input in the current geometric pipeline, although it may serve as contextual or conditioning information in future extensions.

### Visual Quality Control

To verify segmentation quality and assess scar morphology across the cohort, an automated three-dimensional rendering pipeline was implemented using PyVista. For each case in the known set, the LV shell is displayed as a semi-transparent grey surface and the core zone is overlaid in opaque red from an isometric viewpoint. The resulting thumbnails are assembled into a per-patient grid, enabling rapid visual inspection of segmentation consistency, scar extent, and anatomical plausibility. This quality control step confirmed that the 84 known cases produce geometrically coherent scar representations suitable for downstream descriptor extraction.

In summary, the analysis cohort comprises 105 post-myocardial infarction patients with three-dimensional CMR-derived LV scar segmentations, of which 84 are used for exploratory analysis and model development. The data span five acquisition years (2021-2025), include a predominantly male population consistent with the epidemiology of myocardial infarction, and provide multi-class tissue surfaces at a resolution suitable for population-level geometric analysis. The following chapter describes the methods used to extract compact geometric descriptors from these three-dimensional scar representations.

## Annex A. Excluded Patient Cases

The following table lists the 63 patient folders excluded from the geometric analysis because they lacked either a valid LV reconstruction or a usable core zone surface in the working ADAS 3D export.

| Patient ID | Status | Reason |
|---|---|---|
| 00021123 | no LV | no DE-MRI folder found |
| 01041694 | no LV | no DE-MRI folder found |
| 02002319 | no LV | no DE-MRI folder found |
| 02011130 | no LV | no DE-MRI folder found |
| 02020684 | no LV | no DE-MRI folder found |
| 03015023 | no LV | no DE-MRI folder found |
| 03031339 | no LV | no DE-MRI folder found |
| 04020449 | no LV | no DE-MRI folder found |
| 05001858 | no LV | no DE-MRI folder found |
| 05027445 | no LV | no DE-MRI folder found |
| 05044223 | no LV | no DE-MRI folder found |
| 06040109 | no LV | no DE-MRI folder found |
| 06050540 | no LV | no DE-MRI folder found |
| 07037693 | no LV | no DE-MRI folder found |
| 07043391 | no LV | no DE-MRI folder found |
| 07044737 | no LV | no DE-MRI folder found |
| 08022556 | no LV | no DE-MRI folder found |
| 18022861 | no LV | no DE-MRI folder found |
| 18042739 | no LV | no DE-MRI folder found |
| 98027481 | no LV | no DE-MRI folder found |
| 98037676 | no LV | no DE-MRI folder found |
| caso11hsp | no LV | no DE-MRI folder found |
| CASO12HSP | no LV | no DE-MRI folder found |
| caso14hsp | no LV | no DE-MRI folder found |
| caso8hsp | no LV | no DE-MRI folder found |
| caso9hsp | no LV | no DE-MRI folder found |
| 00016504 | no LV | no DE-MRI folder found |
| 00050633 | no LV | no DE-MRI folder found |
| 01004887 | no LV | no DE-MRI folder found |
| 12805382 | no LV | no DE-MRI folder found |
| 19033528 | no LV | no DE-MRI folder found |
| 2025015442 | no LV | no DE-MRI folder found |
| 50105639 | no LV | no DE-MRI folder found |
| 50149701 | no LV | no DE-MRI folder found |
| 50187501 | no LV | no DE-MRI folder found |
| 00044761 | no LV | no `Basal/Data` path |
| 02000042 | no LV | no `Basal/Data` path |
| 02045555 | no LV | no `Basal/Data` path |
| 09019007 | no LV | no `Basal/Data` path |
| 12001780 | no LV | no `Basal/Data` path |
| 13018580 | no LV | no `Basal/Data` path |
| 14008873 | no LV | no `Basal/Data` path |
| 21058337 | no LV | no `Basal/Data` path |
| 972100446 | no LV | no `Basal/Data` path |
| 987012904 | no LV | no `Basal/Data` path |
| 22004656 | no LV | no `Basal/Data` path |
| 20027330 | no LV | no `Basal/Data` path |
| MLF00254 | no LV | no `Basal/Data` path |
| 1111593944 | no LV | no `Basal/Data` path |
| 2491608 | no LV | no `Basal/Data` path |
| 07039239 | no CZ | `TISSUE` exists but contains no VTK files |
| 10025183 | no CZ | `TISSUE` exists but contains no VTK files |
| 13017471 | no CZ | `TISSUE` exists but contains no VTK files |
| 20023806 | no CZ | `TISSUE` exists but contains no VTK files |
| 50033820 | no CZ | `TISSUE` exists but contains no VTK files |
| 50037142 | no CZ | `TISSUE` exists but contains no VTK files |
| 50042306 | no CZ | `TISSUE` exists but contains no VTK files |
| 50045793 | no CZ | `TISSUE` exists but contains no VTK files |
| 50052223 | no CZ | `TISSUE` exists but contains no VTK files |
| 50054988 | no CZ | `TISSUE` exists but contains no VTK files |
| 50056378 | no CZ | `TISSUE` exists but contains no VTK files |
| 95025044 | no CZ | `TISSUE` exists but contains no VTK files |
| 96024662 | no CZ | corridor meshes present, but no CZ surface match |

## References for Cohort Context

1. Albrektsen G, Heuch I, Lochen M-L, Thelle DS, Wilsgaard T, Njolstad I, Bonaa KH. Lifelong gender gap in risk of incident myocardial infarction: the Tromso Study. *JAMA Internal Medicine*. 2016.
2. StatPearls. Myocardial Infarction. NCBI Bookshelf.
