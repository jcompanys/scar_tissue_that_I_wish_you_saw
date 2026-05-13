# Memory Folder Guide

Root entry point: `MEMORY.md`

## Folder Layout

| Folder | Purpose |
|---|---|
| `overview/` | Project-wide context, repo layout, conventions, and high-level status. |
| `notebooks/` | Notebook-specific memory pages, named after the notebook they describe. |
| `pipelines/` | Method or analysis pipeline notes that may span several notebooks. |
| `reference/` | Catalogs, API/module references, and cross-file maps. |

## Current Files

| Path | Contents |
|---|---|
| `overview/project_lv_scar.md` | Main LV scar segmentation project context. |
| `notebooks/M01_scar_exploration.md` | M01 scar exploration notebook notes. |
| `notebooks/M03_polar_diagnostics.md` | M03 polar diagnostics and projection comparison notes. |
| `pipelines/isomap_pipeline.md` | ISOMAP scar-fragment shape analysis pipeline. |
| `reference/M0_notebooks_and_src_py.md` | Catalog for all M0 notebooks and `src/*.py` modules. |

## Naming Convention

- Notebook pages: `notebooks/<notebook_stem>.md`
- Pipeline pages: `pipelines/<short_pipeline_name>.md`
- Project-wide notes: `overview/<topic>.md`
- Cross-file catalogs: `reference/<topic>.md`
