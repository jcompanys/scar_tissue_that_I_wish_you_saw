# Code cleanup plan

Branch: `refactor/code-cleanup` (from tag `thesis-final-2026-10`). Never merged back.
Backups: `backup/thesis-v1` (917f88a), `backup/thesis-final-2026-10` + tag (f6248a6).

## Ground rules (decided)

- Results were already presented. Outputs must stay identical to the reference run.
- Reference outputs: `lv-scar-segmentation/results-reference/`, a copy of the results
  produced by the final thesis code (`backup/thesis-final-2026-10`, tag f6248a6).
- Public names (functions, classes, constants, modules) are kept. Notebooks import them.
  Local variables inside a function may get clearer names.
- Claude does not run notebooks. The user runs them; Claude compares the outputs.
- Bugs or doubtful logic are **reported, not fixed**.

## Change levels

| Level | Allowed | Examples | Check |
|---|---|---|---|
| **A. Docs only** | yes | docstrings, comments, signature type hints | AST check, no run needed |
| **B. Refactor, same outputs** | yes | remove unused imports, split long functions, move duplicated notebook helpers into `src`, merge duplicate imports, rename locals, fix section numbering, delete dead code/cells, clear notebook outputs | user re-runs, `tools/compare_results.py` against reference |
| **C. Changes results** | **no** | bug fixes, new defaults, changed thresholds | reported only |

Level A and level B changes go in separate commits, so an A commit can be verified by the AST
check alone.

## Conventions

1. Docstrings: NumPy style (Parameters / Returns / Notes). Geometry conventions go in Notes.
2. Language: English.
3. Type hints on public `src/` function signatures. Never added inside class bodies
   (dataclass fields are created by annotations).
4. Comments explain *why* (geometry choices, thresholds, orientation flips), not *what*.
5. Module header in every `src/` file: purpose, main API, which notebooks use it.
6. Notebook outputs may be cleared (level B); the backups keep the original outputs.
7. One commit per module or notebook and per level.

## Reference results

`lv-scar-segmentation/results-reference/` (540 files, gitignored) is a copy of the
results from the final thesis code, the same code this branch starts from. Any
difference after a cleanup run is therefore caused by the cleanup.

`results/` is gitignored, so every run on any branch overwrites it. Never overwrite
`results-reference/`. To rebuild it: `git switch backup/thesis-final-2026-10`, run all
notebooks, copy `results/` again, switch back.

## Per-file workflow

1. **Map**: read the file, list functions, find every caller (other modules + notebooks).
2. **Review report** (no edits yet), per function:
   - what it does, inputs/outputs, units, coordinate frame
   - unclear parts, magic numbers, dead code, duplicates elsewhere
   - suspected bugs (flag only)
   - proposed level B changes
3. **User decides** which level B changes to apply.
4. **Level A commit**: header, docstrings, comments, type hints.
   Check: imports cleanly + `tools/check_ast_equiv.py` shows no code difference.
5. **Level B commit(s)**: approved refactors.
   Check: imports cleanly; user re-runs the dependent notebooks;
   `tools/compare_results.py` against the reference.
6. **User reviews diff** before each commit.

## Phase 1: `src/` (dependency order)

| # | Module | Used by |
|---|---|---|
| 1 | `data_loading.py` | all notebooks except M06 |
| 2 | `mesh_utils.py` | F01, M01–M05 |
| 3 | `lv_geometry.py` | F01, M01, M03, M04, M05 (`_hull_apex_candidate` imported by F01, M03) |
| 4 | `clinical_data.py` | M01–M04 |
| 5 | `scar_characterization.py` | F01, M04 |
| 6 | `scar_analysis.py` | M04 |
| 7 | `plot_style.py` | most |
| 8 | `scar_statistics.py`, `cone_bspline_simple.py` | **nothing**: candidates for deletion (level B) |

## Phase 2: notebooks

Order: M01 → M02_clinical_eda → M02_isomaps → M03 → M04 → M04_BZ_SCC → M05 → M06 → F01.

Per notebook:
- intro markdown cell: purpose, inputs, outputs, which later notebooks consume its outputs
- fix section numbering (M04: 2.3 → 2.4 → 2.2; M05: 4c/4d under section 5)
- markdown cells and `#` comments explaining each step
- merge duplicate imports (`scan as _f7_scan`, `long_axis as _la_robust2`, ...)
- replace local copies of `src` helpers (M02_isomaps, M06) when they are identical
- remove dead or debug cells
- clear outputs

## Phase 3: verification

After each level B batch: user re-runs the listed notebooks, then
`python tools/compare_results.py lv-scar-segmentation/results-reference lv-scar-segmentation/results`.

## Phase 4: docs

Refresh `README.md` (tree, notebooks, modules) and `auxiliary/memory`.

## Progress

| Item | Status |
|---|---|
| Phase 0 conventions | agreed |
| Reference results | done (`results-reference/`, matches `results/`) |
| Check scripts (`tools/`) | written |
| `data_loading.py` | next |
