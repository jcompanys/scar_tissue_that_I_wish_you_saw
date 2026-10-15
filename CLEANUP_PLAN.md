# Code cleanup plan

Branch: `refactor/code-cleanup` (from tag `thesis-final-2026-10`). Never merged back.
Backups: `backup/thesis-v1` (917f88a), `backup/thesis-final-2026-10` + tag (f6248a6).

> **PENDING (2026-10-12): update `results-reference/` for F01, M05, M06** with the Method C
> outputs. User re-ran them 2026-10-12 and the differences were checked: only polar-map
> outputs changed (M05 mean polar Dice 0.095 -> 0.172 weighted_euclidean; scar error and
> segment Dice unchanged). The user decides when to update the reference; Claude reminds
> them (at the latest at the end of Phase 1). Until then, `compare_results.py` on these
> three notebooks is expected to differ. Remove this note once done.

## Level C changes (decided by the user)

| Date | Change | Affects |
|---|---|---|
| 2026-10-09 | `lv_geometry`: default `ref_fn` of `prepare_polar_geometry` and `compute_transform` is Method C (LV -> RV centroid, from `rv_mesh`) instead of Method B (`rv_reference_pca`). Raises if `rv_mesh` is missing (every case has an RV, per the user). F01 `_build_pmap` and M05 (2 cells) now also pass `rv_mesh=rv`, which adds RV PC1 as a long-axis candidate. | F01 `_build_pmap`, M05, M06 outputs (each call site keeps the thesis Method B call as a comment, marked THESIS). M03/M04/F01 no-RV branches now raise instead of using B. M01 `_compute_transform` uses the default but is never called. |

## Ground rules (decided)

- Results were already presented. Outputs must stay identical to the reference run.
- Reference outputs: `lv-scar-segmentation/results-reference/`, a copy of the results
  produced by the final thesis code (`backup/thesis-final-2026-10`, tag f6248a6).
- Public names (functions, classes, constants, modules) are kept. Notebooks import them.
  Local variables inside a function may get clearer names. Exception, decided by the user:
  `PatientCase.anatomy` -> `PatientCase.lv_anatomy` (pairs with `rv_anatomy`).
  `mesh_utils.SHELL_PRIORITY` -> `LV_PRIORITY` (pairs with `RV_PRIORITY`; no notebook imports it).
- `src/` is the main implementation. Notebooks call the `src` function and adapt it
  through its parameters (e.g. a different threshold) instead of keeping their own copy.
  Notebook copies of `src` helpers are replaced in Phase 2 (see Phase 2 notes).
- Claude does not run notebooks. The user runs them; Claude compares the outputs.
- Claude asks the user before touching any notebook, even a one-line change, and says
  which notebooks and cells. Notebook edits otherwise wait for Phase 2.
- Bugs or doubtful logic are **reported, not fixed**.
- `src` functions used by no notebook and no other `src` function are deleted.
- Code is simplified when it reads easier (one way of writing a thing, not a comment
  explaining two). Accepted only if a differential check on real data shows identical
  outputs; float noise around 1e-15 is accepted when the derived maps are identical
  (decided by the user 2026-10-09, `lv_geometry`).
- One tolerance per module: `lv_geometry.EPS = 1e-9`.
- Wrong docstrings or comments are corrected.

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
| 8 | `scar_statistics.py` | **nothing**: deleted (level B, 2026-10-13) — unused star-import wrapper over `scar_analysis` + `scar_characterization`. `cone_bspline_simple.py` deleted with `mesh_utils` B |

## Phase 2: notebooks

Order: M01 → M02_clinical_eda → M02_isomaps → M03 → M04 → M04_BZ_SCC → M05 → M06 → F01.

Per notebook:
- intro markdown cell: purpose, inputs, outputs, which later notebooks consume its outputs
- fix section numbering (M04: 2.3 → 2.4 → 2.2; M05: 4c/4d under section 5)
- markdown cells and `#` comments explaining each step
- merge duplicate imports (`scan as _f7_scan`, `long_axis as _la_robust2`, ...)
- replace local copies of `src` helpers with the `src` function, adapted through its
  parameters (M04, M06). Identical behaviour: level B. Different behaviour:
  run a differential check first; level C if any case changes
- remove dead or debug cells
- clear outputs

## Phase 3: verification

Running all notebooks is slow, so per batch use fast checks; full run once at the end.

| When | Check | Needs data |
|---|---|---|
| every A commit | `tools/check_ast_equiv.py --ref HEAD FILE` | no |
| every B commit | `tools/check_notebook_imports.py` (names notebooks import still exist) | no |
| (all tools) | run with `D:/joan/CONDA_ENVS/scar_char/python.exe` (notebook kernel env) | |
| B on `data_loading` | `tools/check_scan_equiv.py` (old vs new `scan`, field by field) | F: |
| B on other modules | differential test: old module from tag vs new, same inputs, compare outputs | F: |
| end of Phase 2 (once) | user re-runs all notebooks, then `python tools/compare_results.py lv-scar-segmentation/results-reference lv-scar-segmentation/results` | F: |

## Phase 4: docs

Refresh `README.md` (tree, notebooks, modules) and `auxiliary/memory`.

## Progress

| Item | Status |
|---|---|
| Phase 0 conventions | agreed |
| Reference results | done (`results-reference/`, matches `results/`) |
| Check scripts (`tools/`) | written |
| `data_loading.py` | done: B + A committed, imports check OK, `check_scan_equiv.py` IDENTICAL |
| `mesh_utils.py` | done: B + A committed (`anatomy` -> `lv_anatomy` also in `data_loading.py`, `tools/check_scan_equiv.py`); `check_scan_equiv.py` IDENTICAL (168 cases), `check_notebook_imports.py` 0 missing |
| `lv_geometry.py` | B + A committed: `long_axis_qc` deleted, `EPS`, `rv_mesh_ref_fn` added; old vs new on real data: masks identical (106 cases 36x72, 81 cases 64x128), max float diff 9e-16. Then level C: default septal reference = Method C (see top). **Waiting for the user's F01/M05/M06 run.** |

## Flagged issues (level C, not fixed)

| Module | Issue |
|---|---|
| `data_loading` | `border_zone` fallback `*Border*Zone*.vtk` also matches `3D Corridor 3D Border Zone.vtk` |
| `data_loading` | `scar` fallback `*Scar*.vtk` matches any file with "Scar" in its name |
| `data_loading` | `rv_endo`/`rv_epi` fallbacks `*Endo*Layer*`/`*Epi*Layer*` match LV names; safe only inside RV folder |
| `data_loading` | which glob pattern matched is not recorded; unknown if fallbacks ever fired |
| `lv_geometry` | `_hull_apex_candidate` picks the end farthest from the centroid; can return the base as apex (M03 noted an RV-proximity fix, never applied) |
| `lv_geometry` | Thesis results of F01 `_build_pmap`, M05 and M06 used Method B (then the default `ref_fn`), not C, although the RV was loaded. Fixed by the level C change above; the presented results still reflect B |
| M05 | `_scatter_cz_on_ax` (cell 22) draws the scar points mirrored top/bottom (only the first theta flip); figures drawn with it (cell 29 `deterministic_retrieval_<vector>_<territory>`, cell 33 `_make_case_figure` case figures) show anterior and inferior swapped. Numbers unaffected. Fixed by the Phase 2 orientation change |
| M03 | Diagnostics mirrored the same way: cell 57 `polar_map_counts` (flip, no column reversal), cell 38 `theta_cw` colour scalar on the 3-D mesh |
| `lv_geometry` | `compute_transform` keeps `R = I` when the axis is ~ -Z instead of a 180 deg turn |
| `lv_geometry` | `rv_reference_pca`/`rv_reference_from_rv_mesh` measure the basal band from the apex, `project_points_to_polar` from `s_min` |

| M04 | `_connected_fragments` merges all pieces into one when every piece is < 50 points (4 cases); `src.mesh_utils` would keep the largest piece |
| `mesh_utils`, M04 | **Unverified.** `connectivity` runs on the mesh as read; `.clean()` (merges duplicate points) only runs per piece afterwards. VTK links cells by point index, not position, so duplicate points along a seam would split one blob into 2 fragments. Unknown whether ADAS3D exports have duplicates. Fix (`.clean()` before `connectivity`) would change results |

## Research ideas (future work, not part of the cleanup)

| Area | Idea |
|---|---|
| `lv_geometry.polar_mask` | A bin is True if **one** vertex falls in it. Study a minimum vertex count per bin (e.g. >= k vertices, or a density relative to bin area) so isolated vertices do not mark a bin as scar. Effect on per-segment burden and fragment tables to be compared against the current 1-vertex rule. |

## Phase 2 notes (found during Phase 1)

`notebooks/old/` and `notebooks/exploratory/` are archaic: not reviewed, may break.

- ~20 notebook lambdas `(m, ax, ap) -> rv_reference_from_rv_mesh(m, rv, ax, ap)`
  (F01, M03, M04): replace with `lv_geometry.rv_mesh_ref_fn(rv)`. If no notebook ends
  up using it, delete it.
- Septal reference: default is now C. Remove the no-RV `else` branches (F01, M03, M04)
  and the C lambdas (use the default). Delete M01's unused `_compute_transform`.
  `rv_reference_pca` stays only for the M03 A/B/C comparison.
- **Polar orientation: remove the double theta flip + rotate RV centroid to 10:30 (45 deg).**
  Decided by the user 2026-10-15; `src` and notebooks change together in Phase 2 (changing
  `src` alone would break the maps). Level C (rotation changes M04 segment/territory results).
  - Why: `project_points_to_polar` theta drawn as is (`set_theta_zero_location('W')`,
    `set_theta_direction(-1)`) is already the AHA apex-to-base view. Notebooks flip it
    (`2*pi - theta`) and `polar_mask` flips it back (`[:, ::-1]`): maps correct, but any path
    that does only the first flip is mirrored. Checked 2026-10-15 in M04 (RV projected: map has
    the basal RV toward 10-11h = anterior, scatter toward 7-8h).
  - `src`: `project_points_to_polar` returns the display theta with the RV centroid at 45 deg
    (one constant, e.g. `RV_CENTROID_DEG = 45`); `polar_mask` drops `[:, ::-1]` and `roll_deg`;
    `compute_transform` applies the same 45 deg for 3-D views. Blue dot = RV centroid, same constant.
  - Notebooks, remove every `2*pi - theta` flip and the matching column reversals:
    F01 cells 14, 20 (flip + `grid[:, ::-1]`), 22, 24; M03 cells 6, 38, 41, 42, 55, 57, 61
    (flip + `grid[:, ::-1]`), 75, 76; M04 cells 14, 60 (flip + `grid[:, ::-1]`), 85 (flip x3 +
    `ti = n_theta - 1 - ti`); M05 cells 8, 22; M06 cell 5 (3 flips; check the flat-map
    flip test still works).
  - Blue-dot helpers (`_mark_septal_point*` in F01, M03, M04, M05, M06): draw at the RV centroid
    angle from `src` instead of the fixed 60 deg.
  - Then re-run F01, M03, M04, M05, M06 and update `results-reference/`.
- Remove every `AHA_DISPLAY_DEG` (config constants, `roll_deg=AHA_DISPLAY_DEG` arguments,
  the M05 scatter `+ np.deg2rad(AHA_DISPLAY_DEG)`): the RV-centroid angle is to live in `src`
  (planned `lv_geometry` constant, 45 deg = 10:30), so the notebook roll becomes redundant.
- M05 cell ~2182 adds `AHA_DISPLAY_DEG` to theta directly; `polar_mask` rolls columns.
  Same while the value is 0.
- M06 has its own `_load_case` (fixed file names, needs flat VTK); compare with `scan`.
- `notebooks/old/05_inspect_patient_50122488.ipynb` uses fields removed from `PatientCase`.
- M04 `_connected_fragments` (cell ~line 8504, `_MIN_FRAGMENT_VERTS = 50`) duplicates
  `src.mesh_utils.connected_fragments` + `filter_fragments`. Replace with
  `filter_fragments(connected_fragments(cz), min_points=_MIN_FRAGMENT_VERTS)`; the `src`
  default stays 30. Same result when at least one piece has >= 50 points. Differs when:
  - all pieces < 50: M04 keeps the whole mesh (or `[]` if the mesh has < 50 points),
    `src` keeps the largest piece;
  - no `RegionId` after `connectivity` (not expected with current pyvista): M04 returns
    the cleaned mesh if >= 50 points else `[]`, `src` returns the raw mesh.
  Differential check (2026-10-09, 105 cases with CZ): 101 identical, 4 differ
  (2023/94006904, 2024/02027088, 2024/50139122, 2024/50141810: all pieces < 50, M04
  keeps the whole mesh, `src` the largest piece). So the swap is **level C**: not done,
  listed in flagged issues. This is why `connected_fragments` and
  `filter_fragments` stay in `mesh_utils` (their only caller, `cone_bspline_simple.py`,
  was deleted).
