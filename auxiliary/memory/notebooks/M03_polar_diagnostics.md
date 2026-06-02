# M03 — Polar Diagnostics and Projection Comparison

**Notebook file:** `M03_polar_diagnostics.ipynb`  
**Output directory:** `results/06_polar_diagnostics/`

---

## Overview

This notebook validates the anatomical reference system used to build LV scar polar maps, then benchmarks three complementary projection methods (vertex, area, surface-normalized coverage) against each other at both the per-patient and population level.

The same coordinate convention is enforced throughout: the long axis runs apex → base, the RV/septal reference defines `theta = 0`, and every 3D display is rolled by 45° so the blue septal vector points toward 10:30 in the apex-to-base view.

| Section | Contents |
|---|---|
| 1 | Imports, data loading, clinical sex labels, and shared helpers |
| 2 | D1–D6 septal-reference diagnostics, raw LV PCA axis vs robust axis, and candidate-axis inspection |
| 3 | Polar-map computational steps, all-case trial maps, and population superposition |
| 4 | Polar-map quality checks |
| 5 | Vertex vs area vs surface-normalized comparison, including AHA and ring-level summaries |

---

## Section 1 — Setup, Data Loading, and Shared Helpers

### 1.1 Configuration constants

| Constant | Value | Purpose |
|---|---|---|
| `N_R` | `64` | Radial bins in the polar grid |
| `N_THETA` | `128` | Angular bins in the polar grid |
| `DIP_Z_LO` / `DIP_Z_HI` | `0.60` / `0.85` | Basal-ring z-fraction window for Method A |
| `DIP_N_BINS` | `180` | Angular bins for radial-envelope estimation |
| `DIP_SMOOTH_BINS` | `7` | Smoothing window for dip detection |
| `AHA_DISPLAY_DEG` | `45.0` | Extra display roll so septal vector → 10:30 |
| `D1_N` / `D4_N` | `16` / `12` | Subset sizes for D1 and D4 visuals |
| `SEED` | `42` | RNG seed for fixed case shuffle |
| `HD_ROOT` | `F:/RM_TEKNON_DEVELOP` | Root path for raw dataset |
| `SPLIT_JSON` | `results/split.json` | Pre-computed 80/20 train/held-out split |
| `OUT` | `results/06_polar_diagnostics/` | All outputs are written here |

### 1.2 Data loading

`scan(HD_ROOT)` discovers all patient cases. The known split (84 cases) is loaded from `split.json`. Cases are filtered to those that have both a Core Zone (CZ) mesh and an LV shell. A fixed shuffle (`SEED=42`) ensures stable per-case subsets across runs.

**Runtime result:** 84 known CZ cases, 84/84 with RV mesh (100%).

Sex labels are loaded from the clinical CSV (`develop-vt.csv`) via `src.clinical_data.load`. The `SEXO` field is mapped to `'M'` (1) or `'F'` (2) and stored in `sex_map` keyed by `patient_id`. **Runtime result:** 94 patients with labels, M=79, F=15 (only those that overlap with the known split are used downstream).

### 1.3 Shared helpers

#### Septal reference methods

Three methods are implemented for finding the RV/septal reference direction. All return a unit vector in the plane perpendicular to the long axis:

| Method | Function | Description |
|---|---|---|
| **A — Dip** | `_rv_reference_dip()` | Extracts the basal ring (60–85% up from apex), builds a radial envelope histogram, and finds the midpoint of the two deepest inward dips. These dips correspond to the RV attachment sites. |
| **B — PCA** | `rv_reference_pca()` | PCA of the basal ring; the half with smaller mean radius is assumed to be the septal/concave side. Fallback when no RV mesh is available. |
| **C — RV mesh** | `rv_reference_from_rv_mesh()` | Computes the vector from LV centroid to RV centroid in the basal cross-section plane. Ground truth when an RV mesh is present. |

`oriented_septal_reference()` canonicalises any of these to a consistent hemisphere.

#### Polar map construction: `_polar_map_for_case()`

The main per-case projection function:

1. Load `shell` (LV surface), `cz` (Core Zone mesh), and optionally `rv` (RV mesh).
2. Compute polar geometry via `prepare_polar_geometry(shell, ref_fn, rv_mesh)` — sets up the `(s, theta)` coordinate frame.
3. Project CZ vertices with `project_points_to_polar(cz.points, geom)`.
4. **Flip for clockwise anatomical orientation:** `theta = (2π − theta) % 2π`. This bakes in `set_theta_zero_location('W') + set_theta_direction(-1)`, giving the standard AHA apical bull's-eye view (SEPT=West, ANT=top, LAT=right, INF=bottom).
5. Bin `(s, theta)` into a `(N_R × N_THETA)` boolean mask via `_polar_mask()`.

#### Polar coordinate system

| Coordinate | Definition | Range |
|---|---|---|
| `s` | Normalised distance along long axis | 0 = apex (centre), 1 = base (outer ring) |
| `theta` | Circumferential angle from septal reference (CW) | 0 = RV/septum (West), π/2 = ANT (top), π = LAT (right), 3π/2 = INF (bottom) |

#### Display helper: `_draw_bullseye_grid()`

Draws concentric rings at `r = 0.25, 0.50, 0.75, 1.0` and radial spokes matching AHA segment boundaries. Always sets `theta_zero_location='W'` and `theta_direction=-1`. Marks the septal reference with a blue dot at `theta=0`.

---

## Section 2 — Septal Reference and Axis Diagnostics

This section validates axis orientation and septal reference quality before any polar maps are computed.

### D1 — Basal ring reference arrows

**What it shows:** For each of 16 cases, a 2D cross-section of the basal ring (60–85% up from apex) viewed from the apex toward the base. The LV wall points are coloured by distance from centroid (yellow = far, dark = close). Three reference arrows are overlaid:

- **Red arrow** — Method A (dip midpoint)
- **Blue arrow** — Method B (PCA concave half)
- **Green arrow** — Method C (RV mesh centroid, ground truth)
- **Red dots** — the two selected dip peaks used by Method A

**What to look for:** All three arrows should point toward the same concave/flat side of the kidney-shaped ring. Consistent agreement across cases indicates a reliable reference. Divergence between Method C and Methods A/B flags cases where the fallback methods may fail.

**Output:** `D1_basal_ring_reference.png`

### D2 — Angular disagreement distribution

Computes pairwise circular angular differences between all three methods for every case and plots a histogram. When the RV mesh is available, Method C is the anatomical reference.

- A vs B: PCA vs dip comparison
- A vs C: dip vs RV ground truth
- B vs C: PCA vs RV ground truth

Low disagreement (< 45°) indicates the automatic methods reliably approximate the RV insertion site. A dashed line at 45° provides a visual threshold.

**Output:** `D2_angular_disagreement.png`

### D3 — Septal direction in raw ADAS coordinates

Plots the RV-based (Method C) and PCA (Method B) septal reference vectors in the original un-aligned mesh coordinate system across three projection planes (X-Y, X-Z, Y-Z).

- **Tight cluster → low circular std:** ADAS acquires patients in a consistent world orientation, meaning the raw coordinates partially encode the septal direction.
- **Scattered cloud → high circular std:** Acquisition orientation varies; the raw coordinate cannot be used as a proxy for anatomy.

**Output:** `D3_septal_world_directions.png`

### D4 — Five-row 3D alignment grid

Renders `D4_N = 12` cases in five different states to show the full alignment pipeline:

| Row | Label | Description |
|---|---|---|
| A | Raw | Raw ADAS coordinates; isometric view |
| B | Raw+RV | Side view showing the LV-RV spatial relationship before any alignment |
| C | Long axis | Long-axis aligned only; +Z = apex-to-base, no in-plane rotation |
| D | Anatomical | Long-axis + RV septal reference aligned; RV positioned at −X for display |
| E | Apex view | Final apex-to-base view; this is the polar map viewpoint |

Rows D and E are the final anatomically-aligned views used for polar map construction.

**Output:** `D4_alignment_grid.png`

### D5 — Long axis and septal vector overlay

Two 3D views per case:

| View | Camera | Check |
|---|---|---|
| Side | From +Y | Apex (red dot) at top, green long axis vertical, RV visible |
| Smashed | From apex toward base | Long axis = red dot at centre; blue septal vector → 10:30 |

Green line = `long_axis()` result. Red sphere = apex. Blue arrow = LV centroid → RV centroid vector, rolled by `AHA_DISPLAY_DEG = 45°`.

**Output:** `D5_axis_septal_overlay.png`

### Axis diagnostic — Raw LV PCA vs robust long axis

Side-by-side comparison of the naive SVD PC1 (orange) vs the robust `long_axis()` result (green) in original un-aligned coordinates for 12 cases. Identifies cases where the raw PC1 is flipped or poorly conditioned relative to the anatomical apex-to-base direction.

**Output:** `axis_diagnostic.png`

### D6 — Long axis candidates

Shows all candidate axes considered by `long_axis()` for 12 cases and highlights the final chosen axis.

| Visual element | Meaning |
|---|---|
| Blue line | RV SVD PC1 candidate, when an RV mesh is available |
| Orange line | LV SVD PC1, the naive long-axis candidate |
| Purple line | LV SVD PC2 |
| Thick green line | Final `long_axis()` output, recomputed as apex-to-centroid after candidate selection |
| Red sphere | Apex selected from the winning candidate |

The candidates are evaluated with `_hull_apex_candidate()`. The winning candidate is the one with the highest apex eccentricity, but the final green vector is not simply the winning raw PC direction: it is recomputed as `normalize(centroid - apex)` to match `long_axis()` exactly.

**Output:** `d6_long_axis_candidates.png`

---

## Section 3 — Polar Map Construction

### 3.1 Computational steps

Each of 4 cases is rendered as a 5-panel pipeline strip:

| Panel | Description |
|---|---|
| 1. Side view | Aligned LV + CZ; green axis, red apex, blue septal vector |
| 2. Smashed | Apex-to-base view; septal vector rolled to 10:30 |
| 3. theta coloring | HSV rainbow on the smashed view — smooth continuous rainbow = correct circumferential parameterisation |
| 4. s coloring | Blue = apex (s=0), red = base (s=1) on the side view |
| 5. Bull's-eye | Final polar map: (s, theta) binned to 64 × 128 |

**Output:** `polar_map_steps.png`

### 3.2 All-case trial polar maps

Computes `pmaps` — the dictionary of per-case binary polar maps — for all 84 known cases. Each map has shape `(N_R, N_THETA) = (64, 128)` and contains `1.0` where at least one CZ vertex falls in that bin, `0.0` elsewhere.

Reference method priority: Method C (RV mesh) when available, else Method B (PCA). Title colour in the grid indicates the method: **blue** = RV mesh, **orange** = PCA fallback.

`pmaps` is the central data structure reused by Sections 3.3, 4, and 5.

**Output:** `polar_maps_all_cases.png`

### 3.3 Population polar-map superposition

**Cell id:** `4fd7ff89`

Computes and displays a pixel-wise mean scar fraction across all patients, split by group:

```
_mean_scar_map(items):
    stack = np.stack([v for _, v in items], axis=0)
    return stack.mean(axis=0)   # pixel-wise mean scar fraction
```

Each pixel value = fraction of patients that have scar at that `(r, theta)` location.

**Layout:** Top row = all patients (n=84); bottom row = Male (n=30) and Female (n=7) subgroups.

**Colour scale note:** `_vmax` is computed from the **99th percentile of all values concatenated across all subgroups** (all + male + female). Because the Female subgroup (n=7) is sparse and near-binary, its ravel contributes many near-zero values and shifts the 99th percentile higher than the vertex-only map's auto-scale would produce. This means the superposition plot and the vertex frequency map in Section 5.3 display the same underlying quantity but with different colour scales. To align them, compute `_vmax` from `_all_map` alone.

**Output:** `polar_maps_superposition.png`

---

## Section 4 — Polar-Map Quality Checks

Checks 1–4 validate axis orientation, coordinate colouring, and method agreement on 12 cases (`_QC_SUBSET`). Check 4f adds population-level scar profile plots.

### 4.1 Shared QC rendering helpers

`_qc_prepare(case)` is the shared setup function for all QC checks. It aligns the mesh and returns a dictionary containing:

- `shell_t` — transformed LV shell
- `cz_t` — transformed CZ mesh (if present)
- `s_vals` — per-vertex normalised longitudinal coordinate
- `theta_vals` — per-vertex circumferential angle
- `R_original` — reference orientation transform
- `R_disp` — display roll transform
- `span` — LV long-axis length in mm

Reference selection: Method C (RV mesh) when available, else Method B (PCA).

### 4.2 Check 1 — Axis and anatomical orientation

Two-row panel per case:
- **Row A:** Tilted axis comparison — orange = raw SVD PC1, green = robust `long_axis()`.
- **Row B:** Side view after alignment — apex should be at top, base at bottom, RV visible on the left.

**What to check:** Green axis should be vertical in row B. Apex (red sphere) should be at the superior tip.

**Output:** `polar_map_qc_axis.png`

### 4.3 Check 2 — Longitudinal s-coordinate sanity

LV shell coloured by the `s` coordinate (blue = apex, red = base) in a side view. A smooth gradient from tip (blue) to base (red) indicates correct long-axis orientation and normalisation.

**Output:** `polar_map_qc_s_coloring.png`

### 4.4 Check 3 — Circumferential theta-coordinate sanity

LV shell coloured by `theta` using the HSV colormap in the apex-to-base view. A smooth, single-revolution rainbow wrapping around the ventricle = correct circumferential parameterisation. The colour seam (HSV reset) marks `theta = 0/2π` — the septal reference direction — and should align with the RV insertion region.

**What breaks it:** Patchy or multi-seam colouring indicates a problem with the reference direction or a topologically inconsistent mesh.

**Output:** `polar_map_qc_theta_coloring.png`

### 4.5 Check 4 — 3D-to-bull's-eye method comparison

Side-by-side comparison of the 3D apex-to-base view and the resulting bull's-eye for each case, rendered separately using Method B (PCA fallback) and Method C (RV mesh ground truth). The two bull's-eyes should show the same scar location; disagreement identifies cases where the PCA reference is unreliable.

**Output:** `polar_map_qc_method_comparison.png`

### 4.6 Complementary scar profiles (Check 4f)

Five supplementary population-level plots derived from the `pmaps` stack:

An additional ring-level vertex coverage plot is computed before the 4f profile plots using `polar_map_counts(case)`. It bins both CZ vertices and LV shell vertices into the same `(N_R × N_THETA)` grid and reports CZ vertices / LV vertices per AHA ring zone.

**Output:** `aha_scar_density_rings.png`

**Runtime values:**

| Ring zone | CZ / LV vertex fraction | Counts |
|---|---:|---:|
| Apex (17) | 0.6196 | 41,654 / 67,228 |
| Apical (13–16) | 0.7199 | 73,290 / 101,803 |
| Mid (7–12) | 0.2957 | 82,900 / 280,392 |
| Basal (1–6) | 0.0526 | 23,644 / 449,545 |

#### 4f-a — Circumferential scar profile by ring zone

Mean scar density as a function of `theta` (0 = septum), broken down by four longitudinal zones:

| Zone | Row range | Colour |
|---|---|---|
| Apex | rows 0–15 (s = 0–0.25) | Purple |
| Apical | rows 16–31 (s = 0.25–0.50) | Red |
| Mid | rows 32–47 (s = 0.50–0.75) | Green |
| Basal | rows 48–63 (s = 0.75–1.0) | Blue |

**Output:** `polar_circumferential_profile.png`

#### 4f-b — Longitudinal scar profile

Mean scar density as a function of `s` (0 = apex, 1 = base), averaged over all theta. Shows the apex-to-base gradient of scar prevalence in the population.

**Output:** `polar_longitudinal_profile.png`

#### 4f-c — Hotspot maps at prevalence thresholds

Three bull's-eyes showing which grid cells have scar in ≥ 2%, ≥ 5%, and ≥ 10% of patients respectively. Cells below threshold are zero. Segment labels are shown.

**Output:** `polar_hotspot_maps.png`

#### 4f-d — Density map with iso-contours

Population scar density map (`hot_r` colormap) with blue iso-contour lines at 5%, 10%, and 20% prevalence overlaid.

**Output:** `polar_density_contours.png`

#### 4f-e — Scar eccentricity (angular entropy per ring row)

For each radial row `r_i`, normalised Shannon entropy is computed over the theta distribution of scar density:

```
H(r_i) = -Σ p_j * log(p_j) / log(N_THETA)
```

`H = 0` means all scar in that ring is concentrated in one angular location (focal). `H = 1` means scar is uniformly distributed around the full circumference (diffuse).

Printed summary: mean eccentricity per zone (Apex / Apical / Mid / Basal).

**Output:** `polar_eccentricity.png`

---

## Section 5 — Vertex, Area, and Surface-Normalized Comparison

Three complementary ways to represent the same scar on a polar grid:

| Method | Quantity | Sensitivity |
|---|---|---|
| **Vertex** | Binary CZ vertex presence per bin | Mesh sampling density; biased by local vertex density |
| **Area** | Raw CZ cell surface area (mm²) per bin | Absolute scar burden; biased toward basal bins (larger physical area) |
| **Coverage** | CZ area / LV shell area per bin (ratio) | Normalised: removes both bin-size and inter-patient LV-size bias |

All three use the same coordinate system and reference method.

### 5.1 Area and coverage helper functions

#### `_surface_cell_centers_and_areas(mesh)`

Extracts the surface of a volumetric mesh via `extract_surface`, computes cell areas using `pyvista.compute_cell_sizes`, and returns `(centers, areas)` for all cells with `area > 0`.

#### `_accumulate_surface_area(mesh, geom)`

Projects cell centers to `(theta, s)` using `project_points_to_polar`, applies the same CW flip as the vertex method, then accumulates areas into the `(N_R × N_THETA)` grid using `np.add.at`. Each cell contributes its entire area to the single bin containing its center.

#### `area_surface_polar_map_06(case)`

Top-level per-case wrapper. Returns a dict with:

| Key | Type | Description |
|---|---|---|
| `cz_area` | `(64, 128) float` | CZ cell area per bin (mm²) |
| `lv_area` | `(64, 128) float` | LV shell cell area per bin (mm²) |
| `coverage` | `(64, 128) float` | `cz_area / lv_area`; `NaN` where `lv_area ≈ 0` |
| `ref_method` | `str` | `'rv'` or `'pca'` |
| `cz_total` | `float` | Total CZ surface area (mm²) |
| `lv_total` | `float` | Total LV shell area (mm²) |

Reference selection mirrors the vertex method: Method C when an RV mesh exists, else Method B.

**Runtime result:** 84/84 area maps computed, 0 failures. All used `ref_method='rv'`.

#### `draw_bullseye_cmp(arr, ax, ...)`

Plotting helper for Section 5 comparison panels. Applies the reference orientation (`theta_zero='W'`, clockwise), draws the grid, adds anatomical labels (SEPT/ANT/LAT/INF), optionally overlays AHA segment numbers, and auto-scales `vmax` to the 99th percentile of finite values when not specified.

### 5.2 Compute area and coverage polar maps

Iterates over all 84 known cases and stores results in `area_results_06` dict keyed by `year/patient_id`. Failed cases are collected in `area_failures_06`.

### 5.3 Population-level projection comparison

**Cell id:** `892324d0`

Computes `_shared` — the intersection of keys present in both `pmaps` and `area_results_06` (84 cases). Builds four population-averaged grids:

```python
_vert_freq = _vert_stack.mean(axis=0)               # fraction with any vertex
_area_freq = (_cz_stack > EPS).astype(float).mean() # fraction with any area
_mean_area = _cz_stack.mean(axis=0)                 # mean CZ area (mm²)
_mean_cov  = np.nanmean(_cov_stack, axis=0)          # mean per-patient coverage
```

**Key interpretation:**

- `_vert_freq` and `_area_freq` should agree on **location** — they measure the same thing differently. Any divergence points to mesh sampling artefacts.
- `_mean_area` is systematically higher in basal bins because basal bins cover more physical area.
- `_mean_cov` corrects for available LV surface per bin and is the most comparable metric across patients and LV sizes.

**Colour scale:** Each of the four panels uses independent auto-scaling (99th percentile of that panel's own values). This is why the "Vertex" panel in this figure and the "All patients" panel in Section 3.3 look different despite computing the same quantity — see the known issue below.

**Output:** `comparison_population_maps.png`

### 5.4 Per-patient projection comparison

Renders a 3-row grid for up to 12 cases (`_CMP_N`): vertex (Reds), area (magma), coverage (viridis). Useful for identifying individual outliers where the three methods disagree strongly — for example, a case where vertex shows focal scar but coverage shows it is actually a small dense patch on a large basal ring.

**Output:** `comparison_per_patient_grid.png`

### 5.5 AHA-17 segment summary

Aggregates the three metrics per AHA segment, producing a bar chart with three panels:

1. **Vertex scar frequency** — fraction of patients with at least one CZ vertex in the segment
2. **Area scar frequency** — fraction of patients with any CZ area in the segment
3. **Pooled CZ/LV coverage** — total CZ area across all patients divided by total LV area in that segment

The AHA segment grid is computed by `_aha_segment_grid_06()`, which assigns each `(r_i, theta_j)` bin to one of the 17 segments based on radial thresholds (0.25, 0.50, 0.75) and angular boundaries.

**Output:** `comparison_aha_segments.png`

### 5.6 CZ area by LV ring

Pools the Section 5 area results by four radial AHA-style ring zones:

| Ring zone | Radial range |
|---|---|
| Apex (17) | `r < 0.25` |
| Apical (13–16) | `0.25 <= r < 0.50` |
| Mid (7–12) | `0.50 <= r < 0.75` |
| Basal (1–6) | `r >= 0.75` |

For each ring, the metric is pooled CZ area divided by pooled LV area across all cases in `_s5_keys`. This is the area-based analogue of the vertex ring coverage plot in Section 4.6 and is less sensitive to mesh sampling density.

**Output:** `s5_ring_area_coverage.png`

### 5.7 CZ area by LV ring — sex split

Repeats the pooled ring area coverage calculation separately for cases whose `patient_id` has a clinical sex label in `sex_map`.

The notebook constructs:

```python
_keys_m = [k for k in _s5_keys if sex_map.get(k.split('/')[1]) == 'M']
_keys_f = [k for k in _s5_keys if sex_map.get(k.split('/')[1]) == 'F']
```

It then plots side-by-side male and female bars for each ring. The metric remains pooled CZ area / pooled LV area per ring, not mean of per-patient ratios.

**Output:** `s5_ring_area_sex_split.png`

---

## Output Files Summary

| File | Section | Description |
|---|---|---|
| `D1_basal_ring_reference.png` | 2 / D1 | Basal ring cross-sections with three reference arrows |
| `D2_angular_disagreement.png` | 2 / D2 | Histogram of pairwise method angular differences |
| `D3_septal_world_directions.png` | 2 / D3 | Septal direction scatter in raw ADAS coordinates |
| `D4_alignment_grid.png` | 2 / D4 | Five-row 3D alignment grid (raw → anatomical → apex view) |
| `D5_axis_septal_overlay.png` | 2 / D5 | Long-axis and septal vector side + smashed views |
| `axis_diagnostic.png` | 2 | Raw SVD PC1 vs robust long axis |
| `d6_long_axis_candidates.png` | 2 / D6 | Candidate axes tested by `long_axis()` plus final chosen axis |
| `polar_map_steps.png` | 3.1 | Pipeline strip: side → smash → theta/s coloring → bull's-eye |
| `polar_maps_all_cases.png` | 3.2 | Per-case bull's-eye grid, 84 cases |
| `polar_maps_superposition.png` | 3.3 | Population mean scar fraction: all / male / female |
| `polar_map_qc_axis.png` | 4.2 | Axis and anatomical orientation check |
| `polar_map_qc_s_coloring.png` | 4.3 | Longitudinal s-coordinate sanity check |
| `polar_map_qc_theta_coloring.png` | 4.4 | Circumferential theta sanity check |
| `polar_map_qc_method_comparison.png` | 4.5 | 3D-to-bull's-eye method B vs C comparison |
| `polar_circumferential_profile.png` | 4.6a | Mean scar density by theta, broken down by ring zone |
| `polar_longitudinal_profile.png` | 4.6b | Mean scar density by s (apex to base) |
| `polar_hotspot_maps.png` | 4.6c | Hotspot maps at 2%, 5%, 10% thresholds |
| `polar_density_contours.png` | 4.6d | Density map with 5%/10%/20% iso-contours |
| `polar_eccentricity.png` | 4.6e | Angular entropy per ring row |
| `aha_scar_density_rings.png` | 4.6 | Vertex-count CZ/LV coverage by AHA ring zone |
| `comparison_population_maps.png` | 5.3 | Four-panel: vertex / area freq / mean area / mean coverage |
| `comparison_per_patient_grid.png` | 5.4 | Per-patient 3-row grid for 12 cases |
| `comparison_aha_segments.png` | 5.5 | AHA-17 bar chart: vertex / area freq / pooled coverage |
| `s5_ring_area_coverage.png` | 5.6 | Pooled CZ area / LV area by ring zone |
| `s5_ring_area_sex_split.png` | 5.7 | Pooled ring area coverage split by male/female labels |

---

## Known Issues and Caveats

### Colour scale mismatch between Section 3.3 and Section 5.3

The "All patients" superposition (Section 3.3) and the "Vertex scar frequency" panel (Section 5.3) compute the **same quantity** (pixel-wise mean of binary `pmaps` across 84 cases), but appear visually different due to different `vmax` computation:

| Plot | vmax source |
|---|---|
| Section 3.3 superposition | 99th percentile of **all subgroups concatenated** (all + male + female). The sparse female group (n=7) shifts this value upward. |
| Section 5.3 vertex panel | 99th percentile of **the vertex map alone** via `draw_bullseye_cmp` auto-scale. |

**Fix:** In Cell `4fd7ff89` (Section 3.3), compute `_vmax` from `_all_map` only:

```python
# Before (current — includes sex subgroups):
_all_vals = np.concatenate([_all_map.ravel()] + [g[2].ravel() for g in _sex_groups])
_vmax = float(np.nanpercentile(_all_vals[np.isfinite(_all_vals)], 99))

# After (consistent with Section 5.3):
_vmax = float(np.nanpercentile(_all_map.ravel(), 99))
```

### Sex-stratified maps are low-n

The female subgroup contains only 7 patients in the known split. The female superposition map is essentially near-binary (most pixels are either 0 or 1 because few patients contribute), which makes it visually distinct from the male and all-patient maps. It should be interpreted cautiously.

### Vertex method sensitivity to mesh density

The vertex map counts CZ vertices, not surface area. Coarsely meshed CZs may miss thin scar patches, while fine meshes may double-count the same physical region. The coverage map (Section 5) is more robust to this effect.

### `draw_bullseye_cmp` auto-scale

When `vmax` is not passed explicitly, `draw_bullseye_cmp` auto-scales each panel independently to its own 99th percentile. This makes individual panels visually comparable but prevents cross-panel comparison of absolute magnitudes. Always pass an explicit `vmax` when comparing multiple panels on the same scale.

---

## Dependencies

### External libraries

| Library | Usage |
|---|---|
| `numpy` | Array operations throughout |
| `matplotlib` | All plotting |
| `pyvista` | 3D mesh loading, rendering, and cell area computation |
| `scipy.signal.find_peaks` | Radial dip detection in Method A |
| `scipy.ndimage.uniform_filter1d` | Envelope smoothing for Method A |
| `tqdm` | Progress bars |

### Project modules (`src/`)

| Module | Functions used |
|---|---|
| `src.data_loading` | `scan()` |
| `src.mesh_utils` | `best_shell()`, `get_cz_mesh()`, `best_rv_mesh()` |
| `src.lv_geometry` | `long_axis()`, `plane_basis()`, `rv_reference_pca()`, `rv_reference_from_rv_mesh()`, `oriented_septal_reference()`, `compute_transform()`, `apply_transform()`, `rotation_about_z()`, `prepare_polar_geometry()`, `project_points_to_polar()` |
| `src.clinical_data` | `load()` — reads the clinical CSV and normalises the `SEXO` field |
