# M02 - ISOMAP Shape Analysis Pipeline

**Notebook file:** `lv-scar-segmentation/notebooks/M02_isomaps.ipynb`  
**Output directory:** `lv-scar-segmentation/results/01_ISOMAPS/`

---

## Overview

ISOMAP is used to analyse the intrinsic geometry of each Core Zone (CZ) scar
fragment. It answers a different question from the polar-map notebooks:

| Method | Question | Coordinate meaning |
|---|---|---|
| ISOMAP | What shape is each scar fragment? | Arbitrary 2D embedding preserving surface/geodesic shape |
| Polar map | Where is scar on the LV? | Anatomical apex-base and circumferential coordinates |

The key design rule is to keep two concepts separate:

1. **Per-fragment shape:** run ISOMAP independently on each connected CZ island.
2. **Whole-scar spatial layout:** store centroid distances between fragments in 3D.

Running ISOMAP on disconnected fragments together is invalid because geodesic
distances between separate islands are infinite. The notebook therefore splits
the CZ mesh first, filters tiny islands, and embeds each retained fragment on its
own.

---

## Current Corrected Notebooks

The corrected notebook series is the `M...` series:

| Notebook | Role |
|---|---|
| `M01_scar_exploration.ipynb` | Dataset scan, EDA, reproducible known/held-out split, and 3D CZ grid |
| `M02_isomaps.ipynb` | ISOMAP scar-fragment shape analysis |
| `M03_polar_diagnostics.ipynb` | Septal-reference diagnostics, polar maps, and vertex/area/coverage comparison |

The older files `01_scar_exploration.ipynb`, `01_ISOMAPS.ipynb`,
`06_polar_diagnostics.ipynb`, and `memory/polar_map_pipeline.md` have been
replaced by the corrected `M...` notebooks and memory notes.

---

## Setup and Inputs

`M02_isomaps.ipynb` runs as a standalone notebook. It:

- finds the project root by walking upward until it sees `src/`, or a nested
  `lv-scar-segmentation/src`
- imports `scan()` from `src.data_loading`
- reads the dataset root from `LV_SCAR_DATA_ROOT`, defaulting to
  `F:/RM_TEKNON_DEVELOP`
- reads the clinical CSV path from `LV_SCAR_CLINICAL_CSV`, defaulting to
  `C:/Users/joan/Desktop/FEINA/UPF/TFG/develop-vt.csv`
- loads `results/split.json` if present, otherwise creates the same seeded
  80/20 split as M01
- writes plots to `results/01_ISOMAPS/`

M02 currently defines `_best_shell()`, `_connected_fragments()`, and
`_filter_fragments()` locally. The same fragment utilities also exist in
`src.mesh_utils` as shared project helpers:

- `best_shell(case)`
- `get_cz_mesh(case)`
- `best_rv_mesh(case)`
- `connected_fragments(mesh)`
- `filter_fragments(fragments, min_points=30)`

If the notebook is refactored later, those `src.mesh_utils` helpers are the
natural source of truth to import instead of duplicating the local functions.

---

## Key Constants

| Constant | Value | Purpose |
|---|---:|---|
| `SEED` | `42` | Reproducible split and ISOMAP subsampling |
| `TRAIN_RATIO` | `0.80` | Known/held-out split ratio if `split.json` is missing |
| `ISOMAP_PTS` | `800` | Maximum sampled vertices per fragment before ISOMAP |
| `ISOMAP_K` | `10` | Neighbours in the ISOMAP k-NN graph |
| `MIN_VERTS_FRAGMENT` | `30` | Minimum vertices to keep a connected fragment |
| `MIN_VERTS_ISOMAP` | `5` | Minimum vertices to attempt a 2D embedding |
| `MAX_PLOT_FRAGMENTS` | `4` | Maximum fragments shown per patient in grids |

---

## Main Result Dataclass

```python
@dataclass
class CZAnalysis:
    key: str                         # "year/patient_id"
    raw_n_fragments: int             # connected islands before filtering
    n_fragments: int                 # retained fragments after filtering
    dropped_fragments: int           # raw_n_fragments - n_fragments
    vertex_counts: list              # vertices per retained fragment
    largest_fraction: float          # largest fragment / total retained vertices
    centroids: np.ndarray            # 3D centroid per retained fragment
    centroid_dist_matrix: np.ndarray # pairwise 3D centroid distances
    fragment_meshes: list            # PyVista mesh per retained fragment
    embeddings: list                 # normalized 2D ISOMAP embedding, or None

    @property
    def is_compact(self):
        return self.n_fragments == 1
```

Fragments are sorted largest to smallest. Fragment 1 (`F1`, index 0) is therefore
always the dominant scar fragment.

---

## Stage 1 - Load the CZ Mesh

```python
cz = pv.read(str(case.tissue_surfaces["core"]))
```

Only Core Zone meshes are analysed. M01 is responsible for scanning the dataset
and creating/loading the known split; M02 reuses that split and runs ISOMAP on
the known cases.

---

## Stage 2 - Split Connected Fragments

```python
labeled = mesh.connectivity(largest=False)
```

`connectivity()` assigns a `RegionId` to each connected mesh island. Depending
on the PyVista mesh, `RegionId` can be stored in `cell_data` or `point_data`, so
the code checks both.

Each region is extracted by thresholding:

```python
piece = labeled.threshold((rid - 0.5, rid + 0.5), scalars="RegionId")
piece = piece.extract_surface().clean()
```

The resulting fragment meshes are sorted by decreasing vertex count. This avoids
point/cell indexing mismatches and guarantees a stable fragment rank.

---

## Stage 3 - Drop Tiny Fragments

```python
kept = [frag for frag in fragments if frag.n_points >= MIN_VERTS_FRAGMENT]
if kept:
    return kept
return fragments[:1]
```

Fragments below 30 vertices are treated as likely segmentation noise. If every
fragment is below the threshold, the largest fragment is still kept so the case
does not disappear silently.

`dropped_fragments` is stored in `CZAnalysis` for diagnostics.

---

## Stage 4 - Run ISOMAP Per Fragment

ISOMAP approximates geodesic distances along the scar surface:

1. Build a k-nearest-neighbour graph in 3D (`k=10`).
2. Compute shortest-path distances through that graph.
3. Use MDS internally to find 2D coordinates that preserve those geodesic
   distances as well as possible.

M02 uses scikit-learn:

```python
def _run_isomap(pts, k=ISOMAP_K, n_pts=ISOMAP_PTS):
    if len(pts) > n_pts:
        idx = np.random.default_rng(SEED).choice(len(pts), n_pts, replace=False)
        pts = pts[idx]
    k_actual = min(k, len(pts) - 1)
    return Isomap(n_neighbors=k_actual, n_components=2).fit_transform(pts)
```

The 800-point cap keeps the pairwise distance matrix tractable. Any fragment
with fewer than `MIN_VERTS_ISOMAP` vertices gets `None` instead of an embedding.
If scikit-learn raises during embedding, the notebook prints the case key and
stores `None` for that fragment.

---

## Stage 5 - Normalize Embeddings

```python
emb = emb - emb.mean(axis=0, keepdims=True)
scale = np.linalg.norm(emb, axis=1).max()
if scale > 0:
    emb = emb / scale
```

Normalization removes translation and scale so embeddings can be compared across
patients and fragments. The embedding axes remain arbitrary: rotation, reflection,
and absolute orientation are not anatomical.

The shape of the point cloud is meaningful; the x/y labels are not.

---

## Stage 6 - Collect Fragment Metadata

For each retained fragment, `analyze_cz(case)` stores:

| Field | Meaning |
|---|---|
| `vertex_counts` | Size proxy for every retained fragment |
| `centroids` | 3D centre of mass for each fragment |
| `centroid_dist_matrix` | Pairwise Euclidean distances between fragment centroids |
| `largest_fraction` | Compactness score; `1.0` means a single retained scar island |
| `fragment_meshes` | Meshes used for 3D diagnostics |
| `embeddings` | Normalized per-fragment ISOMAP embeddings |

`centroid_dist_matrix` is intentionally not part of the ISOMAP embedding. It
describes how far apart fragments are in physical LV space, while ISOMAP
describes the intrinsic shape of each fragment.

---

## Notebook Outputs

| Section | Output | Description |
|---|---|---|
| 1a | console summary | Number of analysed, compact, and fragmented known-set cases |
| 1b | `isomap_grid_known.png` | Known-set grid of normalized embeddings; compact cases have grey borders, fragmented cases orange borders |
| 1c | `isomap_3d_vs_embedding_diagnostic.png` | Side-by-side 3D fragment renderings and ISOMAP plots for selected fragmented cases |
| 1d | `fragmentation_summary.png` | Fragment count, compactness, and top fragmented cases |
| 1d | `top_fragmented_3d.png` | 3D renders of the most fragmented cases |
| 1e | `xyz_candidates.png` | Candidate F1 embeddings with three-arm/XYZ-like radial signatures |

The file names above come from the current `savefig()` calls in
`M02_isomaps.ipynb`.

---

## Interpreting the Embedding

| Embedding shape | Scar interpretation |
|---|---|
| Compact disc | One continuous blob without a dominant elongation |
| Elongated band | Scar with one main geometric axis |
| Star / XYZ-like shape | Three elongated arms in the dominant fragment |
| Multiple plotted lanes | Multiple retained fragments in one patient |
| `x` marker / `None` | Fragment too small or failed to embed |

Fragment colours identify size rank, not anatomy:

| Label | Colour |
|---|---|
| `F1` | blue |
| `F2` | orange |
| `F3` | green |
| `F4` | red |

---

## Relationship to M01 and M03

M01 establishes the dataset inventory and reproducible known/held-out split.
M02 uses that split to quantify CZ shape and fragmentation. M03 then validates
the anatomical reference frame and builds polar maps that describe scar location.

Together:

- **M01:** what data is available and what cases are in the known split
- **M02:** what shape and fragmentation pattern each CZ scar has
- **M03:** where scar lies on the LV, with vertex, area, and coverage projections

---

## Source Modules Checked

| Module | Relevant role |
|---|---|
| `src.data_loading` | Dataset discovery, CZ/LV/RV path resolution, `PatientCase` model |
| `src.mesh_utils` | Shared mesh loading plus connected-fragment splitting/filtering |
| `src.lv_geometry` | Long-axis, septal reference, transforms, and polar projection geometry |
| `src.clinical_data` | Clinical CSV loading and sex-label normalization |
| `src.cone_bspline_simple` | Cone/polar B-spline scar-envelope utilities built on the shared geometry and fragment helpers |

---

## Caveats

- The fragment filter is vertex-count based; an area-based threshold would be
  more physically meaningful.
- ISOMAP embeddings are shape descriptors only. They do not carry anatomical
  orientation.
- The notebook still duplicates some helper functions that now exist in
  `src.mesh_utils`; future cleanup can import the shared helpers directly.
- Full execution requires the project/Jupyter environment with the scientific
  stack installed and access to the external dataset root.
