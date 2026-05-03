# ISOMAP Shape Analysis — Pipeline Reference

## Overview

ISOMAP is used to analyse the **intrinsic geometry** of each CZ scar fragment. Rather than
describing where the scar is (that is the polar map's job), ISOMAP describes **what shape**
the scar has — elongated, compact, branching, etc. — by unfolding the 3D surface into a
meaningful 2D embedding.

The pipeline has two concerns that must be kept separate:
- **Shape** of each individual fragment → ISOMAP embedding per fragment.
- **Spatial layout** of multiple fragments → centroid distance matrix.

---

## Current notebook

The ISOMAP section has been moved out of `notebooks/01_scar_exploration.ipynb` into
`notebooks/01_ISOMAPS.ipynb`.

`01_ISOMAPS.ipynb` is now intended to run standalone:
- imports NumPy, pandas, matplotlib, PyVista, SciPy, scikit-learn, tqdm, and IPython `display`
- detects the project root from either `lv-scar-segmentation/` or `lv-scar-segmentation/notebooks/`
- imports `scan` from `src.data_loading`
- scans `HD_ROOT` (default `F:/RM_TEKNON_DEVELOP`, override with `LV_SCAR_DATA_ROOT`)
- loads the existing canonical split from `results/split.json` when present; otherwise creates it with `SEED=42`
- writes split-out ISOMAP figures to `results/01_ISOMAPS/`
- defines `_best_shell()` locally for the 3D diagnostic cells

The shell Python available during the notebook split did not have `pandas` installed, so the
notebook was syntax-checked but not executed end-to-end from that shell. Use the project/Jupyter
environment with `requirements.txt` installed for a full run.

---

## Key constants

```python
ISOMAP_PTS        = 800   # max vertices subsampled per fragment before ISOMAP
ISOMAP_K          = 10    # number of neighbours in the k-NN graph
MIN_VERTS_FRAGMENT = 30   # minimum vertices to keep a fragment (below = noise)
MIN_VERTS_ISOMAP   = 5    # minimum vertices to attempt embedding (below = None)
```

---

## Output dataclass

```python
@dataclass
class CZAnalysis:
    key:                  str            # "year/patient_id"
    raw_n_fragments:      int            # total fragments before filtering
    n_fragments:          int            # kept after MIN_VERTS_FRAGMENT filter
    dropped_fragments:    int            # raw - kept  (noise islands removed)
    vertex_counts:        list           # number of vertices per kept fragment
    largest_fraction:     float          # vertex_counts[0] / total  → compactness [0,1]
    centroids:            np.ndarray     # 3D centre of mass per fragment  (N×3)
    centroid_dist_matrix: np.ndarray     # pairwise Euclidean dist between centroids
    fragment_meshes:      list           # PyVista mesh per fragment
    embeddings:           list           # normalised 2D ISOMAP per fragment (or None)

    @property
    def is_compact(self):
        return self.n_fragments == 1     # True = single solid scar
```

`embeddings[i]` is parallel to `fragment_meshes[i]`. F1 (index 0) is always the largest
fragment. `None` embedding means the fragment was too small to embed.

---

## Stage 1 — Load raw CZ mesh

```python
cz = pv.read(str(case.tissue_surfaces["core"]))
```

One raw PyVista mesh per patient. May contain multiple disconnected islands — either genuine
separate scar fragments or segmentation noise.

---

## Stage 2 — Split into connected fragments  `_connected_fragments(cz)`

**Goal**: separate the raw mesh into one mesh per disconnected island.

**Why**: geodesic distances between disconnected islands are infinite. Running ISOMAP on the
whole mesh at once would corrupt the distance matrix and produce meaningless embeddings. Each
fragment must be embedded independently.

```python
labeled    = mesh.connectivity(largest=False)   # assigns RegionId to every cell/point
region_ids = np.unique(labeled.cell_data['RegionId'])

fragments = []
for rid in region_ids:
    piece = labeled.threshold((rid - 0.5, rid + 0.5), scalars='RegionId')
    piece = piece.extract_surface().clean()
    if piece.n_points > 0:
        fragments.append(piece)

fragments.sort(key=lambda m: m.n_points, reverse=True)   # largest first → F1, F2, ...
```

`connectivity()` labels regions but RegionId may live on cell data or point data depending
on the mesh — the code checks both. Sorted largest→smallest so F1 is always the dominant
fragment.

**Returns**: list of PyVista meshes, one per connected island.

---

## Stage 3 — Drop noise fragments  `_filter_fragments(fragments)`

**Goal**: remove tiny islands that are almost certainly segmentation artefacts.

```python
MIN_VERTS_FRAGMENT = 30

kept = [frag for frag in fragments if frag.n_points >= MIN_VERTS_FRAGMENT]
if not kept:
    return fragments[:1]    # always keep at least one fragment
```

`dropped_fragments = raw_n_fragments - n_fragments` is stored in `CZAnalysis` for diagnostics.

**Known limitation**: vertex-count filtering is a rough proxy. Area-based filtering would be
more principled but is not currently implemented.

**Returns**: cleaned fragment list.

---

## Stage 4 — Per-fragment ISOMAP embedding  `_run_isomap(pts)`

### What ISOMAP does

ISOMAP unfolds the intrinsic geometry of a 3D surface into 2D coordinates. The shape of the
2D point cloud reflects the true shape of the scar surface — elongated arms in the embedding
mean elongated arms in the real scar, regardless of its orientation in 3D space.

ISOMAP = **k-NN graph → geodesic distances → MDS**.

### Step A — k-NN graph

Each vertex is connected by an edge to its `k=10` nearest neighbours measured by Euclidean
distance. Because neighbours on a surface are close together, each edge is a tiny step *along*
the surface. The result is a graph that hugs the surface geometry.

**Why k=10?**
- Too small (k=2): graph breaks apart, some point pairs have no path → infinite distances.
- Too large (k=50): edges start cutting through 3D space rather than following the surface →
  geodesic approximation degrades.

### Step B — geodesic distances (shortest paths)

Dijkstra's algorithm computes the shortest path between every pair of points through the k-NN
graph. These path lengths approximate the true geodesic distance along the surface — the
distance you would measure if you laid a tape measure along the surface between two points.

**Key contrast with Euclidean distance**: two points on opposite sides of a fold in the scar
may be close in 3D space (small Euclidean distance) but far along the surface (large geodesic
distance). ISOMAP uses the geodesic — plain MDS on Euclidean distances would wrongly place
them near each other.

### Step C — MDS (Multidimensional Scaling)

MDS answers: *given a table of pairwise distances between N points, find 2D coordinates such
that the distances between the dots match the table as closely as possible.*

It does this via eigendecomposition of the double-centred squared distance matrix D². The top
2 eigenvectors give the best 2D layout — the one that preserves the most distance information
in 2 dimensions.

**ISOMAP's only trick**: feed MDS geodesic distances instead of Euclidean ones.

### Code

```python
ISOMAP_PTS = 800
ISOMAP_K   = 10

def _run_isomap(pts, k=ISOMAP_K, n_pts=ISOMAP_PTS):
    # 1. subsample if needed  (ISOMAP is O(N²) in memory)
    if len(pts) > n_pts:
        idx = np.random.default_rng(SEED).choice(len(pts), n_pts, replace=False)
        pts = pts[idx]

    # 2. run ISOMAP  (sklearn handles k-NN graph + geodesic distances + MDS internally)
    k_actual = min(k, len(pts) - 1)
    return Isomap(n_neighbors=k_actual, n_components=2).fit_transform(pts)
```

**Subsampling**: ISOMAP builds an N×N distance matrix — O(N²) memory. Capping at 800 points
keeps it tractable. The random seed is fixed for reproducibility.

**Returns**: raw 2D embedding array (N×2).

---

## Stage 5 — Normalise embedding  `_normalise_embedding(emb)`

**Goal**: remove arbitrary translation and scale so embeddings are comparable across patients
and fragments of different sizes.

```python
def _normalise_embedding(emb):
    emb   = np.asarray(emb, dtype=float)
    emb   = emb - emb.mean(axis=0, keepdims=True)   # centre at origin
    scale = np.linalg.norm(emb, axis=1).max()        # furthest point from centre
    if scale > 0:
        emb = emb / scale                            # scale to unit radius
    return emb
```

After normalisation every embedding fits inside a unit circle. ISOMAP axes are **arbitrary**
(rotation and reflection are not fixed) — the shape of the point cloud is meaningful, not its
absolute orientation.

**Returns**: normalised 2D embedding in unit circle.

---

## Stage 6 — Collect fragment metadata  inside `analyze_cz()`

```python
for frag in fragments:
    pts = np.array(frag.points)
    vertex_counts.append(len(pts))
    centroids.append(pts.mean(axis=0))              # 3D centre of mass

    if len(pts) >= MIN_VERTS_ISOMAP:
        embeddings.append(_normalise_embedding(_run_isomap(pts)))
    else:
        embeddings.append(None)                     # too small to embed

centroids = np.array(centroids)
dist_mat  = cdist(centroids, centroids)             # pairwise Euclidean distances (mm)

largest_fraction = vertex_counts[0] / max(sum(vertex_counts), 1)
```

| Field | Meaning |
|-------|---------|
| `vertex_counts` | Size proxy for each fragment |
| `centroids` | 3D centre of mass — where in the LV each fragment sits |
| `centroid_dist_matrix` | How spatially spread the fragments are (real-space mm distances) |
| `largest_fraction` | 1.0 = single compact scar; <1.0 = fragmented |

**Design**: `centroid_dist_matrix` bridges the per-fragment shape analysis (ISOMAP) and the
global spatial layout of the whole scar. The two together give a complete picture of scar
morphology.

---

## End-to-end call

```python
# run for one patient
analysis = analyze_cz(case)

# key derived fields
analysis.is_compact          # True / False
analysis.n_fragments         # how many islands
analysis.largest_fraction    # compactness score
analysis.embeddings[0]       # normalised 2D embedding of the largest fragment (F1)
analysis.centroid_dist_matrix  # pairwise distances between fragment centres
```

---

## Why ISOMAP and not PCA?

PCA on the raw 3D vertex coordinates finds the axes of maximum variance in 3D space. For a
flat patch this works, but for a curved or folded scar it treats all points as if they were
in a flat cloud — the curvature is lost. ISOMAP respects the surface: points that are far
apart along the surface remain far apart in the embedding, even if they are geometrically
close in 3D.

---

## Interpreting the embedding

| Embedding shape | Scar interpretation |
|----------------|---------------------|
| Compact disc | Single continuous blob with no dominant elongation |
| Elongated band | Scar with one main axis — a stripe along the wall |
| Star / XYZ shape | Three elongated arms — detected by `_f1_xyz_signature()` |
| Multiple separated clouds | Fragmented scar — separate islands in the same plot lane |
| `None` | Fragment too small to embed (< `MIN_VERTS_ISOMAP` vertices) |

The embedding axes are **not anatomical**. They do not correspond to X/Y/Z in the scanner
frame. A star shape means three geometric arms, not three Cartesian directions. The colour
of each fragment (F1=blue, F2=orange, F3=green, F4=red) identifies fragment rank by size, not
anatomy.

---

## Relationship to the Polar Map

| | Polar Map | ISOMAP |
|--|-----------|--------|
| Question answered | Where is the scar? | What shape is the scar? |
| Input | CZ mesh + LV shell | CZ mesh only |
| Output | 64×128 boolean grid | 2D point cloud per fragment |
| Coordinate system | Anatomically grounded (septum=0, apex=centre) | Arbitrary (shape only) |
| Handles fragments? | No — projects all CZ vertices together | Yes — per fragment |
| Limitation | Apical scars distorted by cylindrical projection | Axes not anatomically labelled |
