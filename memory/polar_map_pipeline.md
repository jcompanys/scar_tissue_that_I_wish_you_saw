# Polar Map Bull's Eye - Pipeline Reference

## Current diagnostic preference

For visual 3-D diagnostic grids, prefer a direct anatomical view over the AHA bull's-eye display roll:

- Long axis: `+Z = apex -> base`.
- Camera: view from `-Z`, meaning look from apex toward base.
- Septal/RV reference: keep it on screen-right for diagnostic plots.
- Do not apply the extra AHA display roll unless the plot is explicitly meant to match the bull's-eye layout.
- Putting every apex at the same origin is useful for coordinate computation, but not essential for visual thumbnails; camera centering can use LV bounds or centroid.

After `compute_transform(..., align_septum=True)` has placed the septal/RV reference at `+X`, an apex-to-base diagnostic camera should use a camera position on `-Z` and a `viewup` choice that leaves `+X` on screen-right.

## Current source of truth

Shared code should live in source modules, not copied across notebooks:

- `src.mesh_utils`
  - `best_shell(case)`
  - `best_rv_mesh(case)`
  - `get_cz_mesh(case)`
  - `connected_fragments(mesh)`
  - `filter_fragments(fragments)`
- `src.lv_geometry`
  - `long_axis(mesh)`
  - `plane_basis(axis)`
  - `rv_reference_pca(lv_mesh, axis, apex, ...)`
  - `rv_reference_from_rv_mesh(lv_mesh, rv_mesh, axis, apex, ...)`
  - `oriented_septal_reference(ref, flip=False)`
  - `compute_transform(...)`
  - `apply_transform(mesh, R, apex)`
  - `prepare_polar_geometry(...)`
  - `project_points_to_polar(points, geom)`

The old notebook-private helper names may still appear as wrappers for backward compatibility, but new logic should be moved into these modules.

## What the polar map does

The polar map converts a 3-D Core Zone scar mesh into a 2-D cylindrical bull's-eye grid. The grid can be stacked across patients only if every patient's `theta=0` means the same anatomical direction.

Default grid shape:

```python
N_R = 64       # rows: 0 = apex, 63 = base
N_THETA = 128  # columns: 0 = septal/RV reference
```

Core coordinates:

- `s`: normalized apex-to-base position, clipped to `[0, 1]`.
- `theta`: circumferential angle around the LV, measured relative to the chosen septal/RV reference.

The critical operation is:

```python
theta = (raw_angle - ref_angle) % (2 * np.pi)
```

This subtracts the patient-specific angular offset. If `ref_angle` is inconsistent across patients, scars from the same wall land in different angular bins and population density maps become unreliable.

## Standard 3-D grid vs anatomical aligned 3-D

The standard 3-D grid shows meshes in their exported ADAS/scanner coordinates with a common camera. It can look organized because ADAS often exports cases in a similar world coordinate frame.

The anatomical aligned 3-D grid computes a new frame per patient from the mesh:

1. Estimate LV long axis with SVD.
2. Put the apex at the origin.
3. Rotate the long axis to world `+Z`.
4. Optionally rotate around `+Z` so the inferred septal/RV direction lands on `+X`.

This view can look random when the inferred septal/RV reference is unstable, flipped by 180 degrees, or based on a noisy/ambiguous basal shape. That visual randomness matters for polar maps only if the underlying `theta=0` reference is also inconsistent.

## Reference options

### PCA basal-ring reference

`rv_reference_pca(lv_mesh, axis, apex, z_lo=0.60, z_hi=0.85)` estimates the septal/RV direction from the LV shell alone:

1. Extract a basal band of LV points.
2. Project the band into the plane perpendicular to the long axis.
3. Use PC1 of the 2-D ring as the main asymmetric axis.
4. Choose the side with smaller mean radius as the concave septal/RV side.

This is useful when no RV mesh exists, but it can fail if segmentation artifacts dominate the basal ring.

### RV-mesh reference

`rv_reference_from_rv_mesh(lv_mesh, rv_mesh, axis, apex, z_lo=0.50, z_hi=0.90)` uses the optional RV mesh:

1. Restrict LV and RV points to a basal long-axis band.
2. Compute LV and RV centroids in that band.
3. Use the perpendicular component of `RV_centroid - LV_centroid` as the septal/RV direction.

This is often more anatomically direct than LV-only PCA, but it requires reliable RV surfaces. `best_rv_mesh(case)` now loads `rv`, `rv_endo`, or `rv_epi` when available.

### ADAS/scanner coordinates

ADAS coordinates may visually work better if the exported world frame is already consistent across cases. They are a useful diagnostic baseline. They should not be assumed correct for population polar maps until the angular reference is checked with the diagnostics below.

## Flip support

The code has explicit flip support:

- `oriented_septal_reference(ref, flip=False)` in `src.lv_geometry`
- `compute_transform(..., flip_septal_direction=False)`
- `prepare_polar_geometry(..., flip_septal_direction=False)`

The notebooks may expose this through a `FLIP_SEPTAL_DIRECTION` constant or wrapper. Use it only after diagnostics show a systematic 180-degree inversion.

## How to check whether theta=0 is reliable

Use `notebooks/06_polar_diagnostics.ipynb`.

Diagnostics:

- D1 basal-ring visual: checks whether the selected reference arrow points to the expected basal concavity/RV side.
- D2 angular disagreement: compares reference methods and flags large angular differences.
- D3 ADAS direction clustering: checks whether raw ADAS/scanner directions already cluster consistently.
- D4 three-panel 3-D grid: compares raw, anatomical, and polar-frame views for many cases.
- D5 split-half polar density: splits cases and checks whether population density remains stable.

Generated figures live in `lv-scar-segmentation/results/06_polar_diagnostics/`.

## End-to-end call sequence

```python
shell = best_shell(case)
cz = get_cz_mesh(case)

geom = prepare_polar_geometry(
    shell,
    flip_septal_direction=False,
)

theta, s = project_points_to_polar(cz.points, geom)

grid = np.zeros((N_R, N_THETA), dtype=bool)
ri = np.clip((s * N_R).astype(int), 0, N_R - 1)
ti = np.clip((theta / (2 * np.pi) * N_THETA).astype(int), 0, N_THETA - 1)
grid[ri, ti] = True
```

For an RV-based reference, pass a small wrapper as `ref_fn` that closes over the RV mesh:

```python
rv = best_rv_mesh(case)

def ref_from_rv(shell, axis, apex):
    return rv_reference_from_rv_mesh(shell, rv, axis, apex)

geom = prepare_polar_geometry(shell, ref_fn=ref_from_rv)
```

## Display convention and open caveat

The intended bull's-eye display convention is:

- Septal/RV reference at 9 o'clock.
- Anterior at 12 o'clock.
- Lateral at 3 o'clock.
- Inferior at 6 o'clock.

Open caveat: if `_draw_bullseye()` uses `ax.set_theta_zero_location("W")` together with `ax.set_theta_direction(1)`, then positive `pi/2` is rendered at the bottom, not the top. That conflicts with the intended AHA display text. Before final population figures, verify or patch the display direction so the rendered labels match the mathematical bins.

## Known limitations

- Apical scars spread over the inner disc because cylindrical projection distorts the apex.
- Vertex-count binning can leave empty cells for sparse meshes; area-based binning would be more principled.
- PCA basal-ring reference can fail in noisy basal rings or when the LV-only shape does not clearly encode the septal/RV side.
- RV-mesh reference can fail if RV meshes are missing, incomplete, or in a different coordinate frame from the LV.
