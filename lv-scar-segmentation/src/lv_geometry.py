from __future__ import annotations

import numpy as np


def _hull_apex_candidate(pts, center, axis):
    """Return (apex, eccentricity) using convex hull vertices along *axis*."""
    try:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(pts)
        cand = pts[hull.vertices]
    except Exception:
        cand = pts

    proj = (cand - center) @ axis
    half_span = (((pts - center) @ axis).max() - ((pts - center) @ axis).min()) / 2.0
    # apex = hull vertex furthest from centroid (max |projection|)
    apex = cand[np.argmax(np.abs(proj))]
    apex_dist = float(np.abs((apex - center) @ axis))
    eccentricity = apex_dist / (half_span + 1e-9)
    return apex, eccentricity, half_span


def long_axis(mesh, rv_mesh=None):
    """Return (axis, apex) for an LV mesh using robust convex-hull apex detection.

    Algorithm:
    1. If rv_mesh given: use RV SVD PC1 as primary candidate (cardiac long axis is
       shared between LV and RV — RV PC1 is often a better estimate when the LV mesh
       is truncated/holey and its own PC1 is dominated by circumferential variance).
    2. LV SVD PC1 and PC2 as additional candidates.
    3. Pick the candidate that yields the highest apex eccentricity
       (||apex - center|| / half_span).
    4. Orient axis so it points apex → base (dot(center-apex, axis) > 0).
    """
    pts = np.asarray(mesh.points)
    center = pts.mean(axis=0)
    _, sv, vt = np.linalg.svd(pts - center, full_matrices=False)

    # Build candidate axes: RV PC1 first (highest priority), then LV PC1 / PC2
    candidates = []
    if rv_mesh is not None:
        try:
            rv_pts = np.asarray(rv_mesh.points)
            rv_center = rv_pts.mean(axis=0)
            _, _, rv_vt = np.linalg.svd(rv_pts - rv_center, full_matrices=False)
            candidates.append(rv_vt[0])
        except Exception:
            pass
    candidates.append(vt[0])
    if len(sv) > 1:
        candidates.append(vt[1])

    best_apex, best_ecc = pts[0], -1.0
    for cand in candidates:
        apex, ecc, _ = _hull_apex_candidate(pts, center, cand)
        if ecc > best_ecc:
            best_ecc = ecc
            best_apex = apex

    # Force axis through apex: direction = apex → centroid (guaranteed to pass tip)
    vec = center - best_apex
    norm = np.linalg.norm(vec)
    best_axis = vec / norm if norm > 1e-9 else vt[0]

    return best_axis, best_apex


def long_axis_qc(mesh):
    """Compute QC metrics for long-axis detection.

    Returns dict:
        apex_eccentricity  : ||apex - center|| / half_span  (flag if < 0.35)
        pc1_pc2_ratio      : sv[0] / sv[1]                  (flag if < 1.4)
        mesh_completeness  : fraction of faces with inward normals (flag if > 0.15)
        flag_bad           : True if any threshold triggered
        reason             : human-readable summary
    """
    pts = np.asarray(mesh.points)
    center = pts.mean(axis=0)
    _, sv, vt = np.linalg.svd(pts - center, full_matrices=False)

    axis = vt[0]
    _, apex_ecc, _ = _hull_apex_candidate(pts, center, axis)

    pc1_pc2_ratio = float(sv[0] / (sv[1] + 1e-9)) if len(sv) > 1 else float("inf")

    # mesh_completeness: fraction of faces whose normal points toward centroid
    mesh_completeness = 0.0
    try:
        normals = np.asarray(mesh.face_normals)
        faces_raw = np.asarray(mesh.faces)
        if faces_raw.ndim == 1 and len(faces_raw) > 0:
            # PyVista flat format: [3, i0, i1, i2, 3, i0, ...]
            n_faces = len(normals)
            stride = faces_raw[0] + 1  # assume uniform (triangles → 4)
            face_verts_idx = faces_raw.reshape(n_faces, stride)[:, 1:]
            face_centers = pts[face_verts_idx].mean(axis=1)
        else:
            face_centers = pts[np.asarray(mesh.faces)].mean(axis=1)
        to_center = center - face_centers
        dots = np.einsum("ij,ij->i", normals, to_center)
        mesh_completeness = float((dots > 0).mean())
    except Exception:
        pass

    reasons = []
    if apex_ecc < 0.35:
        reasons.append(f"apex_eccentricity={apex_ecc:.2f}<0.35")
    if pc1_pc2_ratio < 1.4:
        reasons.append(f"pc1_pc2_ratio={pc1_pc2_ratio:.2f}<1.4")
    if mesh_completeness > 0.15:
        reasons.append(f"mesh_completeness={mesh_completeness:.2f}>0.15")

    return {
        "apex_eccentricity": float(apex_ecc),
        "pc1_pc2_ratio": pc1_pc2_ratio,
        "mesh_completeness": mesh_completeness,
        "flag_bad": len(reasons) > 0,
        "reason": "; ".join(reasons) if reasons else "ok",
    }


def plane_basis(axis):
    """Return two orthonormal vectors spanning the plane perpendicular to axis."""
    secondary = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(axis, secondary)) > 0.9:
        secondary = np.array([1.0, 0.0, 0.0])
    u = np.cross(axis, secondary)
    u /= np.linalg.norm(u)
    w = np.cross(axis, u)
    return u, w


def rv_reference_pca(lv_mesh, axis, apex, z_lo=0.60, z_hi=0.85):
    """Estimate the septal/RV direction from basal-ring PCA."""
    u, w = plane_basis(axis)
    pts = np.asarray(lv_mesh.points)
    proj = (pts - apex) @ axis
    span = proj.max() - proj.min()
    mask = (proj >= z_lo * span) & (proj <= z_hi * span)
    ring = pts[mask]
    if len(ring) < 20:
        return u

    ring_2d = np.column_stack([(ring - apex) @ u, (ring - apex) @ w])
    centroid = ring_2d.mean(axis=0)
    diff = ring_2d - centroid
    _, _, vt = np.linalg.svd(diff, full_matrices=False)
    pc1 = vt[0]

    proj_pc1 = diff @ pc1
    mask_pos = proj_pc1 > 0
    mask_neg = proj_pc1 < 0
    r_pos = np.linalg.norm(diff[mask_pos], axis=1).mean() if mask_pos.any() else np.inf
    r_neg = np.linalg.norm(diff[mask_neg], axis=1).mean() if mask_neg.any() else np.inf
    sign = -1.0 if r_pos < r_neg else 1.0
    direction = sign * pc1
    return direction[0] * u + direction[1] * w


def rv_reference_from_rv_mesh(lv_mesh, rv_mesh, axis, apex, z_lo=0.50, z_hi=0.90):
    """Septal direction from RV centroid relative to LV centroid.

    Restricts both meshes to the basal band [z_lo, z_hi] of the LV long axis
    so the full-mesh RV centroid (which extends apically) doesn't bias the direction.
    Falls back to full-mesh centroids if the band yields too few RV points.
    """
    lv_pts = np.asarray(lv_mesh.points)
    rv_pts = np.asarray(rv_mesh.points)

    proj_lv = (lv_pts - apex) @ axis
    span = proj_lv.max() - proj_lv.min()
    lo, hi = z_lo * span, z_hi * span

    lv_band = lv_pts[(proj_lv >= lo) & (proj_lv <= hi)]
    proj_rv = (rv_pts - apex) @ axis
    rv_band = rv_pts[(proj_rv >= lo) & (proj_rv <= hi)]

    lv_centroid = lv_band.mean(axis=0) if len(lv_band) >= 10 else lv_pts.mean(axis=0)
    rv_centroid = rv_band.mean(axis=0) if len(rv_band) >= 10 else rv_pts.mean(axis=0)

    vec = rv_centroid - lv_centroid
    u, w = plane_basis(axis)
    vec_perp = vec - np.dot(vec, axis) * axis
    norm = np.linalg.norm(vec_perp)
    if norm < 1e-9:
        return u
    return vec_perp / norm


def oriented_septal_reference(ref, flip=False):
    """Normalize a septal reference vector, optionally applying a 180-degree flip."""
    ref = np.asarray(ref, dtype=float)
    if flip:
        ref = -ref
    return ref / (np.linalg.norm(ref) + 1e-12)


def rotation_about_z(angle_rad):
    """Return a 3x3 rotation matrix around world +Z."""
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def compute_transform(
    lv_mesh,
    align_septum=False,
    flip_septal_direction=False,
    ref_fn=rv_reference_pca,
    rv_mesh=None,
):
    """Return (R, apex): apex at origin, long axis +Z, optional septal/RV +X."""
    axis, apex = long_axis(lv_mesh, rv_mesh=rv_mesh)

    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(axis, z)
    s = np.linalg.norm(v)
    c = np.dot(axis, z)
    if s < 1e-9:
        R = np.eye(3)
    else:
        vx = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
        R = np.eye(3) + vx + vx @ vx * ((1.0 - c) / (s**2))

    if align_septum:
        ref = oriented_septal_reference(
            ref_fn(lv_mesh, axis, apex),
            flip=flip_septal_direction,
        )
        ref_aligned = ref @ R.T
        ref_aligned = ref_aligned - np.dot(ref_aligned, z) * z
        ref_norm = np.linalg.norm(ref_aligned)
        if ref_norm > 1e-9:
            ref_aligned /= ref_norm
            R = rotation_about_z(-np.arctan2(ref_aligned[1], ref_aligned[0])) @ R

    return R, apex


def apply_transform(mesh, R, apex):
    """Copy a mesh and transform points with apex translation plus rotation R."""
    transformed = mesh.copy()
    transformed.points = (np.asarray(mesh.points) - apex) @ R.T
    return transformed


def prepare_polar_geometry(
    shell,
    reference_points=None,
    flip_septal_direction=False,
    ref_fn=rv_reference_pca,
    rv_mesh=None,
):
    """Build reusable geometry for cylindrical LV polar projections."""
    axis, apex = long_axis(shell, rv_mesh=rv_mesh)
    u, w = plane_basis(axis)
    ref = oriented_septal_reference(
        ref_fn(shell, axis, apex),
        flip=flip_septal_direction,
    )
    ref_angle = np.arctan2(np.dot(ref, w), np.dot(ref, u))

    s_points = np.asarray(shell.points if reference_points is None else reference_points)
    s_raw = (s_points - apex) @ axis
    return {
        "axis": axis,
        "apex": apex,
        "u": u,
        "w": w,
        "ref_angle": ref_angle,
        "s_min": float(s_raw.min()),
        "s_max": float(s_raw.max()),
    }


def project_points_to_polar(points, geom):
    """Project 3-D points to circumferential theta and normalized apex-base s."""
    axis = geom["axis"]
    apex = geom["apex"]
    u = geom["u"]
    w = geom["w"]
    ref_angle = geom["ref_angle"]

    pts = np.asarray(points)
    v = pts - apex
    s_raw = v @ axis
    span = geom["s_max"] - geom["s_min"]
    s = np.clip((s_raw - geom["s_min"]) / (span + 1e-9), 0.0, 1.0)
    perp = v - np.outer(s_raw, axis)
    theta = (np.arctan2(perp @ w, perp @ u) - ref_angle) % (2.0 * np.pi)
    return theta, s


def polar_mask(theta, s, n_r, n_theta, roll_deg=0.0):
    """Binary polar map with X-axis orientation correction and optional AHA roll.

    X-flip ([:, ::-1]) corrects the theta mirroring that arises from the
    clockwise flip applied in callers before binning.
    roll_deg rotates the bullseye CW to align with AHA display convention
    (pass AHA_DISPLAY_DEG from the notebook config).
    """
    mask = np.zeros((n_r, n_theta), dtype=bool)
    ri = np.clip((s * n_r).astype(int), 0, n_r - 1)
    ti = np.clip((theta / (2.0 * np.pi) * n_theta).astype(int), 0, n_theta - 1)
    mask[ri, ti] = True
    mask = mask[:, ::-1]
    if roll_deg:
        k = int(round(roll_deg / 360.0 * n_theta))
        mask = np.roll(mask, k, axis=1)
    return mask
