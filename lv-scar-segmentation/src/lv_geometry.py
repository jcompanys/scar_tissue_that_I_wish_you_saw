from __future__ import annotations

import numpy as np


def long_axis(mesh):
    """Return the apex-to-base unit axis and apex point for an LV mesh."""
    pts = np.asarray(mesh.points)
    center = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - center, full_matrices=False)
    axis = vt[0]
    proj = (pts - center) @ axis
    apex = pts[proj.argmin()]
    if np.dot(center - apex, axis) < 0:
        axis = -axis
    return axis, apex


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
):
    """Return (R, apex): apex at origin, long axis +Z, optional septal/RV +X."""
    axis, apex = long_axis(lv_mesh)

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
):
    """Build reusable geometry for cylindrical LV polar projections."""
    axis, apex = long_axis(shell)
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
