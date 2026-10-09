"""LV long axis, septal reference and polar (bullseye) projection.

All inputs are PyVista meshes in the ADAS3D (scanner) frame, in mm. Every
notebook that builds a polar map or aligns an LV in 3-D goes through these
functions; notebooks only add case-specific wrapping on top.

Pipeline
--------
1. ``long_axis``: apex point and unit axis pointing apex -> base.
2. Septal reference: a unit vector perpendicular to the axis that fixes
   theta = 0.

   - Method C, ``rv_reference_from_rv_mesh``: LV centroid -> RV centroid in
     the basal band. **The method to use.**
   - Method B, ``rv_reference_pca``: LV-only basal-ring PCA. Exploratory.
   - Method A, "dip midpoint" (RV insertion points): exploratory, only in
     M03 (``_rv_reference_dip``), not here.

   C is the default of ``prepare_polar_geometry`` and ``compute_transform``
   (``ref_fn=None``), built from their ``rv_mesh`` argument; every case has
   an RV mesh. Pass ``ref_fn=rv_reference_pca`` to get B explicitly.
   Until 2026-10-09 the default was B, so the thesis results of F01
   ``_build_pmap``, M05 and M06 were made with B (see CLEANUP_PLAN).

3a. Polar maps: ``prepare_polar_geometry`` -> ``project_points_to_polar``
    -> ``polar_mask``.
3b. 3-D alignment (apex at origin, axis on +Z, septum on +X):
    ``compute_transform`` -> ``apply_transform``.

``ref_fn`` chooses the septal method: any callable
``ref_fn(lv_mesh, axis, apex) -> vector``. ``rv_reference_pca`` already has
that shape; ``rv_mesh_ref_fn(rv_mesh)`` gives Method C that shape.

Used by
-------
Polar pipeline: F01, M03, M04_scar_analysis, M05, M06. ``compute_transform``,
``apply_transform``, ``rotation_about_z``, ``plane_basis``,
``oriented_septal_reference``: M01, M03. ``_hull_apex_candidate``: F01, M03
(diagnostics of the long axis). ``rv_mesh_ref_fn``: default ``ref_fn``; notebooks
are to use it instead of their own lambdas (see CLEANUP_PLAN, Phase 2).
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import ConvexHull

# Single tolerance for the whole module: guards divisions by a norm or span
# and decides when a vector counts as zero length.
EPS = 1e-9


def _hull_apex_candidate(
    pts: np.ndarray, center: np.ndarray, axis: np.ndarray
) -> tuple[np.ndarray, float, float]:
    """Find the apex along one candidate axis and score how pointed it is.

    Eccentricity, in basic terms: slide along ``axis`` and look at the two
    ends of the mesh. ``half_span`` is half the mesh length along the axis.
    The apex is the end point farthest from the centroid. Eccentricity is
    that distance divided by ``half_span``.

    - Symmetric direction (e.g. across the LV cavity): the centroid sits in
      the middle, so eccentricity is about 1.
    - The true long axis: the LV is a cone, so most of its mass, and its
      centroid, is near the wide base. The apex is far from the centroid and
      eccentricity is above 1.

    So the candidate with the highest eccentricity is the cone's axis.

    Only convex-hull vertices are tested as apex: the apex is an extreme
    point, so it is always on the hull, and the hull is far fewer points than
    the mesh. ``half_span`` still uses every point, since it only needs the
    min and max projection.

    Parameters
    ----------
    pts : (N, 3) ndarray
        Mesh points, mm.
    center : (3,) ndarray
        Centroid of ``pts``.
    axis : (3,) ndarray
        Unit candidate axis. Its sign does not matter.

    Returns
    -------
    apex : (3,) ndarray
        Hull vertex with the largest ``|projection|`` on ``axis``.
    eccentricity : float
        ``|apex projection| / half_span``.
    half_span : float
        Half the extent of ``pts`` along ``axis``, mm.

    Notes
    -----
    Not a bug fix, flagged only: ``argmax(|proj|)`` takes whichever end is
    farther from the centroid, so it can return the base instead of the apex.
    """
    try:
        cand = pts[ConvexHull(pts).vertices]
    except Exception:
        # Degenerate (e.g. flat) point set: test every point instead.
        cand = pts

    proj_all = (pts - center) @ axis
    half_span = (proj_all.max() - proj_all.min()) / 2.0
    proj = np.abs((cand - center) @ axis)
    i = np.argmax(proj)
    eccentricity = float(proj[i]) / (half_span + EPS)
    return cand[i], eccentricity, half_span


def long_axis(mesh, rv_mesh=None) -> tuple[np.ndarray, np.ndarray]:
    """Find the LV apex and the long axis through it.

    1. Candidate axes: LV PC1 and PC2 (SVD of the LV points) and, if
       ``rv_mesh`` is given, RV PC1. RV PC1 helps when the LV mesh is
       truncated and its own PC1 follows the cavity width instead of its
       length.
    2. For each candidate, ``_hull_apex_candidate`` gives an apex and an
       eccentricity. Keep the apex of the most eccentric candidate. The order
       of candidates only matters for exact ties.
    3. The returned axis is not the winning candidate. It is the unit vector
       from the chosen apex to the LV centroid, so it always passes through
       the apex and points apex -> base.

    Parameters
    ----------
    mesh : pyvista.DataSet
        LV shell.
    rv_mesh : pyvista.DataSet, optional
        RV surface; adds RV PC1 as a candidate.

    Returns
    -------
    axis : (3,) ndarray
        Unit vector, apex -> base.
    apex : (3,) ndarray
        Apex point, mm.

    Notes
    -----
    If the apex coincides with the centroid, the axis falls back to LV PC1,
    whose sign is arbitrary.
    """
    pts = np.asarray(mesh.points)
    center = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - center, full_matrices=False)

    candidates = []
    if rv_mesh is not None:
        try:
            rv_pts = np.asarray(rv_mesh.points)
            rv_center = rv_pts.mean(axis=0)
            _, _, rv_vt = np.linalg.svd(rv_pts - rv_center, full_matrices=False)
            candidates.append(rv_vt[0])
        except Exception:
            pass
    candidates += [vt[0], vt[1]]

    results = [_hull_apex_candidate(pts, center, c) for c in candidates]
    apex = max(results, key=lambda r: r[1])[0]

    vec = center - apex
    norm = np.linalg.norm(vec)
    axis = vec / norm if norm > EPS else vt[0]
    return axis, apex


def plane_basis(axis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Orthonormal basis ``(u, w)`` of the plane perpendicular to ``axis``.

    ``(u, w, axis)`` is right-handed, so angles measured from ``u`` towards
    ``w`` run counter-clockwise when looking from the base down to the apex
    (``axis`` pointing at the viewer).

    Parameters
    ----------
    axis : (3,) ndarray
        Unit vector.

    Returns
    -------
    u, w : (3,) ndarray
        ``u = axis x Z`` (``axis x X`` if the axis is within ~25 deg of Z,
        where the cross product would be unstable), ``w = axis x u``.
    """
    secondary = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(axis, secondary)) > 0.9:
        secondary = np.array([1.0, 0.0, 0.0])
    u = np.cross(axis, secondary)
    u /= np.linalg.norm(u)
    w = np.cross(axis, u)
    return u, w


def rv_reference_pca(
    lv_mesh, axis: np.ndarray, apex: np.ndarray, z_lo: float = 0.60, z_hi: float = 0.85
) -> np.ndarray:
    """Method B (exploratory): septal direction from the LV basal ring alone.

    Takes the LV points between ``z_lo`` and ``z_hi`` of the LV height,
    flattens them into the (u, w) plane and finds the ring's main direction
    (PC1). Of the two halves of the ring split by PC1, the arrow points to
    the half whose points are farther from the ring centroid on average.

    Method C (``rv_reference_from_rv_mesh``) is the method to use and the
    default. B is kept for the M03 A/B/C comparison; pass it as ``ref_fn``.

    Parameters
    ----------
    lv_mesh : pyvista.DataSet
        LV shell.
    axis, apex : (3,) ndarray
        From ``long_axis``.
    z_lo, z_hi : float
        Basal band, as fractions of the LV height measured from the apex.

    Returns
    -------
    (3,) ndarray
        Unit vector perpendicular to ``axis``. ``u`` from ``plane_basis`` if
        the band has fewer than 20 points.

    Notes
    -----
    Flagged only: the band is measured from the apex (``z_lo * span``), not
    from the lowest point like ``project_points_to_polar``.
    """
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

    # Mean radius of each half of the ring; point to the larger one.
    proj_pc1 = diff @ pc1
    mask_pos = proj_pc1 > 0
    mask_neg = proj_pc1 < 0
    r_pos = np.linalg.norm(diff[mask_pos], axis=1).mean() if mask_pos.any() else np.inf
    r_neg = np.linalg.norm(diff[mask_neg], axis=1).mean() if mask_neg.any() else np.inf
    sign = -1.0 if r_pos < r_neg else 1.0
    direction = sign * pc1
    return direction[0] * u + direction[1] * w


def rv_reference_from_rv_mesh(
    lv_mesh, rv_mesh, axis: np.ndarray, apex: np.ndarray, z_lo: float = 0.50, z_hi: float = 0.90
) -> np.ndarray:
    """Method C (the method to use): septal direction as LV -> RV centroid.

    Both centroids are taken in the basal band ``[z_lo, z_hi]`` of the LV
    height, because the whole RV reaches further towards the apex and would
    tilt the vector. A mesh with fewer than 10 points in the band uses its
    full-mesh centroid instead. For ``ref_fn``, wrap it with
    ``rv_mesh_ref_fn(rv_mesh)``.

    Parameters
    ----------
    lv_mesh, rv_mesh : pyvista.DataSet
        LV shell and RV surface.
    axis, apex : (3,) ndarray
        From ``long_axis``.
    z_lo, z_hi : float
        Basal band, as fractions of the LV height measured from the apex.

    Returns
    -------
    (3,) ndarray
        Unit vector perpendicular to ``axis``, pointing towards the RV.
        ``u`` from ``plane_basis`` if the centroids line up along the axis.
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
    vec_perp = vec - np.dot(vec, axis) * axis
    norm = np.linalg.norm(vec_perp)
    if norm < EPS:
        return plane_basis(axis)[0]
    return vec_perp / norm


def rv_mesh_ref_fn(rv_mesh):
    """Turn Method C into a ``ref_fn(lv_mesh, axis, apex)`` for one RV mesh.

    Parameters
    ----------
    rv_mesh : pyvista.DataSet
        RV surface of the case.

    Returns
    -------
    callable
        Pass as ``ref_fn`` to ``compute_transform`` or
        ``prepare_polar_geometry`` (their default when ``ref_fn`` is None).

    Raises
    ------
    ValueError
        If ``rv_mesh`` is None.
    """
    if rv_mesh is None:
        raise ValueError("Method C (default septal reference) needs an RV mesh")
    return lambda lv_mesh, axis, apex: rv_reference_from_rv_mesh(lv_mesh, rv_mesh, axis, apex)


def oriented_septal_reference(ref, flip: bool = False) -> np.ndarray:
    """Normalise a septal reference vector, optionally reversing it.

    Parameters
    ----------
    ref : array_like, shape (3,)
        Output of a ``ref_fn``.
    flip : bool
        Reverse the vector (180 deg turn of the bullseye).

    Returns
    -------
    (3,) ndarray
        Unit vector.
    """
    ref = np.asarray(ref, dtype=float)
    if flip:
        ref = -ref
    return ref / (np.linalg.norm(ref) + EPS)


def rotation_about_z(angle_rad: float) -> np.ndarray:
    """3x3 matrix rotating counter-clockwise by ``angle_rad`` about world +Z."""
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def compute_transform(
    lv_mesh,
    align_septum: bool = False,
    flip_septal_direction: bool = False,
    ref_fn=None,
    rv_mesh=None,
) -> tuple[np.ndarray, np.ndarray]:
    """Rigid transform putting the LV in a standard 3-D pose.

    After ``apply_transform(mesh, R, apex)``: apex at the origin, long axis
    on +Z (apex -> base) and, if ``align_septum``, the septal reference on
    +X.

    Parameters
    ----------
    lv_mesh : pyvista.DataSet
        LV shell.
    align_septum : bool
        Also rotate about Z so the septal reference lands on +X.
    flip_septal_direction : bool
        Reverse the septal reference first (septum on -X).
    ref_fn : callable, optional
        ``ref_fn(lv_mesh, axis, apex)``; see module header. None (default)
        means Method C from ``rv_mesh``.
    rv_mesh : pyvista.DataSet
        RV surface. Used by ``long_axis`` and by the default Method C.

    Returns
    -------
    R : (3, 3) ndarray
        Rotation matrix.
    apex : (3,) ndarray
        Apex point, mm; subtract it before rotating.

    Notes
    -----
    Flagged only: if the axis is already almost -Z, ``R`` stays the identity
    instead of a 180 deg turn, leaving the LV upside down.
    """
    axis, apex = long_axis(lv_mesh, rv_mesh=rv_mesh)

    # Rodrigues' formula: rotation taking ``axis`` onto +Z.
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(axis, z)
    s = np.linalg.norm(v)
    c = np.dot(axis, z)
    if s < EPS:
        R = np.eye(3)
    else:
        vx = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
        R = np.eye(3) + vx + vx @ vx * ((1.0 - c) / (s**2))

    if align_septum:
        if ref_fn is None:
            ref_fn = rv_mesh_ref_fn(rv_mesh)
        ref = oriented_septal_reference(
            ref_fn(lv_mesh, axis, apex),
            flip=flip_septal_direction,
        )
        # Septal reference after the first rotation, in the XY plane.
        ref_xy = (R @ ref)[:2]
        if np.linalg.norm(ref_xy) > EPS:
            R = rotation_about_z(-np.arctan2(ref_xy[1], ref_xy[0])) @ R

    return R, apex


def apply_transform(mesh, R: np.ndarray, apex: np.ndarray):
    """Copy of ``mesh`` with points moved by ``(p - apex)`` then rotated by ``R``."""
    transformed = mesh.copy()
    transformed.points = (np.asarray(mesh.points) - apex) @ R.T
    return transformed


def prepare_polar_geometry(
    shell,
    reference_points=None,
    flip_septal_direction: bool = False,
    ref_fn=None,
    rv_mesh=None,
) -> dict:
    """Everything ``project_points_to_polar`` needs for one case.

    Parameters
    ----------
    shell : pyvista.DataSet
        LV shell.
    reference_points : (M, 3) array_like, optional
        Points whose extent along the axis defines s = 0 and s = 1. Defaults
        to the LV shell. F01 and M05 pass the RV points here.
    flip_septal_direction : bool
        Reverse the septal reference (theta = 0 on the opposite side).
    ref_fn : callable, optional
        ``ref_fn(lv_mesh, axis, apex)``; see module header. None (default)
        means Method C from ``rv_mesh``.
    rv_mesh : pyvista.DataSet
        RV surface. Used by ``long_axis`` and by the default Method C.

    Returns
    -------
    dict
        ``axis``, ``apex``, ``u``, ``w`` as from ``long_axis`` and
        ``plane_basis``; ``ref_angle``, angle of the septal reference in the
        (u, w) plane, rad; ``s_min``, ``s_max``, extent of
        ``reference_points`` along the axis from the apex, mm.
    """
    if ref_fn is None:
        ref_fn = rv_mesh_ref_fn(rv_mesh)
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


def project_points_to_polar(points, geom: dict) -> tuple[np.ndarray, np.ndarray]:
    """Cylindrical coordinates of 3-D points around the LV long axis.

    Parameters
    ----------
    points : (N, 3) array_like
        Points to project (e.g. scar mesh), mm.
    geom : dict
        From ``prepare_polar_geometry``.

    Returns
    -------
    theta : (N,) ndarray
        Angle around the axis from the septal reference, rad in [0, 2*pi),
        counter-clockwise in the (u, w) plane.
    s : (N,) ndarray
        Height along the axis, 0 at ``s_min`` (apex) to 1 at ``s_max``
        (base). Points outside are clipped to 0 or 1.
    """
    axis = geom["axis"]
    apex = geom["apex"]
    u = geom["u"]
    w = geom["w"]
    ref_angle = geom["ref_angle"]

    pts = np.asarray(points)
    v = pts - apex
    s_raw = v @ axis
    span = geom["s_max"] - geom["s_min"]
    s = np.clip((s_raw - geom["s_min"]) / (span + EPS), 0.0, 1.0)
    perp = v - np.outer(s_raw, axis)
    theta = (np.arctan2(perp @ w, perp @ u) - ref_angle) % (2.0 * np.pi)
    return theta, s


def polar_mask(theta, s, n_r: int, n_theta: int, roll_deg: float = 0.0) -> np.ndarray:
    """Bin polar points into a binary (s, theta) grid.

    Parameters
    ----------
    theta, s : (N,) ndarray
        Polar coordinates, theta in rad, s in [0, 1].
    n_r, n_theta : int
        Number of rings (rows, row 0 = apex) and of angular sectors (columns).
    roll_deg : float
        Rotates the map by shifting columns; callers pass ``AHA_DISPLAY_DEG``
        (currently 0, so no shift).

    Returns
    -------
    (n_r, n_theta) bool ndarray
        True where at least one point falls in the bin.

    Notes
    -----
    Every caller first applies ``theta -> 2*pi - theta``, and the column
    reversal ``[:, ::-1]`` here undoes it. Together they equal binning the
    original counter-clockwise theta (only points exactly on a bin edge can
    differ). Kept as is so outputs stay identical.
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
