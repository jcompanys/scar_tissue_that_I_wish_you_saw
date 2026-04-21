from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyvista as pv
from scipy.interpolate import splprep, splev
from scipy.ndimage import binary_closing, binary_erosion, binary_fill_holes
from scipy.spatial import ConvexHull, Delaunay, cKDTree


@dataclass
class SplineFit2D:
    xy: np.ndarray
    hull_xy: np.ndarray | None
    spline_xy: np.ndarray | None
    fit_ok: bool
    fit_error: str | None
    n_points: int


@dataclass
class CaseConeSplineResult:
    key: str
    whole: SplineFit2D
    fragments: list[SplineFit2D]
    raw_n_fragments: int
    kept_n_fragments: int
    dropped_fragments: int
    largest_fraction: float


ALPHA_RADIUS_SCALE = 1.55
BSPLINE_SMOOTHING = 0.0015
BSPLINE_SAMPLES = 300
WHOLE_GRID_SIZE = 220
WHOLE_CLOSE_ITERS = 2


def _best_shell(case):
    for key in ("lv", "myocardium", "endo", "epi"):
        path = case.anatomy.get(key)
        if path is not None and path.exists():
            return pv.read(str(path))
    return None


def get_cz_mesh(case):
    path = case.tissue_surfaces.get("core")
    if path is None or not path.exists():
        return None
    return pv.read(str(path))


def _connected_fragments(mesh):
    labeled = mesh.connectivity(largest=False)
    if "RegionId" in labeled.cell_data:
        region_ids = np.unique(np.asarray(labeled.cell_data["RegionId"]).astype(int))
    elif "RegionId" in labeled.point_data:
        region_ids = np.unique(np.asarray(labeled.point_data["RegionId"]).astype(int))
    else:
        return [mesh]

    fragments = []
    for rid in region_ids:
        piece = labeled.threshold((rid - 0.5, rid + 0.5), scalars="RegionId")
        piece = piece.extract_surface().clean()
        if piece.n_points > 0:
            fragments.append(piece)
    fragments.sort(key=lambda m: m.n_points, reverse=True)
    return fragments


def _filter_fragments(fragments, min_points=30):
    kept = [frag for frag in fragments if frag.n_points >= min_points]
    if kept:
        return kept
    return fragments[:1]


def _long_axis(mesh):
    pts = np.asarray(mesh.points)
    center = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - center, full_matrices=False)
    axis = vt[0]
    proj = (pts - center) @ axis
    apex = pts[proj.argmin()]
    if np.dot(center - apex, axis) < 0:
        axis = -axis
    return axis, apex


def _plane_basis(axis):
    secondary = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(axis, secondary)) > 0.9:
        secondary = np.array([1.0, 0.0, 0.0])
    u = np.cross(axis, secondary)
    u /= np.linalg.norm(u)
    w = np.cross(axis, u)
    return u, w


def _rv_reference_pca(lv_mesh, axis, apex, z_lo=0.60, z_hi=0.85):
    u, w = _plane_basis(axis)
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


def _prepare_polar_geometry(shell, reference_points):
    axis, apex = _long_axis(shell)
    u, w = _plane_basis(axis)
    ref = _rv_reference_pca(shell, axis, apex)
    ref_angle = np.arctan2(np.dot(ref, w), np.dot(ref, u))
    ref_points = np.asarray(reference_points)
    ref_s_raw = (ref_points - apex) @ axis
    return {
        "axis": axis,
        "apex": apex,
        "u": u,
        "w": w,
        "ref_angle": ref_angle,
        "s_min": float(ref_s_raw.min()),
        "s_max": float(ref_s_raw.max()),
    }


def project_points_to_polar(points, geom):
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


def polar_to_cone_xy(theta, s):
    return np.column_stack([s * np.cos(theta), s * np.sin(theta)])


def _triangle_circumradius(a, b, c):
    ab = np.linalg.norm(a - b)
    bc = np.linalg.norm(b - c)
    ca = np.linalg.norm(c - a)
    s = 0.5 * (ab + bc + ca)
    area_sq = max(s * (s - ab) * (s - bc) * (s - ca), 0.0)
    if area_sq <= 1e-12:
        return np.inf
    area = np.sqrt(area_sq)
    return (ab * bc * ca) / (4.0 * area)


def _polygon_area(xy):
    if len(xy) < 3:
        return 0.0
    x = xy[:, 0]
    y = xy[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _extract_boundary_loops(edges):
    adjacency = {}
    for i, j in edges:
        adjacency.setdefault(i, []).append(j)
        adjacency.setdefault(j, []).append(i)
    if not adjacency:
        return []

    unused = {tuple(sorted((i, j))) for i, j in edges}
    loops = []

    while unused:
        start_edge = next(iter(unused))
        start, current = start_edge
        loop = [start, current]
        unused.discard(start_edge)
        prev = start

        for _ in range(len(edges) + 10):
            neighbors = adjacency.get(current, [])
            candidates = [n for n in neighbors if n != prev and tuple(sorted((current, n))) in unused]
            if not candidates:
                break
            if len(candidates) == 1:
                nxt = candidates[0]
            else:
                nxt = min(candidates)
            unused.discard(tuple(sorted((current, nxt))))
            if nxt == start:
                break
            loop.append(nxt)
            prev, current = current, nxt

        if len(loop) >= 3:
            loops.append(loop)

    return loops


def _concave_boundary_xy(unique_xy, alpha_radius_scale=ALPHA_RADIUS_SCALE):
    if len(unique_xy) < 4:
        return None

    tree = cKDTree(unique_xy)
    dists, _ = tree.query(unique_xy, k=min(3, len(unique_xy)))
    if dists.ndim == 1 or dists.shape[1] < 2:
        return None
    nn = dists[:, 1]
    finite_nn = nn[np.isfinite(nn)]
    if finite_nn.size == 0:
        return None
    radius_threshold = float(np.median(finite_nn) * alpha_radius_scale)
    if radius_threshold <= 0:
        return None

    tri = Delaunay(unique_xy)
    edge_counts = {}
    for simplex in tri.simplices:
        a, b, c = unique_xy[simplex]
        if _triangle_circumradius(a, b, c) > radius_threshold:
            continue
        tri_edges = [(simplex[0], simplex[1]), (simplex[1], simplex[2]), (simplex[2], simplex[0])]
        for i, j in tri_edges:
            edge = tuple(sorted((int(i), int(j))))
            edge_counts[edge] = edge_counts.get(edge, 0) + 1

    boundary_edges = [edge for edge, count in edge_counts.items() if count == 1]
    loops = _extract_boundary_loops(boundary_edges)
    if not loops:
        return None

    best_xy = None
    best_area = -1.0
    for loop in loops:
        boundary_xy = unique_xy[np.asarray(loop)]
        if len(boundary_xy) < 4:
            continue
        area = _polygon_area(boundary_xy)
        if area > best_area:
            best_area = area
            best_xy = boundary_xy

    if best_xy is None:
        return None
    return best_xy


def _fit_closed_bspline_on_xy(
    xy,
    smoothing=BSPLINE_SMOOTHING,
    n_samples=BSPLINE_SAMPLES,
    alpha_radius_scale=ALPHA_RADIUS_SCALE,
):
    xy = np.asarray(xy, dtype=float)
    if len(xy) < 4:
        return SplineFit2D(
            xy=xy,
            hull_xy=None,
            spline_xy=None,
            fit_ok=False,
            fit_error="not enough points for hull/spline",
            n_points=len(xy),
        )

    unique_xy = np.unique(np.round(xy, 6), axis=0)
    if len(unique_xy) < 4:
        return SplineFit2D(
            xy=xy,
            hull_xy=None,
            spline_xy=None,
            fit_ok=False,
            fit_error="not enough unique points",
            n_points=len(xy),
        )

    try:
        hull_xy = _concave_boundary_xy(unique_xy, alpha_radius_scale=alpha_radius_scale)
        if hull_xy is None:
            hull = ConvexHull(unique_xy)
            hull_xy = unique_xy[hull.vertices]
        if len(hull_xy) < 4:
            return SplineFit2D(
                xy=xy,
                hull_xy=hull_xy,
                spline_xy=None,
                fit_ok=False,
                fit_error="boundary too small",
                n_points=len(xy),
            )
        closed = np.vstack([hull_xy, hull_xy[0]])
        tck, _ = splprep([closed[:, 0], closed[:, 1]], s=smoothing, per=True)
        u_new = np.linspace(0.0, 1.0, n_samples)
        x_new, y_new = splev(u_new, tck)
        spline_xy = np.column_stack([x_new, y_new])
        return SplineFit2D(
            xy=xy,
            hull_xy=hull_xy,
            spline_xy=spline_xy,
            fit_ok=True,
            fit_error=None,
            n_points=len(xy),
        )
    except Exception as exc:
        return SplineFit2D(
            xy=xy,
            hull_xy=None,
            spline_xy=None,
            fit_ok=False,
            fit_error=str(exc),
            n_points=len(xy),
        )


def _rasterize_xy(xy, grid_size=WHOLE_GRID_SIZE):
    xy = np.asarray(xy, dtype=float)
    mins = xy.min(axis=0)
    maxs = xy.max(axis=0)
    span = np.maximum(maxs - mins, 1e-6)
    uv = (xy - mins) / span
    rows = np.clip(np.rint(uv[:, 1] * (grid_size - 1)).astype(int), 0, grid_size - 1)
    cols = np.clip(np.rint(uv[:, 0] * (grid_size - 1)).astype(int), 0, grid_size - 1)
    mask = np.zeros((grid_size, grid_size), dtype=bool)
    mask[rows, cols] = True
    return mask, mins, span


def _ordered_boundary_pixels(boundary):
    coords = [tuple(rc) for rc in np.argwhere(boundary)]
    if len(coords) < 4:
        return None

    coord_set = set(coords)
    neighbors = {}
    for r, c in coords:
        ngh = []
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                if dr == 0 and dc == 0:
                    continue
                cand = (r + dr, c + dc)
                if cand in coord_set:
                    ngh.append(cand)
        neighbors[(r, c)] = ngh

    start = min(coords, key=lambda rc: (rc[1], rc[0]))
    ordered = [start]
    visited = {start}
    current = start
    previous = None

    while len(visited) < len(coords):
        candidates = [p for p in neighbors[current] if p not in visited]
        if candidates:
            if previous is None:
                next_pt = min(candidates, key=lambda rc: (rc[1], rc[0]))
            else:
                direction = np.array([current[0] - previous[0], current[1] - previous[1]], dtype=float)
                if np.linalg.norm(direction) < 1e-9:
                    direction = np.array([0.0, 1.0])
                def _score(pt):
                    step = np.array([pt[0] - current[0], pt[1] - current[1]], dtype=float)
                    return (-np.dot(direction, step), pt[1], pt[0])
                next_pt = min(candidates, key=_score)
        else:
            remaining = [p for p in coords if p not in visited]
            if not remaining:
                break
            next_pt = min(remaining, key=lambda rc: (rc[0] - current[0]) ** 2 + (rc[1] - current[1]) ** 2)
        ordered.append(next_pt)
        visited.add(next_pt)
        previous, current = current, next_pt

    return np.asarray(ordered, dtype=float)


def _fit_whole_envelope_on_xy(xy, smoothing=BSPLINE_SMOOTHING, n_samples=BSPLINE_SAMPLES):
    xy = np.asarray(xy, dtype=float)
    if len(xy) < 4:
        return SplineFit2D(
            xy=xy,
            hull_xy=None,
            spline_xy=None,
            fit_ok=False,
            fit_error="not enough points for whole-envelope spline",
            n_points=len(xy),
        )

    try:
        mask, mins, span = _rasterize_xy(xy)
        structure = np.ones((5, 5), dtype=bool)
        for _ in range(WHOLE_CLOSE_ITERS):
            mask = binary_closing(mask, structure=structure)
        mask = binary_fill_holes(mask)
        boundary = mask & ~binary_erosion(mask, structure=np.ones((3, 3), dtype=bool))
        ordered = _ordered_boundary_pixels(boundary)
        if ordered is None or len(ordered) < 4:
            raise ValueError("failed to trace whole-envelope boundary")

        hull_xy = np.column_stack([
            mins[0] + (ordered[:, 1] / (mask.shape[1] - 1)) * span[0],
            mins[1] + (ordered[:, 0] / (mask.shape[0] - 1)) * span[1],
        ])
        closed = np.vstack([hull_xy, hull_xy[0]])
        tck, _ = splprep([closed[:, 0], closed[:, 1]], s=smoothing, per=True)
        u_new = np.linspace(0.0, 1.0, n_samples)
        x_new, y_new = splev(u_new, tck)
        spline_xy = np.column_stack([x_new, y_new])
        return SplineFit2D(
            xy=xy,
            hull_xy=hull_xy,
            spline_xy=spline_xy,
            fit_ok=True,
            fit_error=None,
            n_points=len(xy),
        )
    except Exception as exc:
        return SplineFit2D(
            xy=xy,
            hull_xy=None,
            spline_xy=None,
            fit_ok=False,
            fit_error=str(exc),
            n_points=len(xy),
        )


def build_case_cone_splines(case, min_fragment_points=30, alpha_radius_scale=ALPHA_RADIUS_SCALE):
    shell = _best_shell(case)
    if shell is None:
        raise ValueError("no LV shell available")

    cz = get_cz_mesh(case)
    if cz is None:
        raise ValueError("no core-zone mesh available")

    geom = _prepare_polar_geometry(shell, cz.points)

    whole_theta, whole_s = project_points_to_polar(cz.points, geom)
    whole_xy = polar_to_cone_xy(whole_theta, whole_s)
    whole_fit = _fit_whole_envelope_on_xy(whole_xy)

    fragments_all = _connected_fragments(cz)
    fragments = _filter_fragments(fragments_all, min_points=min_fragment_points)

    fragment_fits = []
    for frag in fragments:
        theta, s = project_points_to_polar(frag.points, geom)
        xy = polar_to_cone_xy(theta, s)
        fragment_fits.append(_fit_closed_bspline_on_xy(xy, alpha_radius_scale=alpha_radius_scale))

    vertex_counts = [frag.n_points for frag in fragments]
    total_vertices = max(sum(vertex_counts), 1)
    largest_fraction = max(vertex_counts) / total_vertices if vertex_counts else 0.0

    return CaseConeSplineResult(
        key=f"{case.year}/{case.patient_id}",
        whole=whole_fit,
        fragments=fragment_fits,
        raw_n_fragments=len(fragments_all),
        kept_n_fragments=len(fragments),
        dropped_fragments=len(fragments_all) - len(fragments),
        largest_fraction=largest_fraction,
    )
