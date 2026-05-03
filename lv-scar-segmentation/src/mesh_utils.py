from __future__ import annotations

import numpy as np


SHELL_PRIORITY = ("lv", "myocardium", "endo", "epi")
RV_PRIORITY = ("rv", "rv_endo", "rv_epi")


def best_shell(case, shell_priority=SHELL_PRIORITY):
    """Load the preferred LV shell mesh for a patient case."""
    import pyvista as pv

    for key in shell_priority:
        path = case.anatomy.get(key)
        if path is not None and path.exists():
            return pv.read(str(path))
    return None


def get_cz_mesh(case):
    """Load the core-zone scar mesh for a patient case."""
    import pyvista as pv

    path = case.tissue_surfaces.get("core")
    if path is None or not path.exists():
        return None
    return pv.read(str(path))


def best_rv_mesh(case):
    """Load the preferred RV mesh for a patient case, or None."""
    import pyvista as pv

    rv_anatomy = getattr(case, "rv_anatomy", {})
    if not rv_anatomy:
        anatomy = getattr(case, "anatomy", {})
        rv_anatomy = {key: anatomy.get(key) for key in RV_PRIORITY if key in anatomy}

    for key in RV_PRIORITY:
        path = rv_anatomy.get(key)
        if path is not None and path.exists():
            return pv.read(str(path))
    return None


def connected_fragments(mesh):
    """Split a PyVista mesh into connected surface fragments, largest first."""
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


def filter_fragments(fragments, min_points=30):
    """Drop tiny connected fragments, while always keeping at least one."""
    kept = [frag for frag in fragments if frag.n_points >= min_points]
    if kept:
        return kept
    return fragments[:1]
