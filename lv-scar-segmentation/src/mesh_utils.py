"""Mesh loading helpers for ``PatientCase`` objects.

Reads the ADAS3D meshes whose paths ``src.data_loading.scan`` resolved, and
splits a scar mesh into its connected pieces. Meshes are returned as read by
PyVista, in the ADAS3D (scanner) frame, in mm; nothing is transformed here.

Main API
--------
best_shell(case, shell_priority=LV_PRIORITY) -> mesh or None
    LV surface mesh, first available of ``LV_PRIORITY``.
get_cz_mesh(case) -> mesh or None
    Core-zone scar surface.
best_rv_mesh(case) -> mesh or None
    RV surface mesh, first available of ``RV_PRIORITY``.
connected_fragments(mesh) -> list of meshes
    Connected pieces of a mesh, largest first.
filter_fragments(fragments, min_points=30) -> list of meshes
    Drops tiny pieces, always keeps at least one.

Used by
-------
``best_shell``, ``get_cz_mesh``: F01, M01, M02_clinical_eda, M03,
M04_scar_analysis, M05. ``best_rv_mesh``: F01, M03, M04_scar_analysis, M05.
``connected_fragments``, ``filter_fragments``: no caller yet; M04 is to use
them instead of its own ``_connected_fragments`` (see CLEANUP_PLAN, Phase 2).
"""

from __future__ import annotations

import numpy as np
import pyvista as pv


# Keys of ``PatientCase.lv_anatomy``, most preferred first. "lv" is the full
# Left Ventricle surface; the myocardium, endo and epi layers are fallbacks
# for patients exported without it.
LV_PRIORITY = ("lv", "myocardium", "endo", "epi")

# Keys of ``PatientCase.rv_anatomy``, most preferred first.
RV_PRIORITY = ("rv", "rv_endo", "rv_epi")


def best_shell(case, shell_priority: tuple[str, ...] = LV_PRIORITY):
    """Load the preferred LV shell mesh for a patient case.

    Parameters
    ----------
    case : PatientCase
        Case from ``src.data_loading.scan``.
    shell_priority : tuple of str, optional
        ``case.lv_anatomy`` keys to try, in order. Pass ``("lv",)`` to accept
        only the Left Ventricle mesh (M01 does this).

    Returns
    -------
    pyvista.DataSet or None
        First mesh whose file exists, or ``None`` if none does.

    Notes
    -----
    The fallback is silent: the caller cannot tell which key was used, so a
    patient without ``Left Ventricle.vtk`` gets e.g. the myocardium mesh.
    """
    for key in shell_priority:
        path = case.lv_anatomy.get(key)
        if path is not None and path.exists():
            return pv.read(str(path))
    return None


def get_cz_mesh(case):
    """Load the core-zone scar mesh for a patient case.

    Parameters
    ----------
    case : PatientCase
        Case from ``src.data_loading.scan``.

    Returns
    -------
    pyvista.DataSet or None
        Mesh at ``case.tissue_surfaces["core"]``, or ``None`` if the path is
        missing or the file does not exist.
    """
    path = case.tissue_surfaces.get("core")
    if path is None or not path.exists():
        return None
    return pv.read(str(path))


def best_rv_mesh(case):
    """Load the preferred RV mesh for a patient case, or None.

    Parameters
    ----------
    case : PatientCase
        Case from ``src.data_loading.scan``.

    Returns
    -------
    pyvista.DataSet or None
        First mesh of ``RV_PRIORITY`` whose file exists, or ``None`` (also
        when the patient has no RV folder and ``case.rv_anatomy`` is empty).
    """
    for key in RV_PRIORITY:
        path = case.rv_anatomy.get(key)
        if path is not None and path.exists():
            return pv.read(str(path))
    return None


def connected_fragments(mesh):
    """Split a PyVista mesh into connected surface fragments, largest first.

    Parameters
    ----------
    mesh : pyvista.DataSet
        Mesh to split, e.g. a core-zone scar surface.

    Returns
    -------
    list of pyvista.PolyData
        One cleaned surface per connected region, sorted by number of points
        (descending). Empty regions are dropped. If ``connectivity`` gives no
        ``RegionId`` array, ``[mesh]`` is returned unchanged.

    Notes
    -----
    ``RegionId`` is looked up in cell data first, then point data, because
    PyVista versions differ in where ``connectivity`` stores it. Each region
    is selected with a threshold of +-0.5 around its integer id.
    """
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
        # threshold returns an UnstructuredGrid; back to a surface, with
        # duplicate points merged.
        piece = piece.extract_surface().clean()
        if piece.n_points > 0:
            fragments.append(piece)
    fragments.sort(key=lambda m: m.n_points, reverse=True)
    return fragments


def filter_fragments(fragments, min_points: int = 30):
    """Drop tiny connected fragments, while always keeping at least one.

    Parameters
    ----------
    fragments : list of pyvista.DataSet
        Output of ``connected_fragments`` (largest first).
    min_points : int, optional
        Fragments with fewer points are dropped. Callers may pass their own
        threshold (M04 is to use 50).

    Returns
    -------
    list of pyvista.DataSet
        Fragments with at least ``min_points`` points, in input order. If none
        qualifies, ``fragments[:1]`` (the largest piece when the input is
        sorted), so a non-empty input never gives an empty result.
    """
    kept = [frag for frag in fragments if frag.n_points >= min_points]
    if kept:
        return kept
    return fragments[:1]
