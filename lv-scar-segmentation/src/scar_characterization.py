"""Reusable scar-map characterization helpers.

Characterization turns polar maps into descriptors: AHA/territory/ring burden,
dominant region, summary tables, and low-dimensional pattern scores. It should
not decide whether an association is meaningful; that belongs in
``scar_analysis``.
"""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def patient_id_from_key(key: str) -> str:
    """Return the patient id from the notebook's ``year/patient_id`` key format.

    The notebooks index cases by ``"<year>/<patient_id>"``; the id after the
    last ``/`` is the per-patient identifier used in the output tables.
    """

    return str(key).split("/")[-1]


def region_burden_table(
    maps: Mapping[str, np.ndarray],
    region_grid: np.ndarray,
    labels: Sequence,
    *,
    prefix: str = "",
    denominator_maps: Mapping[str, np.ndarray] | None = None,
    denominator_grid: np.ndarray | None = None,
    normalize_by_region_size: bool = True,
    include_raw: bool = False,
) -> pd.DataFrame:
    """Summarize each polar map inside named regions, one patient per row.

    For every map, each region's cells are reduced to a single burden value:
    the sum of map values in the region divided by a denominator. The
    denominator is chosen in priority order:

    1. ``denominator_maps[key]`` summed over the region (per-patient denominator,
       e.g. normalize scar by that patient's LV surface area);
    2. ``denominator_grid`` summed over the region (one shared denominator grid);
    3. the number of finite cells in the region (plain regional mean).

    The thesis characterized the same scar four ways through this helper, each
    applied at AHA-17 segment, coronary-territory, and radial-ring scales (set by
    the ``region_grid`` argument; the ring scale uses only the vertex and
    coverage forms). The four measures:

    * Vertex occupancy: binary scar-vertex presence per polar bin; with no
      denominator the burden is scar bins / polar bins in the region (priority 3,
      region size in grid cells). On the 64x128 grid an AHA segment spans
      336-2048 bins and a territory ~1.9k-4.4k bins.
    * Raw CZ area: maps hold triangulated scar (CZ) surface area in mm2 per bin,
      summed over the region (``normalize_by_region_size=False``).
    * Scar-normalized area: the CZ-area maps pre-divided by each patient's total
      scar area, so the region value is that patient's share of scar there.
    * LV surface coverage: CZ-area maps divided by per-patient LV-area maps given
      as ``denominator_maps`` (priority 1) -- scar area / available LV area,
      correcting for unequal bin sizes and LV-size differences.

    Args:
        maps: ``key -> polar map`` (2-D array). ``key`` is the ``year/patient``
            key used across the notebooks.
        region_grid: integer/label grid, same shape as the maps, assigning each
            cell to a region.
        labels: region labels to report; each becomes one output column.
        prefix: prepended to region column names (see ``_label_column``).
        denominator_maps: optional per-patient denominator maps (priority 1).
        denominator_grid: optional shared denominator grid (priority 2).
        normalize_by_region_size: if False, store the raw numerator instead of
            the ratio (the denominator is then ignored for the stored value).
        include_raw: also emit ``<col>_raw`` and ``<col>_denom`` audit columns.

    Returns:
        One row per map with columns ``key``, ``pid``, ``total`` (mean over the
        whole map) and one burden column per region label.
    """

    region_grid = np.asarray(region_grid)
    denominator_grid = None if denominator_grid is None else np.asarray(denominator_grid, dtype=float)
    rows: list[dict] = []

    for key, arr in maps.items():
        arr = np.asarray(arr, dtype=float)
        # "total" = whole-map mean, ignoring NaN cells outside the LV.
        row = {"key": key, "pid": patient_id_from_key(key), "total": float(np.nanmean(arr))}

        # Per-patient denominator map wins over the shared grid, when present.
        denom_arr = None
        if denominator_maps is not None and key in denominator_maps:
            denom_arr = np.asarray(denominator_maps[key], dtype=float)

        for label in labels:
            mask = region_grid == label
            values = arr[mask]
            numerator = float(np.nansum(values))

            # Denominator priority: per-patient map > shared grid > cell count.
            if denom_arr is not None:
                denominator = float(np.nansum(denom_arr[mask]))
            elif denominator_grid is not None:
                denominator = float(np.nansum(denominator_grid[mask]))
            else:
                denominator = float(np.isfinite(values).sum())

            col = _label_column(prefix, label)
            if normalize_by_region_size:
                # Guard empty/zero regions so they read NaN, not inf/0-div.
                row[col] = numerator / denominator if denominator > 0 else np.nan
            else:
                row[col] = numerator
            if include_raw:
                row[f"{col}_raw"] = numerator
                row[f"{col}_denom"] = denominator
        rows.append(row)

    return pd.DataFrame(rows)


def region_summary(
    table: pd.DataFrame,
    columns: Sequence[str],
    *,
    labels: Sequence | None = None,
    label_names: Mapping | None = None,
    value_name: str = "value",
) -> pd.DataFrame:
    """Return descriptive cohort summaries for region columns, one row per region.

    Takes the per-patient burden table from :func:`region_burden_table` and
    collapses each region column into cohort statistics: how many patients have
    any scar there, and the all-patient and scar-positive-only distributions.

    Args:
        table: per-patient table (rows = patients, columns include ``columns``).
        columns: region columns to summarize.
        labels: region keys aligned with ``columns`` (defaults to ``columns``);
            used as the ``region`` value and to look up ``label_names``.
        label_names: optional ``label -> human-readable name`` mapping.
        value_name: stem for the statistic column names, e.g. ``"burden"`` gives
            ``mean_burden_all``.

    Returns:
        One row per region with counts and all / positive-only summary stats.
    """

    rows = []
    labels = columns if labels is None else labels
    label_names = {} if label_names is None else label_names
    for label, col in zip(labels, columns):
        vals = pd.Series(table[col]).dropna().astype(float)
        pos = vals[vals > 0]  # scar-positive patients only, for the positive-only mean
        rows.append(
            {
                "region": label,
                "label": label_names.get(label, label),
                "n_patients": int(vals.size),
                "n_with_scar": int((vals > 0).sum()),
                "pct_with_scar": 100 * float((vals > 0).mean()) if vals.size else np.nan,
                f"mean_{value_name}_all": float(vals.mean()) if vals.size else np.nan,
                f"sd_{value_name}_all": float(vals.std(ddof=1)) if vals.size > 1 else np.nan,
                f"median_{value_name}_all": float(vals.median()) if vals.size else np.nan,
                f"mean_{value_name}_positive": float(pos.mean()) if pos.size else np.nan,
            }
        )
    return pd.DataFrame(rows)


def pca_region_table(
    table: pd.DataFrame,
    columns: Sequence[str],
    *,
    n_components: int = 6,
    labels: Sequence | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, PCA, StandardScaler]:
    """Run a standardized PCA over region columns as descriptive pattern scores.

    Patients with any missing region value are dropped, the remaining columns
    are z-scored, and PCA extracts the dominant co-variation patterns across
    regions. Purely descriptive: no inference about what the components mean.

    Args:
        table: per-patient burden table.
        columns: region columns fed into the PCA (the feature space).
        n_components: max components to keep; capped by #patients and #columns.
        labels: region names aligned with ``columns`` for the loadings table
            (defaults to ``columns``).

    Returns:
        ``(pca_df, variance_df, loadings_df, pca, scaler)`` where ``pca_df`` is
        the kept patients plus ``PC*`` score columns, ``variance_df`` holds the
        explained / cumulative variance ratios, ``loadings_df`` the per-region
        component weights, and ``pca`` / ``scaler`` the fitted estimators.
        With fewer than 3 usable patients, returns the (near-empty) frame plus
        empty frames and unfitted estimators.
    """

    labels = columns if labels is None else labels
    pca_df = table.dropna(subset=list(columns)).copy()
    # Need at least 3 patients for PCA to be meaningful; bail out otherwise.
    if len(pca_df) < 3:
        empty = pd.DataFrame()
        return pca_df, empty, empty, PCA(), StandardScaler()

    scaler = StandardScaler()
    x_scaled = scaler.fit_transform(pca_df[list(columns)].to_numpy(dtype=float))
    # Can't ask for more components than samples or features allow.
    n = min(n_components, x_scaled.shape[0], x_scaled.shape[1])
    pca = PCA(n_components=n, random_state=42)  # fixed seed -> reproducible signs/order
    scores = pca.fit_transform(x_scaled)
    pc_cols = [f"PC{i + 1}" for i in range(n)]
    for i, col in enumerate(pc_cols):
        pca_df[col] = scores[:, i]

    variance_df = pd.DataFrame(
        {
            "component": pc_cols,
            "explained_variance_ratio": pca.explained_variance_ratio_,
            "cumulative_variance_ratio": np.cumsum(pca.explained_variance_ratio_),
        }
    )
    loadings_df = pd.DataFrame(pca.components_.T, columns=pc_cols)
    loadings_df.insert(0, "region", list(labels))
    loadings_df.insert(1, "column", list(columns))
    return pca_df, variance_df, loadings_df, pca, scaler


def _label_column(prefix: str, label) -> str:
    """Build the output column name for a region label.

    Integer labels (e.g. AHA segment numbers) are zero-padded to two digits
    when a prefix is given, so ``prefix="seg_", label=4`` -> ``"seg_04"`` and
    columns sort naturally. Without a prefix an integer stays bare (``"4"``).
    String labels (e.g. ``"LAD"``) are just prefixed as-is.
    """

    if isinstance(label, (int, np.integer)):
        return f"{prefix}{int(label):02d}" if prefix else str(int(label))
    return f"{prefix}{label}"
