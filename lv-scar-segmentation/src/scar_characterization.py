"""Reusable scar-map characterization helpers.

Characterization turns polar maps into descriptors: AHA/territory/ring burden,
dominant region, summary tables, and low-dimensional pattern scores. It should
not decide whether an association is meaningful; that belongs in
``scar_analysis``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


@dataclass(frozen=True)
class RegionSpec:
    """Named masks that partition or label a polar map."""

    name: str
    labels: Sequence
    grid: np.ndarray


def patient_id_from_key(key: str) -> str:
    """Return patient id from the notebook's ``year/patient_id`` key format."""

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
    """Summarize each map inside named regions."""

    region_grid = np.asarray(region_grid)
    denominator_grid = None if denominator_grid is None else np.asarray(denominator_grid, dtype=float)
    rows: list[dict] = []

    for key, arr in maps.items():
        arr = np.asarray(arr, dtype=float)
        row = {"key": key, "pid": patient_id_from_key(key), "total": float(np.nanmean(arr))}

        denom_arr = None
        if denominator_maps is not None and key in denominator_maps:
            denom_arr = np.asarray(denominator_maps[key], dtype=float)

        for label in labels:
            mask = region_grid == label
            values = arr[mask]
            numerator = float(np.nansum(values))

            if denom_arr is not None:
                denominator = float(np.nansum(denom_arr[mask]))
            elif denominator_grid is not None:
                denominator = float(np.nansum(denominator_grid[mask]))
            else:
                denominator = float(np.isfinite(values).sum())

            col = _label_column(prefix, label)
            if normalize_by_region_size:
                row[col] = numerator / denominator if denominator > 0 else np.nan
            else:
                row[col] = numerator
            if include_raw:
                row[f"{col}_raw"] = numerator
                row[f"{col}_denom"] = denominator
        rows.append(row)

    return pd.DataFrame(rows)


def dominant_region_table(
    maps: Mapping[str, np.ndarray],
    region_grid: np.ndarray,
    labels: Sequence[str],
    *,
    threshold: float = 0.02,
    denominator_maps: Mapping[str, np.ndarray] | None = None,
    value_prefix: str = "",
) -> pd.DataFrame:
    """Compute per-patient regional burden, dominant region, and involvement."""

    table = region_burden_table(
        maps,
        region_grid,
        labels,
        prefix=value_prefix,
        denominator_maps=denominator_maps,
        normalize_by_region_size=True,
        include_raw=True,
    )
    value_cols = [_label_column(value_prefix, label) for label in labels]

    if len(table) == 0:
        return table

    table["dominant"] = table[value_cols].idxmax(axis=1).str.removeprefix(value_prefix)
    table["n_involved"] = (table[value_cols] >= threshold).sum(axis=1).astype(int)
    return table


def region_summary(
    table: pd.DataFrame,
    columns: Sequence[str],
    *,
    labels: Sequence | None = None,
    label_names: Mapping | None = None,
    value_name: str = "value",
) -> pd.DataFrame:
    """Return descriptive cohort summaries for region columns."""

    rows = []
    labels = columns if labels is None else labels
    label_names = {} if label_names is None else label_names
    for label, col in zip(labels, columns):
        vals = pd.Series(table[col]).dropna().astype(float)
        pos = vals[vals > 0]
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
    """Standardized PCA over region columns as descriptive pattern scores."""

    labels = columns if labels is None else labels
    pca_df = table.dropna(subset=list(columns)).copy()
    if len(pca_df) < 3:
        empty = pd.DataFrame()
        return pca_df, empty, empty, PCA(), StandardScaler()

    scaler = StandardScaler()
    x_scaled = scaler.fit_transform(pca_df[list(columns)].to_numpy(dtype=float))
    n = min(n_components, x_scaled.shape[0], x_scaled.shape[1])
    pca = PCA(n_components=n, random_state=42)
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
    if isinstance(label, (int, np.integer)):
        return f"{prefix}{int(label):02d}" if prefix else str(int(label))
    return f"{prefix}{label}"
