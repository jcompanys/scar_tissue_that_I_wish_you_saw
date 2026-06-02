"""Reusable statistical analysis helpers for scar descriptors.

Analysis asks questions about descriptors created elsewhere: uncertainty,
within-patient differences, pairwise post-hoc tests, and associations with
clinical variables such as age or sex.
"""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, pearsonr, spearmanr, wilcoxon


def bootstrap_samples(data, statistic=np.nanmean, n_boot: int = 2000, rng=None) -> np.ndarray:
    """Bootstrap one-dimensional data after dropping NaNs."""

    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]
    if data.size == 0:
        return np.array([])
    rng = np.random.default_rng(42) if rng is None else rng
    idx = rng.integers(0, data.size, size=(n_boot, data.size))
    return np.apply_along_axis(statistic, 1, data[idx])


def bootstrap_ci(data, statistic=np.nanmean, n_boot: int = 2000, ci: float = 95, rng=None) -> tuple[float, float]:
    """Return percentile bootstrap CI for one-dimensional data."""

    boot = bootstrap_samples(data, statistic=statistic, n_boot=n_boot, rng=rng)
    if boot.size == 0:
        return np.nan, np.nan
    alpha = (100 - ci) / 2
    lo, hi = np.nanpercentile(boot, [alpha, 100 - alpha])
    return float(lo), float(hi)


def bootstrap_region_ci(
    table: pd.DataFrame,
    columns: Sequence[str],
    *,
    labels: Sequence | None = None,
    value_name: str = "mean",
    n_boot: int = 2000,
    rng=None,
) -> pd.DataFrame:
    """Bootstrap mean values for region columns."""

    rng = np.random.default_rng(42) if rng is None else rng
    labels = columns if labels is None else labels
    rows = []
    for label, col in zip(labels, columns):
        vals = table[col].to_numpy(dtype=float)
        lo, hi = bootstrap_ci(vals, np.nanmean, n_boot=n_boot, rng=rng)
        rows.append(
            {
                "region": label,
                value_name: float(np.nanmean(vals)) if np.isfinite(vals).any() else np.nan,
                "ci95_low": lo,
                "ci95_high": hi,
            }
        )
    return pd.DataFrame(rows)


def friedman_regions(table: pd.DataFrame, columns: Sequence[str]) -> dict:
    """Run a Friedman test over paired region columns."""

    sub = table[list(columns)].dropna()
    if len(sub) < 3 or len(columns) < 3:
        return {"n": int(len(sub)), "statistic": np.nan, "p_value": np.nan}
    stat, p_value = friedmanchisquare(*[sub[col].to_numpy(dtype=float) for col in columns])
    return {"n": int(len(sub)), "statistic": float(stat), "p_value": float(p_value)}


def pairwise_wilcoxon_regions(
    table: pd.DataFrame,
    columns: Sequence[str],
    *,
    labels: Sequence | None = None,
) -> pd.DataFrame:
    """Pairwise paired Wilcoxon tests with Bonferroni correction."""

    labels = columns if labels is None else labels
    sub = table[list(columns)].dropna()
    rows = []
    for i, (label_i, col_i) in enumerate(zip(labels, columns)):
        for label_j, col_j in zip(labels[i + 1 :], columns[i + 1 :]):
            x = sub[col_i].to_numpy(dtype=float)
            y = sub[col_j].to_numpy(dtype=float)
            diff = x - y
            if len(diff) == 0:
                stat, p_value = np.nan, np.nan
            elif np.allclose(diff, 0):
                stat, p_value = 0.0, 1.0
            else:
                stat, p_value = wilcoxon(x, y, zero_method="zsplit", alternative="two-sided")
            rows.append(
                {
                    "region_a": label_i,
                    "region_b": label_j,
                    "mean_diff_a_minus_b": float(np.nanmean(diff)) if len(diff) else np.nan,
                    "wilcoxon_stat": float(stat) if np.isfinite(stat) else np.nan,
                    "p_uncorrected": float(p_value) if np.isfinite(p_value) else np.nan,
                }
            )
    out = pd.DataFrame(rows)
    if len(out):
        out["p_bonferroni"] = np.minimum(out["p_uncorrected"] * len(out), 1.0)
        out["significant_bonf_0_05"] = out["p_bonferroni"] < 0.05
    return out


def spearman_by_region(
    table: pd.DataFrame,
    columns: Sequence[str],
    x_col: str,
    *,
    labels: Sequence | None = None,
    bonferroni: bool = True,
) -> pd.DataFrame:
    """Spearman correlation between one scalar column and many region columns."""

    labels = columns if labels is None else labels
    rows = []
    for label, col in zip(labels, columns):
        sub = table[[x_col, col]].dropna()
        if len(sub) >= 3:
            rho, p_value = spearmanr(sub[x_col], sub[col])
        else:
            rho, p_value = np.nan, np.nan
        rows.append(
            {
                "region": label,
                "n": int(len(sub)),
                "spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
                "p_uncorrected": float(p_value) if np.isfinite(p_value) else np.nan,
            }
        )
    out = pd.DataFrame(rows)
    if bonferroni and len(out):
        out["p_bonferroni"] = np.minimum(out["p_uncorrected"] * len(out), 1.0)
        out["significant_bonf_0_05"] = out["p_bonferroni"] < 0.05
    return out


def scalar_associations(
    table: pd.DataFrame,
    y_columns: Sequence[str],
    x_col: str,
    *,
    labels: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    """Spearman and Pearson association between one scalar and many variables."""

    labels = {} if labels is None else labels
    rows = []
    for col in y_columns:
        sub = table[[x_col, col]].dropna()
        if len(sub) >= 3:
            rho, p_s = spearmanr(sub[x_col], sub[col])
        else:
            rho, p_s = np.nan, np.nan
        if len(sub) >= 3 and sub[x_col].nunique() > 1 and sub[col].nunique() > 1:
            r, p_p = pearsonr(sub[x_col], sub[col])
        else:
            r, p_p = np.nan, np.nan
        rows.append(
            {
                "variable": col,
                "label": labels.get(col, col),
                "n": int(len(sub)),
                "spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
                "spearman_p": float(p_s) if np.isfinite(p_s) else np.nan,
                "pearson_r": float(r) if np.isfinite(r) else np.nan,
                "pearson_p": float(p_p) if np.isfinite(p_p) else np.nan,
            }
        )
    return pd.DataFrame(rows)

