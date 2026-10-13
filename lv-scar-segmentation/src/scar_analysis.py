"""Reusable statistical analysis helpers for scar descriptors.

Analysis asks questions about descriptors created elsewhere: uncertainty,
within-patient differences, pairwise post-hoc tests, and associations with
clinical variables such as age or sex.

All helpers take a tidy :class:`pandas.DataFrame` (one patient per row, one
descriptor per column) and return either a plain ``dict`` of test outputs or a
new DataFrame; none mutate their inputs. Random-number helpers seed with a
fixed default (42) so results are reproducible unless the caller passes its own
generator.
"""

from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, pearsonr, spearmanr, wilcoxon


def bootstrap_samples(
    data,
    statistic=np.nanmean,
    n_boot: int = 2000,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Bootstrap a one-dimensional sample after dropping NaNs.

    Draws ``n_boot`` resamples with replacement and evaluates ``statistic`` on
    each, yielding the bootstrap distribution of that statistic.

    Args:
        data: Array-like of values; coerced to float, non-finite entries dropped.
        statistic: Callable reducing a 1-D array to a scalar (default
            ``np.nanmean``), applied to each resample.
        n_boot: Number of bootstrap resamples to draw.
        rng: NumPy random generator; when ``None`` a generator seeded with 42 is
            used so the output is reproducible.

    Returns:
        1-D array of length ``n_boot`` holding ``statistic`` for each resample,
        or an empty array when ``data`` has no finite values.
    """

    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]
    if data.size == 0:
        return np.array([])
    rng = np.random.default_rng(42) if rng is None else rng
    # Each row of idx is one resample: n integer indices drawn with replacement.
    idx = rng.integers(0, data.size, size=(n_boot, data.size))
    return np.apply_along_axis(statistic, 1, data[idx])


def bootstrap_ci(
    data,
    statistic=np.nanmean,
    n_boot: int = 2000,
    ci: float = 95,
    rng: np.random.Generator | None = None,
) -> tuple[float, float]:
    """Return a percentile bootstrap confidence interval for 1-D data.

    Args:
        data: Array-like of values passed through to :func:`bootstrap_samples`.
        statistic: Callable reducing a 1-D array to a scalar (default
            ``np.nanmean``).
        n_boot: Number of bootstrap resamples.
        ci: Central confidence level in percent (e.g. 95 for a 95% CI).
        rng: NumPy random generator; ``None`` seeds with 42 for reproducibility.

    Returns:
        ``(low, high)`` percentile bounds, or ``(nan, nan)`` when ``data`` has
        no finite values.
    """

    boot = bootstrap_samples(data, statistic=statistic, n_boot=n_boot, rng=rng)
    if boot.size == 0:
        return np.nan, np.nan
    # Two-sided interval: split the remaining mass equally into both tails.
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
    rng: np.random.Generator | None = None,
) -> pd.DataFrame:
    """Bootstrap the mean and its CI for each of several region columns.

    Args:
        table: Tidy DataFrame with one row per patient.
        columns: Column names, one per region, to summarise.
        labels: Display names aligned with ``columns``; defaults to the column
            names themselves.
        value_name: Name of the output column holding the per-region mean.
        n_boot: Number of bootstrap resamples per region.
        rng: NumPy random generator; ``None`` seeds with 42 for reproducibility.
            The same generator is reused across regions, so results depend on
            region order.

    Returns:
        One row per region with columns ``region``, ``value_name``, ``ci95_low``
        and ``ci95_high``; the mean is NaN when a region has no finite values.
    """

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
    """Run a Friedman test across paired region columns.

    The Friedman test is the non-parametric analogue of a repeated-measures
    ANOVA: it asks whether the region columns differ, using only patients with
    complete data across every region.

    Args:
        table: Tidy DataFrame with one row per patient.
        columns: Region columns to compare; the test needs at least three.

    Returns:
        Dict with ``n`` (complete-case count), ``statistic`` and ``p_value``.
        Statistic and p-value are NaN when fewer than three complete rows or
        fewer than three columns are available.
    """

    sub = table[list(columns)].dropna()
    # Friedman needs >=3 paired groups and >=3 complete observations.
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
    """Pairwise paired Wilcoxon tests with Bonferroni correction.

    Post-hoc companion to :func:`friedman_regions`: compares every unordered
    pair of region columns on the complete-case rows, then corrects the family
    of p-values for multiple comparisons.

    Args:
        table: Tidy DataFrame with one row per patient.
        columns: Region columns to compare pairwise.
        labels: Display names aligned with ``columns``; defaults to the column
            names themselves.

    Returns:
        One row per region pair with columns ``region_a``, ``region_b``,
        ``mean_diff_a_minus_b``, ``wilcoxon_stat`` and ``p_uncorrected``; when
        any pair was testable, ``p_bonferroni`` and ``significant_bonf_0_05``
        are added. Empty DataFrame when fewer than two columns are given.
    """

    labels = columns if labels is None else labels
    sub = table[list(columns)].dropna()
    rows = []
    # Upper-triangular loop: each unordered pair (i, j>i) is tested once.
    for i, (label_i, col_i) in enumerate(zip(labels, columns)):
        for label_j, col_j in zip(labels[i + 1 :], columns[i + 1 :]):
            x = sub[col_i].to_numpy(dtype=float)
            y = sub[col_j].to_numpy(dtype=float)
            diff = x - y
            if len(diff) == 0:
                stat, p_value = np.nan, np.nan
            elif np.allclose(diff, 0):
                # Wilcoxon is undefined when all differences are zero; the pair
                # is trivially equal, so report a neutral p=1.
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
        # Bonferroni: scale each p by the number of tests, clipped at 1.
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
    """Spearman correlation between one scalar column and many region columns.

    Correlates ``x_col`` (e.g. age) against each region column independently,
    using the complete-case rows of that pair, and optionally Bonferroni-corrects
    across regions.

    Args:
        table: Tidy DataFrame with one row per patient.
        columns: Region columns to correlate against ``x_col``.
        x_col: Scalar column correlated against every region column.
        labels: Display names aligned with ``columns``; defaults to the column
            names themselves.
        bonferroni: When True, add Bonferroni-corrected columns across regions.

    Returns:
        One row per region with columns ``region``, ``n``, ``spearman_rho`` and
        ``p_uncorrected``; when ``bonferroni`` and any rows exist, also
        ``p_bonferroni`` and ``significant_bonf_0_05``. Rho and p are NaN for a
        region with fewer than three complete pairs.
    """

    labels = columns if labels is None else labels
    rows = []
    for label, col in zip(labels, columns):
        sub = table[[x_col, col]].dropna()
        # Spearman is undefined below three paired points.
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
        # Bonferroni: scale each p by the number of regions tested, clipped at 1.
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
    """Spearman and Pearson association between one scalar and many variables.

    Reports both a rank (Spearman) and a linear (Pearson) association between
    ``x_col`` and each column in ``y_columns``, computed on the complete-case
    rows of each pair. No multiple-comparison correction is applied here.

    Args:
        table: Tidy DataFrame with one row per patient.
        y_columns: Variables to associate with ``x_col``.
        x_col: Scalar column correlated against every variable.
        labels: Optional mapping from column name to a display label; columns
            absent from the mapping fall back to their own name.

    Returns:
        One row per variable with columns ``variable``, ``label``, ``n``,
        ``spearman_rho``, ``spearman_p``, ``pearson_r`` and ``pearson_p``.
        Coefficients are NaN when fewer than three complete pairs exist; Pearson
        additionally requires both columns to vary (more than one unique value).
    """

    labels = {} if labels is None else labels
    rows = []
    for col in y_columns:
        sub = table[[x_col, col]].dropna()
        # Spearman needs three paired points.
        if len(sub) >= 3:
            rho, p_s = spearmanr(sub[x_col], sub[col])
        else:
            rho, p_s = np.nan, np.nan
        # Pearson additionally needs variation in both columns, else r is undefined.
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
