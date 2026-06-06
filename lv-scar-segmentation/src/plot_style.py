"""Shared plotting style for the LV scar characterization notebooks.

The palette is derived from the M01/M02 notebook style and kept deliberately
formal for thesis figures: restrained clinical colors, white backgrounds,
light grids, and export-ready DPI defaults.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

TFG_COLORS = {
    "primary": "#256D85",
    "secondary": "#D9822B",
    "positive": "#2E8B57",
    "neutral": "#B8C0C8",
    "missing": "#C44E52",
    "male": "#256D85",
    "female": "#B35C9E",
    "line": "#263238",
    "text": "#263238",
    "background": "#FFFFFF",
    "grid": "#B8C0C8",
}

TFG_PALETTE = [
    TFG_COLORS["primary"],
    TFG_COLORS["secondary"],
    TFG_COLORS["positive"],
    TFG_COLORS["female"],
    TFG_COLORS["missing"],
    TFG_COLORS["neutral"],
]

HIST_COLOR = TFG_COLORS["primary"]
MEDIAN_COLOR = TFG_COLORS["secondary"]
SEG_COLOR = TFG_COLORS["secondary"]
NO_SEG_COLOR = TFG_COLORS["primary"]

SEX_COLORS = {1: TFG_COLORS["male"], 2: TFG_COLORS["female"]}

TERRITORY_COLORS = {
    "LAD": "#C83E4D",
    "RCA": "#2F6B9A",
    "LCX": "#3F8F5F",
}

FRAGMENT_COLORS = [
    TFG_COLORS["primary"],
    TFG_COLORS["secondary"],
    TFG_COLORS["positive"],
    TFG_COLORS["missing"],
]

POLAR_CANDIDATE_COLORS = {
    "RV PC1": TFG_COLORS["primary"],
    "LV PC1": TFG_COLORS["secondary"],
    "LV PC2": "#7A68A6",
}

THRESHOLD_COLORS = {
    "low": "#C44E52",
    "mid": TFG_COLORS["primary"],
    "high": TFG_COLORS["positive"],
}

RING_COLORS = {
    "Apex": "#7A68A6",
    "Apical": "#C83E4D",
    "Mid": TFG_COLORS["positive"],
    "Basal": TFG_COLORS["primary"],
}

SCAR_CMAP = "YlOrRd"
DIVERGING_CMAP = "RdBu_r"


def set_tfg_style(context: str = "notebook", font_scale: float = 1.0) -> None:
    """Apply the shared TFG plotting defaults to matplotlib and seaborn."""
    base_font = 11 if context == "notebook" else 10
    plt.rcParams.update(
        {
            "figure.dpi": 125,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.facecolor": TFG_COLORS["background"],
            "figure.facecolor": TFG_COLORS["background"],
            "axes.facecolor": TFG_COLORS["background"],
            "axes.edgecolor": TFG_COLORS["line"],
            "axes.labelcolor": TFG_COLORS["text"],
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "axes.axisbelow": True,
            "grid.color": TFG_COLORS["grid"],
            "grid.alpha": 0.18,
            "grid.linewidth": 0.7,
            "font.family": "DejaVu Sans",
            "font.size": base_font * font_scale,
            "axes.titleweight": "semibold",
            "axes.titlesize": 11 * font_scale,
            "axes.labelsize": 10 * font_scale,
            "legend.frameon": False,
            "legend.fontsize": 9 * font_scale,
            "xtick.color": TFG_COLORS["text"],
            "ytick.color": TFG_COLORS["text"],
            "xtick.labelsize": 9 * font_scale,
            "ytick.labelsize": 9 * font_scale,
            "lines.linewidth": 1.6,
            "patch.linewidth": 0.8,
        }
    )

    try:
        import seaborn as sns

        sns.set_theme(
            context=context,
            style="whitegrid",
            palette=TFG_PALETTE,
            rc={
                "figure.facecolor": TFG_COLORS["background"],
                "axes.facecolor": TFG_COLORS["background"],
                "axes.edgecolor": TFG_COLORS["line"],
                "axes.labelcolor": TFG_COLORS["text"],
                "axes.spines.top": False,
                "axes.spines.right": False,
                "axes.grid": True,
                "grid.color": TFG_COLORS["grid"],
                "grid.alpha": 0.18,
                "grid.linewidth": 0.7,
                "legend.frameon": False,
            },
        )
    except Exception:
        pass


def apply_tfg_style(ax: Any, xgrid: bool = False, grid_axis: str | None = None) -> Any:
    """Format one matplotlib axis with the shared thesis style."""
    axis = grid_axis or ("x" if xgrid else "y")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(True, axis=axis, alpha=0.18, linewidth=0.7, color=TFG_COLORS["grid"])
    ax.tick_params(labelsize=9, colors=TFG_COLORS["text"])
    ax.xaxis.label.set_color(TFG_COLORS["text"])
    ax.yaxis.label.set_color(TFG_COLORS["text"])
    ax.title.set_color(TFG_COLORS["text"])
    return ax


def annotate_panel(
    ax: Any,
    text: str,
    x: float = 0.03,
    y: float = 0.97,
    ha: str = "left",
    va: str = "top",
    **kwargs: Any,
) -> Any:
    """Add a small thesis-style annotation inside an axis."""
    bbox = {
        "boxstyle": "round,pad=0.25",
        "facecolor": "white",
        "edgecolor": TFG_COLORS["neutral"],
        "alpha": 0.92,
        "linewidth": 0.6,
    }
    bbox.update(kwargs.pop("bbox", {}))
    return ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        ha=ha,
        va=va,
        fontsize=9.5,
        color=TFG_COLORS["text"],
        bbox=bbox,
        **kwargs,
    )


def save_figure(fig: Any, path: str | Path, dpi: int = 300, **kwargs: Any) -> None:
    """Save a figure with consistent export settings."""
    fig.savefig(
        path,
        dpi=dpi,
        bbox_inches=kwargs.pop("bbox_inches", "tight"),
        facecolor=kwargs.pop("facecolor", TFG_COLORS["background"]),
        **kwargs,
    )
