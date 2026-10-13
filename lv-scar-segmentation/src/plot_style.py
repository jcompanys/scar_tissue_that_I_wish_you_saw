"""Shared plotting style for the LV scar characterization notebooks.

The palette is derived from the M01/M02 notebook style and kept deliberately
formal for thesis figures: restrained clinical colors, white backgrounds,
light grids, and export-ready DPI defaults.

The module exposes three kinds of object:

* colour tables (``TFG_COLORS`` and the task-specific maps built from it) so
  every notebook draws the same category with the same hex value;
* colormap names (``SCAR_CMAP``, ``DIVERGING_CMAP``) for image/heatmap plots;
* style helpers (``set_tfg_style`` once per notebook, then ``apply_tfg_style``/
  ``annotate_panel``/``save_figure`` per axis or figure).

Nothing here reads data or produces figures on import; it only configures
matplotlib/seaborn and returns the axes it is handed.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt

# Master palette. Keys are semantic roles (not raw colour names) so a notebook
# asks for "primary"/"missing" rather than a hex string, and every downstream
# map below is built from these entries to stay consistent.
TFG_COLORS = {
    "primary": "#256D85",    # teal: main series / segmented group
    "secondary": "#D9822B",  # amber: second series / medians / highlights
    "positive": "#2E8B57",   # green: agreement / "good" category
    "neutral": "#B8C0C8",    # grey: de-emphasised / "unknown"
    "missing": "#C44E52",    # red: missing data / warning category
    "male": "#256D85",       # sex coding reuses primary (teal)
    "female": "#B35C9E",     # sex coding: magenta
    "line": "#263238",       # near-black for axis spines
    "text": "#263238",       # near-black for labels/ticks/titles
    "background": "#FFFFFF",  # white figure/axes background for print
    "grid": "#B8C0C8",       # light grey gridlines
}

# Ordered categorical cycle handed to seaborn as its default palette; the order
# controls which colour the Nth plotted series receives.
TFG_PALETTE = [
    TFG_COLORS["primary"],
    TFG_COLORS["secondary"],
    TFG_COLORS["positive"],
    TFG_COLORS["female"],
    TFG_COLORS["missing"],
    TFG_COLORS["neutral"],
]

# Semantic aliases used by the EDA notebook so call sites read by meaning.
HIST_COLOR = TFG_COLORS["primary"]       # histogram / scatter fill
MEDIAN_COLOR = TFG_COLORS["secondary"]   # median line / regression fit
SEG_COLOR = TFG_COLORS["secondary"]      # "segmented" patients
NO_SEG_COLOR = TFG_COLORS["primary"]     # "not segmented" patients

# Dataset sex coding: 1 -> male, 2 -> female (as stored in the clinical table).
SEX_COLORS = {1: TFG_COLORS["male"], 2: TFG_COLORS["female"]}

# Coronary artery territories, each given a distinct clinical hue for
# AHA-segment territory overlays.
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

# Colours for the three candidate polar/PCA axes compared in M03.
POLAR_CANDIDATE_COLORS = {
    "RV PC1": TFG_COLORS["primary"],
    "LV PC1": TFG_COLORS["secondary"],
    "LV PC2": "#7A68A6",  # muted purple, distinct from the two PC1 hues
}

# Low/mid/high threshold bands (e.g. retrieval score cut-offs in M05).
THRESHOLD_COLORS = {
    "low": "#C44E52",
    "mid": TFG_COLORS["primary"],
    "high": TFG_COLORS["positive"],
}

# One colour per myocardial level, apex -> base, for ring-coloured plots.
RING_COLORS = {
    "Apex": "#7A68A6",
    "Apical": "#C83E4D",
    "Mid": TFG_COLORS["positive"],
    "Basal": TFG_COLORS["primary"],
}

# Standard AHA-17 segment -> myocardial level mapping (segments 1-6 basal,
# 7-12 mid, 13-16 apical, 17 apex). Built as the default so callers can pass a
# custom mapping but normally rely on the convention.
AHA_SEGMENT_RING = {
    **{s: "Basal" for s in range(1, 7)},
    **{s: "Mid" for s in range(7, 13)},
    **{s: "Apical" for s in range(13, 17)},
    17: "Apex",
}


def color_aha_xticklabels(
    ax: Any, segments: Any, ring_map: dict[int, str] | None = None
) -> None:
    """Colour AHA-17 x-tick labels by myocardial level.

    Each tick is recoloured to its ring (basal/mid/apical/apex) and bolded, so
    a bar/line plot over AHA segments shows level membership without a separate
    axis.

    Args:
        ax: Matplotlib axis whose x-tick labels line up, in order, with
            ``segments``.
        segments: Sequence of AHA segment ids (1-17) in the same left-to-right
            order as the tick labels; each is cast to ``int`` for lookup.
        ring_map: Optional segment-id -> ring-name mapping; defaults to
            :data:`AHA_SEGMENT_RING`.

    Returns:
        None. The axis is modified in place.
    """
    ring_map = ring_map or AHA_SEGMENT_RING
    # Pair each tick with its segment id positionally; skip ids absent from the
    # map (leaves their default colour untouched).
    for tick, seg in zip(ax.get_xticklabels(), segments):
        ring = ring_map.get(int(seg))
        if ring:
            tick.set_color(RING_COLORS[ring])
            tick.set_fontweight("bold")


def add_aha_ring_legend(ax: Any, loc: str = "upper right", **kwargs: Any) -> Any:
    """Add a legend mapping ring colours to Basal/Mid/Apical/Apex.

    Args:
        ax: Matplotlib axis to attach the legend to.
        loc: Matplotlib legend location string.
        **kwargs: Extra keyword arguments forwarded to ``ax.legend``.

    Returns:
        The created :class:`matplotlib.legend.Legend`.
    """
    # Imported locally so the module has no hard dependency on patches at import
    # time (only this helper needs it).
    import matplotlib.patches as mpatches

    # Fixed base-to-apex order so the legend reads anatomically top to bottom.
    handles = [
        mpatches.Patch(color=RING_COLORS[r], label=r)
        for r in ("Basal", "Mid", "Apical", "Apex")
    ]
    return ax.legend(handles=handles, fontsize=8, title="AHA level", loc=loc, **kwargs)


# Sequential colormap for scar magnitude (light -> dark = low -> high burden).
SCAR_CMAP = "YlOrRd"
# Diverging colormap for signed/centred quantities (e.g. differences vs a mean).
DIVERGING_CMAP = "RdBu_r"


def set_tfg_style(context: str = "notebook", font_scale: float = 1.0) -> None:
    """Apply the shared TFG plotting defaults to matplotlib and seaborn.

    Call once near the top of a notebook. It sets global ``plt.rcParams`` and,
    if seaborn is installed, a matching seaborn theme so both backends draw with
    the same fonts, colours and gridlines.

    Args:
        context: ``"notebook"`` for on-screen work (slightly larger base font)
            or anything else for the tighter print base; also passed through to
            seaborn's ``set_theme`` context.
        font_scale: Multiplier applied to every font size, for scaling a figure
            up or down without editing each size individually.

    Returns:
        None.
    """
    # Notebook context gets a 1pt larger base so on-screen text stays legible;
    # print/export uses the smaller base.
    base_font = 11 if context == "notebook" else 10
    plt.rcParams.update(
        {
            "figure.dpi": 125,          # crisp on-screen rendering
            "savefig.dpi": 300,         # print-quality export
            "savefig.bbox": "tight",    # trim surrounding whitespace on save
            "savefig.facecolor": TFG_COLORS["background"],
            "figure.facecolor": TFG_COLORS["background"],
            "axes.facecolor": TFG_COLORS["background"],
            "axes.edgecolor": TFG_COLORS["line"],
            "axes.labelcolor": TFG_COLORS["text"],
            "axes.spines.top": False,   # open "L"-frame: drop top/right spines
            "axes.spines.right": False,
            "axes.grid": True,
            "axes.axisbelow": True,     # draw grid behind data, not over it
            "grid.color": TFG_COLORS["grid"],
            "grid.alpha": 0.18,         # faint grid so data dominates
            "grid.linewidth": 0.7,
            "font.family": "DejaVu Sans",  # bundled with matplotlib -> portable
            "font.size": base_font * font_scale,
            "axes.titleweight": "semibold",
            "axes.titlesize": 11 * font_scale,
            "axes.labelsize": 10 * font_scale,
            "legend.frameon": False,    # borderless legends
            "legend.fontsize": 9 * font_scale,
            "xtick.color": TFG_COLORS["text"],
            "ytick.color": TFG_COLORS["text"],
            "xtick.labelsize": 9 * font_scale,
            "ytick.labelsize": 9 * font_scale,
            "lines.linewidth": 1.6,
            "patch.linewidth": 0.8,
        }
    )

    # Seaborn is optional: mirror the same rc so seaborn-drawn plots match, but
    # degrade silently if seaborn is not installed in the environment.
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
    """Format one matplotlib axis with the shared thesis style.

    Use per-axis to enforce the style on axes created outside
    :func:`set_tfg_style` (e.g. inside subplot grids): removes the top/right
    spines, draws a single faint gridline direction, and recolours ticks,
    labels and title.

    Args:
        ax: Matplotlib axis to restyle in place.
        xgrid: If ``True`` grid the x-axis, else the y-axis. Ignored when
            ``grid_axis`` is given.
        grid_axis: Explicit grid direction (``"x"``, ``"y"`` or ``"both"``);
            overrides ``xgrid`` when set.

    Returns:
        The same ``ax``, for chaining.
    """
    # grid_axis wins if provided; otherwise map the xgrid flag to a direction.
    axis = grid_axis or ("x" if xgrid else "y")
    ax.spines[["top", "right"]].set_visible(False)
    # Match set_tfg_style's faint grid (alpha 0.18) for axes it did not touch.
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
    """Add a small thesis-style annotation inside an axis.

    Draws boxed text in axes-fraction coordinates (so position is independent
    of the data range), by default in the top-left corner.

    Args:
        ax: Target matplotlib axis.
        text: Annotation string.
        x: X position in axes fraction (0 = left, 1 = right).
        y: Y position in axes fraction (0 = bottom, 1 = top).
        ha: Horizontal alignment of the text.
        va: Vertical alignment of the text.
        **kwargs: Forwarded to ``ax.text``; a ``bbox`` dict here is merged into
            (and overrides keys of) the default rounded white box.

    Returns:
        The created :class:`matplotlib.text.Text`.
    """
    # Default rounded, semi-opaque white box; caller's bbox keys override these.
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
        transform=ax.transAxes,  # interpret x/y as axes fractions, not data
        ha=ha,
        va=va,
        fontsize=9.5,
        color=TFG_COLORS["text"],
        bbox=bbox,
        **kwargs,
    )


def save_figure(fig: Any, path: str | Path, dpi: int = 300, **kwargs: Any) -> None:
    """Save a figure with consistent export settings.

    Args:
        fig: Matplotlib figure to save.
        path: Output file path (extension selects the format).
        dpi: Export resolution; defaults to 300 for print quality.
        **kwargs: Forwarded to ``fig.savefig``. ``bbox_inches`` and
            ``facecolor`` fall back to the thesis defaults ("tight" crop, white
            background) when not supplied.

    Returns:
        None.
    """
    fig.savefig(
        path,
        dpi=dpi,
        bbox_inches=kwargs.pop("bbox_inches", "tight"),
        facecolor=kwargs.pop("facecolor", TFG_COLORS["background"]),
        **kwargs,
    )
