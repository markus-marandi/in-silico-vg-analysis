"""
Thesis-wide matplotlib / seaborn style standards.

Usage in any notebook (add near the top, after imports):
    import sys; sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[1]))
    from thesis_style import apply_thesis_style, sig_stars, save_figure, add_panel_letter

Or from the 10_THESIS/ directory:
    import sys, pathlib
    sys.path.insert(0, str(pathlib.Path.cwd().parent))   # project root
    from ten_THESIS.thesis_style import ...

The simplest import from inside a notebook:
    import sys; sys.path.insert(0, '..'); from ten_THESIS import thesis_style; thesis_style.apply()
But the cleanest approach is to resolve the project root (as done in every existing notebook)
and then:
    from thesis_style import apply_thesis_style, sig_stars, save_figure, add_panel_letter
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# ── Paths ─────────────────────────────────────────────────────────────────────
_HERE = Path(__file__).resolve().parent          # 10_THESIS/
PROJECT_ROOT = _HERE.parent                       # in-silico-vg-analysis/
FIGURE_DIR = _HERE / "FIGURES_FOR_THESIS"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

# ── Font stack ────────────────────────────────────────────────────────────────
_FONT_FAMILY = ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"]

# ── Size constants ─────────────────────────────────────────────────────────────
TITLE_SIZE      = 14
LABEL_SIZE      = 12
TICK_SIZE       = 10
LEGEND_SIZE     = 9
ANNOT_SIZE      = 9
LINEWIDTH_BOX   = 1.2   # box/violin outlines
LINEWIDTH_GRID  = 0.8

# ── Save settings ─────────────────────────────────────────────────────────────
SAVE_DPI   = 300         # rasterised export
INLINE_DPI = 150         # Jupyter inline


def apply_thesis_style() -> None:
    """Apply the thesis-wide rcParams and seaborn theme.

    Call once per notebook, after importing matplotlib / seaborn.
    """
    sns.set_theme(
        style="ticks",          # clean base; we despine per-axis
        context="paper",
        font="sans-serif",
        rc={
            # fonts
            "font.family":           "sans-serif",
            "font.sans-serif":       _FONT_FAMILY,
            # sizes
            "axes.titlesize":        TITLE_SIZE,
            "axes.labelsize":        LABEL_SIZE,
            "xtick.labelsize":       TICK_SIZE,
            "ytick.labelsize":       TICK_SIZE,
            "legend.fontsize":       LEGEND_SIZE,
            # lines / patches
            "lines.linewidth":       1.4,
            "patch.linewidth":       LINEWIDTH_BOX,
            "axes.linewidth":        0.8,
            "grid.linewidth":        LINEWIDTH_GRID,
            # grid
            "axes.grid":             True,
            "grid.alpha":            0.35,
            "grid.color":            "#c0c0c0",
            # legend
            "legend.frameon":        False,
            # figure
            "figure.dpi":            INLINE_DPI,
            "savefig.dpi":           SAVE_DPI,
            "figure.constrained_layout.use": True,
            # pdf / font embedding
            "pdf.fonttype":          42,
            "ps.fonttype":           42,
            "svg.fonttype":          "none",
        },
    )


# Call at import time so a bare `import thesis_style` is enough.
apply_thesis_style()


# ── Significance stars ────────────────────────────────────────────────────────

def sig_stars(p: float) -> str:
    """Return significance label including the p-value, e.g. '** (p=3.2e-04)'."""
    if p < 0.001:
        stars = "***"
    elif p < 0.01:
        stars = "**"
    elif p < 0.05:
        stars = "*"
    else:
        stars = "ns"
    return f"{stars} (p={p:.2e})"


def sig_label(p: float) -> str:
    """Return just the star symbol (or 'ns'), without the p-value."""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


# ── Figure saving ─────────────────────────────────────────────────────────────

def save_figure(
    fig: plt.Figure | None = None,
    name: str = "figure",
    outdir: Path | str | None = None,
    formats: Sequence[str] = ("png", "pdf"),
    dpi: int = SAVE_DPI,
    tight: bool = True,
) -> list[Path]:
    """Save *fig* as PNG (raster) and PDF (vector) into *outdir*.

    Returns the list of paths written.
    """
    fig = fig or plt.gcf()
    outdir = Path(outdir) if outdir else FIGURE_DIR
    outdir.mkdir(parents=True, exist_ok=True)

    if tight:
        try:
            fig.tight_layout()
        except Exception:
            pass

    written: list[Path] = []
    for fmt in formats:
        dest = outdir / f"{name}.{fmt}"
        kwargs: dict = {"bbox_inches": "tight", "facecolor": "white"}
        if fmt.lower() in ("png", "jpg", "tiff"):
            kwargs["dpi"] = dpi
        fig.savefig(dest, format=fmt, **kwargs)
        written.append(dest)
        print(f"[thesis_style] saved {dest.name}")

    return written


# ── Panel letters ─────────────────────────────────────────────────────────────

def add_panel_letter(
    ax: plt.Axes,
    letter: str,
    x: float = -0.12,
    y: float = 1.04,
    size: int = TITLE_SIZE,
) -> None:
    """Add a bold panel letter (A, B, C …) in the top-left of *ax*."""
    ax.text(
        x, y, letter,
        transform=ax.transAxes,
        fontsize=size,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def add_panel_letters(
    axes: Iterable[plt.Axes],
    letters: Iterable[str] | None = None,
    **kwargs,
) -> None:
    """Add panel letters to a sequence of axes.

    *letters* defaults to A, B, C, … if not supplied.
    """
    import string
    axes = list(axes)
    letters = list(letters) if letters else list(string.ascii_uppercase[: len(axes)])
    for ax, letter in zip(axes, letters):
        add_panel_letter(ax, letter, **kwargs)


# ── Despine helper ────────────────────────────────────────────────────────────

def despine_all(axes: Iterable[plt.Axes], **kwargs) -> None:
    """Call sns.despine on every axis in *axes*."""
    for ax in axes:
        sns.despine(ax=ax, **kwargs)


# ── Legend helpers ────────────────────────────────────────────────────────────

def outside_legend(ax: plt.Axes, **kwargs) -> None:
    """Place legend to the right of the axes, no frame."""
    defaults = dict(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        fontsize=LEGEND_SIZE,
        borderaxespad=0,
    )
    defaults.update(kwargs)
    ax.legend(**defaults)


# ── Annotation helpers ────────────────────────────────────────────────────────

def annotate_pval(
    ax: plt.Axes,
    p: float,
    x0: float,
    x1: float,
    y: float,
    dy_bar: float = 0.02,
    fontsize: int = ANNOT_SIZE,
    color: str = "black",
) -> None:
    """Draw a significance bracket between x0 and x1 at height y (data coords)."""
    y_bar = y + dy_bar
    ax.plot([x0, x0, x1, x1], [y, y_bar, y_bar, y], lw=0.9, c=color)
    ax.text(
        (x0 + x1) / 2, y_bar,
        sig_stars(p),
        ha="center", va="bottom",
        fontsize=fontsize,
        color=color,
    )


# ── Minor gridlines for log-scale axes ───────────────────────────────────────

def enable_log_minor_grid(ax: plt.Axes) -> None:
    """Show minor gridlines on a log-scale axis (both axes if log/log)."""
    ax.minorticks_on()
    ax.grid(which="minor", linestyle=":", linewidth=0.4, alpha=0.4, color="#c0c0c0")
    ax.grid(which="major", linestyle="-", linewidth=LINEWIDTH_GRID, alpha=0.45, color="#c0c0c0")
