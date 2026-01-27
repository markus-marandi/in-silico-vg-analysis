import os
import re
import datetime
import matplotlib as mpl
import matplotlib.pyplot as plt
from contextlib import contextmanager
from pathlib import Path


UTILS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = UTILS_DIR.parent
DEFAULT_OUTDIR = PROJECT_ROOT / "figures"

DEFAULT_FORMATS = ("pdf", "svg")  # PDF is best for Affinity Designer
EMBED_FONTS = True
PDF_FONTTYPE = 42  # Makes text editable as fonts in Affinity/Illustrator
PNG_DPI = 300
RASTERIZE_POINTCOUNT_THRESHOLD = 10000

mpl.rcParams['pdf.fonttype'] = PDF_FONTTYPE
mpl.rcParams['ps.fonttype'] = PDF_FONTTYPE
mpl.rcParams['svg.fonttype'] = 'none'  # Export text as text, not paths

_fname_sanitize_re = re.compile(r'[^0-9A-Za-z\-_\.]')


def _sanitize(s: str, max_len: int = 150):
    if not s: return "unnamed"
    s = s.strip().replace(" ", "_")
    s = _fname_sanitize_re.sub("", s)
    return s[:max_len]


def _now_timestamp(fmt: str = "%d%m%Y_%H%M"):
    return datetime.datetime.now().strftime(fmt)


def _detect_notebook_name():
    """Best-effort detection of notebook or script name."""
    try:
        # Check if running in a script
        import __main__
        if hasattr(__main__, '__file__'):
            return Path(__main__.__file__).stem
    except Exception:
        pass
    return "interactive"


def _make_basename(notebook_name: str, title: str, ts: str):
    nb = _sanitize(notebook_name)
    ttl = _sanitize(title)
    return f"{nb}_{ttl}_{ts}"


# ---- Main Save Function ----
def save_plot(fig: plt.Figure = None,
              ax: plt.Axes = None,
              title: str = None,
              notebook: str = None,
              outdir: Path | str = None,
              formats: tuple = None,
              rasterize_threshold: int = None,
              force_rasterize_artists: list = None,
              verbose: bool = True):
    """
    Save the provided figure using settings optimized for Affinity Designer.
    """
    fig = fig or plt.gcf()
    ax = ax or (fig.axes[0] if fig.axes else None)

    # Path handling
    outdir = Path(outdir) if outdir else DEFAULT_OUTDIR
    outdir.mkdir(parents=True, exist_ok=True)
    
    formats = formats or DEFAULT_FORMATS
    rasterize_threshold = rasterize_threshold or RASTERIZE_POINTCOUNT_THRESHOLD

    # Metadata
    if title is None:
        if ax is not None:
            title = ax.get_title() or getattr(fig, '_suptitle_text', None) or "plot"
        else:
            title = "plot"

    notebook_name = notebook or _detect_notebook_name()
    ts = _now_timestamp()
    basename = _make_basename(notebook_name, title, ts)

    # Handle heavy data points (Vector-Raster Hybrid)
    if ax is not None and rasterize_threshold:
        for artist in list(ax.collections) + list(ax.lines):
            if hasattr(artist, "get_offsets"):
                try:
                    npts = len(artist.get_offsets())
                    if npts >= rasterize_threshold:
                        artist.set_rasterized(True)
                        if verbose:
                            # FIXED: Removed {fmt} reference here
                            print(f"[save_plot] Rasterizing {npts} points to optimize vector export.")
                except Exception:
                    continue

    if force_rasterize_artists:
        for a in force_rasterize_artists:
            try: a.set_rasterized(True)
            except: pass

    saved_paths = []
    for fmt in formats:
        fname = outdir / f"{basename}.{fmt}"
        try:
            # Affinity Designer prefers PDF for layered vector imports
            save_kwargs = {"bbox_inches": "tight", "transparent": True}
            
            if fmt.lower() in ("png", "jpg", "tiff"):
                save_kwargs["dpi"] = PNG_DPI
            
            fig.savefig(fname, format=fmt, **save_kwargs)
            saved_paths.append(str(fname))
            
            if verbose:
                print(f"[save_plot] Saved: {fname.name}")
        except Exception as e:
            print(f"[save_plot] Error saving {fmt}: {e}")

    return saved_paths


@contextmanager
def autosave(title: str, **kwargs):
    """Context manager for quick plotting blocks."""
    try:
        yield
    finally:
        save_plot(title=title, **kwargs)