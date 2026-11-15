#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
# Stamped copy of the house figure style, version 2026.10.17, made on 2026-11-14 for code that ships
# with a paper. It is left unchanged here so that the paper's figures stay reproducible.
"""House figure style of Dr Dutykh's papers, reports, talks and messages (matplotlib).

Palette, in order of use: navy (the main result), wine (the comparison), teal, ochre, plum,
lilac; neutrals ink, slate, mist, sand, paper. Colour never carries information alone:
series also differ by dash pattern and, when marked, by hollow markers. save() refuses a
figure in which a label, a legend or a colour bar crosses a curve, a field, another label
or a frame.

Usage:
    import manuscript_style as hs
    hs.configure(family="libertine")          # or "latin-modern" for LaTeX articles
    fig, ax = plt.subplots(figsize=(hs.WIDTH, 2.6))
    ax.plot(x, y, **hs.series(0), label="...")
    hs.save(fig, "figures/name.pdf")

Settings read from the environment (the former NLS_ names still work): FIGURE_USETEX=0 avoids
LaTeX, FIGURE_FORCE=1 writes an overlapping figure for inspection instead of refusing it,
FIGURE_PREVIEW_DIR=folder also writes a PNG preview of every figure saved.
"""
from __future__ import annotations

import glob
import os
import shutil
import subprocess
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path.home() / ".cache" / "house-style-matplotlib"))

VERSION = "2026.10.17"


def _env(name: str, default: str | None = None) -> str | None:
    """A setting FIGURE_<name>, or its former name NLS_FIGURE_<name>."""
    return os.environ.get("FIGURE_" + name, os.environ.get("NLS_FIGURE_" + name, default))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from cycler import cycler
from matplotlib import colors as mcolors
from matplotlib import font_manager
from matplotlib.collections import (LineCollection, PathCollection, PolyCollection,
                                    QuadMesh, TriMesh)
from matplotlib.text import Text
from matplotlib.transforms import Bbox, ScaledTranslation

# ---------------------------------------------------------------- palette
NAVY = "#17365D"
WINE = "#8E2A3F"
TEAL = "#007F86"
OCHRE = "#B66A19"
PLUM = "#6E4B9E"
LILAC = "#9C86C8"
INK = "#1C1C1C"
SLATE = "#5C6670"
MIST = "#C7CDD4"
SAND = "#F1EBDD"
PAPER = "#FBF9F4"
SERIES = (NAVY, WINE, TEAL, OCHRE, PLUM, LILAC)
NAMES = {"navy": NAVY, "wine": WINE, "teal": TEAL, "ochre": OCHRE, "plum": PLUM,
         "lilac": LILAC, "ink": INK, "slate": SLATE, "mist": MIST, "sand": SAND,
         "paper": PAPER}
# Dash patterns and hollow markers accompany the series colours in this order.
DASHES = ("-", (0, (5.0, 2.0)), (0, (1.6, 1.6)), (0, (5.0, 1.5, 1.5, 1.5)),
          (0, (3.0, 1.2)), (0, (1.0, 1.2)))
MARKERS = ("o", "s", "D", "^", "v", "p")
# The former role names remain usable by every generator.
BLUE, ORANGE, GREEN, PURPLE, ROSE, GOLD = NAVY, WINE, TEAL, PLUM, LILAC, OCHRE

WIDTH = 6.05   # inches; the text width of the articles
HALF = 2.95    # inches; a half-width figure
# Full and half widths of the text block, in inches, for the page formats in use.
PAGES = {"article": (6.05, 2.95), "a4": (5.748, 2.80), "springer": (4.606, 2.24)}
DENSITY_COLORMAP = "cavity"
SIGNED_COLORMAP = "cavity_diverging"
SECOND_COLORMAP = "ember"
PREAMBLE = r"\usepackage{lmodern}\usepackage[T1]{fontenc}\usepackage{amsmath,amssymb}"
PREAMBLE_LIBERTINE = (r"\usepackage[T1]{fontenc}\usepackage{libertine}"
                      r"\usepackage[libertine]{newtxmath}\usepackage{amsmath}")
METADATA = {"Creator": "House figure style", "CreationDate": None, "ModDate": None}
FIELD_THRESHOLD = 0.02   # normalized field value above which a region counts as occupied


def _register_colormaps() -> None:
    stops = {
        "cavity": [PAPER, "#DCE7E3", "#9CC7C0", "#4FA3A6", "#1E6F86", NAVY, "#0B1A33"],
        "ember": [PAPER, "#F2DFB8", "#DDA153", OCHRE, "#7A2E2E", "#3B1F3A"],
        "cavity_diverging": [WINE, PAPER, NAVY],
    }
    for name, colours in stops.items():
        if name not in matplotlib.colormaps:
            cmap = mcolors.LinearSegmentedColormap.from_list(name, colours)
            matplotlib.colormaps.register(cmap)
            matplotlib.colormaps.register(cmap.reversed())


_register_colormaps()


# ---------------------------------------------------------------- type and layout
def _latex_available() -> bool:
    if _env("USETEX", "1") == "0":
        return False
    if shutil.which("latex") is None or shutil.which("dvipng") is None:
        return False
    try:
        found = subprocess.run(["kpsewhich", "lmodern.sty"], capture_output=True, text=True)
    except OSError:
        return False
    return found.returncode == 0 and bool(found.stdout.strip())


def _register_latin_modern() -> None:
    for pattern in ("/usr/share/texmf/fonts/opentype/public/lm/lmroman*.otf",
                    "/usr/share/texlive/texmf-dist/fonts/opentype/public/lm/lmroman*.otf"):
        for path in glob.glob(pattern):
            try:
                font_manager.fontManager.addfont(path)
            except (OSError, RuntimeError):
                continue


def configure(base: float = 8.5, usetex: bool | None = None,
              family: str = "latin-modern") -> None:
    """Apply the house style. With family="latin-modern" text is set in Latin Modern
    through LaTeX when it is available, and through the matching Computer Modern fonts
    otherwise. With family="libertine" it is set in Linux Libertine (LaTeX packages
    libertine and newtxmath, or the Linux Libertine O system font without LaTeX)."""
    if family not in ("latin-modern", "libertine"):
        raise ValueError(f"unknown family {family!r}")
    if usetex is None:
        usetex = _latex_available()
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": base,
            "axes.titlesize": base + 0.5,
            "axes.titleweight": "normal",
            "axes.titlelocation": "left",
            "axes.labelsize": base,
            "legend.fontsize": base - 1.0,
            "legend.title_fontsize": base - 1.0,
            "xtick.labelsize": base - 1.0,
            "ytick.labelsize": base - 1.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.7,
            "axes.edgecolor": INK,
            "axes.labelcolor": INK,
            "axes.titlecolor": INK,
            "text.color": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.minor.width": 0.5,
            "ytick.minor.width": 0.5,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "xtick.minor.size": 1.8,
            "ytick.minor.size": 1.8,
            "axes.prop_cycle": cycler(color=list(SERIES)),
            "lines.linewidth": 1.3,
            "lines.markersize": 4.0,
            "lines.markeredgewidth": 0.9,
            "axes.grid": False,
            "axes.axisbelow": True,
            "legend.frameon": False,
            "legend.handlelength": 2.4,
            "legend.handletextpad": 0.6,
            "legend.borderaxespad": 0.4,
            "legend.columnspacing": 1.4,
            "legend.labelspacing": 0.35,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 400,
            "figure.dpi": 100,
            "figure.constrained_layout.use": True,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "image.cmap": DENSITY_COLORMAP,
            "axes.unicode_minus": False,
            "axes.formatter.use_mathtext": True,
        }
    )
    if usetex:
        preamble = PREAMBLE_LIBERTINE if family == "libertine" else PREAMBLE
        plt.rcParams.update({"text.usetex": True, "text.latex.preamble": preamble})
    elif family == "libertine":
        plt.rcParams.update(
            {
                "text.usetex": False,
                "font.serif": ["Linux Libertine O", "Libertinus Serif", "DejaVu Serif"],
                "mathtext.fontset": "custom",
                "mathtext.rm": "Linux Libertine O",
                "mathtext.it": "Linux Libertine O:italic",
                "mathtext.bf": "Linux Libertine O:bold",
                "mathtext.fallback": "stix",
            }
        )
    else:
        _register_latin_modern()
        plt.rcParams.update(
            {
                "text.usetex": False,
                "font.serif": ["Latin Modern Roman", "DejaVu Serif"],
                "mathtext.fontset": "cm",
                "mathtext.fallback": "stix",
            }
        )


def grid(axis) -> None:
    """Faint decade guides on a logarithmic ordinate; nothing on a linear one."""
    if axis.get_yscale() == "log":
        axis.grid(True, axis="y", which="major", color=MIST, linewidth=0.5)
        axis.set_axisbelow(True)


def panel(axis, label: str, pad: float = 4.0) -> None:
    axis.set_title(label, loc="left", pad=pad)


def series(index: int) -> dict:
    """Colour, dash pattern and hollow marker of the series with this index."""
    k = index % len(SERIES)
    return {"color": SERIES[k], "linestyle": DASHES[k], "marker": MARKERS[k],
            "markerfacecolor": "white", "markeredgecolor": SERIES[k]}


def legend_outside(target, where: str = "above", ncol: int | None = None,
                   handles=None, labels=None, **keywords):
    """A legend that lies outside the data area: above the axes (beside the
    panel label, which is then raised), to their right, below their abscissa
    label, or, for a figure, outside all panels ("upper", "lower", "right")."""
    figure = target if hasattr(target, "legends") and not hasattr(target, "plot") else target.figure
    if handles is None or labels is None:
        if hasattr(target, "get_legend_handles_labels"):
            handles, labels = target.get_legend_handles_labels()
        else:
            handles, labels = [], []
            for axis in figure.axes:
                got = axis.get_legend_handles_labels()
                for handle, label in zip(*got):
                    if label not in labels:
                        handles.append(handle)
                        labels.append(label)
    if ncol is None:
        ncol = min(len(labels), 3) or 1
    keywords.setdefault("frameon", False)
    if not hasattr(target, "plot"):
        location = {"upper": "outside upper center", "above": "outside upper center",
                    "lower": "outside lower center", "below": "outside lower center",
                    "right": "outside right upper"}[where]
        return figure.legend(handles, labels, loc=location, ncol=ncol, **keywords)
    axis = target
    if where == "right":
        return axis.legend(handles, labels, loc="upper left", bbox_to_anchor=(1.02, 1.0),
                           ncol=1 if ncol == min(len(labels), 3) else ncol,
                           borderaxespad=0.0, **keywords)
    if where == "below":
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        tight = axis.xaxis.get_tightbbox(renderer)
        bottom = axis.get_window_extent(renderer).y0
        drop = 0.0 if tight is None else (bottom - tight.y0) / figure.dpi
        transform = axis.transAxes + ScaledTranslation(0.0, -drop - 5.0 / 72.0,
                                                       figure.dpi_scale_trans)
        return axis.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.0),
                           bbox_transform=transform, ncol=ncol, borderaxespad=0.0, **keywords)
    if where != "above":
        raise ValueError(f"unknown legend position {where!r}")
    title = axis.get_title(loc="left")
    anchor = dict(loc="lower center", bbox_to_anchor=(0.5, 1.0)) if title else dict(loc="lower right", bbox_to_anchor=(1.0, 1.0))
    # The topmost tick label straddles the upper edge of the axes; leave room for it.
    legend = axis.legend(handles, labels, ncol=ncol, borderaxespad=0.8, **anchor, **keywords)
    if title:
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        height = legend.get_window_extent(renderer).height / figure.dpi * 72.0
        offset = 0.8 * plt.rcParams["legend.fontsize"]
        axis.set_title(title, loc="left", pad=height + offset + 5.0)
    return legend


def time_stamp(axis, time: float, occupancy=None, label: str | None = None,
               corner: str | None = None, fontsize: float | None = None):
    """A time stamp in a paper-coloured box, in the panel corner that carries
    the least field (occupancy(x0, x1, y0, y1) -> float, in data coordinates)."""
    corners = {"upper right": (0.99, 0.97, "right", "top"),
               "upper left": (0.01, 0.97, "left", "top"),
               "lower right": (0.99, 0.03, "right", "bottom"),
               "lower left": (0.01, 0.03, "left", "bottom")}
    if corner is None:
        corner = "upper right"
        if occupancy is not None:
            x0, x1 = axis.get_xlim()
            y0, y1 = axis.get_ylim()
            dx, dy = 0.14 * (x1 - x0), 0.34 * (y1 - y0)
            windows = {"upper right": (x1 - dx, x1, y1 - dy, y1),
                       "upper left": (x0, x0 + dx, y1 - dy, y1),
                       "lower right": (x1 - dx, x1, y0, y0 + dy),
                       "lower left": (x0, x0 + dx, y0, y0 + dy)}
            corner = min(windows, key=lambda name: occupancy(*windows[name]))
    x, y, ha, va = corners[corner]
    text = label if label is not None else rf"$t={time:g}$"
    return axis.text(x, y, text, transform=axis.transAxes, ha=ha, va=va,
                     fontsize=fontsize or plt.rcParams["legend.fontsize"], color=INK,
                     bbox=dict(boxstyle="round,pad=0.18", facecolor=PAPER, edgecolor="none"),
                     zorder=6)


def occupancy_of(points: np.ndarray, values: np.ndarray, normalization):
    """The share of samples inside a window whose normalized value exceeds the
    field threshold; used to place time stamps away from the wave."""
    scaled = np.asarray(normalization(np.asarray(values, float)))

    def occupancy(x0, x1, y0, y1):
        inside = ((points[:, 0] >= x0) & (points[:, 0] <= x1)
                  & (points[:, 1] >= y0) & (points[:, 1] <= y1))
        if not inside.any():
            return 0.0
        return float(np.mean(scaled[inside] > FIELD_THRESHOLD))

    return occupancy


def colorbar(figure, image, axes, label: str, **keywords):
    keywords.setdefault("pad", 0.02)
    keywords.setdefault("aspect", 30)
    keywords.setdefault("shrink", 0.9)
    bar = figure.colorbar(image, ax=axes, label=label, **keywords)
    bar.outline.set_linewidth(0.5)
    bar.outline.set_edgecolor(INK)
    bar.ax.tick_params(width=0.5, length=2.5)
    return bar


# ---------------------------------------------------------------- fields
def triangles_from_quads(quads: np.ndarray) -> np.ndarray:
    """Split each bilinear cell into two triangles for Gouraud display."""
    return np.vstack([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]])


def density_panel(axis, points, quads, density, normalization, x_range,
                  cmap: str = DENSITY_COLORMAP, y_range=None):
    """Draw one field on a quadrilateral mesh, rasterized inside an otherwise
    vector figure, with the house colour map."""
    triangles = triangles_from_quads(quads)
    centres = np.mean(points[triangles, 0], axis=1)
    triangulation = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
    triangulation.set_mask(
        (centres < x_range[0] - 1.0e-7) | (centres > x_range[1] + 1.0e-7)
    )
    image = axis.tripcolor(triangulation, density, shading="gouraud", cmap=cmap,
                           norm=normalization, rasterized=True)
    axis.set_xlim(*x_range)
    if y_range is not None:
        axis.set_ylim(*y_range)
    axis.set_aspect("equal")
    return image


def field_panel(axis, points, triangles, values, normalization, x_range=None,
                cmap: str = DENSITY_COLORMAP, y_range=None):
    """Draw one field on a triangulation (points (n, 2), triangles (m, 3))."""
    triangulation = mtri.Triangulation(points[:, 0], points[:, 1], triangles)
    if x_range is not None:
        centres = np.mean(points[triangles, 0], axis=1)
        triangulation.set_mask((centres < x_range[0] - 1.0e-7) | (centres > x_range[1] + 1.0e-7))
    image = axis.tripcolor(triangulation, values, shading="gouraud", cmap=cmap,
                           norm=normalization, rasterized=True)
    if x_range is not None:
        axis.set_xlim(*x_range)
    if y_range is not None:
        axis.set_ylim(*y_range)
    axis.set_aspect("equal")
    return image


# ---------------------------------------------------------------- checks
def _padded(box: Bbox, pad: float) -> Bbox:
    return Bbox.from_extents(box.x0 + pad, box.y0 + pad, box.x1 - pad, box.y1 - pad)


def _visible_text(artist) -> bool:
    return artist is not None and artist.get_visible() and bool(artist.get_text().strip())


def _label_boxes(figure, renderer, pad):
    boxes = []

    def text(artist, name):
        if not _visible_text(artist) or _background(artist):
            return
        # An annotation's own extent would include its leader arrow; the label is the text.
        box = Text.get_window_extent(artist, renderer)
        if box.width > 0 and box.height > 0:
            boxes.append((name, _padded(box, pad)))

    def legend(artist, name):
        if artist is None or not artist.get_visible():
            return
        boxes.append((name, _padded(artist.get_window_extent(renderer), pad)))

    for k, axis in enumerate(figure.axes):
        tag = f"axes {k}"
        legend(axis.get_legend(), f"{tag}: legend")
        for artist in axis.texts:
            text(artist, f"{tag}: text {artist.get_text()[:40]!r}")
        for artist in (axis.title, axis._left_title, axis._right_title):
            text(artist, f"{tag}: title {artist.get_text()[:40]!r}")
        text(axis.xaxis.label, f"{tag}: x label")
        text(axis.yaxis.label, f"{tag}: y label")
        for axis_object, name in ((axis.xaxis, "x"), (axis.yaxis, "y")):
            # compare tick positions in the scaled coordinate of the axis, so that a
            # logarithmic axis is treated like a linear one
            scale = axis_object.get_transform().transform
            low, high = sorted(scale(np.asarray(axis_object.get_view_interval(), float)))
            slack = 1e-9 * max(high - low, 1e-300)
            for tick in axis_object.get_major_ticks() + axis_object.get_minor_ticks():
                if not tick.get_visible():
                    continue
                position = float(scale(np.asarray([tick.get_loc()], float))[0])
                if not (low - slack <= position <= high + slack):
                    continue
                for label in (tick.label1, tick.label2):
                    text(label, f"{tag}: {name} tick label {label.get_text()!r}")
    for artist in figure.legends:
        legend(artist, "figure legend")
    for artist in figure.texts:
        text(artist, f"figure text {artist.get_text()[:40]!r}")
    return boxes


def _background(artist) -> bool:
    gid = artist.get_gid()
    return bool(gid) and gid.startswith("background")


def _field_points(collection, axis):
    """Display-space sample points and normalized values of a field artist."""
    values = collection.get_array()
    if values is None:
        return None, None
    values = np.asarray(values, float).ravel()
    if isinstance(collection, TriMesh):
        triangulation = getattr(collection, "_triangulation", None)
        if triangulation is None:
            return None, None
        xy = np.column_stack([triangulation.x, triangulation.y])
        if len(values) != len(xy):
            xy = xy[triangulation.triangles].mean(axis=1)
    elif isinstance(collection, QuadMesh):
        coordinates = np.asarray(collection.get_coordinates(), float)
        if values.size == coordinates.shape[0] * coordinates.shape[1]:
            xy = coordinates.reshape(-1, 2)
        else:
            centres = 0.25 * (coordinates[:-1, :-1] + coordinates[1:, :-1]
                              + coordinates[:-1, 1:] + coordinates[1:, 1:])
            xy = centres.reshape(-1, 2)
    elif isinstance(collection, PolyCollection):
        xy = np.array([path.vertices[:-1].mean(axis=0) if len(path.vertices) > 1
                       else path.vertices[0] for path in collection.get_paths()])
    else:
        return None, None
    if len(xy) != len(values):
        return None, None
    display = collection.get_transform().transform(xy)
    scaled = np.asarray(collection.norm(values), float)
    return display, scaled


def _elements(figure, renderer):
    """Every drawn element with the box to which it is clipped (its axes)."""
    items = []
    for k, axis in enumerate(figure.axes):
        tag = f"axes {k}"
        clip = axis.get_window_extent(renderer)
        for line in axis.lines:
            if not line.get_visible() or _background(line):
                continue
            path = line.get_transform().transform_path(line.get_path())
            vertices = path.vertices[np.isfinite(path.vertices).all(axis=1)]
            items.append((f"{tag}: line {line.get_label()!r}", "path", (path, vertices), clip if line.get_clip_on() else None))
        for collection in axis.collections:
            if not collection.get_visible() or _background(collection):
                continue
            if isinstance(collection, (TriMesh, QuadMesh)) or (
                    isinstance(collection, PolyCollection) and collection.get_array() is not None):
                items.append((f"{tag}: field", "field", collection, clip))
            elif isinstance(collection, PathCollection):
                offsets = collection.get_offset_transform().transform(collection.get_offsets())
                items.append((f"{tag}: markers", "points", np.asarray(offsets, float), clip if collection.get_clip_on() else None))
            else:
                transform = collection.get_transform()
                for path in collection.get_paths():
                    shown = transform.transform_path(path)
                    items.append((f"{tag}: collection", "path", (shown, shown.vertices), clip if collection.get_clip_on() else None))
        for patch in axis.patches:
            if patch is axis.patch or not patch.get_visible() or _background(patch):
                continue
            shown = patch.get_transform().transform_path(patch.get_path())
            items.append((f"{tag}: patch {patch.get_label()!r}", "path", (shown, shown.vertices), clip if patch.get_clip_on() else None))
        for image in axis.images:
            if image.get_visible():
                items.append((f"{tag}: image", "image", (image, axis), clip))
        for name, spine in axis.spines.items():
            if spine.get_visible():
                shown = spine.get_transform().transform_path(spine.get_path())
                items.append((f"{tag}: frame {name}", "path", (shown, shown.vertices), None))
    return items


def _hits_path(box, path, vertices) -> bool:
    if len(vertices) and np.any((vertices[:, 0] >= box.x0) & (vertices[:, 0] <= box.x1)
                                & (vertices[:, 1] >= box.y0) & (vertices[:, 1] <= box.y1)):
        return True
    try:
        return bool(path.intersects_bbox(box, filled=False))
    except (ValueError, TypeError):
        return False


def _hits_field(box, collection) -> bool:
    display, scaled = _field_points(collection, None)
    if display is None:
        return False
    inside = ((display[:, 0] >= box.x0) & (display[:, 0] <= box.x1)
              & (display[:, 1] >= box.y0) & (display[:, 1] <= box.y1))
    if not inside.any():
        return False
    return bool(np.nanmax(np.abs(scaled[inside] - _field_zero(collection))) > FIELD_THRESHOLD)


def _field_zero(collection) -> float:
    """Normalized value of an empty region: zero for a density map, the centre
    for a signed map."""
    norm = collection.norm
    if isinstance(norm, (mcolors.TwoSlopeNorm, mcolors.CenteredNorm)):
        return 0.5
    try:
        return float(norm(0.0))
    except Exception:
        return 0.0


def _hits_image(box, image, axis) -> bool:
    array = image.get_array()
    if array is None:
        return False
    array = np.ma.filled(np.asarray(array, float), np.nan)
    if array.ndim == 3:
        return True
    x0, x1, y0, y1 = image.get_extent()
    corners = axis.transData.transform([[x0, y0], [x1, y1], [x0, y1], [x1, y0]])
    extent = Bbox.from_extents(corners[:, 0].min(), corners[:, 1].min(),
                               corners[:, 0].max(), corners[:, 1].max())
    if not box.overlaps(extent):
        return False
    inverse = axis.transData.inverted()
    (dx0, dy0), (dx1, dy1) = inverse.transform([[box.x0, box.y0], [box.x1, box.y1]])
    rows, cols = array.shape
    def index(value, low, high, count):
        return int(np.clip(np.floor((value - low) / (high - low) * count), 0, count - 1))
    c0, c1 = sorted((index(dx0, x0, x1, cols), index(dx1, x0, x1, cols)))
    origin_low = image.origin == "lower"
    ys = (y0, y1) if origin_low else (y1, y0)
    r0, r1 = sorted((index(dy0, *ys, rows), index(dy1, *ys, rows)))
    sub = array[r0:r1 + 1, c0:c1 + 1]
    if sub.size == 0 or np.all(np.isnan(sub)):
        return False
    zero = _field_zero(image)
    return bool(np.nanmax(np.abs(np.asarray(image.norm(sub), float) - zero)) > FIELD_THRESHOLD)


def check_overlaps(figure, pad_points: float = 0.5) -> list[str]:
    """Every collision between a label box (legend, text, annotation, title,
    axis label, tick label) and any other element (curve, marker, patch,
    field, image, frame, other label) in display space."""
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    pad = pad_points * figure.dpi / 72.0
    boxes = _label_boxes(figure, renderer, pad)
    elements = _elements(figure, renderer)
    problems = []
    for name, full_box in boxes:
        for element_name, kind, payload, clip in elements:
            box = full_box if clip is None else Bbox.intersection(full_box, clip)
            if box is None or box.width <= 0 or box.height <= 0:
                continue
            if kind == "path" and _hits_path(box, *payload):
                problems.append(f"{name} crosses {element_name}")
            elif kind == "points":
                inside = ((payload[:, 0] >= box.x0) & (payload[:, 0] <= box.x1)
                          & (payload[:, 1] >= box.y0) & (payload[:, 1] <= box.y1))
                if inside.any():
                    problems.append(f"{name} crosses {element_name}")
            elif kind == "field" and _hits_field(box, payload):
                problems.append(f"{name} lies on {element_name}")
            elif kind == "image" and _hits_image(box, *payload):
                problems.append(f"{name} lies on {element_name}")
    for i, (name, box) in enumerate(boxes):
        for other_name, other in boxes[i + 1:]:
            if box.overlaps(other):
                problems.append(f"{name} overlaps {other_name}")
    return problems


def _check_text(figure) -> None:
    if not plt.rcParams.get("text.usetex", False):
        return
    bad = []
    artists = []
    for axis in figure.axes:
        artists += list(axis.texts) + [axis.title, axis._left_title, axis._right_title,
                                       axis.xaxis.label, axis.yaxis.label]
        legend = axis.get_legend()
        if legend is not None:
            artists.extend(legend.get_texts())
    for artist in figure.texts:
        artists.append(artist)
    for legend in figure.legends:
        artists.extend(legend.get_texts())
    for artist in artists:
        string = artist.get_text()
        if any(ord(char) > 127 for char in string):
            bad.append(string)
    if bad:
        raise ValueError("non-ASCII characters are not typeset by LaTeX in figure text: " + repr(bad))


def save(figure, path, allow_overlap: str | None = None) -> Path:
    """Check the figure, write it deterministically, and close it. A
    collision aborts the save unless `allow_overlap` states the reason."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    _check_text(figure)
    problems = check_overlaps(figure)
    if problems and not allow_overlap:
        if _env("FORCE") == "1":
            # Inspection mode: write the figure and report, instead of refusing it.
            print(f"{path.name}: figure elements overlap (written for inspection):\n  " + "\n  ".join(problems))
        else:
            raise RuntimeError(f"{path.name}: figure elements overlap:\n  " + "\n  ".join(problems))
    figure.savefig(path, metadata=METADATA)
    preview = _env("PREVIEW_DIR")
    if preview:
        Path(preview).mkdir(parents=True, exist_ok=True)
        figure.savefig(Path(preview) / (path.stem + ".png"), dpi=150)
    plt.close(figure)
    return path
