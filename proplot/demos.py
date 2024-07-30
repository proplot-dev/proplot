#!/usr/bin/env python3
"""
Functions for displaying colors and fonts.
"""
import os
import re

import cycler
import matplotlib.colors as mcolors
import matplotlib.font_manager as mfonts
import numpy as np

from . import colors as pcolors
from . import constructor, ui
from .config import _get_data_folders, rc
from .internals import ic  # noqa: F401
from .internals import _not_none, _version_mpl, docstring, warnings
from .utils import to_rgb, to_xyz

__all__ = [
    "show_cmaps",
    "show_channels",
    "show_colors",
    "show_colorspaces",
    "show_cycles",
    "show_fonts",
]


# Tables and constants
FAMILY_TEXGYRE = (
    "TeX Gyre Heros",  # sans-serif
    "TeX Gyre Schola",  # serif
    "TeX Gyre Bonum",
    "TeX Gyre Termes",
    "TeX Gyre Pagella",
    "TeX Gyre Chorus",  # cursive
    "TeX Gyre Adventor",  # fantasy
    "TeX Gyre Cursor",  # monospace
)
COLOR_TABLE = {
    # NOTE: Just want the names but point to the dictionaries because
    # they don't get filled until after __init__ imports this module.
    "base": mcolors.BASE_COLORS,
    "css4": mcolors.CSS4_COLORS,
    "opencolor": pcolors.COLORS_OPEN,
    "xkcd": pcolors.COLORS_XKCD,
}
CYCLE_TABLE = {
    "Matplotlib defaults": (
        "default",
        "classic",
    ),
    "Matplotlib stylesheets": (
        # NOTE: Do not include 'solarized' because colors are terrible for
        # colorblind folks.
        "colorblind",
        "colorblind10",
        "tableau",
        "ggplot",
        "538",
        "seaborn",
        "bmh",
    ),
    "ColorBrewer2.0 qualitative": (
        "Accent",
        "Dark2",
        "Paired",
        "Pastel1",
        "Pastel2",
        "Set1",
        "Set2",
        "Set3",
        "tab10",
        "tab20",
        "tab20b",
        "tab20c",
    ),
    "Other qualitative": (
        "FlatUI",
        "Qual1",
        "Qual2",
    ),
}
CMAP_TABLE = {
    # NOTE: No longer rename colorbrewer greys map, just redirect 'grays'
    # to 'greys' in colormap database.
    "Grayscale": (  # assorted origin, but they belong together
        "Greys",
        "Mono",
        "MonoCycle",
    ),
    "Matplotlib sequential": (
        "viridis",
        "plasma",
        "inferno",
        "magma",
        "cividis",
    ),
    "Matplotlib cyclic": ("twilight",),
    "Seaborn sequential": (
        "Rocket",
        "Flare",
        "Mako",
        "Crest",
    ),
    "Seaborn diverging": (
        "IceFire",
        "Vlag",
    ),
    "Proplot sequential": (
        "Fire",
        "Stellar",
        "Glacial",
        "Dusk",
        "Marine",
        "Boreal",
        "Sunrise",
        "Sunset",
    ),
    "Proplot diverging": (
        "Div",
        "NegPos",
        "DryWet",
    ),
    "Other sequential": ("cubehelix", "turbo"),
    "Other diverging": (
        "BR",
        "ColdHot",
        "CoolWarm",
    ),
    "cmOcean sequential": (
        "Oxy",
        "Thermal",
        "Dense",
        "Ice",
        "Haline",
        "Deep",
        "Algae",
        "Tempo",
        "Speed",
        "Turbid",
        "Solar",
        "Matter",
        "Amp",
    ),
    "cmOcean diverging": (
        "Balance",
        "Delta",
        "Curl",
    ),
    "cmOcean cyclic": ("Phase",),
    "Scientific colour maps sequential": (
        "batlow",
        "batlowK",
        "batlowW",
        "devon",
        "davos",
        "oslo",
        "lapaz",
        "acton",
        "lajolla",
        "bilbao",
        "tokyo",
        "turku",
        "bamako",
        "nuuk",
        "hawaii",
        "buda",
        "imola",
        "oleron",
        "bukavu",
        "fes",
    ),
    "Scientific colour maps diverging": (
        "roma",
        "broc",
        "cork",
        "vik",
        "bam",
        "lisbon",
        "tofino",
        "berlin",
        "vanimo",
    ),
    "Scientific colour maps cyclic": (
        "romaO",
        "brocO",
        "corkO",
        "vikO",
        "bamO",
    ),
    "ColorBrewer2.0 sequential": (
        "Purples",
        "Blues",
        "Greens",
        "Oranges",
        "Reds",
        "YlOrBr",
        "YlOrRd",
        "OrRd",
        "PuRd",
        "RdPu",
        "BuPu",
        "PuBu",
        "PuBuGn",
        "BuGn",
        "GnBu",
        "YlGnBu",
        "YlGn",
    ),
    "ColorBrewer2.0 diverging": (
        "Spectral",
        "PiYG",
        "PRGn",
        "BrBG",
        "PuOr",
        "RdGY",
        "RdBu",
        "RdYlBu",
        "RdYlGn",
    ),
    "SciVisColor blues": (
        "Blues1",
        "Blues2",
        "Blues3",
        "Blues4",
        "Blues5",
        "Blues6",
        "Blues7",
        "Blues8",
        "Blues9",
        "Blues10",
        "Blues11",
    ),
    "SciVisColor greens": (
        "Greens1",
        "Greens2",
        "Greens3",
        "Greens4",
        "Greens5",
        "Greens6",
        "Greens7",
        "Greens8",
    ),
    "SciVisColor yellows": (
        "Yellows1",
        "Yellows2",
        "Yellows3",
        "Yellows4",
    ),
    "SciVisColor oranges": (
        "Oranges1",
        "Oranges2",
        "Oranges3",
        "Oranges4",
    ),
    "SciVisColor browns": (
        "Browns1",
        "Browns2",
        "Browns3",
        "Browns4",
        "Browns5",
        "Browns6",
        "Browns7",
        "Browns8",
        "Browns9",
    ),
    "SciVisColor reds": (
        "Reds1",
        "Reds2",
        "Reds3",
        "Reds4",
        "Reds5",
    ),
    "SciVisColor purples": (
        "Purples1",
        "Purples2",
        "Purples3",
    ),
    # Builtin colormaps that re hidden by default. Some are really bad, some
    # are segmented maps that should be cycles, and some are just uninspiring.
    "MATLAB": (
        "bone",
        "cool",
        "copper",
        "autumn",
        "flag",
        "prism",
        "jet",
        "hsv",
        "hot",
        "spring",
        "summer",
        "winter",
        "pink",
        "gray",
    ),
    "GNUplot": (
        "gnuplot",
        "gnuplot2",
        "ocean",
        "afmhot",
        "rainbow",
    ),
    "GIST": (
        "gist_earth",
        "gist_gray",
        "gist_heat",
        "gist_ncar",
        "gist_rainbow",
        "gist_stern",
        "gist_yarg",
    ),
    "Other": (
        "binary",
        "bwr",
        "brg",  # appear to be custom matplotlib
        "Wistia",
        "CMRmap",  # individually released
        "seismic",
        "terrain",
        "nipy_spectral",  # origin ambiguous
        "tab10",
        "tab20",
        "tab20b",
        "tab20c",  # merged colormap cycles
    ),
}

# Docstring snippets
_colorbar_docstring = """
length : unit-spec, optional
    The length of each colorbar.
    %(units.in)s
width : float or str, optional
    The width of each colorbar.
    %(units.in)s
rasterized : bool, default: :rc:`colorbar.rasterized`
    Whether to rasterize the colorbar solids. This increases rendering
    time and decreases file sizes for vector graphics.
"""
docstring._snippet_manager["demos.cmaps"] = ", ".join(f"``{s!r}``" for s in CMAP_TABLE)
docstring._snippet_manager["demos.cycles"] = ", ".join(
    f"``{s!r}``" for s in CYCLE_TABLE
)  # noqa: E501
docstring._snippet_manager["demos.colors"] = ", ".join(
    f"``{s!r}``" for s in COLOR_TABLE
)  # noqa: E501
docstring._snippet_manager["demos.colorbar"] = _colorbar_docstring


def show_channels(
    *args,
    N=100,
    rgb=False,
    saturation=True,
    minhue=0,
    maxsat=500,
    width=100,
    refwidth=1.7,
):
    """
    Show how arbitrary colormap(s) vary with respect to the hue, chroma,
    luminance, HSL saturation, and HPL saturation channels, and optionally
    the red, blue and green channels. Adapted from `this example \
<https://matplotlib.org/stable/tutorials/colors/colormaps.html#lightness-of-matplotlib-colormaps>`__.

    Parameters
    ----------
    *args : colormap-spec, default: :rc:`image.cmap`
        Positional arguments are colormap names or objects.
    N : int, optional
        The number of markers to draw for each colormap.
    rgb : bool, optional
        Whether to also show the red, green, and blue channels in the bottom row.
    saturation : bool, optional
        Whether to show the HSL and HPL saturation channels alongside the raw chroma.
    minhue : float, optional
        The minimum hue. This lets you rotate the hue plot cyclically.
    maxsat : float, optional
        The maximum saturation. Use this to truncate large saturation values.
    width : int, optional
        The width of each colormap line in points.
    refwidth : int or str, optional
        The width of each subplot. Passed to `~proplot.ui.subplots`.

    Returns
    -------
    proplot.figure.Figure
        The figure.
    proplot.gridspec.SubplotGrid
        The subplot grid.

    See also
    --------
    show_cmaps
    show_colorspaces
    """
    # Figure and plot
    if not args:
        raise ValueError("At least one positional argument required.")
    array = [[1, 1, 2, 2, 3, 3]]
    labels = ("Hue", "Chroma", "Luminance")
    if saturation:
        array += [[0, 4, 4, 5, 5, 0]]
        labels += ("HSL saturation", "HPL saturation")
    if rgb:
        array += [np.array([4, 4, 5, 5, 6, 6]) + 2 * int(saturation)]
        labels += ("Red", "Green", "Blue")
    fig, axs = ui.subplots(
        array=array,
        refwidth=refwidth,
        wratios=(1.5, 1, 1, 1, 1, 1.5),
        share="labels",
        span=False,
        innerpad=1,
    )
    # Iterate through colormaps
    mc = ms = mp = 0
    cmaps = []
    for cmap in args:
        # Get colormap and avoid registering new names
        name = cmap if isinstance(cmap, str) else getattr(cmap, "name", None)
        cmap = constructor.Colormap(cmap, N=N)  # arbitrary cmap argument
        if name is not None:
            cmap.name = name
        cmap._init()
        cmaps.append(cmap)

        # Get clipped RGB table
        x = np.linspace(0, 1, N)
        lut = cmap._lut[:-3, :3].copy()
        rgb_data = lut.T  # 3 by N
        hcl_data = np.array([to_xyz(color, space="hcl") for color in lut]).T  # 3 by N
        hsl_data = [to_xyz(color, space="hsl")[1] for color in lut]
        hpl_data = [to_xyz(color, space="hpl")[1] for color in lut]

        # Plot channels
        # If rgb is False, the zip will just truncate the other iterables
        data = tuple(hcl_data)
        if saturation:
            data += (hsl_data, hpl_data)
        if rgb:
            data += tuple(rgb_data)
        for ax, y, label in zip(axs, data, labels):
            ylim, ylocator = None, None
            if label in ("Red", "Green", "Blue"):
                ylim = (0, 1)
                ylocator = 0.2
            elif label == "Luminance":
                ylim = (0, 100)
                ylocator = 20
            elif label == "Hue":
                ylim = (minhue, minhue + 360)
                ylocator = 90
                y = y - 720
                for _ in range(3):  # rotate up to 1080 degrees
                    y[y < minhue] += 360
            else:
                if "HSL" in label:
                    m = ms = max(min(max(ms, max(y)), maxsat), 100)
                elif "HPL" in label:
                    m = mp = max(min(max(mp, max(y)), maxsat), 100)
                else:
                    m = mc = max(min(max(mc, max(y)), maxsat), 100)
                ylim = (0, m)
                ylocator = ("maxn", 5)
            ax.scatter(x, y, c=x, cmap=cmap, s=width, linewidths=0)
            ax.format(title=label, ylim=ylim, ylocator=ylocator)

    # Formatting
    suptitle = (
        ", ".join(repr(cmap.name) for cmap in cmaps[:-1])
        + (", and " if len(cmaps) > 2 else " and " if len(cmaps) == 2 else " ")
        + f"{repr(cmaps[-1].name)} colormap"
        + ("s" if len(cmaps) > 1 else "")
    )
    axs.format(
        xlocator=0.25,
        xformatter="null",
        suptitle=f"{suptitle} by channel",
        ylim=None,
        ytickminor=False,
    )

    # Colorbar on the bottom
    for cmap in cmaps:
        fig.colorbar(
            cmap,
            loc="b",
            span=(2, 5),
            locator="null",
            label=cmap.name,
            labelweight="bold",
        )
    return fig, axs


def show_colorspaces(*, luminance=None, saturation=None, hue=None, refwidth=2):
    """
    Generate hue-saturation, hue-luminance, and luminance-saturation
    cross-sections for the HCL, HSL, and HPL colorspaces.

    Parameters
    ----------
    luminance : float, default: 50
        If passed, saturation-hue cross-sections are drawn for
        this luminance. Must be between ``0`` and ``100``.
    saturation : float, optional
        If passed, luminance-hue cross-sections are drawn for this
        saturation. Must be between ``0`` and ``100``.
    hue : float, optional
        If passed, luminance-saturation cross-sections
        are drawn for this hue. Must be between ``0`` and ``360``.
    refwidth : str or float, optional
        Average width of each subplot. Units are interpreted by
        `~proplot.utils.units`.

    Returns
    -------
    proplot.figure.Figure
        The figure.
    proplot.gridspec.SubplotGrid
        The subplot grid.

    See also
    --------
    show_cmaps
    show_channels
    """
    # Get colorspace properties
    hues = np.linspace(0, 360, 361)
    sats = np.linspace(0, 120, 120)
    lums = np.linspace(0, 99.99, 101)
    if luminance is None and saturation is None and hue is None:
        luminance = 50
    _not_none(luminance=luminance, saturation=saturation, hue=hue)  # warning
    if luminance is not None:
        hsl = np.concatenate(
            (
                np.repeat(hues[:, None], len(sats), axis=1)[..., None],
                np.repeat(sats[None, :], len(hues), axis=0)[..., None],
                np.ones((len(hues), len(sats)))[..., None] * luminance,
            ),
            axis=2,
        )
        suptitle = f"Hue-saturation cross-section for luminance {luminance}"
        xlabel, ylabel = "hue", "saturation"
        xloc, yloc = 60, 20
    elif saturation is not None:
        hsl = np.concatenate(
            (
                np.repeat(hues[:, None], len(lums), axis=1)[..., None],
                np.ones((len(hues), len(lums)))[..., None] * saturation,
                np.repeat(lums[None, :], len(hues), axis=0)[..., None],
            ),
            axis=2,
        )
        suptitle = f"Hue-luminance cross-section for saturation {saturation}"
        xlabel, ylabel = "hue", "luminance"
        xloc, yloc = 60, 20
    elif hue is not None:
        hsl = np.concatenate(
            (
                np.ones((len(lums), len(sats)))[..., None] * hue,
                np.repeat(sats[None, :], len(lums), axis=0)[..., None],
                np.repeat(lums[:, None], len(sats), axis=1)[..., None],
            ),
            axis=2,
        )
        suptitle = "Luminance-saturation cross-section"
        xlabel, ylabel = "luminance", "saturation"
        xloc, yloc = 20, 20

    # Make figure, with black indicating invalid values
    # Note we invert the x-y ordering for imshow
    fig, axs = ui.subplots(refwidth=refwidth, ncols=3, share=False, innerpad=0.5)
    for ax, space in zip(axs, ("hcl", "hsl", "hpl")):
        rgba = np.ones((*hsl.shape[:2][::-1], 4))  # RGBA
        for j in range(hsl.shape[0]):
            for k in range(hsl.shape[1]):
                rgb_jk = to_rgb(hsl[j, k, :], space)
                if not all(0 <= c <= 1 for c in rgb_jk):
                    rgba[k, j, 3] = 0  # black cell
                else:
                    rgba[k, j, :3] = rgb_jk
        ax.imshow(rgba, origin="lower", aspect="auto")
        ax.format(
            xlabel=xlabel,
            ylabel=ylabel,
            suptitle=suptitle,
            grid=False,
            xtickminor=False,
            ytickminor=False,
            xlocator=xloc,
            ylocator=yloc,
            facecolor="k",
            title=space.upper(),
        )
    return fig, axs


@warnings._rename_kwargs("0.8.0", categories="include")
@warnings._rename_kwargs("0.10.0", rasterize="rasterized")
def _draw_bars(
    cmaps,
    *,
    source,
    unknown="User",
    include=None,
    ignore=None,
    length=4.0,
    width=0.2,
    N=None,
    rasterized=None,
):
    """
    Draw colorbars for "colormaps" and "color cycles". This is called by
    `show_cycles` and `show_cmaps`.
    """
    # Categorize the input names
    table = {unknown: []} if unknown else {}
    table.update({cat: [None] * len(names) for cat, names in source.items()})
    for cmap in cmaps:
        cat = None
        name = cmap.name or "_no_name"
        name = name.lower()
        for opt, names in source.items():
            names = list(map(str.lower, names))
            if name in names:
                i, cat = names.index(name), opt
        if cat:
            table[cat][i] = cmap
        elif unknown:
            table[unknown].append(cmap)

    # Filter out certain categories
    options = set(map(str.lower, source))
    if ignore is None:
        ignore = ("matlab", "gnuplot", "gist", "other")
    if isinstance(include, str):
        include = (include,)
    if isinstance(ignore, str):
        ignore = (ignore,)
    if include is None:
        include = options - set(map(str.lower, ignore))
    else:
        include = set(map(str.lower, include))
    if any(cat not in options and cat != unknown for cat in include):
        raise ValueError(
            f"Invalid categories {include!r}. Options are: "
            + ", ".join(map(repr, source))
            + "."
        )
    for cat in tuple(table):
        table[cat][:] = [cmap for cmap in table[cat] if cmap is not None]
        if not table[cat] or cat.lower() not in include and cat != unknown:
            del table[cat]

    # Draw figure
    # Allocate two colorbar widths for each title of sections
    naxs = 2 * len(table) + sum(map(len, table.values()))
    fig, axs = ui.subplots(
        refwidth=length,
        refheight=width,
        nrows=naxs,
        share=False,
        hspace="2pt",
        top="-1em",
    )
    i = -1
    nheads = nbars = 0  # for deciding which axes to plot in
    for cat, cmaps in table.items():
        nheads += 1
        for j, cmap in enumerate(cmaps):
            i += 1
            if j + nheads + nbars > naxs:
                break
            if j == 0:  # allocate this axes for title
                i += 2
                for ax in axs[i - 2 : i]:
                    ax.set_visible(False)
            ax = axs[i]
            if N is not None:
                cmap = cmap.copy(N=N)
            label = cmap.name
            label = re.sub(r"\A_*", "", label)
            label = re.sub(r"(_copy)*\Z", "", label)
            ax.colorbar(
                cmap,
                loc="fill",
                orientation="horizontal",
                locator="null",
                linewidth=0,
                rasterized=rasterized,
            )
            ax.text(
                0 - (rc["axes.labelpad"] / 72) / length,
                0.45,
                label,
                ha="right",
                va="center",
                transform="axes",
            )
            if j == 0:
                ax.set_title(cat, weight="bold")
        nbars += len(cmaps)

    return fig, axs


@docstring._snippet_manager
def show_cmaps(*args, **kwargs):
    """
    Generate a table of the registered colormaps or the input colormaps
    categorized by source. Adapted from `this example \
<http://matplotlib.org/stable/gallery/color/colormap_reference.html>`__.

    Parameters
    ----------
    *args : colormap-spec, optional
        Colormap names or objects.
    N : int, default: :rc:`image.lut`
        The number of levels in each colorbar.
    unknown : str, default: 'User'
        Category name for colormaps that are unknown to proplot.
        Set this to ``False`` to hide unknown colormaps.
    include : str or sequence of str, default: None
        Category names to be shown in the table. Use this to limit
        the table to a subset of categories. Valid categories are
        %(demos.cmaps)s.
    ignore : str or sequence of str, default: 'MATLAB', 'GNUplot', 'GIST', 'Other'
        Used only if `include` was not passed. Category names to be removed from the
        table. Use of the default ignored colormaps is discouraged because they contain
        non-uniform color transitions (see the :ref:`user guide <ug_perceptual>`).
    %(demos.colorbar)s

    Returns
    -------
    proplot.figure.Figure
        The figure.
    proplot.gridspec.SubplotGrid
        The subplot grid.

    See also
    --------
    show_colorspaces
    show_channels
    show_cycles
    show_colors
    show_fonts
    """
    # Get the list of colormaps
    if args:
        cmaps = list(map(constructor.Colormap, args))
        cmaps = [
            (
                cmap
                if isinstance(cmap, mcolors.LinearSegmentedColormap)
                else pcolors._get_cmap_subtype(cmap, "continuous")
            )
            for cmap in args
        ]
        ignore = ()
    else:
        cmaps = [
            cmap
            for cmap in pcolors._cmap_database.values()
            if isinstance(cmap, pcolors.ContinuousColormap)
            and not (cmap.name or "_")[:1] == "_"
        ]
        ignore = None

    # Return figure of colorbars
    kwargs.setdefault("source", CMAP_TABLE)
    kwargs.setdefault("ignore", ignore)
    return _draw_bars(cmaps, **kwargs)


@docstring._snippet_manager
def show_cycles(*args, **kwargs):
    """
    Generate a table of registered color cycles or the input color cycles
    categorized by source. Adapted from `this example \
<http://matplotlib.org/stable/gallery/color/colormap_reference.html>`__.

    Parameters
    ----------
    *args : colormap-spec, optional
        Cycle names or objects.
    unknown : str, default: 'User'
        Category name for cycles that are unknown to proplot.
        Set this to ``False`` to hide unknown colormaps.
    include : str or sequence of str, default: None
        Category names to be shown in the table. Use this to limit
        the table to a subset of categories. Valid categories are
        %(demos.cycles)s.
    ignore : str or sequence of str, default: None
        Used only if `include` was not passed. Category names to be removed
        from the table.
    %(demos.colorbar)s

    Returns
    -------
    proplot.figure.Figure
        The figure.
    proplot.gridspec.SubplotGrid
        The subplot grid.

    See also
    --------
    show_cmaps
    show_colors
    show_fonts
    """
    # Get the list of cycles
    if args:
        cycles = [
            (
                pcolors.DiscreteColormap(
                    cmap.by_key().get("color", ["k"]), name=getattr(cmap, "name", None)
                )
                if isinstance(cmap, cycler.Cycler)
                else (
                    cmap
                    if isinstance(cmap, mcolors.ListedColormap)
                    else pcolors._get_cmap_subtype(cmap, "discrete")
                )
            )
            for cmap in args
        ]
        ignore = ()
    else:
        cycles = [
            cmap
            for cmap in pcolors._cmap_database.values()
            if isinstance(cmap, pcolors.DiscreteColormap)
            and not (cmap.name or "_")[:1] == "_"
        ]
        ignore = None

    # Return figure of colorbars
    kwargs.setdefault("source", CYCLE_TABLE)
    kwargs.setdefault("ignore", ignore)
    return _draw_bars(cycles, **kwargs)


def _filter_colors(hcl, ihue, nhues, minsat):
    """
    Filter colors into categories.

    Parameters
    ----------
    hcl : tuple
        The data.
    ihue : int
        The hue column.
    nhues : int
        The total number of hues.
    minsat : float
        The minimum saturation used for the "grays" column.
    """
    breakpoints = np.linspace(0, 360, nhues)
    gray = hcl[1] <= minsat
    if ihue == 0:
        return gray
    color = breakpoints[ihue - 1] <= hcl[0] < breakpoints[ihue]
    if ihue == nhues - 1:
        color = color or color == breakpoints[ihue]  # endpoint inclusive
    return not gray and color


@docstring._snippet_manager
def show_colors(*, nhues=17, minsat=10, unknown="User", include=None, ignore=None):
    """
    Generate tables of the registered color names. Adapted from
    `this example <https://matplotlib.org/examples/color/named_colors.html>`__.

    Parameters
    ----------
    nhues : int, optional
        The number of breaks between hues for grouping "like colors" in the
        color table.
    minsat : float, optional
        The threshold saturation, between ``0`` and ``100``, for designating
        "gray colors" in the color table.
    unknown : str, default: 'User'
        Category name for color names that are unknown to proplot.
        Set this to ``False`` to hide unknown color names.
    include : str or sequence of str, default: None
        Category names to be shown in the table. Use this to limit
        the table to a subset of categories. Valid categories are
        %(demos.colors)s.
    ignore : str or sequence of str, default: 'CSS4'
        Used only if `include` was not passed. Category names to be removed
        from the colormap table.

    Returns
    -------
    proplot.figure.Figure
        The figure.
    proplot.gridspec.SubplotGrid
        The subplot grid.
    """
    # Tables of known colors to be plotted
    colordict = {}
    if ignore is None:
        ignore = "css4"
    if isinstance(include, str):
        include = (include.lower(),)
    if isinstance(ignore, str):
        ignore = (ignore.lower(),)
    if include is None:
        include = COLOR_TABLE.keys()
        include -= set(map(str.lower, ignore))
    for cat in sorted(include):
        if cat not in COLOR_TABLE:
            raise ValueError(
                f"Invalid categories {include!r}. Options are: "
                + ", ".join(map(repr, COLOR_TABLE))
                + "."
            )
        colordict[cat] = list(COLOR_TABLE[cat])  # copy the names

    # Add "unknown" colors
    if unknown:
        unknown_colors = [
            color
            for color in map(repr, pcolors._color_database)
            if "xkcd:" not in color
            and "tableau:" not in color
            and not any(color in list_ for list_ in COLOR_TABLE)
        ]
        if unknown_colors:
            colordict[unknown] = unknown_colors

    # Divide colors into columns and rows
    # For base and open colors, tables are already organized into like
    # colors, so just reshape them into grids. For other colors, we group
    # them by hue in descending order of luminance.
    namess = {}
    for cat in sorted(include):
        if cat == "base":
            names = np.asarray(colordict[cat])
            ncols, nrows = len(names), 1
        elif cat == "opencolor":
            names = np.asarray(colordict[cat])
            ncols, nrows = 7, 20
        else:
            hclpairs = [(name, to_xyz(name, "hcl")) for name in colordict[cat]]
            hclpairs = [
                sorted(
                    [
                        pair
                        for pair in hclpairs
                        if _filter_colors(pair[1], ihue, nhues, minsat)
                    ],
                    key=lambda x: x[1][2],  # sort by luminance
                )
                for ihue in range(nhues)
            ]
            names = np.array([name for ipairs in hclpairs for name, _ in ipairs])
            ncols, nrows = 4, len(names) // 4 + 1

        names.resize((ncols, nrows))  # fill empty slots with empty string
        namess[cat] = names

    # Draw figures for different groups of colors
    # NOTE: Aspect ratios should be number of columns divided by number
    # of rows, times the aspect ratio of the slot for each swatch-name
    # pair, which we set to 5.
    shape = tuple(namess.values())[0].shape  # sample *first* group
    figwidth = 6.5
    refaspect = (figwidth * 72) / (10 * shape[1])  # points
    maxcols = max(names.shape[0] for names in namess.values())
    hratios = tuple(names.shape[1] for names in namess.values())
    fig, axs = ui.subplots(
        figwidth=figwidth,
        refaspect=refaspect,
        nrows=len(include),
        hratios=hratios,
    )
    title_dict = {
        "css4": "CSS4 colors",
        "base": "Base colors",
        "opencolor": "Open color",
        "xkcd": "XKCD colors",
    }
    for ax, (cat, names) in zip(axs, namess.items()):
        # Format axes
        ax.format(
            title=title_dict.get(cat, cat),
            titleweight="bold",
            xlim=(0, maxcols - 1),
            ylim=(0, names.shape[1]),
            grid=False,
            yloc="neither",
            xloc="neither",
            alpha=0,
        )

        # Draw swatches as lines
        lw = 8  # best to just use trial and error
        swatch = 0.45  # percent of column reserved for swatch
        ncols, nrows = names.shape
        for col, inames in enumerate(names):
            for row, name in enumerate(inames):
                if not name:
                    continue
                y = nrows - row - 1  # start at top
                x1 = col * (maxcols - 1) / ncols  # e.g. idx 3 --> idx 7
                x2 = x1 + swatch  # portion of column
                xtext = x1 + 1.1 * swatch
                ax.text(
                    xtext,
                    y,
                    name,
                    ha="left",
                    va="center",
                    transform="data",
                    clip_on=False,
                )
                ax.plot(
                    [x1, x2],
                    [y, y],
                    color=name,
                    lw=lw,
                    solid_capstyle="butt",  # do not stick out
                    clip_on=False,
                )

    return fig, axs


def show_fonts(
    *args, family=None, user=None, text=None, math=False, fallback=False, **kwargs
):
    """
    Generate a table of fonts. If a glyph for a particular font is unavailable,
    it is replaced with the "¤" dummy character.

    Parameters
    ----------
    *args : str or `~matplotlib.font_manager.FontProperties`
        The font specs, font names, or `~matplotlib.font_manager.FontProperties`\\ s
        to show. If no positional arguments are passed and the `family` argument is
        not passed, then the fonts found in `~proplot.config.Configurator.user_folder`
        and `~proplot.config.Configurator.local_folders` and the *available*
        :rcraw:`font.sans-serif` fonts are shown.
    family \
: {'tex-gyre', 'sans-serif', 'serif', 'monospace', 'cursive', 'fantasy'}, optional
        The family from which *available* fonts are shown. Default is ``'sans-serif'``
        if no arguments were provided. Otherwise the default is to not show family
        fonts. The fonts belonging to each family are listed under :rcraw:`font.serif`,
        :rcraw:`font.sans-serif`, :rcraw:`font.monospace`, :rcraw:`font.cursive`, and
        :rcraw:`font.fantasy`. The special family ``'tex-gyre'`` includes the
        `TeX Gyre <http://www.gust.org.pl/projects/e-foundry/tex-gyre>`__ fonts.
    user : bool, optional
        Whether to include fonts in `~proplot.config.Configurator.user_folder` and
        `~proplot.config.Configurator.local_folders` at the top of the table. Default
        is ``True`` if called without any arguments and ``False`` otherwise.
    text : str, optional
        The sample text shown for each font. If not passed then default math or
        non-math sample text is used.
    math : bool, default: False
        Whether the default sample text should show non-math Latin characters or
        or math equations and Greek letters.
    fallback : bool, default: False
        Whether to use the fallback font :rcraw:`mathtext.fallback` for unavailable
        characters. If ``False`` the dummy glyph "¤" is shown for missing characters.
    **kwargs
        Additional font properties passed to `~matplotlib.font_manager.FontProperties`.
        Default size is ``12`` and default weight, style, and strength are ``'normal'``.

    Other parameters
    ----------------
    size : float, default: 12
        The font size.
    weight : str, default: 'normal'
        The font weight.
    style : str, default: 'normal'
        The font style.
    stretch : str, default: 'normal'
        The font stretch.

    Returns
    -------
    proplot.figure.Figure
        The figure.
    proplot.gridspec.SubplotGrid
        The subplot grid.

    See also
    --------
    show_cmaps
    show_cycles
    show_colors
    """
    # Parse user input fonts and translate into FontProperties.
    s = set()
    props = []  # should be string names
    all_fonts = sorted(mfonts.fontManager.ttflist, key=lambda font: font.name)
    all_fonts = [
        font for font in all_fonts if font.name not in s and not s.add(font.name)
    ]  # noqa: E501
    all_names = [font.name for font in all_fonts]
    for arg in args:
        if isinstance(arg, str):
            arg = mfonts.FontProperties(arg, **kwargs)  # possibly a fontspec
        elif not isinstance(arg, mfonts.FontProperties):
            raise TypeError(f"Expected string or FontProperties but got {type(arg)}.")
        opts = arg.get_family()  # usually a singleton list
        if opts and opts[0] in all_names:
            props.append(arg)
        else:
            warnings._warn_proplot(f"Input font name {opts[:1]!r} not found. Skipping.")

    # Add user and family FontProperties.
    user = _not_none(user, not args and family is None)
    family = _not_none(family, None if args else "sans-serif")
    if user:
        paths = _get_data_folders("fonts", default=False)
        for font in all_fonts:  # fonts sorted by unique name
            if os.path.dirname(font.fname) in paths:
                props.append(mfonts.FontProperties(font.name, **kwargs))
    if family is not None:
        options = ("serif", "sans-serif", "monospace", "cursive", "fantasy", "tex-gyre")
        if family not in options:
            raise ValueError(
                f"Invalid font family {family!r}. Options are: "
                + ", ".join(map(repr, options))
                + "."
            )
        names = FAMILY_TEXGYRE if family == "tex-gyre" else rc["font." + family]
        for name in names:
            if name in all_names:  # valid font name
                props.append(mfonts.FontProperties(name, **kwargs))

    # The default sample text
    linespacing = 0.8 if text is None and math else 1.2
    if text is None:
        if not math:
            text = (
                "the quick brown fox jumps over a lazy dog 01234 ; . , + - * ^ () ||"
                "\n"
                "THE QUICK BROWN FOX JUMPS OVER A LAZY DOG 56789 : ! ? & # % $ [] {}"
            )
        else:
            text = (
                "\n"
                r"$\alpha\beta$ $\Gamma\gamma$ $\Delta\delta$ "
                r"$\epsilon\zeta\eta$ $\Theta\theta$ $\kappa\mu\nu$ "
                r"$\Lambda\lambda$ $\Pi\pi$ $\xi\rho\tau\chi$ $\Sigma\sigma$ "
                r"$\Phi\phi$ $\Psi\psi$ $\Omega\omega$ "
                r"$\{ \; \}^i$ $[ \; ]_j$ $( \; )^k$ $\left<  \right>_n$"
                "\n"
                r"$0^a + 1_b - 2^c \times 3_d = "
                r"4.0^e \equiv 5.0_f \approx 6.0^g \sim 7_h \leq 8^i \geq 9_j"
                r"\ll \prod \, P \gg \sum \, Q \, "
                r"\int \, Y \mathrm{d}y \propto \oint \;\, Z \mathrm{d}z$"
            )

    # Settings for rendering math text
    ctx = {"mathtext.fontset": "custom"}
    if not fallback:
        if _version_mpl < "3.4":
            ctx["mathtext.fallback_to_cm"] = False
        else:
            ctx["mathtext.fallback"] = None
    if "size" not in kwargs:
        for prop in props:
            if prop.get_size() == rc["font.size"]:
                prop.set_size(12)  # only if fontspec did not change the size

    # Create figure
    refsize = props[0].get_size_in_points() if props else rc["font.size"]
    refheight = 1.2 * (text.count("\n") + 2.5) * refsize / 72
    fig, axs = ui.subplots(
        refwidth=4.5,
        refheight=refheight,
        nrows=len(props),
        ncols=1,
        space=0,
    )
    fig._render_context.update(ctx)
    fig.format(
        xloc="neither", yloc="neither", xlocator="null", ylocator="null", alpha=0
    )
    for ax, prop in zip(axs, props):
        name = prop.get_family()[0]
        ax.text(
            0,
            0.5,
            f"{name}:\n{text} ",
            ha="left",
            va="center",
            linespacing=linespacing,
            fontproperties=prop,
        )
    return fig, axs
