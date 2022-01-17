#!/usr/bin/env python3
"""
The first-level axes subclass used for all proplot figures.
Implements basic shared functionality.
"""
import copy
import inspect
import itertools
import re
from numbers import Integral

import matplotlib.axes as maxes
import matplotlib.axis as maxis
import matplotlib.cm as mcm
import matplotlib.colors as mcolors
import matplotlib.container as mcontainer
import matplotlib.contour as mcontour
import matplotlib.legend as mlegend
import matplotlib.offsetbox as moffsetbox
import matplotlib.patches as mpatches
import matplotlib.projections as mprojections
import matplotlib.text as mtext
import matplotlib.ticker as mticker
import matplotlib.transforms as mtransforms
import numpy as np
from matplotlib import cbook

from .. import colors as pcolors
from .. import constructor
from .. import ticker as pticker
from ..config import rc
from ..internals import ic  # noqa: F401
from ..internals import (
    _kwargs_to_args,
    _not_none,
    _pop_kwargs,
    _pop_params,
    _pop_props,
    _pop_rc,
    _translate_loc,
    _version_mpl,
    docstring,
    guides,
    labels,
    rcsetup,
    warnings,
)
from ..utils import _fontsize_to_pt, edges, units

try:
    from cartopy.crs import CRS, PlateCarree
except Exception:
    CRS = PlateCarree = object

__all__ = ['Axes']


# A-b-c label string
ABC_STRING = 'abcdefghijklmnopqrstuvwxyz'

# Legend align options
ALIGN_OPTS = {
    None: {
        'center': 'center',
        'left': 'center left',
        'right': 'center right',
        'top': 'upper center',
        'bottom': 'lower center',
    },
    'left': {
        'top': 'upper right',
        'center': 'center right',
        'bottom': 'lower right',
    },
    'right': {
        'top': 'upper left',
        'center': 'center left',
        'bottom': 'lower left',
    },
    'top': {
        'left': 'lower left',
        'center': 'lower center',
        'right': 'lower right'
    },
    'bottom': {
        'left': 'upper left',
        'center': 'upper center',
        'right': 'upper right'
    },
}


# Projection docstring
_proj_docstring = """
proj, projection : \
str, `cartopy.crs.Projection`, or `~mpl_toolkits.basemap.Basemap`, optional
    The map projection specification(s). If ``'cart'`` or ``'cartesian'``
    (the default), a `~proplot.axes.CartesianAxes` is created. If ``'polar'``,
    a `~proplot.axes.PolarAxes` is created. Otherwise, the argument is
    interpreted by `~proplot.constructor.Proj`, and the result is used
    to make a `~proplot.axes.GeoAxes` (in this case the argument can be
    a `cartopy.crs.Projection` instance, a `~mpl_toolkits.basemap.Basemap`
    instance, or a projection name listed in :ref:`this table <proj_table>`).
"""
_proj_kw_docstring = """
proj_kw, projection_kw : dict-like, optional
    Keyword arguments passed to `~mpl_toolkits.basemap.Basemap` or
    cartopy `~cartopy.crs.Projection` classes on instantiation.
"""
_backend_docstring = """
backend : {'cartopy', 'basemap'}, default: :rc:`geo.backend`
    Whether to use `~mpl_toolkits.basemap.Basemap` or
    `~cartopy.crs.Projection` for map projections.
"""
docstring._snippet_manager['axes.proj'] = _proj_docstring
docstring._snippet_manager['axes.proj_kw'] = _proj_kw_docstring
docstring._snippet_manager['axes.backend'] = _backend_docstring


# Colorbar and legend space
_space_docstring = """
queue : bool, optional
    If ``True`` and `loc` is the same as an existing {name}, the input
    arguments are added to a queue and this function returns ``None``.
    This is used to "update" the same {name} with successive ``ax.{name}(...)``
    calls. If ``False`` (the default) and `loc` is the same as an existing
    *inset* {name}, the old {name} is removed. If ``False`` and `loc` is an
    *outer* {name}, the {name}s are "stacked".
space : unit-spec, default: None
    For outer {name}s only. The fixed space between the {name} and the subplot
    edge. %(units.em)s
    When the :ref:`tight layout algorithm <ug_tight>` is active for the figure,
    `space` is computed automatically (see `pad`). Otherwise, `space` is set to
    a suitable default.
pad : unit-spec, default: :rc:`subplots.panelpad` or :rc:`{default}`
    For outer {name}s, this is the :ref:`tight layout padding <ug_tight>`
    between the {name} and the subplot (default is :rcraw:`subplots.panelpad`).
    For inset {name}s, this is the fixed space between the axes
    edge and the {name} (default is :rcraw:`{default}`).
    %(units.em)s
align : {{'center', 'top', 'bottom', 'left', 'right', 't', 'b', 'l', 'r'}}, optional
    For outer {name}s only. How to align the {name} against the subplot edge.
    The values ``'top'`` and ``'bottom'`` are valid for left and right {name}s
    and ``'left'`` and ``'right'`` are valid for top and bottom {name}s.
    The default is always ``'center'``.
"""
docstring._snippet_manager['axes.legend_space'] = _space_docstring.format(
    name='legend', default='legend.borderaxespad'
)
docstring._snippet_manager['axes.colorbar_space'] = _space_docstring.format(
    name='colorbar', default='colorbar.insetpad'
)


# Transform docstring
# Used for text and add_axes
_transform_docstring = """
transform : {'data', 'axes', 'figure', 'subfigure'} \
or `~matplotlib.transforms.Transform`, optional
    The transform used to interpret the bounds. Can be a
    `~matplotlib.transforms.Transform` instance or a string representing
    the `~matplotlib.axes.Axes.transData`, `~matplotlib.axes.Axes.transAxes`,
    `~matplotlib.figure.Figure.transFigure`, or
    `~matplotlib.figure.Figure.transSubfigure`, transforms.
"""
docstring._snippet_manager['axes.transform'] = _transform_docstring


# Inset docstring
# NOTE: Used by SubplotGrid.inset_axes
_inset_docstring = """
Add an inset axes.
This is similar to `matplotlib.axes.Axes.inset_axes`.

Parameters
----------
bounds : 4-tuple of float
    The (left, bottom, width, height) coordinates for the axes.
%(axes.transform)s
    Default is to use the same projection as the current axes.
%(axes.proj)s
%(axes.proj_kw)s
%(axes.backend)s
zorder : float, default: 4
    The `zorder <https://matplotlib.org/stable/gallery/misc/zorder_demo.html>`__
    of the axes. Should be greater than the zorder of elements in the parent axes.
zoom : bool, default: True or False
    Whether to draw lines indicating the inset zoom using `~Axes.indicate_inset_zoom`.
    The line positions will automatically adjust when the parent or inset axes limits
    change. Default is ``True`` only if both axes are `~proplot.axes.CartesianAxes`.
zoom_kw : dict, optional
    Passed to `~Axes.indicate_inset_zoom`.

Other parameters
----------------
**kwargs
    Passed to `proplot.axes.Axes`.

Returns
-------
proplot.axes.Axes
    The inset axes.

See also
--------
Axes.indicate_inset_zoom
matplotlib.axes.Axes.inset_axes
matplotlib.axes.Axes.indicate_inset
matplotlib.axes.Axes.indicate_inset_zoom
"""
_indicate_inset_docstring = """
Add indicators denoting the zoom range of the inset axes.
This will replace previously drawn zoom indicators.

Parameters
----------
%(artist.patch)s
zorder : float, default: 3.5
    The `zorder <https://matplotlib.org/stable/gallery/misc/zorder_demo.html>`__ of
    the indicators. Should be greater than the zorder of elements in the parent axes.

Other parameters
----------------
**kwargs
    Passed to `~matplotlib.patches.Patch`.

Note
----
This command must be called from the inset axes rather than the parent axes.
It is called automatically when ``zoom=True`` is passed to `~Axes.inset_axes`
and whenever the axes are drawn (so the line positions always track the axis
limits even if they are later changed).

See also
--------
matplotlib.axes.Axes.indicate_inset
matplotlib.axes.Axes.indicate_inset_zoom
"""
docstring._snippet_manager['axes.inset'] = _inset_docstring
docstring._snippet_manager['axes.indicate_inset'] = _indicate_inset_docstring


# Panel docstring
# NOTE: Used by SubplotGrid.panel_axes
_panel_loc_docstring = """
    ==========  =====================
    Location    Valid keys
    ==========  =====================
    left        ``'left'``, ``'l'``
    right       ``'right'``, ``'r'``
    bottom      ``'bottom'``, ``'b'``
    top         ``'top'``, ``'t'``
    ==========  =====================
"""
_panel_docstring = """
Add a panel axes.

Parameters
----------
side : str, optional
    The panel location. Valid location keys are as follows.

%(axes.panel_loc)s

width : unit-spec, default: :rc:`subplots.panelwidth`
    The panel width.
    %(units.in)s
space : unit-spec, default: None
    The fixed space between the panel and the subplot edge.
    %(units.em)s
    When the :ref:`tight layout algorithm <ug_tight>` is active for the figure,
    `space` is computed automatically (see `pad`). Otherwise, `space` is set to
    a suitable default.
pad : unit-spec, default: :rc:`subplots.panelpad`
    The :ref:`tight layout padding <ug_tight>` between the panel and the subplot.
    %(units.em)s
share : bool, default: True
    Whether to enable axis sharing between the *x* and *y* axes of the
    main subplot and the panel long axes for each panel in the "stack".
    Sharing between the panel short axis and other panel short axes
    is determined by figure-wide `sharex` and `sharey` settings.

Other parameters
----------------
**kwargs
    Passed to `proplot.axes.CartesianAxes`.

Returns
-------
proplot.axes.CartesianAxes
    The panel axes.
"""
docstring._snippet_manager['axes.panel_loc'] = _panel_loc_docstring
docstring._snippet_manager['axes.panel'] = _panel_docstring


# Format docstrings
_axes_format_docstring = """
title : str or sequence, optional
    The axes title. Can optionally be a sequence strings, in which case
    the title will be selected from the sequence according to `~Axes.number`.
abc : bool or str or sequence, default: :rc:`abc`
    The "a-b-c" subplot label style. Must contain the character ``a`` or ``A``,
    for example ``'a.'``, or ``'A'``. If ``True`` then the default style of
    ``'a'`` is used. The ``a`` or ``A`` is replaced with the alphabetic character
    matching the `~Axes.number`. If `~Axes.number` is greater than 26, the
    characters loop around to a, ..., z, aa, ..., zz, aaa, ..., zzz, etc.
    Can also be a sequence of strings, in which case the "a-b-c" label
    will simply be selected from the sequence according to `~Axes.number`.
abcloc, titleloc : str, default: :rc:`abc.loc`, :rc:`title.loc`
    Strings indicating the location for the a-b-c label and main title.
    The following locations are valid:

    .. _title_table:

    ========================  ============================
    Location                  Valid keys
    ========================  ============================
    center above axes         ``'center'``, ``'c'``
    left above axes           ``'left'``, ``'l'``
    right above axes          ``'right'``, ``'r'``
    lower center inside axes  ``'lower center'``, ``'lc'``
    upper center inside axes  ``'upper center'``, ``'uc'``
    upper right inside axes   ``'upper right'``, ``'ur'``
    upper left inside axes    ``'upper left'``, ``'ul'``
    lower left inside axes    ``'lower left'``, ``'ll'``
    lower right inside axes   ``'lower right'``, ``'lr'``
    ========================  ============================

abcborder, titleborder : bool, default: :rc:`abc.border` and :rc:`title.border`
    Whether to draw a white border around titles and a-b-c labels positioned
    inside the axes. This can help them stand out on top of artists
    plotted inside the axes.
abcbbox, titlebbox : bool, default: :rc:`abc.bbox` and :rc:`title.bbox`
    Whether to draw a white bbox around titles and a-b-c labels positioned
    inside the axes. This can help them stand out on top of artists plotted
    inside the axes.
abc_kw, title_kw : dict-like, optional
    Additional settings used to update the a-b-c label and title
    with ``text.update()``.
titlepad : float, default: :rc:`title.pad`
    The padding for the inner and outer titles and a-b-c labels.
    %(units.pt)s
titleabove : bool, default: :rc:`title.above`
    Whether to try to put outer titles and a-b-c labels above panels,
    colorbars, or legends that are above the axes.
abctitlepad : float, default: :rc:`abc.titlepad`
    The horizontal padding between a-b-c labels and titles in the same location.
    %(units.pt)s
ltitle, ctitle, rtitle, ultitle, uctitle, urtitle, lltitle, lctitle, lrtitle \
: str or sequence, optional
    Shorthands for the below keywords.
lefttitle, centertitle, righttitle, upperlefttitle, uppercentertitle, upperrighttitle, \
lowerlefttitle, lowercentertitle, lowerrighttitle : str or sequence, optional
    Additional titles in specific positions (see `title` for details). This works as
    an alternative to the ``ax.format(title='Title', titleloc=loc)`` workflow and
    permits adding more than one title-like label for a single axes.
a, alpha, fc, facecolor, ec, edgecolor, lw, linewidth, ls, linestyle : default: \
:rc:`axes.alpha`, :rc:`axes.facecolor`, :rc:`axes.edgecolor`, :rc:`axes.linewidth`, '-'
    Additional settings applied to the background patch, and their
    shorthands. Their defaults values are the ``'axes'`` properties.
"""
_figure_format_docstring = """
rowlabels, collabels, llabels, tlabels, rlabels, blabels
    Aliases for `leftlabels` and `toplabels`, and for `leftlabels`,
    `toplabels`, `rightlabels`, and `bottomlabels`, respectively.
leftlabels, toplabels, rightlabels, bottomlabels : sequence of str, optional
    Labels for the subplots lying along the left, top, right, and
    bottom edges of the figure. The length of each list must match
    the number of subplots along the corresponding edge.
leftlabelpad, toplabelpad, rightlabelpad, bottomlabelpad : float or unit-spec, default\
: :rc:`leftlabel.pad`, :rc:`toplabel.pad`, :rc:`rightlabel.pad`, :rc:`bottomlabel.pad`
    The padding between the labels and the axes content.
    %(units.pt)s
leftlabels_kw, toplabels_kw, rightlabels_kw, bottomlabels_kw : dict-like, optional
    Additional settings used to update the labels with ``text.update()``.
figtitle
    Alias for `suptitle`.
suptitle : str, optional
    The figure "super" title, centered between the left edge of the lefmost
    column of subplots and the right edge of the rightmost column of subplots, and
    automatically offset above figure titles. This is an improvement on matplotlib's
    "super" title, which just centers the text between figure edges.
suptitlepad : float, default: :rc:`suptitle.pad`
    The padding between the super title and the axes content.
    %(units.pt)s
suptitle_kw : optional
    Additional settings used to update the super title with ``text.update()``.
includepanels : bool, default: False
    Whether to include panels when aligning figure "super titles" along the top
    of the subplot grid and when aligning the `spanx` *x* axis labels and
    `spany` *y* axis labels along the sides of the subplot grid.
"""
_rc_init_docstring = """
"""
_rc_format_docstring = """
rc_mode : int, optional
    The context mode passed to `~proplot.config.Configurator.context`.
rc_kw : dict-like, optional
    An alternative to passing extra keyword arguments. See below.
**kwargs
    {}Keyword arguments that match the name of an `~proplot.config.rc` setting are
    passed to `proplot.config.Configurator.context` and used to update the axes.
    If the setting name has "dots" you can simply omit the dots. For example,
    ``abc='A.'`` modifies the :rcraw:`abc` setting, ``titleloc='left'`` modifies the
    :rcraw:`title.loc` setting, ``gridminor=True`` modifies the :rcraw:`gridminor`
    setting, and ``gridbelow=True`` modifies the :rcraw:`grid.below` setting. Many
    of the keyword arguments documented above are internally applied by retrieving
    settings passed to `~proplot.config.Configurator.context`.
"""
docstring._snippet_manager['rc.init'] = _rc_format_docstring.format(
    'Remaining keyword arguments are passed to `matplotlib.axes.Axes`.\n    '
)
docstring._snippet_manager['rc.format'] = _rc_format_docstring.format('')
docstring._snippet_manager['axes.format'] = _axes_format_docstring
docstring._snippet_manager['figure.format'] = _figure_format_docstring


# Colorbar docstrings
_colorbar_args_docstring = """
mappable : mappable, colormap-spec, sequence of color-spec, \
or sequence of `~matplotlib.artist.Artist`
    There are four options here:

    1. A `~matplotlib.cm.ScalarMappable` (e.g., an object returned by
       `~proplot.axes.PlotAxes.contourf` or `~proplot.axes.PlotAxes.pcolormesh`).
    2. A `~matplotlib.colors.Colormap` or registered colormap name used to build a
       `~matplotlib.cm.ScalarMappable` on-the-fly. The colorbar range and ticks depend
       on the arguments `values`, `vmin`, `vmax`, and `norm`. The default for a
       `~proplot.colors.ContinuousColormap` is ``vmin=0`` and ``vmax=1`` (note that
       passing `values` will "discretize" the colormap). The default for a
       `~proplot.colors.DiscreteColormap` is ``values=np.arange(0, cmap.N)``.
    3. A sequence of hex strings, color names, or RGB[A] tuples. A
       `~proplot.colors.DiscreteColormap` will be generated from these colors and
       used to build a `~matplotlib.cm.ScalarMappable` on-the-fly. The colorbar
       range and ticks depend on the arguments `values`, `norm`, and
       `norm_kw`. The default is ``values=np.arange(0, len(mappable))``.
    4. A sequence of `matplotlib.artist.Artist` instances (e.g., a list of
       `~matplotlib.lines.Line2D` instances returned by `~proplot.axes.PlotAxes.plot`).
       A colormap will be generated from the colors of these objects (where the
       color is determined by ``get_color``, if available, or ``get_facecolor``).
       The colorbar range and ticks depend on the arguments `values`, `norm`, and
       `norm_kw`. The default is to infer colorbar ticks and tick labels
       by calling `~matplotlib.artist.Artist.get_label` on each artist.

values : sequence of float or str, optional
    Ignored if `mappable` is a `~matplotlib.cm.ScalarMappable`. This maps the colormap
    colors to numeric values using `~proplot.colors.DiscreteNorm`. If the colormap is
    a `~proplot.colors.ContinuousColormap` then its colors will be "discretized".
    These These can also be strings, in which case the list indices are used for
    tick locations and the strings are applied as tick labels.
"""
_colorbar_kwargs_docstring = """
orientation : {None, 'horizontal', 'vertical'}, optional
    The colorbar orientation. By default this depends on the "side" of the subplot
    or figure where the colorbar is drawn. Inset colorbars are always horizontal.
norm : norm-spec, optional
    Ignored if `mappable` is a `~matplotlib.cm.ScalarMappable`. This is the continuous
    normalizer used to scale the `~proplot.colors.ContinuousColormap` (or passed
    to `~proplot.colors.DiscreteNorm` if `values` was passed). Passed to the
    `~proplot.constructor.Norm` constructor function.
norm_kw : dict-like, optional
    Ignored if `mappable` is a `~matplotlib.cm.ScalarMappable`. These are the
    normalizer keyword arguments. Passed to `~proplot.constructor.Norm`.
vmin, vmax : float, optional
    Ignored if `mappable` is a `~matplotlib.cm.ScalarMappable`. These are the minimum
    and maximum colorbar values. Passed to `~proplot.constructor.Norm`.
label, title : str, optional
    The colorbar label. The `title` keyword is also accepted for
    consistency with `~matplotlib.axes.Axes.legend`.
reverse : bool, optional
    Whether to reverse the direction of the colorbar. This is done automatically
    when descending levels are used with `~proplot.colors.DiscreteNorm`.
rotation : float, default: 0
    The tick label rotation.
grid, edges, drawedges : bool, default: :rc:`colorbar.grid`
    Whether to draw "grid" dividers between each distinct color.
extend : {'neither', 'both', 'min', 'max'}, optional
    Direction for drawing colorbar "extensions" (i.e. color keys for out-of-bounds
    data on the end of the colorbar). Default behavior is to use the value of `extend`
    passed to the plotting command or use ``'neither'`` if the value is unknown.
extendfrac : float, optional
    The length of the colorbar "extensions" relative to the length of the colorbar.
    This is a native matplotlib `~matplotlib.figure.Figure.colorbar` keyword.
extendsize : unit-spec, default: :rc:`colorbar.extend` or :rc:`colorbar.insetextend`
    The length of the colorbar "extensions" in physical units. Default is
    :rcraw:`colorbar.extend` for outer colorbars and :rcraw:`colorbar.insetextend`
    for inset colorbars. %(units.em)s
extendrect : bool, default: False
    Whether to draw colorbar "extensions" as rectangles. If ``False`` then
    the extensions are drawn as triangles.
locator, ticks : locator-spec, optional
    Used to determine the colorbar tick positions. Passed to the
    `~proplot.constructor.Locator` constructor function. By default
    `~matplotlib.ticker.AutoLocator` is used for continuous color levels
    and `~proplot.ticker.DiscreteLocator` is used for discrete color levels.
locator_kw : dict-like, optional
    Keyword arguments passed to `matplotlib.ticker.Locator` class.
minorlocator, minorticks
    As with `locator`, `ticks` but for the minor ticks. By default
    `~matplotlib.ticker.AutoMinorLocator` is used for continuous color levels
    and `~proplot.ticker.DiscreteLocator` is used for discrete color levels.
minorlocator_kw
    As with `locator_kw`, but for the minor ticks.
format, formatter, ticklabels : formatter-spec, optional
    The tick label format. Passed to the `~proplot.constructor.Formatter`
    constructor function.
formatter_kw : dict-like, optional
    Keyword arguments passed to `matplotlib.ticker.Formatter` class.
frame, frameon : bool, default: :rc:`colorbar.frameon`
    For inset colorbars only. Indicates whether to draw a "frame",
    just like `~matplotlib.axes.Axes.legend`.
tickminor : bool, optional
    Whether to add minor ticks using `~matplotlib.colorbar.ColorbarBase.minorticks_on`.
tickloc, ticklocation : {'bottom', 'top', 'left', 'right'}, optional
    Where to draw tick marks on the colorbar. Default is toward the outside
    of the subplot for outer colorbars and ``'bottom'`` for inset colorbars.
tickdir, tickdirection : {'out', 'in', 'inout'}, default: :rc:`tick.dir`
    Direction of major and minor colorbar ticks.
ticklen : unit-spec, default: :rc:`tick.len`
    Major tick lengths for the colorbar ticks.
ticklenratio : float, default: :rc:`tick.lenratio`
    Relative scaling of `ticklen` used to determine minor tick lengths.
tickwidth : unit-spec, default: `linewidth`
    Major tick widths for the colorbar ticks.
    or :rc:`tick.width` if `linewidth` was not passed.
tickwidthratio : float, default: :rc:`tick.widthratio`
    Relative scaling of `tickwidth` used to determine minor tick widths.
ticklabelcolor, ticklabelsize, ticklabelweight \
: default: :rc:`tick.labelcolor`, :rc:`tick.labelsize`, :rc:`tick.labelweight`.
    The font color, size, and weight for colorbar tick labels
labelloc, labellocation : {'bottom', 'top', 'left', 'right'}
    The colorbar label location. Inherits from `tickloc` by default. Default is toward
    the outside of the subplot for outer colorbars and ``'bottom'`` for inset colorbars.
labelcolor, labelsize, labelweight \
: default: :rc:`label.color`, :rc:`label.size`, and :rc:`label.weight`.
    The font color, size, and weight for the colorbar label.
a, alpha, framealpha, fc, facecolor, framecolor, ec, edgecolor, ew, edgewidth : default\
: :rc:`colorbar.framealpha`, :rc:`colorbar.framecolor`
    For inset colorbars only. Controls the transparency and color of
    the background frame.
lw, linewidth, c, color : optional
    Controls the line width and edge color for both the colorbar
    outline and the level dividers.
%(axes.edgefix)s
rasterize : bool, default: :rc:`colorbar.rasterize`
    Whether to rasterize the colorbar solids. The matplotlib default was ``True``
    but proplot changes this to ``False`` since rasterization can cause misalignment
    between the color patches and the colorbar outline.
**kwargs
    Passed to `~matplotlib.figure.Figure.colorbar`.
"""
_edgefix_docstring = """
edgefix : bool or float, default: :rc:`edgefix`
    Whether to fix the common issue where white lines appear between adjacent
    patches in saved vector graphics (this can slow down figure rendering).
    See this `github repo <https://github.com/jklymak/contourfIssues>`__ for a
    demonstration of the problem. If ``True``, a small default linewidth of
    ``0.3`` is used to cover up the white lines. If float (e.g. ``edgefix=0.5``),
    this specific linewidth is used to cover up the white lines. This feature is
    automatically disabled when the patches have transparency.
"""
docstring._snippet_manager['axes.edgefix'] = _edgefix_docstring
docstring._snippet_manager['axes.colorbar_args'] = _colorbar_args_docstring
docstring._snippet_manager['axes.colorbar_kwargs'] = _colorbar_kwargs_docstring


# Legend docstrings
_legend_args_docstring = """
handles : list of artist, optional
    List of matplotlib artists, or a list of lists of artist instances (see the `center`
    keyword). If not passed, artists with valid labels (applied by passing `label` or
    `labels` to a plotting command or calling `~matplotlib.artist.Artist.set_label`)
    are retrieved automatically. If the object is a `~matplotlib.contour.ContourSet`,
    `~matplotlib.contour.ContourSet.legend_elements` is used to select the central
    artist in the list (generally useful for single-color contour plots). Note that
    proplot's `~proplot.axes.PlotAxes.contour` and `~proplot.axes.PlotAxes.contourf`
    accept a legend `label` keyword argument.
labels : list of str, optional
    A matching list of string labels or ``None`` placeholders, or a matching list of
    lists (see the `center` keyword). Wherever ``None`` appears in the list (or
    if no labels were passed at all), labels are retrieved by calling
    `~matplotlib.artist.Artist.get_label` on each `~matplotlib.artist.Artist` in the
    handle list. If a handle consists of a tuple group of artists, labels are inferred
    from the artists in the tuple (if there are multiple unique labels in the tuple
    group of artists, the tuple group is expanded into unique legend entries --
    otherwise, the tuple group elements are drawn on top of eachother). For details
    on matplotlib legend handlers and tuple groups, see the matplotlib `legend guide \
<https://matplotlib.org/stable/tutorials/intermediate/legend_guide.html>`__.
"""
_legend_kwargs_docstring = """
frame, frameon : bool, optional
    Toggles the legend frame. For centered-row legends, a frame
    independent from matplotlib's built-in legend frame is created.
ncol, ncols : int, optional
    The number of columns. `ncols` is an alias, added
    for consistency with `~matplotlib.pyplot.subplots`.
order : {'C', 'F'}, optional
    Whether legend handles are drawn in row-major (``'C'``) or column-major
    (``'F'``) order. Analagous to `numpy.array` ordering. The matplotlib
    default was ``'F'`` but proplot changes this to ``'C'``.
center : bool, optional
    Whether to center each legend row individually. If ``True``, we draw
    successive single-row legends "stacked" on top of each other. If ``None``,
    we infer this setting from `handles`. By default, `center` is set to ``True``
    if `handles` is a list of lists (each sublist is used as a row in the legend).
alphabetize : bool, default: False
    Whether to alphabetize the legend entries according to
    the legend labels.
title, label : str, optional
    The legend title. The `label` keyword is also accepted, for consistency
    with `~matplotlib.figure.Figure.colorbar`.
fontsize, fontweight, fontcolor : optional
    The font size, weight, and color for the legend text. Font size is interpreted
    by `~proplot.utils.units`. The default font size is :rcraw:`legend.fontsize`.
titlefontsize, titlefontweight, titlefontcolor : optional
    The font size, weight, and color for the legend title. Font size is interpreted
    by `~proplot.utils.units`. The default size is `fontsize`.
borderpad, borderaxespad, handlelength, handleheight, handletextpad, \
labelspacing, columnspacing : unit-spec, optional
    Various matplotlib `~matplotlib.axes.Axes.legend` spacing arguments.
    %(units.em)s
a, alpha, framealpha, fc, facecolor, framecolor, ec, edgecolor, ew, edgewidth \
: default: :rc:`legend.framealpha`, :rc:`legend.facecolor`, :rc:`legend.edgecolor`, \
:rc:`axes.linewidth`
    The opacity, face color, edge color, and edge width for the legend frame.
c, color, lw, linewidth, m, marker, ls, linestyle, dashes, ms, markersize : optional
    Properties used to override the legend handles. For example, for a
    legend describing variations in line style ignoring variations
    in color, you might want to use ``color='black'``.
handle_kw : dict-like, optional
    Additional properties used to override legend handles, e.g.
    ``handle_kw={'edgecolor': 'black'}``. Only line properties
    can be passed as keyword arguments.
handler_map : dict-like, optional
    A dictionary mapping instances or types to a legend handler.
    This `handler_map` updates the default handler map found at
    `matplotlib.legend.Legend.get_legend_handler_map`.
**kwargs
    Passed to `~matplotlib.axes.Axes.legend`.
"""
docstring._snippet_manager['axes.legend_args'] = _legend_args_docstring
docstring._snippet_manager['axes.legend_kwargs'] = _legend_kwargs_docstring


def _align_bbox(align, length):
    """
    Return a simple alignment bounding box for intersection calculations.
    """
    if align in ('left', 'bottom'):
        bounds = [[0, 0], [length, 0]]
    elif align in ('top', 'right'):
        bounds = [[1 - length, 0], [1, 0]]
    elif align == 'center':
        bounds = [[0.5 * (1 - length), 0], [0.5 * (1 + length), 0]]
    else:
        raise ValueError(f'Invalid align {align!r}.')
    return mtransforms.Bbox(bounds)


class _TransformedBoundsLocator:
    """
    Axes locator for `~Axes.inset_axes` and other axes.
    """
    def __init__(self, bounds, transform):
        self._bounds = bounds
        self._transform = transform

    def __call__(self, ax, renderer):  # noqa: U100
        transfig = getattr(ax.figure, 'transSubfigure', ax.figure.transFigure)
        bbox = mtransforms.Bbox.from_bounds(*self._bounds)
        bbox = mtransforms.TransformedBbox(bbox, self._transform)
        bbox = mtransforms.TransformedBbox(bbox, transfig.inverted())
        return bbox


class Axes(maxes.Axes):
    """
    The lowest-level `~matplotlib.axes.Axes` subclass used by proplot.
    Implements basic universal features.
    """
    _name = None  # derived must override
    _name_aliases = ()
    _make_inset_locator = _TransformedBoundsLocator

    def __repr__(self):
        # Show the position in the geometry excluding panels. Panels are
        # indicated by showing their parent geometry plus a 'side' argument.
        # WARNING: This will not be used in matplotlib 3.3.0 (and probably next
        # minor releases) because native __repr__ is defined in SubplotBase.
        ax = self._get_topmost_axes()
        name = type(self).__name__
        prefix = '' if ax is self else 'parent_'
        params = {}
        if self._name in ('cartopy', 'basemap'):
            name = name.replace('_' + self._name.title(), 'Geo')
            params['backend'] = self._name
        if self._inset_parent:
            name = re.sub('Axes(Subplot)?', 'AxesInset', name)
            params['bounds'] = tuple(np.round(self._inset_bounds, 2))
        if self._altx_parent or self._alty_parent:
            name = re.sub('Axes(Subplot)?', 'AxesTwin', name)
            params['axis'] = 'x' if self._altx_parent else 'y'
        if self._colorbar_fill:
            name = re.sub('Axes(Subplot)?', 'AxesFill', name)
            params['side'] = self._axes._panel_side
        if self._panel_side:
            name = re.sub('Axes(Subplot)?', 'AxesPanel', name)
            params['side'] = self._panel_side
        try:
            nrows, ncols, num1, num2 = ax.get_subplotspec()._get_geometry()
            params[prefix + 'index'] = (num1, num2)
        except (IndexError, ValueError, AttributeError):  # e.g. a loose axes
            left, bottom, width, height = np.round(self._position.bounds, 2)
            params['left'], params['bottom'], params['size'] = (left, bottom, (width, bottom))  # noqa: E501
        if ax.number:
            params[prefix + 'number'] = ax.number
        params = ', '.join(f'{key}={value!r}' for key, value in params.items())
        return f'{name}({params})'

    def __str__(self):
        return self.__repr__()

    @docstring._snippet_manager
    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        *args
            Passed to `matplotlib.axes.Axes`.
        %(axes.format)s

        Other parameters
        ----------------
        %(rc.init)s

        See also
        --------
        Axes.format
        matplotlib.axes.Axes
        proplot.axes.PlotAxes
        proplot.axes.CartesianAxes
        proplot.axes.PolarAxes
        proplot.axes.GeoAxes
        proplot.figure.Figure.subplot
        proplot.figure.Figure.add_subplot
        """
        # Remove subplot-related args
        # NOTE: These are documented on add_subplot()
        ss = kwargs.pop('_subplot_spec', None)  # see below
        number = kwargs.pop('number', None)
        autoshare = kwargs.pop('autoshare', None)
        autoshare = _not_none(autoshare, True)

        # Remove format-related args and initialize
        rc_kw, rc_mode = _pop_rc(kwargs)
        kw_format = _pop_props(kwargs, 'patch')  # background properties
        if 'zorder' in kw_format:  # special case: refers to the entire axes
            kwargs['zorder'] = kw_format.pop('zorder')
        for cls, sig in self._format_signatures.items():
            if isinstance(self, cls):
                kw_format.update(_pop_params(kwargs, sig))
        super().__init__(*args, **kwargs)

        # Varous scalar properties
        self._active_cycle = rc['axes.prop_cycle']
        self._auto_format = None  # manipulated by wrapper functions
        self._abc_border_kwargs = {}
        self._abc_loc = None
        self._abc_title_pad = rc['abc.titlepad']
        self._title_above = rc['title.above']
        self._title_border_kwargs = {}  # title border properties
        self._title_loc = None
        self._title_pad = rc['title.pad']
        self._title_pad_current = None
        self._altx_parent = None  # for cartesian axes only
        self._alty_parent = None
        self._colorbar_fill = None
        self._inset_parent = None
        self._inset_bounds = None  # for introspection ony
        self._inset_zoom = False
        self._inset_zoom_artists = None
        self._panel_hidden = False  # True when "filled" with cbar/legend
        self._panel_align = {}  # store 'align' and 'length' for "filled" cbar/legend
        self._panel_parent = None
        self._panel_share = False
        self._panel_sharex_group = False
        self._panel_sharey_group = False
        self._panel_side = None
        self._tight_bbox = None  # bounding boxes are saved
        self.xaxis.isDefault_minloc = True  # ensure enabled at start (needed for dual)
        self.yaxis.isDefault_minloc = True

        # Various dictionary properties
        # NOTE: Critical to use self.text() so they are patched with _update_label
        self._legend_dict = {}
        self._colorbar_dict = {}
        d = self._panel_dict = {}
        d['left'] = []  # NOTE: panels will be sorted inside-to-outside
        d['right'] = []
        d['bottom'] = []
        d['top'] = []
        d = self._title_dict = {}
        kw = {'zorder': 3.5, 'transform': self.transAxes}
        d['abc'] = self.text(0, 0, '', **kw)
        d['left'] = self._left_title  # WARNING: track in case mpl changes this
        d['center'] = self.title
        d['right'] = self._right_title
        d['upper left'] = self.text(0, 0, '', va='top', ha='left', **kw)
        d['upper center'] = self.text(0, 0.5, '', va='top', ha='center', **kw)
        d['upper right'] = self.text(0, 1, '', va='top', ha='right', **kw)
        d['lower left'] = self.text(0, 0, '', va='bottom', ha='left', **kw)
        d['lower center'] = self.text(0, 0.5, '', va='bottom', ha='center', **kw)
        d['lower right'] = self.text(0, 1, '', va='bottom', ha='right', **kw)

        # Subplot-specific settings
        # NOTE: Default number for any axes is None (i.e., no a-b-c labels allowed)
        # and for subplots added with add_subplot is incremented automatically
        # WARNING: For mpl>=3.4.0 subplotspec assigned *after* initialization using
        # set_subplotspec. Tried to defer to setter but really messes up both format()
        # and _apply_auto_share(). Instead use workaround: Have Figure.add_subplot pass
        # subplotspec as a hidden keyword arg. Non-subplots don't need this arg.
        # See: https://github.com/matplotlib/matplotlib/pull/18564
        self._number = None
        if number:  # not None or False
            self.number = number
        if ss is not None:  # always passed from add_subplot
            self.set_subplotspec(ss)
        if autoshare:
            self._apply_auto_share()

        # Default formatting
        # NOTE: This ignores user-input rc_mode. Mode '1' applies proplot
        # features which is necessary on first run. Default otherwise is mode '2'
        self.format(rc_kw=rc_kw, rc_mode=1, skip_figure=True, **kw_format)

    def _add_inset_axes(
        self, bounds, transform=None, *, proj=None, projection=None,
        zoom=None, zoom_kw=None, zorder=None, **kwargs
    ):
        """
        Add an inset axes using arbitrary projection.
        """
        # Converting transform to figure-relative coordinates
        transform = self._get_transform(transform, 'axes')
        locator = self._make_inset_locator(bounds, transform)
        bounds = locator(self, None).bounds
        label = kwargs.pop('label', 'inset_axes')
        zorder = _not_none(zorder, 4)

        # Parse projection and inherit from the current axes by default
        # NOTE: The _parse_proj method also accepts axes classes.
        proj = _not_none(proj=proj, projection=projection)
        if proj is None:
            if self._name in ('cartopy', 'basemap'):
                proj = copy.copy(self.projection)
            else:
                proj = self._name
        kwargs = self.figure._parse_proj(proj, **kwargs)

        # Create axes and apply locator. The locator lets the axes adjust
        # automatically if we used data coords. Called by ax.apply_aspect()
        cls = mprojections.get_projection_class(kwargs.pop('projection'))
        ax = cls(self.figure, bounds, zorder=zorder, label=label, **kwargs)
        ax.set_axes_locator(locator)
        ax._inset_parent = self
        ax._inset_bounds = bounds
        self.add_child_axes(ax)

        # Add zoom indicator (NOTE: requires matplotlib >= 3.0)
        zoom = _not_none(zoom, self._name == 'cartesian' and ax._name == 'cartesian')
        ax._inset_zoom = zoom
        if zoom:
            zoom_kw = zoom_kw or {}
            ax.indicate_inset_zoom(**zoom_kw)
        return ax

    def _add_queued_guides(self):
        """
        Draw the queued-up legends and colorbars. Wrapper funcs and legend func let
        user add handles to location lists with successive calls.
        """
        # Draw queued colorbars
        for (loc, align), colorbar in tuple(self._colorbar_dict.items()):
            if not isinstance(colorbar, tuple):
                continue
            handles, labels, kwargs = colorbar
            cb = self._add_colorbar(handles, labels, loc=loc, align=align, **kwargs)
            self._colorbar_dict[(loc, align)] = cb

        # Draw queued legends
        # WARNING: Passing empty list labels=[] to legend causes matplotlib
        # _parse_legend_args to search for everything. Ensure None if empty.
        for (loc, align), legend in tuple(self._legend_dict.items()):
            if not isinstance(legend, tuple) or any(isinstance(_, mlegend.Legend) for _ in legend):  # noqa: E501
                continue
            handles, labels, kwargs = legend
            leg = self._add_legend(handles, labels, loc=loc, align=align, **kwargs)
            self._legend_dict[(loc, align)] = leg

    def _add_guide_frame(
        self, xmin, ymin, width, height, *, fontsize, fancybox=None, **kwargs
    ):
        """
        Add a colorbar or multilegend frame.
        """
        # TODO: Shadow patch does not seem to work. Unsure why.
        # TODO: Add basic 'colorbar' and 'legend' artists with
        # shared control over background frame.
        shadow = kwargs.pop('shadow', None)  # noqa: F841
        renderer = self.figure._get_renderer()
        fontsize = _fontsize_to_pt(fontsize)
        fontsize = (fontsize / 72) / self._get_size_inches()[0]  # axes relative units
        fontsize = renderer.points_to_pixels(fontsize)
        patch = mpatches.FancyBboxPatch(
            (xmin, ymin), width, height,
            snap=True,
            zorder=4.5,
            mutation_scale=fontsize,
            transform=self.transAxes
        )
        patch.set_clip_on(False)
        if fancybox:
            patch.set_boxstyle('round', pad=0, rounding_size=0.2)
        else:
            patch.set_boxstyle('square', pad=0)
        patch.update(kwargs)
        self.add_artist(patch)
        return patch

    def _add_guide_panel(self, loc='fill', align='center', length=0, **kwargs):
        """
        Add a panel to be filled by an "outer" colorbar or legend.
        """
        # NOTE: For colorbars we include 'length' when determining whether to allocate
        # new panel but for legend just test whether that 'align' position was filled.
        # WARNING: Hide content but 1) do not use ax.set_visible(False) so that
        # tight layout will include legend and colorbar and 2) do not use
        # ax.clear() so that top panel title and a-b-c label can remain.
        bbox = _align_bbox(align, length)
        if loc == 'fill':
            ax = self
        elif loc in ('left', 'right', 'top', 'bottom'):
            ax = None
            for pax in self._panel_dict[loc]:
                if not pax._panel_hidden or align in pax._panel_align:
                    continue
                if not any(bbox.overlaps(b) for b in pax._panel_align.values()):
                    ax = pax
                    break
            if ax is None:
                ax = self.panel_axes(loc, filled=True, **kwargs)
        else:
            raise ValueError(f'Invalid filled panel location {loc!r}.')
        for s in ax.spines.values():
            s.set_visible(False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.patch.set_facecolor('none')
        ax._panel_hidden = True
        ax._panel_align[align] = bbox
        return ax

    @warnings._rename_kwargs('0.10', rasterize='rasterized')
    def _add_colorbar(
        self, mappable, values=None, *,
        loc=None, align=None, space=None, pad=None,
        width=None, length=None, shrink=None,
        label=None, title=None, reverse=False,
        rotation=None, grid=None, edges=None, drawedges=None,
        extend=None, extendsize=None, extendfrac=None,
        ticks=None, locator=None, locator_kw=None,
        format=None, formatter=None, ticklabels=None, formatter_kw=None,
        minorticks=None, minorlocator=None, minorlocator_kw=None,
        tickminor=None, ticklen=None, ticklenratio=None,
        tickdir=None, tickdirection=None, tickwidth=None, tickwidthratio=None,
        ticklabelsize=None, ticklabelweight=None, ticklabelcolor=None,
        labelloc=None, labellocation=None, labelsize=None, labelweight=None,
        labelcolor=None, c=None, color=None, lw=None, linewidth=None,
        edgefix=None, rasterized=None, **kwargs
    ):
        """
        The driver function for adding axes colorbars.
        """
        # Parse input args
        # TODO: Get the 'best' inset colorbar location using the legend algorithm
        # and implement inset colorbars the same as inset legends.
        grid = _not_none(grid=grid, edges=edges, drawedges=drawedges, default=rc['colorbar.grid'])  # noqa: E501
        length = _not_none(length=length, shrink=shrink)
        label = _not_none(title=title, label=label)
        labelloc = _not_none(labelloc=labelloc, labellocation=labellocation)
        locator = _not_none(ticks=ticks, locator=locator)
        formatter = _not_none(ticklabels=ticklabels, formatter=formatter, format=format)
        minorlocator = _not_none(minorticks=minorticks, minorlocator=minorlocator)
        color = _not_none(c=c, color=color, default=rc['axes.edgecolor'])
        linewidth = _not_none(lw=lw, linewidth=linewidth)
        ticklen = units(_not_none(ticklen, rc['tick.len']), 'pt')
        tickdir = _not_none(tickdir=tickdir, tickdirection=tickdirection)
        tickwidth = units(_not_none(tickwidth, linewidth, rc['tick.width']), 'pt')
        linewidth = units(_not_none(linewidth, default=rc['axes.linewidth']), 'pt')
        ticklenratio = _not_none(ticklenratio, rc['tick.lenratio'])
        tickwidthratio = _not_none(tickwidthratio, rc['tick.widthratio'])
        rasterized = _not_none(rasterized, rc['colorbar.rasterized'])
        locator_kw = locator_kw or {}
        formatter_kw = formatter_kw or {}
        minorlocator_kw = minorlocator_kw or {}
        for b, kw in enumerate((locator_kw, minorlocator_kw)):
            key = 'maxn_minor' if b else 'maxn'
            name = 'minorlocator' if b else 'locator'
            nbins = kwargs.pop('maxn_minor' if b else 'maxn', None)
            if nbins is not None:
                kw['nbins'] = nbins
                warnings._warn_proplot(
                    f'The colorbar() keyword {key!r} was deprecated in v0.10. To '
                    "achieve the same effect, you can pass 'nbins' to the new default "
                    f"locator DiscreteLocator using {name}_kw={{'nbins': {nbins}}}."
                )
        kw_label = {}
        for key, value in (
            ('size', labelsize),
            ('weight', labelweight),
            ('color', labelcolor),
        ):
            if value is not None:
                kw_label[key] = value
        kw_ticklabels = {}
        for key, value in (
            ('size', ticklabelsize),
            ('weight', ticklabelweight),
            ('color', ticklabelcolor),
            ('rotation', rotation),
        ):
            if value is not None:
                kw_ticklabels[key] = value

        # Generate and prepare the colorbar axes
        # NOTE: The inset axes function needs 'label' to know how to pad the box
        # TODO: Use seperate keywords for frame properties vs. colorbar edge properties?
        if loc in ('fill', 'left', 'right', 'top', 'bottom'):
            length = _not_none(length, rc['colorbar.length'])  # for _add_guide_panel
            kwargs.update({'align': align, 'length': length})
            extendsize = _not_none(extendsize, rc['colorbar.extend'])
            ax = self._add_guide_panel(loc, align, length=length, width=width, space=space, pad=pad)  # noqa: E501
            cax, kwargs = ax._parse_colorbar_filled(**kwargs)
        else:
            kwargs.update({'label': label, 'length': length, 'width': width})
            extendsize = _not_none(extendsize, rc['colorbar.insetextend'])
            cax, kwargs = self._parse_colorbar_inset(loc=loc, pad=pad, **kwargs)  # noqa: E501

        # Parse the colorbar mappable
        # NOTE: Account for special case where auto colorbar is generated from 1D
        # methods that construct an 'artist list' (i.e. colormap scatter object)
        if np.iterable(mappable) and len(mappable) == 1 and isinstance(mappable[0], mcm.ScalarMappable):  # noqa: E501
            mappable = mappable[0]
        if isinstance(mappable, mcm.ScalarMappable):
            pop = _pop_params(kwargs, cax._parse_colorbar_arg, ignore_internal=True)
            locator_default = formatter_default = None
            if pop:
                warnings._warn_proplot(f'Input is already a ScalarMappable. Ignoring unused keyword arg(s): {pop}')  # noqa: E501
        else:
            result = cax._parse_colorbar_arg(mappable, values, **kwargs)
            mappable, locator_default, formatter_default, kwargs = result

        # Parse ticking keyword arguments
        # NOTE: This uses DiscreteLocator for default discrete minor ticks
        name = 'y' if kwargs.get('orientation') == 'vertical' else 'x'
        axis = cax.yaxis if kwargs.get('orientation') == 'vertical' else cax.xaxis
        locator = _not_none(locator, locator_default, None)
        formatter = _not_none(formatter, formatter_default, 'auto')
        formatter = constructor.Formatter(formatter, **formatter_kw)
        categorical = isinstance(formatter, mticker.FixedFormatter)
        discrete = isinstance(mappable.norm, pcolors.DiscreteNorm)
        if locator is not None:
            locator = constructor.Locator(locator, **locator_kw)
        if minorlocator is None and isinstance(locator, pticker.DiscreteLocator):
            if categorical:  # never add default minor ticks
                pass
            elif tickminor or tickminor is None:
                minorlocator = pticker.DiscreteLocator(list(locator.locs), minor=True)
        if tickminor is None:
            if discrete or categorical:  # never use the default minor locator
                tickminor = False
            else:
                tickminor = rc[name + 'tick.minor.visible']
        if minorlocator is not None:
            minorlocator = constructor.Locator(minorlocator, **minorlocator_kw)

        # Prepare colorbar keyword arguments
        # WARNING: Critical to not pass empty major locators in matplotlib < 3.5
        # See: https://github.com/lukelbd/proplot/issues/301
        if isinstance(locator, mticker.NullLocator) or not len(getattr(locator, 'locs', (None,))):  # noqa: E501
            minorlocator, tickminor = None, False  # attempted fix
        for ticker in (locator, formatter, minorlocator):
            if isinstance(ticker, mticker.TickHelper):
                ticker.set_axis(axis)
        if extendsize is not None and extendfrac is not None:
            warnings._warn_proplot(
                f'You cannot specify both an absolute extendsize={extendsize!r} '
                f"and a relative extendfrac={extendfrac!r}. Ignoring 'extendfrac'."
            )
            extendfrac = None
        if extendfrac is None:
            width, height = cax._get_size_inches()
            scale = height if kwargs.get('orientation') == 'vertical' else width
            extendsize = units(extendsize, 'em', 'in')
            extendfrac = extendsize / max(scale - 2 * extendsize, units(1, 'em', 'in'))
        kwargs.update(
            {
                'cax': cax,
                'ticks': locator,
                'format': formatter,
                'extendfrac': extendfrac,
                'drawedges': grid,
            }
        )
        kwargs.setdefault('spacing', 'uniform')
        extend = _not_none(extend, 'neither')
        if isinstance(mappable, mcontour.ContourSet):
            mappable.extend = extend  # required in mpl >= 3.3, else optional
        else:
            kwargs['extend'] = extend

        # Draw and update the colorbar
        # WARNING: Must use colorbar set_label to set text,
        # calling set_text on the axis will do nothing!
        # WARNING: Colorbar _ticker() internally makes dummy axis and updates view
        # limits. Here we apply actual axis rather than dummy, otherwise default nbins
        # of DiscreteLocator will not work. Not sure if this has side effects...
        obj = cax._colorbar_fill = cax.figure.colorbar(mappable, **kwargs)
        obj.minorlocator = minorlocator  # backwards compatibility
        obj.update_ticks = guides._update_ticks.__get__(obj)  # backwards compatibility
        if minorlocator is not None:
            obj.update_ticks()
        elif tickminor:
            obj.minorticks_on()
        else:
            obj.minorticks_off()
        axis.set_tick_params(which='both', color=color, direction=tickdir)
        axis.set_tick_params(which='major', length=ticklen, width=tickwidth)
        axis.set_tick_params(which='minor', length=ticklen * ticklenratio, width=tickwidth * tickwidthratio)  # noqa: E501
        if getattr(mappable.norm, 'descending', None):
            axis.set_inverted(True)
        if reverse:  # potentially double reverse, although that would be weird...
            axis.set_inverted(True)
        if label is not None:
            obj.set_label(label)
        if labelloc is not None:
            axis.set_label_position(labelloc)
        axis.label.update(kw_label)
        for label in axis.get_ticklabels():
            label.update(kw_ticklabels)
        kw_outline = {'edgecolor': color, 'linewidth': linewidth}
        if obj.outline is not None:
            obj.outline.update(kw_outline)
        if obj.dividers is not None:
            obj.dividers.update(kw_outline)
        if obj.solids:
            obj.solids.set_rasterized(rasterized)
            cax._fix_patch_edges(obj.solids, edgefix=edgefix)

        # Return after registering location
        self._register_guide('colorbar', obj, (loc, align))  # possibly replace another
        return obj

    def _add_legend(
        self, handles=None, labels=None, *,
        loc=None, align=None, width=None, pad=None, space=None,
        frame=None, frameon=None, ncol=None, ncols=None,
        alphabetize=False, center=None, order=None, label=None, title=None,
        fontsize=None, fontweight=None, fontcolor=None,
        titlefontsize=None, titlefontweight=None, titlefontcolor=None,
        handle_kw=None, handler_map=None, **kwargs
    ):
        """
        The driver function for adding axes legends.
        """
        # Parse input argument units
        ncol = _not_none(ncols=ncols, ncol=ncol)
        order = _not_none(order, 'C')
        frameon = _not_none(frame=frame, frameon=frameon, default=rc['legend.frameon'])
        fontsize = _not_none(kwargs.pop('fontsize', None), rc['legend.fontsize'])
        titlefontsize = _not_none(
            title_fontsize=kwargs.pop('title_fontsize', None),
            titlefontsize=titlefontsize,
            default=rc['legend.title_fontsize']
        )
        fontsize = _fontsize_to_pt(fontsize)
        titlefontsize = _fontsize_to_pt(titlefontsize)
        if order not in ('F', 'C'):
            raise ValueError(
                f'Invalid order {order!r}. Please choose from '
                "'C' (row-major, default) or 'F' (column-major)."
            )

        # Convert relevant keys to em-widths
        for setting in rcsetup.EM_KEYS:  # em-width keys
            pair = setting.split('legend.', 1)
            if len(pair) == 1:
                continue
            _, key = pair
            value = kwargs.pop(key, None)
            if isinstance(value, str):
                value = units(kwargs[key], 'em', fontsize=fontsize)
            if value is not None:
                kwargs[key] = value

        # Generate and prepare the legend axes
        if loc in ('fill', 'left', 'right', 'top', 'bottom'):
            lax = self._add_guide_panel(loc, align, width=width, space=space, pad=pad)
            kwargs.setdefault('borderaxespad', 0)
            if not frameon:
                kwargs.setdefault('borderpad', 0)
            try:
                kwargs['loc'] = ALIGN_OPTS[lax._panel_side][align]
            except KeyError:
                raise ValueError(f'Invalid align={align!r} for legend loc={loc!r}.')
        else:
            lax = self
            pad = kwargs.pop('borderaxespad', pad)
            kwargs['loc'] = loc  # simply pass to legend
            kwargs['borderaxespad'] = units(pad, 'em', fontsize=fontsize)

        # Handle and text properties that are applied after-the-fact
        # NOTE: Set solid_capstyle to 'butt' so line does not extend past error bounds
        # shading in legend entry. This change is not noticable in other situations.
        kw_frame, kwargs = lax._parse_frame('legend', **kwargs)
        kw_text = {}
        if fontcolor is not None:
            kw_text['color'] = fontcolor
        if fontweight is not None:
            kw_text['weight'] = fontweight
        kw_title = {}
        if titlefontcolor is not None:
            kw_title['color'] = titlefontcolor
        if titlefontweight is not None:
            kw_title['weight'] = titlefontweight
        kw_handle = _pop_props(kwargs, 'line')
        kw_handle.setdefault('solid_capstyle', 'butt')
        kw_handle.update(handle_kw or {})

        # Parse the legend arguments using axes for auto-handle detection
        # TODO: Update this when we no longer use "filled panels" for outer legends
        pairs, multi = lax._parse_legend_handles(
            handles, labels, ncol=ncol, order=order, center=center,
            alphabetize=alphabetize, handler_map=handler_map
        )
        title = _not_none(label=label, title=title)
        kwargs.update(
            {
                'title': title,
                'frameon': frameon,
                'fontsize': fontsize,
                'handler_map': handler_map,
                'title_fontsize': titlefontsize,
            }
        )

        # Add the legend and update patch properties
        if multi:
            objs = lax._parse_legend_centered(pairs, kw_frame=kw_frame, **kwargs)
        else:
            kwargs.update({key: kw_frame.pop(key) for key in ('shadow', 'fancybox')})
            objs = [lax._parse_legend_aligned(pairs, ncol=ncol, order=order, **kwargs)]
            objs[0].legendPatch.update(kw_frame)
        for obj in objs:
            if hasattr(lax, 'legend_') and lax.legend_ is None:
                lax.legend_ = obj  # make first legend accessible with get_legend()
            else:
                lax.add_artist(obj)

        # Update legend patch and elements
        # TODO: Add capacity for categorical labels in a single legend like seaborn
        # rather than manual handle overrides with multiple legends.
        # WARNING: legendHandles only contains the *first* artist per legend because
        # HandlerBase.legend_artist() called in Legend._init_legend_box() only
        # returns the first artist. Instead we try to iterate through offset boxes.
        for obj in objs:
            box = getattr(obj, '_legend_handle_box', None)
            for obj in guides._iter_children(box):
                if isinstance(obj, mtext.Text):
                    kw = kw_text
                else:
                    kw = {key: val for key, val in kw_handle.items() if hasattr(obj, 'set_' + key)}  # noqa: E501
                    if hasattr(obj, 'set_sizes') and 'markersize' in kw_handle:
                        kw['sizes'] = np.atleast_1d(kw_handle['markersize'])
                obj.update(kw)

        # Return after registering location
        for obj in objs:
            obj.set_clip_on(False)  # critical for tight bounding box calcs
        if isinstance(objs[0], mpatches.FancyBboxPatch):
            objs = objs[1:]
        obj = objs[0] if len(objs) == 1 else tuple(objs)
        self._register_guide('legend', obj, (loc, align))  # possibly replace another
        return obj

    def _apply_title_above(self):
        """
        Change assignment of outer titles between main subplot and upper panels.
        This is called when a panel is created or `_update_title` is called.
        """
        # NOTE: Similar to how _apply_axis_sharing() is called in _align_axis_labels()
        # this is called in _align_super_labels() so we get the correct offset.
        paxs = self._panel_dict['top']
        if not paxs:
            return
        pax = paxs[-1]
        names = ('left', 'center', 'right')
        if self._abc_loc in names:
            names += ('abc',)
        if not self._title_above:
            return
        if pax._panel_hidden and self._title_above == 'panels':
            return
        pax._title_pad = self._title_pad
        pax._abc_title_pad = self._abc_title_pad
        for name in names:
            labels._transfer_label(self._title_dict[name], pax._title_dict[name])

    def _apply_auto_share(self):
        """
        Automatically configure axis sharing based on the horizontal and
        vertical extent of subplots in the figure gridspec.
        """
        # Panel axes sharing, between main subplot and its panels
        # NOTE: _panel_share means "include this panel in the axis sharing group"
        # while _panel_sharex_group indicates the group itself and may include main axes
        def shared(paxs):
            return [pax for pax in paxs if not pax._panel_hidden and pax._panel_share]

        # Internal axis sharing, share stacks of panels and main axes with each other
        # NOTE: This is called on the main axes whenver a panel is created.
        # NOTE: This block is why, even though we have figure-wide share[xy], we
        # still need the axes-specific _share[xy]_override attribute.
        if not self._panel_side:  # this is a main axes
            # Top and bottom
            bottom = self
            paxs = shared(self._panel_dict['bottom'])
            if paxs:
                bottom = paxs[-1]
                bottom._panel_sharex_group = False
                for iax in (self, *paxs[:-1]):
                    iax._panel_sharex_group = True
                    iax._sharex_setup(bottom)  # parent is bottom-most
            paxs = shared(self._panel_dict['top'])
            for iax in paxs:
                iax._panel_sharex_group = True
                iax._sharex_setup(bottom)
            # Left and right
            # NOTE: Order of panel lists is always inside-to-outside
            left = self
            paxs = shared(self._panel_dict['left'])
            if paxs:
                left = paxs[-1]
                left._panel_sharey_group = False
                for iax in (self, *paxs[:-1]):
                    iax._panel_sharey_group = True
                    iax._sharey_setup(left)  # parent is left-most
            paxs = shared(self._panel_dict['right'])
            for iax in paxs:
                iax._panel_sharey_group = True
                iax._sharey_setup(left)

        # External axes sharing, sometimes overrides panel axes sharing
        # Share x axes
        parent, *children = self._get_share_axes('x')
        for child in children:
            child._sharex_setup(parent)
        # Share y axes
        parent, *children = self._get_share_axes('y')
        for child in children:
            child._sharey_setup(parent)
        # Global sharing, use the reference subplot because why not
        ref = self.figure._subplot_dict.get(self.figure._refnum, None)
        if self is not ref:
            if self.figure._sharex > 3:
                self._sharex_setup(ref, labels=False)
            if self.figure._sharey > 3:
                self._sharey_setup(ref, labels=False)

    def _artist_fully_clipped(self, artist):
        """
        Return a boolean flag, ``True`` if the artist is clipped to the axes
        and can thus be skipped in layout calculations.
        """
        clip_box = artist.get_clip_box()
        clip_path = artist.get_clip_path()
        types_noclip = (
            maxes.Axes, maxis.Axis, moffsetbox.AnnotationBbox, moffsetbox.OffsetBox
        )
        return not isinstance(artist, types_noclip) and (
            artist.get_clip_on()
            and (clip_box is not None or clip_path is not None)
            and (
                clip_box is None
                or np.all(clip_box.extents == self.bbox.extents)
            )
            and (
                clip_path is None
                or isinstance(clip_path, mtransforms.TransformedPatchPath)
                and clip_path._patch is self.patch
            )
        )

    def _get_share_axes(self, sx, panels=False):
        """
        Return the axes whose horizontal or vertical extent in the main gridspec
        matches the horizontal or vertical extent of this axes.
        """
        # NOTE: The lefmost or bottommost axes are at the start of the list.
        if not isinstance(self, maxes.SubplotBase):
            return [self]
        i = 0 if sx == 'x' else 1
        sy = 'y' if sx == 'x' else 'x'
        argfunc = np.argmax if sx == 'x' else np.argmin
        irange = self._range_subplotspec(sx)
        axs = self.figure._iter_axes(hidden=False, children=False, panels=panels)
        axs = [ax for ax in axs if ax._range_subplotspec(sx) == irange]
        axs = list({self, *axs})  # self may be missing during initialization
        pax = axs.pop(argfunc([ax._range_subplotspec(sy)[i] for ax in axs]))
        return [pax, *axs]  # return with leftmost or bottommost first

    def _get_span_axes(self, side, panels=False):
        """
        Return the axes whose left, right, top, or bottom sides abutt against
        the same row or column as this axes. Deflect to shared panels.
        """
        if side not in ('left', 'right', 'bottom', 'top'):
            raise ValueError(f'Invalid side {side!r}.')
        if not isinstance(self, maxes.SubplotBase):
            return [self]
        x, y = 'xy' if side in ('left', 'right') else 'yx'
        idx = 0 if side in ('left', 'top') else 1  # which side to test
        coord = self._range_subplotspec(x)[idx]  # side for a particular axes
        axs = self.figure._iter_axes(hidden=False, children=False, panels=panels)
        axs = [ax for ax in axs if ax._range_subplotspec(x)[idx] == coord] or [self]
        out = []
        for ax in axs:
            other = getattr(ax, '_share' + y)
            if other and other._panel_parent:  # this is a shared panel
                ax = other
            out.append(ax)
        return out

    def _get_topmost_axes(self):
        """
        Return the topmost axes including panels and parents.
        """
        for _ in range(5):
            self = self._axes or self
            self = self._panel_parent or self
        return self

    def _get_legend_handles(self, handler_map=None):
        """
        Internal implementation of matplotlib's ``get_legend_handles_labels``.
        """
        if not self._panel_hidden:  # this is a normal axes
            axs = [self]
        elif self._panel_parent:  # this is an axes-wide legend
            axs = list(self._panel_parent._iter_axes(hidden=False, children=True))
        else:  # this is a figure-wide legend
            axs = list(self.figure._iter_axes(hidden=False, children=True))
        handles = []
        handler_map_full = mlegend.Legend.get_default_handler_map()
        handler_map_full = handler_map_full.copy()
        handler_map_full.update(handler_map or {})
        for ax in axs:
            for attr in ('lines', 'patches', 'collections', 'containers'):
                for handle in getattr(ax, attr, []):  # guard against API changes
                    label = handle.get_label()
                    handler = mlegend.Legend.get_legend_handler(handler_map_full, handle)  # noqa: E501
                    if handler and label and label[0] != '_':
                        handles.append(handle)
        return handles

    def _get_size_inches(self):
        """
        Return the width and height of the axes in inches.
        """
        width, height = self.figure.get_size_inches()
        bbox = self.get_position()
        width = width * abs(bbox.width)
        height = height * abs(bbox.height)
        return np.array([width, height])

    def _get_transform(self, transform, default='data'):
        """
        Translates user input transform. Also used in an axes method.
        """
        # TODO: Can this support cartopy transforms? Seems not when this
        # is used for inset axes bounds but maybe in other places?
        transform = _not_none(transform, default)
        if isinstance(transform, mtransforms.Transform):
            return transform
        elif CRS is not object and isinstance(transform, CRS):
            return transform
        elif PlateCarree is not object and transform == 'map':
            return PlateCarree()
        elif transform == 'data':
            return self.transData
        elif transform == 'axes':
            return self.transAxes
        elif transform == 'figure':
            return self.figure.transFigure
        elif transform == 'subfigure':
            return self.figure.transSubfigure
        else:
            raise ValueError(f'Unknown transform {transform!r}.')

    def _register_guide(self, guide, obj, key, **kwargs):
        """
        Queue up or replace objects for legends and list-of-artist style colorbars.
        """
        # Initial stuff
        if guide not in ('legend', 'colorbar'):
            raise TypeError(f'Invalid type {guide!r}.')
        dict_ = self._legend_dict if guide == 'legend' else self._colorbar_dict

        # Remove previous instances
        # NOTE: No good way to remove inset colorbars right now until the bounding
        # box and axes are merged into some kind of subclass. Just fine for now.
        if key in dict_ and not isinstance(dict_[key], tuple):
            prev = dict_.pop(key)  # possibly pop a queued object
            if guide == 'colorbar':
                pass
            elif hasattr(self, 'legend_') and prev.axes.legend_ is prev:
                self.legend_ = None  # was never added as artist
            else:
                prev.remove()  # remove legends and inner colorbars

        # Replace with instance or update the queue
        # NOTE: This is valid for both mappable-values pairs and handles-labels pairs
        if not isinstance(obj, tuple) or any(isinstance(_, mlegend.Legend) for _ in obj):  # noqa: E501
            dict_[key] = obj
        else:
            handles, labels = obj
            if not np.iterable(handles) or type(handles) is tuple:
                handles = [handles]
            if not np.iterable(labels) or isinstance(labels, str):
                labels = [labels] * len(handles)
            length = min(len(handles), len(labels))  # mimics 'zip' behavior
            handles_full, labels_full, kwargs_full = dict_.setdefault(key, ([], [], {}))
            handles_full.extend(handles[:length])
            labels_full.extend(labels[:length])
            kwargs_full.update(kwargs)

    def _update_guide(
        self, objs, colorbar=None, colorbar_kw=None, queue_colorbar=True,
        legend=None, legend_kw=None,
    ):
        """
        Update queues for on-the-fly legends and colorbars or track keyword arguments.
        """
        # TODO: Support auto-splitting artists passed to legend into
        # their legend elements. Play with this.
        # WARNING: This should generally be last in the pipeline before calling
        # the plot function or looping over data columns. The colormap parser
        # and standardize functions both modify colorbar_kw and legend_kw.
        if colorbar:
            colorbar_kw = colorbar_kw or {}
            colorbar_kw.setdefault('queue', queue_colorbar)
            self.colorbar(objs, loc=colorbar, **colorbar_kw)
        else:
            guides._guide_kw_to_obj(objs, 'colorbar', colorbar_kw)  # save for later
        if legend:
            legend_kw = legend_kw or {}
            self.legend(objs, loc=legend, queue=True, **legend_kw)
        else:
            guides._guide_kw_to_obj(objs, 'legend', legend_kw)  # save for later

    @staticmethod
    def _parse_frame(guide, fancybox=None, shadow=None, **kwargs):
        """
        Parse frame arguments.
        """
        # NOTE: Here we permit only 'edgewidth' to avoid conflict with
        # 'linewidth' used for legend handles and colorbar edge.
        kw_frame = _pop_kwargs(
            kwargs,
            alpha=('a', 'framealpha', 'facealpha'),
            facecolor=('fc', 'framecolor', 'facecolor'),
            edgecolor=('ec',),
            edgewidth=('ew',),
        )
        _kw_frame_default = {
            'alpha': f'{guide}.framealpha',
            'facecolor': f'{guide}.facecolor',
            'edgecolor': f'{guide}.edgecolor',
            'edgewidth': 'axes.linewidth',
        }
        for key, name in _kw_frame_default.items():
            kw_frame.setdefault(key, rc[name])
        for key in ('facecolor', 'edgecolor'):
            if kw_frame[key] == 'inherit':
                kw_frame[key] = rc['axes.' + key]
        kw_frame['linewidth'] = kw_frame.pop('edgewidth')
        kw_frame['fancybox'] = _not_none(fancybox, rc[f'{guide}.fancybox'])
        kw_frame['shadow'] = _not_none(shadow, rc[f'{guide}.shadow'])
        return kw_frame, kwargs

    @staticmethod
    def _parse_colorbar_arg(
        mappable, values=None, *,
        norm=None, norm_kw=None, vmin=None, vmax=None, **kwargs
    ):
        """
        Generate a mappable from flexible non-mappable input. Useful in bridging
        the gap between legends and colorbars (e.g., creating colorbars from line
        objects whose data values span a natural colormap range).
        """
        # For container objects, we just assume color is the same for every item.
        # Works for ErrorbarContainer, StemContainer, BarContainer.
        if (
            np.iterable(mappable)
            and len(mappable) > 0
            and all(isinstance(obj, mcontainer.Container) for obj in mappable)
        ):
            mappable = [obj[0] for obj in mappable]

        # A colormap instance
        if isinstance(mappable, mcolors.Colormap) or isinstance(mappable, str):
            cmap = constructor.Colormap(mappable)
            if values is None and isinstance(cmap, pcolors.DiscreteColormap):
                values = [None] * cmap.N  # sometimes use discrete norm

        # List of colors
        elif np.iterable(mappable) and all(map(mcolors.is_color_like, mappable)):
            cmap = pcolors.DiscreteColormap(list(mappable), '_no_name')
            if values is None:
                values = [None] * len(mappable)  # always use discrete norm

        # List of artists
        # NOTE: Do not check for isinstance(Artist) in case it is an mpl collection
        elif np.iterable(mappable) and all(
            hasattr(obj, 'get_color') or hasattr(obj, 'get_facecolor') for obj in mappable  # noqa: E501
        ):
            # Generate colormap from colors and infer tick labels
            colors = []
            for obj in mappable:
                if hasattr(obj, 'update_scalarmappable'):  # for e.g. pcolor
                    obj.update_scalarmappable()
                color = obj.get_color() if hasattr(obj, 'get_color') else obj.get_facecolor()  # noqa: E501
                if isinstance(color, np.ndarray):
                    color = color.squeeze()  # e.g. single color scatter plot
                if not mcolors.is_color_like(color):
                    raise ValueError('Cannot make colorbar from artists with more than one color.')  # noqa: E501
                colors.append(color)
            # Try to infer tick values and tick labels from Artist labels
            cmap = pcolors.DiscreteColormap(colors, '_no_name')
            if values is None:
                values = [None] * len(mappable)
            else:
                values = list(values)
            for i, (obj, val) in enumerate(zip(mappable, values)):
                if val is not None:
                    continue
                val = obj.get_label()
                if val and val[0] == '_':
                    continue
                values[i] = val

        else:
            raise ValueError(
                'Input colorbar() argument must be a scalar mappable, colormap name '
                f'or object, list of colors, or list of artists. Got {mappable!r}.'
            )

        # Generate continuous normalizer, and possibly discrete normalizer. Update
        # outgoing locator and formatter if user does not override.
        norm_kw = norm_kw or {}
        norm = norm or 'linear'
        vmin = _not_none(vmin=vmin, norm_kw_vmin=norm_kw.pop('vmin', None), default=0)
        vmax = _not_none(vmax=vmax, norm_kw_vmax=norm_kw.pop('vmax', None), default=1)
        norm = constructor.Norm(norm, vmin=vmin, vmax=vmax, **norm_kw)
        locator = formatter = None
        if values is not None:
            ticks = []
            for i, val in enumerate(values):
                try:
                    val = float(val)
                except (TypeError, ValueError):
                    pass
                if val is None:
                    val = i
                ticks.append(val)
            if not any(isinstance(_, str) for _ in ticks):
                locator = pticker.DiscreteLocator(ticks)
            else:
                formatter = mticker.FixedFormatter(list(map(str, ticks)))
                ticks = np.arange(len(ticks))
                locator = mticker.FixedLocator(ticks)
            if len(ticks) == 1:
                levels = [ticks[0] - 1, ticks[0] + 1]
            else:
                levels = edges(ticks)
            norm, cmap, _ = Axes._parse_discrete_norm(levels, norm, cmap)

        # Return ad hoc ScalarMappable and update locator and formatter
        # NOTE: If value list doesn't match this may cycle over colors.
        mappable = mcm.ScalarMappable(norm, cmap)
        return mappable, locator, formatter, kwargs

    def _parse_colorbar_filled(
        self, length=None, align=None, tickloc=None, ticklocation=None,
        orientation=None, **kwargs
    ):
        """
        Return the axes and adjusted keyword args for a panel-filling colorbar.
        """
        # Parse input arguments
        side = self._panel_side
        side = _not_none(side, 'left' if orientation == 'vertical' else 'bottom')
        align = _not_none(align, 'center')
        length = _not_none(length=length, default=rc['colorbar.length'])
        ticklocation = _not_none(tickloc=tickloc, ticklocation=ticklocation)

        # Calculate inset bounds for the colorbar
        delta = 0.5 * (1 - length)
        if side in ('bottom', 'top'):
            if align == 'left':
                bounds = (0, 0, length, 1)
            elif align == 'center':
                bounds = (delta, 0, length, 1)
            elif align == 'right':
                bounds = (2 * delta, 0, length, 1)
            else:
                raise ValueError(f'Invalid align={align!r} for colorbar loc={side!r}.')
        else:
            if align == 'bottom':
                bounds = (0, 0, 1, length)
            elif align == 'center':
                bounds = (0, delta, 1, length)
            elif align == 'top':
                bounds = (0, 2 * delta, 1, length)
            else:
                raise ValueError(f'Invalid align={align!r} for colorbar loc={side!r}.')

        # Add the axes as a child of the original axes
        cls = mprojections.get_projection_class('proplot_cartesian')
        locator = self._make_inset_locator(bounds, self.transAxes)
        ax = cls(self.figure, locator(self, None).bounds, zorder=5)
        self.add_child_axes(ax)
        ax.set_axes_locator(locator)
        ax.patch.set_facecolor('none')  # ignore axes.alpha application

        # Handle default keyword args
        if side in ('bottom', 'top'):
            outside, inside = 'bottom', 'top'
            if side == 'top':
                outside, inside = inside, outside
            ticklocation = _not_none(ticklocation, outside)
            orientation = _not_none(orientation, 'horizontal')
        else:
            outside, inside = 'left', 'right'
            if side == 'right':
                outside, inside = inside, outside
            ticklocation = _not_none(ticklocation, outside)
            orientation = _not_none(orientation, 'vertical')
        kwargs.update({'orientation': orientation, 'ticklocation': ticklocation})
        return ax, kwargs

    def _parse_colorbar_inset(
        self, loc=None, width=None, length=None, shrink=None,
        frame=None, frameon=None, label=None, pad=None,
        tickloc=None, ticklocation=None, orientation=None, **kwargs,
    ):
        """
        Return the axes and adjusted keyword args for an inset colorbar.
        """
        # Basic colorbar properties
        frame = _not_none(frame=frame, frameon=frameon, default=rc['colorbar.frameon'])
        length = _not_none(length=length, shrink=shrink, default=rc['colorbar.insetlength'])  # noqa: E501
        width = _not_none(width, rc['colorbar.insetwidth'])
        pad = _not_none(pad, rc['colorbar.insetpad'])
        length = units(length, 'em', 'ax', axes=self, width=True)  # x direction
        width = units(width, 'em', 'ax', axes=self, width=False)  # y direction
        xpad = units(pad, 'em', 'ax', axes=self, width=True)
        ypad = units(pad, 'em', 'ax', axes=self, width=False)

        # Extra space accounting for colorbar label and tick labels
        labspace = rc['xtick.major.size'] / 72
        fontsize = rc['xtick.labelsize']
        fontsize = _fontsize_to_pt(fontsize)
        if label is not None:
            labspace += 2.4 * fontsize / 72
        else:
            labspace += 1.2 * fontsize / 72
        labspace /= self._get_size_inches()[1]  # space for labels

        # Location in axes-relative coordinates
        # Bounds are x0, y0, width, height in axes-relative coordinates
        if loc == 'upper right':
            bounds_inset = [1 - xpad - length, 1 - ypad - width]
            bounds_frame = [1 - 2 * xpad - length, 1 - 2 * ypad - width - labspace]
        elif loc == 'upper left':
            bounds_inset = [xpad, 1 - ypad - width]
            bounds_frame = [0, 1 - 2 * ypad - width - labspace]
        elif loc == 'lower left':
            bounds_inset = [xpad, ypad + labspace]
            bounds_frame = [0, 0]
        else:
            bounds_inset = [1 - xpad - length, ypad + labspace]
            bounds_frame = [1 - 2 * xpad - length, 0]
        bounds_inset.extend((length, width))  # inset axes
        bounds_frame.extend((2 * xpad + length, 2 * ypad + width + labspace))

        # Make axes and frame with zorder matching default legend zorder
        cls = mprojections.get_projection_class('proplot_cartesian')
        locator = self._make_inset_locator(bounds_inset, self.transAxes)
        ax = cls(self.figure, locator(self, None).bounds, zorder=5)
        ax.patch.set_facecolor('none')
        ax.set_axes_locator(locator)
        self.add_child_axes(ax)
        kw_frame, kwargs = self._parse_frame('colorbar', **kwargs)
        if frame:
            frame = self._add_guide_frame(*bounds_frame, fontsize=fontsize, **kw_frame)

        # Handle default keyword args
        if orientation is not None and orientation != 'horizontal':
            warnings._warn_proplot(
                f'Orientation for inset colorbars must be horizontal, '
                f'ignoring orientation={orientation!r}.'
            )
        ticklocation = _not_none(tickloc=tickloc, ticklocation=ticklocation)
        if ticklocation is not None and ticklocation != 'bottom':
            warnings._warn_proplot('Inset colorbars can only have ticks on the bottom.')
        kwargs.update({'orientation': 'horizontal', 'ticklocation': 'bottom'})
        return ax, kwargs

    @staticmethod
    def _parse_discrete_norm(
        levels, norm, cmap, *, extend=None, min_levels=None, **kwargs,
    ):
        """
        Create a `~proplot.colors.DiscreteNorm` or `~proplot.colors.BoundaryNorm`
        from the input colormap and normalizer.

        Parameters
        ----------
        levels : sequence of float
            The level boundaries.
        norm : `~matplotlib.colors.Normalize`
            The continuous normalizer.
        cmap : `~matplotlib.colors.Colormap`
            The colormap.
        extend : str, optional
            The extend setting.
        min_levels : int, optional
            The minimum number of levels.

        Returns
        -------
        norm : `~proplot.colors.DiscreteNorm`
            The discrete normalizer.
        cmap : `~matplotlib.colors.Colormap`
            The possibly-modified colormap.
        kwargs
            Unused arguments.
        """
        # Reverse the colormap if input levels or values were descending
        # See _parse_level_list for details
        min_levels = _not_none(min_levels, 2)  # 1 for contour plots
        unique = extend = _not_none(extend, 'neither')
        under = cmap._rgba_under
        over = cmap._rgba_over
        cyclic = getattr(cmap, '_cyclic', None)
        qualitative = isinstance(cmap, pcolors.DiscreteColormap)  # see _parse_cmap
        if len(levels) < min_levels:
            raise ValueError(
                f'Invalid levels={levels!r}. Must be at least length {min_levels}.'
            )

        # Ensure end colors are unique by scaling colors as if extend='both'
        # NOTE: Inside _parse_cmap should have enforced extend='neither'
        if cyclic:
            step = 0.5  # try to allocate space for unique end colors
            unique = 'both'

        # Ensure color list length matches level list length using rotation
        # NOTE: No harm if not enough colors, we just end up with the same
        # color for out-of-bounds extensions. This is a gentle failure
        elif qualitative:
            step = 0.5  # try to sample the central index for safety
            unique = 'both'
            auto_under = under is None and extend in ('min', 'both')
            auto_over = over is None and extend in ('max', 'both')
            ncolors = len(levels) - min_levels + 1 + auto_under + auto_over
            colors = list(itertools.islice(itertools.cycle(cmap.colors), ncolors))
            if auto_under and len(colors) > 1:
                under, *colors = colors
            if auto_over and len(colors) > 1:
                *colors, over = colors
            cmap = cmap.copy(colors, N=len(colors))
            if under is not None:
                cmap.set_under(under)
            if over is not None:
                cmap.set_over(over)

        # Ensure middle colors sample full range when extreme colors are present
        # by scaling colors as if extend='neither'
        else:
            step = 1.0
            if over is not None and under is not None:
                unique = 'neither'
            elif over is not None:  # turn off over-bounds unique bin
                if extend == 'both':
                    unique = 'min'
                elif extend == 'max':
                    unique = 'neither'
            elif under is not None:  # turn off under-bounds unique bin
                if extend == 'both':
                    unique = 'min'
                elif extend == 'max':
                    unique = 'neither'

        # Generate DiscreteNorm and update "child" norm with vmin and vmax from
        # levels. This lets the colorbar set tick locations properly!
        if not isinstance(norm, mcolors.BoundaryNorm) and len(levels) > 1:
            norm = pcolors.DiscreteNorm(levels, norm=norm, unique=unique, step=step)

        return norm, cmap, kwargs

    def _parse_legend_aligned(self, pairs, ncol=None, order=None, **kwargs):
        """
        Draw an individual legend with aligned columns. Includes support
        for switching legend-entries between column-major and row-major.
        """
        # Potentially change the order of handles to column-major
        npairs = len(pairs)
        ncol = _not_none(ncol, 3)
        nrow = npairs // ncol + 1
        array = np.empty((nrow, ncol), dtype=object)
        for i, pair in enumerate(pairs):
            array.flat[i] = pair  # must be assigned individually
        if order == 'C':
            array = array.T

        # Return a legend
        # NOTE: Permit drawing empty legend to catch edge cases
        pairs = [pair for pair in array.flat if isinstance(pair, tuple)]
        args = tuple(zip(*pairs)) or ([], [])
        return mlegend.Legend(self, *args, ncol=ncol, **kwargs)

    def _parse_legend_centered(
        self, pairs, *, fontsize,
        loc=None, title=None, frameon=None, kw_frame=None, **kwargs
    ):
        """
        Draw "legend" with centered rows by creating separate legends for
        each row. The label spacing/border spacing will be exactly replicated.
        """
        # Parse input args
        # NOTE: Main legend() function applies default 'legend.loc' of 'best' when
        # users pass legend=True or call legend without 'loc'. Cannot issue warning.
        kw_frame = kw_frame or {}
        kw_frame['fontsize'] = fontsize
        if loc is None or loc == 'best':  # white lie
            loc = 'upper center'
        if not isinstance(loc, str):
            raise ValueError(
                f'Invalid loc={loc!r} for centered-row legend. Must be string.'
            )
        keys = ('bbox_transform', 'bbox_to_anchor')
        kw_ignore = {key: kwargs.pop(key) for key in keys if key in kwargs}
        if kw_ignore:
            warnings._warn_proplot(
                f'Ignoring invalid centered-row legend keyword args: {kw_ignore!r}'
            )

        # Iterate and draw
        # NOTE: Empirical testing shows spacing fudge factor necessary to
        # exactly replicate the spacing of standard aligned legends.
        # NOTE: We confine possible bounding box in *y*-direction, but do not
        # confine it in *x*-direction. Matplotlib will automatically move
        # left-to-right if you request this.
        legs = []
        kwargs.update({'loc': loc, 'frameon': False})
        space = kwargs.get('labelspacing', None) or rc['legend.labelspacing']
        height = (((1 + space * 0.85) * fontsize) / 72) / self._get_size_inches()[1]
        for i, ipairs in enumerate(pairs):
            extra = int(i > 0 and title is not None)
            if 'upper' in loc:
                base, offset = 1, -extra
            elif 'lower' in loc:
                base, offset = 0, len(pairs)
            else:  # center
                base, offset = 0.5, 0.5 * (len(pairs) - extra)
            y0, y1 = base + (offset - np.array([i + 1, i])) * height
            bb = mtransforms.Bbox([[0, y0], [1, y1]])
            leg = mlegend.Legend(
                self, *zip(*ipairs), bbox_to_anchor=bb, bbox_transform=self.transAxes,
                ncol=len(ipairs), title=title if i == 0 else None, **kwargs
            )
            legs.append(leg)

        # Draw manual fancy bounding box for un-aligned legend
        # WARNING: legendPatch uses the default transform, i.e. universal coordinates
        # in points. Means we have to transform mutation scale into transAxes sizes.
        # WARNING: Tempting to use legendPatch for everything but for some reason
        # coordinates are messed up. In some tests all coordinates were just result
        # of get window extent multiplied by 2 (???). Anyway actual box is found in
        # _legend_box attribute, which is accessed by get_window_extent.
        objs = tuple(legs)
        if frameon and legs:
            rend = self.figure._get_renderer()  # arbitrary renderer
            trans = self.transAxes.inverted()
            bboxes = [leg.get_window_extent(rend).transformed(trans) for leg in legs]
            bb = mtransforms.Bbox.union(bboxes)
            bounds = (bb.xmin, bb.ymin, bb.xmax - bb.xmin, bb.ymax - bb.ymin)
            self._add_guide_frame(*bounds, **kw_frame)
        return objs

    @staticmethod
    def _parse_legend_group(handles, labels=None):
        """
        Parse possibly tuple-grouped input handles.
        """
        # Helper function. Retrieve labels from a tuple group or from objects
        # in a container. Multiple labels lead to multiple legend entries.
        def _legend_label(*objs):  # noqa: E301
            labs = []
            for obj in objs:
                if hasattr(obj, 'get_label'):  # e.g. silent list
                    lab = obj.get_label()
                    if lab is not None and str(lab)[:1] != '_':
                        labs.append(lab)
            return tuple(labs)

        # Helper function. Translate handles in the input tuple group. Extracts
        # legend handles from contour sets and extracts labeled elements from
        # matplotlib containers (important for histogram plots).
        ignore = (mcontainer.ErrorbarContainer,)
        containers = (cbook.silent_list, mcontainer.Container)
        def _legend_tuple(*objs):  # noqa: E306
            handles = []
            for obj in objs:
                if isinstance(obj, ignore) and not _legend_label(obj):
                    continue
                if hasattr(obj, 'update_scalarmappable'):  # for e.g. pcolor
                    obj.update_scalarmappable()
                if isinstance(obj, mcontour.ContourSet):  # extract single element
                    hs, _ = obj.legend_elements()
                    label = getattr(obj, '_legend_label', '_no_label')
                    if hs:  # non-empty
                        obj = hs[len(hs) // 2]
                        obj.set_label(label)
                if isinstance(obj, containers):  # extract labeled elements
                    hs = (obj, *guides._iter_iterables(obj))
                    hs = tuple(filter(_legend_label, hs))
                    if hs:
                        handles.extend(hs)
                    elif obj:  # fallback to first element
                        handles.append(obj[0])
                    else:
                        handles.append(obj)
                elif hasattr(obj, 'get_label'):
                    handles.append(obj)
                else:
                    warnings._warn_proplot(f'Ignoring invalid legend handle {obj!r}.')
            return tuple(handles)

        # Sanitize labels. Ignore e.g. extra hist() or hist2d() return values,
        # auto-detect labels in tuple group, auto-expand tuples with diff labels
        # NOTE: Allow handles and labels of different length like
        # native matplotlib. Just truncate extra values with zip().
        if labels is None:
            labels = [None] * len(handles)
        ihandles, ilabels = [], []
        for hs, label in zip(handles, labels):
            # Filter objects
            if type(hs) is not tuple:  # ignore Containers (tuple subclasses)
                hs = (hs,)
            hs = _legend_tuple(*hs)
            labs = _legend_label(*hs)
            if not hs:
                continue
            # Unfurl tuple of handles
            if label is None and len(labs) > 1:
                hs = tuple(filter(_legend_label, hs))
                ihandles.extend(hs)
                ilabels.extend(_.get_label() for _ in hs)
            # Append this handle with some name
            else:
                hs = hs[0] if len(hs) == 1 else hs  # unfurl for better error messages
                label = label if label is not None else labs[0] if labs else '_no_label'
                ihandles.append(hs)
                ilabels.append(label)
        return ihandles, ilabels

    def _parse_legend_handles(
        self, handles, labels, ncol=None, order=None, center=None,
        alphabetize=None, handler_map=None,
    ):
        """
        Parse input handles and labels.
        """
        # Handle lists of lists
        # TODO: Often desirable to label a "mappable" with one data value. Maybe add a
        # legend option for the *number of samples* or *sample points* when drawing
        # legends for mappables. Look into "legend handlers", might just want to add
        # handlers by passing handler_map to legend() and get_legend_handles_labels().
        is_list = lambda obj: (  # noqa: E731
            np.iterable(obj) and not isinstance(obj, (str, tuple))
        )
        to_list = lambda obj: (  # noqa: E731
            obj.tolist() if isinstance(obj, np.ndarray)
            else obj if obj is None or is_list(obj) else [obj]
        )
        handles, labels = to_list(handles), to_list(labels)
        if handles and not labels and all(isinstance(h, str) for h in handles):
            handles, labels = labels, handles
        multi = any(is_list(h) and len(h) > 1 for h in (handles or ()))
        if multi and order == 'F':
            warnings._warn_proplot(
                'Column-major ordering of legend handles is not supported '
                'for horizontally-centered legends.'
            )
        if multi and ncol is not None:
            warnings._warn_proplot(
                'Detected list of *lists* of legend handles. Ignoring '
                'the user input property "ncol".'
            )
        if labels and not handles:
            warnings._warn_proplot(
                'Passing labels without handles is unsupported in proplot. '
                'Please explicitly pass the handles to legend() or pass labels '
                "to plotting commands with e.g. plot(data_1d, label='label') or "
                "plot(data_2d, labels=['label1', 'label2', ...]). After passing "
                'labels to plotting commands you can call legend() without any '
                'arguments or with the handles as a sole positional argument.'
            )
        ncol = _not_none(ncol, 3)
        center = _not_none(center, multi)

        # Iterate over each sublist and parse independently
        pairs = []
        if not multi:  # temporary
            handles, labels = [handles], [labels]
        elif labels is None:
            labels = [labels] * len(handles)
        for ihandles, ilabels in zip(handles, labels):
            ihandles, ilabels = to_list(ihandles), to_list(ilabels)
            if ihandles is None:
                ihandles = self._get_legend_handles(handler_map)
            ihandles, ilabels = self._parse_legend_group(ihandles, ilabels)
            ipairs = list(zip(ihandles, ilabels))
            if alphabetize:
                ipairs = sorted(ipairs, key=lambda pair: pair[1])
            pairs.append(ipairs)

        # Manage (handle, label) pairs in context of the 'center' option
        if not multi:
            pairs = pairs[0]
            if center:
                multi = True
                pairs = [pairs[i * ncol:(i + 1) * ncol] for i in range(len(pairs))]
        else:
            if not center:  # standardize format based on input
                multi = False  # no longer is list of lists
                pairs = [pair for ipairs in pairs for pair in ipairs]

        if multi:
            pairs = [ipairs for ipairs in pairs if ipairs]
        return pairs, multi

    def _range_subplotspec(self, s):
        """
        Return the column or row range for the subplotspec.
        """
        if not isinstance(self, maxes.SubplotBase):
            raise RuntimeError('Axes must be a subplot.')
        ss = self.get_subplotspec()
        row1, row2, col1, col2 = ss._get_rows_columns()
        if s == 'x':
            return (col1, col2)
        else:
            return (row1, row2)

    def _range_tightbbox(self, s):
        """
        Return the tight bounding box span from the cached bounding box.
        """
        # TODO: Better testing for axes visibility
        bbox = self._tight_bbox
        if bbox is None:
            return np.nan, np.nan
        if s == 'x':
            return bbox.xmin, bbox.xmax
        else:
            return bbox.ymin, bbox.ymax

    def _sharex_setup(self, sharex, **kwargs):
        """
        Configure x-axis sharing for panels. See also `~CartesianAxes._sharex_setup`.
        """
        self._share_short_axis(sharex, 'left', **kwargs)  # x axis of left panels
        self._share_short_axis(sharex, 'right', **kwargs)
        self._share_long_axis(sharex, 'bottom', **kwargs)  # x axis of bottom panels
        self._share_long_axis(sharex, 'top', **kwargs)

    def _sharey_setup(self, sharey, **kwargs):
        """
        Configure y-axis sharing for panels. See also `~CartesianAxes._sharey_setup`.
        """
        self._share_short_axis(sharey, 'bottom', **kwargs)  # y axis of bottom panels
        self._share_short_axis(sharey, 'top', **kwargs)
        self._share_long_axis(sharey, 'left', **kwargs)  # y axis of left panels
        self._share_long_axis(sharey, 'right', **kwargs)

    def _share_short_axis(self, share, side, **kwargs):
        """
        Share the "short" axes of panels in this subplot with other panels.
        """
        if share is None or self._panel_side:
            return  # if this is a panel
        s = 'x' if side in ('left', 'right') else 'y'
        caxs = self._panel_dict[side]
        paxs = share._panel_dict[side]
        caxs = [pax for pax in caxs if not pax._panel_hidden]
        paxs = [pax for pax in paxs if not pax._panel_hidden]
        for cax, pax in zip(caxs, paxs):  # may be uneven
            getattr(cax, f'_share{s}_setup')(pax, **kwargs)

    def _share_long_axis(self, share, side, **kwargs):
        """
        Share the "long" axes of panels in this subplot with other panels.
        """
        # NOTE: We do not check _panel_share because that only controls
        # sharing with main subplot, not other subplots
        if share is None or self._panel_side:
            return  # if this is a panel
        s = 'x' if side in ('top', 'bottom') else 'y'
        paxs = self._panel_dict[side]
        paxs = [pax for pax in paxs if not pax._panel_hidden]
        for pax in paxs:
            getattr(pax, f'_share{s}_setup')(share, **kwargs)

    def _reposition_subplot(self):
        """
        Reposition the subplot axes.
        """
        # WARNING: In later versions self.numRows, self.numCols, and self.figbox
        # are @property definitions that never go stale but in mpl < 3.4 they are
        # attributes that must be updated explicitly with update_params().
        # WARNING: In early versions matplotlib only removes '_layoutbox' and
        # '_poslayoutbox' when calling public set_position but in later versions it
        # calls set_in_layout(False) which removes children from get_tightbbox().
        # Therefore try to use _set_position() even though it is private
        if not isinstance(self, maxes.SubplotBase):
            raise RuntimeError('Axes must be a subplot.')
        setter = getattr(self, '_set_position', self.set_position)
        if _version_mpl >= 3.4:
            setter(self.get_subplotspec().get_position(self.figure))
        else:
            self.update_params()
            setter(self.figbox)  # equivalent to above

    def _update_abc(self, **kwargs):
        """
        Update the a-b-c label.
        """
        # Properties
        # NOTE: Border props only apply for "inner" title locations so we need to
        # store on the axes whenever they are modified in case the current location
        # is an 'outer' location then re-apply in case 'loc' is subsequently changed
        kw = rc.fill(
            {
                'size': 'abc.size',
                'weight': 'abc.weight',
                'color': 'abc.color',
                'family': 'font.family',
            },
            context=True
        )
        kwb = rc.fill(
            {
                'border': 'abc.border',
                'borderwidth': 'abc.borderwidth',
                'bbox': 'abc.bbox',
                'bboxpad': 'abc.bboxpad',
                'bboxcolor': 'abc.bboxcolor',
                'bboxstyle': 'abc.bboxstyle',
                'bboxalpha': 'abc.bboxalpha',
            },
            context=True,
        )
        self._abc_border_kwargs.update(kwb)

        # A-b-c labels. Build as a...z...aa...zz...aaa...zzz
        # NOTE: The abc string should already be validated here
        abc = rc.find('abc', context=True)  # 1st run, or changed
        if abc is True:
            abc = 'a'
        if abc is False:
            abc = ''
        if abc is None or self.number is None:
            pass
        elif isinstance(abc, str):
            nabc, iabc = divmod(self.number - 1, 26)
            if abc:  # should have been validated to contain 'a' or 'A'
                old = re.search('[aA]', abc).group()  # return first occurrence
                new = (nabc + 1) * ABC_STRING[iabc]
                new = new.upper() if old == 'A' else new
                abc = abc.replace(old, new, 1)  # replace first occurrence
            kw['text'] = abc
        else:
            if self.number > len(abc):
                raise ValueError(
                    f'Invalid abc list length {len(abc)} '
                    f'for axes with number {self.number}.'
                )
            else:
                kw['text'] = abc[self._number - 1]

        # Update a-b-c label
        loc = rc.find('abc.loc', context=True)
        loc = self._abc_loc = _translate_loc(loc or self._abc_loc, 'text')
        if loc not in ('left', 'right', 'center'):
            kw.update(self._abc_border_kwargs)
        kw.update(kwargs)
        self._title_dict['abc'].update(kw)

    def _update_title(self, loc, title=None, **kwargs):
        """
        Update the title at the specified location.
        """
        # Titles, with two workflows here:
        # 1. title='name' and titleloc='position'
        # 2. ltitle='name', rtitle='name', etc., arbitrarily many titles
        # NOTE: This always updates the *current* title and deflection to panels
        # is handled later so that titles set with set_title() are deflected too.
        # See notes in _update_super_labels() and _apply_title_above().
        # NOTE: Matplotlib added axes.titlecolor in version 3.2 but we still use
        # custom title.size, title.weight, title.color properties for retroactive
        # support in older matplotlib versions. First get params and update kwargs.
        kw = rc.fill(
            {
                'size': 'title.size',
                'weight': 'title.weight',
                'color': 'title.color',
                'family': 'font.family',
            },
            context=True
        )
        if 'color' in kw and kw['color'] == 'auto':
            del kw['color']  # WARNING: matplotlib permits invalid color here
        kwb = rc.fill(
            {
                'border': 'title.border',
                'borderwidth': 'title.borderwidth',
                'bbox': 'title.bbox',
                'bboxpad': 'title.bboxpad',
                'bboxcolor': 'title.bboxcolor',
                'bboxstyle': 'title.bboxstyle',
                'bboxalpha': 'title.bboxalpha',
            },
            context=True,
        )
        self._title_border_kwargs.update(kwb)

        # Update the padding settings read at drawtime. Make sure to
        # update them on the panel axes if 'title.above' is active.
        pad = rc.find('abc.titlepad', context=True)
        if pad is not None:
            self._abc_title_pad = pad
        pad = rc.find('title.pad', context=True)  # title
        if pad is not None:
            self._title_pad = pad
            self._set_title_offset_trans(pad)

        # Get the title location. If 'titleloc' was used then transfer text
        # from the old location to the new location.
        if loc is not None:
            loc = _translate_loc(loc, 'text')
        else:
            old = self._title_loc
            loc = rc.find('title.loc', context=True)
            loc = self._title_loc = _translate_loc(loc or self._title_loc, 'text')
            if loc != old and old is not None:
                labels._transfer_label(self._title_dict[old], self._title_dict[loc])

        # Update the title text. For outer panels, add text to the panel if
        # necesssary. For inner panels, use the border and bbox settings.
        if loc not in ('left', 'right', 'center'):
            kw.update(self._title_border_kwargs)
        if title is None:
            pass
        elif isinstance(title, str):
            kw['text'] = title
        elif np.iterable(title) and all(isinstance(_, str) for _ in title):
            if self.number is None:
                pass
            elif self.number > len(title):
                raise ValueError(
                    f'Invalid title list length {len(title)} '
                    f'for axes with number {self.number}.'
                )
            else:
                kw['text'] = title[self.number - 1]
        else:
            raise ValueError(f'Invalid title {title!r}. Must be string(s).')
        kw.update(kwargs)
        self._title_dict[loc].update(kw)

    def _update_title_position(self, renderer):
        """
        Update the position of inset titles and outer titles. This is called
        by matplotlib at drawtime.
        """
        # Update title positions
        # NOTE: Critical to do this every time in case padding changes or
        # we added or removed an a-b-c label in the same position as a title
        width, height = self._get_size_inches()
        x_pad = self._title_pad / (72 * width)
        y_pad = self._title_pad / (72 * height)
        for loc, obj in self._title_dict.items():
            x, y = (0, 1)
            if loc == 'abc':  # redirect
                loc = self._abc_loc
            if loc == 'left':
                x = 0
            elif loc == 'center':
                x = 0.5
            elif loc == 'right':
                x = 1
            if loc in ('upper center', 'lower center'):
                x = 0.5
            elif loc in ('upper left', 'lower left'):
                x = x_pad
            elif loc in ('upper right', 'lower right'):
                x = 1 - x_pad
            if loc in ('upper left', 'upper right', 'upper center'):
                y = 1 - y_pad
            elif loc in ('lower left', 'lower right', 'lower center'):
                y = y_pad
            obj.set_position((x, y))

        # Get title padding. Push title above tick marks since matplotlib ignores them.
        # This is known matplotlib problem but especially annoying with top panels.
        # NOTE: See axis.get_ticks_position for inspiration
        pad = self._title_pad
        abcpad = self._abc_title_pad
        if self.xaxis.get_visible() and any(
            tick.tick2line.get_visible() and not tick.label2.get_visible()
            for tick in self.xaxis.majorTicks
        ):
            pad += self.xaxis.get_tick_padding()

        # Avoid applying padding on every draw in case it is expensive to change
        # the title Text transforms every time.
        pad_current = self._title_pad_current
        if pad_current is None or not np.isclose(pad, pad_current):
            self._title_pad_current = pad
            self._set_title_offset_trans(pad)

        # Adjust the above-axes positions with builtin algorithm
        # WARNING: Make sure the name of this private function doesn't change
        super()._update_title_position(renderer)

        # Sync the title position with the a-b-c label position
        aobj = self._title_dict['abc']
        tobj = self._title_dict[self._abc_loc]
        aobj.set_transform(tobj.get_transform())
        aobj.set_position(tobj.get_position())
        aobj.set_ha(tobj.get_ha())
        aobj.set_va(tobj.get_va())

        # Offset title away from a-b-c label
        # NOTE: Title texts all use axes transform in x-direction
        if not tobj.get_text() or not aobj.get_text():
            return
        awidth, twidth = (
            obj.get_window_extent(renderer).transformed(self.transAxes.inverted())
            .width for obj in (aobj, tobj)
        )
        ha = aobj.get_ha()
        pad = (abcpad / 72) / self._get_size_inches()[0]
        aoffset = toffset = 0
        if ha == 'left':
            toffset = awidth + pad
        elif ha == 'right':
            aoffset = -(twidth + pad)
        else:  # guaranteed center, there are others
            toffset = 0.5 * (awidth + pad)
            aoffset = -0.5 * (twidth + pad)
        aobj.set_x(aobj.get_position()[0] + aoffset)
        tobj.set_x(tobj.get_position()[0] + toffset)

    def _update_super_title(self, suptitle=None, **kwargs):
        """
        Update the figure super title.
        """
        # NOTE: This is actually *figure-wide* setting, but that line gets blurred
        # where we have shared axes, spanning labels, etc. May cause redundant
        # assignments if using SubplotGrid.format() but this is fast so nbd.
        if self.number is None:
            # NOTE: Kludge prevents changed *figure-wide* settings from getting
            # overwritten when user makes a new panels or insets. Funky limitation but
            # kind of makes sense to make these inaccessible from panels.
            return
        kw = rc.fill(
            {
                'size': 'suptitle.size',
                'weight': 'suptitle.weight',
                'color': 'suptitle.color',
                'family': 'font.family'
            },
            context=True,
        )
        kw.update(kwargs)
        if suptitle or kw:
            self.figure._update_super_title(suptitle, **kw)

    def _update_super_labels(self, side, labels=None, **kwargs):
        """
        Update the figure super labels.
        """
        fig = self.figure
        if self.number is None:
            return  # NOTE: see above
        kw = rc.fill(
            {
                'color': side + 'label.color',
                'rotation': side + 'label.rotation',
                'size': side + 'label.size',
                'weight': side + 'label.weight',
                'family': 'font.family'
            },
            context=True,
        )
        kw.update(kwargs)
        if labels or kw:
            fig._update_super_labels(side, labels, **kw)

    @docstring._snippet_manager
    def format(
        self, *, title=None, title_kw=None, abc_kw=None,
        ltitle=None, lefttitle=None,
        ctitle=None, centertitle=None,
        rtitle=None, righttitle=None,
        ultitle=None, upperlefttitle=None,
        uctitle=None, uppercentertitle=None,
        urtitle=None, upperrighttitle=None,
        lltitle=None, lowerlefttitle=None,
        lctitle=None, lowercentertitle=None,
        lrtitle=None, lowerrighttitle=None,
        **kwargs
    ):
        """
        Modify the a-b-c label, axes title(s), and background patch,
        and call `proplot.figure.Figure.format` on the axes figure.

        Parameters
        ----------
        %(axes.format)s

        Important
        ---------
        `abc`, `abcloc`, `titleloc`, `titleabove`, `titlepad`, and
        `abctitlepad` are actually :ref:`configuration settings <ug_config>`.
        We explicitly document these arguments here because it is common to
        change them for specific axes. But many :ref:`other configuration
        settings <ug_format>` can be passed to ``format`` too.

        Other parameters
        ----------------
        %(figure.format)s
        %(rc.format)s

        See also
        --------
        proplot.axes.CartesianAxes.format
        proplot.axes.PolarAxes.format
        proplot.axes.GeoAxes.format
        proplot.figure.Figure.format
        proplot.gridspec.SubplotGrid.format
        proplot.config.Configurator.context
        """
        skip_figure = kwargs.pop('skip_figure', False)  # internal keyword arg
        params = _pop_params(kwargs, self.figure._format_signature)

        # Initiate context block
        rc_kw, rc_mode = _pop_rc(kwargs)
        with rc.context(rc_kw, mode=rc_mode):
            # Behavior of titles in presence of panels
            above = rc.find('title.above', context=True)
            if above is not None:
                self._title_above = above  # used for future titles

            # Update a-b-c label and titles
            abc_kw = abc_kw or {}
            title_kw = title_kw or {}
            self._update_abc(**abc_kw)
            self._update_title(
                None,
                title,
                **title_kw
            )
            self._update_title(
                'left',
                _not_none(ltitle=ltitle, lefttitle=lefttitle),
                **title_kw,
            )
            self._update_title(
                'center',
                _not_none(ctitle=ctitle, centertitle=centertitle),
                **title_kw,
            )
            self._update_title(
                'right',
                _not_none(rtitle=rtitle, righttitle=righttitle),
                **title_kw,
            )
            self._update_title(
                'upper left',
                _not_none(ultitle=ultitle, upperlefttitle=upperlefttitle),
                **title_kw,
            )
            self._update_title(
                'upper center',
                _not_none(uctitle=uctitle, uppercentertitle=uppercentertitle),
                **title_kw
            )
            self._update_title(
                'upper right',
                _not_none(urtitle=urtitle, upperrighttitle=upperrighttitle),
                **title_kw
            )
            self._update_title(
                'lower left',
                _not_none(lltitle=lltitle, lowerlefttitle=lowerlefttitle),
                **title_kw
            )
            self._update_title(
                'lower center',
                _not_none(lctitle=lctitle, lowercentertitle=lowercentertitle),
                **title_kw
            )
            self._update_title(
                'lower right',
                _not_none(lrtitle=lrtitle, lowerrighttitle=lowerrighttitle),
                **title_kw
            )

            # Update the axes style
            # NOTE: This will also raise an error if unknown args are encountered
            cycle = rc.find('axes.prop_cycle', context=True)
            if cycle is not None:
                self.set_prop_cycle(cycle)
            self._update_background(**kwargs)

        # Update super labels and super title
        # NOTE: To avoid resetting figure-wide settings when new axes are created
        # we only proceed if using the default context mode. Simliar to geo.py
        if skip_figure:  # avoid recursion
            return
        if rc_mode == 1:  # avoid resetting
            return
        self.figure.format(rc_kw=rc_kw, rc_mode=rc_mode, skip_axes=True, **params)

    def draw(self, renderer=None, *args, **kwargs):
        # Perform extra post-processing steps
        # NOTE: In *principle* these steps go here but should already be complete
        # because auto_layout() (called by figure pre-processor) has to run them
        # before aligning labels. So these are harmless no-ops.
        self._add_queued_guides()
        self._apply_title_above()
        if self._colorbar_fill:
            self._colorbar_fill.update_ticks(manual_only=True)  # only update if needed!
        if self._inset_parent is not None and self._inset_zoom:
            self.indicate_inset_zoom()
        super().draw(renderer, *args, **kwargs)

    def get_tightbbox(self, renderer, *args, **kwargs):
        # Perform extra post-processing steps
        # NOTE: This should be updated alongside draw(). We also cache the resulting
        # bounding box to speed up tight layout calculations (see _range_tightbbox).
        self._add_queued_guides()
        self._apply_title_above()
        if self._colorbar_fill:
            self._colorbar_fill.update_ticks(manual_only=True)  # only update if needed!
        if self._inset_parent is not None and self._inset_zoom:
            self.indicate_inset_zoom()
        self._tight_bbox = super().get_tightbbox(renderer, *args, **kwargs)
        return self._tight_bbox

    def get_default_bbox_extra_artists(self):
        # Further restrict artists to those with disabled clipping
        # or use the axes bounding box or patch path for clipping.
        # NOTE: Critical to ignore x and y axis, spines, and all child axes.
        # For some reason these have clipping 'enabled' but it is not respected.
        # NOTE: Matplotlib already tries to do this inside get_tightbbox() but
        # their approach fails for cartopy axes clipped by paths and not boxes.
        return [
            artist for artist in super().get_default_bbox_extra_artists()
            if not self._artist_fully_clipped(artist)
        ]

    def set_prop_cycle(self, *args, **kwargs):
        # Silent override. This is a strict superset of matplotlib functionality.
        # Includes both proplot syntax with positional arguments interpreted as
        # color arguments and oldschool matplotlib cycler(key, value) syntax.
        if len(args) == 2 and isinstance(args[0], str) and np.iterable(args[1]):
            if _pop_props({args[0]: object()}, 'line'):
                kwargs = {args[0]: args[1]}
                args = ()
        cycle = self._active_cycle = constructor.Cycle(*args, **kwargs)
        return super().set_prop_cycle(cycle)  # set the property cycler after validation

    @docstring._snippet_manager
    def inset(self, *args, **kwargs):
        """
        %(axes.inset)s
        """
        return self._add_inset_axes(*args, **kwargs)

    @docstring._snippet_manager
    def inset_axes(self, *args, **kwargs):
        """
        %(axes.inset)s
        """
        return self._add_inset_axes(*args, **kwargs)

    @docstring._snippet_manager
    def indicate_inset_zoom(self, **kwargs):
        """
        %(axes.indicate_inset)s
        """
        # Add the inset indicators
        parent = self._inset_parent
        if not parent:
            raise ValueError('This command can only be called from an inset axes.')
        kwargs.update(_pop_props(kwargs, 'patch'))  # impose alternative defaults
        if not self._inset_zoom_artists:
            kwargs.setdefault('zorder', 3.5)
            kwargs.setdefault('linewidth', rc['axes.linewidth'])
            kwargs.setdefault('edgecolor', rc['axes.edgecolor'])
        xlim, ylim = self.get_xlim(), self.get_ylim()
        rect = (xlim[0], ylim[0], xlim[1] - xlim[0], ylim[1] - ylim[0])
        rectpatch, connects = parent.indicate_inset(rect, self)

        # Update indicator properties
        # NOTE: Unlike matplotlib we sync zoom box properties with connection lines.
        if self._inset_zoom_artists:
            rectpatch_prev, connects_prev = self._inset_zoom_artists
            rectpatch.update_from(rectpatch_prev)
            rectpatch.set_zorder(rectpatch_prev.get_zorder())
            rectpatch_prev.remove()
            for line, line_prev in zip(connects, connects_prev):
                line.update_from(line_prev)
                line.set_zorder(line_prev.get_zorder())  # not included in update_from
                line_prev.remove()
        rectpatch.update(kwargs)
        for line in connects:
            line.update(kwargs)
        self._inset_zoom_artists = (rectpatch, connects)
        return rectpatch, connects

    @docstring._snippet_manager
    def panel(self, side=None, **kwargs):
        """
        %(axes.panel)s
        """
        return self.figure._add_axes_panel(self, side, **kwargs)

    @docstring._snippet_manager
    def panel_axes(self, side=None, **kwargs):
        """
        %(axes.panel)s
        """
        return self.figure._add_axes_panel(self, side, **kwargs)

    @docstring._obfuscate_params
    @docstring._snippet_manager
    def colorbar(
        self, mappable, values=None,
        loc=None, location=None, align=None, queue=False, **kwargs
    ):
        """
        Add an inset colorbar or an outer colorbar along the edge of the axes.

        Parameters
        ----------
        %(axes.colorbar_args)s
        loc, location : int or str, default: :rc:`colorbar.loc`
            The colorbar location. Valid location keys are shown in the below table.

            .. _colorbar_table:

            ==================  =======================================
            Location            Valid keys
            ==================  =======================================
            outer left          ``'left'``, ``'l'``
            outer right         ``'right'``, ``'r'``
            outer bottom        ``'bottom'``, ``'b'``
            outer top           ``'top'``, ``'t'``
            default inset       ``'best'``, ``'inset'``, ``'i'``, ``0``
            upper right inset   ``'upper right'``, ``'ur'``, ``1``
            upper left inset    ``'upper left'``, ``'ul'``, ``2``
            lower left inset    ``'lower left'``, ``'ll'``, ``3``
            lower right inset   ``'lower right'``, ``'lr'``, ``4``
            "filled"            ``'fill'``
            ==================  =======================================

        shrink
            Alias for `length`. This is included for consistency with
            `matplotlib.figure.Figure.colorbar`.
        length \
: float or unit-spec, default: :rc:`colorbar.length` or :rc:`colorbar.insetlength`
            The colorbar length. For outer colorbars, units are relative to the axes
            width or height (default is :rcraw:`colorbar.length`). For inset
            colorbars, floats interpreted as em-widths and strings interpreted
            by `~proplot.utils.units` (default is :rcraw:`colorbar.insetlength`).
        width : unit-spec, default: :rc:`colorbar.width` or :rc:`colorbar.insetwidth
            The colorbar width. For outer colorbars, floats are interpreted as inches
            (default is :rcraw:`colorbar.width`). For inset colorbars, floats are
            interpreted as em-widths (default is :rcraw:`colorbar.insetwidth`).
            Strings are interpreted by `~proplot.utils.units`.
        %(axes.colorbar_space)s
            Has no visible effect if `length` is ``1``.

        Other parameters
        ----------------
        %(axes.colorbar_kwargs)s

        See also
        --------
        proplot.figure.Figure.colorbar
        matplotlib.figure.Figure.colorbar
        """
        # Either draw right now or queue up for later
        # The queue option lets us successively append objects (e.g. line handles)
        # to a list later used for colorbar levels. Same as legend.
        loc = _not_none(loc=loc, location=location)
        loc = _translate_loc(loc, 'colorbar', default=rc['colorbar.loc'])
        align = _translate_loc(align, 'panel', default='center', c='center', center='center')  # noqa: E501
        kwargs = guides._guide_kw_from_obj(mappable, 'colorbar', kwargs)
        if queue:
            self._register_guide('colorbar', (mappable, values), (loc, align), **kwargs)
        else:
            cb = self._add_colorbar(mappable, values, loc=loc, align=align, **kwargs)
            return cb

    @docstring._concatenate_inherited  # also obfuscates params
    @docstring._snippet_manager
    def legend(
        self, handles=None, labels=None,
        loc=None, location=None, align=None, queue=False, **kwargs
    ):
        """
        Add an *inset* legend or *outer* legend along the edge of the axes.

        Parameters
        ----------
        %(axes.legend_args)s
        loc, location : int or str, default: :rc:`legend.loc`
            The legend location. Valid location keys are shown in the below table.

            .. _legend_table:

            ==================  =======================================
            Location            Valid keys
            ==================  =======================================
            outer left          ``'left'``, ``'l'``
            outer right         ``'right'``, ``'r'``
            outer bottom        ``'bottom'``, ``'b'``
            outer top           ``'top'``, ``'t'``
            "best" inset        ``'best'``, ``'inset'``, ``'i'``, ``0``
            upper right inset   ``'upper right'``, ``'ur'``, ``1``
            upper left inset    ``'upper left'``, ``'ul'``, ``2``
            lower left inset    ``'lower left'``, ``'ll'``, ``3``
            lower right inset   ``'lower right'``, ``'lr'``, ``4``
            center left inset   ``'center left'``, ``'cl'``, ``5``
            center right inset  ``'center right'``, ``'cr'``, ``6``
            lower center inset  ``'lower center'``, ``'lc'``, ``7``
            upper center inset  ``'upper center'``, ``'uc'``, ``8``
            center inset        ``'center'``, ``'c'``, ``9``
            "filled"            ``'fill'``
            ==================  =======================================

        width : unit-spec, optional
            For outer legends only. The space allocated for the legend
            box. This does nothing if the :ref:`tight layout algorithm
            <ug_tight>` is active for the figure.
            %(units.in)s
        %(axes.legend_space)s

        Other parameters
        ----------------
        %(axes.legend_kwargs)s

        See also
        --------
        proplot.figure.Figure.legend
        matplotlib.axes.Axes.legend
        """
        # Either draw right now or queue up for later
        # Handles can be successively added to a single location this way. This
        # is used internally for on-the-fly legends.
        loc = _not_none(loc=loc, location=location)
        loc = _translate_loc(loc, 'legend', default=rc['legend.loc'])
        align = _translate_loc(align, 'panel', default='center', c='center', center='center')  # noqa: E501
        kwargs = guides._guide_kw_from_obj(handles, 'legend', kwargs)
        if queue:
            self._register_guide('legend', (handles, labels), (loc, align), **kwargs)
        else:
            leg = self._add_legend(handles, labels, loc=loc, align=align, **kwargs)
            return leg

    @docstring._concatenate_inherited
    @docstring._snippet_manager
    def text(
        self, *args, border=False, bbox=False,
        bordercolor='w', borderwidth=2, borderinvert=False, borderstyle='miter',
        bboxcolor='w', bboxstyle='round', bboxalpha=0.5, bboxpad=None, **kwargs
    ):
        """
        Add text to the axes.

        Parameters
        ----------
        x, y, [z] : float
            The coordinates for the text. `~proplot.axes.ThreeAxes` accept an
            optional third coordinate. If only two are provided this automatically
            redirects to the `~mpl_toolkits.mplot3d.Axes3D.text2D` method.
        s, text : str
            The string for the text.
        %(axes.transform)s

        Other parameters
        ----------------
        border : bool, default: False
            Whether to draw border around text.
        borderwidth : float, default: 2
            The width of the text border.
        bordercolor : color-spec, default: 'w'
            The color of the text border.
        borderinvert : bool, optional
            If ``True``, the text and border colors are swapped.
        borderstyle : {'miter', 'round', 'bevel'}, optional
            The `line join style \
<https://matplotlib.org/stable/gallery/lines_bars_and_markers/joinstyle.html>`__
            used for the border.
        bbox : bool, default: False
            Whether to draw a bounding box around text.
        bboxcolor : color-spec, default: 'w'
            The color of the text bounding box.
        bboxstyle : boxstyle, default: 'round'
            The style of the bounding box.
        bboxalpha : float, default: 0.5
            The alpha for the bounding box.
        bboxpad : float, default: :rc:`title.bboxpad`
            The padding for the bounding box.
        %(artist.text)s

        **kwargs
            Passed to `matplotlib.axes.Axes.text`.

        See also
        --------
        matplotlib.axes.Axes.text
        """
        # Translate positional args
        # Audo-redirect to text2D for 3D axes if not enough arguments passed
        # NOTE: The transform must be passed positionally for 3D axes with 2D coords
        keys = 'xy'
        func = super().text
        if self._name == 'three':
            if len(args) >= 4 or 'z' in kwargs:
                keys += 'z'
            else:
                func = self.text2D
        keys = (*keys, ('s', 'text'), 'transform')
        args, kwargs = _kwargs_to_args(keys, *args, **kwargs)
        *args, transform = args
        if any(arg is None for arg in args):
            raise TypeError('Missing required positional argument.')
        if transform is None:
            transform = self.transData
        else:
            transform = self._get_transform(transform)
        with warnings.catch_warnings():  # ignore duplicates (internal issues?)
            warnings.simplefilter('ignore', warnings.ProplotWarning)
            kwargs.update(_pop_props(kwargs, 'text'))

        # Update the text object using a monkey patch
        obj = func(*args, transform=transform, **kwargs)
        obj.update = labels._update_label.__get__(obj)
        obj.update(
            {
                'border': border,
                'bordercolor': bordercolor,
                'borderinvert': borderinvert,
                'borderwidth': borderwidth,
                'borderstyle': borderstyle,
                'bbox': bbox,
                'bboxcolor': bboxcolor,
                'bboxstyle': bboxstyle,
                'bboxalpha': bboxalpha,
                'bboxpad': bboxpad,
            }
        )
        return obj

    def _iter_axes(self, hidden=False, children=False, panels=True):
        """
        Return a list of visible axes, panel axes, and child axes of both.

        Parameters
        ----------
        hidden : bool, optional
            Whether to include "hidden" panels.
        children : bool, optional
            Whether to include children. Note this now includes "twin" axes.
        panels : bool or str or sequence of str, optional
            Whether to include panels or the panels to include.
        """
        # Parse panels
        if panels is False:
            panels = ()
        elif panels is True or panels is None:
            panels = ('left', 'right', 'bottom', 'top')
        elif isinstance(panels, str):
            panels = (panels,)
        if not set(panels) <= {'left', 'right', 'bottom', 'top'}:
            raise ValueError(f'Invalid sides {panels!r}.')
        # Iterate
        axs = (self, *(ax for side in panels for ax in self._panel_dict[side]))
        for iax in axs:
            if not hidden and iax._panel_hidden:
                continue  # ignore hidden panel and its colorbar/legend child
            iaxs = (iax, *(iax.child_axes if children else ()))
            for jax in iaxs:
                if not jax.get_visible():
                    continue  # safety first
                yield jax

    @property
    def number(self):
        """
        The axes number. This controls the order of a-b-c labels and the
        order of appearence in the `~proplot.gridspec.SubplotGrid` returned
        by `~proplot.figure.Figure.subplots`.
        """
        return self._number

    @number.setter
    def number(self, num):
        if num is None or isinstance(num, Integral) and num > 0:
            self._number = num
        else:
            raise ValueError(f'Invalid number {num!r}. Must be integer >=1.')


# Apply signature obfuscation after storing previous signature
# NOTE: This is needed for __init__
Axes._format_signatures = {Axes: inspect.signature(Axes.format)}
Axes.format = docstring._obfuscate_kwargs(Axes.format)
