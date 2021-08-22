#!/usr/bin/env python3
"""
The second-level axes subclass used for all ProPlot figures.
Implements plotting method overrides.
"""
import functools
import inspect
import itertools
import re
import sys
from numbers import Integral

import matplotlib.axes as maxes
import matplotlib.cm as mcm
import matplotlib.collections as mcollections
import matplotlib.colors as mcolors
import matplotlib.contour as mcontour
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import numpy as np
import numpy.ma as ma

from .. import colors as pcolors
from .. import constructor
from .. import ticker as pticker
from ..config import rc
from ..internals import ic  # noqa: F401
from ..internals import (
    _guide_kw_to_arg,
    _guide_kw_to_obj,
    _keyword_to_positional,
    _not_none,
    _pop_kwargs,
    _pop_params,
    _pop_props,
    _process_props,
    _snippet_manager,
    _state_context,
    docstring,
    warnings,
)
from ..utils import edges, edges2d, to_xyz
from . import base

try:
    from cartopy.crs import PlateCarree
except ModuleNotFoundError:
    PlateCarree = object

__all__ = ['PlotAxes']


# Constants
EDGEWIDTH = 0.25  # native linewidth used for grid box edges in matplotlib
BASEMAP_FUNCS = (  # default latlon=True
    'barbs', 'contour', 'contourf', 'hexbin',
    'imshow', 'pcolor', 'pcolormesh', 'plot',
    'quiver', 'scatter', 'streamplot', 'step',
)
CARTOPY_FUNCS = (  # default transform=PlateCarree()
    'barbs', 'contour', 'contourf',
    'fill', 'fill_between', 'fill_betweenx',  # NOTE: not sure if these work
    'imshow', 'pcolor', 'pcolormesh', 'plot',
    'quiver', 'scatter', 'streamplot', 'step',
    'tricontour', 'tricontourf', 'tripcolor',  # NOTE: not sure why these work
)

# Data argument docstrings
_args_1d_docstring = """
*args : {y} or {x}, {y}
    The data passed as positional or keyword arguments. Interpreted as follows:

    * If only `{y}` coordinates are passed, try to infer the `{x}` coordinates
      from the `~pandas.Series` or `~pandas.DataFrame` indices or the
      `~xarray.DataArray` coordinates. Otherwise, the `{x}` coordinates
      are ``np.arange(0, {y}.shape[0])``.
    * If the `{y}` coordinates are a 2D array, plot each column of data in succession
      (except where each column of data represents a statistical distribution, as with
      ``boxplot``, ``violinplot``, or when using ``means=True`` or ``medians=True``).
    * If any arguments are `pint.Quantity`, auto-add the pint unit registry
      to matplotlib's unit registry using `~pint.UnitRegistry.setup_matplotlib`.
      A `pint.Quantity` embedded in an `xarray.DataArray` is also supported.
"""
_args_1d_multi_docstring = """
*args : {y}2 or {x}, {y}2, or {x}, {y}1, {y}2
    The data passed as positional or keyword arguments. Interpreted as follows:

    * If only `{y}` coordinates are passed, try to infer the `{x}` coordinates from
      the `~pandas.Series` or `~pandas.DataFrame` indices or the `~xarray.DataArray`
      coordinates. Otherwise, the `{x}` coordinates are ``np.arange(0, {y}2.shape[0])``.
    * If only `{x}` and `{y}2` coordinates are passed, set the `{y}1` coordinates
      to zero. This draws elements originating from the zero line.
    * If both `{y}1` and `{y}2` are provided, draw elements between these points. If
      either are 2D, draw elements by iterating over each column.
    * If any arguments are `pint.Quantity`, auto-add the pint unit registry
      to matplotlib's unit registry using `~pint.UnitRegistry.setup_matplotlib`.
      A `pint.Quantity` embedded in an `xarray.DataArray` is also supported.
"""
_args_2d_docstring = """
*args : {z} or x, y, {z}
    The data passed as positional or keyword arguments. Interpreted as follows:

    * If only {zvar} coordinates are passed, try to infer the `x` and `y` coordinates
      from the `~pandas.DataFrame` indices and columns or the `~xarray.DataArray`
      coordinates. Otherwise, the `y` coordinates are ``np.arange(0, y.shape[0])``
      and the `x` coordinates are ``np.arange(0, y.shape[1])``.
    * For ``pcolor`` and ``pcolormesh``, calculate coordinate *edges* using
      `~proplot.utils.edges` or `~proplot.utils.edges2d` if *centers* were provided.
      For all other methods, calculate coordinate *centers* if *edges* were provided.
    * If the `x` or `y` coordinates are `pint.Quantity`, auto-add the pint unit registry
      to matplotlib's unit registry using `~pint.UnitRegistry.setup_matplotlib`. If the
      {zvar} coordinates are `pint.Quantity`, pass the magnitude to the plotting
      command. A `pint.Quantity` embedded in an `xarray.DataArray` is also supported.
"""
_snippet_manager['plot.args_1d_y'] = _args_1d_docstring.format(x='x', y='y')
_snippet_manager['plot.args_1d_x'] = _args_1d_docstring.format(x='y', y='x')
_snippet_manager['plot.args_1d_multiy'] = _args_1d_multi_docstring.format(x='x', y='y')
_snippet_manager['plot.args_1d_multix'] = _args_1d_multi_docstring.format(x='y', y='x')
_snippet_manager['plot.args_2d'] = _args_2d_docstring.format(z='z', zvar='`z`')
_snippet_manager['plot.args_2d_flow'] = _args_2d_docstring.format(z='u, v', zvar='`u` and `v`')  # noqa: E501


# Shared docstrings
_args_1d_shared_docstring = """
data : dict-like, optional
    A dict-like dataset container (e.g., `~pandas.DataFrame` or
    `~xarray.DataArray`). If passed, positional arguments can optionally
    be string `data` keys and the arrays used for plotting are retrieved
    with ``data[key]``. This is a `native matplotlib feature
    <https://matplotlib.org/stable/gallery/misc/keyword_plotting.html>`__.
autoformat : bool, optional
    Whether the `x` axis labels, `y` axis labels, axis formatters, axes titles,
    legend titles, and colorbar labels are automatically configured when a
    `~pandas.Series`, `~pandas.DataFrame`, `~xarray.DataArray`, or `~pint.Quantity`
    is passed to the plotting command. Default is :rc:`autoformat`.
"""
_args_2d_shared_docstring = """
%(plot.args_1d_shared)s
order : {{'C', 'F'}}, optional
    If ``'C'`` (C-style row-major order), `z` coordinates should be shaped
    ``(y, x)``. If ``'F'`` (Fortran-style column-major order) `z` coordinates
    should be shaped ``(x, y)``. Default is ``'C'``.
globe : bool, optional
    For `proplot.axes.GeoAxes` only. Whether to enforce global coverage.
    Default is ``False``. When set to ``True`` this does the following:

    #. Interpolates input data to the North and South poles by setting the data
       values at the poles to the mean from latitudes nearest each pole.
    #. Makes meridional coverage "circular", i.e. the last longitude coordinate
       equals the first longitude coordinate plus 360\N{DEGREE SIGN}.
    #. When basemap is the backend, cycles 1D longitude vectors to fit within
       the map edges. For example, if the central longitude is 90\N{DEGREE SIGN},
       the data is shifted so that it spans -90\N{DEGREE SIGN} to 270\N{DEGREE SIGN}.
"""
_snippet_manager['plot.args_1d_shared'] = _args_1d_shared_docstring
_snippet_manager['plot.args_2d_shared'] = _args_2d_shared_docstring


# Auto colorbar and legend docstring
_guide_docstring = """
colorbar : bool, int, or str, optional
    If not ``None``, this is a location specifying where to draw an
    *inset* or *panel* colorbar from the resulting object(s). If ``True``,
    the default :rc:`colorbar.loc` is used. Valid locations are described
    in `~proplot.axes.Axes.colorbar`.
colorbar_kw : dict-like, optional
    Extra keyword args for the call to `~proplot.axes.Axes.colorbar`.
legend : bool, int, or str, optional
    If not ``None``, this is a location specifying where to draw an *inset*
    or *panel* legend from the resulting object(s). If ``True``, the default
    :rc:`legend.loc` is used. Valid locations are described in
    `~proplot.axes.Axes.legend`.
legend_kw : dict-like, optional
    Extra keyword args for the call to `~proplot.axes.Axes.legend`.
"""
_snippet_manager['plot.guide'] = _guide_docstring


# Misc shared 1D plotting docstrings
_inbounds_docstring = """
inbounds : bool, optional
    Whether to restrict the default `y` (`x`) axis limits to account for only
    in-bounds data when the `x` (`y`) axis limits have been locked. Default
    is :rc:`axes.inbounds`. See also :rcraw:`cmap.inbounds`.
"""
_error_means_docstring = """
mean, means : bool, optional
    Whether to plot the means of each column for 2D `{y}` coordinates. Means
    are calculated with `numpy.nanmean`. If no other arguments are specified,
    this also sets ``barstd=True`` (and ``boxstd=True`` for violin plots).
median, medians : bool, optional
    Whether to plot the medians of each column for 2D `{y}` coordinates. Medians
    are calculated with `numpy.nanmedian`. If no other arguments arguments are
    specified, this also sets ``barstd=True`` (and ``boxstd=True`` for violin plots).
"""
_error_bars_docstring = """
barstd, barstds : bool, float, or 2-tuple of float, optional
    *Valid only if `mean` or `median` is ``True``*. Standard deviation multiples for
    *thin error bars* with optional whiskers (i.e., caps). If scalar, then +/- that
    multiple is used. If ``True``, the default standard deviation range of +/-3 is used.
barpctile, barpctiles : bool, float, or 2-tuple of float, optional
    *Valid only if `mean` or `median` is ``True``*. As with `barstd`, but instead
    using *percentiles* for the error bars. If scalar, that percentile range is
    used (e.g., ``90`` shows the 5th to 95th percentiles). If ``True``, the default
    percentile range of 0 to 100 is used.
bardata : 2 x N array or 1D array, optional
    *Valid only if `mean` and `median` are ``False``*. If shape is 2 x N, these
    are the lower and upper bounds for the thin error bars. If shape is N, these
    are the absolute, symmetric deviations from the central points.
boxstd, boxstds, boxpctile, boxpctiles, boxdata : optional
    As with `barstd`, `barpctile`, and `bardata`, but for *thicker error bars*
    representing a smaller interval than the thin error bars. If `boxstd` is
    ``True``, the default standard deviation range of +/-1 is used. If `boxpctiles`
    is ``True``, the default percentile range of 25 to 75 is used (i.e., the
    interquartile range). When "boxes" and "bars" are combined, this has the
    effect of drawing miniature box-and-whisker plots.
capsize : float, optional
    The cap size for thin error bars in points. Default is :rc:`errorbar.capsize`.
barz, barzorder, boxz, boxzorder : float, optional
    The "zorder" for the thin and thick error bars. Default is ``2.5``.
barc, barcolor, boxc, boxcolor : color-spec, optional
    Colors for the thin and thick error bars. Default is
    :rc:`boxplot.whiskerprops.color`.
barlw, barlinewidth, boxlw, boxlinewidth : float, optional
    Line widths for the thin and thick error bars, in points. The defaults
    :rc:`boxplot.whiskerprops.linewidth` (bars) and four times that value (boxes).
boxm, boxmarker : bool or marker-spec, optional
    Whether to draw a small marker in the middle of the box denoting the mean or
    median position. Ignored if `boxes` is ``False``. Default is ``'o'``.
boxms, boxmarkersize : size-spec, optional
    The marker size for the `boxmarker` marker in points ** 2. Default size
    is equal to ``(2 * boxlinewidth) ** 2``.
boxmc, boxmarkercolor, boxmec, boxmarkeredgecolor : color-spec, optional
    Color, face color, and edge color for the `boxmarker` marker. Default color
    and edge color are ``'w'``.
"""
_error_shading_docstring = """
shadestd, shadestds, shadepctile, shadepctiles, shadedata : optional
    As with `barstd`, `barpctile`, and `bardata`, but using *shading* to indicate
    the error range. If `shadestds` is ``True``, the default standard deviation
    range of +/-2 is used. If `shadepctiles` is ``True``, the default
    percentile range of 10 to 90 is used.
fadestd, fadestds, fadepctile, fadepctiles, fadedata : optional
    As with `shadestd`, `shadepctile`, and `shadedata`, but for an additional,
    more faded, *secondary* shaded region. If `fadestds` is ``True``, the default
    standard deviation range of +/-3 is used. If `fadepctiles` is ``True``,
    the default percentile range of 0 to 100 is used.
shadec, shadecolor, fadec, fadecolor : color-spec, optional
    Colors for the different shaded regions. Default is to inherit the parent color.
shadez, shadezorder, fadez, fadezorder : float, optional
    The "zorder" for the different shaded regions. Default is ``1.5``.
shadea, shadealpha, fadea, fadealpha : float, optional
    The opacity for the different shaded regions. Defaults are ``0.4`` and ``0.2``.
shadelw, shadelinewidth, fadelw, fadelinewidth : float, optional
    The edge line width for the shading patches. Default is :rc:`patch.linewidth`.
shdeec, shadeedgecolor, fadeec, fadeedgecolor : float, optional
    The edge color for the shading patches. Default is ``'none'``.
shadelabel, fadelabel : bool or str, optional
    Labels for the shaded regions to be used as separate legend entries. To toggle
    labels "on" and apply a *default* label, use e.g. ``shadelabel=True``. To apply
    a *custom* label, use e.g. ``shadelabel='label'``. Otherwise, the shading is
    drawn underneath the line and/or marker in the legend entry.
"""
_snippet_manager['plot.inbounds'] = _inbounds_docstring
_snippet_manager['plot.error_means_y'] = _error_means_docstring.format(y='y')
_snippet_manager['plot.error_means_x'] = _error_means_docstring.format(y='x')
_snippet_manager['plot.error_bars'] = _error_bars_docstring
_snippet_manager['plot.error_shading'] = _error_shading_docstring


# Color docstrings
_cycle_docstring = """
cycle : cycle-spec, optional
    The cycle specifer, passed to the `~proplot.constructor.Cycle` constructor.
    If the returned cycler is unchanged from the current cycler, the axes
    cycler will not be reset to its first position. To disable property cycling
    and just use black for the default color, use ``cycle=False``, ``cycle='none'``,
    or ``cycle=()`` (analogous to disabling ticks with e.g. ``xformatter='none'``).
    To restore the default property cycler, use ``cycle=True``.
cycle_kw : dict-like, optional
    Passed to `~proplot.constructor.Cycle`.
"""
_cmap_norm_docstring = """
cmap : colormap spec, optional
    The colormap specifer, passed to the `~proplot.constructor.Colormap`
    constructor function.
cmap_kw : dict-like, optional
    Passed to `~proplot.constructor.Colormap`.
norm : normalizer spec, optional
    The continuous colormap normalizer, passed to the `~proplot.constructor.Norm`
    constructor function. If `discrete` is ``True`` this is also used to normalize
    values passed to `~proplot.colors.DiscreteNorm` before colors is selected.
norm_kw : dict-like, optional
    Passed to `~proplot.constructor.Norm`.
discrete : bool, optional
    If ``False``, then `~proplot.colors.DiscreteNorm` is not applied to the
    colormap. Instead, for non-contour plots, the number of levels will be
    roughly controlled by :rcraw:`cmap.lut`. This has a similar effect to
    using `levels=large_number` but it may improve rendering speed. Default
    is ``False`` for `~proplot.axes.Axes.imshow`, `~proplot.axes.Axes.matshow`,
    `~proplot.axes.Axes.spy`, `~proplot.axes.Axes.hexbin`, `~proplot.axes.Axes.hist2d`,
    and `~proplot.axes.Axes.heatmap` plots, but ``True`` otherwise.
sequential : bool, optional
    Use :rc:`cmap.sequential` as the default colormap.
diverging : bool, optional
    Use :rc:`cmap.diverging` as the default colormap and use
    `~proplot.colors.DivergingNorm` as the default continuous normalizer.
    This will also ensure auto-generated levels include a value at zero.
cyclic : bool, optional
    Use :rc:`cmap.cyclic` as the default colormap and modify the default
    arguments passed to `~proplot.colors.DiscreteNorm` so that colors
    on either end are distinct.
sequential, diverging, cyclic, qualitative : bool, optional
    Boolean arguments used if `cmap` is not passed. Set these to ``True``
    to use the default :rcraw:`cmap.sequential`, :rcraw:`cmap.diverging`,
    :rcraw:`cmap.cyclic`, and :rcraw:`cmap.qualitative` colormaps. The
    latter three options also change level- and norm-generation behavior.
extend : {{'neither', 'min', 'max', 'both'}}, optional
    Whether to assign unique colors to out-of-bounds data and draw
    colorbar "extensions" when a colorbar is drawn.
"""
_snippet_manager['plot.cycle'] = _cycle_docstring
_snippet_manager['plot.cmap_norm'] = _cmap_norm_docstring


# Levels docstrings
# NOTE: In some functions we only need some components
_vlim_levels_docstring = """
vmin, vmax : float, optional
    Used to determine level locations if `levels` or `values` is an integer.
    Actual levels may not fall exactly on `vmin` and `vmax`, but the minimum
    level will be no smaller than `vmin` and the maximum level will be
    no larger than `vmax`. If `vmin` or `vmax` are not provided, the
    minimum and maximum data values are used.
"""
_manual_levels_docstring = """
N
    Shorthand for `levels`.
levels : int or list of float, optional
    The number of level edges or a list of level edges. If the former,
    `locator` is used to generate this many level edges at "nice" intervals.
    If the latter, the levels should be monotonically increasing or
    decreasing (note that decreasing levels will only work with ``pcolor``
    plots, not ``contour`` plots). Default is :rc:`cmap.levels`.
values : int or list of float, optional
    The number of level centers or a list of level centers. If the former,
    `locator` is used to generate this many level centers at "nice" intervals.
    If the latter, levels are inferred using `~proplot.utils.edges`.
    This will override any `levels` input.
"""
_auto_levels_docstring = """
robust : bool, float, or 2-tuple, optional
    If ``True`` and `vmin` or `vmax` were not provided, they are
    determined from the 2nd and 98th data percentiles rather than the
    minimum and maximum. If float, this percentile range is used (for example,
    ``90`` corresponds to the 5th to 95th percentiles). If 2-tuple of float,
    these specific percentiles should be used. This feature is useful
    when your data has large outliers. Default is :rc:`cmap.robust`.
inbounds : bool, optional
    If ``True`` and `vmin` or `vmax` were not provided, when axis limits
    have been explicitly restricted with `~matplotlib.axes.Axes.set_xlim`
    or `~matplotlib.axes.Axes.set_ylim`, out-of-bounds data is ignored.
    Default is :rc:`cmap.inbounds`. See also :rcraw:`axes.inbounds`.
locator : locator-spec, optional
    The locator used to determine level locations if `levels` or `values`
    is an integer. Passed to the `~proplot.constructor.Locator` constructor.
    Default is `~matplotlib.ticker.MaxNLocator` with ``levels`` integer levels.
locator_kw : dict-like, optional
    Passed to `~proplot.constructor.Locator`.
symmetric : bool, optional
    If ``True``, automatically generated levels are symmetric about zero.
    Default is always ``False``.
positive : bool, optional
    If ``True``, automatically generated levels are positive with a minimum at zero.
    Default is always ``False``.
negative : bool, optional
    If ``True``, automatically generated levels are negative with a maximum at zero.
    Default is always ``False``.
nozero : bool, optional
    If ``True``, ``0`` is removed from the level list. This is mainly useful for
    single-color `~matplotlib.axes.Axes.contour` plots.
"""
_snippet_manager['plot.levels_vlim'] = _vlim_levels_docstring
_snippet_manager['plot.levels_manual'] = _manual_levels_docstring
_snippet_manager['plot.levels_auto'] = _auto_levels_docstring


# Labels docstrings
_labels_1d_docstring = """
label, value : float or str, optional
    The single legend label or colorbar coordinate to be used for this plotted
    element. This is generally used with 1D input coordinates.
labels, values : list of float or list of str, optional
    The legend labels or colorbar coordinates used for each plotted element.
    Can be numeric or string, and must match the number of plotted elements.
    This is generally used with 2D input coordinates.
"""
_labels_2d_docstring = """
labels : bool, optional
    Whether to apply labels to contours and grid boxes. The text will be
    white when the luminance of the underlying filled contour or grid box
    is less than 50 and black otherwise.
labels_kw : dict-like, optional
    Ignored if `labels` is ``False``. Extra keyword args for the labels.
    For contour plots, this is passed to `~matplotlib.axes.Axes.clabel`.
    Otherwise, this is passed to `~matplotlib.axes.Axes.text`.
fmt : format-spec, optional
    Passed to the `~proplot.constructor.Norm` constructor, used to format
    number labels. You can also use the `precision` keyword arg.
precision : int, optional
    Maximum number of decimal places for the number labels. Number labels
    are generated with the `~proplot.ticker.SimpleFormatter` formatter,
    which permits limiting the precision.
label : str, optional
    The legend label to be used for this object. In the case of
    contours, this is paired with the the central artist in the artist
    list returned by `matplotlib.contour.ContourSet.legend_elements`.
"""
_snippet_manager['plot.labels_1d'] = _labels_1d_docstring
_snippet_manager['plot.labels_2d'] = _labels_2d_docstring


# Plot docstring
_plot_docstring = """
Plot standard lines.

Parameters
----------
%(plot.args_1d_{y})s
%(plot.args_1d_shared)s

Other parameters
----------------
%(plot.cycle)s
%(plot.labels_1d)s
%(plot.guide)s
%(plot.error_means_{y})s
%(plot.error_bars)s
%(plot.error_shading)s
%(plot.inbounds)s
**kwargs
    Passed to `~matplotlib.axes.Axes.plot`.

See also
--------
PlotAxes.plot
PlotAxes.plotx
matplotlib.axes.Axes.plot
"""
_snippet_manager['plot.plot'] = _plot_docstring.format(y='y')
_snippet_manager['plot.plotx'] = _plot_docstring.format(y='x')


# Step docstring
# NOTE: Internally matplotlib implements step with thin wrapper of plot
_step_docstring = """
Plot step lines.

Parameters
----------
%(plot.args_1d_{y})s
%(plot.args_1d_shared)s

Other parameters
----------------
%(plot.cycle)s
%(plot.labels_1d)s
%(plot.guide)s
%(plot.inbounds)s
**kwargs
    Passed to `~matplotlib.axes.Axes.step`.

See also
--------
PlotAxes.step
PlotAxes.stepx
matplotlib.axes.Axes.step
"""
_snippet_manager['plot.step'] = _step_docstring.format(y='y')
_snippet_manager['plot.stepx'] = _step_docstring.format(y='x')


# Stem docstring
_stem_docstring = """
Plot stem lines.

Parameters
----------
%(plot.args_1d_{y})s
%(plot.args_1d_shared)s

Other parameters
----------------
%(plot.cycle)s
%(plot.guide)s
%(plot.inbounds)s
**kwargs
    Passed to `~matplotlib.axes.Axes.stem`.
"""
_snippet_manager['plot.stem'] = _stem_docstring.format(y='x')
_snippet_manager['plot.stemx'] = _stem_docstring.format(y='x')


# Lines docstrings
_lines_docstring = """
Plot {orientation} lines.

Parameters
----------
%(plot.args_1d_multi{y})s
%(plot.args_1d_shared)s

Other parameters
----------------
stack, stacked : bool, optional
    Whether to "stack" successive columns of the `{y}1` coordinates. If this
    is ``True`` and `{y}2` was provided, it will be ignored.
negpos : bool, optional
    Whether to color lines greater than zero with `poscolor` and lines less
    than zero with `negcolor`.
negcolor, poscolor : color-spec, optional
    Colors to use for the negative and positive lines. Ignored if `negpos`
    is ``False``. Defaults are :rc:`negcolor` and :rc:`poscolor`.
c, color, colors : color-spec or list thereof, optional
    The line color(s).
ls, linestyle, linestyles : linestyle-spec or list thereof, optional
    The line style(s).
lw, linewidth, linewidths : linewidth-spec or list thereof, optional
    The line width(s).
%(plot.cycle)s
%(plot.labels_1d)s
%(plot.guide)s
%(plot.inbounds)s
**kwargs
    Passed to `~matplotlib.axes.Axes.{prefix}lines`.

See also
--------
PlotAxes.vlines
PlotAxes.hlines
matplotlib.axes.Axes.vlines
matplotlib.axes.Axes.hlines
"""
_snippet_manager['plot.vlines'] = _lines_docstring.format(
    y='y', prefix='v', orientation='vertical'
)
_snippet_manager['plot.hlines'] = _lines_docstring.format(
    y='x', prefix='h', orientation='horizontal'
)


# Scatter docstring
_parametric_docstring = """
Plot a parametric line.

Parameters
----------
%(plot.args_1d_y)s
c, color, colors, values : array-like, optional
    The parametric coordinate. These can be passed as a third
    positional argument or as a keyword argument. They can also
    be string labels instead of numbers and the resulting
    colorbar ticks will be labeled accordingly.
%(plot.args_1d_shared)s

Other parameters
----------------
%(plot.cmap_norm)s
%(plot.levels_vlim)s
%(plot.guide)s
interp : int, optional
    Interpolate to this many additional points between the parametric
    coordinates. Default is ``0``. This can be increased to make the color
    gradations between a small number of coordinates appear "smooth".
scalex, scaley : bool, optional
    Whether the view limits are adapted to the data limits. The values are
    passed on to `~matplotlib.axes.Axes.autoscale_view`.
%(plot.inbounds)s
**kwargs
    Valid `~matplotlib.collections.LineCollection` properties.

Returns
-------
`~matplotlib.collections.LineCollection`
    The parametric line. See `this matplotlib example \
<https://matplotlib.org/stable/gallery/lines_bars_and_markers/multicolored_line>`__.

See also
--------
PlotAxes.plot
PlotAxes.plotx
matplotlib.collections.LineCollection
"""
_snippet_manager['plot.parametric'] = _parametric_docstring


# Scatter function docstring
_scatter_docstring = """
Plot markers with flexible keyword arguments.

Parameters
----------
%(plot.args_1d_{y})s
s, size, ms, markersize : float or list of float, optional
    The marker size(s). If this is an array matching the shape of `x` and `y`,
    the units are scaled by `smin` and `smax`.
c, color, colors, mc, markercolor, markercolors \
: color-spec or list thereof, or array, optional
    The marker fill color(s). If this is an array matching the shape of `x` and `y`,
    the colors are generated using `cmap`, `norm`, `vmin`, and `vmax`.
smin, smax : float, optional
    The minimum and maximum marker size in units ``points**2`` used to scale
    `s`. If not provided, the marker sizes are equivalent to the values in `s`.
%(plot.levels_vlim)s
%(plot.args_1d_shared)s

Other parameters
----------------
lw, linewidth, linewidths, markeredgewidth, markeredgewidths \
: float or list thereof, optional
    The marker edge width.
edgecolors, markeredgecolor, markeredgecolors \
: color-spec or list thereof, optional
    The marker edge color.
%(plot.cmap_norm)s
%(plot.levels_manual)s
%(plot.levels_auto)s
%(plot.cycle)s
%(plot.labels_1d)s
%(plot.guide)s
%(plot.error_means_{y})s
%(plot.error_bars)s
%(plot.error_shading)s
%(plot.inbounds)s
**kwargs
    Passed to `~matplotlib.axes.Axes.scatter`.

See also
--------
PlotAxes.scatter
PlotAxes.scatterx
matplotlib.axes.Axes.scatter
"""
_snippet_manager['plot.scatter'] = _scatter_docstring.format(y='y')
_snippet_manager['plot.scatterx'] = _scatter_docstring.format(y='x')


# Bar function docstring
_bar_docstring = """
Plot individual, grouped, or stacked bars.

Parameters
----------
%(plot.args_1d_{y})s
width : float or array-like, optional
    The width(s) of the bars relative to the {x} coordinate step size.
    Can be passed as a third positional argument.
{bottom} : float or array-like, optional
    The coordinate(s) of the {bottom} edge of the bars. Default is
    ``0``. Can be passed as a fourth positinal argument.
absolute_width : bool, optional
    Whether to make the `width` units *absolute*. If ``True``, this
    restores the default matplotlib behavior. Default is ``False``.
%(plot.args_1d_shared)s

Other parameters
----------------
stack, stacked : bool, optional
    Whether to stack columns of the input array or plot the bars
    side-by-side in groups.
negpos : bool, optional
    Whether to shade bars greater than zero with `poscolor` and bars less
    than zero with `negcolor`.
negcolor, poscolor : color-spec, optional
    Colors to use for the negative and positive bars. Ignored if `negpos`
    is ``False``. Defaults are :rc:`negcolor` and :rc:`poscolor`.
lw, linewidth, linewidths : float, optional
    The edge width for the bar patches.
ec, edgecolor, edgecolors : color-spec, optional
    The edge color for the bar patches.
%(plot.cycle)s
%(plot.labels_1d)s
%(plot.guide)s
%(plot.error_means_{y})s
%(plot.error_bars)s
%(plot.inbounds)s
**kwargs
    Passed to `~matplotlib.axes.Axes.bar{suffix}`.

See also
--------
PlotAxes.bar
PlotAxes.barh
matplotlib.axes.Axes.bar
matplotlib.axes.Axes.barh
"""
_snippet_manager['plot.bar'] = _bar_docstring.format(
    x='x', y='y', bottom='bottom', suffix=''
)
_snippet_manager['plot.barh'] = _bar_docstring.format(
    x='y', y='x', bottom='left', suffix='h'
)


# Area plot docstring
_fill_docstring = """
Plot individual, grouped, or overlaid shading patches.

Parameters
----------
%(plot.args_1d_multi{y})s
%(plot.args_1d_shared)s

Other parameters
----------------
stack, stacked : bool, optional
    Whether to "stack" successive columns of the `{y}1` array. If this is
    ``True`` and `{y}2` was provided, it will be ignored.
negpos : bool, optional
    Whether to shade where ``{y}1 >= {y}2`` with `poscolor` and where ``{y}1 < {y}2``
    with `negcolor`. For example, to shade positive values red and negative values
    blue, simply use ``ax.fill_between{suffix}({x}, {y}, negpos=True)``.
negcolor, poscolor : color-spec, optional
    Colors to use for the negative and positive shaded regions. Ignored if `negpos`
    is ``False``. Defaults are :rc:`negcolor` and :rc:`poscolor`.
where : ndarray, optional
    Boolean ndarray mask for points you want to shade. See `this example \
<https://matplotlib.org/stable/gallery/pyplots/whats_new_98_4_fill_between.html>`__.
lw, linewidth, linewidths : float, optional
    The edge width for the area patches.
ec, edgecolor, edgecolors : color-spec, optional
    The edge color for the area patches.
%(plot.cycle)s
%(plot.labels_1d)s
%(plot.guide)s
%(plot.inbounds)s
**kwargs
    Passed to `~matplotlib.axes.Axes.fill_between{suffix}`.

See also
--------
PlotAxes.area
PlotAxes.areax
PlotAxes.fill_between
PlotAxes.fill_betweenx
matplotlib.axes.Axes.fill_between
matplotlib.axes.Axes.fill_betweenx
"""
_snippet_manager['plot.fill_between'] = _fill_docstring.format(
    x='x', y='y', suffix=''
)
_snippet_manager['plot.fill_betweenx'] = _fill_docstring.format(
    x='y', y='x', suffix='x'
)


# Histogram docstrings
_weight_docstring = """
weights : array-like, optional
    The weights associated with each point. If string this
    can be retrieved from `data` (see below).
"""
_snippet_manager['plot.weights'] = _weight_docstring
_hist_docstring = """
Plot {orientation} histograms.

Parameters
----------
%(plot.args_1d_{y})s
bins : int or list of float, optional
    The bin count or list of bins.
%(plot.weights)s
%(plot.args_1d_shared)s

Other parameters
----------------
%(plot.cycle)s
%(plot.labels_1d)s
%(plot.guide)s
**kwargs
    Passed to `~matplotlib.axes.Axes.hist`.

See also
--------
PlotAxes.hist
PlotAxes.histh
matplotlib.axes.Axes.hist
"""
_snippet_manager['plot.hist'] = _hist_docstring.format(
    y='x', orientation='vertical'
)
_snippet_manager['plot.histh'] = _hist_docstring.format(
    y='x', orientation='horizontal'
)


# Box plot docstrings
_boxplot_docstring = """
Plot {orientation} boxes and whiskers with a nice default style.

Parameters
----------
%(plot.args_1d_{y})s
%(plot.args_1d_shared)s

Other parameters
----------------
mean, means : bool, optional
    If ``True``, this passes ``showmeans=True`` and ``meanline=True`` to
    `~matplotlib.axes.Axes.boxplot`.
fill : bool, optional
    Whether to fill the box with a color.
fc, facecolor, fillcolor : color-spec, list, optional
    The fill color for the boxes. Default is the next color cycler color. If
    a list, it should be the same length as the number of objects.
a, alpha, fa, facealpha, fillalpha : float, optional
    The opacity of the boxes. Default is ``1.0``. If a list,
    should be the same length as the number of objects.
lw, linewidth, linewidths : float, optional
    The linewidth of all objects. Default is :rc:`patch.linewidth`.
c, color, colors, ec, edgecolor, edgecolors : color-spec, list, optional
    The color of all objects. Default is ``'black'``. If a list, it should
    be the same length as the number of objects.
meanls, medianls, meanlinestyle, medianlinestyle, meanlinestyles, medianlinestyles \
: line style-spec, optional
    The line style for the mean and median lines drawn horizontally
    across the box.
boxc, capc, whiskerc, flierc, meanc, medianc, \
boxcolor, capcolor, whiskercolor, fliercolor, meancolor, mediancolor \
boxcolors, capcolors, whiskercolors, fliercolors, meancolors, mediancolors \
: color-spec, list, optional
    The color of various boxplot components. If a list, it should be the
    same length as the number of objects. These are shorthands so you don't
    have to pass e.g. a ``boxprops`` dictionary.
boxlw, caplw, whiskerlw, flierlw, meanlw, medianlw, boxlinewidth, caplinewidth, \
meanlinewidth, medianlinewidth, whiskerlinewidth, flierlinewidth, boxlinewidths, \
caplinewidths, meanlinewidths, medianlinewidths, whiskerlinewidths, flierlinewidths \
: float, optional
    The line width of various boxplot components. These are shorthands so
    you don't have to pass e.g. a ``boxprops`` dictionary.
m, marker : marker-spec, optional
    Marker style for the 'fliers', i.e. outliers.
ms, markersize : float, optional
    Marker size for the 'fliers', i.e. outliers.
%(plot.cycle)s
%(plot.labels_1d)s
**kwargs
    Passed to `matplotlib.axes.Axes.boxplot`.

See also
--------
PlotAxes.boxes
PlotAxes.boxesh
PlotAxes.boxplot
PlotAxes.boxploth
matplotlib.axes.Axes.boxplot
"""
_snippet_manager['plot.boxplot'] = _boxplot_docstring.format(
    y='y', orientation='vertical'
)
_snippet_manager['plot.boxploth'] = _boxplot_docstring.format(
    y='x', orientation='horizontal'
)


# Violin plot docstrings
_violinplot_docstring = """
Plot {orientation} violins with a nice default style matching
`this matplotlib example \
<https://matplotlib.org/stable/gallery/statistics/customized_violin.html>`__.

Parameters
----------
%(plot.args_1d_{y})s
%(plot.args_1d_shared)s

Other parameters
----------------
fc, facecolor, facecolors, fillcolor, fillcolors : color-spec, list, optional
    The violin plot fill color. Default is the next color cycler color. If
    a list, it should be the same length as the number of objects.
a, alpha, fa, facealpha, fillalpha : float, optional
    The opacity of the violins. Default is ``1.0``. If a list,
    it should be the same length as the number of objects.
lw, linewidth, linewidths : float, optional
    The linewidth of the line objects. Default is :rc:`patch.linewidth`.
c, color, colors, ec, edgecolor, edgecolors : color-spec, list, optional
    The edge color for the violin patches. Default is ``'black'``. If a
    list, it should be the same length as the number of objects.
%(plot.cycle)s
%(plot.labels_1d)s
%(plot.error_bars)s
**kwargs
    Passed to `matplotlib.axes.Axes.violinplot`.

Note
----
It is no longer possible to show minima and maxima with whiskers --
while this is useful for `~matplotlib.axes.Axes.boxplot`\\ s it is
redundant for `~matplotlib.axes.Axes.violinplot`\\ s.

See also
--------
PlotAxes.violins
PlotAxes.violinsh
PlotAxes.violinplot
PlotAxes.violinploth
matplotlib.axes.Axes.violinplot
"""
_snippet_manager['plot.violinplot'] = _violinplot_docstring.format(
    y='y', orientation='vertical'
)
_snippet_manager['plot.violinploth'] = _violinplot_docstring.format(
    y='x', orientation='horizontal'
)


# Contour docstrings
_contour_docstring = """
Plot {descrip}.

Parameters
----------
%(plot.args_2d)s

%(plot.args_2d_shared)s

Other parameters
----------------
%(plot.cmap_norm)s
%(plot.levels_manual)s
%(plot.levels_vlim)s
%(plot.levels_auto)s
%(plot.labels_2d)s
%(plot.guide)s
{add}lw, linewidth, linewidths : optional
    The width of the contour lines.
    For `contourf` plots, lines are added between the filled contours.
ls, linestyle, linestyles : optional
    The style of the contour lines.
    For `contourf` plots, lines are added between the filled contours.
ec, edgecolor, edgecolors : optional
    The color for the contour lines.
    For `contourf` plots, lines are added between the filled contours.
c, color, colors : optional
    The color(s) for the contour lines or filled contours. If not passed,
    the color is determined by `cmap` and the `z` data.
**kwargs
    Passed to `matplotlib.axes.Axes.{command}`.

See also
--------
PlotAxes.contour
PlotAxes.contourf
PlotAxes.tricontour
PlotAxes.tricontourf
matplotlib.axes.Axes.{command}
"""
_edgefix_docstring = """
edgefix : bool or float, optional
    Whether to fix an issue where `white lines appear between filled contours
    <https://stackoverflow.com/q/8263769/4970632>`__ in saved vector graphics.
    This can slow down figure rendering. Default is :rc:`cmap.edgefix`.
    If ``True``, a default linewidth is used. If float, this linewidth is used.
"""
_snippet_manager['plot.contour'] = _contour_docstring.format(
    descrip='contour lines', command='contour', add=''
)
_snippet_manager['plot.contourf'] = _contour_docstring.format(
    descrip='filled contours', command='contourf', add=_edgefix_docstring
)
_snippet_manager['plot.tricontour'] = _contour_docstring.format(
    descrip='contour lines on a triangular grid', command='tricontour', add=''
)
_snippet_manager['plot.tricontourf'] = _contour_docstring.format(
    descrip='filled contours on a triangular grid', command='tricontourf', add=_edgefix_docstring  # noqa: E501
)


# Pcolor docstring
_pcolor_docstring = """
Plot {descrip}.

Parameters
----------
%(plot.args_2d)s

%(plot.args_2d_shared)s
{add}

Other parameters
----------------
%(plot.cmap_norm)s
%(plot.levels_manual)s
%(plot.levels_vlim)s
%(plot.levels_auto)s
%(plot.labels_2d)s
%(plot.guide)s
lw, linewidth, linewidths : optional
    The width of lines between grid boxes.
ls, linestyle, linestyles : optional
    The style of lines between grid boxes.
ec, edgecolor, edgecolors : optional
    The color for lines between grid boxes.
c, color, colors : optional
    The color(s) for the grid boxes. If not passed,
    the color is determined by `cmap` and the `z` data.
edgefix : bool, optional
    Whether to fix an issue where `white lines appear between grid boxes
    <https://stackoverflow.com/q/8263769/4970632>`__ in saved vector graphics.
    This can slow down figure rendering. Default is :rc:`cmap.edgefix`.
    If ``True``, a default linewidth is used. If float, this linewidth is used.
**kwargs
    Passed to `matplotlib.axes.Axes.{command}`.

See also
--------
PlotAxes.pcolor
PlotAxes.pcolormesh
PlotAxes.pcolorfast
PlotAxes.heatmap
PlotAxes.tripcolor
matplotlib.axes.Axes.{command}
"""
_heatmap_descrip = """
grid boxes with formatting suitable for heatmaps. Ensures square grid
boxes, adds major ticks to the center of each grid box, disables minor ticks
and gridlines, and sets :rcraw:`cmap.discrete` to ``False`` by default.
"""
_heatmap_aspect = """
aspect : {'equal', 'auto'} or float, optional
    Modify the axes aspect ratio. The aspect ratio is of particular
    relevance for heatmaps since it may lead to non-square grid boxes.
    This parameter is a shortcut for calling `~matplotlib.axes.set_aspect`.
    Default is :rc:`image.aspect`. The options are as follows:

    * Number: The data aspect ratio.
    * ``'equal'``: A data aspect ratio of 1.
    * ``'auto'``: Allows the data aspect ratio to change depending on
      the layout. In general this results in non-square grid boxes.
"""
_snippet_manager['plot.pcolor'] = _pcolor_docstring.format(
    descrip='irregular grid boxes', command='pcolor', add=''
)
_snippet_manager['plot.pcolormesh'] = _pcolor_docstring.format(
    descrip='regular grid boxes', command='pcolormesh', add=''
)
_snippet_manager['plot.pcolorfast'] = _pcolor_docstring.format(
    descrip='grid boxes quickly', command='pcolorfast', add=''
)
_snippet_manager['plot.tripcolor'] = _pcolor_docstring.format(
    descrip='triangular grid boxes', command='tripcolor', add=''
)
_snippet_manager['plot.heatmap'] = _pcolor_docstring.format(
    descrip=_heatmap_descrip.strip(), command='pcolormesh', add=_heatmap_aspect
)


# Flow function docstring
_flow_docstring = """
Plot {descrip}.

Parameters
----------
%(plot.args_2d_flow)s

c, color, colors : array-like or color-spec, optional
    The colors of the {descrip} passed as either a keyword argument
    or a fifth positional argument. This can be a single color or
    a color array to be scaled by `cmap` and `norm`.
%(plot.args_2d_shared)s

Other parameters
----------------
%(plot.cmap_norm)s
%(plot.levels_manual)s
%(plot.levels_vlim)s
%(plot.levels_auto)s
**kwargs
    Passed to `matplotlib.axes.Axes.{command}`

See also
--------
PlotAxes.barbs
PlotAxes.quiver
PlotAxes.stream
PlotAxes.streamplot
matplotlib.axes.Axes.{command}
"""
_snippet_manager['plot.barbs'] = _flow_docstring.format(
    descrip='wind barbs', command='barbs'
)
_snippet_manager['plot.quiver'] = _flow_docstring.format(
    descrip='quiver arrows', command='quiver'
)
_snippet_manager['plot.stream'] = _flow_docstring.format(
    descrip='streamlines', command='streamplot'
)


# Image docstring
_show_docstring = """
Plot {descrip}.

Parameters
----------
z : array-like
    The data passed as a positional argument or keyword argument.
%(plot.args_1d_shared)s

Other parameters
----------------
%(plot.cmap_norm)s
%(plot.levels_manual)s
%(plot.levels_vlim)s
%(plot.levels_auto)s
%(plot.guide)s
**kwargs
    Passed to `matplotlib.axes.Axes.{command}`.

See also
--------
proplot.axes.PlotAxes
matplotlib.axes.Axes.{command}
"""
_snippet_manager['plot.imshow'] = _show_docstring.format(
    descrip='an image', command='imshow'
)
_snippet_manager['plot.matshow'] = _show_docstring.format(
    descrip='a matrix', command='matshow'
)
_snippet_manager['plot.spy'] = _show_docstring.format(
    descrip='a sparcity pattern', command='spy'
)


def _load_objects():
    """
    Load array-like objects.
    """
    # NOTE: We just want to detect if *input arrays* belong to these types -- and if
    # this is the case, it means the module has already been imported! So, we only
    # try loading these classes within autoformat calls. This saves >500ms of import
    # time. We use ndarray as the default value for unimported types and in loops we
    # are careful to check membership to np.ndarray before anything else.
    global ndarray, DataArray, DataFrame, Series, Index, Quantity
    ndarray = np.ndarray
    DataArray = getattr(sys.modules.get('xarray', None), 'DataArray', ndarray)
    DataFrame = getattr(sys.modules.get('pandas', None), 'DataFrame', ndarray)
    Series = getattr(sys.modules.get('pandas', None), 'Series', ndarray)
    Index = getattr(sys.modules.get('pandas', None), 'Index', ndarray)
    Quantity = getattr(sys.modules.get('pint', None), 'Quantity', ndarray)


_load_objects()


# Standardization utilities
def _is_array(data):
    """
    Test whether input is numpy array or pint quantity.
    """
    # NOTE: This is used in _iter_columns to identify 2D matrices that
    # should be iterated over and omit e.g. scalar marker size or marker color.
    _load_objects()
    return isinstance(data, ndarray) or ndarray is not Quantity and isinstance(data, Quantity)  # noqa: E501


def _is_numeric(data):
    """
    Test whether input is numeric array rather than datetime or strings.
    """
    return len(data) and np.issubdtype(_to_numpy_array(data).dtype, np.number)


def _is_categorical(data):
    """
    Test whether input is array of strings.
    """
    return len(data) and isinstance(_to_numpy_array(data).item(0), str)


def _require_centers(x, y, z):
    """
    Enforce that coordinates are centers. Convert from edges if possible.
    """
    xlen, ylen = x.shape[-1], y.shape[0]
    if z.ndim == 2 and z.shape[1] == xlen - 1 and z.shape[0] == ylen - 1:
        # Get centers given edges
        if all(z.ndim == 1 and z.size > 1 and _is_numeric(z) for z in (x, y)):
            x = 0.5 * (x[1:] + x[:-1])
            y = 0.5 * (y[1:] + y[:-1])
        else:
            if x.ndim == 2 and x.shape[0] > 1 and x.shape[1] > 1 and _is_numeric(x):
                x = 0.25 * (x[:-1, :-1] + x[:-1, 1:] + x[1:, :-1] + x[1:, 1:])
            if y.ndim == 2 and y.shape[0] > 1 and y.shape[1] > 1 and _is_numeric(y):
                y = 0.25 * (y[:-1, :-1] + y[:-1, 1:] + y[1:, :-1] + y[1:, 1:])
    elif z.shape[-1] != xlen or z.shape[0] != ylen:
        # Helpful error message
        raise ValueError(
            f'Input shapes x {x.shape} and y {y.shape} '
            f'must match z centers {z.shape} '
            f'or z borders {tuple(i+1 for i in z.shape)}.'
        )
    return x, y


def _require_edges(x, y, z):
    """
    Enforce that coordinates are edges. Convert from centers if possible.
    """
    xlen, ylen = x.shape[-1], y.shape[0]
    if z.ndim == 2 and z.shape[1] == xlen and z.shape[0] == ylen:
        # Get edges given centers
        if all(z.ndim == 1 and z.size > 1 and _is_numeric(z) for z in (x, y)):
            x = edges(x)
            y = edges(y)
        else:
            if x.ndim == 2 and x.shape[0] > 1 and x.shape[1] > 1 and _is_numeric(x):
                x = edges2d(x)
            if y.ndim == 2 and y.shape[0] > 1 and y.shape[1] > 1 and _is_numeric(y):
                y = edges2d(y)
    elif z.shape[-1] != xlen - 1 or z.shape[0] != ylen - 1:
        # Helpful error message
        raise ValueError(
            f'Input shapes x {x.shape} and y {y.shape} must match '
            f'array centers {z.shape} or '
            f'array borders {tuple(i + 1 for i in z.shape)}.'
        )
    return x, y


def _safe_mask(mask, *args):
    """
    Safely apply the mask to the input arrays, accounting for existing masked
    or invalid values. Values matching ``False`` are set to `np.nan`.
    """
    _load_objects()
    invalid = ~mask  # True if invalid
    args_masked = []
    for arg in args:
        units = 1
        if ndarray is not Quantity and isinstance(arg, Quantity):
            arg, units = arg.magnitude, arg.units
        arg = ma.masked_invalid(arg, copy=False)
        arg = arg.astype(np.float).filled(np.nan)
        if arg.size > 1 and arg.shape != invalid.shape:
            raise ValueError(f'Mask shape {mask.shape} incompatible with array shape {arg.shape}.')  # noqa: E501
        if arg.size == 1 or invalid.size == 1:  # NOTE: happens with _restrict_inbounds
            pass
        elif invalid.size == 1:
            arg = np.nan if invalid.item() else arg
        elif arg.size > 1:
            arg[invalid] = np.nan
        args_masked.append(arg * units)
    return args_masked[0] if len(args_masked) == 1 else args_masked


def _safe_range(data, lo=0, hi=100, automin=True, automax=True):
    """
    Safely return the minimum and maximum (default) or percentile range accounting
    for masked values. Use min and max functions when possible for speed. Return
    ``None`` if we faile to get a valid range.
    """
    _load_objects()
    units = 1
    if ndarray is not Quantity and isinstance(data, Quantity):
        data, units = data.magnitude, data.units
    data = ma.masked_invalid(data, copy=False)
    data = data.compressed()  # remove all invalid values
    min_ = max_ = None
    if data.size and automin:
        min_ = float(np.min(data) if lo <= 0 else np.percentile(data, lo))
        if np.isfinite(min_):
            min_ *= units
        else:
            min_ = None
    if data.size and automax:
        max_ = float(np.max(data) if hi >= 100 else np.percentile(data, hi))
        if np.isfinite(max_):
            max_ *= units
        else:
            max_ = None
    return min_, max_


def _to_duck_array(data, strip_units=False):
    """
    Convert arbitrary input to duck array. Preserve array containers with metadata.
    """
    _load_objects()
    if data is None:
        raise ValueError('Invalid data None.')
    if (
        not isinstance(data, (ndarray, DataArray, DataFrame, Series, Index, Quantity))
        or not np.iterable(data)
    ):
        # WARNING: this strips e.g. scalar DataArray metadata
        data = _to_numpy_array(data)
    if strip_units:  # used for z coordinates that cannot have units
        if isinstance(data, (ndarray, Quantity)):
            if Quantity is not ndarray and isinstance(data, Quantity):
                data = data.magnitude
        elif isinstance(data, DataArray):
            if Quantity is not ndarray and isinstance(data.data, Quantity):
                data = data.copy(deep=False)
                data.data = data.data.magnitude
    return data


def _to_numpy_array(data, strip_units=False):
    """
    Convert arbitrary input to numpy array. Preserve masked arrays and unit arrays.
    """
    _load_objects()
    if data is None:
        raise ValueError('Invalid data None.')
    if isinstance(data, ndarray):
        pass
    elif isinstance(data, DataArray):
        data = data.data  # support pint quantities that get unit-stripped later
    elif isinstance(data, (DataFrame, Series, Index)):
        data = data.values
    if Quantity is not ndarray and isinstance(data, Quantity):
        if strip_units:
            return np.atleast_1d(data.magnitude)
        else:
            return np.atleast_1d(data.magnitude) * data.units
    else:
        return np.atleast_1d(data)  # natively preserves masked arrays


# Metadata utilities
def _get_data(data, *args):
    """
    Try to convert positional `key` arguments to `data[key]`. If argument is string
    it could be a valid positional argument like `fmt` so do not raise error.
    """
    if data is None:
        return
    args = list(args)
    for i, arg in enumerate(args):
        if isinstance(arg, str):
            try:
                array = data[arg]
            except KeyError:
                pass
            else:
                args[i] = array
    return args


def _get_coords(*args, which='x', **kwargs):
    """
    Return the index arrays associated with string coordinates and
    keyword arguments updated with index locators and formatters.
    """
    # NOTE: Why FixedLocator and not IndexLocator? The latter requires plotting
    # lines or else an error is raised... very strange.
    # NOTE: Why IndexFormatter and not FixedFormatter? The former ensures labels
    # correspond to indices while the latter can mysteriously truncate labels.
    res = []
    for arg in args:
        arg = _to_duck_array(arg)
        if not _is_categorical(arg):
            res.append(arg)
            continue
        if arg.ndim > 1:
            raise ValueError('Non-1D string coordinate input is unsupported.')
        idx = np.arange(len(arg))
        kwargs.setdefault(which + 'locator', mticker.FixedLocator(idx))
        kwargs.setdefault(which + 'formatter', pticker.IndexFormatter(_to_numpy_array(arg)))  # noqa: E501
        kwargs.setdefault(which + 'minorlocator', mticker.NullLocator())
        res.append(idx)
    return (*res, kwargs)


def _get_labels(data, axis=0, always=True):
    """
    Return the array-like "labels" along axis `axis`. If `always` is ``False``
    we return ``None`` for simple ndarray input.
    """
    # NOTE: Previously inferred 'axis 1' metadata of 1D variable using the
    # data values metadata but that is incorrect. The paradigm for 1D plots
    # is we have row coordinates representing x, data values representing y,
    # and column coordinates representing individual series.
    labels = None
    _load_objects()
    if axis not in (0, 1, 2):
        raise ValueError(f'Invalid axis {axis}.')
    if isinstance(data, (ndarray, Quantity)):
        if not always:
            pass
        elif axis < data.ndim:
            labels = np.arange(data.shape[axis])
        else:  # requesting 'axis 1' on a 1D array
            labels = np.array([0])
    # Xarray object
    # NOTE: Even if coords not present .coords[dim] auto-generates indices
    elif isinstance(data, DataArray):
        if axis < data.ndim:
            labels = data.coords[data.dims[axis]]
        elif not always:
            pass
        else:
            labels = np.array([0])
    # Pandas object
    elif isinstance(data, (DataFrame, Series, Index)):
        if axis == 0 and isinstance(data, (DataFrame, Series)):
            labels = data.index
        elif axis == 1 and isinstance(data, (DataFrame,)):
            labels = data.columns
        elif not always:
            pass
        else:  # beyond dimensionality
            labels = np.array([0])
    # Everything else
    # NOTE: Ensure data is at least 1D in _to_duck_array so this covers everything
    else:
        raise ValueError(f'Unrecognized array type {type(data)}.')
    return labels


def _get_title(data, include_units=True):
    """
    Return the "title" of an array-like object with metadata. Include units in
    the title if `include_units` is ``True``.
    """
    title = units = None
    _load_objects()
    if isinstance(data, ndarray):
        pass
    # Xarray object with possible long_name, standard_name, and units attributes.
    # Output depends on if units is True. Prefer long_name (come last in loop).
    elif isinstance(data, DataArray):
        title = getattr(data, 'name', None)
        for key in ('standard_name', 'long_name'):
            title = data.attrs.get(key, title)
        if include_units:
            units = _get_units(data)
    # Pandas object. DataFrame has no native name attribute but user can add one
    # See: https://github.com/pandas-dev/pandas/issues/447
    elif isinstance(data, (DataFrame, Series, Index)):
        title = getattr(data, 'name', None) or None
    # Pint Quantity
    elif isinstance(data, Quantity):
        if include_units:
            units = _get_units(data)
    # Add units or return units alone
    if title and units:
        title = f'{title} ({units})'
    else:
        title = title or units
    if title is not None:
        title = str(title).strip()
    return title


def _get_units(data):
    """
    Get the unit string from the `xarray.DataArray` attributes or the
    `pint.Quantity`. Format the latter with :rcraw:`unitformat`.
    """
    _load_objects()
    # Get units from the attributes
    if ndarray is not DataArray and isinstance(data, DataArray):
        units = data.attrs.get('units', None)
        data = data.data
        if units is not None:
            return units
    # Get units from the quantity
    if ndarray is not Quantity and isinstance(data, Quantity):
        fmt = rc.unitformat
        try:
            units = format(data.units, fmt)
        except (TypeError, ValueError):
            warnings._warn_proplot(
                f'Failed to format pint quantity with format string {fmt!r}.'
            )
        else:
            if 'L' in fmt:  # auto-apply LaTeX math indicator
                units = '$' + units + '$'
            return units


# Geographic utilties
def _geo_basemap_1d(x, *ys, xmin=-180, xmax=180):
    """
    Fix basemap geographic 1D data arrays.
    """
    # Ensure data is within map bounds
    x = _geo_monotonic(x)
    ys = _geo_clip(ys)
    x_orig, ys_orig, ys = x, ys, []
    for y_orig in ys_orig:
        x, y = _geo_inbounds(x_orig, y_orig, xmin=xmin, xmax=xmax)
        ys.append(y)
    return (x, *ys)


def _geo_basemap_2d(x, y, *zs, globe=False, xmin=-180, xmax=180):
    """
    Fix basemap geographic 2D data arrays.
    """
    x = _geo_monotonic(x)
    y = _geo_clip(y)  # in case e.g. edges() added points beyond poles
    x_orig, y_orig, zs_orig, zs = x, y, zs, []
    for z_orig in zs_orig:
        # Ensure data is within map bounds
        x, y, z = x_orig, y_orig, z_orig
        x, z = _geo_inbounds(x, z, xmin=xmin, xmax=xmax)
        if not globe or x.ndim > 1 or y.ndim > 1:
            zs.append(z)
            continue
        # Fix gaps in coverage
        y, z = _geo_poles(y, z)
        x, z = _geo_seams(x, z, xmin=xmin, modulo=False)
        zs.append(z)
    return (x, y, *zs)


def _geo_cartopy_1d(x, *ys):
    """
    Fix cartopy geographic 1D data arrays.
    """
    # So far not much to do here
    x = _geo_monotonic(x)
    ys = _geo_clip(ys)
    return (x, *ys)


def _geo_cartopy_2d(x, y, *zs, globe=False):
    """
    Fix cartopy geographic 2D data arrays.
    """
    x = _geo_monotonic(x)
    y = _geo_clip(y)  # in case e.g. edges() added points beyond poles
    x_orig, y_orig, zs_orig = x, y, zs
    zs = []
    for z_orig in zs_orig:
        # Bail for 2D coordinates
        x, y, z = x_orig, y_orig, z_orig
        if z is None or not globe or x.ndim > 1 or y.ndim > 1:
            zs.append(z)
            continue
        # Fix gaps in coverage
        y, z = _geo_poles(y, z)
        x, z = _geo_seams(x, z, modulo=True)
        zs.append(z)
    return (x, y, *zs)


def _geo_inbounds(x, y, xmin=-180, xmax=180):
    """
    Fix conflicts with map coordinates by rolling the data to fall between the
    minimum and maximum longitudes and masking out-of-bounds data points.
    """
    # Roll in same direction if some points on right-edge extend
    # more than 360 above min longitude; *they* should be on left side
    if x.ndim != 1:
        return x, y
    lonroll = np.where(x > xmin + 360)[0]  # tuple of ids
    if lonroll.size:  # non-empty
        roll = x.size - lonroll.min()
        x = np.roll(x, roll)
        y = np.roll(y, roll, axis=-1)
        x[:roll] -= 360  # make monotonic
    # Set NaN where data not in range xmin, xmax. Must be done for regional smaller
    # projections or get weird side-effects from valid data outside boundaries
    y = ma.masked_invalid(y, copy=False)
    y = y.astype(np.float).filled(np.nan)
    if x.size - 1 == y.shape[-1]:  # test western/eastern grid cell edges
        mask = (x[1:] < xmin) | (x[:-1] > xmax)
        y[..., mask] = np.nan
    elif x.size == y.shape[-1]:  # test the centers and pad by one for safety
        where, = np.where((x < xmin) | (x > xmax))
        y[..., where[1:-1]] = np.nan
    return x, y


def _geo_clip(*ys):
    """
    Ensure latitudes span only minus 90 to plus 90 degrees.
    """
    ys = tuple(np.clip(y, -90, 90) for y in ys)
    return ys[0] if len(ys) == 1 else ys


def _geo_monotonic(x):
    """
    Ensure longitudes are monotonic without rolling or filtering coordinates.
    """
    # Add 360 until data is greater than base value
    # TODO: Is this necessary for cartopy data? Maybe only with _geo_seams?
    if x.ndim != 1 or all(x < x[0]):  # skip 2D arrays and monotonic backwards data
        return x
    xmin = x[0]
    mask = np.array([True])
    while np.sum(mask):
        mask = x < xmin
        x[mask] += 360
    return x


def _geo_poles(y, z):
    """
    Fix gaps in coverage over the poles by adding data points at the poles
    using averages of the highest latitude data.
    """
    # Get means
    with np.errstate(all='ignore'):
        p1 = np.mean(z[0, :])  # do not ignore NaN if present
        p2 = np.mean(z[-1, :])
    if hasattr(p1, 'item'):
        p1 = np.asscalar(p1)  # happens with DataArrays
    if hasattr(p2, 'item'):
        p2 = np.asscalar(p2)
    # Concatenate
    ps = (-90, 90) if (y[0] < y[-1]) else (90, -90)
    z1 = np.repeat(p1, z.shape[1])[None, :]
    z2 = np.repeat(p2, z.shape[1])[None, :]
    y = ma.concatenate((ps[:1], y, ps[1:]))
    z = ma.concatenate((z1, z, z2), axis=0)
    return y, z


def _geo_seams(x, z, xmin=-180, modulo=False):
    """
    Fix gaps in coverage over longitude seams by adding points to either
    end or both ends of the array.
    """
    # Fix seams by ensuring circular coverage (cartopy can plot over map edges)
    if modulo:
        # Simply test the coverage % 360
        if x[0] % 360 != (x[-1] + 360) % 360:
            x = ma.concatenate((x, [x[0] + 360]))
            z = ma.concatenate((z, z[:, :1]), axis=1)
    else:
        # Ensure exact match between data and seams
        # Edges (e.g. pcolor) fit perfectly against seams. Size is unchanged.
        if x[0] == xmin and x.size - 1 == z.shape[1]:
            pass
        # Edges (e.g. pcolor) do not fit perfectly. Size augmented by 1.
        elif x.size - 1 == z.shape[1]:
            x = ma.append(xmin, x)
            x[-1] = xmin + 360
            z = ma.concatenate((z[:, -1:], z), axis=1)
        # Centers (e.g. contour) interpolated to edge. Size augmented by 2.
        elif x.size == z.shape[1]:
            xi = np.array([x[-1], x[0] + 360])
            if xi[0] == xi[1]:  # impossible to interpolate
                pass
            else:
                zq = ma.concatenate((z[:, -1:], z[:, :1]), axis=1)
                xq = xmin + 360
                zq = (zq[:, :1] * (xi[1] - xq) + zq[:, 1:] * (xq - xi[0])) / (xi[1] - xi[0])  # noqa: E501
                x = ma.concatenate(([xmin], x, [xmin + 360]))
                z = ma.concatenate((zq, z, zq), axis=1)
        else:
            raise ValueError('Unexpected shapes of coordinates or data arrays.')
    return x, z


# Misc utilties
def _get_vert(vert=None, orientation=None, **kwargs):
    """
    Get the orientation specified as either `vert` or `orientation`. This is
    used internally by various helper functions.
    """
    if vert is not None:
        return kwargs, vert
    elif orientation is not None:
        return kwargs, orientation != 'horizontal'  # should already be validated
    else:
        return kwargs, True  # fallback


def _parse_vert(
    vert=None, orientation=None, default_vert=None, default_orientation=None,
    **kwargs
):
    """
    Interpret both 'vert' and 'orientation' and add to outgoing keyword args
    if a default is provided.
    """
    # NOTE: Users should only pass these to hist, boxplot, or violinplot. To change
    # the plot, scatter, area, or bar orientation users should use the differently
    # named functions. Internally, however, they use these keyword args.
    if default_vert is not None:
        kwargs['vert'] = _not_none(
            vert=vert,
            orientation=None if orientation is None else orientation == 'vertical',
            default=default_vert,
        )
    if default_orientation is not None:
        kwargs['orientation'] = _not_none(
            orientation=orientation,
            vert=None if vert is None else 'vertical' if vert else 'horizontal',
            default=default_orientation,
        )
    if kwargs.get('orientation', None) not in (None, 'horizontal', 'vertical'):
        raise ValueError("Orientation must be either 'horizontal' or 'vertical'.")
    return kwargs


def _distribution_reduce(
    y, *, mean=None, means=None, median=None, medians=None, **kwargs
):
    """
    Return distribution columns reduced into means and medians. Tack on a
    distribution keyword argument for processing down the line.
    """
    # TODO: Permit 3D array with error dimension coming first
    data = y
    means = _not_none(mean=mean, means=means)
    medians = _not_none(median=median, medians=medians)
    if means and medians:
        warnings._warn_proplot('Cannot have both means=True and medians=True. Using former.')  # noqa: E501
        medians = None
    if means or medians:
        if data.ndim != 2:
            raise ValueError(f'Expected 2D array for means=True. Got {data.ndim}D.')
        data = ma.masked_invalid(data, copy=False)
        data = data.astype(np.float).filled(np.nan)
        if not data.size:
            raise ValueError('The input data contains all masked or NaN values.')
        elif means:
            y = np.nanmean(data, axis=0)
        elif medians:
            y = np.nanmedian(data, axis=0)
        kwargs['distribution'] = data

    # Save argument passed to _error_bars
    return (y, kwargs)


def _distribution_range(
    y, distribution, *, errdata=None, absolute=False, label=False,
    stds=None, pctiles=None, stds_default=None, pctiles_default=None,
):
    """
    Return a plottable characteristic range for the distribution keyword
    argument relative to the input coordinate means or medians.
    """
    # Parse stds arguments
    # NOTE: Have to guard against "truth value of an array is ambiguous" errors
    if stds is True:
        stds = stds_default
    elif stds is False or stds is None:
        stds = None
    else:
        stds = np.atleast_1d(stds)
        if stds.size == 1:
            stds = sorted((-stds.item(), stds.item()))
        elif stds.size != 2:
            raise ValueError('Expected scalar or length-2 stdev specification.')

    # Parse pctiles arguments
    if pctiles is True:
        pctiles = pctiles_default
    elif pctiles is False or pctiles is None:
        pctiles = None
    else:
        pctiles = np.atleast_1d(pctiles)
        if pctiles.size == 1:
            delta = (100 - pctiles.item()) / 2.0
            pctiles = sorted((delta, 100 - delta))
        elif pctiles.size != 2:
            raise ValueError('Expected scalar or length-2 pctiles specification.')

    # Incompatible settings
    if distribution is None and any(_ is not None for _ in (stds, pctiles)):
        raise ValueError(
            'To automatically compute standard deviations or percentiles on '
            'columns of data you must pass means=True or medians=True.'
        )
    if stds is not None and pctiles is not None:
        warnings._warn_proplot(
            'Got both a standard deviation range and a percentile range for '
            'auto error indicators. Using the standard deviation range.'
        )
        pctiles = None
    if distribution is not None and errdata is not None:
        stds = pctiles = None
        warnings._warn_proplot(
            'You explicitly provided the error bounds but also requested '
            'automatically calculating means or medians on data columns. '
            'It may make more sense to use the "stds" or "pctiles" keyword args '
            'and have *proplot* calculate the error bounds.'
        )

    # Compute error data in format that can be passed to maxes.Axes.errorbar()
    # NOTE: Include option to pass symmetric deviation from central points
    if errdata is not None:
        # Manual error data
        if y.ndim != 1:
            raise ValueError('errdata with 2D y coordinates is not yet supported.')
        label_default = 'uncertainty'
        err = _to_numpy_array(errdata)
        if (
            err.ndim not in (1, 2)
            or err.shape[-1] != y.size
            or err.ndim == 2 and err.shape[0] != 2
        ):
            raise ValueError(f'errdata has shape {err.shape}. Expected (2, {y.shape[-1]}).')  # noqa: E501
        if err.ndim == 1:
            abserr = err
            err = np.empty((2, err.size))
            err[0, :] = y - abserr  # translated back to absolute deviations below
            err[1, :] = y + abserr
    elif stds is not None:
        # Standard deviations
        # NOTE: Invalid values were handled by _distribution_reduce
        label_default = fr'{abs(stds[1])}$\sigma$ range'
        stds = _to_numpy_array(stds)[:, None]
        err = y + stds * np.nanstd(distribution, axis=0)
    elif pctiles is not None:
        # Percentiles
        # NOTE: Invalid values were handled by _distribution_reduce
        label_default = f'{pctiles[1] - pctiles[0]}% range'
        err = np.nanpercentile(distribution, pctiles, axis=0)
    else:
        raise ValueError('You must provide error bounds.')

    # Return data with legend entry
    if not absolute:  # for errorbar() ingestion
        err = err - y
        err[0, :] *= -1  # absolute deviations from central points
    if label is True:
        label = label_default
    elif not label:
        label = None
    return err, label


# Input preprocessor
def _preprocess_data(*keys, keywords=None, allow_extra=True):
    """
    Redirect internal plotting calls to native matplotlib methods. Also perform
    convert keyword args to positional and pass arguments through 'data' dictionary.
    """
    # Keyword arguments processed through 'data'
    # Positional arguments are always processed through data
    keywords = keywords or ()
    if isinstance(keywords, str):
        keywords = (keywords,)

    def decorator(func):
        name = func.__name__

        @functools.wraps(func)
        def _redirect_or_standardize(self, *args, **kwargs):
            if getattr(self, '_internal_call', None):
                # Redirect internal matplotlib call to native function
                func_native = getattr(super(PlotAxes, self), name)
                return func_native(*args, **kwargs)
            else:
                # Impose default coordinate system
                if self.name == 'proplot_basemap' and name in BASEMAP_FUNCS:
                    latlon = kwargs.pop('latlon', None)
                    kwargs['latlon'] = _not_none(latlon, True)
                if self.name == 'proplot_cartopy' and name in CARTOPY_FUNCS:
                    transform = kwargs.pop('transform', None)
                    kwargs['transform'] = _not_none(transform, PlateCarree())

                # Process data args
                # NOTE: Raises error if there are more args than keys
                args, kwargs = _keyword_to_positional(
                    keys, *args, allow_extra=allow_extra, **kwargs
                )
                data = kwargs.pop('data', None)
                if data is not None:
                    args = _get_data(data, *args)
                    for key in set(keywords) & set(kwargs):
                        kwargs[key] = _get_data(data, kwargs[key])

                # Auto-setup matplotlib with the input unit registry
                _load_objects()
                for arg in args:
                    if ndarray is not DataArray and isinstance(arg, DataArray):
                        arg = arg.data
                    if ndarray is not Quantity and isinstance(arg, Quantity):
                        ureg = getattr(arg, '_REGISTRY', None)
                        if hasattr(ureg, 'setup_matplotlib'):
                            ureg.setup_matplotlib(True)

                # Call main function
                return func(self, *args, **kwargs)  # call unbound method

        return _redirect_or_standardize
    return decorator


class PlotAxes(base.Axes):
    """
    The second lowest-level `~matplotlib.axes.Axes` subclass used by proplot.
    Implements all plotting overrides.
    """
    def __init__(self, *args, **kwargs):
        """
        Parameters
        ----------
        *args, **kwargs
            Passed to `~proplot.axes.Axes`.

        See also
        --------
        matplotlib.axes.Axes
        proplot.axes.Axes
        proplot.axes.CartesianAxes
        proplot.axes.PolarAxes
        proplot.axes.GeoAxes
        """
        super().__init__(*args, **kwargs)

    def _plot_safe(self, name, *args, **kwargs):
        """
        Call the plotting method and use context object to redirect internal
        calls to native methods. Finally add attributes to outgoing methods.
        """
        # NOTE: Previously allowed internal matplotlib plotting function calls to run
        # through proplot overrides then avoided awkward conflicts in piecemeal fashion.
        # Now prevent internal calls from running through overrides using preprocessor
        kwargs.pop('distribution', None)  # remove stat distributions
        with _state_context(self, _internal_call=True):
            if getattr(self, 'name', None) == 'proplot_basemap':
                obj = getattr(self.projection, name)(*args, ax=self, **kwargs)
            else:
                obj = getattr(super(), name)(*args, **kwargs)
        return obj

    def _plot_edges(self, method, *args, **kwargs):
        """
        Call the contour method to add "edges" to filled contours.
        """
        # NOTE: This is also used to provide an object that can be used by 'clabel'
        # for auto-labels. Filled contours seem to create weird artifacts.
        # NOTE: Make the default 'line width' identical to one used for pcolor plots
        # rather than rc['contour.linewidth']. See mpl pcolor() source code
        if not any(key in kwargs for key in ('linewidths', 'linestyles', 'edgecolors')):
            kwargs['linewidths'] = 0  # for clabel
        kwargs.setdefault('linewidths', EDGEWIDTH)
        kwargs.pop('cmap', None)
        kwargs['colors'] = kwargs.pop('edgecolors', 'k')
        return self._plot_safe(method, *args, **kwargs)

    def _plot_negpos(
        self, name, x, *ys, negcolor=None, poscolor=None, colorkey='facecolor',
        use_where=False, use_zero=False, **kwargs
    ):
        """
        Call the plot method separately for "negative" and "positive" data.
        """
        if use_where:
            kwargs.setdefault('interpolate', True)  # see fill_between docs
        for key in ('color', 'colors', 'facecolor', 'facecolors', 'where'):
            value = kwargs.pop(key, None)
            if value is not None:
                warnings._warn_proplot(
                    f'{name}() argument {key}={value!r} is incompatible with negpos=True. Ignoring.'  # noqa: E501
                )
        # Negative component
        yneg = list(ys)  # copy
        if use_zero:  # filter bar heights
            yneg[0] = _safe_mask(ys[0] < 0, ys[0])
        elif use_where:  # apply fill_between mask
            kwargs['where'] = ys[1] < ys[0]
        else:
            yneg = _safe_mask(ys[1] < ys[0], *ys)
        kwargs[colorkey] = _not_none(negcolor, rc['negcolor'])
        negobj = self._plot_safe(name, x, *yneg, **kwargs)
        # Positive component
        ypos = list(ys)  # copy
        if use_zero:  # filter bar heights
            ypos[0] = _safe_mask(ys[0] >= 0, ys[0])
        elif use_where:  # apply fill_between mask
            kwargs['where'] = ys[1] >= ys[0]
        else:
            ypos = _safe_mask(ys[1] >= ys[0], *ys)
        kwargs[colorkey] = _not_none(poscolor, rc['poscolor'])
        posobj = self._plot_safe(name, x, *ypos, **kwargs)
        return (negobj, posobj)

    def _add_queued_guide(
        self, objs, colorbar=None, colorbar_kw=None, legend=None, legend_kw=None,
    ):
        """
        Queue the input artist(s) for an automatic legend or colorbar once
        the axes is drawn or auto layout is called.
        """
        # WARNING: This should generally be last in the pipeline before calling
        # the plot function or looping over data columns. The colormap parser
        # and standardize functions both modify colorbar_kw and legend_kw.
        if colorbar:
            colorbar_kw = colorbar_kw or {}
            self.colorbar(objs, loc=colorbar, queue=True, **colorbar_kw)
        else:
            _guide_kw_to_obj(objs, 'colorbar', colorbar_kw)  # save for later
        if legend:
            legend_kw = legend_kw or {}
            self.legend(objs, loc=legend, queue=True, **legend_kw)
        else:
            _guide_kw_to_obj(objs, 'legend', legend_kw)  # save for later

    def _add_sticky_edges(self, objs, axis, *args):
        """
        Add sticky edges to the input artists using the minimum and maximum of the
        input coordinates. This is used to copy `bar` behavior to `area` and `lines`.
        """
        iter_ = list(obj for _ in objs for obj in (_ if isinstance(_, tuple) else (_,)))
        for sides in args:
            sides = np.atleast_1d(sides)
            if not sides.size:
                continue
            min_, max_ = _safe_range(sides)
            if min_ is None or max_ is None:
                continue
            for obj in iter_:
                convert = getattr(self, 'convert_' + axis + 'units')
                edges = getattr(obj.sticky_edges, axis)
                edges.extend(convert((min_, max_)))

    def _auto_format_1d(
        self, x, *ys, zerox=False, autox=True, autoy=True, autoformat=None,
        autoreverse=True, autolabels=True, autovalues=False, autoguide=True,
        label=None, labels=None, value=None, values=None, **kwargs
    ):
        """
        Try to retrieve default coordinates from array-like objects and apply default
        formatting. Also update the keyword arguments.
        """
        # Parse input
        y = max(ys, key=lambda y: y.size)  # try to find a non-scalar y for metadata
        autox = autox and not zerox  # so far just relevant for hist()
        autoformat = _not_none(autoformat, rc['autoformat'])
        kwargs, vert = _get_vert(**kwargs)
        labels = _not_none(
            label=label,
            labels=labels,
            value=value,
            values=values,
            legend_kw_labels=kwargs.get('legend_kw', {}).pop('labels', None),
            colorbar_kw_values=kwargs.get('colorbar_kw', {}).pop('values', None),
        )

        # Retrieve the x coords
        # NOTE: Where columns represent distributions, like for box and violinplot or
        # where we use 'means' or 'medians', columns coords (axis 1) are 'x' coords.
        # Otherwise, columns represent e.g. lines and row coords (axis 0) are 'x'
        # coords. Exception is passing "ragged arrays" to boxplot and violinplot.
        dists = any(kwargs.get(s) for s in ('mean', 'means', 'median', 'medians'))
        raggd = any(getattr(y, 'dtype', None) == 'object' for y in ys)
        xaxis = 0 if raggd else 1 if dists or not autoy else 0
        if autox and x is None:
            x = _get_labels(y, axis=xaxis)  # use the first one

        # Retrieve the labels. We only want default legend labels if this is an
        # object with 'title' metadata and/or the coords are string.
        # WARNING: Confusing terminology differences here -- for box and violin plots
        # labels refer to indices along x axis.
        if autolabels and labels is None:
            laxis = 0 if not autox and not autoy else xaxis if not autoy else xaxis + 1
            if laxis >= y.ndim:
                labels = _get_title(y)
            else:
                labels = _get_labels(y, axis=laxis, always=False)
            notitle = not _get_title(labels)
            if labels is None:
                pass
            elif notitle and not any(isinstance(_, str) for _ in labels):
                labels = None

        # Apply the labels or values
        if labels is not None:
            if autovalues:
                kwargs['values'] = _to_numpy_array(labels)
            elif autolabels:
                kwargs['labels'] = _to_numpy_array(labels)

        # Apply title for legend or colorbar that uses the labels or values
        if autoguide and autoformat:
            title = _get_title(labels)
            if title:  # safely update legend_kw and colorbar_kw
                _guide_kw_to_arg('legend', kwargs, title=title)
                _guide_kw_to_arg('colorbar', kwargs, label=title)

        # Apply the basic x and y settings
        autox = autox and self.name == 'proplot_cartesian'
        autoy = autoy and self.name == 'proplot_cartesian'
        sx, sy = 'xy' if vert else 'yx'
        kw_format = {}
        if autox and autoformat:  # 'x' axis
            title = _get_title(x)
            if title:
                axis = getattr(self, sx + 'axis')
                if axis.isDefault_label:
                    kw_format[sx + 'label'] = title
        if autoy and autoformat:  # 'y' axis
            sy = sx if zerox else sy  # hist() 'y' values are along 'x' axis
            title = _get_title(y)
            if title:
                axis = getattr(self, sy + 'axis')
                if axis.isDefault_label:
                    kw_format[sy + 'label'] = title

        # Convert string-type coordinates
        # NOTE: This should even allow qualitative string input to hist()
        if autox:
            x, kw_format = _get_coords(x, which=sx, **kw_format)
        if autoy:
            *ys, kw_format = _get_coords(*ys, which=sy, **kw_format)
        if autox and autoreverse and x.ndim == 1 and x.size > 1 and x[1] < x[0]:
            kw_format[sx + 'reverse'] = True

        # Apply formatting
        if kw_format:
            self.format(**kw_format)

        # Finally strip metadata
        # WARNING: Most methods that accept 2D arrays use columns of data, but when
        # pandas DataFrame specifically is passed to hist, boxplot, or violinplot, rows
        # of data assumed! Converting to ndarray necessary.
        ys = tuple(map(_to_numpy_array, ys))
        if x is not None:  # pie() and hist()
            x = _to_numpy_array(x)
        return (x, *ys, kwargs)

    def _standardize_1d(self, x, *ys, **kwargs):
        """
        Interpret positional arguments for all "1D" plotting commands.
        """
        # Standardize values
        zerox = not ys
        if zerox or all(y is None for y in ys):  # pad with remaining Nones
            x, *ys = None, x, *ys[1:]
        if len(ys) == 2:  # 'lines' or 'fill_between'
            if ys[1] is None:
                ys = (np.array([0.0]), ys[0])  # user input 1 or 2 positional args
            elif ys[0] is None:
                ys = (np.array([0.0]), ys[1])  # user input keyword 'y2' but no y1
        if any(y is None for y in ys):
            raise ValueError('Missing required data array argument.')
        ys = tuple(map(_to_duck_array, ys))
        if x is not None:
            x = _to_duck_array(x)
        x, *ys, kwargs = self._auto_format_1d(x, *ys, zerox=zerox, **kwargs)

        # Geographic corrections
        if self.name == 'proplot_cartopy' and isinstance(kwargs.get('transform'), PlateCarree):  # noqa: E501
            x, *ys = _geo_cartopy_1d(x, *ys)
        elif self.name == 'proplot_basemap' and kwargs.get('latlon', None):
            xmin, xmax = self.projection.lonmin, self.projection.lonmax
            x, *ys = _geo_basemap_1d(x, *ys, xmin=xmin, xmax=xmax)

        return (x, *ys, kwargs)

    def _auto_format_2d(self, x, y, *zs, autoformat=None, autoguide=True, **kwargs):
        """
        Try to retrieve default coordinates from array-like objects and apply default
        formatting. Also apply optional transpose and update the keyword arguments.
        """
        # Retrieve coordinates
        autoformat = _not_none(autoformat, rc['autoformat'])
        if x is None and y is None:
            z = zs[0]
            if z.ndim == 1:
                x = _get_labels(z, axis=0)
                y = np.zeros(z.shape)  # default barb() and quiver() behavior in mpl
            else:
                x = _get_labels(z, axis=1)
                y = _get_labels(z, axis=0)

        # Apply labels and XY axis settings
        if self.name == 'proplot_cartesian':
            # Apply labels
            # NOTE: Do not overwrite existing labels!
            kw_format = {}
            if autoformat:
                for s, d in zip('xy', (x, y)):
                    title = _get_title(d)
                    if title:
                        axis = getattr(self, s + 'axis')
                        if axis.isDefault_label:
                            kw_format[s + 'label'] = title

            # Handle string-type coordinates
            x, kw_format = _get_coords(x, which='x', **kw_format)
            y, kw_format = _get_coords(y, which='y', **kw_format)
            for s, d in zip('xy', (x, y)):
                if d.size > 1 and d.ndim == 1 and _to_numpy_array(d)[1] < _to_numpy_array(d)[0]:  # noqa: E501
                    kw_format[s + 'reverse'] = True

            # Apply formatting
            if kw_format:
                self.format(**kw_format)

        # Apply title for legend or colorbar
        if autoguide and autoformat:
            title = _get_title(zs[0])
            if title:  # safely update legend_kw and colorbar_kw
                _guide_kw_to_arg('legend', kwargs, title=title)
                _guide_kw_to_arg('colorbar', kwargs, label=title)

        # Finally strip metadata
        x = _to_numpy_array(x)
        y = _to_numpy_array(y)
        zs = tuple(map(_to_numpy_array, zs))
        return (x, y, *zs, kwargs)

    def _standardize_2d(
        self, x, y, *zs, globe=False, edges=False, allow1d=False, order='C',
        **kwargs
    ):
        """
        Interpret positional arguments for all "2D" plotting commands.
        """
        # Standardize values
        # NOTE: Functions pass two 'zs' at most right now
        if all(z is None for z in zs):
            x, y, zs = None, None, (x, y)[:len(zs)]
        if any(z is None for z in zs):
            raise ValueError('Missing required data array argument(s).')
        zs = tuple(_to_duck_array(z, strip_units=True) for z in zs)
        if x is not None:
            x = _to_duck_array(x)
        if y is not None:
            y = _to_duck_array(y)
        if order == 'F':
            zs = tuple(z.T for z in zs)
            if x is not None:
                x = x.T
            if y is not None:
                y = y.T
        x, y, *zs, kwargs = self._auto_format_2d(x, y, *zs, **kwargs)
        if edges:
            # NOTE: These functions quitely pass through 1D inputs, e.g. barb data
            x, y = _require_edges(x, y, zs[0])
        else:
            x, y = _require_centers(x, y, zs[0])

        # Geographic corrections
        if allow1d:
            pass
        elif self.name == 'proplot_cartopy' and isinstance(kwargs.get('transform'), PlateCarree):  # noqa: E501
            x, y, *zs = _geo_cartopy_2d(x, y, *zs, globe=globe)
        elif self.name == 'proplot_basemap' and kwargs.get('latlon', None):
            xmin, xmax = self.projection.lonmin, self.projection.lonmax
            x, y, *zs = _geo_basemap_2d(x, y, *zs, xmin=xmin, xmax=xmax, globe=globe)
            x, y = np.meshgrid(x, y)  # WARNING: required always

        return (x, y, *zs, kwargs)

    def _parse_inbounds(self, *, inbounds=None, **kwargs):
        """
        Capture the `inbounds` keyword arg and return data limit
        extents if it is ``True``. Otherwise return ``None``. When
        ``_restrict_inbounds`` gets ``None`` it will silently exit.
        """
        extents = None
        inbounds = _not_none(inbounds, rc['axes.inbounds'])
        if inbounds:
            extents = list(self.dataLim.extents)  # ensure modifiable
        return kwargs, extents

    def _mask_inbounds(self, x, y, z, *, to_centers=False):
        """
        Restrict the sample data used for automatic `vmin` and `vmax` selection
        based on the existing x and y axis limits.
        """
        # Get masks
        # WARNING: Experimental, seems robust but this is not mission-critical so
        # keep this in a try-except clause for now. However *internally* we should
        # not reach this block unless everything is an array so raise that error.
        xmask = ymask = None
        if self.name != 'proplot_cartesian':
            return z  # TODO: support geographic projections when input is PlateCarree()
        if not all(getattr(a, 'ndim', None) in (1, 2) for a in (x, y, z)):
            raise ValueError('Invalid input coordinates. Must be 1D or 2D arrays.')
        try:
            # Get centers and masks
            if to_centers and z.ndim == 2:
                x, y = _require_centers(x, y, z)
            if not self.get_autoscalex_on():
                xlim = self.get_xlim()
                xmask = (x >= min(xlim)) & (x <= max(xlim))
            if not self.get_autoscaley_on():
                ylim = self.get_ylim()
                ymask = (y >= min(ylim)) & (y <= max(ylim))
            # Get subsample
            if xmask is not None and ymask is not None:
                z = z[np.ix_(ymask, xmask)] if z.ndim == 2 and xmask.ndim == 1 else z[ymask & xmask]  # noqa: E501
            elif xmask is not None:
                z = z[:, xmask] if z.ndim == 2 and xmask.ndim == 1 else z[xmask]
            elif ymask is not None:
                z = z[ymask, :] if z.ndim == 2 and ymask.ndim == 1 else z[ymask]
            return z
        except Exception as err:
            warnings._warn_proplot(
                'Failed to restrict automatic colormap normalization algorithm '
                f'to in-bounds data only. Error message: {err}'
            )
            return z

    def _restrict_inbounds(self, extents, x, y, **kwargs):
        """
        Restrict the `dataLim` to exclude out-of-bounds data when x (y) limits
        are fixed and we are determining default y (x) limits. This modifies
        the mutable input `extents` to support iteration over columns.
        """
        # WARNING: This feature is still experimental. But seems obvious. Matplotlib
        # updates data limits in ad hoc fashion differently for each plotting command
        # but since proplot standardizes inputs we can easily use them for dataLim.
        kwargs, vert = _get_vert(**kwargs)
        if extents is None or self.name != 'proplot_cartesian':
            return
        if not x.size or not y.size:
            return
        if not vert:
            x, y = y, x
        trans = self.dataLim
        autox, autoy = self.get_autoscalex_on(), self.get_autoscaley_on()
        try:
            if autoy and not autox and x.shape == y.shape:
                # Reset the y data limits
                xmin, xmax = sorted(self.get_xlim())
                mask = (x >= xmin) & (x <= xmax)
                ymin, ymax = _safe_range(_safe_mask(mask, y))  # in-bounds y limits
                convert = self.convert_yunits  # handle datetime, pint units
                if ymin is not None:
                    trans.y0 = extents[1] = min(convert(ymin), extents[1])
                if ymax is not None:
                    trans.y1 = extents[3] = max(convert(ymax), extents[3])
                self._request_autoscale_view()
            if autox and not autoy and y.shape == x.shape:
                # Reset the x data limits
                ymin, ymax = sorted(self.get_ylim())
                mask = (y >= ymin) & (y <= ymax)
                xmin, xmax = _safe_range(_safe_mask(mask, x))
                convert = self.convert_xunits  # handle datetime, pint units
                if xmin is not None:
                    trans.x0 = extents[0] = min(convert(xmin), extents[0])
                if xmax is not None:
                    trans.x1 = extents[2] = max(convert(xmax), extents[2])
                    self._request_autoscale_view()
        except Exception as err:
            warnings._warn_proplot(
                'Failed to restrict automatic y (x) axis limit algorithm to '
                f'data within locked x (y) limits only. Error message: {err}'
            )

    def _parse_color(self, x, y, c, *, apply_cycle=True, **kwargs):
        """
        Parse either a colormap or color cycler. Colormap will be discrete and fade
        to subwhite luminance by default. Returns a HEX string if needed so we don't
        get ambiguous color warnings. Used with scatter, streamplot, quiver, barbs.
        """
        # NOTE: This function is positioned right above all _parse_cmap and
        # _parse_cycle functions and helper functions.
        if c is None or mcolors.is_color_like(c):
            if c is not None:
                c = pcolors.to_hex(c)  # avoid scatter() ambiguous color warning
            if apply_cycle:  # False for scatter() so we can wait to get correct 'N'
                kwargs = self._parse_cycle(**kwargs)
            methods = (self._parse_cmap, self._parse_levels, self._parse_autolev, self._parse_vlim)  # noqa: E501
        else:
            kwargs['line_plot'] = True
            kwargs['default_discrete'] = False
            kwargs = self._parse_cmap(x, y, c, **kwargs)
            methods = (self._parse_cycle,)
        pop = _pop_params(kwargs, *methods, ignore_internal=True)
        if pop:
            warnings._warn_proplot(f'Ignoring bad unused keyword arg(s): {pop}')
        return (c, kwargs)

    def _parse_vlim(
        self, *args,
        vmin=None, vmax=None, to_centers=False,
        robust=None, inbounds=None, **kwargs,
    ):
        """
        Return a suitable vmin and vmax based on the input data.

        Parameters
        ----------
        *args
            The sample data.
        vmin, vmax : float, optional
            The user input minimum and maximum.
        robust : bool, optional
            Whether to limit the default range to exclude outliers.
        inbounds : bool, optional
            Whether to filter to in-bounds data.
        to_centers : bool, optional
            Whether to convert coordinates to 'centers'.

        Returns
        -------
        vmin, vmax : float
            The minimum and maximum.
        kwargs
            Unused arguemnts.
        """
        # Parse vmin and vmax
        automin = vmin is None
        automax = vmax is None
        if not automin and not automax:
            return vmin, vmax, kwargs

        # Parse input args
        inbounds = _not_none(inbounds, rc['cmap.inbounds'])
        robust = _not_none(robust, rc['cmap.robust'], False)
        robust = 96 if robust is True else 100 if robust is False else robust
        robust = np.atleast_1d(robust)
        if robust.size == 1:
            pmin, pmax = 50 + 0.5 * np.array([-robust.item(), robust.item()])
        elif robust.size == 2:
            pmin, pmax = robust.flat  # pull out of array
        else:
            raise ValueError(f'Unexpected robust={robust!r}. Must be bool, float, or 2-tuple.')  # noqa: E501

        # Get sample data
        # NOTE: Critical to use _to_duck_array here because some commands
        # are unstandardized.
        # NOTE: Try to get reasonable *count* levels for hexbin/hist2d, but in general
        # have no way to select nice ones a priori (why we disable discretenorm).
        # NOTE: Currently we only ever use this function with *single* array input
        # but in future could make this public as a way for users (me) to get
        # automatic synced contours for a bunch of arrays in a grid.
        vmins, vmaxs = [], []
        if len(args) > 2:
            x, y, *zs = args
        else:
            x, y, *zs = None, None, *args
        for z in zs:
            if z is None:  # e.g. empty scatter color
                continue
            if z.ndim > 2:  # e.g. imshow data
                continue
            z = _to_numpy_array(z)
            if inbounds and x is not None and y is not None:  # ignore if None coords
                z = self._mask_inbounds(x, y, z, to_centers=to_centers)
            vmin, vmax = _safe_range(z, pmin, pmax, automin=automin, automax=automax)
            if vmin is not None:
                vmins.append(vmin)
            if vmax is not None:
                vmaxs.append(vmax)
        return min(vmins, default=0), max(vmaxs, default=1), kwargs

    def _parse_autolev(
        self, *args, levels=None,
        extend='neither', norm=None, norm_kw=None, vmin=None, vmax=None,
        locator=None, locator_kw=None, symmetric=None, **kwargs
    ):
        """
        Return a suitable level list given the input data, normalizer,
        locator, and vmin and vmax.

        Parameters
        ----------
        *args
            The sample data. Passed to `_parse_vlim`.
        levels : int
            The approximate number of levels.
        vmin, vmax : float, optional
            The approximate minimum and maximum level edges. Passed to the locator.
        diverging : bool, optional
            Whether the resulting levels are intended for a diverging normalizer.
        symmetric : bool, optional
            Whether the resulting levels should be symmetric about zero.
        norm, norm_kw : optional
            Passed to `~proplot.constructor.Norm`. Used to change the default
            `locator` (e.g., a `~matplotlib.colors.LogNorm` normalizer will use
            a `~matplotlib.ticker.LogLocator` to generate levels).

        Parameters
        ----------
        levels : list of float
            The level edges.
        kwargs
            Unused arguments.
        """
        # Input args
        # NOTE: Some of this is adapted from the hidden contour.ContourSet._autolev
        # NOTE: We use 'symmetric' with MaxNLocator to ensure boundaries include a
        # zero level but may trim many of these below.
        norm_kw = norm_kw or {}
        locator_kw = locator_kw or {}
        levels = _not_none(levels, rc['cmap.levels'])
        vmin = _not_none(vmin=vmin, norm_kw_vmin=norm_kw.pop('vmin', None))
        vmax = _not_none(vmax=vmax, norm_kw_vmax=norm_kw.pop('vmax', None))
        norm = constructor.Norm(norm or 'linear', **norm_kw)
        symmetric = _not_none(
            symmetric=symmetric,
            locator_kw_symmetric=locator_kw.pop('symmetric', None),
            default=False,
        )

        # Get default locator from input norm
        # NOTE: This normalizer is only temporary for inferring level locs
        norm = constructor.Norm(norm or 'linear', **norm_kw)
        if locator is not None:
            locator = constructor.Locator(locator, **locator_kw)
        elif isinstance(norm, mcolors.LogNorm):
            locator = mticker.LogLocator(**locator_kw)
        elif isinstance(norm, mcolors.SymLogNorm):
            for key, default in (('base', 10), ('linthresh', 1)):
                val = _not_none(getattr(norm, key, None), getattr(norm, '_' + key, None), default)  # noqa: E501
                locator_kw.setdefault(key, val)
            locator = mticker.SymmetricalLogLocator(**locator_kw)
        else:
            locator_kw['symmetric'] = symmetric
            locator = mticker.MaxNLocator(levels, min_n_ticks=1, **locator_kw)

        # Get default level locations
        nlevs = levels
        automin = vmin is None
        automax = vmax is None
        vmin, vmax, kwargs = self._parse_vlim(*args, vmin=vmin, vmax=vmax, **kwargs)
        try:
            levels = locator.tick_values(vmin, vmax)
        except RuntimeError:  # too-many-ticks error
            levels = np.linspace(vmin, vmax, levels)  # TODO: _autolev used N+1

        # Possibly trim levels far outside of 'vmin' and 'vmax'
        # NOTE: This part is mostly copied from matplotlib _autolev
        if not symmetric:
            i0, i1 = 0, len(levels)  # defaults
            under, = np.where(levels < vmin)
            if len(under):
                i0 = under[-1]
                if not automin or extend in ('min', 'both'):
                    i0 += 1  # permit out-of-bounds data
            over, = np.where(levels > vmax)
            if len(over):
                i1 = over[0] + 1 if len(over) else len(levels)
                if not automax or extend in ('max', 'both'):
                    i1 -= 1  # permit out-of-bounds data
            if i1 - i0 < 3:
                i0, i1 = 0, len(levels)  # revert
            levels = levels[i0:i1]

        # Compare the no. of levels we *got* (levels) to what we *wanted* (N)
        # If we wanted more than 2 times the result, then add nn - 1 extra
        # levels in-between the returned levels *in normalized space*.
        # Example: A LogNorm gives too few levels, so we select extra levels
        # here, but use the locator for determining tick locations.
        nn = nlevs // len(levels)
        if nn >= 2:
            olevels = norm(levels)
            nlevels = []
            for i in range(len(levels) - 1):
                l1, l2 = olevels[i], olevels[i + 1]
                nlevels.extend(np.linspace(l1, l2, nn + 1)[:-1])
            nlevels.append(olevels[-1])
            levels = norm.inverse(nlevels)

        return levels, kwargs

    def _parse_levels(
        self, *args, N=None, levels=None, values=None, minlength=2,
        positive=False, negative=False, nozero=False, norm=None, norm_kw=None,
        vmin=None, vmax=None, skip_autolev=False, **kwargs,
    ):
        """
        Return levels resulting from a wide variety of keyword options.

        Parameters
        ----------
        *args
            The sample data. Passed to `_parse_vlim`.
        N
            Shorthand for `levels`.
        levels : int or list of float, optional
            The levels list or (approximate) number of levels to create.
        values : int or list of float, optional
            The level center list or (approximate) number of level centers to create.
        minlength : int, optional
            The minimum number of levels allowed.
        positive, negative, nozero : bool, optional
            Whether to remove out non-positive, non-negative, and zero-valued
            levels. The latter is useful for single-color contour plots.
        norm, norm_kw : optional
            Passed to `~proplot.constructor.Norm`. Used to possbily infer levels
            or to convert values to levels.
        vmin, vmax
            Passed to ``_parse_autolev``.

        Returns
        -------
        levels : list of float
            The level edges.
        kwargs
            Unused arguments.
        """
        # Rigorously check user input levels and values
        # NOTE: Include special case where color levels are referenced
        # by string label values.
        levels = _not_none(
            N=N, levels=levels, norm_kw_levels=norm_kw.pop('levels', None),
        )
        if positive and negative:
            negative = False
            warnings._warn_proplot(
                'Incompatible args positive=True and negative=True. Using former.'
            )
        if levels is not None and values is not None:
            warnings._warn_proplot(
                f'Incompatible args levels={levels!r} and values={values!r}. Using former.'  # noqa: E501
            )
        for key, val in (('levels', levels), ('values', values)):
            if val is None:
                continue
            if isinstance(norm, (mcolors.BoundaryNorm, pcolors.SegmentedNorm)):
                warnings._warn_proplot(
                    f'Ignoring {key}={val}. Instead using norm={norm!r} boundaries.'
                )
            if not np.iterable(val):
                continue
            if len(val) < minlength:
                raise ValueError(
                    f'Invalid {key}={val}. Must be at least length {minlength}.'
                )
            if len(val) >= 2 and np.any(np.sign(np.diff(val)) != np.sign(val[1] - val[0])):  # noqa: E501
                raise ValueError(
                    f'Invalid {key}={val}. Must be monotonically increasing or decreasing.'  # noqa: E501
                )
        if isinstance(norm, (mcolors.BoundaryNorm, pcolors.SegmentedNorm)):
            levels, values = norm.boundaries, None
        else:
            levels = _not_none(levels, rc['cmap.levels'])

        # Infer level edges from level centers if possible
        # NOTE: The only way for user to manually impose BoundaryNorm is by
        # passing one -- users cannot create one using Norm constructor key.
        descending = None
        if values is None:
            pass
        elif isinstance(values, Integral):
            levels = values + 1
        elif np.iterable(values) and len(values) == 1:
            levels = [values[0] - 1, values[0] + 1]  # weird but why not
        elif norm is None or norm in ('segments', 'segmented'):
            # Try to generate levels so SegmentedNorm will place 'values' ticks at the
            # center of each segment. edges() gives wrong result unless spacing is even.
            # Solve: (x1 + x2) / 2 = y --> x2 = 2 * y - x1 with arbitrary starting x1.
            values, descending = pcolors._sanitize_levels(values)
            levels = [values[0] - (values[1] - values[0]) / 2]  # arbitrary x1
            for val in values:
                levels.append(2 * val - levels[-1])
            if any(np.diff(levels) < 0):  # backup plan in event of weird ticks
                levels = edges(values)
            if descending:  # then revert back below
                levels = levels[::-1]
        else:
            # Generate levels by finding in-between points in the
            # normalized numeric space, e.g. LogNorm space.
            norm_kw = norm_kw or {}
            convert = constructor.Norm(norm, **norm_kw)
            levels = convert.inverse(edges(convert(values)))

        # Process level edges and infer defaults
        # NOTE: Matplotlib colorbar algorithm *cannot* handle descending levels so
        # this function reverses them and adds special attribute to the normalizer.
        # Then colorbar() reads this attr and flips the axis and the colormap direction
        if np.iterable(levels) and len(levels) > 2:
            levels, descending = pcolors._sanitize_levels(levels)
        if not np.iterable(levels) and not skip_autolev:
            levels, kwargs = self._parse_autolev(
                *args, levels=levels, vmin=vmin, vmax=vmax, norm=norm, norm_kw=norm_kw, **kwargs  # noqa: E501
            )
        ticks = values if np.iterable(values) else levels
        if descending is not None:
            kwargs.setdefault('descending', descending)  # for _parse_discrete
        if ticks is not None and np.iterable(ticks):
            _guide_kw_to_arg('colorbar', kwargs, locator=ticks)

        # Filter the resulting level boundaries
        if levels is not None and np.iterable(levels):
            if nozero and 0 in levels:
                levels = levels[levels != 0]
            if positive:
                levels = levels[levels >= 0]
            if negative:
                levels = levels[levels <= 0]
        return levels, kwargs

    def _parse_discrete(
        self, levels, norm, cmap, extend='neither', descending=False, **kwargs,
    ):
        """
        Create a `~proplot.colors.DiscreteNorm` or `~proplot.colors.BoundaryNorm`
        from the input colormap and normalizer.

        Parameters
        ----------
        levels : list of float
            The level boundaries.
        norm : `~matplotlib.colors.Normalize`
            The continuous normalizer.
        cmap : `~matplotlib.colors.Colormap`
            The colormap.
        extend : str, optional
            The extend setting.
        descending : bool, optional
            Whether levels are descending.

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
        # See _parse_levels for details
        under = cmap._rgba_under
        over = cmap._rgba_over
        unique = extend  # default behavior
        cyclic = getattr(cmap, '_cyclic', None)
        qualitative = isinstance(cmap, pcolors.DiscreteColormap)  # see _parse_cmap
        if descending:
            cmap = cmap.reversed()

        # Ensure end colors are unique by scaling colors as if extend='both'
        # NOTE: Inside _parse_cmap should have enforced extend='neither'
        if cyclic:
            step = 0.5
            unique = 'both'

        # Ensure color list matches level list
        # NOTE: If user-input levels were integer or omitted then integer levels
        # passed to level will have matched
        elif qualitative:
            # Truncate or wrap color list (see matplotlib.ListedColormap)
            step = 0.5  # try to sample the central color index
            auto_under = under is None and extend in ('min', 'both')
            auto_over = over is None and extend in ('max', 'both')
            ncolors = len(levels) - 1 + auto_under + auto_over
            colors = list(itertools.islice(itertools.cycle(cmap.colors), ncolors))
            # Create new colormap and optionally apply colors to extremes
            if auto_under:
                under, *colors = colors
            if auto_over:
                *colors, over = colors
            cmap = cmap.copy(colors, N=len(colors))
            if under is not None:
                cmap.set_under(under)
            if over is not None:
                cmap.set_over(over)

        # Ensure middle colors sample full range when extreme colors are present
        # by scaling colors as if extend='neither'
        else:
            # Keep unique bins
            step = 1.0
            if over is not None and under is not None:
                unique = 'neither'
            # Turn off unique bin for over-bounds colors
            elif over is not None:
                if extend == 'both':
                    unique = 'min'
                elif extend == 'max':
                    unique = 'neither'
            # Turn off unique bin for under-bounds colors
            elif under is not None:
                if extend == 'both':
                    unique = 'min'
                elif extend == 'max':
                    unique = 'neither'

        # Generate DiscreteNorm and update "child" norm with vmin and vmax from
        # levels. This lets the colorbar set tick locations properly!
        if not isinstance(norm, mcolors.BoundaryNorm) and len(levels) > 1:
            norm = pcolors.DiscreteNorm(
                levels, norm=norm, descending=descending, unique=unique, step=step,
            )

        return norm, cmap, kwargs

    @warnings._rename_kwargs('0.6', centers='values')
    @_snippet_manager
    def _parse_cmap(
        self, *args,
        cmap=None, cmap_kw=None, c=None, color=None, colors=None, default_cmap=None,
        norm=None, norm_kw=None, extend='neither', vmin=None, vmax=None,
        sequential=None, diverging=None, qualitative=None, cyclic=None,
        discrete=None, default_discrete=True, skip_autolev=False,
        line_plot=False, contour_plot=False, **kwargs
    ):
        """
        Parse colormap and normalizer arguments.
        """
        # Parse keyword args
        # NOTE: Always disable 'autodiverging' when an unknown colormap is passed to
        # avoid awkwardly combining 'DivergingNorm' with sequential colormaps.
        # However let people use diverging=False with diverging cmaps because
        # some use them (wrongly IMO but nbd) for increased color contrast.
        cmap_kw = cmap_kw or {}
        norm_kw = norm_kw or {}
        vmin = _not_none(vmin=vmin, norm_kw_vmin=norm_kw.pop('vmin', None))
        vmax = _not_none(vmax=vmax, norm_kw_vmax=norm_kw.pop('vmax', None))
        colors = _not_none(c=c, color=color, colors=colors)  # in case untranslated
        autodiverging = rc['cmap.autodiverging']
        name = getattr(cmap, 'name', cmap)
        if isinstance(name, str) and name not in pcolors.CMAPS_DIVERGING:
            autodiverging = False  # avoid auto-truncation of sequential colormaps

        # Build qualitative colormap using 'colors'
        # NOTE: Try to match number of level centers with number of colors here
        # WARNING: Previously 'colors' set the edgecolors. To avoid all-black
        # colormap make sure to ignore 'colors' if 'cmap' was also passed.
        # WARNING: Previously tried setting number of levels to len(colors) but
        # this would make single-level contour plots and _parse_autolev is designed
        # to only give approximate level count so failed anyway. Users should pass
        # their own levels to avoid truncation/cycling in these very special cases.
        if cmap is not None and colors is not None:
            warnings._warn_proplot(
                f'You specifed both cmap={cmap!r} and the qualitative-colormap '
                f'colors={colors!r}. Ignoring the latter.'
            )
            colors = None
        if colors is not None:
            if mcolors.is_color_like(colors):
                colors = [colors]  # RGB[A] tuple possibly
            cmap = colors = np.atleast_1d(colors)
            cmap_kw['listmode'] = 'discrete'

        # Create the user-input colormap
        # Also force options in special cases
        if line_plot:
            cmap_kw['default_luminance'] = pcolors.CYCLE_LUMINANCE
        if cmap is not None:
            cmap = constructor.Colormap(cmap, **cmap_kw)  # for testing only
        cyclic = _not_none(cyclic, getattr(cmap, '_cyclic', None))
        if cyclic and extend != 'neither':
            warnings._warn_proplot(f"Cyclic colormaps require extend='neither'. Ignoring extend={extend!r}.")  # noqa: E501
            extend = 'neither'
        qualitative = _not_none(qualitative, isinstance(cmap, pcolors.DiscreteColormap))
        if qualitative and discrete is not None and not discrete:
            warnings._warn_proplot('Qualitative colormaps require discrete=True. Ignoring discrete=False.')  # noqa: E501
            discrete = True
        if contour_plot and discrete is not None and not discrete:
            warnings._warn_proplot('Contoured plots require discrete=True. Ignoring discrete=False.')  # noqa: E501
            discrete = True
        keys = ('levels', 'values', 'locator', 'negative', 'positive', 'symmetric')
        if discrete is None and any(key in kwargs for key in keys):
            discrete = True  # override
        else:
            discrete = _not_none(discrete, rc['cmap.discrete'], default_discrete)

        # Determine the appropriate 'vmin', 'vmax', and/or 'levels'
        # NOTE: Unlike xarray, but like matplotlib, vmin and vmax only approximately
        # determine level range. Levels are selected with Locator.tick_values().
        levels = None  # unused
        if discrete:
            levels, kwargs = self._parse_levels(
                *args, vmin=vmin, vmax=vmax, norm=norm, norm_kw=norm_kw,
                extend=extend, skip_autolev=skip_autolev, **kwargs
            )
        elif not skip_autolev:
            vmin, vmax, kwargs = self._parse_vlim(
                *args, vmin=vmin, vmax=vmax, **kwargs
            )
        if autodiverging:
            default_diverging = None
            if levels is not None:
                _, counts = np.unique(np.sign(levels), return_counts=True)
                if counts[counts > 1].size > 1:
                    default_diverging = True
            elif vmin is not None and vmax is not None:
                if np.sign(vmin) != np.sign(vmax):
                    default_diverging = True
            diverging = _not_none(diverging, default_diverging)

        # Create the continuous normalizer. Only use SegmentedNorm if necessary
        # NOTE: We create normalizer here only because auto level generation depends
        # on the normalizer class (e.g. LogNorm). We don't have to worry about vmin
        # and vmax because they get applied to normalizer inside DiscreteNorm.
        if norm is None and levels is not None and len(levels) > 0:
            if len(levels) == 1:  # edge case, use central colormap color
                vmin = _not_none(vmin, levels[0] - 1)
                vmax = _not_none(vmax, levels[0] + 1)
            else:
                vmin, vmax = min(levels), max(levels)
                steps = np.abs(np.diff(levels))
                eps = np.mean(steps) / 1e3
                if np.any(np.abs(np.diff(steps)) >= eps):
                    norm = 'segmented'
        if norm in ('segments', 'segmented'):
            if np.iterable(levels):
                norm_kw['levels'] = levels  # apply levels
            else:
                norm = None  # same result but much faster
        if diverging:
            norm = _not_none(norm, 'div')
        else:
            norm = _not_none(norm, 'linear')
        norm = constructor.Norm(norm, vmin=vmin, vmax=vmax, **norm_kw)
        if autodiverging and isinstance(norm, pcolors.DivergingNorm):
            diverging = _not_none(diverging, True)

        # Create the final colormap
        if cmap is None:
            if default_cmap is not None:  # used internally
                cmap = default_cmap
            elif qualitative:
                cmap = rc['cmap.qualitative']
            elif cyclic:
                cmap = rc['cmap.cyclic']
            elif diverging:
                cmap = rc['cmap.diverging']
            elif sequential:
                cmap = rc['cmap.sequential']
            cmap = _not_none(cmap, rc['image.cmap'])
            cmap = constructor.Colormap(cmap, **cmap_kw)

        # Create the discrete normalizer
        # Then finally warn and remove unused args
        if levels is not None:
            kwargs['extend'] = extend
            norm, cmap, kwargs = self._parse_discrete(levels, norm, cmap, **kwargs)
        methods = (self._parse_levels, self._parse_autolev, self._parse_vlim)
        params = _pop_params(kwargs, *methods, ignore_internal=True)
        if 'N' in params:  # use this for lookup table N instead of levels N
            cmap = cmap.copy(N=params.pop('N'))
        if params:
            warnings._warn_proplot(f'Ignoring unused keyword args(s): {params}')

        # Update outgoing args
        # NOTE: With contour(..., discrete=False, levels=levels) users can bypass
        # proplot's level selection and use native matplotlib level selection
        if contour_plot:
            kwargs['levels'] = levels
            kwargs['extend'] = extend
        kwargs.update({'cmap': cmap, 'norm': norm})
        _guide_kw_to_arg('colorbar', kwargs, extend=extend)

        return kwargs

    def _iter_pairs(self, *args):
        """
        Iterate over ``[x1,] y1, [fmt1,] [x2,] y2, [fmt2,] ...`` input.
        """
        # NOTE: This is copied from _process_plot_var_args.__call__ to avoid relying
        # on private API. We emulate this input style with successive plot() calls.
        args = list(args)
        while args:  # this permits empty input
            x, y, *args = args
            if args and isinstance(args[0], str):  # format string detected!
                fmt, *args = args
            elif isinstance(y, str):  # omits some of matplotlib's rigor but whatevs
                x, y, fmt = None, x, y
            else:
                fmt = None
            yield x, y, fmt

    def _iter_columns(self, *args, label=None, labels=None, values=None, **kwargs):
        """
        Iterate over columns of positional arguments and add successive ``'label'``
        keyword arguments using the input label-list ``'labels'``.
        """
        # Handle cycle args and label lists
        # NOTE: Arrays here should have had metadata stripped by _standardize_1d
        # but could still be pint quantities that get processed by axis converter.
        n = max(1 if not _is_array(a) or a.ndim < 2 else a.shape[-1] for a in args)
        labels = _not_none(label=label, values=values, labels=labels)
        if not np.iterable(labels) or isinstance(labels, str):
            labels = n * [labels]
        if len(labels) != n:
            raise ValueError(f'Array has {n} columns but got {len(labels)} labels.')
        if labels is not None:
            labels = [str(_not_none(label, '')) for label in _to_numpy_array(labels)]
        else:
            labels = n * [None]

        # Yield successive columns
        for i in range(n):
            kw = kwargs.copy()
            kw['label'] = labels[i] or None
            a = tuple(a if not _is_array(a) or a.ndim < 2 else a[..., i] for a in args)
            yield (i, n, *a, kw)

    def _parse_cycle(
        self, ncycle=None, *,
        cycle=None, cycle_kw=None, cycle_manually=None, return_cycle=False, **kwargs
    ):
        """
        Parse property cycle-related arguments.
        """
        # Create the property cycler and update it if necessary
        # NOTE: Matplotlib Cycler() objects have built-in __eq__ operator
        # so really easy to check if the cycler has changed!
        if cycle is not None or cycle_kw:
            cycle_kw = cycle_kw or {}
            if ncycle != 1:  # ignore for column-by-column plotting commands
                cycle_kw.setdefault('N', ncycle)  # if None then filled in Colormap()
            if isinstance(cycle, str) and cycle.lower() == 'none':
                cycle = False
            if not cycle:
                args = ()
            elif cycle is True:  # consistency with 'False' ('reactivate' the cycler)
                args = (rc['axes.prop_cycle'],)
            else:
                args = (cycle,)
            cycle = constructor.Cycle(*args, **cycle_kw)
            with warnings.catch_warnings():  # hide 'elementwise-comparison failed'
                warnings.simplefilter('ignore', FutureWarning)
                if return_cycle:
                    pass
                elif cycle != self._active_cycle:
                    self.set_prop_cycle(cycle)

        # Manually extract and apply settings to outgoing keyword arguments
        # if native matplotlib function does not include desired properties
        cycle_manually = cycle_manually or {}
        parser = self._get_lines  # the _process_plot_var_args instance
        props = {}  # which keys to apply from property cycler
        for prop, key in cycle_manually.items():
            value = kwargs.get(key, None)
            if value is None and prop in parser._prop_keys:
                props[prop] = key
        if props:
            dict_ = next(parser.prop_cycler)
            for prop, key in props.items():
                value = dict_[prop]
                if key == 'c':  # special case: scatter() color must be converted to hex
                    value = pcolors.to_hex(value)
                kwargs[key] = value

        if return_cycle:
            return cycle, kwargs  # needed for stem() to apply in a context()
        else:
            return kwargs

    def _error_bars(
        self, x, y, *_, distribution=None,
        default_bars=True, default_boxes=False,
        barstd=None, barstds=None, barpctile=None, barpctiles=None, bardata=None,
        boxstd=None, boxstds=None, boxpctile=None, boxpctiles=None, boxdata=None,
        capsize=None, **kwargs,
    ):
        """
        Add up to 2 error indicators: thick "boxes" and thin "bars".
        """
        # Parse input args
        # NOTE: Want to keep _error_bars() and _error_shading() separate. But also
        # want default behavior where some default error indicator is shown if user
        # requests means/medians only. Result is the below kludge.
        kwargs, vert = _get_vert(**kwargs)
        barstds = _not_none(barstd=barstd, barstds=barstds)
        boxstds = _not_none(boxstd=boxstd, boxstds=boxstds)
        barpctiles = _not_none(barpctile=barpctile, barpctiles=barpctiles)
        boxpctiles = _not_none(boxpctile=boxpctile, boxpctiles=boxpctiles)
        bars = any(_ is not None for _ in (bardata, barstds, barpctiles))
        boxes = any(_ is not None for _ in (boxdata, boxstds, boxpctiles))
        shade = any(
            prefix + suffix in key for key in kwargs
            for prefix in ('shade', 'fade') for suffix in ('std', 'pctile', 'data')
        )
        if distribution is not None and not bars and not boxes and not shade:
            barstds = bars = default_bars
            boxstds = boxes = default_boxes

        # Error bar properties
        edgecolor = kwargs.get('edgecolor', rc['boxplot.whiskerprops.color'])
        barprops = _pop_props(kwargs, 'line', ignore='marker', prefix='bar')
        barprops['capsize'] = _not_none(capsize, rc['errorbar.capsize'])
        barprops['linestyle'] = 'none'
        barprops.setdefault('color', edgecolor)
        barprops.setdefault('zorder', 2.5)
        barprops.setdefault('linewidth', rc['boxplot.whiskerprops.linewidth'])

        # Error box properties
        # NOTE: Includes 'markerfacecolor' and 'markeredgecolor' props
        boxprops = _pop_props(kwargs, 'line', prefix='box')
        boxprops['capsize'] = 0
        boxprops['linestyle'] = 'none'
        boxprops.setdefault('color', barprops['color'])
        boxprops.setdefault('zorder', barprops['zorder'])
        boxprops.setdefault('linewidth', 4 * barprops['linewidth'])

        # Box marker properties
        boxmarker = {key: boxprops.pop(key) for key in tuple(boxprops) if 'marker' in key}  # noqa: E501
        boxmarker['c'] = _not_none(boxmarker.pop('markerfacecolor', None), 'white')
        boxmarker['s'] = _not_none(boxmarker.pop('markersize', None), boxprops['linewidth'] ** 0.5)  # noqa: E501
        boxmarker['zorder'] = boxprops['zorder']
        boxmarker['edgecolor'] = boxmarker.pop('markeredgecolor', None)
        boxmarker['linewidth'] = boxmarker.pop('markerlinewidth', None)
        if boxmarker.get('marker') is True:
            boxmarker['marker'] = 'o'
        elif default_boxes:  # enable by default
            boxmarker.setdefault('marker', 'o')

        # Draw thin or thick error bars from distributions or explicit errdata
        sy = 'y' if vert else 'x'  # yerr
        ex, ey = (x, y) if vert else (y, x)
        eobjs = []
        if bars:  # now impossible to make thin bar width different from cap width!
            edata, _ = _distribution_range(
                y, distribution,
                stds=barstds, pctiles=barpctiles, errdata=bardata,
                stds_default=(-3, 3), pctiles_default=(0, 100),
            )
            obj = self.errorbar(ex, ey, **barprops, **{sy + 'err': edata})
            eobjs.append(obj)
        if boxes:  # must go after so scatter point can go on top
            edata, _ = _distribution_range(
                y, distribution,
                stds=boxstds, pctiles=boxpctiles, errdata=boxdata,
                stds_default=(-1, 1), pctiles_default=(25, 75),
            )
            obj = self.errorbar(ex, ey, **boxprops, **{sy + 'err': edata})
            if boxmarker.get('marker', None):
                self.scatter(ex, ey, **boxmarker)
            eobjs.append(obj)

        kwargs['distribution'] = distribution
        return (*eobjs, kwargs)

    def _error_shading(
        self, x, y, *_, distribution=None, color_key='color',
        shadestd=None, shadestds=None, shadepctile=None, shadepctiles=None, shadedata=None,  # noqa: E501
        fadestd=None, fadestds=None, fadepctile=None, fadepctiles=None, fadedata=None,
        shadelabel=False, fadelabel=False, **kwargs
    ):
        """
        Add up to 2 error indicators: more opaque "shading" and less opaque "fading".
        """
        kwargs, vert = _get_vert(**kwargs)
        shadestds = _not_none(shadestd=shadestd, shadestds=shadestds)
        fadestds = _not_none(fadestd=fadestd, fadestds=fadestds)
        shadepctiles = _not_none(shadepctile=shadepctile, shadepctiles=shadepctiles)
        fadepctiles = _not_none(fadepctile=fadepctile, fadepctiles=fadepctiles)
        shade = any(_ is not None for _ in (shadedata, shadestds, shadepctiles))
        fade = any(_ is not None for _ in (fadedata, fadestds, fadepctiles))

        # Shading properties
        shadeprops = _pop_props(kwargs, 'patch', prefix='shade')
        shadeprops.setdefault('alpha', 0.4)
        shadeprops.setdefault('zorder', 1.5)
        shadeprops.setdefault('linewidth', rc['patch.linewidth'])
        shadeprops.setdefault('edgecolor', 'none')
        # Fading properties
        fadeprops = _pop_props(kwargs, 'patch', prefix='fade')
        fadeprops.setdefault('zorder', shadeprops['zorder'])
        fadeprops.setdefault('alpha', 0.5 * shadeprops['alpha'])
        fadeprops.setdefault('linewidth', shadeprops['linewidth'])
        fadeprops.setdefault('edgecolor', 'none')
        # Get default color then apply to outgoing keyword args so
        # that plotting function will not advance to next cycler color.
        # TODO: More robust treatment of 'color' vs. 'facecolor'
        if (
            shade and shadeprops.get('facecolor', None) is None
            or fade and fadeprops.get('facecolor', None) is None
        ):
            color = kwargs.get(color_key, None)
            if color is None:  # add to outgoing
                color = kwargs[color_key] = self._get_lines.get_next_color()
            shadeprops.setdefault('facecolor', color)
            fadeprops.setdefault('facecolor', color)

        # Draw dark and light shading from distributions or explicit errdata
        eobjs = []
        fill = self.fill_between if vert else self.fill_betweenx
        if fade:
            edata, label = _distribution_range(
                y, distribution,
                stds=fadestds, pctiles=fadepctiles, errdata=fadedata,
                stds_default=(-3, 3), pctiles_default=(0, 100),
                label=fadelabel, absolute=True,
            )
            eobj = fill(x, *edata, label=label, **fadeprops)
            eobjs.append(eobj)
        if shade:
            edata, label = _distribution_range(
                y, distribution,
                stds=shadestds, pctiles=shadepctiles, errdata=shadedata,
                stds_default=(-2, 2), pctiles_default=(10, 90),
                label=shadelabel, absolute=True,
            )
            eobj = fill(x, *edata, label=label, **shadeprops)
            eobjs.append(eobj)

        kwargs['distribution'] = distribution
        return (*eobjs, kwargs)

    def _label_contours(
        self, obj, cobj, fmt, *, c=None, color=None, colors=None,
        size=None, fontsize=None, inline_spacing=None, **kwargs
    ):
        """
        Add labels to contours with support for shade-dependent filled contour labels.
        Text color is inferred from filled contour object and labels are always drawn
        on unfilled contour object (otherwise errors crop up).
        """
        # Parse input args
        colors = _not_none(c=c, color=color, colors=colors)
        fontsize = _not_none(size=size, fontsize=fontsize, default=rc['font.smallsize'])
        inline_spacing = _not_none(inline_spacing, 2.5)
        text_kw = {}
        clabel_keys = ('levels', 'inline', 'manual', 'rightside_up', 'use_clabeltext')
        for key in tuple(kwargs):  # allow dict to change size
            if key not in clabel_keys:
                text_kw[key] = kwargs.pop(key)

        # Draw hidden additional contour for filled contour labels
        cobj = _not_none(cobj, obj)
        colors = kwargs.pop('colors', None)
        if obj.filled and colors is None:
            colors = []
            for level in obj.levels:
                _, _, lum = to_xyz(obj.cmap(obj.norm(level)))
                colors.append('w' if lum < 50 else 'k')

        # Draw labels
        labs = cobj.clabel(
            fmt=fmt, colors=colors, fontsize=fontsize, inline_spacing=inline_spacing, **kwargs  # noqa: E501
        )
        if labs is not None:  # returns None if no contours
            for lab in labs:
                lab.update(text_kw)

        return labs

    def _label_gridboxes(
        self, obj, fmt, *, c=None, color=None, colors=None, size=None, fontsize=None, **kwargs  # noqa: E501
    ):
        """
        Add labels to pcolor boxes with support for shade-dependent text colors.
        Values are inferred from the unnormalized grid box color.
        """
        # Parse input args
        # NOTE: This function also hides grid boxes filled with NaNs to avoid ugly
        # issue where edge colors surround NaNs. Should maybe move this somewhere else.
        obj.update_scalarmappable()  # update 'edgecolors' list
        color = _not_none(c=c, color=color, colors=colors)
        fontsize = _not_none(size=size, fontsize=fontsize, default=rc['font.smallsize'])
        kwargs.setdefault('ha', 'center')
        kwargs.setdefault('va', 'center')

        # Apply colors and hide edge colors for empty grids
        # NOTE: Could also
        labs = []
        array = obj.get_array()
        paths = obj.get_paths()
        edgecolors = _to_numpy_array(obj.get_edgecolors())
        if len(edgecolors) == 1:
            edgecolors = np.repeat(edgecolors, len(array), axis=0)
        for i, (path, value) in enumerate(zip(paths, array)):
            # Round to the number corresponding to the *color* rather than
            # the exact data value. Similar to contour label numbering.
            if value is ma.masked or not np.isfinite(value):
                edgecolors[i, :] = 0
                continue
            if isinstance(obj.norm, pcolors.DiscreteNorm):
                value = obj.norm._norm.inverse(obj.norm(value))
            icolor = color
            if color is None:
                _, _, lum = to_xyz(obj.cmap(obj.norm(value)), 'hcl')
                icolor = 'w' if lum < 50 else 'k'
            bbox = path.get_extents()
            x = (bbox.xmin + bbox.xmax) / 2
            y = (bbox.ymin + bbox.ymax) / 2
            lab = self.text(x, y, fmt(value), color=icolor, size=fontsize, **kwargs)
            labs.append(lab)

        obj.set_edgecolors(edgecolors)
        return labs

    def _auto_labels(
        self, obj, cobj=None, labels=False, labels_kw=None,
        fmt=None, formatter=None, formatter_kw=None, precision=None,
    ):
        """
        Add number labels. Default formatter is `~proplot.ticker.SimpleFormatter`
        with a default maximum precision of ``3`` decimal places.
        """
        # TODO: Add quiverkey to this!
        if not labels:
            return
        labels_kw = labels_kw or {}
        formatter_kw = formatter_kw or {}
        formatter = _not_none(
            fmt_labels_kw=labels_kw.pop('fmt', None),
            formatter_labels_kw=labels_kw.pop('formatter', None),
            fmt=fmt,
            formatter=formatter,
            default='simple'
        )
        precision = _not_none(
            formatter_kw_precision=formatter_kw.pop('precision', None),
            precision=precision,
            default=3,  # should be lower than the default intended for tick labels
        )
        formatter = constructor.Formatter(formatter, precision=precision, **formatter_kw)  # noqa: E501
        if isinstance(obj, mcontour.ContourSet):
            self._label_contours(obj, cobj, formatter, **labels_kw)
        elif isinstance(obj, mcollections.Collection):
            self._label_gridboxes(obj, formatter, **labels_kw)
        else:
            raise RuntimeError(f'Not possible to add labels to object {obj!r}.')

    def _fix_edges(self, obj, edgefix=None, **kwargs):
        """
        Fix white lines between between filled contours and mesh and fix issues with
        colormaps that are transparent. Also takes collection-translated keyword
        args and if it detects any were passed then we skip this step.
        """
        # See: https://github.com/jklymak/contourfIssues
        # See: https://stackoverflow.com/q/15003353/4970632
        # NOTE: Use default edge width used for pcolor grid box edges. This is thick
        # enough to hide lines but thin enough to not add 'dots' to corners of boxes.
        edgefix = _not_none(edgefix, rc['cmap.edgefix'], True)
        linewidth = EDGEWIDTH if edgefix is True else 0 if edgefix is False else edgefix
        if not linewidth or not isinstance(obj, mcm.ScalarMappable):
            return
        if any(key in kwargs for key in ('linewidths', 'linestyles', 'edgecolors')):
            return

        # Remove edges when cmap has transparency
        cmap = obj.cmap
        if not cmap._isinit:
            cmap._init()
        if all(cmap._lut[:-1, 3] == 1):  # skip for cmaps with transparency
            edgecolor = 'face'
        else:
            edgecolor = 'none'

        # Apply fixes
        # NOTE: This also covers TriContourSet returned by tricontour
        if isinstance(obj, mcollections.Collection):
            obj.set_linewidth(linewidth)
            obj.set_edgecolor(edgecolor)
        if isinstance(obj, mcontour.ContourSet):
            if not obj.filled:
                return
            for contour in obj.collections:
                contour.set_linestyle('-')
                contour.set_linewidth(linewidth)
                contour.set_edgecolor(edgecolor)

    def _apply_plot(self, *pairs, vert=True, **kwargs):
        """
        Plot standard lines.
        """
        # Plot the lines
        objs = []
        kws = kwargs.copy()
        _process_props(kws, 'line')
        kws, extents = self._parse_inbounds(**kws)
        for xs, ys, fmt in self._iter_pairs(*pairs):
            xs, ys, kw = self._standardize_1d(xs, ys, vert=vert, **kws)
            ys, kw = _distribution_reduce(ys, **kw)
            guide_kw = _pop_params(kw, self._add_queued_guide)  # after standardize
            for _, n, x, y, kw in self._iter_columns(xs, ys, **kw):
                *eb, kw = self._error_bars(x, y, vert=vert, **kw)
                *es, kw = self._error_shading(x, y, vert=vert, **kw)
                kw = self._parse_cycle(n, **kw)
                if not vert:
                    x, y = y, x
                a = [x, y]
                if fmt is not None:  # x1, y1, fmt1, x2, y2, fm2... style input
                    a.append(fmt)
                obj, = self._plot_safe('plot', *a, **kw)
                self._restrict_inbounds(extents, x, y)
                objs.append((*eb, *es, obj) if eb or es else obj)

        # Add sticky edges
        axis = 'x' if vert else 'y'
        for obj in objs:
            if not isinstance(obj, mlines.Line2D):
                continue  # TODO: still looks good with error caps???
            data = getattr(obj, 'get_' + axis + 'data')()
            if not data.size:
                continue
            convert = getattr(self, 'convert_' + axis + 'units')
            edges = getattr(obj.sticky_edges, axis)
            min_, max_ = _safe_range(data)
            if min_ is not None:
                edges.append(convert(min_))
            if max_ is not None:
                edges.append(convert(max_))

        self._add_queued_guide(objs, **guide_kw)
        return objs  # always return list to match matplotlib behavior

    @_snippet_manager
    def line(self, *args, **kwargs):
        """
        %(plot.plot)s
        """
        return self.plot(*args, **kwargs)

    @_snippet_manager
    def linex(self, *args, **kwargs):
        """
        %(plot.plotx)s
        """
        return self.plotx(*args, **kwargs)

    @_preprocess_data('x', 'y', allow_extra=True)
    @docstring._concatenate_original
    @_snippet_manager
    def plot(self, *args, **kwargs):
        """
        %(plot.plot)s
        """
        kwargs = _parse_vert(default_vert=True, **kwargs)
        return self._apply_plot(*args, **kwargs)

    @_preprocess_data('y', 'x', allow_extra=True)
    @_snippet_manager
    def plotx(self, *args, **kwargs):
        """
        %(plot.plotx)s
        """
        kwargs = _parse_vert(default_vert=False, **kwargs)
        return self._apply_plot(*args, **kwargs)

    def _apply_step(self, *pairs, vert=True, where='pre', **kwargs):
        """
        Plot the steps.
        """
        # Plot the steps
        # NOTE: Internally matplotlib plot() calls step() so we could use that
        # approach... but instead repeat _apply_plot internals here so we can
        # disable error indications that make no sense for 'step' plots.
        kws = kwargs.copy()
        if where not in ('pre', 'post', 'mid'):
            raise ValueError(f"Invalid where={where!r}. Options are 'pre', 'post', 'mid'.")  # noqa: E501
        _process_props(kws, 'line')
        kws.setdefault('drawstyle', 'steps-' + where)
        kws, extents = self._parse_inbounds(**kws)
        objs = []
        for xs, ys, fmt in self._iter_pairs(*pairs):
            xs, ys, kw = self._standardize_1d(xs, ys, vert=vert, **kws)
            guide_kw = _pop_params(kw, self._add_queued_guide)  # after standardize
            if fmt is not None:
                kw['fmt'] = fmt
            for _, n, x, y, *a, kw in self._iter_columns(xs, ys, **kw):
                kw = self._parse_cycle(n, **kw)
                if not vert:
                    x, y = y, x
                obj, = self._plot_safe('step', x, y, *a, **kw)
                self._restrict_inbounds(extents, x, y)
                objs.append(obj)

        self._add_queued_guide(objs, **guide_kw)
        return objs  # always return list to match matplotlib behavior

    @_preprocess_data('x', 'y', allow_extra=True)
    @docstring._concatenate_original
    @_snippet_manager
    def step(self, *args, **kwargs):
        """
        %(plot.step)s
        """
        kwargs = _parse_vert(default_vert=True, **kwargs)
        return self._apply_step(*args, **kwargs)

    @_preprocess_data('y', 'x', allow_extra=True)
    @_snippet_manager
    def stepx(self, *args, **kwargs):
        """
        %(plot.stepx)s
        """
        kwargs = _parse_vert(default_vert=False, **kwargs)
        return self._apply_step(*args, **kwargs)

    def _apply_stem(
        self, x, y, *,
        linefmt=None, markerfmt=None, basefmt=None, orientation=None, **kwargs
    ):
        """
        Plot stem lines and markers.
        """
        # Parse input
        kw = kwargs.copy()
        kw, extents = self._parse_inbounds(**kw)
        x, y, kw = self._standardize_1d(x, y, **kw)
        guide_kw = _pop_params(kw, self._add_queued_guide)

        # Set default colors
        # NOTE: 'fmt' strings can only be 2 to 3 characters and include color
        # shorthands like 'r' or cycle colors like 'C0'. Cannot use full color names.
        # NOTE: Matplotlib defaults try to make a 'reddish' color the base and 'bluish'
        # color the stems. To make this more robust we temporarily replace the cycler.
        # Bizarrely stem() only reads from the global cycler() so have to update it.
        fmts = (linefmt, basefmt, markerfmt)
        orientation = _not_none(orientation, 'vertical')
        if not any(isinstance(fmt, str) and re.match(r'\AC[0-9]', fmt) for fmt in fmts):
            cycle = constructor.Cycle((rc['negcolor'], rc['poscolor']), name='_no_name')
            kw.setdefault('cycle', cycle)
        kw['basefmt'] = _not_none(basefmt, 'C1-')  # red base
        kw['linefmt'] = linefmt = _not_none(linefmt, 'C0-')  # blue stems
        kw['markerfmt'] = _not_none(markerfmt, linefmt[:-1] + 'o')  # blue marker
        sig = inspect.signature(maxes.Axes.stem)
        if 'use_line_collection' in sig.parameters:
            kw.setdefault('use_line_collection', True)

        # Call function then restore property cycle
        # WARNING: Horizontal stem plots are only supported in recent versions of
        # matplotlib. Let matplotlib raise an error if need be.
        ctx = {}
        cycle, kw = self._parse_cycle(return_cycle=True, **kw)  # allow re-application
        if cycle is not None:
            ctx['axes.prop_cycle'] = cycle
        if orientation == 'horizontal':  # may raise error
            kw['orientation'] = orientation
        with rc.context(ctx):
            obj = self._plot_safe('stem', x, y, **kw)
        self._restrict_inbounds(extents, x, y, orientation=orientation)
        self._add_queued_guide(obj, **guide_kw)
        return obj

    @_preprocess_data('x', 'y')
    @docstring._concatenate_original
    @_snippet_manager
    def stem(self, *args, **kwargs):
        """
        %(plot.stem)s
        """
        kwargs = _parse_vert(default_orientation='vertical', **kwargs)
        return self._apply_stem(*args, **kwargs)

    @_preprocess_data('x', 'y')
    @_snippet_manager
    def stemx(self, *args, **kwargs):
        """
        %(plot.stemx)s
        """
        kwargs = _parse_vert(default_orientation='horizontal', **kwargs)
        return self._apply_stem(*args, **kwargs)

    @_preprocess_data('x', 'y', ('c', 'color', 'colors', 'values'))
    @_snippet_manager
    def parametric(self, x, y, c, *, interp=0, scalex=True, scaley=True, **kwargs):
        """
        %(plot.parametric)s
        """
        # Standardize arguments
        # NOTE: Values are inferred in _auto_format() the same way legend labels are
        # inferred. Will not always return an array like inferred coordinates do.
        kw = kwargs.copy()
        _process_props(kw, 'collection')
        kw, extents = self._parse_inbounds(**kw)
        x, y, kw = self._standardize_1d(x, y, values=c, autovalues=True, autoreverse=False, **kw)  # noqa: E501
        c = kw.pop('values', None)  # permits inferring values e.g. a simple ordinate
        c = np.arange(y.size) if c is None else _to_numpy_array(c)
        c, colorbar_kw = _get_coords(c, which='')
        _guide_kw_to_arg('colorbar', kw, **colorbar_kw)
        _guide_kw_to_arg('colorbar', kw, locator=c)

        # Interpolate values to allow for smooth gradations between values or just
        # to color siwtchover halfway between points (interp True, False respectively)
        if interp > 0:
            x_orig, y_orig, v_orig = x, y, c
            x, y, c = [], [], []
            for j in range(x_orig.shape[0] - 1):
                idx = slice(None)
                if j + 1 < x_orig.shape[0] - 1:
                    idx = slice(None, -1)
                x.extend(np.linspace(x_orig[j], x_orig[j + 1], interp + 2)[idx].flat)
                y.extend(np.linspace(y_orig[j], y_orig[j + 1], interp + 2)[idx].flat)
                c.extend(np.linspace(v_orig[j], v_orig[j + 1], interp + 2)[idx].flat)  # noqa: E501
            x, y, c = np.array(x), np.array(y), np.array(c)

        # Get coordinates and values for points to the 'left' and 'right' of joints
        coords = []
        for i in range(y.shape[0]):
            icoords = np.empty((3, 2))
            for j, arr in enumerate((x, y)):
                icoords[0, j] = arr[0] if i == 0 else 0.5 * (arr[i - 1] + arr[i])
                icoords[1, j] = arr[i]
                icoords[2, j] = arr[-1] if i + 1 == y.shape[0] else 0.5 * (arr[i + 1] + arr[i])  # noqa: E501
            coords.append(icoords)
        coords = np.array(coords)

        # Get the colormap accounting for 'discrete' mode
        discrete = kw.get('discrete', None)
        if discrete is not None and not discrete:
            a = (x, y, c)  # pick levels from vmin and vmax, possibly limiting range
        else:
            a, kw['values'] = (), c
        kw = self._parse_cmap(*a, line_plot=True, **kw)
        cmap, norm = kw.pop('cmap'), kw.pop('norm')

        # Add collection with some custom attributes
        # NOTE: Modern API uses self._request_autoscale_view but this is
        # backwards compatible to earliest matplotlib versions.
        guide_kw = _pop_params(kw, self._add_queued_guide)
        obj = mcollections.LineCollection(
            coords, cmap=cmap, norm=norm,
            linestyles='-', capstyle='butt', joinstyle='miter',
        )
        obj.set_array(c)  # the ScalarMappable method
        obj.update({key: value for key, value in kw.items() if key not in ('color',)})
        self.add_collection(obj)
        self.autoscale_view(scalex=scalex, scaley=scaley)
        self._add_queued_guide(obj, **guide_kw)
        return obj

    def _apply_lines(
        self, xs, ys1, ys2, colors, *,
        vert=True, stack=None, stacked=None, negpos=False, **kwargs
    ):
        """
        Plot vertical or hotizontal lines at each point.
        """
        # Parse input arguments
        kw = kwargs.copy()
        name = 'vlines' if vert else 'hlines'
        if colors is not None:
            kw['colors'] = colors
        _process_props(kw, 'collection')
        kw, extents = self._parse_inbounds(**kw)
        stack = _not_none(stack=stack, stacked=stacked)
        xs, ys1, ys2, kw = self._standardize_1d(xs, ys1, ys2, vert=vert, **kw)
        guide_kw = _pop_params(kw, self._add_queued_guide)

        # Support "negative" and "positive" lines
        # TODO: Ensure 'linewidths' etc. are applied! For some reason
        # previously thought they had to be manually applied.
        y0 = 0
        objs, sides = [], []
        for _, n, x, y1, y2, kw in self._iter_columns(xs, ys1, ys2, **kw):
            kw = self._parse_cycle(n, **kw)
            if stack:
                y1 = y1 + y0  # avoid in-place modification
                y2 = y2 + y0
                y0 = y0 + y2 - y1  # irrelevant that we added y0 to both
            if negpos:
                obj = self._plot_negpos(name, x, y1, y2, colorkey='colors', **kw)
            else:
                obj = self._plot_safe(name, x, y1, y2, **kw)
            for y in (y1, y2):
                self._restrict_inbounds(extents, x, y, vert=vert)
                if y.size == 1:  # add sticky edges if bounds are scalar
                    sides.append(y)
            objs.append(obj)

        # Draw guide and add sticky edges
        self._add_sticky_edges(objs, 'y' if vert else 'x', *sides)
        self._add_queued_guide(objs, **guide_kw)
        return objs[0] if len(objs) == 1 else objs

    # WARNING: breaking change from native 'ymin' and 'ymax'
    @_preprocess_data('x', 'y1', 'y2', ('c', 'color', 'colors'))
    @_snippet_manager
    def vlines(self, *args, **kwargs):
        """
        %(plot.vlines)s
        """
        kwargs = _parse_vert(default_vert=True, **kwargs)
        return self._apply_lines(*args, **kwargs)

    # WARNING: breaking change from native 'xmin' and 'xmax'
    @_preprocess_data('y', 'x1', 'x2', ('c', 'color', 'colors'))
    @_snippet_manager
    def hlines(self, *args, **kwargs):
        """
        %(plot.hlines)s
        """
        kwargs = _parse_vert(default_vert=False, **kwargs)
        return self._apply_lines(*args, **kwargs)

    def _parse_markersize(self, s, *, smin=None, smax=None, **kwargs):
        """
        Scale the marker sizes with optional keyword args.
        """
        if np.atleast_1d(s).size == 1:  # None or scalar
            return s, kwargs
        smin_true, smax_true = _safe_range(s)
        smin_true = _not_none(smin_true, 0)
        smax_true = _not_none(smax_true, rc['lines.markersize'])
        smin = _not_none(smin, smin_true)
        smax = _not_none(smax, smax_true)
        s = smin + (smax - smin) * (s - smin_true) / (smax_true - smin_true)
        return s, kwargs

    def _apply_scatter(self, xs, ys, ss, cc, *, vert=True, **kwargs):
        """
        Apply scatter or scatterx markers.
        """
        # Apply from property cycle. Keys are cycle keys and values are scatter keys
        # NOTE: Matplotlib uses the property cycler in _get_patches_for_fill for
        # scatter() plots. It only ever inherits color from that. We instead use
        # _get_lines to help overarching goal of unifying plot() and scatter().
        cycle_manually = {
            'color': 'c',
            'markersize': 's',
            'linewidth': 'linewidths',
            'linestyle': 'linestyles',
            'markeredgewidth': 'linewidths',
            'markeredgecolor': 'edgecolors',
            'alpha': 'alpha',
            'marker': 'marker',
        }

        # Iterate over the columns
        kw = kwargs.copy()
        _process_props(kw, 'line')
        kw, extents = self._parse_inbounds(**kw)
        xs, ys, kw = self._standardize_1d(xs, ys, vert=vert, autoreverse=False, **kw)
        ss, kw = self._parse_markersize(ss, **kw)  # parse 's'
        cc, kw = self._parse_color(xs, ys, cc, apply_cycle=False, **kw)  # parse 'c'
        ys, kw = _distribution_reduce(ys, **kw)
        guide_kw = _pop_params(kw, self._add_queued_guide)
        objs = []
        for _, n, x, y, s, c, kw in self._iter_columns(xs, ys, ss, cc, **kw):
            kw['s'], kw['c'] = s, c  # make _parse_cycle() detect these
            *eb, kw = self._error_bars(x, y, vert=vert, **kw)
            *es, kw = self._error_shading(x, y, vert=vert, color_key='c', **kw)
            kw = self._parse_cycle(n, cycle_manually=cycle_manually, **kw)
            if not vert:
                x, y = y, x
            obj = self._plot_safe('scatter', x, y, **kw)
            self._restrict_inbounds(extents, x, y)
            objs.append((*eb, *es, obj) if eb or es else obj)

        self._add_queued_guide(objs, **guide_kw)
        return objs[0] if len(objs) == 1 else objs

    @_preprocess_data(
        'x', 'y', ('s', 'ms', 'markersize'), ('c', 'color', 'colors'),
        keywords=('lw', 'linewidth', 'linewidths', 'ec', 'edgecolor', 'edgecolors', 'fc', 'facecolor', 'facecolors')  # noqa: E501
    )
    @docstring._concatenate_original
    @_snippet_manager
    def scatter(self, *args, **kwargs):
        """
        %(plot.scatter)s
        """
        kwargs = _parse_vert(default_vert=True, **kwargs)
        return self._apply_scatter(*args, **kwargs)

    @_preprocess_data(
        'y', 'x', ('s', 'ms', 'markersize'), ('c', 'color', 'colors'),
        keywords=('lw', 'linewidth', 'linewidths', 'ec', 'edgecolor', 'edgecolors', 'fc', 'facecolor', 'facecolors')  # noqa: E501
    )
    @_snippet_manager
    def scatterx(self, *args, **kwargs):
        """
        %(plot.scatterx)s
        """
        kwargs = _parse_vert(default_vert=False, **kwargs)
        return self._apply_scatter(*args, **kwargs)

    def _apply_fill(
        self, xs, ys1, ys2, where, *,
        vert=True, negpos=None, stack=None, stacked=None, **kwargs
    ):
        """
        Apply area shading.
        """
        # Parse input arguments
        kw = kwargs.copy()
        _process_props(kw, 'patch')
        kw, extents = self._parse_inbounds(**kw)
        name = 'fill_between' if vert else 'fill_betweenx'
        stack = _not_none(stack=stack, stacked=stacked)
        xs, ys1, ys2, kw = self._standardize_1d(xs, ys1, ys2, vert=vert, **kw)

        # Draw patches with default edge width zero
        y0 = 0
        objs, xsides, ysides = [], [], []
        guide_kw = _pop_params(kw, self._add_queued_guide)
        for _, n, x, y1, y2, w, kw in self._iter_columns(xs, ys1, ys2, where, **kw):
            kw = self._parse_cycle(n, **kw)
            if stack:
                y1 = y1 + y0  # avoid in-place modification
                y2 = y2 + y0
                y0 = y0 + y2 - y1  # irrelevant that we added y0 to both
            if negpos:
                # NOTE: pass 'where' so plot_negpos can ignore it and issue a warning
                obj = self._plot_negpos(name, x, y1, y2, where=w, use_where=True, **kw)
            else:
                obj = self._plot_safe(name, x, y1, y2, where=w, **kw)
            xsides.append(x)
            for y in (y1, y2):
                self._restrict_inbounds(extents, x, y, vert=vert)
                if y.size == 1:  # add sticky edges if bounds are scalar
                    ysides.append(y)
            objs.append(obj)

        # Draw guide and add sticky edges
        self._add_queued_guide(objs, **guide_kw)
        for axis, sides in zip('xy' if vert else 'yx', (xsides, ysides)):
            self._add_sticky_edges(objs, axis, *sides)
        return objs[0] if len(objs) == 1 else objs

    @_snippet_manager
    def area(self, *args, **kwargs):
        """
        %(plot.fill_between)s
        """
        return self.fill_between(*args, **kwargs)

    @_snippet_manager
    def areax(self, *args, **kwargs):
        """
        %(plot.fill_betweenx)s
        """
        return self.fill_betweenx(*args, **kwargs)

    @_preprocess_data('x', 'y1', 'y2', 'where')
    @docstring._concatenate_original
    @_snippet_manager
    def fill_between(self, *args, **kwargs):
        """
        %(plot.fill_between)s
        """
        kwargs = _parse_vert(default_vert=True, **kwargs)
        return self._apply_fill(*args, **kwargs)

    @_preprocess_data('y', 'x1', 'x2', 'where')
    @docstring._concatenate_original
    @_snippet_manager
    def fill_betweenx(self, *args, **kwargs):
        """
        %(plot.fill_betweenx)s
        """
        # NOTE: The 'horizontal' orientation will be inferred by downstream
        # wrappers using the function name.
        kwargs = _parse_vert(default_vert=False, **kwargs)
        return self._apply_fill(*args, **kwargs)

    @staticmethod
    def _convert_bar_width(x, width=1):
        """
        Convert bar plot widths from relative to coordinate spacing. Relative
        widths are much more convenient for users.
        """
        # WARNING: This will fail for non-numeric non-datetime64 singleton
        # datatypes but this is good enough for vast majority of cases.
        x_test = _to_numpy_array(x)
        if len(x_test) >= 2:
            x_step = x_test[1:] - x_test[:-1]
            x_step = np.concatenate((x_step, x_step[-1:]))
        elif x_test.dtype == np.datetime64:
            x_step = np.timedelta64(1, 'D')
        else:
            x_step = np.array(0.5)
        if np.issubdtype(x_test.dtype, np.datetime64):
            # Avoid integer timedelta truncation
            x_step = x_step.astype('timedelta64[ns]')
        return width * x_step

    def _apply_bar(
        self, xs, hs, ws, bs, *, absolute_width=False,
        stack=None, stacked=None, negpos=False, orientation='vertical', **kwargs
    ):
        """
        Apply bar or barh command. Support default "minima" at zero.
        """
        # Parse args
        kw = kwargs.copy()
        kw, extents = self._parse_inbounds(**kw)
        name = 'barh' if orientation == 'horizontal' else 'bar'
        stack = _not_none(stack=stack, stacked=stacked)
        xs, hs, kw = self._standardize_1d(xs, hs, orientation=orientation, **kw)

        # Call func after converting bar width
        b0 = 0
        objs = []
        _process_props(kw, 'patch')
        kw.setdefault('edgecolor', 'black')
        hs, kw = _distribution_reduce(hs, **kw)
        guide_kw = _pop_params(kw, self._add_queued_guide)
        for i, n, x, h, w, b, kw in self._iter_columns(xs, hs, ws, bs, **kw):
            kw = self._parse_cycle(n, **kw)
            # Adjust x or y coordinates for grouped and stacked bars
            w = _not_none(w, np.array([0.8]))  # same as mpl but in *relative* units
            b = _not_none(b, np.array([0.0]))  # same as mpl
            if not absolute_width:
                w = self._convert_bar_width(x, w)
            if stack:
                b = b + b0
                b0 = b0 + h
            else:  # instead "group" the bars (this is no-op if we have 1 column)
                w = w / n  # rescaled
                o = 0.5 * (n - 1)  # center coordinate
                x = x + w * (i - o)  # += may cause integer/float casting issue
            # Draw simple bars
            *eb, kw = self._error_bars(x, b + h, orientation=orientation, **kw)
            if negpos:
                obj = self._plot_negpos(name, x, h, w, b, use_zero=True, **kw)
            else:
                obj = self._plot_safe(name, x, h, w, b, **kw)
            for y in (b, b + h):
                self._restrict_inbounds(extents, x, y, orientation=orientation)
            objs.append((*eb, obj) if eb else obj)

        self._add_queued_guide(objs, **guide_kw)
        return objs[0] if len(objs) == 1 else objs

    @_preprocess_data('x', 'height', 'width', 'bottom')
    @docstring._concatenate_original
    @_snippet_manager
    def bar(self, *args, **kwargs):
        """
        %(plot.bar)s
        """
        kwargs = _parse_vert(default_orientation='vertical', **kwargs)
        return self._apply_bar(*args, **kwargs)

    # WARNING: Swap 'height' and 'width' here so that they are always relative
    # to the 'tall' axis. This lets people always pass 'width' as keyword
    @_preprocess_data('y', 'height', 'width', 'left')
    @docstring._concatenate_original
    @_snippet_manager
    def barh(self, *args, **kwargs):
        """
        %(plot.barh)s
        """
        kwargs = _parse_vert(default_orientation='horizontal', **kwargs)
        return self._apply_bar(*args, **kwargs)

    def _apply_boxplot(
        self, x, y, *,
        mean=None, means=None, vert=True,
        fill=None, marker=None, markersize=None,
        **kwargs
    ):
        """
        Apply the box plot.
        """
        # Global and fill properties
        kw = kwargs.copy()
        _process_props(kw, 'patch')
        linewidth = kw.pop('linewidth', rc['patch.linewidth'])
        edgecolor = kw.pop('edgecolor', 'black')
        fillcolor = kw.pop('facecolor', None)
        fillalpha = kw.pop('alpha', None)
        fill = fill is True or fillcolor is not None or fillalpha is not None
        if fill and fillcolor is None:  # TODO: support e.g. 'facecolor' cycle?
            parser = self._get_patches_for_fill
            fillcolor = parser.get_next_color()
        fillalpha = _not_none(fillalpha, 1)

        # Arist-specific properties
        # NOTE: Output dict keys are plural but we use singular for keyword args
        props = {}
        for key in ('boxes', 'whiskers', 'caps', 'fliers', 'medians', 'means'):
            prefix = key.rstrip('es')  # singular form
            props[key] = iprops = _pop_props(kw, 'line', prefix=prefix)
            iprops.setdefault('color', edgecolor)
            iprops.setdefault('linewidth', linewidth)
            iprops.setdefault('markeredgecolor', edgecolor)
        means = _not_none(mean=mean, means=means, showmeans=kw.get('showmeans'))
        if means:
            kw['showmeans'] = kw['meanline'] = True

        # Call function
        x, y, kw = self._standardize_1d(x, y, autoy=False, autoguide=False, vert=vert, **kw)  # noqa: E501
        kw.setdefault('positions', x)
        obj = self._plot_safe('boxplot', y, vert=vert, **kw)

        # Modify artist settings
        for key, aprops in props.items():
            if key not in obj:  # possible if not rendered
                continue
            artists = obj[key]
            if not isinstance(fillalpha, list):
                fillalpha = [fillalpha] * len(artists)
            if not isinstance(fillcolor, list):
                fillcolor = [fillcolor] * len(artists)
            for i, artist in enumerate(artists):
                # Update lines used for boxplot components
                # TODO: Test this thoroughly!
                iprops = {
                    key: (
                        value[i // 2 if key in ('caps', 'whiskers') else i]
                        if isinstance(value, (list, ndarray))
                        else value
                    ) for key, value in aprops.items()
                }
                artist.update(iprops)
                # "Filled" boxplot by adding patch beneath line path
                if key == 'boxes':
                    ifillcolor = fillcolor[i]  # must stay within the if statement
                    ifillalpha = fillalpha[i]
                    if ifillcolor is not None or ifillalpha is not None:
                        patch = mpatches.PathPatch(
                            artist.get_path(),
                            linewidth=0, facecolor=ifillcolor, alpha=ifillalpha,
                        )
                        self.add_artist(patch)
                # Outlier markers
                if key == 'fliers':
                    if marker is not None:
                        artist.set_marker(marker)
                    if markersize is not None:
                        artist.set_markersize(markersize)

        return obj

    @_snippet_manager
    def box(self, *args, **kwargs):
        """
        %(plot.boxplot)s
        """
        return self.boxplot(*args, **kwargs)

    @_snippet_manager
    def boxh(self, *args, **kwargs):
        """
        %(plot.boxploth)s
        """
        return self.boxploth(*args, **kwargs)

    @_preprocess_data('positions', 'y')
    @docstring._concatenate_original
    @_snippet_manager
    def boxplot(self, *args, **kwargs):
        """
        %(plot.boxplot)s
        """
        kwargs = _parse_vert(default_vert=True, **kwargs)
        return self._apply_boxplot(*args, **kwargs)

    @_preprocess_data('positions', 'x')
    @_snippet_manager
    def boxploth(self, *args, **kwargs):
        """
        %(plot.boxploth)s
        """
        kwargs = _parse_vert(default_vert=False, **kwargs)
        return self._apply_boxplot(*args, **kwargs)

    def _apply_violinplot(self, x, y, vert=True, **kwargs):
        """
        Apply the violinplot.
        """
        # Parse keyword args
        kw = kwargs.copy()
        _process_props(kw, 'patch')
        linewidth = kw.pop('linewidth', rc['patch.linewidth'])
        edgecolor = kw.pop('edgecolor', 'black')
        fillcolor = kw.pop('facecolor', None)
        fillalpha = kw.pop('alpha', None)
        fillalpha = _not_none(fillalpha, 1)
        kw.setdefault('capsize', 0)  # caps are redundant for violin plots
        kw.setdefault('means', kw.pop('showmeans', None))  # for _indicate_error
        kw.setdefault('medians', kw.pop('showmedians', None))
        if kw.pop('showextrema', None):
            warnings._warn_proplot('Ignoring showextrema=True.')

        # Parse and control error bars
        x, y, kw = self._standardize_1d(x, y, autoy=False, autoguide=False, vert=vert, **kw)  # noqa: E501
        y, kw = _distribution_reduce(y, **kw)
        *eb, kw = self._error_bars(x, y, vert=vert, default_boxes=True, **kw)  # noqa: E501
        kw = self._parse_cycle(**kw)

        # Call function
        kw.pop('labels', None)  # already applied in _standardize_1d
        kw.update({'showmeans': False, 'showmedians': False, 'showextrema': False})
        kw.setdefault('positions', x)
        y = kw.pop('distribution', None)  # 'y' was changes in _distribution_reduce
        obj = self._plot_safe('violinplot', y, vert=vert, **kw)

        # Modify body settings
        artists = (obj or {}).get('bodies', ())
        if not isinstance(fillalpha, list):
            fillalpha = [fillalpha] * len(artists)
        if not isinstance(fillcolor, list):
            fillcolor = [fillcolor] * len(artists)
        if not isinstance(edgecolor, list):
            edgecolor = [edgecolor] * len(artists)
        for i, artist in enumerate(artists):
            artist.set_linewidths(linewidth)
            if fillalpha[i] is not None:
                artist.set_alpha(fillalpha[i])
            if fillcolor[i] is not None:
                artist.set_facecolor(fillcolor[i])
            if edgecolor[i] is not None:
                artist.set_edgecolor(edgecolor[i])

        return obj

    @_snippet_manager
    def violin(self, *args, **kwargs):
        """
        %(plot.violinplot)s
        """
        # WARNING: This disables use of 'violin' by users but
        # probably very few people use this anyway.
        if getattr(self, '_internal_call', None):
            return super().violin(*args, **kwargs)
        else:
            return self.violinplot(*args, **kwargs)

    @_snippet_manager
    def violinh(self, *args, **kwargs):
        """
        %(plot.violinploth)s
        """
        return self.violinploth(*args, **kwargs)

    @_preprocess_data('positions', 'y')
    @docstring._concatenate_original
    @_snippet_manager
    def violinplot(self, *args, **kwargs):
        """
        %(plot.violinplot)s
        """
        kwargs = _parse_vert(default_vert=True, **kwargs)
        return self._apply_violinplot(*args, **kwargs)

    @_preprocess_data('positions', 'x')
    @_snippet_manager
    def violinploth(self, *args, **kwargs):
        """
        %(plot.violinploth)s
        """
        kwargs = _parse_vert(default_vert=False, **kwargs)
        return self._apply_violinplot(*args, **kwargs)

    def _apply_hist(self, xs, bins, *, orientation='vertical', **kwargs):
        """
        Apply the histogram.
        """
        # WARNING: Weirdly while Axes.bar() adds labels to the container
        # Axes.hist() adds them to the first element in the container. The
        # legend handle reader just looks for items with get_label() so we
        # manually apply labels to the container on the result.
        _, xs, kw = self._standardize_1d(xs, orientation=orientation, **kwargs)
        objs = []
        guide_kw = _pop_params(kw, self._add_queued_guide)
        if bins is not None:
            kw['bins'] = bins
        for _, n, x, kw in self._iter_columns(xs, **kw):
            kw = self._parse_cycle(n, **kw)
            obj = self._plot_safe('hist', x, orientation=orientation, **kw)
            if 'label' in kw:
                for arg in obj[2]:
                    arg.set_label(kw['label'])
                if hasattr(obj[2], 'set_label'):  # recent mpl versions
                    obj[2].set_label(kw['label'])
            objs.append(obj)
        self._add_queued_guide(objs, **guide_kw)
        return objs[0] if len(objs) == 1 else objs

    @_preprocess_data('x', 'bins', keywords='weights')
    @docstring._concatenate_original
    @_snippet_manager
    def hist(self, *args, **kwargs):
        """
        %(plot.hist)s
        """
        kwargs = _parse_vert(default_orientation='vertical', **kwargs)
        return self._apply_hist(*args, **kwargs)

    @_preprocess_data('y', 'bins', keywords='weights')
    @_snippet_manager
    def histh(self, *args, **kwargs):
        """
        %(plot.histh)s
        """
        kwargs = _parse_vert(default_orientation='horizontal', **kwargs)
        return self._apply_hist(*args, **kwargs)

    # WARNING: 'labels' and 'colors' no longer passed through `data` (seems like
    # extremely niche usage... `data` variables should be data-like)
    @_preprocess_data('x', 'explode')
    @docstring._concatenate_original
    @_snippet_manager
    def pie(self, x, explode, *, labelpad=None, labeldistance=None, **kwargs):
        """
        Plot a pie chart.

        Parameters
        ----------
        %(plot.args_1d_y)s
        %(plot.args_1d_shared)s

        Other parameters
        ----------------
        %(plot.cycle)s
        %(plot.labels_1d)s
        labelpad, labeldistance : float, optional
            The distance at which labels are drawn in radial coordinates.
        lw, linewidth, linewidths : float, optional
            The edge width for the pie sectors.
        ec, edgecolor, edgecolors : color-spec, optional
            The edge color for the pie sectors.

        See also
        --------
        matplotlib.axes.Axes.pie
        """
        pad = _not_none(labeldistance=labeldistance, labelpad=labelpad, default=1.15)
        props = _pop_props(kwargs, 'patch')
        props.setdefault('edgecolor', 'k')  # sensible default
        _, x, kwargs = self._standardize_1d(
            x, autox=False, autoy=False, **kwargs
        )
        kwargs = self._parse_cycle(**kwargs)
        kwargs['labeldistance'] = pad
        obj = self._plot_safe('pie', x, explode, wedgeprops=props, **kwargs)
        return obj

    @_preprocess_data('x', 'y', 'bins', keywords='weights')
    @docstring._concatenate_original
    @_snippet_manager
    def hist2d(self, x, y, bins, **kwargs):
        """
        Plot a standard 2D histogram.

        Parameters
        ----------
        %(plot.args_1d_y)s
        bins : int or 2-tuple of int, or array or 2-tuple of array, optional
            The bin count or list of bins for each dimension or both dimensions.
        %(plot.weights)s
        %(plot.args_1d_shared)s

        Other parameters
        ----------------
        %(plot.cmap_norm)s
        %(plot.levels_manual)s
        %(plot.levels_vlim)s
        %(plot.levels_auto)s
        %(plot.labels_2d)s
        %(plot.guide)s
        **kwargs
            Passed to `~matplotlib.axes.Axes.hist2d`.

        See also
        --------
        PlotAxes.hist2d
        matplotlib.axes.Axes.hexbin
        """
        # Rely on pcolormesh() override for this.
        if bins is not None:
            kwargs['bins'] = bins
        return super().hist2d(x, y, **kwargs)

    # WARNING: breaking change from native 'C'
    @_preprocess_data('x', 'y', 'weights')
    @docstring._concatenate_original
    @_snippet_manager
    def hexbin(self, x, y, weights, **kwargs):
        """
        Plot a 2D hexagonally binned histogram.

        Parameters
        ----------
        %(plot.args_1d_y)s
        %(plot.weights)s
        %(plot.args_1d_shared)s

        Other parameters
        ----------------
        %(plot.cmap_norm)s
        %(plot.levels_manual)s
        %(plot.levels_vlim)s
        %(plot.labels_2d)s
        %(plot.guide)s
        **kwargs
            Passed to `~matplotlib.axes.Axes.hexbin`.

        See also
        --------
        PlotAxes.hist2d
        matplotlib.axes.Axes.hexbin
        """
        # WARNING: Cannot use automatic level generation here until counts are
        # estimated. Inside _parse_levels if no manual levels were provided then
        # _parse_autolev is skipped and args like levels=10 or locator=5 are ignored
        x, y, kw = self._standardize_1d(x, y, autovalues=True, **kwargs)
        _process_props(kw, 'collection')  # takes LineCollection props
        kw = self._parse_cmap(x, y, y, skip_autolev=True, default_discrete=False, **kw)
        norm = kw.get('norm', None)
        if norm is not None and not isinstance(norm, pcolors.DiscreteNorm):
            norm.vmin = norm.vmax = None  # remove nonsense values
        labels_kw = _pop_params(kw, self._auto_labels)
        guide_kw = _pop_params(kw, self._add_queued_guide)
        m = self._plot_safe('hexbin', x, y, weights, **kw)
        self._auto_labels(m, **labels_kw)
        self._add_queued_guide(m, **guide_kw)
        return m

    @_preprocess_data('x', 'y', 'z')
    @docstring._concatenate_original
    @_snippet_manager
    def contour(self, x, y, z, **kwargs):
        """
        %(plot.contour)s
        """
        x, y, z, kw = self._standardize_2d(x, y, z, **kwargs)
        _process_props(kw, 'collection')
        kw = self._parse_cmap(x, y, z, minlength=1, contour_plot=True, **kw)
        cmap = kw.pop('cmap', None)
        if isinstance(cmap, pcolors.DiscreteColormap) and len(set(cmap.colors)) == 1:
            kw['colors'] = cmap.colors[0]  # otherwise negative linestyle fails
        else:
            kw['cmap'] = cmap
        labels_kw = _pop_params(kw, self._auto_labels)
        guide_kw = _pop_params(kw, self._add_queued_guide)
        label = kw.pop('label', None)
        m = self._plot_safe('contour', x, y, z, **kw)
        m._legend_label = label
        self._auto_labels(m, **labels_kw)
        self._add_queued_guide(m, **guide_kw)
        return m

    @_preprocess_data('x', 'y', 'z')
    @docstring._concatenate_original
    @_snippet_manager
    def contourf(self, x, y, z, **kwargs):
        """
        %(plot.contourf)s
        """
        x, y, z, kw = self._standardize_2d(x, y, z, **kwargs)
        _process_props(kw, 'collection')
        kw = self._parse_cmap(x, y, z, contour_plot=True, **kw)
        contour_kw = _pop_kwargs(kw, 'edgecolors', 'linewidths', 'linestyles')
        edgefix_kw = _pop_params(kw, self._fix_edges)
        labels_kw = _pop_params(kw, self._auto_labels)
        guide_kw = _pop_params(kw, self._add_queued_guide)
        label = kw.pop('label', None)
        m = cm = self._plot_safe('contourf', x, y, z, **kw)
        m._legend_label = label
        self._fix_edges(m, **edgefix_kw, **contour_kw)  # skipped if bool(contour_kw)
        if contour_kw or labels_kw:
            cm = self._plot_edges('contour', x, y, z, **kw, **contour_kw)
        self._auto_labels(m, cm, **labels_kw)
        self._add_queued_guide(m, **guide_kw)
        return m

    @_preprocess_data('x', 'y', 'z')
    @docstring._concatenate_original
    @_snippet_manager
    def pcolor(self, x, y, z, **kwargs):
        """
        %(plot.pcolor)s
        """
        x, y, z, kw = self._standardize_2d(x, y, z, edges=True, **kwargs)
        _process_props(kw, 'collection')
        kw = self._parse_cmap(x, y, z, to_centers=True, **kw)
        edgefix_kw = _pop_params(kw, self._fix_edges)
        labels_kw = _pop_params(kw, self._auto_labels)
        guide_kw = _pop_params(kw, self._add_queued_guide)
        m = self._plot_safe('pcolor', x, y, z, **kw)
        self._fix_edges(m, **edgefix_kw, **kw)
        self._auto_labels(m, **labels_kw)
        self._add_queued_guide(m, **guide_kw)
        return m

    @_preprocess_data('x', 'y', 'z')
    @docstring._concatenate_original
    @_snippet_manager
    def pcolormesh(self, x, y, z, **kwargs):
        """
        %(plot.pcolormesh)s
        """
        x, y, z, kw = self._standardize_2d(x, y, z, edges=True, **kwargs)
        _process_props(kw, 'collection')
        kw = self._parse_cmap(x, y, z, to_centers=True, **kw)
        edgefix_kw = _pop_params(kw, self._fix_edges)
        labels_kw = _pop_params(kw, self._auto_labels)
        guide_kw = _pop_params(kw, self._add_queued_guide)
        m = self._plot_safe('pcolormesh', x, y, z, **kw)
        self._fix_edges(m, **edgefix_kw, **kw)
        self._auto_labels(m, **labels_kw)
        self._add_queued_guide(m, **guide_kw)
        return m

    @_preprocess_data('x', 'y', 'z')
    @docstring._concatenate_original
    @_snippet_manager
    def pcolorfast(self, x, y, z, **kwargs):
        """
        %(plot.pcolorfast)s
        """
        x, y, z, kw = self._standardize_2d(x, y, z, edges=True, **kwargs)
        _process_props(kw, 'collection')
        kw = self._parse_cmap(x, y, z, to_centers=True, **kw)
        edgefix_kw = _pop_params(kw, self._fix_edges)
        labels_kw = _pop_params(kw, self._auto_labels)
        guide_kw = _pop_params(kw, self._add_queued_guide)
        m = self._plot_safe('pcolorfast', x, y, z, **kw)
        self._fix_edges(m, **edgefix_kw, **kw)
        self._auto_labels(m, **labels_kw)
        self._add_queued_guide(m, **guide_kw)
        return m

    @_snippet_manager
    def heatmap(self, *args, aspect=None, **kwargs):
        """
        %(plot.heatmap)s
        """
        obj = self.pcolormesh(*args, default_discrete=False, **kwargs)
        aspect = _not_none(aspect, rc['image.aspect'])
        if self.name != 'proplot_cartesian':
            warnings._warn_proplot(
                'The heatmap() command is meant for CartesianAxes. '
                'Please use pcolor() or pcolormesh() instead.'
            )
        else:
            coords = getattr(obj, '_coordinates', None)
            xlocator = ylocator = None
            if coords is not None:
                coords = 0.5 * (coords[1:, ...] + coords[:-1, ...])
                coords = 0.5 * (coords[:, 1:, :] + coords[:, :-1, :])
                xlocator, ylocator = coords[0, :, 0], coords[:, 0, 1]
            kw = {'aspect': aspect, 'xgrid': False, 'ygrid': False}
            if xlocator is not None and self.xaxis.isDefault_majloc:
                kw['xlocator'] = xlocator
            if ylocator is not None and self.yaxis.isDefault_majloc:
                kw['ylocator'] = ylocator
            if self.xaxis.isDefault_minloc:
                kw['xtickminor'] = False
            if self.yaxis.isDefault_minloc:
                kw['ytickminor'] = False
            self.format(**kw)
        return obj

    @_preprocess_data('x', 'y', 'u', 'v', ('c', 'color', 'colors'))
    @docstring._concatenate_original
    @_snippet_manager
    def barbs(self, x, y, u, v, c, **kwargs):
        """
        %(plot.barbs)s
        """
        x, y, u, v, kw = self._standardize_2d(x, y, u, v, allow1d=True, autoguide=False, **kwargs)  # noqa: E501
        _process_props(kw, 'line')  # applied to barbs
        c, kw = self._parse_color(x, y, c, **kw)
        if mcolors.is_color_like(c):
            kw['barbcolor'], c = c, None
        a = [x, y, u, v]
        if c is not None:
            a.append(c)
        kw.pop('colorbar_kw', None)  # added by _parse_cmap
        m = self._plot_safe('barbs', *a, **kw)
        return m

    @_preprocess_data('x', 'y', 'u', 'v', ('c', 'color', 'colors'))
    @docstring._concatenate_original
    @_snippet_manager
    def quiver(self, x, y, u, v, c, **kwargs):
        """
        %(plot.quiver)s
        """
        x, y, u, v, kw = self._standardize_2d(x, y, u, v, allow1d=True, autoguide=False, **kwargs)  # noqa: E501
        _process_props(kw, 'line')  # applied to arrow outline
        c, kw = self._parse_color(x, y, c, **kw)
        color = None
        if mcolors.is_color_like(c):
            color, c = c, None
        if color is not None:
            kw['color'] = color
        a = [x, y, u, v]
        if c is not None:
            a.append(c)
        kw.pop('colorbar_kw', None)  # added by _parse_cmap
        m = self._plot_safe('quiver', *a, **kw)
        return m

    @_snippet_manager
    def stream(self, *args, **kwargs):
        """
        %(plot.stream)s
        """
        return self.streamplot(*args, **kwargs)

    # WARNING: breaking change from native streamplot() fifth positional arg 'density'
    @_preprocess_data('x', 'y', 'u', 'v', ('c', 'color', 'colors'), keywords='start_points')  # noqa: E501
    @docstring._concatenate_original
    @_snippet_manager
    def streamplot(self, x, y, u, v, c, **kwargs):
        """
        %(plot.stream)s
        """
        x, y, u, v, kw = self._standardize_2d(x, y, u, v, **kwargs)
        _process_props(kw, 'line')  # applied to lines
        c, kw = self._parse_color(x, y, c, **kw)
        if c is None:  # throws an error if color not provided
            c = pcolors.to_hex(self._get_lines.get_next_color())
        kw['color'] = c  # always pass this
        guide_kw = _pop_params(kw, self._add_queued_guide)
        label = kw.pop('label', None)
        m = self._plot_safe('streamplot', x, y, u, v, **kw)
        m.lines.set_label(label)  # the collection label
        self._add_queued_guide(m.lines, **guide_kw)  # lines inside StreamplotSet
        return m

    @_preprocess_data('x', 'y', 'z')
    @docstring._concatenate_original
    @_snippet_manager
    def tricontour(self, x, y, z, **kwargs):
        """
        %(plot.tricontour)s
        """
        kw = kwargs.copy()
        if x is None or y is None or z is None:
            raise ValueError('Three input arguments are required.')
        _process_props(kw, 'collection')
        kw = self._parse_cmap(x, y, z, minlength=1, contour_plot=True, **kw)
        cmap = kw.pop('cmap', None)
        if isinstance(cmap, pcolors.DiscreteColormap) and len(set(cmap.colors)) == 1:
            kw['colors'] = cmap.colors[0]  # otherwise negative linestyle fails
        else:
            kw['cmap'] = cmap
        labels_kw = _pop_params(kw, self._auto_labels)
        guide_kw = _pop_params(kw, self._add_queued_guide)
        label = kw.pop('label', None)
        m = self._plot_safe('tricontour', x, y, z, **kw)
        m._legend_label = label
        self._auto_labels(m, **labels_kw)
        self._add_queued_guide(m, **guide_kw)
        return m

    @_preprocess_data('x', 'y', 'z')
    @docstring._concatenate_original
    @_snippet_manager
    def tricontourf(self, x, y, z, **kwargs):
        """
        %(plot.tricontourf)s
        """
        kw = kwargs.copy()
        if x is None or y is None or z is None:
            raise ValueError('Three input arguments are required.')
        _process_props(kw, 'collection')
        contour_kw = _pop_kwargs(kw, 'edgecolors', 'linewidths', 'linestyles')
        kw = self._parse_cmap(x, y, z, contour_plot=True, **kw)
        edgefix_kw = _pop_params(kw, self._fix_edges)
        labels_kw = _pop_params(kw, self._auto_labels)
        guide_kw = _pop_params(kw, self._add_queued_guide)
        label = kw.pop('label', None)
        m = cm = self._plot_safe('tricontourf', x, y, z, **kw)
        m._legend_label = label
        self._fix_edges(m, **edgefix_kw, **contour_kw)  # skipped if bool(contour_kw)
        if contour_kw or labels_kw:
            cm = self._plot_edges('tricontour', x, y, z, **kw, **contour_kw)
        self._auto_labels(m, cm, **labels_kw)
        self._add_queued_guide(m, **guide_kw)
        return m

    @_preprocess_data('x', 'y', 'z')
    @docstring._concatenate_original
    @_snippet_manager
    def tripcolor(self, x, y, z, **kwargs):
        """
        %(plot.tripcolor)s
        """
        kw = kwargs.copy()
        if x is None or y is None or z is None:
            raise ValueError('Three input arguments are required.')
        _process_props(kw, 'collection')
        kw = self._parse_cmap(x, y, z, **kw)
        edgefix_kw = _pop_params(kw, self._fix_edges)
        labels_kw = _pop_params(kw, self._auto_labels)
        guide_kw = _pop_params(kw, self._add_queued_guide)
        m = self._plot_safe('tripcolor', x, y, z, **kw)
        self._fix_edges(m, **edgefix_kw, **kw)
        self._auto_labels(m, **labels_kw)
        self._add_queued_guide(m, **guide_kw)
        return m

    # WARNING: breaking change from native 'X'
    @_preprocess_data('z')
    @docstring._concatenate_original
    @_snippet_manager
    def imshow(self, z, **kwargs):
        """
        %(plot.imshow)s
        """
        kw = kwargs.copy()
        kw = self._parse_cmap(z, default_discrete=False, **kw)
        guide_kw = _pop_params(kw, self._add_queued_guide)
        m = self._plot_safe('imshow', z, **kw)
        self._add_queued_guide(m, **guide_kw)
        return m

    # WARNING: breaking change from native 'Z'
    @_preprocess_data('z')
    @docstring._concatenate_original
    @_snippet_manager
    def matshow(self, z, **kwargs):
        """
        %(plot.matshow)s
        """
        kw = kwargs.copy()
        kw = self._parse_cmap(z, **kw)
        guide_kw = _pop_params(kw, self._add_queued_guide)
        m = self._plot_safe('matshow', z, **kw)
        self._add_queued_guide(m, **guide_kw)
        return m

    # WARNING: breaking change from native 'Z'
    @_preprocess_data('z')
    @docstring._concatenate_original
    @_snippet_manager
    def spy(self, z, **kwargs):
        """
        %(plot.spy)s
        """
        kw = kwargs.copy()
        _process_props(kw, 'line')  # takes valid Line2D properties
        default_cmap = pcolors.DiscreteColormap(['w', 'k'], '_no_name')
        kw = self._parse_cmap(z, default_cmap=default_cmap, **kw)
        guide_kw = _pop_params(kw, self._add_queued_guide)
        m = self._plot_safe('spy', z, **kw)
        self._add_queued_guide(m, **guide_kw)
        return m

    def set_prop_cycle(self, *args, **kwargs):
        # Silent override. This is a strict superset of matplotlib functionality
        # with one exception: you cannot use e.g. set_prop_cycle('color', color_list).
        # Instead keyword args are required (but note naked positional arguments
        # are assumed color arguments). Cycles are still validated in rcsetup.cycler()
        cycle = self._active_cycle = constructor.Cycle(*args, **kwargs)
        return super().set_prop_cycle(cycle)  # set the property cycler after validation

    # Rename the shorthands
    boxes = warnings._rename_objs('0.8', boxes=box)
    violins = warnings._rename_objs('0.8', violins=violin)
