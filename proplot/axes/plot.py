#!/usr/bin/env python3
"""
The plotting wrappers that add functionality to various `~matplotlib.axes.Axes`
methods. "Wrapped" `~matplotlib.axes.Axes` methods accept the additional
arguments documented in the wrapper function.
"""
import functools
import inspect
import re
import sys
from numbers import Integral, Real

import matplotlib.artist as martist
import matplotlib.axes as maxes
import matplotlib.cm as mcm
import matplotlib.collections as mcollections
import matplotlib.colors as mcolors
import matplotlib.container as mcontainer
import matplotlib.contour as mcontour
import matplotlib.legend as mlegend
import matplotlib.patches as mpatches
import matplotlib.patheffects as mpatheffects
import matplotlib.text as mtext
import matplotlib.ticker as mticker
import matplotlib.transforms as mtransforms
import numpy as np
import numpy.ma as ma

from .. import colors as pcolors
from .. import constructor
from .. import ticker as pticker
from ..config import rc
from ..internals import ic  # noqa: F401
from ..internals import (
    _dummy_context,
    _flexible_getattr,
    _not_none,
    _state_context,
    docstring,
    warnings,
)
from ..utils import edges, edges2d, to_rgb, to_xyz, units

try:
    from cartopy.crs import PlateCarree
except ModuleNotFoundError:
    PlateCarree = object

__all__ = [
    'default_latlon',
    'default_transform',
    'standardize_1d',
    'standardize_2d',
    'indicate_error',
    'apply_cmap',
    'apply_cycle',
    'colorbar_extras',
    'legend_extras',
    'text_extras',
    'bar_extras',
    'barh_extras',
    'fill_between_extras',
    'fill_betweenx_extras',
    'boxplot_extras',
    'violinplot_extras',
    'vlines_extras',
    'hlines_extras',
    'scatter_extras',
]


# Consistent keywords for styling apply_cmap-overridden plots
# TODO: Deprecate linewidth and linestyle interpretation. Think these
# already have flexible interpretation for all plotting funcs.
STYLE_ARGS_TRANSLATE = {
    'contour': {
        'colors': 'colors',
        'linewidths': 'linewidths',
        'linestyles': 'linestyles',
    },
    'tricontour': {
        'colors': 'colors',
        'linewidths': 'linewidths',
        'linestyles': 'linestyles',
    },
    'pcolor': {
        'colors': 'edgecolors',
        'linewidths': 'linewidth',
        'linestyles': 'linestyle',
    },
    'pcolormesh': {
        'colors': 'edgecolors',
        'linewidths': 'linewidth',
        'linestyles': 'linestyle',
    },
    'pcolorfast': {
        'colors': 'edgecolors',
        'linewidths': 'linewidth',
        'linestyles': 'linestyle',
    },
    'tripcolor': {
        'colors': 'edgecolors',
        'linewidths': 'linewidth',
        'linestyles': 'linestyle',
    },
    'parametric': {
        'colors': 'color',
        'linewidths': 'linewidth',
        'linestyles': 'linestyle',
    },
    'hexbin': {
        'colors': 'edgecolors',
        'linewidths': 'linewidths',
        'linestyles': 'linestyles',
    },
    'hist2d': {
        'colors': 'edgecolors',
        'linewidths': 'linewidths',
        'linestyles': 'linestyles',
    },
    'barbs': {
        'colors': 'barbcolor',
        'linewidths': 'linewidth',
        'linestyles': 'linestyle',
    },
    'quiver': {  # NOTE: linewidth/linestyle apply to *arrow outline*
        'colors': 'color',
        'linewidths': 'linewidth',
        'linestyles': 'linestyle',
    },
    'streamplot': {
        'colors': 'color',
        'linewidths': 'linewidth',
        'linestyles': 'linestyle',
    },
    'spy': {
        'colors': 'color',
        'linewidths': 'linewidth',
        'linestyles': 'linestyle',
    },
    'matshow': {
        'colors': 'color',
        'linewidths': 'linewidth',
        'linestyles': 'linestyle',
    },
    'imshow': None,
}


docstring.snippets['standardize.autoformat'] = """
autoformat : bool, optional
    Whether *x* axis labels, *y* axis labels, axis formatters, axes titles,
    legend labels, and colorbar labels are automatically configured when
    a `~pandas.Series`, `~pandas.DataFrame` or `~xarray.DataArray` is passed
    to the plotting command. Default is :rc:`autoformat`.
"""

docstring.snippets['axes.apply_cmap'] = """
cmap : colormap spec, optional
    The colormap specifer, passed to the `~proplot.constructor.Colormap`
    constructor.
cmap_kw : dict-like, optional
    Passed to `~proplot.constructor.Colormap`.
norm : normalizer spec, optional
    The colormap normalizer, used to warp data before passing it
    to `~proplot.colors.DiscreteNorm`. This is passed to the
    `~proplot.constructor.Norm` constructor.
norm_kw : dict-like, optional
    Passed to `~proplot.constructor.Norm`.
extend : {{'neither', 'min', 'max', 'both'}}, optional
    Whether to assign unique colors to out-of-bounds data and draw
    "extensions" (triangles, by default) on the colorbar.
vmin, vmax : float, optional
    Used to determine level locations if `levels` is an integer. Actual
    levels may not fall exactly on `vmin` and `vmax`, but the minimum
    level will be no smaller than `vmin` and the maximum level will be
    no larger than `vmax`. If `vmin` or `vmax` is not provided, the
    minimum and maximum data values are used.
N, levels : int or list of float, optional
    The number of level edges, or a list of level edges. If the former,
    `locator` is used to generate this many levels at "nice" intervals.
    If the latter, the levels should be monotonically increasing or
    decreasing (note that decreasing levels will only work with ``pcolor``
    plots, not ``contour`` plots). Default is :rc:`image.levels`.
    Note this means you can now discretize your colormap colors in a
    ``pcolor`` plot just like with ``contourf``.
values : int or list of float, optional
    The number of level centers, or a list of level centers. If provided,
    levels are inferred using `~proplot.utils.edges`. This will override
    any `levels` input.
symmetric : bool, optional
    If ``True``, automatically generated levels are symmetric about zero.
positive : bool, optional
    If ``True``, automatically generated levels are positive with a minimum at zero.
negative : bool, optional
    If ``True``, automatically generated levels are negative with a maximum at zero.
nozero : bool, optional
    If ``True``, ``0`` is removed from the level list.
locator : locator-spec, optional
    The locator used to determine level locations if `levels` or `values`
    is an integer and `vmin` and `vmax` were not provided. Passed to the
    `~proplot.constructor.Locator` constructor. Default is
    `~matplotlib.ticker.MaxNLocator` with ``levels`` integer levels.
locator_kw : dict-like, optional
    Passed to `~proplot.constructor.Locator`.
"""

_fill_between_docstring = """
Supports overlaying and stacking successive columns of data, and permits
using different colors for "negative" and "positive" regions.

Important
---------
This function wraps `~matplotlib.axes.Axes.fill_between{suffix}` and
`~proplot.axes.Axes.area{suffix}`.

Parameters
----------
*args : ({y}1,), ({x}, {y}1), or ({x}, {y}1, {y}2)
    The *{x}* and *{y}* coordinates. If `{x}` is not provided, it will be inferred
    from `{y}1`. If `{y}1` and `{y}2` are provided, this function will shade between
    respective columns of the arrays. The default value for `{y}2` is ``0``.
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
<https://matplotlib.org/3.1.0/gallery/pyplots/whats_new_98_4_fill_between.html#sphx-glr-gallery-pyplots-whats-new-98-4-fill-between-py>`__.
lw, linewidth : float, optional
    The edge width for the area patches.
edgecolor : color-spec, optional
    The edge color for the area patches.

Other parameters
----------------
**kwargs
    Passed to `~matplotlib.axes.Axes.fill_between`.

See also
--------
proplot.axes.Axes.area{suffix}
standardize_1d
apply_cycle
"""
docstring.snippets['axes.fill_between'] = _fill_between_docstring.format(
    x='x', y='y', suffix='',
)
docstring.snippets['axes.fill_betweenx'] = _fill_between_docstring.format(
    x='y', y='x', suffix='x',
)

_bar_docstring = """
Supports grouping and stacking successive columns of data, and changes
the default bar style.

Important
---------
This function wraps `~matplotlib.axes.Axes.bar{suffix}`.

Parameters
----------
{x}, {height}, width, {bottom} : float or list of float, optional
    The dimensions of the bars. If the *{x}* coordinates are not provided,
    they are set to ``np.arange(0, len(height))``. Note that the units
    for `width` are now *relative*.
orientation : {{None, 'vertical', 'horizontal'}}, optional
    The orientation of the bars. If ``'horizontal'``, bars are drawn horizontally
    rather than vertically. This is the default for `~matplotlib.axes.Axes.barh`.
vert : bool, optional
    Alternative to the `orientation` keyword arg. If ``False``, bars are drawn
    horizontally. Added for consistency with `~matplotlib.axes.Axes.boxplot`.
stack, stacked : bool, optional
    Whether to stack columns of the input array or plot the bars
    side-by-side in groups.
negpos : bool, optional
    Whether to shade bars greater than zero with `poscolor` and bars less
    than zero with `negcolor`.
negcolor, poscolor : color-spec, optional
    Colors to use for the negative and positive bars. Ignored if `negpos`
    is ``False``. Defaults are :rc:`negcolor` and :rc:`poscolor`.
lw, linewidth : float, optional
    The edge width for the bar patches.
edgecolor : color-spec, optional
    The edge color for the bar patches.

Other parameters
----------------
**kwargs
    Passed to `~matplotlib.axes.Axes.bar{suffix}`.

See also
--------
standardize_1d
indicate_error
apply_cycle
"""
docstring.snippets['axes.bar'] = _bar_docstring.format(
    x='x', height='height', bottom='bottom', suffix='',
)
docstring.snippets['axes.barh'] = _bar_docstring.format(
    x='y', height='right', bottom='left', suffix='h',
)

docstring.snippets['axes.lines'] = """
negpos : bool, optional
    Whether to color lines greater than zero with `poscolor` and lines less
    than zero with `negcolor`.
negcolor, poscolor : color-spec, optional
    Colors to use for the negative and positive lines. Ignored if `negpos`
    is ``False``. Defaults are :rc:`negcolor` and :rc:`poscolor`.
"""


def _load_objects():
    """
    Delay loading expensive modules. We just want to detect if *input
    arrays* belong to these types -- and if this is the case, it means the
    module has already been imported! So, we only try loading these classes
    within autoformat calls. This saves >~500ms of import time.
    """
    global DataArray, DataFrame, Series, Index, ndarray
    ndarray = np.ndarray
    DataArray = getattr(sys.modules.get('xarray', None), 'DataArray', ndarray)
    DataFrame = getattr(sys.modules.get('pandas', None), 'DataFrame', ndarray)
    Series = getattr(sys.modules.get('pandas', None), 'Series', ndarray)
    Index = getattr(sys.modules.get('pandas', None), 'Index', ndarray)


_load_objects()


def _is_number(data):
    """
    Test whether input is numeric array rather than datetime or strings.
    """
    return len(data) and np.issubdtype(_to_ndarray(data).dtype, np.number)


def _is_string(data):
    """
    Test whether input is array of strings.
    """
    return len(data) and isinstance(_to_ndarray(data).flat[0], str)


def _to_arraylike(data):
    """
    Convert list of lists to array-like type.
    """
    _load_objects()
    if not isinstance(data, (ndarray, DataArray, DataFrame, Series, Index)):
        data = np.asarray(data)
    if not np.iterable(data):
        data = np.atleast_1d(data)
    return data


def _to_indexer(data):
    """
    Return indexible attribute of array-like type.
    """
    return getattr(data, 'iloc', data)


def _to_ndarray(data):
    """
    Convert arbitrary input to ndarray cleanly. Returns a masked
    array if input is a masked array.
    """
    return np.atleast_1d(getattr(data, 'values', data))


def default_latlon(self, func, *args, latlon=True, **kwargs):
    """
    Makes ``latlon=True`` the default for basemap plots.
    This means you no longer have to pass ``latlon=True`` if your data
    coordinates are longitude and latitude.

    Important
    ---------
    This function wraps {methods} for `~proplot.axes.BasemapAxes`.
    """
    return func(self, *args, latlon=latlon, **kwargs)


def default_transform(self, func, *args, transform=None, **kwargs):
    """
    Makes ``transform=cartopy.crs.PlateCarree()`` the default
    for cartopy plots. This means you no longer have to
    pass ``transform=cartopy.crs.PlateCarree()`` if your data
    coordinates are longitude and latitude.

    Important
    ---------
    This function wraps {methods} for `~proplot.axes.CartopyAxes`.
    """
    # Apply default transform
    # TODO: Do some cartopy methods reset backgroundpatch or outlinepatch?
    # Deleted comment reported this issue
    if transform is None:
        transform = PlateCarree()
    result = func(self, *args, transform=transform, **kwargs)
    return result


def _get_label(obj):
    """
    Return a valid non-placeholder artist label from the artist or a tuple of
    artists destined for a legend. Prefer final artist (drawn last and on top).
    """
    # NOTE: BarContainer and StemContainer are instances of tuple
    while not hasattr(obj, 'get_label') and isinstance(obj, tuple) and len(obj) > 1:
        obj = obj[-1]
    label = getattr(obj, 'get_label', lambda: None)()
    return label if label and label[:1] != '_' else None


def _get_labels(data, axis=0, always=True):
    """
    Return the array-like "labels" along axis `axis` from an array-like
    object. These might be an xarray `DataArray` or pandas `Index`. If
    `always` is ``False`` we return ``None`` for simple ndarray input.
    """
    # NOTE: Previously inferred 'axis 1' metadata of 1D variable using the
    # data values metadata but that is incorrect. The paradigm for 1D plots
    # is we have row coordinates representing x, data values representing y,
    # and column coordinates representing individual series.
    if axis not in (0, 1):
        raise ValueError(f'Invalid axis {axis}.')
    _load_objects()
    labels = None
    if isinstance(data, ndarray):
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
        elif always:
            labels = np.array([0])
    # Pandas object
    elif axis == 0 and isinstance(data, (DataFrame, Series)):
        labels = data.index
    elif axis == 1 and isinstance(data, (DataFrame,)):
        labels = data.columns
    elif axis == 1 and isinstance(data, (Series, Index)) and always:
        labels = np.array([0])
    return labels


def _get_title(data, units=True):
    """
    Return the "title" associated with an array-like object with metadata. This
    might be a pandas `DataFrame` `name` or a name constructed from xarray `DataArray`
    attributes. In the latter case we search for `long_name` and `standard_name`,
    preferring the former, and append `(units)` if `units` is ``True``. If no
    names are available but units are available we just use the units string.
    """
    _load_objects()
    title = None
    if isinstance(data, ndarray):
        pass
    # Xarray object with possible long_name, standard_name, and units attributes.
    # Output depends on if units is True
    elif isinstance(data, DataArray):
        title = getattr(data, 'name', None)
        for key in ('standard_name', 'long_name'):
            title = data.attrs.get(key, title)
        if units:
            units = data.attrs.get('units', None)
            if title and units:
                title = f'{title} ({units})'
            elif units:
                title = units
    # Pandas object. Note DataFrame has no native name attribute but user can add one
    # See: https://github.com/pandas-dev/pandas/issues/447
    elif isinstance(data, (DataFrame, Series, Index)):
        title = getattr(data, 'name', None) or None
    # Standardize result
    if title is not None:
        title = str(title).strip()
    return title


def _parse_string_coords(*args, which='x', **kwargs):
    """
    Convert string arrays and lists to index coordinates.
    """
    # NOTE: Why FixedLocator and not IndexLocator? The latter requires plotting
    # lines or else error is raised... very strange.
    # NOTE: Why IndexFormatter and not FixedFormatter? The former ensures labels
    # correspond to indices while the latter can mysteriously truncate labels.
    res = []
    for arg in args:
        if _is_string(arg) and arg.ndim > 1:
            raise ValueError('Non-1D string coordinate input is unsupported.')
        if not _is_string(arg):
            res.append(arg)
            continue
        idx = np.arange(len(arg))
        kwargs.setdefault(which + 'locator', mticker.FixedLocator(idx))
        kwargs.setdefault(which + 'formatter', pticker._IndexFormatter(_to_ndarray(arg)))  # noqa: E501
        kwargs.setdefault(which + 'minorlocator', mticker.NullLocator())
        res.append(idx)
    return *res, kwargs


def _auto_format_1d(
    self, x, *ys, name='plot', vert=True, autoformat=False,
    label=None, values=None, labels=None, **kwargs
):
    """
    Try to retrieve default coordinates from array-like objects and apply default
    formatting. Also update the keyword arguments.
    """
    # Parse input
    vert = kwargs.get('vert', kwargs.get('orientation', 'vertical') == 'vertical')
    projection = hasattr(self, 'projection')
    parametric = name in ('parametric',)
    scatter = name in ('scatter',)
    hist = name in ('hist',)
    box = name in ('boxplot', 'violinplot')
    nocycle = name in ('stem', 'hlines', 'vlines', 'hexbin', 'parametric')  # WARNING
    labels = _not_none(
        label=label,
        values=values,
        labels=labels,
        colorbar_kw_values=kwargs.get('colorbar_kw', {}).pop('values', None),
        legend_kw_labels=kwargs.get('legend_kw', {}).pop('labels', None),
    )

    # The x coords or histogram labels
    # NOTE: Where columns represent distributions, like for box and violin plots or
    # where we use 'means' or 'medians', columns coords (axis 1) are 'x' coords.
    # Otherwise, columns represent e.g. lines, and row coords (axis 0) are 'x' coords.
    reduce = any(kwargs.get(s) for s in ('mean', 'means', 'median', 'medians'))
    sx = int(box or hist or reduce or False)
    if hist:
        if labels is None:
            labels = _get_labels(ys[0], axis=sx, always=False)  # hist labels
    else:
        if x is None:
            x = _get_labels(ys[0], axis=sx)  # infer from rows or columns
        x = _to_arraylike(x)

    # The labels and XY axis settings
    if not projection:
        # Apply label
        # NOTE: Do not overwrite existing labels!
        sx, sy = 'xy' if vert else 'yx'
        sy = sx if hist else sy  # histogram 'y' values end up along 'x' axis
        kw_format = {}
        if autoformat:  # 'y' axis
            title = _get_title(ys[0])
            if title and not getattr(self, f'get_{sy}label')():
                kw_format[sy + 'label'] = title
        if autoformat and not hist:  # 'x' axis
            title = _get_title(x)
            if title and not getattr(self, f'get_{sx}label')():
                kw_format[sx + 'label'] = title

        # Handle string-type coordinates
        if not hist:
            x, kw_format = _parse_string_coords(x, which=sx, **kw_format)
        if not hist and not box:
            *ys, kw_format = _parse_string_coords(*ys, which=sy, **kw_format)
        if not scatter and x.ndim == 1 and x.size > 1 and _to_ndarray(x)[1] < _to_ndarray(x)[0]:  # noqa: E501
            kw_format[sx + 'reverse'] = True  # auto reverse

        # Appply
        if kw_format:
            self.format(**kw_format)

    # Default legend or colorbar labels and title
    # WARNING: This will fail for any funcs not wrapped by apply_cycle. Keep
    # the 'nocycle' list updated until 'wrapper' functions disbanded.
    if labels is None and not nocycle and autoformat:
        labels = _get_labels(ys[0] if vert else x, axis=1)
        title = _get_title(labels)
        if not title and not any(isinstance(_, str) for _ in labels):
            labels = None
        else:
            kwargs.setdefault('legend_kw', {})
            kwargs.setdefault('colorbar_kw', {})
            kwargs['legend_kw'].setdefault('title', title)
            kwargs['colorbar_kw'].setdefault('label', title)

    # Finally strip metadata
    # WARNING: Most methods that accept 2D arrays use columns of data, but when
    # pandas DataFrame specifically is passed to hist, boxplot, or violinplot, rows
    # of data assumed! Converting to ndarray necessary.
    if labels is not None:
        if not nocycle:
            kwargs['labels'] = _to_ndarray(labels)
        elif parametric:
            kwargs['values'] = _to_ndarray(labels)
    return _to_ndarray(x), *map(_to_ndarray, ys), kwargs


@docstring.add_snippets
def standardize_1d(self, func, *args, autoformat=None, **kwargs):
    """
    Interpret positional arguments for the "1D" plotting methods so usage is
    consistent. Positional arguments are standardized as follows:

    * If a 2D array is passed, the corresponding plot command is called for
      each column of data (except for ``boxplot`` and ``violinplot``, in which
      case each column is interpreted as a distribution).
    * If *x* and *y* or *latitude* and *longitude* coordinates were not provided,
      and a `~pandas.DataFrame` or `~xarray.DataArray`, we try to infer them from
      the metadata. Otherwise, ``np.arange(0, data.shape[0])`` is used.

    Important
    ---------
    This function wraps {methods}

    Parameters
    ----------
    %(standardize.autoformat)s

    See also
    --------
    apply_cycle
    indicate_error
    """
    name = func.__name__
    fill = name in ('fill_between', 'fill_betweenx')
    box = name in ('boxplot', 'violinplot')
    autoformat = _not_none(autoformat, rc['autoformat'])
    _load_objects()
    if not args:
        return func(self, *args, **kwargs)

    # Parse input args
    # WARNING: This will temporarily add pseudo-x coordinates to hist, boxplot, and
    # pie but they are ignored when we reach apply_cycle.
    if len(args) == 1:
        x = None
        y, *args = args
    elif len(args) <= 4:  # max signature is x, y1, y2, color
        x, y, *args = args
    else:
        raise ValueError(f'Expected 1-4 positional arguments, got {len(args)}.')
    if fill and len(args) >= 1:
        *ys, args = y, args[0], args[1:]
    else:
        ys = (y,)
    ys = tuple(_to_arraylike(y) for y in ys)

    # Automatic formatting and coordinates
    x, *ys, kwargs = _auto_format_1d(
        self, x, *ys, name=name, autoformat=autoformat, **kwargs
    )

    # Ensure data is monotonic and falls within map bounds
    if getattr(self, 'name', None) == 'proplot_basemap' and kwargs.get('latlon', None):
        xmin, xmax = self.projection.lonmin, self.projection.lonmax
        x_orig, ys_orig = x, ys
        ys = []
        for y_orig in ys_orig:
            x, y = _fix_bounds(*_fix_latlon(x_orig, y_orig), xmin, xmax)
            ys.append(y)

    # Call function
    if box:
        kwargs.setdefault('positions', x)  # *this* is how 'x' is passed to boxplot
    return func(self, x, *ys, *args, **kwargs)


def _auto_format_2d(self, x, y, *Zs, order='C', autoformat=False, **kwargs):
    """
    Try to retrieve default coordinates from array-like objects and apply default
    formatting. Also apply optional transpose and update the keyword arguments.
    """
    # Retrieve coordinates
    projection = hasattr(self, 'projection')
    if x is None and y is None:
        Z = Zs[0]
        if isinstance(Z, ndarray):
            x = np.arange(Z.shape[-1])  # in case 1D
            y = np.zeros(Z.shape) if Z.ndim == 1 else np.arange(Z.shape[0])
        elif isinstance(Z, DataArray):
            x = Z.coords[Z.dims[-1]]  # in case 1D
            y = np.zeros(Z.shape) if Z.ndim == 1 else Z.coords[Z.dims[0]]
        elif isinstance(Z, DataFrame):
            x, y = Z.columns, Z.index
        elif isinstance(Z, Series):  # e.g. barbs or quiver
            x, y = Z.index, np.zeros(Z.shape)
        elif isinstance(Z, Index):  # e.g. barbs or quiver
            x, y = np.arange(Z.size), np.zeros(Z.shape)
        else:  # snould be unreachable due to _to_arraylike
            raise ValueError(f'Unrecognized array type {type(Z)}.')
        if order == 'F':
            x, y = y, x

    # Check coordinates
    x, y = _to_arraylike(x), _to_arraylike(y)
    if x.ndim != y.ndim:
        raise ValueError(
            f'x coordinates are {x.ndim}-dimensional, '
            f'but y coordinates are {y.ndim}-dimensional.'
        )
    for s, array in zip(('x', 'y'), (x, y)):
        if array.ndim not in (1, 2):
            raise ValueError(
                f'{s} coordinates are {array.ndim}-dimensional, '
                f'but must be 1 or 2-dimensional.'
            )
    if order == 'F':  # TODO: double check this
        x, y = x.T, y.T  # in case they are 2-dimensional
        Zs = tuple(Z.T for Z in Zs)

    # The labels and XY axis settings
    if not projection:
        # Apply labels
        # NOTE: Do not overwrite existing labels!
        kw_format = {}
        if autoformat:
            for s, d in zip('xy', (x, y)):
                title = _get_title(d)
                if title and not getattr(self, f'get_{s}label')():
                    kw_format[s + 'label'] = title

        # Handle string-type coordinates
        x, kw_format = _parse_string_coords(x, which='x', **kw_format)
        y, kw_format = _parse_string_coords(y, which='y', **kw_format)
        for s, d in zip('xy', (x, y)):
            if d.size > 1 and d.ndim == 1 and _to_ndarray(d)[1] < _to_ndarray(d)[0]:
                kw_format[s + 'reverse'] = True

        # Apply formatting
        if kw_format:
            self.format(**kw_format)

    # Default colorbar label
    # WARNING: This will fail for any funcs wrapped by standardize_2d but not
    # wrapped by apply_cmap. So far there are none.
    if autoformat:
        kwargs.setdefault('colorbar_kw', {})
        title = _get_title(Zs[0])
        if title and True:
            kwargs['colorbar_kw'].setdefault('label', title)

    # Finally strip metadata
    return _to_ndarray(x), _to_ndarray(y), *map(_to_ndarray, Zs), kwargs


def _enforce_centers(x, y, Z):
    """
    Enforce that coordinates are centers. Convert from edges if possible.
    """
    xlen, ylen = x.shape[-1], y.shape[0]
    if Z.ndim == 2 and Z.shape[1] == xlen - 1 and Z.shape[0] == ylen - 1:
        # Get centers given edges
        if all(z.ndim == 1 and z.size > 1 and _is_number(z) for z in (x, y)):
            x = 0.5 * (x[1:] + x[:-1])
            y = 0.5 * (y[1:] + y[:-1])
        else:
            if (
                x.ndim == 2 and x.shape[0] > 1 and x.shape[1] > 1
                and _is_number(x)
            ):
                x = 0.25 * (x[:-1, :-1] + x[:-1, 1:] + x[1:, :-1] + x[1:, 1:])
            if (
                y.ndim == 2 and y.shape[0] > 1 and y.shape[1] > 1
                and _is_number(y)
            ):
                y = 0.25 * (y[:-1, :-1] + y[:-1, 1:] + y[1:, :-1] + y[1:, 1:])
    elif Z.shape[-1] != xlen or Z.shape[0] != ylen:
        # Helpful error message
        raise ValueError(
            f'Input shapes x {x.shape} and y {y.shape} '
            f'must match Z centers {Z.shape} '
            f'or Z borders {tuple(i+1 for i in Z.shape)}.'
        )
    return x, y


def _enforce_edges(x, y, Z):
    """
    Enforce that coordinates are edges. Convert from centers if possible.
    """
    xlen, ylen = x.shape[-1], y.shape[0]
    if Z.ndim == 2 and Z.shape[1] == xlen and Z.shape[0] == ylen:
        # Get edges given centers
        if all(z.ndim == 1 and z.size > 1 and _is_number(z) for z in (x, y)):
            x = edges(x)
            y = edges(y)
        else:
            if (
                x.ndim == 2 and x.shape[0] > 1 and x.shape[1] > 1
                and _is_number(x)
            ):
                x = edges2d(x)
            if (
                y.ndim == 2 and y.shape[0] > 1 and y.shape[1] > 1
                and _is_number(y)
            ):
                y = edges2d(y)
    elif Z.shape[-1] != xlen - 1 or Z.shape[0] != ylen - 1:
        # Helpful error message
        raise ValueError(
            f'Input shapes x {x.shape} and y {y.shape} must match '
            f'array centers {Z.shape} or '
            f'array borders {tuple(i + 1 for i in Z.shape)}.'
        )
    return x, y


def _fix_bounds(x, y, xmin, xmax):
    """
    Ensure data for basemap plots is restricted between the minimum and
    maximum longitude of the projection. Input is the ``x`` and ``y``
    coordinates. The ``y`` coordinates are rolled along the rightmost axis.
    """
    if x.ndim != 1:
        return x, y

    # Roll in same direction if some points on right-edge extend
    # more than 360 above min longitude; *they* should be on left side
    lonroll = np.where(x > xmin + 360)[0]  # tuple of ids
    if lonroll.size:  # non-empty
        roll = x.size - lonroll.min()
        x = np.roll(x, roll)
        y = np.roll(y, roll, axis=-1)
        x[:roll] -= 360  # make monotonic

    # Set NaN where data not in range xmin, xmax. Must be done
    # for regional smaller projections or get weird side-effects due
    # to having valid data way outside of the map boundaries
    y = y.copy()
    if x.size - 1 == y.shape[-1]:  # test western/eastern grid cell edges
        y[..., (x[1:] < xmin) | (x[:-1] > xmax)] = np.nan
    elif x.size == y.shape[-1]:  # test the centers and pad by one for safety
        where = np.where((x < xmin) | (x > xmax))[0]
        y[..., where[1:-1]] = np.nan

    return x, y


def _fix_latlon(x, y):
    """
    Ensure longitudes are monotonic and make `~numpy.ndarray` copies so the
    contents can be modified. Ignores 2D coordinate arrays.
    """
    if x.ndim != 1 or all(x < x[0]):  # skip 2D arrays and monotonic backwards data
        return x, y
    lon1 = x[0]
    while True:
        filter_ = (x < lon1)
        if filter_.sum() == 0:
            break
        x[filter_] += 360
    return x, y


def _fix_poles(y, Z):
    """
    Add data points on the poles as the average of highest latitude data.
    """
    # Get means
    with np.errstate(all='ignore'):
        p1 = Z[0, :].mean()  # pole 1, make sure is not 0D DataArray!
        p2 = Z[-1, :].mean()  # pole 2
    if hasattr(p1, 'item'):
        p1 = np.asscalar(p1)  # happens with DataArrays
    if hasattr(p2, 'item'):
        p2 = np.asscalar(p2)
    # Concatenate
    ps = (-90, 90) if (y[0] < y[-1]) else (90, -90)
    Z1 = np.repeat(p1, Z.shape[1])[None, :]
    Z2 = np.repeat(p2, Z.shape[1])[None, :]
    y = ma.concatenate((ps[:1], y, ps[1:]))
    Z = ma.concatenate((Z1, Z, Z2), axis=0)
    return y, Z


def _fix_cartopy(x, y, *Zs, globe=False):
    """
    Fix cartopy geographic data arrays.
    """
    # Fix coordinates
    x, y = _fix_latlon(x, y)

    # Fix data
    x_orig, y_orig, Zs_orig = x, y, Zs
    Zs = []
    for Z_orig in Zs_orig:
        # Bail for 2D coordinates
        if not globe or x_orig.ndim > 1 or y_orig.ndim > 1:
            Zs.append(Z_orig)
            continue
        # Fix holes over poles by *interpolating* there
        y, Z = _fix_poles(y_orig, Z_orig)
        # Fix seams by ensuring circular coverage (cartopy can plot over map edges)
        if x_orig[0] % 360 != (x_orig[-1] + 360) % 360:
            x = ma.concatenate((x_orig, [x_orig[0] + 360]))
            Z = ma.concatenate((Z, Z[:, :1]), axis=1)
        Zs.append(Z)

    return x, y, *Zs


def _fix_basemap(x, y, *Zs, globe=False, projection=None):
    """
    Fix basemap geographic data arrays.
    """
    # Fix coordinates
    x, y = _fix_latlon(x, y)

    # Fix data
    xmin, xmax = projection.lonmin, projection.lonmax
    x_orig, y_orig, Zs_orig = x, y, Zs
    Zs = []
    for Z_orig in Zs_orig:
        # Ensure data is within map bounds
        x, Z_orig = _fix_bounds(x_orig, Z_orig, xmin, xmax)
        # Bail for 2D coordinates
        if not globe or x_orig.ndim > 1 or y_orig.ndim > 1:
            Zs.append(Z_orig)
            continue
        # Fix holes over poles by *interpolating* there
        y, Z = _fix_poles(y_orig, Z_orig)
        # Fix seams at map boundary
        if x[0] == xmin and x.size - 1 == Z.shape[1]:  # scenario 1
            # Edges (e.g. pcolor) fit perfectly against seams. Size is unchanged.
            pass
        elif x.size - 1 == Z.shape[1]:  # scenario 2
            # Edges (e.g. pcolor) do not fit perfectly. Size augmented by 1.
            x = ma.append(xmin, x)
            x[-1] = xmin + 360
            Z = ma.concatenate((Z[:, -1:], Z), axis=1)
        elif x.size == Z.shape[1]:  # scenario 3
            # Centers (e.g. contour) must be interpolated to edge. Size augmented by 2.
            xi = np.array([x[-1], x[0] + 360])
            if xi[0] == xi[1]:  # impossible to interpolate
                pass
            else:
                Zq = ma.concatenate((Z[:, -1:], Z[:, :1]), axis=1)
                xq = xmin + 360
                Zq = (Zq[:, :1] * (xi[1] - xq) + Zq[:, 1:] * (xq - xi[0])) / (xi[1] - xi[0])  # noqa: E501
                x = ma.concatenate(([xmin], x, [xmin + 360]))
                Z = ma.concatenate((Zq, Z, Zq), axis=1)
        else:
            raise ValueError('Unexpected shapes of coordinates or data arrays.')
        Zs.append(Z)

    # Convert coordinates
    if x.ndim == 1 and y.ndim == 1:
        x, y = np.meshgrid(x, y)
    x, y = projection(x, y)

    return x, y, *Zs


@docstring.add_snippets
def standardize_2d(
    self, func, *args, autoformat=None, order='C', globe=False, **kwargs
):
    """
    Interpret positional arguments for the "2D" plotting methods so usage is
    consistent. Positional arguments are standardized as follows:

    * If *x* and *y* or *latitude* and *longitude* coordinates were not
      provided, and a `~pandas.DataFrame` or `~xarray.DataArray` is passed, we
      try to infer them from the metadata. Otherwise, ``np.arange(0, data.shape[0])``
      and ``np.arange(0, data.shape[1])`` are used.
    * For ``pcolor`` and ``pcolormesh``, coordinate *edges* are calculated
      if *centers* were provided. For all other methods, coordinate *centers*
      are calculated if *edges* were provided.

    Important
    ---------
    This function wraps {methods}

    Parameters
    ----------
    %(standardize.autoformat)s
    order : {{'C', 'F'}}, optional
        If ``'C'``, arrays should be shaped ``(y, x)``. If ``'F'``, arrays
        should be shaped ``(x, y)``. Default is ``'C'``.
    globe : bool, optional
        Whether to ensure global coverage for `~proplot.axes.GeoAxes` plots.
        Default is ``False``. When set to ``True`` this does the following:

        #. Interpolates input data to the North and South poles by setting the data
           values at the poles to the mean from latitudes nearest each pole.
        #. Makes meridional coverage "circular", i.e. the last longitude coordinate
           equals the first longitude coordinate plus 360\N{DEGREE SIGN}.
        #. For `~proplot.axes.BasemapAxes`, 1D longitude vectors are also cycled to
           fit within the map edges. For example, if the projection central longitude
           is 90\N{DEGREE SIGN}, the data is shifted so that it spans
           -90\N{DEGREE SIGN} to 270\N{DEGREE SIGN}.

    See also
    --------
    apply_cmap
    """
    name = func.__name__
    pcolor = name in ('pcolor', 'pcolormesh', 'pcolorfast')
    allow_1d = name in ('barbs', 'quiver')  # these also allow 1D data
    autoformat = _not_none(autoformat, rc['autoformat'])
    _load_objects()
    if not args:
        return func(self, *args, **kwargs)

    # Parse input args
    if len(args) > 5:
        raise TypeError(f'Expected 1-5 positional arguments, got {len(args)}.')
    if len(args) > 2:
        x, y, *args = args
    else:
        x = y = None
    if order not in ('C', 'F'):
        raise ValueError(f"Expected order 'C' (row-major) or 'F' (column-major). Got {order!r}.")  # noqa: E501
    Zs = []
    for Z in args:
        Z = _to_arraylike(Z)
        if Z.ndim != 2 and not (allow_1d and Z.ndim == 1):
            raise ValueError(f'Array must be 2-dimensional, but got shape {Z.shape}.')
        Zs.append(Z)
    shapes = tuple(Z.shape for Z in Zs)
    if len(set(shapes)) > 1:
        raise ValueError(f'Arrays must have same shape, but got shapes {shapes}.')

    # Automatic formatting
    x, y, *Zs, kwargs = _auto_format_2d(
        self, x, y, *Zs, autoformat=autoformat, **kwargs
    )

    # Standardize coordinates
    if pcolor:
        x, y = _enforce_edges(x, y, Zs[0])
    else:
        x, y = _enforce_centers(x, y, Zs[0])

    # Cartopy projection axes
    if (
        not allow_1d and getattr(self, 'name', None) == 'proplot_cartopy'
        and isinstance(kwargs.get('transform', None), PlateCarree)
    ):
        x, y, *Zs = _fix_cartopy(x, y, *Zs, globe=globe)

    # Basemap projection axes
    elif (
        not allow_1d and getattr(self, 'name', None) == 'proplot_basemap'
        and kwargs.get('latlon', None)
    ):
        x, y, *Zs = _fix_basemap(x, y, *Zs, globe=globe, projection=self.projection)
        kwargs['latlon'] = False

    # Call function
    return func(self, x, y, *Zs, **kwargs)


def _deprecate_add_errorbars(func):
    """
    Translate old-style keyword arguments to new-style in way that is too complex
    for _rename_kwargs. Use a decorator to avoid call signature pollution.
    """
    @functools.wraps(func)
    def wrapper(
        *args, bars=None, boxes=None,
        barstd=None, boxstd=None, barrange=None, boxrange=None, **kwargs
    ):
        for (prefix, b, std, range_) in zip(
            ('bar', 'box'), (bars, boxes), (barstd, boxstd), (barrange, boxrange),
        ):
            if b is not None or std is not None or range_ is not None:
                warnings._warn_proplot(
                    f"Keyword args '{prefix}s', '{prefix}std', and '{prefix}range' "
                    'are deprecated and will be removed in a future release. '
                    f"Please use '{prefix}stds' or '{prefix}pctiles' instead."
                )
            if range_ is None and b:  # means 'use the default range'
                range_ = b
            if std:
                kwargs.setdefault(prefix + 'stds', range_)
            else:
                kwargs.setdefault(prefix + 'pctiles', range_)
        return func(*args, **kwargs)
    return wrapper


def _get_error_data(
    data, y, errdata=None, stds=None, pctiles=None,
    stds_default=None, pctiles_default=None,
    reduced=True, absolute=False, label=False,
):
    """
    Return values that can be passed to the `~matplotlib.axes.Axes.errorbar`
    `xerr` and `yerr` keyword args.
    """
    # Parse arguments
    # NOTE: Have to guard against "truth value of an array is ambiguous" errors
    if stds is True:
        stds = stds_default
    elif stds is False or stds is None:
        stds = None
    elif isinstance(stds, Real):
        stds = sorted((-stds, stds))
    elif not np.iterable(stds) or len(stds) != 2:
        raise ValueError('Expected scalar or length-2 tuple stdev specification.')
    if pctiles is True:
        pctiles = pctiles_default
    elif pctiles is False or pctiles is None:
        pctiles = None
    elif isinstance(pctiles, Real):
        delta = (100 - pctiles) / 2.0
        pctiles = sorted((delta, 100 - delta))
    elif not np.iterable(pctiles) or len(pctiles) != 2:
        raise ValueError('Expected length-2 tuple pctiles specification.')

    # Incompatible settings
    if stds is not None and pctiles is not None:
        warnings._warn_proplot(
            'You passed both a standard deviation range and a percentile range for '
            'drawing error indicators. Using the former.'
        )
        pctiles = None
    if not reduced and (stds is not None or pctiles is not None):
        raise ValueError(
            'To automatically compute standard deviations or percentiles on columns '
            'of data you must pass means=True or medians=True.'
        )
    if reduced and errdata is not None:
        stds = pctiles = None
        warnings._warn_proplot(
            'You explicitly provided the error bounds but also requested '
            'automatically calculating means or medians on data columns. '
            'It may make more sense to use the "stds" or "pctiles" keyword args '
            'and have *proplot* calculate the error bounds.'
        )

    # Compute error data in format that can be passed to matplotlib.axes.Axes.errorbar()
    # NOTE: Include option to pass symmetric deviation from central points
    if errdata is not None:
        label_default = 'uncertainty'
        err = _to_ndarray(errdata)
        if (
            err.ndim not in (1, 2)
            or err.shape[-1] != y.shape[-1]
            or err.ndim == 2 and err.shape[0] != 2
        ):
            raise ValueError(
                f'errdata must have shape (2, {y.shape[-1]}), but got {err.shape}.'
            )
        if err.ndim == 1:
            abserr = err
            err = np.empty((2, err.size))
            err[0, :] = y - abserr  # translated back to absolute deviations below
            err[1, :] = y + abserr
    elif stds is not None:
        label_default = fr'{abs(stds[1])}$\sigma$ range'
        err = y + np.std(data, axis=0)[None, :] * _to_ndarray(stds)[:, None]
    elif pctiles is not None:
        label_default = f'{pctiles[1] - pctiles[0]}% range'
        err = np.percentile(data, pctiles, axis=0)
    else:
        raise ValueError('You must provide error bounds.')
    if label is True:
        label = label_default
    elif not label:
        label = None
    if not absolute:
        err = err - y
        err[0, :] *= -1  # absolute deviations from central points

    # Return data with legend entry
    return err, label


@_deprecate_add_errorbars
def indicate_error(
    self, func, *args,
    mean=None, means=None, median=None, medians=None,
    barstd=None, barstds=None, barpctile=None, barpctiles=None, bardata=None,
    boxstd=None, boxstds=None, boxpctile=None, boxpctiles=None, boxdata=None,
    shadestd=None, shadestds=None, shadepctile=None, shadepctiles=None, shadedata=None,
    fadestd=None, fadestds=None, fadepctile=None, fadepctiles=None, fadedata=None,
    boxmarker=None, boxmarkercolor='white',
    boxcolor=None, barcolor=None, shadecolor=None, fadecolor=None,
    shadelabel=False, fadelabel=False, shadealpha=0.4, fadealpha=0.2,
    boxlinewidth=None, boxlw=None, barlinewidth=None, barlw=None, capsize=None,
    boxzorder=2.5, barzorder=2.5, shadezorder=1.5, fadezorder=1.5,
    **kwargs
):
    """
    Adds support for drawing error bars and error shading on-the-fly. Includes
    options for interpreting columns of data as *samples*, representing the mean
    or median of each sample with lines, points, or bars, and drawing error bars
    representing percentile ranges or standard deviation multiples for each sample.
    Also supports specifying error bar data explicitly.

    Important
    ---------
    This function wraps {methods}

    Parameters
    ----------
    *args
        The input data.
    mean, means : bool, optional
        Whether to plot the means of each column in the input data. If no other
        arguments specified, this also sets ``barstd=True`` (and ``boxstd=True``
        for violin plots).
    median, medians : bool, optional
        Whether to plot the medians of each column in the input data. If no other
        arguments specified, this also sets ``barstd=True`` (and ``boxstd=True``
        for violin plots).
    vert : bool, optional
        If ``False``, error data is drawn horizontally rather than vertially. Set
        automatically by methods like `bar`, `barh`, `area`, and `areax`.
    orientation : {{None, 'vertical', 'horizontal'}}, optional
        Alternative to the `vert` keyword arg. If ``'horizontal'``, error data is
        drawn horizontally rather than vertically.
    barstd, barstds : float, (float, float), or bool, optional
        Standard deviation multiples for *thin error bars* with optional whiskers
        (i.e. caps). If scalar, then +/- that number is used. If ``True``, the
        default of +/-3 standard deviations is used. This argument is only valid
        if `means` or `medians` is ``True``.
    barpctile, barpctiles : float, (float, float) or bool, optional
        As with `barstd`, but instead using *percentiles* for the error bars. The
        percentiles are calculated with `numpy.percentile`. If scalar, that width
        surrounding the 50th percentile is used (e.g. ``90`` shows the 5th to 95th
        percentiles). If ``True``, the default percentile range of 0 to 100 is
        used. This argument is only valid if `means` or `medians` is ``True``.
    bardata : 2 x N array or 1D array, optional
        If shape is 2 x N these are the lower and upper bounds for the thin error bars.
        If array is 1D these are the absolute, symmetric deviations from the central
        points.  This should be used if `means` and `medians` are both ``False`` (i.e.
        you did not provide dataset columns from which statistical properties can be
        calculated automatically).
    boxstd, boxstds, boxpctile, boxpctiles, boxdata : optional
        As with `barstd`, `barpctile`, and `bardata`, but for *thicker error bars*
        representing a smaller interval than the thin error bars. If `boxstds` is
        ``True``, the default standard deviation range of +/-1 is used. If `boxpctiles`
        is ``True``, the default percentile range of 25 to 75 is used (i.e. the
        interquartile range). When "boxes" and "bars" are combined, this has the effect
        of drawing miniature box-and-whisker plots.
    shadestd, shadestds, shadepctile, shadepctiles, shadedata : optional
        As with `barstd`, `barpctile`, and `bardata`, but using *shading* to indicate
        the error range. If `shadestds` is ``True``, the default standard deviation
        range of +/-2 is used. If `shadepctiles` is ``True``, the default
        percentile range of 10 to 90 is used. Shading is generally useful for
        `~matplotlib.axes.Axes.plot` plots.
    fadestd, fadestds, fadepctile, fadepctiles, fadedata : optional
        As with `shadestd`, `shadepctile`, and `shadedata`, but for an additional,
        more faded, *secondary* shaded region. If `fadestds` is ``True``, the default
        standard deviation range of +/-3 is used. If `fadepctiles` is ``True``,
        the default percentile range of 0 to 100 is used.
    barcolor, boxcolor, shadecolor, fadecolor : color-spec, optional
        Colors for the different error indicators. For error bars, the default is
        ``'k'``. For shading, the default behavior is to inherit color from the
        primary `~matplotlib.artist.Artist`.
    shadelabel, fadelabel : bool or str, optional
        Labels for the shaded regions to be used as separate legend entries. To toggle
        labels "on" and apply a *default* label, use e.g. ``shadelabel=True``. To apply
        a *custom* label, use e.g. ``shadelabel='label'``. Otherwise, the shading is
        drawn underneath the line and/or marker in the legend entry.
    barlinewidth, boxlinewidth, barlw, boxlw : float, optional
        Line widths for the thin and thick error bars, in points. The defaults
        are ``barlw=0.8`` and ``boxlw=4 * barlw``.
    boxmarker : bool, optional
        Whether to draw a small marker in the middle of the box denoting the mean or
        median position. Ignored if `boxes` is ``False``.
    boxmarkercolor : color-spec, optional
        Color for the `boxmarker` marker. Default is ``'w'``.
    capsize : float, optional
        The cap size for thin error bars in points.
    barzorder, boxzorder, shadezorder, fadezorder : float, optional
        The "zorder" for the thin error bars, thick error bars, and shading.

    Returns
    -------
    h, err1, err2, ...
        The original plot object and the error bar or shading objects.
    """
    name = func.__name__
    bar = name in ('bar',)
    violin = name in ('violinplot',)
    plot = name in ('plot', 'scatter')
    means = _not_none(mean=mean, means=means)
    medians = _not_none(median=median, medians=medians)
    barstds = _not_none(barstd=barstd, barstds=barstds)
    boxstds = _not_none(boxstd=boxstd, boxstds=boxstds)
    shadestds = _not_none(shadestd=shadestd, shadestds=shadestds)
    fadestds = _not_none(fadestd=fadestd, fadestds=fadestds)
    barpctiles = _not_none(barpctile=barpctile, barpctiles=barpctiles)
    boxpctiles = _not_none(boxpctile=boxpctile, boxpctiles=boxpctiles)
    shadepctiles = _not_none(shadepctile=shadepctile, shadepctiles=shadepctiles)
    fadepctiles = _not_none(fadepctile=fadepctile, fadepctiles=fadepctiles)
    bars = any(_ is not None for _ in (bardata, barstds, barpctiles))
    boxes = any(_ is not None for _ in (boxdata, boxstds, boxpctiles))
    shade = any(_ is not None for _ in (shadedata, shadestds, shadepctiles))
    fade = any(_ is not None for _ in (fadedata, fadestds, fadepctiles))
    if means and medians:
        warnings._warn_proplot('Cannot have both means=True and medians=True. Using former.')  # noqa: E501

    # Get means or medians while preserving metadata for autoformat
    # TODO: Permit 3D array with error dimension coming first
    # NOTE: Previously went to great pains to preserve metadata but now retrieval
    # of default legend handles moved to _auto_format_1d so can strip.
    x, data, *args = args
    y = data
    if means or medians:
        if data.ndim != 2:
            raise ValueError(f'Expected 2D array for means=True. Got {data.ndim}D.')
        if not any((bars, boxes, shade, fade)):
            bars = barstds = True
            if violin:
                boxes = boxstds = True
        if means:
            y = np.mean(data, axis=0)
        elif medians:
            y = np.percentile(data, 50, axis=0)

    # Parse keyword args and apply defaults
    # NOTE: Should not use plot() 'linewidth' for bar elements
    # NOTE: violinplot_extras passes some invalid keyword args with expectation
    # that indicate_error pops them and uses them for error bars.
    method = kwargs.pop if violin else kwargs.get if bar else lambda *args: None
    boxmarker = _not_none(boxmarker, True if violin else False)
    capsize = _not_none(capsize, 3.0)
    linewidth = _not_none(method('linewidth', None), method('lw', None), 1.0)
    barlinewidth = _not_none(barlinewidth=barlinewidth, barlw=barlw, default=linewidth)
    boxlinewidth = _not_none(boxlinewidth=boxlinewidth, boxlw=boxlw, default=4 * barlinewidth)  # noqa: E501
    edgecolor = _not_none(method('edgecolor', None), 'k')
    barcolor = _not_none(barcolor, edgecolor)
    boxcolor = _not_none(boxcolor, barcolor)
    shadecolor_infer = shadecolor is None
    fadecolor_infer = fadecolor is None
    shadecolor = _not_none(shadecolor, kwargs.get('color'), kwargs.get('facecolor'), edgecolor)  # noqa: E501
    fadecolor = _not_none(fadecolor, shadecolor)

    # Draw dark and light shading
    method = kwargs.pop if plot else kwargs.get
    vert = method('vert', method('orientation', 'vertical') == 'vertical')
    eobjs = []
    method = self.fill_between if vert else self.fill_betweenx
    if fade:
        edata, label = _get_error_data(
            data, y, errdata=fadedata, stds=fadestds, pctiles=fadepctiles,
            stds_default=(-3, 3), pctiles_default=(0, 100), absolute=True,
            reduced=means or medians, label=fadelabel,
        )
        eobj = method(
            x, *edata, linewidth=0, label=label,
            color=fadecolor, alpha=fadealpha, zorder=fadezorder,
        )
        eobjs.append(eobj)
    if shade:
        edata, label = _get_error_data(
            data, y, errdata=shadedata, stds=shadestds, pctiles=shadepctiles,
            stds_default=(-2, 2), pctiles_default=(10, 90), absolute=True,
            reduced=means or medians, label=shadelabel,
        )
        eobj = method(
            x, *edata, linewidth=0, label=label,
            color=shadecolor, alpha=shadealpha, zorder=shadezorder,
        )
        eobjs.append(eobj)

    # Draw thin error bars and thick error boxes
    sy = 'y' if vert else 'x'  # yerr
    ex, ey = (x, y) if vert else (y, x)
    if boxes:
        edata, _ = _get_error_data(
            data, y, errdata=boxdata, stds=boxstds, pctiles=boxpctiles,
            stds_default=(-1, 1), pctiles_default=(25, 75),
            reduced=means or medians,
        )
        if boxmarker:
            self.scatter(
                ex, ey, s=boxlinewidth, marker='o', color=boxmarkercolor, zorder=5
            )
        eobj = self.errorbar(
            ex, ey, color=boxcolor, linewidth=boxlinewidth, linestyle='none',
            capsize=0, zorder=boxzorder, **{sy + 'err': edata}
        )
        eobjs.append(eobj)
    if bars:  # now impossible to make thin bar width different from cap width!
        edata, _ = _get_error_data(
            data, y, errdata=bardata, stds=barstds, pctiles=barpctiles,
            stds_default=(-3, 3), pctiles_default=(0, 100),
            reduced=means or medians,
        )
        eobj = self.errorbar(
            ex, ey, color=barcolor, linewidth=barlinewidth, linestyle='none',
            markeredgecolor=barcolor, markeredgewidth=barlinewidth,
            capsize=capsize, zorder=barzorder, **{sy + 'err': edata}
        )
        eobjs.append(eobj)

    # Call main function
    # NOTE: Provide error objects for inclusion in legend, but *only* provide
    # the shading. Never want legend entries for error bars.
    xy = (x, data) if name == 'violinplot' else (x, y)
    kwargs.setdefault('_errobjs', eobjs[:int(shade + fade)])
    res = obj = func(self, *xy, *args, **kwargs)

    # Apply inferrred colors to objects
    i = 0
    if isinstance(res, (list, tuple)):  # avoid BarContainer
        obj = res[0]
    for b, infer in zip((fade, shade), (fadecolor_infer, shadecolor_infer)):
        if not b or not infer:
            continue
        if hasattr(obj, 'get_facecolor'):
            color = obj.get_facecolor()
        elif hasattr(obj, 'get_color'):
            color = obj.get_color()
        else:
            color = None
        if color is not None:
            eobjs[i].set_facecolor(color)
        i += 1

    # Return objects
    # NOTE: This should not affect internal matplotlib calls to these funcs
    # NOTE: Avoid expanding matplolib collections that are list subclasses here
    if eobjs:
        return (*res, *eobjs) if isinstance(res, (list, tuple)) else (res, *eobjs)
    else:
        return res


def _plot_extras(self, func, *args, cmap=None, values=None, **kwargs):
    """
    Adds the option `orientation` to change the default orientation of the
    lines. See also `~proplot.axes.Axes.plotx`.
    """
    if len(args) > 3:  # e.g. with fmt string
        raise ValueError(f'Expected 1-3 positional args, got {len(args)}.')
    if cmap is not None:
        warnings._warn_proplot(
            'Drawing "parametric" plots with ax.plot(x, y, values=values, cmap=cmap) '
            'is deprecated and will be removed in a future release. Please use '
            'ax.parametric(x, y, values, cmap=cmap) instead.'
        )
        return self.parametric(*args, cmap=cmap, values=values, **kwargs)

    # Call function
    result = func(self, *args, values=values, **kwargs)

    # Add sticky edges? No because there is no way to check whether "dependent variable"
    # is x or y axis like with area/areax and bar/barh. Better to always have margin.
    # for obj in result:
    #     xdata = obj.get_xdata()
    #     obj.sticky_edges.x.append(self.convert_xunits(min(xdata)))
    #     obj.sticky_edges.x.append(self.convert_yunits(max(xdata)))

    return result


def _parametric_extras(self, func, *args, interp=0, **kwargs):
    """
    Calls `~proplot.axes.Axes.parametric` and optionally interpolates values before
    they get passed to `apply_cmap` and the colormap boundaries are drawn. Full
    documentation is on public axes method.
    """
    # Parse input arguments
    # WARNING: This is separated from parametric because its task is analogous to
    # the other wrappers -- standardizing arguments before ingestion by other
    # functions. So far this only works for 1D *x* and *y* coordinates.
    if len(args) == 3:
        x, y, values = args
    elif 'values' in kwargs:
        values = kwargs.pop('values')
        if len(args) == 1:
            y = args[0]
            x = np.arange(y.shape[-1])
        elif len(args) == 2:
            x, y = args
        else:
            raise ValueError(f'Expected 1-3 positional arguments, got {len(args)}.')
    else:
        raise ValueError('Missing required keyword argument "values".')
    x = _to_arraylike(x)
    y = _to_arraylike(y)
    values = _to_arraylike(values)
    if any(_.ndim != 1 or _.size != x.size for _ in (x, y, values)):
        raise ValueError(
            f'x {x.shape}, y {y.shape}, and values {values.shape} '
            'must be 1-dimensional and have the same size.'
        )

    # Interpolate values to allow for smooth gradations between values
    # (interp=False) or color switchover halfway between points
    # (interp=True). Then optionally interpolate the colormap values.
    if interp > 0:
        x_orig, y_orig, v_orig = x, y, values
        x, y, values = [], [], []
        for j in range(x_orig.shape[0] - 1):
            idx = slice(None)
            if j + 1 < x_orig.shape[0] - 1:
                idx = slice(None, -1)
            x.extend(np.linspace(x_orig[j], x_orig[j + 1], interp + 2)[idx].flat)
            y.extend(np.linspace(y_orig[j], y_orig[j + 1], interp + 2)[idx].flat)
            values.extend(np.linspace(v_orig[j], v_orig[j + 1], interp + 2)[idx].flat)
        x, y, values = np.array(x), np.array(y), np.array(values)

    # Call main function
    return func(self, x, y, values=values, **kwargs)


def _stem_extras(
    self, func, *args, linefmt=None, basefmt=None, markerfmt=None, **kwargs
):
    """
    Make `use_line_collection` the default to suppress warning message.
    """
    # Set default colors
    # NOTE: 'fmt' strings can only be 2 to 3 characters and include color shorthands
    # like 'r' or cycle colors like 'C0'. Cannot use full color names.
    # NOTE: Matplotlib defaults try to make a 'reddish' color the base and 'bluish'
    # color the stems. To make this more robust we temporarily replace the cycler
    # with a negcolor/poscolor cycler, otherwise try to point default colors to the
    # blush 'C0' and reddish 'C1' from the new default 'colorblind' cycler.
    fmts = (linefmt, basefmt, markerfmt)
    if not any(isinstance(fmt, str) and re.match(r'\AC[0-9]', fmt) for fmt in fmts):
        cycle = constructor.Cycle((rc['negcolor'], rc['poscolor']), name='_neg_pos')
        context = rc.context({'axes.prop_cycle': cycle})
    else:
        context = _dummy_context()

    # Add stem lines with bluish stem color and reddish base color
    with context:
        kwargs['linefmt'] = linefmt = _not_none(linefmt, 'C0-')
        kwargs['basefmt'] = _not_none(basefmt, 'C1-')
        kwargs['markerfmt'] = _not_none(markerfmt, linefmt[:-1] + 'o')
        kwargs.setdefault('use_line_collection', True)
        try:
            return func(self, *args, **kwargs)
        except TypeError:
            kwargs.pop('use_line_collection')  # old version
            return func(self, *args, **kwargs)


def _lines_apply(
    self, func, *args, negpos=False, negcolor=None, poscolor=None, **kwargs
):
    """
    Parse lines arguments. Support automatic *x* coordinates and default
    "minima" at zero.
    """
    # Parse positional arguments, use default "base" position of zero
    x = 'x' if func.__name__ == 'vlines' else 'y'
    y = 'y' if x == 'x' else 'x'
    args = list(args)
    if x in kwargs:
        args.insert(0, kwargs.pop(x))
    for suffix in ('min', 'max'):
        key = y + suffix
        if key in kwargs:
            args.append(kwargs.pop(key))
    if len(args) == 1:
        x = np.arange(len(np.atleast_1d(args[0])))
        args.insert(0, x)
    if len(args) == 2:
        args.insert(1, 0.0)
    elif len(args) != 3:
        raise TypeError(f'Expected 1-3 positional arguments, got {len(args)}.')

    # Support "negative" and "positive" lines
    x, y1, y2 = args
    if negpos and kwargs.get('color', None) is None:
        y1 = _to_arraylike(y1)
        y2 = _to_arraylike(y2)
        y1array = _to_ndarray(y1)
        y2array = _to_ndarray(y2)

        # Negative colors
        mask = y2array >= y1array  # positive
        y1neg = y1.copy()
        y2neg = y2.copy()
        if mask.size == 1:
            if mask.item():
                y1neg = y2neg = np.nan
        else:
            if y1.size > 1:
                _to_indexer(y1neg)[mask] = np.nan
            if y2.size > 1:
                _to_indexer(y2neg)[mask] = np.nan
        color = _not_none(negcolor, rc['negcolor'])
        negobj = func(self, x, y1neg, y2neg, color=color, **kwargs)

        # Positive colors
        mask = y2array < y1array  # negative
        y1pos = y1.copy()
        y2pos = y2.copy()
        if mask.size == 1:
            if mask.item():
                y1pos = y2pos = np.nan
        else:
            if y1.size > 1:
                _to_indexer(y1pos)[mask] = np.nan
            if y2.size > 1:
                _to_indexer(y2pos)[mask] = np.nan
        color = _not_none(poscolor, rc['poscolor'])
        posobj = func(self, x, y1pos, y2pos, color=color, **kwargs)

        # Return both objects
        return (negobj, posobj)
    else:
        return func(self, x, y1, y2, **kwargs)


@docstring.add_snippets
def hlines_extras(self, func, *args, **kwargs):
    """
    Plot horizontal lines with flexible positional arguments and optionally
    use different colors for "negative" and "positive" lines.

    Important
    ---------
    This function wraps {methods}

    Parameters
    ----------
    %(axes.lines)s

    See also
    --------
    standardize_1d
    """
    return _lines_apply(self, func, *args, **kwargs)


@docstring.add_snippets
def vlines_extras(self, func, *args, **kwargs):
    """
    Plot vertical lines with flexible positional arguments and optionally
    use different colors for "negative" and "positive" lines.

    Important
    ---------
    This function wraps {methods}

    Parameters
    ----------
    %(axes.lines)s

    See also
    --------
    standardize_1d
    """
    return _lines_apply(self, func, *args, **kwargs)


@docstring.add_snippets
def scatter_extras(
    self, func, *args,
    s=None, size=None, markersize=None,
    c=None, color=None, markercolor=None, smin=None, smax=None,
    cmap=None, cmap_kw=None, norm=None, norm_kw=None,
    vmin=None, vmax=None, extend='neither', N=None, levels=None, values=None,
    symmetric=False, locator=None, locator_kw=None,
    lw=None, linewidth=None, linewidths=None,
    markeredgewidth=None, markeredgewidths=None,
    edgecolor=None, edgecolors=None,
    markeredgecolor=None, markeredgecolors=None,
    **kwargs
):
    """
    Adds keyword arguments to `~matplotlib.axes.Axes.scatter` that are more
    consistent with the `~matplotlib.axes.Axes.plot` keyword arguments and
    supports `apply_cmap` features.

    Important
    ---------
    This function wraps {methods}

    Parameters
    ----------
    s, size, markersize : float or list of float, optional
        The marker size(s). The units are optionally scaled by
        `smin` and `smax`.
    smin, smax : float, optional
        The minimum and maximum marker size in units ``points^2`` used to scale
        `s`. If not provided, the marker sizes are equivalent to the values in `s`.
    c, color, markercolor : color-spec or list thereof, or array, optional
        The marker fill color(s). If this is an array of scalar values, colors
        will be generated using the colormap `cmap` and normalizer `norm`.
    %(axes.apply_cmap)s
    lw, linewidth, linewidths, markeredgewidth, markeredgewidths : \
float or list thereof, optional
        The marker edge width.
    edgecolors, markeredgecolor, markeredgecolors : \
color-spec or list thereof, optional
        The marker edge color.

    Other parameters
    ----------------
    **kwargs
        Passed to `~matplotlib.axes.Axes.scatter`.

    See also
    --------
    standardize_1d
    indicate_error
    apply_cycle
    """
    # Manage input arguments
    # NOTE: Parse 1d must come before this
    if len(args) > 4:
        raise ValueError(f'Expected 1-4 positional arguments, got {len(args)}.')
    args = list(args)
    if len(args) == 4:
        c = args.pop(1)
    if len(args) == 3:
        s = args.pop(0)

    # Apply some aliases for keyword arguments
    c = _not_none(c=c, color=color, markercolor=markercolor)
    s = _not_none(s=s, size=size, markersize=markersize)
    lw = _not_none(
        lw=lw, linewidth=linewidth, linewidths=linewidths,
        markeredgewidth=markeredgewidth, markeredgewidths=markeredgewidths,
    )
    ec = _not_none(
        edgecolor=edgecolor, edgecolors=edgecolors,
        markeredgecolor=markeredgecolor, markeredgecolors=markeredgecolors,
    )

    # Get colormap
    cmap_kw = cmap_kw or {}
    if cmap is not None:
        cmap = constructor.Colormap(cmap, **cmap_kw)

    # Get normalizer and levels
    ticks = None
    carray = np.atleast_1d(c)
    if (
        np.issubdtype(carray.dtype, np.number)
        and not (carray.ndim == 2 and carray.shape[1] in (3, 4))
    ):
        carray = carray.ravel()
        norm, cmap, _, ticks = _build_discrete_norm(
            carray,  # sample data for getting suitable levels
            N=N, levels=levels, values=values,
            cmap=cmap, norm=norm, norm_kw=norm_kw, vmin=vmin, vmax=vmax, extend=extend,
            symmetric=symmetric, locator=locator, locator_kw=locator_kw,
        )

    # Fix 2D arguments but still support scatter(x_vector, y_2d) usage
    # NOTE: Since we are flattening vectors the coordinate metadata is meaningless,
    # so converting to ndarray and stripping metadata is no problem.
    # NOTE: numpy.ravel() preserves masked arrays
    if len(args) == 2 and all(_to_ndarray(arg).squeeze().ndim > 1 for arg in args):
        args = tuple(np.ravel(arg) for arg in args)

    # Scale s array
    if np.iterable(s) and (smin is not None or smax is not None):
        smin_true, smax_true = min(s), max(s)
        if smin is None:
            smin = smin_true
        if smax is None:
            smax = smax_true
        s = smin + (smax - smin) * (np.array(s) - smin_true) / (smax_true - smin_true)
    obj = objs = func(
        self, *args, c=c, s=s, cmap=cmap, norm=norm,
        linewidths=lw, edgecolors=ec, **kwargs
    )
    if not isinstance(objs, tuple):
        objs = (obj,)
    for iobj in objs:
        iobj._colorbar_extend = extend
        iobj._colorbar_ticks = ticks
    return obj


def _fill_between_apply(
    self, func, *args, negcolor=None, poscolor=None, negpos=None,
    lw=None, linewidth=None, stack=None, stacked=None, where=None, **kwargs
):
    """
    Helper function that powers `fill_between` and `fill_betweenx`.
    """
    # Parse input arguments as follows:
    # * Permit using 'x', 'y1', and 'y2' or 'y', 'x1', and 'x2' as keyword
    #   arguments.
    # * When negpos is True, use fill_between(x, y1=0, y2) as the default
    #   instead of fill_between(x, y1, y2=0).
    name = func.__name__
    sx = 'y' if 'x' in name else 'x'  # i.e. fill_betweenx
    sy = 'x' if sx == 'y' else 'y'
    stack = _not_none(stack=stack, stacked=stacked)
    args = list(args)
    if sx in kwargs:  # keyword 'x'
        args.insert(0, kwargs.pop(sx))
    if len(args) == 1:
        args.insert(0, np.arange(len(args[0])))
    for yi in (sy + '1', sy + '2'):
        if yi in kwargs:  # keyword 'y'
            args.append(kwargs.pop(yi))
    if len(args) == 2:
        args.append(0)
    elif len(args) == 3 and stack:  # ignore argument 3 down-the-line
        warnings._warn_proplot(
            f'{name}() cannot have three positional arguments with stack=True. '
            'Ignoring third argument.'
        )
    elif len(args) != 3:
        raise ValueError(f'Expected 2-3 positional args, got {len(args)}.')

    # Draw patches with default edge width zero
    # TODO: Test 'negpos' with 3 arguments (i.e. zero location is not at zero)
    x, y1, y2 = args
    x = _to_arraylike(x)
    y1 = _to_arraylike(y1)
    y2 = _to_arraylike(y2)
    kwargs.update({'linewidth': _not_none(lw=lw, linewidth=linewidth, default=0)})
    if not negpos or kwargs.get('color') is not None:
        # Plot basic patches
        kwargs.update({'stack': stack, 'where': where})
        result = func(self, x, y1, y2, **kwargs)
        objs = (result,)
    else:
        # Plot negative and positive patches
        kwargs.setdefault('interpolate', True)
        message = (
            f'{name}() argument {{}}={{!r}} is incompatible with negpos=True. Ignoring.'
        )
        if where is not None:
            warnings._warn_proplot(message.format('where', where))
        if stack:
            warnings._warn_proplot(message.format('stack', stack))
        if y1.ndim > 1 or y2.ndim > 1:
            raise ValueError(f'{name}() arguments with negpos=True must be 1D.')
        negcolor = _not_none(negcolor, rc['negcolor'])
        poscolor = _not_none(poscolor, rc['poscolor'])
        obj1 = func(self, x, y1, y2, where=(y1 < y2), color=negcolor, **kwargs)
        obj2 = func(self, x, y1, y2, where=(y1 >= y2), color=poscolor, **kwargs)
        result = objs = (obj1, obj2)  # may be tuple of tuples due to apply_cycle

    # Add sticky edges in x-direction, and sticky edges in y-direction
    # *only* if one of the y limits is scalar. This should satisfy most users.
    # NOTE: Could also retrieve data from PolyCollection but that's tricky.
    xsides = (np.min(_to_ndarray(x)), np.max(_to_ndarray(x)))
    ysides = []
    if y1.size == 1:
        ysides.append(_to_ndarray(y1).item())
    if y2.size == 1:
        ysides.append(_to_ndarray(y2).item())
    objs = tuple(obj for _ in objs for obj in (_ if isinstance(_, tuple) else (_,)))
    for obj in objs:
        for s, sides in zip((sx, sy), (xsides, ysides)):
            convert = getattr(self, 'convert_' + s + 'units')
            edges = getattr(obj.sticky_edges, s)
            edges.extend(convert(sides))

    return result


@docstring.add_snippets
def fill_between_extras(self, func, *args, **kwargs):
    """
    %(axes.fill_between)s
    """
    return _fill_between_apply(self, func, *args, **kwargs)


@docstring.add_snippets
def fill_betweenx_extras(self, func, *args, **kwargs):
    """
    %(axes.fill_betweenx)s
    """
    return _fill_between_apply(self, func, *args, **kwargs)


def _hist_extras(self, func, x, bins=None, **kwargs):
    """
    Forces `bar_extras` to interpret `width` as literal rather than relative
    to step size and enforces all arguments after `bins` are keyword-only.
    """
    with _state_context(self, _absolute_bar_width=True):
        return func(self, x, bins=bins, **kwargs)


@docstring.add_snippets
def bar_extras(
    self, func, x=None, height=None, width=0.8, bottom=None, *,
    vert=None, orientation='vertical', stack=None, stacked=None,
    lw=None, linewidth=None, edgecolor='black',  # default edge color black
    negpos=False, negcolor=None, poscolor=None,
    **kwargs
):
    """
    %(axes.bar)s
    """
    # Parse arguments
    # WARNING: Implementation is really weird... we flip around arguments for horizontal
    # plots only to flip them back in apply_cycle when iterating through columns.
    if vert is not None:
        orientation = 'vertical' if vert else 'horizontal'
    if orientation == 'horizontal':
        x, bottom, width, height = bottom, x, height, width

    # Parse args
    # TODO: Stacked feature is implemented in `apply_cycle`, but makes more
    # sense do document here; figure out way to move it here?
    if kwargs.get('left', None) is not None:
        warnings._warn_proplot('bar() keyword "left" is deprecated. Use "x" instead.')
        x = kwargs.pop('left')
    if x is None and height is None:
        raise ValueError('bar() requires at least 1 positional argument, got 0.')
    elif height is None:
        x, height = None, x
    args = (x, height)
    stack = _not_none(stack=stack, stacked=stacked)
    linewidth = _not_none(lw=lw, linewidth=linewidth, default=rc['patch.linewidth'])
    kwargs.update({'width': width, 'bottom': bottom, 'orientation': orientation})
    kwargs.update({'linewidth': linewidth, 'edgecolor': edgecolor})

    # Call func
    # NOTE: This *must* also be wrapped by apply_cycle, which ultimately
    # permutes back the x/bottom args for horizontal bars! Need to clean up.
    if not negpos or kwargs.get('color') is not None:
        # Draw simple bars
        kwargs['stack'] = stack
        result = func(self, *args, **kwargs)
    else:
        # Draw negative and positive bars
        # NOTE: apply_cycle makes bar widths *relative* to step size between
        # x coordinates so cannot just omit data. Instead make some heights nan.
        message = 'bar() argument {}={!r} is incompatible with negpos=True. Ignoring.'
        stack = kwargs.pop('stack', None)
        if stack:
            warnings._warn_proplot(message.format('stack', stack))
        height = _to_ndarray(height)
        if height.ndim > 1:
            raise ValueError('bar() heights with negpos=True must be 1D.')
        height1 = height.astype(np.float64)
        height1[height >= 0] = np.nan
        height2 = height.astype(np.float64)
        height2[height < 0] = np.nan
        negcolor = _not_none(negcolor, rc['negcolor'])
        poscolor = _not_none(poscolor, rc['poscolor'])
        obj1 = func(self, x, height1, color=negcolor, **kwargs)
        obj2 = func(self, x, height2, color=poscolor, **kwargs)
        result = (obj1, obj2)
    return result


@docstring.add_snippets
def barh_extras(self, func, y=None, right=None, width=0.8, left=None, **kwargs):
    """
    %(axes.barh)s
    """
    # NOTE: Also the second positional argument is called 'right' so that 'width'
    # always means the width of *bars*.
    # NOTE: This reverses the input arguemnts by setting y-->bottom, left-->x,
    # width-->height, height-->width. Arguments then passed through standardize_1d
    # indicate_error and reversed back by apply_cycle.
    # NOTE: You *must* do juggling of barh keyword order --> bar keyword order
    # --> barh keyword order, because horizontal hist passes arguments to bar
    # directly and will not use a 'barh' method with overridden argument order!
    func  # avoid U100 error
    height = _not_none(height=kwargs.pop('height', None), width=width, default=0.8)
    kwargs.setdefault('orientation', 'horizontal')
    if y is None and width is None:
        raise ValueError('barh() requires at least 1 positional argument, got 0.')
    return self.bar(x=left, width=right, height=height, bottom=y, **kwargs)


def boxplot_extras(
    self, func, *args,
    orientation=None, means=None,
    fill=True, fillcolor=None, fillalpha=None,
    lw=None, linewidth=None,
    color=None, edgecolor=None,
    boxcolor=None, boxlw=None, boxlinewidth=None,
    capcolor=None, caplw=None, caplinewidth=None,
    whiskercolor=None, whiskerlw=None, whiskerlinewidth=None,
    fliercolor=None, flierlw=None, flierlinewidth=None,  # fliers have no line width
    meancolor=None, meanlw=None, meanlinewidth=None,
    mediancolor=None, medianlw=None, medianlinewidth=None,
    meanls=None, meanlinestyle=None, medianls=None, medianlinestyle=None,
    marker=None, markersize=None,
    **kwargs
):
    """
    Adds convenient keyword arguments and changes the default boxplot style.

    Important
    ---------
    This function wraps {methods}

    Parameters
    ----------
    *args : 1D or 2D ndarray
        The data array.
    vert : bool, optional
        If ``False``, box plots are drawn horizontally.
    orientation : {{None, 'vertical', 'horizontal'}}, optional
        Alternative to the native `vert` keyword arg. Added for
        consistency with `~matplotlib.axes.Axes.bar`.
    means : bool, optional
        If ``True``, this passes ``showmeans=True`` and ``meanline=True`` to
        `~matplotlib.axes.Axes.boxplot`.
    fill : bool, optional
        Whether to fill the box with a color.
    fillcolor : color-spec, list, optional
        The fill color for the boxes. Default is the next color cycler color. If
        a list, it should be the same length as the number of objects.
    fillalpha : float, optional
        The opacity of the boxes. Default is ``0.7``. If a list, should be
        the same length as the number of objects.
    lw, linewidth : float, optional
        The linewidth of all objects. Default is ``0.8``.
    color, edgecolor : color-spec, list, optional
        The color of all objects. Defalut is ``'black'``. If a list, it should
        be the same length as the number of objects.
    meanls, medianls, meanlinestyle, medianlinestyle : line style-spec, optional
        The line style for the mean and median lines drawn horizontally
        across the box.
    boxcolor, capcolor, whiskercolor, fliercolor, meancolor, mediancolor : \
color-spec, list, optional
        The color of various boxplot components. If a list, it should be the
        same length as the number of objects. These are shorthands so you don't
        have to pass e.g. a ``boxprops`` dictionary.
    boxlw, caplw, whiskerlw, flierlw, meanlw, medianlw, boxlinewidth, caplinewidth, \
meanlinewidth, medianlinewidth, whiskerlinewidth, flierlinewidth : float, optional
        The line width of various boxplot components. These are shorthands so
        you don't have to pass e.g. a ``boxprops`` dictionary.
    marker : marker-spec, optional
        Marker style for the 'fliers', i.e. outliers.
    markersize : float, optional
        Marker size for the 'fliers', i.e. outliers.

    Other parameters
    ----------------
    **kwargs
        Passed to the matplotlib plotting method.

    See also
    --------
    proplot.axes.Axes.boxes
    standardize_1d
    indicate_error
    apply_cycle
    """
    # Parse keyword args
    fill = fill is True or fillcolor is not None or fillalpha is not None
    fillalpha = _not_none(fillalpha, default=0.7)
    if fill and fillcolor is None:
        cycler = next(self._get_lines.prop_cycler)
        fillcolor = cycler.get('color', None)
    color = _not_none(color=color, edgecolor=edgecolor, default='black')
    linewidth = _not_none(lw=lw, linewidth=linewidth, default=0.8)
    boxlinewidth = _not_none(boxlw=boxlw, boxlinewidth=boxlinewidth)
    caplinewidth = _not_none(caplw=caplw, caplinewidth=caplinewidth)
    whiskerlinewidth = _not_none(whiskerlw=whiskerlw, whiskerlinewidth=whiskerlinewidth)
    flierlinewidth = _not_none(flierlw=flierlw, flierlinewidth=flierlinewidth)
    meanlinewidth = _not_none(meanlw=meanlw, meanlinewidth=meanlinewidth)
    medianlinewidth = _not_none(medianlw=medianlw, medianlinewidth=medianlinewidth)
    meanlinestyle = _not_none(meanls=meanls, meanlinestyle=meanlinestyle)
    medianlinestyle = _not_none(medianls=medianls, medianlinestyle=medianlinestyle)
    if means or kwargs.get('showmeans'):
        kwargs['showmeans'] = kwargs['meanline'] = True
    if orientation == 'horizontal':
        kwargs['vert'] = False
    elif orientation == 'vertical':
        kwargs['vert'] = True
    elif orientation is not None:
        raise ValueError("Orientation must be 'horizontal' or 'vertical', got {orientation!r}.")  # noqa: E501

    # Call function
    if len(args) > 2:
        raise ValueError(f'Expected 1-2 positional args, got {len(args)}.')
    obj = func(self, *args, **kwargs)
    if not args:
        return obj

    # Modify artist settings
    # TODO: Pass props keyword args instead? Maybe does not matter.
    for key, icolor, ilinewidth, ilinestyle in (
        ('boxes', boxcolor, boxlinewidth, None),
        ('caps', capcolor, caplinewidth, None),
        ('whiskers', whiskercolor, whiskerlinewidth, None),
        ('fliers', fliercolor, flierlinewidth, None),
        ('means', meancolor, meanlinewidth, meanlinestyle),
        ('medians', mediancolor, medianlinewidth, medianlinestyle),
    ):
        if key not in obj:  # possible if not rendered
            continue
        artists = obj[key]
        icolor = _not_none(icolor, color)
        ilinewidth = _not_none(ilinewidth, linewidth)
        if not isinstance(fillalpha, list):
            fillalpha = [fillalpha] * len(artists)
        if not isinstance(fillcolor, list):
            fillcolor = [fillcolor] * len(artists)
        for i, artist in enumerate(artists):
            # Lines used for boxplot components
            jcolor = icolor
            if isinstance(icolor, list):
                jcolor = icolor[i // 2 if key in ('caps', 'whiskers') else i]
            if ilinestyle is not None:
                artist.set_linestyle(ilinestyle)
            if ilinewidth is not None:
                artist.set_linewidth(ilinewidth)
                artist.set_markeredgewidth(ilinewidth)
            if jcolor is not None:
                artist.set_color(jcolor)
                artist.set_markeredgecolor(jcolor)
            # "Filled" boxplot by adding patch beneath line path
            if fill and key == 'boxes':
                patch = mpatches.PathPatch(
                    artist.get_path(), linewidth=0,
                    facecolor=fillcolor[i], alpha=fillalpha[i]
                )
                self.add_artist(patch)
            # Outlier markers
            if key == 'fliers':
                if marker is not None:
                    artist.set_marker(marker)
                if markersize is not None:
                    artist.set_markersize(markersize)

    return obj


def violinplot_extras(
    self, func, *args, orientation=None,
    fillcolor=None, fillalpha=None,
    lw=None, linewidth=None,
    color=None, edgecolor=None,
    **kwargs
):
    """
    Adds convenient keyword arguments and changes the default violinplot style
    to match `this matplotlib example \
<https://matplotlib.org/3.1.0/gallery/statistics/customized_violin.html>`__.
    It is also no longer possible to show minima and maxima with whiskers --
    while this is useful for `~matplotlib.axes.Axes.boxplot`\\ s it is
    redundant for `~matplotlib.axes.Axes.violinplot`\\ s.

    Important
    ---------
    This function wraps {methods}

    Parameters
    ----------
    *args : 1D or 2D ndarray
        The data array.
    vert : bool, optional
        If ``False``, box plots are drawn horizontally.
    orientation : {{None, 'vertical', 'horizontal'}}, optional
        Alternative to the native `vert` keyword arg.
        Added for consistency with `~matplotlib.axes.Axes.bar`.
    fillcolor : color-spec, list, optional
        The violin plot fill color. Default is the next color cycler color. If
        a list, it should be the same length as the number of objects.
    fillalpha : float, optional
        The opacity of the violins. Default is ``0.7``. If a list, it
        should be the same length as the number of objects.
    lw, linewidth : float, optional
        The linewidth of the line objects. Default is ``0.8``.
    color, edgecolor : color-spec, list, optional
        The edge color for the violin patches. Default is ``'black'``. If a
        list, it should be the same length as the number of objects.

    Other parameters
    ----------------
    **kwargs
        Passed to `~matplotlib.axes.Axes.violinplot`.

    See also
    --------
    proplot.axes.Axes.violins
    standardize_1d
    indicate_error
    apply_cycle
    """
    # Parse keyword args
    # NOTE: Some of these are caught by indicate_error for drawing error bars
    if orientation == 'horizontal':
        kwargs['vert'] = False
    elif orientation == 'vertical':
        kwargs['vert'] = True
    elif orientation is not None:
        raise ValueError("Orientation must be 'horizontal' or 'vertical', got {orientation!r}.")  # noqa: E501
    color = _not_none(color=color, edgecolor=edgecolor, default='black')
    linewidth = _not_none(lw=lw, linewidth=linewidth, default=0.8)
    fillalpha = _not_none(fillalpha, default=0.7)
    if kwargs.pop('showextrema', None):
        warnings._warn_proplot('Ignoring showextrema=True.')
    if 'showmeans' in kwargs:  # native argument overridden by proplot handling
        kwargs.setdefault('mean', kwargs.pop('showmeans'))
    if 'showmedians' in kwargs:  # native argument overridden by proplot handling
        kwargs.setdefault('median', kwargs.pop('showmedians'))
    kwargs.update({'showmeans': False, 'showmedians': False, 'showextrema': False})
    kwargs.setdefault('capsize', 0)  # caps are redundant for violin plots

    # Call function
    if len(args) > 2:
        raise ValueError(f'Expected 1-2 positional args, got {len(args)}.')
    obj = result = func(self, *args, linewidth=linewidth, **kwargs)
    if not args:
        return result

    # Modify body settings
    if isinstance(obj, (list, tuple)):
        obj = obj[0]
    artists = obj['bodies']
    if not isinstance(fillalpha, list):
        fillalpha = [fillalpha] * len(artists)
    if not isinstance(fillcolor, list):
        fillcolor = [fillcolor] * len(artists)
    if not isinstance(color, list):
        color = [color] * len(artists)
    for i, artist in enumerate(artists):
        artist.set_linewidths(linewidth)
        if fillalpha[i] is not None:
            artist.set_alpha(fillalpha[i])
        if fillcolor[i] is not None:
            artist.set_facecolor(fillcolor[i])
        if color[i] is not None:
            artist.set_edgecolor(color[i])

    return result


def _get_transform(self, transform):
    """
    Translates user input transform. Also used in an axes method.
    """
    try:
        from cartopy.crs import CRS
    except ModuleNotFoundError:
        CRS = None
    cartopy = getattr(self, 'name', None) == 'proplot_cartopy'
    if (
        isinstance(transform, mtransforms.Transform)
        or CRS and isinstance(transform, CRS)
    ):
        return transform
    elif transform == 'figure':
        return self.figure.transFigure
    elif transform == 'axes':
        return self.transAxes
    elif transform == 'data':
        return PlateCarree() if cartopy else self.transData
    elif cartopy and transform == 'map':
        return self.transData
    else:
        raise ValueError(f'Unknown transform {transform!r}.')


def _update_text(self, props):
    """
    Monkey patch that adds pseudo "border" and "bbox" properties to text objects without
    wrapping the entire class. Overrides update to facilitate updating inset titles.
    """
    props = props.copy()  # shallow copy

    # Update border
    border = props.pop('border', None)
    bordercolor = props.pop('bordercolor', 'w')
    borderinvert = props.pop('borderinvert', False)
    borderwidth = props.pop('borderwidth', 2)
    if border:
        facecolor, bgcolor = self.get_color(), bordercolor
        if borderinvert:
            facecolor, bgcolor = bgcolor, facecolor
        kwargs = {
            'linewidth': borderwidth,
            'foreground': bgcolor,
            'joinstyle': 'miter',
        }
        self.update({
            'color': facecolor,
            'path_effects': [mpatheffects.Stroke(**kwargs), mpatheffects.Normal()],
        })
    elif border is False:
        self.update({
            'path_effects': None,
        })

    # Update bounding box
    # NOTE: We use '_title_pad' and '_title_above' for both titles and a-b-c labels
    # because always want to keep them aligned.
    # NOTE: For some reason using pad / 10 results in perfect alignment. Matplotlib
    # docs are vague about bounding box units, maybe they are tens of points?
    bbox = props.pop('bbox', None)
    bboxcolor = props.pop('bboxcolor', 'w')
    bboxstyle = props.pop('bboxstyle', 'round')
    bboxalpha = props.pop('bboxalpha', 0.5)
    bboxpad = _not_none(props.pop('bboxpad', None), self.axes._title_pad / 10)
    if isinstance(bbox, dict):  # *native* matplotlib usage
        props['bbox'] = bbox
    elif bbox:
        self.set_bbox({
            'edgecolor': 'black',
            'facecolor': bboxcolor,
            'boxstyle': bboxstyle,
            'alpha': bboxalpha,
            'pad': bboxpad,
        })
    elif bbox is False:
        self.set_bbox(None)  # disables the bbox

    return type(self).update(self, props)


def text_extras(
    self, func,
    x=0, y=0, text='', transform='data',
    family=None, fontfamily=None, fontname=None, fontsize=None, size=None,
    border=False, bordercolor='w', borderwidth=2, borderinvert=False,
    bbox=False, bboxcolor='w', bboxstyle='round', bboxalpha=0.5, bboxpad=None,
    **kwargs
):
    """
    Enables specifying `tranform` with a string name and adds a feature for
    drawing borders and bbox around text.

    Important
    ---------
    This function wraps {methods}

    Parameters
    ----------
    x, y : float
        The *x* and *y* coordinates for the text.
    text : str
        The text string.
    transform : {{'data', 'axes', 'figure'}} \
or `~matplotlib.transforms.Transform`, optional
        The transform used to interpret `x` and `y`. Can be a
        `~matplotlib.transforms.Transform` object or a string corresponding to
        `~matplotlib.axes.Axes.transData`, `~matplotlib.axes.Axes.transAxes`,
        or `~matplotlib.figure.Figure.transFigure` transforms. Default is
        ``'data'``, i.e. the text is positioned in data coordinates.
    fontsize, size : float or str, optional
        The font size. If float, units are inches. If string, units are
        interpreted by `~proplot.utils.units`.
    fontname, fontfamily, family : str, optional
        The font name (e.g. ``'Fira Math'``) or font family name (e.g. ``'serif'``).
        Matplotlib falls back to the system default if not found.
    fontweight, weight, fontstyle, style, fontvariant, variant : str, optional
        Additional font properties. See `~matplotlib.text.Text` for details.
    border : bool, optional
        Whether to draw border around text.
    borderwidth : float, optional
        The width of the text border. Default is ``2`` points.
    bordercolor : color-spec, optional
        The color of the text border. Default is ``'w'``.
    borderinvert : bool, optional
        If ``True``, the text and border colors are swapped.
    bbox : bool, optional
        Whether to draw a bounding box around text.
    bboxcolor : color-spec, optional
        The color of the text bounding box. Default is ``'w'``.
    bboxstyle : boxstyle, optional
        The style of the bounding box. Default is ``'round'``.
    bboxalpha : float, optional
        The alpha for the bounding box. Default is ``'0.5'``.
    bboxpad : float, optional
        The padding for the bounding box. Default is :rc:`title.bboxpad`.

    Other parameters
    ----------------
    **kwargs
        Passed to `~matplotlib.axes.Axes.text`.
    """
    # Parse input args
    # NOTE: Previously issued warning if fontname did not match any of names
    # in ttflist but this would result in warning for e.g. family='sans-serif'.
    # Matplotlib font API makes it very difficult to inject warning in
    # correct place. Simpler to just
    # NOTE: Do not emit warning if user supplied conflicting properties
    # because matplotlib has like 100 conflicting text properties for which
    # it doesn't emit warnings. Prefer not to fix all of them.
    fontsize = _not_none(fontsize, size)
    fontfamily = _not_none(fontname, fontfamily, family)
    if fontsize is not None:
        try:
            rc._scale_font(fontsize)  # *validate* but do not translate
        except KeyError:
            fontsize = units(fontsize, 'pt')
        kwargs['fontsize'] = fontsize
    if fontfamily is not None:
        kwargs['fontfamily'] = fontfamily
    if not transform:
        transform = self.transData
    else:
        transform = _get_transform(self, transform)

    # Apply monkey patch to text object
    # TODO: Why only support this here, and not in arbitrary places throughout
    # rest of matplotlib API? Units engine needs better implementation.
    obj = func(self, x, y, text, transform=transform, **kwargs)
    obj.update = _update_text.__get__(obj)
    obj.update({
        'border': border,
        'bordercolor': bordercolor,
        'borderinvert': borderinvert,
        'borderwidth': borderwidth,
        'bbox': bbox,
        'bboxcolor': bboxcolor,
        'bboxstyle': bboxstyle,
        'bboxalpha': bboxalpha,
        'bboxpad': bboxpad,
    })
    return obj


def _convert_bar_width(x, width=1, ncols=1):
    """
    Convert bar plot widths from relative to coordinate spacing. Relative
    widths are much more convenient for users.
    """
    # WARNING: This will fail for non-numeric non-datetime64 singleton
    # datatypes but this is good enough for vast majority of cases.
    x_test = np.atleast_1d(_to_ndarray(x))
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
    return width * x_step / ncols


def _iter_objs_labels(objs):
    """
    Retrieve the (object, label) pairs for objects with actual labels
    from nested lists and tuples of objects.
    """
    # Account for (1) multiple columns of data, (2) functions that return
    # multiple values (e.g. hist() returns (bins, values, patches)), and
    # (3) matplotlib.Collection list subclasses.
    label = _get_label(objs)
    if label:
        yield (objs, label)
    elif isinstance(objs, list):
        for obj in objs:
            yield from _iter_objs_labels(obj)


def _update_cycle(self, cycle, scatter=False, **kwargs):
    """
    Try to update the `~cycler.Cycler` without resetting it if it has not changed.
    Also return keys that should be explicitly iterated over for commands that
    otherwise don't use the property cycler (currently just scatter).
    """
    # Get the original property cycle
    # NOTE: Matplotlib saves itertools.cycle(cycler), not the original
    # cycler object, so we must build up the keys again.
    # NOTE: Axes cycle has no getter, only set_prop_cycle, which sets a
    # prop_cycler attribute on the hidden _get_lines and _get_patches_for_fill
    # objects. This is the only way to query current axes cycler! Should not
    # wrap set_prop_cycle because would get messy and fragile.
    # NOTE: The _get_lines cycler is an *itertools cycler*. Has no length, so
    # we must cycle over it with next(). We try calling next() the same number
    # of times as the length of input cycle. If the input cycle *is* in fact
    # the same, below does not reset the color position, cycles us to start!
    i = 0
    by_key = {}
    cycle_orig = self._get_lines.prop_cycler
    for i in range(len(cycle)):  # use the cycler object length as a guess
        prop = next(cycle_orig)
        for key, value in prop.items():
            if key not in by_key:
                by_key[key] = set()
            if isinstance(value, (list, np.ndarray)):
                value = tuple(value)
            by_key[key].add(value)

    # Reset property cycler if it differs
    reset = set(by_key) != set(cycle.by_key())
    if not reset:  # test individual entries
        for key, value in cycle.by_key().items():
            if by_key[key] != set(value):
                reset = True
                break
    if reset:
        self.set_prop_cycle(cycle)

    # Psuedo-expansion of matplotlib's native property cycling for scatter(). Check out
    # the current property cycle for extra keys not interpreted by matplotlib (all keys
    # other than color, linestyle, and dashes). If they are present, and the user didn't
    # explicitly override them with some other keyword arg, then take note of that.
    # NOTE: By default matplotlib uses _get_patches_for_fill.get_next_color
    # for scatter next scatter color, but cannot get anything else! We simultaneously
    # iterate through the _get_lines property cycler and apply relevant properties.
    apply_manually = {}  # which keys to apply from property cycler
    if scatter:
        prop_keys = set(self._get_lines._prop_keys) - {'color', 'linestyle', 'dashes'}
        for prop, key in (
            ('markersize', 's'),
            ('linewidth', 'linewidths'),
            ('markeredgewidth', 'linewidths'),
            ('markeredgecolor', 'edgecolors'),
            ('alpha', 'alpha'),
            ('marker', 'marker'),
        ):
            value = kwargs.get(key, None)  # a apply_cycle argument
            if prop in prop_keys and value is None:  # if key in cycler and prop unset
                apply_manually[prop] = key

    return apply_manually  # set indicating additional keys we cycle through


def apply_cycle(
    self, func, *args,
    cycle=None, cycle_kw=None,
    label=None, labels=None, values=None,
    legend=None, legend_kw=None,
    colorbar=None, colorbar_kw=None,
    _errobjs=None, **kwargs
):
    """
    Adds features for controlling colors in the property cycler and drawing
    legends or colorbars in one go.

    Important
    ---------
    This function wraps {methods}

    This wrapper also *standardizes acceptable input* -- these methods now all
    accept 2D arrays holding columns of data, and *x*-coordinates are always
    optional. Note this alters the behavior of `~matplotlib.axes.Axes.boxplot`
    and `~matplotlib.axes.Axes.violinplot`, which now compile statistics on
    *columns* of data instead of *rows*.

    Parameters
    ----------
    cycle : cycle-spec, optional
        The cycle specifer, passed to the `~proplot.constructor.Cycle`
        constructor. If the returned list of colors is unchanged from the
        current axes color cycler, the axes cycle will **not** be reset to the
        first position.
    cycle_kw : dict-like, optional
        Passed to `~proplot.constructor.Cycle`.
    label : float or str, optional
        The legend label to be used for this plotted element.
    labels, values : list of float or list of str, optional
        Used with 2D input arrays. The legend labels or colorbar coordinates for
        each column in the array. Can be numeric or string, and must match
        the number of columns in the 2D array.
    legend : bool, int, or str, optional
        If not ``None``, this is a location specifying where to draw an *inset*
        or *panel* legend from the resulting handle(s). If ``True``, the
        default location is used. Valid locations are described in
        `~proplot.axes.Axes.legend`.
    legend_kw : dict-like, optional
        Ignored if `legend` is ``None``. Extra keyword args for our call
        to `~proplot.axes.Axes.legend`.
    colorbar : bool, int, or str, optional
        If not ``None``, this is a location specifying where to draw an *inset*
        or *panel* colorbar from the resulting handle(s). If ``True``, the
        default location is used. Valid locations are described in
        `~proplot.axes.Axes.colorbar`.
    colorbar_kw : dict-like, optional
        Ignored if `colorbar` is ``None``. Extra keyword args for our call
        to `~proplot.axes.Axes.colorbar`.

    Other parameters
    ----------------
    *args, **kwargs
        Passed to the matplotlib plotting method.

    See also
    --------
    standardize_1d
    indicate_error
    proplot.constructor.Cycle
    proplot.constructor.Colors
    """
    # NOTE: Requires standardize_1d wrapper before reaching this. Also note
    # that the 'x' coordinates are sometimes ignored below.
    name = func.__name__
    plot = name in ('plot',)
    scatter = name in ('scatter',)
    fill = name in ('fill_between', 'fill_betweenx')
    hist = name in ('hist',)
    bar = name in ('bar',)
    pie = name in ('pie',)
    box = name in ('boxplot', 'violinplot')
    violin = name in ('violinplot',)
    if not args:
        return func(self, *args, **kwargs)
    x, y, *args = args
    ys = (y,)
    if fill and len(args) >= 1:
        *ys, args = y, args[0], args[1:]
    ncols = 1 if pie or box or y.ndim == 1 else y.shape[1]

    # Parse keyword args
    # NOTE: Already pull out labels/values from legend_kw/colorbar_kw in _auto_format_1d
    cycle_kw = cycle_kw or {}
    legend_kw = legend_kw or {}
    colorbar_kw = colorbar_kw or {}
    labels = _not_none(label=label, values=values, labels=labels)
    if pie:  # add x coordinates as default pie chart labels
        labels = _not_none(labels, x)  # TODO: move to pie wrapper?

    # Bar plot width and origin
    barh = stack = False
    if bar or fill:
        stack = kwargs.pop('stack', False)
    if bar:
        barh = kwargs.get('orientation', None) == 'horizontal'
        width = kwargs.pop('width', 0.8)  # 'width' for bar *and* barh (see bar_extras)
        if not stack and not getattr(self, '_absolute_bar_width', None):
            width = _convert_bar_width(x, width, ncols)
        kwargs['height' if barh else 'width'] = width
        kwargs.setdefault('x' if barh else 'bottom', 0)

    # Update the property cycler
    apply_manually = {}
    if cycle is not None or cycle_kw:
        if y.ndim > 1 and y.shape[1] > 1:  # default samples count
            cycle_kw.setdefault('N', y.shape[1])
        cycle_args = () if cycle is None else (cycle,)
        cycle = constructor.Cycle(*cycle_args, **cycle_kw)
        apply_manually = _update_cycle(self, cycle, scatter=scatter, **kwargs)

    # Handle legend labels. Several scenarios:
    # 1. Always prefer input labels
    # 2. Always add labels if this is a *named* dimension.
    # 3. Even if not *named* dimension add labels if labels are string
    # Then convert labels to string because e.g. scatter() applies default label
    # if input is False-ey. So numeric '0' would be overridden.
    if pie or box:
        # Functions handle multiple labels on their own
        # NOTE: Using boxplot() without this will overwrite labels previously
        # set along the x-axis by _auto_format_1d.
        if not violin and labels is not None:
            kwargs['labels'] = _to_ndarray(labels)  # error raised down the line
    else:
        # Check and standardize labels
        # NOTE: Must convert to ndarray or can get singleton DataArrays
        if not np.iterable(labels) or isinstance(labels, str):
            labels = [labels] * ncols
        if len(labels) != ncols:
            raise ValueError(f'Array has {ncols} columns but got {len(labels)} labels.')
        if labels is not None:
            labels = [str(_not_none(label, '')) for label in _to_ndarray(labels)]
        else:
            labels = [None] * ncols

    # Plot successive columns
    objs = []
    for i in range(ncols):
        # Property cycling for scatter plots
        kw = kwargs.copy()
        if apply_manually:
            props = next(self._get_lines.prop_cycler)
        for prop, key in apply_manually.items():
            kw[key] = props[prop]

        # The x coordinates for bar plots
        ix, iy, *_ = x, *ys  # samples
        if bar and stack and iy.ndim > 1:
            kw['x' if barh else 'bottom'] = _to_indexer(iy)[:, :i].sum(axis=1)
        if bar and not stack:
            offset = width * (i - 0.5 * (ncols - 1))
            ix = x + offset

        # The y coordinates and labels
        # WARNING: If stack=True then we always *ignore* second argument passed to
        # fill_between. Warning should be issued by fill_between_extras in this case.
        if pie or box:  # only ever have one y value, cannot have legend labels
            iys = ys[:1]
        elif fill and stack:  # ignore argument 3 as warned in _fill_between_apply
            iys = tuple(iy if iy.ndim == 1 else _to_indexer(iy)[:, :ii].sum(axis=1) for ii in (i, i + 1))  # noqa: E501
            kw['label'] = labels[i] or None
        else:
            iys = tuple(iy if iy.ndim == 1 else _to_indexer(iy)[:, i] for iy in ys)
            kw['label'] = labels[i] or None

        # Call function for relevant column
        # TODO: Why pull result out of singleton list?
        iargs = ()
        if barh:  # special case, use kwargs only!
            kw.update({'bottom': ix, 'width': iys[0]})  # always single y vector
        elif pie or hist or box:
            iargs = iys
        else:  # has x-coordinates, and maybe more than one y
            iargs = ix, *iys
        obj = func(self, *iargs, *args, **kw)
        if isinstance(obj, (list, tuple)) and len(obj) == 1:
            obj = obj[0]
        objs.append(obj)

    # Add colorbar
    # NOTE: Colorbar will get the labels from the artists. Don't need to extract
    # them because can't have multiple-artist entries like for legend()
    if colorbar:
        self.colorbar(objs, loc=colorbar, queue=True, **colorbar_kw)

    # Add legend
    # NOTE: Put error bounds objects *before* line objects in the tuple
    # so that line gets drawn on top of bounds.
    # NOTE: If error objects have separate label, allocate separate legend entry.
    # If they do not, try to combine with current legend entry.
    if legend:
        if not isinstance(_errobjs, (list, tuple)):
            _errobjs = (_errobjs,)
        eobjs = [obj for obj in _errobjs if obj and not _get_label(obj)]
        hobjs = [(*eobjs[::-1], *objs)] if eobjs else objs.copy()
        hobjs.extend(obj for obj in _errobjs if obj and _get_label(obj))
        try:
            hobjs, labels = list(zip(*_iter_objs_labels(hobjs)))
        except ValueError:
            hobjs = labels = ()
        self.legend(hobjs, labels, loc=legend, queue=True, **legend_kw)

    # Return
    # WARNING: Make sure plot always returns tuple of objects, and bar always
    # returns singleton unless we have bulk drawn bar plots! Other matplotlib
    # methods call these internally and expect a certain output format!
    if plot:
        return tuple(objs)  # always return tuple of objects
    elif box:
        return objs[0]  # always return singleton
    else:
        return objs[0] if len(objs) == 1 else tuple(objs)


def _auto_levels_locator(
    *args, N=None, norm=None, norm_kw=None,
    vmin=None, vmax=None, extend='neither', locator=None, locator_kw=None,
    symmetric=False, positive=False, negative=False, nozero=False,
):
    """
    Automatically generate level locations based on the input data, the
    input locator, and the input normalizer.

    Parameters
    ----------
    *args
        The sample dataset(s).
    N : int, optional
        The (approximate) number of levels to create.
    norm, norm_kw
        Passed to `~proplot.constructor.Norm`. Used to determine suitable
        level locations if `locator` is not passed.
    vmin, vmax : float, optional
        The data limits.
    extend : str, optional
        The extend setting.
    locator, locator_kw
        Passed to `~proplot.constructor.Locator`. Used to determine suitable
        level locations.
    symmetric, positive, negative : bool, optional
        Whether the automatic levels should be symmetric, should be all positive
        with a minimum at zero, or should be all negative with a maximum at zero.
    nozero : bool, optional
        Whether zero should be excluded from automatic levels. This is also
        implemented in `apply_cmap` so that `nozero` can be used to remove user
        input levels (e.g. ``ax.contour(..., levels=plot.arange(-5, 5), nozero=True)``),
        but is replecated here so power users can use this function in isolation.

    Returns
    -------
    levels : ndarray
        The levels.
    locator : ndarray or `matplotlib.ticker.Locator`
        The locator used for colorbar tick locations.
    """
    if np.iterable(N):
        return N, N
    if N is None:
        N = 11
    norm_kw = norm_kw or {}
    locator_kw = locator_kw or {}
    norm = constructor.Norm(norm or 'linear', **norm_kw)
    if positive and negative:
        raise ValueError('Incompatible options: positive=True and negative=True.')
    if locator is not None:
        level_locator = tick_locator = constructor.Locator(locator, **locator_kw)
    elif isinstance(norm, mcolors.LogNorm):
        level_locator = tick_locator = mticker.LogLocator(**locator_kw)
    elif isinstance(norm, mcolors.SymLogNorm):
        locator_kw.setdefault('base', _flexible_getattr(norm, 'base', 10))
        locator_kw.setdefault('linthresh', _flexible_getattr(norm, 'linthresh', 1))
        level_locator = tick_locator = mticker.SymmetricalLogLocator(**locator_kw)
    else:
        nbins = N * 2 if positive or negative else N
        locator_kw.setdefault('symmetric', symmetric or positive or negative)
        level_locator = mticker.MaxNLocator(nbins, min_n_ticks=1, **locator_kw)
        tick_locator = None

    # Get locations
    automin = vmin is None
    automax = vmax is None
    if automin or automax:
        vmins = []
        vmaxs = []
        for data in args:
            data = ma.masked_invalid(data, copy=False)
            if automin:
                vmin = float(data.min())
            if automax:
                vmax = float(data.max())
            if vmin == vmax or ma.is_masked(vmin) or ma.is_masked(vmax):
                vmin, vmax = 0, 1
            vmins.append(vmin)
            vmaxs.append(vmax)
        vmin = min(vmins)
        vmax = max(vmaxs)
    try:
        levels = level_locator.tick_values(vmin, vmax)
    except RuntimeError:  # too-many-ticks error
        levels = np.linspace(vmin, vmax, N)  # TODO: _autolev used N+1

    # Trim excess levels the locator may have supplied
    # NOTE: This part is mostly copied from matplotlib _autolev
    if not locator_kw.get('symmetric', None):
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
    nn = N // len(levels)
    if nn >= 2:
        olevels = norm(levels)
        nlevels = []
        for i in range(len(levels) - 1):
            l1, l2 = olevels[i], olevels[i + 1]
            nlevels.extend(np.linspace(l1, l2, nn + 1)[:-1])
        nlevels.append(olevels[-1])
        levels = norm.inverse(nlevels)

    # Filter the remaining contours
    if nozero and 0 in levels:
        levels = levels[levels != 0]
    if positive:
        levels = levels[levels >= 0]
    if negative:
        levels = levels[levels <= 0]

    # Use auto-generated levels for ticks if still None
    locator = tick_locator or levels
    return levels, locator


def _build_discrete_norm(
    data=None, N=None, levels=None, values=None,
    cmap=None, norm=None, norm_kw=None, vmin=None, vmax=None, extend='neither',
    minlength=2,
    **kwargs,
):
    """
    Build a `~proplot.colors.DiscreteNorm` or `~proplot.colors.BoundaryNorm`
    from the input arguments. This automatically calculates "nice" level
    boundaries if they were not provided.

    Parameters
    ----------
    data : ndarray, optional
        The data.
    levels, values : ndarray, optional
        The explicit boundaries.
    norm, norm_kw
        Passed to `~proplot.constructor.Norm` and then to `DiscreteNorm`.
    vmin, vmax : float, optional
        The minimum and maximum values for the normalizer.
    cmap : `matplotlib.colors.Colormap`, optional
        The colormap. Passed to `DiscreteNorm`.
    extend : str, optional
        The extend setting.
    minlength : int, optional
        The minimum length for level lists.
    **kwargs
        Passed to `_auto_levels_locator`.

    Returns
    -------
    norm : `matplotlib.colors.Normalize`
        The normalizer.
    ticks : `numpy.ndarray` or `matplotlib.locator.Locator`
        The axis locator or the tick location candidates.
    """
    # Parse flexible keyword args
    norm_kw = norm_kw or {}
    levels = _not_none(
        N=N, levels=levels, norm_kw_levels=norm_kw.pop('levels', None),
        default=rc['image.levels']
    )
    vmin = _not_none(vmin=vmin, norm_kw_vmin=norm_kw.pop('vmin', None))
    vmax = _not_none(vmax=vmax, norm_kw_vmax=norm_kw.pop('vmax', None))
    if norm == 'segments':  # TODO: remove
        norm = 'segmented'

    # NOTE: Matplotlib colorbar algorithm *cannot* handle descending levels
    # so this function reverses them and adds special attribute to the
    # normalizer. Then colorbar_extras reads this attribute and flips the
    # axis and the colormap direction.
    # Check input levels and values
    for key, val in (('levels', levels), ('values', values)):
        if not np.iterable(val):
            continue
        if len(val) < minlength or len(val) >= 2 and any(
            np.sign(np.diff(val)) != np.sign(val[1] - val[0])
        ):
            raise ValueError(
                f'{key!r} must be monotonically increasing or decreasing '
                f'and at least length {minlength}, got {val}.'
            )

    # Get level edges from level centers
    locator = None
    if isinstance(values, Integral):
        levels = values + 1
    elif np.iterable(values) and len(values) == 1:
        levels = [values[0] - 1, values[0] + 1]  # weird but why not
    elif np.iterable(values) and len(values) > 1:
        # Try to generate levels such that a LinearSegmentedNorm will
        # place values ticks at the center of each colorbar level.
        # utils.edges works only for evenly spaced values arrays.
        # We solve for: (x1 + x2)/2 = y --> x2 = 2*y - x1
        # with arbitrary starting point x1. We also start the algorithm
        # on the end with *smaller* differences.
        if norm is None or norm == 'segmented':
            reverse = abs(values[-1] - values[-2]) < abs(values[1] - values[0])
            if reverse:
                values = values[::-1]
            levels = [values[0] - (values[1] - values[0]) / 2]
            for val in values:
                levels.append(2 * val - levels[-1])
            if reverse:
                levels = levels[::-1]
            if any(np.sign(np.diff(levels)) != np.sign(levels[1] - levels[0])):
                levels = edges(values)  # backup plan, weird tick locations
        # Generate levels by finding in-between points in the
        # normalized numeric space, e.g. LogNorm space.
        else:
            inorm = constructor.Norm(norm, **norm_kw)
            levels = inorm.inverse(edges(inorm(values)))
    elif values is not None:
        raise ValueError(
            f'Unexpected input values={values!r}. '
            'Must be integer or list of numbers.'
        )

    # Get default normalizer
    # Only use LinearSegmentedNorm if necessary, because it is slow
    descending = False
    if np.iterable(levels):
        if len(levels) == 1:
            norm = mcolors.Normalize(vmin=levels[0] - 1, vmax=levels[0] + 1)
        else:
            levels, descending = pcolors._check_levels(levels)
            if len(levels) > 2 and norm is None:
                steps = np.abs(np.diff(levels))
                eps = np.mean(steps) / 1e3
                if np.any(np.abs(np.diff(steps)) >= eps):
                    norm = 'segmented'
    if norm == 'segmented':
        if not np.iterable(levels):
            norm = 'linear'  # same result with improved speed
        else:
            norm_kw['levels'] = levels
    norm = constructor.Norm(norm or 'linear', **norm_kw)

    # Use the locator to determine levels
    # Mostly copied from the hidden contour.ContourSet._autolev
    # NOTE: Subsequently, we *only* use the locator to determine ticks if
    # *levels* and *values* were not passed.
    if isinstance(norm, mcolors.BoundaryNorm):
        # Get levels from bounds
        # TODO: Test this feature?
        # NOTE: No warning because we get here internally?
        levels = norm.boundaries
    elif np.iterable(values):
        # Prefer ticks in center, but subsample if necessary
        locator = _to_ndarray(values)
    elif np.iterable(levels):
        # Prefer ticks on level edges, but subsample if necessary
        locator = _to_ndarray(levels)
    else:
        # Determine levels automatically
        levels, locator = _auto_levels_locator(
            data, N=levels, norm=norm, vmin=vmin, vmax=vmax, extend=extend, **kwargs
        )

    # Generate DiscreteNorm and update "child" norm with vmin and vmax from
    # levels. This lets the colorbar set tick locations properly!
    # TODO: Move these to DiscreteNorm?
    if not isinstance(norm, mcolors.BoundaryNorm) and len(levels) > 1:
        norm = pcolors.DiscreteNorm(
            levels, cmap=cmap, norm=norm, descending=descending, unique=extend,
        )
    if descending:
        cmap = cmap.reversed()
    return norm, cmap, levels, locator


def _fix_white_lines(obj):
    """
    Fix white lines between between filled contours and mesh and fix issues with
    colormaps that are transparent.
    """
    # See: https://github.com/jklymak/contourfIssues
    # See: https://stackoverflow.com/q/15003353/4970632
    # 0.4pt is thick enough to hide lines but thin enough to not add "dots"
    # in corner of pcolor plots so good compromise.
    cmap = obj.cmap
    if not cmap._isinit:
        cmap._init()
    if all(cmap._lut[:-1, 3] == 1):  # skip for cmaps with transparency
        edgecolor = 'face'
    else:
        edgecolor = 'none'
    # Contour fixes
    # NOTE: This also covers TriContourSet returned by tricontour
    if isinstance(obj, mcontour.ContourSet):
        for contour in obj.collections:
            contour.set_edgecolor(edgecolor)
            contour.set_linewidth(0.4)
            contour.set_linestyle('-')
    # Pcolor fixes
    # NOTE: This ignores AxesImage and PcolorImage sometimes returned by pcolorfast
    elif isinstance(obj, (mcollections.PolyCollection, mcollections.QuadMesh)):
        if hasattr(obj, 'set_linewidth'):  # not always true for pcolorfast
            obj.set_linewidth(0.4)
        if hasattr(obj, 'set_edgecolor'):  # not always true for pcolorfast
            obj.set_edgecolor(edgecolor)


def _labels_contour(self, obj, *args, fmt=None, **kwargs):
    """
    Add labels to contours with support for shade-dependent filled contour labels
    flexible keyword arguments. Requires original arguments passed to contour function.
    """
    # Parse input args
    text_kw = {}
    for key in (*kwargs,):  # allow dict to change size
        if key not in (
            'levels', 'colors', 'fontsize', 'inline', 'inline_spacing',
            'manual', 'rightside_up', 'use_clabeltext',
        ):
            text_kw[key] = kwargs.pop(key)
    kwargs.setdefault('inline_spacing', 3)
    kwargs.setdefault('fontsize', rc['text.labelsize'])
    fmt = _not_none(fmt, pticker.SimpleFormatter())

    # Draw hidden additional contour for filled contour labels
    colors = None
    if obj.filled:  # guard against changes?
        obj = self.contour(*args, levels=obj.levels, linewidths=0)
        lums = [to_xyz(obj.cmap(obj.norm(level)), 'hcl')[2] for level in obj.levels]
        colors = ['w' if lum < 50 else 'k' for lum in lums]
    kwargs.setdefault('colors', colors)

    # Draw labels
    labs = obj.clabel(fmt=fmt, **kwargs)
    if labs is not None:  # returns None if no contours
        for lab in labs:
            lab.update(text_kw)

    return labs


def _labels_pcolor(self, obj, fmt=None, **kwargs):
    """
    Add labels to pcolor boxes with support for shade-dependent text colors.
    """
    # Parse input args and populate _facecolors, which is initially unfilled
    # See: https://stackoverflow.com/a/20998634/4970632
    fmt = _not_none(fmt, pticker.SimpleFormatter())
    labels_kw = {'size': rc['text.labelsize'], 'ha': 'center', 'va': 'center'}
    labels_kw.update(kwargs)
    obj.update_scalarmappable()  # populate _facecolors

    # Get positions and contour colors
    array = obj.get_array()
    paths = obj.get_paths()
    colors = _to_ndarray(obj.get_facecolors())
    edgecolors = _to_ndarray(obj.get_edgecolors())
    if len(colors) == 1:  # weird flex but okay
        colors = np.repeat(colors, len(array), axis=0)
    if len(edgecolors) == 1:
        edgecolors = np.repeat(edgecolors, len(array), axis=0)

    # Apply colors
    labs = []
    for i, (color, path, num) in enumerate(zip(colors, paths, array)):
        if not np.isfinite(num):
            edgecolors[i, :] = 0
            continue
        bbox = path.get_extents()
        x = (bbox.xmin + bbox.xmax) / 2
        y = (bbox.ymin + bbox.ymax) / 2
        if 'color' not in kwargs:
            _, _, lum = to_xyz(color, 'hcl')
            if lum < 50:
                color = 'w'
            else:
                color = 'k'
            labels_kw['color'] = color
        lab = self.text(x, y, fmt(num), **labels_kw)
        labs.append(lab)
    obj.set_edgecolors(edgecolors)

    return labs


@warnings._rename_kwargs('0.6', centers='values')
@docstring.add_snippets
def apply_cmap(
    self, func, *args,
    cmap=None, cmap_kw=None, norm=None, norm_kw=None,
    vmin=None, vmax=None, extend='neither', N=None, levels=None, values=None,
    symmetric=False, positive=False, negative=False, nozero=False,
    locator=None, locator_kw=None,
    edgefix=None, labels=False, labels_kw=None, fmt=None, precision=2,
    colorbar=False, colorbar_kw=None,
    lw=None, linewidth=None, linewidths=None,
    ls=None, linestyle=None, linestyles=None,
    color=None, colors=None, edgecolor=None, edgecolors=None,
    **kwargs
):
    """
    Adds several new keyword args and features for specifying the colormap,
    levels, and normalizers. Uses the `~proplot.colors.DiscreteNorm`
    normalizer to bin data into discrete color levels (see notes).

    Important
    ---------
    This function wraps {methods}

    Parameters
    ----------
    %(axes.apply_cmap)s
    edgefix : bool, optional
        Whether to fix the the `white-lines-between-filled-contours \
<https://stackoverflow.com/q/8263769/4970632>`__
        and `white-lines-between-pcolor-rectangles \
<https://stackoverflow.com/q/27092991/4970632>`__
        issues. This slows down figure rendering by a bit. Default is
        :rc:`image.edgefix`.
    labels : bool, optional
        For `~matplotlib.axes.Axes.contour`, whether to add contour labels
        with `~matplotlib.axes.Axes.clabel`. For `~matplotlib.axes.Axes.pcolor`
        or `~matplotlib.axes.Axes.pcolormesh`, whether to add labels to the
        center of grid boxes. In the latter case, the text will be black
        when the luminance of the underlying grid box color is >50%%, and
        white otherwise.
    labels_kw : dict-like, optional
        Ignored if `labels` is ``False``. Extra keyword args for the labels.
        For `~matplotlib.axes.Axes.contour`, passed to
        `~matplotlib.axes.Axes.clabel`.  For `~matplotlib.axes.Axes.pcolor`
        or `~matplotlib.axes.Axes.pcolormesh`, passed to
        `~matplotlib.axes.Axes.text`.
    fmt : format-spec, optional
        Passed to the `~proplot.constructor.Norm` constructor, used to format
        number labels. You can also use the `precision` keyword arg.
    precision : int, optional
        Maximum number of decimal places for the number labels.
        Number labels are generated with the
        `~proplot.ticker.SimpleFormatter` formatter, which allows us to
        limit the precision.
    colorbar : bool, int, or str, optional
        If not ``None``, this is a location specifying where to draw an *inset*
        or *panel* colorbar from the resulting mappable. If ``True``, the
        default location is used. Valid locations are described in
        `~proplot.axes.Axes.colorbar`.
    colorbar_kw : dict-like, optional
        Ignored if `colorbar` is ``None``. Extra keyword args for our call
        to `~proplot.axes.Axes.colorbar`.

    Other parameters
    ----------------
    lw, linewidth, linewidths
        The width of `~matplotlib.axes.Axes.contour` lines and
        `~proplot.axes.Axes.parametric` lines. Also the width of lines
        *between* `~matplotlib.axes.Axes.pcolor` boxes,
        `~matplotlib.axes.Axes.pcolormesh` boxes, and
        `~matplotlib.axes.Axes.contourf` filled contours.
    ls, linestyle, linestyles
        As above, but for the line style.
    color, colors, edgecolor, edgecolors
        As above, but for the line color. For `~matplotlib.axes.Axes.contourf`
        plots, if you provide `colors` without specifying the `linewidths`
        or `linestyles`, this argument is used to manually specify the *fill
        colors*. See the `~matplotlib.axes.Axes.contourf` documentation for
        details.
    *args, **kwargs
        Passed to the matplotlib plotting method.

    See also
    --------
    standardize_2d
    proplot.constructor.Colormap
    proplot.constructor.Norm
    proplot.colors.DiscreteNorm
    """
    name = func.__name__
    contour = name in ('contour', 'tricontour')
    contourf = name in ('contourf', 'tricontourf')
    pcolor = name in ('pcolor', 'pcolormesh', 'pcolorfast')
    hexbin = name in ('hexbin',)
    parametric = name in ('parametric',)
    no_discrete_norm = hexbin  # TODO: this should be global setting!!!
    if not args:
        return func(self, *args, **kwargs)
    sample = args[-1]  # used for labels

    # Parse keyword args
    # NOTE: For now when drawing contour or contourf plots with no colormap,
    # cannot use 'values' to specify level centers or level center count.
    # NOTE: For now need to duplicate 'levels' parsing here and in
    # _build_discrete_norm so that it works with contour plots with no cmap.
    cmap_kw = cmap_kw or {}
    norm_kw = norm_kw or {}
    labels_kw = labels_kw or {}
    locator_kw = locator_kw or {}
    colorbar_kw = colorbar_kw or {}
    norm_kw = norm_kw or {}
    edgefix = _not_none(edgefix, rc['image.edgefix'])
    linewidths = _not_none(lw=lw, linewidth=linewidth, linewidths=linewidths)
    linestyles = _not_none(ls=ls, linestyle=linestyle, linestyles=linestyles)
    colors = _not_none(
        color=color, colors=colors, edgecolor=edgecolor, edgecolors=edgecolors,
    )
    levels = _not_none(
        N=N, levels=levels, norm_kw_levels=norm_kw.pop('levels', None), default=rc['image.levels']  # noqa: E501
    )

    # Get colormap, but do not use cmap when 'colors' are passed to contour()
    # or to contourf() -- the latter only when 'linewidths' and 'linestyles'
    # are also *not* passed. This wrapper lets us add "edges" to contourf
    # plots by calling contour() after contourf() if 'linewidths' or
    # 'linestyles' are explicitly passed, but do not want to disable the
    # native matplotlib feature for manually coloring filled contours.
    # https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.contourf
    add_contours = contourf and (linewidths is not None or linestyles is not None)
    use_cmap = colors is None or (not contour and (not contourf or add_contours))
    if not use_cmap:
        if cmap is not None:
            warnings._warn_proplot(
                f'Ignoring input colormap cmap={cmap!r}, using input colors '
                f'colors={colors!r} instead.'
            )
            cmap = None
        if contourf:
            kwargs['colors'] = colors  # this was not done above
            colors = None
    else:
        if cmap is None:
            if name == 'spy':
                cmap = pcolors.ListedColormap(['w', 'k'], '_binary')
            else:
                cmap = rc['image.cmap']
        cmap = constructor.Colormap(cmap, **cmap_kw)
        if getattr(cmap, '_cyclic', None) and extend != 'neither':
            warnings._warn_proplot(
                'Cyclic colormap requires extend="neither". '
                f'Overriding user input extend={extend!r}.'
            )
            extend = 'neither'

    # Translate standardized keyword arguments back into the keyword args
    # accepted by native matplotlib methods. Also disable edgefix if user want
    # to customize the "edges".
    style_kw = STYLE_ARGS_TRANSLATE.get(name, None)
    for key, value in (
        ('colors', colors),
        ('linewidths', linewidths),
        ('linestyles', linestyles)
    ):
        if value is None or add_contours:
            continue
        if not style_kw or key not in style_kw:  # no known conversion table
            raise TypeError(f'{name}() got an unexpected keyword argument {key!r}')
        edgefix = False  # disable edgefix when specifying borders!
        kwargs[style_kw[key]] = value

    # Build colormap normalizer and update keyword args
    # NOTE: Standard algorithm for obtaining default levels does not work
    # for hexbin, because it colors *counts*, not data values!
    ticks = None
    if not use_cmap and not np.iterable(levels):
        levels, _ = _auto_levels_locator(
            sample, N=levels,
            norm=norm, norm_kw=norm_kw, locator=locator, locator_kw=locator_kw,
            vmin=vmin, vmax=vmax, extend=extend, symmetric=symmetric,
        )
    if use_cmap and not no_discrete_norm:
        norm, cmap, levels, ticks = _build_discrete_norm(
            sample, levels=levels, values=values,
            cmap=cmap, norm=norm, norm_kw=norm_kw, vmin=vmin, vmax=vmax, extend=extend,
            symmetric=symmetric, positive=positive, negative=negative, nozero=nozero,
            locator=locator, locator_kw=locator_kw,
            minlength=(2 - int(contour)),  # allow single contours
        )
    if nozero and np.iterable(levels) and 0 in levels:
        levels = _to_ndarray(levels)
        levels = levels[levels != 0]
    if cmap is not None:
        kwargs['cmap'] = cmap
    if norm is not None:
        kwargs['norm'] = norm
    if no_discrete_norm:
        kwargs['vmin'] = vmin
        kwargs['vmax'] = vmax
    if contour or contourf:
        kwargs['levels'] = levels
        kwargs['extend'] = extend
    if parametric:
        kwargs['values'] = values

    # Call function and possibly add solid contours between filled ones or
    # fix common "white lines" issues with vector graphic output
    obj = func(self, *args, **kwargs)
    if not isinstance(obj, tuple):  # hist2d
        obj._colorbar_extend = extend  # used by proplot colorbar
        obj._colorbar_ticks = ticks  # used by proplot colorbar
    if add_contours:
        colors = _not_none(colors, 'k')
        self.contour(
            *args, levels=levels, linewidths=linewidths,
            linestyles=linestyles, colors=colors
        )
    if edgefix:
        _fix_white_lines(obj)

    # Apply labels
    # TODO: Add quiverkey to this!
    if labels:
        fmt = _not_none(labels_kw.pop('fmt', None), fmt, 'simple')
        fmt = constructor.Formatter(fmt, precision=precision)
        if contour or contourf:
            _labels_contour(self, obj, *args, fmt=fmt, **labels_kw)
        elif pcolor:
            _labels_pcolor(self, obj, fmt=fmt, **labels_kw)
        else:
            raise RuntimeError(f'Not possible to add labels to {name!r} plot.')

    # Optionally add colorbar
    if colorbar:
        loc = self._loc_translate(colorbar, 'colorbar')
        if parametric and values is not None:
            colorbar_kw.setdefault('values', values)
        if loc != 'fill':
            colorbar_kw.setdefault('loc', loc)
        self.colorbar(obj, **colorbar_kw)

    return obj


def _generate_mappable(
    mappable, values=None, *, orientation='horizontal',
    locator=None, formatter=None, norm=None, norm_kw=None, rotation=None,
):
    """
    Generate a mappable from flexible non-mappable input. Useful in bridging
    the gap between legends and colorbars (e.g., creating colorbars from line
    objects whose data values span a natural colormap range).
    """
    # A colormap instance
    # TODO: Pass remaining arguments through Colormap()? This is really
    # niche usage so maybe not necessary.
    if isinstance(mappable, mcolors.Colormap):
        # NOTE: 'Values' makes no sense if this is just a colormap. Just
        # use unique color for every segmentdata / colors color.
        cmap = mappable
        values = np.linspace(0, 1, cmap.N)

    # List of colors
    elif np.iterable(mappable) and all(
        isinstance(obj, str) or (np.iterable(obj) and len(obj) in (3, 4))
        for obj in mappable
    ):
        cmap = mcolors.ListedColormap(list(mappable), '_no_name')
        if values is None:
            values = np.arange(len(mappable))
        locator = _not_none(locator, values)  # tick *all* values by default

    # List of artists
    # NOTE: Do not check for isinstance(Artist) in case it is an mpl collection
    elif np.iterable(mappable) and all(
        hasattr(obj, 'get_color') or hasattr(obj, 'get_facecolor')
        for obj in mappable
    ):
        # Generate colormap from colors and infer tick labels
        colors = []
        for obj in mappable:
            if hasattr(obj, 'get_color'):
                color = obj.get_color()
            else:
                color = obj.get_facecolor()
            if isinstance(color, np.ndarray):
                color = color.squeeze()  # e.g. scatter plot
                if color.ndim != 1:
                    raise ValueError(
                        'Cannot make colorbar from list of artists '
                        f'with more than one color: {color!r}.'
                    )
            colors.append(to_rgb(color))

        # Try to infer tick values and tick labels from Artist labels
        cmap = mcolors.ListedColormap(colors, '_no_name')
        if values is None:
            # Get object labels and values
            labels = []
            values = []
            for obj in mappable:
                label = _get_label(obj)  # could be None
                try:
                    value = float(label)  # could be float(None)
                except (TypeError, ValueError):
                    value = None
                labels.append(label)
                values.append(value)

            # Use default values if labels are non-numeric (numeric labels are
            # common when making on-the-fly colorbars). Try to use object labels
            # for ticks with default vertical rotation, like datetime axes.
            if any(value is None for value in values):
                values = np.arange(len(mappable))
                if formatter is None and any(label is not None for label in labels):
                    formatter = labels  # use these fixed values for ticks
                    if orientation == 'horizontal':
                        rotation = _not_none(rotation, 90)
        locator = _not_none(locator, values)  # tick *all* values by default

    else:
        raise ValueError(
            'Input mappable must be a matplotlib artist, '
            'list of objects, list of colors, or colormap. '
            f'Got {mappable!r}.'
        )

    # Build ad hoc ScalarMappable object from colors
    if np.iterable(mappable) and len(values) != len(mappable):
        raise ValueError(
            f'Passed {len(values)} values, but only {len(mappable)} '
            f'objects or colors.'
        )
    norm, *_ = _build_discrete_norm(
        values=values,
        cmap=cmap,
        norm=norm,
        norm_kw=norm_kw,
        extend='neither',
    )
    mappable = mcm.ScalarMappable(norm, cmap)

    return mappable, rotation


def colorbar_extras(
    self, mappable, values=None, *,  # analogous to handles and labels
    extend=None, extendsize=None,
    title=None, label=None,
    grid=None, tickminor=None,
    reverse=False, tickloc=None, ticklocation=None, tickdir=None, tickdirection=None,
    locator=None, ticks=None, maxn=None, maxn_minor=None,
    minorlocator=None, minorticks=None,
    locator_kw=None, minorlocator_kw=None,
    formatter=None, ticklabels=None, formatter_kw=None, rotation=None,
    norm=None, norm_kw=None,  # normalizer to use when passing colors/lines
    orientation='horizontal',
    edgecolor=None, linewidth=None,
    labelsize=None, labelweight=None, labelcolor=None,
    ticklabelsize=None, ticklabelweight=None, ticklabelcolor=None,
    **kwargs
):
    """
    Adds useful features for controlling colorbars.

    Important
    ---------
    This function wraps `proplot.axes.Axes.colorbar`
    and `proplot.figure.Figure.colorbar`.

    Parameters
    ----------
    mappable : mappable, list of plot handles, list of color-spec, \
or colormap-spec
        There are four options here:

        1. A mappable object. Basically, any object with a ``get_cmap`` method,
           like the objects returned by `~matplotlib.axes.Axes.contourf` and
           `~matplotlib.axes.Axes.pcolormesh`.
        2. A list of "plot handles". Basically, any object with a ``get_color``
           method, like `~matplotlib.lines.Line2D` instances. A colormap will
           be generated from the colors of these objects, and colorbar levels
           will be selected using `values`.  If `values` is ``None``, we try
           to infer them by converting the handle labels returned by
           `~matplotlib.artist.Artist.get_label` to `float`. Otherwise, it is
           set to ``np.linspace(0, 1, len(mappable))``.
        3. A list of hex strings, color string names, or RGB tuples. A colormap
           will be generated from these colors, and colorbar levels will be
           selected using `values`. If `values` is ``None``, it is set to
           ``np.linspace(0, 1, len(mappable))``.
        4. A `~matplotlib.colors.Colormap` instance. In this case, a colorbar
           will be drawn using this colormap and with levels determined by
           `values`. If `values` is ``None``, it is set to
           ``np.linspace(0, 1, cmap.N)``.

    values : list of float, optional
        Ignored if `mappable` is a mappable object. This maps each color or
        plot handle in the `mappable` list to numeric values, from which a
        colormap and normalizer are constructed.
    norm : normalizer spec, optional
        Ignored if `values` is ``None``. The normalizer for converting `values`
        to colormap colors. Passed to `~proplot.constructor.Norm`.
    norm_kw : dict-like, optional
        The normalizer settings. Passed to `~proplot.constructor.Norm`.
    extend : {None, 'neither', 'both', 'min', 'max'}, optional
        Direction for drawing colorbar "extensions" (i.e. references to
        out-of-bounds data with a unique color). These are triangles by
        default. If ``None``, we try to use the ``extend`` attribute on the
        mappable object. If the attribute is unavailable, we use ``'neither'``.
    extendsize : float or str, optional
        The length of the colorbar "extensions" in *physical units*.
        If float, units are inches. If string, units are interpreted
        by `~proplot.utils.units`. Default is :rc:`colorbar.insetextend`
        for inset colorbars and :rc:`colorbar.extend` for outer colorbars.
    reverse : bool, optional
        Whether to reverse the direction of the colorbar.
    tickloc, ticklocation : {'bottom', 'top', 'left', 'right'}, optional
        Where to draw tick marks on the colorbar.
    tickdir, tickdirection : {'out', 'in', 'inout'}, optional
        Direction that major and minor tick marks point.
    tickminor : bool, optional
        Whether to add minor ticks to the colorbar with
        `~matplotlib.colorbar.ColorbarBase.minorticks_on`.
    grid : bool, optional
        Whether to draw "gridlines" between each level of the colorbar.
        Default is :rc:`colorbar.grid`.
    label, title : str, optional
        The colorbar label. The `title` keyword is also accepted for
        consistency with `legend`.
    locator, ticks : locator spec, optional
        Used to determine the colorbar tick positions. Passed to the
        `~proplot.constructor.Locator` constructor.
    maxn : int, optional
        Used if `locator` is ``None``. Determines the maximum number of levels
        that are ticked. Default depends on the colorbar length relative
        to the font size. The keyword name "maxn" is meant to mimic
        the `~matplotlib.ticker.MaxNLocator` class name.
    locator_kw : dict-like, optional
        The locator settings. Passed to `~proplot.constructor.Locator`.
    minorlocator, minorticks, maxn_minor, minorlocator_kw
        As with `locator`, `maxn`, and `locator_kw`, but for the minor ticks.
    formatter, ticklabels : formatter spec, optional
        The tick label format. Passed to the `~proplot.constructor.Formatter`
        constructor.
    formatter_kw : dict-like, optional
        The formatter settings. Passed to `~proplot.constructor.Formatter`.
    rotation : float, optional
        The tick label rotation. Default is ``0``.
    edgecolor, linewidth : optional
        The edge color and line width for the colorbar outline.
    labelsize, labelweight, labelcolor : optional
        The font size, weight, and color for colorbar label text.
    ticklabelsize, ticklabelweight, ticklabelcolor : optional
        The font size, weight, and color for colorbar tick labels.
    orientation : {{None, 'horizontal', 'vertical'}}, optional
        The colorbar orientation. By default this depends on the "side" of the subplot
        or figure where the colorbar is drawn. Inset colorbars are always horizontal.

    Other parameters
    ----------------
    **kwargs
        Passed to `matplotlib.figure.Figure.colorbar`.

    See also
    --------
    proplot.axes.Axes.colorbar
    proplot.figure.Figure.colorbar
    """
    # NOTE: There is a weird problem with colorbars when simultaneously
    # passing levels and norm object to a mappable; fixed by passing vmin/vmax
    # instead of levels. (see: https://stackoverflow.com/q/40116968/4970632).
    # NOTE: Often want levels instead of vmin/vmax, while simultaneously
    # using a Normalize (for example) to determine colors between the levels
    # (see: https://stackoverflow.com/q/42723538/4970632). Workaround makes
    # sure locators are in vmin/vmax range exclusively; cannot match values.
    # NOTE: In legend_extras() we try to add to the objects accepted by
    # legend() using handler_map. We can't really do anything similar for
    # colorbars; input must just be insnace of mixin class cm.ScalarMappable
    # Mutable args
    norm_kw = norm_kw or {}
    formatter_kw = formatter_kw or {}
    locator_kw = locator_kw or {}
    minorlocator_kw = minorlocator_kw or {}

    # Parse input args
    label = _not_none(title=title, label=label)
    locator = _not_none(ticks=ticks, locator=locator)
    minorlocator = _not_none(minorticks=minorticks, minorlocator=minorlocator)
    ticklocation = _not_none(tickloc=tickloc, ticklocation=ticklocation)
    tickdirection = _not_none(tickdir=tickdir, tickdirection=tickdirection)
    formatter = _not_none(ticklabels=ticklabels, formatter=formatter)

    # Colorbar kwargs
    grid = _not_none(grid, rc['colorbar.grid'])
    kwargs.update({
        'cax': self,
        'use_gridspec': True,
        'orientation': orientation,
        'spacing': 'uniform',
    })
    kwargs.setdefault('drawedges', grid)

    # Special case where auto colorbar is generated from 1d methods, a list is
    # always passed, but some 1d methods (scatter) do have colormaps.
    if (
        np.iterable(mappable)
        and len(mappable) == 1
        and hasattr(mappable[0], 'get_cmap')
    ):
        mappable = mappable[0]

    # For container objects, we just assume color is the same for every item.
    # Works for ErrorbarContainer, StemContainer, BarContainer.
    if (
        np.iterable(mappable)
        and len(mappable) > 0
        and all(isinstance(obj, mcontainer.Container) for obj in mappable)
    ):
        mappable = [obj[0] for obj in mappable]

    # Test if we were given a mappable, or iterable of stuff; note Container
    # and PolyCollection matplotlib classes are iterable.
    if not isinstance(mappable, (martist.Artist, mcontour.ContourSet)):
        mappable, rotation = _generate_mappable(
            mappable, values, locator=locator, formatter=formatter,
            norm=norm, norm_kw=norm_kw, orientation=orientation, rotation=rotation
        )

    # Define text property keyword args
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

    # Try to get tick locations from *levels* or from *values* rather than
    # random points along the axis.
    # NOTE: Do not necessarily want e.g. minor tick locations at logminor for LogNorm!
    # In _build_discrete_norm we sometimes select evenly spaced levels in log-space
    # *between* powers of 10, so logminor ticks would be misaligned with levels.
    if locator is None:
        locator = getattr(mappable, '_colorbar_ticks', None)
        if locator is None:
            # This should only happen if user calls plotting method on native
            # matplotlib axes.
            if isinstance(norm, mcolors.LogNorm):
                locator = 'log'
            elif isinstance(norm, mcolors.SymLogNorm):
                locator = 'symlog'
                locator_kw.setdefault('linthresh', norm.linthresh)
            else:
                locator = 'auto'

        elif not isinstance(locator, mticker.Locator):
            # Get default maxn, try to allot 2em squares per label maybe?
            # NOTE: Cannot use Axes.get_size_inches because this is a
            # native matplotlib axes
            width, height = self.figure.get_size_inches()
            if orientation == 'horizontal':
                scale = 3  # em squares alotted for labels
                length = width * abs(self.get_position().width)
                fontsize = kw_ticklabels.get('size', rc['xtick.labelsize'])
            else:
                scale = 1
                length = height * abs(self.get_position().height)
                fontsize = kw_ticklabels.get('size', rc['ytick.labelsize'])
            fontsize = rc._scale_font(fontsize)
            maxn = _not_none(maxn, int(length / (scale * fontsize / 72)))
            maxn_minor = _not_none(maxn_minor, int(length / (0.5 * fontsize / 72)))

            # Get locator
            if tickminor and minorlocator is None:
                step = 1 + len(locator) // max(1, maxn_minor)
                minorlocator = locator[::step]
            step = 1 + len(locator) // max(1, maxn)
            locator = locator[::step]

    # Get extend triangles in physical units
    width, height = self.figure.get_size_inches()
    if orientation == 'horizontal':
        scale = width * abs(self.get_position().width)
    else:
        scale = height * abs(self.get_position().height)
    extendsize = units(_not_none(extendsize, rc['colorbar.extend']))
    extendsize = extendsize / (scale - 2 * extendsize)

    # Draw the colorbar
    # NOTE: Set default formatter here because we optionally apply a FixedFormatter
    # using *labels* from handle input.
    extend = _not_none(extend, getattr(mappable, '_colorbar_extend', 'neither'))
    locator = constructor.Locator(locator, **locator_kw)
    formatter = constructor.Formatter(_not_none(formatter, 'auto'), **formatter_kw)
    kwargs.update({
        'ticks': locator,
        'format': formatter,
        'ticklocation': ticklocation,
        'extendfrac': extendsize,
    })
    if isinstance(mappable, mcontour.ContourSet):
        mappable.extend = extend  # required in mpl >= 3.3, else optional
    else:
        kwargs['extend'] = extend
    cb = self.figure.colorbar(mappable, **kwargs)
    axis = self.xaxis if orientation == 'horizontal' else self.yaxis

    # The minor locator
    # TODO: Document the improved minor locator functionality!
    # NOTE: Colorbar._use_auto_colorbar_locator() is never True because we use
    # the custom DiscreteNorm normalizer. Colorbar._ticks() always called.
    if minorlocator is None:
        if tickminor:
            cb.minorticks_on()
        else:
            cb.minorticks_off()
    elif not hasattr(cb, '_ticker'):
        warnings._warn_proplot(
            'Matplotlib colorbar API has changed. '
            f'Cannot use custom minor tick locator {minorlocator!r}.'
        )
        cb.minorticks_on()  # at least turn them on
    else:
        # Set the minor ticks just like matplotlib internally sets the
        # major ticks. Private API is the only way!
        minorlocator = constructor.Locator(minorlocator, **minorlocator_kw)
        ticks, *_ = cb._ticker(minorlocator, mticker.NullFormatter())
        axis.set_ticks(ticks, minor=True)
        axis.set_ticklabels([], minor=True)

    # Label and tick label settings
    # WARNING: Must use colorbar set_label to set text, calling set_text on
    # the axis will do nothing!
    if label is not None:
        cb.set_label(label)
    axis.label.update(kw_label)
    for obj in axis.get_ticklabels():
        obj.update(kw_ticklabels)

    # Ticks consistent with rc settings and overrides
    s = axis.axis_name
    for which in ('minor', 'major'):
        kw = rc.category(s + 'tick.' + which)
        kw.pop('visible', None)
        if tickdirection:
            kw['direction'] = tickdirection
        if edgecolor:
            kw['color'] = edgecolor
        if linewidth:
            kw['width'] = linewidth
        axis.set_tick_params(which=which, **kw)
    axis.set_ticks_position(ticklocation)

    # Fix alpha-blending issues. Cannot set edgecolor to 'face' because blending will
    # occur, end up with colored lines instead of white ones. Need manual blending.
    # NOTE: For some reason cb solids uses listed colormap with always 1.0 alpha,
    # then alpha is applied after. See: https://stackoverflow.com/a/35672224/4970632
    cmap = cb.cmap
    blend = 'pcolormesh.snap' not in rc or not rc['pcolormesh.snap']
    if not cmap._isinit:
        cmap._init()
    if blend and any(cmap._lut[:-1, 3] < 1):
        warnings._warn_proplot(
            f'Using manual alpha-blending for {cmap.name!r} colorbar solids.'
        )
        # Generate "secret" copy of the colormap!
        lut = cmap._lut.copy()
        cmap = mcolors.Colormap('_cbar_fix', N=cmap.N)
        cmap._isinit = True
        cmap._init = lambda: None
        # Manually fill lookup table with alpha-blended RGB colors!
        for i in range(lut.shape[0] - 1):
            alpha = lut[i, 3]
            lut[i, :3] = (1 - alpha) * 1 + alpha * lut[i, :3]  # blend *white*
            lut[i, 3] = 1
        cmap._lut = lut
        # Update colorbar
        cb.cmap = cmap
        cb.draw_all()

    # Fix colorbar outline
    kw_outline = {
        'edgecolor': _not_none(edgecolor, rc['axes.edgecolor']),
        'linewidth': _not_none(linewidth, rc['axes.linewidth']),
    }
    if cb.outline is not None:
        cb.outline.update(kw_outline)
    if cb.dividers is not None:
        cb.dividers.update(kw_outline)

    # *Never* rasterize because it causes misalignment with border lines
    if cb.solids:
        cb.solids.set_rasterized(False)
        cb.solids.set_linewidth(0.4)
        cb.solids.set_edgecolor('face')

    # Invert the axis if descending DiscreteNorm
    norm = mappable.norm
    if getattr(norm, '_descending', None):
        axis.set_inverted(True)
    if reverse:  # potentially double reverse, although that would be weird...
        axis.set_inverted(True)
    return cb


def _iter_legend_children(children):
    """
    Iterate recursively through `_children` attributes of various `HPacker`,
    `VPacker`, and `DrawingArea` classes.
    """
    for obj in children:
        if hasattr(obj, '_children'):
            yield from _iter_legend_children(obj._children)
        else:
            yield obj


def _multiple_legend(self, pairs, loc=None, ncol=None, order=None, **kwargs):
    """
    Draw "legend" with centered rows by creating separate legends for
    each row. The label spacing/border spacing will be exactly replicated.
    """
    # Message when overriding some properties
    legs = []
    overridden = []
    frameon = kwargs.pop('frameon', None)  # then add back later!
    for override in ('bbox_transform', 'bbox_to_anchor'):
        prop = kwargs.pop(override, None)
        if prop is not None:
            overridden.append(override)
    if ncol is not None:
        warnings._warn_proplot(
            'Detected list of *lists* of legend handles. '
            'Ignoring user input property "ncol".'
        )
    if overridden:
        warnings._warn_proplot(
            'Ignoring user input properties '
            + ', '.join(map(repr, overridden))
            + ' for centered-row legend.'
        )

    # Determine space we want sub-legend to occupy as fraction of height
    # NOTE: Empirical testing shows spacing fudge factor necessary to
    # exactly replicate the spacing of standard aligned legends.
    width, height = self.get_size_inches()
    fontsize = kwargs.get('fontsize', None) or rc['legend.fontsize']
    fontsize = rc._scale_font(fontsize)
    spacing = kwargs.get('labelspacing', None) or rc['legend.labelspacing']
    if pairs:
        interval = 1 / len(pairs)  # split up axes
        interval = (((1 + spacing * 0.85) * fontsize) / 72) / height

    # Iterate and draw
    # NOTE: We confine possible bounding box in *y*-direction, but do not
    # confine it in *x*-direction. Matplotlib will automatically move
    # left-to-right if you request this.
    ymin = ymax = None
    if order == 'F':
        raise NotImplementedError(
            'When center=True, ProPlot vertically stacks successive '
            "single-row legends. Column-major (order='F') ordering "
            'is un-supported.'
        )
    loc = _not_none(loc, 'upper center')
    if not isinstance(loc, str):
        raise ValueError(
            f'Invalid location {loc!r} for legend with center=True. '
            'Must be a location *string*.'
        )
    elif loc == 'best':
        warnings._warn_proplot(
            'For centered-row legends, cannot use "best" location. '
            'Using "upper center" instead.'
        )

    # Iterate through sublists
    for i, ipairs in enumerate(pairs):
        if i == 1:
            title = kwargs.pop('title', None)
        if i >= 1 and title is not None:
            i += 1  # add extra space!

        # Legend position
        if 'upper' in loc:
            y1 = 1 - (i + 1) * interval
            y2 = 1 - i * interval
        elif 'lower' in loc:
            y1 = (len(pairs) + i - 2) * interval
            y2 = (len(pairs) + i - 1) * interval
        else:  # center
            y1 = 0.5 + interval * len(pairs) / 2 - (i + 1) * interval
            y2 = 0.5 + interval * len(pairs) / 2 - i * interval
        ymin = min(y1, _not_none(ymin, y1))
        ymax = max(y2, _not_none(ymax, y2))

        # Draw legend
        bbox = mtransforms.Bbox([[0, y1], [1, y2]])
        leg = mlegend.Legend(
            self, *zip(*ipairs), loc=loc, ncol=len(ipairs),
            bbox_transform=self.transAxes, bbox_to_anchor=bbox,
            frameon=False, **kwargs
        )
        legs.append(leg)

    # Simple cases
    if not frameon:
        return legs
    if len(legs) == 1:
        legs[0].set_frame_on(True)
        return legs

    # Draw manual fancy bounding box for un-aligned legend
    # WARNING: The matplotlib legendPatch transform is the default transform, i.e.
    # universal coordinates in points. Means we have to transform mutation scale
    # into transAxes sizes.
    # WARNING: Tempting to use legendPatch for everything but for some reason
    # coordinates are messed up. In some tests all coordinates were just result
    # of get window extent multiplied by 2 (???). Anyway actual box is found in
    # _legend_box attribute, which is accessed by get_window_extent.
    width, height = self.get_size_inches()
    renderer = self.figure._get_renderer()
    bboxs = [
        leg.get_window_extent(renderer).transformed(self.transAxes.inverted())
        for leg in legs
    ]
    xmin = min(bbox.xmin for bbox in bboxs)
    xmax = max(bbox.xmax for bbox in bboxs)
    ymin = min(bbox.ymin for bbox in bboxs)
    ymax = max(bbox.ymax for bbox in bboxs)
    fontsize = (fontsize / 72) / width  # axes relative units
    fontsize = renderer.points_to_pixels(fontsize)

    # Draw and format patch
    patch = mpatches.FancyBboxPatch(
        (xmin, ymin), xmax - xmin, ymax - ymin,
        snap=True, zorder=4.5,
        mutation_scale=fontsize,
        transform=self.transAxes
    )
    if kwargs.get('fancybox', rc['legend.fancybox']):
        patch.set_boxstyle('round', pad=0, rounding_size=0.2)
    else:
        patch.set_boxstyle('square', pad=0)
    patch.set_clip_on(False)
    self.add_artist(patch)

    # Add shadow
    # TODO: This does not work, figure out
    if kwargs.get('shadow', rc['legend.shadow']):
        shadow = mpatches.Shadow(patch, 20, -20)
        self.add_artist(shadow)

    # Add patch to list
    return patch, *legs


def _single_legend(self, pairs, ncol=None, order=None, **kwargs):
    """
    Draw an individual legend with support for changing legend-entries
    between column-major and row-major.
    """
    # Optionally change order
    # See: https://stackoverflow.com/q/10101141/4970632
    # Example: If 5 columns, but final row length 3, columns 0-2 have
    # N rows but 3-4 have N-1 rows.
    ncol = _not_none(ncol, 3)
    if order == 'C':
        split = [pairs[i * ncol:(i + 1) * ncol] for i in range(len(pairs) // ncol + 1)]
        pairs = []
        nrows_max = len(split)  # max possible row count
        ncols_final = len(split[-1])  # columns in final row
        nrows = [nrows_max] * ncols_final + [nrows_max - 1] * (ncol - ncols_final)
        for col, nrow in enumerate(nrows):  # iterate through cols
            pairs.extend(split[row][col] for row in range(nrow))

    # Draw legend
    return mlegend.Legend(self, *zip(*pairs), ncol=ncol, **kwargs)


def legend_extras(
    self, handles=None, labels=None, *, loc=None, ncol=None, ncols=None,
    center=None, order='C', label=None, title=None,
    fontsize=None, fontweight=None, fontcolor=None,
    color=None, marker=None, lw=None, linewidth=None,
    dashes=None, linestyle=None, markersize=None, frameon=None, frame=None,
    **kwargs
):
    """
    Adds useful features for controlling legends, including "centered-row"
    legends.

    Important
    ---------
    This function wraps `proplot.axes.Axes.legend`
    and `proplot.figure.Figure.legend`.

    Parameters
    ----------
    handles : list of `~matplotlib.artist.Artist`, optional
        List of artists instances, or list of lists of artist instances (see
        the `center` keyword). If ``None``, the artists are retrieved with
        `~matplotlib.axes.Axes.get_legend_handles_labels`.
    labels : list of str, optional
        Matching list of string labels, or list of lists of string labels (see
        the `center` keywod). If ``None``, the labels are retrieved by calling
        `~matplotlib.artist.Artist.get_label` on each
        `~matplotlib.artist.Artist` in `handles`.
    loc : int or str, optional
        The legend location. The following location keys are valid:

        ==================  ================================================
        Location            Valid keys
        ==================  ================================================
        "best" possible     ``0``, ``'best'``, ``'b'``, ``'i'``, ``'inset'``
        upper right         ``1``, ``'upper right'``, ``'ur'``
        upper left          ``2``, ``'upper left'``, ``'ul'``
        lower left          ``3``, ``'lower left'``, ``'ll'``
        lower right         ``4``, ``'lower right'``, ``'lr'``
        center left         ``5``, ``'center left'``, ``'cl'``
        center right        ``6``, ``'center right'``, ``'cr'``
        lower center        ``7``, ``'lower center'``, ``'lc'``
        upper center        ``8``, ``'upper center'``, ``'uc'``
        center              ``9``, ``'center'``, ``'c'``
        ==================  ================================================

    ncol, ncols : int, optional
        The number of columns. `ncols` is an alias, added
        for consistency with `~matplotlib.pyplot.subplots`.
    order : {'C', 'F'}, optional
        Whether legend handles are drawn in row-major (``'C'``) or column-major
        (``'F'``) order. Analagous to `numpy.array` ordering. For some reason
        ``'F'`` was the original matplotlib default. Default is ``'C'``.
    center : bool, optional
        Whether to center each legend row individually. If ``True``, we
        actually draw successive single-row legends stacked on top of each
        other. If ``None``, we infer this setting from `handles`. Default is
        ``True`` if `handles` is a list of lists (each sublist is used as a *row*
        in the legend). Otherwise, default is ``False``.
    title, label : str, optional
        The legend title. The `label` keyword is also accepted, for consistency
        with `colorbar`.
    fontsize, fontweight, fontcolor : optional
        The font size, weight, and color for legend text.
    color, lw, linewidth, marker, linestyle, dashes, markersize : \
property-spec, optional
        Properties used to override the legend handles. For example, for a
        legend describing variations in line style ignoring variations in color, you
        might want to use ``color='k'``. For now this does not include `facecolor`,
        `edgecolor`, and `alpha`, because `~matplotlib.axes.Axes.legend` uses these
        keyword args to modify the frame properties.

    Other parameters
    ----------------
    **kwargs
        Passed to `matplotlib.axes.Axes.legend`.

    See also
    --------
    proplot.axes.Axes.legend
    proplot.figure.Figure.legend
    """
    # Parse input args
    # TODO: Legend entries for colormap or scatterplot objects! Idea is we
    # pass a scatter plot or contourf or whatever, and legend is generated by
    # drawing patch rectangles or markers using data values and their
    # corresponding cmap colors! For scatterplots just test get_facecolor()
    # to see if it contains more than one color.
    # TODO: It is *also* often desirable to label a colormap object with
    # one data value. Maybe add a legend option for the *number of samples*
    # or the *sample points* when drawing legends for colormap objects.
    # Look into "legend handlers", might just want to add own handlers by
    # passing handler_map to legend() and get_legend_handles_labels().
    if order not in ('F', 'C'):
        raise ValueError(
            f'Invalid order {order!r}. Choose from '
            '"C" (row-major, default) and "F" (column-major).'
        )
    ncol = _not_none(ncols=ncols, ncol=ncol)
    title = _not_none(label=label, title=title)
    frameon = _not_none(frame=frame, frameon=frameon, default=rc['legend.frameon'])
    if handles is not None and not np.iterable(handles):  # e.g. a mappable object
        handles = [handles]
    if labels is not None and (not np.iterable(labels) or isinstance(labels, str)):
        labels = [labels]
    if title is not None:
        kwargs['title'] = title
    if frameon is not None:
        kwargs['frameon'] = frameon
    if fontsize is not None:
        kwargs['fontsize'] = rc._scale_font(fontsize)

    # Handle and text properties that are applied after-the-fact
    # NOTE: Set solid_capstyle to 'butt' so line does not extend past error bounds
    # shading in legend entry. This change is not noticable in other situations.
    kw_text = {}
    for key, value in (('color', fontcolor), ('weight', fontweight)):
        if value is not None:
            kw_text[key] = value
    kw_handle = {'solid_capstyle': 'butt'}
    for key, value in (
        ('color', color),
        ('marker', marker),
        ('linewidth', lw),
        ('linewidth', linewidth),
        ('markersize', markersize),
        ('linestyle', linestyle),
        ('dashes', dashes),
    ):
        if value is not None:
            kw_handle[key] = value

    # Get axes for legend handle detection
    # TODO: Update this when no longer use "filled panels" for outer legends
    axs = [self]
    if self._panel_hidden:
        if self._panel_parent:  # axes panel
            axs = list(self._panel_parent._iter_axes(hidden=False, children=True))
        else:
            axs = list(self.figure._iter_axes(hidden=False, children=True))

    # Handle list of lists (centered row legends)
    # NOTE: Avoid very common plot() error where users draw individual lines
    # with plot() and add singleton tuples to a list of handles. If matplotlib
    # gets a list like this but gets no 'labels' argument, it raises error.
    list_of_lists = False
    if handles is not None:
        handles = [h[0] if isinstance(h, tuple) and len(h) == 1 else h for h in handles]
        list_of_lists = any(isinstance(h, (list, np.ndarray)) for h in handles)
    if list_of_lists:
        if any(not np.iterable(_) for _ in handles):
            raise ValueError(f'Invalid handles={handles!r}.')
        if not labels:
            labels = [None] * len(handles)
        elif not all(np.iterable(_) and not isinstance(_, str) for _ in labels):
            # e.g. handles=[obj1, [obj2, obj3]] requires labels=[lab1, [lab2, lab3]]
            raise ValueError(f'Invalid labels={labels!r} for handles={handles!r}.')

    # Parse handles and legends with native matplotlib parser
    if not list_of_lists:
        if isinstance(handles, np.ndarray):
            handles = handles.tolist()
        if isinstance(labels, np.ndarray):
            labels = labels.tolist()
        handles, labels, *_ = mlegend._parse_legend_args(
            axs, handles=handles, labels=labels,
        )
        pairs = list(zip(handles, labels))
    else:
        pairs = []
        for ihandles, ilabels in zip(handles, labels):
            if isinstance(ihandles, np.ndarray):
                ihandles = ihandles.tolist()
            if isinstance(ilabels, np.ndarray):
                ilabels = ilabels.tolist()
            ihandles, ilabels, *_ = mlegend._parse_legend_args(
                axs, handles=ihandles, labels=ilabels,
            )
            pairs.append(list(zip(ihandles, ilabels)))

    # Manage pairs in context of 'center' option
    center = _not_none(center, list_of_lists)
    if not center and list_of_lists:  # standardize format based on input
        list_of_lists = False  # no longer is list of lists
        pairs = [pair for ipairs in pairs for pair in ipairs]
    elif center and not list_of_lists:
        list_of_lists = True
        ncol = _not_none(ncol, 3)
        pairs = [pairs[i * ncol:(i + 1) * ncol] for i in range(len(pairs))]
        ncol = None
    if list_of_lists:  # remove empty lists, pops up in some examples
        pairs = [ipairs for ipairs in pairs if ipairs]

    # Bail if no pairs
    if not pairs:
        return mlegend.Legend(self, [], [], loc=loc, ncol=ncol, **kwargs)
    # Multiple-legend pseudo-legend
    elif center:
        objs = _multiple_legend(self, pairs, loc=loc, ncol=ncol, order=order, **kwargs)
    # Individual legend
    else:
        objs = [_single_legend(self, pairs, loc=loc, ncol=ncol, order=order, **kwargs)]

    # Add legends manually so matplotlib does not remove old ones
    for obj in objs:
        if isinstance(obj, mpatches.FancyBboxPatch):
            continue
        if hasattr(self, 'legend_') and self.legend_ is None:
            self.legend_ = obj  # set *first* legend accessible with get_legend()
        else:
            self.add_artist(obj)

    # Apply legend box properties
    outline = rc.fill({
        'linewidth': 'axes.linewidth',
        'edgecolor': 'axes.edgecolor',
        'facecolor': 'axes.facecolor',
        'alpha': 'legend.framealpha',
    })
    for key in (*outline,):
        if key != 'linewidth':
            if kwargs.get(key, None):
                outline.pop(key, None)
    for obj in objs:
        if isinstance(obj, mpatches.FancyBboxPatch):
            obj.update(outline)  # the multiple-legend bounding box
        else:
            obj.legendPatch.update(outline)  # no-op if frame is off

    # Apply *overrides* to legend elements
    # WARNING: legendHandles only contains the *first* artist per legend because
    # HandlerBase.legend_artist() called in Legend._init_legend_box() only
    # returns the first artist. Instead we try to iterate through offset boxes.
    # TODO: Remove this feature? Idea was this lets users create *categorical*
    # legends in clunky way, e.g. entries denoting *colors* and entries denoting
    # *markers*. But would be better to add capacity for categorical labels in a
    # *single* legend like seaborn rather than multiple legends.
    for obj in objs:
        try:
            children = obj._legend_handle_box._children
        except AttributeError:  # older versions maybe?
            children = []
        for obj in _iter_legend_children(children):
            # Account for mixed legends, e.g. line on top of error bounds shading
            if isinstance(obj, mtext.Text):
                obj.update(kw_text)
            else:
                for key, value in kw_handle.items():
                    getattr(obj, f'set_{key}', lambda value: None)(value)

    # Append attributes and return, and set clip property!!! This is critical
    # for tight bounding box calcs!
    for obj in objs:
        obj.set_clip_on(False)
    if isinstance(objs[0], mpatches.FancyBboxPatch):
        objs = objs[1:]
    return objs[0] if len(objs) == 1 else tuple(objs)


def _basemap_redirect(func):
    """
    Docorator that calls the basemap version of the function of the
    same name. This must be applied as the innermost decorator.
    """
    name = func.__name__

    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        if getattr(self, 'name', None) == 'proplot_basemap':
            return getattr(self.projection, name)(*args, ax=self, **kwargs)
        else:
            return func(self, *args, **kwargs)
    wrapper.__doc__ = None
    return wrapper


def _basemap_norecurse(func):
    """
    Decorator to prevent recursion in basemap method overrides.
    See `this post https://stackoverflow.com/a/37675810/4970632`__.
    """
    name = func.__name__
    func._called_from_basemap = False

    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        if func._called_from_basemap:
            result = getattr(maxes.Axes, name)(self, *args, **kwargs)
        else:
            with _state_context(func, _called_from_basemap=True):
                result = func(self, *args, **kwargs)
        return result
    return wrapper


def _process_wrapper(driver):
    """
    Generate generic wrapper decorator and dynamically modify the docstring
    to list methods wrapped by this function. Also set `__doc__` to ``None`` so
    that ProPlot fork of automodapi doesn't add these methods to the website
    documentation. Users can still call help(ax.method) because python looks
    for superclass method docstrings if a docstring is empty.
    """
    driver._docstring_orig = driver.__doc__ or ''
    driver._methods_wrapped = []
    proplot_methods = ('parametric', 'heatmap', 'area', 'areax')

    def decorator(func):
        # Define wrapper and suppress documentation
        # We only document wrapper functions, not the methods they wrap
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            return driver(self, func, *args, **kwargs)
        name = func.__name__
        if name not in proplot_methods:
            wrapper.__doc__ = None

        # List wrapped methods in the driver function docstring
        # Prevents us from having to both explicitly apply decorators in
        # axes.py and explicitly list functions *again* in this file
        docstring = driver._docstring_orig
        if '{methods}' in docstring:
            if name in proplot_methods:
                link = f'`~proplot.axes.Axes.{name}`'
            else:
                link = f'`~matplotlib.axes.Axes.{name}`'
            methods = driver._methods_wrapped
            if link not in methods:
                methods.append(link)
                string = (
                    ', '.join(methods[:-1])
                    + ',' * int(len(methods) > 2)  # Oxford comma bitches
                    + ' and ' * int(len(methods) > 1)
                    + methods[-1]
                )
                driver.__doc__ = docstring.format(methods=string)
        return wrapper
    return decorator


def _concatenate_docstrings(func):
    """
    Concatenate docstrings from a matplotlib axes method with a ProPlot axes
    method and obfuscate the call signature to avoid misleading users.

    Warning
    -------
    This is not yet used but will be in the future.
    """
    # NOTE: Originally had idea to use numpydoc.docscrape.NumpyDocString to
    # interpolate docstrings but *enormous* number of assupmtions would go into
    # this. And simple is better than complex.
    # Get matplotlib axes func
    name = func.__name__
    orig = getattr(maxes.Axes, name)
    odoc = inspect.getdoc(orig)
    if not odoc:  # should never happen
        return func

    # Prepend summary and potentially bail
    # TODO: Does this break anything on sphinx website?
    fdoc = inspect.getdoc(func) or ''  # also dedents
    regex = re.search(r'\.( | *\n|\Z)', odoc)
    if regex:
        fdoc = odoc[:regex.start() + 1] + '\n\n' + fdoc
    if rc['docstring.hardcopy']:  # True when running sphinx
        func.__doc__ = fdoc
        return func

    # Obfuscate signature by converting to *args **kwargs. Note this does
    # not change behavior of function! Copy parameters from a dummy function
    # because I'm too lazy to figure out inspect.Parameters API
    # See: https://stackoverflow.com/a/33112180/4970632
    dsig = inspect.signature(lambda *args, **kwargs: None)
    fsig = inspect.signature(func)
    func.__signature__ = fsig.replace(parameters=tuple(dsig.parameters.values()))

    # Concatenate docstrings and copy summary
    # Make sure different sections are very visible
    pad = '=' * len(name)
    doc = f"""
    ================================{pad}
    proplot.axes.Axes.{name} documentation
    ================================{pad}
    {fdoc}
    ==================================={pad}
    matplotlib.axes.Axes.{name} documentation
    ==================================={pad}
    {odoc}
    """
    func.__doc__ = inspect.cleandoc(doc)  # dedents and trims whitespace

    # Return
    return func


# Generate decorators and fill wrapper function docstrings. Each wrapper
# function should call function(self, ...) somewhere.
# Hidden wrapper functions providing only internal functionality
_hist_extras = _process_wrapper(_hist_extras)
_stem_extras = _process_wrapper(_stem_extras)
_plot_extras = _process_wrapper(_plot_extras)
_parametric_extras = _process_wrapper(_parametric_extras)
# Public wrapper functions providing important functionality
_apply_cmap = _process_wrapper(apply_cmap)
_apply_cycle = _process_wrapper(apply_cycle)
_bar_extras = _process_wrapper(bar_extras)
_barh_extras = _process_wrapper(barh_extras)
_boxplot_extras = _process_wrapper(boxplot_extras)
_default_latlon = _process_wrapper(default_latlon)
_default_transform = _process_wrapper(default_transform)
_fill_between_extras = _process_wrapper(fill_between_extras)
_fill_betweenx_extras = _process_wrapper(fill_betweenx_extras)
_hlines_extras = _process_wrapper(hlines_extras)
_indicate_error = _process_wrapper(indicate_error)
_scatter_extras = _process_wrapper(scatter_extras)
_standardize_1d = _process_wrapper(standardize_1d)
_standardize_2d = _process_wrapper(standardize_2d)
_text_extras = _process_wrapper(text_extras)
_violinplot_extras = _process_wrapper(violinplot_extras)
_vlines_extras = _process_wrapper(vlines_extras)

# Deprecated
# NOTE: Unlike other deprecated functions don't bother importing them to the
# top-level namespace. Will clutter things and it would have been impossible
# for any user to work with these outside of the Axes interface.
(
    cmap_changer,
    cycle_changer,
    colorbar_wrapper,
    legend_wrapper,
    text_wrapper,
    bar_wrapper,
    barh_wrapper,
    fill_between_wrapper,
    fill_betweenx_wrapper,
    boxplot_wrapper,
    violinplot_wrapper,
    vlines_wrapper,
    hlines_wrapper,
    scatter_wrapper,
) = warnings._rename_objs(
    '0.7',
    cmap_changer=apply_cmap,
    cycle_changer=apply_cycle,
    colorbar_wrapper=colorbar_extras,
    legend_wrapper=legend_extras,
    text_wrapper=text_extras,
    bar_wrapper=bar_extras,
    barh_wrapper=barh_extras,
    fill_between_wrapper=fill_between_extras,
    fill_betweenx_wrapper=fill_betweenx_extras,
    boxplot_wrapper=boxplot_extras,
    violinplot_wrapper=violinplot_extras,
    vlines_wrapper=vlines_extras,
    hlines_wrapper=hlines_extras,
    scatter_wrapper=scatter_extras,
)
