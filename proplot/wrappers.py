#!/usr/bin/env python3
"""
Imported by axes.py, wrappers for various plotting functions.
"""
import re
import numpy as np
import warnings
import functools
from . import utils, colortools, fonttools, axistools
from .utils import _default, ic
import cycler
import matplotlib.axes as maxes
import matplotlib.contour as mcontour
import matplotlib.ticker as mticker
import matplotlib.transforms as mtransforms
import matplotlib.patheffects as mpatheffects
import matplotlib.collections as mcollections
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.artist as martist
import matplotlib.legend as mlegend
from numbers import Number
from .rcmod import rc

# Xarray and pandas integration
ndarray = np.ndarray
try:
    from xarray import DataArray
except ModuleNotFoundError:
    DataArray = ndarray
try:
    from pandas import DataFrame, Series, Index
except ModuleNotFoundError:
    DataFrame, Series, Index = ndarray, ndarray, ndarray

# Cartopy
try:
    from cartopy.crs import PlateCarree
    from cartopy.mpl.geoaxes import GeoAxes
except ModuleNotFoundError:
    PlateCarree, GeoAxes = object, object

# Methods for wrapping
# TODO: 'quiver', 'streamplot' for cmap?
_edges_methods = ('pcolor', 'pcolormesh',)
_centers_methods = ('contour', 'contourf', 'quiver', 'streamplot', 'barbs')
_2d_methods = (*_centers_methods, *_edges_methods)
_1d_methods = ('plot', 'scatter', 'bar', 'hist', 'boxplot', 'violinplot', 'pie', 'fill_between', 'fill_betweenx', 'hexbin')
_cycle_methods = ('plot', 'scatter', 'bar', 'hist', 'boxplot', 'violinplot', 'pie', 'fill_between', 'fill_betweenx')
_cmap_methods = ('contour', 'contourf', 'pcolor', 'pcolormesh',
    'tripcolor', 'tricontour', 'tricontourf', 'cmapline',
    'hexbin', 'matshow', 'imshow', 'spy', 'hist2d')
_crs_methods = ('get_extent', 'set_extent', 'set_xticks', 'set_yticks',) # adds crs=PlateCarree()
_latlon_methods = ('plot', 'scatter', *_edges_methods, *_centers_methods) # adds latlon=True
_transform_methods = ('plot', 'scatter', *_edges_methods, *_centers_methods, 'tripcolor', 'tricontour', 'tricontourf',) # adds transform=PlateCarree()

# Disabled methods; keys are error messages
# TODO: rigorous support for violin plots, bar, barh, streamline and quiver
# TODO: 'table', 'eventplot', 'pie', 'xcorr', 'acorr', 'psd', 'csd',
# 'magnitude_spectrum', 'angle_spectrum', 'phase_spectrum', 'cohere', 'specgram'
_disabled_methods = {
    "Redundant function {} has been disabled. Control axis scale with format(xscale='scale', yscale='scale').":
        ('semilogx', 'semilogy', 'loglog'),
    "Redundant function {} has been disabled. Date formatters will be used automatically when x/y coordinates are python datetime or numpy datetime64.":
        ('plot_date',),
    "Redundant function {} has been disabled. Use proj='polar' in subplots() call, then use the angle as your *x* coordinate and radius as your *y* coordinate.":
        ('polar',)
    }
_map_disabled_methods = (
    # These are obvious
    'matshow', 'imshow', 'spy', # don't disable 'bar' or 'barh', can be used in polar plots
    'hist', 'hist2d', 'errorbar', 'boxplot', 'violinplot', 'step', 'stem',
    'hlines', 'vlines', 'axhline', 'axvline', 'axhspan', 'axvspan',
    # Look into these
    'stackplot', 'table', 'eventplot', 'pie',
    'xcorr', 'acorr', 'psd', 'csd', 'magnitude_spectrum',
    'angle_spectrum', 'phase_spectrum', 'cohere', 'specgram',
    )

# Keywords
# These may automatically override the 'fix' option!
_cmap_options = {
    'contour$': {'colors':'colors', 'linewidths':'linewidths', 'linestyles':'linestyles'},
    '^cmapline': {'colors':'color',  'linewidths':'linewidth', 'linestyles':'linestyle'},
    '^pcolor':  {'colors':'edgecolors', 'linewidths':'linewidth', 'linestyles':'linestyle'},
    }
# Translator for inset colorbars and legends
_loc_translate = {
    # Numbers
    0:'best',
    1:'upper right',
    2:'upper left',
    3:'lower left',
    4:'lower right',
    5:'center left',
    6:'center right',
    7:'lower center',
    8:'upper center',
    9:'center',
    # Simple
    'i':'best', # for inset
    'b':'best',
    True:'best', # i.e. "just turn this thing on"
    'inset':'best',
    'ur':'upper right',
    'ul':'upper left',
    'll':'lower left',
    'lr':'lower right',
    # Centered
    'cr':'center right',
    'cl':'center left',
    'uc':'upper center',
    'lc':'lower center',
    }


#------------------------------------------------------------------------------#
# For documentation
#------------------------------------------------------------------------------#
def _sphinx_name(name):
    """Gets sphinx name."""
    if name=='cmapline':
        return f'`~proplot.axes.BaseAxes.{name}`'
    else:
        return f'`~matplotlib.axes.Axes.{name}`'

def _expand_methods_list(func):
    """Fills `_method_list` with a list of methods that link to matplotlib
    documentation. Previously had documented public tuples storing the method
    names, this is much cleaner."""
    doc = func.__doc__
    for name,methods in (
        ('_centers_methods',       _centers_methods),
        ('_edges_methods',         _edges_methods),
        ('_centers_edges_methods', (*_centers_methods, *_edges_methods)),
        ('_latlon_methods',        _latlon_methods),
        ('_crs_methods',           _crs_methods),
        ('_transform_methods',     _transform_methods),
        ('_cycle_methods',         _cycle_methods),
        ('_cmap_methods',          _cmap_methods),
        ('_disabled_methods',      (*(method for methods in _disabled_methods.values() for method in methods),)),
        ('_map_disabled_methods',  _map_disabled_methods),
        ):
        if f'`{name}`' not in doc:
            continue
        doc = re.sub(f'`{name}`',
            ', '.join(_sphinx_name(method) for method in methods[:-1])
            + ','*min((len(methods)-2, 1)) + f' and {_sphinx_name(methods[-1])}', doc)
    func.__doc__ = doc
    return func

#------------------------------------------------------------------------------#
# Standardized inputs and automatic formatting
# NOTE: These do not have to be used explicitly, they are called by wrappers
#------------------------------------------------------------------------------#
def _to_iloc(data):
    """Get indexible attribute of array, so we can perform axis wise operations."""
    return getattr(data, 'iloc', data)

def _to_array(data):
    """Convert to ndarray cleanly."""
    data = getattr(data, 'values', data)
    return np.array(data)

def _array_std(data):
    """Converts list of lists to array."""
    if not isinstance(data, (ndarray, DataArray, DataFrame, Series, Index)):
        data = np.array(data)
    if not np.iterable(data):
        data = np.atleast_1d(data)
    return data

def _auto_label(data, axis=None, units=True):
    """Gets data and label for pandas or xarray objects or their coordinates."""
    label = ''
    if isinstance(data, ndarray):
        if axis is not None and data.ndim>axis:
            data = np.arange(data.shape[axis])
    # Xarray with common NetCDF attribute names
    elif isinstance(data, DataArray):
        if axis is not None and data.ndim>axis:
            data = data.coords[data.dims[axis]]
        label = getattr(data, 'name', '') or ''
        for key in ('long_name', 'standard_name'):
            label = data.attrs.get(key, label)
        if units:
            units = data.attrs.get('units', '')
            if label and units:
                label = f'{label} ({units})'
            elif units:
                label = units
    # Pandas object with name attribute
    # if not label and isinstance(data, DataFrame) and data.columns.size==1:
    elif isinstance(data, (DataFrame, Series, Index)):
        if axis==0 and isinstance(data, (DataFrame, Series)):
            data = data.index
        elif axis==1 and isinstance(data, DataFrame):
            data = data.columns
        elif axis is not None:
            data = np.arange(len(data)) # e.g. for Index
        label = getattr(data, 'name', '') or '' # DataFrame has no native name attribute but user can add one: https://github.com/pandas-dev/pandas/issues/447
    return data, str(label).strip()

def _autoformat_1d(self, func, *args, **kwargs):
    """Accepts 1d DataArray or Series, or
    2D DataArray or DataFrame, in which case list of lines or points
    are drawn. Used by `plot_wrapper` and `scatter_wrapper`."""
    # Sanitize input
    name = func.__name__
    if len(args)==0:
        return func(*args, **kwargs)
    elif len(args)==1:
        x = None
        y, *args = args
    elif len(args) in (2,3,4):
        x, y, *args = args # same
    else:
        raise ValueError(f'Too many arguments passed to {name}. Max is 4.')

    # Iterate through list of ys that we assume are identical
    # Standardize based on the first y input
    if len(args)>=1 and 'fill_between' in name:
        ys, args = (y, args[0]), args[1:]
    else:
        ys = (y,)
    ys = [_array_std(y) for y in ys]

    # Auto x coords
    y = ys[0] # test the first y input
    if x is None:
        axis = 1 if (name in ('boxplot','violinplot') or kwargs.get('means', None) or kwargs.get('medians', None)) else 0
        x, _ = _auto_label(y, axis=axis)
    x = _array_std(x)
    if x.ndim!=1:
        raise ValueError(f'x coordinates must be 1-dimensional, but got {x.ndim}.')

    # Auto formatting
    if self.figure._autoformat and not self._is_map:
        xaxis = 'y' if (kwargs.get('orientation', None)=='horizontal') else 'x'
        yaxis = 'x' if xaxis=='y' else 'y'
        # Ylabel
        kw = {}
        y, label = _auto_label(y)
        if label:
            axis = xaxis if name in ('hist',) else yaxis # for histogram, this indicates x coordinate
            kw[axis + 'label'] = label
        # Xlabel
        x, label = _auto_label(x)
        if label and name not in ('hist',):
            kw[xaxis + 'label'] = label
        if name!='scatter' and len(x)>1 and all(isinstance(x, Number) for x in x[:2]) and x[1]<x[0]:
            kw[xaxis + 'reverse'] = True
        # For bar plots, so we can offset by widths, need to convert to
        # array prematurely and apply the IndexFormatter
        if name in ('bar',):
            x = _to_array(x)
            if x.dtype=='object':
                kw[xaxis + 'formatter'] = mticker.IndexFormatter(x)
                kw[xaxis + 'minorlocator'] = mticker.NullLocator()
                x = np.arange(len(x))
        # Boxplot accepts a 'labels' argument but violinplot does not, so
        # we try to set the locators and formatters here. Note x is ignored later.
        if name in ('boxplot','violinplot'):
            x, label = _auto_label(y, axis=1)
            ys[0] = _to_array(y) # store naked array
            if label:
                kw[xaxis + 'label'] = label
            if x is not None: # TODO: support for e.g. date axes?
                x = _to_array(x)
                if x.dtype=='object':
                    kw[xaxis + 'formatter'] = mticker.IndexFormatter(x)
                    kw[xaxis + 'minorlocator'] = mticker.NullLocator()
                    if name=='boxplot':
                        kwargs['labels'] = x # requires or get overwritten
                if x.dtype!='object' or name=='violinplot':
                    kwargs['positions'] = np.arange(len(x))
        if kw:
            self.format(**kw)
    return func(x, *ys, *args, **kwargs)

def _autoformat_2d(self, func, *args, order='C', **kwargs):
    """Gets 2d data. Accepts ndarray and DataArray. Used by `enforce_centers`
    and `enforce_edges`, which are used for all 2d plot methods."""
    # Sanitize input
    name = func.__name__
    if len(args)==0:
        return func(*args, **kwargs)
    elif len(args)>4:
        raise ValueError(f'Too many arguments passed to {name}. Max is 4.')
    x, y = None, None
    if len(args)>2:
        x, y, *args = args

    # Ensure DataArray, DataFrame or ndarray
    # WARNING: Why is DataFrame always column major? Is this best behavior?
    Zs = []
    for Z in args:
        Z = _array_std(Z)
        if Z.ndim!=2:
            raise ValueError(f'Z must be 2-dimensional, got shape {Z.shape}.')
        Zs.append(Z)
    if not all(Zs[0].shape==Z.shape for Z in Zs):
        raise ValueError(f'Zs must be same shape, got shapes {[Z.shape for Z in Zs]}.')

    # Retrieve coordinates
    if x is None and y is None:
        Z = Zs[0]
        if order=='C': # TODO: check order stuff works
            idx, idy = 1, 0
        else:
            idx, idy = 0, 1
        if isinstance(Z, ndarray):
            x = np.arange(Z.shape[idx])
            y = np.arange(Z.shape[idy])
        elif isinstance(Z, DataArray): # DataArray
            x = Z.coords[Z.dims[idx]]
            y = Z.coords[Z.dims[idy]]
        else: # DataFrame; never Series or Index because these are 1d
            x = Z.index
            y = Z.columns

    # Check coordinates
    x, y = _array_std(x), _array_std(y)
    if x.ndim != y.ndim:
        raise ValueError(f'x coordinates are {x.ndim}-dimensional, but y coordinates are {y.ndim}-dimensional.')
    for name,array in zip(('x','y'), (x,y)):
        if array.ndim not in (1,2):
            raise ValueError(f'{name} coordinates are {array.ndim}-dimensional, but must be 1 or 2-dimensional.')

    # Auto formatting
    if self.figure._autoformat:
        kw = {}
        # Xlabel and Ylabel
        if not self._is_map:
            for key,xy in zip(('xlabel','ylabel'), (x,y)):
                _, label = _auto_label(xy)
                if label:
                    kw[key] = label
                if len(xy)>1 and all(isinstance(xy, Number) for xy in xy[:2]) and xy[1]<xy[0]:
                    kw[key[0] + 'reverse'] = True
        # Title
        _, title = _auto_label(Zs[0], units=False)
        if title:
            kw['title'] = title
        if kw:
            self.format(**kw)
    return func(x, y, *Zs, **kwargs)

#------------------------------------------------------------------------------
# 2D plot wrappers
#------------------------------------------------------------------------------
@_expand_methods_list
def enforce_centers(self, func, *args, order='C', **kwargs):
    """Wraps 2D plotting functions that take coordinate *centers* (`_centers_methods`),
    calculates centers if graticule *edges* were provided."""
    # Checks whether sizes match up, checks whether graticule was input
    x, y, *Zs = args
    xlen, ylen = x.shape[-1], y.shape[0]
    for Z in Zs:
        if Z.ndim!=2:
            raise ValueError(f'Input arrays must be 2D, instead got shape {Z.shape}.')
        elif Z.shape[1]==xlen-1 and Z.shape[0]==ylen-1 and x.ndim==1 and y.ndim==1:
            # Get centers given edges
            # NOTE: Do not print warning message because this use case is super
            # common, will just be annoying for user.
            if x.ndim==1 and y.ndim==1:
                x = (x[1:] + x[:-1])/2
                y = (y[1:] + y[:-1])/2
        elif Z.shape[1]!=xlen or Z.shape[0]!=ylen:
            raise ValueError(f'Input shapes x {x.shape} and y {y.shape} must match Z centers {Z.shape} or Z borders {tuple(i+1 for i in Z.shape)}.')
    # Optionally re-order
    # TODO: Double check this
    if order=='F':
        x, y = x.T, y.T # in case they are 2-dimensional
        Zs = (Z.T for Z in Zs)
    elif order!='C':
        raise ValueError(f'Invalid order "{order}". Choose from "C" (row-major, default) and "F" (column-major).')
    result = func(x, y, *Zs, **kwargs)
    return result

@_expand_methods_list
def enforce_edges(self, func, *args, order='C', **kwargs):
    """Wraps 2D plotting functions that take graticule *edges* (`_edges_methods`),
    calculates edges if coordinate *centers* were provided."""
    # Checks that sizes match up, checks whether graticule was input
    x, y, *Zs = args
    xlen, ylen = x.shape[-1], y.shape[0]
    for Z in Zs:
        if Z.ndim!=2:
            raise ValueError(f'Input arrays must be 2D, instead got shape {Z.shape}.')
        elif Z.shape[1]==xlen and Z.shape[0]==ylen:
            # If 2D, don't raise error, but don't fix either, because
            # matplotlib pcolor accepts grid center inputs.
            # NOTE: Do not print warning message because this use case is super
            # common, will just be annoying for user.
            if x.ndim==1 and y.ndim==1:
                x, y = utils.edges(x), utils.edges(y)
        elif Z.shape[1]!=xlen-1 or Z.shape[0]!=ylen-1:
            raise ValueError(f'Input shapes x {x.shape} and y {y.shape} must match Z centers {Z.shape} or Z borders {tuple(i+1 for i in Z.shape)}.')
    # Optionally re-order
    # TODO: Double check this
    if order=='F':
        x, y = x.T, y.T # in case they are 2-dimensional
        Zs = (Z.T for Z in Zs)
    elif order!='C':
        raise ValueError(f'Invalid order "{order}". Choose from "C" (row-major, default) and "F" (column-major).')
    result = func(x, y, *Zs, **kwargs)
    return result

#------------------------------------------------------------------------------#
# 1D plot wrappers
#------------------------------------------------------------------------------#
def plot_wrapper(self, func, *args, cmap=None, values=None, **kwargs):
    """
    Wraps `~matplotlib.axes.Axes.plot`, calls `~proplot.axes.BaseAxes.cmapline`
    if ``cmap`` is passed by the user.

    Parameters
    ----------
    *args
        Passed to `~matplotlib.axes.Axes.plot`.
    cmap, values
        Passed to `~proplot.axes.BaseAxes.cmapline`.
    **kwargs
        `~matplotlib.lines.Line2D` properties.
    """
    if len(args) not in (2,3): # e.g. with fmt string
        raise ValueError(f'Expected 1-3 plot args, got {len(args)}.')
    if cmap is None:
        lines = func(*args, **kwargs)
    else:
        lines = self.cmapline(*args, cmap=cmap, values=values, **kwargs)
    return lines

def scatter_wrapper(self, func, *args,
    s=None, size=None, markersize=None,
    c=None, color=None, markercolor=None,
    smin=None, smax=None,
    cmap=None, cmap_kw={}, vmin=None, vmax=None, norm=None, norm_kw={},
    lw=None, linewidth=None, linewidths=None, markeredgewidth=None, markeredgewidths=None,
    edgecolor=None, edgecolors=None, markeredgecolor=None, markeredgecolors=None,
    **kwargs):
    """
    Wraps `~matplotlib.axes.Axes.scatter`, adds optional keyword args
    more consistent with the `~matplotlib.axes.Axes.plot` keywords.

    Parameters
    ----------
    s, size, markersize : None, float, or list of float, optional
        Aliases for the marker size.
    smin, smax : None or float, optional
        Used to scale the `s` array. These are the minimum and maximum marker
        sizes. Defaults to the minimum and maximum of the `s` array.
    c, color, markercolor : None, color-spec, sequency of color-spec, or array, optional
        Aliases for the marker fill color. If just an array of values, colors
        will be generated by passing the values through the `norm` normalizer
        and drawing from the `cmap` colormap.
    cmap : None or colormap-spec, optional
        The colormap specifer, passed to the `~proplot.colortools.Colormap`
        constructor.
    cmap_kw : dict-like, optional
        Passed to `~proplot.colortools.Colormap`.
    vmin, vmax : None or float, optional
        Used to generate a `norm` for scaling the `c` array. These are the
        values corresponding to the leftmost and rightmost colors in the
        colormap. Defaults to the minimum and maximum values of the `c` array.
    norm : None or normalizer spec, optional
        The colormap normalizer, passed to the `~proplot.colortools.Norm`
        constructor.
    norm_kw : None or norm-spec, optional
        Passed to `~proplot.colortools.Norm`.
    lw, linewidth, linewidths, markeredgewidth, markeredgewidths : None or float, or list thereof, optional
        Aliases for the marker edge width.
    edgecolors, markeredgecolor, markeredgecolors : None or str or (R,G,B) tuple, or list thereof, optional
        Aliases for the marker edge color.
    **kwargs
        Passed to `~matplotlib.axes.Axes.scatter`.
    """
    # Manage input arguments
    # NOTE: Parse 1d must come before this
    x, y, *args = args
    if len(args)==2:
        c = args.pop(1)
    if len(args)==1:
        s = args.pop(0)
    if args:
        raise ValueError(f'Expected 1-4 scatter args, got {len(args)}.')
    # Format cmap and norm
    if cmap is not None:
        if isinstance(cmap, (str, dict, mcolors.Colormap)):
            cmap = cmap, # make a tuple
        cmap = colortools.Colormap(*cmap, N=None, **cmap_kw)
    if norm is not None:
        if isinstance(norm, (str, dict, mcolors.Normalize)):
            norm = norm,
        norm = colortools.Norm(*norm, N=None, **norm_kw)
    # Apply some aliases for keyword arguments
    c = _default(c, color, markercolor)
    s = _default(s, size, markersize)
    lw = _default(lw, linewidth, linewidths, markeredgewidth, markeredgewidths)
    ec = _default(edgecolor, edgecolors, markeredgecolor, markeredgecolors)
    # Scale s array
    if np.iterable(s):
        smin_true, smax_true = min(s), max(s)
        if smin is None:
            smin = smin_true
        if smax is None:
            smax = smax_true
        s = smin + (smax - smin)*(np.array(s) - smin_true)/(smax_true - smin_true)
    # Call function
    return func(x, y, c=c, s=s,
        cmap=cmap, vmin=vmin, vmax=vmax,
        norm=norm, linewidths=lw, edgecolors=ec,
        **kwargs)

def _fill_between_parse(func, *args, negcolor='blue', poscolor='red', negpos=False, **kwargs):
    """Parse args and call function."""
    # Allow common keyword usage
    xy = 'y' if 'x' in func.__name__ else 'y'
    yx = 'x' if xy=='y' else 'y'
    if xy in kwargs:
        args = (kwargs.pop(xy), *args)
    for yx in (yx + '1', yx + '2'):
        if yx in kwargs:
            args = (*args, kwargs.pop(yx))
    if len(args)==1:
        args = (np.arange(len(args[0])), *args)
    if len(args)==2:
        if kwargs.get('stacked', False):
            args = (*args, 0)
        else:
            args = (args[0], 0, args[1]) # default behavior
    if len(args)!=3:
        raise ValueError(f'Expected 2-3 positional args, got {len(args)}.')
    if not negpos:
        obj = func(*args, **kwargs)
        return obj
    # Get zero points
    objs = []
    kwargs.setdefault('interpolate', True)
    y1, y2 = np.atleast_1d(args[-2]).squeeze(), np.atleast_1d(args[-1]).squeeze()
    if y1.ndim>1 or y2.ndim>1:
        raise ValueError(f'When "negpos" is True, y must be 1-dimensional.')
    if kwargs.get('where', None) is not None:
        raise ValueError('When "negpos" is True, you cannot set the "where" keyword.')
    for i in range(2):
        kw = {**kwargs}
        kw.setdefault('color', negcolor if i==0 else poscolor)
        where = (y2<y1) if i==0 else (y2>=y1)
        obj = func(*args, where=where, **kw)
        objs.append(obj)
    return (*objs,)

def fill_between_wrapper(self, func, *args, **kwargs):
    """
    Wraps `~matplotlib.axes.Axes.fill_between`, also accessible via the
    `~proplot.axes.BaseAxes.area` alias.

    Parameters
    ----------
    *args : (y1,), (x,y1), or (x,y1,y2)
        The *x* and *y* coordinates. If `x` is not provided, it will be
        inferred from `y1`. If `y1` and `y2` are provided, their shapes
        must be identical, and we fill between respective columns of these
        arrays.
    stacked : bool, optional
        If `y2` is ``None``, this indicates whether to "stack" successive
        columns of the `y1` array.
    negpos : bool, optional
        Whether to shade where `y2` is greater than `y1` with the color `poscolor`,
        and where `y1` is greater than `y2` with the color `negcolor`. For
        example, to shade positive values red and negtive blue, use
        ``ax.fill_between(x, 0, y)``.
    negcolor, poscolor : color-spec, optional
        Colors to use for the negative and positive values. Ignored if `negpos`
        is ``False``.
    where : ndarray, optional
        Boolean ndarray mask for points you want to shade. See
        `this matplotlib example <https://matplotlib.org/3.1.0/gallery/pyplots/whats_new_98_4_fill_between.html#sphx-glr-gallery-pyplots-whats-new-98-4-fill-between-py>`__.
    **kwargs
        Passed to `~matplotlib.axes.Axes.fill_between`.
    """
    # WARNING: Unlike others, this wrapper is applied *before* parse_1d, so
    # we can handle common keyword arg usage.
    return _fill_between_parse(func, *args, **kwargs)

def fill_betweenx_wrapper(self, func, *args, **kwargs):
    """Wraps `~matplotlib.axes.Axes.fill_betweenx`, usage is same as `fill_between_wrapper`."""
    return _fill_between_parse(func, *args, **kwargs)

def barh_wrapper(self, func, *args, **kwargs):
    """Wraps `~matplotlib.axes.Axes.barh`, usage is same as `bar_wrapper`."""
    kwargs['orientation'] = 'horizontal'
    return self.bar(*args, **kwargs)

def bar_wrapper(self, func, *args, edgecolor=None, lw=None, linewidth=None,
    vert=True, orientation='vertical',
    stacked=False, medians=False, means=False,
    box=True, bar=True, boxrange=(25, 75), barrange=(5, 95),
    boxcolor=None, barcolor=None,
    boxlw=None, barlw=None,
    capsize=None,
    **kwargs):
    """
    Wraps `~matplotlib.axes.Axes.bar`, applies sensible defaults.

    Parameters
    ----------
    edgecolor : None or color-spec, optional
        The default edge color.
    lw, linewidth : None or float, optional
        The default edge width.
    orientation : {'vertical', 'horizontal'}, optional
        The orientation of the bars.
    vert : bool, optional
        Alternative to the `orientation` keyword arg. If ``False``, horizontal
        bars are drawn. This is for consistency with `~matplotlib.axes.Axes.boxplot`
        and `~matplotlib.axes.Axes.violinplot`.
    stacked : bool, optional
        Whether to stack columns of input data, or plot the bars side-by-side.
    medians : bool, optional
        Whether to plot the mean of columns of input data, instead of stacking
        or grouping the bars.
    means : bool, optional
        Whether to plot the median of columns of input data, instead of stacking
        or grouping the bars.
    box, bar : bool, optional
        Toggles thick and thin error bars when either of `means` or `medians`
        is ``True``.
    boxrange, barrange : (float, float), optional
        Percentile ranges for drawing thick and thin central error bars.
        The defaults are ``(25, 75)`` and ``(5, 95)``, respectively.
        Ignored if `medians` and `means` are both ``False``.
    boxcolor, barcolor : None or float, optional
        Colors for the thick and thin error bars.
    boxlw, barlw : None or float, optional
        Line widths for the thick and thin error bars.
    capsize : None or float, optional
        The cap size for thin error bars.
    """
    # Bail
    # TODO: Stacked feature is implemented in `cycle_wrapper`, but makes more
    # sense do document here; figure out way to move it here?
    # TODO: The 1d parser picks column index for x-coordinates if it detects
    # the means or medians keyword arg; move that behavior to this function?
    if not args:
        return func(*args, **kwargs)
    if len(args) not in (1,2,3,4):
        raise ValueError(f'Expected 1-4 positional args, got {len(args)}.')
    if len(args)==4: # note these are translated by cycle_wrapper for horizontal bars
        *args, kwargs['bottom'] = args
    if len(args)==3:
        *args, kwargs['width'] = args
    if 'left' in kwargs: # 'left' does not exist
        kwargs['bottom'] = kwargs.pop('left')
    # Defaults
    lw = _default(lw, linewidth, 0.7)
    edgecolor = _default(edgecolor, 'k')
    x, y, *args = args
    if not vert:
        orientation = 'horizontal'
    # Function
    iy = y
    if means:
        iy = y.mean(axis=0) # should work with DataFrame and DataArray
    if medians:
        if hasattr(y, 'median'): # true for DataFrame and DataArray
            iy = y.median(axis=0)
        else:
            iy = np.percentile(y, 50, axis=0)
    obj = func(x, iy, *args, linewidth=lw, edgecolor=edgecolor, stacked=stacked, orientation=orientation, **kwargs)
    # Add error bars
    if means or medians:
        if orientation=='horizontal':
            axis = 'x' # xerr
            xy = (iy,x)
        else:
            axis = 'y' # yerr
            xy = (x,iy)
        if box:
            err = np.percentile(y, boxrange, axis=0)
            err = err - np.array(iy)
            err[0,:] *= -1 # array now represents error bar sizes
            boxlw = _default(boxlw, 4*lw)
            boxcolor = _default(boxcolor, edgecolor)
            self.scatter(*xy, marker='o', color='white', s=boxlw, zorder=5)
            self.errorbar(*xy, **{axis+'err': err, 'capsize':0,
                'color':boxcolor, 'linestyle':'none', 'linewidth':boxlw})
        if bar: # note it is now impossible to make thin bar width different from cap width!
            err = np.percentile(y, barrange, axis=0)
            err = err - np.array(iy)
            err[0,:] *= -1 # array now represents error bar sizes
            barlw = _default(barlw, lw)
            barcolor = _default(barcolor, edgecolor)
            capsize = _default(capsize, 3) # better default than 5
            self.errorbar(*xy, **{axis+'err': err, 'capsize':capsize,
                'color':barcolor, 'linewidth':barlw, 'linestyle':'none',
                'markeredgecolor':barcolor, 'markeredgewidth':barlw})
    # Return object
    return obj

def boxplot_wrapper(self, func, *args,
    color=None, fill=True, fillcolor=None, fillalpha=None,
    lw=None, linewidth=None, orientation=None,
    marker=None, markersize=None,
    boxcolor=None, boxlw=None,
    capcolor=None, caplw=None,
    meancolor=None, meanlw=None,
    mediancolor=None, medianlw=None,
    whiskercolor=None, whiskerlw=None,
    fliercolor=None, flierlw=None,
    **kwargs):
    """
    Wraps `~matplotlib.axes.Axes.boxplot`, adds convenient keyword args.
    Fills the objects with a cycle color by default.

    Parameters
    ----------
    color : None or color-spec, optional
        The color of all objects.
    fill : bool, optional
        Whether to fill the box with a color.
    fillcolor : None or color-spec, optional
        The fill color for the boxes. Defaults to the next color cycler color.
    fillalpha : None or float, optional
        The opacity of the boxes. Defaults to ``1``.
    lw, linewidth : None or float, optional
        The linewidth of all objects.
    orientation : {None, 'horizontal', 'vertical'}, optional
        Alternative to the native `vert` keyword arg. Controls orientation.
    marker : None or marker-spec, optional
        Marker style for the 'fliers', i.e. outliers.
    markersize : None or float, optional
        Marker size for the 'fliers', i.e. outliers.
    boxcolor, capcolor, meancolor, mediancolor, whiskercolor : None or color-spec, optional
        The color of various boxplot components. These are shorthands so you
        don't have to pass e.g. a ``boxprops`` dictionary.
    boxlw, caplw, meanlw, medianlw, whiskerlw : None or float, optional
        The line width of various boxplot components. These are shorthands so you
        don't have to pass e.g. a ``boxprops`` dictionary.
    """
    # Call function
    if not args:
        return func(*args, **kwargs)
    if orientation is not None:
        if orientation=='horizontal':
            kwargs['vert'] = False
        elif orientation!='vertical':
            raise ValueError('Orientation must be "horizontal" or "vertical", got "{orientation}".')
    obj = func(*args, **kwargs)
    # Modify results
    # TODO: Pass props keyword args instead? Maybe does not matter.
    lw = _default(lw, linewidth)
    fillalpha = _default(fillalpha, 0.7)
    if fillcolor is None:
        cycler = next(self._get_lines.prop_cycler)
        fillcolor = cycler.get('color', None)
    for key,icolor,ilw in (
        ('boxes',boxcolor,boxlw),
        ('caps',capcolor,caplw),
        ('whiskers',whiskercolor,whiskerlw),
        ('means',meancolor,meanlw),
        ('medians',mediancolor,medianlw),
        ('fliers',fliercolor,flierlw),
        ):
        if key not in obj: # possible if not rendered
            continue
        artists = obj[key]
        ilw = _default(ilw, lw)
        icolor = _default(icolor, color)
        for artist in artists:
            if icolor is not None:
                artist.set_color(icolor)
                artist.set_markeredgecolor(icolor)
            if ilw is not None:
                artist.set_linewidth(ilw)
                artist.set_markeredgewidth(ilw)
            if key=='boxes' and fill:
                patch = mpatches.PathPatch(artist.get_path(), color=fillcolor, alpha=fillalpha, linewidth=0)
                self.add_artist(patch)
            if key=='fliers':
                if marker is not None:
                    artist.set_marker(marker)
                if markersize is not None:
                    artist.set_markersize(markersize)
    return obj

def violinplot_wrapper(self, func, *args,
    color=None, fillcolor=None, fillalpha=None,
    lw=None, linewidth=None, orientation=None,
    boxrange=(25, 75), barrange=(5, 95),
    showboxes=True, showbars=True, showmedians=None, showmeans=None,
    bodycolor=None, bodylw=None,
    boxcolor=None, boxlw=None,
    barcolor=None, barlw=None,
    meancolor=None, meanlw=None,
    mediancolor=None, medianlw=None,
    mincolor=None, minlw=None,
    maxcolor=None, maxlw=None,
    **kwargs):
    """
    Wraps `~matplotlib.axes.Axes.violinplot`, adds convenient keyword args.
    Makes the style shown in right plot of `this matplotlib example
    <https://matplotlib.org/3.1.0/gallery/statistics/customized_violin.html>`__
    the default. Among other things, it is no longer possible to show minima
    and maxima with whiskers, because that is redundant anyway.

    Parameters
    ----------
    color : None or color-spec, optional
        The color of all line objects. Defaults to ``'k'``.
    fillcolor : None or color-spec, optional
        The violin plot fill color. Defaults to the next color cycler color.
    fillalpha : None or float, optional
        The opacity of the violins. Defaults to ``1``.
    lw, linewidth : None or float, optional
        The linewidth of all line objects. Defaults to ``1``.
    orientation : {None, 'horizontal', 'vertical'}, optional
        Alternative to the native `vert` keyword arg. Controls orientation.
    boxrange, barrange : (float, float), optional
        Percentile ranges for the thick and thin central bars. The defaults
        are ``(25, 75)`` and ``(5, 95)``, respectively.
    showboxes, showbars, showmedians, showmeans : bool, optional
        Toggles showing the central thick vertical bar, thin vertical bar,
        median position marker, mean position marker, and the minimum maximum
        datum markers. These defaults are different from the matplotlib defaults.
    bodycolor, boxcolor, barcolor, meancolor, mediancolor, mincolor, maxcolor : None or color-spec, optional
        The color of various violin plot components. Applies to violin body,
        thick central bar, thin central bar, mean and median markers, and
        minimum and maximum markers.
    bodylw, boxlw, barlw, meanlw, medianlw, minlw, maxlw : None or float, optional
        The line width of various violin plot components. Applies to violin body,
        thick central bar, thin central bar, mean and median markers, and
        minimum and maximum markers.
    """
    # Call function
    if not args:
        return func(*args, **kwargs)
    if orientation is not None:
        if orientation=='horizontal':
            kwargs['vert'] = False
        elif orientation!='vertical':
            raise ValueError('Orientation must be "horizontal" or "vertical", got "{orientation}".')
    obj = func(*args, showmeans=False, showmedians=False, showextrema=False, **kwargs)
    # Defaults
    lw = _default(lw, linewidth)
    color = _default(color, 'k') # use black for edges
    fillalpha = _default(fillalpha, 0.7)
    if boxlw is None: # use a multiplier
        boxlw = 4*_default(lw, rc['lines.linewidth'])
    # Add custom thin and thick lines
    # WARNING: This wrapper comes after the parse 1d wrapper, which always
    # passes an (x,y) pair as args; but in this case the x is nonsense, just
    # indicates each datum index. We want the positions keyword arg.
    if 'positions' not in kwargs:
        raise ValueError('Could not find violinplot positions, but 1d parser should have supplied them.')
    x, y = kwargs['positions'], args[1]
    if showboxes:
        box1, box2 = np.percentile(y, boxrange, axis=0)
        obj['cboxes'] = self.vlines(x, box1, box2, linestyle='-')
    if showbars:
        bar1, bar2 = np.percentile(y, barrange, axis=0)
        obj['cbars'] = self.vlines(x, bar1, bar2, linestyle='-')
    if showmeans is None and showmedians is None:
        showmeans = True
    if showmeans:
        means = np.mean(y, axis=0)
        obj['cmeans'] = self.scatter(x, means, marker='o', color='white', s=boxlw, zorder=5)
    if showmedians:
        medians = np.percentile(y, 50, axis=0)
        obj['cmedians'] = self.scatter(x, medians, marker='o', color='white', s=boxlw, zorder=5)
    # Modify results
    # NOTE: In this case we can *only* modify after instantiation
    for key,icolor,ilw in (
        ('bodies',bodycolor,bodylw), # poly collection
        ('cboxes',boxcolor,boxlw), # line collections
        ('cbars',barcolor,barlw), # line collections
        ('cmeans',meancolor,meanlw),
        ('cmedians',mediancolor,medianlw),
        ):
        if key not in obj: # possible if not rendered
            continue
        artists = obj[key]
        ilw = _default(ilw, lw)
        icolor = _default(icolor, color)
        if not isinstance(artists, (list,tuple)): # only for bodies
            artists = [artists]
        for artist in artists:
            if key=='bodies':
                if isinstance(artist, tuple):
                    continue
                artist.set_alpha(fillalpha)
                if fillcolor is not None:
                    artist.set_facecolor(fillcolor)
                if icolor is not None:
                    artist.set_edgecolor(icolor)
                if ilw is not None:
                    artist.set_linewidths(ilw)
            else:
                if icolor is not None:
                    artist.set_color(icolor)
                if ilw is not None:
                    artist.set_linewidth(ilw)
    return obj

#------------------------------------------------------------------------------#
# Text wrapper
#------------------------------------------------------------------------------#
def text_wrapper(self, func, x, y, text, transform='data', fontname=None,
    border=False, bordercolor='w', invert=False, lw=None, linewidth=2,
    **kwargs):
    """
    Wraps `~matplotlib.axes.Axes.text`, enables specifying `tranform` with
    a string name and adds feature for drawing borders around text.

    Parameters
    ----------
    x, y : float
        The *x* and *y* coordinates for the text.
    text : str
        The text.
    transform : {'data', 'axes', 'figure'} or `~matplotlib.transforms.Transform`, optional
        The transform, or a string pointing to either of the
        `~matplotlib.axes.Axes.transData`, `~matplotlib.axes.Axes.transAxes`,
        or `~matplotlib.figure.Figure.transFigure` transforms. Default is
        ``'data'`` (unchanged).
    fontname : None or str, optional
        Alias for the ``fontfamily`` `~matplotlib.text.Text` property.
    border : bool, optional
        Whether to draw border around text.
    bordercolor : color-spec, optional
        The color of the border. Defaults to ``'w'``.
    invert : bool, optional
        If ``False``, ``'color'`` is used for the text and ``bordercolor``
        for the border. If ``True``, this is inverted.
    lw, linewidth : float, optional
        Ignored if `border` is ``False``. The width of the text border.

    Other parameters
    ----------------
    **kwargs
        Passed to `~matplotlib.text.Text` instantiator.
    """
    # Default transform by string name
    if not transform:
        transform = self.transData
    elif isinstance(transform, mtransforms.Transform):
        pass # do nothing
    elif transform=='figure':
        transform = self.figure.transFigure
    elif transform=='axes':
        transform = self.transAxes
    elif transform=='data':
        transform = self.transData
    else:
        raise ValueError(f"Unknown transform {transform}. Use string \"axes\" or \"data\".")
    # Font name strings
    if 'family' in kwargs: # builtin matplotlib alias
        kwargs.setdefault('fontfamily', kwargs.pop('family'))
    if fontname is not None:
        kwargs.setdefault('fontfamily', fontname)
    # Raise helpful error message if font unavailable
    font = kwargs.get('fontfamily', None)
    if font and font not in fonttools.fonts:
        warnings.warn(f'Font "{font}" unavailable. Available fonts are {", ".join(fonttools.fonts)}.')
        fonttools._missing_fonts.add(font)
        kwargs.pop('fontfamily')
    # Apply color default which is sometimes ignored?
    kwargs.setdefault('color', rc.get('text.color'))
    obj = func(x, y, text, transform=transform, **kwargs)
    # Draw border around text
    if border:
        linewidth = lw or linewidth
        facecolor, bgcolor = kwargs['color'], bordercolor
        if invert:
            facecolor, bgcolor = bgcolor, facecolor
        kwargs = {'linewidth':linewidth, 'foreground':bgcolor, 'joinstyle':'miter'}
        obj.update({
            'color':facecolor, 'zorder':100,
            'path_effects': [mpatheffects.Stroke(**kwargs), mpatheffects.Normal()]
            })
    return obj

#------------------------------------------------------------------------------#
# Geographic wrappers
#------------------------------------------------------------------------------#
# First basemap recursion fix
# Normally we *cannot* modify the underlying *axes* pcolormesh etc. because this
# this will cause basemap's self.m.pcolormesh etc. to use *custom* version and
# cause suite of weird errors. Prevent this recursion with the below decorator.
def _general_norecurse(self, func):
    """Decorator to prevent recursion in certain method overrides.
    See `this post https://stackoverflow.com/a/37675810/4970632`__."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        name = getattr(func, '__name__')
        if self._hasrecurred:
            # Return the *original* version of the matplotlib method (i.e.
            # the one we have not wrapped by overriding the __getattribute__
            # method). We reach this block e.g. when pcolormesh calls pcolor
            # or when basemap.Basemap tries to call something.
            self._hasrecurred = False
            result = object.__getattribute__(self, name)(*args, **kwargs)
        else:
            # Return the version we have wrapped
            self._hasrecurred = True
            result = func(*args, **kwargs)
        self._hasrecurred = False # cleanup, in case recursion never occurred
        return result
    return wrapper

def _basemap_call(self, func):
    """Docorator that calls the basemap version of the function of the same name."""
    name = func.__name__
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return self.m.__getattribute__(name)(ax=self, *args, **kwargs)
    return wrapper

@_expand_methods_list
def cartopy_transform(self, func, *args, transform=PlateCarree, **kwargs):
    """
    Wraps plotting functions for `~proplot.axes.CartopyAxes` (`_transform_methods`).

    With the default `~cartopy.mpl.geoaxes.GeoAxes` API, you need to pass
    ``transform=cartopy.crs.PlateCarree()`` if your data coordinates are
    longitude and latitude instead of map projection units. Now,
    ``transform=cartopy.crs.PlateCarree()`` is the default.
    """
    # Simple
    if isinstance(transform, type):
        transform = transform() # instantiate
    result = func(*args, transform=transform, **kwargs)
    # Re-enforce settings because some plot functions seem to reset the
    # outlinepatch or backgroundpatch (TODO: Double check this)
    self.format()
    return result

@_expand_methods_list
def cartopy_crs(self, func, *args, crs=PlateCarree, **kwargs):
    """
    Wraps axes functions for `~proplot.axes.CartopyAxes` (`_crs_methods`).

    As in `cartopy_transform`, but sets ``crs=cartopy.crs.PlateCarree()``
    as the default.
    """
    # Simple
    name = func.__name__
    if isinstance(crs, type):
        crs = crs() # instantiate
    try:
        result = func(*args, crs=crs, **kwargs)
    except TypeError as err:
        if not args:
            raise err
        args, crs = args[:-1], args[-1]
        result = func(*args, crs=crs, **kwargs)
    # Fix extent, so axes tight bounding box gets correct box!
    # From this issue: https://github.com/SciTools/cartopy/issues/1207#issuecomment-439975083
    # NOTE: May still get weird positioning because ProPlot assumes aspect
    # ratio from the full size projection, not the zoomed in version
    if name=='set_extent':
        clipped_path = self.outline_patch.orig_path.clip_to_bbox(self.viewLim)
        self.outline_patch._path = clipped_path
        self.background_patch._path = clipped_path
    return result

@_expand_methods_list
def basemap_latlon(self, func, *args, latlon=True, **kwargs):
    """
    Wraps plotting functions for `~proplot.axes.BasemapAxes` (`_latlon_methods`).

    With the default `~mpl_toolkits.basemap` API, you need to pass
    ``latlon=True`` if your data coordinates are longitude and latitude
    instead of map projection units. Now, ``latlon=True`` is the default.
    """
    return func(*args, latlon=latlon, **kwargs)

@_expand_methods_list
def cartopy_gridfix(self, func, lon, lat, *Zs, globe=False, **kwargs):
    """
    Wraps 2D plotting functions for `~proplot.axes.CartopyAxes` (`_centers_edges_methods`).

    Makes 1D longitude vectors monotonic and adds the `globe` keyword arg to
    optionally make data coverage *global*. Passing ``globe=True`` does the
    following:

    1. Makes longitudinal coverage *circular* (i.e. the last longitude coordinate
       equals the first longitude coordinate plus 360 degrees).
    2. Interpolates data to the North and South poles.

    If latitude and longitude arrays are 2D, `globe` is set to ``False``.
    """
    # Ensure monotonic lons or things get messed up
    # WARNING: In first geophysical data example, got persistent
    # changes to input data arrays without below line, resulting in non
    # monotonic coordinates and discovery of this error
    # WARNING: Cartopy contouring methods have issues with circularly wrapped data.
    # Triggers annoying ``TopologyException`` statements, which we suppress
    # with the IPython `~IPython.utils.io.capture_output` tool. See `this issue
    # <https://github.com/SciTools/cartopy/issues/946>`__.
    lon, lat = np.array(lon), np.array(lat) # no need to retain metadata on e.g. DataArray
    if not isinstance(kwargs.get('transform', None), PlateCarree):
        return func(lon, lat, *Zs, **kwargs)
    if lon.ndim==1 and not (lon<lon[0]).all(): # skip backwards data
        lonmin = lon[0]
        while True:
            filter_ = (lon<lonmin)
            if filter_.sum()==0:
                break
            lon[filter_] += 360
    # Below only works for vector data
    Zss = []
    for Z in Zs:
        if not globe or lon.ndim!=1 or lat.ndim!=1:
            Zss.append(Z)
            continue
        # 1) Fix holes over poles by *interpolating* there (equivalent to
        # simple mean of highest/lowest latitude points)
        Z_south = np.repeat(Z[0,:].mean(),  Z.shape[1])[None,:]
        Z_north = np.repeat(Z[-1,:].mean(), Z.shape[1])[None,:]
        lat = np.concatenate(([-90], lat, [90]))
        Z = np.concatenate((Z_south, Z, Z_north), axis=0)
        # 2) Fix seams at map boundary; by ensuring circular coverage
        if (lon[0] % 360) != ((lon[-1] + 360) % 360):
            lon = np.array((*lon, lon[0] + 360)) # make longitudes circular
            Z = np.concatenate((Z, Z[:,:1]), axis=1) # make data circular
        # Append
        Zss.append(Z)

    # Call function
    return func(lon, lat, *Zss, **kwargs)

@_expand_methods_list
def basemap_gridfix(self, func, lon, lat, *Zs, globe=False, **kwargs):
    """
    Wraps 2D plotting functions for `~proplot.axes.BasemapAxes` (`_centers_edges_methods`).

    Makes 1D longitude vectors monotonic and cycles them to fit within the map
    edges (i.e. if the projection central longitude is 90 degrees, will permute
    data to span from -90 degrees to 270 degrees longitude).

    Also adds the `globe` keyword arg to optionally make data coverage *global*.
    Passing ``globe=True`` does the following:

    1. Makes longitudinal coverage *circular* (i.e. the last longitude coordinate
       equals the first longitude coordinate plus 360 degrees).
    2. Interpolates data to the North and South poles.

    If latitude and longitude arrays are 2D, `globe` is set to ``False``.
    """
    # Bail out if map coordinates already provided
    lon, lat = np.array(lon), np.array(lat) # no need to retain metadata on e.g. DataArray
    if not kwargs.get('latlon', None):
        return func(lon, lat, *Zs, **kwargs)
    # Raise errors
    eps = 1e-3
    lonmin, lonmax = self.m.lonmin, self.m.lonmax
    if lon.max() > lon.min() + 360 + eps:
        raise ValueError(f'Longitudes span {lon.min()} to {lon.max()}. Can only span 360 degrees at most.')
    if lon.min() < -360 or lon.max() > 360:
        raise ValueError(f'Longitudes span {lon.min()} to {lon.max()}. Must fall in range [-360, 360].')
    # Establish 360-degree range
    lon -= 720
    while True:
        filter_ = (lon<lonmin)
        if filter_.sum()==0:
            break
        lon[filter_] += 360
    # Below only works with vectors
    Zss = []
    for Z in Zs:
        if lon.ndim!=1 or lat.ndim!=1:
            Zss.append(Z)
            continue
        # 1. Roll, accounting for whether ends are identical
        # If go from 0,1,...,359,0 (i.e. the borders), returns id of first zero
        roll = -np.argmin(lon) # always returns *first* value
        if lon[0]==lon[-1]:
            lon = np.roll(lon[:-1], roll)
            lon = np.append(lon, lon[0]+360)
        else:
            lon = np.roll(lon, roll)
        Z = np.roll(Z, roll, axis=1)
        # 2. Roll in same direction some more, if some points on right-edge
        # extend more than 360 above the minimum longitude; *they* should be the
        # ones on west/left-hand-side of map
        lonroll = np.where(lon>lonmin+360)[0] # tuple of ids
        if lonroll: # non-empty
            roll = lon.size-min(lonroll) # e.g. if 10 lons, lonmax id is 9, we want to roll once
            lon = np.roll(lon, roll) # need to roll foreward
            Z = np.roll(Z, roll, axis=1) # roll again
            lon[:roll] -= 360 # retains monotonicity
        # 3. Set NaN where data not in range lonmin, lonmax
        # This needs to be done for some regional smaller projections or otherwise
        # might get weird side-effects due to having valid data way outside of the
        # map boundaries -- e.g. strange polygons inside an NaN region
        Z = Z.copy()
        if lon.size-1==Z.shape[1]: # test western/eastern grid cell edges
            # remove data where east boundary is east of min longitude or west
            # boundary is west of max longitude
            Z[:,(lon[1:]<lonmin) | (lon[:-1]>lonmax)] = np.nan
        elif lon.size==Z.shape[1]: # test the centers
            # this just tests centers and pads by one for safety
            # remember that a *slice* with no valid range just returns empty array
            where = np.where((lon<lonmin) | (lon>lonmax))[0]
            Z[:,where[1:-1]] = np.nan
        # Global coverage
        if globe:
            # 4. Fix holes over poles by interpolating there (equivalent to
            # simple mean of highest/lowest latitude points)
            Z_south = np.repeat(Z[0,:].mean(),  Z.shape[1])[None,:]
            Z_north = np.repeat(Z[-1,:].mean(), Z.shape[1])[None,:]
            lat = np.concatenate(([-90], lat, [90]))
            Z = np.concatenate((Z_south, Z, Z_north), axis=0)
            # 5. Fix seams at map boundary; 3 scenarios here:
            # a. Have edges (e.g. for pcolor), and they fit perfectly against basemap seams.
            # This does not augment size
            if lon[0]==lonmin and lon.size-1==Z.shape[1]: # borders fit perfectly
                pass # do nothing
            # b. Have edges (e.g. for pcolor), and the projection edge is in-between grid cell boundaries.
            # This augments size by 1.
            elif lon.size-1==Z.shape[1]: # no interpolation necessary; just make a new grid cell
                lon = np.append(lonmin, lon) # append way easier than concatenate
                lon[-1] = lonmin + 360 # we've added a new tiny cell to the end
                Z = np.concatenate((Z[:,-1:], Z), axis=1) # don't use pad; it messes up masked arrays
            # c. Have centers (e.g. for contourf), and we need to interpolate to the
            # left/right edges of the map boundary.
            # This augments size by 2.
            elif lon.size==Z.shape[1]: # linearly interpolate to the edges
                x = np.array([lon[-1], lon[0]+360]) # x
                if x[0] != x[1]:
                    y = np.concatenate((Z[:,-1:], Z[:,:1]), axis=1)
                    xq = lonmin+360
                    yq = (y[:,:1]*(x[1]-xq) + y[:,1:]*(xq-x[0]))/(x[1]-x[0]) # simple linear interp formula
                    Z = np.concatenate((yq, Z, yq), axis=1)
                    lon = np.append(np.append(lonmin, lon), lonmin+360)
            else:
                raise ValueError('Unexpected shape of longitude, latitude, data arrays.')
        # Add
        Zss.append(Z)

    # Prevent error where old boundary, drawn on a different axes, remains
    # on the Basemap instance, which means it is not in self.patches, which
    # means Basemap tries to draw it again so it can clip the *contours* by the
    # resulting path, which raises error because you can't draw on Artist on multiple axes
    self.m._mapboundarydrawn = self.boundary # stored the axes-specific boundary here

    # Convert to projection coordinates and call function
    lat[lat>90], lat[lat<-90] = 90, -90 # otherwise, weird stuff happens
    if lon.ndim==1 and lat.ndim==1:
        lon, lat = np.meshgrid(lon, lat)
    x, y = self.m(lon, lat)
    kwargs['latlon'] = False
    return func(x, y, *Zss, **kwargs)

#------------------------------------------------------------------------------#
# Colormaps and color cycles
#------------------------------------------------------------------------------#
def _get_panel(self, loc):
    """Used to interpret colorbar=x and legend=x keyword args.
    And get inset location."""
    # Panel index
    if np.iterable(loc) and not isinstance(loc, str) and len(loc)==2:
        idx, loc = loc
    else:
        idx = 0
    # Try getting the panel
    ax = None
    if isinstance(loc, str):
        ax = getattr(self, loc + 'panel', None)
    if ax is not None:
        # Verify panel is available
        try:
            ax = ax[idx]
        except IndexError:
            raise ValueError(f'Stack index {idx} for panel "{loc}" is invalid. You must make room for it in your call to subplots() with e.g. axcolorbars_kw={{"{loc}stack":2}}.')
        if not ax or not ax.get_visible():
            raise ValueError(f'Panel "{loc}" does not exist. You must make room for it in your call to subplots() with e.g. axcolorbars="{loc}".')
        loc = 'fill'
    else:
        # Translate to a location
        ax = self
        loc = _loc_translate.get(loc, loc)
    return ax, loc

@_expand_methods_list
def cycle_wrapper(self, func, *args,
    cycle=None, cycle_kw={},
    markers=None, linestyles=None,
    label=None, labels=None, values=None,
    legend=None, legend_kw={},
    colorbar=None, colorbar_kw={},
    **kwargs):
    """
    Wraps methods that use the property cycler (`_cycle_methods`),
    adds features for controlling colors in the property cycler and drawing
    legends or colorbars in one go. Critically, this wrapper also **standardizes
    acceptable input** -- these methods now all accept 2D arrays holding columns
    of data, and *x*-coordinates are always optional. Note this alters the
    behavior of `~matplotlib.axes.Axes.boxplot` and `~matplotlib.axes.Axes.violinplot`,
    which now compile statistics on *columns* of data instead of *rows*.

    Parameters
    ----------
    cycle : None or cycle-spec, optional
        The cycle specifer, passed to the `~proplot.colortools.Cycle`
        constructor. If the returned list of colors is unchanged from the
        current axes color cycler, the axes cycle will **not** be reset to the
        first position.
    cycle_kw : dict-like, optional
        Passed to `~proplot.colortools.Cycle`.
    labels, values : None or list, optional
        The legend labels or colorbar coordinates for each line in the
        input array. Can be numeric or string.
    legend : bool or str, optional
        Whether to draw a legend from the resulting handle list.
        If ``'l'``, ``'r'``, ``'b'``, ``'left'``, ``'right'``, or ``'bottom'``,
        an axes panel is filled with a legend. Note in this case that the
        panel must already exist (i.e. it was generated by your call to
        `~proplot.subplots.subplots`)!

        Otherwise, an *inset* legend is drawn, and this sets the position.
        The following locations and location aliases are valid:

        ==================  ==============================================
        Location            Valid aliases
        ==================  ==============================================
        ``'best'``          ``0``, ``'b'``, ``'i'``, ``'inset'``, ``True``
        ``'upper right'``   ``1``, ``'ur'``
        ``'upper left'``    ``2``, ``'ul'``
        ``'lower left'``    ``3``, ``'ll'``
        ``'lower right'``   ``4``, ``'lr'``
        ``'center left'``   ``5``, ``'cl'``
        ``'center right'``  ``6``, ``'cr'``
        ``'lower center'``  ``7``, ``'lc'``
        ``'upper center'``  ``8``, ``'uc'``
        ``'center'``        ``9``, ``'c'``
        ==================  ==============================================

    legend_kw : dict-like, optional
        Ignored if `legend` is ``None``. Extra keyword args for our call
        to `~proplot.axes.BaseAxes` `~proplot.axes.BaseAxes.legend` or
        `~proplot.axes.PanelAxes` `~proplot.axes.PanelAxes.legend`.
    colorbar : bool or str, optional
        As with `legend`, but draws a panel or inset *colorbar*. For valid
        inset colorbar locations, see the `~proplot.axes.BaseAxes.colorbar`
        method.
    colorbar_kw : dict-like, optional
        Ignored if `colorbar` is ``None``. Extra keyword args for our call
        to the `~proplot.axes.BaseAxes` `~proplot.axes.BaseAxes.colorbar` or
        `~proplot.axes.PanelAxes` `~proplot.axes.PanelAxes.colorbar` methods.

    Other parameters
    ----------------
    *args, **kwargs
        Passed to the matplotlib plotting method.

    See also
    --------
    `~proplot.axes.BaseAxes`, `~proplot.colortools.Cycle`

    Note
    ----
    See the `matplotlib source 
    <https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/axes/_base.py>`_.
    The `set_prop_cycle` command modifies underlying 
    `_get_lines` and `_get_patches_for_fill`.
    """
    # Test input
    # NOTE: Requires _autoformat_1d wrapper before reaching this. Also note that
    # the 'x' coordinates are sometimes ignored below.
    name = func.__name__
    if not args:
        return func(*args, **kwargs)
    barh = (name=='bar' and kwargs.get('orientation', None)=='horizontal')
    x, y, *args = args
    if len(args)>=1 and 'fill_between' in name:
        ys = (y, args[0])
        args = args[1:]
    else:
        ys = (y,)
    is1d = (y.ndim==1)

    # Determine and temporarily set cycler
    # NOTE: Axes cycle has no getter, only set_prop_cycle, which sets a
    # prop_cycler attribute on the hidden _get_lines and _get_patches_for_fill
    # objects. This is the only way to query current axes cycler! Should not
    # wrap set_prop_cycle because would get messy and fragile.
    # NOTE: The _get_lines cycler is an *itertools cycler*. Has no length, so
    # we must cycle over it with next(). We try calling next() the same number
    # of times as the length of user input cycle. If the input cycle *is* in
    # fact the same, below does not reset the color position, cycles us to start!
    if cycle is not None or cycle_kw:
        # Get the new cycler
        cycle_kw = {**cycle_kw} # copy
        if not is1d and y.shape[1]>1: # default samples count
            cycle_kw.setdefault('samples', y.shape[1])
        cycle = colortools.Cycle(cycle, **cycle_kw)
        # Get the original property cycle
        # WARNING: Matplotlib saves itertools.cycle(cycler), not the original
        # cycler object, so we must build up the keys again.
        i = 0
        by_key = {}
        cycle_orig = self._get_lines.prop_cycler
        for i in range(len(cycle)): # use the cycler object length as a guess
            prop = next(cycle_orig)
            for key,value in prop.items():
                if key not in by_key:
                    by_key[key] = {*()} # set
                by_key[key].add(value)
        # Reset property cycler if it differs
        reset = ({*by_key} != {*cycle.by_key()}) # reset if keys are different
        if not reset: # test individual entries
            for key,value in cycle.by_key().items():
                if by_key[key] != {*value}:
                    reset = True
                    break
        if reset:
            self.set_prop_cycle(cycle)

    # Custom property cycler additions
    # NOTE: By default matplotlib uses _get_patches_for_fill.get_next_color
    # for scatter properties! So we simultaneously iterate through the
    # _get_lines property cycler and apply them.
    apply = {*()} # which keys to apply from property cycler
    if name=='scatter':
        # Figure out which props should be updated
        keys = {*self._get_lines._prop_keys} - {'color','linestyle','dashes'} # color already applied, linestyle ignored
        for key,prop in (
            ('markersize','s'),
            ('linewidth','linewidths'),
            ('markeredgewidth','linewidths'),
            ('markeredgecolor','edgecolors'),
            ('alpha','alpha'),
            ('marker','marker'),
            ):
            prop = kwargs.get(prop,None)
            if key in keys and prop is None:
                apply.add(key)

    # Plot susccessive columns
    # WARNING: Most methods that accept 2d arrays use columns of data, but when
    # pandas DataFrame passed to hist, boxplot, or violinplot, rows of data assumed!
    # This is fixed in parse_1d by converting to values.
    objs = []
    ncols = 1
    label_leg = None # for colorbar or legend
    labels = _default(values, labels, label, None)
    stacked = kwargs.pop('stacked', False) # for 'bar', 'barh', 'area', 'areax' plots
    if name in ('pie','boxplot','violinplot'):
        if labels is not None:
            kwargs['labels'] = labels
    else:
        ncols = (1 if is1d else y.shape[1])
        if labels is None or isinstance(labels, str):
            labels = [labels]*ncols
    if name in ('bar',):
        width = kwargs.pop('width', 0.8) # for bar plots; 0.8 is matplotlib default
        kwargs['height' if barh else 'width'] = width if stacked else width/ncols
    for i in range(ncols):
        # Prop cycle properties
        kw = {**kwargs} # copy
        if apply:
            props = next(self._get_lines.prop_cycler)
            for key in apply:
                value = props[key]
                if key in ('size','markersize'):
                    key = 's'
                elif key in ('linewidth','markeredgewidth'): # translate
                    key = 'linewidths'
                elif key=='markeredgecolor':
                    key = 'edgecolors'
                kw[key] = value
        # Get x coordinates
        ix, iy = x, ys[0] # samples
        if name in ('pie',):
            kw['labels'] = _default(labels, ix) # TODO: move to pie wrapper?
        if name in ('bar',): # adjust
            if not stacked:
                ix = x + (i - ncols/2 + 0.5)*width/ncols
            elif stacked and not is1d:
                key = 'x' if barh else 'bottom'
                # iy = _to_iloc(iy)[:,:i].sum(axis=1) # sum of empty slice will be zero
                # kw[key] = iy
                kw[key] = _to_iloc(iy)[:,:i].sum(axis=1) # sum of empty slice will be zero
        # Get y coordinates and labels
        if name in ('pie','boxplot','violinplot'):
            iys = (iy,) # only ever have one y value, and cannot have legend labels
        else:
            # The coordinates
            if stacked and 'fill_between' in name:
                iys = tuple(iy if is1d else _to_iloc(iy)[:,:j].sum(axis=1) for j in (i,i+1))
            else:
                iys = tuple(iy if is1d else _to_iloc(iy)[:,i] for iy in ys)
            # Possible legend labels
            label = labels[i]
            values, label_leg = _auto_label(iy, axis=1) # _auto_label(iy) # e.g. a pd.Series name
            if label_leg and label is None:
                label = _to_array(values)[i]
            if label is not None:
                kw['label'] = label
        # Call with correct args
        xy = ()
        if barh: # special, use kwargs only!
            kw.update({'bottom':ix, 'width':iys[0]})
            kw.setdefault('x', kwargs.get('bottom', 0)) # must always be provided
        elif name in ('pie','hist','boxplot','violinplot'): # no x-coordinate
            xy = (*iys,)
        else: # has x-coordinates, and maybe more than one y
            xy = (ix,*iys)
        obj = func(*xy, *args, **kw)
        if isinstance(obj, (list,tuple)) and len(obj)==1: # plot always returns list or tuple
            obj = obj[0]
        objs.append(obj)

    # Add colorbar and/or legend
    if colorbar:
        # Add handles
        ax, loc = _get_panel(self, colorbar)
        if loc not in ax._auto_colorbar:
            ax._auto_colorbar[loc] = []
            ax._auto_colorbar_kw[loc] = {}
        ax._auto_colorbar[loc].extend(objs)
        # Add keywords
        kw = {**colorbar_kw} # copy
        if loc!='fill':
            kw['loc'] = loc
        if label_leg:
            kw['label'] = label_leg
        ax._auto_colorbar_kw[loc].update(kw)
    if legend:
        # Add handles
        ax, loc = _get_panel(self, legend)
        if loc not in ax._auto_legend:
            ax._auto_legend[loc] = []
            ax._auto_legend_kw[loc] = {}
        ax._auto_legend[loc].extend(objs)
        # Add keywords
        kw = {**legend_kw}
        if loc!='fill':
            kw['loc'] = loc
        if label_leg:
            kw['label'] = label_leg
        ax._auto_legend_kw[loc].update(kw)

    # Return
    # WARNING: Make sure plot always returns tuple of objects, and bar always
    # returns singleton unless we have bulk drawn bar plots! Other matplotlib
    # methods call these internally!
    if name in ('hist',):
        objs = tuple(obj[-1] for obj in objs) # just the patch objects
    if name=='plot':
        return (*objs,) # always return tuple of objects
    elif name in ('boxplot', 'violinplot'):
        return objs[0] # always singleton, because these methods accept the whole 2d object
    else:
        return objs[0] if is1d else (*objs,) # sensible default behavior

@_expand_methods_list
def cmap_wrapper(self, func, *args, cmap=None, cmap_kw={},
    norm=None, norm_kw={},
    extend='neither',
    values=None, levels=None, zero=False, vmin=None, vmax=None, # override levels to be *centered* on zero
    edgefix=None, labels=False, labels_kw={}, precision=2,
    colorbar=False, colorbar_kw={},
    lw=None, linewidth=None, linewidths=None,
    ls=None, linestyle=None, linestyles=None,
    color=None, colors=None, edgecolor=None, edgecolors=None,
    **kwargs):
    """
    Wraps methods that take a ``cmap`` argument (`_cmap_methods`),
    adds several new keyword args and features.
    Uses the `~proplot.colortools.BinNorm` normalizer to bin data into
    discrete color levels (see notes).

    Parameters
    ----------
    cmap : None or colormap spec, optional
        The colormap specifer, passed to the `~proplot.colortools.Colormap`
        constructor.
    cmap_kw : dict-like, optional
        Passed to `~proplot.colortools.Colormap`.
    norm : None or normalizer spec, optional
        The colormap normalizer, used to warp data before passing it
        to `~proplot.colortools.BinNorm`. This is passed to the
        `~proplot.colortools.Norm` constructor.
    norm_kw : dict-like, optional
        Passed to `~proplot.colortools.Norm`.
    extend : {'neither', 'both', 'min', 'max'}, optional
        Whether to assign a unique color to out-of-bounds data. Also means
        when the colorbar is drawn, colorbar "extensions" will be drawn
        (triangles by default).
    levels : None or int or list of float, optional
        If list of float, level *edges*. If integer, the number of level
        edges, and boundaries are chosen by matplotlib based on the data
        range. Default is ``rc['image.levels']``.

        Since this function also wraps `~matplotlib.axes.Axes.pcolor` and
        `~matplotlib.axes.Axes.pcolormesh`, this means they now
        accept the `levels` keyword arg. You can now discretize your
        colors in a ``pcolor`` plot just like with ``contourf``.
    values : None or list of float, optional
        The level *centers*. If not ``None``, infers levels using
        `~proplot.utils.edges` and overrides the `levels` keyword arg.
    zero : bool, optional
        Ignored if `values` or `levels` are lists.
        If colormap levels were automatically selected by matplotlib, toggle
        this to modify the levels to be **symmetric about zero**.
    edgefix : None or bool, optional
        Whether to fix the the `white-lines-between-filled-contours
        <https://stackoverflow.com/q/8263769/4970632>`__
        and `white-lines-between-pcolor-rectangles
        <https://stackoverflow.com/q/27092991/4970632>`__
        issues. This slows down figure rendering by a bit. Default is
        ``rc['image.edgefix']``.
    labels : bool, optional
        For `~matplotlib.axes.Axes.contour`, whether to add contour labels
        with `~matplotlib.axes.Axes.clabel`. For `~matplotlib.axes.Axes.pcolor`
        or `~matplotlib.axes.Axes.pcolormesh`, whether to add labels to the
        center of grid boxes. In the latter case, the text will be black
        when the luminance of the underlying grid box color is >50%, and
        white otherwise (see the `~proplot.colortools` documentation).
    labels_kw : dict-like, optional
        Ignored if `labels` is ``False``. Extra keyword args for the labels.
        For `~matplotlib.axes.Axes.contour`, passed to `~matplotlib.axes.Axes.clabel`.
        For `~matplotlib.axes.Axes.pcolor` or `~matplotlib.axes.Axes.pcolormesh`,
        passed to `~matplotlib.axes.Axes.text`.
    precision : int, optional
        Maximum precision (in decimal places) for the number labels.
        Number labels are generated with the `~proplot.axistools.SimpleFormatter`
        formatter, which allows us to limit the precision.
    colorbar : bool or str, optional
        Whether to draw a colorbar from the resulting mappable object.
        If ``'l'``, ``'r'``, ``'b'``, ``'left'``, ``'right'``, or ``'bottom'``,
        an axes panel is filled with a colorbar. Note in this case that the
        panel must already exist (i.e. it was generated by your call to
        `~proplot.subplots.subplots`)!

        Otherwise, an *inset* colorbar is drawn, and this sets the position
        (e.g. `'upper right'`). For valid inset colorbar locations, see the
        `~proplot.axes.BaseAxes.colorbar` method.
    colorbar_kw : dict-like, optional
        Ignored if `colorbar` is ``None``. Extra keyword args for our call
        to `~proplot.axes.BaseAxes` `~proplot.axes.BaseAxes.colorbar` or
        `~proplot.axes.PanelAxes` `~proplot.axes.PanelAxes.colorbar`.

    Other parameters
    ----------------
    lw, linewidth, linewidths
        Aliases. Refers to `linewidths` for `~matplotlib.axes.Axes.contour`,
        `linewidth` for `~matplotlib.axes.Axes.pcolor` and
        `~matplotlib.axes.Axes.pcolormesh`, and `linewidth` for `~proplot.axes.BaseAxes.cmapline`.
    ls, linestyle, linestyles
        Aliases. Refers to `linestyles` for `~matplotlib.axes.Axes.contour`,
        `linestyle` for `~matplotlib.axes.Axes.pcolor` and
        `~matplotlib.axes.Axes.pcolormesh`, and `linestyle` for `~proplot.axes.BaseAxes.cmapline`.
    color, colors, edgecolor, edgecolors
        Aliases. Refers to `colors` for `~matplotlib.axes.Axes.contour`,
        `edgecolors` for `~matplotlib.axes.Axes.pcolor` and
        `~matplotlib.axes.Axes.pcolormesh`, and `color` for `~proplot.axes.BaseAxes.cmapline`.
    *args, **kwargs
        Passed to the matplotlib plotting method.

    Note
    ----
    The `~proplot.colortools.BinNorm` normalizer, used with all colormap
    plots, makes sure that your "levels" always span the full range of colors
    in the colormap, whether you are extending max, min, neither, or both. By
    default, when you select `extend` not ``'both'``, matplotlib seems to just
    cut off the most intense colors (reserved for coloring "out of bounds"
    data), even though they are not being used.

    This could also be done by limiting the number of colors in the colormap lookup
    table by selecting a smaller ``N`` (see `~matplotlib.colors.LinearSegmentedColormap`).
    But I prefer the approach of always building colormaps with hi-res lookup
    tables, and leaving the job of normalizing data values to colormap locations
    to the `~matplotlib.colors.Normalize` object.

    See also
    --------
    `~proplot.axes.BaseAxes`, `~proplot.colortools.Colormap`,
    `~proplot.colortools.Norm`, `~proplot.colortools.BinNorm`,
    `~matplotlib.colors.Colormap`, `~matplotlib.colors.Normalize`
    """
    # Input levels
    # See this post: https://stackoverflow.com/a/48614231/4970632
    name = func.__name__
    if not args:
        return func(*args, **kwargs)
    if np.iterable(values):
        if name in ('cmapline',):
            kwargs['values'] = values
            levels = utils.edges(values) # special case, used by colorbar factory
        elif norm is None or isinstance(norm, str) and 'segment' in norm:
            levels = utils.edges(values) # special case, used by colorbar factory
        else:
            norm_tmp = colortools.Norm(norm, **norm_kw)
            levels = norm_tmp.inverse(utils.edges(norm_tmp(values)))
    # Input colormap and other args
    # TODO: Add similar support for quiver, and automatic quiver keys
    cyclic = False
    colors = _default(color, colors, edgecolor, edgecolors)
    levels = _default(levels, rc['image.levels']) # e.g. pcolormesh can auto-determine levels if you input a number
    if vmin is not None and vmax is not None and not np.iterable(levels):
        levels = np.linspace(vmin, vmax, levels)
    linewidths = _default(lw, linewidth, linewidths)
    linestyles = _default(ls, linestyle, linestyles)
    if not re.match('contour$', name): # contour, tricontour, i.e. not a method where cmap is optional
        cmap = cmap or rc['image.cmap']
    if cmap is not None:
        cmap = colortools.Colormap(cmap, N=None, **cmap_kw)
        cyclic = cmap._cyclic
        if cyclic and extend!='neither':
            warnings.warn(f'Cyclic colormap requires extend="neither". Overriding user input extend="{extend}".')
            extend = 'neither'
        kwargs['cmap'] = cmap
    if 'contour' in name: # contour, contourf, tricontour, tricontourf
        kwargs.update({'levels': levels, 'extend': extend})

    # Aspect ratio handling
    # NOTE: For some bizarre reason, in first pass to draw, aspect ratio calculated
    # from bbox will be inverted after using set_aspect and apply_aspect. Need
    # to store custom attribute.
    aspect = kwargs.get('aspect', rc['image.aspect'])
    if name in ('imshow', 'matshow', 'spy', 'hist2d') and aspect=='equal': # aspect ratio fix
        self._aspect_equal = args[-1].shape[1] / args[-1].shape[0]

    # Disable edgefix=True for certain keyword combinations, e.g. if user wants
    # white lines around their pcolor mesh.
    # TODO: Allow calling contourf with linewidth?
    for regex,names in _cmap_options.items():
        if not re.search(regex, name):
            continue
        # Different cmap functions may or may not accept 'colors', 'linewidths',
        # or 'linestyles' as arguments
        for key,value in (('colors',colors), ('linewidths',linewidths), ('linestyles',linestyles)):
            if value is None:
                continue
            if key not in names:
                if value:
                    raise ValueError(f'Unknown keyword arg {key} for function {name}.')
                continue
            edgefix = False # override!
            kwargs[names[key]] = value

    # Call function with custom kwargs
    obj = func(*args, **kwargs)

    # Colormap features
    if cmap:
        # Get levels automatically determined by contourf, or make them
        # from the automatically chosen pcolor/imshow clims.
        # TODO: This still is not respected for hexbin 'log' norm
        if not getattr(obj, 'extend', None):
            obj.extend = extend # will already be on there for some funcs
        if not np.iterable(levels): # i.e. was an integer
            # Some tools automatically generate levels, like contourf.
            # Others will just automatically impose clims, like pcolor.
            levels = getattr(obj, 'levels', np.linspace(*obj.get_clim(), levels))
            if all(levels==levels[0]):
                warnings.warn('Failed to infer colormap levels. Something is probably wrong.')
                levels[1:] = levels[0] + 1
            # Default linear normalizer
            if not zero:
                norm = _default(norm, 'linear') # matplotlib will have chosen *linearly* spaced levels, this helps us reduce BinNorm time
            # Get centered levels (when contourf has already rendered contours,
            # they cannot be changed or updated; must be redrawn)
            else: # non-linearly spaced, need
                abs_max = max([abs(max(levels)), abs(min(levels))])
                levels = np.linspace(-abs_max, abs_max, len(levels))
                if not hasattr(obj, 'levels'): # e.g. contourf, Artist must be re-drawn!
                    obj.set_clim(-abs_max, abs_max)
                else:
                    if hasattr(obj, 'collections'):
                        for artist in obj.collections:
                            artist.set_visible(False)
                    elif hasattr(obj, 'set_visible'):
                        obj.set_visible(False)
                    else:
                        raise ValueError(f'Unknown object {obj}. Cannot center colormap levels.')
                    kwargs['levels'] = levels
                    obj = func(*args, **kwargs)
        # Always add 'levels' as attribute, even for pcolor
        obj.levels = levels
        # Get 'pre-processor' norm -- e.g. maybe user wants colormap scaled
        # in logarithmic space, or warped to diverge from center from a midpoint
        # NOTE: LinearSegmentedNorm is slow! Test for linearity of input levels.
        diff = np.diff(levels)
        if (diff<0).any() or diff.size<2:
            raise ValueError(f'Need at least 3 monotonically increasing levels, got {levels}.')
        if norm is None:
            eps = diff.mean()/1e3
            if (np.abs(np.diff(diff)) >= eps).any():
                norm = 'segmented'
            else:
                norm = 'linear'
        # Set the normalizer
        step = 1.0
        if cyclic:
            step = 0.5
            extend = 'both'
        norm_preprocess = colortools.Norm(norm, levels=levels, clip=False, **norm_kw)
        norm_main = colortools.BinNorm(norm=norm_preprocess, levels=levels, step=step, extend=extend)
        obj.set_norm(norm_main)

    # Set labels
    if labels:
        # Very simple, use clabel args
        fmt = axistools.Formatter('simple', precision=precision)
        if name=='contour': # TODO: document alternate keyword args!
            labels_kw_ = {'fmt':fmt, 'inline_spacing':3, 'fontsize':rc['small']} # for rest, we keep the defaults
            for key1,key2 in (('size','fontsize'),):
                value = labels_kw.pop(key1, None)
                if value:
                    labels_kw[key2] = value
            labels_kw_.update(labels_kw)
            self.clabel(obj, **labels_kw_)
        # Label each box manually
        # See: https://stackoverflow.com/a/20998634/4970632
        elif 'pcolor' in name:
            obj.update_scalarmappable() # populates the _facecolors attribute, initially filled with just a single color
            labels_kw_ = {'size':rc['small'], 'ha':'center', 'va':'center'}
            labels_kw_.update(labels_kw)
            array = obj.get_array()
            paths = obj.get_paths()
            colors = obj.get_facecolors() # *flattened* list of objects
            for color,path,num in zip(colors,paths,array):
                if not np.isfinite(num):
                    continue
                bbox = path.get_extents()
                x = (bbox.xmin + bbox.xmax)/2
                y = (bbox.ymin + bbox.ymax)/2
                if 'color' not in labels_kw:
                    _, _, lum = colortools.to_xyz(color, 'hcl')
                    if lum<50:
                        color = 'w'
                    else:
                        color = 'k'
                    labels_kw_['color'] = color
                self.text(x, y, fmt(num), **labels_kw_)
        else:
            raise RuntimeError(f'Not possible to add labels to {name} plot.')

    # Fix white lines between filled contours/mesh, allow user to override!
    edgefix = _default(edgefix, rc['image.edgefix'])
    if edgefix:
        color = 'face'
        linewidth = 0.4 # seems to be lowest threshold where white lines disappear
        linestyle = '-'
        if 'pcolor' in name: # 'pcolor', 'pcolormesh', 'tripcolor'
            obj.set_edgecolor(color)
            obj.set_linewidth(linewidth) # seems to do the trick, without dots in corner being visible
        elif 'contourf' in name: # 'contourf', 'tricontourf'
            for contour in obj.collections:
                contour.set_edgecolor(color)
                contour.set_linewidth(linewidth)
                contour.set_linestyle(linestyle)

    # Add colorbar
    if colorbar:
        ax, loc = _get_panel(self, colorbar)
        colorbar_kw = {**colorbar_kw} # make copy of mutable default object!
        if 'label' not in colorbar_kw and self.figure._autoformat:
            _, label = _auto_label(args[-1]) # last one is data, we assume
            if label:
                colorbar_kw['label'] = label
        if name in ('cmapline',) and values is not None:
            colorbar_kw['values'] = values
        if loc!='fill':
            colorbar_kw['loc'] = loc
        ax.colorbar(obj, **colorbar_kw)
    return obj

#------------------------------------------------------------------------------#
# Legends and colorbars
#------------------------------------------------------------------------------#
def legend_wrapper(self, handles=None, labels=None, ncol=None, ncols=None,
    center=None, order='C', loc=None, label=None, title=None,
    color=None, marker=None, lw=None, linewidth=None,
    dashes=None, linestyle=None, markersize=None, frameon=None, frame=None,
    **kwargs):
    """
    Wraps `~matplotlib.axes.Axes` `~matplotlib.axes.Axes.legend` and
    `~proplot.axes.PanelAxes` `~proplot.axes.PanelAxes.legend`, adds some
    handy features.

    Parameters
    ----------
    handles : None or list of `~matplotlib.artist.Artist`, optional
        List of artists instances, or list of lists of artist instances (see
        the `center` keyword). If ``None``, the artists are retrieved with
        `~matplotlib.axes.Axes.get_legend_handles_labels`.
    labels : None or list of str, optional
        Matching list of string labels, or list of lists of string labels (see
        the `center` keywod). If ``None``, the labels are retrieved by calling
        `~matplotlib.artist.Artist.get_label` on each `~matplotlib.artist.Artist`
        in `handles`.
    ncol, ncols : int, optional
        The number of columns. `ncols` is an alias, added
        for consistency with `~matplotlib.pyplot.subplots`.
    order : {'C', 'F'}, optional
        Whether legend handles are drawn in row-major (``'C'``) or column-major
        (``'F'``) order. Analagous to `numpy.array` ordering. For some reason
        ``'F'`` was the original matplotlib default; the default is now ``'C'``.
    center : None or bool, optional
        Whether to center each legend row individually. If ``True``, we
        actually draw successive single-row legends stacked on top of each
        other.

        If ``None``, we infer this setting from `handles`. Defaults to ``True``
        if `handles` is a list of lists; each sublist is used as a *row*
        in the legend. Otherwise, defaults to ``False``.
    loc : None or str, optional
        The legend location. The following locations and location aliases are
        valid.

        ==================  ==============================================
        Location            Valid aliases
        ==================  ==============================================
        ``'best'``          ``0``, ``'b'``, ``'i'``, ``'inset'``, ``True``
        ``'upper right'``   ``1``, ``'ur'``
        ``'upper left'``    ``2``, ``'ul'``
        ``'lower left'``    ``3``, ``'ll'``
        ``'lower right'``   ``4``, ``'lr'``
        ``'center left'``   ``5``, ``'cl'``
        ``'center right'``  ``6``, ``'cr'``
        ``'lower center'``  ``7``, ``'lc'``
        ``'upper center'``  ``8``, ``'uc'``
        ``'center'``        ``9``, ``'c'``
        ==================  ==============================================

    label, title : None or str, optional
        The legend title. The `label` keyword is also accepted, for consistency
        with `colorbar_wrapper`.
    color, lw, linewidth, marker, linestyle, dashes, markersize : None or property-spec, optional
        Properties used to override the legend handles. For example, if you
        want a legend that describes variations in line style ignoring variations
        in color, you might want to use ``color='k'``. For now this does not
        include `facecolor`, `edgecolor`, and `alpha`, because
        `~matplotlib.axes.Axes.legend` uses these keyword args to modify the
        frame properties.

    Other parameters
    ----------------
    **kwargs
        Passed to `~matplotlib.axes.Axes.legend`.

    See also
    --------
    `~proplot.axes.BaseAxes.colorbar`, `~proplot.axes.PanelAxes.colorbar`, `~matplotlib.axes.Axes.legend`
    """
    # First get legend settings and interpret kwargs.
    if order not in ('F','C'):
        raise ValueError(f'Invalid order "{order}". Choose from "C" (row-major, default) and "F" (column-major).')
    ncol = _default(ncols, ncol) # may still be None, wait till later
    title = _default(label, title)
    frameon = _default(frame, frameon, rc['legend.frameon'])
    if title is not None:
        kwargs['title'] = title
    if frameon is not None:
        kwargs['frameon'] = frameon
    kwargs['prop'] = {'family': rc['fontname']} # 'prop' can be a FontProperties object or a dict for the kwargs

    # Automatically get labels and handles
    # Also accept non-list input
    if handles is None:
        if self._filled:
            raise ValueError('You must pass a handles list for panel axes "filled" with a legend.')
        else:
            handles, labels_default = self.get_legend_handles_labels() # ignores artists with labels '_nolegend_'
            if labels is None:
                labels = labels_default
            if not handles:
                raise ValueError('No labeled artists found. To generate a legend without providing the artists explicitly, pass label="label" in your plotting commands.')
    if not np.iterable(handles): # e.g. a mappable object
        handles = [handles]
    if labels is not None and (not np.iterable(labels) or isinstance(labels, str)):
        labels = [labels]

    # Legend entry for colormap or scatterplot object
    # TODO: Idea is we pass a scatter plot or contourf or whatever, and legend
    # is generating by drawing patch rectangles or markers with different colors.
    if any(not hasattr(handle, 'get_facecolor') and hasattr(handle, 'get_cmap') for handle in handles) and len(handles)>1:
        raise ValueError(f'Handles must be objects with get_facecolor attributes or a single mappable object from which we can draw colors.')

    # Build pairs of handles and labels
    # This allows alternative workflow where user specifies labels when
    # creating the legend.
    pairs = []
    list_of_lists = (not hasattr(handles[0], 'get_label')) # e.g. not including BarContainer
    if labels is None:
        for handle in handles:
            if list_of_lists:
                ipairs = []
                for ihandle in handle:
                    if not hasattr(ihandle, 'get_label'):
                        raise ValueError(f'Object {ihandle} must have a "get_label" attribute.')
                    ipairs.append((ihandle, ihandle.get_label()))
                pairs.append(ipairs)
            else:
                if not hasattr(handle, 'get_label'):
                    raise ValueError(f'Object {handle} must have a "get_label" attribute.')
                pairs.append((handle, handle.get_label()))
    else:
        if len(labels)!=len(handles):
            raise ValueError(f'Got {len(labels)} labels, but {len(handles)} handles.')
        for label,handle in zip(labels,handles):
            if list_of_lists:
                ipairs = []
                if not np.iterable(label) or isinstance(label, str):
                    raise ValueError(f'Got list of lists of handles, but just list of labels.')
                elif len(label)!=len(handle):
                    raise ValueError(f'Got {len(label)} labels in sublist, but {len(handle)} handles.')
                for ilabel,ihandle in zip(label,handle):
                    ipairs.append((ihandle, ilabel))
                pairs.append(ipairs)
            else:
                if not isinstance(label, str) and np.iterable(label):
                    raise ValueError(f'Got list of lists of labels, but just list of handles.')
                pairs.append((handle, label))

    # Manage pairs in context of 'center' option
    if center is None: # automatically guess
        center = list_of_lists
    elif center and list_of_lists and ncol is not None:
        warnings.warn('Detected list of *lists* of legend handles. Ignoring user input property "ncol".')
    elif not center and list_of_lists: # standardize format based on input
        list_of_lists = False # no longer is list of lists
        pairs = [pair for ipairs in pairs for pair in ipairs]
    elif center and not list_of_lists:
        list_of_lists = True
        ncol = _default(ncol, 3)
        pairs = [pairs[i*ncol:(i+1)*ncol] for i in range(len(pairs))] # to list of iterables
    if list_of_lists: # remove empty lists, pops up in some examples
        pairs = [ipairs for ipairs in pairs if ipairs]

    # Now draw legend(s)
    legs = []
    width, height = self.figure.get_size_inches()
    width, height = width*abs(self._position.width), height*abs(self._position.height)
    # Individual classic legend
    if not center:
        # Optionally change order
        # See: https://stackoverflow.com/q/10101141/4970632
        # Example: If 5 columns, but final row length 3, columns 0-2 have
        # N rows but 3-4 have N-1 rows.
        ncol = _default(ncol, 3)
        if order=='C':
            fpairs = []
            split = [pairs[i*ncol:(i+1)*ncol] for i in range(len(pairs)//ncol+1)] # split into rows
            nrowsmax, nfinalrow = len(split), len(split[-1]) # max possible row count, and columns in final row
            nrows = [nrowsmax]*nfinalrow + [nrowsmax-1]*(ncol-nfinalrow)
            for col,nrow in enumerate(nrows): # iterate through cols
                fpairs.extend(split[row][col] for row in range(nrow))
            pairs = fpairs
        if loc is not None:
            kwargs['loc'] = _loc_translate.get(loc, loc)
        # Make legend object
        leg = mlegend.Legend(self, *zip(*pairs), ncol=ncol, **kwargs)
        legs.append(leg)
    # Separate legend for each row. The label spacing/border spacing will be
    # exactly replicated, as if we were using the original legend command.
    else:
        # Message when overriding some properties
        overridden = []
        kwargs.pop('frameon', None) # then add back later!
        for override in ('bbox_transform', 'bbox_to_anchor'):
            prop = kwargs.pop(override, None)
            if prop is not None:
                overridden.append(override)
        if overridden:
            warnings.warn(f'For centered-row legends, must override user input properties "' + '", "'.join(prop for prop in overridden) + '".')
        # Default location
        loc = _loc_translate.get(_default(loc, 'upper center'), loc)
        if loc=='best':
            warnings.warn('For centered-row legends, cannot use "best" location. Defaulting to "upper center".')
            loc = 'upper center'
        # Determine space we want sub-legend to occupy as fraction of height
        # NOTE: Empirical testing shows spacing fudge factor necessary to exactly
        # replicate the spacing of standard aligned legends.
        fontsize = kwargs.get('fontsize', None) or rc['legend.fontsize']
        spacing  = kwargs.get('labelspacing', None) or rc['legend.labelspacing']
        interval = 1/len(pairs) # split up axes
        interval = (((1 + spacing*0.85)*fontsize)/72)/height
        # Iterate and draw
        # NOTE: We confine possible bounding box in *y*-direction, but do not
        # confine it in *x*-direction. Matplotlib will automatically move
        # left-to-right if you request this.
        ymin, ymax = None, None
        if order=='F':
            raise NotImplementedError(f'When center=True, ProPlot vertically stacks successive single-row legends. Column-major (order="F") ordering is un-supported.')
        for i,ipairs in enumerate(pairs):
            if i==1:
                kwargs.pop('title', None)
            if i>=1 and title is not None:
                i += 1 # extra space!
            # Legend position
            if 'upper' in loc:
                y1 = 1 - (i+1)*interval
                y2 = 1 - i*interval
            elif 'lower' in loc:
                y1 = (len(pairs) + i - 2)*interval
                y2 = (len(pairs) + i - 1)*interval
            else: # center
                y1 = 0.5 + interval*len(pairs)/2 - (i+1)*interval
                y2 = 0.5 + interval*len(pairs)/2 - i*interval
            ymin = min(y1, _default(ymin, y1))
            ymax = max(y2, _default(ymax, y2))
            # Draw legend
            bbox = mtransforms.Bbox([[0, y1], [1, y2]])
            leg = mlegend.Legend(self, *zip(*ipairs), loc=loc, ncol=len(ipairs),
                bbox_transform=self.transAxes, bbox_to_anchor=bbox, frameon=False,
                **kwargs) # _format_legend is overriding original legend Method
            legs.append(leg)

    # Add legends manually so matplotlib does not remove old ones
    # Also apply override settings
    override = {}
    outline = rc.fill({
        'linewidth':'axes.linewidth',
        'edgecolor':'axes.edgecolor',
        'facecolor':'axes.facecolor',
        'alpha':'legend.framealpha',
        }, cache=False)
    for key in (*outline,):
        if key!='linewidth':
            if kwargs.get(key, None):
                outline.pop(key, None)
    for key,value in (
        ('color',color),
        ('marker',marker),
        ('linewidth',lw),
        ('linewidth',linewidth),
        ('markersize',markersize),
        ('linestyle',linestyle),
        ('dashes',dashes),
        ):
        if value is not None:
            override[key] = value
    for leg in legs:
        self.add_artist(leg)
        leg.legendPatch.update(outline) # or get_frame()
        for obj in leg.legendHandles:
            obj.update(override)
    # Draw manual fancy bounding box for un-aligned legend
    # WARNING: The matplotlib legendPatch transform is the default transform,
    # i.e. universal coordinates in points. Means we have to transform
    # mutation scale into transAxes sizes.
    # WARNING: Tempting to use legendPatch for everything but for some reason
    # coordinates are messed up. In some tests all coordinates were just result
    # of get window extent multiplied by 2 (???). Anyway actual box is found in
    # _legend_box attribute, which is accessed by get_window_extent.
    if center and frameon:
        if len(legs)==1:
            legs[0].set_frame_on(True) # easy!
        else:
            # Get coordinates
            renderer = self.figure.canvas.get_renderer()
            bboxs = [leg.get_window_extent(renderer).transformed(self.transAxes.inverted()) for leg in legs]
            xmin, xmax = min(bbox.xmin for bbox in bboxs), max(bbox.xmax for bbox in bboxs)
            ymin, ymax = min(bbox.ymin for bbox in bboxs), max(bbox.ymax for bbox in bboxs)
            fontsize = (fontsize/72)/width # axes relative units
            fontsize = renderer.points_to_pixels(fontsize)
            # Draw and format patch
            patch = mpatches.FancyBboxPatch((xmin,ymin), xmax-xmin, ymax-ymin,
                    snap=True, zorder=4.5,
                    mutation_scale=fontsize, transform=self.transAxes) # fontsize defined in if statement
            if kwargs.get('fancybox', rc['legend.fancybox']):
                patch.set_boxstyle('round', pad=0, rounding_size=0.2)
            else:
                patch.set_boxstyle('square', pad=0)
            patch.set_clip_on(False)
            patch.update(outline)
            self.add_artist(patch)
            # Add shadow
            # TODO: This does not work, figure out
            if kwargs.get('shadow', rc['legend.shadow']):
                shadow = mpatches.Shadow(patch, 20, -20)
                self.add_artist(patch)
            # Add patch to list
            legs = (patch, *legs)
    # Return legend(s)
    return legs[0] if len(legs)==1 else (*legs,)

def colorbar_wrapper(self, mappable, values=None,
    extend=None, extendsize=None,
    title=None, label=None,
    grid=None, tickminor=None,
    tickloc=None, ticklocation=None,
    locator=None, ticks=None, minorlocator=None, minorticks=None, locator_kw={}, minorlocator_kw={},
    formatter=None, ticklabels=None, formatter_kw={},
    fixticks=False, norm=None, norm_kw={}, # normalizer to use when passing colors/lines
    orientation='horizontal',
    **kwargs):
    """
    Wraps `~proplot.axes.BaseAxes` `~proplot.axes.BaseAxes.colorbar` and
    `~proplot.axes.PanelAxes` `~proplot.axes.PanelAxes.colorbar`, adds some
    handy features.

    Parameters
    ----------
    mappable : mappable or list of str or list of plot handles
        There are three options here:

        1. A mappable object. Basically, any object with a ``get_cmap`` method,
           like the objects returned by `~matplotlib.axes.Axes.contourf` and
           `~matplotlib.axes.Axes.pcolormesh`.
        2. A list of hex strings, color string names, or RGB tuples. From this,
           a colormap will be generated and used with the colorbar. Requires
           `values` is not ``None``.
        3. A list of "plot handles". Basically, any object with a ``get_color``
           method, like `~matplotlib.lines.Line2D` instances. From this,
           a colormap will be generated and used with the colorbar. If `values`
           is ``None``, they will try to be inferred by converting the handle
           labels returned by `~matplotlib.artist.Artist.get_label` to `float`.

    values : None or list of float, optional
        Ignored if `mappable` is a mappable object. Maps each color or plot
        handle in the `mappable` list to numeric values. From this, a
        colormap and normalizer are constructed.
    extend : {None, 'neither', 'both', 'min', 'max'}, optional
        Direction for drawing colorbar "extensions" (i.e. references to
        out-of-bounds data with a unique color). These are triangles by
        default. If ``None``, we try to use the ``extend`` attribute on the
        mappable object. If the attribute is unavailable, we use ``'neither'``.
    extendsize : None or float or str, optional
        The length of the colorbar "extensions" in *physical units*.
        If float, units are inches. If string, units are interpreted
        by `~proplot.utils.units`. Default is ``rc['colorbar.extend']``.

        This is handy if you have multiple colorbars in one figure.
        With the matplotlib API, it is really hard to get triangle
        sizes to match, because the `extendsize` units are *relative*.
    tickloc, ticklocation : {'bottom', 'top', 'left', 'right'}, optional
        Where to draw tick marks on the colorbar.
    label, title : None or str, optional
        The colorbar label. The `title` keyword is also accepted for
        consistency with `legend_wrapper`.
    grid : None or bool, optional
        Whether to draw "gridlines" between each level of the colorbar.
        Default is ``rc['colorbar.grid']``.
    tickminor : bool, optional
        Whether to put minor ticks on the colorbar. Default is ``False``.
    locator : None or locator spec, optional
        The colorbar tick mark positions. Passed to the
        `~proplot.axistools.Locator` constructor.
    locator_kw : dict-like, optional
        The locator settings. Passed to `~proplot.axistools.Locator`.
    minorlocator
        As with `locator`, but for the minor tick marks.
    minorlocator_kw
        As for `locator_kw`, but for the minor locator.
    formatter : None or formatter spec, optional
        The tick label format. Passed to the `~proplot.axistools.Formatter`
        constructor.
    formatter_kw : dict-like, optional
        The formatter settings. Passed to `~proplot.axistools.Formatter`.
    fixticks : bool, optional
        For complicated normalizers (e.g. `~matplotlib.colors.LogNorm`), the
        colorbar minor and major ticks can appear misaligned. When `fixticks`
        is ``True``, this misalignment is fixed. The default is ``False``.

        This will give incorrect positions when the colormap index does not
        appear to vary "linearly" from left-to-right across the colorbar (for
        example, when the leftmost colormap colors seem to be "pulled" to the
        right farther than normal). In this case, you should stick with
        ``fixticks=False``.
    norm : None or normalizer spec, optional
        Ignored if `values` is ``None``. The normalizer
        for converting `values` to colormap colors. Passed to the
        `~proplot.colortools.Norm` constructor. As an example, if your
        values are logarithmically spaced but you want the level boundaries
        to appear halfway in-between the colorbar tick marks, try
        ``norm='log'``.
    norm_kw : dict-like, optional
        The normalizer settings. Passed to `~proplot.colortools.Norm`.
    orientation : {'horizontal', 'vertical'}, optional
        The colorbar orientation. Generally, you shouldn't have to explicitly
        set this.

    Other parameters
    ----------------
    **kwargs
        Passed to `~matplotlib.figure.Figure.colorbar`.

    See also
    --------
    `~proplot.axes.BaseAxes.colorbar`, `~proplot.axes.PanelAxes.colorbar`, `~matplotlib.figure.Figure.colorbar`,
    `~proplot.axistools.Locator`, `~proplot.axistools.Formatter`, `~proplot.colortools.Norm`
    """
    # Developer notes
    # * Colorbar axes must be of type `matplotlib.axes.Axes`,
    #   not `~proplot.axes.BaseAxes`, because colorbar uses some internal methods
    #   that are wrapped by `~proplot.axes.BaseAxes`.
    # * There is an insanely weird problem with colorbars when simultaneously
    #   passing levels and norm object to a mappable; fixed by passing
    #   vmin/vmax instead of levels.
    #   (see: https://stackoverflow.com/q/40116968/4970632).
    # * Problem is, often want levels instead of vmin/vmax, while simultaneously
    #   using a Normalize (for example) to determine colors between the levels
    #   (see: https://stackoverflow.com/q/42723538/4970632). Workaround is to
    #   make sure locators are in vmin/vmax range exclusively; cannot match/exceed values.
    # Parse flexible input
    label = _default(title, label)
    locator = _default(ticks, locator)
    formatter = _default(ticklabels, formatter, 'default')
    minorlocator = _default(minorticks, minorlocator)
    ticklocation = _default(tickloc, ticklocation)
    # Apply user-kwargs
    # WARNING: PathCollection scatter objects have an extend method!
    if extend is None:
        if isinstance(getattr(mappable, 'extend', None), str):
            extend = mappable.extend or 'neither'
        else:
            extend = 'neither'
    grid = _default(grid, rc['colorbar.grid'])
    kwdefault = {'cax':self, 'orientation':orientation, 'use_gridspec':True, # use space afforded by entire axes
                 'spacing':'uniform', 'extend':extend, 'drawedges':grid} # this is default case unless mappable has special props
    kwdefault.update(kwargs)
    kwargs = kwdefault

    # Special case where auto colorbar is generated from 1d methods, a list is
    # always passed but some 1d methods (scatter) do have colormaps.
    if np.iterable(mappable) and len(mappable)==1 and hasattr(mappable[0], 'get_cmap'):
        mappable = mappable[0]
    # Test if we were given a mappable, or iterable of stuff; note Container and
    # PolyCollection matplotlib classes are iterable.
    cmap = None
    if not isinstance(mappable, martist.Artist) and not isinstance(mappable, mcontour.ContourSet):
        # Get object for testing
        if not np.iterable(mappable):
            mappable = [None] # raises error below
        obj = mappable[0]
        try:
            obj = obj[0] # e.g. for BarContainer, which is not numpy.iterable
        except (TypeError,KeyError):
            pass
        # Draw from handles
        if hasattr(obj, 'get_color') or hasattr(obj, 'get_facecolor'): # simplest approach
            colors = []
            for obj in mappable:
                if np.iterable(obj):
                    obj = obj[0]
                color = getattr(obj, 'get_color', None) or getattr(obj, 'get_facecolor')
                colors.append(color())
            cmap = colortools.Colormap(colors)
            if values is None:
                values = []
                for obj in mappable:
                    val = obj.get_label()
                    try:
                        val = float(val)
                    except ValueError:
                        raise ValueError(f'To generate colorbars from line handles or other objects, pass the "values" keyword arg to colorbar(), or give your handles numeric values with e.g. plot(..., label=123) or line.set_label(123).')
                    values.append(val)
        # List of colors
        elif all(isinstance(obj, str) for obj in mappable) or all(np.iterable(obj) and len(obj) in (3,4) for obj in mappable):
            cmap = colortools.Colormap(mappable)
            if values is None:
                raise ValueError(f'To generate colorbars from lists of colors, pass the "values" keyword arg to colorbar().')
        # Invalid
        else:
            raise ValueError(f'Input mappable must be a matplotlib artist, list of objects, or list of colors. Got {mappable}.')
    # Build new ad hoc mappable object from handles
    if cmap is not None:
        func = _cmap_wrapper(self, self.contourf)
        mappable = func([[0,0],[0,0]],
            cmap=cmap, extend='neither', values=np.array(values),
            norm=norm, norm_kw=norm_kw) # workaround
        if locator is None:
            nstep = 1 + len(values)//20
            locator = values[::nstep]

    # By default, label the discretization levels (if there aren't too many)
    # Prefer centers (i.e. 'values') to edges (i.e. 'levels')
    if locator is None:
        locator = getattr(mappable, 'values', getattr(mappable, 'levels', None))
        if locator is not None:
            step = 1 + len(locator)//20
            locator = locator[::step]

    # Determine major formatters and major/minor tick locators
    # Can pass locator/minorlocator as the *jump values* between the mappables
    # vmin/vmax if desired
    ivalues = None # so linter doesn't detect error in if i==1 block
    normfix = False # whether we need to modify the norm object
    locators = [] # put them here
    for i,(ilocator,ilocator_kw) in enumerate(zip((locator,minorlocator), (locator_kw,minorlocator_kw))):
        # Get the locator values
        # Need to use tick_values instead of accessing 'locs' attribute because
        # many locators don't have these attributes; require norm.vmin/vmax as input
        if i==1 and (not tickminor and ilocator is None): # means we never wanted minor ticks
            locators.append(axistools.Locator('null'))
            continue
        ilocator = _default(ilocator, 'auto')
        jvalues = np.array(axistools.Locator(ilocator, **ilocator_kw).tick_values(mappable.norm.vmin, mappable.norm.vmax)) # get the current values
        # Modify ticks to work around mysterious error, and to prevent annoyance
        # where minor ticks extend beyond extendsize.
        # We need to figure out the numbers that will eventually be rendered to
        # solve the error, so we will always use a fixedlocator.
        values_min = np.where(jvalues>=mappable.norm.vmin)[0]
        values_max = np.where(jvalues<=mappable.norm.vmax)[0]
        if len(values_min)==0 or len(values_max)==0:
            locators.append(axistools.Locator('null'))
            continue
        values_min, values_max = values_min[0], values_max[-1]
        jvalues = jvalues[values_min:values_max+1]
        if jvalues[0]==mappable.norm.vmin:
            normfix = True
        # Prevent annoying major/minor overlaps where one is slightly shifted left/right
        # Consider floating point weirdness too
        if i==1:
            eps = 1e-10
            jvalues = [v for v in jvalues if not any(o+eps >= v >= o-eps for o in ivalues)]
        ivalues = jvalues # record as new variable
        locators.append(axistools.Locator(ivalues)) # final locator object

    # Fix the norm object; get weird error without this block
    # * The error is triggered when a *major* tick sits exactly on vmin, but
    #   the actual error is due to processing of *minor* ticks, even if the 
    #   minor locator was set to NullLocator; very weird. Happens when we call
    #   get_ticklabels(which='both') below. Can be prevented by just calling
    #   which='major'. Minor ticklabels are never drawn anyway.
    # * We can eliminate the normfix below, but that actually causes an annoying
    #   warning to be printed (related to same issue I guess). So we keep this.
    #   The culprit for all of this seems to be the colorbar API line:
    #        z = np.take(y, i0) + (xn - np.take(b, i0)) * dy / db
    #   Also strange that minorticks extending *below* the minimum
    #   don't raise the error. It is only when they are exactly on the minimum.
    # * When changing the levels attribute, need to make sure the levels
    #   datatype is float; otherwise division will be truncated and bottom
    #   level will still lie on same location, so error will occur
    if normfix:
        mappable.norm.vmin -= (mappable.norm.vmax-mappable.norm.vmin)/10000
    if hasattr(mappable.norm, 'levels'):
        mappable.norm.levels = np.atleast_1d(mappable.norm.levels).astype(np.float)
        if normfix:
            mappable.norm.levels[0] -= np.diff(mappable.norm.levels[:2])[0]/10000

    # Final settings
    # NOTE: The only way to avoid bugs seems to be to pass the major formatter
    # and locator to colorbar commmand directly, but edit the minor locators
    # and formatters manually; set_locator methods are completely ignored.
    width, height = self.figure.get_size_inches()
    formatter = axistools.Formatter(formatter, **formatter_kw)
    if orientation=='horizontal':
        scale = width*abs(self.get_position().width)
    else:
        scale = height*abs(self.get_position().height)
    extendsize = utils.units(_default(extendsize, rc['colorbar.extend']))
    extendsize = extendsize/(scale - 2*extendsize)
    kwargs.update({'ticks':locators[0], # WARNING: without this, set_ticks screws up number labels for some reason
                   'format':formatter,
                   'ticklocation':ticklocation,
                   'extendfrac':extendsize})
    # Draw the colorbar
    cb = self.figure.colorbar(mappable, **kwargs)
    # Make edges/dividers style consistent with gridline style
    if cb.dividers is not None:
        cb.dividers.update(rc['grid'])

    # The minor locators and formatters
    # NOTE: Re-apply major ticks here because for some reason minor ticks don't
    # align with major ones for LogNorm. When we call set_ticks, labels and
    # numbers are not changed; just re-adjust existing ticks to proper locations.
    minorvals = np.array(locators[1].tick_values(mappable.norm.vmin, mappable.norm.vmax))
    majorvals = np.array(locators[0].tick_values(mappable.norm.vmin, mappable.norm.vmax))
    if isinstance(mappable.norm, colortools.BinNorm):
        minorvals = mappable.norm._norm(minorvals) # use *child* normalizer
        majorvals = mappable.norm._norm(majorvals)
    else:
        minorvals = mappable.norm(minorvals)
        majorvals = mappable.norm(majorvals) # use *child* normalizer
    minorvals = [tick for tick in minorvals if 0<=tick<=1]
    majorvals = [tick for tick in majorvals if 0<=tick<=1]
    if orientation=='horizontal':
        axis = self.xaxis
    else:
        axis = self.yaxis
    if fixticks:
        axis.set_ticks(majorvals, minor=False)
    axis.set_ticks(minorvals, minor=True)
    axis.set_minor_formatter(mticker.NullFormatter()) # to make sure
    if label is not None:
        axis.label.update({'text':label})

    # Fix alpha issues. Cannot set edgecolor to 'face' if alpha non-zero
    # because blending will occur, will get colored lines instead of white ones;
    # need to perform manual alpha blending.
    # NOTE: For some reason cb solids uses listed colormap with always 1.0
    # alpha, then alpha is applied after.
    # See: https://stackoverflow.com/a/35672224/4970632
    alpha = None
    if cb.solids: # for e.g. contours with colormap, colorbar will just be lines
        alpha = cb.solids.get_alpha()
    if alpha is not None and alpha<1:
        # First get reference color
        warnings.warn('Performing manual alpha-blending for colorbar solids.')
        reference = mappable.axes.get_facecolor() # the axes facecolor
        reference = [(1 - reference[-1]) + reference[-1]*color for color in reference[:3]]
        # Next get solids
        reference = [1,1,1] # override?
        alpha = 1 - (1 - alpha)**2 # make more colorful
        colors = cb.solids.get_cmap().colors
        colors = np.array(colors)
        for i in range(3): # Do not include the last column!
            colors[:,i] = (reference[i] - alpha) + alpha*colors[:,i]
        cmap = mcolors.ListedColormap(colors, name='colorbar-fix')
        cb.solids.set_cmap(cmap)
        cb.solids.set_alpha(1.0)

    # Fix pesky white lines between levels + misalignment with border due
    # to rasterized blocks.
    if cb.solids:
        cb.solids.set_linewidth(0.4) # lowest size that works
        cb.solids.set_edgecolor('face')
        cb.solids.set_rasterized(False)
    axis.set_ticks_position(ticklocation)
    return cb

#------------------------------------------------------------------------------#
# Now that we have defined all wrappers in a way that will make nice documentation/
# minimize the amount of new code that has to be generated when a method is
# wrapped inside __getattribute__, construct the *actual* wrappers.
#------------------------------------------------------------------------------#
# Helper funcs
def _wrapper(driver):
    def decorator(self, func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return driver(self, func, *args, **kwargs)
        return wrapper
    return decorator
def _simple_wrapper(driver):
    def decorator(self, func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return driver(self, *args, **kwargs)
        return wrapper
    return decorator
# Hidden wrappers
# Also _basemap_call and _general_norecurse
_autoformat_1d_ = _wrapper(_autoformat_1d)
_autoformat_2d_ = _wrapper(_autoformat_2d)
# Documented
_enforce_centers       = _wrapper(enforce_centers)
_enforce_edges         = _wrapper(enforce_edges)
_basemap_gridfix       = _wrapper(basemap_gridfix)
_basemap_latlon        = _wrapper(basemap_latlon)
_cartopy_gridfix       = _wrapper(cartopy_gridfix)
_cartopy_transform     = _wrapper(cartopy_transform)
_cartopy_crs           = _wrapper(cartopy_crs)
_cmap_wrapper          = _wrapper(cmap_wrapper)
_cycle_wrapper         = _wrapper(cycle_wrapper)
_bar_wrapper           = _wrapper(bar_wrapper)
_barh_wrapper          = _wrapper(barh_wrapper)
_plot_wrapper          = _wrapper(plot_wrapper)
_scatter_wrapper       = _wrapper(scatter_wrapper)
_boxplot_wrapper       = _wrapper(boxplot_wrapper)
_violinplot_wrapper    = _wrapper(violinplot_wrapper)
_fill_between_wrapper  = _wrapper(fill_between_wrapper)
_fill_betweenx_wrapper = _wrapper(fill_betweenx_wrapper)
_text_wrapper          = _wrapper(text_wrapper)
_legend_wrapper        = _simple_wrapper(legend_wrapper) # special

