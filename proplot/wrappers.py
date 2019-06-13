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
import matplotlib.axes as maxes
import matplotlib.contour as mcontour
import matplotlib.ticker as mticker
import matplotlib.transforms as mtransforms
import matplotlib.patheffects as mpatheffects
import matplotlib.colors as mcolors
import matplotlib.artist as martist
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

# Keywords
# These may automatically override the 'fix' option!
_options = {
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

# Methods for wrapping
# TODO: 'quiver', 'streamplot' for cmap?
_centers_methods = ('contour', 'contourf', 'quiver', 'streamplot', 'barbs')
_edges_methods = ('pcolor', 'pcolormesh',)
_2d_methods = (*_centers_methods, *_edges_methods)
_1d_methods = ('plot', 'scatter', 'bar', 'barh')
_cycle_methods  = ('plot', 'scatter', 'bar', 'barh')
_cmap_methods = ('contour', 'contourf', 'pcolor', 'pcolormesh',
    'tripcolor', 'tricontour', 'tricontourf',
    'cmapline', 'hexbin', 'matshow', 'imshow', 'spy', 'hist2d',)
_latlon_methods = ('plot', 'scatter',) # adds latlon=True
_transform_methods = ('plot', 'scatter', 'tripcolor', 'tricontour', 'tricontourf',) # adds transform=PlateCarree()
_crs_methods = ('get_extent', 'set_extent', 'set_xticks', 'set_yticks',) # adds crs=PlateCarree()

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

#------------------------------------------------------------------------------#
# For documentation
#------------------------------------------------------------------------------#
def _sphinx_name(name):
    """Gets sphinx name."""
    if name=='cmapline':
        return f'`~BaseAxes.{name}`'
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
        ('_latlon_methods',       _latlon_methods),
        ('_crs_methods',          _crs_methods),
        ('_transform_methods',    _transform_methods),
        ('_cycle_methods',        _cycle_methods),
        ('_cmap_methods',         _cmap_methods),
        ('_disabled_methods',     (*(method for methods in _disabled_methods.values() for method in methods),)),
        ('_map_disabled_methods', _map_disabled_methods),
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
def _array_std(data):
    """Converts list of lists to array, but no other input."""
    # First convert to array
    if not isinstance(data, (ndarray, DataArray, DataFrame, Series, Index)):
        data = np.array(data)
    return data

def _auto_label(data, units=True):
    """Gets label from pandas or xarray objects."""
    if isinstance(data, ndarray):
        return ''
    # Xarray with common NetCDF attribute names
    elif isinstance(data, DataArray):
        label = getattr(data, 'name', '') or ''
        for key in ('standard_name', 'long_name'):
            label = data.attrs.get(key, label)
        if units:
            units = data.attrs.get('units', '')
            if label and units:
                label = f'{label} ({units})'
            elif units:
                label = units
    # Pandas, account for common situation with 1-column DataFrame where the
    # column label is descriptor for the data
    elif isinstance(data, (DataFrame, Series, Index)):
        label = getattr(data, 'name', '') or '' # DataFrame has no native name attribute but user can add one: https://github.com/pandas-dev/pandas/issues/447
        if not label and isinstance(data, DataFrame) and data.columns.size==1:
            label = str(df.columns[0])
    return str(label).strip()

def _parse_1d(self, func, *args, **kwargs):
    """Accepts 1d DataArray or Series, or
    2D DataArray or DataFrame, in which case list of lines or points
    are drawn. Used by `plot_wrapper` and `scatter_wrapper`."""
    # Sanitize input
    if len(args)==1:
        x = None
        y, *args = args
    elif len(args) in (2,3,4):
        x, y, *args = args # same
    else:
        raise ValueError(f'Passed {len(args)} arguments to plotting command. Only 1-4 are valid.')
    # Detect 1d
    is1d = True # i.e. args is a 1d vector
    if not np.iterable(y):
        raise ValueError(f'Invalid y data {y}.')
    elif getattr(y, 'ndim', None)==2: # e.g. DataFrame[0] can raise error, indexing is by column name, so this test is best
        is1d = False
    elif np.iterable(y[0]): # e.g. list of lists
        is1d = False
    # Ensure 2d and draw sample column
    y = _array_std(y)
    if y.ndim==2:
        iy = getattr(y, 'iloc', y)
        iy = iy[:,0]
    elif y.ndim==1:
        iy = y
    else:
        raise ValueError(f'y must be 1 or 2-dimensional, got shape {y.shape}.')
    # Auto coords
    if x is None:
        if isinstance(iy, ndarray):
            x = np.arange(iy.size)
        elif isinstance(iy, DataArray): # DataArray
            x = iy.coords[iy.dims[0]]
        elif isinstance(iy, Series): # Series
            x = iy.index
        else: # Index
            raise ValueError(f'Unable to infer x coordinates from pandas.Index-type y coordinates.')
    # Check coordinates
    x = _array_std(x)
    if x.ndim!=1:
        raise ValueError(f'x coordinates must be 1-dimensional, but got {x.ndim}.')
    # Auto formatting
    if self.figure._autoformat:
        kw = {}
        iy = None
        if not self._is_map:
            # Xlabel
            label = _auto_label(x)
            if label:
                kw['xlabel'] = label
            if all(isinstance(x, Number) for x in x[:2]) and x[1]<x[0]:
                kw['xreverse'] = True
            # Ylabel; only try if the input was 1d, otherwise info is title
            if y.ndim==2 and y.shape[1]==1:
                iy = getattr(y, 'iloc', y)
                iy = iy[:,0]
            elif y.ndim==1:
                iy = y
            if iy is not None:
                label = _auto_label(iy)
                if label:
                    kw['ylabel'] = label
        # Title
        if iy is None:
            label = _auto_label(y)
            if label:
                kw['title'] = label
        self.format(**kw)
    return func(x, y, *args, **kwargs)

def _parse_2d(self, func, *args, order='C', **kwargs):
    """Gets 2d data. Accepts ndarray and DataArray. Used by `check_centers`
    and `check_edges`, which are used for all 2d plot methods."""
    # Sanitize input
    if len(args)>4:
        raise ValueError(f'Too many arguments passed to plotting command. Max is 4.')
    elif len(args)<1:
        raise ValueError(f'Too few arguments passed to plotting command. Min is 1.')
    x, y = None, None
    if len(args)>2:
        x, y, *args = args
    # Ensure DataArray, DataFrame or ndarray
    # NOTE: All of these have shape and ndim attributes
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
        if order=='C':
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
        # Labels
        if not self._is_map:
            for key,z in zip(('xlabel','ylabel'), (x,y)):
                label = _auto_label(z)
                if label:
                    kw[key] = label
                if z[1]<z[0]:
                    kw[key[0] + 'reverse'] = True
        # Title
        title = _auto_label(Zs[0], units=False)
        if title:
            kw['title'] = title
        self.format(**kw)
    return func(x, y, *Zs, **kwargs)

#------------------------------------------------------------------------------
# 2D plot wrappers
#------------------------------------------------------------------------------
@_expand_methods_list
def check_centers(self, func, *args, order='C', **kwargs):
    """
    Wraps 2D plotting functions that take coordinate *centers* (`_centers_methods`),
    calculates centers if graticule *edges* were provided.

    Note
    ----
    Optional numbers of arguments:

    * Z
    * U, V
    * x, y, Z
    * x, y, U, V
    """
    # Checks whether sizes match up, checks whether graticule was input
    x, y, *Zs = args
    xlen, ylen = x.shape[-1], y.shape[0]
    for Z in Zs:
        if Z.ndim!=2:
            raise ValueError(f'Input arrays must be 2D, instead got shape {Z.shape}.')
        elif Z.shape[1]==xlen-1 and Z.shape[0]==ylen-1 and x.ndim==1 and y.ndim==1:
            x = (x[1:] + x[:-1])/2
            y = (y[1:] + y[:-1])/2 # get centers, given edges
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
def check_edges(self, func, *args, order='C', **kwargs):
    """
    Wraps 2D plotting functions that take graticule *edges* (`_edges_methods`),
    calculates edges if coordinate *centers* were provided.

    Note
    ----
    Optional numbers of arguments:

    * Z
    * U, V
    * x, y, Z
    * x, y, U, V
    """
    # Checks that sizes match up, checks whether graticule was input
    # return func(*args, **kwargs)
    x, y, *Zs = args
    xlen, ylen = x.shape[-1], y.shape[0]
    for Z in Zs:
        if Z.ndim!=2:
            raise ValueError(f'Input arrays must be 2D, instead got shape {Z.shape}.')
        elif Z.shape[1]==xlen and Z.shape[0]==ylen:
            # If 2D, don't raise error, but don't fix either, because
            # matplotlib pcolor accepts grid center inputs.
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
    Wraps `~matplotlib.axes.Axes.plot`, calls `~BaseAxes.cmapline`
    if ``cmap`` is passed by the user.

    Parameters
    ----------
    *args
        Passed to `~matplotlib.axes.Axes.plot`.
    cmap, values
        Passed to `~BaseAxes.cmapline`.
    **kwargs
        `~matplotlib.lines.Line2D` properties.
    """
    # Parse input
    if len(args) not in (2,3): # e.g. with fmt string
        raise ValueError(f'Expected 1-3 plot args, got {len(args)}')
    # Make normal boring lines
    if cmap is None:
        lines = func(*args, **kwargs)
    # Make special colormap lines
    else:
        lines = self.cmapline(*args, cmap=cmap, values=values, **kwargs)
    return lines

def scatter_wrapper(self, func, *args,
    c=None, color=None, markercolor=None,
    s=None, size=None, markersize=None,
    lw=None, linewidth=None, linewidths=None, markeredgewidth=None, markeredgewidths=None,
    edgecolor=None, edgecolors=None, markeredgecolor=None, markeredgecolors=None,
    **kwargs):
    """
    Wraps `~matplotlib.axes.Axes.scatter`, adds optional keyword args
    more consistent with the `~matplotlib.axes.Axes.plot` keywords.

    Parameters
    ----------
    c, color, markercolor : None or str or (R,G,B) tuple, or list thereof, optional
        Aliases for the marker fill color.
    s, size, markersize : None or float, or list thereof, optional
        Aliases for the marker size.
    lw, linewidth, linewidths, markeredgewidth, markeredgewidths : None or float, or list thereof, optional
        Aliases for the marker edge width.
    edgecolors, markeredgecolor, markeredgecolors : None or str or (R,G,B) tuple, or list thereof, optional
        Aliases for the marker edge color.
    **kwargs
        Passed to `~matplotlib.axes.Axes.scatter`.
    """
    # Manage input arguments
    args = [*args] # convert to list
    if len(args)>3:
        c = args.pop(3)
    if len(args)>2:
        s = args.pop(2)
    # Apply some aliases for keyword arguments
    c = _default(c, color, markercolor)
    s = _default(s, size, markersize)
    lws = _default(lw, linewidths, linewidth, markeredgewidths, markeredgewidth)
    ecs = _default(edgecolors, edgecolor, markeredgecolors, markeredgecolor)
    return func(*args, c=c, s=s, linewidths=lws, edgecolors=ecs, **kwargs)

#------------------------------------------------------------------------------#
# Text wrapper
#------------------------------------------------------------------------------#
def text_wrapper(self, func, x, y, text,
    transform='data',
    border=False, border_kw={},
    invert=False,
    linewidth=2, lw=None,
    **kwargs): # linewidth is for the border
    """
    Wraps `~matplotlib.axes.Axes.text`, enables specifying `tranform` with
    a string name and adds feature for drawing borders around text.

    Parameters
    ----------
    x, y : float
        The *x* and *y* coordinates for the text.
    text : str
        The text.
    linewidth, lw : float, optional
        Ignored if `border` is ``False``. The width of the text border.
    border : bool, optional
        Whether to draw border around text.
    border_kw : dict-like, optional
        Passed to `~matplotlib.path_effects.Stroke` if drawing a border.
    invert : bool, optional
        Ignored if `border` is ``False``. Whether to draw black text
        with a white border (``False``), or white text on a black
        border (``True``).
    transform : {'data', 'axes', 'figure'} or `~matplotlib.transforms.Transform`, optional
        The transform, or a string pointing to either of the
        `~matplotlib.axes.Axes.transData`, `~matplotlib.axes.Axes.transAxes`,
        or `~matplotlib.figure.Figure.transFigure` transforms. Default is
        ``'data'`` (unchanged).

    Other parameters
    ----------------
    **kwargs
        Passed to `~matplotlib.text.Text` instantiator.
    """
    # Get default transform by string name
    # Note basemap gridlining methods call text, so if you change the
    # default transform, you will not be able to draw latitude and
    # longitude labels! Leave it alone.
    linewidth = lw or linewidth
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
    # Raise more helpful error message if font unavailable
    name = kwargs.pop('fontname', rc['fontname']) # is actually font.sans-serif
    if name not in fonttools.fonts:
        suffix = ''
        if name not in fonttools._missing_fonts:
            suffix = f' Available fonts are: {", ".join(fonttools.fonts)}.'
        warnings.warn(f'Font "{name}" unavailable, falling back to DejaVu Sans. Consider running proplot.install_fonts().' + suffix)
        fonttools._missing_fonts.append(name)
        name = 'DejaVu Sans'
    # Call parent, with custom rc settings
    # These seem to sometimes not get used by default
    size   = kwargs.pop('fontsize', rc.get('font.size'))
    color  = kwargs.pop('color', rc.get('text.color'))
    weight = kwargs.pop('font', rc.get('font.weight'))
    obj = func(x, y, text,
        transform=transform, fontname=name,
        fontsize=size, color=color, fontweight=weight, **kwargs)
    # Optionally draw border around text
    if border:
        facecolor, bgcolor = color, 'w'
        if invert:
            facecolor, bgcolor = bgcolor, facecolor
        kwargs = {'linewidth':linewidth, 'foreground':bgcolor, 'joinstyle':'miter'}
        kwargs.update(border_kw)
        obj.update({'color':facecolor, 'zorder':1e10, # have to update after-the-fact for path effects
            'path_effects': [mpatheffects.Stroke(**kwargs), mpatheffects.Normal()]})
    return obj

#------------------------------------------------------------------------------#
# Geographic wrappers
#------------------------------------------------------------------------------#
# First basemap recursion fix
# Normally we *cannot* modify the underlying *axes* pcolormesh etc. because this
# this will cause basemap's self.m.pcolormesh etc. to use *custom* version and
# cause suite of weird errors. Prevent this recursion with the below decorator.
def _m_call(self, func):
    """Docorator that calls the basemap version of the function of the same name."""
    name = func.__name__
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return self.m.__getattribute__(name)(ax=self, *args, **kwargs)
    return wrapper

def _m_norecurse(self, func):
    """Decorator to prevent recursion in Basemap method overrides.
    See `this post https://stackoverflow.com/a/37675810/4970632`__."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        name = getattr(func, '__name__')
        if self._recurred:
            # Return the *original* version of the matplotlib method (i.e.
            # the one we have not wrapped by overriding the __getattribute__
            # method). We reach this block when basemap.Basemap tries to
            # call the original method internally.
            self._recurred = False
            result = object.__getattribute__(self, name)(*args, **kwargs)
        else:
            # Return the version we have wrapped, which itself will call the
            # basemap.Basemap method thanks to the _m_call wrapper.
            self._recurred = True
            result = func(*args, **kwargs)
        self._recurred = False # cleanup, in case recursion never occurred
        return result
    return wrapper

@_expand_methods_list
def basemap_latlon(self, func, *args, latlon=True, **kwargs):
    """
    Wraps plotting functions for `BasemapAxes` (`_latlon_methods`).

    With the default `~mpl_toolkits.basemap` API, you need to pass
    ``latlon=True`` if your data coordinates are longitude and latitude,
    instead of map projection coordinates. Now, ``latlon=True`` is always
    used.
    """
    return func(*args, latlon=latlon, **kwargs)

@_expand_methods_list
def cartopy_transform(self, func, *args, transform=PlateCarree, **kwargs):
    """
    Wraps plotting functions for `CartopyAxes` (`_transform_methods`).

    With the default `~cartopy.mpl.geoaxes.GeoAxes` API, you need to pass
    ``transform=cartopy.crs.PlateCarree()`` if your data coordinates are
    longitude and latitude, instead of map projection coordinates.
    Now, ``transform=cartopy.crs.PlateCarree()`` is the default behavior.
    """
    # Simple
    if isinstance(transform, type):
        transform = transform() # instantiate
    result = func(*args, transform=transform, **kwargs)
    # Re-enforce settings because some plot functions seem to reset the
    # outlinepatch or backgroundpatch (???)
    # TODO: Double check this
    self.format()
    return result

@_expand_methods_list
def cartopy_crs(self, func, *args, crs=PlateCarree, **kwargs):
    """
    Wraps axes functions for `CartopyAxes` (`_crs_methods`).

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
def cartopy_gridfix(self, func, lon, lat, *Zs, transform=PlateCarree, globe=False, **kwargs):
    """
    Wraps 2D plotting functions for `CartopyAxes` (`_centers_edges_methods`).

    As in `cartopy_transform`, but adds the `globe` keyword arg
    to optionally make data coverage *global*. Passing ``globe=True``
    does the following:

    1. Makes longitudinal coverage *circular* (i.e. the last
       longitude coordinate equals the first longitude coordinate plus 360
       degrees).
    2. Interpolates data to the North and South poles.

    Warning
    -------
    Cartopy contouring methods have issues with circularly wrapped data.
    Triggers annoying ``TopologyException`` statements, which we suppress
    with the IPython `~IPython.utils.io.capture_output` tool.
    This is a workaround. See `this issue
    <https://github.com/SciTools/cartopy/issues/946>`__.
    """
    # Ensure monotonic lons or things get messed up
    # Unlike basemap can be monotonic from any starting point
    # WARNING: In first geophysical data example, got persistent
    # changes to input data arrays without below line, resulting in non
    # monotonic coordinates and discovery of this error
    lon, lat = np.array(lon), np.array(lat) # no need to retain metadata on e.g. DataArray
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
    # with io.capture_output() as captured:
    if isinstance(transform, type):
        transform = transform() # instantiate
    result = func(lon, lat, *Zss, transform=transform, **kwargs)
    # Re-enforce settings because some plot functions seem to reset the
    # outlinepatch or backgroundpatch (???)
    # TODO: Double check this
    self.format()
    return result

@_expand_methods_list
def basemap_gridfix(self, func, lon, lat, *Zs, globe=False, latlon=True, **kwargs):
    """
    Wraps 2D plotting functions for `BasemapAxes` (`_centers_edges_methods`).

    As in `basemap_latlon`, but cycles longitudes to fit within the
    map edges (i.e. if the projection central longitude is 90 degrees, will
    permute data to span from -90 degrees to 270 degrees longitude).

    Also adds the `globe` keyword arg to optionally make data coverage *global*.
    Passing ``globe=True`` does the following:

    1. Makes longitudinal coverage *circular* (i.e. the last
       longitude coordinate equals the first longitude coordinate plus 360
       degrees).
    2. Interpolates data to the North and South poles.
    """
    # Bail out if map coordinates already provided
    lon, lat = np.array(lon), np.array(lat) # no need to retain metadata on e.g. DataArray
    if not latlon:
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
            # if self.m.projection[:4] != 'merc': # did not fix the problem where Mercator goes way too far
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
    x, y = self.m(*np.meshgrid(lon, lat))
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
        label=None, labels=None, values=None,
        legend=None, legend_kw={},
        colorbar=None, colorbar_kw={},
        **kwargs):
    """
    Wraps methods that use the property cycler (`_cycle_methods`),
    adds features for controlling colors in the property cycler.

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

        Otherwise, an inset legend is drawn, and this sets the position.
        The following position keys and position key aliases are valid:

        ==================  ==============================================
        Position            Valid aliases
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
    `BaseAxes`, `~proplot.colortools.Cycle`

    Note
    ----
    See the `matplotlib source 
    <https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/axes/_base.py>`_.
    The `set_prop_cycle` command modifies underlying 
    `_get_lines` and `_get_patches_for_fill`.
    """
    # Test input
    # NOTE: Since _plot_wrapper and _scatter_wrapper come before this, input
    # will be standardized, so below is valid.
    x, y, *args = args
    is1d = (y.ndim==1)

    # Determine and temporarily set cycler
    # WARNING: Axes cycle has no getter, only setter (set_prop_cycle), which
    # sets a 'prop_cycler' attribute on the hidden _get_lines and
    # _get_patches_for_fill objects. Only way to query current axes cycler!
    # NOTE: Should not wrap set_prop_cycle because it calls several secondary
    # methods, would get messy and fragile.
    if cycle is not None:
        # Get the new group of colors
        if not np.iterable(cycle) or isinstance(cycle, str):
            cycle = cycle,
        if not is1d and y.shape[1]>1: # apply new default
            cycle_kw['samples'] = y.shape[1]
        cycle = colortools.Cycle(*cycle, **cycle_kw)
        # Compare to the original group of colors, reset if different
        # NOTE: The _get_lines cycler is an *itertools cycler*. Has no length,
        # so we must cycle over it with next(). We try calling next() the same
        # number of times as the length of user input cycle. If the input cycle
        # *is* in fact the same, below does not reset the color position, cycles
        # us back to start!
        i = 0
        cycler = self._get_lines.prop_cycler
        cycle_orig = []
        while i<len(cycle):
            next_ = next(cycler)
            if 'color' in next_:
                cycle_orig.append(next_['color'])
            i += 1
        if {*cycle_orig} != {*cycle} or cycle_kw.get('shift', None): # order is immaterial
            self.set_prop_cycle(color=cycle)

    # Iterate
    labels = _default(values, labels, label, None)
    if isinstance(labels, str) or labels is None:
        labels = [labels]*(1 if is1d else y.shape[1])
    clabel = None
    objs = []
    y = getattr(y, 'iloc', y) # for indexing
    for i,label in enumerate(labels):
        iy = y if is1d else y[:,i]
        if not label:
            # Try to get it from coordiantes
            if isinstance(y, ndarray):
                pass
            elif isinstance(y, DataArray):
                label = y.coords[y.dims[1]].values[i]
                clabel = _auto_label(y.coords[y.dims[1]])
            elif isinstance(y, DataFrame):
                label = y.columns[i]
                clabel = _auto_label(y.columns)
            # Try just getting e.g. a pd.Series name
            else:
                label = _auto_label(iy)
        obj = func(x, iy, *args, label=label, **kwargs)
        if isinstance(obj, (list,tuple)): # plot always returns list or tuple
            obj = obj[0]
        objs.append(obj)

    # Add colorbar and/or legend
    if colorbar:
        ax, loc = _get_panel(self, colorbar)
        if loc not in ax._auto_colorbar:
            ax._auto_colorbar[loc] = []
            ax._auto_colorbar_kw[loc] = {}
        ax._auto_colorbar[loc].extend(objs)
        ax._auto_colorbar_kw[loc].update(colorbar_kw)
        if loc!='fill':
            ax._auto_colorbar_kw[loc].update({'loc':loc})
        if clabel:
            ax._auto_colorbar_kw[loc].update({'label':clabel})
    if legend:
        ax, loc = _get_panel(self, legend)
        if loc not in ax._auto_legend:
            ax._auto_legend[loc] = []
            ax._auto_legend_kw[loc] = {}
        ax._auto_legend[loc].extend(objs)
        ax._auto_legend_kw[loc].update(legend_kw)
        if loc!='fill':
            ax._auto_legend_kw[loc].update({'loc':loc})

    return objs[0] if is1d else objs # list of PathCollection or Line2D

@_expand_methods_list
def cmap_wrapper(self, func, *args, cmap=None, cmap_kw={},
    norm=None, norm_kw={},
    extend='neither',
    values=None, levels=None, zero=False, vmin=None, vmax=None, # override levels to be *centered* on zero
    edgefix=True, labels=False, labels_kw={}, precision=2,
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
        constructor. See `~proplot.colortools.Colormap` for options.
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
        range. Defaults to 11.

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
    edgefix : bool, optional
        Whether to fix the the `white-lines-between-filled-contours
        <https://stackoverflow.com/q/8263769/4970632>`__
        and `white-lines-between-pcolor-rectangles
        <https://stackoverflow.com/q/27092991/4970632>`__
        issues. Note this will slow down the figure rendering by a bit.
        Defaults to ``True``.
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
        (e.g. `'upper right'`). See `~proplot.axes.BaseAxes.colorbar` for
        details.
    colorbar_kw : dict-like, optional
        Ignored if `colorbar` is ``None``. Extra keyword args for our call
        to `~proplot.axes.BaseAxes` `~proplot.axes.BaseAxes.colorbar` or
        `~proplot.axes.PanelAxes` `~proplot.axes.PanelAxes.colorbar`.

    Other parameters
    ----------------
    lw, linewidth, linewidths
        Aliases. Refers to `linewidths` for `~matplotlib.axes.Axes.contour`,
        `linewidth` for `~matplotlib.axes.Axes.pcolor` and
        `~matplotlib.axes.Axes.pcolormesh`, and `linewidth` for `BaseAxes.cmapline`.
    ls, linestyle, linestyles
        Aliases. Refers to `linestyles` for `~matplotlib.axes.Axes.contour`,
        `linestyle` for `~matplotlib.axes.Axes.pcolor` and
        `~matplotlib.axes.Axes.pcolormesh`, and `linestyle` for `BaseAxes.cmapline`.
    color, colors, edgecolor, edgecolors
        Aliases. Refers to `colors` for `~matplotlib.axes.Axes.contour`,
        `edgecolors` for `~matplotlib.axes.Axes.pcolor` and
        `~matplotlib.axes.Axes.pcolormesh`, and `color` for `BaseAxes.cmapline`.
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

    This could also be done by limiting the number of colors in the
    colormap lookup table by selecting a smaller ``N`` (see
    `~matplotlib.colors.LinearSegmentedColormap`) -- but I prefer the approach
    of always building colormaps with hi-res lookup tables, and leaving the
    job of normalizing data values to colormap locations to the
    `~matplotlib.colors.Normalize` object.

    See also
    --------
    `BaseAxes`, `~proplot.colortools.Colormap`,
    `~proplot.colortools.Norm`, `~proplot.colortools.BinNorm`,
    `~matplotlib.colors.Colormap`, `~matplotlib.colors.Normalize`
    """
    # Input levels
    # See this post: https://stackoverflow.com/a/48614231/4970632
    name = func.__name__
    values_as_keyword = (name=='cmapline') # functions for which 'values' is a native keyword
    if np.iterable(values):
        if values_as_keyword:
            kwargs['values'] = values
            levels = utils.edges(values) # special case, used by colorbar factory
        else:
            if norm is None or isinstance(norm, str) and 'segment' in norm:
                levels = utils.edges(values) # special case, used by colorbar factory
            else:
                norm_tmp = colortools.Norm(norm, **norm_kw)
                levels = norm_tmp.inverse(utils.edges(norm_tmp(values)))
    # Input colormap and other args
    # TODO: Add similar support for quiver, and automatic quiver keys
    # TODO: Make sure matshow, imshow, work! And just like aspect ratios can
    # change for cartopy plots, allow aspect ratio to change for these
    cyclic = False
    colors = _default(color, colors, edgecolor, edgecolors)
    levels = _default(levels, 11) # e.g. pcolormesh can auto-determine levels if you input a number
    if vmin is not None and vmax is not None and not np.iterable(levels):
        levels = np.linspace(vmin, vmax, levels)
    linewidths = _default(lw, linewidth, linewidths)
    linestyles = _default(ls, linestyle, linestyles)
    if not re.match('contour$', name): # contour, tricontour, i.e. not a method where cmap is optional
        cmap = cmap or rc['image.cmap']
    if cmap is not None:
        if isinstance(cmap, (str, dict, mcolors.Colormap)):
            cmap = cmap, # make a tuple
        cmap = colortools.Colormap(*cmap, N=None, **cmap_kw)
        cyclic = cmap._cyclic
        if cyclic and extend!='neither':
            warnings.warn(f'Cyclic colormap selected. Overriding user input extend "{extend}".')
            extend = 'neither'
        kwargs['cmap'] = cmap
    if 'contour' in name: # contour, contourf, tricontour, tricontourf
        kwargs.update({'levels': levels, 'extend': extend})
    elif name in ('imshow', 'matshow', 'spy', 'hist2d'): # aspect ratio fix
        kwargs['aspect'] = 'auto'
    # Disable fix=True for certain keyword combinations, e.g. if user wants
    # white lines around their pcolor mesh.
    # TODO: Allow calling contourf with linewidth?
    for regex,names in _options.items():
        if not re.search(regex, name):
            continue
        # Different cmap functions may or may not accept 'colors', 'linewidths',
        # or 'linestyles' as arguments
        for key,value in (('colors',colors), ('linewidths',linewidths), ('linestyles',linestyles)):
            if not value:
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
                x, y = bbox.intervalx.mean(), bbox.intervaly.mean()
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
            label = _auto_label(args[-1]) # last one is data, we assume
            if label:
                colorbar_kw['label'] = label
        if values_as_keyword and values is not None:
            colorbar_kw['values'] = values
        if loc!='fill':
            colorbar_kw['loc'] = loc
        ax.colorbar(obj, **colorbar_kw)
    return obj

#------------------------------------------------------------------------------#
# Legends and colorbars
#------------------------------------------------------------------------------#
def legend_wrapper(self, handles=None, align=None, order='C',
        loc=None, color=False, linewidth=False, linestyle=False,
        ncol=None, ncols=None, frameon=None, frame=None,
        **kwargs):
    """
    Function for drawing a legend, with some handy added features.

    Parameters
    ----------
    self : `~matplotlib.axes.Axes`
        The axes.
    handles : None or list of `~matplotlib.artist.Artist`, optional
        List of artists instances -- for example, `~matplotlib.lines.Line2D`.
    loc : None or str, optional
        The legend location.
    order : {'C', 'F'}, optional
        Whether legend handles are drawn in column-major (``'C'``) or row-major
        (``'F'``) order. Analagous to `numpy.array` ordering. For some reason
        ``'F'`` was the original matplotlib default; the default is now ``'C'``.
    align : None or bool, optional
        Whether to align rows of legend handles. If ``False``, we actually
        draw successive single-row legends stacked on top of each other,
        and you cannot have a "legend box".

        If ``None``, we try to infer this setting from the `handles` keyword
        arg. Becomes ``True`` if `handles` is a list of lists (implies
        each sublist is a *row* in the legend), ``False`` if not.
    ncol, ncols : int, optional
        The number of columns. `ncols` is an alias, added
        for consistency with `~matplotlib.pyplot.subplots`.

    Other parameters
    ----------------
    **kwargs
        Passed to `~matplotlib.axes.Axes.legend`.

    See also
    --------
    `BaseAxes.colorbar`, `PanelAxes.colorbar`, `~matplotlib.axes.Axes.legend`
    """
    # Parse input
    if order not in ('F','C'):
        raise ValueError(f'Invalid order "{order}". Choose from "C" (row-major, default) and "F" (column-major).')
    # First get legend settings and interpret kwargs.
    ncol = _default(ncols, ncol) # may still be None, wait till later
    frameon = _default(frame, frameon)
    if frameon is not None:
        kwargs.update({'frameon':frameon})
    kwargs.update({'prop': {'family': rc['fontname']}}) # 'prop' can be a FontProperties object or a dict for the kwargs
    # Legend text and handle properties
    hsettings = {}
    for candidate in ['linewidth', 'color']: # candidates for modifying legend objects
        if candidate in kwargs:
            hsettings[candidate] = kwargs.pop(candidate)
    # Legend entry for colormap object
    for i,handle in enumerate(handles):
        if hasattr(handle, 'get_facecolor') or not hasattr(handle, 'get_cmap'): # latter is for scatter (TODO: add cmap_wrapper for scatter?)
            continue
        warnings.warn('Getting legend entry from colormap.')
        size = np.mean(handle.get_sizes())
        cmap = handle.get_cmap()
        handles[i] = self.scatter([0], [0], markersize=size,
                                color=[cmap(0.5)],
                                label=handle.get_label())

    # Detect if user wants to specify rows manually
    # Gives huge latitude for user input:
    #   1) user can specify nothing and align will be inferred (list of iterables
    #      will always be False, i.e. we draw consecutive legends, and list of handles is always true)
    #   2) user can specify align (needs list of handles for True, list of handles or list
    #      of iterables for False and if the former, will turn into list of iterables)
    if not handles and not self._filled:
        handles = self.get_legend_handles_labels()[0]
        if not handles:
            raise ValueError('No axes artists with labels were found.')
    elif not handles:
        raise ValueError('You must pass a handles list for panel axes "filled" with a legend.')
    # Mange input
    list_of_lists = not isinstance(handles[0], martist.Artist)
    if align is None: # automatically guess
        align = not list_of_lists
    elif align and list_of_lists: # standardize format based on input
        handles = [handle for sublist in handles for handle in sublist]
        list_of_lists = False # no longer is list of lists
    elif not align and not list_of_lists:
        list_of_lists = True
        ncol = _default(ncol, 3)
        handles = [handles[i*ncol:(i+1)*ncol] for i in range(len(handles))] # to list of iterables
    elif not align and list_of_lists and ncol is not None:
        warnings.warn('Detected list of *lists* of legend handles. Ignoring user input property "ncol".')
    # Remove empty lists; pops up in some examples, not sure how
    handles = [sublist for sublist in handles if sublist]

    # Now draw legend, with two options
    # 1) Normal legend, just draw everything like normal and columns
    # will be aligned; we re-order handles to be row-major, is only difference
    if align:
        # Optionally change order
        # See: https://stackoverflow.com/q/10101141/4970632
        ncol = _default(ncol, 3)
        if order=='C':
            newhandles = []
            handlesplit = [handles[i*ncol:(i+1)*ncol] for i in range(len(handles)//ncol+1)] # split into rows
            nrowsmax, nfinalrow = len(handlesplit), len(handlesplit[-1]) # max possible row count, and columns in final row
            # e.g. if 5 columns, but final row length 3, columns 0-2 have N rows but 3-4 have N-1 rows
            nrows = [nrowsmax]*nfinalrow + [nrowsmax-1]*(ncol-nfinalrow)
            for col,nrow in enumerate(nrows): # iterate through cols
                newhandles.extend(handlesplit[row][col] for row in range(nrow))
            handles = newhandles
        # Finally draw legend, mimicking row-major ordering
        if loc is not None:
            kwargs['loc'] = _loc_translate.get(loc, loc)
        leg = maxes.Axes.legend(self, ncol=ncol, handles=handles, **kwargs)
        legends = [leg]

    # 2) Separate legend for each row
    # The label spacing/border spacing will be exactly replicated, as if we were
    # using the original legend command
    # Means we also have to overhaul some settings
    else:
        # Warn when user input props are overridden
        overridden = []
        loc = _default(loc, 'upper center')
        loc = _loc_translate.get(loc, loc)
        if loc=='best':
            warnings.warn('Cannot use "best" location for un-aligned legend. Defaulting to "upper center".')
            overridden.append('loc')
            loc = 'upper center'
        for override in ['bbox_transform', 'bbox_to_anchor', 'frameon']:
            prop = kwargs.pop(override, None)
            if prop is not None:
                overridden.append(override)
        if overridden:
            warnings.warn(f'Creating unaligned legend. Overriding user input legend properties "' + '", "'.join(prop for prop in overridden) + '".')
        # Determine space we want sub-legend to occupy, as fraction of height
        # Don't normally save "height" and "width" of axes so keep here
        fontsize = kwargs.get('fontsize', None)     or rc['legend.fontsize']
        spacing  = kwargs.get('labelspacing', None) or rc['legend.labelspacing']
        interval = 1/len(handles) # split up axes
        interval = (((1 + spacing)*fontsize)/72) / \
                (self.figure.get_figheight() * np.diff(self._position.intervaly))
        # Iterate and draw
        # NOTE: We confine possible bounding box (within which legend position
        # is allowed to vary) in *y*-direction, but do not confine it in
        # *x*-direction (notice the *x*-coordinate range 0-1). Matplotlib will
        # automatically move left-to-right if you request this.
        legends = []
        if order=='F':
            raise NotImplementedError(f'When align=False, ProPlot vertically stacks successive single-row legends. Column-major (order="F") ordering is un-supported.')
        for h,hs in enumerate(handles):
            if 'upper' in loc:
                y1 = 1 - (h+1)*interval
                y2 = 1 - h*interval
            elif 'lower' in loc:
                y1 = (len(handles) + h - 2)*interval
                y2 = (len(handles) + h - 1)*interval
            else: # center
                y1 = 0.5 + interval*len(handles)/2 - (h+1)*interval
                y2 = 0.5 + interval*len(handles)/2 - h*interval
            bbox = mtransforms.Bbox([[0, y1], [1, y2]])
            leg = maxes.Axes.legend(self,
                handles=hs, ncol=len(hs), loc=loc,
                frameon=False,
                bbox_transform=self.transAxes,
                bbox_to_anchor=bbox,
                **kwargs) # _format_legend is overriding original legend Method
            legends.append(leg)
        for l in legends[:-1]:
            self.add_artist(l) # because matplotlib deletes previous ones

    # Properties for legends
    outline = rc.fill({
        'linewidth':'axes.linewidth',
        'edgecolor':'axes.edgecolor',
        'facecolor':'axes.facecolor',
        }, cache=False)
    for leg in legends:
        # for t in leg.texts:
        leg.legendPatch.update(outline) # or get_frame()
        for obj in leg.legendHandles:
            obj.update(hsettings)
    return legends[0] if len(legends)==1 else legends

def colorbar_wrapper(self, mappable, values=None,
        extend=None, extendsize=None, label=None,
        grid=None, tickminor=None,
        tickloc=None, ticklocation=None,
        locator=None, ticks=None, minorlocator=None, minorticks=None, locator_kw={}, minorlocator_kw={},
        formatter=None, ticklabels=None, formatter_kw={},
        fixticks=False, norm=None, norm_kw={}, # normalizer to use when passing colors/lines
        orientation='horizontal',
        **kwargs):
    """
    Function for filling an axes with a colorbar, with some handy added
    features.

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
           a colormap will be generated and used with the colorbar. Requires
           `values` is not ``None``.

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
        If float, units are inches. If string,
        units are interpreted by `~proplot.utils.units`.

        This is handy if you have multiple colorbars in one figure.
        With the matplotlib API, it is really hard to get triangle
        sizes to match, because the `extendsize` units are *relative*.
    tickloc
        Alias for `ticklocation`.
    ticklocation : {'bottom', 'top', 'left', 'right'}, optional
        Where to draw tick marks on the colorbar.
    label : None or str, optional
        The colorbar label.
    grid : bool, optional
        Whether to draw "gridlines" (i.e. separators) between each level
        across the colorbar. Default is ``False``.
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
    `BaseAxes.colorbar`, `PanelAxes.colorbar`, `~matplotlib.figure.Figure.colorbar`,
    `~proplot.axistools.Locator`, `~proplot.axistools.Formatter`, `~proplot.colortools.Norm`

    Warning
    -------
    Colorbar axes must be of type `matplotlib.axes.Axes`,
    not `~proplot.axes.BaseAxes`, because colorbar uses some internal methods
    that are wrapped by `~proplot.axes.BaseAxes`.
    """
    # Developer notes
    # * There are options on the colorbar object (cb.locator,
    #   cb.formatter with cb.update_ticks) and by passing kwargs (ticks=x,
    #   format=y) that allow user to not reference the underlying "axes"
    #   when fixing ticks. Don't use this functionality because not necessary
    #   for us and is missing many features, e.g. minorlocators/minorformatters.
    #   Also is different syntax.
    # * There is an insanely weird problem with colorbars when simultaneously
    #   passing levels and norm object to a mappable; fixed by passing
    #   vmin/vmax instead of levels.
    #   (see: https://stackoverflow.com/q/40116968/4970632).
    # * Problem is, often want levels instead of vmin/vmax, while simultaneously
    #   using a Normalize (for example) to determine colors between the levels
    #   (see: https://stackoverflow.com/q/42723538/4970632). Workaround is to
    #   make sure locators are in vmin/vmax range exclusively; cannot match/exceed values.
    # * The 'extend' kwarg is used for the case when you are manufacturing
    #   colorbar from list of colors or lines. Most of the time want 'neither'.
    # Parse flexible input
    ticklocation = _default(tickloc, ticklocation)
    locator = _default(ticks, locator)
    minorlocator = _default(minorticks, minorlocator)
    formatter = _default(ticklabels, formatter, 'default')

    # Test if we were given a mappable, or iterable of stuff; note Container and
    # PolyCollection matplotlib classes are iterable.
    fromlines, fromcolors = False, False
    if np.iterable(mappable) and len(mappable)==2:
        mappable, values = mappable
    if not isinstance(mappable, martist.Artist) and not isinstance(mappable, mcontour.ContourSet):
        if isinstance(mappable[0], martist.Artist):
            fromlines = True # we passed a bunch of line handles; just use their colors
        else:
            fromcolors = True # we passed a bunch of color strings or tuples
    # Update with user-kwargs
    if extend is None:
        if hasattr(mappable, 'extend'):
            extend = mappable.extend or 'neither'
        else:
            extend = 'neither'
    kwdefault = {'cax':self, 'orientation':orientation, 'use_gridspec':True, # use space afforded by entire axes
                 'spacing':'uniform', 'extend':extend, 'drawedges':grid} # this is default case unless mappable has special props
    kwdefault.update(kwargs)
    kwargs = kwdefault

    # Option to generate colorbar/colormap from line handles
    # * Note the colors are perfect if we don't extend them by dummy color on either side,
    #   but for some reason labels for edge colors appear offset from everything
    #   Too tired to figure out why so just use this workaround
    # * Note contourf will not be overridden for colorbar axes! Need to
    #   manually wrap with cmap_wrapper.
    if fromlines or fromcolors:
        # Colors
        if fromcolors:
            colors = mappable
        else:
            colors = [h.get_color() for h in mappable]
            if values is None:
                values = []
                for obj in mappable:
                    val = obj.get_label()
                    values.append(val)
        # Verify values are numeric
        ivalues = []
        for val in values:
            try:
                val = float(val)
            except ValueError:
                raise ValueError(f'To generate colorbars from line handles or other objects, pass the "values" keyword arg to colorbar(), or give your handles numeric values with e.g. plot(..., label=123) or line.set_label(123).')
            ivalues.append(val)
        values = np.array(ivalues)
        # Make plot
        cmap = colortools.Colormap(colors)
        func = _cmap_wrapper(self, self.contourf)
        mappable = func([[0,0],[0,0]],
            values=values, cmap=cmap, extend='neither',
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

    # Fix the norm object
    # Check out the *insanely weird error* that occurs when you comment out this block!
    # * The error is triggered when a *major* tick sits exactly on vmin, but
    #   the actual error is due to processing of *minor* ticks, even if the 
    #   minor locator was set to NullLocator; very weird
    # * Happens when we call get_ticklabels(which='both') below. Can be prevented
    #   by just calling which='major'. Minor ticklabels are never drawn anyway.
    # * We can eliminate the normfix below, but that actually causes an annoying
    #   warning to be printed (related to same issue I guess). So we keep this.
    #   The culprit for all of this seems to be the colorbar API line:
    #        z = np.take(y, i0) + (xn - np.take(b, i0)) * dy / db
    # * Also strange that minorticks extending *below* the minimum
    #   don't raise the error. It is only when they are exaclty on the minimum.
    # * Note that when changing the levels attribute, need to make sure the
    #   levels datatype is float; otherwise division will be truncated and bottom
    #   level will still lie on same location, so error will occur
    if normfix:
        mappable.norm.vmin -= (mappable.norm.vmax-mappable.norm.vmin)/10000
    if hasattr(mappable.norm, 'levels'):
        mappable.norm.levels = np.atleast_1d(mappable.norm.levels).astype(np.float)
        if normfix:
            mappable.norm.levels[0] -= np.diff(mappable.norm.levels[:2])[0]/10000

    # Draw the colorbar
    # NOTE: Only way to avoid bugs seems to be to pass the major formatter/locator
    # to colorbar commmand and directly edit the minor locators/formatters;
    # update_ticks after the fact ignores the major formatter.
    # axis.set_major_locator(locators[0]) # does absolutely nothing
    # axis.set_major_formatter(formatter)
    width, height = self.figure.get_size_inches()
    if orientation=='horizontal':
        axis = self.xaxis
        scale = width*np.diff(getattr(self.get_position(),'intervalx'))[0]
    else:
        axis = self.yaxis
        scale = height*np.diff(getattr(self.get_position(),'intervaly'))[0]
    extendsize = utils.units(_default(extendsize, rc.get('colorbar.extend')))
    extendsize = extendsize/(scale - 2*extendsize)
    formatter    = axistools.Formatter(formatter, **formatter_kw)
    kwargs.update({'ticks':locators[0], # WARNING: without this, set_ticks screws up number labels for some reason
                   'format':formatter,
                   'ticklocation':ticklocation,
                   'extendfrac':extendsize})
    cb = self.figure.colorbar(mappable, **kwargs)
    # Make edges/dividers consistent with axis edges
    if cb.dividers is not None:
        cb.dividers.update(rc['grid'])

    # The minor locators and formatters
    # * The minor locator must be set with set_ticks after transforming an array
    #   using the mappable norm object; see: https://stackoverflow.com/a/20079644/4970632
    # * The set_minor_locator seems to be completely ignored depending on the colorbar
    #   in question, for whatever reason, and cb.minorticks_on() gives no control.
    # NOTE: Re-apply major ticks here because for some reason minor ticks don't
    # align with major ones for LogNorm. When we call set_ticks, labels (and
    # numbers) are not changed; just re-adjust existing ticks to proper locations.
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
    axis.set_ticks(minorvals, minor=True)
    if fixticks:
        axis.set_ticks(majorvals, minor=False)
    axis.set_minor_formatter(mticker.NullFormatter()) # to make sure
    # The label
    if label is not None:
        axis.label.update({'text':label})

    # Fix alpha issues (cannot set edgecolor to 'face' if alpha non-zero
    # because blending will occur, will get colored lines instead of white ones;
    # need to perform manual alpha blending)
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
    # to rasterized blocks
    if cb.solids:
        cb.solids.set_linewidth(0.2) # something small
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
# Also _m_call and _m_norecurse
_parse_1d_ = _wrapper(_parse_1d)
_parse_2d_ = _wrapper(_parse_2d)
# Documented
_check_centers = _wrapper(check_centers)
_check_edges   = _wrapper(check_edges)
_basemap_gridfix   = _wrapper(basemap_gridfix)
_basemap_latlon    = _wrapper(basemap_latlon)
_cartopy_gridfix   = _wrapper(cartopy_gridfix)
_cartopy_transform = _wrapper(cartopy_transform)
_cartopy_crs       = _wrapper(cartopy_crs)
_cmap_wrapper    = _wrapper(cmap_wrapper)
_cycle_wrapper   = _wrapper(cycle_wrapper)
_plot_wrapper    = _wrapper(plot_wrapper)
_scatter_wrapper = _wrapper(scatter_wrapper)
_text_wrapper    = _wrapper(text_wrapper)
_legend_wrapper  = _simple_wrapper(legend_wrapper) # special

