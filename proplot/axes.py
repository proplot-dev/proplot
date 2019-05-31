#!/usr/bin/env python3
"""
The axes subclasses central to this library, plus the enhanced
colorbar and legend functions. You should focus on the ``format`` and
``smart_update`` methods:

* `BaseAxes.format`
* `BaseAxes.smart_update`
* `XYAxes.smart_update`
* `MapAxes.smart_update`

`BaseAxes.format` calls the ``smart_update`` methods in turn.
This is your one-stop-shop for changing axes settings like
*x* and *y* axis limits, axis labels, tick locations, tick label
format, grid lines, scales, titles, a-b-c labelling, adding
geographic features, and much more.


ProPlot also adds various new settings by wrapping native
`~matplotlib.axes.Axes` plotting methods with the ``wrapper`` functions.
See the notes below for details.

.. raw:: html

   <h1>Developer notes</h1>

The wrapping is done by dynamically overriding the
`~proplot.axes.BaseAxes.__getattribute__` methods on
`~proplot.axes.BaseAxes` and its subclasses. You may be wondering: why did I
choose to do this, rather than using **decorators**? Two reasons:

1. Brevity. For example: `wrapper_cmap` overrides **a dozen** different
   methods. This lets me override these methods in *one* line, instead of 50
   lines. To see which methods are overriden, the user can simply check the
   documentation.
2. Documentation. If I wrapped every method, the sphinx `autodoc <http://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html>`_
   documentation generator would inherit docstrings from the parent methods.
   In other words, the plotting method docstrings would get duplicated on
   my website from the matplotlib website, generally with a bunch of errors.

   I could also override these methods with my own docstrings, but that would
   mean when the user tries e.g. ``help(ax.contourf)``, they would see my
   own, very brief docstring instead of the comprehensive matplotlib reference
   they were probably looking for! I only add a handful of features to these
   functions, but the native methods generally have way more options.

It should be noted that dynamically wrapping every time the user accesses
the corresponding method will be slower than "decoration", which just wraps them
once on declaration of the axes. But this was not found to significantly affect
performance. And anyway, `Premature Optimization is the Root of All Evil
<http://wiki.c2.com/?PrematureOptimization>`__.
"""
# Note that even if not in IPython notebook, io capture output still works
# functools.wraps preserves __name__ metadata; see comment:
# https://stackoverflow.com/a/739665/4970632
import os
import re
import numpy as np
import warnings
import functools
from numbers import Number
from IPython.utils import io
from matplotlib.cbook import mplDeprecation
# from matplotlib.lines import _get_dash_pattern, _scale_dashes
import matplotlib.projections as mproj
import matplotlib.figure as mfigure
import matplotlib.axes as maxes
import matplotlib.scale as mscale
import matplotlib.contour as mcontour
import matplotlib.patheffects as mpatheffects
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import matplotlib.text as mtext
import matplotlib.ticker as mticker
import matplotlib.artist as martist
import matplotlib.gridspec as mgridspec
import matplotlib.transforms as mtransforms
import matplotlib.collections as mcollections

# Local modules, projection sand formatters and stuff
# TODO: Import matplotlib docstring func, use it?
from .rcmod import rc, _rc_names_nodots
from . import utils, projs, colortools, fonttools, axistools
from .utils import _default, _timer, _counter, ic, units
from .gridspec import FlexibleGridSpec, FlexibleGridSpecFromSubplotSpec

# Silly recursive function, returns a...z...aa...zz...aaa...zzz
_abc_string = 'abcdefghijklmnopqrstuvwxyz'
def _abc(i):
    if i < 26:
        return _abc_string[i]
    else:
        return _abc(i - 26) + _abc_string[i % 26]

# Filter warnings, seems to be necessary before drawing stuff for first time,
# otherwise this has no effect (e.g. if you stick it in a function)
warnings.filterwarnings('ignore', category=mplDeprecation)
# Optionally import mapping toolboxes
# Main conda distro says they are incompatible, so make sure not required!
try:
    from cartopy.mpl.geoaxes import GeoAxes
    from cartopy.crs import PlateCarree
except ModuleNotFoundError:
    GeoAxes = PlateCarree = object

#------------------------------------------------------------------------------#
# Decorators
# To bulk decorate lists of methods, instead of wrapping methods explicitly,
# we do some hacky bullshit with __getattribute__ and apply 'wrapper' there
# Do this because it's cleaner, not really any performance issues since even
# if we look up attributes thousands of times, testing membership in a length-4
# list is nanoseconds level
#------------------------------------------------------------------------------#
# 2D plot functions that require coordinate centers and edges
_centers_methods = (
    'contour', 'contourf', 'quiver', 'streamplot', 'barbs'
    )
_edges_methods = (
    'pcolor', 'pcolormesh',
    )
_latlon_methods = ( # adds latlon=True, methods not in edges or centers methods
    'plot', 'scatter',
    )
_transform_methods = ( # adds transform=PlateCarree, methods not in edges or centers methods
    'plot', 'scatter', 'tripcolor', 'tricontour', 'tricontourf',
    )
_crs_methods = ( # special extras
    'get_extent', 'set_extent', 'set_xticks', 'set_yticks',
    )
# Whether to wrap plot functions with cycle features or cmap features
_cycle_methods  = (
    'plot', 'scatter', 'bar', 'barh', 'hist', 'boxplot', 'errorbar'
    )
_cmap_methods = (
    # 'quiver', 'streamplot',
    'contour', 'contourf', 'pcolor', 'pcolormesh',
    'cmapline', 'hexbin', # special
    'matshow', 'imshow', 'spy', 'hist2d',
    'tripcolor', 'tricontour', 'tricontourf',
    )
# Disable some stuff for all axes, and just for map projection axes
# The keys in below dictionary are error messages
_disabled_methods = {
    # "Unsupported plotting function {}. May be added soon.":
    #     ('table', 'eventplot', # pie?
    #     'xcorr', 'acorr', 'psd', 'csd', 'magnitude_spectrum',
    #     'angle_spectrum', 'phase_spectrum', 'cohere', 'specgram'),
    "Redundant function {} has been disabled. Control axis scale with format(xscale='scale', yscale='scale').":
        ('semilogx', 'semilogy', 'loglog'),
    "Redundant function {} has been disabled. Date formatters will be used automatically when x/y coordinates are python datetime or numpy datetime64.":
        ('plot_date',),
    "Redundant function {} has been disabled. Use proj='polar' in subplots() call, then use the angle as your *x* coordinate and radius as your *y* coordinate.":
        ('polar',)
    }
_map_disabled_methods = (
    'matshow', 'imshow', 'spy', # don't disable 'bar' or 'barh', can be used in polar plots
    'hist', 'hist2d', 'errorbar', 'boxplot', 'violinplot', 'step', 'stem',
    'hlines', 'vlines', 'axhline', 'axvline', 'axhspan', 'axvspan',
    # 'triplot', 'tricontour', 'tricontourf', 'tripcolor',
    # 'fill_between', 'fill_betweenx', 'fill', # can be used to fill between lons/lats
    'stackplot',
    'table', 'eventplot', 'pie',
    'xcorr', 'acorr', 'psd', 'csd', 'magnitude_spectrum',
    'angle_spectrum', 'phase_spectrum', 'cohere', 'specgram',
    )
# Aliases for panel names
_aliases = {
    'bpanel': 'bottompanel',
    'rpanel': 'rightpanel',
    'tpanel': 'toppanel',
    'lpanel': 'leftpanel'
    }
# Keywords
# These may automatically override the 'fix' option!
# NOTE: Pcolor edgecolors is really just more *flexible* version of 'color';
# and if the latter is specified will overwrite the former by default, so we
# just do the same thing.
_options = {
    'contour$': {'colors':'colors', 'linewidths':'linewidths', 'linestyles':'linestyles'},
    'cmapline': {'colors':'color',  'linewidths':'linewidth', 'linestyles':'linestyle'},
    '^pcolor':  {'colors':'edgecolors', 'linewidths':'linewidth', 'linestyles':'linestyle'},
    }

#------------------------------------------------------------------------------#
# For documentation
#------------------------------------------------------------------------------#
def _sphinx_name(name):
    """Gets sphinx name."""
    if name=='cmapline':
        return f'`~BaseAxes.{name}`'
    else:
        return f'`~matplotlib.axes.Axes.{name}`'
def _expand_wrapped_methods(func):
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
        doc = re.sub(f'`{name}`', ', '.join(_sphinx_name(method) for method in methods[:-1])
                + ','*min((len(methods)-2, 1)) + f' and {_sphinx_name(methods[-1])}', doc)
    func.__doc__ = doc
    return func

#------------------------------------------------------------------------------
# Wrappers for standardizing 2D grid inputs for contour, etc.
# NOTE: All wrappers are actully wrapper *pairs*. The hidden just call the
# public ones of the same name.
# used, the public ones do the work. 
#------------------------------------------------------------------------------
def _parse_args(args):
    """Helper function for `wrapper_check_edges` and `wrapper_check_centers`."""
    # Sanitize input
    if len(args)>2:
        Zs = args[2:]
    else:
        Zs = args
    Zs = [np.array(Z) for Z in Zs] # ensure array
    if any(Z.ndim!=2 for Z in Zs):
        raise ValueError(f'Zs must be 2-dimensional, got {[Z.ndim for Z in Zs]}.')
    if len(args)>2:
        x, y = args[:2]
        x, y = np.array(x), np.array(y)
    else:
        x = np.arange(Zs[0].shape[0])
        y = np.arange(Zs[0].shape[1])
    # Custom error messages, more useful than default ones that would arise
    if any(Z.ndim != 2 for Z in Zs):
        raise ValueError(f'Data array dimensionalities are {[Z.ndim for Z in Zs]}, but should be 2D.')
    for string,array in zip(('X','Y'), (x,y)):
        if array.ndim not in (1,2):
            raise ValueError(f'{string} coordinates are {array.ndim}D, but should be 1D or 2D.')
    if x.ndim != y.ndim:
        raise ValueError(f'X coordinates are {x.ndim}D, but Y coordinates are {y.ndim}D.')
    return x, y, Zs

@_expand_wrapped_methods
def wrapper_check_centers(func, *args, order='C', **kwargs):
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
    x, y, Zs = _parse_args(args)
    xlen, ylen = x.shape[-1], y.shape[0]
    for Z in Zs:
        if Z.ndim!=2:
            raise ValueError(f'Input arrays must be 2D, instead got shape {Z.shape}.')
        elif Z.shape[1]==xlen-1 and Z.shape[0]==ylen-1 and x.ndim==1 and y.ndim==1:
            x = (x[1:] + x[:-1])/2
            y = (y[1:] + y[:-1])/2 # get centers, given edges
        elif Z.shape[1]!=xlen or Z.shape[0]!=ylen:
            raise ValueError(f'X ({"x".join(str(i) for i in x.shape)}) '
                    f'and Y ({"x".join(str(i) for i in y.shape)}) must correspond to '
                    f'nrows ({Z.shape[0]}) and ncolumns ({Z.shape[1]}) of Z, or its borders.')
    # Optionally re-order
    if order=='F':
        x, y = x.T, y.T # in case they are 2-dimensional
        Zs = (Z.T for Z in Zs)
    elif order!='C':
        raise ValueError(f'Invalid order "{order}". Choose from "C" (row-major, default) and "F" (column-major).')
    result = func(x, y, *Zs, **kwargs)
    return result

@_expand_wrapped_methods
def wrapper_check_edges(func, *args, order='C', **kwargs):
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
    x, y, Zs = _parse_args(args)
    xlen, ylen = x.shape[-1], y.shape[0]
    for Z in Zs:
        if Z.ndim!=2:
            raise ValueError(f'Input arrays must be 2D, instead got shape {Z.shape}.')
        elif Z.shape[1]==xlen and Z.shape[0]==ylen:
            # If 2D, don't raise error, but don't fix either, because
            # matplotlib pcolor accepts grid center inputs.
            # TODO: Fix 2D inputs?
            if x.ndim==1 and y.ndim==1:
                x, y = utils.edges(x), utils.edges(y)
        elif Z.shape[1]!=xlen-1 or Z.shape[0]!=ylen-1:
            raise ValueError(f'X ({"x".join(str(i) for i in x.shape)}) '
                    f'and Y ({"x".join(str(i) for i in y.shape)}) must correspond to '
                    f'nrows ({Z.shape[0]}) and ncolumns ({Z.shape[1]}) of Z, or its borders.')
    # Optionally re-order
    if order=='F':
        x, y = x.T, y.T # in case they are 2-dimensional
        Zs = (Z.T for Z in Zs)
    elif order!='C':
        raise ValueError(f'Invalid order "{order}". Choose from "C" (row-major, default) and "F" (column-major).')
    result = func(x, y, *Zs, **kwargs)
    return result

#------------------------------------------------------------------------------#
# Simple geographic wrappers
#------------------------------------------------------------------------------#
# First basemap recursion fix
# Normally we *cannot* modify the underlying *axes* pcolormesh etc. because this
# this will cause basemap's self.m.pcolormesh etc. to use *custom* version and
# cause suite of weird errors. Prevent this recursion with the below decorator.
def _wrapper_m_call(self, func):
    """Docorator that calls the basemap version of the function of the same name."""
    name = func.__name__
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return self.m.__getattribute__(name)(ax=self, *args, **kwargs)
    return wrapper

def _wrapper_m_norecurse(self, func):
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
            # basemap.Basemap method thanks to the _wrapper_m_call wrapper.
            self._recurred = True
            result = func(*args, **kwargs)
        self._recurred = False # cleanup, in case recursion never occurred
        return result
    return wrapper

@_expand_wrapped_methods
def wrapper_cartopy_transform(self, func, *args, transform=PlateCarree, **kwargs):
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
    # Some plot functions seem to reset the outlinepatch or
    # backgroundpatch (???), so need to re-enforce settings.
    # TODO: Double check this
    self.format()
    return result

@_expand_wrapped_methods
def wrapper_cartopy_crs(self, func, *args, crs=PlateCarree, **kwargs):
    """
    Wraps axes functions for `CartopyAxes` (`_crs_methods`).

    As in `wrapper_cartopy_transform`, but sets ``crs=cartopy.crs.PlateCarree()``
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

@_expand_wrapped_methods
def wrapper_basemap_latlon(self, func, *args, latlon=True, **kwargs):
    """
    Wraps plotting functions for `BasemapAxes` (`_latlon_methods`).

    With the default `~mpl_toolkits.basemap` API, you need to pass
    ``latlon=True`` if your data coordinates are longitude and latitude,
    instead of map projection coordinates. Now, ``latlon=True`` is always
    used.
    """
    return func(*args, latlon=latlon, **kwargs)

#------------------------------------------------------------------------------#
# Geographic grid fixes
#------------------------------------------------------------------------------#
@_expand_wrapped_methods
def wrapper_cartopy_gridfix(self, func, lon, lat, Z, transform=PlateCarree, globe=False, **kwargs):
    """
    Wraps 2D plotting functions for `CartopyAxes` (`_centers_edges_methods`).

    As in `wrapper_cartopy_transform`, but adds the `globe` keyword arg
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
    # Below only works for vector data
    if lon.ndim==1 and lat.ndim==1:
        # 1) Fix holes over poles by *interpolating* there (equivalent to
        # simple mean of highest/lowest latitude points)
        if globe:
            Z_south = np.repeat(Z[0,:].mean(),  Z.shape[1])[None,:]
            Z_north = np.repeat(Z[-1,:].mean(), Z.shape[1])[None,:]
            lat = np.concatenate(([-90], lat, [90]))
            Z = np.concatenate((Z_south, Z, Z_north), axis=0)
        # 2) Fix seams at map boundary; by ensuring circular coverage
        if globe and (lon[0] % 360) != ((lon[-1] + 360) % 360):
            lon = np.array((*lon, lon[0] + 360)) # make longitudes circular
            Z = np.concatenate((Z, Z[:,:1]), axis=1) # make data circular
    # Instantiate transform
    if isinstance(transform, type):
        transform = transform() # instantiate
    # Call function
    # with io.capture_output() as captured:
    result = func(lon, lat, Z, transform=transform, **kwargs)
    # Some plot functions seem to reset the outlinepatch or
    # backgroundpatch (???), so need to re-enforce settings.
    # TODO: Double check this
    self.format()
    return result

@_expand_wrapped_methods
def wrapper_basemap_gridfix(self, func, x, y, Z, globe=False, latlon=True, **kwargs):
    """
    Wraps 2D plotting functions for `BasemapAxes` (`_centers_edges_methods`).

    As in `wrapper_basemap_latlon`, but cycles longitudes to fit within the
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
    if not latlon:
        return func(x, y, Z, **kwargs)
    # Raise errors
    eps = 1e-3
    lon, lat = x, y
    lonmin, lonmax = self.m.lonmin, self.m.lonmax
    if lon.max() > lon.min() + 360 + eps:
        raise ValueError(f'Longitudes span {lon.min()} to {lon.max()}. Can only span 360 degrees at most.')
    if lon.min() < -360 or lon.max() > 360:
        raise ValueError(f'Longitudes span {lon.min()} to {lon.max()}. Must fall in range [-360, 360].')
    # Establish 360-degree range
    lon -= 720
    while True:
        filter_ = lon<lonmin
        if filter_.sum()==0:
            break
        lon[filter_] += 360
    # Below only works with vectors
    if lon.ndim==1 and lat.ndim==1:
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
    # Finally get grid of x/y map projection coordinates
    lat[lat>90], lat[lat<-90] = 90, -90 # otherwise, weird stuff happens
    x, y = self.m(*np.meshgrid(lon, lat))
    # Prevent error where old boundary, drawn on a different axes, remains
    # to the Basemap instance, which means it is not in self.patches, which
    # means Basemap tries to draw it again so it can clip the contours by the
    # resulting path, which raises error because you can't draw on Artist on multiple axes
    self.m._mapboundarydrawn = self.boundary # stored the axes-specific boundary here
    # Call function
    return func(x, y, Z, **kwargs)

#------------------------------------------------------------------------------#
# Colormaps and color cycles
#------------------------------------------------------------------------------#
@_expand_wrapped_methods
def wrapper_cmap(self, func, *args, cmap=None, cmap_kw={},
    norm=None, norm_kw={},
    extend='neither',
    values=None, levels=None, zero=False, # override levels to be *centered* on zero
    labels=False, labels_kw={}, precision=2,
    edgefix=True,
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
    labels_kw
        Ignored if `labels` is ``False``. Extra keyword args for the labels.
        For `~matplotlib.axes.Axes.contour`, passed to `~matplotlib.axes.Axes.clabel`.
        For `~matplotlib.axes.Axes.pcolor` or `~matplotlib.axes.Axes.pcolormesh`,
        passed to `~matplotlib.axes.Axes.text`.
    precision : int, optional
        Maximum precision (in decimal places) for the number labels.
        Number labels are generated with the `~proplot.axistools.SimpleFormatter`
        formatter, which allows us to limit the precision.

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
    # Input colormap
    cyclic = False
    colors = _default(color, colors, edgecolor, edgecolors)
    levels = _default(levels, 11) # e.g. pcolormesh can auto-determine levels if you input a number
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
    elif name in ('imshow', 'matshow', 'spy', 'hist2d'): # aspect ratio settings
        kwargs['aspect'] = 'auto'
    # Disable fix=True for certain keyword combinations, e.g. if user wants
    # white lines around their pcolor mesh.
    # TODO: Double check matshow, etc.!
    for regex,names in _options.items():
        if not re.search(regex, name):
            continue
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
        # TODO: This still is not respected for hexbin 'log' norm for some reason, fix.
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
        # TODO: Should use normalizer and colormap to get colors instead of
        # update scalar mappable? Also gets called at drawtime, maybe redundant.
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

    # Fix white lines between filled contours/mesh
    # Allow user to override properties!
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
    return obj

@_expand_wrapped_methods
def wrapper_cycle(self, func, *args, cycle=None, cycle_kw={}, **kwargs):
    """
    Wraps methods that use the property cycler (`_cycle_methods`),
    adds features for controlling colors in the property cycler.

    Parameters
    ----------
    cycle : None or cycle spec, optional
        The cycle specifer, passed to the `~proplot.colortools.Cycle`
        constructor. If the returned list of colors is unchanged from the
        current axes color cycler, the axes cycle will not be reset to the first
        position.
    cycle_kw : dict-like, optional
        Passed to `~proplot.colortools.Cycle`.

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
    # Determine and temporarily set cycler
    # WARNING: Axes cycle has no getter, only setter (set_prop_cycle), which
    # sets a 'prop_cycler' attribute on the hidden _get_lines and
    # _get_patches_for_fill objects. This is the only way to query the "current"
    # axes cycler! Could also have wrapped set_prop_cycle but that calls
    # several secondary methods, may get messy.
    if cycle is not None:
        # Get the new group of colors
        if not np.iterable(cycle) or isinstance(cycle, str):
            cycle = cycle,
        cycle = colortools.Cycle(*cycle, **cycle_kw)
        # Compare to the original group of colors, and only reset the
        # cycle if this group of colors is new
        # NOTE: The _get_lines cycler is an *itertools cycler*. Has no length,
        # so we must cycle over it with next(). We try calling next() the same
        # number of times as the length of user input cycle.
        # NOTE: If the input cycle *is* in fact the same, below does not reset
        # the color position, as it cycles us back to the original position!
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
    return func(*args, **kwargs)

#------------------------------------------------------------------------------#
# Simple wrappers that apply to just one function
#------------------------------------------------------------------------------#
def wrapper_plot(self, func, *args, cmap=None, values=None, **kwargs):
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
    # Make normal boring lines
    if cmap is None:
        lines = func(*args, **kwargs)
    # Make special colormap lines
    else:
        lines = self.cmapline(*args, cmap=cmap, values=values, **kwargs)
    return lines

def wrapper_scatter(self, func, *args,
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
    if len(args)>4:
        raise ValueError(f'Function accepts up to 4 args, received {len(args)}.')
    args = [*args]
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

def _update_text(obj, kwargs):
    """Allow updating new text properties introduced by override."""
    # Allow updating properties introduced by the BaseAxes.text() override.
    # Don't really want to subclass mtext.Text; only have a few features
    try:
        obj.update(kwargs)
        return obj
    except Exception:
        pass
    # Destroy original text instance and get its properties
    obj.set_visible(False)
    text = kwargs.pop('text', obj.get_text())
    if 'color' not in kwargs:
        kwargs['color'] = obj.get_color()
    if 'weight' not in kwargs:
        kwargs['weight'] = obj.get_weight()
    if 'fontsize' not in kwargs:
        kwargs['fontsize'] = obj.get_fontsize()
    # Position
    pos = obj.get_position()
    x, y = None, None
    if 'position' in kwargs:
        x, y = kwargs.pop('position')
    x = kwargs.pop('x', x)
    y = kwargs.pop('y', y)
    if x is None:
        x = pos[0]
    if y is None:
        y = pos[1]
    if np.iterable(x):
        x = x[0]
    if np.iterable(y):
        y = y[0]
    return obj.axes.text(x, y, text, **kwargs)

def wrapper_text(self, func, x, y, text,
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
# Legends and colorbars
#------------------------------------------------------------------------------#
def legend_factory(ax, handles=None, align=None, order='C', **kwargs):
    """
    Function for drawing a legend, with some handy added features.

    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        The axes.
    handles : None or list of `~matplotlib.artist.Artist`, optional
        List of artists instances -- for example, `~matplotlib.lines.Line2D`.
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
    ncol : int, optional
        The number of columns.
    ncols
        Alias for `ncol`. Added for consistency with
        `~matplotlib.pyplot.subplots`.

    Other parameters
    ----------------
    **kwargs
        Passed to `~matplotlib.axes.Axes.legend`.

    See also
    --------
    `BaseAxes.colorbar`, `PanelAxes.colorbar`, `~matplotlib.axes.Axes.legend`
    """
    # First get legend settings (usually just one per plot so don't need to declare
    # this dynamically/globally), and interpret kwargs.
    # TODO: Should update this function to clip the legend box when it goes
    # outside axes area, so the legend width and bottom or right widths can be
    # chosen propertly/separately.
    for name,alias in [('ncol', 'ncols'), ('frameon', 'frame')]:
        if alias in kwargs:
            kwargs[name] = kwargs.pop(alias)
    if order not in ('F','C'):
        raise ValueError(f'Invalid order "{order}". Choose from "C" (row-major, default) and "F" (column-major).')
    # Setup legend text and handle properties
    hsettings = {}
    for candidate in ['linewidth', 'color']: # candidates for modifying legend objects
        if candidate in kwargs:
            hsettings[candidate] = kwargs.pop(candidate)
    # Font name; 'prop' can be a FontProperties object or a dict for the kwargs
    # to instantiate one.
    kwargs.update({'prop': {'family': rc['fontname']}})

    # Detect if user wants to specify rows manually
    # Gives huge latitude for user input:
    #   1) user can specify nothing and align will be inferred (list of iterables
    #      will always be False, i.e. we draw consecutive legends, and list of handles is always true)
    #   2) user can specify align (needs list of handles for True, list of handles or list
    #      of iterables for False and if the former, will turn into list of iterables)
    if not handles and not ax._filled:
        handles = ax.get_legend_handles_labels()[0]
        if not handles:
            raise ValueError('No axes artists with labels were found.')
    elif not handles:
        raise ValueError('You must pass a handles list for panel axes "filled" with a legend.')
    for i,handle in enumerate(handles):
        if not hasattr(handle, 'get_cmap') or hasattr(handle, 'get_facecolor'): # latter is for scatter (TODO: add wrapper_cmap for scatter?)
            continue
        # Make sure we sample the *center* of the colormap
        warnings.warn('Getting legend entry from colormap.')
        size = np.mean(handle.get_sizes())
        cmap = handle.get_cmap()
        handles[i] = ax.scatter([0], [0], markersize=size,
                                color=[cmap(0.5)],
                                label=handle.get_label())
    list_of_lists = not isinstance(handles[0], martist.Artist)
    if align is None: # automatically guess
        align = not list_of_lists
    else: # standardize format based on input
        if not align and not list_of_lists: # separate into columns
            ncol = kwargs.pop('ncol', 3)
            handles = [handles[i*ncol:(i+1)*ncol]
                        for i in range(len(handles))] # to list of iterables
            list_of_lists = True
        elif not align and list_of_lists and 'ncol' in kwargs:
            kwargs.pop('ncol')
            warnings.warn('Detected list of *lists* of legend handles. Ignoring user input property "ncol".')
        if align and list_of_lists: # unfurl, because we just want one legend!
            handles = [handle for sublist in handles for handle in sublist]
            list_of_lists = False # no longer is list of lists
    # Remove empty lists... pops up in some examples, not sure how
    handles = [sublist for sublist in handles if sublist]

    # Now draw legend, with two options
    # 1) Normal legend, just draw everything like normal and columns
    # will be aligned; we re-order handles to be row-major, is only difference
    if align:
        # Optionally change order
        # See: https://stackoverflow.com/q/10101141/4970632
        if 'ncol' not in kwargs:
            kwargs['ncol'] = 3
        if order=='C':
            newhandles = []
            ncol = kwargs['ncol'] # number of columns
            handlesplit = [handles[i*ncol:(i+1)*ncol] for i in range(len(handles)//ncol+1)] # split into rows
            nrowsmax, nfinalrow = len(handlesplit), len(handlesplit[-1]) # max possible row count, and columns in final row
            # e.g. if 5 columns, but final row length 3, columns 0-2 have N rows but 3-4 have N-1 rows
            nrows = [nrowsmax]*nfinalrow + [nrowsmax-1]*(kwargs['ncol']-nfinalrow)
            for col,nrow in enumerate(nrows): # iterate through cols
                newhandles.extend(handlesplit[row][col] for row in range(nrow))
            handles = newhandles
        # Finally draw legend, mimicking row-major ordering
        leg = super(BaseAxes, ax).legend(handles=handles, **kwargs)
        legends = [leg]

    # 2) Separate legend for each row
    # The label spacing/border spacing will be exactly replicated, as if we were
    # using the original legend command
    # Means we also have to overhaul some settings
    else:
        # Warn when user input props are overridden
        overridden = []
        loc = kwargs.pop('loc', 'upper center')
        loc = {0:'best',
            1:'upper right',
            2:'upper left',
            3:'lower left',
            4:'lower right',
            5:'right',
            6:'center left',
            7:'center right',
            8:'lower center',
            9:'upper center',
            10:'center'}.get(loc, loc)
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
                (ax.figure.get_figheight() * np.diff(ax._position.intervaly))
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
            leg = super(BaseAxes, ax).legend(handles=hs, ncol=len(hs), loc=loc,
                frameon=False,
                bbox_transform=ax.transAxes,
                bbox_to_anchor=bbox,
                **kwargs) # _format_legend is overriding original legend Method
            legends.append(leg)
        for l in legends[:-1]:
            ax.add_artist(l) # because matplotlib deletes previous ones

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

def colorbar_factory(ax, mappable, values=None,
        extend=None, extendlength=None, label=None,
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
    ax : `~matplotlib.axes.Axes`
        The axes to fill with a colorbar.
    mappable : mappable or list of str or list of plot handles
        There are three options here:

        1. A mappable object. Basically, any object with a `get_cmap` method, like
           the objects returned by `~matplotlib.axes.Axes.contourf` and
           `~matplotlib.axes.Axes.pcolormesh`.
        2. A list of hex strings, color string names, or RGB tuples. From this,
           a colormap will be generated and used with the colorbar. Requires
           `values` is not ``None``.
        3. A list of "plot handles". Basically, any object with a `get_color`
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
    extendlength : None or float or str, optional
        The length of the colorbar "extensions" in *physical units*.
        If float, units are inches. If string,
        units are interpreted by `~proplot.utils.units`.

        This is handy if you have multiple colorbars in one figure.
        With the matplotlib API, it is really hard to get triangle
        sizes to match, because the `extendlength` units are *relative*.
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
    that are wrapped by `~proplot.axes.BaseAxes` (causing errors due to new
    usage).
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
    # See comment under colorbar() method def for PanelAxes class. Will get
    # weird results if axes is a special BaseAxes.
    if isinstance(ax, BaseAxes):
        raise ValueError('The colorbar axes cannot be an instance of proplot.BaseAxes. Must be native matplotlib axes.Axes class.')
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
    if not isinstance(mappable, martist.Artist) and \
        not isinstance(mappable, mcontour.ContourSet):
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
    kwdefault = {'cax':ax, 'orientation':orientation, 'use_gridspec':True, # use space afforded by entire axes
                 'spacing':'uniform', 'extend':extend, 'drawedges':grid} # this is default case unless mappable has special props
    kwdefault.update(kwargs)
    kwargs = kwdefault

    # Option to generate colorbar/colormap from line handles
    # * Note the colors are perfect if we don't extend them by dummy color on either side,
    #   but for some reason labels for edge colors appear offset from everything
    # * Too tired to figure out why so just use this workaround
    if fromcolors: # we passed the colors directly
        colors = mappable
        if values is None:
            raise ValueError('Must pass "values", corresponding to list of colors.')
    if fromlines: # the lines
        if values is None:
            raise ValueError('Must pass "values", corresponding to list of handles.')
        if len(mappable)!=len(values):
            raise ValueError('Number of "values" should equal number of handles.')
        colors = [h.get_color() for h in mappable]
    # Get colors, and by default, label each value directly
    # Note contourf will not be overridden for colorbar axes! Need to
    # manually wrap with wrapper_cmap.
    if fromlines or fromcolors:
        cmap   = colortools.Colormap(colors)
        func = _wrapper_cmap(ax, ax.contourf)
        mappable = func([[0,0],[0,0]],
            values=np.array(values), cmap=cmap, extend='neither',
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
        # where minor ticks extend beyond extendlength.
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
    width, height = ax.figure.get_size_inches()
    if orientation=='horizontal':
        axis = ax.xaxis
        scale = width*np.diff(getattr(ax.get_position(),'intervalx'))[0]
    else:
        axis = ax.yaxis
        scale = height*np.diff(getattr(ax.get_position(),'intervaly'))[0]
    extendlength = utils.units(_default(extendlength, rc.get('colorbar.extend')))
    extendlength = extendlength/(scale - 2*extendlength)
    formatter    = axistools.Formatter(formatter, **formatter_kw)
    kwargs.update({'ticks':locators[0], # WARNING: without this, set_ticks screws up number labels for some reason
                   'format':formatter,
                   'ticklocation':ticklocation,
                   'extendfrac':extendlength})
    cb = ax.figure.colorbar(mappable, **kwargs)
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
def _self_func_wrapper(driver):
    def decorator(self, func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return driver(self, func, *args, **kwargs)
        return wrapper
    return decorator
def _self_wrapper(driver):
    def decorator(self, func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return driver(self, *args, **kwargs)
        return wrapper
    return decorator
def _func_wrapper(driver):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return driver(func, *args, **kwargs)
        return wrapper
    return decorator
# Build wrappers
_wrapper_check_centers = _func_wrapper(wrapper_check_centers)
_wrapper_check_edges   = _func_wrapper(wrapper_check_edges)
_wrapper_basemap_gridfix   = _self_func_wrapper(wrapper_basemap_gridfix)
_wrapper_basemap_latlon    = _self_func_wrapper(wrapper_basemap_latlon)
_wrapper_cartopy_gridfix   = _self_func_wrapper(wrapper_cartopy_gridfix)
_wrapper_cartopy_transform = _self_func_wrapper(wrapper_cartopy_transform)
_wrapper_cartopy_crs       = _self_func_wrapper(wrapper_cartopy_crs)
_wrapper_legend  = _self_wrapper(legend_factory) # special
_wrapper_cmap    = _self_func_wrapper(wrapper_cmap)
_wrapper_cycle   = _self_func_wrapper(wrapper_cycle)
_wrapper_plot    = _self_func_wrapper(wrapper_plot)
_wrapper_scatter = _self_func_wrapper(wrapper_scatter)
_wrapper_text    = _self_func_wrapper(wrapper_text)

#------------------------------------------------------------------------------#
# Generalized custom axes class
# WARNING: Wanted to bulk wrap methods using __new__ on a *metaclass*, since
# this would just wrap method *once* and not every single time user accesses
# object. More elegant, but __new__ *does not receive inherited methods* (that
# comes later down the line), so we can't wrap them. Anyway overriding
# __getattribute__ is fine, and premature optimiztaion is root of all evil!
#------------------------------------------------------------------------------#
class BaseAxes(maxes.Axes):
    """Lowest-level axes subclass. Handles titles and axis
    sharing. Adds several new methods and overrides existing ones."""
    name = 'base'
    """The registered projection name."""
    def __init__(self, *args, number=None,
            sharex=None, sharey=None, spanx=None, spany=None,
            sharex_level=0, sharey_level=0,
            **kwargs):
        """
        Parameters
        ----------
        number : None or int
            The subplot number, used for a-b-c labelling (see
            `~BaseAxes.format`).
        sharex_level, sharey_level : {3, 2, 1, 0}, optional
            The "axis sharing level" for the *x* axis, *y* axis, or both
            axes.
        sharex, sharey : None or `BaseAxes`, optional
            Axes to use for *x* and *y* axis sharing. Should correspond
            to the subplot in the bottommost row, leftmost column.
        spanx, spany : None or `BaseAxes`, optional
            Axes to use for the "spanning" *x* and *y* axis labels. Should
            correspond to the subplot in the leftmost column, bottommost row.

        See also
        --------
        `~proplot.subplots.subplots`,
        `XYAxes`, `CartopyAxes`, `BasemapAxes`,
        `wrapper_cmap`, `wrapper_cycle`
        """
        # Properties
        self.number = number # for abc numbering
        self._hatch = None # background hatching
        self._spanx = spanx # boolean toggles, whether we want to span axes labels
        self._spany = spany
        self._yrange = None # geometry, filled later
        self._xrange = None
        self._nrows = None
        self._ncols = None
        # Ugly but necessary
        self._abc_inside = False
        self._title_inside = False # toggle this to figure out whether we need to push 'super title' up
        self._gridliner_on = False # whether cartopy gridliners are enabled
        # Children and related properties
        self.bottompanel = EmptyPanel()
        self.toppanel    = EmptyPanel()
        self.leftpanel   = EmptyPanel()
        self.rightpanel  = EmptyPanel()
        self._panels_main_gridspec = None # filled with gridspec used for axes subplot and its panels
        self._panels_stack_gridspec = None # filled with gridspec used for 'stacked' panels
        self._tight_bbox = None # save these
        self._zoom = None # the 'zoom lines' for inset zoom-in axes
        self._panel_parent = None
        self._inset_parent = None # filled later
        self._inset_children = [] # arbitrary number of insets possible
        self._colorbar_parent = None
        self._colorbar_child = None # the *actual* axes, with content and whatnot; may be needed for tight subplot stuff
        self._filled = False # turned off when panels filled with colorbar or legend
        self._alty_child = None
        self._altx_child = None
        self._alty_parent = None
        self._altx_parent = None
        self._dualy_scale = None # for scaling units on opposite side of ax, and holding data limits fixed
        self._dualx_scale = None

        # Call parent
        super().__init__(*args, **kwargs)

        # Set up axis sharing, save geometry
        width, height = self.figure.get_size_inches()
        self.width  = np.diff(self._position.intervalx)*width # position is in figure units
        self.height = np.diff(self._position.intervaly)*height
        if isinstance(self, maxes.SubplotBase):
            nrows, ncols, subspec = self._topmost_subspec()
            self._yrange = ((subspec.num1 // ncols) // 2, (subspec.num2 // ncols) // 2)
            self._xrange = ((subspec.num1 % ncols) // 2,  (subspec.num2 % ncols) // 2)
            self._nrows = 1 + nrows // 2 # remember, we have rows masquerading as empty spaces!
            self._ncols = 1 + ncols // 2
        # Axis sharing, title stuff, new text attributes
        self._sharex_setup(sharex, sharex_level)
        self._sharey_setup(sharey, sharey_level)
        self._title_transform = self.title.get_transform() # save in case it changes
        self.abc = self.text(0, 0, '') # position tbd
        self.collabel = self.text(0, 0, '', va='bottom', ha='center', transform=self._title_transform)
        self.rowlabel = self.text(0, 0, '', va='center', ha='right', transform=self.transAxes)
        # Apply 'special' props
        rc._getitem_mode = 0 # might still be non-zero if had error
        self.format(mode=1)

    @_expand_wrapped_methods
    def __getattribute__(self, attr, *args):
        """
        Wraps methods when they are requested by the user. See `wrapper_cmap`,
        `wrapper_cycle`, `wrapper_plot`, `wrapper_scatter`, `wrapper_text`,
        and `legend_factory` (which wraps `~matplotlib.axes.Axes.legend`)
        for more info.

        Disables the redundant methods `_disabled_methods`.

        Enables the attribute aliases ``bpanel`` for ``bottompanel``,
        ``tpanel`` for ``toppanel``, ``lpanel`` for ``leftpanel``, and
        ``rpanel`` for ``rightpanel``.
        """
        for message,attrs in _disabled_methods.items():
            if attr in attrs:
                raise NotImplementedError(message.format(attr))
        attr = _aliases.get(attr, attr)
        obj = super().__getattribute__(attr, *args)
        if attr=='text':
            obj = _wrapper_text(self, obj)
        elif attr=='plot':
            obj = _wrapper_plot(self, obj)
        elif attr=='scatter':
            obj = _wrapper_scatter(self, obj)
        elif attr=='legend' and not isinstance(self, PanelAxes):
            obj = _wrapper_legend(self, obj)
        return obj

    def _topmost_subspec(self):
        """Get the top-level SubplotSpec (i.e. the one encompassed by an
        axes and all its panels, if any are present)."""
        subspec = self.get_subplotspec()
        gridspec = subspec.get_gridspec()
        while isinstance(gridspec, mgridspec.GridSpecFromSubplotSpec):
            try:
                subspec = gridspec._subplot_spec
            except AttributeError:
                raise ValueError('The _subplot_spec attribute is missing from this GridSpecFromSubplotSpec. Cannot determine the parent GridSpec rows/columns occupied by this slot.')
            gridspec = subspec.get_gridspec()
        nrows, ncols = gridspec.get_geometry()
        return nrows, ncols, subspec

    def _share_short_axis(self, share, side, level):
        """When sharing main subplots, shares the short axes of their side panels."""
        if isinstance(self, PanelAxes):
            return
        axis = 'x' if side[0] in 'lr' else 'y'
        paxs1 = getattr(self, side + 'panel') # calling this means, share properties on this axes with input 'share' axes
        paxs2 = getattr(share, side + 'panel')
        if not all(pax and pax.get_visible() and not pax._filled for pax in paxs1) or \
           not all(pax and pax.get_visible() and not pax._filled for pax in paxs2):
            return
        if len(paxs1) != len(paxs2):
            raise RuntimeError('Sync error. Different number of stacked panels along axes on like column/row of figure.')
        for pax1,pax2 in zip(paxs1,paxs2):
            getattr(pax1, '_share' + axis + '_setup')(pax2, level)

    def _share_long_axis(self, share, side, level):
        """When sharing main subplots, shares the long axes of their side panels,
        assuming long axis sharing is enabled for that panel."""
        if isinstance(self, PanelAxes):
            return
        axis = 'x' if side[0] in 'tb' else 'y'
        paxs = getattr(self, side + 'panel') # calling this means, share properties on this axes with input 'share' axes
        if not all(pax and pax.get_visible() and not pax._filled for pax in paxs) or \
           not all(pax._share for pax in paxs):
            return
        for pax in paxs:
            getattr(pax, '_share' + axis + '_setup')(share, level)

    def _sharex_setup(self, sharex, level):
        """Sets up shared axes. The input is the 'parent' axes, from which
        this one will draw its properties."""
        if sharex is None or sharex is self:
            return
        if isinstance(self, MapAxes) or isinstance(sharex, MapAxes):
            return
        if level not in range(4):
            raise ValueError('Level can be 1 (do not share limits, just hide axis labels), 2 (share limits, but do not hide tick labels), or 3 (share limits and hide tick labels).')
        # Account for side panels 
        self._share_short_axis(sharex, 'left',   level)
        self._share_short_axis(sharex, 'right',  level)
        self._share_long_axis(sharex,  'bottom', level)
        self._share_long_axis(sharex,  'top',    level)
        # Builtin features
        self._sharex = sharex
        if level>1:
            self._shared_x_axes.join(self, sharex)
        # "Shared" axis and tick labels
        # WARNING: Assigning *another* axis label to this axes will raise error,
        # because matplotlib tries to draw same Artist twice. Just make it invisible.
        if level>2:
            for t in self.xaxis.get_ticklabels():
                t.set_visible(False)
        self.xaxis.label.set_visible(False)

    def _sharey_setup(self, sharey, level):
        """Sets up shared axes. The input is the 'parent' axes, from which
        this one will draw its properties."""
        if sharey is None or sharey is self:
            return
        if isinstance(self, MapAxes) or isinstance(sharey, MapAxes):
            return
        if level not in range(4):
            raise ValueError('Level can be 1 (do not share limits, just hide axis labels), 2 (share limits, but do not hide tick labels), or 3 (share limits and hide tick labels).')
        # Account for side panels
        self._share_short_axis(sharey, 'bottom', level)
        self._share_short_axis(sharey, 'top',    level)
        self._share_long_axis(sharey,  'left',   level)
        self._share_long_axis(sharey,  'right',  level)
        # Builtin features
        self._sharey = sharey
        if level>1:
            self._shared_y_axes.join(self, sharey)
        # "Shared" axis and tick labels
        if level>2:
            for t in self.yaxis.get_ticklabels():
                t.set_visible(False)
        self.yaxis.label.set_visible(False)

    def _parse_title_args(self, abc=False, pos=None, border=None, linewidth=None, **kwargs):
        """Position title text to the left, center, or right and either
        inside or outside the axes (default is center, outside)."""
        # Maybe all this crap was done already in first call to format
        if pos is None:
            return kwargs
        # Determine position
        ypad = rc.get('axes.titlepad')/(72*self.height) # to inches --> to axes relative
        xpad = rc.get('axes.titlepad')/(72*self.width)  # why not use the same for x?
        xpad_i, ypad_i = xpad*1.5, ypad*1.5 # inside labels need a bit more room
        if not isinstance(pos, str):
            # Coordinate position
            ha = va = 'center'
            x, y = pos
            transform = self.transAxes
        else:
            # Get horizontal position
            if not any(c in pos for c in 'lcr'):
                pos += 'c'
            if not any(c in pos for c in 'oi'):
                pos += 'o'
            if 'c' in pos:
                x = 0.5
                ha = 'center'
            elif 'l' in pos:
                x = 0 + xpad_i*('i' in pos)
                ha = 'left'
            elif 'r' in pos:
                x = 1 - xpad_i*('i' in pos)
                ha = 'right'
            # Get vertical position
            if 'o' in pos:
                y = 1 # leave it alone, may be adjusted during draw-time to account for axis label (fails to adjust for tick labels; see notebook)
                va = 'bottom' # matches title alignment maybe
                transform = self._title_transform
                kwargs['border'] = False
            elif 'i' in pos:
                y = 1 - ypad_i
                va = 'top'
                transform = self.transAxes
                if border is not None:
                    kwargs['border'] = border
                if border and linewidth:
                    kwargs['linewidth'] = linewidth
        # Record _title_inside so we can automatically deflect suptitle
        # If *any* object is outside (title or abc), want to deflect it up
        inside = (not isinstance(pos, str) or 'i' in pos)
        if abc:
            self._abc_inside = inside
        else:
            self._title_inside = inside
        return {'x':x, 'y':y, 'ha':ha, 'va':va, 'transform':transform, **kwargs}

    def format(self, rc_kw=None, mode=2, **kwargs):
        """
        Sets up temporary rc settings and calls `~BaseAxes.smart_update`.

        Parameters
        ----------
        rc_kw : None or dict, optional
            A dictionary containing "rc" configuration settings that will
            be applied to this axes. Temporarily updates the
            `~proplot.rcmod.rc` object. See `~proplot.rcmod` for details.
        **kwargs
            Any of three options:

            * A keyword arg for `BaseAxes.smart_update`,
              `XYAxes.smart_update`, or `MapAxes.smart_update`.
            * A global "rc" keyword arg, like ``linewidth`` or ``color``.
            * A standard "rc" keyword arg **with the dots omitted**.
              For example, ``land.color`` becomes ``landcolor``.

            The latter two options update the `~proplot.rcmod.rc`
            object, just like `rc_kw`.

        Other parameters
        ----------------
        mode : int, optional
            The "getitem mode". This is used under-the-hood; you shouldn't
            have to use it directly. Determines whether queries to the
            `~proplot.rcmod.rc` object will ignore :ref:`rcParams`.
            This can help prevent a massive number of unnecessary lookups
            when the settings haven't been changed by the user.
            See `~proplot.rcmod.rc_configurator` for details.

        See also
        --------
        `~proplot.rcmod`,
        `BaseAxes.smart_update`, `XYAxes.smart_update`, `MapAxes.smart_update`
        """
        # Figure out which kwargs are valid rc settings
        # WARNING: First part will fail horribly if mode is not zero!
        # WARNING: If declare {} as kwarg, can get filled by user, so default
        # changes! Need to make default None.
        rc_kw = rc_kw or {}
        kw_extra = {}
        for key,value in kwargs.items():
            key_fixed = _rc_names_nodots.get(key, None)
            if key_fixed is None:
                kw_extra[key] = value
            else:
                rc_kw[key_fixed] = value
        rc._getitem_mode = 0 # might still be non-zero if had error
        # Apply custom settings, potentially
        # TODO: Put this somewhere else?
        kw = {}
        if isinstance(self, PanelAxes) and self._flush:
            if self._side=='right':
                kw.update({'yloc':'right', 'ylabelloc':'right', 'yticklabelloc':'right'})
            elif self._side=='top':
                kw.update({'xloc':'top', 'xlabelloc':'top', 'xticklabelloc':'top'})
        kw.update(kw_extra)
        with rc.context(rc_kw, mode=mode):
            self.smart_update(**kw)

    def smart_update(self, title=None, abc=None,
        figtitle=None, suptitle=None, collabels=None, rowlabels=None, # label rows and columns
        top=True, # nopanel optionally puts title and abc label in main axes
        ):
        """
        Function for formatting axes titles and labels.

        Parameters
        ----------
        title : None or str, optional
            The axes title.
        abc : None or bool, optional
            Whether to apply "a-b-c" subplot labelling based on the
            `number` attribute.

            If `number` is >26, labels will loop around to a, ..., z, aa,
            ..., zz, aaa, ..., zzz, ... God help you if you ever need that
            many labels.
        top : bool, optional
            Whether to try to put title and a-b-c label above the top
            axes panel, if it exists (``True``), or to always put them on
            the main subplot (``False``). Defaults to ``True``.
        rowlabels, colllabels : None or list of str, optional
            The subplot row and column labels. If list, length must match
            the number of subplot rows, columns.
        figtitle, suptitle : None or str, optional
            The figure "super" title, centered between the left edge of
            the lefmost column of subplots and the right edge of the rightmost
            column of subplots, and automatically offset above figure titles.

            This is more sophisticated than matplotlib's builtin "super title",
            which is just centered between the figure edges and offset from
            the top edge.

        See also
        --------
        `~proplot.subplots.subplots`, `~proplot.rcmod`,
        `XYAxes.smart_update`, `MapAxes.smart_update`
        """
        # Figure patch (for some reason needs to be re-asserted even if
        # declared before figure is drawn)
        # Look into `~matplotlib.axes.SubplotBase.is_last_row` and
        # `~matplotlib.axes.SubplotBase.is_first_column` methods.
        kw = rc.fill({'facecolor':'figure.facecolor'})
        self.figure.patch.update(kw)

        # Super title and labels
        # NOTE: These are actually *figure-wide* settings, but that line seems
        # to get blurred -- where we have shared axes, spanning labels, and
        # whatnot. May result in redundant assignments if formatting more than
        # one axes, but operations are fast so some redundancy is nbd.
        fig = self.figure # the figure
        suptitle = figtitle or suptitle
        if suptitle is not None:
            kw = rc.fill({
                'fontsize': 'suptitle.fontsize',
                'weight':   'suptitle.weight',
                'color':    'suptitle.color',
                'fontname': 'fontname'
                }, cache=False)
            fig._suptitle_setup(suptitle, **kw)
        if rowlabels is not None:
            kw = rc.fill({
                'fontsize': 'rowlabel.fontsize',
                'weight':   'rowlabel.weight',
                'color':    'rowlabel.color',
                'fontname': 'fontname'
                }, cache=False)
            fig._rowlabels(rowlabels, **kw)
        if collabels is not None:
            kw = rc.fill({
                'fontsize': 'collabel.fontsize',
                'weight':   'collabel.weight',
                'color':    'collabel.color',
                'fontname': 'fontname'
                }, cache=False)
            fig._collabels(collabels, **kw)

        # Create axes title
        # WARNING: Aligning title flush against left or right of axes is
        # actually already a matplotlib feature! Use set_title(loc=loc), and
        # it assigns text to a different hidden object. My version is just
        # more flexible, allows specifying arbitrary postiion.
        if title is not None:
            kw = rc.fill({
                'pos':       'title.pos',
                'border':    'title.border',
                'linewidth': 'title.linewidth',
                'fontsize':  'title.fontsize',
                'weight':    'title.weight',
                'fontname':  'fontname'
                }, cache=False)
            pax = self.toppanel[0]
            if top and pax and pax.get_visible() and not pax._filled:
                ax = self.toppanel[0]
            else:
                ax = self
            ax = ax._altx_child or ax # always on top!
            kw = ax._parse_title_args(**kw)
            ax.title = _update_text(ax.title, {'text':title, 'visible':True, **kw})

        # Create axes numbering
        if self.number is not None and abc:
            # Position and format
            abcformat = rc.get('abc.format')
            if not ('a' in abcformat or 'A' in abcformat):
                raise ValueError(f'Invalid abcformat {abcformat}.')
            abcedges = abcformat.split('a' if 'a' in abcformat else 'A')
            text = abcedges[0] + _abc(self.number-1) + abcedges[-1]
            if 'A' in abcformat:
                text = text.upper()
            # Settings, make
            kw = rc.fill({
                'pos':       'abc.pos',
                'border':    'abc.border',
                'linewidth': 'abc.linewidth',
                'fontsize':  'abc.fontsize',
                'weight':    'abc.weight',
                'color':     'abc.color',
                'fontname':  'fontname'
                }, cache=False)
            pax = self.toppanel[0]
            if top and pax and pax.get_visible() and not pax._filled:
                ax = self.toppanel[0]
            else:
                ax = self
            ax = ax._altx_child or ax # always on top!
            kw = ax._parse_title_args(abc=True, **kw)
            ax.abc = _update_text(ax.abc, {'text':text, 'visible':True, **kw})
        else:
            for ax in (self, *self.toppanel):
                if ax and ax.abc and not abc and abc is not None:
                    ax.abc.set_visible(False)

    def colorbar(self, *args, loc=None, xspace=None, pad=None, width=None,
        length=None, extendlength=None, label=None, clabel=None, 
        **kwargs):
        """
        Adds an *inset* colorbar, sort of like `~matplotlib.axes.Axes.legend`.

        Parameters
        ----------
        loc : None or str or int, optional
            The colorbar location. Just like ``loc`` for the native matplotlib
            `~matplotlib.axes.Axes.legend`, but filtered to only corner
            positions. Options are:

                * ``1`` or ``'upper right'``
                * ``2`` or ``'upper left'``
                * ``3`` or ``'lower left'``
                * ``4`` or ``'lower right'``

            Default is from the rc configuration.
        xspace : None or str or float, optional
            Space allocated for the bottom x-label of the colorbar.
            If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`. If ``None``,
            read from `~proplot.rcmod.rc` configuration.
        pad : None or str or float, optional
            Space between the axes edge and the colorbar.
            If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`. If ``None``,
            read from `~proplot.rcmod.rc` configuration.
        width : None or str or float, optional
            Colorbar width.
            If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`. If ``None``,
            read from `~proplot.rcmod.rc` configuration.
        length : None or str or float, optional
            Colorbar length.
            If float, units are inches. If string,
            units are interpreted by `~proplot.utils.units`. If ``None``,
            read from `~proplot.rcmod.rc` configuration.
        **kwargs, label, clabel, extendlength
            Passed to `colorbar_factory`.
        """
        # Default props
        loc = _default(loc, rc.get('colorbar.loc'))
        extend = units(_default(extendlength, rc.get('colorbar.extendinset')))
        length = units(_default(pad, rc.get('colorbar.length')))/self.width
        width = units(_default(pad, rc.get('colorbar.width')))/self.height
        pad = units(_default(pad, rc.get('colorbar.axespad')))
        xpad = pad/self.width
        ypad = pad/self.height
        # Space for labels
        xspace = units(_default(xspace, rc.get('colorbar.xspace'))) # for x label
        if not label and not clabel:
            xspace -= 1.2*rc.get('font.size')/72
        xspace /= self.height
        # Get location in axes-relative coordinates
        if loc in (1,'upper right'):
            bounds = (1-xpad-length, 1-ypad-width, length, width)
        elif loc in (2,'upper left'):
            bounds = (xpad, 1-ypad-width, length, width) # x0, y0, width, height in axes-relative coordinate to start
        elif loc in (3,'lower left'):
            bounds = (xpad, ypad+xspace, length, width)
        elif loc in (4,'lower right'):
            bounds = (1-xpad-length, ypad+xspace, length, width)
        else:
            raise ValueError(f'Invalid location {loc}.')
        # Test
        bbox = mtransforms.Bbox.from_bounds(*bounds)
        bb = mtransforms.TransformedBbox(bbox, self.transAxes)
        tr = self.figure.transFigure.inverted()
        bb = mtransforms.TransformedBbox(bb, tr)
        # Axes
        locator = self._make_inset_locator(bounds, self.transAxes)
        bb = locator(None, None)
        ax = maxes.Axes(self.figure, bb.bounds, zorder=5)
        ax.set_axes_locator(locator)
        self.add_child_axes(ax)
        kwargs.update({'ticklocation':'bottom', 'extendlength':extend, 'label':label, 'clabel':clabel})
        return colorbar_factory(ax, *args, **kwargs)

    def cmapline(self, *args, values=None,
            cmap=None, norm=None,
            interp=0, **kwargs):
        """
        Creates a "colormap line", i.e. a line whose color changes as a function
        of the coordinate `values`. This is actually a collection of lines,
        added as a `~matplotlib.collections.LineCollection` instance. See `this
        matplotlib example <https://matplotlib.org/gallery/lines_bars_and_markers/multicolored_line.html>`__.
        This method is invoked if you call `~matplotlib.axes.Axes.plot` with
        the ``cmap`` keyword arg.

        Parameters
        ----------
        values : list of float
            Values corresponding to points on the line.
        cmap : None or colormap spec, optional
            Colormap specifier, passed to `~proplot.colortools.Colormap`.
        norm : None or `~matplotlib.colors.Normalize`, optional
            The normalizer, used for mapping `values` to colormap colors.
        interp : int, optional
            Number of values between each line joint and each *halfway* point
            between line joints to which you want to interpolate.

        """
        # WARNING: So far this only works for 1D *x* and *y* coordinates. Cannot
        # draw multiple colormap lines at once, unlike `~matplotlib.axes.Axes.plot`.
        # First error check
        if values is None:
            raise ValueError('Requires a "values" keyword arg.')
        if len(args) not in (1,2):
            raise ValueError(f'Requires 1-2 arguments, got {len(args)}.')
        y = np.array(args[-1]).squeeze()
        x = np.arange(y.shape[-1]) if len(args)==1 else np.array(args[0]).squeeze()
        values = np.array(values).squeeze()
        if x.ndim!=1 or y.ndim!=1 or values.ndim!=1:
            raise ValueError(f'x ({x.ndim}-d), y ({y.ndim}-d), and values ({values.ndim}-d) must be 1-dimensional.')
        if len(x)!=len(y) or len(x)!=len(values) or len(y)!=len(values):
            raise ValueError(f'{len(x)} xs, {len(y)} ys, but {len(values)} colormap values.')

        # Next draw the line
        # Interpolate values to optionally allow for smooth gradations between
        # values (bins=False) or color switchover halfway between points (bins=True)
        # Next optionally interpolate the corresponding colormap values
        # NOTE: We linearly interpolate here, but user might use a normalizer that
        # e.g. performs log before selecting linear color range; don't need to
        # implement that here
        if interp>0:
            xorig, yorig, vorig = x, y, values
            x, y, values = [], [], []
            for j in range(xorig.shape[0]-1):
                idx = (slice(None, -1) if j+1<xorig.shape[0]-1 else slice(None))
                x.extend(np.linspace(xorig[j], xorig[j+1], interp + 2)[idx].flat)
                y.extend(np.linspace(yorig[j], yorig[j+1], interp + 2)[idx].flat)
                values.extend(np.linspace(vorig[j], vorig[j+1], interp + 2)[idx].flat)
            x, y, values = np.array(x), np.array(y), np.array(values)
        coords, vals = [], []
        levels = utils.edges(values)
        for j in range(y.shape[0]):
            # Get x/y coordinates and values for points to the 'left' and
            # 'right' of each joint. Also prevent duplicates.
            if j==0:
                xleft, yleft = [], []
            else:
                xleft = [(x[j-1] + x[j])/2, x[j]]
                yleft = [(y[j-1] + y[j])/2, y[j]]
            if j+1==y.shape[0]:
                xright, yright = [], []
            else:
                xleft  = xleft[:-1] # prevent repetition when joined with xright/yright
                yleft  = yleft[:-1] # actually need numbers of x/y coordinates to be same for each segment
                xright = [x[j], (x[j+1] + x[j])/2]
                yright = [y[j], (y[j+1] + y[j])/2]
            pleft  = np.stack((xleft,  yleft), axis=1)
            pright = np.stack((xright, yright), axis=1)
            coords.append(np.concatenate((pleft, pright), axis=0))

        # Create LineCollection and update with values
        # TODO: Why not just pass kwargs to class?
        collection = mcollections.LineCollection(np.array(coords),
                cmap=cmap, norm=norm, linestyles='-', capstyle='butt', joinstyle='miter')
        collection.set_array(np.array(values))
        collection.update({key:value for key,value in kwargs.items() if key not in ('color',)})

        # Add collection, with some custom attributes
        self.add_collection(collection)
        self.autoscale_view() # data limits not updated otherwise
        collection.values = values
        collection.levels = levels # needed for other functions some
        return collection

#------------------------------------------------------------------------------#
# Specific classes, which subclass the base one
#------------------------------------------------------------------------------#
def _rcloc_to_stringloc(xy, string): # figures out string location
    """Gets location string from the boolean rc settings, for a given string
    prefix like ``'axes.spines'`` or ``'xtick'``."""
    if xy=='x':
        top = rc[f'{string}.top']
        bottom = rc[f'{string}.bottom']
        if top is None and bottom is None:
            return None
        elif top and bottom:
            return 'both'
        elif top:
            return 'top'
        elif bottom:
            return 'bottom'
        else:
            return 'neither'
    elif xy=='y':
        left = rc[f'{string}.left']
        right = rc[f'{string}.right']
        if left is None and right is None:
            return None
        elif left and right:
            return 'both'
        elif left:
            return 'left'
        elif right:
            return 'right'
        else:
            return 'neither'
    else:
        raise ValueError(f'"xy" must equal "x" or "y".')

class XYAxes(BaseAxes):
    """
    Axes subclass for ordinary Cartesian axes. Adds several new methods and
    overrides existing ones.

    See also
    --------
    `~proplot.subplots.subplots`, `BaseAxes`, `PanelAxes`
    """
    name = 'xy'
    """The registered projection name."""
    def __init__(self, *args, **kwargs):
        # Create simple x by y subplot.
        super().__init__(*args, **kwargs)
        # Change the default formatter
        # Mine trims trailing zeros, but numbers no longer aligned. Matter
        # of taste really; will see if others like it.
        formatter = axistools.Formatter('default')
        self.xaxis.set_major_formatter(formatter)
        formatter = axistools.Formatter('default')
        self.yaxis.set_major_formatter(formatter)
        # Reset this; otherwise matplotlib won't automatically change
        # formatter when it encounters certain types of data, like datetime.
        self.xaxis.isDefault_majfmt = True
        self.yaxis.isDefault_majfmt = True
        # For tight layout; matplotlib draws tick bbox around ticks, even
        # if they are actually all turned off.
        self._xtick_pad_error = (0,0)
        self._ytick_pad_error = (0,0)

    @_expand_wrapped_methods
    def __getattribute__(self, attr, *args):
        """Wraps methods when they are requested by the user. See
        `wrapper_cmap`, `wrapper_cycle`, `wrapper_check_centers`, and
        `wrapper_check_edges` for more info."""
        obj = super().__getattribute__(attr, *args)
        if attr in _cmap_methods: # must come first!
            obj = _wrapper_cmap(self, obj)
        elif attr in _cycle_methods:
            obj = _wrapper_cycle(self, obj)
        if attr in _centers_methods:
            obj = _wrapper_check_centers(obj)
        elif attr in _edges_methods:
            obj = _wrapper_check_edges(obj)
        return obj

    def smart_update(self,
        xloc=None,          yloc=None,          # aliases for 'where to put spine'
        xmargin=None,       ymargin=None,
        xcolor=None,        ycolor=None,        # separate color for x or y axis spines, ticks, tick labels, and axis labels; useful for twin axes
        xspineloc=None,     yspineloc=None,     # deals with spine options
        xtickloc=None,      ytickloc=None,      # which spines to draw ticks on
        fixticks=False, # whether to always transform locator to FixedLocator
        xlabelloc=None,     ylabelloc=None,
        xticklabelloc=None, yticklabelloc=None, # where to put tick labels
        xtickdir=None,      ytickdir=None,      # which direction ('in', 'out', or 'inout')
        xgrid=None,         ygrid=None,         # gridline toggle
        xgridminor=None,    ygridminor=None,    # minor grids on/off (if ticks off, grid will always be off)
        xtickminor=True,    ytickminor=True,    # minor ticks on/off
        xticklabeldir=None, yticklabeldir=None, # which direction to draw labels
        xtickrange=None,    ytickrange=None,    # limit regions where we assign ticklabels to major-ticks
        xreverse=False,     yreverse=False,     # special properties
        xlabel=None,        ylabel=None,        # axis labels
        xlim=None,          ylim=None,
        xbounds=None,       ybounds=None,       # limit spine bounds?
        xscale=None,        yscale=None,
        xformatter=None, yformatter=None, xticklabels=None, yticklabels=None,
        xticks=None, xminorticks=None, xlocator=None, xminorlocator=None,
        yticks=None, yminorticks=None, ylocator=None, yminorlocator=None, # locators, or derivatives that are passed to locators
        xlabel_kw={}, ylabel_kw={},
        xscale_kw={}, yscale_kw={},
        xlocator_kw={}, ylocator_kw={},
        xformatter_kw={}, yformatter_kw={},
        xminorlocator_kw={}, yminorlocator_kw={},
        **kwargs): # formatter
        """
        Format the *x* and *y* axis labels, tick locations, tick label
        format, scales, and more.

        Parameters
        ----------
        xlabel, ylabel : None or str, optional
            The *x* and *y* axis labels. Applied with
            `~matplotlib.axes.Axes.set_xlabel`
            and `~matplotlib.axes.Axes.set_ylabel`.
        xlabel_kw, ylabel_kw : dict-like, optional
            The *x* and *y* axis label settings. Applied with the
            `~matplotlib.artist.Artist.update` method on the
            `~matplotlib.text.Text` instance. Options include ``'color'``,
            ``'size'``, and ``'weight'``.
        xlim, ylim : None or length-2 list of float or None, optional
            The *x* and *y* axis data limits. Applied with
            `~matplotlib.axes.Axes.set_xlim` and `~matplotlib.axes.Axes.set_ylim`.
        xreverse, yreverse : bool, optional
            Simply reverses the order of the `xlim` or `ylim` tuples if
            they were provided.
        xbounds, ybounds : None or length-2 list of float, optional
            The *x* and *y* axis data bounds within which to draw the spines.
            For example, this can be used to separate the *x*
            and *y* axis spines, so they don't meet at the corner.
        xtickrange, ytickrange : None or length-2 list of float, optional
            The *x* and *y* axis data ranges within which major tick marks
            are labelled. For example, the tick range ``[-1,1]`` with
            axis range ``[-5,5]`` and a tick interval of 1 will only
            label the ticks marks at -1, 0, and 1.
        fixticks : bool, optional
            Whether to always transform the tick locators to a
            `~matplotlib.ticker.FixedLocator` instance. Defaults to ``False``.
            If your axis ticks are doing weird things (for example, ticks
            drawn outside of the axis spine), try setting this to ``True``.
        xscale, yscale : None or scale spec, optional
            The *x* and *y* axis scales. Arguments are passed to and interpreted
            by `~proplot.axistools.Scale`. Examples include ``'linear'``
            and ``('cutoff', 0.5, 2)``.
        xscale_kw, yscale_kw : dict-like, optional
            The *x* and *y* axis scale settings. Passed to
            `~proplot.axistools.Scale`.
        xloc, yloc
            Aliases for `xspineloc`, `yspineloc`.
        xspineloc, yspineloc : {'both', 'bottom', 'top', 'left', 'right', 'neither', 'center', 'zero'}, optional
            The *x* and *y* axis spine locations.
        xtickloc, ytickloc : {'both', 'bottom', 'top', 'left', 'right', 'neither'}, optional
            Which *x* and *y* axis spines should have major and minor tick marks.
        xtickminor, ytickminor : None or bool, optional
            Whether to draw minor ticks on the *x* and *y* axes.
        xtickdir, ytickdir : {'out', 'in', 'inout'}
            Direction that major and minor tick marks point for the *x* and
            *y* axis.
        xgrid, ygrid : None or bool, optional
            Whether to draw major gridlines on the *x* and *y* axis.
        xgridminor, ygridminor : None or bool, optional
            Whether to draw minor gridlines for the *x* and *y* axis.
        xticklabeldir, yticklabeldir : {'in', 'out'}
            Whether to place *x* and *y* axis tick label text inside
            or outside the axes.
        xticks, yticks
            Aliases for `xlocator`, `ylocator`.
        xlocator, ylocator : None or locator spec, optional
            The *x* and *y* axis tick mark positions. Passed to the
            `~proplot.axistools.Locator` constructor.
        xlocator_kw, ylocator_kw : dict-like, optional
            The *x* and *y* axis locator settings. Passed to
            `~proplot.axistools.Locator`.
        xminorticks, yminorticks
            Aliases for `xminorlocator`, `yminorlocator`.
        xminorlocator, yminorlocator
            As for `xlocator`, `ylocator`, but for the minor ticks.
        xminorlocator_kw, yminorlocator_kw
            As for `xlocator_kw`, `ylocator_kw`, but for the minor locator.
        xticklabels, yticklabels
            Aliases for `xformatter`, `yformatter`.
        xformatter, yformatter : None or format spec, optional
            The *x* and *y* axis tick label format. Passed to the
            `~proplot.axistools.Formatter` constructor.
        xformatter_kw, yformatter_kw : dict-like, optional
            The *x* and *y* axis formatter settings. Passed to
            `~proplot.axistools.Formatter`.

        Note
        ----
        If you plot something with a `numpy`
        `datetime64 <https://docs.scipy.org/doc/numpy/reference/arrays.datetime.html>`__,
        `pandas.Timestamp`, `pandas.DatetimeIndex`, `datetime.date`,
        `datetime.time`, or `datetime.datetime` array as the *x* or *y*-axis
        coordinate, the axis ticks and tick labels will be formatted as dates
        by default.

        See also
        --------
        `BaseAxes.format`, `BaseAxes.smart_update`, `~proplot.axistools.Scale`,
        `~proplot.axistools.Locator`, `~proplot.axistools.Formatter`
        """
        # Background patch basics
        self.patch.set_clip_on(False)
        self.patch.set_zorder(-1)
        self.patch.update(rc.fill({'facecolor': 'axes.facecolor', 'alpha': 'axes.alpha'}))

        # Background hatching (useful where we want to highlight invalid data)
        # TODO: Implement for PolarAxes and map axes.
        hatch = rc['axes.hatch']
        if not self._hatch and hatch: # non-empty
            self._hatch = self.fill_between([0,1], 0, 1, zorder=0, # put in back
                linewidth=0, hatch=hatch, # linewidth affects only patch edge; hatch width controlled by hatch.linewidth
                alpha=rc.get('hatch.alpha'), edgecolor=rc.get('hatch.color'),
                facecolor='none', transform=self.transAxes)
        if self._hatch:
            kw = rc.fill({
                'alpha':     'hatch.alpha',
                'edgecolor': 'hatch.color',
                })
            self._hatch.update(kw)
            self._hatch.set_hatch(hatch)

        # Flexible keyword args, declare defaults
        xmargin       = _default(xmargin, rc['axes.xmargin'])
        ymargin       = _default(ymargin, rc['axes.ymargin'])
        xtickdir      = _default(xtickdir, rc['xtick.direction'])
        ytickdir      = _default(ytickdir, rc['ytick.direction'])
        xtickminor    = _default(xtickminor, rc['xtick.minor.visible'])
        ytickminor    = _default(ytickminor, rc['ytick.minor.visible'])
        xspineloc     = _default(xloc, xspineloc, _rcloc_to_stringloc('x', 'axes.spines'))
        yspineloc     = _default(yloc, yspineloc, _rcloc_to_stringloc('y', 'axes.spines'))
        xformatter    = _default(xticklabels, xformatter) # default is just matplotlib version
        yformatter    = _default(yticklabels, yformatter)
        xlocator      = _default(xticks, xlocator) # default is AutoLocator, no setting
        ylocator      = _default(yticks, ylocator)
        xminorlocator = _default(xminorticks, xminorlocator) # default is AutoMinorLocator, no setting
        yminorlocator = _default(yminorticks, yminorlocator)
        # Grid defaults are more complicated
        # NOTE: rcmod will always change *grid* and *grid.which* at the same time.
        # raise ValueError('Setting unavailable. Older matplotlib version?')
        # TODO: Make xgridminor, ygridminor, xtickloc, ytickloc all global rc
        # settings that modify xtick.left, axes.grid.axis, etc., instead of
        # these keyword args. Results in duplicate behavior maybe.
        which = rc['axes.grid.which']
        grid = rc.get('axes.grid', cache=True)
        axis = rc.get('axes.grid.axis') # always need this property
        if which is not None or grid is not None: # only if *one* was changed recently!
            # But no matter what we need *both* to figure out proper xgrid, ygrid arguments
            # NOTE: Should we try to make xgridminor, ygridminor part of thing?
            if grid is None:
                grid = rc.get('axes.grid')
            elif which is None:
                which = rc.get('axes.grid.which')
            xgrid      = _default(xgrid, grid and axis in ('x','both') and which in ('major','both'))
            ygrid      = _default(ygrid, grid and axis in ('y','both') and which in ('major','both'))
            xgridminor = _default(xgridminor, grid and axis in ('x','both') and which in ('minor','both'))
            ygridminor = _default(ygridminor, grid and axis in ('y','both') and which in ('minor','both'))

        # Override for weird bug where title doesn't get automatically offset
        # from ticklabels in certain circumstance; check out notebook
        xlabelloc = _default(xlabelloc, xticklabelloc)
        ylabelloc = _default(ylabelloc, yticklabelloc)
        xtickloc = _default(xtickloc, xticklabelloc, _rcloc_to_stringloc('x', 'xtick'))
        ytickloc = _default(ytickloc, yticklabelloc, _rcloc_to_stringloc('x', 'ytick'))
        if xlabelloc=='both':
            xlabelloc = 'bottom'
        if ylabelloc=='both':
            ylabelloc = 'left'
        if xticklabelloc=='both' and xtickloc!='both':
            xticklabelloc = xtickloc
        if yticklabelloc=='both' and ytickloc!='both':
            yticklabelloc = ytickloc
        if xticklabelloc in ('both','top') and (xlabelloc!='top' or not xlabel): # xtickloc *cannot* be 'top', *only* appears for 'both'
            warnings.warn('This keyword combo causes matplotlib bug where title is not offset from tick labels. Try again with xticklabelloc="bottom" or xlabelloc="top". Defaulting to the former.')
            xticklabelloc = xlabelloc = 'bottom'

        # Begin loop
        for (axis, label, color, margin,
            tickloc, spineloc, ticklabelloc, labelloc,
            bounds, grid, gridminor, tickminor, tickminorlocator,
            lim, reverse, scale, locator, formatter, tickrange, tickdir, ticklabeldir,
            scale_kw, label_kw, formatter_kw, locator_kw, minorlocator_kw) in zip(
            (self.xaxis, self.yaxis), (xlabel, ylabel), (xcolor, ycolor), (xmargin, ymargin),
            (xtickloc, ytickloc), (xspineloc, yspineloc), (xticklabelloc, yticklabelloc), (xlabelloc, ylabelloc),
            (xbounds, ybounds), (xgrid, ygrid), (xgridminor, ygridminor), (xtickminor, ytickminor), (xminorlocator, yminorlocator), # minor ticks
            (xlim, ylim), (xreverse, yreverse), (xscale, yscale), (xlocator, ylocator), (xformatter, yformatter), (xtickrange, ytickrange), (xtickdir, ytickdir), (xticklabeldir, yticklabeldir),
            (xscale_kw, yscale_kw), (xlabel_kw, ylabel_kw), (xformatter_kw, yformatter_kw), (xlocator_kw, ylocator_kw), (xminorlocator_kw, yminorlocator_kw),
            ):
            # Axis label properties
            # Redirect user request to the correct *shared* axes, then
            # to the correct *spanning* axis label.
            name = axis.axis_name
            if label is not None:
                kw = rc.fill({
                    'color':    'axes.edgecolor',
                    'fontname': 'fontname',
                    'fontsize': 'axes.labelsize',
                    'weight':   'axes.labelweight'
                    })
                if axis.get_label_position() == 'top':
                    kw['va'] = 'bottom' # baseline was cramped if no ticklabels present
                if color:
                    kw['color'] = color
                kw.update(label_kw)
                self.figure._axis_label_update(axis, text=label, **kw)

            # Axis scale and limits. These don't have axis-specific setters.
            # If user specified xlocator or ylocator and scale is log, enforce
            # custom formatter; this generally means we want specific tick labels
            # on a log-scale plot, but log formatter overrides this and only shows powers of 10.
            if scale is not None:
                if hasattr(scale, 'name'): # class was passed
                    scale = scale.name
                getattr(self, f'set_{name}scale')(axistools.Scale(scale, **scale_kw))
                if scale in ('log','inverse') and locator is not None and formatter is None:
                    formatter = 'default' # override
            if lim is not None:
                if reverse:
                    lim = lim[::-1]
                getattr(self, f'set_{name}lim')(lim)

            # Fix spines
            kw = rc.fill({
                'linewidth': 'axes.linewidth',
                'color':     'axes.edgecolor',
                })
            if color is not None:
                kw['color'] = color
            sides = ('bottom','top') if name=='x' else ('left','right')
            spines = [self.spines[s] for s in sides]
            for spine,side in zip(spines,sides):
                # Line properties
                # Override if we're settings spine bounds
                spineloc = getattr(self, f'{name}spine_override', spineloc) # optionally override; necessary for twinx/twiny situation
                if bounds is not None and spineloc not in sides:
                    spineloc = sides[0] # by default, should just have spines on edges in this case
                # Eliminate sides
                if spineloc=='neither':
                    spine.set_visible(False)
                elif spineloc=='both':
                    spine.set_visible(True)
                elif spineloc in sides: # make relevant spine visible
                    b = True if side==spineloc else False
                    spine.set_visible(b)
                elif spineloc is not None:
                    # Special spine location
                    # Note special 'spine location' options include 'zero', 'center',
                    # and tuple with (units, location) where units can be axes, data, or outward
                    if side==sides[0]: # move the left/semabottom spine onto the specified location, with set_position
                        spine.set_visible(True)
                        spine.set_position(spineloc)
                    else:
                        spine.set_visible(False)
                # Apply spine bounds
                if bounds is not None and spine.get_visible():
                    spine.set_bounds(*bounds)
                spine.update(kw)
            # Get which spines are visible; needed for setting tick locations
            spines = [side for side,spine in zip(sides,spines) if spine.get_visible()]

            # Tick and ticklabel properties
            # * Weird issue seems to cause set_tick_params to reset/forget that the grid
            #   is turned on if you access tick.gridOn directly, instead of passing through tick_params.
            #   Since gridOn is undocumented feature, don't use it. So calling _format_axes() a second time will remove the lines
            # * Can specify whether the left/right/bottom/top spines get ticks; sides will be 
            #   group of left/right or top/bottom
            # * Includes option to draw spines but not draw ticks on that spine, e.g.
            #   on the left/right edges
            # First determine tick sides
            kw_sides = {}
            translate = {None: None, 'both': sides, 'neither': (), 'none': ()}
            if bounds is not None and tickloc not in sides:
                tickloc = sides[0] # override to just one side
            ticklocs = translate.get(tickloc, (tickloc,))
            if ticklocs is not None:
                kw_sides.update({side: (side in ticklocs) for side in sides})
            kw_sides.update({side: False for side in sides if side not in spines}) # override
            # Next the tick label sides
            # Will override to make sure only appear where ticks are
            ticklabellocs = translate.get(ticklabelloc, (ticklabelloc,))
            if ticklabellocs is not None:
                kw_sides.update({f'label{side}': (side in ticklabellocs) for side in sides})
            kw_sides.update({f'label{side}': False for side in sides
                if (side not in spines or (ticklocs is not None and side not in ticklocs))}) # override
            # The axis label side
            if labelloc is None:
                if ticklocs is not None:
                    options = [side for side in sides if (side in ticklocs and side in spines)]
                    if len(options)==1:
                        labelloc = options[0]
            elif labelloc not in sides:
                raise ValueError(f'Got labelloc "{labelloc}", valid options are {sides}.')
            if labelloc is not None:
                axis.set_label_position(labelloc)
            # Style of ticks
            kw_shared = rc.fill({
                # 'labelcolor': f'{name}tick.color',
                # 'labelsize':  f'{name}tick.labelsize'
                'labelcolor': f'tick.labelcolor', # new props
                'labelsize':  f'tick.labelsize',
                'color':      f'{name}tick.color',
                })
            if color:
                kw_shared['color'] = color
                kw_shared['labelcolor'] = color
            if tickdir=='in':
                kw_shared['pad'] = 1 # ticklabels should be much closer
            if ticklabeldir=='in': # put tick labels inside the plot; TODO check this
                tickdir = 'in'
                pad = rc.get(f'{name}tick.major.size') + rc.get(f'{name}tick.major.pad') + rc.get(f'{name}tick.labelsize')
                kw_shared['pad'] = -pad
            if tickdir is not None:
                kw_shared['direction'] = tickdir
            # Apply settings. Also apply gridline settings and major
            # and minor tick settings.
            dict_ = lambda prefix: {
                'grid_color':     f'{prefix}.color',
                'grid_alpha':     f'{prefix}.alpha',
                'grid_linewidth': f'{prefix}.linewidth',
                'grid_linestyle': f'{prefix}.linestyle',
                }
            override = getattr(self, 'grid_override', None)
            # The override should just be a "new default", but user should
            # be able to override that.
            # grid = _default(override, grid)
            # gridminor = _default(override, gridminor)
            grid = _default(grid, override)
            gridminor = _default(gridminor, override)
            for which,igrid in zip(('major', 'minor'), (grid,gridminor)):
                # Turn grid on
                if igrid is not None:
                    axis.grid(igrid, which=which) # toggle with special global props
                # Grid and tick style
                if which=='major':
                    kw_grid = rc.fill(dict_('grid'))
                else:
                    kw_grid_major = kw_grid
                    kw_grid = rc.fill(dict_('gridminor'))
                    kw_grid.update({key:value for key,value in kw_grid_major.items() if key not in kw_grid})
                kw_tick = rc[f'{name}tick.{which}']
                kw_tick.pop('visible', None) # invalid setting
                axis.set_tick_params(which=which, **kw_tick, **kw_grid)
            axis.set_tick_params(which='both', **kw_sides, **kw_shared)
            # Settings that can't be controlled by set_tick_params
            kw = rc.fill({'fontname': 'fontname', 'weight':'tick.labelweight'})
            for t in axis.get_ticklabels():
                t.update(kw)

            # Update margins
            if margin is not None:
                self.margins(**{name:margin})
            # Major and minor locator
            # Also automatically detect whether axis is a 'time axis' (i.e.
            # whether user has plotted something with x/y as datetime/date/np.datetime64
            # objects, and matplotlib automatically set the unit converter)
            time = isinstance(axis.converter, mdates.DateConverter)
            if locator is not None:
                locator = axistools.Locator(locator, time=time, **locator_kw)
                axis.set_major_locator(locator)
                if isinstance(locator, mticker.IndexLocator):
                    tickminor = False # minor ticks make no sense for 'index' data
            if not tickminor and tickminorlocator is None:
                axis.set_minor_locator(axistools.Locator('null'))
            elif tickminorlocator is not None:
                axis.set_minor_locator(axistools.Locator(tickminorlocator,
                    minor=True, time=time, **minorlocator_kw))
            # Major and minor formatter
            fixedformatfix = False
            if formatter is not None or tickrange is not None:
                # Override
                if tickrange is not None:
                    if formatter not in (None, 'default'):
                        warnings.warn('The tickrange feature requires proplot.ScalarFormatter formatter. Overriding input formatter.')
                    formatter = 'default'
                # Set the formatter
                formatter = axistools.Formatter(formatter, tickrange=tickrange, time=time, **formatter_kw)
                axis.set_major_formatter(formatter)
                if isinstance(formatter, mticker.FixedFormatter): # if locator is MultipleLocator, first tick gets cut off!
                    fixedformatfix = True
            axis.set_minor_formatter(mticker.NullFormatter())

            # Ensure no out-of-bounds ticks! Even set_smart_bounds() fails sometimes.
            # Notes
            # * Using set_bounds also failed, and fancy method overrides did
            #   not work, so instead just turn locators into fixed version
            # * Most locators take no arguments in __call__, and some do not
            #   have tick_values method, so we just call them.
            # TODO: Add optional override to do this every time
            if fixticks or fixedformatfix or bounds is not None or axis.get_scale()=='cutoff':
                if bounds is None:
                    bounds = getattr(self, f'get_{name}lim')()
                locator = axistools.Locator([x for x in axis.get_major_locator()() if bounds[0] <= x <= bounds[1]])
                axis.set_major_locator(locator)
                locator = axistools.Locator([x for x in axis.get_minor_locator()() if bounds[0] <= x <= bounds[1]])
                axis.set_minor_locator(locator)

        # Pass stuff to parent formatter, e.g. title and abc labeling
        if (xlim is not None or ylim is not None) and self._inset_parent:
            self.indicate_inset_zoom()
        super().smart_update(**kwargs)

    def dualx(self, offset=0, scale=1, xscale='linear', xlabel=None, **kwargs):
        """As with `~XYAxes.dualy`, but for the *x*-axis. See `~XYAxes.dualy`."""
        parent = self.get_xscale()
        if parent!='linear':
            warnings.warn(f'Parent axis scale must be linear. Overriding current "{parent}" scale.')
            self.set_xscale('linear')
        ax = self.twiny()
        if xlabel is None:
            warnings.warn('Axis label is highly recommended for "alternate units" axis. Use the "xlabel" keyword argument.')
        xscale = axistools.InvertedScaleFactory(xscale)
        ax.format(xscale=xscale, xlabel=xlabel, **kwargs)
        self._dualx_scale = (offset, scale)

    def dualy(self, offset=0, scale=1, yscale='linear', ylabel=None, **kwargs):
        """
        Makes a secondary *y*-axis for denoting equivalent *y*
        coordinates in **alternate units**. Returns nothing.

        Parameters
        ----------
        scale : float, optional
            The constant multiple applied after scaling data with `transform`.
            Defaults to ``1``.
            For example, if your *y*-axis is meters and you
            want kilometers on the other side, use ``scale=1e-3``.
        offset : float, optional
            The constant offset added after multipyling by `scale`.
            Defaults to ``0``.
            For example, if your *y*-axis is Kelvin and you want degrees
            Celsius on the opposite side, use ``offset=-273.15``.
        yscale : str, optional
            The registered scale name used to transform data to the alternate
            units.  Defaults to ``'linear'``.
            For example, if your *y*-axis is wavenumber and you want wavelength on
            the opposite side, use ``yscale='inverse'``. If your-*y* axis
            is height and you want pressure on the opposite side, use
            ``yscale='pressure'`` (and vice versa).
        ylabel : None or str, optional
            The axis label (highly recommended). A warning will be issued if
            this is not supplied.
        **kwargs
            Formats the new axis. Passed to `BaseAxes.format`.

        Note
        ----
        The axis scale `yscale` is used to transform units on the left axis,
        linearly spaced, to units on the right axis. This means the right
        'axis scale' must scale its data with the *inverse* of this transform.
        We make this inverted scale with `~proplot.axistools.InvertedScaleFactory`.
        """
        # Notes:
        # For some reason, when scale is applied, it can change the default
        # formatter. For example, a linear scale will change default formatter
        # to original matplotlib version instead of my custom override. Need
        # to apply it explicitly.
        parent = self.get_yscale()
        if parent!='linear':
            warnings.warn(f'Parent axis scale must be linear. Overriding current "{parent}" scale.')
            self.set_yscale('linear')
        ax = self.twinx()
        if ylabel is None:
            warnings.warn('Axis label is highly recommended for "alternate units" axis. Use the "ylabel" keyword argument.')
        yscale = axistools.InvertedScaleFactory(yscale)
        ax.format(yscale=yscale, ylabel=ylabel, **kwargs)
        self._dualy_scale = (offset, scale)

    def altx(self, *args, **kwargs):
        """Alias (and more intuitive name) for `~XYAxes.twiny`.
        The matplotlib `~matplotlib.axes.Axes.twiny` function
        actually generates two *x*-axes with a shared ("twin") *y*-axis."""
        return self.twiny(*args, **kwargs)

    def alty(self, *args, **kwargs):
        """Alias (and more intuitive name) for `~XYAxes.twinx`.
        The matplotlib `~matplotlib.axes.Axes.twinx` function
        actually generates two *y*-axes with a shared ("twin") *x*-axis."""
        return self.twinx(*args, **kwargs)

    def twinx(self):
        """Mimics matplotlib's `~matplotlib.axes.Axes.twinx` and intelligently handles
        axis ticks, gridlines, axis tick labels, axis labels, and axis sharing.
        Returns an `XYAxes` instance."""
        # Note: Cannot wrap twinx() because then the axes created will be
        # instantiated from the parent class, which doesn't have format method.
        # Instead, use hidden method _make_twin_axes.
        if self._alty_child:
            raise ValueError('No more than two twin axes!')
        if self._alty_parent:
            raise ValueError('This *is* a twin axes!')
        ax = self._make_twin_axes(sharex=self, projection=self.name)
        # Setup
        self.yaxis.tick_left()
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')
        ax.yaxis.set_offset_position('right')
        ax.set_autoscalex_on(self.get_autoscalex_on())
        ax.xaxis.set_visible(False)
        ax.patch.set_visible(False)
        ax.grid(False)
        # Special settings, force spine locations when format called
        self.yspine_override = 'left' # original axis ticks on left
        ax.yspine_override = 'right' # new axis ticks on right
        ax.xspine_override = 'neither'
        ax.grid_override = False
        # Return
        ax._alty_parent = self
        self._alty_child = ax
        return ax

    def twiny(self):
        """Mimics matplotlib's `~matplotlib.axes.Axes.twiny` and intelligently handles
        axis ticks, gridlines, axis tick labels, axis labels, and axis sharing.
        Returns an `XYAxes` instance."""
        # Note: Cannot wrap twiny() because we want to use our own XYAxes,
        # not the matplotlib Axes. Instead use hidden method _make_twin_axes.
        # See https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/axes/_subplots.py
        if self._altx_child:
            raise ValueError('No more than two twin axes!')
        if self._altx_parent:
            raise ValueError('This *is* a twin axes!')
        ax = self._make_twin_axes(sharey=self, projection=self.name)
        # Setup
        self.xaxis.tick_bottom()
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_autoscaley_on(self.get_autoscaley_on())
        ax.yaxis.set_visible(False)
        ax.patch.set_visible(False)
        ax.grid(False)
        # Special settings, force spine locations when format called
        self.xspine_override = 'bottom' # original axis ticks on bottom
        ax.xspine_override = 'top' # new axis ticks on top
        ax.yspine_override = 'neither'
        ax.grid_override = False
        # Return
        ax._altx_parent = self
        self._altx_child = ax
        return ax

    def inset(self, *args, **kwargs):
        """Alias for `~XYAxes.inset_axes`."""
        return self.inset_axes(*args, **kwargs)

    def inset_axes(self, bounds, *, transform=None, zorder=5, zoom=True, zoom_kw={}, **kwargs):
        """Draws an inset `XYAxes` axes. Otherwise, this is a carbon copy
        of the `~matplotlib.axes.Axes.inset_axes` method."""
        # Carbon copy, but use my custom axes
        # Defaults
        if transform is None:
            transform = self.transAxes
        label = kwargs.pop('label', 'inset_axes')
        # This puts the rectangle into figure-relative coordinates.
        # TODO: Use default matplotlib attributes, instead of custom _inset_children
        # and _inset_parent attribute?
        locator = self._make_inset_locator(bounds, transform)
        bb = locator(None, None)
        ax = maxes.Axes(self.figure, bb.bounds, zorder=zorder, label=label, **kwargs)
        # The following locator lets the axes move if in data coordinates, gets called in ax.apply_aspect()
        ax.set_axes_locator(locator)
        self.add_child_axes(ax)
        self._inset_children.append(ax)
        ax._inset_parent = self
        # Finally add zoom (NOTE: Requires version >=3.0)
        if zoom:
            ax.indicate_inset_zoom(**zoom_kw)
        return ax

    def inset_zoom(self, *args, **kwargs):
        """Alias for `~XYAxes.indicate_inset_zoom`."""
        # Just an alias
        return self.indicate_inset_zoom(*args, **kwargs)

    def indicate_inset_zoom(self, alpha=None, linewidth=None, color=None, edgecolor=None, **kwargs):
        """Custom inset zoom indicator that can be *refreshed* as the
        parent axis limits change."""
        # Makes more sense to be defined on the inset axes, since parent
        # could have multiple insets
        parent = self._inset_parent
        alpha = alpha or 1.0
        linewidth = linewidth or rc['axes.linewidth']
        edgecolor = color or edgecolor or rc['axes.edgecolor']
        if not parent:
            raise ValueError(f'{self} is not an inset axes.')
        xlim = self.get_xlim()
        ylim = self.get_ylim()
        rect = [xlim[0], ylim[0], xlim[1] - xlim[0], ylim[1] - ylim[0]]
        kwargs.update({'linewidth':linewidth, 'edgecolor':edgecolor, 'alpha':alpha})
        rectpatch, connects = parent.indicate_inset(rect, self, **kwargs)
        # Adopt properties from old one
        if self._zoom:
            rectpatch_old, connects_old = self._zoom
            rectpatch.update_from(rectpatch_old)
            rectpatch_old.set_visible(False)
            for line,line_old in zip(connects,connects_old):
                # Actually want to *preserve* whether line is visible! This
                # is automatically determined!
                visible = line.get_visible()
                line.update_from(line_old)
                line.set_visible(visible)
                line_old.set_visible(False)
        # By default linewidth is only applied to box
        else:
            for line in connects:
                line.set_linewidth(linewidth)
                line.set_color(edgecolor)
                line.set_alpha(alpha)
        self._zoom = (rectpatch, connects)
        return (rectpatch, connects)

    def _make_inset_locator(self, bounds, trans):
        """Helper function, copied from private matplotlib version."""
        def inset_locator(ax, renderer):
            bbox = mtransforms.Bbox.from_bounds(*bounds)
            bb = mtransforms.TransformedBbox(bbox, trans)
            tr = self.figure.transFigure.inverted()
            bb = mtransforms.TransformedBbox(bb, tr)
            return bb
        return inset_locator

class EmptyPanel(object):
    """
    Replaces `PanelAxes` when the axes or figure panel does not exist.
    This gives a nicer error message than if we had just ``None``, and
    permits indexing to mimick the behavior of a singleton
    `~proplot.subplots.axes_list`.

    Note
    ----
    `__getattr__` is invoked only when `__getattribute__` fails, i.e.
    when the user requests anything that isn't a builtin method.
    """
    def __bool__(self):
        """Returns False. Provides shorthand way to check whether panel
        attribute is specifically EmptyPanel."""
        return False # it's empty, so this is 'falsey'

    def __len__(self, key):
        """Returns 1. This allows us to treat `EmptyPanel` like an
        `~proplot.subplots.axes_list` of stacked panels."""
        return 1

    def __getitem__(self, key):
        """Returns itself. This allows us to treat `EmptyPanel` like an
        `~proplot.subplots.axes_list` of stacked panels."""
        # See: https://stackoverflow.com/a/26611639/4970632
        if key>0:
            raise IndexError('End of panel list.')
        return self

    def __getattr__(self, attr, *args):
        """Raises RuntimeError."""
        raise RuntimeError('Panel does not exist.')

class PanelAxes(XYAxes):
    """`~proplot.axes.XYAxes` subclass, adds `~PanelAxes.legend` and
    `~PanelAxes.colorbar` methods that "fill" the entire axes."""
    # Notes:
    # See `this post <https://stackoverflow.com/a/52121237/4970632>`_
    # and `this example <https://stackoverflow.com/q/26236380/4970632>`_.
    name = 'panel'
    """The registered projection name."""
    def __init__(self, *args,
        side=None, share=False, flush=False,
        visible=True, parent=None, **kwargs):
        """
        Parameters
        ----------
        side : {'left', 'right', 'bottom', 'top'}
            The side on which the panel is drawn.
        share : bool, optional
            Whether to share panel *x* and *y* axes with "parent" axes.
            Irrelevant for figure panels, and defaults to ``False``.
        flush : bool, optional
            Whether panel should always be "flush" against its parent
            subplot or not.
        visible : bool, optional
            Whether to make the axes visible or not.
        parent : None or `~matplotlib.axes.Axes`
            The "parent" of the panel. Not relevant for "outer panel"
            axes.
        *args, **kwargs
            Passed to the `XYAxes` initializer.

        See also
        --------
        `~proplot.subplots.subplots`,
        `~proplot.subplots.Figure.add_subplot_and_panels`,
        `XYAxes`, `BaseAxes`
        """
        # Misc props
        # WARNING: Need to set flush property before calling super init!
        self._share = share
        self._side = side
        self._flush = flush # should panel always be flush against subplot?
        self._parent = parent # used when declaring parent
        # Initialize
        super().__init__(*args, **kwargs)
        if side not in ('left', 'right', 'bottom', 'top'):
            raise ValueError(f'Invalid panel side "{side}".')
        if not visible:
            self.set_visible(False)

    def legend(self, *args, fill=True, **kwargs):
        """
        Use ``fill=True`` (the default) to "fill the panel" with a legend.
        Use ``fill=False`` to add an inset legend with `matplotlib.axes.Axes.legend`.
        If the former, a centered legend is drawn and the axes spines, ticks,
        etc.  are made invisible.

        See `~matplotlib.axes.Axes.legend` and `legend_factory` for details.
        """
        # Regular old inset legend
        if not fill:
            return super().legend(*args, **kwargs)
        # Hide content
        self._filled = True
        for s in self.spines.values():
            s.set_visible(False)
        self.xaxis.set_visible(False)
        self.yaxis.set_visible(False)
        self.patch.set_alpha(0)
        # Allocate invisible axes for drawing legend; by default try to
        # make handles and stuff flush against the axes edge
        kwdefault = {'borderaxespad': 0}
        if not kwargs.get('frameon', rc.get('legend.frameon')):
            kwdefault['borderpad'] = 0
        kwdefault.update(kwargs)
        kwargs = kwdefault
        # Location determined by panel side
        # WARNING: center left and center right also turn off horizontal
        # center alignment, so not an option.
        if 'loc' in kwargs:
            warnings.warn(f'Overriding user input legend property "loc".')
        kwargs['loc'] = {'bottom':'upper center', 'right':'center',
                          'left':'center',  'top':'lower center'}[self._side]
        return legend_factory(self, *args, **kwargs)

    def colorbar(self, *args, fill=True, length=1, **kwargs):
        """
        Use ``fill=True`` (the default) to fill the panel with a colorbar.
        Use ``fill=False`` to add an inset colorbar with `BaseAxes.colorbar`.
        If the former, `length` denotes the fractional extent of the panel
        that is filled.

        See `~matplotlib.figure.Figure.colorbar` and `colorbar_factory`
        for details.
        """
        # Inset 'legend-style' colorbar
        if not fill:
            return super().colorbar(*args, **kwargs)
        # Hide content
        self._filled = True
        for s in self.spines.values():
            s.set_visible(False)
        self.xaxis.set_visible(False)
        self.yaxis.set_visible(False)
        self.patch.set_alpha(0)
        # Draw colorbar with arbitrary length relative to full length of panel
        # TODO: Require the stacked colorbar thing to be declared right away! Will
        # then have panels accessible with the slice panel[i,j] instead of panel[i].
        # space = _default(hspace, wspace, space) # flexible arguments
        fig = self.figure
        side = self._side
        subspec = self.get_subplotspec()
        if side=='top': # this is ugly, and hard to implement with title, super title, and stuff
            raise NotImplementedError('Colorbars in upper panels are not allowed.')
        if length!=1:
            if side in ['bottom']:
                gridspec = FlexibleGridSpecFromSubplotSpec(
                        nrows=1, ncols=3, wspace=0, #hspace=space,
                        subplot_spec=subspec,
                        width_ratios=((1-length)/2, length, (1-length)/2),
                        )
                subspec = gridspec[1]
            elif side in ['left','right']:
                gridspec = FlexibleGridSpecFromSubplotSpec(
                        nrows=3, ncols=1, hspace=0, #wspace=space,
                        subplot_spec=subspec,
                        height_ratios=((1-length)/2, length, (1-length)/2),
                        )
                subspec = gridspec[1]
        # Get properties
        ax = fig.add_subplot(subspec, projection=None)
        if side in ['bottom','top']:
            outside, inside = 'bottom', 'top'
            if side=='top':
                outside, inside = inside, outside
            ticklocation = outside
            orientation  = 'horizontal'
        elif side in ['left','right']:
            outside, inside = 'left', 'right'
            if side=='right':
                outside, inside = inside, outside
            ticklocation = outside
            orientation  = 'vertical'
        # Make colorbar
        self._colorbar_child = ax # these are so far unused
        ax._colorbar_parent = self
        kwargs.update({'orientation':orientation, 'ticklocation':ticklocation})
        return colorbar_factory(ax, *args, **kwargs)

class MapAxes(BaseAxes):
    """Intermediate class, shared by `CartopyAxes` and `BasemapAxes`.
    Disables methods that are inappropriate for map projections and adds
    `MapAxes.smart_update`, so that arguments passed to `~BaseAxes.format` are
    identical for `CartopyAxes` and `BasemapAxes`."""
    def __init__(self, *args, **kwargs): # just to disable docstring inheritence
        """
        See also
        --------
        `~proplot.subplots.subplots`, `BaseAxes`, `CartopyAxes`, `BasemapAxes`
        """
        super().__init__(*args, **kwargs)

    # Disable some methods to prevent weird shit from happening
    # Originally used property decorators for this but way too verbose
    # See: https://stackoverflow.com/a/23126260/4970632
    @_expand_wrapped_methods
    def __getattribute__(self, attr, *args):
        """Disables methods inappropriate for map projections: `_map_disabled_methods`."""
        if attr in _map_disabled_methods:
            raise NotImplementedError('Invalid plotting function {} for map projection axes.'.format(attr))
        return super().__getattribute__(attr, *args)

    # Note this *actually* just returns some standardized arguments
    # to the CartopyAxes.smart_update and BasemapAxes.smart_update methods; they
    # both jump over this intermediate class and call BaseAxes.smart_update
    def smart_update(self, labels=None, latlabels=None, lonlabels=None,
        latmax=None, lonlim=None, latlim=None, grid=None,
        lonlocator=None, lonlines=None, lonticks=None,
        latlocator=None, latlines=None, latticks=None,
        **kwargs,
        ):
        """
        Format the map projection.

        Parameters
        ----------
        labels : None or bool, optional
            Whether to draw longitude and latitude labels. If ``None``, read
            from `~proplot.rcmod.rc` configuration.
        lonlabels, latlabels
            Whether to label longitudes and latitudes, and on which sides
            of the map. There are four different options:

            1. Boolean ``True``. Indicates left side for latitudes,
               bottom for longitudes.
            2. A string, e.g. ``'lr'`` or ``'bt'``.
            3. A boolean ``(left,right)`` tuple for longitudes,
               ``(bottom,top)`` for latitudes.
            4. A boolean ``(n1,n2,n3,n4)`` tuple as in the
               `~mpl_toolkits.basemap.Basemap.drawmeridians` and
               `~mpl_toolkits.basemap.Basemap.drawparallels` methods.
               The boolean values indicate whether to label gridlines intersecting
               the left, right, top, and bottom sides, respectively.

        latmax : None or float, optional
            Meridian gridlines are cut off poleward of this latitude. If
            ``None``, read from the configuration.
        lonlim, latlim : None or length-2 list of float, optional
            Longitude and latitude limits of projection, applied
            with `~cartopy.mpl.geoaxes.GeoAxes.set_extent`.
        grid : None or bool, optional
            Whether to add gridlines.
        lonlines, latlines, lonticks, latticks
            Aliases for `lonlocator`, `latlocator`.
        lonlocator, latlocator : None or list of float, optional
            List of longitudes and latitudes for drawing gridlines.
        **kwargs
            Passed to `BaseAxes.smart_update`.

        See also
        --------
        `BaseAxes.format`, `BaseAxes.smart_update`, `~proplot.subplots.subplots`, `~proplot.rcmod`
        """
        # Parse alternative keyword args
        # NOTE: If labels keyword args were passed, automatically turn grid on
        grid = _default(grid, rc.get('geogrid'))
        labels = _default(labels, rc.get('geogrid.labels')) or bool(lonlabels or latlabels)
        lonlocator = _default(lonlines, lonticks, lonlocator)
        latlocator = _default(latlines, latticks, latlocator)
        latlocator = _default(latlocator, rc.get('geogrid.latlines')) # gridlines by default
        lonlocator = _default(lonlocator, rc.get('geogrid.lonlines'))

        # Interptet latitude
        latmax = _default(latmax, rc.get('geogrid.latmax'))
        if isinstance(self, CartopyAxes):
            lon_0 = self.projection.proj4_params.get('lon_0', 0)
        else:
            lon_0 = self.m.lonmin
        if lonlocator is not None:
            if not np.iterable(lonlocator):
                lonlocator = utils.arange(lon_0 - 180, lon_0 + 180, lonlocator)
            lonlocator = [*lonlocator]
        if latlocator is not None:
            if not np.iterable(latlocator):
                latlocator = utils.arange(-latmax, latmax, latlocator)
            latlocator = [*latlocator]

        # Length-4 boolean arrays of whether and where to goggle labels
        if lonlabels or latlabels:
            labels = True # toggle them at all?
        ilabels = [lonlabels, latlabels]
        for i,(mode,jlabels) in enumerate(zip(('x', 'y'), tuple(ilabels))):
            if jlabels is False:
                return [0]*4
            if jlabels is None:
                jlabels = 1
                # jlabels = False
                # jlabels = True # will label lons on bottom, lats on left
            if isinstance(jlabels, str):
                string = jlabels
                jlabels = [0]*4
                for idx,char in zip([0,1,2,3],'lrbt'):
                    if char in string:
                        jlabels[idx] = 1
            if isinstance(jlabels, Number): # e.g. *boolean*
                jlabels = np.atleast_1d(jlabels)
            if len(jlabels)==1:
                jlabels = [*jlabels, 0] # default is to label bottom/left
            if len(jlabels)==2:
                if mode=='x':
                    jlabels = [0, 0, *jlabels]
                elif mode=='y':
                    jlabels = [*jlabels, 0, 0]
            elif len(jlabels)!=4:
                raise ValueError(f'Invalid {mode} labels: {jlabels}.')
            ilabels[i] = jlabels
        lonlabels, latlabels = ilabels
        return grid, latmax, lonlim, latlim, lonlocator, latlocator, labels, lonlabels, latlabels, kwargs

class PolarAxes(MapAxes, mproj.PolarAxes):
    """Intermediate class, mixes `~matplotlib.projections.polar.PolarAxes`
    with `MapAxes`."""
    def __init__(self, *args, **kwargs):
        """
        See also
        --------
        `~proplot.subplots.subplots`, `BaseAxes`, `MapAxes`

        Warning
        -------
        Polar axes have not been tested yet!
        """
        super().__init__(*args, **kwargs)

    def smart_update(self, *args, **kwargs):
        """Calls `BaseAxes.smart_update`."""
        super(MapAxes, self).smart_update(*args, **kwargs)

    name = 'newpolar'
    """The registered projection name."""

# Cartopy takes advantage of documented feature where any class with method
# named _as_mpl_axes can be passed as 'projection' object.
# Feature documented here: https://matplotlib.org/devel/add_new_projection.html
class CartopyAxes(MapAxes, GeoAxes):
    """
    Axes subclass for plotting `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_
    projections. Initializes the `cartopy.crs.Projection` instance. Also
    allows for *partial* coverage of azimuthal projections by zooming into
    the full projection, then drawing a circle boundary around some latitude
    away from the center (this is surprisingly difficult to do).
    """
    name = 'cartopy'
    """The registered projection name."""
    # Helper
    _n_bounds = 100 # number of points for drawing circle map boundary
    _proj_np = ('npstere',)
    _proj_sp = ('spstere',)
    _proj_circles = ('laea', 'aeqd', 'stere', 'npstere', 'spstere')
    def __init__(self, *args, map_projection=None, centerlat=90, boundinglat=0, **kwargs):
        """
        Parameters
        ----------
        map_projection : `~mpl_toolkits.basemap.Basemap`
            The `~mpl_toolkits.basemap.Basemap` instance.
        centerlat : {90, -90}, optional
            For polar projections, the center latitude of the circle.
        boundinglat : float, optional
            For polar projections, the edge latitude of the circle.
        *args, **kwargs
            Passed to `BaseAxes.__init__`.

        Note
        ----
        The circle stuff for polar projection was developed from `this example
        <https://scitools.org.uk/cartopy/docs/v0.15/examples/always_circular_stereo.html>`_.

        See also
        --------
        `~proplot.proj`, `BaseAxes`, `MapAxes`, `~proplot.subplots.subplots`
        """
        # Dependencies
        import cartopy.crs as ccrs # verify package is available

        # Do the GeoAxes initialization steps manually (there are very few)
        # WARNING: If _hold is set False or None, cartopy will call cla() on
        # axes before plotting stuff, which will wipe out row and column labels
        # even though they appear to stick around; maybe the artist persists
        # but is no longer associated with the axes. Does not matter whether
        # attribute is hidden.
        # self._hold = None # dunno
        if not isinstance(map_projection, ccrs.Projection):
            raise ValueError('You must initialize CartopyAxes with map_projection=<cartopy.crs.Projection>.')
        self._cartopy_gl = None # gridliner
        self.projection = map_projection # attribute used extensively by GeoAxes methods, and by builtin one

        # Call BaseAxes
        super().__init__(*args, map_projection=map_projection, **kwargs)
        # Apply circle boundary
        # self.projection.threshold = kwargs.pop('threshold', self.projection.threshold) # optionally modify threshold
        name = map_projection.proj4_params['proj']
        if name not in self._proj_circles:
            self.set_global() # see: https://stackoverflow.com/a/48956844/4970632
        else:
            if name in self._proj_np:
                centerlat = 90
            elif name in self._proj_sp:
                centerlat = -90
            eps = 1e-3 # had bug with full -180, 180 range when lon_0 was not 0
            center = self.projection.proj4_params['lon_0']
            self.set_extent([center - 180 + eps, center + 180 - eps, boundinglat, centerlat], PlateCarree()) # use platecarree transform
            self.set_boundary(projs.Circle(self._n_bounds), transform=self.transAxes)

    def __getattribute__(self, attr, *args):
        """Wraps methods when they are requested by the user. See
        `wrapper_cmap`, `wrapper_cycle`,
        `wrapper_cartopy_transform`, `wrapper_check_centers`,
        `wrapper_check_edges`, and `wrapper_cartopy_gridfix` for more info."""
        # Wrappers
        obj = super().__getattribute__(attr, *args)
        if attr in _cmap_methods:
            obj = _wrapper_cmap(self, obj)
        elif attr in _cycle_methods:
            obj = _wrapper_cycle(self, obj)
        if attr in _transform_methods:
            obj = _wrapper_cartopy_transform(self, obj)
        elif attr in _crs_methods:
            obj = _wrapper_cartopy_crs(self, obj)
        elif attr in _edges_methods or attr in _centers_methods:
            obj = _wrapper_cartopy_gridfix(self, obj)
            if attr in _edges_methods:
                obj = _wrapper_check_edges(obj)
            else:
                obj = _wrapper_check_centers(obj)
        return obj

    def smart_update(self, grid=None, **kwargs):
        # Documentation inherited from MapAxes
        # Dependencies
        import cartopy.feature as cfeature
        import cartopy.crs as ccrs

        # Parse flexible input
        grid, _, lonlim, latlim, lonlocator, latlocator, labels, lonlabels, latlabels, kwargs = \
                super().smart_update(**kwargs)

        # Projection extent
        # NOTE: They may add this as part of set_xlim and set_ylim in the
        # near future; see: https://github.com/SciTools/cartopy/blob/master/lib/cartopy/mpl/geoaxes.py#L638
        # WARNING: The set_extent method tries to set a *rectangle* between
        # the *4* (x,y) coordinate pairs (each corner), so something like
        # (-180,180,-90,90) will result in *line*, causing error!
        if lonlim is not None or latlim is not None:
            lonlim = lonlim or [None, None]
            latlim = latlim or [None, None]
            lonlim, latlim = [*lonlim], [*latlim]
            lon_0 = self.projection.proj4_params.get('lon_0', 0)
            if lonlim[0] is None:
                lonlim[0] = lon_0 - 180
            if lonlim[1] is None:
                lonlim[1] = lon_0 + 180
            eps = 1e-3 # had bug with full -180, 180 range when lon_0 was not 0
            lonlim[0] += eps
            if latlim[0] is None:
                latlim[0] = -90
            if latlim[1] is None:
                latlim[1] = 90
            self.set_extent([*lonlim, *latlim], PlateCarree())

        # Draw gridlines
        # Make old one invisible
        if self._cartopy_gl is not None:
            gl = self._cartopy_gl
            gl.xlines         = False
            gl.xlabels_top    = False
            gl.xlabels_bottom = False
            gl.ylines         = False
            gl.ylabels_left   = False
            gl.ylabels_right  = False
        # Make new gridliner
        # WARNING: For some reason very weird side effects happen if you try
        # to call gridlines twice on same axes. So impossible to do 'major'
        # and 'minor' gridlines.
        if grid:
            kw = rc.fill({
                'alpha':     'geogrid.alpha',
                'color':     'geogrid.color',
                'linewidth': 'geogrid.linewidth',
                'linestyle': 'geogrid.linestyle',
                }, cache=False)
            labels = labels and isinstance(self.projection, (ccrs.Mercator, ccrs.PlateCarree))
            gl = self.gridlines(draw_labels=labels, zorder=100, **kw)
            self._cartopy_gl = gl

            # Grid locations
            eps = 1e-3
            if latlocator[0]==-90:
                latlocator[0] += eps
            if latlocator[-1]==90:
                latlocator[-1] -= eps
            gl.ylocator = mticker.FixedLocator(latlocator)
            gl.xlocator = mticker.FixedLocator(lonlocator)

            # Grid labels
            if labels:
                from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
                self._gridliner_on = True
                gl.xformatter = LONGITUDE_FORMATTER
                gl.yformatter = LATITUDE_FORMATTER
                gl.xlabels_bottom, gl.xlabels_top = lonlabels[2:]
                gl.ylabels_left, gl.ylabels_right = latlabels[:2]
            else:
                self._gridliner_on = False

        # Geographic features
        # WARNING: Seems cartopy features can't be updated!
        # See: https://scitools.org.uk/cartopy/docs/v0.14/_modules/cartopy/feature.html#Feature
        # Change the _kwargs property also does *nothing*
        # WARNING: Changing linewidth is evidently impossible with cfeature. Bug?
        # See: https://stackoverflow.com/questions/43671240/changing-line-width-of-cartopy-borders
        # NOTE: The natural_earth_shp method is deprecated, use add_feature instead.
        # See: https://cartopy-pelson.readthedocs.io/en/readthedocs/whats_new.html
        # NOTE: The e.g. cfeature.COASTLINE features are just for convenience,
        # hi res versions. Use cfeature.COASTLINE.name to see how it can be looked
        # up with NaturalEarthFeature.
        reso = rc.reso
        if reso not in ('lo','med','hi'):
            raise ValueError(f'Invalid resolution {reso}.')
        reso = {
            'lo':  '110m',
            'med': '50m',
            'hi':  '10m',
            }.get(reso)
        features = {
            'land':         ('physical', 'land'),
            'ocean':        ('physical', 'ocean'),
            'lakes':        ('physical', 'lakes'),
            'coast':        ('physical', 'coastline'),
            'rivers':       ('physical', 'rivers_lake_centerlines'),
            'borders':      ('cultural', 'admin_0_boundary_lines_land'),
            'innerborders': ('cultural', 'admin_1_states_provinces_lakes'),
            }
        for name, args in features.items():
            # Get feature
            if not rc.get(name): # toggled
                continue
            if getattr(self, f'_{name}', None): # already drawn
                continue
            feat = cfeature.NaturalEarthFeature(*args, reso)
            # Customize
            # For 'lines', need to specify edgecolor and facecolor individually
            # See: https://github.com/SciTools/cartopy/issues/803
            kw = {**rc.get(name, nodict=False)}
            if name in ('coast', 'rivers', 'borders', 'innerborders'):
                kw['edgecolor'] = kw.pop('color')
                kw['facecolor'] = 'none'
            else:
                kw['linewidth'] = 0
            self.add_feature(feat, **kw)
            setattr(self, f'_{name}', feat)

        # Update patch
        kw = rc.fill({
            'facecolor': 'axes.facecolor'
            }, cache=False)
        self.background_patch.update(kw)
        kw = rc.fill({
            'edgecolor': 'axes.edgecolor',
            'linewidth': 'axes.linewidth'
            }, cache=False)
        self.outline_patch.update(kw)

        # Pass stuff to parent formatter, e.g. title and abc labeling
        BaseAxes.smart_update(self, **kwargs)

def _ls_translate(obj, style):
    """Make basemap gridlines look like cartopy lines using the `dashes`
    tuple."""
    # Notes:
    # * For some reason basemap gridlines look different from cartopy ones.
    #   Have absolutely **no idea** why. Cartopy seems to do something weird because
    #   there is no ``_dashSeq`` attribute on the lines, and the line styles
    #   are always "officially" ``'-'``, even though they are dashed.
    # * Dots actually look better for cartopy so we mimick them. After testing,
    #   below works *really* well... maybe it is hardcoded in cartopy somewhere?
    # See `this reference <https://matplotlib.org/gallery/lines_bars_and_markers/line_styles_reference.html>`_.
    dashes = [None, None]
    if style==':':
        dashes[0] = 0.05
        dashes[1] = 2.5
    elif style=='--':
        dashes[0] = 2.5
        dashes[1] = 2.5
    elif style!='-':
        raise ValueError(f'Invalid line style {style}.')
    return dashes

class BasemapAxes(MapAxes):
    """
    Axes subclass for plotting `basemap <https://matplotlib.org/basemap/>`_
    projections. The `~mpl_toolkits.basemap.Basemap` projection instance is added as
    the `m` attribute, but this is all abstracted away -- you can use
    `~matplotlib.axes.Axes` methods like `~matplotlib.axes.Axes.plot` and
    `~matplotlib.axes.Axes.contour` with
    your raw longitude-latitude data.
    """
    name = 'basemap'
    """The registered projection name."""
    # Note non-rectangular projections; for rectnagular ones, axes spines are
    # used as boundaries, but for these, have different boundary.
    _proj_non_rectangular = (
            # Always non-rectangular
            'ortho', 'geos', 'nsper',
            'moll', 'hammer', 'robin',
            'eck4', 'kav7', 'mbtfpq', # last one is McBryde-Thomas flat polar quartic
            'sinu', 'vandg', # last one is van der Grinten
            # Only non-rectangular if we pass 'round' kwarg
            # This is done by default, currently no way to change it
            'npstere', 'spstere', 'nplaea',
            'splaea', 'npaeqd', 'spaeqd',
            )
    def __init__(self, *args, map_projection=None, **kwargs):
        """
        Parameters
        ----------
        map_projection : `~mpl_toolkits.basemap.Basemap`
            The `~mpl_toolkits.basemap.Basemap` instance.
        **kwargs
            Passed to `BaseAxes.__init__`.

        See also
        --------
        `~proplot.proj`, `BaseAxes`, `MapAxes`, `~proplot.subplots.subplots`
        """
        # Some notes
        # * Must set boundary before-hand, otherwise the set_axes_limits method called
        #   by mcontourf/mpcolormesh/etc draws two mapboundary Patch objects called "limb1" and
        #   "limb2" automatically: one for fill and the other for the edges
        # * Then, since the patch object in _mapboundarydrawn is only the fill-version, calling
        #   drawmapboundary again will replace only *that one*, but the original visible edges
        #   are still drawn -- so e.g. you can't change the color
        # * If you instead call drawmapboundary right away, _mapboundarydrawn will contain
        #   both the edges and the fill; so calling it again will replace *both*
        import mpl_toolkits.basemap as mbasemap # verify package is available
        if not isinstance(map_projection, mbasemap.Basemap):
            raise ValueError('You must initialize BasemapAxes with map_projection=(basemap.Basemap instance).')
        self.m = map_projection
        self.boundary = None
        self._mapboundarydrawn = None
        self._recurred = False # use this so we can override plotting methods
        self._parallels = None
        self._meridians = None

        # Initialize
        super().__init__(*args, **kwargs)

    # Basemap overrides
    # WARNING: Never ever try to just make blanket methods on the Basemap
    # instance accessible from axes instance! Can of worms and had bunch of
    # weird errors! Just pick the ones you think user will want to use.
    # NOTE: Need to again use cmap and cycle wrapper here, because
    # Basemap internally bypasses my BaseAxes superclass.
    def __getattribute__(self, attr, *args):
        """
        Wraps methods when they are requested by the user. See
        `wrapper_basemap_latlon`, `wrapper_check_centers`,
        `wrapper_check_edges`, and `wrapper_basemap_gridfix` for more info.

        Wraps all plotting methods with ``_wrapper_m_call`` and
        ``_wrapper_m_norecurse``, which (former) calls methods on the
        `~mpl_toolkits.basemap.Basemap` instance and (latter) prevents
        recursion issues arising from internal calls to axes methods by
        the `~mpl_toolkits.basemap.Basemap` instance.
        """
        obj = super().__getattribute__(attr, *args)
        if attr in _latlon_methods or attr in _edges_methods or attr in _centers_methods:
            obj = _wrapper_m_call(self, obj) # this must be the *last* step!
            if attr in _cmap_methods:
                obj = _wrapper_cmap(self, obj)
            elif attr in _cycle_methods:
                obj = _wrapper_cycle(self, obj)
            if attr in _latlon_methods:
                obj = _wrapper_basemap_latlon(self, obj)
            elif attr in _edges_methods or attr in _centers_methods:
                obj = _wrapper_basemap_gridfix(self, obj)
                if attr in _edges_methods:
                    obj = _wrapper_check_edges(obj)
                else:
                    obj = _wrapper_check_centers(obj)
            obj = _wrapper_m_norecurse(self, obj)
        return obj

    def smart_update(self, grid=None, **kwargs):
        # Documentation inherited from MapAxes
        # Parse flexible input
        grid, latmax, lonlim, latlim, lonlocator, latlocator, labels, lonlabels, latlabels, kwargs = \
                super().smart_update(**kwargs)
        if lonlim is not None or latlim is not None:
            warnings.warn('You cannot "zoom into" a basemap projection after creating it. Pass a proj_kw dictionary in your call to subplots, with any of the following basemap keywords: llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, llcrnrx, llcrnry, urcrnrx, urcrnry, width, or height.')

        # Map boundary
        # * First have to *manually replace* the old boundary by just deleting
        #   the original one
        # * If boundary is drawn successfully should be able to call
        #   self.m._mapboundarydrawn.set_visible(False) and edges/fill color disappear
        # * For now will enforce that map plots *always* have background whereas
        #   axes plots can have transparent background
        kw_face = rc.fill({
            'facecolor': 'axes.facecolor'
            }, cache=False)
        kw_edge = rc.fill({
            'linewidth': 'axes.linewidth',
            'edgecolor': 'axes.edgecolor'
            }, cache=False)
        self.axesPatch = self.patch # bugfix or something
        if self.m.projection in self._proj_non_rectangular:
            self.patch.set_alpha(0) # make patch invisible
            if not self.m._mapboundarydrawn:
                p = self.m.drawmapboundary(ax=self) # set fill_color to 'none' to make transparent
            else:
                p = self.m._mapboundarydrawn
            p.update({**kw_face, **kw_edge})
            p.set_rasterized(False) # not sure about this; might be rasterized
            p.set_clip_on(False)    # so edges of *line* denoting boundary aren't cut off
            self.boundary = p       # not sure why this one
        else:
            self.patch.update({**kw_face, 'edgecolor':'none'})
            for spine in self.spines.values():
                spine.update(kw_edge)

        # Longitude/latitude lines
        # Make sure to turn off clipping by invisible axes boundary; otherwise
        # get these weird flat edges where map boundaries, parallel/meridian markers come up to the axes bbox
        if grid:
            lkw = rc.fill({
                'alpha':     'geogrid.alpha',
                'color':     'geogrid.color',
                'linewidth': 'geogrid.linewidth',
                'linestyle': 'geogrid.linestyle',
                }, cache=False)
            tkw = rc.fill({
                'color':    'geogrid.color',
                'fontsize': 'geogrid.labelsize',
                }, cache=False)
            # Latitudes
            # TODO: latlocator and lonlocator are always Truthy, because we
            # if latlocator:
            # if lonlocator:
            # Remove old ones
            if self._parallels:
                for pi in self._parallels.values():
                    for obj in [i for j in pi for i in j]: # magic
                        obj.set_visible(False)
            # Change from left/right/bottom/top to left/right/top/bottom
            if labels:
                latlabels[2:] = latlabels[2:][::-1]
            else:
                latlabels = 4*[0]
            p = self.m.drawparallels(latlocator, latmax=latmax, labels=latlabels, ax=self)
            for pi in p.values(): # returns dict, where each one is tuple
                # Tried passing clip_on to the below, but it does nothing; must set
                # for lines created after the fact
                for obj in [i for j in pi for i in j]: # magic
                    if isinstance(obj, mtext.Text):
                        obj.update(tkw)
                    else:
                        obj.update(lkw)
                        obj.set_dashes(_ls_translate(obj, lkw['linestyle']))
            self._parallels = p

            # Longitudes
            # Remove old ones
            if self._meridians:
                for pi in self._meridians.values():
                    for obj in [i for j in pi for i in j]: # magic
                        obj.set_visible(False)
            # Draw new ones
            if labels:
                lonlabels[2:] = lonlabels[2:][::-1]
            else:
                lonlabels = 4*[0]
            p = self.m.drawmeridians(lonlocator, latmax=latmax, labels=lonlabels, ax=self)
            for pi in p.values():
                for obj in [i for j in pi for i in j]: # magic
                    if isinstance(obj, mtext.Text):
                        obj.update(tkw)
                    else:
                        obj.update(lkw)
                        obj.set_dashes(_ls_translate(obj, lkw['linestyle']))
            self._meridians = p

        # Geography
        # NOTE: Also notable are drawcounties, blumarble, drawlsmask,
        # shadedrelief, and etopo methods.
        features = {
            'land':      'fillcontinents',
            'coast':     'drawcoastlines',
            'rivers':    'drawrivers',
            'borders':   'drawcountries',
            'innerborders': 'drawstates',
            }
        for name, method in features.items():
            if not rc.get(name): # toggled
                continue
            if getattr(self, f'_{name}', None): # already drawn
                continue
            feat = getattr(self.m, method)(ax=self, **rc.get(name, nodict=False))
            setattr(self, f'_{name}', feat)

        # Pass stuff to parent formatter, e.g. title and abc labeling
        BaseAxes.smart_update(self, **kwargs)

# Register the projections
mproj.register_projection(BaseAxes)
mproj.register_projection(XYAxes)
mproj.register_projection(PanelAxes)
mproj.register_projection(PolarAxes)
mproj.register_projection(BasemapAxes)
mproj.register_projection(CartopyAxes)

