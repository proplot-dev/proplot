#!/usr/bin/env python3
"""
The axes subclasses central to this library, plus the enhanced
colorbar and legend functions.

You should focus on the `format` methods: `BaseAxes.format`,
`XYAxes.format`, and `CartopyAxes.format` and `BasemapAxes.format`.
These methods are your one-stop-shop for changing axes settings like
*x* and *y* axis limits, axis labels, tick locations, tick label
format, grid lines, scales, titles, a-b-c labelling, adding
geographic features, and much more.
"""
# Note that even if not in IPython notebook, io capture output still works
# functools.wraps preserves __name__ metadata; see comment:
# https://stackoverflow.com/a/739665/4970632
import os
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
from matplotlib import docstring
from .rcmod import rc
from . import utils, colortools, fonttools, axistools
from .utils import _dot_dict, _default, _timer, _counter, ic, units
from .gridspec import FlexibleGridSpec, FlexibleGridSpecFromSubplotSpec

# Silly recursive function, returns a...z...aa...zz...aaa...zzz
# God help you if you ever need that many labels
_abc = 'abcdefghijklmnopqrstuvwxyz'
def _ascii(i, prefix=''):
    if i < 26:
        return prefix + _abc[i]
    else:
        return _ascii(i - 26, prefix) + _abc[i % 26]

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
# Aliases for axes names
# NOTE: In this package, we only create a new method 'pcolorpoly', don't override
# the 'pcolor' method, so don't run into recursion issues due to internal
# matplotlib pcolormesh() calls to pcolor()
_aliases = {
    'pcolorpoly': 'pcolor',
    'bpanel': 'bottompanel',
    'rpanel': 'rightpanel',
    'tpanel': 'toppanel',
    'lpanel': 'leftpanel'
    }

# Basic plotting tool categories, used in various places
_line_methods = (
    'plot', 'scatter', 'tripcolor', 'tricontour', 'tricontourf'
    )
_contour_methods = (
    'contour', 'tricontour',
    )
_pcolor_methods = (
    'pcolor', 'pcolormesh', 'pcolorpoly', 'tripcolor'
    )
_contourf_methods = (
    'contourf', 'tricontourf',
    )
_show_methods = (
    'imshow', 'matshow', 'spy', 'hist2d',
    )

# 2D plot functions that require coordinate centers and edges
_center_methods = (
    'contour', 'contourf', 'quiver', 'streamplot', 'barbs'
    )
_edge_methods = (
    'pcolor', 'pcolormesh', 'pcolorpoly',
    )

# Whether to wrap plot functions with cycle features or cmap features
cycle_methods  = (
    'plot', 'scatter', 'bar', 'barh', 'hist', 'boxplot', 'errorbar'
    )
"""List of plotting methods wrapped by `cycle_features`."""
cmap_methods = (
    'cmapline', 'hexbin', # special
    'contour', 'contourf', 'pcolor', 'pcolormesh',
    'matshow', 'imshow', 'spy', 'hist2d',
    'tripcolor', 'tricontour', 'tricontourf',
    )
"""List of plotting methods wrapped by `cmap_features`."""

# Finally disable some stuff for all axes, and just for map projection axes
# The keys in below dictionary are error messages
_disabled_methods = {
    "Unsupported plotting function {}. May be added soon.":
        ('pie', 'table', 'eventplot',
        'xcorr', 'acorr', 'psd', 'csd', 'magnitude_spectrum',
        'angle_spectrum', 'phase_spectrum', 'cohere', 'specgram'),
    "Redundant function {} has been disabled. Control axis scale with format(xscale='scale', yscale='scale').":
        ('semilogx', 'semilogy', 'loglog'),
    "Redundant function {} has been disabled. Date formatters will be used automatically when x/y coordinates are python datetime or numpy datetime64.":
        ('plot_date',),
    "Redundant function {} has been disabled. Use proj='polar' in subplots() call, then use angle as 'x' and radius as 'y'.":
        ('polar',)
    }
_map_disabled_methods = (
    'matshow', 'imshow', 'spy', # don't disable 'bar' or 'barh', can be used in polar plots
    # 'triplot', 'tricontour', 'tricontourf', 'tripcolor',
    'hist', 'hist2d', 'errorbar', 'boxplot', 'violinplot', 'step', 'stem',
    'hlines', 'vlines', 'axhline', 'axvline', 'axhspan', 'axvspan',
    # 'fill_between', 'fill_betweenx', 'fill',
    'stackplot'
    )

#------------------------------------------------------------------------------
# Helper functions for plot overrides
# List of stuff in pcolor/contourf that need to be fixed:
#   * White lines between the edges; cover them by changing edgecolors to 'face'.
#   * Determination of whether we are using graticule edges/centers; not sure
#       what default behavior is but harder to debug. My decorator is nicer.
#   * Pcolor can't take an extend argument, and colorbar can take an extend argument
#       but it is ignored when the mappable is a contourf. Make our pcolor decorator
#       add an "extend" attribute on the mappable that our colorbar decorator detects.
#   * Extend used in contourf causes color-change between in-range values and
#       out-of-range values, while extend used in colorbar on pcolor has no such
#       color change. Standardize by messing with the colormap.
#------------------------------------------------------------------------------
def _parse_args(args):
    """
    Parse arguments for checking 2D data centers/edges.
    """
    # Sanitize input
    if len(args)>2:
        Zs = args[2:]
    else:
        Zs = args
    Zs = [np.array(Z) for Z in Zs] # ensure array
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

def _check_centers(func):
    """
    Check shape of arguments passed to contour, and fix result.

    Notes
    -----
    Optional numbers of arguments:

    * Z
    * U, V
    * x, y, Z
    * x, y, U, V
    """
    @functools.wraps(func)
    def decorator(*args, order='C', **kwargs):
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
    return decorator

def _check_edges(func):
    """
    Check shape of arguments passed to pcolor, and fix result.
    """
    @functools.wraps(func)
    def decorator(*args, order='C', **kwargs):
        # Checks that sizes match up, checks whether graticule was input
        x, y, Zs = _parse_args(args)
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
    return decorator

def cycle_features(self, func):
    """
    Wrapper generator that allows specifying the "color cycler" at plot-time.
    This simply sets the axes property cycler before calling the plot method,
    and sets the cycle

    Parameters
    ----------
    cycle : None or cycle spec, optional
        The cycle specifer, passed to the `~proplot.colortools.Cycle`
        constructor.
    cycle_kw : dict-like, optional
        Passed to `~proplot.colortools.Cycle`.

    Other parameters
    ----------------
    *args, **kwargs
        Passed to the matplotlib plotting method.

    See also
    --------
    `cycle_methods`, `BaseAxes`, `~proplot.colortools.Cycle`

    Notes
    -----
    See the `matplotlib source 
    <https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/axes/_base.py>`_.
    The `set_prop_cycle` command modifies underlying 
    `_get_lines` and `_get_patches_for_fill`.
    """
    @functools.wraps(func)
    def decorator(*args, cycle=None, cycle_kw={}, **kwargs):
        # Determine and temporarily set cycler
        if not np.iterable(cycle) or isinstance(cycle, str):
            cycle = cycle,
        if cycle[0] is not None and not (isinstance(cycle[0], str) and cycle == rc.cycle):
            cycle = colortools.Cycle(*cycle, **cycle_kw)
            self.set_prop_cycle(color=cycle)
        return func(*args, **kwargs)
    return decorator

def cmap_features(self, func):
    """
    Wrapper generator that adds a bunch of new features, including
    flexible and on-the-fly colormap generation.

    Parameters
    ----------
    cmap : None or colormap spec, optional
        The colormap specifer, passed to the `~proplot.colortools.Colormap`
        constructor.
    cmap_kw : dict-like, optional
        Passed to `~proplot.colortools.Colormap`.
    extend : {'neither', 'both', 'min', 'max'}, optional
        Whether to assign a unique color to out-of-bounds data. Also means
        when the colorbar is drawn, colorbar "extensions" will be drawn
        (triangles by default).
    levels : None or int or list of float, optional
        If list of float, level *edges*. If integer, the number of level
        edges, and boundaries are chosen by matplotlib based on the data
        range. Defaults to 11.
    values : None or list of float, optional
        The level *centers*. If not ``None``, infers levels using
        `~proplot.utils.edges` and overrides the `levels` keyword arg.
    zero : bool, optional
        Ignored if `values` or `levels` are lists.
        If colormap levels were automatically selected, toggle this to
        modify the levels to be **symmetric about zero**.
    norm : None or normalizer spec, optional
        The colormap normalizer, used to map data values to colormap colors.
        This is passed to the `~proplot.colortools.Norm` constructor.
    norm_kw : dict-like, optional
        Passed to `~proplot.colortools.Norm`.
    values_as_levels : bool, optional
        Used internally. Toggles whether to infer `values` from `levels`, or
        to bypass them and add `values` as a keyword arg to the main function.

    * Create new colormaps on the fly, and merge arbitrary named
      or created colormaps.
    * Always use full range of colormap, whether you are extending
      max, min, neither, or both. For the first three, will reconstruct
      colormap so 'out-of-bounds' have same color as edge colors
      from 'in-bounds' region.

    Also see `this post <https://stackoverflow.com/a/48614231/4970632>`.

    Other parameters
    ----------------
    *args, **kwargs
        Passed to the matplotlib plotting method.

    See also
    --------
    `cmap_methods`, `BaseAxes`, `~proplot.colortools.Colormap`,
    `~proplot.colortools.Norm`, `~matplotlib.colors.Colormap`,
    `~matplotlib.colors.Normalizer`

    Notes
    -----
    The `bins` argument lets you choose between:

    1. (True) Use a *discrete* normalizer with a *continuous* (i.e. very
        high resolution) color table.
    2. (False) Use a *continuous* normalizer with a *discrete* (containing
        the number of colors you want) color table.
    """
    @functools.wraps(func)
    def decorator(*args, cmap=None, cmap_kw={}, extend='neither',
                values=None, levels=None, zero=False, # override levels to be *centered* on zero
                norm=None, norm_kw={},
                values_as_levels=True, # if values are passed, treat them as levels? or just use them for e.g. cmapline, then do whatever?
                **kwargs):
        # Optionally pass level centers instead of level edges
        # NOTE: See the norm_preprocessor section below for why we funnel
        # results through a normalizer here
        if kwargs.get('interp', 0): # e.g. for cmapline, we want to *interpolate*
            values_as_levels = False # get levels later down the line
        if np.iterable(values) and values_as_levels:
            # Special case of LinearSegmentedNorm, just get edges
            if isinstance(norm, colortools.LinearSegmentedNorm) or \
                (isinstance(norm, str) and 'segment' in norm):
                levels = utils.edges(values)
            # Else see what pops out
            else:
                norm_tmp = colortools.Norm(norm, **norm_kw)
                if norm_tmp: # is not None
                    levels = norm_tmp.inverse(utils.edges(norm_tmp(values)))
                else:
                    levels = utils.edges(values)
        levels = _default(levels, 11) # e.g. pcolormesh can auto-determine levels if you input a number

        # Call function with custom kwargs
        name = func.__name__
        if name in _contour_methods or name in _contourf_methods: # only valid kwargs for contouring
            kwargs.update({'levels': levels, 'extend': extend})
        if name == 'cmapline':
            kwargs.update({'values': values}) # implement this directly
        if name in _show_methods: # *do not* auto-adjust aspect ratio! messes up subplots!
            kwargs.update({'aspect': 'auto'})
        result = func(*args, **kwargs)
        if not getattr(result, 'extend', None):
            result.extend = extend # will already be on there for some funcs

        # Get levels automatically determined by contourf, or make them
        # from the automatically chosen pcolor/imshow clims.
        # TODO: This still is not respected for hexbin 'log' norm for
        # some reason, figure out fix.
        if not np.iterable(levels): # i.e. was an integer
            # Some tools automatically generate levels, like contourf.
            # Others will just automatically impose clims, like pcolor.
            # Below accounts for both options.
            if hasattr(result, 'levels'):
                levels = result.levels
            else:
                levels = np.linspace(*result.get_clim(), levels)
            # Center the levels
            # NOTE: When contourf has already rendered contours, they *cannot* be
            # changed/updated -- must be redrawn!
            if zero:
                abs_max = max([abs(max(levels)), abs(min(levels))])
                levels = np.linspace(-abs_max, abs_max, len(levels))
                if hasattr(result, 'levels'): # e.g. contourf, Artist must be re-drawn!
                    if hasattr(result, 'collections'):
                        for artist in result.collections:
                            artist.set_visible(False)
                    elif hasattr(result, 'set_visible'):
                        result.set_visible(False)
                    else:
                        raise ValueError(f'Unknown object {result}. Cannot center colormap levels.')
                    kwargs['levels'] = levels
                    result = func(*args, **kwargs)
        result.levels = levels # make sure they are on there!

        # Contour *lines* can be colormapped, but this should not be
        # default if user did not input a cmap
        if name in _contour_methods and cmap is None:
            return result

        # Get 'pre-processor' norm -- e.g. maybe user wants colormap scaled
        # in logarithmic space, or warped to diverge from center from a
        # given midpoint like zero
        norm_preprocess = colortools.Norm(norm, levels=levels, **norm_kw)
        if hasattr(result, 'norm') and norm_preprocess is None:
            norm_preprocess = result.norm

        # Create colors
        N = None # will be ignored
        norm = colortools.BinNorm(norm=norm_preprocess, levels=levels, extend=extend)
        result.set_norm(norm)

        # Specify colormap
        cmap = cmap or rc['image.cmap']
        if isinstance(cmap, (str, dict, mcolors.Colormap)):
            cmap = cmap, # make a tuple
        cmap = colortools.Colormap(*cmap, N=N, extend=extend, **cmap_kw)
        if not cmap._isinit:
            cmap._init()
        result.set_cmap(cmap)

        # Fix white lines between filled contours/mesh
        linewidth = 0.4 # seems to be lowest threshold where white lines disappear
        if name in _contourf_methods:
            for contour in result.collections:
                contour.set_edgecolor('face')
                contour.set_linewidth(linewidth)
        if name in _pcolor_methods:
            result.set_edgecolor('face')
            result.set_linewidth(linewidth) # seems to do the trick, without dots in corner being visible
        return result

    return decorator

#------------------------------------------------------------------------------#
# Helper functions for basemap and cartopy plot overrides
# NOTE: These wrappers should be invoked *after* _check_centers and _check_edges,
# which perform basic shape checking and permute the data array, so the data
# will now be y by x (or lat by lon) instead of lon by lat.
#------------------------------------------------------------------------------#
# Normally we *cannot* modify the underlying *axes* pcolormesh etc. because this
# this will cause basemap's self.m.pcolormesh etc. to use my *custom* version and
# cause a suite of weird errors. Prevent this recursion with the below decorator.
def _m_call(self, func):
    """
    Call the basemap version of the function of the same name.
    """
    name = func.__name__
    @functools.wraps(func)
    def decorator(*args, **kwargs):
        return self.m.__getattribute__(name)(ax=self, *args, **kwargs)
    return decorator

def _no_recurse(self, func):
    """
    Decorator to prevent recursion in Basemap method overrides.
    See: https://stackoverflow.com/a/37675810/4970632
    """
    @functools.wraps(func)
    # def decorator(self, *args, **kwargs):
    def decorator(*args, **kwargs):
        name = getattr(func, '__name__')
        if self._recurred:
            # Don't call func again, now we want to call the original
            # matplotlib version of this function
            self._recurred = False
            result = object.__getattribute__(self, name)(*args, **kwargs)
        else:
            # Actually return the basemap version
            self._recurred = True
            result = func(*args, **kwargs)
        self._recurred = False # cleanup, in case recursion never occurred
        return result
    return decorator

def _linefix_basemap(self, func):
    """
    Simply add an additional kwarg. Needs whole function because we
    want to @wrap it to preserve documentation.
    """
    @functools.wraps(func)
    # def decorator(self, *args, **kwargs):
    def decorator(*args, **kwargs):
        kwargs.update(latlon=True)
        return func(*args, **kwargs)
    return decorator

def _gridfix_basemap(self, func):
    """
    Interpret coordinates and fix discontinuities in grid.
    """
    @functools.wraps(func)
    def decorator(lon, lat, Z, fix_poles=True, **kwargs):
        # Raise errors
        eps = 1e-3
        lonmin, lonmax = self.m.lonmin, self.m.lonmax
        if lon.max() > lon.min() + 360 + eps:
            raise ValueError(f'Longitudes span {lon.min()} to {lon.max()}. Can only span 360 degrees at most.')
        if lon.min() < -360 or lon.max() > 360:
            raise ValueError(f'Longitudes span {lon.min()} to {lon.max()}. Must fall in range [-360, 360].')
        # if lonmin < -360 or lonmin > 0:
        #     print(f'Warning: Minimum longitude is {lonmin}, not in range [-360,0].')
        #     raise ValueError('Minimum longitude must fall in range [-360, 0].')
        # Establish 360-degree range
        lon -= 720
        while True:
            filter_ = lon<lonmin
            if filter_.sum()==0:
                break
            lon[filter_] += 360
        # Below only works with vectors
        if lon.ndim==1 and lat.ndim==1:
            # 1) Roll, accounting for whether ends are identical
            # If go from 0,1,-->,359,0 (borders), returns id of first zero
            roll = -np.argmin(lon) # always returns *first* value
            if lon[0]==lon[-1]:
                lon = np.roll(lon[:-1], roll)
                lon = np.append(lon, lon[0]+360)
            else:
                lon = np.roll(lon, roll)
            Z = np.roll(Z, roll, axis=1)
            # 2) Roll in same direction some more, if some points on right-edge
            # extend more than 360 above the minimum longitude; THEY should be the
            # ones on west/left-hand-side of map
            lonroll = np.where(lon>lonmin+360)[0] # tuple of ids
            if lonroll: # non-empty
                roll = lon.size-min(lonroll) # e.g. if 10 lons, lonmax id is 9, we want to roll once
                lon = np.roll(lon, roll) # need to roll foreward
                Z = np.roll(Z, roll, axis=1) # roll again
                lon[:roll] -= 360 # retains monotonicity
            # 3) Set NaN where data not in range lonmin, lonmax
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
            # 4) Fix holes over poles by interpolating there (equivalent to
            # simple mean of highest/lowest latitude points)
            # if self.m.projection[:4] != 'merc': # did not fix the problem where Mercator goes way too far
            if fix_poles:
                Z_south = np.repeat(Z[0,:].mean(),  Z.shape[1])[None,:]
                Z_north = np.repeat(Z[-1,:].mean(), Z.shape[1])[None,:]
                lat = np.concatenate(([-90], lat, [90]))
                Z = np.concatenate((Z_south, Z, Z_north), axis=0)
            # 5) Fix seams at map boundary; 3 scenarios here:
            # Have edges (e.g. for pcolor), and they fit perfectly against basemap seams
            # this does not augment size
            if lon[0]==lonmin and lon.size-1==Z.shape[1]: # borders fit perfectly
                pass # do nothing
            # Have edges (e.g. for pcolor), and the projection edge is in-between grid cell boundaries
            # this augments size by 1
            elif lon.size-1==Z.shape[1]: # no interpolation necessary; just make a new grid cell
                lon = np.append(lonmin, lon) # append way easier than concatenate
                lon[-1] = lonmin + 360 # we've added a new tiny cell to the end
                Z = np.concatenate((Z[:,-1:], Z), axis=1) # don't use pad; it messes up masked arrays
            # Have centers (e.g. for contourf), and we need to interpolate to the
            # left/right edges of the map boundary
            # this augments size by 2
            elif lon.size==Z.shape[1]: # linearly interpolate to the edges
                x = np.array([lon[-1], lon[0]+360]) # x
                if x[0] != x[1]:
                    y = np.concatenate((Z[:,-1:], Z[:,:1]), axis=1)
                    xq = lonmin+360
                    yq = (y[:,:1]*(x[1]-xq) + y[:,1:]*(xq-x[0]))/(x[1]-x[0]) # simple linear interp formula
                    Z = np.concatenate((yq, Z, yq), axis=1)
                    lon = np.append(np.append(lonmin, lon), lonmin+360)
            else:
                raise ValueError()
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
    return decorator

def _linefix_cartopy(self, func):
    """
    Simply add an additional kwarg. Needs whole function because we
    want to @wrap it to preserve documentation.
    """
    @functools.wraps(func)
    def decorator(*args, transform=PlateCarree, **kwargs):
        # Simple
        if isinstance(transform, type):
            transform = transform() # instantiate
        result = func(*args, transform=transform, **kwargs)
        # Some plot functions seem to reset the outlinepatch or
        # backgroundpatch (???), so need to re-enforce settings.
        self._rcupdate()
        return result
    return decorator

def _gridfix_cartopy(self, func):
    """
    Apply default transform and fix discontinuities in grid.
    Note for cartopy, we don't have to worry about meridian at which longitude
    wraps around; projection handles all that.

    Todo
    ----
    Contouring methods for some reason have issues with circularly wrapped
    data. Triggers annoying TopologyException statements, which we suppress
    with IPython capture_output() tool, like in nbsetup().
    See: https://github.com/SciTools/cartopy/issues/946
    """
    @functools.wraps(func)
    def decorator(lon, lat, Z, transform=PlateCarree, fix_poles=True, **kwargs):
        # Below only works for vector data
        if lon.ndim==1 and lat.ndim==1:
            # 1) Fix holes over poles by *interpolating* there (equivalent to
            # simple mean of highest/lowest latitude points)
            if fix_poles:
                Z_south = np.repeat(Z[0,:].mean(),  Z.shape[1])[None,:]
                Z_north = np.repeat(Z[-1,:].mean(), Z.shape[1])[None,:]
                lat = np.concatenate(([-90], lat, [90]))
                Z = np.concatenate((Z_south, Z, Z_north), axis=0)
            # 2) Fix seams at map boundary; by ensuring circular coverage
            if (lon[0] % 360) != ((lon[-1] + 360) % 360):
                lon = np.array((*lon, lon[0] + 360)) # make longitudes circular
                Z = np.concatenate((Z, Z[:,:1]), axis=1) # make data circular
        # Instantiate transform
        if isinstance(transform, type):
            transform = transform() # instantiate
        # Call function
        with io.capture_output() as captured:
            result = func(lon, lat, Z, transform=transform, **kwargs)
        # Some plot functions seem to reset the outlinepatch or
        # backgroundpatch (???), so need to re-enforce settings.
        self._rcupdate()
        return result
    return decorator

#------------------------------------------------------------------------------#
# Generalized custom axes class
#------------------------------------------------------------------------------#
class BaseAxes(maxes.Axes):
    """
    Lowest-level `~matplotlib.axes.Axes` override. Handles titles and axis
    sharing. Overrides the legend, colorbar, and plot methods.

    Notes
    -----
    It is impossible to subclass `~matplotlib.axes.SubplotBase` directly.
    The subplot subclass is made automatically by
    `~matplotlib.axes.subplot_class_factory`
    after calling `~matplotlib.figure.Figure.add_subplot`.

    See also
    --------
    `~proplot.subplots.subplots`, `XYAxes`, `PanelAxes`,
    `MapAxes`, `CartopyAxes`, `BasemapAxes`

    Todo
    ----
    Consider moving the share axis stuff to `XYAxes` class.
    """
    # Initial stuff
    name = 'base'
    """The registered projection name."""
    def __init__(self, *args, number=None,
            sharex=None, sharey=None, spanx=None, spany=None,
            sharex_level=0, sharey_level=0,
            map_name=None,
            panel_parent=None, panel_side=None,
            **kwargs):
        # Copied sharex stuff from subplots documentation
        """
        Parameters
        ----------
        number : None or int
            The subplot number, used for a-b-c labelling (see
            `~BaseAxes.format`).
        sharex_level, sharey_level : {3, 2, 1, 0}, optional
            The "axis sharing level" for the *x* axis, *y* axis, or both
            axes.

            See `~proplot.subplots.subplots` for details.
        sharex, sharey : None or `BaseAxes`, optional
            Axes to use for *x* and *y* axis sharing. Should correspond
            to the subplot in the bottommost row, leftmost column.
        spanx, spany : None or `BaseAxes`, optional
            Axes to use for the "spanning" *x* and *y* axis labels. Should
            correspond to the subplot in the leftmost column, bottommost row.

            See `~proplot.subplots.subplots` for details.
        panel_parent : None or `BaseAxes`, optional
            The parent of the panel. Used only with `PanelAxes`.
        panel_side : {None, 'left', 'right', 'top', 'bottom'}, optional
            The side the panel was drawn on. Used only with `PanelAxes`.
        map_name : None or str, optional
            The name of the map projection. Used only with `BasemapAxes`.
        """
        # Initialize
        self._spanx = spanx # boolean toggles, whether we want to span axes labels
        self._spany = spany
        self._title_inside = True # toggle this to figure out whether we need to push 'super title' up
        self._zoom = None # if non-empty, will make invisible
        self._inset_parent = None # change this later
        self._insets = [] # add to these later
        self._map_name = map_name # consider conditionally allowing 'shared axes' for certain projections
        self._gridliner_on = False # whether cartopy gridliners are plotted here; matplotlib tight bounding box does not detect them! so disable smart_tight_layout if True
        super().__init__(*args, **kwargs)
        self._title_pos_init = self.title.get_position() # position of title outside axes
        self._title_pos_transform = self.title.get_transform()

        # Panels
        if panel_side not in (None, 'left','right','bottom','top'):
            raise ValueError(f'Invalid panel side "{panel_side}".')
        self.panel_side = panel_side
        self.panel_parent = panel_parent # used when declaring parent
        self.bottompanel = EmptyPanel()
        self.toppanel    = EmptyPanel()
        self.leftpanel   = EmptyPanel()
        self.rightpanel  = EmptyPanel()

        # Number and size
        if isinstance(self, maxes.SubplotBase):
            nrows, ncols, subspec = self._topmost_subspec()
            self._row_span = ((subspec.num1 // ncols) // 2, (subspec.num2 // ncols) // 2)
            self._col_span = ((subspec.num1 % ncols) // 2,  (subspec.num2 % ncols) // 2)
        else:
            self._row_span = None
            self._col_span = None
        self.number = number # for abc numbering
        self.width  = np.diff(self._position.intervalx)*self.figure.width # position is in figure units
        self.height = np.diff(self._position.intervaly)*self.figure.height

        # Turn off tick labels and axis label for shared axes
        # Want to do this ***manually*** because want to have the ability to
        # add shared axes ***after the fact in general***. If the API changes,
        # will modify the below methods.
        self._sharex_setup(sharex, sharex_level)
        self._sharey_setup(sharey, sharey_level)

        # Add extra text properties for abc labeling, rows/columns labels
        # (can only be filled with text if axes is on leftmost column/topmost row)
        self.abc = self.text(0, 0, '') # position tbd
        self.collabel = self.text(*self.title.get_position(), '',
                va='baseline', ha='center', transform=self.title.get_transform())
        self.rowlabel = self.text(*self.yaxis.label.get_position(), '',
                va='center', ha='right', transform=self.transAxes)

        # Enforce custom rc settings! And only look for rcSpecial settings,
        # because the 'global' settings will be applied on initialization/do
        # not have to be explicitly re-applied.
        with rc._context(mode=1):
            self._rcupdate()

    # Apply some simple featueres, and disable spectral and triangular features
    # See: https://stackoverflow.com/a/23126260/4970632
    # Also see: https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/axes/_axes.py
    # for all Axes methods ordered logically in class declaration.
    def __getattribute__(self, attr, *args):
        for message,attrs in _disabled_methods.items():
            if attr in attrs:
                raise NotImplementedError(message.format(attr))
        attr = _aliases.get(attr, attr)
        obj = super().__getattribute__(attr, *args)
        if attr in cmap_methods:
            obj = cmap_features(self, obj)
        elif attr in cycle_methods:
            obj = cycle_features(self, obj)
        return obj

    def _topmost_subspec(self):
        # Needed for e.g. getting the top-level SubplotSpec (i.e. the one
        # encompassed by an axes and all its panels, if any are present)
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

    def _sharex_setup(self, sharex, level):
        if sharex is None:
            return
        if self is sharex:
            return
        if isinstance(self, MapAxes) or isinstance(sharex, MapAxes):
            return
        if level not in range(4):
            raise ValueError('Level can be 1 (do not share limits, just hide axis labels), 2 (share limits, but do not hide tick labels), or 3 (share limits and hide tick labels).')
        # Share vertical panel x-axes with *eachother*
        if self.leftpanel and sharex.leftpanel:
            self.leftpanel._sharex_setup(sharex.leftpanel, level)
        if self.rightpanel and sharex.rightpanel:
            self.rightpanel._sharex_setup(sharex.rightpanel, level)
        # Share horizontal panel x-axes with *sharex*
        if self.bottompanel and sharex is not self.bottompanel:
            self.bottompanel._sharex_setup(sharex, level)
        if self.toppanel and sharex is not self.toppanel:
            self.toppanel._sharex_setup(sharex, level)
        # Builtin features
        self._sharex = sharex
        if level>1:
            self._shared_x_axes.join(self, sharex)
        # Simple method for setting up shared axes
        # WARNING: It turned out setting *another axes' axis label* as
        # this attribute caused error, because matplotlib tried to add
        # the same artist instance twice. Can only make it invisible.
        if level>2:
            for t in self.xaxis.get_ticklabels():
                t.set_visible(False)
        self.xaxis.label.set_visible(False)

    def _sharey_setup(self, sharey, level):
        if sharey is None:
            return
        if self is sharey:
            return
        if isinstance(self, MapAxes) or isinstance(sharey, MapAxes):
            return
        if level not in range(4):
            raise ValueError('Level can be 1 (do not share limits, just hide axis labels), 2 (share limits, but do not hide tick labels), or 3 (share limits and hide tick labels).')
        # Share horizontal panel y-axes with *eachother*
        if self.bottompanel and sharey.bottompanel:
            self.bottompanel._sharey_setup(sharey.bottompanel, level)
        if self.toppanel and sharey.toppanel:
            self.toppanel._sharey_setup(sharey.toppanel, level)
        # Share vertical panel y-axes with *sharey*
        if self.leftpanel:
            self.leftpanel._sharey_setup(sharey, level)
        if self.rightpanel:
            self.rightpanel._sharey_setup(sharey, level)
            # sharey = self.leftpanel._sharey or self.leftpanel
        # Builtin features
        self._sharey = sharey
        if level>1:
            self._shared_y_axes.join(self, sharey)
        # Simple method for setting up shared axes
        if level>2:
            for t in self.yaxis.get_ticklabels():
                t.set_visible(False)
        self.yaxis.label.set_visible(False)

    def _sharex_panels(self):
        # Call this once panels are all declared
        if self.bottompanel:
            self._sharex_setup(self.bottompanel, 3)
        bottom = self.bottompanel or self
        if self.toppanel:
            self.toppanel._sharex_setup(bottom, 3)

    def _sharey_panels(self):
        # Same but for y
        if self.leftpanel:
            self._sharey_setup(self.leftpanel, 3)
        left = self.leftpanel or self
        if self.rightpanel:
            self.rightpanel._sharey_setup(left, 3)

    def _rcupdate(self):
        # Figure patch (for some reason needs to be re-asserted even if declared before figure drawn)
        kw = rc.fill({'facecolor':'figure.facecolor'})
        self.figure.patch.update(kw)

        # Axes, figure title (builtin settings)
        kw = rc.fill({'fontsize':'axes.titlesize', 'weight':'axes.titleweight', 'fontname':'fontname'})
        self.title.update(kw)
        kw = rc.fill({'fontsize':'figure.titlesize', 'weight':'figure.titleweight', 'fontname':'fontname'})
        self.figure._suptitle.update(kw)

        # Row and column labels, ABC labels
        kw = rc.fill({'fontsize':'abc.fontsize', 'weight':'abc.weight', 'color':'abc.color', 'fontname':'fontname'})
        self.abc.update(kw)
        kw = rc.fill({'fontsize':'rowlabel.fontsize', 'weight':'rowlabel.weight', 'color':'rowlabel.color', 'fontname':'fontname'})
        self.rowlabel.update(kw)
        kw = rc.fill({'fontsize':'collabel.fontsize', 'weight':'collabel.weight', 'color':'collabel.color', 'fontname':'fontname'})
        self.collabel.update(kw)

        # Hatching options (useful where we want to highlight invalid data)
        # NOTE: So that we can keep re-accessing hatches between multiple calls,
        # we will add hatches to the patch object directly.
        # TODO: This only works for XYAxes! Need to do something else for
        # PolarAxes! Or just add a _bg_hatch property or something?
        # self.patch.update(kw)
        kw = rc.fill({'hatch':'facehatch', 'edgecolor':'axes.hatchcolor', 'alpha':'hatchalpha'})
        if kw and kw.get('hatch',None): # non-empty
            self.fill_between([0,1], 0, 1, zorder=0, # put in back
                facecolor='none', transform=self.transAxes, **kw)


    def _text_update(self, obj, kwargs):
        """
        Update text, and allow for border.
        """
        # Allow updating properties introduced by the BaseAxes.text() override.
        # Don't really want to subclass mtext.Text; only have a few features
        # NOTE: Don't use kwargs because want this to look like standard
        # artist self.update.
        try:
            obj.update(kwargs)
        except Exception:
            obj.set_visible(False)
            text     = kwargs.pop('text', obj.get_text())
            color    = kwargs.pop('color', obj.get_color())
            weight   = kwargs.pop('weight', obj.get_weight())
            fontsize = kwargs.pop('fontsize', obj.get_fontsize())
            x, y = None, None
            pos = obj.get_position()
            if 'position' in kwargs:
                x, y = kwargs.pop('position')
            x = kwargs.pop('x', x)
            y = kwargs.pop('y', y)
            if x is None:
                x = pos[0]
            if y is None:
                y = pos[1]
            try:
                x = x[0]
            except TypeError:
                pass
            try:
                y = y[0]
            except TypeError:
                pass
            obj = self.text(x, y, text, color=color, weight=weight, fontsize=fontsize, **kwargs)
        return obj

    def _title_pos(self, pos, **kwargs):
        """
        Position title text to the left, center, or right and either inside or
        outside the axes (default is center, outside).
        """
        ypad = (rc['axes.titlepad']/72)/self.height # to inches --> to axes relative
        xpad = (rc['axes.titlepad']/72)/self.width # why not use the same for x?
        xpad_i, ypad_i = xpad*1.5, ypad*1.5 # inside labels need a bit more room
        pos = pos or 'oc'
        extra = {}
        if not isinstance(pos, str):
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
            # Record _title_inside so we can automatically deflect suptitle
            # If *any* object is outside (title or abc), want to deflect it up
            if 'o' in pos:
                y = 1 # 1 + ypad # leave it alone, may be adjusted during draw-time to account for axis label (fails to adjust for tick labels; see notebook)
                va = 'bottom'
                self._title_inside = False
                transform = self.title.get_transform()
            elif 'i' in pos:
                y = 1 - ypad_i
                va = 'top'
                self._title_inside = True
                transform = self.transAxes
                extra['border'] = _default(kwargs.pop('border', None), True) # by default
        return {'x':x, 'y':y, 'transform':transform, 'ha':ha, 'va':va, **extra}

    def invisible(self):
        """
        Make axes invisible.
        """
        for s in self.spines.values():
            s.set_visible(False)
        self.xaxis.set_visible(False)
        self.yaxis.set_visible(False)
        self.patch.set_alpha(0)

    # New convenience feature
    # The title position can be a mix of 'l/c/r' and 'i/o'
    def format(self,
        suptitle=None, suptitle_kw={},
        collabels=None, collabels_kw={},
        rowlabels=None, rowlabels_kw={}, # label rows and columns
        title=None, titlepos=None, title_kw={},
        abc=None, abcpos=None, abcformat=None, abc_kw={},
        rc_kw={}, **kwargs,
        ):
        """
        Function for formatting axes titles and labels.

        Parameters
        ----------
        suptitle : None or str, optional
            The figure "super" title, centered between the left edge of
            the lefmost column of subplots and the right edge of the rightmost
            column of subplots, and automatically offset above figure titles.

            This is more sophisticated than matplotlib's builtin "super title",
            which is just centered between the figure edges and offset from
            the top edge.
        suptitle_kw : dict-like, optional
            The super title settings. Applied with the
            `~matplotlib.artist.Artist.update` method on the
            `~matplotlib.text.Text` instance. Options include ``'color'``,
            ``'size'``, and ``'weight'``.
        title : None or str, optional
            The axes title.
        titlepos : None or str, optional
            The position of the axes title, set by a string consisting of up
            to 2 characters. The default is ``'co'``.

            There are 2 options for the vertical position:

            * ``'o'``: Outside the axes.
            * ``'i'``: Inside the axes. In this case, a white border
              will be drawn around the title text by default.
              Disable with ``title_kw={'border':False}``.

            There are 3 options for the horizontal position:

            * ``'c'``: Centered, the default.
            * ``'l'``: Left aligned.
            * ``'r'``: Right aligned.

        title_kw : dict-like, optional
            The title settings. See `suptitle_kw`.
        rowlabels, colllabels : None or list of str, optional
            The subplot row and column labels. If list, length must match
            the number of subplot rows, columns.
        rowlabels_kw, collabels_kw : dict-like, optional
            The row and column label settings. See `suptitle_kw`.
        abc : None or bool, optional
            Whether to apply "a-b-c" subplot labelling based on the
            `number` attribute.

            If `number` is >26, labels will loop around to a, ..., z, aa,
            ..., zz, aaa, ..., zzz, ... God help you if you ever need that
            many labels.
        abcpos : None or str, optional
            The position of the subplot label. See `titlepos`.
            The default is ``'li'``.
        abcformat : None or str, optional
            A string containing the character ``'a'``, specifying the
            format of the a-b-c labelling. ``'a.'`` is the default, but
            ``'(a)'`` or ``'a)'`` might also be desirable.
        abc_kw : dict-like, optional
            The subplot label settings. See `suptitle_kw`.
        rc_kw, kwargs
            "rc" configuration settings specific to this axes, used to
            temporarily update the `~proplot.rcmod.rc` object.

            See `~proplot.rcmod` for details.

        Todo
        ----
        * Add options for datetime handling. Note possible date axes handles
          are `pandas.TimeStamp`, `numpy.datetime64`, and `pandas.DateTimeIndex`.
        * Can fix with `~matplotlib.figure.Figure.autofmt_xdate` or manually
          set options. The problem is there is no `~matplotib.figure.Figure.autofmt_ydate`,
          so really should implement my own version of this.
        * Look into `~matplotlib.axes.SubplotBase.is_last_row` and
          `~matplotlib.axes.SubplotBase.is_first_column` methods.

        Notes
        -----
        `pandas.TimeStamp`, `numpy.datetime64`, and `pandas.DateTimeIndex`.

        See also
        --------
        `~proplot.subplots.subplots`, `~proplot.rcmod`,
        `XYAxes.format`, `BasemapAxes.format`, `CartopyAxes.format`
        """
        # NOTE: These next two are actually *figure-wide* settings, but that
        # line seems to get blurred -- where we have shared axes, spanning
        # labels, and whatnot. May result in redundant assignments if formatting
        # more than one axes, but operations are fast so some redundancy is nbd.
        # Create figure title
        fig = self.figure # the figure
        if suptitle is not None:
            fig._suptitle_setup(text=suptitle, auto=False, **suptitle_kw)
        if rowlabels is not None:
            fig._rowlabels(rowlabels, **rowlabels_kw)
        if collabels is not None:
            fig._collabels(collabels, **collabels_kw)

        # Create axes title
        # Input needs to be emptys string
        if title is not None:
            pos_kw = self._title_pos(titlepos or 'oc', **title_kw)
            self.title = self._text_update(self.title, {'text':title, 'visible':True, **pos_kw, **title_kw})

        # Create axes numbering
        if self.number is not None and abc:
            # Get text
            abcformat = abcformat or 'a'
            if 'a' not in abcformat:
                raise ValueError(f'Invalid abcformat {abcformat}.')
            abcedges = abcformat.split('a')
            text = abcedges[0] + _ascii(self.number-1) + abcedges[-1]
            pos_kw = self._title_pos(abcpos or 'il')
            self.abc = self._text_update(self.abc, {'text':text, **abc_kw, **pos_kw})
        elif hasattr(self, 'abc') and abc is not None and not abc:
            # Hide
            self.abc.set_visible(False)

        # First update (note that this will call _rcupdate overridden by child
        # classes, which can in turn call the parent class version, so we only
        # need to call this from the base class, and all settings will be applied)
        with rc._context(rc_kw, mode=2, **kwargs):
            self._rcupdate()

    # Create legend creation method
    def legend(self, *args, **kwargs):
        """
        Add legend.
        """
        # Call custom legend() function.
        return legend_factory(self, *args, **kwargs)

    # Fill entire axes with colorbar
    # TODO: Make the default behavior draw a tiny colorbar inset
    # in the axes itself with InsetAxes! Then panel axes overrides this
    # default behavior, so that panel.colorbar *fills* the axes with a colorbar.
    def colorbar(self, *args, **kwargs):
        """
        Add *inset* colorbar, sort of like a legend.
        """
        raise NotImplementedError('Inset colorbars for non-panel axes is not yet implemented.')
        # WARNING: Below will have bugs because 'self' is a BaseAxes! See comment
        # under colorbar method on PanelAxes class.
        # return colorbar_factory(self, *args, **kwargs)

    # Fancy wrappers
    def text(self, x, y, text,
            transform='data',
            border=False, border_kw={},
            invert=False,
            linewidth=2, lw=None, **kwargs): # linewidth is for the border
        """
        Wrapper around original text method. Mainly adds feature for drawing
        "borders" around text.

        Parameters
        ----------
        x, y : float
            The *x* and *y* coordinates for the text.
        text : str
            The text.
        transform : {'data', 'axes', 'figure'} or `~matplotlib.transforms.Transform`, optional
            The coordinate system (i.e. "transform"), or a string pointing
            to a transform that we will look up. Default is ``'data'``.
        border : bool, optional
            Whether to draw border around text.
        border_kw : dict-like, optional
            Passed to `~matplotlib.path_effects.Stroke` if drawing a border.
        invert : bool, optional
            Ignored if `border` is ``False``. Whether to draw black text
            with a white border (``False``), or white text on a black
            border (``True``).
        linewidth : float, optional
            Ignored if `border` is ``False``. The width of the text border.
        lw
            Alias for `linewidth`.

        Other parameters
        ----------------
        **kwargs
            Passed to `~matplotlib.text.Text` instantiator.

        Warning
        -------
        Basemap gridlining methods call text, so if you change the default
        transform, you will not be able to draw latitude and longitude
        labels! Leave it alone.
        """
        # Get default transform by string name
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
            print(f'Warning: Font "{name}" unavailable, falling back to DejaVu Sans.' + suffix)
            fonttools._missing_fonts.append(name)
            name = 'DejaVu Sans'
        # Call parent, with custom rc settings
        # These seem to sometimes not get used by default
        size   = kwargs.pop('fontsize', rc['font.size'])
        color  = kwargs.pop('color', rc['text.color'])
        weight = kwargs.pop('font', rc['font.weight'])
        t = super().text(x, y, text, transform=transform, fontname=name,
            fontsize=size, color=color, fontweight=weight, **kwargs)
        # Optionally draw border around text
        if border:
            facecolor, bgcolor = color, 'w'
            if invert:
                facecolor, bgcolor = bgcolor, facecolor
            # facecolor, bgcolor = ('wk' if invert else 'kw')
            kwargs = {'linewidth':linewidth, 'foreground':bgcolor, 'joinstyle':'miter'}
            kwargs.update(border_kw)
            t.update({'color':facecolor, 'zorder':1e10, # have to update after-the-fact for path effects
                'path_effects': [mpatheffects.Stroke(**kwargs), mpatheffects.Normal()]})
        return t

    def plot(self, *args, cmap=None, values=None, **kwargs):
        """
        As in `~matplotlib.axes.Axes.plot`, but adds functionality
        for making `~matplotlib.collections.LineCollection` lines. These
        are lines whose colors change as a function of the coordinates
        `values`.

        Parameters
        ----------
        *args
            Passed to `~matplotlib.axes.Axes.plot`.
        cmap, values
            Passed to `~BaseAxes.cmapline`.
        **kwargs
            `~matplotlib.lines.Line2D` properties.

        """
        if cmap is None and values is None:
            # Make normal boring lines
            lines = super().plot(*args, **kwargs)
        elif cmap is not None and values is not None:
            # Make special colormap lines
            lines = self.cmapline(*args, cmap=cmap, values=values, **kwargs)
        else:
            # Error
            raise ValueError('To draw colormap line, must provide kwargs "values" and "cmap".')
        return lines

    def scatter(self, *args, **kwargs):
        """
        Just makes keyword arg conventions consistent with `plot`. This is
        something that always bothered me.
        """
        # Manage input arguments
        if len(args)>4:
            raise ValueError(f'Function accepts up to 4 args, received {len(args)}.')
        args = [*args]
        if len(args)>3:
            kwargs['c'] = args.pop(3)
        if len(args)>2:
            kwargs['s'] = args.pop(2)
        # Apply some aliases for keyword arguments
        aliases = {
            'c':          ['color', 'markercolor'],
            's':          ['size', 'markersize'],
            'linewidths': ['lw', 'linewidth', 'markeredgewidth', 'markeredgewidths'],
            'edgecolors': ['markeredgecolor', 'markeredgecolors'],
            }
        for name,options in aliases.items():
            for option in options:
                if option in kwargs:
                    kwargs[name] = kwargs.pop(option)
        return super().scatter(*args, **kwargs)

    def cmapline(self, *args, values=None,
            cmap=None, norm=None,
            interp=0, **kwargs):
        """
        Create lines with colormap. See `this matplotlib example
        <https://matplotlib.org/gallery/lines_bars_and_markers/multicolored_line.html>`_.
        Calling `BaseAxes.plot` with kwargs `cmap` and `values` will
        invoke this method.

        Parameters
        ----------
        values : list of float
            Values corresponding to points on the line.
        cmap : None or colormap spec, optional
            Colormap specifier, passed to `~proplot.colortools.colormap`.
        norm : None or `~matplotlib.colors.Normalizer`, optional
            The normalizer, used for mapping `values` to colormap colors.
        interp : int, optional
            Number of values between each line joint and each *halfway* point
            between line joints to which you want to interpolate.

        Warning
        -------
        So far this only works for **1D** *x* and *y* coordinates. Cannot draw
        multiple colormap lines at once, unlike `~matplotlib.axes.Axes.plot`.
        """
        # First error check
        if values is None:
            raise ValueError('For line with a "colormap", must input values=<iterable> to which colors will be mapped.')
        if len(args) not in (1,2):
            raise ValueError(f'Function requires 1-2 arguments, got {len(args)}.')
        y = np.array(args[-1]).squeeze()
        x = np.arange(y.shape[-1]) if len(args)==1 else np.array(args[0]).squeeze()
        values = np.array(values).squeeze()
        if x.ndim!=1 or y.ndim!=1 or values.ndim!=1:
            raise ValueError(f'Input x ({x.ndim}-d), y ({y.ndim}-d), and values ({values.ndim}-d) must be 1-dimensional.')
        if len(x)!=len(y) or len(x)!=len(values) or len(y)!=len(values):
            raise ValueError(f'Got {len(x)} xs, {len(y)} ys, but {len(values)} colormap values.')

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
        collection = mcollections.LineCollection(np.array(coords), cmap=cmap, norm=norm, linestyles='-', joinstyle='miter')
        collection.set_array(np.array(values))
        collection.update({key:value for key,value in kwargs.items() if key not in ('color',)})

        # Add collection, with some custom attributes
        self.add_collection(collection)
        collection.values = values
        collection.levels = levels # needed for other functions some
        return collection

#------------------------------------------------------------------------------#
# Specific classes, which subclass the base one
#------------------------------------------------------------------------------#
class XYAxes(BaseAxes):
    """
    Subclass for ordinary Cartesian axes.

    See also
    --------
    `~proplot.subplots.subplots`, `BaseAxes`, `PanelAxes`
    """
    # Initialize
    name = 'xy'
    """The registered projection name."""
    def __init__(self, *args, **kwargs):
        # Create simple x by y subplot.
        super().__init__(*args, **kwargs)
        # Change the default formatter
        # My trims trailing zeros, but numbers no longer aligned. Matter
        # of taste really; will see if others like it.
        formatter = axistools.Formatter('default')
        self.xaxis.set_major_formatter(formatter)
        formatter = axistools.Formatter('default')
        self.yaxis.set_major_formatter(formatter)
        # Reset this; otherwise matplotlib won't automatically change
        # formatter when it encounters certain types of data, like
        # datetime.
        self.xaxis.isDefault_majfmt = True
        self.yaxis.isDefault_majfmt = True

    def __getattribute__(self, attr, *args):
        # Attribute
        obj = super().__getattribute__(attr, *args)
        if attr in _center_methods:
            obj = _check_centers(obj)
        elif attr in _edge_methods:
            obj = _check_edges(obj)
        return obj

    def _share_span_label(self, axis):
        # Bail
        # TODO: This works, but for spanning y-axes there is danger that we
        # pick the y-label where y ticks are really narrow and spanning
        # label crashes into tick labels on another axes.
        name = axis.axis_name
        base = self
        base = getattr(base, '_share' + name, None) or base
        base = getattr(base, '_share' + name, None) or base
        if not getattr(base, '_span'  + name):
            return getattr(base, name + 'axis').label
        # Get the 'edge' we want to share (bottom row, or leftmost column),
        # and then finding the coordinates for the spanning axes along that edge
        axs = []
        span = lambda ax: getattr(ax, '_col_span') if name=='x' else getattr(ax, '_row_span')
        edge = lambda ax: getattr(ax, '_row_span')[1] if name=='x' else getattr(ax, '_col_span')[0]
        # Identify the *main* axes spanning this edge, and if those axes have
        # a panel and are shared with it (i.e. has a _sharex/_sharey attribute
        # declared with _sharex_panels), point to the panel label
        axs = [ax for ax in self.figure.axes if isinstance(ax, BaseAxes)
            and not isinstance(ax, PanelAxes) and edge(ax)==edge(base)]
        span_all = np.array([span(ax) for ax in axs])

        # Build the transform object
        idx = slice(span_all.min(), span_all.max() + 1)
        if name=='x': # span columns
            subspec = self.figure._gridspec[0,idx]
        else: # spans rows
            subspec = self.figure._gridspec[idx,0]
        bbox = subspec.get_position(self.figure) # in figure-relative coordinates
        x0, y0, width, height = bbox.bounds
        if name=='x':
            transform = mtransforms.blended_transform_factory(self.figure.transFigure, mtransforms.IdentityTransform())
            position = (x0 + width/2, 1)
        else:
            transform = mtransforms.blended_transform_factory(mtransforms.IdentityTransform(), self.figure.transFigure)
            position = (1, y0 + height/2)

        # Update the label we selected
        if axis not in self.figure._span_labels:
            self.figure._span_labels.append(axis)
        ax = axs[np.argmin(span_all[:,0])]
        if name=='x':
            axis = (ax._sharex or ax).xaxis
        else:
            axis = (ax._sharey or ax).yaxis
        axis.label.update({'visible':True, 'position':position, 'transform':transform})
        return axis.label

    def _rcupdate(self):
        # Axis settings
        for name,axis in zip('xy', (self.xaxis, self.yaxis)):
            # Optionally apply an x/y axis specific color
            label_side = axis.get_label_position()
            axis_color = rc[name + 'color'] # the special 'xcolor' and 'ycolor' property, for changing color on one spine
            sides = ('bottom','top') if name=='x' else ('left','right')

            # Update the rcParams according to user input.
            # The ticks/spines can be on both sides or just one, while the
            # tick labels and axis label on just one side
            kw = rc.fill({'lw':'axes.linewidth', 'color':'axes.edgecolor'})
            if axis_color:
                kw['color'] = axis_color
            for side in sides:
                self.spines[side].update(kw)

            # Tick marks
            # NOTE: We decide that tick location should be controlled only
            # by format(), so don't override that here.
            kw = rc.fill({
                'color': name + 'tick.color',
                'labelcolor': 'axes.edgecolor',
                'labelsize':  'axes.labelsize'
                })
            if axis_color:
                kw['color'] = axis_color
            axis.set_tick_params(which='both', **kw)
            for which in ('major','minor'):
                axis.set_tick_params(which=which, **rc[name + 'tick.' + which])

            # Tick labels
            # TODO: Figure out how to change fontname for all ticks, like in
            # set_tick_params. But fontname is global.
            kw = rc.fill({'fontname': 'fontname'})
            if axis_color:
                kw['color'] = axis_color
            for t in axis.get_ticklabels():
                t.update(kw)

            # Axis label
            kw = rc.fill({
                'color':    'axes.edgecolor',
                'fontname': 'fontname',
                'fontsize': 'axes.labelsize',
                'weight':   'axes.labelweight'})
            if axis_color:
                kw['color'] = axis_color
            axis.label.update(kw)

            # Manually update gridlines
            for grid,ticks in zip(['grid','gridminor'],[axis.get_major_ticks(), axis.get_minor_ticks()]):
                kw = rc[grid]
                for tick in ticks:
                    tick.gridline.update(kw)

        # Background patch basics
        self.patch.set_clip_on(False)
        self.patch.set_zorder(-1)
        self.patch.update(rc.fill({'facecolor': 'axes.facecolor'}))

        # Call parent
        super()._rcupdate()

    # Cool overrides
    def format(self,
        xloc=None, yloc=None, # aliases for 'where to put spine'
        xspineloc=None,  yspineloc=None,  # deals with spine options
        xtickloc=None,   ytickloc=None,   # which spines to draw ticks on
        xlabelloc=None,  ylabelloc=None,
        xticklabelloc=None, yticklabelloc=None, # where to put tick labels
        xtickdir=None,  ytickdir=None,   tickdir=None,  # which direction ('in', 'out', or 'inout')
        tickminor=None, xtickminor=True, ytickminor=True, # minor ticks on/off
        grid=None,      xgrid=None,      ygrid=None,      # gridline toggle
        gridminor=None, xgridminor=None, ygridminor=None, # minor grids on/off (if ticks off, grid will always be off)
        xticklabeldir=None, yticklabeldir=None, ticklabeldir=None, # which direction to draw labels
        xtickrange=None,    ytickrange=None,    # limit regions where we assign ticklabels to major-ticks
        xreverse=False, yreverse=False, # special properties
        xlabel=None,    ylabel=None,    # axis labels
        xlim=None,      ylim=None,
        xbounds=None,   ybounds=None, # limit spine bounds?
        xscale=None,    yscale=None,
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
            Whether to reverse the *x* and *y* axis limits.
        xbounds, ybounds : None or length-2 list of float, optional
            The *x* and *y* axis data bounds within which to draw the spines.
            For example, this can be used to separate the *x*
            and *y* axis spines, so they don't meet at the corner.
        xtickrange, ytickrange : None or length-2 list of float, optional
            The *x* and *y* axis data ranges within which major tick marks
            are labelled. For example, the tick range ``[-1,1]`` with
            axis range ``[-5,5]`` and a tick interval of 1 will only
            label the ticks marks at -1, 0, and 1.
        xscale, yscale : None or scale spec, optional
            The *x* and *y* axis scales. Arguments are passed to and interpreted
            by `~proplot.axistools.Scale`. Examples include ``'linear'``
            and ``('cutoff', 0.5, 2)``.
        xscale_kw, yscale_kw : dict-like, optional
            The *x* and *y* axis scale settings. Passed to
            `~proplot.axistools.Scale`.
        xspineloc, yspineloc : {'both', 'bottom', 'top', 'left', 'right', 'neither', 'zero'}, optional
            The *x* and *y* axis spine locations.
        xloc, yloc
            Aliases for `xspineloc`, `yspineloc`.
        xtickloc, ytickloc : {'both', 'bottom', 'top', 'left', 'right', 'neither'}, optional
            Which *x* and *y* axis spines should have major and minor tick marks.
        xtickminor, ytickminor, tickminor : None or bool, optional
            Whether to draw minor ticks on the *x* axis, *y* axis, and
            both axes.
        xtickdir, ytickdir, tickdir : {'out', 'in', 'inout'}
            Direction that major and minor tick marks point for the *x* axis,
            *y* axis, and both axes.
        xgrid, ygrid, grid : None or bool, optional
            Whether to draw major gridlines on the *x* axes, *y* axes,
            and both axes.
        xgridminor, ygridminor, gridminor : None or bool, optional
            Whether to draw minor gridlines on the *x* axis, *y* axis, and
            both axes.
        xticklabeldir, yticklabeldir : {'in', 'out'}
            Whether to place *x* and *y* axis tick label text inside
            or outside the axes.
        xlocator, ylocator : None or locator spec, optional
            The *x* and *y* axis tick mark positions. Passed to the
            `~proplot.axistools.Locator` constructor.
        xlocator_kw, ylocator_kw : dict-like, optional
            The *x* and *y* axis locator settings. Passed to
            `~proplot.axistools.Locator`.
        xticks, yticks
            Aliases for `xlocator`, `ylocator`.
        xminorlocator, yminorlocator
            As for `xlocator`, `ylocator`, but for the minor ticks.
        xminorticks, yminorticks
            Aliases for `xminorlocator`, `yminorlocator`.
        xminorlocator_kw, yminorlocator_kw
            As for `xlocator_kw`, `ylocator_kw`, but for the minor locator.
        xformatter, yformatter : None or format spec, optional
            The *x* and *y* axis tick label format. Passed to the
            `~proplot.axistools.Formatter` constructor.
        xticklabels, yticklabels
            Aliases for `xformatter`, `yformatter`.
        xformatter_kw, yformatter_kw : dict-like, optional
            The *x* and *y* axis formatter settings. Passed to
            `~proplot.axistools.Formatter`.

        Todo
        ----
        * Consider redirecting user to another label rather than making
          this one invisible for spanning and shared axes.
        * More intelligent axis sharing with map axes. Consider optionally
          anchoring them by the locator and limits like the default API, or just
          disable ticklabels and labels for some.

        See also
        --------
        `BaseAxes.format`, `~proplot.axistools.Scale`,
        `~proplot.axistools.Locator`, `~proplot.axistools.Formatter`
        """
        # Set axis scaling and limits
        # These do not seem to have their own axes-specific public methods,
        # so do the x/y one by one here
        # WARNING: Special override here, force custom formatter when scale
        # is changed to log and user uses 'xlocator'/'ylocator'! Generally
        # means want specific tick labels on log-scale plot, but log formatter
        # will *override* and only show powers of 10.
        rc._getitem_mode = 0 # might still be non-zero if had error
        if xscale is not None:
            if hasattr(xscale,'name'):
                xscale = xscale.name
            self.set_xscale(axistools.Scale(xscale, **xscale_kw))
            if xscale in ('log','inverse') and xlocator is not None and xformatter is None:
                xformatter = 'default'
        if yscale is not None:
            if hasattr(yscale,'name'):
                yscale = yscale.name
            self.set_yscale(axistools.Scale(yscale, **yscale_kw))
            if yscale in ('log','inverse') and ylocator is not None and yformatter is None:
                yformatter = 'default'
        if xlim is not None:
            if xreverse:
                xlim = xlim[::-1]
            self.set_xlim(xlim)
        if ylim is not None:
            if yreverse:
                ylim = ylim[::-1]
            self.set_ylim(ylim)
        if (xlim is not None or ylim is not None) and self._inset_parent:
            self.indicate_inset_zoom()

        # Control axis ticks and labels and stuff
        # Allow for flexible input
        xspineloc = _default(xloc, xspineloc)
        yspineloc = _default(yloc, yspineloc)
        xformatter = _default(xticklabels, xformatter)
        yformatter = _default(yticklabels, yformatter)
        xlocator = _default(xticks, xlocator)
        ylocator = _default(yticks, ylocator)
        xminorlocator = _default(xminorticks, xminorlocator)
        yminorlocator = _default(yminorticks, yminorlocator)
        xtickminor = _default(tickminor, xtickminor)
        ytickminor = _default(tickminor, ytickminor)
        xgrid = _default(grid, xgrid)
        ygrid = _default(grid, ygrid)
        xgridminor = _default(gridminor, xgridminor)
        ygridminor = _default(gridminor, ygridminor)
        xtickdir = _default(tickdir, xtickdir)
        ytickdir = _default(tickdir, ytickdir)
        xticklabeldir = _default(ticklabeldir, xticklabeldir)
        yticklabeldir = _default(ticklabeldir, yticklabeldir)
        # Override for weird bug where title doesn't get automatically offset
        # from ticklabels in certain circumstance; check out notebook
        xtickloc = _default(xtickloc, xticklabelloc) # if user specified labels somewhere, make sure to put ticks there by default!
        ytickloc = _default(ytickloc, yticklabelloc)
        if xtickloc=='both' and xticklabelloc in ('both','top') and not xlabel: # xtickloc *cannot* be 'top', *only* appears for 'both'
            print('Warning: This keyword combination causes matplotlib bug where title is not offset from tick labels. Try adding an x-axis label or ticking only the top axis.')
        # Begin loop
        for axis, label, tickloc, spineloc, ticklabelloc, labelloc, bounds, gridminor, tickminor, tickminorlocator, \
                grid, ticklocator, tickformatter, tickrange, tickdir, ticklabeldir, \
                label_kw, formatter_kw, locator_kw, minorlocator_kw in \
            zip((self.xaxis, self.yaxis), (xlabel, ylabel),
                (xtickloc,ytickloc), (xspineloc, yspineloc), # other stuff
                (xticklabelloc, yticklabelloc), (xlabelloc, ylabelloc),
                (xbounds, ybounds),
                (xgridminor, ygridminor), (xtickminor, ytickminor), (xminorlocator, yminorlocator), # minor ticks
                (xgrid, ygrid),
                (xlocator, ylocator), (xformatter, yformatter), # major ticks
                (xtickrange, ytickrange), # range in which we label major ticks
                (xtickdir, ytickdir), (xticklabeldir, yticklabeldir), # tick direction
                (xlabel_kw, ylabel_kw), (xformatter_kw, yformatter_kw), (xlocator_kw, ylocator_kw), (xminorlocator_kw, yminorlocator_kw),
                ):
            # NOTE: Some of these settings are also rc settings, but I think a
            # good rule of thumb is format() methods should control toggling of
            # features, while _rcupdate() controls the look of those features.
            # Example: Set spine/tick locations with this func, but control
            # color/linewidth through _rcupdate().
            # TODO: Maybe the '_kw' stuff should be done in _rcupdate()?
            # Axis spine visibility and location
            sides = ('bottom','top') if axis.axis_name=='x' else ('left','right')
            spines = [self.spines[s] for s in sides]
            for spine, side in zip(spines, sides):
                # Line properties
                spineloc = getattr(self, f'twin_{axis.axis_name}spine_override', spineloc) # optionally override; necessary for twinx/twiny situation
                # Override if we're settings spine bounds
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
            spines = [side for side,spine in zip(sides,spines) if spine.get_visible()]

            # Set the major and minor locators and formatters
            # Also automatically detect whether axis is a 'time axis' (i.e.
            # whether user has plotted something with x/y as datetime/date/np.datetime64
            # objects, and matplotlib automatically set the unit converter)
            time = isinstance(axis.converter, mdates.DateConverter)
            if ticklocator is not None:
                axis.set_major_locator(axistools.Locator(ticklocator, time=time, **locator_kw))
            if tickformatter is not None or tickrange is not None:
                axis.set_major_formatter(axistools.Formatter(tickformatter, tickrange=tickrange, time=time, **formatter_kw))
            if not tickminor and tickminorlocator is None:
                axis.set_minor_locator(axistools.Locator('null'))
            elif tickminorlocator is not None:
                locator = axistools.Locator(tickminorlocator, minor=True, time=time, **minorlocator_kw)
                axis.set_minor_locator(locator)
            axis.set_minor_formatter(mticker.NullFormatter())

            # Tick properties
            # * Weird issue seems to cause set_tick_params to reset/forget that the grid
            #   is turned on if you access tick.gridOn directly, instead of passing through tick_params.
            #   Since gridOn is undocumented feature, don't use it. So calling _format_axes() a second time will remove the lines
            # * Can specify whether the left/right/bottom/top spines get ticks; sides will be 
            #   group of left/right or top/bottom
            # * Includes option to draw spines but not draw ticks on that spine, e.g.
            #   on the left/right edges
            # First determine tick sides
            ticklocs_kw = {None: None, 'both': sides, 'neither': (), 'none': ()}
            if bounds is not None and tickloc not in sides:
                tickloc = sides[0] # override to just one side
            ticklocs = ticklocs_kw.get(tickloc, (tickloc,))
            if ticklocs is None:
                ticks_sides = {}
            else:
                ticks_sides = {side: (side in ticklocs) for side in sides}
            ticks_sides.update({side: False for side in sides if side not in spines}) # override
            # Next the tick label sides
            # Will override to make sure sides match
            ticklabellocs = ticklocs_kw.get(ticklabelloc, (ticklabelloc,))
            if ticklabellocs is None:
                ticklabels_sides = {}
            else:
                ticklabels_sides = {'label' + side: (side in ticklabellocs) for side in sides}
            ticklabels_sides.update({'label' + side: False for side in sides
                if (side not in spines or (ticklocs is not None and side not in ticklocs))}) # override
            # Finally the label side
            if labelloc is None:
                if ticklocs is not None:
                    options = [side for side in sides if (side in ticklocs and side in spines)]
                    if len(options)==1:
                        labelloc = options[0]
            elif labelloc not in sides:
                raise ValueError('Got labelloc "{labelloc}", valid options are {sides}.')
            if labelloc is not None:
                axis.set_label_position(labelloc)
            # Apply settings to ticks
            ticks_major, ticks_minor = {}, {}
            if tickdir is not None:
                ticks_major.update({'direction':tickdir})
                ticks_minor.update({'direction':tickdir})
            if tickdir=='in':
                ticks_major.update({'pad':1}) # ticklabels should be much closer
                ticks_minor.update({'pad':1})
            if ticklabeldir=='in': # put tick labels inside the plot; sometimes might actually want this
                pad = rc['xtick.major.size'] + rc['xtick.major.pad'] + rc['xtick.labelsize']
                ticks_major.update({'pad':-pad})
                ticks_minor.update({'pad':-pad})
            axis.set_tick_params(which='major', **ticks_sides, **ticklabels_sides, **ticks_major)
            axis.set_tick_params(which='minor', **ticks_sides, **ticklabels_sides, **ticks_minor) # have length

            # Ensure no out-of-bounds ticks! Even set_smart_bounds() does not
            # always fix this! Need to try manual approach.
            # NOTE: set_bounds also failed, and fancy method overrides did
            # not work, so instead just turn locators into fixed version
            # NOTE: most locators take no arguments in call(), and some have
            # no tick_values method; so do the following
            # TODO: add optional override to do this every time
            if bounds is not None or axis.get_scale()=='cutoff':
                if bounds is None: # no API for this on axis
                    bounds = getattr(self, 'get_' + axis.axis_name + 'lim')()
                locator = axistools.Locator([x for x in axis.get_major_locator()() if bounds[0] <= x <= bounds[1]])
                axis.set_major_locator(locator)
                locator = axistools.Locator([x for x in axis.get_minor_locator()() if bounds[0] <= x <= bounds[1]])
                axis.set_minor_locator(locator)

            # Axis label properties
            # First redirect user request to the correct *shared* axes, then
            # redirect to the correct *spanning* axes if the label is meant
            # to span multiple subplot
            if label is not None:
                # Shared and spanning axes; try going a few layers deep
                # The _span_label method changes label position so it spans axes
                # If axis spanning not enabled, will just return the shared axis
                label_text = label
                label = self._share_span_label(axis)
                label.update({'text':label_text, **label_kw})
                if axis.get_label_position() == 'top':
                    label.set_va('bottom') # baseline was cramped if no ticklabels present

            # Gridline activation and setting (necessary because rcParams has no 'minorgrid'
            # property, must be set in rcSpecial settings)
            # NOTE: Inexplicably, for a twinx axis, could only get the minor gridlines
            # to disappear if we changed the 'visible' property on each one.
            # for tick in axis.get_major_ticks():
            #     if grid is not None:
            #         tick.gridline.set_visible(grid)
            #     tick.gridline.update(rc['grid']) # already set but why not, for symmetry
            # # for tick in axis.minorTicks:
            # for tick in axis.get_minor_ticks():
            #     if gridminor is not None:
            #         tick.gridline.set_visible(gridminor)
            #     tick.gridline.update(rc['gridminor'])
            # For some insane reasion, these are ***both*** needed
            # Without this below stuff, e.g. gridminor=True doesn't draw gridlines
            if grid is not None: # grid changes must be after tick
                axis.grid(grid, which='major')
            if gridminor is not None:
                axis.grid(gridminor, which='minor') # ignore if no minor ticks

        # Pass stuff to parent formatter, e.g. title and abc labeling
        super().format(**kwargs)

    def twiny(self, **kwargs):
        """
        Create second *x* axis extending from a shared ("twin") *y*
        axis. Adds features for intelligently moving around tick and axis
        labels, handling axis sharing.
        """
        # Create second x-axis extending from shared ("twin") y-axis
        # Note: Cannot wrap twiny() because then the axes created will be
        # instantiated from the parent class, which doesn't have format() method.
        # Instead, use hidden method _make_twin_axes.
        # See https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/axes/_subplots.py
        ax = self._make_twin_axes(sharey=self, projection=self.name)
        self.xaxis.tick_bottom()
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_autoscaley_on(self.get_autoscaley_on())
        ax.yaxis.set_visible(False)
        ax.patch.set_visible(False)
        ax.grid(False)
        # Special settings, force spine locations when format() called
        self.twin_xspine_override = 'bottom' # original axis ticks on bottom
        ax.twin_xspine_override   = 'top' # new axis ticks on top
        ax.twin_yspine_override   = 'neither'
        return ax

    def twinx(self, yscale=None, **kwargs):
        """
        Create second *y* axis extending from a shared ("twin") *x*
        axis. Adds features for intelligently moving around tick and axis
        labels, handling axis sharing.
        """
        # Create second y-axis extending from shared ("twin") x-axis
        # Note: Cannot wrap twinx() because then the axes created will be
        # instantiated from the parent class, which doesn't have format() method.
        # Instead, use hidden method _make_twin_axes.
        ax = self._make_twin_axes(sharex=self, projection=self.name)
        self.yaxis.tick_left()
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')
        ax.yaxis.set_offset_position('right')
        ax.set_autoscalex_on(self.get_autoscalex_on())
        ax.xaxis.set_visible(False)
        ax.patch.set_visible(False)
        ax.grid(False)
        # Apply scale and axes sharing
        # NOTE: Forget about this because often (for height/pressure scales)
        # the units won't match because we didn't use p0
        # if yscale:
        #     transform = mscale.scale_factory(yscale, self.yaxis).get_transform()
        #     lims = self.get_ylim()
        #     lims = transform.transform(np.array(lims))
        #     ax.set_ylim(lims)
        # Special settings, force spine locations when format() called
        self.twin_yspine_override = 'left' # original axis ticks on left
        ax.twin_yspine_override   = 'right' # new axis ticks on right
        ax.twin_xspine_override   = 'neither'
        return ax

    def _make_inset_locator(self, bounds, trans):
        # Helper function, had to be copied from private matplotlib version.
        def inset_locator(ax, renderer):
            bbox = mtransforms.Bbox.from_bounds(*bounds)
            bb = mtransforms.TransformedBbox(bbox, trans)
            tr = self.figure.transFigure.inverted()
            bb = mtransforms.TransformedBbox(bb, tr)
            return bb
        return inset_locator

    def inset_axes(self, bounds, *, transform=None, zorder=5, zoom=True, zoom_kw={}, **kwargs):
        """
        Draw an inset `XYAxes` axes. Otherwise, this is a carbon copy
        of the `~matplotlib.axes.Axes.inset_axes` method.
        """
        # Carbon copy, but use my custom axes
        # Defaults
        if transform is None:
            transform = self.transAxes
        label = kwargs.pop('label', 'inset_axes')
        # This puts the rectangle into figure-relative coordinates.
        locator = self._make_inset_locator(bounds, transform)
        bb = locator(None, None)
        ax = XYAxes(self.figure, bb.bounds, zorder=zorder, label=label, **kwargs)
        # The following locator lets the axes move if in data coordinates, gets called in ax.apply_aspect()
        ax.set_axes_locator(locator)
        self.add_child_axes(ax)
        self._insets += [ax]
        ax._inset_parent = self
        # Finally add zoom
        # NOTE: Requirs version >=3.0
        if zoom:
            ax.indicate_inset_zoom(**zoom_kw)
        return ax

    def indicate_inset_zoom(self, alpha=None, linewidth=None, color=None, edgecolor=None, **kwargs):
        """
        Custom inset zoom indicator that can be *refreshed* as the
        parent axis limits change.
        """
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
        kwargs.update({'linewidth': linewidth, 'edgecolor':edgecolor, 'alpha':alpha})
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

    def inset(self, *args, **kwargs):
        """
        Alias for `~XYAxes.inset_axes`.
        """
        # Just an alias
        return self.inset_axes(*args, **kwargs)

    def inset_zoom(self, *args, **kwargs):
        """
        Alias for `~XYAxes.indicate_inset_zoom`.
        """
        # Just an alias
        return self.indicate_inset_zoom(*args, **kwargs)

class EmptyPanel(object):
    """
    Dummy object to put in place when an axes or figure panel does not exist.
    This gives a nicer error message than if we had just used ``None`` or put
    nothing there at all.

    Notes
    -----
    `__getattr__` is invoked only when `__getattribute__` fails, i.e.
    when the user requests anything that isn't a hidden `object` method.
    """
    def __bool__(self):
        return False # it's empty, so this is 'falsey'

    def __getattr__(self, attr, *args):
        raise AttributeError('Panel does not exist.')

class PanelAxes(XYAxes):
    """
    An `~proplot.axes.XYAxes` with `legend` and `colorbar` methods
    overridden. Calling these will "fill" the entire axes with a legend
    or colorbar.

    This is suitable for axes meant to reference content in several
    other subplots at once.

    Notes
    -----
    See `this post <https://stackoverflow.com/a/52121237/4970632>`_
    and `this example <https://stackoverflow.com/q/26236380/4970632>`_.

    See also
    --------
    `~proplot.subplots.subplots`, `~proplot.subplots.FigureBase.panel_factory`, `XYAxes`, `BaseAxes`
    """
    # Name
    name = 'panel'
    """The registered projection name."""
    def __init__(self, *args, panel_side=None, invisible=False, **kwargs):
        """
        Parameters
        ----------
        panel_side : {'left', 'right', 'bottom', 'top'}
            The side on which the panel is drawn.
        invisible : bool, optional
            Whether to make the axes invisible at the start. This is the
            default when using e.g. the `bottomlegend` keyword arg
            in `~proplot.subplots.subplots`.
        *args, **kwargs
            Passed to the `XYAxes` initializer.
        """
        # Initiate
        if panel_side is None:
            raise ValueError('Must specify side.')
        super().__init__(*args, panel_side=panel_side, **kwargs)
        # Make everything invisible
        if invisible:
            self.invisible()

    # TODO: Should have two versions of both of these
    # 1) Adds a legend or *small* colorbar to the axes
    # 2) *Fills* the entire axes with a colorbar or legend box
    # TODO: Get rid of legend and colorbar factories, implement them as
    # direct overrides!
    def legend(self, handles, entire=True, **kwargs_override):
        """
        "Fill the panel" with a legend. That is, draw a centered legend
        and make the axes spines, ticks, etc. invisible.
        """
        if not entire:
            # Do the normal thing
            kwargs = kwargs_override
        else:
            # Allocate invisible axes for drawing legend.
            # Returns the axes and the output of legend_factory().
            self.invisible()
            kwargs = {'borderaxespad':  0,
                    'frameon':        False,
                    'loc':            'upper center',
                    'bbox_transform': self.transAxes}
            kwargs.update(kwargs_override)
        return legend_factory(self, handles, **kwargs)

    def colorbar(self, *args, i=0, n=1, length=1,
        space=0, hspace=None, wspace=None,
        **kwargs):
        """
        Fill the panel with colorbar.
        """
        # Draw colorbar with arbitrary length relative to full length of the
        # panel, and optionally *stacking* multiple colorbars
        # Will always redraw an axes with new subspec
        self.invisible()
        side = self.panel_side
        space = _default(hspace, _default(wspace, space)) # flexible arguments
        figure = self.figure
        subspec = self.get_subplotspec()

        # Create colorbar axes
        # if n>2:
        #     raise ValueError('I strongly advise against drawing more than 2 stacked colorbars.')
        if length!=1 or n!=1:
            # First get gridspec
            # Note formula: total width = n*<colorbar width> + (n-1)*<space width>
            if side in ['bottom','top']:
                hwidth = (self.height - (n-1)*space)/n # express height ratios in physical units
                if hwidth<0:
                    raise ValueError(f'Space {space} too big for {n} colorbars on panel with width {self.height}.')
                gridspec = FlexibleGridSpecFromSubplotSpec(
                        nrows=n,  ncols=3,
                        wspace=0, hspace=space,
                        subplot_spec=subspec,
                        width_ratios=((1-length)/2, length, (1-length)/2),
                        height_ratios=hwidth,
                        )
                subspec = gridspec[i,1]
            elif side in ['left','right']:
                wwidth = (self.width - (n-1)*space)/n
                if wwidth<0:
                    raise ValueError(f'Space {space} too big for {n} colorbars on panel with width {self.width}.')
                gridspec = FlexibleGridSpecFromSubplotSpec(
                        nrows=3,  ncols=n,
                        wspace=wspace, hspace=hspace,
                        subplot_spec=subspec,
                        height_ratios=((1-length)/2, length, (1-length)/2),
                        width_ratios=wwidth,
                        )
                subspec = gridspec[1,i]
            # Next redraw axes
            # self.remove() # save memory
            self.set_visible(False)

        # Allocate axes for drawing colorbar.
        # Returns the axes and the output of colorbar_factory().
        # WARNING: Using BaseAxes for the colorbar axes seems to cause bugs,
        # because colorbar uses some internal methods that are wrapped by
        # cmap_features! Not exactly sure why, but bottom line, do not use
        # custom projection.
        # ax = figure.add_subplot(subspec, projection='base')
        ax = figure.add_subplot(subspec, projection=None)
        if side in ['bottom','top']:
            outside, inside = 'bottom', 'top'
            if side=='top':
                outside, inside = inside, outside
            # ticklocation = outside if i==n-1 else inside
            ticklocation = outside
            orientation  = 'horizontal'
        elif side in ['left','right']:
            outside, inside = 'left', 'right'
            if side=='right':
                outside, inside = inside, outside
            # ticklocation = outside if i==n-1 else inside
            ticklocation = outside
            orientation  = 'vertical'
        kwargs.update({'orientation':orientation, 'ticklocation':ticklocation})
        return colorbar_factory(ax, *args, **kwargs)

class MapAxes(BaseAxes):
    """
    Dummy intermediate class. Disables a bunch of methods that are
    inappropriate for map projections.

    See also
    --------
    `~proplot.subplots.subplots`, `BaseAxes`, `CartopyAxes`, `BasemapAxes`
    """
    # Disable some methods to prevent weird shit from happening
    # Originally used property decorators for this but way too verbose
    # See: https://stackoverflow.com/a/23126260/4970632
    def __getattribute__(self, attr, *args):
        if attr in _map_disabled_methods:
            raise NotImplementedError('Invalid plotting function {} for map projection axes.'.format(attr))
        return super().__getattribute__(attr, *args)

    @staticmethod
    def ls_translate(obj, style):
        """
        Make basemap gridlines look like cartopy lines using the `dashes`
        tuple.

        Notes
        -----
        For some reason basemap gridlines look different from cartopy ones.
        Have absolutely **no idea** why. Cartopy seems to do something weird because
        there is no ``_dashSeq`` attribute on the lines, and the line styles
        are always "officially" ``'-'``, even though they are dashed.
        See `this reference <https://matplotlib.org/gallery/lines_bars_and_markers/line_styles_reference.html>`_.

        The dots ``':'`` actually look better on cartopy so we try to mimick them
        below.
        """
        if style=='-':
            dashes = [None,None]
        else:
            dashes = [*obj._dashSeq]
            if style==':':
                dashes[0] /= 10
                dashes[1] *= 1.5
            elif style=='--':
                dashes[0] /= 1.5
                dashes[1] *= 1.5
            else:
                raise ValueError(f'Invalid style {style}.')
        return dashes

    @staticmethod
    def parse_labels(labels, mode):
        """
        Parses complex ``lonlabels`` and ``latlabels`` arguments. There are
        four different options:

        1. A string, e.g. ``'lr'`` or ``'bt'``.
        2. Boolean ``True``. Indicates left side for latitudes,
           bottom for longitudes.
        3. A boolean ``(left,right)`` tuple for latitudes,
           ``(bottom,top)`` for longitudes.
        4. A boolean ``(n1,n2,n3,n4)`` tuple as in the
           `~mpl_toolkits.basemap.Basemap.drawmeridians` and
           `~mpl_toolkits.basemap.Basemap.drawparallels` methods.
           Indicates whether to label left, right, top, and bottom
           sides, respectively.
        """
        if labels is False:
            return [0]*4
        if labels is None:
            labels = True # use the default
        if isinstance(labels, str):
            string = labels
            labels = [0]*4
            for idx,char in zip([0,1,2,3],'lrbt'):
                if char in string:
                    labels[idx] = 1
        if isinstance(labels, Number): # e.g. *boolean*
            labels = np.atleast_1d(labels)
        if len(labels)==1:
            labels = [*labels, 0] # default is to label bottom/left
        if len(labels)==2:
            if mode=='x':
                labels = [0, 0, *labels]
            elif mode=='y':
                labels = [*labels, 0, 0]
        elif len(labels)!=4:
            raise ValueError(f'Invalid labels: {labels}.')
        return labels

class BasemapAxes(MapAxes):
    """
    Axes subclass for plotting `basemap <https://matplotlib.org/basemap/>`_
    projections. The `~mpl_toolkits.basemap.Basemap` projection instance is added as
    the `m` attribute, but this is all abstracted away -- you can use
    `~matplotlib.axes.Axes` methods like `~matplotlib.axes.Axes.plot` and
    `~matplotlib.axes.Axes.contour` with
    your raw longitude-latitude data.

    See also
    --------
    `~proplot.subplots.subplots`, `~proplot.proj`, `BaseAxes`, `MapAxes`
    """
    name = 'basemap'
    # Note non-rectangular projections; for rectnagular ones, axes spines are
    # used as boundaries, but for these, have different boundary.
    _proj_non_rectangular = [
            'ortho', 'geos', 'nsper',
            'moll', 'hammer', 'robin',
            'eck4', 'kav7', 'mbtfpq', # last one is McBryde-Thomas flat polar quartic
            'sinu', 'vandg', # last one is van der Grinten
            ]
    """The registered projection name."""
    def __init__(self, *args, map_projection=None, **kwargs):
        # * Must set boundary before-hand, otherwise the set_axes_limits method called
        #   by mcontourf/mpcolormesh/etc draws two mapboundary Patch objects called "limb1" and
        #   "limb2" automatically: one for fill and the other for the edges
        # * Then, since the patch object in _mapboundarydrawn is only the fill-version, calling
        #   drawmapboundary() again will replace only *that one*, but the original visible edges
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
        self._land = None
        self._ocean = None
        self._coastline = None
        # Initialize
        super().__init__(*args, map_name=self.m.projection, **kwargs)

    def _rcupdate(self):
        # Map boundary
        # * First have to *manually replace* the old boundary by just deleting
        #   the original one
        # * If boundary is drawn successfully should be able to call
        #   self.m._mapboundarydrawn.set_visible(False) and edges/fill color disappear
        # * For now will enforce that map plots *always* have background whereas
        #   axes plots can have transparent background
        kw_face = rc.fill({'facecolor': 'map.facecolor'})
        self.axesPatch = self.patch # bugfix or something
        if self.m.projection in self._proj_non_rectangular:
            self.patch.set_alpha(0) # make patch invisible
            kw_edge = rc.fill({'linewidth': 'map.linewidth', 'edgecolor': 'map.edgecolor'}) # necessary to update after drawn, because patch 'color' is the fill but kwarg for edge color is 'color'
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
            kw_edge = rc.fill({'linewidth': 'map.linewidth', 'color': 'map.edgecolor'})
            for spine in self.spines.values():
                spine.update(kw_edge)

        # Call parent
        super()._rcupdate()

    # Basemap overrides
    # WARNING: Never ever try to just make blanket methods on the Basemap
    # instance accessible from axes instance! Can of worms and had bunch of
    # weird errors! Just pick the ones you think user will want to use.
    def __getattribute__(self, attr, *args):
        if attr=='pcolorpoly': # need to specify this again to access the .m method
            attr = 'pcolor' # use alias so don't run into recursion issues due to internal pcolormesh calls to pcolor()
        obj = super().__getattribute__(attr, *args)
        if attr in _line_methods or attr in _edge_methods or attr in _center_methods:
            obj = _m_call(self, obj) # this must be the *last* step!
            if attr in _line_methods:
                if attr[:3] != 'tri':
                    obj = cycle_features(self, obj)
                obj = _linefix_basemap(self, obj)
            elif attr in _edge_methods or attr in _center_methods:
                obj = cmap_features(self, obj)
                obj = _gridfix_basemap(self, obj)
                if attr in _edge_methods:
                    obj = _check_edges(obj)
                else:
                    obj = _check_centers(obj)
            obj = _no_recurse(self, obj)
        return obj

    # Format basemap axes
    # Add documentation here.
    def format(self,
        lonlim=None, latlim=None, xlim=None, ylim=None,
        lonticks=None, lonminorticks=None, lonlocator=None, lonminorlocator=None,
        latticks=None, latminorticks=None, latlocator=None, latminorlocator=None,
        xticks=None, xminorticks=None, xlocator=None, xminorlocator=None,
        yticks=None, yminorticks=None, ylocator=None, yminorlocator=None,
        land=False, ocean=False, coastline=False, # coastlines and land
        latlabels=None, lonlabels=None, # sides for labels [left, right, bottom, top]
        xlabels=None, ylabels=None,
        landcolor=None, oceancolor=None, coastcolor=None, coastlinewidth=None,
        **kwargs):
        """
        Format the map projection.

        Parameters
        ----------
        lonlim, latlim : None or length-2 list of float, optional
            Longitude and latitude limits of projection.
        xlim, ylim
            Aliases for `lonlim`, `latlim`.
        lonlocator, latlocator : None or list of float, optional
            List of longitudes and latitudes for drawing gridlines.
        xlocator, ylocator, lonticks, latticks, xticks, yticks
            Aliases for `lonlocator`, `latlocator`.
        lonminorlocator, latminorlocator : None or list of float, optional
            As with `lonlocator` and `latlocator`, but for minor gridlines.
        xminorlocator, yminorlocator, lonminorticks, latminorticks, xminorticks, yminorticks
            Aliases for `lonminorlocator`, `latminorlocator`.
        lonlabels, latlabels
            Positions for drawing longitude and latitude coordinate labels
            for the gridlines. Interpreted by `~MapAxes.parse_labels`.
        xlabels, ylabels
            Aliases for `lonlabels`, `latlabels`.

        Other parameters
        ----------------
        **kwargs
            Passed to the parent `BaseAxes.format`.

        See also
        --------
        `~proplot.subplots.subplots`, `~proplot.rcmod`, `BaseAxes.format`
        """
        # Parse flexible input
        xlim = _default(lonlim, xlim)
        ylim = _default(latlim, ylim)
        lonlocator = _default(lonlocator, _default(lonticks, _default(xlocator, xticks)))
        latlocator = _default(latlocator, _default(latticks, _default(ylocator, yticks)))
        lonminorlocator = _default(lonminorlocator, _default(lonminorticks, _default(xminorlocator, xminorticks)))
        latminorlocator = _default(latminorlocator, _default(latminorticks, _default(yminorlocator, yminorticks)))
        lonlabels = self.parse_labels(_default(xlabels, lonlabels), 'x')
        latlabels = self.parse_labels(_default(ylabels, latlabels), 'y')

        # Basemap axes setup
        # Coastlines, parallels, meridians
        if land and not self._land:
            color = _default(landcolor, rc.landcolor)
            self._land = self.m.fillcontinents(ax=self, color=color)
        if coastline and not self._coastline:
            lw = _default(coastlinewidth, rc.coastlinewidth)
            color = _default(coastcolor, rc.coastcolor)
            self._coastline = self.m.drawcoastlines(ax=self, color=color, lw=lw)

        # Longitude/latitude lines
        # Make sure to turn off clipping by invisible axes boundary; otherwise
        # get these weird flat edges where map boundaries, parallel/meridian markers come up to the axes bbox
        tsettings = {'color':rc['xtick.color'], 'fontsize':rc['xtick.labelsize']}
        latlabels[2:] = latlabels[2:][::-1] # default is left/right/top/bottom which is dumb
        lonlabels[2:] = lonlabels[2:][::-1] # change to left/right/bottom/top
        lsettings = rc['lonlatlines']
        linestyle = lsettings['linestyle']
        latlocator = _default(latlocator, 20) # gridlines by default
        lonlocator = _default(lonlocator, 60)
        if latlocator is not None:
            if isinstance(latlocator, Number):
                latlocator = utils.arange(self.m.latmin+latlocator, self.m.latmax-latlocator, latlocator)
            p = self.m.drawparallels(latlocator, labels=latlabels, ax=self)
            for pi in p.values(): # returns dict, where each one is tuple
                # Tried passing clip_on to the below, but it does nothing; must set
                # for lines created after the fact
                for obj in [i for j in pi for i in j]: # magic
                    if isinstance(obj, mtext.Text):
                        obj.update(tsettings)
                    else:
                        obj.update(lsettings)
                        obj.set_dashes(self.ls_translate(obj, linestyle))
        if lonlocator is not None:
            if isinstance(lonlocator, Number):
                lonlocator = utils.arange(self.m.lonmin+lonlocator, self.m.lonmax-lonlocator, lonlocator)
            p = self.m.drawmeridians(lonlocator, labels=lonlabels, ax=self)
            for pi in p.values():
                for obj in [i for j in pi for i in j]: # magic
                    if isinstance(obj, mtext.Text):
                        obj.update(tsettings)
                    else:
                        obj.update(lsettings)
                        obj.set_dashes(self.ls_translate(obj, linestyle))

        # Pass stuff to parent formatter, e.g. title and abc labeling
        super().format(**kwargs)

# Cartopy takes advantage of documented feature where any class with method
# named _as_mpl_axes can be passed as 'projection' object.
class CartopyAxes(MapAxes, GeoAxes): # custom one has to be higher priority, so the methods can overwrite stuff
# Feature documented here: https://matplotlib.org/devel/add_new_projection.html
    """
    Axes subclass for plotting `cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_
    projections. Initializes the `cartopy.crs.Projection` instance. Also
    allows for *partial* coverage of azimuthal projections by zooming into
    the full projection, then drawing a circle boundary around some latitude
    away from the center (this is surprisingly difficult to do).

    Parameters
    ----------
    map_projection : str
        String name for projection.
    circle_center : float, optional
        For polar projections, the center latitude of the circle (``-90``
        or ``90``).
    circle_edge : float, optional
        For polar projections, the edge latitude of the circle.

    Notes
    -----
    The circle stuff for polar projection was developed from `this example
    <https://scitools.org.uk/cartopy/docs/v0.15/examples/always_circular_stereo.html>`_.

    See also
    --------
    `~proplot.subplots.subplots`, `~proplot.proj`, `BaseAxes`, `MapAxes`
    """
    # Used in Projection parent class here: https://scitools.org.uk/cartopy/docs/v0.13/_modules/cartopy/crs
    name = 'cartopy'
    """The registered projection name."""
    def __init__(self, *args, map_projection=None, circle_center=90, circle_edge=0, **kwargs):
        # Dependencies
        import cartopy.crs as ccrs # verify package is available

        # Do the GeoAxes initialization steps manually (there are very few)
        if not isinstance(map_projection, ccrs.Projection):
            raise ValueError('You must initialize CartopyAxes with map_projection=(cartopy.crs.Projection instance).')
        self._hold = None # dunno
        self._land = None
        self._ocean = None
        self._coastline = None
        self.projection = map_projection # attribute used extensively by GeoAxes methods, and by builtin one

        # Below will call BaseAxes, which will call GeoAxes as the superclass
        # NOTE: Previously did stuff in __init__ manually, and called self._boundary,
        # which hides existing border patch and rewrites as None. Don't do that again.
        try:
            map_name = map_projection.name
        except AttributeError:
            map_name = map_projection.proj4_params['proj']
        super().__init__(*args, map_projection=map_projection, map_name=map_name, **kwargs)

        # Apply circle boundary
        crs_circles = (ccrs.LambertAzimuthalEqualArea, ccrs.AzimuthalEquidistant)
        if any(isinstance(map_projection, cp) for cp in crs_circles):
            self.set_extent([-180, 180, circle_edge, circle_center], PlateCarree()) # use platecarree transform
            self.set_boundary(proj.Circle(100), transform=self.transAxes)
            # self.projection.threshold = kwargs.pop('threshold', self.projection.threshold) # optionally modify threshold
        self.set_global() # see: https://stackoverflow.com/a/48956844/4970632

    def __getattribute__(self, attr, *args):
        obj = super().__getattribute__(attr, *args)
        if attr in _line_methods:
            obj = _linefix_cartopy(self, obj)
        elif attr in _edge_methods or attr in _center_methods:
            obj = _gridfix_cartopy(self, obj)
            if attr in _edge_methods:
                obj = _check_edges(obj)
            else:
                obj = _check_centers(obj)
        return obj

    def _rcupdate(self):
        # Update properties controlled by custom rc settings
        # WARNING: Seems cartopy features can't be updated!
        # See: https://scitools.org.uk/cartopy/docs/v0.14/_modules/cartopy/feature.html#Feature
        # Change the _kwargs property also does *nothing*
        self.set_global() # see: https://stackoverflow.com/a/48956844/4970632
        kw = rc.fill({'facecolor': 'map.facecolor'})
        self.background_patch.update(kw)
        kw = rc.fill({'edgecolor': 'map.edgecolor', 'linewidth': 'map.linewidth'})
        self.outline_patch.update(kw)

        # Call parent
        super()._rcupdate()

    # Format cartopy GeoAxes.
    # Add documentation here.
    def format(self,
        xlim=None, ylim=None, lonlim=None, latlim=None,
        xticks=None, xminorticks=None, xlocator=None, xminorlocator=None,
        yticks=None, yminorticks=None, ylocator=None, yminorlocator=None,
        latticks=None, latminorticks=None, latlocator=None, latminorlocator=None,
        lonticks=None, lonminorticks=None, lonlocator=None, lonminorlocator=None,
        land=False, ocean=False, coastline=False, # coastlines and continents
        landcolor=None, oceancolor=None, coastcolor=None, coastlinewidth=None,
        xlabels=None, ylabels=None,
        latlabels=None, lonlabels=None, # sides for labels [left, right, bottom, top]
        reso=None,
        **kwargs):
        """
        Format the map projection.

        Parameters
        ----------
        lonlim, latlim : None or length-2 list of float, optional
            Longitude and latitude limits of projection.
        xlim, ylim
            Aliases for `lonlim`, `latlim`.
        lonlocator, latlocator : None or list of float, optional
            List of longitudes and latitudes for drawing gridlines.
        xlocator, ylocator, lonticks, latticks, xticks, yticks
            Aliases for `lonlocator`, `latlocator`.
        lonminorlocator, latminorlocator : None or list of float, optional
            As with `lonlocator` and `latlocator`, but for minor gridlines.
        xminorlocator, yminorlocator, lonminorticks, latminorticks, xminorticks, yminorticks
            Aliases for `lonminorlocator`, `latminorlocator`.
        lonlabels, latlabels
            Positions for drawing longitude and latitude coordinate labels
            for the gridlines. Interpreted by `~MapAxes.parse_labels`.
        xlabels, ylabels
            Aliases for `lonlabels`, `latlabels`.
        reso : {None, 'lo', 'med', 'hi'}, optional
            Resolution for geographic features. Inferred from configuration
            if ``None``. See the `~proplot.rcmod` documentation.
        landcolor, oceancolor, coastcolor, coastlinewidth
            Geographic feature settings. If ``None``, read from the
            configuration. See the `~proplot.rcmod` documentation.

        Other parameters
        ----------------
        **kwargs
            Passed to the parent `BaseAxes.format`.

        See also
        --------
        `~proplot.subplots.subplots`, `~proplot.rcmod`, `BaseAxes.format`
        """
        # Dependencies
        import cartopy.feature as cfeature
        import cartopy.crs as ccrs # verify package is available
        from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

        # Parse flexible input
        xlim = _default(lonlim, xlim)
        ylim = _default(latlim, ylim)
        lonlocator = _default(lonlocator, _default(lonticks, _default(xlocator, xticks)))
        latlocator = _default(latlocator, _default(latticks, _default(ylocator, yticks)))
        lonminorlocator = _default(lonminorlocator, _default(lonminorticks, _default(xminorlocator, xminorticks)))
        latminorlocator = _default(latminorlocator, _default(latminorticks, _default(yminorlocator, yminorticks)))
        lonlabels = self.parse_labels(_default(xlabels, lonlabels), 'x')
        latlabels = self.parse_labels(_default(ylabels, latlabels), 'y')

        # Configure extents?
        # WARNING: The set extents method tries to set a *rectangle* between
        # the *4* (x,y) coordinate pairs (each corner), so something like
        # (-180,180,-90,90) will result in *vertical line*, causing error!
        # NOTE: proj4_params stores keyword-arg pairs, proj4_init stores
        # the shell string passed
        # NOTE: They may add this in set_xlim and set_ylim in the
        # near future; see:
        # https://github.com/SciTools/cartopy/blob/master/lib/cartopy/mpl/geoaxes.py#L638
        if xlim is not None or ylim is not None:
            xlim = xlim or [None,None]
            ylim = ylim or [None,None]
            xlim, ylim = [*xlim], [*ylim]
            lon_0 = self.projection.proj4_params.get('lon_0', 0)
            if xlim[0] is None:
                xlim[0] = lon_0 - 180
            if xlim[1] is None:
                xlim[1] = lon_0 + 180
            if ylim[0] is None:
                ylim[0] = -90
            if ylim[1] is None:
                ylim[1] = 90
            self.set_extent([*xlim, *ylim], PlateCarree())

        # Add geographic features
        # Use the NaturalEarthFeature to get more configurable resolution; can choose
        # between 10m, 50m, and 110m (scales 1:10mil, 1:50mil, and 1:110mil)
        reso = _default(reso, rc['map.reso'])
        if reso not in ('lo','med','hi'):
            raise ValueError(f'Invalid resolution {reso}.')
        reso = {'lo':'110m', 'med':'50m', 'hi':'10m'}.get(reso)
        if coastline and not self._coastline:
            # self.add_feature(cfeature.COASTLINE, **rc['coastlines'])
            feat  = cfeature.NaturalEarthFeature('physical', 'coastline', reso)
            color = _default(coastcolor, rc.coastcolor)
            lw    = _default(coastlinewidth, rc.coastlinewidth)
            self.add_feature(feat, color=color, lw=lw)
            self._coastline = feat
        if land and not self._land:
            # self.add_feature(cfeature.LAND, **rc['continents'])
            feat  = cfeature.NaturalEarthFeature('physical', 'land', reso)
            color = _default(landcolor, rc.landcolor)
            self.add_feature(feat, color=color)
            self._land = feat
        if ocean and not self._ocean:
            # self.add_feature(cfeature.OCEAN, **rc['oceans'])
            feat = cfeature.NaturalEarthFeature('physical', 'ocean', reso)
            color = _default(oceancolor, rc.oceancolor)
            self.add_feature(feat, color=oceancolor)
            self._ocean = feat

        # Draw gridlines
        # WARNING: For some reason very weird side effects happen if you try
        # to call gridlines() twice on same axes. Can't do it. Which is why
        # we do this nonsense with the formatter below, instead of drawing 'major'
        # grid lines and 'minor' grid lines.
        lonvec = lambda v: [] if v is None else [*v] if np.iterable(v) else [*utils.arange(-180,180,v)]
        latvec = lambda v: [] if v is None else [*v] if np.iterable(v) else [*utils.arange(-90,90,v)]
        lonminorlocator, latminorlocator = lonvec(lonminorlocator), latvec(latminorlocator)
        lonlocator, latlocator = lonvec(lonlocator), latvec(latlocator)
        lonlines = lonminorlocator or lonlocator # where we draw gridlines
        latlines = latminorlocator or latlocator

        # First take care of gridlines
        eps = 1e-3
        draw_labels = (isinstance(self.projection, ccrs.Mercator) or isinstance(self.projection, ccrs.PlateCarree))
        if latlines and latlines[0]==-90:
            latlines[0] += eps
        if lonlines and lonlines[0]==-90:
            lonlines[0] -= eps
        gl = self.gridlines(**rc['lonlatlines'], draw_labels=draw_labels)
        if lonlines: # NOTE: using mticker.NullLocator results in error!
            gl.xlocator = mticker.FixedLocator(lonlines)
        if latlines:
            gl.ylocator = mticker.FixedLocator(latlines)

        # Now take care of labels
        if draw_labels:
            self._gridliner_on = True
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlabels_bottom, gl.xlabels_top = lonlabels[2:]
            gl.ylabels_left, gl.ylabels_right = latlabels[:2]

        # Pass stuff to parent formatter, e.g. title and abc labeling
        super().format(**kwargs)

class PolarAxes(MapAxes, mproj.PolarAxes):
    """
    Thin decorator around `~matplotlib.projections.PolarAxes` with my
    new plotting features. So far, just mixes the two classes.

    Warning
    -------
    Polar axes have not been tested yet!

    See also
    --------
    `~proplot.subplots.subplots`, `BaseAxes`, `MapAxes`
    """
    name = 'newpolar'
    """The registered projection name."""

# Register the projections
mproj.register_projection(BaseAxes)
mproj.register_projection(XYAxes)
mproj.register_projection(PanelAxes)
mproj.register_projection(PolarAxes)
mproj.register_projection(BasemapAxes)
mproj.register_projection(CartopyAxes)

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
        Alias for ncol. Added for consistency with
        `~matplotlib.pyplot.subplots`.

    Other parameters
    ----------------
    **kwargs
        Passed to `~matplotlib.axes.Axes.legend`.

    Todo
    ----
    Should update this function to clip the legend box when it goes
    outside axes area, so the legend width and bottom or right widths can be
    chosen propertly/separately.

    See also
    --------
    `BaseAxes.colorbar`, `PanelAxes.colorbar`, `~matplotlib.axes.Axes.legend`
    """
    # First get legend settings (usually just one per plot so don't need to declare
    # this dynamically/globally), and interpret kwargs.
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
    # Overwrite alpha? Bad idea
    # hsettings.update({'alpha':1.0})
    # Font name; 'prop' can be a FontProperties object or a dict for the kwargs
    # to instantiate one.
    kwargs.update({'prop': {'family':rc['fontname']}})

    # Detect if user wants to specify rows manually
    # Gives huge latitude for user input:
    #   1) user can specify nothing and align will be inferred (list of iterables
    #      will always be False, i.e. we draw consecutive legends, and list of handles is always true)
    #   2) user can specify align (needs list of handles for True, list of handles or list
    #      of iterables for False and if the former, will turn into list of iterables)
    if handles is None:
        handles = ax.get_legend_handles_labels()[0]
    for i,handle in enumerate(handles):
        if hasattr(handle, 'cmap'):
            # Make sure we sample the *center* of the colormap
            print('Warning: Getting legend entry from colormap.')
            size = np.mean(handle.get_sizes())
            handles[i] = ax.scatter([0], [0],
                                 markersize=size,
                                 color=[handle.cmap(0.5)],
                                 label=handle.get_label())
    # handles = np.array(handles).squeeze().tolist()
    list_of_lists = not isinstance(handles[0], martist.Artist)
    if align is None: # automatically guess
        align = not list_of_lists
    else: # standardize format based on input
        if not align and not list_of_lists: # separate into columns
            if 'ncol' not in kwargs:
                kwargs['ncol'] = 3
            ncol = kwargs['ncol']
            list_of_lists = True
            handles = [handles[i*ncol:(i+1)*ncol]
                        for i in range(len(handles))] # to list of iterables
        if align and list_of_lists: # unfurl, because we just want one legend!
            list_of_lists = False
            handles = [handle for isiterable in handles for handle in isiterable]
            list_of_lists = False # no longer is list of lists

    # Now draw legend, with two options
    # 1) Normal legend, just draw everything like normal and columns
    # will be aligned; we re-order handles to be row-major, is only difference
    if align:
        # Prepare settings
        if list_of_lists:
            kwargs['ncol'] = len(handles[0]) # choose this for column length
        elif 'ncol' not in kwargs:
            kwargs['ncol'] = 3
        # Optionally change order
        # See: https://stackoverflow.com/q/10101141/4970632
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
    # The labelspacing/borderspacing will be exactly replicated, as if we were
    # using the original legend command
    # Means we also have to overhaul some settings
    else:
        legends = []
        overridden = []
        for override in ['loc','ncol','bbox_to_anchor','borderpad',
                         'borderaxespad','frameon','framealpha']:
            prop = kwargs.pop(override, None)
            if prop is not None:
                overridden.append(override)
        if overridden:
            print(f'Warning: Overriding legend properties {", ".join(prop for prop in overridden)}.')
        # Determine space we want sub-legend to occupy, as fraction of height
        # Don't normally save "height" and "width" of axes so keep here
        fontsize = kwargs.get('fontsize', None)     or rc['legend.fontsize']
        spacing  = kwargs.get('labelspacing', None) or rc['legend.labelspacing']
        interval = 1/len(handles) # split up axes
        interval = (((1 + spacing)*fontsize)/72) / \
                (ax.figure.get_figheight() * np.diff(ax._position.intervaly))
        # Iterate and draw
        if order=='F':
            raise NotImplementedError(f'When align=False, proplot vertically stacks successive single-row legends. Column-major (order="F") ordering is un-supported.')
        for h,hs in enumerate(handles):
            bbox = mtransforms.Bbox([[0,1-(h+1)*interval],[1,1-h*interval]])
            leg = super(BaseAxes, ax).legend(handles=hs, ncol=len(hs),
                loc='center',
                frameon=False,
                borderpad=0,
                bbox_to_anchor=bbox,
                **kwargs) # _format_legend is overriding original legend Method
            legends.append(leg)
        for l in legends[:-1]:
            ax.add_artist(l) # because matplotlib deletes previous ones

    # Properties for legends
    outline = {'linewidth': rc['axes.linewidth'],
               'edgecolor': rc['axes.edgecolor'],
               'facecolor': rc['axes.facecolor']}
    for leg in legends:
        leg.legendPatch.update(outline) # or get_frame()
        for obj in leg.legendHandles:
            obj.update(hsettings)
        # for t in leg.texts:
        #     t.update(tsettings) # or get_texts()
    return legends

def colorbar_factory(ax, mappable, values=None,
        orientation='horizontal', extend=None, extendlength=0.2,
        clabel=None, label=None,
        ctickminor=False, tickminor=None,
        cgrid=False, grid=None,
        ticklocation=None, cticklocation=None, ctickdir=None, tickdir='out',
        clocator=None, locator=None, cminorlocator=None, minorlocator=None,
        clocator_kw={}, locator_kw=None, cminorlocator_kw={}, minorlocator_kw=None,
        cformatter=None, formatter=None,
        cticklabels=None, ticklabels=None,
        norm=None, norm_kw={}, # normalizer to use when passing colors/lines
        **kwargs): #, settings=None):
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
    orientation : {'horizontal', 'vertical'}, optional
        The colorbar orientation.
    extend : {None, 'neither', 'both', 'min', 'max'}, optional
        Direction for drawing colorbar "extensions" (i.e. references to
        out-of-bounds data with a unique color). These are triangles by
        default. If ``None``, we try to use the ``extend`` attribute on the
        mappable object. If the attribute is unavailable, we use ``'neither'``.
    extendlength : float or str, optional
        The length of the colorbar "extensions" in *physical units*.
        If float, units are inches. If string,
        units are interpreted by `~proplot.utils.units`.

        This is handy if you have multiple colorbars in one figure.
        With the matplotlib API, it is really hard to get triangle
        sizes to match, because the `extendlength` units are *relative*.
    tickdir : {'out', 'in'}, optional
        Whether to draw tick marks pointing inside or outside.
    ctickdir, ticklocation, cticklocation
        Aliases of `tickdir`.
    label : None or str, optional
        The colorbar label.
    tickminor : bool, optional
        Whether to put minor ticks on the colorbar. Default is ``False``.
    grid : bool, optional
        Whether to draw "gridlines" (i.e. separators) between each level
        across the colorbar. Default is ``False``.
    clabel, ctickminor, cgrid
        Aliases for `label`, `tickminor`, `grid`.
    locator : None or locator spec, optional
        The colorbar tick mark positions. Passed to the
        `~proplot.axistools.Locator` constructor.
    locator_kw : dict-like, optional
        The locator settings. Passed to `~proplot.axistools.Locator`.
    minorlocator
        As with `locator`, but for the minor tick marks.
    minorlocator_kw
        As for `locator_kw`, but for the minor locator.
    clocator, cminorlocator, clocator_kw, cminorlocator_kw
        Aliases for `locator`, `minorlocator`, `locator_kw`, `minorlocator_kw`
    formatter : None or formatter spec, optional
        The tick label format. Passed to the `~proplot.axistools.Formatter`
        constructor.
    cformatter, ticklabels, cticklabels
        Aliases for `formatter`.
    norm : None or normalizer spec, optional
        Ignored if `mappable` has the ``norm`` attribute. The normalizer
        for converting `values` to colormap colors. Passed to the
        `~proplot.colortools.Norm` constructor.
    norm_kw : dict-like, optional
        The normalizer settings. Passed to `~proplot.colortools.Norm`.

    Other parameters
    ----------------
    **kwargs
        Passed to `~matplotlib.figure.Figure.colorbar`.

    See also
    --------
    `BaseAxes.colorbar`, `PanelAxes.colorbar`, `~matplotlib.figure.Figure.colorbar`,
    `~proplot.axistools.Locator`, `~proplot.axistools.Formatter`, `~proplot.colortools.Norm`
    """
    # Developer notes
    # * There are options on the colorbar object (cb.locator,
    #   cb.formatter with cb.update_ticks) and by passing kwargs (ticks=x,
    #   format=y) that allow user to not reference the underlying "axes"
    #   when fixing ticks. Don't use this functionality because not necessary
    #   for us and is missing many features, e.g. minorlocators/minorformatters.
    #   Also is different syntax.
    # * There is an INSANELY WEIRD problem with colorbars when simultaneously
    #   passing levels and norm object to a mappable; fixed by passing
    #   vmin/vmax INSTEAD OF levels.
    #   (see: https://stackoverflow.com/q/40116968/4970632).
    # * Problem is, often WANT levels instead of vmin/vmax, while simultaneously
    #   using a Normalize (for example) to determine colors between the levels
    #   (see: https://stackoverflow.com/q/42723538/4970632).
    # * Workaround is to make sure locators are in vmin/vmax range exclusively;
    #   cannot match/exceed values.
    # * The 'extend' kwarg is used for the case when you are manufacturing
    #   colorbar from list of colors or lines. Most of the time want 'neither'.
    # See comment under colorbar() method def for PanelAxes class. Will get
    # weird results if axes is a special BaseAxes.
    if isinstance(ax, BaseAxes):
        raise ValueError('The colorbar axes cannot be an instance of proplot.BaseAxes. Must be native matplotlib axes.Axes class.')
    # Parse flexible input
    clocator      = _default(locator, clocator)
    cgrid         = _default(grid, cgrid)
    ctickminor    = _default(tickminor, ctickminor)
    cminorlocator = _default(minorlocator, cminorlocator)
    cformatter    = _default(ticklabels, _default(cticklabels, _default(formatter, cformatter)))
    clabel        = _default(label, clabel)
    clocator_kw = _default(locator_kw, clocator_kw)
    cminorlocator_kw = _default(minorlocator_kw, cminorlocator_kw)

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
    csettings = {'cax':ax, 'orientation':orientation, 'use_gridspec':True, # use space afforded by entire axes
                 'spacing':'uniform', 'extend':extend, 'drawedges':cgrid} # this is default case unless mappable has special props
    csettings.update(kwargs)

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
    # manually wrap with cmap_features.
    if fromlines or fromcolors:
        cmap   = colortools.Colormap(colors)
        func = cmap_features(ax, ax.contourf)
        mappable = func([[0,0],[0,0]],
            values=np.array(values), cmap=cmap, extend='neither',
            norm=(norm or 'segmented')
            ) # workaround
        if clocator is None:
            nstep = 1 + len(values)//20
            clocator = values[::nstep]
    # By default, label the discretization levels (if there aren't too many)
    # Prefer centers (i.e. 'values') to edges (i.e. 'levels')
    if clocator is None:
        clocator = getattr(mappable, 'values', getattr(mappable, 'levels', None))
        if clocator is not None:
            step = 1 + len(clocator)//20
            clocator = clocator[::step]

    # Determine major formatters and major/minor tick locators
    # Can pass clocator/cminorlocator as the *jump values* between the mappables
    # vmin/vmax if desired
    fixed = None # so linter doesn't detect error in if i==1 block
    normfix = False # whether we need to modify the norm object
    locators = [] # put them here
    for i,(locator,locator_kw) in enumerate(zip((clocator,cminorlocator),(clocator_kw,cminorlocator_kw))):
        # Get the locator values
        # Need to use tick_values instead of accessing 'locs' attribute because
        # many locators don't have these attributes; require norm.vmin/vmax as input
        if i==1 and not ctickminor and locator is None: # means we never wanted minor ticks
            locators.append(axistools.Locator('null'))
            continue
        values = np.array(axistools.Locator(locator, **locator_kw).tick_values(mappable.norm.vmin, mappable.norm.vmax)) # get the current values
        # Modify ticks to work around mysterious error, and to prevent annoyance
        # where minor ticks extend beyond extendlength.
        # We need to figure out the numbers that will eventually be rendered to
        # solve the error, so we will always use a fixedlocator.
        values_min = np.where(values>=mappable.norm.vmin)[0]
        values_max = np.where(values<=mappable.norm.vmax)[0]
        if len(values_min)==0 or len(values_max)==0:
            locators.append(axistools.Locator('null'))
            continue
        values_min, values_max = values_min[0], values_max[-1]
        values = values[values_min:values_max+1]
        if values[0]==mappable.norm.vmin:
            normfix = True
        # Prevent annoying major/minor overlaps where one is slightly shifted left/right
        # Consider floating point weirdness too
        if i==1:
            eps = 1e-10
            values = [v for v in values if not any(o+eps >= v >= o-eps for o in fixed)]
        fixed = values # record as new variable
        locators.append(axistools.Locator(fixed)) # final locator object
    # Next the formatter
    cformatter = axistools.Formatter(cformatter)

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
    # NOTE: For whatever reason the only way to avoid bugs seems to be to pass
    # the major formatter/locator to colorbar commmand and directly edit the
    # minor locators/formatters; update_ticks after the fact ignores the major formatter
    # axis.set_major_locator(locators[0]) # does absolutely nothing
    # axis.set_major_formatter(cformatter)
    if orientation=='horizontal':
        axis = ax.xaxis
        scale = ax.figure.width*np.diff(getattr(ax.get_position(),'intervalx'))[0]
    else:
        axis = ax.yaxis
        scale = ax.figure.height*np.diff(getattr(ax.get_position(),'intervaly'))[0]
    extendlength = utils.units(extendlength)
    extendlength = extendlength/(scale - 2*extendlength)
    csettings.update({'extendfrac':extendlength}) # width of colorbar axes and stuff
    ticklocation = ticklocation or ctickdir or tickdir
    ticklocation = {'out':'outer', 'in':'inner'}.get(ticklocation, ticklocation)
    cb = ax.figure.colorbar(mappable,
            ticklocation=ticklocation,
            ticks=locators[0],
            format=cformatter,
            **csettings)

    # Make edges/dividers consistent with axis edges
    if cb.dividers is not None:
        cb.dividers.update(rc['grid'])

    # The minor locators and formatters
    # * The minor locator must be set with set_ticks after transforming an array
    #   using the mappable norm object; see: https://stackoverflow.com/a/20079644/4970632
    # * The set_minor_locator seems to be completely ignored depending on the colorbar
    #   in question, for whatever reason
    # * The major locator and formatter settings here are also not ideal since we'd have to
    #   update_ticks which might throw off the minor ticks again
    # WARNING: If functionality of BoundaryNorm is modified so data is transformed
    # by some linear transformation before getting binned, below may fail.
    # cb.minorticks_on() # alternative, but can't control the god damn spacing/set our own version
    # axis.set_minor_locator(locators[1]) # does absolutely nothing
    # WARNING: For some reason, pcolor mappables need to take *un-normalized
    # ticks* when set_ticks is called, while contourf mappables need to
    # take *normalized* data (verify by printing)
    minorvals = np.array(locators[1].tick_values(mappable.norm.vmin, mappable.norm.vmax))
    # axis.set_ticks(mappable.norm(majorvals), minor=False)
    if isinstance(mappable.norm, mcolors.BoundaryNorm): # including my own version
        vmin, vmax = mappable.norm.vmin, mappable.norm.vmax
        minorvals = (minorvals-vmin)/(vmax-vmin)
    elif hasattr(mappable, 'levels'):
        minorvals = mappable.norm(minorvals)
    axis.set_ticks(minorvals, minor=True)
    axis.set_minor_formatter(mticker.NullFormatter()) # to make sure

    # Set up the label
    if clabel is not None:
        axis.label.update({'text':clabel})

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
        print('Warning: Performing manual alpha-blending for colorbar solids.')
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
        # cb.solids.set_cmap()

    # Fix pesky white lines between levels + misalignment with border due to rasterized blocks
    if cb.solids:
        cb.solids.set_linewidth(0.2) # something small
        cb.solids.set_edgecolor('face')
        cb.solids.set_rasterized(False)
    return cb

