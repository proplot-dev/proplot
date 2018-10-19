#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Imports, all
# Note that because of how matplotlib imports work (a module is only imported
# for a first time; import statements elsewhere in script during process just point
# to the already loaded module), settings rcParams here will change properties
# everywhere else until you start new instance
# TODO: Add inner panels as attributes to axes, instead of as seperate
# return arguments
# TODO: Standardize the order that panels are delivered (e.g. in a list
# [None, axh, None, axh]), then *fix axis sharing across multiple panels*!
# TODO: Note severe change; twiny now means "share the x-axis but
# now make two y-axes"; this makes way more sense to me than default behavior
# perhaps should undo this
# TODO: Add feature for multiple (stacked) bottompanels and rightpanels, useful for having
# different colorbars on same figure, or legend and colorbar
# TODO: Add feature to draw axes in *row-major* or *column-major* mode
# TODO: Add options to axes formatter to *disable* the 'quiet' formatting steps
# that normally take place, e.g. colorizing the label and title, formatters, setting
# linewidths, etc. This way can, for example, change color of spines only, then call
# format again with new rc() settings and xlabels/ylabels/stuff are different
# colors. In general *allow overwriting rc() settings* perhaps. Should make
# globals a *convenience feature* but *not necessary* to change certain properties, 
# everything should *also* be accessible with format function.
#------------------------------------------------------------------------------#
# Decorators used a lot here; below is very simple example that demonstrates
# how simple wrapper decorators work
# def decorator1(func):
#     def wrapper():
#         print('decorator 1 called')
#         func()
#         print('decorator 1 finished')
#     return wrapper
# def decorator2(func):
#     def wrapper():
#         print('decorator 2 called')
#         func()
#         print('decorator 2 finished')
#     return wrapper
# @decorator1
# @decorator2
# def hello():
#     print('hello world!')
# hello()
#------------------------------------------------------------------------------
# Global module dependencies
# Recommended using functools.wraps from comment:
# https://stackoverflow.com/a/739665/4970632
# This tool preserve __name__ metadata.
# from .basics import Dict, arange
import os
import numpy as np
import warnings
import time
# import io
# from contextlib import redirect_stdout
from matplotlib.cbook import mplDeprecation
from string import ascii_lowercase
from types import FunctionType
from functools import wraps
from inspect import cleandoc
import matplotlib.figure as mfigure
import matplotlib.axes as maxes
import matplotlib.projections as mprojections
import matplotlib.path as mpath
import matplotlib.contour as mcontour
import matplotlib.patheffects as mpatheffects
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import matplotlib.text as mtext
import matplotlib.ticker as mticker
import matplotlib.gridspec as mgridspec
import matplotlib.artist as martist
import matplotlib.transforms as mtransforms
import matplotlib.collections as mcollections
import matplotlib.cm as mcm
import matplotlib.pyplot as plt
# Local modules, projection sand formatters and stuff
from .rcmod import rc
from .axis import Locator, Formatter # default axis norm and formatter
from .proj import Aitoff, Hammer, KavrayskiyVII, WinkelTripel
from . import colortools
from . import utils
# Filter warnings, seems to be necessary before drawing stuff for first time,
# otherwise this has no effect (e.g. if you stick it in a function)
warnings.filterwarnings('ignore', category=mplDeprecation)
# Optionally import mapping toolboxes
# Main conda distro says they are incompatible, so make sure not required!
try:
    import mpl_toolkits.basemap as mbasemap
except ModuleNotFoundError:
    print('Warning: Basemap is not available.')
try: # crs stands for "coordinate reference system", leading c is "cartopy"
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.geoaxes import GeoAxes
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    PlateCarree = ccrs.PlateCarree() # global variable
except ModuleNotFoundError:
    GeoAxes = PlateCarree = object
    print('Warning: cartopy is not available.')

#------------------------------------------------------------------------------#
# Create matplotlib.Path objects
#------------------------------------------------------------------------------#
def circle(N):
    """
    Draw a circle.
    """
    theta = np.linspace(0, 2*np.pi, N)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    return mpath.Path(verts * radius + center)

#------------------------------------------------------------------------------#
# Misc tools
#------------------------------------------------------------------------------#
def fancy_decorator(decorator):
    """
    Normally to make a decorator that accepts arguments, you have to create
    3 nested function definitions. This abstracts that away -- if you decorate
    your decorator-function declaration with this, the decorator will now accept arguments.
    See: https://stackoverflow.com/a/1594484/4970632
    """
    @wraps(decorator)
    def decorator_maker(*args, **kwargs):
        def decorator_wrapper(func):
            return decorator(func, *args, **kwargs)
        return decorator_wrapper
    return decorator_maker

def timer(func):
    """
    A decorator that prints the time a function takes to execute.
    See: https://stackoverflow.com/a/1594484/4970632
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        t = time.clock()
        res = func(*args, **kwargs)
        print(f'{func.__name__} time: {time.clock()-t}s')
        return res
    return wrapper

def logger(func):
    """
    A decorator that logs the activity of the script (it actually just prints it,
    but it could be logging!)
    See: https://stackoverflow.com/a/1594484/4970632
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        res = func(*args, **kwargs)
        print(f'{func.__name__} called with: {args} {kwargs}')
        return res
    return wrapper

def counter(func):
    """
    A decorator that counts and prints the number of times a function
    has been executed.
    See: https://stackoverflow.com/a/1594484/4970632
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        wrapper.count = wrapper.count + 1
        res = func(*args, **kwargs)
        print(f'{func.__name__} has been used: {wrapper.count}x')
        return res
    wrapper.count = 0
    return wrapper

#------------------------------------------------------------------------------
# Helper functions for plot overrides
# List of stuff in pcolor/contourf that need to be fixed:
#   * White lines between the edges; cover them by changing edgecolors to 'face'.
#   * Determination of whether we are using graticule edges/centers; not sure
#       what default behavior is but harder to debug. My wrapper is nicer.
#   * Pcolor can't take an extend argument, and colorbar can take an extend argument
#       but it is ignored when the mappable is a contourf. Make our pcolor wrapper
#       add an "extend" attribute on the mappable that our colorbar wrapper detects.
#   * Extend used in contourf causes color-change between in-range values and
#       out-of-range values, while extend used in colorbar on pcolor has no such
#       color change. Standardize by messing with the colormap.
#------------------------------------------------------------------------------
def check_centers(func):
    """
    Check shape of arguments passed to contour, and fix result.
    """
    @wraps(func)
    def wrapper(self, x, y, *args, **kwargs):
        # Checks whether sizes match up, checks whether graticule was input
        x, y = np.array(x), np.array(y)
        args = [np.array(Z) for Z in args]
        xlen, ylen = x.shape[0], y.shape[-1]
        for Z in args:
            if Z.shape[0]==xlen-1 and Z.shape[1]==ylen-1:
                x, y = (x[1:]+x[:-1])/2, (y[1:]+y[:-1])/2 # get centers, given edges
            elif Z.shape[0]!=xlen or Z.shape[1]!=ylen:
                raise ValueError(f'X ({"x".join(str(i) for i in x.shape)}) '
                        f'and Y ({"x".join(str(i) for i in y.shape)}) must correspond to '
                        f'nrows ({Z.shape[0]}) and ncolumns ({Z.shape[1]}) of Z, or its borders.')
        args = [Z.T for Z in args]
        return func(self, x, y, *args, **kwargs)
    return wrapper

def check_edges(func):
    """
    Check shape of arguments passed to pcolor, and fix result.
    """
    @wraps(func)
    def wrapper(self, x, y, *args, **kwargs):
        # Checks that sizes match up, checks whether graticule was input
        x, y = np.array(x), np.array(y)
        args = [np.array(Z) for Z in args]
        xlen, ylen = x.shape[0], y.shape[-1]
        for Z in args:
            if Z.shape[0]==xlen and Z.shape[1]==ylen:
                x, y = utils.edges(x), utils.edges(y)
            elif Z.shape[0]!=xlen-1 or Z.shape[1]!=ylen-1:
                raise ValueError(f'X ({"x".join(str(i) for i in x.shape)}) '
                        f'and Y ({"x".join(str(i) for i in y.shape)}) must correspond to '
                        f'nrows ({Z.shape[0]}) and ncolumns ({Z.shape[1]}) of Z, or its borders.')
        args = [Z.T for Z in args]
        return func(self, x, y, *args, **kwargs)
    return wrapper

def cmap_features(func):
    """
    Manage output of contour and pcolor functions.
    New features:
        * Create new colormaps on the fly, and merge arbitrary named
          or created colormaps.
        * Always use full range of colormap, whether you are extending
          max, min, neither, or both. For the first three, will reconstruct
          colormap so 'out-of-bounds' have same color as edge colors
          from 'in-bounds' region.
    Also see: https://stackoverflow.com/a/48614231/4970632
    """
    @wraps(func)
    def wrapper(self, *args, cmap=None, norm=None, levels=None, extend='neither', extremes=True, **kwargs):
        # Specify normalizer
        name = func.__name__
        default = ('continuous' if 'contour' in name else 'discrete')
        norm = colortools.Norm(norm, default=default, levels=levels)
        if levels is None:
            levels = getattr(norm, 'levels', None)
        # Specify colormap
        if any(isinstance(cmap,t) for t in (str,dict,mcolors.Colormap)):
            cmap = cmap, # make a tuple
        cmap_extend = 'neither' if not extremes else extend
        cmap = colortools.Colormap(*cmap, levels=levels, extend=cmap_extend) # pass as arguments to generalized colormap constructor
        # Call function with special args removed
        kwextra = dict()
        if 'contour' in name:
            kwextra = dict(levels=levels, extend=extend)
        result = func(self, *args, cmap=cmap, norm=norm, **kwextra, **kwargs)
        if 'pcolor' in name:
            result.extend = extend
        # Optionally use same color for data in 'edge bins' as 'out of bounds' data
        # NOTE: This shouldn't mess up normal contour() calls because if we
        # specify color=<color>, should override cmap instance anyway
        # Fix white lines between filled contours/mesh
        linewidth = 0.2
        if name=='contourf':
            for contour in result.collections:
                contour.set_edgecolor('face')
                contour.set_linewidth(linewidth)
        elif not utils.isiterable(result):
            result.set_edgecolor('face')
            result.set_linewidth(linewidth) # seems to do the trick, without dots in corner being visible
        return result
    return wrapper

def cycle_features(func):
    """
    Allow specification of color cycler at plot-time. Will simply set the axes
    property cycler, and if it differs from user input, update it.
    See: https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/axes/_base.py
    The set_prop_cycle command modifies underlying _get_lines and _get_patches_for_fill.
    """
    @wraps(func)
    def wrapper(self, *args, cycle=None, **kwargs):
        # Determine and temporarily set cycler
        if cycle is not None:
            if not utils.isiterable(cycle) or type(cycle) is str:
                cycle = cycle,
            cycle = colortools.Cycle(*cycle)
            self.set_prop_cycle(color=cycle)
        return func(self, *args, **kwargs)
    return wrapper

#------------------------------------------------------------------------------#
# Helper functions for basemap and cartopy plot overrides
# NOTE: These wrappers should be invoked *after* check_centers and check_edges,
# which perform basic shape checking and permute the data array, so the data
# will now be y by x (or lat by lon) instead of lon by lat.
#------------------------------------------------------------------------------#
# Normally we *cannot* modify the underlying *axes* pcolormesh etc. because this
# this will cause basemap's self.m.pcolormesh etc. to use my *custom* version and
# cause a suite of weird errors. Prevent this recursion with the below decorator.
def norecurse(func):
    """
    Decorator to prevent recursion in Basemap method overrides.
    This will call the function 'func' on the first run (which does some stuff,
    then calls the basemap method), and will call the super-class method of the
    same name on the second run. See: https://stackoverflow.com/a/37675810/4970632
    """
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        name = getattr(func, '__name__')
        if self.recurred:
            self.recurred = False
            call = getattr(super(BasemapAxes,self), name) # don't call func again, now we want to call the parent function
            obj = () # this time 'self' is repeated in position args[0]
        else:
            self.recurred = True
            call = func
            obj = self,
        result = call(*obj, *args, **kwargs)
        self.recurred = False # cleanup
        return result
    return wrapper

def gridfix_cartopy(func):
    """
    Apply default transform and fix discontinuities in grid.
    """
    @wraps(func)
    def wrapper(self, lon, lat, Z, transform=PlateCarree, **kwargs):
        # The todo list for cartopy is much shorter, as we don't have to worry
        # about longitudes going from 0 to 360, -180 to 180, etc.; projection handles all that
        # 1) Fix holes over poles by *interpolating* there (equivalent to
        # simple mean of highest/lowest latitude points)
        Z_south = np.repeat(Z[0,:].mean(),  Z.shape[1])[None,:]
        Z_north = np.repeat(Z[-1,:].mean(), Z.shape[1])[None,:]
        lat = np.concatenate(([-90], lat, [90]))
        Z = np.concatenate((Z_south, Z, Z_north), axis=0)
        # 2) Fix seams at map boundary; this time the fancy projection will
        # handle issues with pcolor seams that we'd have with basemap, but still
        # have to ensure *circular* coverage if doing e.g. contourf
        lon = np.array((*lon, lon[0]+360)) # make longitudes circular
        Z = np.concatenate((Z, Z[:,:1]), axis=1) # make data circular
        # Call function
        # _ = io.StringIO() # message has a bunch of unnecessary newlines; will modify it
        # with redirect_stdout(_):
        # with warnings.catch_warnings():
        #     warnings.simplefilter('ignore')
        result = func(self, lon, lat, Z, transform=transform, **kwargs)
        # Call function
        return result
    return wrapper

def gridfix_basemap(func):
    """
    Interpret coordinates and fix discontinuities in grid.
    """
    @wraps(func)
    def wrapper(self, lon, lat, Z, **kwargs):
        # Raise errors
        # print('lon', lon, 'lat', lat, 'Z', Z)
        lonmin, lonmax = self.m.lonmin, self.m.lonmax
        if lon.max()>lon.min()+360:
            print(lon.max()-lon.min())
            raise ValueError(f'Longitudes span {lon.min()} to {lon.max()}. Can only span 360 degrees at most.')
        if lon.min()<-360 or lon.max()>360:
            raise ValueError(f'Longitudes span {lon.min()} to {lon.max()}. Must fall in range [-360, 360].')
        if lonmin<-360 or lonmin>0:
            print(f'Warning: Minimum longitude is {lonmin}, not in range [-360,0].')
            # raise ValueError('Minimum longitude must fall in range [-360, 0].')
        # 1) Establish 360-degree range
        lon -= 720
        while True:
            filter_ = lon<lonmin
            if filter_.sum()==0:
                break
            lon[filter_] += 360
        # 2) Roll, accounting for whether ends are identical
        # If go from 0,1,-->,359,0 (borders), returns id of first zero
        roll = -np.argmin(lon) # always returns *first* value
        if lon[0]==lon[-1]:
            lon = np.roll(lon[:-1], roll)
            lon = np.append(lon, lon[0]+360)
        else:
            lon = np.roll(lon, roll)
        Z = np.roll(Z, roll, axis=1)
        # 3) Roll in same direction some more, if some points on right-edge
        # extend more than 360 above the minimum longitude; THEY should be the
        # ones on west/left-hand-side of map
        lonroll = np.where(lon>lonmin+360)[0] # tuple of ids
        if lonroll: # non-empty
            roll = lon.size-min(lonroll) # e.g. if 10 lons, lonmax id is 9, we want to roll once
            lon = np.roll(lon, roll) # need to roll foreward
            Z = np.roll(Z, roll, axis=1) # roll again
            lon[:roll] -= 360 # retains monotonicity
        # 4) Set NaN where data not in range lonmin, lonmax
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
        # 5) Fix holes over poles by interpolating there (equivalent to
        # simple mean of highest/lowest latitude points)
        if lon.size==Z.shape[1]: # have centers, not grid cell edges
            Z_south = np.repeat(Z[0,:].mean(),  Z.shape[1])[None,:]
            Z_north = np.repeat(Z[-1,:].mean(), Z.shape[1])[None,:]
            lat = np.concatenate(([-90], lat, [90]))
            Z = np.concatenate((Z_south, Z, Z_north), axis=0)
        # 6) Fix seams at map boundary; 3 scenarios here:
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
        return func(self, x, y, Z, **kwargs)
    return wrapper

#------------------------------------------------------------------------------
# Custom classes class
#------------------------------------------------------------------------------
def docstring_fix(child):
    """
    Decorator function for appending documentation from overridden method
    onto the overriding method docstring.
    Adapted from: https://stackoverflow.com/a/8101598/4970632
    """
    for name,chfunc in vars(child).items(): # returns __dict__ object
        if not isinstance(chfunc, FunctionType):
            continue
        for parent in getattr(child, '__bases__', ()):
            parfunc = getattr(parent, name, None)
            if not getattr(parfunc, '__doc__', None):
                continue
            if not getattr(chfunc, '__doc__', None):
                chfunc.__doc__ = '' # in case it's None
            cmessage = f'Full name: {parfunc.__qualname__}()'
            pmessage = f'Parent method (documentation below): {chfunc.__qualname__}()'
            chfunc.__doc__ = f'\n{cmessage}\n{cleandoc(chfunc.__doc__)}\n\n{pmessage}\n{cleandoc(parfunc.__doc__)}'
            break # only do this for the first parent class
    return child

@docstring_fix
class Figure(mfigure.Figure):
    """
    Subclass of the mfigure.Figure class, with lots of special formatting
    options. Can be called by using pyplot.figure(FigureClass=Figure) kwargument
    in my subplots function.
    """
    def __init__(self, left=None, bottom=None, right=None, top=None,
        bwidth=None, bspace=None, rwidth=None, rspace=None,
        width=None, height=None, gridspec=None,
        **kwargs):
        """
        Add a few extra attributes, accessed by other methods.
        """
        self.left     = left
        self.bottom   = bottom
        self.right    = right
        self.top      = top
        self.bwidth   = bwidth
        self.bspace   = bspace
        self.rwidth   = rwidth
        self.rspace   = rspace
        self.width    = width
        self.height   = height
        self.gridspec = gridspec
        super().__init__(**kwargs) # python 3 only

    def panel_factory(self, subspec, width, height,
            whichpanels=None, hspace=None, wspace=None, hwidth=None, wwidth=None,
            **kwargs):
        """
        Helper function for creating panelled axes
        Will take in some arguments from parent function so don't have to pass them
        every time.
        """
        # Format the plot to have tiny spaces between subplots and narrow panels by default
        if whichpanels is None:
            whichpanels = 'b' # default is a bottom panel, cuz why not
        if hspace is None:
            hspace = 0.15
        if wspace is None:
            wspace = 0.15
        if hwidth is None:
            hwidth = 0.5
        if wwidth is None:
            wwidth = 0.5
        if any(s.lower() not in 'tblr' for s in whichpanels) or whichpanels=='':
            raise ValueError('Whichpanels argument can contain characters l (left), r (right), b (bottom), or t (top).')
        # Determine number of rows/columns, position of the
        # main/primary axes, and potential position of corners that
        # we want to ignore for now
        nrows_i, ncols_i = 1, 1
        for s in whichpanels:
            if s in ('b','t'):
                nrows_i += 1
            if s in ('l','r'):
                ncols_i += 1
        empty_pos = []
        main_pos = [0,0]
        if 't' in whichpanels:
            main_pos[0] += 1
        if 'l' in whichpanels:
            main_pos[1] += 1
        corners = {'tl':(0,0),             'tr':(0,main_pos[1]+1),
                'bl':(main_pos[0]+1,0), 'br':(main_pos[0]+1,main_pos[1]+1)}
        for corner,position in corners.items():
            if all(s in whichpanels for s in corner):
                empty_pos.append(position)
        # Fix wspace/hspace in inches, using the Bbox from get_postition
        # on the subspec object to determine physical width of axes to be created
        # * Consider writing some convenience funcs to automate this unit conversion
        bbox_i = subspec.get_position(self) # valid since axes not drawn yet
        if hspace is not None:
            height_i = np.diff(bbox_i.intervaly)[0]*height - hspace*(nrows_i-1)
            hspace = hspace/(height_i/nrows_i)
        if wspace is not None:
            width_i = np.diff(bbox_i.intervalx)[0]*width - wspace*(ncols_i-1)
            wspace = wspace/(width_i/ncols_i)
        # Figure out hratios/wratios
        # Will enforce (main_width + panel_width)/total_width = 1
        wwidth_ratios = [width_i-wwidth*(ncols_i-1)]*ncols_i
        if wwidth_ratios[0]<0:
            raise ValueError(f'Panel wwidth is too large. Must be less than {width_i/(nrows_i-1):.3f}.')
        for i in range(ncols_i):
            if i!=main_pos[1]: # this is a panel entry
                wwidth_ratios[i] = wwidth
        hwidth_ratios = [height_i-hwidth*(nrows_i-1)]*nrows_i
        if hwidth_ratios[0]<0:
            raise ValueError(f'Panel hwidth is too large. Must be less than {height_i/(ncols_i-1):.3f}.')
        for i in range(nrows_i):
            if i!=main_pos[0]: # this is a panel entry
                hwidth_ratios[i] = hwidth
        # Determine axes to be shared
        if ncols_i==1:
            ybase = ()
            yshare = []
        elif 'l' in whichpanels:
            ybase = (main_pos[0], main_pos[1]-1) # the base axes
            yshare = [(main_pos[0], main_pos[1]+i) for i in range(ncols_i-1)]
        else:
            ybase = (main_pos[0], main_pos[1]) # the base axes
            yshare = [(main_pos[0], main_pos[1]+i+1) for i in range(ncols_i-1)]
        if nrows_i==1:
            xbase = ()
            xshare = []
        elif 'b' in whichpanels:
            xbase = (main_pos[0]+1, main_pos[1]) # the base axes
            xshare = [(main_pos[0]-i, main_pos[1]) for i in range(nrows_i-1)]
        else:
            xbase = (main_pos[0], main_pos[1]) # the base axes
            xshare = [(main_pos[0]-i-1, main_pos[1]) for i in range(nrows_i-1)]
        # Create subplotspec and draw the axes
        # Will create axes in order of rows/columns so that the "base" axes
        # are always built before the axes to be "shared" with them
        panels = []
        gs_i = mgridspec.GridSpecFromSubplotSpec(
                nrows         = nrows_i,
                ncols         = ncols_i,
                subplot_spec  = subspec,
                wspace        = wspace,
                hspace        = hspace,
                width_ratios  = wwidth_ratios,
                height_ratios = hwidth_ratios,
                )
        ax_ybase, ax_xbase = None, None
        if 'sharex' in kwargs:
            kwargs.pop('sharex')
        if 'sharey' in kwargs:
            kwargs.pop('sharey')
        # Create panel axes
        for r in range(nrows_i)[::-1]: # order is bottom-top
            for c in range(ncols_i):   # order is left-right
                if (r,c) in empty_pos:
                    continue
                main = (r==main_pos[0] and c==main_pos[1])
                # Add the subplot first
                # Also, we optionally allow for *projection axes* with *panels*, which
                # would make sense for pseudo-cylindrical or rectangular projections
                if main:
                    ax_kwargs = kwargs
                else:
                    ax_kwargs = kwargs.copy()
                    ax_kwargs['projection'] = 'xy'
                    ax_kwargs.pop('transform', None)
                ax = self.add_subplot(gs_i[r,c], **ax_kwargs) # make this my custom XY axes
                if r==main_pos[0] and c==main_pos[1]:
                    axmain = ax
                else:
                    panels += [ax]
                # Configure shared axes
                if (r,c) == ybase:
                    ax_ybase = ax
                if (r,c) == xbase:
                    ax_xbase = ax
                if (r,c) in yshare:
                    if ax_ybase is None:
                        raise ValueError("Axes with base for x-axis not drawn.")
                    ax._sharey = ax_ybase
                    for t in ax.yaxis.get_ticklabels(): t.set_visible(False)
                    ax.yaxis.label.set_visible(False)
                if (r,c) in xshare:
                    if ax_xbase is None:
                        raise ValueError("Axes with base for x-axis not drawn.")
                    ax._sharex = ax_xbase
                    for t in ax.xaxis.get_ticklabels(): t.set_visible(False)
                    ax.xaxis.label.set_visible(False)
        # Add panels as children of main axes, and return main axes
        if len(panels)!=len(whichpanels):
            raise ValueError(f'Created {len(panels)} panels but {len(whichpanels)} panels ("{whichpanels}") were requested.')
        translate = {'b':'bottom', 't':'top', 'l':'left', 'r':'right'}
        # axmain.panel = {} # alternative
        for panel,string in zip(panels,whichpanels):
            # axmain.panel[string] = panel # access panel['b'], panel['t'], et cetera
            setattr(axmain, f'{translate[string]}panel', panel) # access 'bottompanel', 'toppanel', 'leftpanel', and 'rightpanel'
        return axmain

    # @timer
    # def draw(self, *args, **kwargs):
    #     """
    #     Wrapper so we can time stuff.
    #     """
    #     return super().draw(*args, **kwargs)

    @timer
    def save(self, filename, silent=False, pad=None, **kwargs):
        """
        Echo some helpful information before saving.
        Note that the gridspec object must be updated before figure is printed to screen
        in interactive environment; will fail to update after that. Seems to be glitch,
        should open thread on GitHub.
        Note to color axes patches, you will have to explicitly pass the
        transparent=False kwarg.
        """
        # Some kwarg translations, to pass to savefig
        # Note almost never want an 'edgecolor' for the figure frame, though it can be set
        if 'alpha' in kwargs:
            kwargs['transparent'] = not bool(kwargs.pop('alpha')) # 1 is non-transparent
        if 'color' in kwargs:
            kwargs['facecolor'] = kwargs.pop('color') # the color
        # Get bounding box that encompasses *all artists*, compare to bounding
        # box used for saving *figure*
        obbox = self.bbox_inches # original bbox
        bbox = self.get_tightbbox(self.canvas.get_renderer())
        ox, oy, x, y = obbox.intervalx, obbox.intervaly, bbox.intervalx, bbox.intervaly
        x1, y1, x2, y2 = x[0], y[0], ox[1]-x[1], oy[1]-y[1] # deltas
        width, height = ox[1], oy[1] # desired save-width
        # Echo some information
        if not silent:
            x1fix, y2fix = self.left-x1, self.top-y2
            if self.rightpanel is not None:
                x2fix = self.rspace - x2
            else:
                x2fix = self.right  - x2
            if self.bottompanel is not None:
                y1fix = self.bspace - y1
            else:
                y1fix = self.bottom - y1
            formatter = lambda x,y: f'{x:.2f} ({-y:+.2f})'
            print(f'Graphics bounded by L{formatter(x1fix,x1)} R{formatter(x2fix,x2)} '
                  f'B{formatter(y1fix,y1)} T{formatter(y2fix,y2)}.')
        # Perform gridspec adjustment with some padding
        if pad is not None:
            self.gridspec.left   -= (x1-pad)/width
            self.gridspec.bottom -= (y1-pad)/height
            self.gridspec.right  += (x2-pad)/width
            self.gridspec.top    += (y2-pad)/height
            self.gridspec.update()
        # Finally, save
        if not silent:
            print(f'Saving to "{filename}".')
        self.savefig(os.path.expanduser(filename), **kwargs) # specify DPI for embedded raster objects

@docstring_fix
class Axes(maxes.Axes):
    """
    Subclass the default Axes class. Then register it as the 'base' projection,
    and you will get a subclass Subplot by calling fig.add_subplot(projection='pubplot').
    Notes:
    * You cannot subclass SubplotBase directly, should only be done with
      maxes.subplot_class_factory, which is called automatically when using add_subplot.
    * Cartopy projections should use same methods as for ordinary 'cartesian'
      plot, so we put a bunch of definition overrides in here.
    """
    # Initial stuff
    name = 'base'
    def __init__(self, *args, **kwargs):
        """
        Initialize and add some custom attributes.
        """
        super().__init__(*args, **kwargs)
        self.number = None
        self.width  = np.diff(self._position.intervalx)*self.figure.width # position is in figure units
        self.height = np.diff(self._position.intervaly)*self.figure.height

    # Helpful override, to prevent annoying IndexError when user accidentally
    # tries to index single axes instead of list of axes
    def __getitem__(self, key):
        if key==0:
            return self
        else:
            raise IndexError('Figure contains only one subplot.')

    # New convenience feature
    def format(self,
        hatch=None, color=None, # control figure/axes background; hatch just applies to axes
        suptitle=None, suptitlepos=None, title=None, titlepos=None, titlepad=0.1, titledict={},
        abc=False, abcpos=None, abcformat='', abcpad=0.1, abcdict={},
        ):
        """
        Function for formatting axes of all kinds; some arguments are only relevant to special axes, like 
        colorbar or basemap axes. By default, simply applies the dictionary values from settings() above,
        but can supply many kwargs to further modify things.
        TODO:
        * Add options for datetime handling; note possible date axes handles are TimeStamp (pandas),
          np.datetime64, DateTimeIndex; can fix with fig.autofmt_xdate() or manually set options; uses
          ax.is_last_row() or ax.is_first_column(), so should look that up.
        * Problem is there is no autofmt_ydate(), so really should implement my own
          version of this.
        """
        # Create figure title
        fig = self.figure # the figure
        if suptitle is not None and getattr(self,'number',0)==1:
            xpos = fig.left/fig.width + .5*(fig.width - fig.left - fig.right)/fig.width
            ypos = (self.title._transform.transform(self.title.get_position())[1]/fig.dpi \
                    + self.title._linespacing*self.title.get_size()/72)/fig.height
            if not title: # just place along axes if no title here
                ypos -= (self.title._linespacing*self.title.get_size()/72)/fig.height
            # default linespacing is 1.2; it has no getter, only a setter; see:
            # https://matplotlib.org/api/text_api.html#matplotlib.text.Text.set_linespacing
            # print(self.title.get_size()/72, ypos/fig.height)
            fig.suptitle = self.text(xpos,ypos,suptitle,
                    transform=fig.transFigure) # title
            fig.suptitle.update({'ha':'center', 'va':'baseline', **rc('suptitle')})
            if suptitlepos=='title': # not elevated
                ypos = self.title._transform.transform(self.title.get_position())[1]/fig.dpi/fig.height
                fig.suptitle.update({'position':(xpos,ypos)})
            elif isinstance(suptitlepos,str):
                raise ValueError(f"Unknown suptitle position: {suptitlepos}.")
            elif suptitlepos is not None:
                fig.suptitle.update({'position':suptitlepos})

        # Create axes title
        # Input needs to be emptys string
        self.title.update({**rc('title'), 'text':title or ''})
        self.title.update(titledict)
        if titlepos=='left':
            self.title.update({'position':(0,1), 'ha':'left'})
        elif titlepos=='right':
            self.title.update({'position':(1,1), 'ha':'right'})
        elif titlepos=='inside':
            self.title.update({'position':(0.5,1-titlepad/self.height),
                'transform':self.transAxes, 'va':'top'})
        elif isinstance(titlepos,str):
            raise ValueError(f'Unknown title position: {titlepos}.')
        elif titlepos is not None: 
            self.title.update({'position':titlepos, 'ha':'center', 'va':'center'})

        # Create axes numbering
        if self.number is not None and abc:
            abcedges = abcformat.split('a')
            self.abc = self.text(0, 1, abcedges[0] + ascii_lowercase[self.number-1] + abcedges[-1],
                    transform=self.title._transform, **abcdict) # optionally include paren
            self.abc.update({'ha':'left', 'va':'baseline', **rc('abc')})
            if abcpos=='inside':
                self.abc.update({'position':(abcpad/self.width, 1-abcpad/self.height),
                    'transform':self.transAxes, 'ha':'left', 'va':'top'})
            elif isinstance(abcpos,str):
                raise ValueError(f'Unknown abc position: {abcpos}.')
            elif abcpos is not None:
                self.abc.update({'position':abcpos, 'ha':'left', 'va':'top'})

        # Color setup, optional hatching in background of axes
        # You should control transparency by passing transparent=True or False
        # to the savefig command
        if color is not None:
            self.patch.set_color(color)
        self.patch.set_zorder(-1)
        self.patch.set_clip_on(False)
        if hatch: # non-empty string or not none
            self.fill_between([0,1], 0, 1, hatch=hatch, zorder=0, # put in back
                facecolor='none', edgecolor='k', transform=self.transAxes)
        return self

    # Create legend creation method
    def legend(self, *args, **kwargs):
        """
        Call custom legend() function.
        """
        return legend_factory(self, *args, **kwargs)

    # Fill entire axes with colorbar
    def colorbar(self, *args, **kwargs):
        """
        Call colorbar() function.
        """
        return colorbar_factory(self, *args, **kwargs)

    # Fancy wrappers
    def text(self, x, y, text, transform='axes', fancy=False, black=True,
            linewidth=2, lw=None, **kwarg): # linewidth is for the border
        """
        Wrapper around original text method.
        If requested, can return specially created black-on-white text.
        So far no other features implemented.
        """
        linewidth = lw or linewidth
        if type(transform) is not str:
            pass # leave alone
            # raise ValueError("Just name the transform with string \"axes\" or \"data\".")
        elif transform=='figure':
            transform = self.figure.transFigure
        elif transform=='axes':
            transform = self.transAxes
        elif transform=='data':
            transform = self.transData
        else:
            raise ValueError("Unknown transform name. Use string \"axes\" or \"data\".")
        t = super().text(x, y, text, transform=transform, **kwarg)
        if fancy:
            fcolor, bcolor = 'wk'[black], 'kw'[black]
            t.update({'color':fcolor, 'zorder':1e10, # have to update after-the-fact for path effects
                'path_effects': [mpatheffects.Stroke(linewidth=linewidth, foreground=bcolor), mpatheffects.Normal()]})
            # t.update({'size':11, 'zorder':1e10,
            #     'path_effects':[mpatheffects.PathPatchEffect(edgecolor=bcolor,linewidth=.6,facecolor=fcolor)]})
        return t

    @cycle_features
    def scatter(self, *args, **kwargs):
        """
        Just add some more consistent keyword argument options.
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
        aliases = {'c':   ['markercolor', 'color'],
                   's':   ['markersize', 'size'],
            'linewidths': ['lw','linewidth','markeredgewidth', 'markeredgewidths'],
            'edgecolors': ['markeredgecolor', 'markeredgecolors']}
        for name,options in aliases.items():
            for option in options:
                if option in kwargs:
                    kwargs[name] = kwargs.pop(option)
        return super().scatter(*args, **kwargs)

    @cmap_features
    def cmapline(*args, cmap=None, values=None, norm=None, bins=True, nsample=1):
        """
        Create lines with colormap.
        See: https://matplotlib.org/gallery/lines_bars_and_markers/multicolored_line.html
        Will manage input more strictly, this is harder to generalize.
        """
        # First error check
        if len(args) not in (1,2):
            raise ValueError(f'Function requires 1-2 arguments, got {len(args)}.')
        if cmap is None or values is None:
            raise ValueError(f'Function requires "cmap" and "values" kwargs.')
        y = np.array(args[-1]).squeeze()
        x = np.arange(y.shape[-1]) if len(args)==1 else np.array(args[0]).squeeze()
        values = np.array(values).squeeze()
        if x.ndim!=1 or y.ndim!=1 or values.ndim!=1:
            raise ValueError(f'Input x ({x.ndim}D), y ({y.ndim}D), and values ({values.ndim}D) must be 1-dimensional.')
        # Next draw the line
        # Interpolate values to optionally allow for smooth gradations between
        # values (bins=True) or color switchover halfway between points (bins=False)
        newx, newy, newvalues = [], [], []
        edges = utils.edges(values)
        if bins:
            norm = colortools.DiscreteNorm(edges)
        else:
            norm = colortools.ContinuousNorm(edges)
        for j in range(y.shape[0]-1):
            newx.extend(np.linspace(x[j], x[j+1], nsample+2))
            newy.extend(np.linspace(y[j], y[j+1], nsample+2))
            interp = np.linspace(np.asscalar(norm(values[j])), np.asscalar(norm(values[j+1])), nsample+2)
            newvalues.extend(norm.inverse(interp))
        # Create LineCollection and update with values
        newvalues  = np.array(newvalues)
        points     = np.array([newx, newy]).T.reshape(-1, 1, 2) # -1 means value is inferred from size of array, neat!
        segments   = np.concatenate([points[:-1], points[1:]], axis=1)
        collection = mcollections.LineCollection(segments, cmap=cmap, norm=norm, linestyles='-')
        collection.set_array(newvalues)
        collection.update({k:v for k,v in kwargs.items() if k not in ['color']})
        line = self.add_collection(collection) # FIXME: for some reason using 'line' as the mappable results in colorbar with color *cutoffs* at values, instead of centered levels at values
        line = colortools.mappable(cmap, values, norm=norm) # use hacky mchackerson instead
        return line,

    @cycle_features
    def plot(self, *args, cmap=None, values=None, **kwargs):
        """
        Expand functionality of plot to also make LineCollection lines, i.e. lines
        whose colors change as a function of some key/indicator.
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

    # Apply some simple enhancements
    # These just fix the white edges between patch objects, and optionally
    # enable using same color for data in 'edge bins' as 'out of bounds' data
    @cmap_features
    def contour(self, *args, **kwargs):
        return super().contour(*args, **kwargs)

    @cmap_features
    def contourf(self, *args, **kwargs):
        return super().contourf(*args, **kwargs)

    @cmap_features
    def pcolorpoly(self, *args, **kwargs):
        return super().pcolor(*args, **kwargs)

    @cmap_features
    def pcolormesh(self, *args, **kwargs):
        return super().pcolor(*args, **kwargs)

    # Matrix visualizations
    @cmap_features
    def matshow(self, *args, **kwargs):
        return super().matshow(*args, **kwargs)

    @cmap_features
    def imshow(self, *args, **kwargs):
        return super().imshow(*args, **kwargs)

    @cmap_features
    def spy(self, *args, **kwargs):
        return super().spy(*args, **kwargs)

    @cmap_features
    def hist2d(self, *args, **kwargs): # expand to allow 2D kernel
        return super().hist2d(*args, **kwargs)

    # Stuff that needs to be worked on
    @cycle_features
    def polar(self, *args, **kwargs): # dunno
        return super().polar(*args, **kwargs)

    @cycle_features
    def bar(self, *args, **kwargs): # have bar plot cycle through default color cycle
        return super().bar(*args, **kwargs)

    @cycle_features
    def barh(self, *args, **kwargs):
        return super().barh(*args, **kwargs)

    @cycle_features
    def hist(self, *args, **kwargs): # expand so we can apply transparency to histograms, so they can be overlaid
        return super().hist(*args, **kwargs)

    @cycle_features
    def boxplot(self, *args, **kwargs): # box whisker plot
        return super().boxplot(*args, **kwargs)

    @cycle_features
    def violinplot(self, *args, **kwargs): # ditto with box plot function
        return super().violinplot(*args, **kwargs)

    def errorbar(self, *args, **kwargs):
        return super().errorbar(*args, **kwargs)

    # Redundant stuff that want to cancel
    message1 = 'Redundant function has been disabled.'
    @property
    def plot_date(self): # use xdates=True or ydates=True
        raise NotImplementedError(self.message1)
    @property
    def semilogx(self): # use xscale='log' instead, this is dumb!
        raise NotImplementedError(self.message1)
    @property
    def semilogy(self):
        raise NotImplementedError(self.message1)
    @property
    def loglog(self):
        raise NotImplementedError(self.message1)

    # Weird stuff that probably will not wrap
    message2 = 'Unsupported plotting function.'
    @property
    def pie(self): # dunno
        raise NotImplementedError(self.message2)
    @property
    def table(self): # dude... fugly
        raise NotImplementedError(self.message2)
    @property
    def hexbin(self): # expand to allow making grouped boxes
        raise NotImplementedError(self.message2)
    @property
    def eventplot(self):
        raise NotImplementedError(self.message2)
    @property
    def triplot(self):
        raise NotImplementedError(self.message2)
    @property
    def tricontour(self):
        raise NotImplementedError(self.message2)
    @property
    def tricontourf(self):
        raise NotImplementedError(self.message2)
    @property
    def tripcolor(self):
        raise NotImplementedError(self.message2)

    # Disable spectral and triangular features
    # See: https://stackoverflow.com/a/23126260/4970632
    # Also see: https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/axes/_axes.py
    # for all Axes methods ordered logically in class declaration.
    @property
    def xcorr(self):
        raise NotImplementedError(self.message2)
    @property
    def acorr(self):
        raise NotImplementedError(self.message2)
    @property
    def psd(self):
        raise NotImplementedError(self.message2)
    @property
    def csd(self):
        raise NotImplementedError(self.message2)
    @property
    def magnitude_spectrum(self):
        raise NotImplementedError(self.message2)
    @property
    def angle_spectrum(self):
        raise NotImplementedError(self.message2)
    @property
    def phase_spectrum(self):
        raise NotImplementedError(self.message2)
    @property
    def cohere(self):
        raise NotImplementedError(self.message2)
    @property
    def specgram(self):
        raise NotImplementedError(self.message2)

@docstring_fix
class XYAxes(Axes):
    """
    Subclass for ordinary Cartesian-grid axes.
    """
    # Initialize
    name = 'xy'
    def __init__(self, *args, **kwargs):
        """
        Create simple x by y subplot.
        """
        super().__init__(*args, **kwargs)

    # Cool overrides
    def format(self,
        xgrid=None,      ygrid=None,      # gridline toggle
        xdates=False,    ydates=False,    # whether to format axis labels as long datetime strings; the formatter should be a date %-style string
        xspineloc=None,  yspineloc=None,  # deals with spine options
        xtickminor=True, ytickminor=True, # minor ticks on/off
        xgridminor=None, ygridminor=None, # minor grids on/off (if ticks off, grid will always be off)
        xtickloc=None,   ytickloc=None,   # which spines to draw ticks on
        xtickdir=None,   ytickdir=None,   # which direction ('in', 'our', or 'inout')
        xticklabeldir=None, yticklabeldir=None, # which direction to draw labels
        xtickrange=None,    ytickrange=None,    # limit regions where we assign ticklabels to major-ticks
        xreverse=False, yreverse=False, # special properties
        xlabel=None,    ylabel=None,    # axis labels
        xlim=None,      ylim=None,
        xscale=None,    yscale=None,    xscale_kwargs={},    yscale_kwargs={},
        xlocator=None,   xminorlocator=None, ylocator=None, yminorlocator=None, # locators, or derivatives that are passed to locators
        xformatter=None, yformatter=None,
        **kwargs): # formatter
        """
        Format the x/y labels, tick locators, tick formatters, and more.
        Needs more documentation.
        """
        # Pass stuff to parent formatter, e.g. title and abc labeling
        super().format(**kwargs)

        # Set axis scaling and limits
        if xscale is not None:
            if hasattr(xscale,'name'):
                xscale = xscale.name
            self.set_xscale(xscale, **xscale_kwargs)
        if yscale is not None:
            if hasattr(yscale,'name'):
                yscale = yscale.name
            self.set_yscale(yscale, **yscale_kwargs)
        if xlim is None:
            xlim = self.get_xlim()
        if ylim is None:
            ylim = self.get_ylim()
        if xreverse:
            xlim = xlim[::-1]
        if yreverse:
            ylim = ylim[::-1]
        self.set_xlim(xlim)
        self.set_ylim(ylim)

        # Control axis ticks and labels and stuff
        for xy, axis, label, tickloc, spineloc, gridminor, tickminor, tickminorlocator, \
                grid, ticklocator, tickformatter, tickrange, tickdir, ticklabeldir in \
            zip('xy', (self.xaxis, self.yaxis), (xlabel, ylabel), \
                (xtickloc,ytickloc), (xspineloc, yspineloc), # other stuff
                (xgridminor, ygridminor), (xtickminor, ytickminor), (xminorlocator, yminorlocator), # minor ticks
                (xgrid, ygrid), (xlocator, ylocator), (xformatter, yformatter), # major ticks
                (xtickrange, ytickrange), # range in which we label major ticks
                (xtickdir, ytickdir), (xticklabeldir, yticklabeldir)): # tick direction
            # Axis spine visibility and location
            sides = ('bottom','top') if xy=='x' else ('left','right')
            for spine, side in zip((self.spines[s] for s in sides), sides):
                # Line properties
                spine.update(rc('spine'))
                spineloc = getattr(self, f'{xy}spine_override', spineloc) # optionally override; necessary for twinx/twiny situation
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

            # Set the major and minor locators and formatters
            # Also automatically detect whether axis is a 'time axis' (i.e.
            # whether user has plotted something with x/y as datetime/date/np.datetime64
            # objects, and matplotlib automatically set the unit converter)
            time = isinstance(axis.converter, mdates.DateConverter)
            axis.set_major_locator(Locator(ticklocator, time=time))
            axis.set_major_formatter(Formatter(tickformatter, tickrange=tickrange, time=time))
            if not tickminor and tickminorlocator is None:
                axis.set_minor_locator(Locator('null'))
            else:
                axis.set_minor_locator(Locator(tickminorlocator, minor=True, time=time))
            axis.set_minor_formatter(mticker.NullFormatter())
            # FIXME: Is this necessary, or shouldn't rcParams have taken care of it?
            for t in axis.get_ticklabels():
                t.update(rc('ticklabels'))

            # Tick properties
            # * Weird issue seems to cause set_tick_params to reset/forget that the grid
            #   is turned on if you access tick.gridOn directly, instead of passing through tick_params.
            #   Since gridOn is undocumented feature, don't use it. So calling _format_axes() a second time will remove the lines
            # * Can specify whether the left/right/bottom/top spines get ticks; sides will be 
            #   group of left/right or top/bottom
            # * Includes option to draw spines but not draw ticks on that spine, e.g.
            #   on the left/right edges
            ticklocs = sides if tickloc=='both' else () if tickloc=='neither' else tickloc
            if ticklocs is None: # only turn off ticks if the spine is invisible; don't explicitly turn on
                ticks_sides = {side: False for side in sides if not self.spines[side].get_visible()}
            else:
                ticks_sides = {side: self.spines[side].get_visible() and side in ticklocs for side in sides}
            ticks_major, ticks_minor = rc('tick'), rc('tickminor')
            if tickdir is not None:
                ticks_major.update({'direction':tickdir})
                ticks_minor.update({'direction':tickdir})
            if tickdir=='in':
                ticks_major.update({'pad':1}) # should be much closer
                ticks_minor.update({'pad':1})
            if ticklabeldir=='in': # put tick labels inside the plot; sometimes might actually want this
                pad = rc('globals', 'tickpad') + rc('globals', 'small') + rc('globals', 'ticklen')
                ticks_major.update({'pad':-pad})
                ticks_minor.update({'pad':-pad})
            axis.set_tick_params(which='major', **ticks_sides, **ticks_major)
            axis.set_tick_params(which='minor', **ticks_sides, **ticks_minor) # have length

            # Gridline activation and setting (necessary because rcParams has no 'minorgrid'
            # property, must be set in rcExtras)
            # NOTE: Inexplicably, for a twinx axis, could only get the minor gridlines
            # to disappear if we changed the 'visible' property on each one.
            # Fugly way, which lets us not explicitly turn on/off gridlines if user did
            # not specify (by updating objects directly) for tick in axis.majorTicks
            for tick in axis.get_major_ticks():
                if grid is not None:
                    tick.gridline.set_visible(grid)
                tick.gridline.update(rc('grid'))
            # for tick in axis.minorTicks:
            for tick in axis.get_minor_ticks():
                if gridminor is not None:
                    tick.gridline.set_visible(gridminor)
                tick.gridline.update(rc('gridminor'))
            # For some insane reasion, these are ***both*** needed
            # Without this below stuff, e.g. gridminor=True doesn't draw gridlines
            if grid is not None: # grid changes must be after tick
                axis.grid(grid, which='major', **rc('grid'))
            if gridminor is not None:
                axis.grid(gridminor, which='minor', **rc('gridminor')) # ignore if no minor ticks

            # Label properties
            axis.label.update(rc('label'))
            if label is not None:
                axis.label.set_text(label)

        # Finished
        return self

    def twiny(self, **kwargs):
        """
        Create second x-axis extending from shared ("twin") y-axis
        Note: Cannot wrap twiny() because then the axes created will be
        instantiated from the parent class, which doesn't have format() method.
        Instead, use hidden method _make_twin_axes.
        """
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
        self.xspine_override = 'bottom' # original axis ticks on bottom
        ax.xspine_override   = 'top' # new axis ticks on top
        ax.yspine_override   = 'neither'
        return ax

    def twinx(self, **kwargs):
        """
        Create second y-axis extending from shared ("twin") x-axis
        Note: Cannot wrap twinx() because then the axes created will be
        instantiated from the parent class, which doesn't have format() method.
        Instead, use hidden method _make_twin_axes.
        """
        # As above
        print(self.name)
        ax = self._make_twin_axes(sharex=self, projection=self.name)
        self.yaxis.tick_left()
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')
        ax.yaxis.set_offset_position('right')
        ax.set_autoscalex_on(self.get_autoscalex_on())
        ax.xaxis.set_visible(False)
        ax.patch.set_visible(False)
        ax.grid(False)
        # Special settings, force spine locations when format() called
        self.yspine_override = 'left' # original axis ticks on left
        ax.yspine_override   = 'right' # new axis ticks on right
        ax.xspine_override   = 'neither'
        return ax

    def make_inset_locator(self, bounds, trans):
        """
        Helper function, had to be copied from private matplotlib version.
        """
        def inset_locator(ax, renderer):
            bbox = mtransforms.Bbox.from_bounds(*bounds)
            bb = mtransforms.TransformedBbox(bbox, trans)
            tr = self.figure.transFigure.inverted()
            bb = mtransforms.TransformedBbox(bb, tr)
        return bb

    def inset_axes(self, bounds, *, transform=None, zorder=5,
            **kwargs):
        """
        Create inset of same type.
        """
        # Defaults
        if transform is None:
            transform = self.transAxes
        label = kwargs.pop('label', 'inset_axes')
        # This puts the rectangle into figure-relative coordinates.
        locator = self.make_inset_locator(bounds, transform)
        bb = locator(None, None)
        ax = AxesXY(self.figure, bb.bounds, zorder=zorder, label=label, **kwargs)
        # The following locator lets the axes move if in data coordinates, gets called in ax.apply_aspect()
        ax.set_axes_locator(locator)
        self.add_child_axes(ax)
        return ax


    # Wrappers that check dimensionality, give better error messages, and change
    # the convention so dim1 is always x-axis, dim2 is always y-axis
    @check_centers
    def contour(self, *args, **kwargs):
        return super().contour(*args, **kwargs)

    @check_centers
    def contourf(self, *args, **kwargs):
        return super().contourf(*args, **kwargs)

    @check_centers
    def quiver(self, *args, **kwargs):
        return super().quiver(*args, **kwargs)

    @check_centers
    def streamplot(self, *args, **kwargs):
        return super().streamplot(*args, **kwargs)

    @check_centers
    def barbs(self, *args, **kwargs):
        return super().barbs(*args, **kwargs)

    @check_edges
    def pcolorpoly(self, *args, **kwargs):
        return super().pcolor(*args, **kwargs)

    @check_edges
    def pcolormesh(self, *args, **kwargs):
        return super().pcolormesh(*args, **kwargs)

class MapAxes(Axes):
    """
    Dummy intermediate class that just disables a bunch of methods
    inappropriate for map projections.
    """
    # This one we want to work, but it doesn't
    # Property decorators mean error will be raised right after method invoked,
    # will not raise error on account of args/kwargs, see: https://stackoverflow.com/a/23126260/4970632
    @property
    def pcolormesh(self):
        raise NotImplementedError('Mesh version of pcolor fails for map projections. Use pcolorpoly instead.')
    # Disable some methods to prevent weird shit from happening
    message = 'Invalid plotting function for map projection axes.'
    @property
    def matshow(self):
        raise NotImplementedError(self.message)
    @property
    def imshow(self):
        raise NotImplementedError(self.message)
    @property
    def spy(self):
        raise NotImplementedError(self.message)
    @property
    def polar(self):
        raise NotImplementedError(self.message)
    @property
    def bar(self):
        raise NotImplementedError(self.message)
    @property
    def barh(self):
        raise NotImplementedError(self.message)
    @property
    def hist(self):
        raise NotImplementedError(self.message)
    @property
    def hist2d(self):
        raise NotImplementedError(self.message)
    @property
    def errorbar(self):
        raise NotImplementedError(self.message)
    @property
    def boxplot(self):
        raise NotImplementedError(self.message)
    @property
    def violinplot(self):
        raise NotImplementedError(self.message)
    @property
    def step(self):
        raise NotImplementedError(self.message)
    @property
    def stem(self):
        raise NotImplementedError(self.message)
    @property
    def hlines(self):
        raise NotImplementedError(self.message)
    @property
    def vlines(self):
        raise NotImplementedError(self.message)
    @property
    def axhline(self):
        raise NotImplementedError(self.message)
    @property
    def axvline(self):
        raise NotImplementedError(self.message)
    @property
    def axhspan(self):
        raise NotImplementedError(self.message)
    @property
    def axvspan(self):
        raise NotImplementedError(self.message)
    @property
    def fill_between(self):
        raise NotImplementedError(self.message)
    @property
    def fill_betweenx(self):
        raise NotImplementedError(self.message)
    @property
    def fill(self):
        raise NotImplementedError(self.message)
    @property
    def stackplot(self):
        raise NotImplementedError(self.message)

@docstring_fix
class BasemapAxes(MapAxes):
    """
    Axes subclass for basemap plotting.
    """
    name = 'basemap'
    pseudocyl = ['moll','robin','eck4','kav7','sinu','mbtfpq','vandg','hammer']
    def __init__(self, *args, map_projection=None, **kwargs):
        """
        Declare basemap projection instance, add it as the 'm' attribute.
        The 'map_projection' argument sets projection, because this axes itself
        is called from add_subplot using projection='basemap'.
        """
        # * Must set boundary before-hand, otherwise the set_axes_limits method called
        #   by mcontourf/mpcolormesh/etc draws two mapboundary Patch objects called "limb1" and
        #   "limb2" automatically: one for fill and the other for the edges
        # * Then, since the patch object in _mapboundarydrawn is only the fill-version, calling
        #   drawmapboundary() again will replace only *that one*, but the original visible edges
        #   are still drawn -- so e.g. you can't change the color
        # * If you instead call drawmapboundary right away, _mapboundarydrawn will contain
        #   both the edges and the fill; so calling it again will replace *both*
        if not isinstance(map_projection, mbasemap.Basemap):
            raise ValueError('You must initialize BasemapAxes with map_projection=(basemap.Basemap instance).')
        super().__init__(*args, **kwargs)
        self.m = map_projection
        self.boundary = None
        self.recurred = False # use this so we can override plotting methods
        if map_projection.projection in self.pseudocyl: # otherwise the spines are map boundary
            self.boundary = self.m.drawmapboundary(ax=self)

    def format(self, color=None,
        oceans=False, coastlines=False, continents=False, # coastlines and continents
        latlabels=[0,0,0,0], lonlabels=[0,0,0,0], # sides for labels [left, right, bottom, top]
        latlocator=None, latminorlocator=None,
        lonlocator=None, lonminorlocator=None,
        **kwargs):
        """
        Format basemap axes.
        Add documentation here.
        """
        # Pass stuff to parent formatter, e.g. title and abc labeling
        super().format(**kwargs)
        # Basemap axes setup
        # Coastlines, parallels, meridians
        if coastlines:
            p = self.m.drawcoastlines(**rc('coastlines'), ax=self)
        if continents:
            p = self.m.fillcontinents(**rc('continents'), ax=self)
        # Longitude/latitude lines
        # Make sure to turn off clipping by invisible axes boundary; otherwise
        # get these weird flat edges where map boundaries, parallel/meridian markers come up to the axes bbox
        if latlocator is not None:
            if not utils.isiterable(latlocator):
                latlocator = AutoLocate(self.m.latmin+latlocator, self.m.latmax-latlocator, latlocator)
            p = self.m.drawparallels(latlocator, labels=latlabels, ax=self)
            for pi in p.values(): # returns dict, where each one is tuple
                for _ in [i for j in pi for i in j]: # magic
                    if isinstance(_, mtext.Text):
                        _.update(rc('ticklabels'))
                    else:
                        _.set_clip_on(True) # no gridlines past boundary
                        _.update(rc('lonlatlines'))
                        # _.set_linestyle(linestyle)
                # tried passing clip_on to the above, but it does nothing; must set
                # for lines created after the fact
        if lonlocator is not None:
            latlabels[2:] = latlabels[2:][::-1] # default is left/right/top/bottom which is dumb
            lonlabels[2:] = lonlabels[2:][::-1] # change to left/right/bottom/top
            if not utils.isiterable(lonlocator):
                lonlocator = AutoLocate(self.m.lonmin+lonlocator, self.m.lonmax-lonlocator, lonlocator)
            p = self.m.drawmeridians(lonlocator, labels=lonlabels, ax=self)
            for pi in p.values():
                for _ in [i for j in pi for i in j]: # magic
                    if isinstance(_, mtext.Text):
                        _.update(rc('ticklabels'))
                    else:
                        _.set_clip_on(True) # no gridlines past boundary
                        _.update(rc('lonlatlines'))
                        # _.set_linestyle(linestyle)
        # Map boundary
        # * First have to *manually replace* the old boundary by just deleting
        #   the original one; note this requires drawmapboundary() was called
        #   when the basemap was first instantiated; see notes in subplots() command.
        # * If boundary is drawn successfully should be able to call
        #   self.m._mapboundarydrawn.set_visible(False) and edges/fill color disappear
        # * For now will enforce that map plots *always* have background whereas
        #   axes plots can have transparent background
        if self.m._mapboundarydrawn:
            self.m._mapboundarydrawn.remove()
        if self.m.projection in self.pseudocyl:
            self.patch.set_alpha(0) # make patch invisible
            p = self.m.drawmapboundary(fill_color=color, ax=self, **rc('spine')) # set fill_color to 'none' to make transparent
            p.set_rasterized(False) # not sure about this; might be rasterized
            p.set_clip_on(False) # so edges of LINE denoting boundary aren't cut off
        else: # use the settings to apply to Axes patch; the basemap API fails here
            self.patch.set_facecolor(color)
            self.patch.set_edgecolor('none')
            for spine in self.spines.values():
                spine.update(rc('spine'))
        return self

    # Basemap overrides
    # The decorators assure that, when method is called the second time by
    # the Basemap copy, we will then call the superclass Axes version instead.
    # For plot/scatter we can just use latlon=True to naively convert to map
    # coordinates, but want to be safe for contour/pcolor methods.
    @norecurse
    def plot(self, *args, **kwargs):
        return self.m.plot(*args, ax=self, latlon=True, **kwargs) # this will call *itself*

    @norecurse
    def scatter(self, *args, **kwargs):
        return self.m.scatter(*args, ax=self, latlon=True, **kwargs)

    @norecurse
    @check_centers
    @gridfix_basemap
    def contour(self, *args, **kwargs):
        return self.m.contour(*args, ax=self, **kwargs)

    @norecurse
    @check_centers
    @gridfix_basemap
    def contourf(self, *args, **kwargs):
        return self.m.contourf(*args, ax=self, **kwargs)

    @norecurse
    @check_edges
    @gridfix_basemap
    def pcolorpoly(self, *args, **kwargs):
        return self.m.pcolor(*args, ax=self, **kwargs)

    @norecurse
    @check_centers
    @gridfix_basemap
    def barbs(self, *args, **kwargs):
        return super().barbs(*args, **kwargs)

    @norecurse
    @check_centers
    @gridfix_basemap
    def quiver(self, *args, **kwargs):
        return super().quiver(*args, **kwargs)

    @norecurse
    @check_centers
    @gridfix_basemap
    def streamplot(self, *args, **kwargs):
        return super().streamplot(*args, **kwargs)

@docstring_fix
# class CartopyAxes(GeoAxes, MapAxes):
class CartopyAxes(MapAxes, GeoAxes): # custom one has to be higher priority, so the methods can overwrite stuff
    """
    Cartopy subclass.
    Some notes:
    * Cartopy takes advantage of documented feature where any class with method
      named _as_mpl_axes can be passed as 'projection' object.
      Feature documented here: https://matplotlib.org/devel/add_new_projection.html
      Used in Projection parent class here: https://scitools.org.uk/cartopy/docs/v0.13/_modules/cartopy/crs
    """
    name = 'cartopy'
    def __init__(self, *args, map_projection=None, circle_center=90, circle_edge=0, **kwargs):
        """
        Initialize cartopy projection, and allow for *partial* (i.e. not global) coverage for
        azimuthal projections by zooming into the full projection, then drawing a circle boundary
        around some latitude away from the center.
        * The 'map_projection' argument sets projection, because this axes itself
          is called from add_subplot using projection='basemap'.
        * Number 'ncircle' controls number of points for drawing circular projection boundary.
          For more info see: https://scitools.org.uk/cartopy/docs/v0.15/examples/always_circular_stereo.html
        """
        # Initialize
        # Do the GeoAxes initialization steps manually, because there are very few, and since GeoAxes and Axes subclass
        # same base classes, might get some repitition or weird errors if try to do both.
        if not isinstance(map_projection, ccrs.Projection):
            raise ValueError('You must initialize CartopyAxes with map_projection=(cartopy.crs.Projection instance).')
        # super().__init__(*args, map_projection=map_projection, **kwargs)
        self.projection = map_projection # attribute used extensively by GeoAxes methods, and by builtin one
        super(GeoAxes, self).__init__(*args, **kwargs)
        self.outline_patch = None    # patch that provides the line bordering the projection.
        self.background_patch = None # patch that provides the filled background of the projection
        self.img_factories = []
        self._gridliners = []
        self._done_img_factory = False
        self._boundary()
        # Apply circle boundary
        crs_circles = (ccrs.LambertAzimuthalEqualArea, ccrs.AzimuthalEquidistant)
        if any(isinstance(map_projection, cp) for cp in crs_circles):
            # self.projection.threshold = kwargs.pop('threshold', self.projection.threshold) # optionally modify threshold
            self.set_extent([-180, 180, circle_edge, circle_center], PlateCarree) # use platecarree transform
            self.set_boundary(circle(100), transform=self.transAxes)

    def _as_mpl_axes(self):
        """
        See above for details. Since this class is *itself* a projection, we cannot
        instantiate it with another projection, will get weird error. Instead manually
        declare this method just like in ccrs.
        """
        return geoaxes.GeoAxes, {'map_projection': self}

    def format(self, color=None,
        oceans=False, coastlines=False, continents=False, # coastlines and continents
        latlabels=[0,0,0,0], lonlabels=[0,0,0,0], # sides for labels [left, right, bottom, top]
        latlocator=None, latminorlocator=None, lonlocator=None, lonminorlocator=None,
        **kwargs):
        """
        Format cartopy GeoAxes.
        Add documentation here.
        """
        # Pass stuff to parent formatter, e.g. title and abc labeling
        super().format(**kwargs)
        # Boundary for projection regoin
        self.outline_patch.update(rc('outline'))
        # Make background patch invisible, so can color axes patch instead
        # See: https://stackoverflow.com/a/32208700/4970632
        # self.background_patch.set_fill(False)
        # This doesn't seem to work because axes patch is always invisible and has
        # been overridden by background patch, so will just show underlying figure color
        if color is not None:
            self.background_patch.set_facecolor(color)
        else:
            self.background_patch.set_facecolor('w')
        # Draw gridlines
        # WARNING: For some reason very weird side effects happen if you try
        # to call gridlines() twice on same axes. Can't do it. Which is why
        # we do this nonsense with the formatter below, instead of drawing 'major'
        # grid lines and 'minor' grid lines.
        lonvec = lambda v: [] if v is None else [*v] if utils.isiterable(v) else [*AutoLocate(-180,180,v)]
        latvec = lambda v: [] if v is None else [*v] if utils.isiterable(v) else [*AutoLocate(-90,90,v)]
        lonminorlocator, latminorlocator = lonvec(lonminorlocator), latvec(latminorlocator)
        lonlocator, latlocator = lonvec(lonlocator), latvec(latlocator)
        lonlines = lonminorlocator or lonlocator # where we draw gridlines
        latlines = latminorlocator or latlocator
        # First take care of gridlines
        draw_labels = (isinstance(self.projection,ccrs.Mercator) or isinstance(self.projection,ccrs.PlateCarree))
        if latlines and latlines[0]==-90:
            latlines[0] += 0.001
        if lonlines and lonlines[0]==-90:
            lonlines[0] -= 0.001
        gl = self.gridlines(**rc('lonlatlines'), draw_labels=draw_labels)
        gl.xlocator = mticker.FixedLocator(lonlines)
        gl.ylocator = mticker.FixedLocator(latlines)
        # Now take care of labels
        if draw_labels:
            lonfunc = lambda x,y: LONGITUDE_FORMATTER(x) if x in lonlocator else ''
            latfunc = lambda x,y: LATITUDE_FORMATTER(x) if x in latlocator else ''
            gl.xformatter = mticker.FuncFormatter(lonfunc)
            gl.yformatter = mticker.FuncFormatter(latfunc)
            gl.xlabels_bottom, gl.xlabels_top = latlabels[2:]
            gl.ylabels_left, gl.ylabels_right = lonlabels[:2]
        # Add geographic features
        # Use the NaturalEarthFeature to get more configurable resolution; can choose
        # between 10m, 50m, and 110m (scales 1:10mil, 1:50mil, and 1:110mil)
        if coastlines:
            # self.add_feature(cfeature.COASTLINE, **rc('coastlines'))
            print('Add coastlines.')
            self.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m'), **rc('coastlines'))
        if continents:
            print('Add continents.')
            # self.add_feature(cfeature.LAND, **rc('continents'))
            self.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m'), **rc('continents'))
        if oceans:
            print('Add oceans.')
            # self.add_feature(cfeature.OCEAN, **rc('oceans'))
            self.add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', '50m'), **rc('oceans'))
        return self

    # Cartopy overrides
    # Does some simple manipulation and adds the default transform PlateCarree if not declared.
    def plot(self, *args, transform=PlateCarree, **kwargs):
        return super().plot(*args, transform=transform, **kwargs) # this will call *itself*

    def scatter(self, *args, transform=PlateCarree, **kwargs):
        return super().scatter(*args, transform=transform, **kwargs)

    @check_centers
    @gridfix_cartopy
    def contour(self, *args, **kwargs):
        return super().contour(*args, **kwargs)

    @check_centers
    @gridfix_cartopy
    def contourf(self, *args, **kwargs):
        return super().contourf(*args, **kwargs)

    @check_edges
    @gridfix_cartopy
    def pcolorpoly(self, *args, **kwargs):
        return super().pcolorpoly(*args, **kwargs)

    @check_centers
    @gridfix_cartopy
    def barbs(self, *args, **kwargs):
        return super().barbs(*args, **kwargs)

    @check_centers
    @gridfix_cartopy
    def quiver(self, *args, **kwargs):
        return super().quiver(*args, **kwargs)

    @check_centers
    @gridfix_cartopy
    def streamplot(self, *args, **kwargs):
        return super().streamplot(*args, **kwargs)

class Panel(object):
    """
    Helper class for 'dummy axes' that defers drawing until user attempts
    to access its attributes, e.g. in calling a method.
    For inheriting from instantiated object see: https://stackoverflow.com/a/33468799/4970632
    but this only works for copying attributes, not methods.
    Goal of this is to have generalized panel class suitable for both 'inner'
    and 'outer' panels.
    """
    def __init__(self, figure, subspec, side, n=1, length=1):
    # def __init__(self, figure, subspec, n=1, length=1, side='bottom', projection=None):
        """
        Initialize abstract panel object.
        Behavior depends on which attribute you access first:
            * Try to access colorbar or legend, and this is a default
              axes.Axes instance.
            * Try to access anything else, and this is 
        """
        # Assignments
        npanel = n
        self.npanel  = npanel # number of panels stacked together
        self.drawn   = False # defer drawing these panels
        self.figure  = figure
        self.side    = side
        self.index   = 0 # axes index that user has queried
        self.length  = length # length relative to full length
        self.axs     = [None for i in range(npanel)]
        # Generate list of subspecs
        # Will optionally shrink in the lengthwise direction and stack
        # multiple panels in the crosswise dimension
        if side in ('bottom','top'):
            gridspec = mgridspec.GridSpecFromSubplotSpec(
                    nrows        = npanel,
                    ncols        = 3,
                    wspace       = 0,
                    hspace       = 0,
                    subplot_spec = subspec,
                    width_ratios = ((1-length)/2, length, (1-length)/2)
                    )
            self.subspecs = [gridspec[i,1] for i in range(npanel)]
        elif side in ('left','right'):
            gridspec = mgridspec.GridSpecFromSubplotSpec(
                    nrows         = 3,
                    ncols         = npanel,
                    hspace        = 0,
                    wspace        = 0,
                    subplot_spec  = subspec,
                    height_ratios = ((1-length)/2, length, (1-length)/2)
                    )
            self.subspecs = [gridspec[1,i] for i in range(npanel)]
        else:
            raise ValueError(f'Invalid panel side "{side}".')

    def __getitem__(self, key):
        """
        Set the axes index user is querying. Also don't implement __delitem__
        and __setitem__ because we'd want to forbid that anyway!
        This is done so user can, before doing anything else, call (e.g.)
        ax.bottompanel[i].colorbar, and under-the-hood, axes will be drawn
        then the i'th axes accessed.
        """
        if key<0:
            key = self.npanel + key # e.g. -1
        if key<0 or key>self.npanel-1:
            raise ValueError(f'Key {key} invalid, only {self.npanel} panels present.')
        self.index = key

    def __getattribute__(self, attribute, *args):
        """
        When method is first invoked, ensure axes are drawn.
        Depending on which property user accesses first, axes can be
        allocated as:
          * A legend box (just dead space for storing a legend)
          * A colorbar (optionally shrunken relative to subplotspec length)
          * A normal axes (for plotting and stuff)
        """
        # First, if trying to access attribute (not method) of the panel,
        # do not do anything
        getter = super().__getattribute__
        if attribute in getter('__dict__'):
            return getter(attribute)
        # Next, instantiate axes
        if not getter('axs')[getter('index')]: # use super-class method
            getter('instantiate')(attribute)
        # Now that axes are initiated, optionally return a couple
        # overridden functions and methods here here
        if attribute in ('colorbar','legend'):
            return getter(attribute, *args)
        else:
            return getter('axs')[getter('index')].__getattribute__(attribute, *args)

    # def instantiate(self):
    def instantiate(self, attribute):
        """
        Function for instantiating axes belonging to this panel.
        Will read from the 'n' attribute to figure out how many
        panels to draw.
        """
        # This function is invoked from within __getattribute__
        # so need to make sure we don't trigger infinite loops
        getter = super().__getattribute__
        axs        = getter('axs')
        index      = getter('index')
        figure     = getter('figure')
        subspecs   = getter('subspecs')
        # projection = getter('projection')
        projection = None if attribute in ('colorbar','legend') else 'xy'
        axs[index] = figure.add_subplot(subspecs[index], projection=projection)

    def legend(self, handles, **kwargs):
        """
        Allocate invisible axes for drawing legend.
        Returns the axes and the output of legend_factory().
        """
        # Simply call legend on axes, but also work some magic
        # to make axes invisible
        ax = self.axs[self.index] # get axes
        for s in ax.spines.values():
            s.set_visible(False)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.patch.set_alpha(0)
        kwnew = {'frameon':False, 'loc':'center', 'borderaxespad':0,
                 'bbox_transform':ax.transAxes}
        kwnew.update(kwargs)
        return ax, legend_factory(ax, handles, **kwnew)

    def colorbar(self, *args, inner=False, outer=False, length=1, **kwargs):
        """
        Allocate axes for drawing colorbar.
        Returns the axes and the output of colorbar_factory().
        """
        # Simply call the 
        ax = self.axs[self.index] # get axes
        side = self.side
        index = self.index
        npanel = self.npanel
        if side in ('bottom','top'):
            ticklocation = 'bottom' if index==npanel-1 else 'top'
            orientation  = 'horizontal'
        elif side in ('left','right'):
            ticklocation = 'right' if index==npanel-1 else 'right'
            orientation = 'vertical'
        kwargs.update({'orientation':orientation, 'ticklocation':ticklocation})
        return ax, colorbar_factory(ax, *args, **kwargs)

# Register the projection
mprojections.register_projection(XYAxes)
mprojections.register_projection(BasemapAxes)
mprojections.register_projection(CartopyAxes)

#------------------------------------------------------------------------------#
# Custom legend and colorbar factories
#------------------------------------------------------------------------------#
def projection_factory(package, projection, **kwargs):
    """
    Returns Basemap object or cartopy ccrs instance.
    """
    # Initial stuff
    silent = kwargs.pop('silent', False)
    crs_translate = { # less verbose keywords, actually match proj4 keywords and are similar to basemap
        **{k:'central_latitude'  for k in ['lat0','lat_0']},
        **{k:'central_longitude' for k in ['lon0', 'lon_0']},
        }
    crs_dict = { # interpret string, create cartopy projection
        **{key: ccrs.PlateCarree   for key in ['cyl','rectilinear','pcarree','platecarree']},
        **{key: ccrs.Mollweide     for key in ['moll','mollweide']},
        **{key: ccrs.Stereographic for key in ['stereo','stereographic']},
        'aeqd': ccrs.AzimuthalEquidistant, 'aeqa': ccrs.LambertAzimuthalEqualArea,
        'mercator': ccrs.Mercator, 'robinson': ccrs.Robinson, 'ortho': ccrs.Orthographic,
        'hammer': Hammer,       'aitoff': Aitoff,
        'wintri': WinkelTripel, 'kav7':   KavrayskiyVII,
        }
    # Create projection and determine required aspect ratio
    if package=='basemap':
        projection = mbasemap.Basemap(projection=(projection or 'cyl'), **{**kwargs, 'fix_aspect':True}) # cylindrical by default
        aspect = (projection.urcrnrx - projection.llcrnrx)/(projection.urcrnry - projection.llcrnry)
    # Get the projection instance from a string and determine required aspect ratio
    elif package=='cartopy':
        projection = projection or 'cyl'
        if projection not in crs_dict:
            raise ValueError(f'For cartopy, projection must be one of the following: {", ".join(crs_dict.keys())}.')
        projection = crs_dict[projection](**{crs_translate.get(key,key):value for key,value in kwargs.items()})
        aspect = (np.diff(projection.x_limits)/np.diff(projection.y_limits))[0]
    if not silent:
        print(f'Forcing aspect ratio: {aspect:.3g}')
    return projection, aspect

def legend_factory(ax, handles=None, align=None, handlefix=False, **kwargs): #, settings=None): # can be updated
    """
    Function for formatting legend-axes (invisible axes with centered legends on them).
    Should update my legend function to CLIP the legend box when it goes outside axes area, so
    the legend-width and bottom/right widths can be chosen propertly/separately.
    """
    # First get legend settings (usually just one per plot so don't need to declare
    # this dynamically/globally), and interpret kwargs
    lsettings = rc('legend')
    if 'ncols' in kwargs:
        kwargs['ncol'] = kwargs.pop('ncols') # pyplot subplot uses 'ncols', but legend uses 'ncol'... annoying!
    if 'frame' in kwargs: # again, confusing choice
        kwargs['frameon'] = kwargs.pop('frame')
    lsettings.update(**kwargs)
    # Setup legend text and handle properties
    tsettings = rc('ticklabels')
    hsettings = {}
    for candidate in ['linewidth', 'color']: # candidates for modifying legend objects
        if candidate in lsettings:
            hsettings[candidate] = lsettings.pop(candidate)
    for candidate,translate in [('fontsize','size'),('fontweight','weight'),('va','va'),('ha','ha')]:
        if candidate in lsettings:
            tsettings[translate] = lsettings.pop(candidate)

    # Detect if user wants to specify rows manually
    # Gives huge latitude for user input:
    #   1) user can specify nothing and align will be inferred (list of iterables
    #      will always be False, i.e. we draw consecutive legends, and list of handles is always true)
    #   2) user can specify align (needs list of handles for True, list of handles or list
    #      of iterables for False and if the former, will turn into list of iterables)
    list_of_lists = not isinstance(handles[0], martist.Artist)
    if align is None: # automatically guess
        align = not list_of_lists
    else: # standardize format based on input
        if not align and not list_of_lists: # separate into columns
            # raise ValueError("Need to specify number of columns with ncol.")
            list_of_lists = True
            lsettings['ncol'] = lsettings.get('ncol',3)
            handles = [handles[i*lsettings['ncol']:(i+1)*lsettings['ncol']]
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
            lsettings['ncol'] = len(handles[0]) # choose this for column length
        elif 'ncol' not in lsettings:
            lsettings['ncol'] = 3
        # Split up into rows and columns -- by default matplotlib will
        # sort them in ***column-major*** order but that's dumb, we want row-major!
        # See: https://stackoverflow.com/q/10101141/4970632
        newhandles = []
        ncol = lsettings['ncol'] # number of columns
        handlesplit = [handles[i*ncol:(i+1)*ncol] for i in range(len(handles)//ncol+1)] # split into rows
        nrowsmax, nfinalrow = len(handlesplit), len(handlesplit[-1]) # max possible row count, and columns in final row
        nrows = [nrowsmax]*nfinalrow + [nrowsmax-1]*(lsettings['ncol']-nfinalrow)
            # e.g. if 5 columns, but final row length 3, columns 0-2 have N rows but 3-4 have N-1 rows
        for col,nrow in enumerate(nrows): # iterate through cols
            newhandles.extend(handlesplit[row][col] for row in range(nrow))
        handles = newhandles
        # Finally draw legend, mimicking row-major ordering
        leg = super(Axes, ax).legend(handles=handles, **lsettings)
        # Format handles, text, and legend patch
        leg.legendPatch.update(rc('outline')) # or get_frame()
        for obj in leg.legendHandles:
            obj.update(hsettings)
        for t in leg.texts:
            t.update(tsettings) # or get_texts()
        legends = leg
    # 2) Separate legend for each row
    # The labelspacing/borderspacing will be exactly replicated, as if we were
    # using the original legend command
    # Means we also have to overhaul some settings
    else:
        legends = []
        if 'labelspacing' not in lsettings:
            raise ValueError('Need to know labelspacing before drawing unaligned-column legends. Add it to settings.')
        if 'size' not in tsettings:
            raise ValueError('Need to know fontsize before drawing unaligned-column legends. Add it to settings.')
        for override in ['loc','ncol','bbox_to_anchor','borderpad','borderaxespad','frameon','framealpha']:
            lsettings.pop(override, None)
        # Determine space we want sub-legend to occupy, as fraction of height
        # Don't normally save "height" and "width" of axes so keep here
        interval = 1/len(handles) # split up axes
        interval = ((tsettings['size'] + lsettings['labelspacing']*tsettings['size'])/72) / \
                (ax.figure.get_figheight() * np.diff(ax._position.intervaly))
        # Iterate and draw
        for h,hs in enumerate(handles):
            bbox = mtransforms.Bbox([[0,1-(h+1)*interval],[1,1-h*interval]])
            leg = super(Axes, ax).legend(handles=hs, ncol=len(hs),
                bbox_to_anchor=bbox,
                frameon=False, borderpad=0, loc='center', **lsettings) # _format_legend is overriding original legend Method
            leg.legendPatch.update(rc('outline')) # or get_frame()
            for obj in leg.legendHandles:
                obj.update(hsettings)
            for t in leg.texts:
                t.update(tsettings) # or get_texts()
            legends.append(leg)
        for l in legends[:-1]:
            ax.add_artist(l) # because matplotlib deletes previous ones...
    return legends

def colorbar_factory(ax, mappable, cgrid=False, clocator=None,
        ctickminor=False, cminorlocator=None, cformatter=None, clabel=None,
        errfix=True, extend='neither', extendlength=0.2, # in inches
        values=None, orientation='horizontal', ticklocation='outer', **kwargs): #, settings=None):
    """
    Function for formatting colorbar-axes (axes that are "filled" by a colorbar).
    * There are options on the colorbar object (cb.locator, cb.formatter with cb.update_ticks)
    and by passing kwargs (ticks=x, format=y) that allow uer to not reference the underlying
    "axes" when fixing ticks. Don't use this functionality because not necessary for us and
    is missing many features, e.g. minorlocators/minorformatters. Also is different syntax.
    * There is an INSANELY WEIRD problem with colorbars when simultaneously passing levels
    and norm object to a mappable; fixed by passing vmin/vmax INSTEAD OF levels 
    (see: https://stackoverflow.com/q/40116968/4970632).
    * Problem is, often WANT levels instead of vmin/vmax, while simultaneously
    using a Normalize (for example) to determine colors between the levels
    (see: https://stackoverflow.com/q/42723538/4970632).
    * Workaround is to make sure locators are in vmin/vmax range exclusively;
    cannot match/exceed values.
    * The 'extend' kwarg is used for the case when you are manufacturing colorbar
    from list of colors or lines. Most of the time want 'neither'.
    TODO: Issue appears where the normalization vmin/vmax are outside of explicitly declared "levels"
    minima and maxima but that is probaby appropriate. If your levels are all within vmin/vmax,
    you will get discrete jumps outside of range and the extendlength at ends of colorbars will be weird.
    """
    # Make sure to eliminate ticks
    # cax.xaxis.set_tick_params(which='both', bottom=False, top=False)
    # cax.yaxis.set_tick_params(which='both', bottom=False, top=False)
    # Test if we were given a mappable, or iterable of stuff; note Container and
    # PolyCollection matplotlib classes are iterable.
    fromlines, fromcolors = False, False
    if not isinstance(mappable, martist.Artist) and not isinstance(mappable, mcontour.ContourSet):
        if isinstance(mappable[0], martist.Artist):
            fromlines = True # we passed a bunch of line handles; just use their colors
        else:
            fromcolors = True # we passed a bunch of color strings or tuples
    csettings = {'cax':ax, 'orientation':orientation, 'use_gridspec':True, # use space afforded by entire axes
                 'spacing':'uniform', 'extend':extend, 'drawedges':cgrid} # this is default case unless mappable has special props
    # Update with user-kwargs
    csettings.update(**kwargs)
    if hasattr(mappable, 'extend') and mappable.extend is not None:
        csettings.update({'extend':mappable.extend})

    # Option to generate colorbar/colormap from line handles
    # * Note the colors are perfect if we don't extend them by dummy color on either side,
    #   but for some reason labels for edge colors appear offset from everything
    # * Too tired to figure out why so just use this workaround
    # TODO TODO TODO: This should be abstracted away into the colormap
    # retrieval routines in colortools.py
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
    if fromlines or fromcolors:
        # colors = ['#ffffff'] + colors + ['#ffffff']
        values = np.array(values) # needed for below
        # colors = colors[:1] + colors + colors[-1:]
        colormap = mcolors.LinearSegmentedColormap.from_list('tmp', colors, N=len(colors)) # note that
            # the 'N' is critical; default 'N' is otherwise 256, and can get weird artifacts due to
            # unintentionally sampling some 'new' colormap colors; very bad!
        levels = utils.edges(values) # get "edge" values between centers desired
        values = utils.edges(values) # this works empirically; otherwise get weird situation with edge-labels appearing on either side
        mappable = plt.contourf([[0,0],[0,0]], levels=levels, cmap=colormap,
                extend='neither', norm=colortools.DiscreteNorm(values)) # workaround
        if clocator is None: # in this case, easy to assume user wants to label each value
            clocator = values
    if clocator is None:
        clocator = getattr(mappable.norm, 'levels', None)

    # Determine major formatters and major/minor tick locators
    # Can pass clocator/cminorlocator as the *jump values* between the mappables
    # vmin/vmax if desired
    normfix = False # whether we need to modify the norm object
    locators = [] # put them here
    for i,locator in enumerate((clocator,cminorlocator)):
        # Get the locator values
        # Need to use tick_values instead of accessing 'locs' attribute because
        # many locators don't have these attributes; require norm.vmin/vmax as input
        if i==1 and not ctickminor and locator is None: # means we never wanted minor ticks
            locators.append(Locator('null'))
            continue
        values = np.array(Locator(locator).tick_values(mappable.norm.vmin, mappable.norm.vmax)) # get the current values
        # Modify ticks to work around mysterious error, and to prevent annoyance
        # where minor ticks extend beyond extendlength.
        # We need to figure out the numbers that will eventually be rendered to
        # solve the error, so we will always use a fixedlocator.
        values_min = np.where(values>=mappable.norm.vmin)[0]
        values_max = np.where(values<=mappable.norm.vmax)[0]
        if len(values_min)==0 or len(values_max)==0:
            raise ValueError(f'No ticks are within the colorbar range {mappable.norm.vmin:.3g} to {mappable.norm.vmax:.3g}.')
        values_min, values_max = values_min[0], values_max[-1]
        values = values[values_min:values_max+1]
        if values[0]==mappable.norm.vmin:
            normfix = True
        if i==1: # prevent annoying major/minor overlaps where one is slightly shifted left/right
            values = [v for v in values if not any(v>=o-abs(v)/1000 and v<=o+abs(v)/1000 for o in fixed)] # floating point weirdness is fixed below
        fixed = values # record as new variable
        locators.append(Locator(fixed)) # final locator object
    # Next the formatter
    formatter = Formatter(cformatter)

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
    if hasattr(mappable.norm,'levels'):
        mappable.norm.levels = np.atleast_1d(mappable.norm.levels).astype(np.float)
        if normfix:
            mappable.norm.levels[0] -= np.diff(mappable.norm.levels[:2])[0]/1000

    # Draw the colorbar
    # TODO: For whatever reason the only way to avoid bugs seems to be to pass
    # the major formatter/locator to colorbar commmand and directly edit the
    # minor locators/formatters; update_ticks after the fact ignores the major formatter
    if orientation=='horizontal':
        axis = ax.xaxis
        scale = ax.figure.width*np.diff(getattr(ax.get_position(),'intervalx'))[0]
    else:
        axis = ax.yaxis
        scale = ax.figure.height*np.diff(getattr(ax.get_position(),'intervaly'))[0]
    extendlength = extendlength/(scale - 2*extendlength)
    csettings.update({'extendfrac':extendlength}) # width of colorbar axes and stuff
    # cb = ax.figure.colorbar(mappable, **csettings)
    cb = ax.figure.colorbar(mappable, ticks=locators[0], ticklocation=ticklocation,
            format=cformatter, **csettings)

    # The ticks/ticklabels basic properties
    for t in axis.get_ticklabels(which='major'):
        t.update(rc('ticklabels'))
    axis.set_tick_params(which='major', **rc('ctick'))
    axis.set_tick_params(which='minor', **rc('ctickminor'))
        # properties are obscure; would have to use hidden methods to do this manually; 
        # using update method (inhereted from Artist) doesn't work
    # The minor locators and formatters
    # * The minor locator must be set with set_ticks after transforming an array
    #   using the mappable norm object; see: https://stackoverflow.com/a/20079644/4970632
    # * The set_minor_locator seems to be completely ignored depending on the colorbar
    #   in question, for whatever reason
    # * The major locator and formatter settings here are also not ideal since we'd have to
    #   update_ticks which might throw off the minor ticks again
    # axis.set_major_locator(locators[0])
    # axis.set_major_formatter(cformatter)
    # axis.set_ticks(locators[1].locs, minor=True)
    # axis.set_minor_locator(locators[1])
    axis.set_ticks(mappable.norm(np.array(locators[1].tick_values(mappable.norm.vmin, mappable.norm.vmax))), minor=True)
    axis.set_minor_formatter(mticker.NullFormatter()) # to make sure
    # Set up the label
    axis.label.update(rc('label'))
    if clabel is not None:
        axis.label.update({'text':clabel})
    # Fix pesky white lines between levels + misalignment with border due to rasterized blocks
    cb.solids.set_edgecolor('face')
    cb.solids.set_rasterized(False)
    # Make edges/dividers consistent with axis edges
    if cb.dividers is not None:
        cb.dividers.update(rc('cgrid'))
    cb.outline.update(rc('outline'))
    # Update the tickers colorbar
    # cb.formatter = cformatter # fucks shit up -- why?
    # cb.update_ticks() # doesn't actually update the formatter
    # axis.set_major_formatter(cformatter) # fucks shit up -- why?
    return cb

#------------------------------------------------------------------------------#
# Custom settings for various journals
# Add to this throughout your career, or as standards change
# PNAS info: http://www.pnas.org/page/authors/submission
# AMS info: https://www.ametsoc.org/ams/index.cfm/publications/authors/journal-and-bams-authors/figure-information-for-authors/
# AGU info: https://publications.agu.org/author-resource-center/figures-faq/
#------------------------------------------------------------------------------#
def journalsize(width, height):
    # User wants to define their own size
    if type(width) is not str:
        return width, height
    # Determine automatically
    cm2in = 0.3937
    mm2in = 0.03937
    table = {
        'pnas1': 8.7*cm2in,
        'pnas2': 11.4*cm2in,
        'pnas3': 17.8*cm2in,
        'ams1': 3.2,
        'ams2': 4.5,
        'ams3': 5.5,
        'ams4': 6.5,
        'agu1': (95*mm2in, 115*mm2in),
        'agu2': (190*mm2in, 115*mm2in),
        'agu3': (95*mm2in, 230*mm2in),
        'agu4': (190*mm2in, 230*mm2in),
        }
    value = table.get(width, None)
    if value is None:
        raise ValueError(f'Unknown journal figure width specifier "{width}". ' +
                          'Options are: ' + ', '.join(table.keys()))
    # Output width, and optionally also specify height
    if not utils.isiterable(value):
        width = value
    else:
        width, height = value
    return width, height

#-------------------------------------------------------------------------------
# Primary plotting function; must be used to create figure/axes if user wants
# to use the other features
#-------------------------------------------------------------------------------
def subplots(array=None, rowmajor=True, # mode 1: specify everything with array
        tight=False,  # whether to set up tight bbox from gridspec object
        silent=False, # how much stuff to print
        nrows=1, ncols=1, emptycols=None, emptyrows=None, # mode 2: use convenient kwargs for simple grids
        sharex=True, sharey=True, # for sharing x/y axis limits/scales/locators for axes with matching GridSpec extents, and making ticklabels/labels invisible
        spanx=True, spany=True,   # custom setting, optionally share axis labels for axes with same xmin/ymin extents
        aspect=1, height=None, width=None,                  # for controlling aspect ratio, default is control for width
        hspace=0.4, wspace=0.2, hratios=None, wratios=None, # spacing between axes, in inches (hspace should be bigger, allowed room for title)
        left=0.5, bottom=0.5, right=0.15, top=0.3,          # spaces around edge of main plotting area, in inches
        bottompanel=False, bottompanels=False, rightpanel=False, rightpanels=False, bottompanelrows=1, rightpanelcols=1, # optionally draw extra rows
        bottomcolorbar=False, bottomcolorbars=False, rightcolorbar=False, rightcolorbars=False, bottomlegend=False, bottomlegends=False,
        bwidth=0.17, bspace=0.5, rwidth=0.17, rspace=0.5, # default to no space between panels
        innerpanels=None, whichpanels=None, # same as below; list of numbers where we want subplotspecs
        ihspace=None, iwspace=None, ihwidth=None, iwwidth=None,
        maps=None, # set maps to True for maps everywhere, or to list of numbers
        package='basemap', projection=None, projection_dict={}, **projection_kwargs): # for projections; can be 'basemap' or 'cartopy'
    """
    Special creation of subplots grids, allowing for arbitrarily overlapping 
    axes objects. Will return figure handle and axes objects. Need to finish
    the documentation of this; extremely feature rich.
    * Easiest way to create subplots is with nrows=1 and ncols=1. If you want extra space
      between a row or column, specify the row/column number that you want to be 'empty' with
      emptyrows=row/emptycolumn=column, and adjust wratios/hratios for the desired width of that space.
    * For more complicated plots, can pass e.g. array=[[1,2,3,4],[0,5,5,0]] to create a grid
      of 4 plots on the top, single plot spanning the middle 2-columns on the bottom, and empty
      spaces where the 0 appears.
    * Use bottompanel/bottompanels to make several or multiple panels on the bottom
      that can be populated with multiple colorbars/legend; bottompanels=True will
      just make one 'space' for every column, and bottompanels=[1,1,2] for example will
      make a panel spanning the first two columns, then a single panel for the final column.
      This will add a bottompanel attribute to the figure; can index that attribute if there
      are multiple places for colorbars/legend.
    * Use rightpanel/rightpanels in the same way.
    * Create extra panels *within* a grid of subplots (e.g. a 2x2 grid, each of which has a teeny
      panel attached) innerpanels=True, and optionally filter to subplot numbers with whichpanels=[list];
      then use ihspace/iwspace/ihwidth/iwwidth to control the separation and widths of the subpanels.
    * Initialize cartopy plots with package='basemap' or package='cartopy'. Can control which plots
      we want to be maps with maps=True (everything) or maps=[numbers] (the specified subplot numbers).
    Notes:
        * Matplotlib set_aspect option seems to behave strangely on some plots (trend-plots from
        SST paper); for this reason we override the fix_aspect option provided by basemap and
        just draw figure with appropriate aspect ratio to begin with. Otherwise get weird
        differently-shaped subplots that seem to make no sense.
        * Shared axes will generally end up with the same axis limits/scaling/majorlocators/minorlocators;
        the sharex and sharey detection algorithm really is just to get instructions to make the
        ticklabels/axis labels invisible for certain axes.
        * We create the cartopy axes on.
    Todo:
        * For spanning axes labels, right now only detect **x labels on bottom**
        and **ylabels on top**; generalize for all subplot edges.
        * Figure size should be constrained by the dimensions of the axes, not vice
        versa; might make things easier.
    """
    # Translate some arguments
    # This is convenience feature to get more sensible default spacing
    if rightcolorbar or rightcolorbars:
        rwidth = 0.17
        rspace = 0.5
        rightpanel = rightcolorbar
        rightpanels = rightcolorbars
    if bottomcolorbar or bottomcolorbars:
        bwidth = 0.17
        bspace = 0.5
        bottompanel = bottomcolorbar
        bottompanels = bottomcolorbars
    if bottomlegend or bottomlegends:
        bwidth = 0.25
        bspace = 0
        bottompanel = bottomlegend
        bottompanels = bottomlegends
    # Fill in some basic required settings
    if package not in ['basemap','cartopy']:
        raise ValueError("Plotting package must be one of basemap or cartopy.")
    width, height = journalsize(width, height) # if user passed width=<string>, will use that journal size
    if width is None and height is None: # at least one must have value
        width = 5
    if width is not None and height is not None: # specify exact size
        aspect = width/height
    if wratios is not None:
        aspect = aspect/(wratios[0]/np.mean(wratios)) # e.g. if 2 columns, 5:1 width ratio, change the 'average' aspect ratio
    if hratios is not None:
        aspect = aspect*(hratios[0]/np.mean(hratios))
    # Automatically generate array of first *arg not provided, or use nrows/ncols
    if array is None:
        array = np.arange(1,nrows*ncols+1)[...,None]
        array = array.reshape((nrows, ncols)) # numpy is row-major, remember
    array = np.array(array) # enforce array type
    if array.ndim==1:
        array = array[None,:] if rowmajor else array[:,None] # interpret as single row or column
    elif not rowmajor:
        array = np.reshape(array.flatten(), (array.shape[0],array.shape[1]), order='F') # make column major
    array[array==None] = 0 # use zero for placeholder; otherwise have issues
    # Include functionality to declare certian rows/columns empty, if user requested it
    # Useful if user wants to include extra space between some columns greater than the
    # spaces between axes in other columns
    if emptycols is not None:
        if not utils.isiterable(emptycols):
            emptycols = [emptycols]
        for col in emptycols:
            array[:,col-1] = 0
    if emptyrows is not None:
        if not utils.isiterable(emptyrows):
            emptyrows = [emptyrows]
        for row in emptyrows:
            array[row-1,:] = 0
    nrows = array.shape[0]
    ncols = array.shape[1]
    # # Enforce consistent numbering; row-major increasing from 1 every time
    # # a new axes is encountered; e.g. [[1,2],[1,3],[1,4],[1,4]]
    # # Maybe ignore for now, but this makes sure the output axes list is
    # # always in same row-major order, not random user-input order
    # number = 1
    # newarray = np.zeros(array.shape)
    # for row in newarray.shape[0]:
    #     for col in newarray.shape[1]:
    #         if array[row,col] not in array.flat: # not already declared
    #             newarray[array==array[row,col]] = number
    #             number += 1

    #--------------------------------------------------------------------------
    # Projection setup
    #--------------------------------------------------------------------------
    # Get basemap.Basemap or cartopy.CRS instances
    map_kwargs = {}
    if maps and (wratios is not None or hratios is not None):
        raise NotImplementedError('Not yet possible.')
    if not maps and projection_kwargs:
        raise ValueError(f'Unknown kwargs: {", ".join(projection_kwargs.keys())}. If you want to create maps, you must change the "maps" kwarg.')
    if maps:
        map_projection, aspect = projection_factory(package, projection, **{**projection_kwargs, **projection_dict})
        map_kwargs = {'projection':package, 'map_projection':map_projection}

    #--------------------------------------------------------------------------
    # Panel considerations; keep things general for future imporvements
    #--------------------------------------------------------------------------
    # Spacings due to panel considerations
    if bottompanel or bottompanels:
        bottom_extra_axes = bwidth
        bottom_extra_hspace = bspace
    else:
        bottom_extra_axes, bottom_extra_hspace = 0, 0
    bottom_extra = bottom_extra_axes + bottom_extra_hspace

    # And right spacings
    if rightpanel or rightpanels:
        right_extra_axes   = rwidth
        right_extra_wspace = rspace
    else:
        right_extra_axes, right_extra_wspace = 0, 0
    right_extra = right_extra_axes + right_extra_wspace

    #--------------------------------------------------------------------------
    # Apply aspect ratio to axes, infer hspaces/wspaces/bottom/top/left/right
    #--------------------------------------------------------------------------
    # Fixed aspect, infer required figure size [try wspace=0.2 for gs with cols=4, axes gs[0,:2] and gs[0,2:] 
    # vs. gs with cols=2, axes gs[0,0] and gs[0,1], and spacing should differ])
    # FORMULA: aspect = ((width - left - right - (ncol-1)*wspace)/ncol) / ((height - top - bottom - (nrow-1)*hspace)/nrow)
    # Fixing height, determine figure aspect ratio
    if height is not None:
        axheight_ave_nopanel = (height - top - bottom - (nrows-1)*hspace - bottom_extra)/nrows
        if axheight_ave_nopanel<0:
            raise ValueError('Not enough room for axes. Reduce left/bottom/wspace.')
        axwidth_ave_nopanel = axheight_ave_nopanel*aspect
        width               = axwidth_ave_nopanel*ncols + left + right + (ncols-1)*wspace + right_extra

    # Fixing width, determine aspect ratio
    else:
        # print(width, left, right, ncols, wspace, right_extra, ncols)
        # print(width, left, right, ncols, wspace, right_extra)
        axwidth_ave_nopanel = (width - left - right - (ncols-1)*wspace - right_extra)/ncols
        if axwidth_ave_nopanel<0:
            raise ValueError('Not enough room for axes. Reduce left/bottom/wspace.')
        axheight_ave_nopanel = axwidth_ave_nopanel/aspect
        height               = axheight_ave_nopanel*nrows + top + bottom + (nrows-1)*hspace + bottom_extra

    # Figure size, and some other stuff
    # The "total" axes width include hspace/wspace 
    axwidth_total  = ncols*axwidth_ave_nopanel + (ncols-1)*wspace
    axheight_total = nrows*axheight_ave_nopanel + (nrows-1)*hspace
    figsize        = (width, height)

    # Properties for outer GridSpec object, need borders in fractional units
    ileft, ibottom, iright, itop = left, bottom, right, top
    left = left/width
    top  = 1-top/height
    if bottom_extra_axes:
        hspace_outer        = bottom/((axheight_total + bottom_extra_axes)/2) # same for main axes and bottom panel, but 'bottom'
        height_ratios_outer = np.array([axheight_total, bottom_extra_axes])/(axheight_total + bottom_extra_axes)
        bottom              = bottom_extra_hspace/height
    else:
        hspace_outer        = 0
        height_ratios_outer = np.array([1])
        bottom              = bottom/height
    if right_extra_axes:
        width_ratios_outer = np.array([axwidth_total, right_extra_axes])/(axwidth_total + right_extra_axes)
        wspace_outer       = right/((axwidth_total + right_extra_axes)/2) # want space between main axes and panel to be 'right'
        right              = 1-right_extra_wspace/width
    else:
        width_ratios_outer = np.array([1])
        wspace_outer       = 0
        right              = 1-right/width
    # Properties for inner GridSpec object
    wspace_orig = wspace
    hspace_orig = hspace
    wspace = wspace/axwidth_ave_nopanel
    hspace = hspace/axheight_ave_nopanel
    if wratios is not None:
        wratios = np.array(wratios)/sum(wratios)
    else:
        wratios = np.ones(ncols)/ncols
    if hratios is not None:
        hratios = np.array(hratios)/sum(hratios)
    else:
        hratios = np.ones(nrows)/nrows

    # Create gridspec for outer plotting regions (divides 'main area' from side panels)
    GS = mgridspec.GridSpec(
            nrows         = 1+int(bottom_extra_axes>0),
            ncols         = 1+int(right_extra_axes>0),
            left          = left,
            bottom        = bottom,
            right         = right, # unique spacing considerations
            top           = top, # so far no panels allowed here
            wspace        = wspace_outer,
            hspace        = hspace_outer,
            width_ratios  = width_ratios_outer,
            height_ratios = height_ratios_outer,
            ) # set wspace/hspace to match the top/bottom spaces
    # Create axes within the 'main' plotting area
    # Will draw individual axes using this GridSpec object later
    gs = mgridspec.GridSpecFromSubplotSpec(
            nrows         = nrows,
            ncols         = ncols,
            subplot_spec  = GS[0,0],
            wspace        = wspace,
            hspace        = hspace,
            width_ratios  = wratios,
            height_ratios = hratios,
            )
    # Create figure
    fig = plt.figure(figsize=figsize, FigureClass=Figure,
        left=ileft,    bottom=ibottom, right=iright,  top=itop,
        bwidth=bwidth, bspace=bspace,  rwidth=rwidth, rspace=rspace,
        height=height, width=width,    gridspec=GS
        )

    #--------------------------------------------------------------------------
    # Manage shared axes/axes with spanning labels
    #--------------------------------------------------------------------------
    # Find shared axes; get sets with identical spans in x or y (ignore panels)
    # Preliminary stuff
    axes_ids = [np.where(array==i) for i in np.unique(array) if i>0] # 0 stands for empty, -1 for colorbar, -2 for legend
    # note that these locations should be **sorted** by axes id
    yrange = np.array([[xy[0].min(), xy[0].max()+1] for xy in axes_ids]) # yrange is shared columns
    xrange = np.array([[xy[1].min(), xy[1].max()+1] for xy in axes_ids])
    xmin   = np.array([xy[0].min() for xy in axes_ids])
    xmax   = np.array([xy[0].max() for xy in axes_ids])
    ymin   = np.array([xy[1].min() for xy in axes_ids])
    ymax   = np.array([xy[1].max() for xy in axes_ids])
    num_axes = len(axes_ids)
    # Find axes that use map projections
    # TODO: Implement this, not currently used! Implementation will be similar to
    # innerpanels_ids below.
    if maps is not None:
        if maps is True:
            maps = [*np.unique(array)]
        if not utils.isiterable(maps):
            maps = [maps] # just want a single map
        else:
            maps = [*maps] # force into list, not array
        maps_ids = [i for i,a in enumerate(np.unique(array).flat) if a in maps]
    else:
        maps_ids = []
    # Find axes that have inner panels
    panel_kwargs = {'whichpanels':whichpanels,
        'hspace':ihspace, 'wspace':iwspace,
        'hwidth':ihwidth, 'wwidth':iwwidth,
        }
    if innerpanels is not None:
        if innerpanels is True:
            innerpanels = [*np.unique(array)]
        if not utils.isiterable(innerpanels):
            innerpanels = [innerpanels]
        else:
            innerpanels = [*innerpanels] # force into list, not array
        innerpanels_ids = [i for i,a in enumerate(np.unique(array)) if a in innerpanels]
    else:
        innerpanels_ids = []
    # Find pairs with edges on same gridspec
    if spanx:
        xgroups_span_base, xgroups_span, grouped = [], [], []
        for i in range(num_axes):
            matching_axes = np.where(xmin[i]==xmin)[0]
            if i not in grouped and matching_axes.size>1:
                xgroups_span      += [matching_axes] # add ndarray of ids to list
                xgroups_span_base += [matching_axes[np.argmin(xmin[matching_axes])]]
                    # add the axes that is farthest down
            grouped += [*matching_axes] # bookkeeping; record ids that have been grouped already
    if spany:
        ygroups_span_base, ygroups_span, grouped = [], [], []
        for i in range(num_axes):
            matching_axes = np.where(ymax[i]==ymax)[0]
            if i not in grouped and matching_axes.size>1:
                ygroups_span      += [matching_axes] # add ndarray of ids to list
                ygroups_span_base += [matching_axes[np.argmax(ymax[matching_axes])]]
                    # add the axes that is farthest left
            grouped += [*matching_axes] # bookkeeping; record ids that have been grouped already
    # Get shared x axes
    if sharex:
        xgroups_base, xgroups, grouped = [], [], []
        for i in range(num_axes): # axes now have pseudo-numbers from 0 to num_axes-1
            matches       = (xrange[i:i+1,:]==xrange).all(axis=1) # broadcasting rules work here
            matching_axes = np.where(matches)[0] # gives ID NUMBER of matching_axes, from 0 to num_axes-1
            if i not in grouped and matching_axes.size>1:
                xgroups      += [matching_axes] # add ndarray of ids to list
                xgroups_base += [matching_axes[np.argmax(yrange[matching_axes,1])]]
                    # bottom-most axis with shared x; should be single number
            grouped += [*matching_axes] # bookkeeping; record ids that have been grouped already
    # Get shared y axes
    if sharey:
        ygroups_base, ygroups, grouped = [], [], []
        for i in range(num_axes):
            matches       = (yrange[i:i+1,:]==yrange).all(axis=1) # broadcasting rules work here
            matching_axes = np.where(matches)[0]
            if i not in grouped and matching_axes.size>1:
                ygroups      += [matching_axes] # add ndarray of ids to list
                ygroups_base += [matching_axes[np.argmin(xrange[matching_axes,0])]] # left-most axis with shared y, for matching_axes
            grouped += [*matching_axes] # bookkeeping; record ids that have been grouped already

    #--------------------------------------------------------------------------
    # Draw axes
    #--------------------------------------------------------------------------
    # Base axes; to be shared with other axes as ._sharex, ._sharey attributes
    axs = num_axes*[None] # list of axes
    allgroups_base = []
    if sharex:
        allgroups_base += xgroups_base
    if sharey:
        allgroups_base += ygroups_base
    for i in allgroups_base: # this is just list of indices in axes_ids, yrange, etc.
        ax_kwargs = map_kwargs if i in maps_ids else {'projection':'xy'}
        if axs[i] is None: # not created as a x-base already, for example
            if i in innerpanels_ids:
                axs[i] = fig.panel_factory(gs[slice(*yrange[i,:]), slice(*xrange[i,:])], width, height,
                        **ax_kwargs, **panel_kwargs) # main axes handle
            else:
                axs[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        **ax_kwargs) # main axes can be a cartopy projection

    # Dependent axes
    for i in range(num_axes):
        # Detect if we want to share this axis with another
        # If so, get that axis
        sharex_ax, sharey_ax = None, None # by default, don't share with other axes objects
        ax_kwargs = map_kwargs if i in maps_ids else {'projection':'xy'}
        if sharex:
            igroup = np.where([i in g for g in xgroups])[0] # np.where works on lists
            if igroup.size==1:
                sharex_ax = axs[xgroups_base[igroup[0]]]
                if sharex_ax is None:
                    raise ValueError('Something went wrong; shared x axes was not already drawn.')
            elif igroup.size>1:
                raise ValueError(f'Something went wrong; axis {i:d} belongs to multiple groups.')
        if sharey:
            igroup = np.where([i in g for g in ygroups])[0] # np.where works on lists
            if igroup.size==1:
                sharey_ax = axs[ygroups_base[igroup[0]]]
                if sharey_ax is None:
                    raise ValueError('Something went wrong; shared x axes was not already drawn.')
            elif igroup.size>1:
                raise ValueError(f'Something went wrong; axis {i:d} belongs to multiple groups.')

        # Draw axes, and add to list
        if axs[i] is not None:
            # Axes is a BASE and has already been drawn, but might still need to add
            # a _sharex or _sharey property; e.g. is bottom-axes of three-column plot
            # and we want it to share a y-axis
            if axs[i] is not sharex_ax:
                axs[i]._sharex = sharex_ax
            else:
                sharex_ax = None
            if axs[i] is not sharey_ax:
                axs[i]._sharey = sharey_ax
            else:
                sharey_ax = None
        else:
            # Virgin axes; these are not an x base or a y base
            # We draw them now and pass the sharex/sharey as kwargs
            if i in innerpanels_ids:
                axs[i] = fig.panel_factory(gs[slice(*yrange[i,:]), slice(*xrange[i,:])], width, height,
                        sharex=sharex_ax, sharey=sharey_ax, **ax_kwargs, **panel_kwargs)
            else:
                axs[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        sharex=sharex_ax, sharey=sharey_ax, **ax_kwargs) # main axes can be a cartopy projection

        # Hide tick labels (not default behavior for manual sharex, sharey use)
        if sharex_ax is not None:
            for t in axs[i].xaxis.get_ticklabels(): t.set_visible(False)
            axs[i].xaxis.label.set_visible(False)
        if sharey_ax is not None:
            for t in axs[i].yaxis.get_ticklabels(): t.set_visible(False)
            axs[i].yaxis.label.set_visible(False)

    # Spanning axes; allow xlabels/ylabels to span them
    if spanx and len(xgroups_span)>0:
        for g, b in zip(xgroups_span, xgroups_span_base):
            # Specify x, y transform in Figure coordinates
            axs[b].xaxis.label.set_transform(mtransforms.blended_transform_factory(
                    fig.transFigure, mtransforms.IdentityTransform()
                    ))
            # Get min/max positions, in figure coordinates, of spanning axes
            xmin = min(axs[i].get_position().xmin for i in g)
            xmax = max(axs[i].get_position().xmax for i in g)
            # This is the shared xlabel
            # print('Group:', g, 'Base:', b, 'Span:', xmin, xmax)
            axs[b].xaxis.label.set_position(((xmin+xmax)/2, 0))
            for i in g:
                if i!=b: axs[i].xaxis.label.set_visible(False)
    if spany and len(ygroups_span)>0:
        for g, b in zip(ygroups_span, ygroups_span_base):
            axs[b].yaxis.label.set_transform(mtransforms.blended_transform_factory(
                    mtransforms.IdentityTransform(), fig.transFigure # specify x, y transform
                    ))
            ymin = min(axs[i].get_position().ymin for i in g)
            ymax = max(axs[i].get_position().ymax for i in g)
            # print('Group:', g, 'Base:', b, 'Span:', ymin, ymax)
            axs[b].yaxis.label.set_position((0, (ymin+ymax)/2))
                # this is the shared ylabel
            for i in g:
                if i!=b: axs[i].yaxis.label.set_visible(False)

    #---------------------------------------------------------------------------
    # Create panel axes
    # TODO: Fix lformat and cformat methods to apply to SubplotSpec objects
    # so they can be drawn immediately; can pass the figure handle needed to draw
    # the axes with add_subplot using a defaulted kwarg fig=fig
    #---------------------------------------------------------------------------
    # First the bottompanel options
    if bottompanel: # one spanning panel
        bottompanels = [1]*ncols
    elif bottompanels: # more flexible spec
        if not utils.isiterable(bottompanels):
            bottompanels = [*range(ncols)] # pass True to make panel for each column
    if bottompanels: # non-empty
        # First get the number of columns requested and their width 
        # ratios to align with main plotting area
        num = bottompanels[0]
        npanel = len(bottompanels)
        widths_orig = [axwidth_ave_nopanel*ncols*wratio for wratio in wratios] # assumes wratios normalized
        widths_new = [widths_orig[0]] # start with this
        for p,panel in enumerate(bottompanels):
            if p==0:
                continue
            newnum = panel
            if newnum==num: # same as previous number
                widths_new[-1] += widths_orig[p] + wspace_orig
                npanel -= 1 # we want one less panel object
            else:
                widths_new += [widths_orig[p]] # add to width
            num = newnum
        rratios = [width/sum(widths_new) for width in widths_new] # just normalize it
        rspace = wspace_orig/(sum(widths_new)/npanel) # divide by new average width
        # Now create new GridSpec and add each of its
        # SubplotSpecs to the figure instance
        P = mgridspec.GridSpecFromSubplotSpec(
                nrows         = bottompanelrows,
                ncols         = npanel,
                subplot_spec  = GS[1,0],
                wspace        = rspace, # same as above
                width_ratios  = rratios,
                )
        panels = [Panel(fig, P[0,i], 'bottom') for i in range(npanel)]
        if bottompanel:
            panels = panels[0] # no indexing if user specified single panel, but does mean indexing if specified bottompanels=True with single column grid
        fig.bottompanel = panels

    # Next the rightpanel options (very similar)
    if rightpanel:
        rightpanels = [1]*nrows
    elif rightpanels:
        if not utils.isiterable(rightpanels):
            rightpanels = [*range(nrows)] # pass True to make panel for each row
    if rightpanels: # non-empty
        # First get the number of columns requested and their width 
        # ratios to align with main plotting area
        num = rightpanels[0]
        npanel = len(rightpanels)
        heights_orig = [axheight_ave_nopanel*nrows*hratio for hratio in hratios] # assumes wratios normalized
        heights_new = [heights_orig[0]] # start with this
        for p,panel in enumerate(rightpanels):
            if p==0:
                continue
            newnum = panel
            if newnum==num: # same as previous number
                heights_new[-1] += heights_orig[p] + hspace_orig
                npanel -= 1 # we want one less panel object
            else:
                heights_new += [heights_orig[p]] # add to width
            num = newnum
        rratios = [height/sum(heights_new) for height in heights_new] # just normalize it
        rspace = hspace_orig/(sum(heights_new)/npanel) # divide by new average width
        # Now create new GridSpec and add each of its
        # SubplotSpecs to the figure instance
        P = mgridspec.GridSpecFromSubplotSpec(
                ncols         = rightpanelcols,
                nrows         = npanel,
                subplot_spec  = GS[0,1],
                hspace        = rspace,
                height_ratios = rratios,
                )
        panels = [Panel(fig, P[i,0], 'right') for i in range(npanel)] # pass the SubplotSpec objects
        if rightpanel:
            panels = panels[0]
        fig.rightpanel = panels

    #--------------------------------------------------------------------------
    # Return results
    # Will square singleton arrays
    #--------------------------------------------------------------------------
    if not silent:
        print('Figure setup complete.')
    for i,ax in enumerate(axs): # add this dynamically because it's easier
        ax.number = i+1 # know axes number ahead of time; start at 1
    if len(axs)==1:
        axs = axs[0]
    return fig, axs

#------------------------------------------------------------------------------#
# Miscellaneous
#------------------------------------------------------------------------------#
def figure(*args, **kwargs):
    """
    Simple alias for 'subplots', perhaps more intuitive.
    """
    return subplots(*args, **kwargs)
def close():
    """
    Close all figures 'open' in memory. This does not delete images printed
    in an ipython notebook; those are rendered versions of the abstract figure objects.
    """
    plt.close('all') # easy peasy
def show():
    """
    Show all figures.
    """
    plt.show()

