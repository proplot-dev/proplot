#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Imports, all
# Note that because of how matplotlib imports work (a module is only imported
# for a first time; import statements elsewhere in script during process just point
# to the already loaded module), settings rcParams here will change properties
# everywhere else until you start new instance
#------------------------------------------------------------------------------
# Global module dependencies
# from .basics import Dict, arange
import os
import numpy as np
# from string import ascii_uppercase
from string import ascii_lowercase
from types import MethodType
import matplotlib.figure as mfigure
# import matplotlib.axes as maxes
# import matplotlib.cm as mcm
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
import matplotlib.pyplot as plt
# Optionally import mapping toolboxes
# Main conda distro says they are incompatible, so make sure not required!
try:
    import mpl_toolkits.basemap as mbasemap
except ModuleNotFoundError:
    print("Warning: Basemap is not available.")
try: # crs stands for "coordinate reference system", leading c is "cartopy"
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.geoaxes import GeoAxes
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    Transform = ccrs.PlateCarree() # global variable
except ModuleNotFoundError:
    GeoAxes = None
    Transform = None # note, never pass transform=None! default is mtransforms.IdentityTransform()
    print("Warning: cartopy is not available.")
# Local modules, projection sand formatters and stuff
from .rc import globals
from .axis import Norm, Formatter, LatFormatter, LonFormatter, AutoLocate # default axis norm and formatter
from .proj import Aitoff, Hammer, KavrayskiyVII, WinkelTripel

#------------------------------------------------------------------------------
# Custom figure class
#------------------------------------------------------------------------------
def close():
    """
    Close all figures 'open' in memory. This does not delete images printed
    in an ipython notebook; those are rendered versions of the abstract figure objects.
    This prevents us having to load pyplot explicitly in ipython notebook.
    """
    plt.close('all') # easy peasy
class Figure(mfigure.Figure):
    """
    Subclass of the mfigure.Figure class, with lots of special formatting
    options. Can be called by using pyplot.figure(FigureClass=Figure) kwargument
    in my subplots function.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs) # python 3 only

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
            if self.rightcolorbar is not None:
                x2fix = self.cspace - x2
            else:
                x2fix = self.right  - x2
            if self.bottomlegend is not None:
                y1fix = self.lwidth - y1
            elif self.bottomcolorbar is not None:
                y1fix = self.cspace - y1
            elif self.bottompanel is not None:
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
def _graticule_fix(x, y):
    # Get graticule; want it to work with e.g. Gaussian-grid model output
    # dx, dy = x[1]-x[0], y[1]-y[0]
    # xb, yb = np.arange(x[0]-dx/2, x[-1]+dx, dx), np.arange(y[0]-dy/2, y[-1]+dy, dy)
    #     # can just use arange, because don't care about propagating float errors
    # return xb, yb
    edges = lambda vec: np.concatenate((vec[:1]-(vec[1]-vec[0])/2,
        (vec[1:]+vec[:-1])/2,
        vec[-1:]+(vec[-1]-vec[-2])/2))
    return edges(x), edges(y)
def _cartopy_seams(self, lon, lat, data):
    # The todo list for cartopy is much shorter, as we don't have to worry
    # about longitudes going from 0 to 360, -180 to 180, etc.; projection handles all that
    if isinstance(self, GeoAxes): # the main GeoAxes class; others like GeoAxesSubplot subclass this
        # Never need to fix 'seams' for pcolor/pcolormesh data!
        # So don't call this function from inside pcolor methods
        if lon.size!=data.shape[0]:
            raise ValueError('Longitude length does not match data length along dimension 0.')
        # print('Fixing cartopy seams!')
        # 1) Fix holes over poles by *interpolating* there (equivalent to
        # simple mean of highest/lowest latitude points)
        dataS = np.repeat(data[:,0].mean(), data.shape[0])[:,None]
        dataN = np.repeat(data[:,-1].mean(), data.shape[0])[:,None]
        lat = np.concatenate(([-90], lat, [90]))
        data = np.concatenate((dataS, data, dataN), axis=1)
        # 2) Fix seams at map boundary; this time the fancy projection will
        # handle issues with pcolor seams that we'd hvae with basemap, but still
        # have to ensure *circular* coverage if doing e.g. contourf
        lon = np.array((*lon, lon[0]+360)) # make longitudes circular
        data = np.concatenate((data, data[:1,:]), axis=0) # make data circular
    return lon, lat, data
def _m_seams(basemap, lon, lat, data):
    # Get some params
    lonmin, lonmax = basemap.lonmin, basemap.lonmax
    # Raise errors
    if lon.max()>lon.min()+360:
        raise ValueError('Longitudes must span 360 degrees at most.')
    if lon.min()<-360 or lon.max()>360:
        raise ValueError('Longitudes must fall in range [-360, 360].')
    if lonmin<-360 or lonmin>0:
        print("Warning: Minimum longitude not in range [-360,0].")
        # raise ValueError('Minimum longitude must fall in range [-360, 0].')
    # 1) Establish 360-degree range
    lon -= 720
    while True:
        filter_ = lon<lonmin
        if filter_.sum()==0:
            break
        lon[filter_] += 360
    # 2) Roll, accounting for whether ends are identical
    # If go from 0,1,...,359,0 (borders), returns id of first zero
    roll = -np.argmin(lon) # always returns FIRST value
    if lon[0]==lon[-1]:
        lon = np.roll(lon[:-1], roll)
        lon = np.append(lon, lon[0]+360)
    else:
        lon = np.roll(lon, roll)
    data = np.roll(data, roll, axis=0)
    # print("After roll 1:", lon)
    # 3) Roll in same direction some more, if some points on right-edge
    # extend more than 360 above the minimum longitude; THEY should be the
    # ones on west/left-hand-side of map
    lonroll = np.where(lon>lonmin+360)[0] # tuple of ids
    if lonroll: # non-empty
        roll = lon.size-min(lonroll) # e.g. if 10 lons, lonmax id is 9, we want to roll once
        lon = np.roll(lon, roll) # need to roll foreward
        data = np.roll(data, roll, axis=0) # roll again
        lon[:roll] -= 360 # retains monotonicity
    # print("After roll 2:", lon)
    # 4) Set NaN where data not in range lonmin, lonmax
    # This needs to be done for some regional smaller projections or otherwise
    # might get weird side-effects due to having valid data way outside of the
    # map boundaries -- e.g. strange polygons inside an NaN region
    data = data.copy()
    if lon.size-1==data.shape[0]: # test western/eastern grid cell edges
        # remove data where east boundary is east of min longitude or west
        # boundary is west of max longitude
        data[(lon[1:]<lonmin) | (lon[:-1]>lonmax),:] = np.nan
    elif lon.size==data.shape[0]: # test the centers
        # this just tests centers and pads by one for safety
        # remember that a SLICE with no valid range just returns empty array
        dwhere = np.where((lon<lonmin) | (lon>lonmax))[0]
        data[dwhere[1:-1],:] = np.nan
    # 5) Fix holes over poles by interpolating there (equivalent to
    # simple mean of highest/lowest latitude points)
    if lon.size==data.shape[0]: # have centers, not grid cell edges
        dataS = np.repeat(data[:,0].mean(), data.shape[0])[:,None]
        dataN = np.repeat(data[:,-1].mean(), data.shape[0])[:,None]
        lat = np.concatenate(([-90], lat, [90]))
        data = np.concatenate((dataS, data, dataN), axis=1)
    # 6) Fix seams at map boundary; 3 scenarios here:
    # a) have edges (e.g. for pcolor), and they fit perfectly against basemap seams
    # this does not augment size
    if lon[0]==lonmin and lon.size-1==data.shape[0]: # borders fit perfectly
        pass # do nothing
    # b) have edges (e.g. for pcolor), and the projection edge is in-between grid cell boundaries
    # this augments size by 1
    elif lon.size-1==data.shape[0]: # no interpolation necessary; just make a new grid cell
        lon = np.append(lonmin, lon) # append way easier than concatenate
        lon[-1] = lonmin+360 # we've added a new tiny cell to the end
        data = np.concatenate((data[-1:,:], data), axis=0) # don't use pad, because messes up masked arrays
        # pad_tuple = ((1,0),) + ((0,0),)*(data.ndim-1)
        # data = np.pad(data, pad_tuple, 'wrap') # pad before
    # c) have centers (e.g. for contourf), and we need to interpolate to the
    # left/right edges of the map boundary
    # this augments size by 2
    elif lon.size==data.shape[0]: # linearly interpolate to the edges
        x = np.array([lon[-1], lon[0]+360]) # x
        y = np.concatenate((data[-1:,:], data[:1,:]), axis=0)
        xq = lonmin+360
        yq = (y[:1,:]*(x[1]-xq) + y[1:,:]*(xq-x[0]))/(x[1]-x[0]) # simple linear interp formula
        data = np.concatenate((yq, data, yq), axis=0)
        lon = np.append(np.append(lonmin, lon), lonmin+360)
    else:
        raise ValueError('Longitude length does not match data length along dimension 0.')
    return lon, lat, data
def _m_coords(m, lon, lat):
    # Convert to 2d x/y by calling basemap object
    lat[lat>90], lat[lat<-90] = 90, -90 # otherwise, weird stuff happens
    X, Y = m(*np.meshgrid(lon, lat))
    return X, Y
def _pcolor_fix(p, linewidth=0.2):
    # Fill in white lines
    p.set_linewidth(linewidth) # seems to do the trick, without dots in corner being visible
    p.set_edgecolor('face')
    return p
def _contour_fix(c, linewidth=0.2):
    # Fill in white lines
    for _ in c.collections:
        _.set_edgecolor('face')
        _.set_linewidth(linewidth)
    return c
def _contourf_levels(kwargs):
    # Processes keyword-arguments to allow levels for pcolor/pcolormesh
    # Using "lut" argument prevents colors in the extention-triangles/rectangles
    # from being different from the deepest colors with in the range of min(levels) max(levels)
    if 'levels' in kwargs:
        if 'cmap' not in kwargs: kwargs['cmap'] = 'viridis'
        kwargs['cmap'] = plt.cm.get_cmap(kwargs['cmap'], lut=len(kwargs['levels'])-1)
        if 'norm' not in kwargs:
            kwargs['norm'] = Norm(kwargs['levels']) # what is wanted 99% of time
            # kwargs['norm'] = Norm(levels, ncolors=kwargs['cmap'].N, clip=True)
    return kwargs
def _pcolor_levels(kwargs):
    # Processes keyword-arguments to allow levels for pcolor/pcolormesh
    if 'levels' in kwargs:
        if 'cmap' not in kwargs: kwargs['cmap'] = 'viridis'
        levels = kwargs.pop('levels') # pcolor can't actually have 'levels'
        kwargs['cmap'] = plt.cm.get_cmap(kwargs['cmap'], lut=len(levels)-1)
        if 'norm' not in kwargs:
            kwargs['norm'] = Norm(levels)
            # kwargs['norm'] = Norm(levels, ncolors=kwargs['cmap'].N, clip=True)
    return kwargs
def _pcolor_check(x, y, Z):
    # Checks that sizes match up, checks whether graticule was input
    x, y, Z = np.array(x), np.array(y), np.array(Z)
    xlen, ylen = x.shape[0], y.shape[-1]
    if Z.shape[0]==xlen and Z.shape[1]==ylen:
        x, y = _graticule_fix(x, y) # guess edges, given centers
    elif Z.shape[0]!=xlen-1 or Z.shape[1]!=ylen-1:
        raise ValueError(f'X ({"x".join(str(i) for i in x.shape)}) '
                f'and Y ({"x".join(str(i) for i in y.shape)}) must correspond to '
                f'nrows ({Z.shape[0]}) and ncolumns ({Z.shape[1]}) of Z, or its borders.')
    return x, y
def _contour_check(x, y, Z):
    # Checks whether sizes match up, checks whether graticule was input
    x, y, Z = np.array(x), np.array(y), np.array(Z)
    xlen, ylen = x.shape[0], y.shape[-1]
    if Z.shape[0]==xlen-1 and Z.shape[1]==ylen-1:
        x, y = (x[1:]+x[:-1])/2, (y[1:]+y[:-1])/2 # guess centers, given edges
    elif Z.shape[0]!=xlen or Z.shape[1]!=ylen:
        raise ValueError(f'X ({"x".join(str(i) for i in x.shape)}) '
                f'and Y ({"x".join(str(i) for i in y.shape)}) must correspond to '
                f'nrows ({Z.shape[0]}) and ncolumns ({Z.shape[1]}) of Z, or its borders.')
    return x, y

#-------------------------------------------------------------------------------
# Cartesian plot overrides; they each reference common helper functions
# TODO: Fix for pcolormesh/pcolor with set of levels
#-------------------------------------------------------------------------------
# Dummy ones that just detect cartopy projection
def _cartopy_kwargs(self, kwargs):
    if isinstance(self, GeoAxes): # the main GeoAxes class; others like GeoAxesSubplot subclass this
        kwargs['transform'] = Transform # default is PlateCarree (what you want 99% of time)
    return kwargs
def _plot(self, *args, **kwargs):
    lines = self._plot(*args, **_cartopy_kwargs(self, kwargs))
    return lines
def _scatter(self, *args, **kwargs):
    points = self._scatter(*args, **_cartopy_kwargs(self, kwargs))
    return points

# TODO TODO TODO special wrapper functions to write!
def _quiver(self, x, y, Z, **kwargs): # quiver plot
    # * Options confusing, so needs wrapper.
    # * Function can plot unpredictable vectors if your array orientation is wrong, and will
    #   not raise an error.
    pass
def _histogram(self, *args, **kwargs): # histogram
    # * Should be able to apply transparency so histograms can be overlaid
    pass
def _box(self, *args, **kwargs): # box whisker plot\
    # * Should be able to make grouped boxes
    pass
def _violin(self, *args, **kwargs): # violin plot
    # * Ditto with box plot function
    pass
def _bar(self, *args, **kwargs): # bar plot with different color options
    # * Perhaps by default cycle through color cycle, otherwise use input color
    pass

# Wrappers that check dimensionality, give better error messages, and change
# the convention so dim1 is always x-axis, dim2 is always y-axis
# WARNING: m.pcolormesh and ax.pcolormesh call m.pcolor and ax.pcolor
# respectively, so ***cannot*** override pcolor! Get weird errors. Since
# pcolor alone returns PolyCollection objects, call it pcolorpoly.
def _contour(self, x, y, Z, **kwargs):
    x, y = _contour_check(x, y, Z)
    x, y, Z = _cartopy_seams(self, x, y, Z)
    c = self._contour(x, y, Z.T, **_cartopy_kwargs(self, kwargs))
    return c
def _contourf(self, x, y, Z, **kwargs):
    kwargs = _contourf_levels(kwargs)
    x, y = _contour_check(x, y, Z)
    x, y, Z = _cartopy_seams(self, x, y, Z)
    c = self._contourf(x, y, Z.T, **_cartopy_kwargs(self, kwargs))
    return _contour_fix(c)
def _pcolorpoly(self, x, y, Z, **kwargs):
    kwargs = _pcolor_levels(kwargs)
    x, y = _pcolor_check(x, y, Z)
    extend = kwargs.pop('extend',None)
    p = self.pcolor(x, y, Z.T, **_cartopy_kwargs(self, kwargs))
    p.extend = extend # add attribute to be used in colorbar creation
    return _pcolor_fix(p)
def _pcolormesh(self, x, y, Z, **kwargs):
    # Not allowed if cartopy GeoAxes
    if isinstance(self, GeoAxes):
        raise ValueError('Mesh version of pcolor fails for map projections. Use pcolorpoly instead.')
    kwargs = _pcolor_levels(kwargs)
    x, y = _pcolor_check(x, y, Z)
    extend = kwargs.pop('extend',None)
    p = self._pcolormesh(x, y, Z.T, **_cartopy_kwargs(self, kwargs))
    p.extend = extend # add attribute to be used in colorbar creation
    return _pcolor_fix(p)

#-------------------------------------------------------------------------------
# Basemap overrides
# Will convert longitude latitude coordinates to projected coordinates
#-------------------------------------------------------------------------------
# Dummy ones that just set latlon True by default
def _m_plot(self, *args, **kwargs):
    return self._plot(*args, **kwargs, latlon=True)
def _m_scatter(self, *args, **kwargs):
    return self._scatter(*args, **kwargs, latlon=True)

# More complex ones; assume regularly spaced data
def _m_contour(self, lon, lat, Z, **kwargs):
    lon, lat = _contour_check(lon, lat, Z)
    lon, lat, Z = _m_seams(self, lon, lat, Z)
    X, Y = _m_coords(self, lon, lat)
    c = self._contour(X, Y, Z.T, **kwargs)
    return c
def _m_contourf(self, lon, lat, Z, **kwargs):
    kwargs = _contourf_levels(kwargs)
    lon, lat = _contour_check(lon, lat, Z)
    lon, lat, Z = _m_seams(self, lon, lat, Z)
    X, Y = _m_coords(self, lon, lat)
    c = self._contourf(X, Y, Z.T, **kwargs)
    return _contour_fix(c)
def _m_pcolorpoly(self, lon, lat, Z, **kwargs):
    kwargs = _pcolor_levels(kwargs)
    lon, lat = _pcolor_check(lon, lat, Z)
    lon, lat, Z = _m_seams(self, lon, lat, Z)
    X, Y = _m_coords(self, lon, lat)
    extend = kwargs.pop('extend',None)
    p = self.pcolor(X, Y, Z.T, **kwargs)
    p.extend = extend # add attribute to be used in colorbar creation
    return _pcolor_fix(p)
def _m_pcolormesh(self, lon, lat, Z, **kwargs):
    # Dummy function, this is not allowed
    raise ValueError('Mesh version of pcolor fails for map projections. Use pcolorpoly instead.')

#------------------------------------------------------------------------------
# Formatting functions, assigned using MethodType onto various Axes
# instances and GridSpec instances
#------------------------------------------------------------------------------
def _text(self, x, y, text, transform='axes', fancy=False, black=True,
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
    t = self._text(x, y, text, transform=transform, **kwarg)
    if fancy:
        fcolor, bcolor = 'wk'[black], 'kw'[black]
        t.update({'color':fcolor, 'zorder':1e10, # have to update after-the-fact for path effects
            'path_effects': [mpatheffects.Stroke(linewidth=linewidth, foreground=bcolor), mpatheffects.Normal()]})
        # t.update({'size':11, 'zorder':1e10,
        #     'path_effects':[mpatheffects.PathPatchEffect(edgecolor=bcolor,linewidth=.6,facecolor=fcolor)]})
    return t

def _format_colorbar(self, mappable, cgrid=False, clocator=None, cminorlocator=None,
        cformatter=None, clabel=None, errfix=True, extend='neither', triangles=.15, # in inches
        length=0.8, values=None, **kwargs): #, settings=None):
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
      using a Norm (for example) to determine colors between the levels
      (see: https://stackoverflow.com/q/42723538/4970632).
    * Workaround is to make sure locators are in vmin/vmax range EXCLUSIVELY;
      cannot match/exceed values.
    * The 'extend' kwarg is used for the case when you are manufacturing colorbar
      from list of colors or lines. Most of the time want 'neither'.
    TODO:
    Issue appears where the normalization vmin/vmax are outside of explicitly declared "levels"
    minima and maxima but that is probaby appropriate. If your levels are all within vmin/vmax,
    you will get discrete jumps outside of range and the triangles at ends of colorbars will be weird.
    """
    # TODO: Still experimental, consider fixing
    # This is a different idea altogether for colorbar
    if isinstance(self, mgridspec.SubplotSpec): # self is a SubplotSpec object
        if self.bottompanel:
            C = mgridspec.GridSpecFromSubplotSpec(
                    nrows        = 1,
                    ncols        = 3,
                    wspace       = 0,
                    subplot_spec = self,
                    width_ratios = ((1-length)/2, length, (1-length)/2)
                    )
            cax = self.figure.add_subplot(C[0,1])
            orientation = 'horizontal'
        elif self.rightpanel:
            C = mgridspec.GridSpecFromSubplotSpec(
                    nrows         = 3,
                    ncols         = 1,
                    hspace        = 0,
                    subplot_spec  = self,
                    height_ratios = ((1-length)/2, length, (1-length)/2)
                    )
            cax = self.figure.add_subplot(C[1,0])
            orientation = 'vertical'
        else:
            raise ValueError("Object does not appear to have been created by subplots method.")
    else:
        cax = self # the axes for drawing
        if self.bottomcolorbar:
            orientation = 'horizontal'
        elif self.rightcolorbar:
            orientation = 'vertical'
        else:
            raise ValueError('Axes passed to colorbar-formatting function is not colorbar-dedicated.')
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
    csettings = {'spacing':'uniform', 'cax':cax, 'use_gridspec':True, # use space afforded by entire axes
        'extend':extend, 'orientation':orientation, 'drawedges':cgrid} # this is default case unless mappable has special props
    # Update with user-kwargs
    csettings.update(**kwargs)
    if hasattr(mappable, 'extend') and mappable.extend is not None:
        csettings.update({'extend':mappable.extend})
    # Option to generate colorbar/colormap from line handles
    # * Note the colors are perfect if we don't extend them by dummy color on either side,
    #   but for some reason labels for edge colors appear offset from everything
    # * Too tired to figure out why so just use this workaround
    if fromcolors: # we passed the colors directly
        colors = mappable
        if values is None:
            raise ValueError("Must pass \"values\", corresponding to list of colors.")
    if fromlines: # the lines
        if values is None:
            raise ValueError("Must pass \"values\", corresponding to list of handles.")
        if len(mappable)!=len(values):
            raise ValueError("Number of \"values\" should equal number of handles.")
        colors = [h.get_color() for h in mappable]
    if fromlines or fromcolors:
        values = np.array(values) # needed for below
        # colors = ['#ffffff'] + colors + ['#ffffff']
        colors = colors[:1] + colors + colors[-1:]
        colormap = mcolors.LinearSegmentedColormap.from_list('tmp', colors, N=len(colors)) # note that
            # the 'N' is critical; default 'N' is otherwise 256, and can get weird artifacts due to
            # unintentionally sampling some 'new' colormap colors; very bad!
        levels = np.concatenate((values[0]-np.diff(values[:2])/2, # get "edge" values between centers desired
            (values[1:]+values[:-1])/2, values[-1]+np.diff(values[-2:])/2))
        values = np.concatenate((values[:1]-np.diff(levels[:2]), values,
            values[-1:]+np.diff(levels[-2:]))) # this works empirically; otherwise get
            # weird situation with edge-labels appearing on either side
        mappable = plt.contourf([[0,0],[0,0]], levels=levels, cmap=colormap,
                extend=extend, norm=Norm(values)) # workaround
        if clocator is None: # in this case, easy to assume user wants to label each value
            clocator = values[1:-1]
            # clocator = values # restore this if values weren't padded
    # Determine tick locators for major/minor ticks
    # Can pass clocator/cminorlocator as the *jump values* between the mappables
    # vmin/vmax if desired
    normfix = False # whether we need to modify the norm object
    locators = [] # put them here
    for i,locator in enumerate((clocator,cminorlocator)):
        # Create a preliminary object first
        if isinstance(locator, mticker.Locator):
            pass # nothing to do
        elif locator is None and i==1: # means we never wanted minor ticks
            locators.append(mticker.NullLocator())
            continue
        elif locator is None:
            locator = mticker.AutoLocator() # determine automatically
        else: # set based on user input
            try: iter(locator)
            except TypeError:
                locator = AutoLocate(mappable.norm.vmin, mappable.norm.vmax, locator)
            locator = mticker.FixedLocator(np.array(locator)) # basic
        # Then modify to work around that mysterious error, and to prevent annoyance
        # where minor ticks extend beyond the "extend" triangles
        # * Need to use tick_values instead of accessing "locs" attribute because
        #   many locators don't have these attributes; require norm.vmin/vmax as input
        # * Previously the below block was unused; consider changing?
        values = np.array(locator.tick_values(mappable.norm.vmin, mappable.norm.vmax)) # get the current values
        values_min = np.where(values>=mappable.norm.vmin)[0]
        values_max = np.where(values<=mappable.norm.vmax)[0]
        if len(values_min)==0 or len(values_max)==0:
            raise ValueError(f"Ticks are not within the colorbar range {mappable.norm.vmin:.3g} to {mappable.norm.vmax:.3g}.")
        values_min, values_max = values_min[0], values_max[-1]
        values = values[values_min:values_max+1]
        if values[0]==mappable.norm.vmin:
            normfix = True
        locators.append(mticker.FixedLocator(values)) # final locator object
    # Figure out the tickformatter
    if isinstance(cformatter, mticker.Formatter):
        pass # just use that
    elif cformatter is None:
        cformatter = Formatter() # default formatter is my special one
    else:
        if isinstance(cformatter, str):
            cformatter = mticker.FormatStrFormatter(cformatter) # %-style, numbers
        else:
            cformatter = mticker.FixedFormatter(cformatter) # manually set list of strings
    # Fix the norm object
    # Check out the INSANELY WEIRD ERROR that occurs when you comment out
    # this block!
    # * The error is triggered when a MAJOR tick sits exactly on vmin, but
    #   the actual error is due to processing of MINOR ticks, even if the 
    #   minor locator was set to NullLocator; very weird
    # * Happens when we call get_ticklabels(which='both') below. Can be prevented
    #   by just calling which='major'. Minor ticklabels are never drawn anyway.
    # * We can eliminate the normfix below, but that actually causes an annoying
    #   warning to be printed (related to same issue I guess). So we keep this.
    #   The culprit for all of this seems to be the colorbar API line:
    #        z = np.take(y, i0) + (xn - np.take(b, i0)) * dy / db
    # * Also strange that minorticks extending BELOW the minimum
    #   don't raise the error. It is only when they are exaclty on the minimum.
    # * Note that when changing the levels attribute, need to make sure the
    #   levels datatype is float; otherwise division will be truncated and bottom
    #   level will still lie on same location, so error will occur
    if normfix:
        # else: # need to just use the vmax/vmin
        if hasattr(mappable.norm,'levels'): # can use the levels difference
            mappable.norm.levels = np.array(mappable.norm.levels).astype(np.float)
            mappable.norm.levels[0] -= np.diff(mappable.norm.levels[:2])[0]/1000
        mappable.norm.vmin -= (mappable.norm.vmax-mappable.norm.vmin)/1000
        # print(mappable.norm.levels, mappable.norm.vmin, mappable.norm.vmax)
    # Figure out axis stuff
    if orientation=='horizontal':
        axis = cax.xaxis
        interval = 'intervalx'
    else:
        axis = cax.yaxis
        interval = 'intervaly'
    # Draw the colorbar
    # TODO: For WHATEVER REASON the only way to avoid bugs seems to be to PASS
    # THE MAJOR FORMATTER/LOCATOR TO COLORBAR COMMMAND and DIRETLY EDIT the 
    # minor locators/formatters; update_ticks after the fact ignores the major formatter
    # cb = self.figure.colorbar(mappable, **csettings)
    triangles = triangles/(cax.figure.width*np.diff(getattr(cax.get_position(),interval))[0]-2*triangles)
    csettings.update({'extendfrac':triangles}) # width of colorbar axes and stuff
    cb = self.figure.colorbar(mappable, ticks=locators[0], format=cformatter, **csettings)
    # The ticks/ticklabels basic properties
    for t in axis.get_ticklabels(which='major'):
        t.update(globals('ticklabels'))
    axis.set_tick_params(which='major', **globals('ctick'))
    axis.set_tick_params(which='minor', **globals('ctickminor'))
        # properties are obscure; would have to use hidden methods to do this manually; 
        # using update method (inhereted from Artist) doesn't work
    # The locators and formatters
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
    axis.label.update(globals('label'))
    if clabel is not None:
        axis.label.update({'text':clabel})
    # Fix pesky white lines between levels + misalignment with border due to rasterized blocks
    cb.solids.set_rasterized(False)
    cb.solids.set_edgecolor('face')
    # Make edges/dividers consistent with axis edges
    if cb.dividers is not None:
        cb.dividers.update(globals('cgrid'))
    cb.outline.update(globals('outline'))
    # Update the tickers colorbar
    # cb.formatter = cformatter # fucks shit up -- why?
    # cb.update_ticks() # doesn't actually update the formatter
    # axis.set_major_formatter(cformatter) # fucks shit up -- why?
    return cb

def _format_legend(self, handles=None, align=None, handlefix=False, **kwargs): #, settings=None): # can be updated
    """
    Function for formatting legend-axes (invisible axes with centered legends on them).
    Should update my legend function to CLIP the legend box when it goes outside axes area, so
    the legend-width and bottom/right widths can be chosen propertly/separately.
    """
    # First get legend settings (usually just one per plot so don't need to declare
    # this dynamically/globally)
    lsettings = globals('legend')
    # TODO: Still experimental, consider fixing
    # This is a different idea altogether for bottom and right PANELS; in this case
    # the lformat function is called on a subplotspec object and axes drawn after the fact
    if isinstance(self, mgridspec.SubplotSpec): # self is a SubplotSpec object
        self = self.figure.add_subplot(self) # easy peasy
        for s in self.spines.values():
            s.set_visible(False)
        self.xaxis.set_visible(False)
        self.yaxis.set_visible(False)
        self.patch.set_alpha(0)
        self.bottomlegend = True # for reading below
    if self.bottomlegend: # this axes is assigned to hold the bottomlegend?
        if not handles: # must specify
            raise ValueError('Must input list of handles.')
        lsettings.update( # need to input handles manually, because for separate axes
                bbox_transform=self.transAxes, # in case user passes bbox_to_anchor
                borderaxespad=0, loc='upper center', # this aligns top of legend box with top of axes
                frameon=False) # turn off frame by default
    else: # not much to do if this is just an axes
        if not handles: # test if any user input
            handles, _ = self.get_legend_handles_labels()
        if not handles: # just test if Falsey
            raise ValueError('Must input list of handles or pass "label" attributes to your plot calls.')
    # Interpret kwargs, and apply them
    if 'ncols' in kwargs:
        kwargs['ncol'] = kwargs.pop('ncols') # pyplot subplot uses 'ncols', but legend uses 'ncol'... annoying!
    if 'frame' in kwargs: # again, confusing choice
        kwargs['frameon'] = kwargs.pop('frame')
    lsettings.update(**kwargs)
    # Setup legend text properties
    tsettings = globals('ticklabels')
    # Setup legend handle properties
    hsettings = {}
    candidates = ['linewidth','color'] # candidates for modifying legend objects
    for candidate in candidates:
        if candidate in kwargs:
            hsettings[candidate] = lsettings.pop(candidate)
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
            handles = [handle for iterable in handles for handle in iterable]

    # Now draw legend, with two options
    # 1) Normal legend, just draw everything like normal and columns
    # will be aligned; we re-order handles to be row-major, is only difference
    if align:
        # Prepare settings
        if list_of_lists:
            lsettings['ncol'] = len(handles[0]) # choose this for column length
        elif 'ncol' not in lsettings:
            # raise ValueError("Need to specify number of columns with ncol.")
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
        if hasattr(self, '_legend'): # this Method was renamed to "legend" on an Axes; old version is in _legend
            leg = self._legend(handles=handles, **lsettings) # includes number columns
        else: # this Method did not override the original "legend" method
            leg = self.legend(handles=handles, **lsettings)
        # Format handles, text, and legend patch
        leg.legendPatch.update(globals('outline')) # or get_frame()
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
            raise ValueError("Need to know labelspacing before drawing unaligned-column legends. Add it to settings.")
        if 'size' not in tsettings:
            raise ValueError("Need to know fontsize before drawing unaligned-column legends. Add it to settings.")
        for override in ['loc','ncol','bbox_to_anchor','borderpad','borderaxespad','frameon','framealpha']:
            if override in lsettings:
                lsettings.pop(override)
        interval = 1/len(handles) # split up axes
        interval = ((tsettings['size'] + lsettings['labelspacing']*tsettings['size'])/72) / \
                (self.figure.get_figheight() * np.diff(self._position.intervaly))
            # space we want sub-legend to occupy, as fraction of height
            # don't normally save "height" and "width" of axes so keep here
        for h,hs in enumerate(handles):
            bbox = mtransforms.Bbox([[0,1-(h+1)*interval],[1,1-h*interval]])
            if hasattr(self, '_legend'):
                leg = self._legend(handles=hs, ncol=len(hs), bbox_to_anchor=bbox,
                    frameon=False, borderpad=0, loc='center', **lsettings) # _format_legend is overriding original legend Method
            else:
                leg = self.legend(handles=hs, ncol=len(hs), bbox_to_anchor=bbox,
                    frameon=False, borderpad=0, loc='center', **lsettings) # not overriding original Method
            leg.legendPatch.update(globals('outline')) # or get_frame()
            for obj in leg.legendHandles:
                obj.update(hsettings)
            for t in leg.texts:
                t.update(tsettings) # or get_texts()
            legends.append(leg)
        for l in legends[:-1]:
            self.add_artist(l) # because matplotlib deletes previous ones...
    return legends

def _format_axes(self,
    hatch=None, color=None, # control figure/axes background; hatch just applies to axes
    oceans=False, coastlines=True, continents=False, # coastlines and continents
    latlabels=[0,0,0,0], lonlabels=[0,0,0,0], # sides for labels (left, right, bottom, top)
    latlocator=None, latminorlocator=None, lonlocator=None, lonminorlocator=None,
    xgrid=None, ygrid=None, # gridline toggle
    xdates=False, ydates=False, # whether to format axis labels as long datetime strings; the formatter should be a date %-style string
    xtickminor=None, ytickminor=None, xgridminor=None, ygridminor=None, # minor ticks/grids; if ticks off, grid will be off
    xspineloc=None, yspineloc=None, # deals with spine options
    xtickloc=None, ytickloc=None, # tick location
    xtickdir=None, ytickdir=None, # change ytick/xtick location; can be in, out, or inout (left-right-up-down depends on spine to which this applied)
    xtickrange=None, ytickrange=None, # limit regions where we assign ticklabels to major-ticks
    xlim=None, ylim=None, xscale=None, yscale=None, xscale_kwargs={}, yscale_kwargs={},
    xreverse=False, yreverse=False, # special properties
    xlabel=None, ylabel=None, # axis labels
    suptitle=None, suptitlepos=None, title=None, titlepos=None, titlepad=0.1, titledict={},
    abc=False, abcpos=None, abcformat='', abcpad=0.1, abcdict={},
    # TODO: add options for allowing UNLABELLED major ticklines -- maybe pass a special Formatter?
    xlocator=None, xminorlocator=None, ylocator=None, yminorlocator=None, # locators, or derivatives that are passed to locators
    xformatter=None, yformatter=None): # formatter
    # legend=False, handles=None, # legend options; if declared and have no handles with 'label', does nothing
    # settings=None): # remaining are passed to Settings()
    """
    Function for formatting axes of all kinds; some arguments are only relevant to special axes, like 
    colorbar or basemap axes. By default, simply applies the dictionary values from settings() above,
    but can supply many kwargs to further modify things.
    TODO:
    * Add options for datetime handling; note possible date axes handles are TimeStamp (pandas),
      np.datetime64, DateTimeIndex; can fix with fig.autofmt_xdate() or manually set options; uses
      ax.is_last_row() or ax.is_first_column(), so should look that up. Problem is there
      is no autofmt_ydate(), so really should implement my own version of this.
    """
    #---------------------------------------------------------------------------
    # Preliminary stuff, stuff that also applies to basemap
    #---------------------------------------------------------------------------
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
        fig.suptitle.update({'ha':'center', 'va':'baseline', **globals('suptitle')})
        if suptitlepos=='title': # not elevated
            ypos = self.title._transform.transform(self.title.get_position())[1]/fig.dpi/fig.height
            fig.suptitle.update({'position':(xpos,ypos)})
        elif isinstance(suptitlepos,str):
            raise ValueError(f"Unknown suptitle position: {suptitlepos}.")
        elif suptitlepos is not None:
            fig.suptitle.update({'position':suptitlepos})
    # Create axes title
    # Input needs to be emptys string
    self.title.update({**globals('title'), 'text':title or ''})
    self.title.update(titledict)
    if titlepos=='left':
        self.title.update({'position':(0,1), 'ha':'left'})
    elif titlepos=='right':
        self.title.update({'position':(1,1), 'ha':'right'})
    elif titlepos=='inside':
        self.title.update({'position':(0.5,1-titlepad/self.height),
            'transform':self.transAxes, 'va':'top'})
    elif isinstance(titlepos,str):
        raise ValueError(f"Unknown title position: {titlepos}.")
    elif titlepos is not None: 
        self.title.update({'position':titlepos, 'ha':'center', 'va':'center'})
    # Create axes numbering
    if self.number is not None and abc:
        abcedges = abcformat.split('a')
        self.abc = self.text(0, 1, abcedges[0] + ascii_lowercase[self.number-1] + abcedges[-1],
                transform=self.title._transform, **abcdict) # optionally include paren
        self.abc.update({'ha':'left', 'va':'baseline', **globals('abc')})
        if abcpos=='inside':
            self.abc.update({'position':(abcpad/self.width, 1-abcpad/self.height),
                'transform':self.transAxes, 'ha':'left', 'va':'top'})
        elif isinstance(abcpos,str):
            raise ValueError(f"Unknown abc position: {abcpos}.")
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

    #--------------------------------------------------------------------------
    # Basemap axes setup
    #--------------------------------------------------------------------------
    # For now force background to never be transparent, logic being that the background
    # itself contains information (e.g. "this is ocean" or "this is land")
    if getattr(self, 'm', None) is not None: # the attribute exists and isn't None
        # Coastlines, parallels, meridians
        m = self.m
        if coastlines:
            p = m.drawcoastlines(**globals('coastlines'))
        if continents:
            p = m.fillcontinents(**globals('continents'))
        # Longitude/latitude lines
        # Make sure to turn off clipping by invisible axes boundary; otherwise
        # get these weird flat edges where map boundaries, parallel/meridian markers come up to the axes bbox
        if ylocator is not None:
            if not hasattr(latlocator,'__iter__'):
                latlocator = AutoLocate(m.latmin+latlocator, m.latmax-latlocator, latlocator)
            p = m.drawparallels(latlocator, labels=latlabels)
            for pi in p.values(): # returns dict, where each one is tuple
                for _ in [i for j in pi for i in j]: # magic
                    if isinstance(_, mtext.Text):
                        _.update(globals('ticklabels'))
                    else:
                        _.set_clip_on(True) # no gridlines past boundary
                        _.update(globals('lonlatlines'))
                        # _.set_linestyle(linestyle)
                # tried passing clip_on to the above, but it does nothing; must set
                # for lines created after the fact
        if lonlocator is not None:
            latlabels[2:] = lablabels[2:][::-1] # default is left/right/top/bottom which is dumb
            lonlabels[2:] = lonlabels[2:][::-1] # change to left/right/bottom/top
            if not hasattr(lonlocator,'__iter__'):
                lonlocator = AutoLocate(m.lonmin+lonlocator, m.lonmax-lonlocator, lonlocator)
            p = m.drawmeridians(lonlocator, labels=lonlabels)
            for pi in p.values():
                for _ in [i for j in pi for i in j]: # magic
                    if isinstance(_, mtext.Text):
                        _.update(globals('ticklabels'))
                    else:
                        _.set_clip_on(True) # no gridlines past boundary
                        _.update(globals('lonlatlines'))
                        # _.set_linestyle(linestyle)
        # Map boundary
        # * First have to MANUALLY REPLACE the old boundary by just deleting
        #   the original one; note this requires drawmapboundary() was called
        #   when the basemap was first instantiated; see notes in subplots() command
        # * If boundary is drawn successfully should be able to call
        #   m._mapboundarydrawn.set_visible(False) and edges/fill color disappear
        # * For now will enforce that map plots ALWAYS have background whereas
        #   axes plots can have transparent background
        if m._mapboundarydrawn:
            m._mapboundarydrawn.remove()
        if m.projection in m._pseudocyl:
            self.patch.set_alpha(0) # make patch invisible
            p = m.drawmapboundary(fill_color=color, **globals('spine')) # set fill_color to 'none' to make transparent
            p.set_rasterized(False) # not sure about this; might be rasterized
            p.set_clip_on(False) # so edges of LINE denoting boundary aren't cut off
        else: # use the settings to apply to Axes patch; the basemap API fails here
            self.patch.set_facecolor(color)
            self.patch.set_edgecolor('none')
            for spine in self.spines.values():
                spine.update(globals('spine'))
        return # skip everything else

    #--------------------------------------------------------------------------
    # Cartopy axes setup
    #--------------------------------------------------------------------------
    if isinstance(self, GeoAxes): # the main GeoAxes class; others like GeoAxesSubplot subclass this
        # Boundary
        self.outline_patch.update(globals('outline'))
        # Make background patch invisible, so can color axes patch instead
        # See: https://stackoverflow.com/a/32208700/4970632
        # self.background_patch.set_fill(False)
        # This doesn't seem to work because axes patch is always invisible and has
        # been overridden by background patch, so will just show underlying figure color
        if color is not None:
            self.background_patch.set_facecolor(color)
        # Draw gridlines
        # WARNING: For some reason very weird side effects happen if you try
        # to call gridlines() twice on same axes. Can't do it. Which is why
        # we do this nonsense with the formatter below, instead of drawing 'major'
        # grid lines and 'minor' grid lines.
        lonvec = lambda v: [] if v is None else [*v] if hasattr(v,'__iter__') else [*AutoLocate(-180,180,v)]
        latvec = lambda v: [] if v is None else [*v] if hasattr(v,'__iter__') else [*AutoLocate(-90,90,v)]
        lonminorlocator, latminorlocator = lonvec(lonminorlocator), latvec(latminorlocator)
        lonlocator, latlocator = lonvec(lonlocator), latvec(latlocator)
        lonlines = lonminorlocator or lonlocator # where we draw gridlines
        latlines = latminorlocator or latlocator
        # First take care of gridlines
        draw_labels = (isinstance(self.projection,ccrs.Mercator) or isinstance(self.projection,ccrs.PlateCarree))
        if latlines[0]==-90:
            latlines[0] += 0.001
        if lonlines[0]==-90:
            lonlines[0] -= 0.001
        gl = self.gridlines(**globals('lonlatlines'), draw_labels=draw_labels)
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
            # self.add_feature(cfeature.COASTLINE, **globals('coastlines'))
            self.add_feature(cfeature.NaturalEarthFeature('physical', 'coastline', '50m'), **globals('coastlines'))
        if continents:
            # self.add_feature(cfeature.LAND, **globals('continents'))
            self.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m'), **globals('continents'))
        if oceans:
            # self.add_feature(cfeature.OCEAN, **globals('oceans'))
            self.add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', '50m'), **globals('oceans'))
        # return asdfs
        return # skip everything else

    #--------------------------------------------------------------------------
    # Process normal axes, various x/y settings individually
    #--------------------------------------------------------------------------
    # Axes scaling, limits, and reversal options (alternatively, supply
    # your own xlim/ylim that go from high to low)
    if xscale is not None: self.set_xscale(xscale, **xscale_kwargs)
    if yscale is not None: self.set_yscale(yscale, **yscale_kwargs)
    if xlim is None: xlim = self.get_xlim()
    if ylim is None: ylim = self.get_ylim()
    if xreverse: xlim = xlim[::-1]
    if yreverse: ylim = ylim[::-1]
    self.set_xlim(xlim)
    self.set_ylim(ylim)
    for axis, label, dates, sides, tickloc, spineloc, gridminor, tickminor, tickminorlocator, \
            grid, ticklocator, tickformatter, tickrange, tickdir in \
        zip((self.xaxis, self.yaxis), (xlabel, ylabel), (xdates, ydates), \
            (('bottom','top'),('left','right')), (xtickloc,ytickloc), (xspineloc, yspineloc), # other stuff
            (xgridminor, ygridminor), (xtickminor, ytickminor), (xminorlocator, yminorlocator), # minor ticks
            (xgrid, ygrid), (xlocator, ylocator), (xformatter, yformatter), # major ticks
            (xtickrange, ytickrange), # range in which we label major ticks
            (xtickdir, ytickdir)): # tick direction
        # Axis spine visibility and location
        for spine, side in zip((self.spines[s] for s in sides), sides):
            # Line properties
            spine.update(globals('spine'))
            if axis==self.xaxis and self.xspine_override is not None:
                spineloc = self.xspine_override
            if axis==self.yaxis and self.yspine_override is not None:
                spineloc = self.yspine_override
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
                if side==sides[0]: # move the left/bottom spine onto the specified location, with set_position
                    spine.set_visible(True)
                    spine.set_position(spineloc)
                else:
                    spine.set_visible(False)

        # Tick properties
        # * Some weird issue seems to cause set_tick_params to reset/forget that the grid
        #   is turned on if you access tick.gridOn directly, instead of passing through tick_params.
        #   Since gridOn is undocumented feature, don't use it.
        #   So calling _format_axes() a second time will remove the lines
        # * Can specify whether the left/right/bottom/top spines get ticks; sides will be 
        #   group of left/right or top/bottom
        # * Includes option to draw spines but not draw ticks on that spine, e.g.
        #   on the left/right edges
        ticklocs = sides if tickloc=='both' else () if tickloc=='neither' else tickloc
        if ticklocs is None: # only turn off ticks if the spine is invisible; don't explicitly turn on
            ticks_sides = {side: False for side in sides if not self.spines[side].get_visible()}
        else:
            ticks_sides = {side: self.spines[side].get_visible() and side in ticklocs for side in sides}
        ticks_major = globals('tick') if tickdir is None else dict(globals('tick'), direction=tickdir)
        ticks_minor = globals('tickminor') if tickdir is None else dict(globals('tickminor'), direction=tickdir)
        axis.set_tick_params(which='major', **ticks_sides, **ticks_major)
        axis.set_tick_params(which='minor', **ticks_sides, **ticks_minor) # have length

        # Gridline activation and setting (necessary because rcParams has no 'minorgrid'
        # property, must be set in rcExtras)
        if grid is not None: # grid changes must be after tick
            axis.grid(grid, which='major', **globals('grid'))
        if gridminor is not None:
            axis.grid(gridminor, which='minor', **globals('gridminor')) # ignore if no minor ticks
        # if type(grid) is bool: # grid changes must be after tick
        #     axis.grid(grid, which='major')
        # for tick in axis.majorTicks:
        #     tick.gridline.update(globals('grid'))
        # if type(gridminor) is bool:
        #     axis.grid(gridminor, which='minor', alpha=1) # ignore if no minor ticks
        # for tick in axis.minorTicks:
        #     tick.gridline.update(globals('gridminor'))

        # Label properties
        axis.label.update(globals('label'))
        if label is not None:
            axis.label.set_text(label)

        # First, major tick locators (should not affect limits)
        lim = axis.get_view_interval() # to be used, automatically
        if isinstance(ticklocator, mticker.Locator):
            axis.set_major_locator(ticklocator)
        elif ticklocator is not None:
            try: iter(ticklocator)
            except TypeError: ticklocator = AutoLocate(lim[0], lim[1], ticklocator)
            # axis.set_ticks(ticklocator)
            axis.set_major_locator(mticker.FixedLocator(ticklocator))

        # Next, minor tick locators (toggle visibility later)
        if tickminor is not None and not tickminor:
            axis.set_minor_locator(mticker.NullLocator())
        if isinstance(tickminorlocator, mticker.Locator):
            axis.set_minor_locator(tickminorlocator) # pass locator class/subclass (all documented ones subclass the parent class, mticker.Locator)
        elif tickminorlocator is not None:
            try: iter(tickminorlocator)
            except TypeError: tickminorlocator = AutoLocate(lim[0], lim[1], tickminorlocator) # use **spacing/setp** as specified
            axis.set_minor_locator(mticker.FixedLocator(tickminorlocator))

        # Next, major tick formatters (enforce Null, always, for minor ticks), and text styling
        # Includes option for %-formatting of numbers and dates, passing a list of strings
        # for explicitly overwriting the text
        axis.set_minor_formatter(mticker.NullFormatter())
        if tickformatter in ['lat']:
            axis.set_major_formatter(LatFormatter(sine=False))
        elif tickformatter in ['sine']:
            axis.set_major_formatter(LatFormatter(sine=True))
        elif tickformatter is None: # default use my special super cool formatter
            axis.set_major_formatter(Formatter(2, tickrange))
        elif tickformatter in ['none','None','NONE']: # eliminate labels
            axis.set_major_formatter(mticker.NullFormatter())
        elif isinstance(tickformatter, mticker.Formatter): # formatter object
            axis.set_major_formatter(tickformatter)
        elif isinstance(tickformatter, str): # a %-style formatter
            if dates: axis.set_major_formatter(mdates.DateFormatter(tickformatter)) # %-style, dates
            else: axis.set_major_formatter(mticker.FormatStrFormatter(tickformatter)) # %-style, numbers
        else: # list of strings
            axis.set_major_formatter(mticker.FixedFormatter(tickformatter)) # list of strings
            # axis.set_ticklabels(tickformatter) # no issues with FixedFormatter so far
        for t in axis.get_ticklabels():
            t.update(globals('ticklabels'))
    return # we're done

#------------------------------------------------------------------------------#
# Other special methods and attributes for our custom axes
#------------------------------------------------------------------------------#
# Helper function; generates 'circle' path, to be applied in axes coordinates
# as the new GeoAxes boundary.
def _circle(N):
    theta = np.linspace(0, 2*np.pi, N)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    return mpath.Path(verts * radius + center)

# Some functions that need to be declared here
# * TODO WARNING: If drawing more than x lines and don't explicitly set every item below
#   its value will revert to the cycler; this may give unexpected results
# * To print current cycle, use list(next(ax._get_lines.prop_cycler)['color'] for i in range(10))
def _atts_global(self):
    # Set up some properties
    # self.set_prop_cycle(propcycle)
    self.m = None # optional basemap instance
    self.bottomlegend, self.bottomcolorbar, self.rightcolorbar = False, False, False # identifiers
    self.xspine_override, self.yspine_override = None, None
    self.number = None # default
    self.width  = np.diff(self._position.intervalx)*self.figure.width # position is in figure units
    self.height = np.diff(self._position.intervaly)*self.figure.height
    # Set up some methods
    # TODO Note severe change; twiny now means "share the x-axis but
    # now make two y-axes"; this makes way more sense to me than default behavior
    self._legend = self.legend # custom legend drawing on any axes
    self._twinx  = self.twinx
    self._twiny  = self.twiny
    self._text   = self.text
    self.legend  = MethodType(_format_legend, self) # MethodType approach
    self.twinx   = MethodType(_twinx, self)
    self.twiny   = MethodType(_twiny, self)
    self.text    = MethodType(_text, self)
    self.format  = MethodType(_format_axes, self)

def _atts_special(self, maps, package, projection, **kwargs):
    # Instantiate the Basemap object and add contouring methods
    # Note we CANNOT modify the underlying AXES pcolormesh etc. because this
    # this will cause basemap's m.pcolormesh etc. to use my CUSTOM version and
    # cause a suite of weird errors
    # * Must set boundary before-hand, otherwise the set_axes_limits method called
    #   by mcontourf/mpcolormesh/etc draws two mapboundary Patch objects called "limb1" and
    #   "limb2" automatically: one for fill and the other for the edges
    # * Then, since the patch object in _mapboundarydrawn is only the fill-version, calling
    #   drawmapboundary() again will replace only *that one*, but the original visible edges
    #   are still drawn -- so e.g. you can't change the color
    # * If you instead call drawmapboundary right away, _mapboundarydrawn will contain
    #   both the edges and the fill; so calling it again will replace *both*
    if package!='basemap' or not maps:
        # Set up cartopy stuff
        # For some reason trying to do circle in projection coordinates, e.g.
        # line from 0 to 360 degrees longitude, fails; gotta use axes coordiantes
        if isinstance(self, GeoAxes) and any(isinstance(self.projection, projection) \
            for projection in (ccrs.LambertAzimuthalEqualArea, ccrs.AzimuthalEquidistant)):
            # self.projection.threshold = \
            #     kwargs.pop('threshold', self.projection.threshold) # optionally modify threshold
            latmin = kwargs.pop('latmin', 0)
            latmin = kwargs.pop('lat_min', latmin) # also try this guy
            self.set_extent([-180,180,latmin,90], Transform) # use platecarree transform
            circle = _circle(100)
            self.set_boundary(circle, transform=self.transAxes)
        # Employ some simple overrides
        # Will save the *original* versions in _<method name>
        # WARNING: m.pcolormesh and ax.pcolormesh call m.pcolor and ax.pcolor
        # respectively, so ***cannot*** override pcolor! Get weird errors. Since
        # pcolor alone returns PolyCollection objects, call it pcolorpoly.
        self._pcolormesh = self.pcolormesh
        self._contourf   = self.contourf
        self._contour    = self.contour
        self._scatter    = self.scatter
        self._plot       = self.plot
        self.pcolorpoly  = MethodType(_pcolorpoly, self)
        self.pcolormesh  = MethodType(_pcolormesh, self)
        self.contourf    = MethodType(_contourf, self)
        self.contour     = MethodType(_contour, self)
        self.scatter     = MethodType(_scatter, self)
        self.plot        = MethodType(_plot, self)
    else:
        # Always need to plot with the basemap versions
        # These cannot be passed directly as elements of a, for some reason
        # Can't override original name of axes object, because m.pcolor calls m.ax.pcolor
        self.m             = mbasemap.Basemap(ax=self, projection=projection, **kwargs)
        self.m.projection  = projection # save projection here
        self.m._pseudocyl  = ['moll','robin','eck4','kav7','sinu','mbtfpq','vandg','hammer']
        self.m._pcolormesh = self.m.pcolormesh
        self.m._contour    = self.m.contour
        self.m._contourf   = self.m.contourf # save old methods
        self.m._scatter    = self.m.scatter
        self.m._plot       = self.m.plot
        self.m.pcolormesh  = MethodType(_m_pcolormesh, self.m) # function calls the old methods
        self.m.pcolorpoly  = MethodType(_m_pcolorpoly, self.m)
        self.m.contour     = MethodType(_m_contour, self.m)
        self.m.contourf    = MethodType(_m_contourf, self.m)
        self.m.scatter     = MethodType(_m_scatter, self.m)
        self.m.plot        = MethodType(_m_plot, self.m)
        # Initialize map boundary if pseudocylindrical projection
        # Initialization has to happen here (see above)
        if projection in self.m._pseudocyl:
            self.m.drawmapboundary()

def _twinx(self, **kwargs):
    # Create secondary x-axes
    # Format function will read extra properties and *enforce* (ignoring
    # user settings) the spine locations.
    ax = self._twiny(**kwargs)
    _atts_global(ax) # basic setup
    self.xspine_override = 'bottom' # original axis ticks on bottom
    ax.yspine_override = 'neither'
    ax.xspine_override = 'top' # new axis ticks on top
    return ax

def _twiny(self, **kwargs):
    # Create secondary y-axes
    # Same as above
    ax = self._twinx(**kwargs)
    _atts_global(ax) # basic setup
    self.yspine_override = 'left' # original axis ticks on left
    ax.xspine_override = 'neither'
    ax.yspine_override = 'right' # new axis ticks on right
    return ax

# Helper function for creating panelled axes
# Will take in some arguments from parent function so don't have to pass them
# every time
def _panel_factory(fig, subspec, width, height,
        whichpanels=None, hspace=None, wspace=None, hwidth=None, wwidth=None,
        **kwargs):
    # TODO: Standardize the order that panels are delivered (e.g. in a list
    # [None, axh, None, axh]), then *fix axis sharing across multiple panels*!
    # Checks and defaults
    # Defaults format the plot to have tiny spaces between subplots
    # and narrow panels
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
        raise ValueError("Whichpanels argument can contain characters l (left), r (right), "
                "b (bottom), or t (top).")
    # Determine number of rows/columns, position of the
    # main/primary axes, and potential position of corners that
    # we want to ignore for now
    nrows_i, ncols_i = 1, 1
    for s in whichpanels:
        if s in ('b','t'):
            nrows_i += 1
        if s in ('l','r'):
            ncols_i += 1
    bad_pos = []
    main_pos = [0,0]
    if 't' in whichpanels:
        main_pos[0] += 1
    if 'l' in whichpanels:
        main_pos[1] += 1
    corners = {'tl':(0,0),             'tr':(0,main_pos[1]+1),
               'bl':(main_pos[0]+1,0), 'br':(main_pos[0]+1,main_pos[1]+1)}
    for corner,position in corners.items():
        if all(s in whichpanels for s in corner):
            bad_pos.append(position)
    # Fix wspace/hspace in inches, using the Bbox from get_postition
    # on the subspec object to determine physical width of axes to be created
    # * Consider writing some convenience funcs to automate this unit conversion
    bbox_i = subspec.get_position(fig) # valid since axes not drawn yet
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
        raise ValueError(f"Panel wwidth is too large. Must be less than {width_i/(nrows_i-1):.3f}.")
    for i in range(ncols_i):
        if i!=main_pos[1]: # this is a panel entry
            wwidth_ratios[i] = wwidth
    hwidth_ratios = [height_i-hwidth*(nrows_i-1)]*nrows_i
    if hwidth_ratios[0]<0:
        raise ValueError(f"Panel hwidth is too large. Must be less than {height_i/(ncols_i-1):.3f}.")
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
    axlist = []
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
    for r in range(nrows_i)[::-1]: # order is bottom-top
        for c in range(ncols_i): # order is left-right
            if (r,c) not in bad_pos:
                # Add the subplot first
                iax = fig.add_subplot(gs_i[r,c], **kwargs)
                if r==main_pos[0] and c==main_pos[1]:
                    axmain = iax
                else:
                    axlist.append(iax) # add axes
                # Configure shared axes
                if (r,c)==ybase:
                    ax_ybase = iax
                if (r,c)==xbase:
                    ax_xbase = iax
                if (r,c) in yshare:
                    if ax_ybase is None:
                        raise ValueError("Axes with base for x-axis not drawn.")
                    iax._sharey = ax_ybase
                    for t in iax.yaxis.get_ticklabels(): t.set_visible(False)
                    iax.yaxis.label.set_visible(False)
                if (r,c) in xshare:
                    if ax_xbase is None:
                        raise ValueError("Axes with base for x-axis not drawn.")
                    iax._sharex = ax_xbase
                    for t in iax.xaxis.get_ticklabels(): t.set_visible(False)
                    iax.xaxis.label.set_visible(False)
    # Then should add panels as children of the main axes
    # * Can gain access to panel
    # * Unsure how the hell to do this
    return axmain, axlist

#------------------------------------------------------------------------------#
# Custom settings for various journals
# Add to this throughout your career, or as standards change
# PNAS info: http://www.pnas.org/page/authors/submission
# AMS info: https://www.ametsoc.org/ams/index.cfm/publications/authors/journal-and-bams-authors/figure-information-for-authors/
# AGU info: https://publications.agu.org/author-resource-center/figures-faq/
#------------------------------------------------------------------------------#
def _journal_sizes(width, height):
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
    try: iter(value)
    except TypeError:
        width = value
    else:
        width, height = value
    return width, height

#-------------------------------------------------------------------------------
# Primary plotting function; must be used to create figure/axes if user wants
# to use the other features
#-------------------------------------------------------------------------------
def subplots(array=None, nrows=1, ncols=1, emptycols=None, emptyrows=None, silent=False,
        tight=False, # whether to set up tight bbox from gridspec object
        sharex=True, sharey=True, # for sharing x/y axis limits/scales/locators for axes with matching GridSpec extents, and making ticklabels/labels invisible
        spanx=True, spany=True, # custom setting; share axis labels for axes with same xmin/ymin extents
            # what about bottom vs. top xlabels, left vs. right ylabels?
        aspect=None, height=None, width=None, # for controlling aspect ratio; default is control for width
        hspace=0.4, wspace=0.2, hratios=None, wratios=None, # options for gridspec (instead of gridpsec_kw)
            # spacing betwen axes, in inches; generally horizontal-spacing small, but vertical-spacing big to allow for title
        left=0.5, bottom=0.5, right=0.1, top=0.3,
            # edges of plot for GridSpec object, in inches; left needs a lot of space
        bottomcolorbar=False, rightcolorbar=False, bottomlegend=False,
            # for legend/colorbar referencing multiple subplots
        bottompanel=False, bottompanels=False, rightpanel=False, rightpanels=False,
        bwidth=0.25, bspace=0.5, rwidth=0.25, rspace=0.5, # default to no space between panels
            # also need to consider drawing legend *centered* so can put e.g. colorbar next to
            # legend in a bottompanel and have it not look super weird
        lwidth=0.15, cwidth=0.25, cspace=0.5, cshrink=0.9,
            # spacing for colorbar text, colorbar axes width, legend axes with, and padding for interior ABC label
        innerpanels=None, # same as below; list of numbers where we want subplotspecs
        whichpanels=None, ihspace=None, iwspace=None, ihwidth=None, iwwidth=None,
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
    * Use bottomcolorbar/bottomlegend/rightcolorbar convenience features to create empty
      axes, from which you can create objects with ax.bottomcolorbar.format and ax.bottomlegend.format;
      this was the old API before idea for 'panels' was hatched.
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
    Todo:
        * For spanning axes labels, right now only detect **x labels on bottom**
        and **ylabels on top**; generalize for all subplot edges.
        * Figure size should be constrained by the dimensions of the axes, not vice
        versa; might make things easier.
    """
    # Fill in some basic required settings
    if package not in ['basemap','cartopy']:
        raise ValueError("Plotting package must be one of basemap or cartopy.")
    width, height = _journal_sizes(width, height) # if user passed width=<string>, will use that journal size
    if width is None and height is None: # at least one must have value
        width = 5
    if aspect is None: # aspect is width to height ratio
        aspect = 1.5   # try this as a default; square plots (1) look too tall to me
    if wratios is not None:
        aspect = aspect/(wratios[0]/np.mean(wratios)) # e.g. if 2 columns, 5:1 width ratio, change the 'average' aspect ratio
    if hratios is not None:
        aspect = aspect*(hratios[0]/np.mean(hratios))
    # Automatically generate array of first *arg not provided, or use nrows/ncols
    if array is None:
        array = np.arange(1,nrows*ncols+1)[...,None]
        array = array.reshape((nrows, ncols)) # numpy is row-major, remember
    array = np.array(array) # enforce array type
    array[array==None] = 0 # use zero for placeholder; otherwise have issues
    # Include functionality to declare certian rows/columns empty, if user requested it
    # Useful if user wants to include extra space between some columns greater than the
    # spaces between axes in other columns
    if emptycols is not None:
        try: iter(emptycols)
        except TypeError:
            emptycols = [emptycols]
        for col in emptycols:
            array[:,col-1] = 0
    if emptyrows is not None:
        try: iter(emptyrows)
        except TypeError:
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
    # Allow user to pass projection dictionary, cause sometimes kwargs overlap
    # First make sure user didn't mess up
    projection_kwargs.update(projection_dict)
    if projection_kwargs and not maps:
        raise ValueError(f"Unknown keyword args: {', '.join(projection_kwargs.keys())}. "
            "If you want to create maps, you must change the \"maps\" kwarg.")

    # Basemap stuff
    if maps and package=='basemap':
        # Just determine the correct aspect ratio
        # To do this need to instantiate basemap object which has dual utility of
        # telling us before draw-time whether any kwpair in projection_kwargs is bad
        projection_kwargs.update({'fix_aspect':True})
        projection = projection or 'cyl' # cylindrical by default
        mexample = mbasemap.Basemap(projection=projection, **projection_kwargs)
        # Override aspect ratio
        aspect = (mexample.urcrnrx-mexample.llcrnrx)/(mexample.urcrnry-mexample.llcrnry)
        if not silent: print(f"Forcing aspect ratio: {aspect:.3g}")

    # Cartopy stuff; get crs instance, and create dictionary to add to add_subplot calls
    cartopy_kwargs = {}
    if maps and package=='cartopy':
        # Get the projection instance from a string and determine the
        # correct aspect ratio (TODO) also prevent auto-scaling
        projection = projection or 'cyl'
        crs_dict = {
            **{key: ccrs.PlateCarree   for key in ('cyl','rectilinear','pcarree','platecarree')},
            **{key: ccrs.Mollweide     for key in ('moll','mollweide')},
            **{key: ccrs.Stereographic for key in ('stereo','stereographic')},
            'mercator': ccrs.Mercator,
            'robinson': ccrs.Robinson,
            'ortho':    ccrs.Orthographic,
            'aeqd':     ccrs.AzimuthalEquidistant,
            'aeqa':     ccrs.LambertAzimuthalEqualArea,
            'wintri': WinkelTripel,
            'hammer': Hammer,
            'aitoff': Aitoff,
            'kav7':   KavrayskiyVII,
            }
        crs_translate = { # less verbose keywords, actually match proj4 keywords and similar to basemap
            'lat0':  'central_latitude',
            'lat_0': 'central_latitude',
            'lon0':  'central_longitude',
            'lon_0': 'central_longitude',
            }
        postprocess_keys = ('latmin', 'lat_min', 'threshold') # will be processed down the line
        if projection not in crs_dict:
            raise ValueError(f"For cartopy, projection must be one of the following: {', '.join(crs_dict.keys())}.")
        init_kwargs = {crs_translate.get(key,key):value for key,value in projection_kwargs.items() if key not in postprocess_keys}
        projection_kwargs = {key:value for key,value in projection_kwargs.items() if key in postprocess_keys}
        projection = crs_dict[projection](**init_kwargs)
        cartopy_kwargs = {'projection': projection}
        # Override aspect ratio
        aspect = (np.diff(projection.x_limits)/np.diff(projection.y_limits))[0]
        if not silent: print(f"Forcing aspect ratio: {aspect:.3g}")

    # Aspect ratio test
    if maps:
        wtest = [0] if wratios is None else wratios
        htest = [0] if hratios is None else hratios
        if any(w!=wtest[0] for w in wtest) or any(h!=htest[0] for h in htest):
           raise ValueError("Cannot set wratios or hratios when plotting map projections.")

    #--------------------------------------------------------------------------
    # Panel considerations; keep things general for future imporvements
    #--------------------------------------------------------------------------
    # Detect conflicts
    # Or Ignore dat shit.
    # if bottomcolorbar and rightcolorbar:
    #     raise ValueError('You can only specify one global colorbar.')
    if bottomcolorbar and bottomlegend:
        raise ValueError('Can only place a global colorbar **or** global legend on bottom, not both.')

    # Spacings due to panel considerations
    if bottompanel or bottompanels:
        bottom_extra_axes = bwidth
        bottom_extra_hspace = bspace
    elif bottomcolorbar:
        bottom_extra_axes   = cwidth
        bottom_extra_hspace = cspace
    elif bottomlegend:
        bottom_extra_axes   = lwidth
        bottom_extra_hspace = 0 # don't need space
    else:
        bottom_extra_axes, bottom_extra_hspace = 0, 0
    bottom_extra = bottom_extra_axes + bottom_extra_hspace

    # And right spacings
    if rightpanel or rightpanels:
        right_extra_axes   = rwidth
        right_extra_wspace = rspace
    elif rightcolorbar:
        right_extra_axes   = cwidth
        right_extra_wspace = cspace
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
    if bottom_extra_axes>0:
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

    # Create figure
    fig = plt.figure(FigureClass=Figure, figsize=figsize)
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
            height_ratios = height_ratios_outer
            ) # set wspace/hspace to match the top/bottom spaces
    # Initialize some stuff
    fig.left     = ileft
    fig.bottom   = ibottom
    fig.right    = iright
    fig.top      = itop
    fig.lwidth   = lwidth
    fig.cwidth   = lwidth
    fig.cspace   = cspace
    fig.bwidth   = bwidth
    fig.bspace   = bspace
    fig.rwidth   = rwidth
    fig.rspace   = rspace
    fig.height   = height
    fig.width    = width
    fig.gridspec = GS # add to figure, for later reference
    # Initialize some more stuff
    fig.bottomlegend   = None
    fig.bottomcolorbar = None
    fig.rightcolorbar  = None
    fig.bottompanel    = None
    fig.rightpanel     = None
    # Create axes within the 'main' plotting area
    # Will draw individual axes using this GridSpec object later
    gs = mgridspec.GridSpecFromSubplotSpec(
            nrows         = nrows,
            ncols         = ncols,
            subplot_spec  = GS[0,0],
            wspace        = wspace,
            hspace        = hspace,
            width_ratios  = wratios,
            height_ratios = hratios
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
        try: iter(maps)
        except TypeError:
            maps = [maps] # just want a single map
        else:
            maps = [*maps] # force into list, not array
        maps_ids = [i for i,a in enumerate(np.unique(array).flat) if a in maps]
    else:
        maps_ids = []
    # Find axes that have inner panels
    if innerpanels is not None:
        if innerpanels is True:
            innerpanels = [*np.unique(array)]
        try: iter(innerpanels)
        except TypeError:
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
    # For panels
    panel_kwargs = {'whichpanels':whichpanels,
            'hspace':ihspace, 'wspace':iwspace,
            'hwidth':ihwidth, 'wwidth':iwwidth}

    # Base axes; to be shared with other axes as ._sharex, ._sharey attributes
    axps = num_axes*[None] # empty for now
    axs = num_axes*[None] # list of axes
    allgroups_base = []
    if sharex: allgroups_base += xgroups_base
    if sharey: allgroups_base += ygroups_base
    for i in allgroups_base: # this is just list of indices in axes_ids, yrange, etc.
        if axs[i] is None: # not created as a x-base already, for example
            if i in innerpanels_ids:
                axs[i], axps[i] = _panel_factory(fig, gs[slice(*yrange[i,:]), slice(*xrange[i,:])], width, height,
                        **panel_kwargs) # main axes handle
            else:
                axs[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        **cartopy_kwargs) # main axes can be a cartopy projection

    # Dependent axes
    for i in range(num_axes):
        # Detect if we want to share this axis with another
        # If so, get that axis
        sharex_ax, sharey_ax = None, None # by default, don't share with other axes objects
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
                axs[i], axps[i] = _panel_factory(fig, gs[slice(*yrange[i,:]), slice(*xrange[i,:])], width, height,
                        sharex=sharex_ax, sharey=sharey_ax, **panel_kwargs)
            else:
                axs[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        sharex=sharex_ax, sharey=sharey_ax, **cartopy_kwargs) # main axes can be a cartopy projection

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
            axs[b].xaxis.label.set_transform(mtransforms.blended_transform_factory(
                    fig.transFigure, mtransforms.IdentityTransform()
                    ))
                # specify x, y transform in Figure coordinates
            xmin = min(axs[i].get_position().xmin for i in g)
            xmax = max(axs[i].get_position().xmax for i in g)
                # get min/max positions, in figure coordinates, of spanning axes
            # print('Group:', g, 'Base:', b, 'Span:', xmin, xmax)
            axs[b].xaxis.label.set_position(((xmin+xmax)/2, 0))
                # this is the shared xlabel
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
        try: iter(bottompanels)
        except TypeError:
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
                nrows         = 1,
                ncols         = npanel,
                subplot_spec  = GS[1,0],
                wspace        = rspace, # same as above
                width_ratios  = rratios,
                )
        specs = [P[0,i] for i in range(npanel)] # pass the SubplotSpec objects
        for s in specs:
            s.figure = fig
            s.rightpanel = False
            s.bottompanel = True
            s.colorbar = MethodType(_format_colorbar, s)
            s.legend = MethodType(_format_legend, s)
        if len(specs)==1: specs = specs[0] # no indexing of singleton panels
        fig.bottompanel = specs
            # don't necessarily want to draw axes (e.g. for colorbar
            # need to make another SubplotSpec from each element)
    # Next the rightpanel options (very similar)
    if rightpanel:
        rightpanels = [1]*nrows
    elif rightpanels:
        try: iter(rightpanels)
        except TypeError:
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
                ncols         = 1,
                nrows         = npanel,
                subplot_spec  = GS[0,1],
                hspace        = rspace,
                height_ratios = rratios,
                )
        specs = [P[i,0] for i in range(npanel)] # pass the SubplotSpec objects
        for s in specs:
            s.figure = fig
            s.rightpanel = True
            s.bottompanel = False
            s.colorbar = MethodType(_format_colorbar, s)
            s.legend = MethodType(_format_legend, s)
        if len(specs)==1: specs = specs[0]
        fig.rightpanel = specs
            # don't necessarily want to draw axes (e.g. for colorbar
            # need to make another SubplotSpec from each element)

    #--------------------------------------------------------------------------
    # Create colorbar and legend axes
    #--------------------------------------------------------------------------
    if bottomcolorbar:
        C = mgridspec.GridSpecFromSubplotSpec(
                nrows        = 1,
                ncols        = 3,
                wspace       = 0,
                subplot_spec = GS[1,0],
                width_ratios = ((1-cshrink)/2, cshrink, (1-cshrink)/2)
                )
        c = fig.add_subplot(C[0,1])
        c.bottomlegend, c.bottomcolorbar, c.rightcolorbar = False, True, False
        c.format = MethodType(_format_colorbar, c) # MethodType approach
        fig.bottomcolorbar = c
        # fig.bottomcolorbar = MethodType(_format_colorbar, c)
        # fig.colorbar = c
    if rightcolorbar:
        C = mgridspec.GridSpecFromSubplotSpec(
                nrows         = 3,
                ncols         = 1,
                hspace        = 0,
                subplot_spec  = GS[0,1],
                height_ratios = ((1-cshrink)/2, cshrink, (1-cshrink)/2)
                )
        c = fig.add_subplot(C[1,0])
        c.bottomlegend, c.bottomcolorbar, c.rightcolorbar = False, False, True
        c.format = MethodType(_format_colorbar, c) # MethodType approach
        fig.rightcolorbar = c
        # fig.rightcolorbar = MethodType(_format_colorbar, c)
        # fig.colorbar = c
    # Bottom legend panel
    if bottomlegend:
        l = fig.add_subplot(GS[1,0])
        l.bottomlegend, l.bottomcolorbar, l.rightcolorbar = True, False, False
        for s in l.spines.values():
            s.set_visible(False)
        l.xaxis.set_visible(False)
        l.yaxis.set_visible(False)
        l.patch.set_alpha(0)
        l.format = MethodType(_format_legend, l) # MethodType approach
        fig.bottomlegend = l # had to set some properties first
        # fig.bottomlegend = MethodType(_format_legend, l)
        # fig.legend = l

    # Set up attributes and methods
    # or is empty string, methods in _format_axes will make them visible
    for i,ax in enumerate(axs):
        _atts_global(ax) # default methods and stuff
        _atts_special(ax, maps, package, projection, **projection_kwargs) # default methods and stuff
        ax.number = i+1 # know axes number ahead of time; start at 1

    # Repeat some of the above for the panel axes
    # Will only need the format method really
    # Also eliminate from the list axes for which no panels are drawn
    axps = [axp for axp in axps if axp is not None]
    for axp in axps:
        for ax in axp:
            _atts_global(ax)

    #--------------------------------------------------------------------------
    # Return results
    # Will square singleton arrays
    #--------------------------------------------------------------------------
    if not silent:
        print('Figure setup complete.')
    if len(axs)==1:
        axs = axs[0]
    axps = [axp[0] if len(axp)==1 else axp for axp in axps]
    if len(axps)==1: # will be list if have panels on more than one side
        axps = axps[0]
    if innerpanels:
        return fig, axs, axps
    else:
        return fig, axs

