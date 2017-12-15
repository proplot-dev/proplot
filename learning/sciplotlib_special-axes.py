#------------------------------------------------------------------------------
# Imports, all
#------------------------------------------------------------------------------
# from .basics import Dict, arange
import os # need for saving figures
from types import MethodType
from fractions import Fraction
from string import ascii_uppercase
from cycler import cycler
from glob import glob
import matplotlib.figure as mfigure
# import matplotlib.axes as maxes
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import matplotlib.text as mtext
import matplotlib.ticker as mticker
import matplotlib.gridspec as mgridspec
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as mbasemap
# import string # for converting number to a/b/c
# import matplotlib.pyplot as plt # will attach a suitable backend, make available the figure/axes modules
import numpy as np # of course
import cartopy.crs as ccrs # crs stands for "coordinate reference system", leading c is "cartopy"
import cartopy.mpl as cmpl
import cartopy.feature as cfeature
# __all__ = [
#     'Figure',
#     'subplots', 'cmapload', 'settings', # basic setup
#     'Normalize', 'Formatter', 'FracFormatter', # custom colormap normalization
#     ]
# Sample fixes for white lines between plots
# c = ax.contourf()
# for _ in c.collections: _.set_edgecolor('face')
# p = ax.pcolormesh(linewidth=0)
# p.set_edgecolor('face')

#------------------------------------------------------------------------------
# Colormap stuff
#------------------------------------------------------------------------------
# Register colormaps immediately on import
_announcement = False
for _file in glob(f'{os.path.dirname(__file__)}/cmaps/*.rgb'):
    # Load each colormap
    _name = os.path.basename(_file).replace('.rgb','')
    if _name not in plt.colormaps(): # don't want to re-register every time
        _load = {'skiprows':1, 'delimiter':','} if _name.startswith('hcl') else {}
        try: _cmap = np.loadtxt(_file, **_load)
        except:
            print(f'Failed to load {_name}.')
            continue
        if (_cmap>1).any(): _cmap = _cmap/255
        plt.register_cmap(name=_name, cmap=mcolors.LinearSegmentedColormap.from_list(name=_name, colors=_cmap, N=256))
        plt.register_cmap(name=_name+'_r', cmap=mcolors.LinearSegmentedColormap.from_list(name=_name+'_r', colors=_cmap[::-1], N=256))
        if not _announcement:
            _announcement = True
            print(f'Registered colormaps.')
# Plot all currently registered colormaps
def cmapshow(ignore=('Qualitative','Miscellaneous','Sequential Alt')):
    """
    Plot all current colormaps, along with their catgories.
    This example comes from the Cookbook on www.scipy.org.  According to the
    history, Andrew Straw did the conversion from an old page, but it is
    unclear who the original author is.
    """
    # Have colormaps separated into categories:
    # http://matplotlib.org/examples/color/colormaps_reference.html
    categories = { 'HCL': [], 'Custom': [], # will add to these lists
            'Perceptually Uniform Sequential': ('viridis', 'plasma', 'inferno', 'magma'),
            'Diverging': ('PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
                'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'),
            'Sequential': ('Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn'),
            'Sequential Alt': ('binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
                'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
                'hot', 'afmhot', 'gist_heat', 'copper'),
            'Qualitative': ('Pastel1', 'Pastel2', 'Paired', 'Accent',
                'Dark2', 'Set1', 'Set2', 'Set3', 'Vega10', 'Vega20', 'Vega20b', 'Vega20c'),
            'Miscellaneous': ('flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
                    'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv', 'spectral',
                    'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar') }
    cmaps = [m for m in plt.colormaps() if not m.endswith('_r')]
    cmaps_known = []
    for v in categories.values(): # big list of all colormaps
        cmaps_known += v # add to this existing list
    for cm in cmaps: # add to 'Custom' if not in above dictionary
        if cm not in cmaps_known:
            if cm.startswith('hcl'):
                categories['HCL'].append(cm)
            else:
                categories['Custom'].append(cm)
    # Array for producing visualization with imshow
    a = np.linspace(0, 1, 257).reshape(1,-1)
    a = np.vstack((a,a))
    # Figure
    twidth = 2 # number of axes-widths to allocate for titles
    nmaps = len(cmaps) + len(categories)*twidth # title for each category
    for cat in ignore:
        nmaps -= (twidth + len(categories[cat])) # reduce size
    fig = plt.figure(figsize=(5,.3*nmaps))
    fig.subplots_adjust(top=0.99, bottom=0.01, left=0.2, right=0.99)
    # Make plot
    ntitles, nplots = 0, 0 # for deciding which axes to plot in
    cats = [cat for cat in categories if cat not in ignore]
    for cat in cats:
        # Space for title
        ntitles += twidth # two axes-widths
        for i,m in enumerate(categories[cat]):
            # Draw, and make invisible
            ax = plt.subplot(nmaps,1,i+ntitles+nplots)
            for s in ax.spines.values():
                s.set_visible(False)
            ax.patch.set_alpha(0)
            ax.set_xticks([])
            ax.set_yticks([])
            # Draw colormap
            ax.imshow(a, aspect='auto', cmap=m, origin='lower')
            if cat not in ('Custom','HCL'):
                cmaps_known.remove(m)
            # Category title
            if i==0:
                t = ax.title
                t.set_text(cat) # category name
                t.set_visible(True)
            # Label
            yl = ax.yaxis.label
            yl.set_text(m) # map name
            yl.set_visible(True)
            yl.set_rotation(0)
            yl.set_ha('right')
            yl.set_va('center')
        # Space for plots
        nplots += len(categories[cat])
    # Check
    for cat in ignore:
        for m in categories[cat]:
            if cat not in ('Custom','HCL'):
                cmaps_known.remove(m)
    if len(cmaps_known)>0:
        print(f'Colormaps in dictionary, but not found: {", ".join(cmaps_known)}')
    # Save
    fig.savefig(f'{os.path.dirname(__file__)}/colormaps.pdf', format='pdf')
    return

#------------------------------------------------------------------------------
# Important class
#------------------------------------------------------------------------------
class Figure(mfigure.Figure):
    """
    Subclass of the mfigure.Figure class, with lots of special formatting
    options. Can be called by using pyplot.figure(FigureClass=Figure) kwargument
    in my subplots function.
    """
    def __init__(self, *args, **kwargs):
        # Initialize
        # mfigure.Figure.__init__(self, *args, **kwargs)
        # super(Figure, self).__init__(*args, **kwargs)
        super().__init__(*args, **kwargs) # python 3 only
        # Special axes
        self.axes_main = []
        self.bottomcolorbar = None # overall colorbar; replace with axis object
        self.rightcolorbar = None
        self.bottomlegend = None # overall legend; replace with axis object
        self.title = None
    
    def save(self, filename, desktop=True):
        """
        Saving figures.
        """
        if '.' not in filename:
            filename = filename + '.pdf'
        if desktop: 
            filename = os.environ['HOME'] + '/Desktop/figures/' + filename
        print('Saving figure, filename: %s' % filename)
        self.savefig(filename, format='pdf', pad_inches=0)
        # self.savefig(filename, format='pdf', bbox_inches='tight', pad_inches=0.05, dpi=300,
        #         frameon=False)
        #     # pad_inches is padding around bbox, and frameon makes figure window transparent
        #     # even though saving as pdf, think dpi will affect any embedded raster objects

#------------------------------------------------------------------------------
# Helper functions
#------------------------------------------------------------------------------
def _round(x, base=5):
    # Round to nearest N
    return base*round(float(x)/base)

def _autolocate(min_, max_, base=5):
    # Return auto-generated levels by provided interval
    return np.arange(_round(min_,base), _round(max_,base)+base/2, base)

def _graticule_fix(x, y):
    # Get graticule; want it to work with e.g. Gaussian-grid model output
    # dx, dy = x[1]-x[0], y[1]-y[0]
    # xb, yb = np.arange(x[0]-dx/2, x[-1]+dx, dx), np.arange(y[0]-dy/2, y[-1]+dy, dy)
    #     # can just use arange, because don't care about propagating float errors
    # return xb, yb
    b = lambda vec: np.concatenate([vec[0]-np.diff(vec[:2])/2,
                (vec[:1]+vec[:-1])/2,
                vec[-1]+np.diff(vec[-2:])/2])
    return b(x), b(y)

def _seam_fix(lon, data, lonmin=-180):
    # Raise errors
    if lon.max()>lon.min()+360:
        print(lon.min(), lon.max())
        print(lon.max() - lon.min()-360)
        raise ValueError('Longitudes must span 360 degrees at most.')
    if lon.min()<-360 or lon.max()>360:
        raise ValueError('Longitudes must fall in range [-360, 360].')
    if lonmin<-360 or lonmin>0:
        raise ValueError('Minimum longitude must fall in range [-360, 0].')

    # Establish 360-degree range
    lon -= 720
    while True:
        filter_ = lon<lonmin
        if filter_.sum()==0:
            break
        lon[filter_] += 360
    
    # Roll, accounting for whether ends are identical
    roll = -np.argmin(lon) # always returns FIRST value; if go from 0...359,0 (borders), will
        # return initial value
    if lon[0]==lon[-1]:
        lon = np.roll(lon[:-1], roll)
        lon = np.append(lon, lon[0]+360)
    else:
        lon = np.roll(lon, roll)
    data = np.roll(data, roll, axis=0)

    # Fix seams
    if lon[0]==lonmin:
        print('No seam fix necessary. Size augmented by 0.')
        return lon, data
    elif lon.size-1==data.shape[0]: # want to avoid testing floats
        print('Fixing seam using lon borders. Size augmented by 1.')
        lon = np.append(lonmin, lon) # append way easier than concatenate
        lon[-1] = lonmin+360 # we've added a new tiny cell to the end
        if data is not None:
            data = np.concatenate((data[-1:,...], data), axis=0) # don't use pad, because messes up masked arrays
            # pad_tuple = ((1,0),) + ((0,0),)*(data.ndim-1)
            # data = np.pad(data, pad_tuple, 'wrap') # pad before
    elif lon.size==data.shape[0]:
        print('Fixing seam using lon centers/interpolation. Size augmented by 2.')
        if data is not None:
            x = np.array([lon[-1], lon[0]+360]) # x
            y = np.concatenate((data[-1:,...], data[:1,...]), axis=0)
            xq = lonmin+360
            yq = (y[:1,...]*(x[1]-xq) + y[1:,...]*(xq-x[0]))/(x[1]-x[0]) # simple linear interp formula
            data = np.concatenate((yq, data, yq), axis=0)
        lon = np.append(np.append(lonmin, lon), lonmin+360)
    else:
        raise ValueError('Longitude length does not match data length along dimension 0.')

    # Return
    return lon, data

def _basemap_fix(m, lon, lat):
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

#------------------------------------------------------------------------------
# Plot method overrides
#------------------------------------------------------------------------------
def _pcolorcheck(x, y, Z):
    x, y, z = np.array(x), np.array(y), np.array(z)
    if Z.shape[0]==x.size and Z.shape[1]==y.size:
        x, y = _graticule_fix(x, y)
    elif Z.shape[0]!=x.size-1 or Z.shape[1]!=y.size-1:
        raise ValueError(f'X (size {x.size}) and Y (size {y.size}) must correspond to '
                f'nrows ({Z.shape[0]}) and ncolumns ({Z.shape[1]}) of Z, or its borders.')
    return x, y
def _contourcheck(x, y, Z):
    x, y, z = np.array(x), np.array(y), np.array(z)
    if Z.shape[0]==x.size-1 and Z.shape[1]==y.size-1:
        x, y = (x[1:]+x[:-1])/2, (y[1:]+y[:-1])/2
        print(x,y)
    elif Z.shape[0]!=x.size or Z.shape[1]!=y.size:
        raise ValueError(f'X (size {x.size}) and Y (size {y.size}) must correspond to '
                f'nrows ({Z.shape[0]}) and ncolumns ({Z.shape[1]}) of Z, or its borders.')
    return x, y

# Ordinary cartesian plot overrides
def _pcolor(self, x, y, Z, **kwargs):
    x, y = _pcolorcheck(x, y, Z)
    p = self._pcolor(x, y, Z.T, **kwargs)
    return _pcolor_fix(p)
def _pcolormesh(self, x, y, Z, **kwargs):
    x, y = _pcolorcheck(x, y, Z)
    p = self._pcolormesh(x, y, Z.T, **kwargs)
    return _pcolor_fix(p)
def _contourf(self, x, y, Z, **kwargs):
    x, y = _contourcheck(x, y, Z)
    c = self._contourf(x, y, Z.T, **kwargs)
    return _contour_fix(c)
def _contour(self, x, y, Z, **kwargs):
    x, y = _contourcheck(x, y, Z)
    c = self._contour(x, y, Z.T, **kwargs)
    return c

# Basemap overrides; assumes regularly-spaced data
def _pcolor_basemap(self, lon, lat, Z, **kwargs):
    # lon, lat = _graticule_fix(lon, lat)
    # lon, lat = _graticule_fix(lon, lat)
    lon, lat = _pcolorcheck(lon, lat, Z)
    lon, Z = _seam_fix(lon, Z, lonmin=self.lonmin)
    X, Y = _basemap_fix(self, lon, lat)
    p = self._pcolor(X, Y, Z.T, **kwargs)
    return _pcolor_fix(p)
def _pcolormesh_basemap(self, lon, lat, Z, **kwargs):
    # lon, lat = _graticule_fix(lon, lat)
    lon, lat = _pcolorcheck(lon, lat, Z)
    lon, Z = _seam_fix(lon, Z, lonmin=self.lonmin)
    X, Y = _basemap_fix(self, lon, lat)
    p = self._pcolormesh(X, Y, Z.T, **kwargs)
    return _pcolor_fix(p)
def _contourf_basemap(self, lon, lat, Z, **kwargs):
    lon, lat = _contourcheck(lon, lat, Z)
    lon, Z = _seam_fix(lon, Z, lonmin=self.lonmin)
    X, Y = _basemap_fix(self, lon, lat)
    c = self._contourf(X, Y, Z.T, **kwargs)
    return _contour_fix(c)
def _contour_basemap(self, lon, lat, Z, **kwargs):
    lon, lat = _contourcheck(lon, lat, Z)
    lon, Z = _seam_fix(lon, Z, lonmin=self.lonmin)
    X, Y = _basemap_fix(self, lon, lat)
    c = self._contour(X, Y, Z.T, **kwargs)
    return c

#------------------------------------------------------------------------------
# Most important stuff
#------------------------------------------------------------------------------
class Settings():
    """
    This function initializes DEFAULT kwargs for various plot elements; user can
    UPDATE the defaults by passing dictionaries as kwargs when creating instance.
    """
    def __init__(self, color='k', linewidth=0.7, transparent=False, **kwargs):
        # Axes look (updated)
        self.tick      = {'length':4, 'direction':'in', 'width':linewidth, 'color':color}
        self.tickminor = {'length':2, 'direction':'in', 'width':linewidth, 'color':color}
        self.spine     = {'linewidth':linewidth, 'color':color} # spine keyword
        self.box       = {'linewidth':linewidth, 'edgecolor':color}
        self.ctick     = {'length':4, 'direction':'out', 'width':linewidth, 'color':color}
        self.cgrid    = {'color':color, 'linewidth':linewidth}
            # for outlines of colorbar, legend, cartopy frame, and cartopy features
        # Text (updated, or set on creation)
        self.ticklabels = {'size':8, 'weight':'normal', 'color':color}
        self.label      = {'size':8, 'weight':'normal', 'color':color}
        self.title      = {'size':10, 'weight':'normal', 'color':color}
        self.abc        = {'size':10, 'weight':'bold', 'color':color}
        # Grid lines (set on creation; should make half-as-wide as axes)
        self.grid      = {'linestyle':'-', 'linewidth':linewidth/2, 'color':color, 'alpha':0.1}
        self.gridminor = {'linestyle':':', 'linewidth':linewidth/2, 'color':color, 'alpha':0.05}
        # Axes scale special properties (set on creation)
        self.xscale = {}
        self.yscale = {}
        # Geographic features look (set on creation)
        self.coastlines  = {'linewidth':linewidth, 'color':'#888888'}
        self.continents  = {'color':'#CCCCCC'}
        self.lonlatlines = {'linewidth':linewidth, 'color':color, 'alpha':0.2, 'dashes':(linewidth,linewidth)}
        self.boundary    = {'linewidth':2*linewidth, 'color':color}
            # boundary should be thick... for some reason looks skinnier than parallels unless multiply by 2
            # even though (checking after the fact) reveals the linewidths are held
        # Auxilliary elements (set on creation)
        self.legend   = {'framealpha':0.6, 'fancybox':False, 'frameon':True, 'fontsize':8}
        self.colorbar = {'extend':'both', 'spacing':'uniform'}
        # User-specified updated properties (provided as 'attribute_kw')
        # for attr in (attr for attr in dir(self) if not attr.startswith('__')):
        for attr, settings in self.__dict__.items():
            settings.update(kwargs.pop(attr+'_kw', {}))
        if kwargs: # non-empty
            raise ValueError(f'Unknown kwarg-settings: {", ".join(kwargs.keys())}.')

def _lformat(self, 
        ):
    """
    Function for formatting colorbar-axes (axes that are "filled" by a colorbar).
    """

def _lformat(self,
        ):
    """
    Function for formatting legend-axes (invisible axes with centered legends on them).
    Should update my legend function to CLIP the legend box when it goes outside axes area, so
    the legend-width and bottom/right widths can be chosen propertly/separately.
    """

def _format(self, transparent=True, 
    legend=False, handles=None, # legend options; if declared and have no handles with 'label', does nothing
    mappable=None, # for colorbar
    coastlines=True, continents=False, latlocator=None, lonlocator=None, # latlocator/lonlocator work just like xlocator/ylocator
    cgrid=None, xgrid=None, ygrid=None, # grid minor; cgrid means, add dividers to colorbars, gets passed as drawedges kwarg
    xdates=False, ydates=False, # whether to format axis labels as long datetime strings
    xtickminor=None, ytickminor=None, xgridminor=None, ygridminor=None, # minor ticks/grids; if ticks off, grid will be off
    xspineloc=None, yspineloc=None, # deals with spine options
    xtickdir=None, ytickdir=None, # change ytick/xtick location; can be in, out, or inout (left-right-up-down depends on spine to which this applied)
    xlim=None, ylim=None, xscale=None, yscale=None, xreverse=False, yreverse=False, # special properties
    xlabel=None, ylabel=None, clabel=None, # axis/colorbar labels
    title=None, titleinside=None, titleleft=None, titleright=None, abcinside=False, abc=False, # titles
    xlocator=None, xminorlocator=None, clocator=None, ylocator=None, yminorlocator=None, # locators, or derivatives that are passed to locators
    xformatter=None, cformatter=None, yformatter=None, # formatter
    **kwargs): # remaining are passed to Settings()
    """
    Function for formatting axes of all kinds; some arguments are only relevant to special axes, like 
    colorbar or basemap axes. By default, simply applies the dictionary values from settings() above,
    but can supply many kwargs to further modify things.

    TODO...
    * Add options for special legend handling; allow overwriting the label linewidth, etc.
    if we want to differentiate colors, for example.
    * Add options for datetime handling; note possible date axes handles are TimeStamp (pandas),
    np.datetime64, DateTimeIndex; can fix with fig.autofmt_xdate() or manually set options; uses
    ax.is_last_row() or ax.is_first_column(), so should look that up...problem is there is no 
    autofmt_ydate(), so really should implement my own version of this...
    * Add options/considerations for twinx/twiny; since this iterates through all axes children in the figure object, 
    should be valid. Should be **added to my.subplots**, actually, and then my.subplots can handle their placement in the array.
    """
    
    # Load default settings, and override with any user-provided dictionaries
    if any(type(kw) is not dict for kw in kwargs.values()):
        bad_args = (name for name,kw in kwargs.items() if type(kw) is not dict)
        raise ValueError(f'Unknown arguments: {", ".join(bad_args)}. '
                'Remember, can change settings by passing dictionaries of new settings.')
    props = Settings(**kwargs)

    #--------------------------------------------------------------------------
    # Special axes setup
    #--------------------------------------------------------------------------
    # Legend axes setup
    if self.bottomlegend: # e.g. can input -1, or axes.Axes object
        props.legend.update(handles=handles, loc='upper center', framealpha=1) # make opaque; no lines underneath
    elif handles is not None:
        props.legend.update(handles=handles)
    if legend or self.bottomlegend:
        # Draw legend
        handleprops = dict()
        candidates = ('linewidth',) # candidates for modifying legend objects
        for candidate in candidates:
            if candidate in props.legend: handleprops[candidate] = props.legend.pop(candidate)
        leg = self.legend(**props.legend)

        # Change properties of legend box (colorbar gets set to same properties),
        # legend lines, and legend text
        leg.legendPatch.update(props.box) # or get_frame()
        for obj in leg.legendHandles:
            obj.update(handleprops)
        for t in leg.texts: t.update(props.ticklabels) # or get_texts()
        if self.bottomlegend:
            return

    #--------------------------------------------------------------------------
    # Titles (should be done before projection stuff)
    #--------------------------------------------------------------------------
    # Main title
    self.title.update(props.title)
    if title is not None:
        self.title.set_text(title)

    # Inner title options (left, right, and centered)
    self.titleinside.update(props.title)
    self.titleright.update(props.title)
    self.titleleft.update(props.title)
    if titleinside is not None:
        self.titleinside.set_text(titleinside)
    if titleleft is not None:
        self.titleleft.set_text(titleleft)
    if titleright is not None:
        self.titleright.set_text(titleright)

    # ABC labelling
    self.abc.update(props.abc) # should have an 'abc' property from subplots method
    self.abcinside.update(props.abc)
    if abc:
        self.abc.set_visible(True) # turn on
    if abcinside:
        self.abcinside.set_visible(True) # turn on

    #--------------------------------------------------------------------------
    # Projections
    #--------------------------------------------------------------------------
    # Basemap axes setup
    if self.m is not None:
        # Coastlines, parallels, meridians
        m = self.m
        if coastlines:
            p = m.drawcoastlines(**props.coastlines)
        if continents:
            p = m.fillcontinents(**props.continents)

        # Longitude/latitude lines
        # Make sure to turn off clipping by invisible axes boundary; otherwise
        # get these weird flat edges where map boundaries, parallel/meridian markers come up to the axes bbox
        # ...setup
        labels = [1,1,1,1]
        labels = [0,0,0,0]
        if latlocator is not None:
            try: iter(latlocator)
            except TypeError:
                latlocator = _autolocate(m.latmin+latlocator, m.latmax-latlocator, latlocator)
            p = m.drawparallels(latlocator, labels=labels) 
            for pi in p.values(): # returns dict, where each one is tuple
                for _ in [i for j in pi for i in j]: # magic
                    if isinstance(_, mtext.Text):
                        _.update(props.ticklabels)
                    else:
                        _.set_clip_on(False)
                        _.update(props.lonlatlines)
                        # _.set_linestyle(linestyle)
                # tried passing clip_on to the above, but it does nothing; must set
                # for lines created after the fact
        if lonlocator is not None:
            try: iter(lonlocator)
            except TypeError:
                lonlocator = _autolocate(m.lonmin+lonlocator, m.lonmax-lonlocator, lonlocator)
            p = m.drawmeridians(lonlocator, labels=labels)
            for pi in p.values():
                for _ in [i for j in pi for i in j]: # magic
                    if isinstance(_, mtext.Text):
                        _.update(props.ticklabels)
                    else:
                        _.set_clip_on(False)
                        _.update(props.lonlatlines)
                        # _.set_linestyle(linestyle)

        # Map boundary
        p = m.drawmapboundary(**props.boundary)
        p.set_facecolor(None)
        p.set_zorder(-1)
        p.set_clip_on(False)
        p.set_rasterized(False) # not sure about this... might be rasterized
        return # skip everything else

    # # Cartopy axes setup
    # if isinstance(self, cmpl.geoaxes.GeoAxes): # the main GeoAxes class; others like GeoAxesSubplot subclass this
    #     self.add_feature(cfeature.COASTLINE, **props.box)
    #     plt.setp(self.outline_patch, **props.box)
    #     return

    #--------------------------------------------------------------------------
    # X/Y axis stuff, label management
    #--------------------------------------------------------------------------
    # Axes scaling, limits, and reversal options (alternatively, supply
    # ...your own xlim/ylim that go from high to low)
    if xscale is not None: self.set_xscale(xscale, **props.xscale)
    if yscale is not None: self.set_yscale(yscale, **props.yscale)
    if xlim is None: xlim = self.get_xlim()
    if ylim is None: ylim = self.get_ylim()
    if xreverse: xlim = xlim[::-1]
    if yreverse: ylim = ylim[::-1]
    self.set_xlim(xlim)
    self.set_ylim(ylim)

    for axis, label, dates, sides, spineloc, gridminor, tickminor, tickminorlocator, grid, ticklocator, tickformatter, tickdir in zip(
            (self.xaxis, self.yaxis), (xlabel, ylabel), (xdates, ydates), (('bottom','top'),('left','right')), (xspineloc, yspineloc), # other stuff
            (xgridminor, ygridminor), (xtickminor, ytickminor), (xminorlocator, yminorlocator), # minor ticks
            (xgrid, ygrid), (xlocator, ylocator), (xformatter, yformatter), # major ticks
            (xtickdir, ytickdir), # tick direction
            ):
        # Axis spine and tick locations
        for spine, side in zip((self.spines[s] for s in sides), sides):
            # Line properties
            spine.update(props.spine)
            # Standard
            if spineloc=='both':
                spine.set_visible(True)
            elif spineloc in sides: # make relevant spine visible
                b = True if side==spineloc else False
                spine.set_visible(b)
            elif spineloc is not None:
                # Special location
                if side==sides[0]: # move the left/bottom spine onto the specified location, with set_position
                    spine.set_visible(True)
                    spine.set_position(spineloc) # just gets passed to the function; options include 
                        # 'zero', 'center', and tuple with (units, location) where units can be axes, data, or outward
                else:
                    spine.set_visible(False)

        # Tick/tick gridline properties
        # Some weird issue seems to cause set_tick_params to reset/forget that the grid
        # is turned on if you access tick.gridOn directly, instead of passing through tick_params.
        # ...so calling _format() a second time will remove the lines
        sidesOn = {side: self.spines[side].get_visible() for side in sides}
        gridOn = dict() if grid is None else dict(gridOn=grid)
        gridminorOn = dict() if gridminor is None else dict(gridOn=gridminor and tickminor)
        majorSet = props.tick if tickdir is None else dict(props.tick, direction=tickdir)
        minorSet = props.tickminor if tickdir is None else dict(props.tickminor, direction=tickdir)
        # Apply
        axis.set_tick_params(which='major', **sidesOn, **gridOn, **majorSet)
        axis.set_tick_params(which='minor', **sidesOn, **gridminorOn, **minorSet) # have length
        # print(axis.minorTicks)
        # print(sidesOn, gridminorOn, minorSet)
        # print(tickminor, sidesOn)
        for tick in axis.majorTicks:
            tick.gridline.update(props.grid)
        for tick in axis.minorTicks:
            tick.gridline.update(props.gridminor)
            if tickminor is not None:
                tick.set_visible(tickminor)

        # Label properties
        axis.label.update(props.label)
        if label is not None:
            axis.label.set_text(label)

        # First, major tick locators (should not affect limits)
        lim = axis.get_view_interval() # to be used, automatically
        if isinstance(ticklocator, mticker.Locator):
            axis.set_major_locator(ticklocator)
        elif ticklocator is not None:
            try: iter(ticklocator)
            except TypeError: ticklocator = _autolocate(lim[0], lim[1], ticklocator)
            # axis.set_ticks(ticklocator)
            axis.set_major_locator(mticker.FixedLocator(ticklocator))

        # Next, minor tick locators (toggle visibility later)
        if isinstance(tickminorlocator, mticker.Locator):
            axis.set_minor_locator(tickminorlocator) # pass locator class/subclass (all documented ones subclass the parent class, mticker.Locator)
        elif tickminorlocator is not None:
            try: iter(tickminorlocator)
            except TypeError: tickminorlocator = _autolocate(lim[0], lim[1], tickminorlocator) # use **spacing/setp** as specified
            axis.set_minor_locator(mticker.FixedLocator(tickminorlocator))

        # Next, major tick formatters (enforce Null, always, for minor ticks), and text styling
        axis.set_minor_formatter(mticker.NullFormatter())
        if isinstance(tickformatter, mticker.Formatter): # use isinstance, because r"stuff" or f"stuff" **do not work**
            axis.set_major_formatter(tickformatter)
        elif tickformatter is not None:
            if isinstance(tickformatter, str):
                if dates: axis.set_major_formatter(mdates.DateFormatter(tickformatter)) # %-style, dates
                else: axis.set_major_formatter(mticker.FormatStrFormatter(tickformatter)) # %-style, numbers
            else:
                # axis.set_ticklabels(tickformatter) # ...FixedFormatter alone has issues
                axis.set_major_formatter(mticker.FixedFormatter(tickformatter)) # list of strings
        for t in axis.get_ticklabels():
            t.update(props.ticklabels)

    # # ...try again?
    # if xscale is not None: self.set_xscale(xscale, **props.xscale)
    # if yscale is not None: self.set_yscale(yscale, **props.yscale)
    # if xlim is None: xlim = self.get_xlim()
    # if ylim is None: ylim = self.get_ylim()
    # if xreverse: xlim = xlim[::-1]
    # if yreverse: ylim = ylim[::-1]
    # self.set_xlim(xlim)
    # self.set_ylim(ylim)
    return # we're done

def subplots(array=None, nrows=1, ncols=1,
        tight=False, # whether to set up tight bbox from gridspec object
        transparent=True, # whether figure should be transparent
        sharex=True, sharey=True, # for sharing x/y axis limits/scales/locators for axes with matching GridSpec extents, and making ticklabels/labels invisible
        spanx=True, spany=True, # custom setting; share axis labels for axes with same xmin/ymin extents
            # what about bottom vs. top xlabels, left vs. right ylabels?
        aspect=None, height=None, width=5, # for controlling aspect ratio; default is control for width
        hspace=0.4, wspace=0.2, height_ratios=None, width_ratios=None, # options for gridspec (instead of gridpsec_kw)
            # spacing betwen axes, in inches; generally horizontal-spacing small, but vertical-spacing big to allow for title
        left=0.8, bottom=0.5, right=0.1, top=0.3,
            # edges of plot for GridSpec object, in inches; left needs a lot of space
        bottomcolorbar=False, rightcolorbar=False, bottomlegend=False,
            # for legend/colorbar referencing multiple subplots
        cspace=0.8, shrink=0.9, colorbar_width=0.2, legend_width=0.3, abcpad=0.1, titlepad=0.1,
            # spacing for colorbar text, colorbar axes width, legend axes with, and padding for interior ABC label
        package='basemap', projection=None, **projection_kw): # for projections; can be 'basemap' or 'cartopy'
    """
    Special creation of subplots grids, allowing for arbitrarily overlapping 
    axes objects. Will return figure handle and axes object. Can also use exaclty as
    plt.subplot, with some added convenience features. Use panelx, panely as 
    templates for figure-wide legends/colorbars; if you actually want ax smallish panel
    to plot data, just use width/height ratios. don't bother with kwarg dictionaries because
    it is really, really customized; we do use them in formatting function though.
    
    TODO:
        * if want **differing spacing** between certain subplots, must use 
        GridSpecFromSubplotSpec with its **own** hspace/wspace properties...requires
        some consideration, but don't worry until you come across it
        * figure out how to automate fixed aspect ratios for basemap subplots, so
        don't have to guess/redraw a bunch of times
        * figure size should be constrained by WIDTH/HEIGHT/ASPECT RATIO OF AXES, 
        with everything else provided; must pick two of these; default, is to pick
        width and aspect, but if height provided, width is ignored instead. no way
        right now to determine aspect by constraining width/height
        * Start stepping away from SubplotSpec/GridSpec, and make axes from scratch.
        ...should then be able to seamlessly make differing wspaes/hspaces.
        * For spanning axes labels, right now only dectect **x labels on bottom**
        and **ylabels on top**; generalize for all subplot edges
    
    NOTE that, because I usually format all subplots with .format() method, shared axes should already
    end up with the same axis limits/scaling/majorlocators/minorlocators... the sharex and sharey detection algorithm
    really is just to get INSTRUCTIONS to make the ticklabels and axis labels INVISIBLE for certain axes.
    """

    #--------------------------------------------------------------------------
    # Panel considerations; keep things general for future imporvements
    #--------------------------------------------------------------------------
    # Detect conflicts
    if bottomcolorbar and rightcolorbar:
        raise ValueError('You can only specify one global colorbar.')
    if bottomcolorbar and bottomlegend:
        raise ValueError('Can only place a global colorbar **or** global legend on bottom, not both.')

    # Spacings due to panel considerations
    if bottomcolorbar:
        bottom_extra_axes   = colorbar_width
        bottom_extra_hspace = cspace
    elif bottomlegend:
        bottom_extra_axes   = legend_width
        bottom_extra_hspace = 0 # don't need space
    else:
        bottom_extra_axes, bottom_extra_hspace = 0, 0
    bottom_extra = bottom_extra_axes + bottom_extra_hspace

    # And right spacings
    if rightcolorbar:
        right_extra_axes   = colorbar_width
        right_extra_wspace = cspace
    else:
        right_extra_axes, right_extra_wspace = 0, 0
    right_extra = right_extra_axes + right_extra_wspace
     
    #--------------------------------------------------------------------------
    # Array setup
    #--------------------------------------------------------------------------
    # Automatically generate array of first *arg not provided, or use nrows/ncols
    if array is None:
        array = np.arange(1,nrows*ncols+1)[...,None]
        array = array.reshape((nrows, ncols))
    array = np.array(array)
    nrows = array.shape[0]
    ncols = array.shape[1]
    
    # Projection setup stuff
    # First make sure user didn't mess up
    if len(projection_kw)>0 and projection is None:
        raise ValueError(f"Unknown arguments: {', '.join(projection_kw.keys())}.")
    # For both basemap/cartopy, determine default aspect ratio
    if projection is not None and aspect is None:
        aspect = 2 # try this
        # aspects = (yrange[:,1]-yrange[:,0])/(xrange[:,1]-xrange[:,0])
        # if not np.all(aspects[0]==aspects): # should ammend to accept width/height ratios
        #     raise ValueError('If each subplot is map projection, their aspect ratios should match.')
        # aspect = .5*aspects[0] # .6 is good bet
    elif aspect is None:
        aspect = 1.5 # try this as a default; square plots (1) look too tall to me
    # Cartopy stuff; get crs instance, and create dictionary to add to add_subplot calls
    if projection is not None and package=='cartopy':
        crs_dict = {
                'rectilinear':    ccrs.PlateCarree,
                'pcarree':        ccrs.PlateCarree,
                'cyl':            ccrs.PlateCarree, # matches basemap
                'platecarree':    ccrs.PlateCarree,
                'robinson':       ccrs.Robinson,
                'stereo':         ccrs.Stereographic,
                'stereographic':  ccrs.Stereographic,
                'moll':           ccrs.Mollweide,
                'mollweide':      ccrs.Mollweide,
                'aeqd':           ccrs.AzimuthalEquidistant,
                'np':             ccrs.NorthPolarStereo,
                'sp':             ccrs.SouthPolarStereo,
                }
        cartopy_kw = {'projection': crs_dict[projection](**projection_kw)}
    else:
        cartopy_kw = {}

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
        axwidth_ave_nopanel = (width - left - right - (ncols-1)*wspace - right_extra)/ncols
        if axwidth_ave_nopanel<0:
            raise ValueError('Not enough room for axes. Reduce left/bottom/wspace.')
        axheight_ave_nopanel = axwidth_ave_nopanel/aspect
        height               = axheight_ave_nopanel*nrows + top + bottom + (nrows-1)*hspace + bottom_extra

    # Figure size, and some other stuff
    axwidth_total  = ncols*axwidth_ave_nopanel + (ncols-1)*wspace
    axheight_total = nrows*axheight_ave_nopanel + (nrows-1)*hspace
    figsize        = (width, height)

    # Convert hspace/wspace/edges to INCHES units
    # ...properties for outer GridSpec object
    left = left/width
    top  = 1-top/height
    if bottomlegend or bottomcolorbar:
        hspace_outer        = bottom/((axheight_total + bottom_extra_axes)/2) # same for main axes and bottom panel, but 'bottom'
        height_ratios_outer = np.array([axheight_total, bottom_extra_axes])/(axheight_total + bottom_extra_axes)
        bottom              = bottom_extra_hspace/height
    else:
        hspace_outer        = 0
        height_ratios_outer = np.array([1])
        bottom              = bottom/height
    if rightcolorbar:
        width_ratios_outer = np.array([axwidth_total, right_extra_axes])/(axwidth_total + right_extra_axes)
        wspace_outer       = right/((axwidth_total + right_extra_axes)/2) # want space between main axes and panel to be 'right'
        right              = 1-right_extra_wspace/width
    else:
        width_ratios_outer = np.array([1])
        wspace_outer       = 0
        right              = 1-right/width
    # ...properties for inner GridSpec object
    wspace = wspace/axwidth_ave_nopanel
    hspace = hspace/axheight_ave_nopanel

    # Get default width/height ratios
    # ...width
    if width_ratios is None:
        width_ratios = np.ones(ncols)/ncols
    else:
        width_ratios = np.array(width_ratios)/sum(width_ratios)
    # ...height
    if height_ratios is None:
        height_ratios = np.ones(nrows)/nrows
    else:
        height_ratios = np.array(height_ratios)/sum(height_ratios)

    # Setup up figure
    fig = plt.figure(FigureClass=Figure, figsize=figsize)

    # Outer area; includes regions for special panels
    GS = mgridspec.GridSpec(
            nrows         = 1+int(bottomcolorbar or bottomlegend),
            ncols         = 1+int(rightcolorbar),
            left          = left,
            bottom        = bottom,
            right         = right, # unique spacing considerations
            top           = top, # so far no panels allowed here
            wspace        = wspace_outer,
            hspace        = hspace_outer,
            width_ratios  = width_ratios_outer,
            height_ratios = height_ratios_outer
            ) # set wspace/hspace to match the top/bottom spaces

    # Main plotting area
    gs = mgridspec.GridSpecFromSubplotSpec(
            nrows         = nrows,
            ncols         = ncols,
            subplot_spec  = GS[0,0],
            wspace        = wspace,
            hspace        = hspace,
            width_ratios  = width_ratios,
            height_ratios = height_ratios
            )

    #--------------------------------------------------------------------------
    # Manage shared axes/axes with spanning labels
    #--------------------------------------------------------------------------
    # Find shared axes; get sets with identical spans in x or y (ignore panels)
    # ...preliminary stuff
    axes_ids = [np.where(array==i) for i in np.unique(array) if i>0] # 0 stands for empty, -1 for colorbar, -2 for legend
        # note that these locations should be **sorted** by axes id
    yrange = np.array([[xy[0].min(), xy[0].max()+1] for xy in axes_ids]) # yrange is shared columns
    xrange = np.array([[xy[1].min(), xy[1].max()+1] for xy in axes_ids])
    xmin   = np.array([xy[0].min() for xy in axes_ids])
    xmax   = np.array([xy[0].max() for xy in axes_ids])
    ymin   = np.array([xy[1].min() for xy in axes_ids])
    ymax   = np.array([xy[1].max() for xy in axes_ids])
    num_axes = len(axes_ids)
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
    ax = num_axes*[None] # list of axes
    if sharex:
        for i in xgroups_base: # sorts locations implicitly
            ax[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                    **cartopy_kw)
    if sharey:
        for i in ygroups_base:
            if ax[i] is None:
                ax[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        **cartopy_kw)

    # Dependent axes
    for i in range(num_axes):
        # Detect if we want to share this axis with another
        sharex_ax, sharey_ax = None, None # by default, don't share with other axes objects
        if sharex:
            igroup = np.where([i in g for g in xgroups])[0] # np.where works on lists
            if igroup.size==1:
                sharex_ax = ax[xgroups_base[igroup[0]]]
                if sharex_ax is None:
                    raise ValueError('Something went wrong; shared x axes was not already drawn.')
            elif igroup.size>1:
                raise ValueError('Something went wrong; axis %d belongs to multiple groups.' % i)
        if sharey:
            igroup = np.where([i in g for g in ygroups])[0] # np.where works on lists
            if igroup.size==1:
                sharey_ax = ax[ygroups_base[igroup[0]]]
                if sharey_ax is None:
                    raise ValueError('Something went wrong; shared x axes was not already drawn.')
            elif igroup.size>1:
                raise ValueError('Something went wrong; axis %d belongs to multiple groups.' % i)

        # Draw axes, and add to list
        if ax[i] is not None:
            # ...special considerations for axes that might, for example, be an x base, but share a y-axis
            if ax[i] is not sharex_ax:
                ax[i]._sharex = sharex_ax
            else:
                sharex_ax = None
            if ax[i] is not sharey_ax:
                ax[i]._sharey = sharey_ax
            else:
                sharey_ax = None
        else:
            # ...virgin axes... these are not an x base or a y base
            ax[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])], sharex=sharex_ax, sharey=sharey_ax,
                    **cartopy_kw)

        # Hide tick labels (not default behavior for manual sharex, sharey use)
        if sharex_ax is not None:
            for t in ax[i].xaxis.get_ticklabels(): t.set_visible(False)
            ax[i].xaxis.label.set_visible(False)
        if sharey_ax is not None:
            for t in ax[i].yaxis.get_ticklabels(): t.set_visible(False)
            ax[i].yaxis.label.set_visible(False)

    # Spanning axes; allow xlabels/ylabels to span them
    if spanx and len(xgroups_span)>0:
        for g, b in zip(xgroups_span, xgroups_span_base):
            ax[b].xaxis.label.set_transform(mtransforms.blended_transform_factory(
                    fig.transFigure, mtransforms.IdentityTransform()
                    ))
                # specify x, y transform in Figure coordinates
            xmin = min(ax[i].get_position().xmin for i in g)
            xmax = max(ax[i].get_position().xmax for i in g)
                # get min/max positions, in figure coordinates, of spanning axes
            # print('Group:', g, 'Base:', b, 'Span:', xmin, xmax)
            ax[b].xaxis.label.set_position(((xmin+xmax)/2, 0))
                # this is the shared xlabel
            for i in g:
                if i!=b: ax[i].xaxis.label.set_visible(False)
    if spany and len(ygroups_span)>0:
        for g, b in zip(ygroups_span, ygroups_span_base):
            ax[b].yaxis.label.set_transform(mtransforms.blended_transform_factory(
                    mtransforms.IdentityTransform(), fig.transFigure # specify x, y transform
                    ))
            ymin = min(ax[i].get_position().ymin for i in g)
            ymax = max(ax[i].get_position().ymax for i in g)
            # print('Group:', g, 'Base:', b, 'Span:', ymin, ymax)
            ax[b].yaxis.label.set_position((0, (ymin+ymax)/2))
                # this is the shared ylabel
            for i in g:
                if i!=b: ax[i].yaxis.label.set_visible(False)
    # Save
    fig.axes_main = ax

    #--------------------------------------------------------------------------
    # Create panel axes
    #--------------------------------------------------------------------------
    # Colorbar panels
    if bottomcolorbar:
        C = mgridspec.GridSpecFromSubplotSpec(
                nrows        = 1,
                ncols        = 3,
                wspace       = 0,
                subplot_spec = GS[1,0],
                width_ratios = ((1-shrink)/2, shrink, (1-shrink)/2)
                )
        c = fig.add_subplot(C[0,1])
        c.bottomlegend, c.bottomcolorbar, c.rightcolorbar = False, True, False
        c.format = MethodType(_format, c) # MethodType approach
        fig.bottomcolorbar = c
    if rightcolorbar:
        C = mgridspec.GridSpecFromSubplotSpec(
                nrows         = 3,
                ncols         = 1,
                hspace        = 0,
                subplot_spec  = GS[0,1],
                height_ratios = ((1-shrink)/2, shrink, (1-shrink)/2)
                )
        c = fig.add_subplot(C[1,0])
        c.bottomlegend, c.bottomcolorbar, c.rightcolorbar = False, False, True
        c.format = MethodType(_format, c) # MethodType approach
        fig.rightcolorbar = c
    # Bottom legend panel
    if bottomlegend:
        l = fig.add_subplot(GS[1,0])
        l.bottomlegend, l.bottomcolorbar, l.rightcolorbar = True, False, False
        for s in l.spines.values():
            s.set_visible(False)
        l.xaxis.set_visible(False)
        l.yaxis.set_visible(False)
        l.patch.set_alpha(0)
        l.format = MethodType(_format, l) # MethodType approach
        fig.bottomlegend = l # had to set some properties first

    #--------------------------------------------------------------------------
    # Dynamically add properties and methods to Axes objects, and meethods
    #--------------------------------------------------------------------------
    # Create some attributes; analagous to default behavior where title 'exists' but is not visible
    # or is empty string, methods in _format will make them visible
    for i, a in enumerate(ax):
        # Identification
        a.bottomlegend, a.bottomcolorbar, a.rightcolorbar = False, False, False # identifiers
        # Helper properties (user can investigate theses)
        a.number = i # know axes number ahead of time
        a.width  = np.diff(a._position.intervalx)*width
        a.height = np.diff(a._position.intervaly)*height
        # Alternative title placement
        # Inside-center, on top/left-aligned, on top/right-aligned
        # a.titleinside = a.text(titlepad/a.width, 1-titlepad/a.height, '',
        #         zorder=-1, transform=a.transAxes, ha='left', va='top')
        a.titleinside = a.text(0.5, 1-titlepad/a.height, '',
                # bbox=dict(facecolor='w', edgecolor=None, linewidth=0, alpha=0.9, pad=2),
                zorder=-1, transform=a.transAxes, ha='center', va='top')
        a.titleleft = a.text(0+titlepad/a.width, 1-titlepad/a.height, '',
                zorder=-1, transform=a.transAxes, ha='left', va='top')
        a.titleright = a.text(1-titlepad/a.width, 1-titlepad/a.height, '',
                zorder=-1, transform=a.transAxes, ha='right', va='top')
        # ABC labeling
        a.abcinside = a.text(abcpad/a.width, 1-abcpad/a.height, ascii_uppercase[i], 
                # bbox=dict(facecolor='w', edgecolor=None, alpha=1),
                transform=a.transAxes, ha='left', va='top', visible=False)
        a.abc = a.text(0, 1, ascii_uppercase[i],
                transform=a.title._transform, ha='left', va='baseline', visible=False) # copies the a.title transform
        # Property cycling
        # To print current cycle, use list(next(ax._get_lines.prop_cycler)['color'] for i in range(10))
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        a.set_prop_cycle(cycler('color', [colors[i%10] for i in range(40)])
                + cycler('linestyle', [i for i in ('-','--','--','--') for n in range(10)])
                + cycler('dashes', [i for i in (tuple(), (1,1), (3,2), (6,3)) for n in range(10)]))

    # Method overrides; plotting (e.g. to fix white lines between polygons), and formatting
    for a in ax:
        # Formatting function
        a.format = MethodType(_format, a) # MethodType approach
        # Mapped overrides
        if projection is not None and package=='basemap':
            # Setup basemap
            a.m = mbasemap.Basemap(projection=projection, ax=a, **projection_kw)
            a.m.drawmapboundary() # **settings().line)
                # must be set before-hand, otherwise this wonky set_axes method automatically draws two mapboundary Patch objects
                # when you plot something, one for fill and the other for the edges; then, since the patch object in 
                # _mapboundarydrawn is only the fill-version, calling drawmapboundar() again will replace only ***that one***, 
                # but the original visible edges is still drawn; if you instead call drawmapboundary right away, calling it again will ***replace*** the object
            # First, set up new methods that fix white edges
            a.m._pcolormesh = a.m.pcolormesh
            a.m._pcolor     = a.m.pcolor
            a.m._contour    = a.m.contour
            a.m._contourf   = a.m.contourf # save old methods
            # Next, add fixed methods to basemap instances
            a.m.pcolormesh  = MethodType(_pcolormesh_basemap, a.m) # function calls the old methods
            a.m.pcolor      = MethodType(_pcolor_basemap, a.m)
            a.m.contour     = MethodType(_contour_basemap, a.m)
            a.m.contourf    = MethodType(_contourf_basemap, a.m)
                # these cannot be passed directly as elements of a, for some reason
                # can't override original name of axes object, because m.pcolor calls m.ax.pcolor
        # General polygon stuff
        else:
            # Same procecure as for basemap
            a.m = None
            a._pcolormesh = a.pcolormesh
            a._pcolor     = a.pcolor
            a._contourf   = a.contourf # copy the old methods
            a._contour   = a.contour # copy the old methods
            a.pcolor      = MethodType(_pcolor, a)
            a.pcolormesh  = MethodType(_pcolormesh, a)
            a.contourf    = MethodType(_contourf, a)
            a.contour    = MethodType(_contour, a)

    #--------------------------------------------------------------------------
    # Some final settings, and return result
    #--------------------------------------------------------------------------
    fig.patch.set_alpha(int(not transparent))
    for a in ax:
        a.patch.set_alpha(int(not transparent))
        a.format() # apply default properties
    if tight:
        GS.tight_layout(fig)
    print('Axes grid constructed.')
    if len(ax)==1:
        ax = ax[0]
    return fig, ax

#------------------------------------------------------------------------------
# Classes and Formatters
#------------------------------------------------------------------------------
class Normalize(mcolors.Normalize):
    """
    Pass as norm=<instance>, when declaring new pcolor or contourf objects.
    Base class for warping either or both sides of colormap.
    TODO: FIX THIS! IDEALLY, COLORS ON COLORBAR HAVE SAME GRADIATION BUT XTICKLABELS ARE WARPED; 
    INSTEAD, COLORBAR COLORS SEEM TO STRETCH AND TICKS ARE STILL LINEAR. OR... DO WE WANT THAT?
    """
    def __init__(self, exp=0, extend='neither', midpoint=None, vmin=None, vmax=None, clip=None):
        # User will use -10 to 10 scale; converted to value used in equation
        if abs(exp) > 10: raise ValueError('Warping scale must be between -10 and 10.')
        super().__init__(vmin, vmax, clip)
        self.midpoint = midpoint
        self.exp = exp
        self.extend = extend
        # mcolors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # Function
        def warp(x, exp, exp_max=4):
            # Returns indices stretched so neutral/low values are sampled more heavily
            # ...will artifically use exp to signify stretching away from neutral vals,
            # ...or compressing toward neutral vals
            if exp > 0:
                invert = True
            else:
                invert, exp = False, -exp
            exp = exp*(exp_max/10)
            # ...apply function; approaches x=1 as a-->Inf, x=x as a-->0
            if invert: x = 1-x
            value =  (x-1+(np.exp(x)-x)**exp)/(np.e-1)**exp
            if invert: value = 1-value # flip on y-axis
            # ...and return
            return value
        # Initial stuff
        if self.midpoint is None:
            midpoint = self.vmin
        else:
            midpoint = self.midpoint
        # Get middle point in 0-1 coords, and value
        midpoint_scaled = (midpoint - self.vmin)/(self.vmax - self.vmin)
        value_scaled    = (value - self.vmin)/(self.vmax - self.vmin)
        try:
            iter(value_scaled)
        except TypeError:
            value_scaled = np.arange(value_scaled)
        value_cmap = np.ma.empty(value_scaled.size)
        for i,v in enumerate(value_scaled):
            # ...get values, accounting for midpoints
            if v < 0: v = 0
            if v > 1: v = 1
            if v >= midpoint_scaled:
                block_width = 1 - midpoint_scaled
                value_cmap[i] = (midpoint_scaled + 
                        block_width*warp((v - midpoint_scaled)/block_width, self.exp)
                        )
            else:
                block_width = midpoint_scaled
                value_cmap[i] = (midpoint_scaled - 
                        block_width*warp((midpoint_scaled - v)/block_width, self.exp)
                        )
        if self.extend=='both' or self.extend=='max':
            value_cmap[value_cmap>1] = 1
        if self.extend=='both' or self.extend=='min':
            value_cmap[value_cmap<0] = 0
        return value_cmap

class Midpoint(mcolors.Normalize):
    """
    Pass as norm=<instance>, when declaring new pcolor or contourf objects.
    Creates new normalization of existing registered cmap by changing midpoint
    away from (vmin+vmax)/2.
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        # super(Midpoint,self).__init__(self, vmin, vmax, clip)
        self.midpoint = midpoint
        mcolors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        # print(value)
        # print(np.ma.masked_array(np.interp(value,x,y)))
        return np.ma.masked_array(np.interp(value, x, y))

def Formatter(sigfig=3):
    """
    Format as a number, with N sigfigs, and trimming trailing zeros.
    Recall, must pass function in terms of n (number) and loc.
    """
    # Format definition
    def f(n, loc, sigfig=sigfig):
        formatstr = '%%.%.ff' % (sigfig,)
        string    = formatstr % (n,)
        if '.' in string:
            string = string.rstrip('0').rstrip('.')
        return string
    # And create object
    return mticker.FuncFormatter(f)

def LatFormatter(sigfig=0, sine=False, NS=False):
    """
    Format latitude labels; can convert sine-lats back into lats (for
    areal weighting) and can apply N/S instead of postiive/negative.
    """
    def f(n, loc, sine=sine, NS=NS, sigfig=sigfig):
        # Convert from sine to latitude number
        if sine:
            n = np.arcsin(n)*180/np.pi
        # Suffix to apply
        if NS and n<0:
            s = 'S'
            n *= -1
        elif NS:
            s = 'N'
        else:
            s = ''
        formatstr = f'%.{sigfig:d}f'
        string = formatstr % (n,)
        return string
    # And create object
    return mticker.FuncFormatter(f)

def FracFormatter(fact=np.pi, symbol=r'\pi'):
    """
    Format as fractions, multiples of some value.
    """
    # Start with fraction definition
    def f(n, loc, fact=fact, symbol=symbol): # must accept location argument
        frac = Fraction(n/fact).limit_denominator()
        if n==0: # zero
            return '0'
        elif frac.denominator==1: # denominator is one
            if frac.numerator==1:
                return r'$%s$' % (symbol,)
            elif frac.numerator==-1:
                return r'${-}%s$' % (symbol,)
            else:
                return r'$%d%s$' % (frac.numerator, symbol)
        elif frac.numerator==1: # numerator is +/-1
            return r'$%s/%d$' % (symbol, frac.denominator)
        elif frac.numerator==-1:
            return r'${-}%s/%d$' % (symbol, frac.denominator)
        else: # otherwise
            return r'$%d%s/%d$' % (frac.numerator, symbol, frac.denominator)
    # And create FuncFormatter class
    return mticker.FuncFormatter(f)

