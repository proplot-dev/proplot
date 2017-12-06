#------------------------------------------------------------------------------
# Imports, all
#------------------------------------------------------------------------------
# from .basics import Dict, arange
import os # need for saving figures
from types import MethodType
from fractions import Fraction
# from string import ascii_uppercase
from string import ascii_lowercase
from cycler import cycler
from glob import glob
import matplotlib.figure as mfigure
# import matplotlib.axes as maxes
# import matplotlib.cm as mcm
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import matplotlib.text as mtext
import matplotlib.ticker as mticker
import matplotlib.gridspec as mgridspec
import matplotlib.container as mcontainer
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as mbasemap
# import string # for converting number to a/b/c
# import matplotlib.pyplot as plt # will attach a suitable backend, make available the figure/axes modules
import numpy as np # of course
import cartopy.crs as ccrs # crs stands for "coordinate reference system", leading c is "cartopy"
# import cartopy.mpl as cmpl
# import cartopy.feature as cfeature
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
# Colormap stuff and initial stuff
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
        if not _announcement: # only do this if register at least one new map
            _announcement = True
            print("Registered colormaps.")
# Create class for storing plot settings, and wrapper function that adds class
# instance as a global variable for this whole module; also declare default settings
class Settings():
    """
    This function initializes DEFAULT kwargs for various plot elements; user can
    UPDATE the defaults by passing dictionaries as kwargs when creating instance.
    """
    def __init__(self, color='k', linewidth=0.7, size1=8, size2=9, size3=9, **kwargs):
        # Axes look (updated)
        if any(type(kwarg) is not dict for kwarg in kwargs.values()):
            for kwarg in kwargs.values():
                print(f"{type(kwarg)}: {kwarg}.")
            raise ValueError("Arguments passed to settings must be dictionary objects.")
        self.tickminor = {'length':2, 'direction':'in', 'width':linewidth, 'color':color}
        self.tick      = {'length':4, 'direction':'in', 'width':linewidth, 'color':color}
        self.ctick     = {'length':4, 'direction':'out', 'width':linewidth, 'color':color}
            # lengths here are in points (i.e. 1/72in)
        self.spine     = {'linewidth':linewidth, 'color':color} # spine keyword
        self.box       = {'linewidth':linewidth, 'edgecolor':color}
        self.cgrid    = {'color':color, 'linewidth':linewidth}
            # for outlines of colorbar, legend, cartopy frame, and cartopy features
        # Text (updated, or set on creation)
        self.ticklabels = {'size':size1, 'weight':'normal', 'color':color}
        self.label      = {'size':size1, 'weight':'normal', 'color':color}
        self.title      = {'size':size2, 'weight':'normal', 'color':color}
        self.abc        = {'size':size3, 'weight':'bold', 'color':color}
        # Grid lines (set on creation; should make half-as-wide as axes)
        self.grid      = {'linestyle':'-', 'linewidth':linewidth/2, 'color':color, 'alpha':0.1}
        self.gridminor = {'linestyle':':', 'linewidth':linewidth/2, 'color':color, 'alpha':0.05}
        # Axes scale special properties (set on creation)
        self.xscale = {}
        self.yscale = {}
        # Geographic features look (set on creation)
        # self.continents  = {'color':'#CCCCCC'}
        # self.continents  = {'color':'moccasin'}
        self.continents  = {'color':'#eeeeee'} # best looking
        self.coastlines  = {'linewidth':linewidth, 'color':'#888888'}
        self.lonlatlines = {'linewidth':linewidth, 'color':color, 'alpha':0.2, 'dashes':(linewidth,linewidth)}
        self.boundary    = {'linewidth':linewidth, 'color':color, 'fill_color':'none'}
            # make sure to set fill_color to string 'none'; otherwise chooses some default,
            # accidentally made a bunch of plots blue this way
        # self.boundary    = self.spine # ignore that comment below; dumb
            # boundary should be thick... for some reason looks skinnier than parallels unless multiply by 2
            # even though (checking after the fact) reveals the linewidths are held
        # Auxilliary elements (set on creation)
        self.legend   = {'framealpha':0.6, 'fancybox':False, 'frameon':False, 'fontsize':size1}
        self.colorbar = {'extend':'both', 'spacing':'uniform'}
        # User-specified updated properties (provided as 'attribute_kw')
        # for attr in (attr for attr in dir(self) if not attr.startswith('__')):
        for attr, settings in self.__dict__.items():
            settings.update(kwargs.pop(attr+'_kw', {}))
            settings.update(kwargs.pop(attr, {}))
        if kwargs: # non-empty
            raise ValueError(f'Unknown kwarg-settings: {", ".join(kwargs.keys())}.')
    def __repr__(self):
        string = ['Current settings...']
        for n,v in self.__dict__.items():
            intro = f'{n}:'.ljust(15) # add whitespace to right
            info = ', '.join(f'{n}={v}' for n,v in v.items())
            string.append(intro + info)
        return '\n'.join(string)
# Now the wrapper function, and start
def setup(message="Plot settings updated.", **kwargs):
    if type(message) is not str: # because user called it manually, probably didn't declare
            # message, and passed something as an *arg without keyword
        raise ValueError("Plot settings can only be updated by keyword-dictionary pairs.")
    global settings # global w.r.t. this module
    settings = Settings(**kwargs)
    print(message)
setup("Default settings enabled.")

#------------------------------------------------------------------------------
# Colormap display
#------------------------------------------------------------------------------
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
    fig.savefig(f'{os.path.dirname(__file__)}/colormaps.pdf',
            bbox_inches='tight', format='pdf')
    return

#------------------------------------------------------------------------------
# Useful mainly when working with plots
#------------------------------------------------------------------------------
def arange(min_, *args):
    """
    Duplicate behavior of np.arange, except with inclusive endpoints; dtype is
    controlled very carefully, so should be 'most precise' among min/max/step args.
    Input...
        stop
        start, stop, [step]
        ...just like np.arange
    Output...
        the array sequence
    """
    # Optional arguments just like np.arange
    if len(args)==0:
        max_ = min_
        min_ = 0 # this re-assignes the NAME "min_" to 0
        step = 1
    elif len(args)==1:
        max_ = args[0]
        step = 1
    elif len(args)==2:
        max_ = args[0]
        step = args[1]
    else:
        raise ValueError('Function takes from one to three arguments.')
    # All input is integer? Get new "max"
    if min_//1==min_ and max_//1==max_ and step//1==step:
        min_, max_, step = np.int64(min_), np.int64(max_), np.int64(step)
        max_ += 1
    # Input is float or mixed; cast all to float64, then get new "max"
    else:
        min_, max_, step = np.float64(min_), np.float64(max_), np.float64(step)
        max_ += step/2
        # max_ = np.nextafter(max_, np.finfo(np.dtype(np.float64)).max)
            # gives the next FLOATING POINT, in direction of the second argument
            # ...forget this; round-off errors from continually adding step to min mess this up
    return np.arange(min_, max_, step)

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
        self.left = None
        self.bottom = None
        self.right = None
        self.top = None
        self.lwidth = None
        self.cspace = None
        self.bottomcolorbar = None # overall colorbar; replace with axis object
        self.bottomlegend = None # overall legend; replace with axis object
        self.rightcolorbar = None
        self.gridspec = None
        # General text
        # self.title = None
            # add this; functionality for title across multiple axes

    def save(self, filename, tight=True, pad=0.05): #, desktop=True):
        """
        Echo some helpful information before saving.
        Note that the gridspec object must be updated before figure is printed to screen
        in interactive environment... will fail to update after that. Seems to be glitch,
        should open thread on GitHub.
        """
        # Get bounding box, from axes
        # Unfortunately this does not find annotations either
        xs, ys = [], []
        for ax in self.axes: # try this method
            abbox = ax.get_tightbbox(self.canvas.get_renderer())
            abbox = mtransforms.TransformedBbox(abbox, mtransforms.Affine2D().scale(1./self.dpi))
            xs.append(abbox.intervalx)
            ys.append(abbox.intervaly)
        xa = [min(x[0] for x in xs), max(x[1] for x in xs)]
        ya = [min(y[0] for y in ys), max(y[1] for y in ys)]
        # Get bounding box
        obbox = self.bbox_inches # original bbox
        bbox = self.get_tightbbox(self.canvas.get_renderer())
        ox, oy, x, y = obbox.intervalx, obbox.intervaly, bbox.intervalx, bbox.intervaly
        x1, y1, x2, y2 = x[0], y[0], ox[1]-x[1], oy[1]-y[1] # deltas
        width, height = ox[1], oy[1] # desired save-width
        print(f'Extra border space: left {x1:.2f}, bottom {y1:.2f}, right {x2:.2f}, top {y2:.2f}.')
        # Echo some information
        # print('Axes-inferred intervals: ', xa, ya)
        # print('Figure intervals: ', x, y)
        left, top = f'left {self.left-x1:.2f}', f'top {self.top-y2:.2f}'
        if self.bottomlegend is not None:
            bottom = f'lwidth {self.lwidth-y1:.2f}'
        elif self.bottomcolorbar is not None:
            bottom = f'cspace {self.cspace-y1:.2f}'
        else:
            bottom = f'bottom {self.bottom-y1:.2f}'
        if self.rightcolorbar is not None:
            right = f'cspace {self.cspace-x2:.2f}'
        else:
            right = f'right {self.right-x2:.2f}'
        print(f'Try these params to avoid adjustment: {left}, {bottom}, {right}, {top}.')
        # Apply adjustment
        if tight:
            self.gridspec.left -= (x1-pad)/width
            self.gridspec.bottom -= (y1-pad)/height
            self.gridspec.right += (x2-pad)/width
            self.gridspec.top += (y2-pad)/height
            self.gridspec.update()
        # Finally, save
        if filename.split('.pdf')[-1]!='':
            filename = f'{filename}.pdf'
        print(f'Saving to {filename}.')
        self.savefig(filename, dpi=300, pad_inches=0, format='pdf') # specify DPI for embedded raster objects
        # print(f'Original sizing: left {self.left:.2f}, bottom {self.bottom:.2f}, right {self.right:.2f}, top {self.top:.2f}, '\
        #     + f'legend {self.lwidth:.2f}, colorbar {self.cwidth:.2f}.')

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
    # THIS FUNCTION IS BROKEN -- NEED TO FIX IT
    # Get graticule; want it to work with e.g. Gaussian-grid model output
    # dx, dy = x[1]-x[0], y[1]-y[0]
    # xb, yb = np.arange(x[0]-dx/2, x[-1]+dx, dx), np.arange(y[0]-dy/2, y[-1]+dy, dy)
    #     # can just use arange, because don't care about propagating float errors
    # return xb, yb
    edges = lambda vec: np.concatenate((vec[:1]-(vec[1]-vec[0])/2,
        (vec[1:]+vec[:-1])/2,
        vec[-1:]+(vec[-1]-vec[-2])/2))
    return edges(x), edges(y)

def _seam_fix(lon, data, lonmin=-180):
    # Raise errors
    if lon.max()>lon.min()+360:
        # print(lon.min(), lon.max())
        # print(lon.max()-lon.min()-360)
        raise ValueError('Longitudes must span 360 degrees at most.')
    if lon.min()<-360 or lon.max()>360:
        raise ValueError('Longitudes must fall in range [-360, 360].')
    if lonmin<-360 or lonmin>0:
        raise ValueError('Minimum longitude must fall in range [-360, 0].')

    # Establish 360-degree range
    # print(lon)
    lon -= 720
    while True:
        # print(lon, lonmin)
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
    # if lon[0]==lonmin:
    if lon[0]==lonmin and lon.size-1==data.shape[0]: # used longitude borders
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
# Twin axes override; add our special method to resulting generated axes
#------------------------------------------------------------------------------
def _twinx(self, **kwargs):
    ax = self._twinx(**kwargs)
    ax.format = MethodType(_format, ax)
    return ax
def _twiny(self, **kwargs):
    ax = self._twiny(**kwargs)
    ax.format = MethodType(_format, ax)
    return ax

#------------------------------------------------------------------------------
# Plot method overrides
#------------------------------------------------------------------------------
def _pcolorcheck(x, y, Z):
    x, y, Z = np.array(x), np.array(y), np.array(Z)
    if Z.shape[0]==x.size and Z.shape[1]==y.size:
        print("Guessing graticule edges from input coordinate centers.")
        x, y = _graticule_fix(x, y)
    elif Z.shape[0]!=x.size-1 or Z.shape[1]!=y.size-1:
        raise ValueError(f'X (size {x.size}) and Y (size {y.size}) must correspond to '
                f'nrows ({Z.shape[0]}) and ncolumns ({Z.shape[1]}) of Z, or its borders.')
    return x, y
def _contourcheck(x, y, Z):
    x, y, Z = np.array(x), np.array(y), np.array(Z)
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
    # print('initial:', lon, lat)
    # print('after shape check:',lon, lat)
    # print('after seam fix:', lon, lat)
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
def _cformat(self, mappable, cgrid=False, clocator=None, cformatter=None, clabel=None, old=False, **kwargs): #, settings=None):
    """
    Function for formatting colorbar-axes (axes that are "filled" by a colorbar).
    There is an INSANELY WEIRD problem with colorbars when simultaneously passing levels
    and norm object to a mappable fixed by passing vmin/vmax INSTEAD OF levels 
    (see: https://stackoverflow.com/q/40116968/4970632); problem is, often WANT levels instead
    of vmin/vmax to get this (https://stackoverflow.com/q/42723538/4970632) behavior.
    Solution is to make sure locators are in vmin/vmax range EXCLUSIVELY; cannot match/exceed values.
    TODO: Issue appears where the normalization vmin/vmax are outside of explicitly declared "levels" minima and maxima
    ...but that is probaby appropriate. If your levels are all within vmin/vmax, you will get discrete jumps outside of range
    ...and the triangles at ends of colorbars will be weird.
    """
    # Some settings
    if self.bottomcolorbar:
        orientation = 'horizontal'
    elif self.rightcolorbar:
        orientation = 'vertical'
    else:
        raise ValueError('Axes passed to colorbar-formatting function is not colorbar-dedicated.')
    csettings = settings.colorbar.copy()
    csettings.update(cax=self, use_gridspec=True, # use space afforded by entire axes
            orientation=orientation, drawedges=cgrid) # use ***ticklocation*** arg to set ticks on bottom/top, left/right
    # Update with user-kwargs
    csettings.update(**kwargs)
    # Draw colorbar
    if not old:
        # Simple method; just pass locator to ticks object
        # The other method can get mysterious errors where colorbar reverses directions for no reason, 
        # and the triangle extensions go white, and get some weird Warning; changing vmin/vmax can fix that
        if isinstance(clocator, mticker.Locator):
            print("WARNING: If this locator includes points matching or outside vmin/vmax, you may"\
                    "trigger a myserious error.")
        elif clocator is not None and not isinstance(clocator, mticker.AutoLocator):
            try: iter(clocator)
            except TypeError:
                clocator = _autolocate(mappable.vmin, mappable.vmax, clocator)
            # Option A: modify locations
            # lodiff, hidiff = clocator[1]-clocator[0], clocator[-1]-clocator[-2]
            # for i,c in enumerate(clocator.flat):
            #     if c<=mappable.vmin: # this one necessary
            #         clocator[i] = mappable.vmin + lodiff/10000
            #     if c>=mappable.vmax: # ...this one wasn't!!!
            #         clocator[i] = mappable.vmax + hidiff/10000
            # Option B: modify norm
            # clocator = [c for c in clocator if c>mappable.vmin]
            # Option C: modify vmin in original object; interestingly, changing mappable.vmin
            # does nothing, but changing the NORM vmin fixes the entire problem
            if mappable.norm is not None:
                # Adjust mappable parameters so the "actual" minimum is a bit below every clocator
                # First for my special BoundaryNorm object
                if hasattr(mappable.norm,'levels'):
                    mappable.norm.levels[0] -= np.diff(mappable.norm.levels[:2])[0]/10
                else:
                    mappable.norm.vmin -= (mappable.norm.vmax-mappable.norm.vmin)/1000
            clocator = mticker.FixedLocator(clocator)
            # print(mappable.norm.vmin, mappable.norm.vmax)
        else:
            raise ValueError("For safety, you must declare the clocator tick locator explicitly.")
        if not isinstance(cformatter, mticker.Formatter) and cformatter is not None:
            if isinstance(cformatter, str):
                cformatter = mticker.FormatStrFormatter(cformatter) # %-style, numbers
            else:
                cformatter = mticker.FixedFormatter(cformatter)
        cb = self.figure.colorbar(mappable, ticks=clocator, format=cformatter, **csettings)
    else:
        # First declare
        cb = self.figure.colorbar(mappable, **csettings)
        # cb = self.figure.colorbar(mappable, ticks=[-100,-5,5], **csettings)

        # Next, ticks formatters/locators (different from axis_setup version)
        # Account for weird error issue discussed in docstring
        if isinstance(clocator, mticker.Locator):
            print("WARNING: If this locator includes points matching or outside vmin/vmax, you may"\
                    "trigger a myserious error.")
            cb.locator = clocator
        elif clocator is not None:
            try: iter(clocator)
            except TypeError: clocator = _autolocate(cb.vmin, cb.vmax, clocator)
            lodiff, hidiff = clocator[1]-clocator[0], clocator[-1]-clocator[-2]
            for i,c in enumerate(clocator.flat):
                if c<=cb.vmin: # this one necessary
                    clocator[i] = cb.vmin + lodiff/10000
                # if c>=cb.vmax: # ...this one wasn't!!!
                #     clocator[i] = cb.vmax + hidiff/10000
            cb.locator = mticker.FixedLocator(clocator)
        if isinstance(cformatter, mticker.Formatter):
            cb.formatter = cformatter
        elif cformatter is not None:
            if isinstance(cformatter, str):
                cb.formatter = mticker.FormatStrFormatter(cformatter) # %-style, numbers
            else:
                cb.formatter = mticker.FixedFormatter(cformatter)
        # And update (use this instead of update_bruteforce)
        cb.update_ticks() # updates formatters/locators; colorbar uses weird combo
            # of its own patch objects, and underlying axis, so have to update manually

    # The ticks/ticklabels basic properties
    if cb.orientation=='horizontal':
        axis = cb.ax.xaxis
    else:
        axis = cb.ax.yaxis
    for t in axis.get_ticklabels(which='both'):
        t.update(settings.ticklabels)
    axis.set_tick_params(which='both', **settings.ctick)
        # properties are obscure; would have to use hidden methods to do this manually; 
        # using update method (inhereted from Artist) doesn't work
    # And label
    axis.label.update(settings.label)
    if clabel is not None:
        axis.label.set(text=clabel) # before, used set_label_text

    # Fix pesky white lines between levels + misalignment with border due to rasterized blocks
    cb.solids.set_rasterized(False)
    cb.solids.set_edgecolor('face')

    # Make edges/dividers consistent with axis edges
    if cb.dividers is not None:
        cb.dividers.update(settings.cgrid)
    cb.outline.update(settings.box)

    return cb

def _lformat(self, handles=None, **kwargs): #, settings=None): # can be updated
    """
    Function for formatting legend-axes (invisible axes with centered legends on them).
    Should update my legend function to CLIP the legend box when it goes outside axes area, so
    the legend-width and bottom/right widths can be chosen propertly/separately.
    """
    # Setup legend properties
    candidates = ['linewidth','color'] # candidates for modifying legend objects
    lsettings = settings.legend.copy()
    if self.bottomlegend: # is this axes assigned to hold the bottomlegend?
        if handles is None: # must specify
            raise ValueError('Must input list of handles.')
        lsettings.update( # need to input handles manually, because for separate axes
                bbox_transform=self.transAxes, # in case user passes bbox_to_anchor
                borderaxespad=0, loc='upper center', # this aligns top of legend box with top of axes
                framealpha=0) # make opaque; there will never be lines underneath
    elif handles is None:
        handles, _ = self.get_legend_handles_labels()
    lsettings.update(**kwargs)
    # Setup legend text properties
    tsettings = settings.ticklabels.copy()
    if 'fontsize' in lsettings:
        tsettings.pop('size')
    # Setup legend handle properties
    hsettings = {}
    for candidate in candidates:
        if candidate in kwargs:
            hsettings[candidate] = lsettings.pop(candidate)
    # Detect if user wants to specify rows manually
    try: iter(handles[0])
    except TypeError:
        multirow = False
        handles = [handles]
    else: # catch exception: only Container objects can be iterated
        if isinstance(handles[0], mcontainer.Container):
            multirow = False
            handles = [handles]
        else:
            multirow = True
    # Next draw legend over each row
    legends = []
    interval = 1/len(handles) # split up axes
    for h,hs in enumerate(handles):
        bbox = mtransforms.Bbox([[0,1-(h+1)*interval],[1,1-h*interval]])
        if 'ncol' in lsettings:
            lsettings.pop('ncol')
        if 'bbox_to_anchor' in lsettings:
            lsettings.pop('bbox_to_anchor')
        leg = self.legend(handles=hs, ncol=len(hs), bbox_to_anchor=bbox, **lsettings)
        leg.legendPatch.update(settings.box) # or get_frame()
        for obj in leg.legendHandles:
            obj.update(hsettings)
        for t in leg.texts:
            t.update(tsettings) # or get_texts()
        legends.append(leg)
    for l in legends[:-1]:
        self.add_artist(l) # because matplotlib deletes them...
    return legends

def _format(self, transparent=True, 
    coastlines=True, continents=False, # coastlines and continents
    latlabels=[0,0,0,0], lonlabels=[0,0,0,0], latlocator=None, lonlocator=None, # latlocator/lonlocator work just like xlocator/ylocator
    xgrid=None, ygrid=None, # gridline toggle
    xdates=False, ydates=False, # whether to format axis labels as long datetime strings
    xtickminor=None, ytickminor=None, xgridminor=None, ygridminor=None, # minor ticks/grids; if ticks off, grid will be off
    xspineloc=None, yspineloc=None, # deals with spine options
    xtickdir=None, ytickdir=None, # change ytick/xtick location; can be in, out, or inout (left-right-up-down depends on spine to which this applied)
    xlim=None, ylim=None, xscale=None, yscale=None, xreverse=False, yreverse=False, # special properties
    xlabel=None, ylabel=None, # axis labels
    title=None, titleleft=None, titleright=None, titleinside=None, abc=False, abcinside=False, # titles
    # TODO: add options for allowing UNLABELLED major ticklines -- maybe pass a special Formatter?
    xlocator=None, xminorlocator=None, ylocator=None, yminorlocator=None, # locators, or derivatives that are passed to locators
    xformatter=None, yformatter=None): # formatter
    # legend=False, handles=None, # legend options; if declared and have no handles with 'label', does nothing
    # settings=None): # remaining are passed to Settings()
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
    #--------------------------------------------------------------------------
    # Titles (should be done before projection stuff)
    #--------------------------------------------------------------------------
    # Title options
    padding = 0.1 # as fractino of axes
    self.title.update(settings.title)
    if title is not None:
        self.title.set_text(title)
    if titleleft is not None:
        t = self.text(0, 1, titleleft,
                transform=self.title._transform, ha='left', va='baseline')
        t.update(settings.title)
    if titleright is not None:
        t = self.text(1, 1, titleright,
                transform=self.title._transform, ha='right', va='baseline')
        t.update(settings.title)
    if titleinside is not None:
        t = self.text(0.5, 1-padding/self.height, titleinside, 
                transform=self.transAxes, ha='center', va='top')
        t.update(settings.title)
    if abcinside:
        if not hasattr(self, "number"):
            raise ValueError("Axes must have number property.")
        t = self.text(padding/self.width, 1-padding/self.height, ascii_lowercase[self.number], 
                transform=self.transAxes, ha='left', va='top')
        t.update(settings.abc)
    if abc:
        if not hasattr(self, "number"):
            raise ValueError("Axes must have number property.")
        t = self.text(0, 1, ascii_lowercase[self.number],
                transform=self.title._transform, ha='left', va='baseline')
        t.update(settings.abc)

    #--------------------------------------------------------------------------
    # Process projection/map axes
    #--------------------------------------------------------------------------
    # Basemap axes setup
    if hasattr(self,'m'):
        if self.m is not None:
            # Coastlines, parallels, meridians
            m = self.m
            if coastlines:
                p = m.drawcoastlines(**settings.coastlines)
            if continents:
                p = m.fillcontinents(**settings.continents)
            # Longitude/latitude lines
            # Make sure to turn off clipping by invisible axes boundary; otherwise
            # get these weird flat edges where map boundaries, parallel/meridian markers come up to the axes bbox
            if latlocator is not None:
                try: iter(latlocator)
                except TypeError:
                    latlocator = _autolocate(m.latmin+latlocator, m.latmax-latlocator, latlocator)
                p = m.drawparallels(latlocator, labels=latlabels) 
                for pi in p.values(): # returns dict, where each one is tuple
                    for _ in [i for j in pi for i in j]: # magic
                        if isinstance(_, mtext.Text):
                            _.update(settings.ticklabels)
                        else:
                            _.set_clip_on(False)
                            _.update(settings.lonlatlines)
                            # _.set_linestyle(linestyle)
                    # tried passing clip_on to the above, but it does nothing; must set
                    # for lines created after the fact
            if lonlocator is not None:
                try: iter(lonlocator)
                except TypeError:
                    lonlocator = _autolocate(m.lonmin+lonlocator, m.lonmax-lonlocator, lonlocator)
                p = m.drawmeridians(lonlocator, labels=lonlabels)
                for pi in p.values():
                    for _ in [i for j in pi for i in j]: # magic
                        if isinstance(_, mtext.Text):
                            _.update(settings.ticklabels)
                        else:
                            _.set_clip_on(False)
                            _.update(settings.lonlatlines)
                            # _.set_linestyle(linestyle)
            # Map boundary
            # p = m.drawmapboundary(**settings.boundary)
            p = m.drawmapboundary(**settings.boundary)
            p.set_facecolor(None)
            p.set_zorder(-1)
            p.set_clip_on(False)
            p.set_rasterized(False) # not sure about this... might be rasterized
            return # skip everything else
    # # Cartopy axes setup
    # if isinstance(self, cmpl.geoaxes.GeoAxes): # the main GeoAxes class; others like GeoAxesSubplot subclass this
    #     self.add_feature(cfeature.COASTLINE, **settings.box)
    #     plt.setp(self.outline_patch, **settings.box)
    #     return
    # # Legend
    # if legend:
    #     kwarg = ({} if handles is None else {'handles':handles})
    #     _lformat(self, **kwarg)

    #--------------------------------------------------------------------------
    # Process normal axes, various x/y settings individually
    #--------------------------------------------------------------------------
    # Axes scaling, limits, and reversal options (alternatively, supply
    # ...your own xlim/ylim that go from high to low)
    if xscale is not None: self.set_xscale(xscale, **settings.xscale)
    if yscale is not None: self.set_yscale(yscale, **settings.yscale)
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
            spine.update(settings.spine)
            # Standard
            if spineloc=='neither':
                spine.set_visible(False)
            elif spineloc=='both':
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
        majorSet = settings.tick if tickdir is None else dict(settings.tick, direction=tickdir)
        minorSet = settings.tickminor if tickdir is None else dict(settings.tickminor, direction=tickdir)
        # Apply
        axis.set_tick_params(which='major', **sidesOn, **gridOn, **majorSet)
        axis.set_tick_params(which='minor', **sidesOn, **gridminorOn, **minorSet) # have length
        for tick in axis.majorTicks:
            tick.gridline.update(settings.grid)
        for tick in axis.minorTicks:
            tick.gridline.update(settings.gridminor)
            if tickminor is not None:
                tick.set_visible(tickminor)

        # Label properties
        axis.label.update(settings.label)
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
        if tickformatter in ['none','None','NONE']:
            axis.set_major_formatter(mticker.NullFormatter())
        elif isinstance(tickformatter, mticker.Formatter): # use isinstance, because r"stuff" or f"stuff" **do not work**
            axis.set_major_formatter(tickformatter)
        elif tickformatter is not None:
            if isinstance(tickformatter, str):
                if dates: axis.set_major_formatter(mdates.DateFormatter(tickformatter)) # %-style, dates
                else: axis.set_major_formatter(mticker.FormatStrFormatter(tickformatter)) # %-style, numbers
            else:
                # axis.set_ticklabels(tickformatter) # ...FixedFormatter alone has issues
                axis.set_major_formatter(mticker.FixedFormatter(tickformatter)) # list of strings
        for t in axis.get_ticklabels():
            t.update(settings.ticklabels)
    return # we're done

def subplots(array=None, nrows=1, ncols=1,
        tight=False, # whether to set up tight bbox from gridspec object
        hatch=None, # hatching throughout background
        transparent=True, # whether figure should be transparent
        sharex=True, sharey=True, # for sharing x/y axis limits/scales/locators for axes with matching GridSpec extents, and making ticklabels/labels invisible
        spanx=True, spany=True, # custom setting; share axis labels for axes with same xmin/ymin extents
            # what about bottom vs. top xlabels, left vs. right ylabels?
        aspect=None, height=None, width=None, # for controlling aspect ratio; default is control for width
        hspace=0.4, wspace=0.2, hratios=None, wratios=None, # options for gridspec (instead of gridpsec_kw)
            # spacing betwen axes, in inches; generally horizontal-spacing small, but vertical-spacing big to allow for title
        left=0.8, bottom=0.5, right=0.1, top=0.3,
            # edges of plot for GridSpec object, in inches; left needs a lot of space
        bottomcolorbar=False, rightcolorbar=False, bottomlegend=False,
            # for legend/colorbar referencing multiple subplots
        bottomcolorbars=False, rightcolorbars=False, 
            # for colorbar labelling of EACH ROW/EACH COLUMN (have not implemented yet)
        lwidth=0.15, cwidth=0.2, cspace=0.8, cshrink=0.9,
            # spacing for colorbar text, colorbar axes width, legend axes with, and padding for interior ABC label
        # abcpad=0.1, titlepad=0.1,
        #     # will delete this stuff
        package='basemap', projection=None, **projectionkw): # for projections; can be 'basemap' or 'cartopy'
    """
    Special creation of subplots grids, allowing for arbitrarily overlapping 
    axes objects. Will return figure handle and axes object. Can also use exaclty as
    plt.subplot, with some added convenience features. Use panelx, panely as 
    templates for figure-wide legends/colorbars; if you actually want ax smallish panel
    to plot data, just use width/height ratios. don't bother with kwarg dictionaries because
    it is really, really customized; we do use them in formatting function though.
    
    TODO:
        * for rightcolorbar/bottomcolorbar, include kwarg option "row=<x>" and we can optionally
        create spare axes for a colorbar to reference EVERY SINGLE ROW/COLUMN
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

    # Patch, make sure width never None
    if width is None:
        width = 5
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
        bottom_extra_axes   = cwidth
        bottom_extra_hspace = cspace
    elif bottomlegend:
        bottom_extra_axes   = lwidth
        bottom_extra_hspace = 0 # don't need space
    else:
        bottom_extra_axes, bottom_extra_hspace = 0, 0
    bottom_extra = bottom_extra_axes + bottom_extra_hspace

    # And right spacings
    if rightcolorbar:
        right_extra_axes   = cwidth
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
    if len(projectionkw)>0 and projection is None:
        raise ValueError(f"Unknown arguments: {', '.join(projectionkw.keys())}.")
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
        cartopy_kw = {'projection': crs_dict[projection](**projectionkw)}
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
        # print(width, left, right, ncols, wspace, right_extra, ncols)
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
    # Properties for outer GridSpec object, need borders in fractional units
    ileft, ibottom, iright, itop = left, bottom, right, top
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
    # Properties for inner GridSpec object
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

    # Setup up figure and plot areas
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
    fig.left = ileft
    fig.bottom = ibottom
    fig.right = iright
    fig.top = itop
    fig.lwidth = lwidth
    fig.cspace = cspace
    fig.gridspec = GS # add to figure, for later reference
    # Main plotting area
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
                width_ratios = ((1-cshrink)/2, cshrink, (1-cshrink)/2)
                )
        c = fig.add_subplot(C[0,1])
        c.bottomlegend, c.bottomcolorbar, c.rightcolorbar = False, True, False
        c.format = MethodType(_cformat, c) # MethodType approach
        fig.bottomcolorbar = c
        # fig.bottomcolorbar = MethodType(_cformat, c)
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
        c.format = MethodType(_cformat, c) # MethodType approach
        fig.rightcolorbar = c
        # fig.rightcolorbar = MethodType(_cformat, c)
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
        l.format = MethodType(_lformat, l) # MethodType approach
        fig.bottomlegend = l # had to set some properties first
        # fig.bottomlegend = MethodType(_lformat, l)
        # fig.legend = l

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
        #----------------------------------------------------------------------
        # DELETE THIS STUFF; has proven useless
        # MAYBE JUST ADD **METHODS** INSTEAD... or move this stuff to the FORMAT 
        # function because never really want to pad by inches, usually will want to pad by 
        # FRACTION OF AXES SIZE... and is AWKWARD to have to add this tertiary stuff right away
        # Alternative title placement
        # Inside-center, on top/left-aligned, on top/right-aligned
        # a.titleinside = a.text(0.5, 1-titlepad/a.height, '',
        #         # bbox=dict(facecolor='w', edgecolor=None, linewidth=0, alpha=0.9, pad=2),
        #         zorder=-1, transform=a.transAxes, ha='center', va='top')
        # a.titleleft = a.text(0+titlepad/a.width, 1-titlepad/a.height, '',
        #         zorder=-1, transform=a.transAxes, ha='left', va='top')
        # a.titleright = a.text(1-titlepad/a.width, 1-titlepad/a.height, '',
        #         zorder=-1, transform=a.transAxes, ha='right', va='top')
        # # ABC labeling
        # a.abcinside = a.text(abcpad/a.width, 1-abcpad/a.height, ascii_uppercase[i],
        #         # bbox=dict(facecolor='w', edgecolor=None, alpha=1),
        #         transform=a.transAxes, ha='left', va='top', visible=False)
        # a.abc = a.text(0, 1, ascii_uppercase[i],
        #         transform=a.title._transform, ha='left', va='baseline', visible=False) # copies the a.title transform
        #----------------------------------------------------------------------
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
        # Legend-drawing on any axes
        a.legend = MethodType(_lformat, a) # MethodType approach
        # Mapped overrides
        if projection is not None and package=='basemap':
            # Setup basemap
            a.m = mbasemap.Basemap(projection=projection, ax=a, **projectionkw)
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
        # Simple changes; cover white lines and change input shape
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
        # Change twinx/twiny properties
        a._twinx = a.twinx
        a.twinx = MethodType(_twinx, a)
        a._twiny = a.twiny
        a.twiny = MethodType(_twiny, a)

    #--------------------------------------------------------------------------
    # Some final settings, and return result
    #--------------------------------------------------------------------------
    fig.patch.set_alpha(int(not transparent))
    for a in ax:
        a.patch.set_alpha(int(not transparent))
        # a.format() # apply default properties
        # or not... messes up projections
        if hatch is not None:
            a.fill_between([0,1], 0, 1, hatch=hatch, facecolor='none', edgecolor='k', 
                    transform=a.transAxes) # fill hatching in background
    if len(ax)==1:
        ax = ax[0]
    print('Axes grid constructed.')
    return fig, ax

#------------------------------------------------------------------------------
# Classes and Formatters
#------------------------------------------------------------------------------
class WarpNorm(mcolors.Normalize):
    """
    Pass as norm=<instance>, when declaring new pcolor or contourf objects.
    Base class for warping either or both sides of colormap.
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

class MidpointNorm(mcolors.Normalize):
    """
    Pass as norm=<instance>, when declaring new pcolor or contourf objects.
    Creates new normalization of existing registered cmap by changing midpoint
    away from (vmin+vmax)/2.
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        # Default None values so can declare by name in any order
        # super(Midpoint,self).__init__(self, vmin, vmax, clip)
        self.midpoint = midpoint
        if any(x is None for x in [vmin,vmax,midpoint]):
            raise ValueError("Must declare vmin, vmax, and midpoint explicitly.")
        mcolors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # How to map data values in range (vmin,vmax) to color indices in colormap
        # if self.midpoint<self.vmin or self.midpoint>self.vmax:
        #     raise ValueError("Midpoint {self.midpoint} is not between"\
        #             "vmin {self.vmin} and vmax {self.vmax}.")
        # print(x, np.interp(value, x, y))
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

class BoundaryNorm(mcolors.Normalize):
    """
    Like the default BoundaryNorm, except instead of declaring level integers
    from the exact RGB file, this interpolates from between 0 and 1 for each level.
    Example: Your levels edges are weirdly spaced [-1000, 100, 0, 100, 1000] or
    even [0, 10, 12, 20, 22], but center "colors" are always at colormap
    coordinates [.2, .4, .6, .8] no matter the spacing; levels just must be monotonic.
    """
    def __init__(self, levels, midpoint=None, clip=False):
        # Very simple
        # super(Midpoint,self).__init__(self, vmin, vmax, clip)
        try: iter(levels)
        except TypeError:
            raise ValueError("Must call BoundaryNorm with your boundary vaues.")
        self.midpoint = midpoint
        self.levels = np.array(levels)
        mcolors.Normalize.__init__(self, min(levels), max(levels), clip)

    def __call__(self, value, clip=None):
        # TOTO: Add optional midpoint; this class will probably end up being one of
        # my most used if so; midpoint would just ensure <value> corresponds to 0.5 in cmap
        # Some checks (maybe not necessary)
        try: iter(value)
        except TypeError:
            value = np.array([value])
        if value.ndim>1:
            raise ValueError("Array is multi-dimensional... not sure what to do.")
        # Map data values in range (vmin,vmax) to color indices in colormap
        nvalues = np.empty(value.shape)
        for i,v in enumerate(value.flat):
            if np.isnan(v):
                continue
            locs = np.where(v>=self.levels)[0]
            if locs.size==0:
                nvalues[i] = 0
            elif locs.size==self.levels.size:
                nvalues[i] = 1
            else:
                interpolee = self.levels[[locs[-1],locs[-1]+1]] # the boundary level values
                interpolant = np.array([locs[-1],locs[-1]+1])/(self.levels.size-1) # the boundary level integers
                    # so if 11 levels and between ids 9 and 10, interpolant is between .9 and 1
                nvalues[i] = np.interp(v, interpolee, interpolant)
                # nvalues[i] = max(0, nvalues[i]-.5/(self.levels.size-1))
                # print(self.vmin, self.vmax, min(self.levels), max(self.levels))
        return np.ma.masked_array(nvalues, np.isnan(value))

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

