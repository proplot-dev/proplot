#!/usr/bin/env python3
#------------------------------------------------------------------------------
# Imports, all
# Note that because of how matplotlib imports work (a module is only imported
# for a first time; import statements elsewhere in script during process just point
# to the already loaded module), settings rcParams here will change properties
# everywhere else until you start new instance
#------------------------------------------------------------------------------
# from .basics import Dict, arange
import os # need for saving figures
import copy # for copying Settings
from types import MethodType
from fractions import Fraction
# from string import ascii_uppercase
from string import ascii_lowercase
from cycler import cycler
from glob import glob
import matplotlib as mpl
import matplotlib.figure as mfigure
# import matplotlib.axes as maxes
# import matplotlib.cm as mcm
import matplotlib.patheffects as mpatheffects
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import matplotlib.text as mtext
import matplotlib.font_manager as mfonts
import matplotlib.ticker as mticker
import matplotlib.gridspec as mgridspec
import matplotlib.container as mcontainer
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as mbasemap
from matplotlib import matplotlib_fname
# import string # for converting number to a/b/c
# import matplotlib.pyplot as plt # will attach a suitable backend, make available the figure/axes modules
import numpy as np # of course
import cartopy.crs as ccrs # crs stands for "coordinate reference system", leading c is "cartopy"
import cartopy.feature as cfeature
import cartopy.mpl.geoaxes as cgeoaxes
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
# Initialization; stuff called on import
# Adds colormap names and lists available font names
#------------------------------------------------------------------------------
# List the system font names
# See: https://olgabotvinnik.com/blog/2012-11-15-how-to-set-helvetica-as-the-default-sans-serif-font-in/
fonts = [font.split('/')[-1].split('.')[0] for font in # system fonts
            mfonts.findSystemFonts(fontpaths=None, fontext='ttf')] + \
        [os.path.basename(font.rstrip('.ttf')) for font in # hack-installed fonts
            glob(f"{matplotlib_fname().rstrip('matplotlibrc')}/fonts/ttf/*.ttf")]
fonts = sorted(set(fonts))
# Alternate version, works on Linux but not Mac, because can't find mac system fonts?
# fonts = mfonts.get_fontconfig_fonts()
# fonts = [mfonts.FontProperties(fname=fname).get_name() for fname in flist]
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
        # Register colormap and a reversed version
        plt.register_cmap(name=_name, cmap=mcolors.LinearSegmentedColormap.from_list(name=_name, colors=_cmap, N=256))
        plt.register_cmap(name=_name+'_r', cmap=mcolors.LinearSegmentedColormap.from_list(name=_name+'_r', colors=_cmap[::-1], N=256))
        if not _announcement: # only do this if register at least one new map
            _announcement = True
            print("Registered colormaps.")
# For retrieving colormap
def cmapcolors(name, N):
    """
    Just spits out colors to be used for e.g. consecutive lines. First argument
    is the colormap name, second is the number of colors needed.
    """
    cmap = plt.cm.get_cmap(name) # the cmap object itself
    return [cmap(i/(N-1)) for i in range(N)] # from smallest to biggest numbers
# Used to create this neat class but it was a stupid idea. Work with existing
# API, not against it. Just create a function that sets rcParam defaults with
# optional override. Then, better to make the format function add actual information
# to the plot and do nothing to change its style/color/etc.
def rc(name, *args, silent=False, **kwargs):
    for arg in args:
        try:
            kwargs = {**kwargs, **arg}
        except TypeError:
            raise ValueError("Extra arguments must be dictionaries.")
    if kwargs: # if non-empty
        for key,value in kwargs.items():
            if f'{name}.{key}' in mpl.rcParams:
                mpl.rcParams[f'{name}.{key}'] = value
            else:
                # if f'{name}.{key}' not in mpl.rcExtras and not silent:
                #     print(f"Adding {name}.{key} to rcExtras.")
                if not hasattr(mpl, 'rcExtras'):
                    mpl.rcExtras = {}
                mpl.rcExtras[f'{name}.{key}'] = value
        return None
    else: # if empty
        dictionary = {}
        for key,value in {**mpl.rcParams, **mpl.rcExtras}.items():
            if key.split('.')[0]==name:
                dictionary[key.split('.')[1]] = value
        if not dictionary:
            raise ValueError(f"Could not find settings for {name}.")
        return dictionary
def setup(everything=True, **kwargs): # fontname is matplotlib default
    """
    See this page: https://matplotlib.org/users/customizing.html
    Quick list of rcParam categories:
        "lines", "patch", "hatch", "legend"
        "font", "text", "mathtext"
            includes usetex for making all matplotlib fonts latex
            and includes settings options for mathmode math, e.g. sf=sans
        "axes", "figure"
            includes tons of options
        "date", "xtick", "ytick", "grid"
            includes autoformatter
        "contour", "image"
            include lut for size of colormap table
        "boxplot", "errorbar", "hist", "scatter"
        "path", "savefig", "ps", "tk", "pdf", "svg"
        "debug", "keymap", "examples"
        "animation"
    * Problem is the settings in rcParams are scattershot and missing some important
      ones we want, but should use it when we can; options won't be rescinded in future versions
    * Currently features that can't be rcParam'ed are *gridminor* properties, the
      *cgrid* idea (lines between colorbar colors), the continent/coastlines/lonlatlines
      settings, and some legend settings.
    * One idea might be to **derive custom settings from existing matplotlib
      rcParam settings**.
    * Use this function to temporarily change settings, but only the specified settings.
    """
    # Create dictionary of settings; revert to default if user requested
    d = {'color':'k', 'linewidth':0.7, 'ssize':8, 'bsize':9, 'fontname':'DejaVu Sans'}
    for name in d.keys():
        d[name] = kwargs.pop(name,None) or d[name]
    if everything:
        # Make a cycler
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
            '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        propcycle = cycler('color', [colors[i%10] for i in range(40)]) \
                + cycler('linestyle', [i for i in ('-', '--', ':', '-.') for n in range(10)])
            # + cycler('dashes', [i for i in ((), (1,1), (3,2), (6,3)) for n in range(10)])
        # Reset the rcParams
        mpl.rcdefaults() # reset defaults
        mpl.rcExtras = {}
        # Included settings
        rc('axes', xmargin=0, ymargin=0, labelsize=d['ssize'], titlesize=d['bsize'],
                edgecolor=d['color'], labelcolor=d['color'],
                grid=True, linewidth=d['linewidth'],
                labelpad=3, prop_cycle=propcycle)
        rc('font', {'size':d['ssize'], 'family':'sans-serif', 'sans-serif':d['fontname']}) # when user calls .text()
        rc('text', color=d['color']) # when user calls .text()
        rc('grid', linewidth=d['linewidth']/2, alpha=0.1, color=d['color'])
        for xy in 'xy':
            rc(f'{xy}tick', labelsize=d['ssize'], color=d['color'], direction='out')
            rc(f'{xy}tick.major', pad=2, width=d['linewidth'], size=4, # tick length is size
                **{'x':{'bottom':True,'top':False},'y':{'left':True,'right':False}}[xy])
            rc(f'{xy}tick.minor', pad=2, width=d['linewidth'], visible=True, size=2,
                **{'x':{'bottom':True,'top':False},'y':{'left':True,'right':False}}[xy])
        rc('legend', framealpha=0.6, fancybox=False, frameon=False,
            labelspacing=0.5, columnspacing=1, handletextpad=0.5, borderpad=0.25)
        rc('patch', linewidth=d['linewidth'],   edgecolor=d['color']) # e.g. bars?
        rc('hatch', linewidth=d['linewidth'],   color=d['color'])
        rc('lines', linewidth=d['linewidth']*2, markersize=d['linewidth']*4, markeredgewidth=0)
        rc('markers', fillstyle='full')
        rc('scatter', marker='o')
        # Custom settings
        # Should begin to delete these and control things with rcParams whenever possible
        # Remember original motivation was realization that some rcParams can't be changes
        # by passing a kwarg (don't remember any examples)
        rc('abc', size=d['bsize'], weight='bold', color=d['color'])
        rc('title', size=d['bsize'], weight='normal', color=d['color'], fontname=d['fontname'])
        rc('label', size=d['ssize'], weight='normal', color=d['color'], fontname=d['fontname'])
        rc('ticklabels', size=d['ssize'], weight='normal', color=d['color'], fontname=d['fontname'])
        rc('gridminor', linestyle=':', linewidth=d['linewidth']/2, color=d['color'], alpha=0.05)
        rc('cgrid', color=d['color'], linewidth=d['linewidth'])
        # rc('xscale')
        # rc('yscale')
        rc('continents', color=d['color'])
        rc('tickminor', length=2, width=d['linewidth'], color=d['color'])
        rc('tick', length=4, width=d['linewidth'], color=d['color'])
        rc('ctickminor', length=2, width=d['linewidth'], color=d['color'])
        rc('ctick', length=4, width=d['linewidth'], color=d['color'])
        # self.continents  = {'color':'moccasin'}
        rc('coastlines', linewidth=d['linewidth'], color=d['color'])
        rc('lonlatlines', linewidth=d['linewidth'], linestyle=':', color=d['color'], alpha=0.2)
        rc('spine', color=d['color'], linewidth=d['linewidth'])
        rc('outline', edgecolor=d['color'], linewidth=d['linewidth'])
    # Overrides
    for key,value in kwargs.items():
        if not isinstance(value, dict):
            raise ValueError("Can only pass dictionaries.")
        rc(key, **value) # update dictionary
# Now call the function to configure params
setup(True)

#------------------------------------------------------------------------------
# Colormap display
#------------------------------------------------------------------------------
def cmapshow(N=11, ignore=['Qualitative','Miscellaneous','Sequential Alt']):
    """
    Plot all current colormaps, along with their catgories.
    This example comes from the Cookbook on www.scipy.org. According to the
    history, Andrew Straw did the conversion from an old page, but it is
    unclear who the original author is.
    See: http://matplotlib.org/examples/color/colormaps_reference.html
    """
    # Have colormaps separated into categories:
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
            # Get object
            cmap = plt.get_cmap(m, lut=N)
            # Draw, and make invisible
            ax = plt.subplot(nmaps,1,i+ntitles+nplots)
            for s in ax.spines.values():
                s.set_visible(False)
            ax.patch.set_alpha(0)
            ax.set_xticks([])
            ax.set_yticks([])
            # Draw colormap
            ax.imshow(a, aspect='auto', cmap=cmap, origin='lower')
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
    return fig

#------------------------------------------------------------------------------
# Function used often in plotting context but kind of unrelated
#------------------------------------------------------------------------------
def arange(min_, *args):
    """
    Duplicate behavior of np.arange, except with inclusive endpoints; dtype is
    controlled very carefully, so should be 'most precise' among min/max/step args.
    Input...
        stop
        start, stop, [step]
        Just like np.arange!
    Output...
        The array sequence.
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
        super().__init__(*args, **kwargs) # python 3 only

    def save(self, filename, squeeze=None, tight=None, pad=0.05): #, desktop=True):
        """
        Echo some helpful information before saving.
        Note that the gridspec object must be updated before figure is printed to screen
        in interactive environment... will fail to update after that. Seems to be glitch,
        should open thread on GitHub.
        """
        # Parse input
        if squeeze is None and tight is None:
            tight = True
        if squeeze is None and tight is not None:
            tight = squeeze
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
def _seam_fix(basemap, lon, lat, data, globe=True):
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
        data[(lon[1:]<lonmin) | (lon[:-1]>lonmax),...] = np.nan
    elif lon.size==data.shape[0]: # test the centers
        # this just tests centers and pads by one for safety
        # remember that a SLICE with no valid range just returns empty array
        dwhere = np.where((lon<lonmin) | (lon>lonmax))[0]
        data[dwhere[1:-1],...] = np.nan
    # 5) Fix holes over poles by interpolating there (equivalent to
    # simple mean of highest/lowest latitude points)
    if lon.size==data.shape[0]: # have centers, not grid cell edges
        dataS = np.repeat(data[:,0].mean(), data.shape[0])[:,None]
        dataN = np.repeat(data[:,-1].mean(), data.shape[0])[:,None]
        lat = np.concatenate(([-90], lat, [90]))
        data = np.concatenate((dataS, data, dataN), axis=1)
    # 6) Fix seams at map boundary
    # print(lon.shape, lat.shape, data.shape)
    if lon[0]==lonmin and lon.size-1==data.shape[0]: # borders fit perfectly
        # print('No seam fix necessary. Size augmented by 0.')
        pass # do nothing
    elif lon.size-1==data.shape[0]: # no interpolation necessary; just make a new grid cell
        # print('Fixing seam using lon borders. Size augmented by 1.')
        lon = np.append(lonmin, lon) # append way easier than concatenate
        lon[-1] = lonmin+360 # we've added a new tiny cell to the end
        data = np.concatenate((data[-1:,...], data), axis=0) # don't use pad, because messes up masked arrays
        # pad_tuple = ((1,0),) + ((0,0),)*(data.ndim-1)
        # data = np.pad(data, pad_tuple, 'wrap') # pad before
    elif lon.size==data.shape[0]: # linearly interpolate to the edges
        # print('Fixing seam using lon centers/interpolation. Size augmented by 2.')
        x = np.array([lon[-1], lon[0]+360]) # x
        y = np.concatenate((data[-1:,...], data[:1,...]), axis=0)
        xq = lonmin+360
        yq = (y[:1,...]*(x[1]-xq) + y[1:,...]*(xq-x[0]))/(x[1]-x[0]) # simple linear interp formula
        data = np.concatenate((yq, data, yq), axis=0)
        lon = np.append(np.append(lonmin, lon), lonmin+360)
    else:
        raise ValueError('Longitude length does not match data length along dimension 0.')
    return lon, lat, data
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
def _contourflevels(kwargs):
    # Processes keyword-arguments to allow levels for pcolor/pcolormesh
    # Using "lut" argument prevents colors in the extention-triangles/rectangles
    # from being different from the deepest colors with in the range of min(levels) max(levels)
    if 'levels' in kwargs:
        if 'cmap' not in kwargs: kwargs['cmap'] = 'viridis'
        kwargs['cmap'] = plt.cm.get_cmap(kwargs['cmap'], lut=len(kwargs['levels'])-1)
        kwargs['norm'] = BoundaryNorm(kwargs['levels']) # what is wanted 99% of time
        # kwargs['norm'] = BoundaryNorm(levels, ncolors=kwargs['cmap'].N, clip=True)
    return kwargs
def _pcolorlevels(kwargs):
    # Processes keyword-arguments to allow levels for pcolor/pcolormesh
    if 'levels' in kwargs:
        if 'cmap' not in kwargs: kwargs['cmap'] = 'viridis'
        levels = kwargs.pop('levels')
        kwargs['cmap'] = plt.cm.get_cmap(kwargs['cmap'], lut=len(levels)-1)
        kwargs['norm'] = BoundaryNorm(levels)
        # kwargs['norm'] = BoundaryNorm(levels, ncolors=kwargs['cmap'].N, clip=True)
    return kwargs
def _pcolorcheck(x, y, Z):
    # Checks that sizes match up, checks whether graticule was input
    x, y, Z = np.array(x), np.array(y), np.array(Z)
    if Z.shape[0]==x.size and Z.shape[1]==y.size:
        # print("Guessing graticule edges from input coordinate centers.")
        x, y = _graticule_fix(x, y)
    elif Z.shape[0]!=x.size-1 or Z.shape[1]!=y.size-1:
        raise ValueError(f'X ({"x".join(str(i) for i in x.shape)}) '
                f'and Y ({"x".join(str(i) for i in y.shape)}) must correspond to '
                f'nrows ({Z.shape[0]}) and ncolumns ({Z.shape[1]}) of Z, or its borders.')
    return x, y
def _contourcheck(x, y, Z):
    # Checks whether sizes match up, checks whether graticule was input
    x, y, Z = np.array(x), np.array(y), np.array(Z)
    if Z.shape[0]==x.size-1 and Z.shape[1]==y.size-1:
        x, y = (x[1:]+x[:-1])/2, (y[1:]+y[:-1])/2
    elif Z.shape[0]!=x.size or Z.shape[1]!=y.size:
        raise ValueError(f'X ({"x".join(str(i) for i in x.shape)}) '
                f'and Y ({"x".join(str(i) for i in y.shape)}) must correspond to '
                f'nrows ({Z.shape[0]}) and ncolumns ({Z.shape[1]}) of Z, or its borders.')
    return x, y
# The linecolorbar method is now self-contained in _cformat method
# def _linecolorbar(self, handles, values):
#     # Creates legend from a bunch of lines
#     if len(handles)!=len(values):
#         raise ValueError("Number of handles should equal number of levels.")
#     colors = [h.get_color() for h in handles]
#     colormap = mcolors.LinearSegmentedColormap.from_list('tmp', colors)
#     mappable = plt.contourf([[0,0],[0,0]], levels=values, cmap=colormap)
#     return mappable

#-------------------------------------------------------------------------------
# Cartesian plot overrides overrides
#-------------------------------------------------------------------------------
# TODO: Fix for pcolormesh/pcolor with set of levels
def _contour(self, x, y, Z, **kwargs):
    x, y = _contourcheck(x, y, Z)
    c = self._contour(x, y, Z.T, **kwargs)
    return c
def _contourf(self, x, y, Z, **kwargs):
    kwargs = _contourflevels(kwargs)
    x, y = _contourcheck(x, y, Z)
    c = self._contourf(x, y, Z.T, **kwargs)
    return _contour_fix(c)
def _pcolor(self, x, y, Z, **kwargs):
    kwargs = _pcolorlevels(kwargs)
    x, y = _pcolorcheck(x, y, Z)
    if 'extend' in kwargs:
        extend = kwargs.pop('extend')
    else:
        extend = None
    p = self._pcolor(x, y, Z.T, **kwargs)
    if extend is not None:
        p.extend = extend # add attribute to be used in colorbar creation
    return _pcolor_fix(p)
def _pcolormesh(self, x, y, Z, **kwargs):
    kwargs = _pcolorlevels(kwargs)
    x, y = _pcolorcheck(x, y, Z)
    if 'extend' in kwargs:
        extend = kwargs.pop('extend')
    else:
        extend = None
    p = self._pcolormesh(x, y, Z.T, **kwargs)
    if extend is not None:
        p.extend = extend # add attribute to be used in colorbar creation
    return _pcolor_fix(p)

#-------------------------------------------------------------------------------
# Basemap overrides
# Will convert coordinates to longitude latitude
#-------------------------------------------------------------------------------
# Dummy ones that just set latlon True by default
def _plot_basemap(self, *args, **kwargs):
    return self._plot(*args, **kwargs, latlon=True)
def _scatter_basemap(self, *args, **kwargs):
    return self._scatter(*args, **kwargs, latlon=True)
# Assumes regularly spaced data
def _contour_basemap(self, lon, lat, Z, **kwargs):
    lon, lat = _contourcheck(lon, lat, Z)
    lon, lat, Z = _seam_fix(self, lon, lat, Z)
    X, Y = _basemap_fix(self, lon, lat)
    c = self._contour(X, Y, Z.T, **kwargs)
    return c
def _contourf_basemap(self, lon, lat, Z, **kwargs):
    kwargs = _contourflevels(kwargs)
    lon, lat = _contourcheck(lon, lat, Z)
    lon, lat, Z = _seam_fix(self, lon, lat, Z)
    X, Y = _basemap_fix(self, lon, lat)
    c = self._contourf(X, Y, Z.T, **kwargs)
    return _contour_fix(c)
def _pcolor_basemap(self, lon, lat, Z, **kwargs):
    # lon, lat = _graticule_fix(lon, lat)
    kwargs = _pcolorlevels(kwargs)
    lon, lat = _pcolorcheck(lon, lat, Z)
    lon, lat, Z = _seam_fix(self, lon, lat, Z)
    X, Y = _basemap_fix(self, lon, lat)
    if 'extend' in kwargs:
        extend = kwargs.pop('extend')
    else:
        extend = None
    p = self._pcolor(X, Y, Z.T, **kwargs)
    if extend is not None:
        p.extend = extend # add attribute to be used in colorbar creation
    return _pcolor_fix(p)
def _pcolormesh_basemap(self, lon, lat, Z, **kwargs):
    # lon, lat = _graticule_fix(lon, lat)
    kwargs = _pcolorlevels(kwargs)
    lon, lat = _pcolorcheck(lon, lat, Z)
    lon, lat, Z = _seam_fix(self, lon, lat, Z)
    X, Y = _basemap_fix(self, lon, lat)
    if 'extend' in kwargs:
        extend = kwargs.pop('extend')
    else:
        extend = None
    p = self._pcolormesh(X, Y, Z.T, **kwargs)
    if extend is not None:
        p.extend = extend # add attribute to be used in colorbar creation
    return _pcolor_fix(p)

#------------------------------------------------------------------------------
# Formatting functions, assigned using MethodType onto various Axes
# instances and GirdSpec instances
#------------------------------------------------------------------------------
def _text(self, x, y, text, transform='axes', fancy=False, black=True, edgewidth=2, **kwarg):
    """
    Wrapper around original text method.
    If requested, can return specially created black-on-white text.
    So far no other features implemented.
    """
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
            'path_effects': [mpatheffects.Stroke(linewidth=2, foreground=bcolor), mpatheffects.Normal()]})
        # t.update({'size':11, 'zorder':1e10,
        #     'path_effects':[mpatheffects.PathPatchEffect(edgecolor=bcolor,linewidth=.6,facecolor=fcolor)]})
    return t

def _autolocate(min_, max_, base=5):
    """
    Return auto-generated levels in the provided interval. Minimum and maximum
    values will be rounded to the nearest number divisible by the specified base.
    """
    _round = lambda x: base*round(float(x)/base)
    return np.arange(_round(min_), _round(max_)+base/2, base)

def _cformat(self, mappable, cgrid=False, clocator=None, cminorlocator=None,
        cformatter=None, clabel=None, errfix=True, extend=.15, # in inches
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
      using a BoundaryNorm (for example) to determine colors between the levels
      (see: https://stackoverflow.com/q/42723538/4970632).
    * Workaround is to make sure locators are in vmin/vmax range EXCLUSIVELY;
      cannot match/exceed values.
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
    # Test if mappable is iterable; if so, user wants to use colorbar to label lines
    try: iter(mappable)
    except TypeError:
        fromlines = False
    else: # catch exception: the only iterable matplotlib objects are Container objects
        if isinstance(mappable, mcontainer.Container):
            fromlines = False
        else:
            fromlines = True
    csettings = {'spacing':'uniform', 'cax':cax, 'use_gridspec':True, # use space afforded by entire axes
        'extend':'both', 'orientation':orientation, 'drawedges':cgrid} # this is default case unless mappable has special props
    # Update with user-kwargs
    csettings.update(**kwargs)
    if hasattr(mappable, 'extend'):
        csettings.update({'extend':mappable.extend})
    # Option to generate colorbar/colormap from line handles
    # * Note the colors are perfect if we don't extend them by dummy color on either side,
    #   but for some reason labels for edge colors appear offset from everything
    # * Too tired to figure out why so just use this workaround
    if fromlines:
        if values is None:
            raise ValueError("Must pass values corresponding to list of handle objects.")
        if len(mappable)!=len(values):
            raise ValueError("Number of values should equal number of handles.")
        values = np.array(values) # needed for below
        colors = [h.get_color() for h in mappable]
        colors = ['#ffffff'] + colors + ['#ffffff']
        colormap = mcolors.LinearSegmentedColormap.from_list('tmp', colors)
        # levels = np.concatenate((values[0]-np.diff(values[:2]), # get "edge" values
        #     (values[1:]+values[:-1])/2, values[-1]+np.diff(values[-2:])))
        levels = np.concatenate((values[0]-np.diff(values[:2])/2, # get "edge" values
            (values[1:]+values[:-1])/2, values[-1]+np.diff(values[-2:])/2))
        values = np.concatenate((values[:1]-np.diff(levels[:2]), values,
            values[-1:]+np.diff(levels[-2:]))) # this works empirically; otherwise get
            # weird situation with edge-labels appearing on either side
        mappable = plt.contourf([[0,0],[0,0]], levels=levels, cmap=colormap, extend='neither',
                norm=BoundaryNorm(values)) # workaround
        if clocator is None: # in this case, easy to assume user wants to label each value
            # clocator = values # restore this if values weren't padded
            clocator = values[1:-1]
    # Determine tick locators for major/minor ticks
    # Can pass clocator/cminorlocator as the *jump values* between the mappables
    # vmin/vmax if desired
    normfix = False # whether we need to modify the norm object
    locators = [] # put them here
    for locator in (clocator,cminorlocator):
        # Create a preliminary object first
        if isinstance(locator, mticker.Locator):
            pass # nothing to do
        elif locator is None: # final; if cminorlocator wasn't set, don't bother at all
            if locator is cminorlocator:
                locators.append(mticker.NullLocator()) # don't set it
                continue
            locator = mticker.AutoLocator() # really basic
        else: # set based on user input
            try: iter(locator)
            except TypeError:
                locator = _autolocate(mappable.norm.vmin, mappable.norm.vmax, locator)
            locator = mticker.FixedLocator(np.array(locator)) # basic
        # Then modify to work around that mysterious error, and to prevent annoyance
        # where minor ticks extend beyond the "extend" triangles
        values = np.array(locator.tick_values(mappable.norm.vmin, mappable.norm.vmax)) # get the current values
            # need to use tick_values instead of accessing "locs" attribute because
            # many locators don't have these attributes; require norm.vmin/vmax as input
        values_min = np.where(values>=mappable.norm.vmin)[0]
        values_max = np.where(values<=mappable.norm.vmax)[0]
        if len(values_min)==0 or len(values_max)==0:
            raise ValueError(f"Ticks are not within the colorbar range {mapable.norm.vmin:.3g} to {mappable.norm.vmax:.3g}.")
        values_min, values_max = values_min[0], values_max[-1]
        values_inrange = values[values_min:values_max+1]
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
    extend = extend/(cax.figure.width*np.diff(getattr(cax.get_position(),interval))[0]-2*extend)
    csettings.update({'extendfrac':extend}) # width of colorbar axes and stuff
    cb = self.figure.colorbar(mappable, ticks=locators[0], format=cformatter, **csettings)
    # The ticks/ticklabels basic properties
    for t in axis.get_ticklabels(which='major'):
        t.update(rc('ticklabels'))
    axis.set_tick_params(which='major', **rc('ctick'))
    axis.set_tick_params(which='minor', **rc('ctickminor'))
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
    axis.label.update(rc('label'))
    if clabel is not None:
        axis.label.update({'text':clabel})
    # Fix pesky white lines between levels + misalignment with border due to rasterized blocks
    cb.solids.set_rasterized(False)
    cb.solids.set_edgecolor('face')
    # Make edges/dividers consistent with axis edges
    if cb.dividers is not None:
        cb.dividers.update(rc('cgrid'))
    cb.outline.update(rc('outline'))
    # Update the tickers colorbar
    # cb.formatter = cformatter # fucks shit up -- why?
    # cb.update_ticks() # doesn't actually update the formatter
    # axis.set_major_formatter(cformatter) # fucks shit up -- why?
    return cb

def _lformat(self, handles=None, multi=None, handlefix=False, **kwargs): #, settings=None): # can be updated
    """
    Function for formatting legend-axes (invisible axes with centered legends on them).
    Should update my legend function to CLIP the legend box when it goes outside axes area, so
    the legend-width and bottom/right widths can be chosen propertly/separately.
    """
    # First get legend settings (usually just one per plot so don't need to declare
    # this dynamically/globally)
    lsettings = rc('legend')
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
    lsettings.update(**kwargs)
    # Setup legend text properties
    tsettings = rc('ticklabels')
    if 'fontsize' in lsettings:
        tsettings['size'] = lsettings['fontsize']
    # Setup legend handle properties
    hsettings = {}
    candidates = ['linewidth','color'] # candidates for modifying legend objects
    for candidate in candidates:
        if candidate in kwargs:
            hsettings[candidate] = lsettings.pop(candidate)
    if not handlefix:
        hsettings = {}
    # Detect if user wants to specify rows manually
    # Gives huge latitude for user input:
    #   1) user can specify nothing and align will be inferred (list of iterables is False,
    #      use consecutive legends)
    #   2) user can specify align (needs list of handles for True, list of handles or list
    #      of iterables for False and if the former, will turn into list of iterables)
    if multi is None: # automatically guess
        try: iter(handles[0])
        except TypeError:
            multi = False
        else: # catch exception: the only iterable matplotlib objects are Container objects
            if isinstance(handles[0], mcontainer.Container):
                multi = False
            else:
                multi = True
    else: # standardize format based on "multi" input
        try: iter(handles[0])
        except TypeError:
            implied = False # no need to fix
        else:
            if not isinstance(handles[0], mcontainer.Container):
                implied = False
            else:
                implied = True # no need to fix
        if multi and not implied: # apply columns
            if 'ncol' not in lsettings:
                raise ValueError("Need to specify number of columns with ncol.")
            handles = [handles[i*lsettings['ncol']:(i+1)*lsettings['ncol']]
                    for i in range(len(handles))] # to list of iterables
        elif not multi and implied:
            handles = [handle for iterable in handles for handle in iterable]

    # Now draw legend, with two options
    # 1) Normal legend, just draw everything like normal
    # Re-orders handles to be row-major, is only difference
    if not multi:
        if 'ncol' not in lsettings:
            raise ValueError("Need to specify number of columns with ncol.")
        newhandles = []
        ncol = lsettings['ncol'] # number of columns
        handlesplit = [handles[i*ncol:(i+1)*ncol] for i in range(len(handles)//ncol+1)] # split into rows
        nrowsmax, nfinalrow = len(handlesplit), len(handlesplit[-1]) # max possible row count, and columns in final row
        nrows = [nrowsmax]*nfinalrow + [nrowsmax-1]*(lsettings['ncol']-nfinalrow)
            # e.g. if 5 columns, but final row length 3, columns 0-2 have N rows but 3-4 have N-1 rows
        for col,nrow in enumerate(nrows): # iterate through cols
            newhandles.extend(handlesplit[row][col] for row in range(nrow))
        handles = newhandles
        if hasattr(self, '_legend'): # this Method was renamed to "legend" on an Axes; old version is in _legend
            leg = self._legend(handles=handles, **lsettings) # includes number columns
        else: # this Method did not override the original "legend" method
            leg = self.legend(handles=handles, **lsettings) # includes number columns
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
            raise ValueError("Need to know labelspacing before drawing multi-row legends. Add it to settings.")
        if 'size' not in tsettings:
            raise ValueError("Need to know fontsize before drawing multi-row legends. Add it to settings.")
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
                    frameon=False, borderpad=0, loc='center', **lsettings) # _lformat is overriding original legend Method
            else:
                leg = self.legend(handles=hs, ncol=len(hs), bbox_to_anchor=bbox,
                    frameon=False, borderpad=0, loc='center', **lsettings) # not overriding original Method
            leg.legendPatch.update(rc('outline')) # or get_frame()
            for obj in leg.legendHandles:
                obj.update(hsettings)
            for t in leg.texts:
                t.update(tsettings) # or get_texts()
            legends.append(leg)
        for l in legends[:-1]:
            self.add_artist(l) # because matplotlib deletes previous ones...
    return legends

def _format(self,
    transparent=True, hatch=None, color='w',
    coastlines=True, continents=False, # coastlines and continents
    latlabels=[0,0,0,0], lonlabels=[0,0,0,0], latlocator=None, lonlocator=None, # latlocator/lonlocator work just like xlocator/ylocator
    xgrid=None, ygrid=None, # gridline toggle
    xdates=False, ydates=False, # whether to format axis labels as long datetime strings
    xtickminor=None, ytickminor=None, xgridminor=None, ygridminor=None, # minor ticks/grids; if ticks off, grid will be off
    xspineloc=None, yspineloc=None, # deals with spine options
    xtickloc=None, ytickloc=None, # tick location
    xtickdir=None, ytickdir=None, # change ytick/xtick location; can be in, out, or inout (left-right-up-down depends on spine to which this applied)
    xtickrange=None, ytickrange=None, # limit regions where we assign ticklabels to major-ticks
    xlim=None, ylim=None, xscale=None, yscale=None, xreverse=False, yreverse=False, # special properties
    xlabel=None, ylabel=None, # axis labels
    suptitle=None, suptitlepos=None, title=None, titlepos=None, abc=False, abcpos=None, abcformat='', padding=0.1,
    # TODO: add options for allowing UNLABELLED major ticklines -- maybe pass a special Formatter?
    xlocator=None, xminorlocator=None, ylocator=None, yminorlocator=None, # locators, or derivatives that are passed to locators
    xformatter=None, yformatter=None): # formatter
    # legend=False, handles=None, # legend options; if declared and have no handles with 'label', does nothing
    # settings=None): # remaining are passed to Settings()
    """
    Function for formatting axes of all kinds; some arguments are only relevant to special axes, like 
    colorbar or basemap axes. By default, simply applies the dictionary values from settings() above,
    but can supply many kwargs to further modify things.
    * Add options for datetime handling; note possible date axes handles are TimeStamp (pandas),
      np.datetime64, DateTimeIndex; can fix with fig.autofmt_xdate() or manually set options; uses
      ax.is_last_row() or ax.is_first_column(), so should look that up.
    * Problem is there is no autofmt_ydate(), so really should implement my own version of this.
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
        fig.suptitle.update({'ha':'center', 'va':'baseline', **rc('title')})
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
    if titlepos=='left':
        self.title.update({'position':(0,1), 'ha':'left'})
    elif titlepos=='right':
        self.title.update({'position':(1,1), 'ha':'right'})
    elif titlepos=='inside':
        self.title.update({'position':(0.5,1-padding/self.height),
            'transform':self.transAxes, 'va':'top'})
    elif isinstance(titlepos,str):
        raise ValueError(f"Unknown title position: {titlepos}.")
    elif titlepos is not None: 
        self.title.update({'position':titlepos, 'ha':'center', 'va':'center'})
    # Create axes numbering
    if self.number is not None and abc:
        abcedges = abcformat.split('a')
        self.abc = self.text(0, 1, abcedges[0] + ascii_lowercase[self.number-1] + abcedges[-1],
                transform=self.title._transform) # optionally include paren
        self.abc.update({'ha':'left', 'va':'baseline', **rc('abc')})
        if abcpos=='inside':
            self.abc.update({'position':(padding/self.width, 1-padding/self.height),
                'transform':self.transAxes, 'ha':'left', 'va':'top'})
        elif isinstance(abcpos,str):
            raise ValueError(f"Unknown abc position: {abcpos}.")
        elif abcpos is not None:
            self.abc.update({'position':abcpos, 'ha':'left', 'va':'top'})
    # Color setup, optional hatching in background of axes
    self.figure.patch.set_alpha(0) # make transparent by default
    self.patch.set_alpha(1) # make not transparent
    if hatch: # non-empty string or not none
        self.fill_between([0,1], 0, 1, hatch=hatch,
                facecolor='none', edgecolor='k', 
                transform=self.transAxes) # fill hatching in background
    else: # fill with color instead
        self.patch.set_color(color)
        self.patch.set_zorder(-1)
        self.patch.set_clip_on(False)

    #--------------------------------------------------------------------------
    # Process projection/map axes
    #--------------------------------------------------------------------------
    # Basemap axes setup
    # For now force background to never be transparent, logic being that the background
    # itself contains information (e.g. "this is ocean" or "this is land")
    if getattr(self, 'm', None) is not None: # the attribute exists and isn't None
        # Coastlines, parallels, meridians
        m = self.m
        if coastlines:
            p = m.drawcoastlines(**rc('coastlines'))
        if continents:
            p = m.fillcontinents(**rc('continents'))
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
                        _.update(rc('ticklabels'))
                    else:
                        _.set_clip_on(True) # no gridlines past boundary
                        _.update(rc('lonlatlines'))
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
                        _.update(rc('ticklabels'))
                    else:
                        _.set_clip_on(True) # no gridlines past boundary
                        _.update(rc('lonlatlines'))
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
            p = m.drawmapboundary(fill_color=color, **rc('spine')) # set fill_color to 'none' to make transparent
            p.set_rasterized(False) # not sure about this; might be rasterized
            p.set_clip_on(False) # so edges of LINE denoting boundary aren't cut off
        else: # use the settings to apply to Axes patch; the basemap API fails here
            self.patch.set_facecolor(color)
            self.patch.set_edgecolor('none')
            for spine in self.spines.values():
                spine.update(rc('spine'))
        # p.set_zorder(-1) # is already zero so this is dumb
        # p.set_facecolor(None) # don't do this as it will change the color
        return # skip everything else
    # Cartopy axes setup
    if isinstance(self, cgeoaxes.GeoAxes): # the main GeoAxes class; others like GeoAxesSubplot subclass this
        print("WARNING: Cartopy axes setup not yet implemented.")
        # self.add_feature(cfeature.COASTLINE, **rc('outline'))
        # self.outline_patch.update(rc('outline'))
        return

    #--------------------------------------------------------------------------
    # Process normal axes, various x/y settings individually
    #--------------------------------------------------------------------------
    # Axes scaling, limits, and reversal options (alternatively, supply
    # your own xlim/ylim that go from high to low)
    # if xscale is not None: self.set_xscale(xscale, **rc('xscale'))
    # if yscale is not None: self.set_yscale(yscale, **rc('yscale'))
    if xscale is not None: self.set_xscale(xscale)
    if yscale is not None: self.set_yscale(yscale)
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
        # Axis spine and tick locations
        # May change settings if twin axes is present
        for spine, side in zip((self.spines[s] for s in sides), sides):
            # Line properties
            spine.update(rc('spine'))
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
                # Special location
                if side==sides[0]: # move the left/bottom spine onto the specified location, with set_position
                    spine.set_visible(True)
                    spine.set_position(spineloc) # just gets passed to the function; options include 
                        # 'zero', 'center', and tuple with (units, location) where units can be axes, data, or outward
                else:
                    spine.set_visible(False)

        # Tick/tick gridline properties
        # * Some weird issue seems to cause set_tick_params to reset/forget that the grid
        #   is turned on if you access tick.gridOn directly, instead of passing through tick_params.
        #   So calling _format() a second time will remove the lines
        # * Can specify whether the left/right/bottom/top spines get ticks; sides will be 
        #   group of left/right or top/bottom
        # * Includes option to draw spines but not draw ticks on that spine, e.g.
        #   on the left/right edges
        # * Weirdly gridOn is undocumented feature
        ticklocs = sides if tickloc=='both' else () if tickloc=='neither' else tickloc
        if ticklocs is None: # only turn off ticks if the spine is invisible; don't explicitly turn on
            sides_kw = {side: False for side in sides if not self.spines[side].get_visible()}
        else:
            sides_kw = {side: self.spines[side].get_visible() and side in ticklocs for side in sides}
        grid_kw = {} if grid is None else {'gridOn': grid}
        gridminor_kw = {} if gridminor is None else {'gridOn': gridminor and tickminor}
        major_kw = rc('tick') if tickdir is None else dict(rc('tick'), direction=tickdir)
        minor_kw = rc('tickminor') if tickdir is None else dict(rc('tickminor'), direction=tickdir)
        axis.set_tick_params(which='major', **sides_kw, **grid_kw, **major_kw)
        axis.set_tick_params(which='minor', **sides_kw, **gridminor_kw, **minor_kw) # have length
        for tick in axis.majorTicks:
            tick.gridline.update(rc('grid'))
        for tick in axis.minorTicks:
            tick.gridline.update(rc('gridminor'))
            if tickminor is not None:
                tick.set_visible(tickminor)

        # Label properties
        axis.label.update(rc('label'))
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
            # axis.set_ticklabels(tickformatter) # ...FixedFormatter alone has issues
            axis.set_major_formatter(mticker.FixedFormatter(tickformatter)) # list of strings
        for t in axis.get_ticklabels():
            t.update(rc('ticklabels'))
    return # we're done

#-------------------------------------------------------------------------------
# Primary plotting function; must be used to create figure/axes if user wants
# to use the other features
#-------------------------------------------------------------------------------
def subplots(array=None, nrows=1, ncols=1, emptycols=None, emptyrows=None,
        tight=False, # whether to set up tight bbox from gridspec object
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
        bottompanel=False, bottompanels=False, rightpanel=False, rightpanels=False,
        bwidth=0.25, bspace=0, rwidth=0.25, rspace=0, # default to no space between panels
            # also need to consider drawing legend *centered* so can put e.g. colorbar next to
            # legend in a bottompanel and have it not look super weird
        lwidth=0.15, cwidth=0.25, cspace=0.8, cshrink=0.9,
            # spacing for colorbar text, colorbar axes width, legend axes with, and padding for interior ABC label
        # abcpad=0.1, titlepad=0.1,
        #     # will delete this stuff
        innerpanels=None, # same as below; list of numbers where we want subplotspecs
        whichpanels=None, hsep=None, wsep=None, hwidth=None, wwidth=None,
        maps=None, # set maps to True for maps everywhere, or to list of numbers
        package='basemap', projection='cyl', projectiondict={}, **projectionkw): # for projections; can be 'basemap' or 'cartopy'
    """
    Special creation of subplots grids, allowing for arbitrarily overlapping 
    axes objects. Will return figure handle and axes object. Can also use exaclty as
    plt.subplot, with some added convenience features. Use panelx, panely as 
    templates for figure-wide legends/colorbars; if you actually want ax smallish panel
    to plot data, just use width/height ratios. don't bother with kwarg dictionaries because
    it is really, really customized; we do use them in formatting function though.

    NOTES:
        * Matplotlib set_aspect option seems to behave strangely on some plots (trend-plots from
        SST paper); for this reason we override the fix_aspect option provided by basemap and
        just draw figure with appropriate aspect ratio to begin with.
        * Otherwise get weird differently-shaped subplots that seem to make no sense.

    TODO:
        * If want **differing spacing** between certain subplots, must use 
        GridSpecFromSubplotSpec with its **own** hspace/wspace properties...requires
        some consideration, but don't worry until you come across it. Could make axes
        all from scratch instead of with SubplotSpec/GridSpec.
        * Figure size should be constrained by WIDTH/HEIGHT/ASPECT RATIO OF AXES, 
        with everything else provided; must pick two of these; default, is to pick
        width and aspect, but if height provided, width is ignored instead. no way
        right now to determine aspect by constraining width/height
        * For spanning axes labels, right now only detect **x labels on bottom**
        and **ylabels on top**; generalize for all subplot edges.

    NOTE that, because I usually format all subplots with .format() method, shared axes should already
    end up with the same axis limits/scaling/majorlocators/minorlocators... the sharex and sharey detection algorithm
    really is just to get INSTRUCTIONS to make the ticklabels and axis labels INVISIBLE for certain axes.
    """
    # Fill in some basic required settings
    if package not in ['basemap','cartopy']:
        raise ValueError("Plotting package must be one of basemap or cartopy.")
    if width is None and height is None: # at least one must have value
        width = 5
    if aspect is None:
        aspect = 1.5 # try this as a default; square plots (1) look too tall to me
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
    # First make sure user didn't mess up
    projectionkw.update(projectiondict)
    if projectionkw and not maps:
        raise ValueError(f"Unknown arguments: {', '.join(projectionkw.keys())}. If you want to "
                "create maps, you must change the \"maps\" kwarg.")
    # Basemap stuff
    if maps and package=='basemap':
        # Just determine the correct aspect ratio
        # To do this need to instantiate basemap object which has dual utility of
        # telling us before draw-time whether any kwpair in projectionkw is bad
        projectionkw.update({'fix_aspect':True})
        mexample = mbasemap.Basemap(projection=projection, **projectionkw)
        aspect = (mexample.urcrnrx-mexample.llcrnrx)/(mexample.urcrnry-mexample.llcrnry)
        print(f"Forcing aspect ratio: {aspect:.3g}")
    # Cartopy stuff; get crs instance, and create dictionary to add to add_subplot calls
    cartopy_kw = {}
    if maps and package=='cartopy':
        # Get the projection instance from a string and determine the
        # correct aspect ratio (TODO) also prevent auto-scaling
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
        if projection not in crs_dict:
            raise ValueError(f"For cartopy, projection must be one of the following: "
                + f"{', '.join(crs_dict.keys())}.")
        cartopy_kw = {'projection': crs_dict[projection](**projectionkw)}
    elif not maps:
        projection = None
    # Aspect ratio test
    if maps:
        wtest = [0] if wratios is None else wratios
        htest = [0] if hratios is None else hratios
        if any(w!=wtest[0] for w in wtest) or any(h!=htest[0] for h in htest):
           raise ValueError("Aspect ratios of Axes boxes must be identical. "
                   "Improve this library to remove this restriction.")

    #--------------------------------------------------------------------------
    # Panel considerations; keep things general for future imporvements
    #--------------------------------------------------------------------------
    # Detect conflicts
    if bottomcolorbar and rightcolorbar:
        raise ValueError('You can only specify one global colorbar.')
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

    # Setup up figure and plot areas
    fig = plt.figure(FigureClass=Figure, figsize=figsize)
    fig.height = height
    fig.width = width
    # Outer area; includes regions for special panels
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
    fig.left = ileft
    fig.bottom = ibottom
    fig.right = iright
    fig.top = itop
    fig.lwidth = lwidth
    fig.cwidth = lwidth
    fig.cspace = cspace
    fig.bwidth = bwidth
    fig.bspace = bspace
    fig.rwidth = rwidth
    fig.rspace = rspace
    fig.gridspec = GS # add to figure, for later reference
    # Initialize some more stuff
    fig.bottomlegend = None
    fig.bottomcolorbar = None
    fig.rightcolorbar = None
    fig.bottompanel = None
    fig.rightpanel = None
    # Main plotting area
    # Will draw individual axes within this GridSpec object later
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
    # Find axes that have special properties
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
    # Helper function for creating panelled axes
    # Will take in some arguments from parent function so don't have to pass them
    # every time
    def panelfactory(subspec, whichpanels=whichpanels,
            hsep=hsep, wsep=wsep,
            hwidth=hwidth, wwidth=wwidth,
            **kwarg):
        # Checks and defaults
        # Defaults format the plot to have tiny spaces between subplots
        # and narrow panels
        if whichpanels is None:
            whichpanels = 'b' # default is a bottom panel, cuz why not
        if hsep is None:
            hsep = 0.1
        if wsep is None:
            wsep = 0.1
        if hwidth is None:
            hwidth = 0.5
        if wwidth is None:
            wwidth = 0.5
        if any(s not in 'tblr' for s in whichpanels) or whichpanels=='':
            raise ValueError("Whichpanels argument can contain characters l (left), r (right), "
                    "b (bottom), or t (top).")
        # Determine number of rows/columns, position of the
        # main/primary axes, and potential position of corners that
        # we want to ignore for now
        nrows_i, ncols_i = 1, 1
        for s in whichpanels:
            if s in ('b','t'):
                ncols_i += 1
            if s in ('l','r'):
                nrows_i += 1
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
        if hsep is not None:
            height_i = np.diff(bbox_i.intervaly)[0]*height - hsep*(nrows_i-1)
            hsep = hsep/(height_i/nrows_i)
        if wsep is not None:
            width_i = np.diff(bbox_i.intervalx)[0]*width - wsep*(ncols_i-1)
            wsep = wsep/(width_i/ncols_i)
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
                wspace        = wsep,
                hspace        = hsep,
                width_ratios  = wwidth_ratios,
                height_ratios = hwidth_ratios,
                )
        ax_ybase, ax_xbase = None, None
        if 'sharex' in kwarg:
            kwarg.pop('sharex')
        if 'sharey' in kwarg:
            kwarg.pop('sharey')
        for r in range(nrows_i)[::-1]:
            for c in range(ncols_i):
                if (r,c) not in bad_pos:
                    # Add the subplot first
                    iax = fig.add_subplot(gs_i[r,c], **kwarg)
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
        if len(axlist)==1:
            axlist = axlist[0] # so user doesn't have to unnecessarily index
        # Then should add panels as children of the main axes
        # * Can gain access to panel
        # * Unsure how the hell to do this
        return axmain, axlist

    # Base axes; to be shared with other axes as ._sharex, ._sharey attributes
    ax_panels = [] # empty for now
    ax = num_axes*[None] # list of axes
    allgroups_base = []
    if sharex: allgroups_base += xgroups_base
    if sharey: allgroups_base += ygroups_base
    for i in allgroups_base: # this is just list of indices in axes_ids, yrange, etc.
        if ax[i] is None: # not created as a x-base already, for example
            if i in innerpanels_ids:
                ax[i], axp = panelfactory(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        **cartopy_kw) # main axes handle
                ax_panels.append(axp) # panels in list; items are lists themselves
                    # if user wanted more than one panels, e.g. left and bottom
            else:
                ax[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        **cartopy_kw)

    # Dependent axes
    for i in range(num_axes):
        # Detect if we want to share this axis with another
        # If so, get that axis
        sharex_ax, sharey_ax = None, None # by default, don't share with other axes objects
        if sharex:
            igroup = np.where([i in g for g in xgroups])[0] # np.where works on lists
            if igroup.size==1:
                sharex_ax = ax[xgroups_base[igroup[0]]]
                if sharex_ax is None:
                    raise ValueError('Something went wrong; shared x axes was not already drawn.')
            elif igroup.size>1:
                raise ValueError(f'Something went wrong; axis {i:d} belongs to multiple groups.')
        if sharey:
            igroup = np.where([i in g for g in ygroups])[0] # np.where works on lists
            if igroup.size==1:
                sharey_ax = ax[ygroups_base[igroup[0]]]
                if sharey_ax is None:
                    raise ValueError('Something went wrong; shared x axes was not already drawn.')
            elif igroup.size>1:
                raise ValueError(f'Something went wrong; axis {i:d} belongs to multiple groups.')

        # Draw axes, and add to list
        if ax[i] is not None:
            # Axes is a BASE and has already been drawn, but might still need to add
            # a _sharex or _sharey property; e.g. is bottom-axes of three-column plot
            # and we want it to share a y-axis
            if ax[i] is not sharex_ax:
                ax[i]._sharex = sharex_ax
            else:
                sharex_ax = None
            if ax[i] is not sharey_ax:
                ax[i]._sharey = sharey_ax
            else:
                sharey_ax = None
        else:
            # Virgin axes; these are not an x base or a y base
            # We draw them now and pass the sharex/sharey as kwargs
            if i in innerpanels_ids:
                ax[i], axp = panelfactory(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        sharex=sharex_ax, sharey=sharey_ax,
                        **cartopy_kw)
                ax_panels.append(axp) # panels in list; items are lists themselves
                    # if user wanted more than one panels, e.g. left and bottom
            else:
                ax[i] = fig.add_subplot(gs[slice(*yrange[i,:]), slice(*xrange[i,:])],
                        sharex=sharex_ax, sharey=sharey_ax,
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
            s.colorbar = MethodType(_cformat, s)
            s.legend = MethodType(_lformat, s)
        if len(specs)==1: specs = specs[0]
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
            s.colorbar = MethodType(_cformat, s)
            s.legend = MethodType(_lformat, s)
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
    # Dynamically add properties and methods to Axes objects, other setup
    #--------------------------------------------------------------------------
    # Some functions that need to be declared here
    # * TODO WARNING: If drawing more than x lines and don't explicitly set every item below
    #   its value will revert to the cycler; this may give unexpected results
    # * To print current cycle, use list(next(ax._get_lines.prop_cycler)['color'] for i in range(10))
    def _setup(self):
        # Set up some properties
        # self.set_prop_cycle(propcycle)
        self.bottomlegend, self.bottomcolorbar, self.rightcolorbar = False, False, False # identifiers
        self.xspine_override, self.yspine_override = None, None
        self.number = None # default
        self.projection = None # default
        self.width  = np.diff(self._position.intervalx)*width # position is in figure units
        self.height = np.diff(self._position.intervaly)*height
        self.m = None # optional basemap instance
        # Set up some methods
        # TODO Note severe change; change twiny to mean "share the x-axis but
        # now make two y-axes"; this makes way more sense to me than default
        self._legend = self.legend # custom legend drawing on any axes
        self._twinx  = self.twinx
        self._twiny  = self.twiny
        self._text   = self.text
        self.legend  = MethodType(_lformat, self) # MethodType approach
        self.twinx   = MethodType(_twinx, self)
        self.twiny   = MethodType(_twiny, self)
        self.text    = MethodType(_text, self)
        self.format  = MethodType(_format, self)
    def _twinx(self, **kwargs):
        # Create secondary x-axes
        a = self._twiny(**kwargs)
        _setup(a) # basic setup
        self.xspine_override = 'top'
        a.xspine_override = 'bottom'
        a.yspine_override = 'neither'
        return a
    def _twiny(self, **kwargs):
        # Create secondary y-axes
        a = self._twinx(**kwargs)
        _setup(a) # basic setup
        self.yspine_override = 'left'
        a.yspine_override = 'right'
        a.xspine_override = 'neither'
        return a
    # Create some attributes; analagous to default behavior where title 'exists' but is not visible
    # or is empty string, methods in _format will make them visible
    for i,a in enumerate(ax):
        # Basic setup
        _setup(a) # default methods and stuff
        a.number = i+1 # know axes number ahead of time; start at 1
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
        if not maps:
            a._pcolormesh = a.pcolormesh
            a._pcolor     = a.pcolor
            a._contourf   = a.contourf # copy the old methods
            a._contour    = a.contour # copy the old methods
            a.pcolor     = MethodType(_pcolor, a)
            a.pcolormesh = MethodType(_pcolormesh, a)
            a.contourf   = MethodType(_contourf, a)
            a.contour    = MethodType(_contour, a)
        elif package=='basemap':
            a.m = mbasemap.Basemap(ax=a, projection=projection, **projectionkw)
            a.m.projection = projection # save projection here
            a.m._pseudocyl = ['moll','robin','eck4','kav7','sinu','mbtfpq','vandg','hammer']
            a.m._pcolormesh = a.m.pcolormesh
            a.m._pcolor     = a.m.pcolor
            a.m._contour    = a.m.contour
            a.m._contourf   = a.m.contourf # save old methods
            a.m._scatter   = a.m.scatter
            a.m._plot   = a.m.plot
            a.m.pcolormesh  = MethodType(_pcolormesh_basemap, a.m) # function calls the old methods
            a.m.pcolor      = MethodType(_pcolor_basemap, a.m)
            a.m.contour     = MethodType(_contour_basemap, a.m)
            a.m.contourf    = MethodType(_contourf_basemap, a.m)
            a.m.scatter    = MethodType(_scatter_basemap, a.m)
            a.m.plot     = MethodType(_plot_basemap, a.m)
            # these cannot be passed directly as elements of a, for some reason
            # can't override original name of axes object, because m.pcolor calls m.ax.pcolor
            if projection in a.m._pseudocyl:
                a.m.drawmapboundary() # just initialize (see above)
            # b = a.m.drawmapboundary() # **rc('line')
            # b.set_zorder(-10) # make sure is in back
        elif package=='cartopy':
            print("WARNING: Use of this library with cartopy projections untested. May need work.")
            pass # TODO nothing so far
    # Repeat some of the above for the panel axes
    # Will only need the format method really
    for axp in ax_panels:
        try: iter(axp)
        except TypeError:
            axp = [axp]
        for a in axp:
            _setup(a)

    #--------------------------------------------------------------------------
    # Return results
    #--------------------------------------------------------------------------
    print('Figure setup complete.')
    if len(ax)==1:
        ax = ax[0]
    if len(ax_panels)==1:
        ax_panels = ax_panels[0] # this might itself be a list, if user e.g.
            # used the panel-API to create a single axes with bordering left/right panels
    if innerpanels:
        return fig, ax, ax_panels
    else:
        return fig, ax

#------------------------------------------------------------------------------
# Classes and Formatters
#------------------------------------------------------------------------------
class WarpNorm(mcolors.Normalize):
    """
    Pass as norm=<instance>, when declaring new pcolor or contourf objects.
    Base class for warping either or both sides of colormap.
    In most cases it may be more suitable to use BoundaryNorm.
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
        try: iter(value_scaled)
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
        # mcolors.Normalize.__init__(self, vmin, vmax, clip)
        super().__init__(vmin, vmax, clip)
        self.midpoint = midpoint
        if any(x is None for x in [vmin,vmax,midpoint]):
            raise ValueError("Must declare vmin, vmax, and midpoint explicitly.")

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
    def __init__(self, levels, midpoint=None, clip=False, **kwargs):
        # Very simple
        try: iter(levels)
        except TypeError:
            raise ValueError("Must call BoundaryNorm with your boundary vaues.")
        # mcolors.Normalize.__init__(self, min(levels), max(levels), clip, **kwargs)
        super().__init__(min(levels), max(levels), clip)
        self.midpoint = midpoint
        self.levels = np.array(levels)

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

def Formatter(sigfig=3, tickrange=None):
    """
    Format as a number, with N sigfigs, and trimming trailing zeros.
    Recall, must pass function in terms of n (number) and loc.
    """
    # Format definition
    def f(value, location, sigfig=sigfig, tickrange=tickrange):
        # Exit if not in tickrange
        if tickrange is not None:
            eps = value/1000
            if value<tickrange[0]-eps or value>tickrange[1]+eps:
                return '' # the data
        # Return special string
        # * Note CANNOT use "g" because "g" precision operator specifies count of
        #   significant digits, not places after decimal place.
        # * There is no format that specifies digits after decimal place AND trims trailing zeros.
        # print(string)
        string = f'{{{0}:.{sigfig:d}f}}'.format(value) # f-string compiled, then format run
        if '.' in string: # g-style trimming
            string = string.rstrip('0').rstrip('.')
        if string=='-0': # special case
            string = '0'
        return string
    # And create object
    return mticker.FuncFormatter(f)

def LatFormatter(sigfig=0, sine=False, cardinal=True):
    """
    Format latitude labels; can convert sine-lats back into lats (for
    areal weighting) and can apply N/S instead of postiive/negative.
    """
    def f(value, location, sine=sine, cardinal=cardinal, sigfig=sigfig):
        # Convert from sine to latitude number
        if sine:
            value = np.arcsin(value)*180/np.pi
        # Suffix to apply
        if cardinal and value<0:
            value *= -1
            suffix = 'S'
        elif cardinal:
            suffix = 'N'
        else:
            suffix = ''
        # string = formatstr % (value, string)
        # Return special string, as in Formatter method
        string = f'{{{0}:.{sigfig:d}f}}'.format(value) # f-string compiled, then call format
        if '.' in string: # g-style trimming
            string = string.rstrip('0').rstrip('.')
        if string=='-0': # special case
            string = '0'
        return string+suffix
    # And create object
    return mticker.FuncFormatter(f)

def LonFormatter(sigfig=0, cardinal=True):
    """
    Format latitude labels; can convert sine-lats back into lats (for
    areal weighting) and can apply N/S instead of postiive/negative.
    """
    def f(value, location, cardinal=cardinal, sigfig=sigfig):
        # Suffix to apply
        if cardinal and value<0:
            value *= -1
            suffix = 'W'
        elif cardinal:
            suffix = 'E'
        else:
            suffix = ''
        # string = formatstr % (value, string)
        # Return special string, as in Formatter method
        string = f'{{{0}:.{sigfig:d}f}}'.format(value) # f-string compiled, then call format
        if '.' in string: # g-style trimming
            string = string.rstrip('0').rstrip('.')
        if string=='-0': # special case
            string = '0'
        return string+suffix
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

