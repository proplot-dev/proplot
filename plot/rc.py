#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# This configures the global working environment
#------------------------------------------------------------------------------#
# First just make sure some dependencies are loaded
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.font_manager as mfonts
from glob import glob
from cycler import cycler
from matplotlib import matplotlib_fname
from types import MethodType
# Will add our own dictionary to the top-level matplotlib module, to go
# alongside rcParams
import matplotlib as mpl
# Local dependencies
from   .colors import shade, cmapcolors
# Default settings to be loaded when globals() is
# called without any arguments
rcDefaults = {
    'linewidth': 0.7,
    'color':     'k',
    'small':     8,
    'medium':    8,
    'large':     9,
    'ticklen':   4,
    'fontname':  'DejaVu Sans',
    'cycle':     'colorblind',
    }

#-------------------------------------------------------------------------------
# Settings management function
#-------------------------------------------------------------------------------
# I used to create this neat class but it was a stupid idea. Work with existing
# API, not against it. Just create a function that sets rcParam defaults with
# optional override. Then, better to make the format function add actual information
# to the plot and do nothing to change its style/color/etc.
# def globals(**kwargs): # fontname is matplotlib default
def globals(*args, verbose=False, **kwargs):
    """
    This has multiple uses, all rolled up into one function.
    1. *Initialize* everything with default settings. Creates special rcExtras
       dictionary assigned to 'mpl' module, just like rcParams, except the values
       in rcExtras have my own special naming and are read by my _format() function.
    2. *Set* rcParams and rcExtras parameters belonging to some category with a single dictionary,
       list of dictionaries, kwarg pairs, or all of the above.
    3. *Set* special global params that are applied to a bunch of different
       rcParams and rcExtras (e.g. 'color'), while leaving others (e.g. 'linewidth') alone.
    4. *Retrieve* a single value belong to a category.subcategory, or a dictionary
       of *all* or *filtered/selected* subcategory=value pairs for subcategories
       belonging to 'category'.
    Here's a quick list of rcParam categories (see: https://matplotlib.org/users/customizing.html)
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
    Some other notes:
    * Note the figure settings are used when printing interactively or just making
      the figure object, but the savefig ones are used when calling savefig.
    * Note that if *autoreload* is triggered/global defauls are reset, it seems that
      any options set with InlineBackend in ipython are ***ignored***. But if you query
      settings, options still present -- you just need to call nbsetup again.
    * Problem is the settings in rcParams are scattershot and missing some important
      ones we want, but should use it when we can; options won't be rescinded in future versions
    * Currently features that can't be rcParam'ed are *gridminor* properties, the
      *cgrid* idea (lines between colorbar colors), the continent/coastlines/lonlatlines
      settings, and some legend settings.
    TODO: Change default background axes color to 0072B2, from seaborn colorblind, just
    a bit darker than figure background; is nice.
    TODO: Expand idea of default 'color', 'linewidth', whatever and create more
    dictionaries with defaults settings, then update them. No more hardcoded values
    in the main text below.
    """
    # Helper function; adds stuff to rcParams or rcExtras in the
    # stylesheet style of category.subcategory = value, or very
    # occasionally category.subcategory.subsubcategory = value, in which case
    # user should input a dictionary {'subcategory.subsubcategory':value}
    def add(category, kwargs): # pass
        if not hasattr(mpl, 'rcExtras'):
            mpl.rcExtras = {} # initialize empty dictionary
        for subcategory,value in kwargs.items():
            if f'{category}.{subcategory}' in mpl.rcParams:
                mpl.rcParams[f'{category}.{subcategory}'] = value
            else:
                mpl.rcExtras[f'{category}.{subcategory}'] = value
    # Manage input, and intialization
    category = None # not necessarily anything
    if args and type(args[0]) is str: # second part of 'and' only tested if first part true
        category, *args = args # pull out category; but args might be a bunch of dictionaries
    newargs = [] # new list
    for arg in args:
        if isinstance(arg,dict):
            kwargs = {**kwargs, **arg} # add dictionary
        else:
            newargs += [arg] # just means we want to retrieve individual arguments
    args = newargs # retain old ones
    # *Retrieve* settings without changing any, if user passed a string 'category' name
    # and did not pass any kwargs for assignment as new settings
    if category is not None and not kwargs:
        if category=='globals':
            return {key.split('.')[-1]:val for key,val in mpl.rcExtras.items() if 'globals.' in key}
        dictionary = {}
        for catstring,value in {**mpl.rcParams, **mpl.rcExtras}.items():
            if catstring.split('.')[0]==category:
                dictionary[catstring.split('.')[1]] = value
        if args: # filter them
            dictionary = {k:v for k,v in dictionary.items() if k in args}
            if len(dictionary)==1:
                value, = dictionary.values() # turn a single value
                return value # early return, ugly
        if not dictionary:
            raise ValueError(f"Could not find settings for \"{category}\".")
        return dictionary
    # Now the section that applies settings
    if args:
        raise ValueError(f"Improper use of globals(). Only supply extra *args without any **kwargs.")
    #--------------------------------------------------------------------------#
    # *Initialize* default settings; that is, both rcParams and rcExtras
    # Requires processing in lines below
    category = category or 'globals' # None defaults to globals in this section
    if category=='globals' and not kwargs: # default settings
        if verbose: print('Resetting to defaults.')
        mpl.rcdefaults() # apply *builtin* default settings
        mpl.rcExtras = {} # reset extras
        add('globals', rcDefaults) # apply *custom* default settings
    #--------------------------------------------------------------------------#
    # *Update* rcParam settings or global settings; just add them to rcParams or rcExtras
    # If one of the 'global' settings was changed, requires processing in lines below
    # Also if any 'value' is a dictionary, update properties with that key as a
    # 'category' and each key-value air in the dictionary as kwargs
    else:
        updated = False
        for key,value in kwargs.items():
            if isinstance(value,dict):
                add(key, value) # add the 'value' as a kwarg dictionary
                if verbose: print(f"Added settings to category '{key}':", value)
            elif category=='globals' and key not in rcDefaults:
                raise ValueError(f"Key \"{key}\" unknown. Not a global property like \"color\" or \"linewidth\".")
            else:
                value = value if value!='default' or category!='globals' else rcDefaults[key]
                add(category, {key:value}) # easy peasy
                updated = True
                if verbose: print(f"Added setting to category '{category}':", {key:value})
        if category!='globals' or (category=='globals' and not updated):
            if verbose: print('No global properties changed.')
            return None # don't need to re-apply the globals
    #--------------------------------------------------------------------------#
    # *Apply* global settings; if this function was not called without arguments, then
    # only settings specifically requested to be changed, will be changed
    current = globals('globals') # the dictionary
    # Make a cycler for drawing a bunch of lines
    # Will add the color cycles from all the 
    dashes = ('-', '--', ':', '-.') # dash cycles; these succeed color changes
    colors = cycles.get(current['cycle'],None) # load colors from cyclers
    if colors is None:
        raise ValueError(f"Unknown color cycler designator \"{current['cycle']}\". Options are: "
            + ', '.join(f'"{k}"' for k in cycles.keys()) + '.') # cat strings
    if isinstance(colors[0],str) and colors[0][0]!='#': # fix; absolutely necessary (try without)
        colors = [f'#{color}' for color in colors]
    propcycle = cycler('color', [colors[i%len(colors)] for i in range(len(colors)*len(dashes))]) \
        + cycler('linestyle', [i for i in dashes for n in range(len(colors))])
    # if not isinstance(colors[0], str) and any(v>1 for v in colors[0]): # convert to hex
    #--------------------------------------------------------------------------#
    # In this *extra special instance*, search for all currently active
    # figures and apply set_prop_cycle on every single axes
    # * Calling plt.figure with 'num' as first argument or 'num=<num>' just
    #   brings that figure to foreground
    # * Tested this and it takes 1ms; no large deal
    figs = list(map(plt.figure, plt.get_fignums())) # from: https://stackoverflow.com/a/26485683/4970632
    for fig in figs:
        for ax in fig.axes:
            ax.set_prop_cycle(propcycle)
    #--------------------------------------------------------------------------#
    # First the rcParam settings
    # Here are ones related to axes and figure properties
    color, linewidth, ticklen, small, large, fontname = \
        current['color'], current['linewidth'], current['ticklen'], current['small'], current['large'], current['fontname']
    add('savefig', {'transparent':True, 'facecolor':'w', 'dpi':300,
        'directory':'', 'pad_inches':0, 'bbox':'standard', 'format':'pdf'}) # empty means current directory
    add('figure', {'facecolor':(.95,.95,.95,1), 'dpi':90, # matches nbsetup
        'max_open_warning':0, 'constrained_layout':False, 'autolayout':False}) # zero disables max open warning
    add('axes', {'xmargin':0, 'ymargin':0.05, 'labelsize':small, 'titlesize':large,
        'edgecolor':color, 'labelcolor':color, 'grid':True, 'linewidth':linewidth,
        'labelpad':3, 'prop_cycle':propcycle})
    add('font', {'size':small, 'family':'sans-serif', 'sans-serif':fontname}) # when user calls .text()
    add('text', {'color':color}) # when user calls .text()
    add('grid', {'linewidth':linewidth/2, 'alpha':0.1, 'color':color})
    for xy in 'xy':
        tickloc = {'x':{'bottom':True,'top':False},'y':{'left':True,'right':False}}[xy]
        add(f'{xy}tick',       {'labelsize':small, 'color':color, 'direction':'out'})
        add(f'{xy}tick.major', {'pad':2, 'width':linewidth, 'size':ticklen, **tickloc}) # size==length
        add(f'{xy}tick.minor', {'pad':2, 'width':linewidth, 'visible':True, 'size':ticklen/2, **tickloc})
    #--------------------------------------------------------------------------#
    # Ones related to stuff plotted inside axes
    # For styles, see: https://matplotlib.org/examples/api/joinstyle.html
    add('patch',    {'linewidth':linewidth,   'edgecolor':color}) # e.g. bars?
    add('hatch',    {'linewidth':linewidth,   'color':color})
    add('lines',    {'linewidth':linewidth*2,
        'markersize':linewidth*4, 'markeredgewidth':0,
        'dash_joinstyle':'miter', 'dash_capstyle':'projecting', # joinstyle opts: miter, round, bevel
        'solid_joinstyle':'miter', 'solid_capstyle':'projecting'}) # capstyle opts: butt, round, projecting
    add('markers',  {'fillstyle':'full'})
    add('scatter',  {'marker':'o'})
    add('legend',   {'framealpha':1, 'fancybox':False, 'frameon':False, # see: https://matplotlib.org/api/legend_api.html
        'labelspacing':0.5, 'handletextpad':0.5, 'handlelength':2, 'columnspacing':1,
        'borderpad':0.3, 'borderaxespad':0, 'numpoints':1, 'facecolor':'inherit'})
    # add('legend',   {'framealpha':0.6, 'fancybox':False, 'frameon':False,
    #     'labelspacing':0.5, 'columnspacing':1, 'handletextpad':0.5, 'borderpad':0.25})
    add('mathtext', {'default':'regular', 'bf':'sans:bold', 'it':'sans:it'}) # no italicization
        # this can only be accomplished with rcParams; impossible to specify with API!
    #--------------------------------------------------------------------------#
    # Next the settings with custom names
    # Some might be redundant, and should consider eliminating
    # Should begin to delete these and control things with rcParams whenever possible
    # Remember original motivation was realization that some rcParams can't be changes
    # by passing a kwarg (don't remember any examples)
    add('abc',         {'size':large, 'weight':'bold', 'color':color})
    add('title',       {'size':large, 'weight':'normal', 'color':color, 'fontname':fontname})
    add('suptitle',    {'size':large, 'weight':'normal', 'color':color, 'fontname':fontname})
    add('label',       {'size':small, 'weight':'normal', 'color':color, 'fontname':fontname})
    add('ticklabels',  {'size':small, 'weight':'normal', 'color':color, 'fontname':fontname})
    add('gridminor',   {'linestyle':'-', 'linewidth':linewidth/2, 'color':color, 'alpha':0.1})
    add('cgrid',       {'color':color, 'linewidth':linewidth})
    add('continents',  {'color':color, 'linewidth':0}) # make sure no lines!
    add('tickminor',   {'length':ticklen/2, 'width':linewidth, 'color':color})
    add('tick',        {'length':ticklen, 'width':linewidth, 'color':color})
    add('ctickminor',  {'length':ticklen/2, 'width':linewidth, 'color':color})
    add('ctick',       {'length':ticklen, 'width':linewidth, 'color':color})
    add('coastlines',  {'linewidth':linewidth, 'color':color})
    add('lonlatlines', {'linewidth':linewidth, 'linestyle':':', 'color':color, 'alpha':0.2})
    add('spine',       {'color':color, 'linewidth':linewidth})
    add('outline',     {'edgecolor':color, 'linewidth':linewidth})
    # add('xscale'); add('yscale'); add('contents', color='moccasin')
    if verbose: print('Applied global properties.')
    return None

#------------------------------------------------------------------------------
# Initialization; stuff called on import
# Adds colormap names and lists available font names
# * If leave 'name' empty in register_cmap, name will be taken from the
#   Colormap instance. So do that.
# * Note that **calls to cmap instance do not interpolate values**; this is only
#   done by specifying levels in contourf call, specifying lut in get_cmap,
#   and using LinearSegmentedColormap.from_list with some N.
# * The cmap object itself only **picks colors closest to the "correct" one
#   in a "lookup table**; using lut in get_cmap interpolates lookup table.
#   See LinearSegmentedColormap doc: https://matplotlib.org/api/_as_gen/matplotlib.colors.LinearSegmentedColormap.html#matplotlib.colors.LinearSegmentedColormap
# * If you want to always disable interpolation, use ListedColormap. This type
#   of colormap instance will choose nearest-neighbors when using get_cmap, levels, etc.
#------------------------------------------------------------------------------
# List the current font names, original version; works on Linux but not Mac, because can't find mac system fonts?
# fonts = mfonts.get_fontconfig_fonts()
# fonts = [mfonts.FontProperties(fname=fname).get_name() for fname in flist]
# List the system font names, smarter version
# See: https://olgabotvinnik.com/blog/2012-11-15-how-to-set-helvetica-as-the-default-sans-serif-font-in/
fonts = [font.split('/')[-1].split('.')[0] for font in # system fonts
            mfonts.findSystemFonts(fontpaths=None, fontext='ttf')] + \
        [os.path.basename(font.rstrip('.ttf')) for font in # hack-installed fonts
            glob(f"{matplotlib_fname().rstrip('matplotlibrc')}/fonts/ttf/*.ttf")]
fonts = sorted(set(fonts)) # unique ones only

#-------------------------------------------------------------------------------
# Register new colormaps; must come before registering the color cycles
_announcement = False
for _file in glob(f'{os.path.dirname(__file__)}/cmaps/*'):
    if '.rgb' in _file or '.hex' in _file:
        _name = os.path.basename(_file)[:-4]
        if _name not in plt.colormaps(): # don't want to re-register every time
            if '.rgb' in _file: # table or RGB values
                _load = {'hc':{'skiprows':1, 'delimiter':','},
                         'cb':{'delimiter':','},
                         'nc':{'skiprows':2}}.get(_name[:2],{}) # default empty
                try: _cmap = np.loadtxt(_file, **_load)
                except:
                    print(f'Failed to load {_name}.')
                    continue
                if (_cmap>1).any(): _cmap = _cmap/255
            else: # list of hex strings
                _cmap = [*open(_file)][0] # just a single line
                _cmap = _cmap.strip().split(',') # csv hex strings
                _cmap = np.array([mcolors.to_rgb(_) for _ in _cmap]) # from list of tuples
            _N = len(_cmap) # simple as that; number of rows of colors
            if 'lines' not in _name.lower():
                _N = 256-len(_cmap)%1 # do this until figure out why colors get segmented
            plt.register_cmap(cmap=mcolors.LinearSegmentedColormap.from_list(_name, _cmap, _N))
            plt.register_cmap(cmap=mcolors.LinearSegmentedColormap.from_list(_name+'_r', _cmap[::-1], _N))
            if not _announcement: # only do this if register at least one new map
                _announcement = True
                print("Registered colormaps.")

#-------------------------------------------------------------------------------
# Register colors by adding them to _colors_full_map
# * So far register opencolors and the "N most popular" xkcd colors; downloaded
#   them directly from the txt file.
_announcement = False
_opencolors = { "gray":   ["#f8f9fa", "#f1f3f5", "#e9ecef", "#dee2e6", "#ced4da", "#adb5bd", "#868e96", "#495057", "#343a40", "#212529"],
  "red":    ["#fff5f5", "#ffe3e3", "#ffc9c9", "#ffa8a8", "#ff8787", "#ff6b6b", "#fa5252", "#f03e3e", "#e03131", "#c92a2a"],
  "pink":   ["#fff0f6", "#ffdeeb", "#fcc2d7", "#faa2c1", "#f783ac", "#f06595", "#e64980", "#d6336c", "#c2255c", "#a61e4d"],
  "grape":  ["#f8f0fc", "#f3d9fa", "#eebefa", "#e599f7", "#da77f2", "#cc5de8", "#be4bdb", "#ae3ec9", "#9c36b5", "#862e9c"],
  "violet": ["#f3f0ff", "#e5dbff", "#d0bfff", "#b197fc", "#9775fa", "#845ef7", "#7950f2", "#7048e8", "#6741d9", "#5f3dc4"],
  "indigo": ["#edf2ff", "#dbe4ff", "#bac8ff", "#91a7ff", "#748ffc", "#5c7cfa", "#4c6ef5", "#4263eb", "#3b5bdb", "#364fc7"],
  "blue":   ["#e7f5ff", "#d0ebff", "#a5d8ff", "#74c0fc", "#4dabf7", "#339af0", "#228be6", "#1c7ed6", "#1971c2", "#1864ab"],
  "cyan":   ["#e3fafc", "#c5f6fa", "#99e9f2", "#66d9e8", "#3bc9db", "#22b8cf", "#15aabf", "#1098ad", "#0c8599", "#0b7285"],
  "teal":   ["#e6fcf5", "#c3fae8", "#96f2d7", "#63e6be", "#38d9a9", "#20c997", "#12b886", "#0ca678", "#099268", "#087f5b"],
  "green":  ["#ebfbee", "#d3f9d8", "#b2f2bb", "#8ce99a", "#69db7c", "#51cf66", "#40c057", "#37b24d", "#2f9e44", "#2b8a3e"],
  "lime":   ["#f4fce3", "#e9fac8", "#d8f5a2", "#c0eb75", "#a9e34b", "#94d82d", "#82c91e", "#74b816", "#66a80f", "#5c940d"],
  "yellow": ["#fff9db", "#fff3bf", "#ffec99", "#ffe066", "#ffd43b", "#fcc419", "#fab005", "#f59f00", "#f08c00", "#e67700"],
  "orange": ["#fff4e6", "#ffe8cc", "#ffd8a8", "#ffc078", "#ffa94d", "#ff922b", "#fd7e14", "#f76707", "#e8590c", "#d9480f"]}
mcolors.OPEN_COLORS = {} # create separate dictionary for them
for _name,_colors in _opencolors.items(): # iterate through json values
    for _i,_color in enumerate(_colors):
        mcolors.OPEN_COLORS[_name+str(_i)] = _color
_xkcdcolors = pd.read_csv(f"{os.path.dirname(__file__)}/xkcd.txt",
    delimiter='\t', header=None, usecols=(0,1)).values.tolist() # note numpy loadtxt fails
mcolors.XKCD_SORTED = {tup[0]:tup[1] for i,tup in enumerate(_xkcdcolors[::-1])
    if i<256} # filter to N most popular
colors = {**mcolors.XKCD_SORTED, **mcolors.OPEN_COLORS} # initialize the color dictionary
for name,value in colors.items():
    _color = name.split('xkcd:')[-1] # for now this is it
    if _color not in mcolors._colors_full_map:
        mcolors._colors_full_map[_color] = value # no more xkcd: now call them directly
        if not _announcement: # only do this if adding at least one new color
            _announcement = True
            print("Registered colors.")

#-------------------------------------------------------------------------------
# Helper function that returns global cycles
# Append it to pyplot module
def _gen_cycles():
    return cycles
plt.cycles = _gen_cycles

#-------------------------------------------------------------------------------
# Register color new color styles
cycles = {'default':['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'], # default V2 matplotlib
    # copied from stylesheets; stylesheets just add color themese from every
    # possible tool, not already present as a colormap
    'ggplot':['#E24A33', '#348ABD', '#988ED5', '#777777', '#FBC15E', '#8EBA42', '#FFB5B8'],
    'bmh': ['#348ABD', '#A60628', '#7A68A6', '#467821', '#D55E00', '#CC79A7', '#56B4E9', '#009E73', '#F0E442', '#0072B2'],
    'solarized': ['#268BD2', '#2AA198', '#859900', '#B58900', '#CB4B16', '#DC322F', '#D33682', '#6C71C4'],
    '538': ['#008fd5', '#fc4f30', '#e5ae38', '#6d904f', '#8b8b8b', '#810f7c'],
    'seaborn': ['#4C72B0', '#55A868', '#C44E52', '#8172B2', '#CCB974', '#64B5CD'],
    'pastel': ['#92C6FF', '#97F0AA', '#FF9F9A', '#D0BBFF', '#FFFEA3', '#B0E0E6'],
    'deep': ['#4C72B0', '#55A868', '#C44E52', '#8172B2', '#CCB974', '#64B5CD'],
    'muted': ['#4878CF', '#6ACC65', '#D65F5F', '#B47CC7', '#C4AD66', '#77BEDB'],
    'colorblind': ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#F0E442', '#56B4E9'],
    # copied using digital color meter from papers with pretty plots
    # note shade is a function from .colors
    'colorblind2':[shade(color, saturation=1.3) for color in # appears to just be pale colorblind sheme
        [(68,139,177), (200,126,72), (68,163,137), (229,220,124), (205,154,182)]],
    # created with online tools
    'cinematic': [(51,92,103), (255,243,176), (224,159,62), (158,42,43), (84,11,14)],
    'cinematic2': [[1,116,152], [231,80,0], [123,65,75], [197,207,255], [241,255,47]],
    }
cycles = {**cycles, **{_cycle.lower():cmapcolors(_cycle) for _cycle in
    ['Accent','Dark2','Set2','Set3','Paired','cbLines1','cbLines2']}} # add in the colormap cycles
for _cycle in cycles.values():
    if isinstance(_cycle[0],str) and _cycle[0][0]!='#': # fix; absolutely necessary (try without)
        _cycle[:] = [f'#{_}' for _ in _cycle] # modify contents; super cool trick
    if not isinstance(_cycle[0],str) and any(c>1 for tup in _cycle for c in tup):
        _cycle[:] = [tuple(np.array(_)/255) for _ in _cycle] # tuple of decimal RGB values
    _cycle[:] = [mcolors.to_hex(_, keep_alpha=False) for _ in _cycle] # standardize; will also convert hes to lower case

#-------------------------------------------------------------------------------
# Now call the function to configure params with default values
globals(verbose=True)

