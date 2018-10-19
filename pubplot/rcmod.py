#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# This configures the global working environment
# TODO: Allow ***seamless*** transitioning between color line cycles and
# colormaps. Can just register all colormaps with 10 'shades' that sample points
# alont the map, call with e.g. color='hclBlue0' or switch to that cycler with
# plot.globals('globals', cycle='hclBlue'), and maybe optionally use tuple-named
# cycles to configure sampling N times along that cycle. Work on this while on
# the plane! For colormaps/cycles without natural gradations, will just number
# along each index.
# TODO: Why don't I just use seaborn? Because that is meant to make *quick* and
# *pretty* plots, but my philosophy is exact opposite -- make tweaking stuff to
# perfection as painless as possible. Also they try to sort of replicate R, not
# exactly the same demographic as us -- just want easy statistics plots where you
# pass a dataset and everything else is abstracted away. I just need to make
# a simple violin plot/split violin plot API, a 1-D and 2-D kernel density API,
# a histogram API, and a boxplot API. Or, actually, could steal Seaborn source
# code and just make my own functions.
#------------------------------------------------------------------------------#
# First just make sure some dependencies are loaded
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
from glob import glob
from cycler import cycler
from functools import wraps
from . import colortools
from . import utils
# Will add our own dictionary to the top-level matplotlib module, to go
# alongside rcParams
import matplotlib as mpl
# Default settings to be loaded when globals() is
# called without any arguments
rcDefaults = {
    'linewidth': 0.7, # points
    'small':     8, # inches
    'medium':    8, # inches
    'large':     9, # inches
    'ticklen':   4, # points
    'tickpad':   2, # points (distance between ticks and labels)
    'color':     'k',
    'cycle':     'colorblind',
    'fontname':  'DejaVu Sans',
    }
# Set up dictionary of default parameters
base_context = {
    'font.size':         12,
    'axes.labelsize':    12,
    'axes.titlesize':    12,
    'xtick.labelsize':   11,
    'ytick.labelsize':   11,
    'legend.fontsize':   11,
    'axes.linewidth':    1.25,
    'grid.linewidth':    1,
    'lines.linewidth':   1.5,
    'lines.markersize':  6,
    'patch.linewidth':   1,
    'xtick.major.width': 1.25,
    'ytick.major.width': 1.25,
    'xtick.minor.width': 1,
    'ytick.minor.width': 1,
    'xtick.major.size':  6,
    'ytick.major.size':  6,
    'xtick.minor.size':  4,
    'ytick.minor.size':  4,
    }
# Name keys
context_keys = [*base_context.keys()]
font_keys = ['axes.labelsize',  'axes.titlesize',  'legend.fontsize',
             'xtick.labelsize', 'ytick.labelsize', 'font.size']

#------------------------------------------------------------------------------#
# Contextual settings management
# Adapted from seaborn; see: https://github.com/mwaskom/seaborn/blob/master/seaborn/rcmod.py
#------------------------------------------------------------------------------#
class PlottingContext(dict):
    """
    Plotting context object.
    """
    # Helper functions
    _keys = context_keys
    @staticmethod
    def _set(context=None, font_scale=1, rc=None):
        context_object = plotting_context(context, font_scale, rc)
        mpl.rcParams.update(context_object)
    # Necessary methods for creating 'with' block
    def __enter__(self):
        rc = mpl.rcParams
        self._orig = {k: rc[k] for k in self._keys}
        self._set(self)
    def __exit__(self, exc_type, exc_value, exc_tb):
        self._set(self._orig)
    def __call__(self, func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            with self:
                return func(*args, **kwargs)
        return wrapper

def plotting_context(context=None, font_scale=1, rc=None):
    """
    Returns plotting context object.
    """
    # Create context
    if context is None: # this is the PlottingContext instance
        context_dict = {k: mpl.rcParams[k] for k in context_keys}
    elif isinstance(context, dict):
        context_dict = context
    else:
        contexts = ['paper', 'notebook', 'talk', 'poster']
        if context not in contexts:
            raise ValueError(f'Context must be in {join(contexts)}.')
        # Scale all the parameters by the same factor depending on the context
        scaling = dict(paper=.8, notebook=1, talk=1.5, poster=2)[context]
        context_dict = {k:v*scaling for k, v in base_context.items()}
        # Now independently scale the fonts
        font_dict = {k: context_dict[k]*font_scale for k in font_keys}
        context_dict.update(font_dict)
    # Override these settings with the provided rc dictionary
    if rc is not None:
        rc = {k: v for k, v in rc.items() if k in context_keys}
        context_dict.update(rc)
    # Wrap in a PlottingContext object so this can be used in a with statement
    context_object = PlottingContext(context_dict)
    return context_object

#-------------------------------------------------------------------------------
# Settings management function
#-------------------------------------------------------------------------------
# I used to create this neat class but it was a stupid idea. Work with existing
# API, not against it. Just create a function that sets rcParam defaults with
# optional override. Then, better to make the format function add actual information
# to the plot and do nothing to change its style/color/etc.
# def rc(**kwargs): # fontname is matplotlib default
def rc(*args, verbose=False, **kwargs):
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
    See: https://matplotlib.org/users/customizing.html
    Here's a quick list of rcParam categories:
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
    #--------------------------------------------------------------------------#
    # Helper function; adds stuff to rcParams or rcExtras in the
    # stylesheet style of category.subcategory = value, or very
    # occasionally category.subcategory.subsubcategory = value, in which case
    # user should input a dictionary {'subcategory.subsubcategory':value}
    #--------------------------------------------------------------------------#
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
        # Get dictionary
        if category=='globals':
            dictionary = {key.split('.')[-1]:val for key,val in mpl.rcExtras.items() if 'globals.' in key}
        else:
            dictionary = {}
            for catstring,value in {**mpl.rcParams, **mpl.rcExtras}.items():
                if catstring.split('.')[0]==category:
                    dictionary[catstring.split('.')[1]] = value
        # Optionally filter results, or just return one value
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
        raise ValueError(f'Improper use of rc(). Only supply extra *args without any **kwargs.')
    #--------------------------------------------------------------------------#
    # *Initialize* default settings; that is, both rcParams and rcExtras
    #--------------------------------------------------------------------------#
    # Requires processing in lines below
    # WARNING: rcdefaults() changes the backend! Inline plotting will fail for
    # rest of notebook session if you call rcdefaults before drawing a figure!!!!
    # After first figure made, backend property is 'sticky', never changes!
    category = category or 'globals' # None defaults to globals in this section
    if category=='globals' and not kwargs: # default settings
        # See: https://stackoverflow.com/a/48322150/4970632
        # backend_save = mpl.rcParams['backend']
        # mpl.rcdefaults() # apply *builtin* default settings
        # mpl.rcParams['backend'] = backend_save
        mpl.style.use('default') # mpl.style does not change the backend
        mpl.rcExtras = {} # reset extra settings
        add('globals', rcDefaults) # apply *global* settings
    #--------------------------------------------------------------------------#
    # *Update* rcParam settings or global settings; just add them to rcParams or rcExtras
    #--------------------------------------------------------------------------#
    # If one of the 'global' settings was changed, requires processing in lines below
    # Also if any 'value' is a dictionary, update properties with that key as a
    # 'category' and each key-value air in the dictionary as kwargs
    else:
        updated = False
        for key,value in kwargs.items():
            if isinstance(value,dict):
                add(key, value) # add the 'value' as a kwarg dictionary
                if verbose:
                    print(f"Added settings to category '{key}':", value)
            elif category=='globals' and key not in rcDefaults:
                raise ValueError(f"Key \"{key}\" unknown. Not a global property like \"color\" or \"linewidth\".")
            else:
                value = value if value!='default' or category!='globals' else rcDefaults[key]
                add(category, {key:value}) # easy peasy
                updated = True
                if verbose:
                    print(f"Added setting to category '{category}':", {key:value})
        if category!='globals' or (category=='globals' and not updated):
            if verbose:
                print('No global properties changed.')
            return None # don't need to re-apply the globals
    #--------------------------------------------------------------------------#
    # *Apply* global settings; if this function was not called without arguments, then
    # only settings specifically requested to be changed, will be changed
    #--------------------------------------------------------------------------#
    current = rc('globals') # the dictionary
    # Make a cycler for drawing a bunch of lines
    # dashes = ('-', '--', ':', '-.') # dash cycles; these succeed color changes
    # cycle = cycler('color', [colors[i%len(colors)] for i in range(len(colors)*len(dashes))]) \
    #       + cycler('linestyle', [i for i in dashes for n in range(len(colors))])
    # [ax.set_prop_cycle(cycle) for fig in map(plt.figure, plt.get_fignums()) for ax in fig.axes]
    cycle = current['cycle']
    if not utils.isiterable(cycle) or type(cycle) is str:
        cycle = cycle,
    colors = colortools.Cycle(*cycle)
    cycle  = cycler('color', colors)
    # First the rcParam settings
    # Here are ones related to axes and figure properties
    # NOTE the figure and axes colors will not be reset on saving if
    # you specify them explicitly, even if you use transparent True.
    color, linewidth, ticklen, tickpad, small, large, fontname = \
        current['color'], current['linewidth'], current['ticklen'], current['tickpad'], current['small'], current['large'], current['fontname']
    add('savefig', {'transparent':True, 'facecolor':'w', 'dpi':300,
        'directory':'', 'pad_inches':0, 'bbox':'standard', 'format':'pdf'}) # empty means current directory
    add('figure', {'facecolor':(.95,.95,.95,1), 'dpi':90, # matches nbsetup
        'max_open_warning':0, 'constrained_layout':False, 'autolayout':False}) # zero disables max open warning
    add('axes', {'xmargin':0, 'ymargin':0.05, 'labelsize':small, 'titlesize':large,
        'edgecolor':color, 'labelcolor':color, 'grid':True, 'linewidth':linewidth,
        'labelpad':3, 'prop_cycle':cycle})
    add('font', {'size':small, 'family':'sans-serif', 'sans-serif':fontname}) # when user calls .text()
    add('text', {'color':color}) # when user calls .text()
    add('grid', {'linewidth':linewidth/2, 'alpha':0.1, 'color':color})
    for xy in 'xy':
        tickloc = {'x':{'bottom':True,'top':False},'y':{'left':True,'right':False}}[xy]
        add(f'{xy}tick',       {'labelsize':small, 'color':color, 'direction':'out'})
        add(f'{xy}tick.major', {'pad':tickpad, 'width':linewidth, 'size':ticklen, **tickloc}) # size==length
        add(f'{xy}tick.minor', {'pad':tickpad, 'width':linewidth, 'size':ticklen/2, **tickloc, 'visible':True})
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
        'labelspacing':0.5, 'handletextpad':0.5, 'handlelength':1.5, 'columnspacing':1,
        'borderpad':0.5,    'borderaxespad':0,   'numpoints':1,    'facecolor':'inherit'})
    # add('legend',   {'framealpha':0.6, 'fancybox':False, 'frameon':False,
    #     'labelspacing':0.5, 'columnspacing':1, 'handletextpad':0.5, 'borderpad':0.25})
    # This can only be accomplished with rcParams; impossible to specify with API!
    add('mathtext', {'default':'regular', 'bf':'sans:bold', 'it':'sans:it'}) # no italicization
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
    add('oceans',      {'color':'w', 'linewidth':0}) # make sure no lines!
    add('tickminor',   {'length':ticklen/2, 'width':linewidth, 'color':color})
    add('tick',        {'length':ticklen, 'width':linewidth, 'color':color})
    add('ctickminor',  {'length':ticklen/2, 'width':linewidth, 'color':color})
    add('ctick',       {'length':ticklen, 'width':linewidth, 'color':color})
    add('coastlines',  {'linewidth':linewidth, 'color':color})
    add('lonlatlines', {'linewidth':linewidth, 'linestyle':':', 'color':color, 'alpha':0.2})
    add('spine',       {'color':color, 'linewidth':linewidth})
    add('outline',     {'edgecolor':color, 'linewidth':linewidth})
    if verbose:
        print('Global properties set.')
    return None

# Now call the function to configure params with default values
rc(verbose=True)

