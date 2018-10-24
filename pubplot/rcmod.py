#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# This configures the global working environment
"""
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

Notes
-----
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

Todo
-----
* Change default background axes color to 0072B2, from seaborn colorblind, just
    a bit darker than figure background; is nice.
* Expand idea of default 'color', 'linewidth', whatever and create more
    dictionaries with defaults settings, then update them. No more hardcoded values
    in the main text below.
"""
#------------------------------------------------------------------------------#
# First just make sure some dependencies are loaded
import os
import re
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
# Default settings
# Roughly speaking, these tend to apply across multiple rc settings
rcDefaults = {
    'axes_color':          'k',
    'grid_color':          'k',
    'grid_alpha':          0.1,
    'patch_color':         'k',
    'patch_linewidth':     1.0,
    'figure_facecolor': (0.95,0.95,0.95,1),
    'axes_facecolor':   (1,1,1,1),
    'cmap':                'sunset',
    'cycle':               'colorblind',
    'fontname':            'DejaVu Sans',
    'fontsize_small':      8, # inches
    'fontsize_large':      9, # inches
    'axes_linewidth':      0.8, # points
    'grid_linewidth':      0.6, # a bit thinner
    'gridminor_linewidth': 0.4, # a bit thinner still
    'plot_linewidth':      1.3,
    'major_ticklen':       4.0,
    'minor_ticklen':       2.0, # points
    'tickpad' :            2.0, # points (distance between ticks and labels)
    'labelpad' :           3.0, # points (distance between ticks and labels)
    'markersize' :         3.0,
    }
# rcDefaults = {
#     'axes_color':          'b',
#     'grid_color':          'b',
#     'grid_alpha':          0.7,
#     'patch_color':         'green',
#     'cmap':                'Greys',
#     'cycle':               'cinematic1',
#     'fontname':            'DejaVu Sans',
#     'fontsize_small':      16, # inches
#     'fontsize_large':      25, # inches
#     'axes_linewidth':      2, # points
#     'grid_linewidth':      4, # a bit thinner
#     'gridminor_linewidth': 6, # a bit thinner still
#     'plot_linewidth':      8,
#     'patch_linewidth':     10,
#     'major_ticklen':       12,
#     'minor_ticklen':       14, # points
#     'tickpad' :            5, # points (distance between ticks and labels)
#     'labelpad' :           8, # points (distance between ticks and labels)
#     'markersize' :         8,
#     }

#------------------------------------------------------------------------------#
# Contextual settings management
# Adapted from seaborn; see: https://github.com/mwaskom/seaborn/blob/master/seaborn/rcmod.py
#------------------------------------------------------------------------------#
class rc_context(object):
    """
    Context object.
    """
    def __init__(self, *args, **kwargs):
        """
        Save user input.
        """
        self._args = args
        self._kwargs = kwargs
        self._rcSpecial = rc._rcSpecial.copy()
        self._rcParams = mpl.rcParams.copy()

    def __enter__(self):
        """
        Apply them.
        """
        for key,value in self._kwargs.items():
            rc.__setattr__(key, value)

    def __exit__(self, _type, _value, _traceback):
        """
        Return to previous state.
        """
        for key,value in self._rcParams.items():
            mpl.rcParams[key] = value
        for key,value in self._rcSpecial.items():
            rc._rcSpecial[key] = value

#-------------------------------------------------------------------------------
# Settings management
#-------------------------------------------------------------------------------
class rcConfigurator(object):
    """
    Abstract class for handling settings.
    """
    def __init__(self):
        """
        Initialize.

        Notes
        -----
        Note rcdefaults() changes the backend! Inline plotting will fail for
        rest of notebook session if you call rcdefaults before drawing a figure!
        After first figure made, backend property is 'sticky', never changes!
        See: https://stackoverflow.com/a/48322150/4970632
        """
        # First initialize matplotlib
        mpl.style.use('default') # mpl.style does not change the backend
        mpl.rcPubPlot = {} # reset extra settings
        # Next add global attributes
        self._settings = rcDefaults.copy()
        self._rcSpecial  = {}
        # Finally add these attributes to rcParams
        self._apply()

    def reset(self): # user can use this
        self.__init__()

    def __str__(self, item):
        print(self._settings)

    def __repr__(self, item):
        print(repr(self._settings))

    def _set_rcparam(self, category, kwargs):
        """
        Input
        -----
            category : the rcParams category
            kwargs : ditionary of subcategory-value pairsk

        Description
        -----------
        Adds values rcParams or rcPubPlot in the style:
            category.subcategory = value
        or very occasionally
            category.subcategory.subsubcategory = value
        in which case the 'category' should be 'category.subcategory'.
        """
        for subcategory,value in kwargs.items():
            key = category+'.'+subcategory
            if key in mpl.rcParams:
                mpl.rcParams[key] = value
            else:
                raise ValueError(f'Invalid rcParam: "{key}".')

    def _set_special(self, category, kwargs):
        """
        Input
        -----
            category : settings

        Description
        -----------
        Add stuff to special dictionary.
        """
        for subcategory,value in kwargs.items():
            key = category+'.'+subcategory
            self._rcSpecial[key] = value

    def __getitem__(self, item):
        """
        Gets arbitrary stuff.
        """
        # Get dictionary
        params = {}
        for d in (mpl.rcParams, self._rcSpecial):
            for category,value in d.items():
                if re.match(f'^{item}$', category):
                    return mpl.rcParams[item]
                elif re.search(f'^{item}\.', category):
                    subcategory = re.sub(f'^{item}.', '', category)
                    if subcategory and '.' not in subcategory:
                        params[subcategory] = d[category]
        if not params:
            raise ValueError(f'Invalid key "{item}".')
        return params

    def __setitem__(self, item, kwargs):
        """
        Apply to rcParams.
        """
        if isinstance(kwargs, dict):
            updated = False
            for key,value in kwargs.items():
                key = f'{item}.{key}'
                if key in rcParams:
                    rcParams[key] = item
                elif key in self._rcSpecial:
                    self._rcSpecial[key] = item
                else:
                    raise ValueError(f'Invalid key "{key}" for parameter "{item}".')

    def __getattr__(self, attr):
        """
        Get an rc setting.
        """
        if attr.startswith('_') or attr=='reset':
            return super().__getattr__(attr)
        elif attr not in rcDefaults:
            raise ValueError(f'Invalid global rcparam. Options are: {", ".join(rcDefaults.keys())}.')
        else:
            return self._settings[attr]

    def __setattr__(self, attr, value):
        """
        Set an rc setting.
        """
        if attr.startswith('_') or attr=='reset':
            super().__setattr__(attr, value)
        elif attr not in rcDefaults:
            print(attr)
            raise ValueError(f'Invalid global rcparam. Options are: {", ".join(rcDefaults.keys())}.')
        else:
            self._settings[attr] = value
            self._apply()

    def _apply(self, verbose=False):
        """
        Apply global settings.

        Note
        ----
        The figure colors will not be reset on saving if you specify them
        explicitly, even if you use transparent True! Only works if they
        remain 'None' until saving.
        """
        # Figure settings
        self._set_rcparam('savefig', dict(transparent=True, facecolor=self.axes_facecolor,
            dpi=300, pad_inches=0, directory='', # empty means current directory
            bbox='standard', format='pdf'))
        self._set_rcparam('figure', dict(facecolor=self.figure_facecolor, dpi=90,
            titlesize=self.fontsize_large, titleweight='bold', # for suptitle
            max_open_warning=0, # zero disables max open warning
            autolayout=False))

        # Axes settings (note the gridminor one is special)
        # First get cycler, then apply
        cycle = self.cycle
        if utils.isscalar(cycle):
            cycle = cycle,
        colors = colortools.Cycle(*cycle)
        prop_cycle = cycler('color', colors)
        self._set_special('gridminor', dict(linestyle='-', linewidth=self.gridminor_linewidth,
            color=self.grid_color, alpha=self.grid_alpha))
        self._set_rcparam('grid', dict(linestyle='-', linewidth=self.grid_linewidth,
            color=self.grid_color, alpha=self.grid_alpha))
        self._set_special('abc',  dict(size=self.fontsize_large, weight='bold',
            color=self.axes_color))
        self._set_rcparam('axes', dict(xmargin=0, ymargin=0.05,
            titlesize=self.fontsize_large, titleweight='normal',
            labelsize=self.fontsize_small, labelweight='bold', labelcolor=self.axes_color,
            linewidth=self.axes_linewidth, edgecolor=self.axes_color,
            labelpad=self.labelpad, grid=True,
            prop_cycle=prop_cycle)) # consider adding linestyle cycle

        # Text settings
        # The mathtext stuff can only be accomplished with rcParams; impossible
        # to specify with API!
        self._set_rcparam('font', dict({'sans-serif':self.fontname}, size=self.fontsize_small, family='sans-serif'))
        self._set_rcparam('text', dict(color=self.axes_color)) # when user calls .text()
        self._set_rcparam('mathtext', dict(default='regular', bf='sans:bold', it='sans:it')) # no italicization

        # Tick marks
        tickloc = {'x':dict(bottom=True, top=False), 'y':dict(left=True,right=False)}
        for xy in 'xy':
            self._set_rcparam(xy+'tick', dict(labelsize=self.fontsize_small,
                color=self.axes_color, direction='out'))
            self._set_rcparam(xy+'tick.minor', dict(tickloc[xy], visible=True,
                pad=self.tickpad,
                width=self.gridminor_linewidth, size=self.minor_ticklen))
            self._set_rcparam(xy+'tick.major', dict(tickloc[xy],
                pad=self.tickpad,
                width=self.grid_linewidth, size=self.major_ticklen))

        # Hatching and colorsmaps and stuff
        self._set_rcparam('image', dict(cmap=self.cmap, lut=256)) # colormap stuff
        self._set_rcparam('patch', dict(linewidth=self.patch_linewidth, edgecolor=self.patch_color))
        self._set_rcparam('hatch', dict(linewidth=self.patch_linewidth, color=self.patch_color))

        # Lines and markers
        # For join/cap styles, see: https://matplotlib.org/examples/api/joinstyle.html
        self._set_rcparam('markers', dict(fillstyle='full'))
        self._set_rcparam('scatter', dict(marker='o'))
        self._set_rcparam('lines', dict(linewidth=self.plot_linewidth, color='C0',
            markeredgewidth=0, markersize=self.markersize,
            dash_joinstyle='miter',  dash_capstyle='projecting',   # joinstyle opts= miter, round, bevel
            solid_joinstyle='miter', solid_capstyle='projecting')) # capstyle opts= butt, round, projecting

        # Legend settings
        # See: https://matplotlib.org/api/legend_api.html
        self._set_rcparam('legend',   dict(framealpha=1, fancybox=False, frameon=False,
            fontsize=self.fontsize_small,
            labelspacing=0.5, handletextpad=0.5,   handlelength=1.5,
            columnspacing=1,  facecolor='inherit', numpoints=1,
            borderpad=0.5,    borderaxespad=0))

        # Geographic stuff (these are custom)
        self._set_special('coastlines',  dict(linewidth=self.axes_linewidth, color=self.axes_color))
        self._set_special('land',        dict(linewidth=0, color=self.axes_color)) # no lines!
        self._set_special('oceans',      dict(linewidth=0, color='w')) # no lines!
        self._set_special('lonlatlines', dict(linewidth=self.grid_linewidth, color=self.axes_color,
            linestyle='=', alpha=self.grid_alpha))
        if verbose:
            print('Updated global settings.')

# Instantiate object
rc = rcConfigurator()

