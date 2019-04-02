#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# This configures the global working environment
"""
See: https://matplotlib.org/users/customizing.html
Here's a quick list of rcParam categories:
    "lines", "patch", "hatch", "legend"
    "font", "text", "mathtext"
    "axes", "figure"
    "date", "xtick", "ytick", "grid"
    "contour", "image"
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
"""
#------------------------------------------------------------------------------#
# First just make sure some dependencies are loaded
import re
from matplotlib.pyplot import figure, get_fignums
from cycler import cycler
from . import colortools
from . import utils
from .utils import timer, counter, ic
from matplotlib import rcParams, style
# Get default font
# WARNING: Had issues with Helvetica Neue on Linux, weirdly some characters
# failed to render/printed nonsense, but Helvetica fine
import sys
_default_font = 'Helvetica' if sys.platform=='linux' else 'Helvetica Neue' # says 'darwin' on mac
# Will add our own dictionary to the top-level matplotlib module, to go
# alongside rcParams
# Default settings
# List of linked settings
rcGlobals = {
    # Apply these ones to list of rcParams
    'color':     'k',
    'xcolor':    None, # these are special; can be used to set particular spine colors
    'ycolor':    None,
    'cycle':     'colorblind',
    # 'facecolor':  '#0072b2', # 0072B2
    'facecolor':  'w', # 0072B2
    'facehatch':  None, # hatching on background, useful for indicating invalid data
    'gridalpha':  0.1,
    'small':      8,
    'large':      9,
    'linewidth':  0.6,
    'gridwidth':  0.6,
    'bottom':     True,
    'top':        False,
    'left':       True,
    'right':      False,
    'ticklen':    4.0,
    'tickpad':    2.0,
    'tickdir' :   'out',
    # Convenient aliases (i.e. they do not bulk apply to a bunch of settings, just shorter names)
    'fontname':       _default_font, # best one; and less crammed than Helvetica
    'margin':         0.0,
    'xmargin':        0.0, # found I wanted to change these *a lot*
    'ymargin':        0.0, # found I wanted to change these *a lot*
    'abcweight':      'bold',
    'titleweight':    'normal',
    'suptitleweight': 'bold',
    # Special ones
    'tickratio':  0.5, # ratio of major-to-minor tick size
    'minorwidth': 0.8, # ratio of major-to-minor tick width
    'gridratio':  0.5, # ratio of major-to-minor grid line widths
    }
rcGlobals_children = {
    # Most important ones, expect these to be used a lot
    # The xcolor/ycolor we don't use 'special' props (since we'd be duplicating ones
    # that already exist for all spines/labels). Instead just manually use the
    # global property in the format script.
    # NOTE: Hatches are weird: https://stackoverflow.com/questions/29549530/how-to-change-the-linewidth-of-hatch-in-matplotlib
    # Linewidths can only be controlled with a global property!
    # NOTE: Should I even bother setting these? Yes: Idea is maybe we change
    # underlying global keywords, but have rcupdate always refer to
    # corresponding builtin values.
    'xcolor':     [],
    'ycolor':     [],
    'color':      ['axes.labelcolor', 'axes.edgecolor', 'axes.hatchcolor', 'map.color', 'map.hatchcolor', 'xtick.color', 'ytick.color'], # change the 'color' of an axes
    'facecolor':  ['axes.facecolor', 'map.facecolor'], # simple alias
    'facehatch':  ['axes.facehatch', 'map.facehatch'], # optionally apply background hatching
    'small':      ['font.size', 'xtick.labelsize', 'ytick.labelsize', 'axes.labelsize', 'legend.fontsize'], # the 'small' fonts
    'large':      ['abc.fontsize', 'figure.titlesize', 'axes.titlesize'], # the 'large' fonts
    'linewidth':  ['axes.linewidth', 'map.linewidth', 'hatch.linewidth', 'axes.hatchlw',
                   'map.hatchlw', 'xtick.major.width', 'ytick.major.width'], # gridline widths same as tick widths
                   # 'grid.linewidth', # should not be coupled, looks ugly
    'gridalpha':  ['grid.alpha',     'gridminor.alpha'],
    'gridcolor':  ['grid.color',     'gridminor.color'],
    'gridstyle':  ['grid.linestyle', 'gridminor.linestyle'],
    # Aliases
    'margin':         ['axes.xmargin', 'axes.ymargin'],
    'xmargin':        ['axes.xmargin'], # found I wanted to change these *a lot*
    'ymargin':        ['axes.ymargin'], # found I wanted to change these *a lot*
    'fontname':       ['font.family'], # specify family directly, so we can easily switch between serif/sans-serif; requires text.usetex = False; see below
    'abcweight':      ['abc.weight'],
    'titleweight':    ['axes.titleweight'],
    'suptitleweight': ['figure.titleweight'],
    # Less important ones
    'bottom':     ['xtick.major.bottom',  'xtick.minor.bottom'], # major and minor ticks should always be in the same place
    'top':        ['xtick.major.top',     'xtick.minor.top'],
    'left':       ['ytick.major.left',    'ytick.minor.left'],
    'right':      ['ytick.major.right',   'ytick.minor.right'],
    'ticklen' :   ['xtick.major.size',    'ytick.major.size'],
    'tickdir':    ['xtick.direction',     'ytick.direction'],
    'tickpad':    ['xtick.major.pad', 'xtick.minor.pad', 'ytick.major.pad', 'ytick.minor.pad'],
    }
# Settings that apply to just one thing, and are
# already implemented by matplotlib
rcDefaults = {
    # Some of these will be overwritten by global alises
    # Figure stuff
    'figure.dpi':              90, # save ipython notebook space
    'figure.facecolor':        (0.95,0.95,0.95,1),
    'figure.max_open_warning': 0,
    'figure.autolayout':       False,
    'figure.titleweight':      'bold',
    'savefig.facecolor':       (1,1,1,1),
    'savefig.transparent':     True,
    'savefig.dpi':             300,
    'savefig.pad_inches':      0.0,
    'savefig.directory':       '',
    'savefig.bbox':            'standard',
    'savefig.format':          'pdf',
    # Axes elements
    'axes.xmargin':            0.0,
    'axes.ymargin':            0.0,
    'axes.titleweight':        'normal',
    'axes.grid':               True,
    'axes.labelweight':        'normal',
    'axes.labelpad':           3.0,
    'axes.titlepad':           3.0,
    'axes.axisbelow':          'lines', # for ticks/gridlines *above* patches, *below* lines, use 'lines'
    'xtick.minor.visible' :    True,
    'ytick.minor.visible' :    True,
    'grid.color':              'k',
    'grid.alpha':              0.1,
    'grid.linestyle':          '-',
    'grid.linewidth':          0.6, # a bit thinner
    # Plot elements
    'image.cmap':              'sunset',
    'image.lut':               256,
    'patch.facecolor':         'C0',
    'patch.edgecolor':         'k',
    'patch.linewidth':         1.0,
    'hatch.color':             'k',
    'hatch.linewidth':         0.7,
    'markers.fillstyle':       'full',
    'scatter.marker':          'o',
    'lines.linewidth' :        1.3,
    'lines.color' :            'C0',
    'lines.markeredgewidth' :  0,
    'lines.markersize' :       3.0,
    'lines.dash_joinstyle' :   'miter',
    'lines.dash_capstyle' :    'projecting',
    'lines.solid_joinstyle' :  'miter', # joinstyle opts= miter, round, bevel
    'lines.solid_capstyle' :   'projecting', # capstyle opts= butt, round, projecting
    'legend.fancybox' :        False,
    'legend.frameon' :         False,
    'legend.labelspacing' :    0.5,
    'legend.handletextpad' :   0.5,
    'legend.handlelength' :    1.5,
    'legend.columnspacing' :   1,
    'legend.facecolor' :       'w',
    'legend.numpoints' :       1,
    'legend.borderpad' :       0.5,
    'legend.borderaxespad' :   0,
    # Text
    # Tried enabling multi-color from: https://stackoverflow.com/a/42768093/4970632, did not work
    # 'text.usetex': True,      # use LaTeX to write all text
    # 'pgf.rcfonts': False,     # Ignore Matplotlibrc
    # 'pgf.preamble': [ r'\usepackage{color}' ], # xcolor for colours
    # 'font.family': 'DejaVu Sans', # allowed to be concrete name(s) when usetex is False
    # 'font.family':       'sans-serif',
    # 'font.sans-serif':   'DejaVu Sans',
    'text.latex.preamble': r'\usepackage{cmbright}', # https://stackoverflow.com/a/16345065/4970632
    'text.usetex':         False, # use TeX for *all* font handling (limits available fonts)
    'mathtext.default':    'regular', # no italics
    'mathtext.bf' :        'sans:bold',
    'mathtext.it' :        'sans:it',
    }
# Special settings, should be thought of as extension of rcParams
rcDefaults_sp = {
    # These ones just need to be present, will get reset by globals
    'map.facecolor':       None,
    'map.color':           None,
    'map.linewidth':       None,
    'abc.fontsize':        None,
    'rowlabel.fontsize':   None,
    'collabel.fontsize':   None,
    'gridminor.alpha':     None,
    'axes.facehatch':      None,
    'axes.hatchcolor':     None,
    'axes.hatchlw':        None,
    'map.facehatch':       None,
    'map.hatchcolor':      None,
    'map.hatchlw' :        None,
    # The rest can be applied as-is
    'abc.weight':            'bold',
    'abc.color':             'k',
    'rowlabel.weight':       'bold',
    'rowlabel.color':        'k',
    'collabel.weight':       'bold',
    'collabel.color':        'k',
    'gridminor.color':       'k',
    'gridminor.linestyle':   '-',
    'gridminor.linewidth':   0.1,
    'land.linewidth':        0, # no boundary for patch object
    'land.color':            'k',
    'ocean.linewidth':       0, # no boundary for patch object
    'ocean.color':           'w',
    'coastline.linewidth' :  1.0,
    'coastline.color' :      'k',
    'lonlatlines.linewidth': 1.0,
    'lonlatlines.linestyle': ':',
    # 'lonlatlines.linestyle': '--',
    'lonlatlines.alpha':     0.4,
    'lonlatlines.color':     'k',
    'gridspec.title':        0.2, # extra space for title/suptitle
    'gridspec.inner':        0.2, # just have ticks, no labeels
    'gridspec.legend':       0.25, # default legend space (bottom of figure)
    'gridspec.cbar':         0.17, # default colorbar width
    'gridspec.ylab':         0.7, # default space wherever we expect tick and axis labels (a bit large if axis has no negative numbers/minus sign tick labels)
    'gridspec.xlab':         0.55, # for horizontal text should have more space
    'gridspec.nolab':        0.15, # only ticks
    }
rcParams_sp = rcDefaults_sp.copy()
n_sp = len(rcParams_sp) # instead of fancy rcGlobals-style key verification, just check if dictionary has grown when user sets an item -- in this case, property was invalid, so we throw error
# Generate list of valid names, and names with subcategories
rc_names = {
    *rcParams.keys(),
    *rcParams_sp.keys(),
    }
rc_categories = {
    *(re.sub('\.[^.]*$', '', name) for name in rc_names),
    *(re.sub('\..*$', '', name) for name in rc_names)
    }

#-------------------------------------------------------------------------------
# Contextual settings management
# Adapted from seaborn; see: https://github.com/mwaskom/seaborn/blob/master/seaborn/rcmod.py
#-------------------------------------------------------------------------------
class AttributeDict(dict):
    # Dictionary elements are attributes too.
    def __getattr__(self, attr): # invoked only if __getattribute__ fails
        return self[attr]
    def __setattr__(self, attr, value):
        self[attr] = value

class rc_configurator(object):
    _public_api = ('reset', 'update', 'fill') # getattr and setattr will not look for these items on underlying dictionary
    def __init__(self):
        """
        Magical abstract class for handling custom settings, builtin rcParams
        settings, and artificial 'global' params that keep certain groups
        of settings synced. Also includes context manager.
        """
        # First initialize matplotlib
        # Note rcdefaults() changes the backend! Inline plotting will fail for
        # rest of notebook session if you call rcdefaults before drawing a figure!
        # After first figure made, backend property is 'sticky', never changes!
        # See: https://stackoverflow.com/a/48322150/4970632
        style.use('default') # mpl.style function does not change the backend
        # Add simple attributes to rcParams
        self._rcCache = {}
        self._rcGlobals = rcGlobals.copy()
        for key,value in rcDefaults.items():
            rcParams[key] = value
        for key,value in rcDefaults_sp.items():
            rcParams_sp[key] = value
        # Apply linked attributes to rcParams
        self._set_cycler('colorblind')
        rc, rc_sp = self._get_globals()
        rcParams.update(rc)
        rcParams_sp.update(rc_sp)
        # Settings
        # TODO: Looks like _getitem_mode can get stuck on a higher, more
        # restrictive value (e.g. 1 or 2) when cell fails to execute. Should
        # consider improving this.
        self._init = True
        self._cache_orig = {}
        self._cache_added = {}
        self._getitem_mode = 0

    def __enter__(self):
        # Apply new settings (will get added to _rcCache)
        for key,value in self._cache_added.items():
            self[key] = value # applies globally linked and individual settings

    def __exit__(self, _type, _value, _traceback):
        # Restore configurator cache to its previous state.
        self._rcCache = self._cache_orig
        self._cache_orig = {}
        self._cache_added = {}
        self._getitem_mode = 0

    # @counter
    def __getitem__(self, key):
        # Can get a whole bunch of different things
        # Get full dictionary e.g. for rc[None]
        if not key:
            return {**rcParams, **rcParams_sp}
        # Allow for special time-saving modes where we *ignore rcParams*
        # or even *ignore rcParams_sp*.
        mode = self._getitem_mode
        if mode==0:
            kws = (self._rcCache, rcParams_sp, rcParams)
        elif mode==1:
            kws = (self._rcCache, rcParams_sp)
        elif mode==2:
            kws = (self._rcCache,)
        else:
            raise ValueError(f'Invalid _getitem_mode {mode}.')
        # If it is available, return the values corresponding to names in
        # user dictionary; e.g. {'color':'axes.facecolor'} becomes {'color':'w'}
        # NOTE: Got weird bugs here. Dunno why. Use self.fill method instead.
        if key in rc_categories:
            params = {}
            for kw in kws:
                for category,value in kw.items():
                    if re.search(f'^{key}\.', category):
                        subcategory = re.sub(f'^{key}\.', '', category)
                        if subcategory and '.' not in subcategory:
                            params[subcategory] = value
            if mode==0 and not params:
                raise ValueError(f'Invalid category "{key}".')
            else:
                return params
        # Get individual property. Will successively index a few different dicts
        # Try to return the value
        for kw in (*kws[:1], self._rcGlobals, *kws[1:]):
            try:
                return kw[key]
            except KeyError:
                continue
        # If we were in one of the exlusive modes, return None
        if mode==0:
            raise ValueError(f'Invalid prop name "{key}".')
        else:
            return None

    # @counter
    def __setitem__(self, key, value):
        # First the special cycler
        # NOTE: No matter the 'setitem mode' this will always set the axes
        # prop_cycle rc settings
        if key=='cycle':
            self._set_cycler(value)
            self._rcCache['cycle'] = value
        # Apply global settings
        elif key in rcGlobals:
            if value=='default':
                value = rcGlobals[key]
            rc, rc_sp = self._get_globals(key, value)
            self._rcCache.update(rc)
            self._rcCache.update(rc_sp)
            self._rcCache[key] = value # also update cached global property itself
            self._rcGlobals[key] = value
            rcParams.update(rc)
            rcParams_sp.update(rc_sp)
        # Directly modify single parameter
        # NOTE: If 'setitem mode' is 0, this means user has directly set
        # something (we are not in a with..as context in format()), so we
        # want to directly modify rcParams.
        elif key in rc_names:
            self._rcCache[key] = value
            try:
                rcParams[key] = value
            except KeyError:
                pass
        else:
            raise ValueError(f'Invalid key "{key}".')
        self._init = False # no longer in initial state

    def __getattribute__(self, attr):
        # Alias to getitem
        if attr[:1]=='_' or attr in self._public_api: # no recursion since second comparison won't be evaluated if first comparison evaluates True
            return super().__getattribute__(attr)
        else:
            return self.__getitem__(attr)

    def __setattr__(self, attr, value):
        # Alias to setitem
        if attr[:1]=='_' or attr in self._public_api:
            super().__setattr__(attr, value)
        else:
            self.__setitem__(attr, value)

    def __str__(self):
        # Alias to __repr__
        return self.__repr__()

    def __repr__(self):
        # Nice string representation
        length = 1 + max(len(key) for key in self._rcGlobals.keys())
        string = '\n'.join(f'{key}: {" "*(length-len(key))}{value}'
                                    for key,value in self._rcGlobals.items())
        return string

    def _set_cycler(self, value):
        # Set the color cycler.
        # NOTE: Generally if user uses 'C0', et cetera, assume they want to
        # refer to the *default* cycler colors; so first reset
        if isinstance(value, str) or utils.isnumber(value):
            value = value,
        colors = colortools.colors('colorblind')
        rcParams['axes.prop_cycle'] = cycler('color', colors)
        colors = colortools.colors(*value)
        rcParams['axes.prop_cycle'] = cycler('color', colors)
        figs = list(map(figure, get_fignums()))
        for fig in figs:
            for ax in fig.axes:
                ax.set_prop_cycle(cycler('color', colors))

    def _get_globals(self, key=None, value=None):
        # Apply all properties in some group.
        kw = {}
        kw_sp = {}
        if key is not None and value is not None:
            items = [(key,value)]
        else:
            items = self._rcGlobals.items()
        for key,value in items:
            # Tick length/major-minor tick length ratio
            if key in ('ticklen','tickratio'):
                if key=='tickratio':
                    ticklen = self._rcGlobals['ticklen']
                    ratio = value
                else:
                    ticklen = value
                    ratio = self._rcGlobals['tickratio']
                kw['xtick.minor.size'] = ticklen*ratio
                kw['ytick.minor.size'] = ticklen*ratio
            # Spine width/major-minor tick width ratio
            if key in ('linewidth','minorwidth'):
                if key=='linewidth':
                    tickwidth = value
                    ratio = self._rcGlobals['minorwidth']
                else:
                    tickwidth = self._rcGlobals['linewidth']
                    ratio = value
                kw['xtick.minor.width'] = tickwidth*ratio
                kw['ytick.minor.width'] = tickwidth*ratio
                # kw_sp['gridminor.linewidth'] = tickwidth*ratio # special
            # Grid line
            if key in ('gridwidth', 'gridratio'):
                if key=='gridwidth':
                    gridwidth = value
                    ratio = self._rcGlobals['gridratio']
                else:
                    gridwidth = self._rcGlobals['gridwidth']
                    ratio = value
                kw_sp['gridminor.linewidth'] = gridwidth*ratio
            # Now update linked settings
            for name in rcGlobals_children.get(key, []):
                if name in rcParams_sp:
                    kw_sp[name] = value
                else:
                    kw[name] = value
        return kw, kw_sp

    def _context(self, *args, mode=0, **kwargs):
        """
        Temporarily modify rc configuration. Do this by simply
        saving the cache, allowing modification of the cache, then
        restoring the old cache.
        Three modes:
            0) __getitem__ searches everything, the default.
            1) __getitem__ ignores rcParams (assumption is these have already
               been set). Used during Axes __init__ calls to _rcupdate.
            2) __getitem__ ignores rcParams and rcParams_sp; only read from
                cache, i.e. settings that user has manually changed.
               Used during Axes format() calls to _rcupdate.
        Notes
        -----
        This is kept private, because it's only mean to be used within
        the 'format()' method automatically! Instead of user having to use
        with..as, they should just pass rc_kw dict or kwargs to format().
        """
        # Apply mode
        if mode not in range(3):
            raise ValueError(f'Invalid _getitem_mode {mode}.')
        for arg in args:
            if not isinstance(arg, dict):
                raise ValueError('rc_context() only accepts dictionary args and kwarg pairs.')
            kwargs.update(arg)
        self._getitem_mode = mode
        self._cache_orig   = rc._rcCache.copy()
        self._cache_added  = kwargs # could be empty
        return self

    def update(self, *args, **kwargs):
        """
        Same as setting rc['axes'] = {'name':value}, but do this for bunch
        of different propts.
        """
        if len(args)==0:
            args = [{}]
        kw = args[-1]
        kw.update(kwargs)
        if len(args)==1:
            for key,value in kw.items():
                self[key] = value
        elif len(args)==2:
            category = args[0]
            for key,value in kw.items():
                self[category + '.' + key] = value
        else:
            raise ValueError('rc.update() accepts 1-2 positional arguments. Use rc.update(kw) to update a bunch of names, or rc.update(category, kw) to update subcategories belonging to single category e.g. axes. All kwargs will be added to the dict.')

    def fill(self, props):
        """
        Function that only updates a property if self.__getitem__ returns not None.
        Meant for optimization; hundreds of 200-item dictionary lookups over several
        subplots end up taking toll, almost 1s runtime.
        """
        props_out = {}
        for key,value in props.items():
            value = self[value]
            if value is not None:
                props_out[key] = value
        return props_out

    def reset(self):
        """
        Restore settings to default.
        """
        return self.__init__()

# Instantiate object
rc = rc_configurator()

