#!/usr/bin/env python3
"""
This module manges global "rc" settings with a bunch of new convenience
features, and defines some brand new settings. There are now three different
types of settings:

* **Builtin** :ref:`rcParams` settings. These have the format
  `category.subcategory`, so approaches 1 and 2 are invalid (see below).
* **Custom** :ref:`rcSpecial` settings. These also have the format
  `category.subcategory`.
* **Global** :ref:`rcGlobal` settings. These are simple, short names
  used to change **multiple** "builtin" and "custom" settings at once,
  or as shorthands for settings with longer names.

Your one-stop-shop for changing settings is the `rc` object. To change
the setting named ``name`` to ``value``, use any of the following 4
approaches:

1. `rc.name = value`
2. `rc.update(name=value)`
3. `rc['name'] = value`
4. `rc.update({'name':value})`

The "builtin", "custom", and "global" settings are described in detail
below.

rcParams
--------
These are the builtin matplotlib settings. See `this page
<https://matplotlib.org/users/customizing.html>`_ for more info.
The `rcParams` categories are as follows:

* Axes: ``axes``.
* Text: ``font``, ``text``, ``mathtext``.
* Plot elements: ``lines``, ``patch``, ``hatch``, ``legend``, ``contour``, ``image``,
  ``boxplot``, ``errorbar``, ``hist``, ``scatter``, ``animation``.
* Axis elements: ``date``, ``xtick``, ``ytick``, ``grid``.
* Printing and saving: ``path``, ``figure``, ``savefig``, ``ps``, ``tk``, ``pdf``, ``svg``.
* Other: ``keymap``, ``examples``, ``debug``.

.. ``figure`` settings are used when printing interactively or
   just making the figure object, while ``savefig`` settings are used when calling
   `~matplotlib.figure.Figure.savefig`.

rcSpecial
---------
My brand new settings, meant to configure special ProPlot featues. The
`rcSpecial` categories are as follows:

* Background hatching: ``axes``, ``map``.
* Subplots: ``gridspec``.
* New labels: ``abc``, ``rowlabel``, ``collabel``
* Gridlines: ``gridminor``, ``lonlatlines``
* Geographic features: ``land``, ``ocean``, ``coastline``.

The ``map`` and ``axes`` subcategories:

==============  ==================================================================
Key             Description
==============  ==================================================================
``facehatch``   Background hatching string pattern, if not ``None`` [1]_.
``hatchcolor``  Color of background hatching pattern. Default is same as spines.
``hatchlw``     Line width for background hatching. Default is same as spines.
==============  ==================================================================

The ``gridspec`` subcategories (all values are in inches):

==========  ==================================================================
Key         Description
==========  ==================================================================
``title``   Vertical space for titles.
``legend``  Width of "legend" panels.
``cbar``    Width of "colorbar" panels.
``ylab``    Horizontal space between subplots alotted for *y*-labels.
``xlab``    Vertical space between subplots alotted for *x*-labels.
``nolab``   Space between subplots alotted for tick marks.
``inner``   Totally empty space between subplots.
==========  ==================================================================

The ``abc``, ``rowlabel``, and ``collabel`` subcategories:

============  ==================================================================
Key           Description
============  ==================================================================
``fontsize``  The font size.
``weight``    The font weight.
``color``     The font color.
============  ==================================================================

The ``lonlatlines`` and ``gridminor`` subcategories:

=============  ==================================================================
Key            Description
=============  ==================================================================
``linewidth``  The line width.
``linestyle``  The line style.
``alpha``      The line transparency.
``color``      The line color.
=============  ==================================================================

The ``land``, ``ocean``, and ``coastlines`` subcategories:

=============  ==================================================================
Key            Description
=============  ==================================================================
``linewidth``  The line width or patch edge width.
``color``      The line color or patch color.
=============  ==================================================================

rcGlobal
--------
Global settings are used to change :ref:`rcParams` and :ref:`rcSpecial` settings
in bulk, or as shorthands for common settings with longer names.

==================  ==============================================================================================================
Key                 Description
==================  ==============================================================================================================
``cycle``           The default color cycle name, used e.g. for lines.
``color``           The color of axis spines, tick marks, tick labels, and labels.
``xcolor``          As with `color`, but specific to *x*-axes.
``ycolor``          As with `color`, but specific to *y*-axes.
``margin``          The margin of space around subplot `~matplotlib.artist.Artist` instances, if ``xlim`` and ``ylim`` are unset.
``xmargin``         As with `margin`, but specific to the *x* direction.
``ymargin``         As with `margin`, but specific to the *y* direction.
``facecolor``       The axes background color.
``facehatch``       The background hatching pattern [1]_. Useful for highlighting "NaN" data in a ``contourf`` plot.
``small``           Font size for legend text, tick labels, axis labels, and text generated with `~proplot.axes.BaseAxes.text`.
``large``           Font size for titles, "super" titles, and a-b-c subplot labels.
``fontname``        Name of font used for all text in the figure [2]_.
``linewidth``       Thickness of axes spines and major tick lines.
``minorwidth``      Ratio of minor to major tick line thickness.
``gridwidth``       Thickness of gridlines.
``gridminor``       Ratio of minor to major gridline thickness.
``gridalpha``       Transparency of major and minor gridlines.
``gridcolor``       Color of major and minor gridlines.
``gridstyle``       Linestyle of major and minor gridlines.
``ticklen``         Length of major ticks.
``tickratio``       Ratio of minor to major tick lengths.
``tickdir``         Major and minor tick direction; one of ``out``, ``in``, or ``inout``.
``abcweight``       Font weight for a-b-c labels [3]_.
``titleweight``     Font weight for titles [3]_.
``suptitleweight``  Font weight for "super" titles [3]_.
==================  ==============================================================================================================

.. [1] For example, ``'xxx'`` or ``'..'``. See `this demo <https://matplotlib.org/gallery/shapes_and_collections/hatch_demo.html>`__.
.. [2] This module *changes the default* from DejaVu Sans (or Bitstream Vera) to
       Helvetica Neue or Helvetica. Install these with `~proplot.fonttools.install_fonts`
       after installing ProPlot for the first time, and after updating
       matplotlib.
.. [3] Valid weights are ``'ultralight'``, ``'light'``, ``'normal'``,
       ``'medium'``, ``'demi'``, ``'bold'``, ``'very bold'``, or ``'black'``.
       Note that many fonts only have ``normal`` or ``bold`` available.
       If you request another weight, the “closest” availble weight will
       be selected.

Todo
----
Disable the ``gridspec`` stuff by automating the inter-subplot spacing?

"""
# First just make sure some dependencies are loaded
import re
import cycler
import matplotlib.pyplot as plt
from . import colortools
from . import utils
from .utils import _timer, _counter, ic
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
_rcGlobal = {
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

_rcGlobal_children = {
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
_rcDefaults = {
    # Some of these will be overwritten by global alises
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
    'font.family':             'DejaVu Sans', # allowed to be concrete name(s) when usetex is False
    # 'font.family':             'sans-serif',
    # 'font.sans-serif':         'DejaVu Sans',
    'text.latex.preamble':      r'\usepackage{cmbright}', # https://stackoverflow.com/a/16345065/4970632
    'text.usetex':             False, # use TeX for *all* font handling (limits available fonts)
    'mathtext.default':        'regular', # no italics
    'mathtext.bf' :            'sans:bold',
    'mathtext.it' :            'sans:it',
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
    }

# Special settings, should be thought of as extension of rcParams
_rcDefaults_sp = {
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
    # the rest can be applied as-is
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
    'gridspec.inner':        0.2, # just have ticks, no labels
    'gridspec.legend':       0.25, # default legend space (bottom of figure)
    'gridspec.cbar':         0.17, # default colorbar width
    'gridspec.ylab':         0.7, # default space wherever we expect tick and axis labels (a bit large if axis has no negative numbers/minus sign tick labels)
    'gridspec.xlab':         0.55, # for horizontal text should have more space
    'gridspec.nolab':        0.15, # only ticks
    }
_rcParams_sp = _rcDefaults_sp.copy()

# Generate list of valid names, and names with subcategories
# TODO: Display this somewhere? Maybe in repr of rc?
_rc_names = {
    *rcParams.keys(),
    *_rcParams_sp.keys(),
    }
_rc_categories = {
    *(re.sub('\.[^.]*$', '', name) for name in _rc_names),
    *(re.sub('\..*$', '', name) for name in _rc_names)
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
        Magical abstract class for managing builtin `rcParams` settings, 
        our artificial `_rcSpecial` settings, and new "global" settings
        that keep certain groups of settings synced.

        See the `rc` documentation for details.
        """
        # First initialize matplotlib
        # Note rcdefaults() changes the backend! Inline plotting will fail for
        # rest of notebook session if you call rcdefaults before drawing a figure!
        # After first figure made, backend property is 'sticky', never changes!
        # See: https://stackoverflow.com/a/48322150/4970632
        style.use('default') # mpl.style function does not change the backend
        # Add simple attributes to rcParams
        self._rcCache = {}
        self._rcGlobal = _rcGlobal.copy()
        for key,value in _rcDefaults.items():
            rcParams[key] = value
        for key,value in _rcDefaults_sp.items():
            _rcParams_sp[key] = value
        # Apply linked attributes to rcParams
        self._set_cycler('colorblind')
        rc, rc_sp = self._get_globals()
        rcParams.update(rc)
        _rcParams_sp.update(rc_sp)
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

    def __getitem__(self, key):
        # Can get a whole bunch of different things
        # Get full dictionary e.g. for rc[None]
        if not key:
            return {**rcParams, **_rcParams_sp}
        # Allow for special time-saving modes where we *ignore rcParams*
        # or even *ignore _rcParams_sp*.
        mode = self._getitem_mode
        if mode==0:
            kws = (self._rcCache, _rcParams_sp, rcParams)
        elif mode==1:
            kws = (self._rcCache, _rcParams_sp)
        elif mode==2:
            kws = (self._rcCache,)
        else:
            raise ValueError(f'Invalid _getitem_mode {mode}.')
        # If it is available, return the values corresponding to names in
        # user dictionary; e.g. {'color':'axes.facecolor'} becomes {'color':'w'}
        # NOTE: Got weird bugs here. Dunno why. Use self.fill method instead.
        if key in _rc_categories:
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
        for kw in (*kws[:1], self._rcGlobal, *kws[1:]):
            try:
                return kw[key]
            except KeyError:
                continue
        # If we were in one of the exlusive modes, return None
        if mode==0:
            raise ValueError(f'Invalid prop name "{key}".')
        else:
            return None

    def __setitem__(self, key, value):
        # First the special cycler
        # NOTE: No matter the 'setitem mode' this will always set the axes
        # prop_cycle rc settings
        if key=='cycle':
            self._set_cycler(value)
            self._rcCache['cycle'] = value
        # Apply global settings
        elif key in _rcGlobal:
            if value=='default':
                value = _rcGlobal[key]
            rc, rc_sp = self._get_globals(key, value)
            self._rcCache.update(rc)
            self._rcCache.update(rc_sp)
            self._rcCache[key] = value # also update cached global property itself
            self._rcGlobal[key] = value
            rcParams.update(rc)
            _rcParams_sp.update(rc_sp)
        # Directly modify single parameter
        # NOTE: If 'setitem mode' is 0, this means user has directly set
        # something (we are not in a with..as context in format()), so we
        # want to directly modify rcParams.
        elif key in _rc_names:
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
        # Nice string representation
        length = 1 + max(len(key) for key in self._rcGlobal.keys())
        string = 'Global settings\n---------------\n' + '\n'.join(f'{key}: {" "*(length-len(key))}{value}'
                                    for key,value in self._rcGlobal.items())
        string += '\n\nAll settings\n------------\n' + ', '.join(_rc_names)
        return string

    def __repr__(self):
        # Simple one, good for auto docs
        length = 1 + max(len(key) for key in self._rcGlobal.keys())
        string = ', '.join(f'{key}: {" "*(length-len(key))}{value}'
                                    for key,value in self._rcGlobal.items())
        return string

    def _set_cycler(self, value):
        # Set the color cycler.
        # NOTE: Generally if user uses 'C0', et cetera, assume they want to
        # refer to the *default* cycler colors; so first reset
        if isinstance(value, str) or isinstance(value, Number):
            value = value,
        colors = colortools.colors('colorblind')
        rcParams['axes.prop_cycle'] = cycler.cycler('color', colors)
        colors = colortools.colors(*value)
        rcParams['axes.prop_cycle'] = cycler.cycler('color', colors)
        figs = list(map(plt.figure, plt.get_fignums()))
        for fig in figs:
            for ax in fig.axes:
                ax.set_prop_cycle(cycler.cycler('color', colors))

    def _get_globals(self, key=None, value=None):
        # Apply all properties in some group.
        kw = {}
        kw_sp = {}
        if key is not None and value is not None:
            items = [(key,value)]
        else:
            items = self._rcGlobal.items()
        for key,value in items:
            # Tick length/major-minor tick length ratio
            if key in ('ticklen','tickratio'):
                if key=='tickratio':
                    ticklen = self._rcGlobal['ticklen']
                    ratio = value
                else:
                    ticklen = value
                    ratio = self._rcGlobal['tickratio']
                kw['xtick.minor.size'] = ticklen*ratio
                kw['ytick.minor.size'] = ticklen*ratio
            # Spine width/major-minor tick width ratio
            if key in ('linewidth','minorwidth'):
                if key=='linewidth':
                    tickwidth = value
                    ratio = self._rcGlobal['minorwidth']
                else:
                    tickwidth = self._rcGlobal['linewidth']
                    ratio = value
                kw['xtick.minor.width'] = tickwidth*ratio
                kw['ytick.minor.width'] = tickwidth*ratio
                # kw_sp['gridminor.linewidth'] = tickwidth*ratio # special
            # Grid line
            if key in ('gridwidth', 'gridratio'):
                if key=='gridwidth':
                    gridwidth = value
                    ratio = self._rcGlobal['gridratio']
                else:
                    gridwidth = self._rcGlobal['gridwidth']
                    ratio = value
                kw_sp['gridminor.linewidth'] = gridwidth*ratio
            # Now update linked settings
            for name in _rcGlobal_children.get(key, []):
                if name in _rcParams_sp:
                    kw_sp[name] = value
                else:
                    kw[name] = value
        return kw, kw_sp

    def _context(self, *args, mode=0, **kwargs):
        """
        Temporarily modify the rc settings. Do this by simply
        saving the old cache, allowing modification of a new cache,
        then restoring the old one.

        This has three modes:

            0. `__getitem__` searches everything, the default.
            1. `__getitem__` ignores `_rcParams` (assumption is these
                have already been set). Used during axes `__init__`
                calls to `_rcupdate`.
            2. `__getitem__` ignores `_rcParams` and `_rcParams_sp`; only
                reads from cache, i.e. settings that user has manually changed.
                Used during `format()` calls to `_rcupdate`.

        Notes
        -----
        This is kept private, because it's only meant to be used within
        the '~proplot.axes.BaseAxes.format' automatically! Instead of the
        user having to use the `with x as y` construct, they should just
        pass an `rc_kw` dict or extra `kwargs` to `~proplot.axes.BaseAxes.format`.
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

    def fill(self, props):
        """
        Fill a dictionary containing rc property string names with the
        corresponding property. This is mostly used internally by
        `~proplot.axes.BaseAxes` and its subclasses.

        Parameters
        ----------
        props : dict-like
            Dictionary whose values are names of rc properties. The values
            are replaced with the corresponding property only if
            `~rc_configurator.__getitem__` does not return ``None``. Otherwise,
            that key, value pair is omitted from the output dictionary.

        Returns
        -------
        dict
            Dictionary filled with rc properties.

        `~rc_configurator.fill` is used to build dictionaries for updating
        `~matplotlib.artist.Artist` instances. Of course, the artist property
        won't need updating unless an rc setting has changed since it was
        drawn (i.e. when the figure and axes were created).

        With this in mind, `~rc_configurator.__getitem__` returns ``None``
        when a setting has not been changed (changed settings are cached
        in a separate dictionary). This prevents hundreds of 1000-item
        dictionary lookups; without caching, runtime increases by 1s even
        for relatively simple plots.
        """
        props_out = {}
        for key,value in props.items():
            value = self[value]
            if value is not None:
                props_out[key] = value
        return props_out

    def update(self, *args, **kwargs):
        """
        Update global settings. Also can be done with `self.name = value`
        or `self['name'] = value`.

        Parameters
        ----------
        category : str, optional
            Category of rc settings to update, passed as the first positional arg.
            If used, all keys are prefixed with the string ``category + '.'``.
        kw : dict-like, optional
            Dictionary of new rc settings, passed as the last positional arg.
        kwargs
            Dictionary of new rc settings, passed as keyword args.
        """
        if len(args)==0:
            args = [{}]
        kw = args[-1]
        kw.update(kwargs)
        if len(args)==2:
            prefix = args[0] + '.'
        elif len(args)==1:
            prefix = ''
        else:
            raise ValueError('Accepts 1-2 positional arguments only. Use rc.update(kw) to update a bunch of names, or rc.update(category, kw) to update subcategories belonging to single category e.g. axes. All kwargs will be added to the dict.')
        for key,value in kw.items():
            self[prefix + key] = value

    def reset(self):
        """
        Restore settings to the default.
        """
        return self.__init__()

# Instantiate object
rc = rc_configurator()
"""
Instance of `rc_configurator`. Use this to change rc settings.

See the `~proplot.rcmod` root page for details.
"""

