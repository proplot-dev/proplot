#!/usr/bin/env python3
"""
This module manges global "rc" settings with a bunch of new convenience
features, and defines some brand new settings. There are now three different
types of settings:

* **Builtin** :ref:`rcParams` settings. These have the format
  ``category.subcategory``, so approaches 1 and 2 are invalid (see below).
* **Custom** :ref:`rcParams_new` settings. These also have the format
  ``category.subcategory``.
* **Global** :ref:`rcGlobals` settings. These are simple, short names
  used to change **multiple** "builtin" and "custom" settings at once,
  or as shorthands for settings with longer names.

Your one-stop-shop for changing settings is the `rc` object. To change
the setting named ``name`` to ``value``, use any of the following 4
approaches:

1. ``rc.name = value``
2. ``rc.update(name=value)``
3. ``rc['name'] = value``
4. ``rc.update({'name':value})``

The "builtin", "custom", and "global" settings are described in detail
below.

rcParams
--------

These are the builtin matplotlib settings. See `this page
<https://matplotlib.org/users/customizing.html>`_ for more info.
The ``rcParams`` categories are as follows:

* Axes: ``axes``.
* Text: ``font``, ``text``, ``mathtext``.
* Plot elements: ``lines``, ``patch``, ``hatch``, ``legend``, ``contour``, ``image``,
  ``boxplot``, ``errorbar``, ``hist``, ``scatter``, ``animation``.
* Axis elements: ``date``, ``xtick``, ``ytick``, ``grid``.
* Printing and saving: ``path``, ``figure``, ``savefig``, ``ps``, ``tk``, ``pdf``, ``svg``.
* Other: ``keymap``, ``examples``, ``debug``.

rcParams_new
------------

My brand new settings, meant to configure special ProPlot featues. The
`rcParams_new` categories are as follows:

* Subplots: ``gridspec``.
* Map settings: ``map``
* Background hatching: ``axes``, ``map``.
* New labels: ``abc``, ``rowlabel``, ``collabel``
* Gridlines: ``gridminor``, ``lonlatlines``

A miscellaneous setting is the boolean ``axes.formatter.zerotrim``; use this
to trim trailing zeros on tick labels. Default is ``True``.

The ``gridspec`` subcategories (all values are in inches):

===================  ==================================================================
Key                  Description
===================  ==================================================================
``gridspec.title``   Vertical space for titles.
``gridspec.legend``  Width of "legend" panels.
``gridspec.cbar``    Width of "colorbar" panels.
``gridspec.ylab``    Horizontal space between subplots alotted for *y*-labels.
``gridspec.xlab``    Vertical space between subplots alotted for *x*-labels.
``gridspec.nolab``   Space between subplots alotted for tick marks.
``gridspec.inner``   Totally empty space between subplots.
===================  ==================================================================

The ``map`` subcategory :

=================  ==================================================================
Key                Description
=================  ==================================================================
``map.reso``       Resolution of geographic features, one of ``'lo'``, ``'med'``, or ``'hi'``
``map.facecolor``  Background color for the map projection.
``map.linewidth``  Line width of map boundary.
``map.edgecolor``  Edge color of map boundary.
=================  ==================================================================

Background hatching with the ``map`` and ``axes`` subcategories :

==================  ==================================================================
Key                 Description
==================  ==================================================================
``xxx.facehatch``   Background hatching string pattern, if not ``None`` [1]_.
``xxx.hatchcolor``  Color of background hatching pattern. Default is same as spines.
==================  ==================================================================

The ``abc``, ``rowlabel``, and ``collabel`` subcategories:

================  ==================================================================
Key               Description
================  ==================================================================
``xxx.fontsize``  The font size.
``xxx.weight``    The font weight.
``xxx.color``     The font color.
================  ==================================================================

The ``lonlatlines`` and ``gridminor`` subcategories:

=================  ==================================================================
Key                Description
=================  ==================================================================
``xxx.linewidth``  The line width.
``xxx.linestyle``  The line style.
``xxx.alpha``      The line transparency.
``xxx.color``      The line color.
=================  ==================================================================

rcGlobals
---------

These settings are used to change :ref:`rcParams` and :ref:`rcParams_new` settings
in bulk, or as shorthands for common settings with longer names.

==================  ==============================================================================================================
Key                 Description
==================  ==============================================================================================================
``cycle``           The default color cycle name, used e.g. for lines.
``cmap``            The default colormap.
``lut``             The number of colors to put in the colormap lookup table.
``color``           The color of axis spines, tick marks, tick labels, and labels.
``xcolor``          As with ``'color'``, but specific to *x*-axes.
``ycolor``          As with ``'color'``, but specific to *y*-axes.
``margin``          The margin of space around subplot `~matplotlib.artist.Artist` instances, if ``xlim`` and ``ylim`` are unset.
``xmargin``         As with ``'margin'``, but specific to the *x* direction.
``ymargin``         As with ``'margin'``, but specific to the *y* direction.
``facecolor``       The axes background color.
``facehatch``       The background hatching pattern [1]_. Useful for highlighting "NaN" data in a ``contourf`` plot.
``small``           Font size for legend text, tick labels, axis labels, and text generated with `~proplot.axes.BaseAxes.text`.
``large``           Font size for titles, "super" titles, and a-b-c subplot labels.
``fontname``        Name of font used for all text in the figure [2]_.
``linewidth``       Thickness of axes spines and major tick lines.
``gridwidth``       Thickness of gridlines.
``gridratio``       Ratio of minor to major gridline thickness.
``gridalpha``       Transparency of major and minor gridlines.
``gridcolor``       Color of major and minor gridlines.
``gridstyle``       Linestyle of major and minor gridlines.
``ticklen``         Length of major ticks.
``tickdir``         Major and minor tick direction; one of ``out``, ``in``, or ``inout``.
``tickratio``       Ratio of minor to major tick line thickness.
``ticklenratio``    Ratio of minor to major tick lengths.
``abcweight``       Font weight for a-b-c labels [3]_.
``titleweight``     Font weight for titles [3]_.
``suptitleweight``  Font weight for "super" titles [3]_.
==================  ==============================================================================================================

.. [1] For example, ``'xxx'`` or ``'..'``. See `this demo <https://matplotlib.org/gallery/shapes_and_collections/hatch_demo.html>`__.
.. [2] This module *changes the default* from DejaVu Sans (or Bitstream Vera) to
       Helvetica Neue or Helvetica. Run `~proplot.fonttools.install_fonts`
       to install them when you download ProPlot for the first time, and
       whenever you update matplotlib.
.. [3] Valid weights are ``'ultralight'``, ``'light'``, ``'normal'``,
       ``'medium'``, ``'demi'``, ``'bold'``, ``'very bold'``, or ``'black'``.
       Note that many fonts only have ``normal`` or ``bold`` available.
       If you request another weight, the “closest” availble weight will
       be selected.

Todo
----
Disable the ``gridspec`` stuff by automating the inter-subplot spacing?

"""
# First import stuff
# try:
#     from icecream import ic
# except ImportError:  # graceful fallback if IceCream isn't installed.
#     ic = lambda *a: None if not a else (a[0] if len(a) == 1 else a) # noqa
import re
import os
import yaml
import cycler
import matplotlib.pyplot as plt
from . import colortools
import numpy as np
import matplotlib as mpl
_rcParams = mpl.rcParams
_rcGlobals = {}
_rcParams_new = {}
_rcGlobals_default = {} # container for ***default*** settings

# Timer
# Cannot import utils because it imports this module
# import time
# import functools
# def _timer(func):
#     """
#     A decorator that prints the time a function takes to execute.
#     See: https://stackoverflow.com/a/1594484/4970632
#     """
#     @functools.wraps(func)
#     def decorator(*args, **kwargs):
#         t = time.clock()
#         print(f'{func.__name__}()')
#         res = func(*args, **kwargs)
#         print(f'{func.__name__}() time: {time.clock()-t}s')
#         return res
#     return decorator

# Get default font
# WARNING: Had issues with Helvetica Neue on Linux, weirdly some characters
# failed to render/printed nonsense, but Helvetica fine
import sys
_user_rc = os.path.join(os.path.expanduser("~"), '.proplotrc')
_default_rc = os.path.join(os.path.dirname(__file__), '.proplotrc') # or parent, but that makes pip install distribution hard
_default_font = 'Helvetica' if sys.platform=='linux' else 'Helvetica Neue' # says 'darwin' on mac
if not os.path.exists(_default_rc):
    raise ValueError('Default configuration file does not exist.')

# "Global" settings and the lower-level settings they change
_rcGlobals_children = {
    # Most important ones, expect these to be used a lot
    # For xcolor/ycolor we just manually use the
    # global property in the format script.
    'cycle':          [], # special handling, passed through Cycle
    'cmap':           [], # special handling, passed through Colormap
    'lut':            ['image.lut'],
    'color':          ['axes.labelcolor', 'axes.edgecolor', 'axes.hatchcolor', 'map.edgecolor', 'map.hatchcolor',
                       'xtick.color', 'ytick.color'], # change the 'color' of an axes
    'hatchlw':        ['hatch.linewidth'],
    'hatchalpha':     [],
    'xcolor':         [], # specially used in the `~matplotlib.axes.XYAxes._rcupdate` function
    'ycolor':         [],
    'landcolor':      [],
    'oceancolor':     [],
    'coastlinewidth': [],
    'coastcolor':     [],
    'facecolor':      ['axes.facecolor', 'map.facecolor'], # simple alias
    'facehatch':      ['axes.facehatch', 'map.facehatch'], # optionally apply background hatching
    'hatchcolor':     ['axes.hatchcolor', 'map.hatchcolor'],
    'gridalpha':      ['grid.alpha', 'gridminor.alpha'],
    'small':          ['font.size', 'xtick.labelsize', 'ytick.labelsize', 'axes.labelsize', 'legend.fontsize'], # the 'small' fonts
    'large':          ['abc.fontsize', 'figure.titlesize', 'axes.titlesize', 'rowlabel.fontsize', 'collabel.fontsize'], # the 'large' fonts
    'linewidth':      ['axes.linewidth', 'map.linewidth', 'hatch.linewidth', 'xtick.major.width', 'ytick.major.width'],
    'gridwidth':      ['grid.linewidth'],
    'gridcolor':      ['grid.color',     'gridminor.color'],
    'gridstyle':      ['grid.linestyle', 'gridminor.linestyle'],
    'ticklen' :       ['xtick.major.size',    'ytick.major.size'],
    'tickpad':        ['xtick.major.pad', 'xtick.minor.pad', 'ytick.major.pad', 'ytick.minor.pad'],
    'tickdir':        ['xtick.direction',     'ytick.direction'],
    'bottom':         ['xtick.major.bottom',  'xtick.minor.bottom'], # major and minor ticks should always be in the same place
    'top':            ['xtick.major.top',     'xtick.minor.top'],
    'left':           ['ytick.major.left',    'ytick.minor.left'],
    'right':          ['ytick.major.right',   'ytick.minor.right'],

    # Simple aliases
    'fontname':       ['font.family'], # specify family directly, so we can easily switch between serif/sans-serif; requires text.usetex = False; see below
    'margin':         ['axes.xmargin', 'axes.ymargin'],
    'xmargin':        ['axes.xmargin'], # found I wanted to change these *a lot*
    'ymargin':        ['axes.ymargin'], # found I wanted to change these *a lot*
    'abcweight':      ['abc.weight'],
    'titleweight':    ['axes.titleweight'],
    'suptitleweight': ['figure.titleweight'],

    # Special ones
    'ticklenratio': [],
    'tickratio': [],
    'gridratio': [],
    }

# Names of the new settings
_rcGlobals_keys = {*_rcGlobals_children.keys()}
_rcParams_new_keys = {
    'axes.formatter.zerotrim',
    'axes.facehatch',
    'axes.hatchcolor',
    'map.facecolor',
    'map.facehatch',
    'map.hatchcolor',
    'map.edgecolor',
    'map.linewidth',
    'map.reso',
    'abc.fontsize',
    'abc.weight',
    'abc.color',
    'rowlabel.fontsize',
    'rowlabel.weight',
    'rowlabel.color',
    'collabel.fontsize',
    'collabel.weight',
    'collabel.color',
    'gridminor.alpha',
    'gridminor.color',
    'gridminor.linestyle',
    'gridminor.linewidth',
    'lonlatlines.alpha',
    'lonlatlines.color',
    'lonlatlines.linewidth',
    'lonlatlines.linestyle',
    'gridspec.title',
    'gridspec.inner',
    'gridspec.legend',
    'gridspec.cbar',
    'gridspec.ylab',
    'gridspec.xlab',
    'gridspec.nolab',
    }

# Generate list of valid names, and names with subcategories
# TODO: Display this somewhere? Maybe in repr of rc?
_rc_names = {
    *_rcParams.keys(),
    *_rcParams_new_keys,
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
        Magical abstract class for managing builtin :ref:`rcParams` settings, 
        our artificial :ref:`rcParams_new` settings, and new "global" settings
        that keep certain groups of settings synced.

        See the module documentation for details.
        """
        # First initialize matplotlib
        # Note rcdefaults() changes the backend! Inline plotting will fail for
        # rest of notebook session if you call rcdefaults before drawing a figure!
        # After first figure made, backend property is 'sticky', never changes!
        # See: https://stackoverflow.com/a/48322150/4970632
        mpl.style.use('default') # mpl.style function does not change the backend

        # Load the defaults from file
        for i,file in enumerate((_default_rc, _user_rc)):
            # Load
            if not os.path.exists(file):
                continue
            with open(file) as f:
                try:
                    data = yaml.safe_load(f)
                except yaml.YAMLError as err:
                    print('Error: Invalid .proplotrc file.')
                    raise err
            # Test
            keys = {*data.keys()}
            if i==0:
                # Check file
                if keys != {'rcGlobals', 'rcParams', 'rcParams_new'}:
                    raise RuntimeError(f'Default .proplotrc file has unexpected sections.')
                # Check contents of each sub dictionary
                _dict = data['rcGlobals']
                if {*_dict.keys()} != _rcGlobals_keys:
                    raise RuntimeError(f'Default .proplotrc file has incomplete or invalid rcGlobals keys.')
                _rcGlobals.update(_dict)
                _rcGlobals_default.update(_dict)
                _dict = data['rcParams_new']
                if {*_dict} != _rcParams_new_keys:
                    raise RuntimeError(f'Default .proplotrc file has incomplete or invalid rcParams_new keys.')
                _rcParams_new.update(_dict)
            else:
                # Check file
                if keys > {'rcGlobals', 'rcParams', 'rcParams_new'}:
                    raise RuntimeError(f'User .proplotrc file has unexpected sections.')
                # Check contents of each sub dictionary
                _dict = data.get('rcGlobals', {})
                if {*_dict.keys()} > _rcGlobals_keys:
                    raise RuntimeError(f'User .proplotrc file has invalid rcGlobals keys.')
                _rcGlobals.update(_dict)
                _dict = data.get('rcParams_new', {})
                if {*_dict.keys()} > _rcGlobals_keys:
                    raise RuntimeError(f'User .proplotrc file has invalid rcParams_new keys.')
                _rcParams_new.update(_dict)

            # Update (this one already checks against invalid keys)
            for key,value in data.get('rcParams', {}).items():
                _rcParams[key] = value

        # Set default fontname and cycler
        # These ones are special. I had issues with Helvetica Neue
        # on UNIX, and Helvetica looked better; but not on macOS.
        if _rcGlobals.get('fontname', None) is None:
            _rcGlobals['fontname'] = _default_font

        # Set the default cycler and colormap (they must
        # be passed through constructor)
        self._set_cmap(_rcGlobals['cmap'])
        self._set_cycler(_rcGlobals['cycle'])

        # Apply *global settings* to children settings
        rc, rc_new = self._get_globals()
        for key,value in rc.items():
            _rcParams[key] = value
        for key,value in rc_new.items():
            _rcParams_new[key] = value

        # Caching stuff
        # TODO: Looks like _getitem_mode can get stuck on a higher, more
        # restrictive value (e.g. 1 or 2) when cell fails to execute. Should
        # consider improving this.
        self._init = True
        self._getitem_mode = 0
        self._setitem_cache = False
        self._cache = {}
        self._cache_restore = {}
        self._context_kwargs = {}
        self._context_cache_backup = {}

    def __enter__(self):
        """
        Apply temporary user global settings.
        """
        self._setitem_cache = True # cache the originals when they are changed?
        self._context_cache_backup = rc._cache.copy()
        for key,value in self._context_kwargs.items():
            self[key] = value # applies globally linked and individual settings

    def __exit__(self, _type, _value, _traceback):
        """
        Restore configurator cache to initial state.
        """
        self._getitem_mode = 0
        self._setitem_cache = False
        for key,value in self._cache_restore.items():
            self[key] = value
        self._cache = self._context_cache_backup
        self._cache_restore = {}
        self._context_kwargs = {}
        self._context_cache_backup = {}

    def __getitem__(self, key):
        """
        Retrieve property. If we are in a `_context` block, will only
        return cached properties (i.e. properties that user wants to
        temporarily change). If not cached, returns None.
        """
        # Can get a whole bunch of different things
        # Get full dictionary e.g. for rc[None]
        if not key:
            return {**_rcParams, **_rcParams_new}
        # Allow for special time-saving modes where we *ignore _rcParams*
        # or even *ignore _rcParams_new*.
        mode = self._getitem_mode
        if mode==0:
            kws = (self._cache, _rcParams_new, _rcParams)
        elif mode==1:
            kws = (self._cache, _rcParams_new)
        elif mode==2:
            kws = (self._cache,)
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
        for kw in (*kws[:1], _rcGlobals, *kws[1:]):
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
        """
        Set :ref:`rcGlobals`, :ref:`rcParams`, or :ref:`rcParams_new` settings.
        """
        # Save changed properties?
        cache = self._cache
        cache[key] = value
        restore = self._setitem_cache
        if restore:
            cache_restore = self._cache_restore
        # First the special cycler and colormaps
        # NOTE: No matter the 'setitem mode' this will always set the axes
        # prop_cycle rc settings
        if key=='cycle':
            if restore:
                for ikey in ('cycle', 'axes.prop_cycle'):
                    cache_restore[ikey] = _rcGlobals[ikey]
            self._set_cycler(value)
        elif key=='cmap':
            if restore:
                for ikey in ('cycle', 'image.cmap'):
                    cache_restore[ikey] = _rcGlobals[ikey]
                cache_restore[key] = _rcGlobals[key]
            self._set_cmap(value)
        # Global settings
        elif key in _rcGlobals:
            # Update globals
            if restore:
                cache_restore[key] = _rcGlobals[key]
            if value=='default':
                value = _rcGlobals_default[key]
            if key=='color': # recursively update global children of global settings
                self['xcolor'] = value
                self['ycolor'] = value
            _rcGlobals[key] = value
            # Update children
            rc, rc_new = self._get_globals(key, value)
            cache.update(rc)
            cache.update(rc_new)
            if restore:
                cache_restore.update({key:_rcParams[key] for key in rc})
                cache_restore.update({key:_rcParams_new[key] for key in rc_new})
            _rcParams.update(rc)
            _rcParams_new.update(rc_new)
        # Directly modify single parameter
        # NOTE: If 'setitem mode' is 0, this means user has directly set
        # something (we are not in a with..as context in format()), so we
        # want to directly modify _rcParams.
        elif key in _rc_names:
            try:
                if restore:
                    cache_restore[key] = _rcParams[key]
                _rcParams[key] = value # rcParams dict has key validation
            except KeyError:
                if restore:
                    cache_restore[key] = _rcParams_new[key]
                _rcParams_new[key] = value
        else:
            raise ValueError(f'Invalid key "{key}".')
        self._init = False # no longer in initial state

    def __getattribute__(self, attr):
        """
        Alias to getitem.
        """
        if attr[:1]=='_' or attr in self._public_api: # no recursion since second comparison won't be evaluated if first comparison evaluates True
            return super().__getattribute__(attr)
        else:
            return self.__getitem__(attr)

    def __setattr__(self, attr, value):
        """
        Alias to setitem.
        """
        if attr[:1]=='_' or attr in self._public_api:
            super().__setattr__(attr, value)
        else:
            self.__setitem__(attr, value)

    def __repr__(self):
        """
        Nice string representation.
        """
        length = 1 + max(len(key) for key in _rcGlobals.keys())
        string = 'rc = {\n' + '\n'.join(f'  {key}: {" "*(length-len(key))}{value}'
                                    for key,value in _rcGlobals.items()) + '\n}'
        return string

    def __str__(self):
        """
        Short string representation.
        """
        length = 1 + max(len(key) for key in _rcGlobals.keys())
        string = ', '.join(f'{key}: {value}' for key,value in _rcGlobals.items())
        return string

    def _set_cmap(self, value):
        """
        Set the default colormap. Value is passed through
        `~proplot.colortools.Colormap`.
        """
        kw = {}
        if np.iterable(value) and len(value)==2 and isinstance(value[-1], dict):
            value, kw = value[0], value[-1]
        if isinstance(value, str) or not np.iterable(value):
            value = value,
        cmap = colortools.Colormap(*value, **kw)
        _rcParams['image.cmap'] = cmap.name

    def _set_cycler(self, value):
        """
        Set the default color cycler. Value is passed through
        `~proplot.colortools.Cycle`.
        """
        # Generally if user uses 'C0', et cetera, assume they want to
        # refer to the *default* cycler colors; so we reset that first
        current = _rcGlobals['cycle']
        _rcParams['axes.prop_cycle'] = cycler.cycler('color', colortools.Cycle(current))
        # Set cycler
        kw = {}
        if np.iterable(value) and len(value)==2 and isinstance(value[-1], dict):
            value, kw = value[0], value[-1]
        if isinstance(value, str) or not np.iterable(value):
            value = value,
        colors = colortools.Cycle(*value, **kw)
        _rcParams['axes.prop_cycle'] = cycler.cycler('color', colors)
        figs = list(map(plt.figure, plt.get_fignums()))
        for fig in figs:
            for ax in fig.axes:
                ax.set_prop_cycle(cycler.cycler('color', colors))

    def _get_globals(self, key=None, value=None):
        """
        Return dictionaries for updating "child" properties in
        `rcParams` and `rcParams_new` with global property.
        """
        kw = {}
        kw_new = {}
        if key is not None and value is not None:
            items = [(key,value)]
        else:
            items = _rcGlobals.items()
        for key,value in items:
            # Tick length/major-minor tick length ratio
            if key in ('ticklen', 'ticklenratio'):
                if key=='ticklenratio':
                    ticklen = _rcGlobals['ticklen']
                    ratio = value
                else:
                    ticklen = value
                    ratio = _rcGlobals['ticklenratio']
                kw['xtick.minor.size'] = ticklen*ratio
                kw['ytick.minor.size'] = ticklen*ratio
            # Spine width/major-minor tick width ratio
            if key in ('linewidth', 'tickratio'):
                if key=='linewidth':
                    tickwidth = value
                    ratio = _rcGlobals['tickratio']
                else:
                    tickwidth = _rcGlobals['linewidth']
                    ratio = value
                kw['xtick.minor.width'] = tickwidth*ratio
                kw['ytick.minor.width'] = tickwidth*ratio
                # kw_new['gridminor.linewidth'] = tickwidth*ratio # special
            # Grid line
            if key in ('gridwidth', 'gridratio'):
                if key=='gridwidth':
                    gridwidth = value
                    ratio = _rcGlobals['gridratio']
                else:
                    gridwidth = _rcGlobals['gridwidth']
                    ratio = value
                kw_new['gridminor.linewidth'] = gridwidth*ratio
            # Now update linked settings
            for name in _rcGlobals_children[key]:
                if name in _rcParams_new:
                    kw_new[name] = value
                else:
                    kw[name] = value
        return kw, kw_new

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
        2. `__getitem__` ignores `_rcParams` and `_rcParams_new`; only
            reads from cache, i.e. settings that user has manually changed.
            Used during `~proplot.BaseAxes.format` calls to
            `~proplot.BaseAxes._rcupdate`.

        Notes
        -----
        This is kept private, because it's only meant to be used within
        the '~proplot.axes.BaseAxes.format' automatically! Instead of the
        user having to use the ``with x as y`` construct, they should just
        pass an `rc_kw` dict or extra `kwargs` to `~proplot.axes.BaseAxes.format`.
        """
        # Apply mode
        if mode not in range(3):
            raise ValueError(f'Invalid _getitem_mode {mode}.')
        for arg in args:
            if not isinstance(arg, dict):
                raise ValueError('rc_context() only accepts dictionary args and kwarg pairs.')
            kwargs.update(arg)
        self._context_kwargs = kwargs # could be empty
        self._getitem_mode = mode
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

        Notes
        -----
        `~rc_configurator.fill` is used to build dictionaries for updating
        `~matplotlib.artist.Artist` instances. Of course, the artist property
        won't need updating unless an rc setting has changed since it was
        instantiated.

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
        Update global settings. Also can be done with ``self.name = value``
        or ``self['name'] = value``.

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

