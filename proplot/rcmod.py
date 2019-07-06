#!/usr/bin/env python3
"""
A special object named `~proplot.rcmod.rc`, belonging to the
`~proplot.rcmod.rc_configurator` class, is created on import.
This is your **one-stop shop for changing global settings** belonging to any of the
following three categories.

1. Builtin matplotlib `rcParams <https://matplotlib.org/users/customizing.html>`__
   settings. These have the format ``x.y`` or ``x.y.z``.
2. ProPlot :ref:`rcCustom` settings. These also have the format ``x.y`` (see below).
3. ProPlot :ref:`rcGlobals` settings. These have no dots (see below). They are **simple,
   short** names used to change multiple matplotlib and ProPlot settings at once,
   as shorthands for settings with longer names, or for special options. For example,
   ``ticklen`` changes the tick length for the *x* and *y* axes in one go.

You can change settings with the `~proplot.rcmod.rc` object as follows.

* ``plot.rc.name = value``
* ``plot.rc['name'] = value``
* ``plot.rc.update(name1=value1, name2=value2)``
* ``plot.rc.update({'name1':value1, 'name2':value2})``

To temporarily change settings on a particular axes, use either of the following.

* ``ax.format(name=value)``
* ``ax.format(rc_kw={'name':value})``

In all of these examples, if the setting name ``name`` contains
any dots, you can simply **omit the dots**. For example, to change the
:ref:`rcCustom` property ``title.loc``, use ``plot.rc.titleloc = value``,
``plot.rc.update(titleloc=value)``, or ``ax.format(titleloc=value)``.

..
    The ``rcParams`` categories can be grouped as follows:

    * Axes: ``axes``.
    * Text: ``font``, ``text``, ``mathtext``.
    * Plot elements: ``lines``, ``patch``, ``hatch``, ``legend``, ``contour``, ``image``,
    ``boxplot``, ``errorbar``, ``hist``, ``scatter``, ``animation``.
    * Axis elements: ``date``, ``xtick``, ``ytick``, ``grid``.
    * Printing and saving: ``path``, ``figure``, ``savefig``, ``ps``, ``tk``, ``pdf``, ``svg``.
    * Other: ``keymap``, ``examples``, ``debug``.

#########
rcGlobals
#########

==================  ==================================================================================================================================================================
Key                 Description
==================  ==================================================================================================================================================================
``nbsetup``         Whether to run `nb_setup` on import.
``autosave``        If non-empty and ``nbsetup`` is ``True``, passed to `%autosave <https://www.webucator.com/blog/2016/03/change-default-autosave-interval-in-ipython-notebook/>`__.
``autoreload``      If non-empty and ``nbsetup`` is ``True``, passed to `%autoreload <https://ipython.readthedocs.io/en/stable/config/extensions/autoreload.html#magic-autoreload>`__.
``tight``           Whether to auto-adjust figure bounds and subplot spacings.
``fontname``        Name of font used for all text in the figure. The default is ``Helvetica`` for Linux and ``Helvetica Neue`` for Windows/Mac. See `~proplot.fonttools` for details.
``cmap``            The default colormap.
``lut``             The number of colors to put in the colormap lookup table.
``cycle``           The default color cycle name, used e.g. for lines.
``rgbcycle``        Whether to register cycles names as ``'r'``, ``'b'``, ``'g'``, etc., like in `seaborn <https://seaborn.pydata.org/tutorial/color_palettes.html>`__.
``color``           The color of axis spines, tick marks, tick labels, and labels.
``alpha``           The opacity of the background axes patch.
``facecolor``       The color of the background axes patch.
``small``           Font size for legend text, tick labels, axis labels, and text generated with `~matplotlib.axes.Axes.text`.
``large``           Font size for titles, "super" titles, and a-b-c subplot labels.
``linewidth``       Thickness of axes spines and major tick lines.
``margin``          The margin of space between axes edges and objects plotted inside the axes, if ``xlim`` and ``ylim`` are unset.
``ticklen``         Length of major ticks in points.
``tickdir``         Major and minor tick direction. Must be one of ``out``, ``in``, or ``inout``.
``tickpad``         Padding between ticks and tick labels in points.
``grid``            Boolean, toggles major grid lines on and off.
``gridminor``       Boolean, toggles minor grid lines on and off.
``tickratio``       Ratio of minor tickline width to major tickline width.
``gridratio``       Ratio of minor gridline width to major gridline width.
``ticklenratio``    Ratio of minor tickline length to major tickline length.
``reso``            Resolution of geographic features, one of ``'lo'``, ``'med'``, or ``'hi'``
``geogrid``         Boolean, toggles meridian and parallel gridlines on and off.
``land``            Boolean, toggles land patches on and off.
``ocean``           Boolean, toggles ocean patches on and off.
``lakes``           Boolean, toggles lake patches on and off.
``coast``           Boolean, toggles coastline lines on and off.
``borders``         Boolean, toggles country border lines on and off.
``innerborders``    Boolean, toggles internal border lines on and off, e.g. for states and provinces.
``rivers``          Boolean, toggles river lines on and off.
==================  ==================================================================================================================================================================

#############
rcCustom
#############
For a-b-c labels and axes title settings, use the ``abc`` and ``title``
categories. For axis tick label settings, use the ``tick`` category.
For figure title, row label, and column label settings, use the new
``suptitle``, ``rowlabel``, and ``collabel`` categories.

======================================  =================================================================================================
Key                                     Description
======================================  =================================================================================================
``abc.format``                          a-b-c label format (for options, see `~proplot.axes.BaseAxes.format_partial`).
``abc.loc``, ``title.loc``              a-b-c label or title position (for options, see `~proplot.axes.BaseAxes.format_partial`).
``abc.border``, ``title.border``        Boolean, indicates whether to draw labels inside the axes with a white border.
``abc.linewidth``, ``title.linewidth``  Width of the (optional) white border.
``xxxx.color``                          The font color, valid for ``abc``, ``title``, ``rowlabel``, ``collabel``, and ``suptitle``.
``xxxx.fontsize``                       The font size, valid for ``abc``, ``title``, ``rowlabel``, ``collabel``, and ``suptitle``.
``xxxx.weight``                         The font weight [1]_, valid for ``abc``, ``title``, ``rowlabel``, ``collabel``, and ``suptitle``.
``tick.labelweight``                    Weight [1]_ of axis tick labels, mirrors ``axes.labelweight``.
``tick.labelcolor``                     Color of axis tick labels, mirrors ``axes.labelcolor``.
``tick.labelsize``                      Font size of axis tick labels, mirrors ``axes.labelsize``.
``axes.formatter.zerotrim``             Boolean, indicates whether trailing decimal zeros are trimmed on tick labels.
``axes.formatter.timerotation``         Float, indicates the default *x* axis tick label rotation for datetime tick labels.
======================================  =================================================================================================

.. [1] Valid font weights are ``'ultralight'``, ``'light'``, ``'normal'``,
       ``'medium'``, ``'demi'``, ``'bold'``, ``'very bold'``, or ``'black'``.
       Many fonts only have ``normal`` or ``bold``. If you request an
       unavailable weight, matplotlib picks the closest availble weight.

For minor gridlines, use the ``gridminor`` category. For meridian and parallel
gridlines on `~proplot.axes.MapAxes`, use the ``geogrid`` category``.
If a property is empty, the corresponding property from the buildin ``grid``
category is used.

==============================================  ==================================================================
Key                                             Description
==============================================  ==================================================================
``gridminor.linewidth``, ``geogrid.linewidth``  The line width.
``gridminor.linestyle``, ``geogrid.linestyle``  The line style.
``gridminor.alpha``, ``geogrid.alpha``          The line transparency.
``gridminor.color``, ``geogrid.color``          The line color.
==============================================  ==================================================================

The below properties are particular to `~proplot.axes.MapAxes`. The properties
for geographic elements like ``land`` are used when the corresponding
:ref:`rcGlobals` geographic feature is toggled on.

==================================================  ==================================================================
Key                                                 Description
==================================================  ==================================================================
``geogrid.labels``                                  Boolean, indicates whether to label the parallels and meridians.
``geogrid.latmax``                                  Meridian gridlines are cut off poleward of this latitude.
``geogrid.lonstep``                                 Interval for meridian gridlines, in degrees.
``geogrid.latstep``                                 Interval for parallel gridlines, in degrees.
``geogrid.labelsize``                               Font size for latitide and longitude labels.
``land.color``                                      Face color for continents and island patches.
``ocean.color``                                     Face color for ocean patches.
``lakes.color``                                     Face color for lake patches.
``rivers.color``, ``rivers.linewidth``              Line color and linewidth for river lines.
``borders.color``, ``borders.linewidth``            Line color and linewidth for country border lines.
``innerborders.color``, ``innerborders.linewidth``  Line color and linewidth for internal border lines.
==================================================  ==================================================================

The ``subplot`` category has settings that control the default figure
layout. As with all sizing arguments, if specified as a number, the units
are inches. If string, the units are interpreted by `~proplot.utils.units`.

=======================  ==================================================================
Key                      Description
=======================  ==================================================================
``subplot.axwidth``      Default width of each axes.
``subplot.panelwidth``   Width of side panels.
``subplot.cbarwidth``    Width of "colorbar" panels.
``subplot.legwidth``     Width of "legend" panels.
``subplot.borderpad``    Padding around edges for tight subplot.
``subplot.subplotpad``   Padding between main subplots.
``subplot.panelpad``     Padding between axes panels and their parents.
``subplot.titlespace``   Vertical space for titles.
``subplot.ylabspace``    Horizontal space between subplots alotted for *y*-labels.
``subplot.xlabspace``    Vertical space between subplots alotted for *x*-labels.
``subplot.nolabspace``   Space between subplots alotted for tick marks.
``subplot.innerspace``   Purely empty space between subplots.
``subplot.panelspace``   Purely empty space between main axes and side panels.
=======================  ==================================================================

The ``colorbar`` category has settings that control the default colorbar
layout. It is analogous to the builtin ``legend`` category, but configures
both *inset* and *panel* colorbars (see the `~proplot.axes.BaseAxes`
`~proplot.axes.BaseAxes.colorbar` and `~proplot.axes.PanelAxes`
`~proplot.axes.PanelAxes.colorbar` methods for details).

========================  =========================================================================================================================
Key                       Description
========================  =========================================================================================================================
``colorbar.frameon``      Boolean, indicates whether to draw a frame behind inset colorbars.
``colorbar.framealpha``   The opacity of colorbar frame.
``colorbar.loc``          Default inset colorbar location, one of "upper right", "upper left", "lower left", or "lower right", or "center" options.
``colorbar.length``       Default length of inset colorbars.
``colorbar.width``        Default width of inset colorbars.
``colorbar.extend``       Length of rectangular or triangular "extensions" for panel colorbars.
``colorbar.extendinset``  Length of rectangular or triangular "extensions" for inset colorbars.
``colorbar.pad``          Default padding between figure edge of rectangular or triangular "extensions" for inset colorbars.
``colorbar.xspace``       Extra space for *x* label of inset colorbars.
========================  =========================================================================================================================
"""
# First import stuff
# WARNING: Must import pyplot here, because otherwise 'style' attribute
# is not added to matplotlib top-level module!
import re
import os
import yaml
import cycler
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
import numpy as np
import matplotlib as mpl
import warnings
try:
    import IPython
    get_ipython = IPython.get_ipython
except ModuleNotFoundError:
    get_ipython = lambda: None
from .utils import ic, units, _timer, _counter
_rcParams = mpl.rcParams
_rcGlobals = {}
_rcCustom = {}

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
# NOTE: This whole section, declaring dictionaries and sets, takes 1ms
_rcGlobals_children = {
    # Notebooks
    'nbsetup':      [],
    'autosave':     [],
    'autoreload':   [],
    # Subplots
    'tight':        [],
    # Style
    'fontname':     ['font.family'], # specify family directly, so we can easily switch between serif/sans-serif; requires text.usetex = False; see below
    'cmap':         [],
    'lut':          ['image.lut'],
    'cycle':        [],
    'rgbcycle':     [],
    'color':        ['axes.labelcolor', 'axes.edgecolor', 'tick.labelcolor', 'hatch.color', 'xtick.color', 'ytick.color'], # change the 'color' of an axes
    'alpha':        ['axes.alpha'], # this is a custom setting
    'facecolor':    ['axes.facecolor'],
    'small':        ['font.size', 'tick.labelsize', 'xtick.labelsize', 'ytick.labelsize', 'axes.labelsize', 'legend.fontsize', 'geogrid.labelsize'], # the 'small' fonts
    'large':        ['abc.fontsize', 'figure.titlesize', 'axes.titlesize', 'suptitle.fontsize', 'title.fontsize', 'rowlabel.fontsize', 'collabel.fontsize'], # the 'large' fonts
    'linewidth':    ['axes.linewidth', 'hatch.linewidth', 'xtick.major.width', 'ytick.major.width'],
    'margin':       ['axes.xmargin', 'axes.ymargin'],
    'grid':         ['axes.grid'],
    'gridminor':    ['axes.gridminor'],
    'ticklen' :     ['xtick.major.size', 'ytick.major.size'],
    'tickdir':      ['xtick.direction',  'ytick.direction'],
    'tickpad':      ['xtick.major.pad', 'xtick.minor.pad', 'ytick.major.pad', 'ytick.minor.pad'],
    'tickratio':    [],
    'ticklenratio': [],
    'gridratio':    [],
    # Geography
    'reso':         [],
    'geogrid':      ['axes.geogrid'],
    'land':         [],
    'ocean':        [],
    'lakes':        [],
    'coast':        [],
    'borders':      [],
    'innerborders': [],
    'rivers':       [],
    }

# Names of the new settings
_rc_names_global = {*_rcGlobals_children.keys()}
_rc_names_old = {*_rcParams.keys()}
_rc_names_new = {
    'axes.formatter.zerotrim', 'axes.formatter.timerotation',
    'axes.gridminor', 'axes.geogrid', 'axes.alpha',
    'land.color', 'ocean.color', 'lakes.color', 'coast.color', 'coast.linewidth',
    'borders.color', 'borders.linewidth', 'innerborders.color', 'innerborders.linewidth', 'rivers.color', 'rivers.linewidth',
    'abc.fontsize', 'abc.weight', 'abc.color', 'abc.loc', 'abc.format', 'abc.border', 'abc.linewidth',
    'title.loc', 'title.color', 'title.border', 'title.linewidth', 'title.weight', 'title.fontsize',
    'suptitle.fontsize', 'suptitle.weight', 'suptitle.color',
    'rowlabel.fontsize', 'rowlabel.weight', 'rowlabel.color',
    'collabel.fontsize', 'collabel.weight', 'collabel.color',
    'gridminor.alpha', 'gridminor.color', 'gridminor.linestyle', 'gridminor.linewidth',
    'geogrid.labels', 'geogrid.alpha', 'geogrid.color', 'geogrid.labelsize', 'geogrid.linewidth', 'geogrid.linestyle', 'geogrid.latmax', 'geogrid.lonstep', 'geogrid.latstep',
    'subplot.subplotpad', 'subplot.panelpad', 'subplot.borderpad', 'subplot.titlespace',
    'subplot.innerspace',
    'tick.labelweight', 'tick.labelcolor', 'tick.labelsize',
    'subplot.legwidth', 'subplot.cbarwidth', 'subplot.ylabspace', 'subplot.xlabspace', 'subplot.nolabspace',
    'subplot.axwidth', 'subplot.panelwidth', 'subplot.panelspace',
    'colorbar.frameon', 'colorbar.framealpha', 'colorbar.length', 'colorbar.width', 'colorbar.loc', 'colorbar.extend', 'colorbar.extendinset', 'colorbar.axespad', 'colorbar.xspace',
    }
# Used by BaseAxes.format, allows user to pass rc settings as keyword args,
# way less verbose. For example, compare landcolor='b' to
# rc_kw={'land.color':'b'}.
# Not used by getitem, because that would be way way too many lookups; is
# only looked up if user *manually* passes something to BaseAxes.format.
_rc_names_nodots = { # useful for passing these as kwargs
    name.replace('.', ''):name for names in
    (_rc_names_new, _rc_names_old, _rc_names_global)
    for name in names
    }
# Categories for returning dict of subcategory properties
_rc_categories = {
    *(re.sub('\.[^.]*$', '', name) for names in (_rc_names_new, _rc_names_old) for name in names),
    *(re.sub('\..*$', '', name) for names in (_rc_names_new, _rc_names_old) for name in names)
    }

#-------------------------------------------------------------------------------
# Contextual settings management
# Adapted from seaborn; see: https://github.com/mwaskom/seaborn/blob/master/seaborn/rcmod.py
#-------------------------------------------------------------------------------
def rc_defaults():
    """Resets all settings to the matplotlib defaults."""
    mpl.style.use('default') # mpl.style function does not change the backend

class _locked(dict):
    """Locked dictionary."""
    __getattr__ = dict.get
    def __repr__(self):
        """Wrap notebook pretty print."""
        string = ''
        prefix = '{'
        if len(self)<=5:
            suffix = ', '
        else:
            suffix = ',\n '
        for key,value in self.items():
            string += (prefix + repr(key) + ': ' + repr(value))
            prefix = suffix
        return '<_locked {' + string + '}>'
    def __setattr__(self, attr, value):
        """Locked."""
        raise NotImplementedError
    def __setitem__(self, attr, value):
        """Locked."""
        raise NotImplementedError
    def __delattr__(self, attr, value):
        """Locked."""
        raise NotImplementedError
    def __delitem__(self, attr, value):
        """Locked."""
        raise NotImplementedError

class rc_configurator(object):
    _public_api = ('reset', 'context', 'get', 'update', 'fill') # getattr and setattr will not look for these items on underlying dictionary
    # @_counter # about 0.05s
    def __init__(self):
        """Magical abstract class for managing matplotlib `rcParams
        <https://matplotlib.org/users/customizing.html>`__ settings, ProPlot
        :ref:`rcCustom` settings, and :ref:`rcGlobals` "global" settings.
        See the `~proplot.rcmod` documentation for details."""
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
                if keys != {'rcGlobals', 'rcParams', 'rcCustom'}:
                    raise RuntimeError(f'Default .proplotrc file has unexpected sections.')
                # Check contents of each sub dictionary
                _dict = data['rcGlobals']
                if {*_dict.keys()} != _rc_names_global:
                    print({*_dict.keys()} - _rc_names_global)
                    print(_rc_names_global - {*_dict.keys()})
                    raise RuntimeError(f'Default .proplotrc file has incomplete or invalid rcGlobals keys.')
                _rcGlobals.update(_dict)
                _dict = data['rcCustom']
                if {*_dict} != _rc_names_new:
                    print({*_dict.keys()} - _rc_names_new)
                    print(_rc_names_new - {*_dict.keys()})
                    raise RuntimeError(f'Default .proplotrc file has incomplete or invalid rcCustom keys.')
                _rcCustom.update(_dict)
            else:
                # Check file
                if keys > {'rcGlobals', 'rcParams', 'rcCustom'}:
                    raise RuntimeError(f'User .proplotrc file has unexpected sections.')
                # Check contents of each sub dictionary
                _dict = data.get('rcGlobals', {})
                if {*_dict.keys()} > _rc_names_global:
                    raise RuntimeError(f'User .proplotrc file has invalid rcGlobals keys.')
                _rcGlobals.update(_dict)
                _dict = data.get('rcCustom', {})
                if {*_dict.keys()} > _rc_names_global:
                    raise RuntimeError(f'User .proplotrc file has invalid rcCustom keys.')
                _rcCustom.update(_dict)

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
            _rcCustom[key] = value

        # Caching stuff
        self._init = True
        self._no_dict = False # whether to prevent returning dictionary of category values
        self._getitem_mode = 0
        self._setitem_cache = False
        self._cache = {}
        self._cache_restore = {}
        self._context_kwargs = {}
        self._context_cache_backup = {}

    def __getitem__(self, key):
        """Retrieves property. If we are in a `~rc_configurator.context`
        block, will return ``None`` if property is not cached (i.e.
        property was changed by user)."""
        # Can get a whole bunch of different things
        # Get full dictionary e.g. for rc[None]
        if not key:
            return {**_rcParams, **_rcCustom}

        # Standardize
        # NOTE: If key is invalid, and caching disabled, will get error
        # down the line. Caching should only be enabled for internal use.
        if '.' not in key and key not in _rcGlobals:
            key = _rc_names_nodots.get(key, key)

        # Allow for special time-saving modes where we *ignore _rcParams*
        # or even *ignore _rcCustom*.
        mode = self._getitem_mode
        if mode==0:
            kws = (self._cache, _rcGlobals, _rcCustom, _rcParams)
        elif mode==1:
            kws = (self._cache, _rcGlobals, _rcCustom)
        elif mode==2:
            kws = (self._cache, _rcGlobals)
        else:
            raise KeyError(f'Invalid _getitem_mode {mode}.')

        # If it is available, return the values corresponding to names in
        # user dictionary; e.g. {'color':'axes.facecolor'} becomes {'color':'w'}
        # Just use this to see what is available
        # NOTE: If have global property name that matches a category name (e.g.
        # 'land' and 'land.color'), use get() to temporarily change mode to
        # prevent looking for categories.
        if not self._no_dict and key in _rc_categories:
            params = {}
            for kw in kws:
                for category,value in kw.items():
                    if re.search(f'^{key}\.', category):
                        subcategory = re.sub(f'^{key}\.', '', category)
                        if subcategory and '.' not in subcategory:
                            params[subcategory] = value
            return _locked(params)

        # Get individual property. Will successively index a few different dicts
        # Try to return the value
        for kw in kws:
            try:
                return kw[key]
            except KeyError:
                continue
        # If we were in one of the exclusive modes, return None
        if mode==0:
            raise KeyError(f'Invalid prop name "{key}".')
        else:
            return None

    def __setitem__(self, key, value):
        """Sets `rcParams <https://matplotlib.org/users/customizing.html>`__,
        :ref:`rcCustom`, and :ref:`rcGlobals` settings."""
        # Standardize name
        # Save changed properties?
        cache = self._cache
        cache[key] = value
        restore = self._setitem_cache
        if restore:
            cache_restore = self._cache_restore

        # Standardize
        # NOTE: If key is invalid, raise error down the line.
        if '.' not in key and key not in _rcGlobals:
            key = _rc_names_nodots.get(key, key)

        # First the special cycler and colormaps
        # NOTE: We don't need to cache these because props aren't looked up
        # by format function; plot methods will read global properties
        if key=='rgbcolors':
            _rcGlobals[key] = value
            key = 'cycle'
            value = _rcGlobals[key]
        if key=='cycle':
            if restore:
                cache_restore[key] = _rcGlobals[key]
                cache_restore['patch.facecolor'] = _rcParams['patch.facecolor']
                cache_restore['axes.prop_cycle'] = _rcParams['axes.prop_cycle']
            self._set_cycler(value)
        elif key=='cmap':
            if restore:
                cache_restore[key] = _rcGlobals[key]
                cache_restore['image.cmap'] = _rcParams['image.cmap']
            self._set_cmap(value)

        # Grid properties, not sure why it has to be this complicated
        elif key in ('grid', 'gridminor'):
            value_orig = _rcParams['axes.grid']
            which_orig = _rcParams['axes.grid.which']
            if restore:
                cache_restore[key] = _rcGlobals[key]
                cache_restore['axes.grid'] = value_orig
                cache_restore['axes.grid.which'] = which_orig
            if not value:
                # Turn off gridlines for just major or minor ticks
                if not value_orig or (which_orig=='major' and key=='grid') or (which_orig=='minor' and key=='gridminor'):
                    which = 'both' # disable both sides
                # Some gridlines were already on! Does this new setting disable
                # everything, or just one of them?
                elif which_orig=='both': # and value_orig is True, as we already tested
                    value = True
                    which = 'major' if key=='gridminor' else 'minor' # if gridminor=False, enable major, and vice versa
                # Do nothing; e.g. if original was which=minor, user input grid=False,
                # nothing should happen.
                else:
                    value = True
                    which = which_orig
            else:
                # Turn on gridlines for major and minor, if was just one before
                if which_orig=='both' or (key=='grid' and which_orig=='minor') or (key=='gridminor' and which_orig=='major'):
                    which = 'both'
                # Do nothing; e.g. if original was which=major, user input was grid=True
                else:
                    value = True
                    which = which_orig
            cache.update({'axes.grid':value, 'axes.grid.which':which})
            _rcParams.update({'axes.grid':value, 'axes.grid.which':which})

        # Global settings
        elif key in _rcGlobals:
            # Update globals
            if restore:
                cache_restore[key] = _rcGlobals[key]
            _rcGlobals[key] = value
            # Update children
            rc, rc_new = self._get_globals(key, value)
            cache.update(rc)
            cache.update(rc_new)
            if restore:
                cache_restore.update({key:_rcParams[key] for key in rc})
                cache_restore.update({key:_rcCustom[key] for key in rc_new})
            _rcParams.update(rc)
            _rcCustom.update(rc_new)

        # Directly modify single parameter
        elif key in _rc_names_new:
            if restore:
                cache_restore[key] = _rcCustom[key]
            _rcCustom[key] = value
        elif key in _rc_names_old:
            if restore:
                cache_restore[key] = _rcParams[key]
            _rcParams[key] = value # rcParams dict has key validation
        else:
            raise KeyError(f'Invalid key "{key}".')
        self._init = False # no longer in initial state

    def __getattribute__(self, attr):
        """Alias to getitem."""
        # No recursion since second comparison won't be evaluated if first comparison evaluates True
        if attr[:1]=='_' or attr in self._public_api:
            return object.__getattribute__(self, attr)
        else:
            return self.__getitem__(attr)

    def __setattr__(self, attr, value):
        """Alias to setitem."""
        if attr[:1]=='_' or attr in self._public_api:
            object.__setattr__(self, attr, value)
        else:
            self.__setitem__(attr, value)

    def __enter__(self):
        """Apply temporary user global settings."""
        self._setitem_cache = True # cache the originals when they are changed?
        self._context_cache_backup = rc._cache.copy()
        for key,value in self._context_kwargs.items():
            self[key] = value # applies globally linked and individual settings

    def __exit__(self, _type, _value, _traceback):
        """Restore configurator cache to initial state."""
        self._getitem_mode = 0
        self._setitem_cache = False
        for key,value in self._cache_restore.items():
            self[key] = value
        self._cache = self._context_cache_backup
        self._cache_restore = {}
        self._context_kwargs = {}
        self._context_cache_backup = {}

    def __delitem__(self, *args):
        """Disable."""
        raise NotImplementedError

    def __delattr__(self, *args):
        """Disable."""
        raise NotImplementedError

    def __repr__(self):
        """Nice string representation."""
        length = 1 + max(len(key) for key in _rcGlobals.keys())
        string = 'rc = {\n' + '\n'.join(f'  {key}: {" "*(length-len(key))}{value}'
                                    for key,value in _rcGlobals.items()) + '\n}'
        return string

    def __str__(self):
        """Short string representation."""
        length = 1 + max(len(key) for key in _rcGlobals.keys())
        string = ', '.join(f'{key}: {value}' for key,value in _rcGlobals.items())
        return string

    def _set_cmap(self, name):
        """Sets the default colormap."""
        # Draw from dictionary
        cmap = mcm.cmap_d[name] # colortools.Colormap(*value)
        _rcParams['image.cmap'] = cmap.name

    def _set_cycler(self, name):
        """Sets the default color cycler."""
        # Draw from dictionary
        colors = mcm.cmap_d[name].colors # colortools.Cycle(*value)
        # Apply new color name definitions
        if _rcGlobals['rgbcycle']:
            if name.lower()=='colorblind': # apply
                regcolors = colors + [(0.1, 0.1, 0.1)]
            else: # reset
                regcolors = [(0.0, 0.0, 1.0), (0.0, .50, 0.0), (1.0, 0.0, 0.0), (.75, .75, 0.0), (.75, .75, 0.0), (0.0, .75, .75), (0.0, 0.0, 0.0)]
            for code,color in zip('brgmyck', regcolors):
                rgb = mcolors.colorConverter.to_rgb(color)
                mcolors.ColorConverter.colors[code] = rgb
                mcolors.ColorConverter.cache[code]  = rgb
        # Pass to cycle constructor
        # WARNING: Used to apply to every axes: list(map(plt.figure, plt.get_fignums()))
        # This was dumb because should not be expected in general that properties
        # apply to figures that already exist, right?
        _rcParams['patch.facecolor'] = colors[0]
        _rcParams['axes.prop_cycle'] = cycler.cycler('color', colors)

    def _get_globals(self, key=None, value=None):
        """Returns dictionaries for updating "child" properties in
        `rcParams` and `rcCustom` with global property."""
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
                kw['xtick.major.size'] = ticklen
                kw['ytick.major.size'] = ticklen
                kw['xtick.minor.size'] = ticklen*ratio
                kw['ytick.minor.size'] = ticklen*ratio
            # Spine width/major-minor tick width ratio
            elif key in ('linewidth', 'tickratio'):
                if key=='linewidth':
                    tickwidth = value
                    ratio = _rcGlobals['tickratio']
                else:
                    tickwidth = _rcGlobals['linewidth']
                    ratio = value
                kw['axes.linewidth'] = tickwidth
                kw['hatch.linewidth'] = tickwidth
                kw['xtick.major.width'] = tickwidth
                kw['ytick.major.width'] = tickwidth
                kw['xtick.minor.width'] = tickwidth*ratio
                kw['ytick.minor.width'] = tickwidth*ratio
            # Grid line
            # NOTE: Changing minor width does not affect major width
            elif key in ('grid.linewidth', 'gridratio'):
                if key=='gridratio':
                    gridwidth = _rcParams['grid.linewidth']
                    ratio = value
                else:
                    gridwidth = value
                    ratio = _rcGlobals['gridratio']
                kw['grid.linewidth'] = gridwidth
                kw_new['gridminor.linewidth'] = gridwidth*ratio
            # Now update linked settings
            else:
                for name in _rcGlobals_children[key]:
                    if name in _rcCustom:
                        kw_new[name] = value
                    else:
                        kw[name] = value
        return kw, kw_new

    # Internally used, but public methods.
    def context(self, *args, mode=0, **kwargs):
        """
        Temporarily modify the rc settings. Do this by simply
        saving the old cache, allowing modification of a new cache,
        then restoring the old one.

        Possible arguments for `mode` are ``0``, ``1``, and ``2``:

        0. `~rc_configurator.__getitem__` searches everything, the default.
        1. `~rc_configurator.__getitem__` ignores `rcParams <https://matplotlib.org/users/customizing.html>`__
           keys. The assumption is these have already been applied. This is used
           during the `~proplot.axes.BaseAxes.__init__` calls to
           `~proplot.axes.BaseAxes.format`.
        2. `~rc_configurator.__getitem__` ignores `rcParams <https://matplotlib.org/users/customizing.html>`__
           and :ref:`rcCustom` keys. It only reads from cache, i.e. settings that
           user has manually changed. This is used during user calls to
           `~proplot.axes.BaseAxes.format`.

        Warning
        -------
        You should not have to use this directly! To temporarily change
        settings, just pass them to `~proplot.axes.BaseAxes.format` (which
        calls this method internally and applies the new properties).
        """
        # Apply mode
        if mode not in range(3):
            raise ValueError(f'Invalid _getitem_mode {mode}.')
        for arg in args:
            if not isinstance(arg, dict):
                raise ValueError('Non-dictionary argument.')
            kwargs.update(arg)
        self._context_kwargs = kwargs # could be empty
        self._getitem_mode = mode
        return self

    def _get(self, key):
        """Powers `get` and `fill` methods."""
        try:
            item = self[key]
        except KeyError as err:
            self._no_dict = False
            self._getitem_mode = orig
            raise err
        return item

    def get(self, key, cache=False, nodict=True):
        """Getter with some special options. This is for internal use by
        `~proplot.axes.BaseAxes` and its subclasses. Defaults to looking
        for **all** properties (i.e. not just the cached ones).

        Parameters
        ----------
        key : str
            The key name.
        cache : bool, optional
            Whether to look for all properties or just cached (i.e.
            recently changed) properties.
        nodict : bool, optional
            Whether to disable the `~rc_configurator.__getitem__` feature
            that returns a **dictionary** of properties belonging to
            some category. This can be helpful when a category and a "global"
            setting have the same name, like "land".
        """
        self._no_dict = nodict
        if not cache:
            orig = self._getitem_mode
            self._getitem_mode = 0
        item = self._get(key)
        self._no_dict = False
        if not cache:
            self._getitem_mode = orig
        return item

    def fill(self, props, cache=True, nodict=True):
        """
        Fill a dictionary containing rc property string names with the
        corresponding property. Defaults to **only** looking for cached
        (i.e. changed) properties.

        Parameters
        ----------
        props : dict-like
            Dictionary whose values are names of rc properties. The values
            are replaced with the corresponding property only if
            `~rc_configurator.__getitem__` does not return ``None``. Otherwise,
            that key, value pair is omitted from the output dictionary.
        cache : bool, optional
            Whether to look for all properties or just cached (i.e.
            recently changed) properties. In the later case, keys that
            return ``None`` will be omitted from the output dictionary.

        Returns
        -------
        dict
            Dictionary filled with rc properties.

        Note
        ----
        `~rc_configurator.fill` is used to build dictionaries for updating
        `~matplotlib.artist.Artist` instances. Of course, the artist property
        won't need updating unless an rc setting has changed since it was
        instantiated.

        With this in mind, `~rc_configurator.__getitem__` returns ``None``
        when a setting has not been changed (changed settings are cached
        in a separate dictionary). This prevents hundreds of 1000-item
        dictionary lookups. Without caching, runtime increases by 1s even
        for relatively simple plots.
        """
        props_out = {}
        self._no_dict = nodict
        if not cache:
            orig = self._getitem_mode
            self._getitem_mode = 0
        for key,value in props.items():
            item = self._get(value)
            if item is not None:
                props_out[key] = item
        self._no_dict = False
        if not cache:
            self._getitem_mode = orig
        return props_out

    def reset(self):
        """Restores settings to the default."""
        return self.__init__()

    def update(self, *args, **kwargs):
        """
        Bulk updates global settings. Usage is just like with python `dict`
        objects.

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
            raise ValueError('Accepts 1-2 positional arguments only. Use plot.rc.update(kw) to update a bunch of names, or plot.rc.update(category, kw) to update subcategories belonging to single category e.g. axes. All kwargs will be added to the dict.')
        for key,value in kw.items():
            self[prefix + key] = value

# Rc object
rc = rc_configurator()
"""Instance of `rc_configurator`. This is used to change rc settings.
See the `~proplot.rcmod` documentation for details."""

#------------------------------------------------------------------------------#
# Ipython notebook behavior
#------------------------------------------------------------------------------#
# Unbelievably weird problem:
#   * Apparently you can change the backend until the first plot is drawn, then
#     it stays the same. So the rcdefault() command changes the backend to a
#     non-inline version. See: https://stackoverflow.com/q/48320804/4970632
#
# Notes on python figures:
#   * Can use InlineBackend rc configuration to make inline figure properties
#     different from figure.<subproperty> settings in rcParams. Problem is,
#     whenever rcParams are reset/pyfuncs module is reloaded, the previous
#     InlineBackend properties disappear.
#   * It is *also* necessary to maintain separate savefig options, including 'transparent'
#     and 'facecolor' -- cannot just set these to use the figure properties.
#     If transparent set to False, saved figure will have no transparency *even if* the
#     default figure.facecolor has zero alpha. Will *only* be transparent if alpha explicitly
#     changed by user command. Try playing with settings in plot.globals to see.
#
# Notes on jupyter configuration:
#   * In .jupyter, the jupyter_nbconvert_config.json sets up locations of stuff; templates
#       for nbextensions and formatting files for markdown/code cells.
#   * In .jupyter, the jupyter_notebook_config.json installs the configurator extension
#       for managing extra plugins.
#   * In .jupyter, not sure yet how to successfully use jupyter_console_config.py and
#       jupyter_notebook_config.py; couldn't get it to do what this function does on startup.
#   * In .jupyter/custom, current_theme.txt lists the current jupyterthemes theme, custom.css
#       contains CSS formatting for it, and fonts should contain font files -- note that there
#       are not font files on my Mac, even though jupyterthemes works; sometimes may be empty
#   * In .jupyter/nbconfig, tree.json loads the extra tab for the NBconfigurator, and
#       common.json gives option to hide incompatible plugs, and notebook.json contains all
#       the new settings; just copy it over to current notebook to update
#------------------------------------------------------------------------------#
def nb_setup():
    """
    Optionally called on import, results in higher-quality iPython notebook
    inline figures. Also enables the useful `autoreload
    <https://ipython.readthedocs.io/en/stable/config/extensions/autoreload.html>`__
    and autosave extensions. The latter automatically saves your notebook
    every `autosave` seconds.

    To disable running this on import, use ``nbsetup: False`` in your
    ``.proplotrc`` file. See the `~proplot.rcmod` documentation for details.
    """
    # Make sure we are in session
    ipython = get_ipython() # save session
    if ipython is None:
        warnings.warn("ProPlot should be used within IPython.")
        return

    # Only do this if not already loaded -- otherwise will get *recursive* 
    # reloading, even with unload_ext command!
    if _rcGlobals['autoreload']:
        if 'autoreload' not in ipython.magics_manager.magics['line']:
            ipython.magic("reload_ext autoreload") # reload instead of load, to avoid annoying message
        ipython.magic("autoreload " + str(_rcGlobals['autoreload'])) # turn on expensive autoreloading

    # Initialize with default 'inline' settings
    # Reset rc object afterwards
    try:
        # For notebooks
        ipython.magic("matplotlib inline") # change print_figure_kwargs to see edges
        rc.reset()
    except Exception:
        # For ipython sessions
        ipython.magic("matplotlib qt") # change print_figure_kwargs to see edges
        rc.reset()
    else:
        # Autosaving, capture the annoying message + 2 line breaks
        if _rcGlobals['autosave']:
            with IPython.utils.io.capture_output() as _:
                ipython.magic("autosave " + str(_rcGlobals['autosave'])) # autosave every minute
        # Retina probably more space efficient (high-res bitmap), but svg is prettiest
        # and is only one preserving vector graphics
        ipython.magic("config InlineBackend.figure_formats = ['retina','svg']")
        # Control all settings with 'rc' object, *no* notebook-specific overrides
        ipython.magic("config InlineBackend.rc = {}")
        # Disable matplotlib tight layout, use proplot instead
        ipython.magic("config InlineBackend.print_figure_kwargs = {'bbox_inches':None}")
        # So don't have memory issues/have to keep re-closing them
        ipython.magic("config InlineBackend.close_figures = True")

# Run nbsetup
if _rcGlobals['nbsetup']:
    nb_setup()
