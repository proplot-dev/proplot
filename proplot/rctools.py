#!/usr/bin/env python3
"""
A special object named `~proplot.rctools.rc`, belonging to the
`~proplot.rctools.rc_configurator` class, is created on import.
This is your **one-stop shop for changing global settings** belonging to any of the
following three categories.

1. Builtin matplotlib `rcParams <https://matplotlib.org/users/customizing.html>`__
   settings. These have the format ``x.y`` or ``x.y.z``.
2. ProPlot :ref:`rcExtraParams` settings. These also have the format ``x.y``
   (see below).
3. ProPlot :ref:`rcGlobals` settings. These have no dots (see below). They are **simple,
   short** names used to change multiple matplotlib and ProPlot settings at once,
   as shorthands for settings with longer names, or for special options. For example,
   ``ticklen`` changes the tick length for the *x* and *y* axes in one go.

You can change settings with the `~proplot.rctools.rc` object as follows.

* ``plot.rc.name = value``
* ``plot.rc['name'] = value``
* ``plot.rc.update(name1=value1, name2=value2)``
* ``plot.rc.update({'name1':value1, 'name2':value2})``

To temporarily change settings on a particular axes, use either of the following.

* ``ax.format(name=value)``
* ``ax.format(rc_kw={'name':value})``

In all of these examples, if the setting name ``name`` contains
any dots, you can simply **omit the dots**. For example, to change the
:ref:`rcExtraParams` property ``title.loc``, use ``plot.rc.titleloc = value``,
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
``nbsetup``         Whether to run `setup` on import.
``backend``         The matplotlib backend, passed to `%matplotlib <https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-matplotlib>`__.
``autosave``        If non-empty and ``nbsetup`` is ``True``, passed to `%autosave <https://www.webucator.com/blog/2016/03/change-default-autosave-interval-in-ipython-notebook/>`__.
``autoreload``      If non-empty and ``nbsetup`` is ``True``, passed to `%autoreload <https://ipython.readthedocs.io/en/stable/config/extensions/autoreload.html#magic-autoreload>`__.
``abc``             Boolean, indicates whether to draw a-b-c labels by default.
``tight``           Boolean, indicates whether to auto-adjust figure bounds and subplot spacings.
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
rcExtraParams
#############
For a-b-c labels and axes title settings, use the ``abc`` and ``title``
categories. For axis tick label settings, use the ``tick`` category.
For figure title, row label, and column label settings, use the new
``suptitle``, ``rowlabel``, and ``collabel`` categories.

==============================================================  ================================================================================================================================================================================================
Key(s)                                                          Description
==============================================================  ================================================================================================================================================================================================
``abc.format``                                                  a-b-c label format (for options, see `~proplot.axes.BaseAxes.format_partial`).
``abc.loc``, ``title.loc``                                      a-b-c label or title position (for options, see `~proplot.axes.BaseAxes.format_partial`).
``abc.border``, ``title.border``                                Boolean, indicates whether to draw a white border around a-b-c labels or titles located inside an axes.
``abc.linewidth``, ``title.linewidth``                          Width of the white border around a-b-c labels or titles.
``abc.color``, ``abc.fontsize``, ``abc.weight``                 Font color, size, and weight for a-b-c labels.
``title.color``, ``title.fontsize``, ``title.weight``           Font color, size, and weight for subplot titles.
``rowlabel.color``, ``rowlabel.fontsize``, ``rowlabel.weight``  Font color, size, and weight for row labels.
``collabel.color``, ``collabel.fontsize``, ``collabel.weight``  Font color, size, and weight for column labels.
``suptitle.color``, ``suptitle.fontsize``, ``suptitle.weight``  Font color, size, and weight for the figure title.
``tick.labelcolor``, ``tick.labelsize``, ``tick.labelweight``   Font color, size, and weight for axis tick labels. These mirror the ``axes.labelcolor``, ``axes.labelsize``, and ``axes.labelweight`` `~matplotlib.rcParams` settings used for axes labels.
``axes.formatter.zerotrim``                                     Boolean, indicates whether trailing decimal zeros are trimmed on tick labels.
``axes.formatter.timerotation``                                 Float, indicates the default *x* axis tick label rotation for datetime tick labels.
==============================================================  ================================================================================================================================================================================================

For minor gridlines, use the ``gridminor`` category. For meridian and parallel
gridlines on `~proplot.axes.MapAxes`, use the ``geogrid`` category.
If a property is empty, the corresponding property from the buildin ``grid``
category is used.

==============================================  ===============================================================================
Key                                             Description
==============================================  ===============================================================================
``gridminor.linewidth``, ``geogrid.linewidth``  The line width.
``gridminor.linestyle``, ``geogrid.linestyle``  The line style.
``gridminor.alpha``, ``geogrid.alpha``          The line transparency.
``gridminor.color``, ``geogrid.color``          The line color.
``geogrid.labels``                              Boolean, indicates whether to label the parallels and meridians.
``geogrid.labelsize``                           Font size for latitide and longitude labels. Inherits from ``small``.
``geogrid.latmax``                              Absolute latitude in degrees, poleward of which meridian gridlines are cut off.
``geogrid.lonstep``, ``geogrid.latstep``        Interval for meridian and parallel gridlines, in degrees.
==============================================  ===============================================================================

The below properties are particular to `~proplot.axes.MapAxes`. The properties
for geographic elements like ``land`` are used when the corresponding
:ref:`rcGlobals` geographic feature is toggled on.

=======================================================  =============================================================================
Key                                                      Description
=======================================================  =============================================================================
``map.facecolor``, ``map.edgecolor``, ``map.linewidth``  Face color, edge color, and edge width for the map outline patch.
``land.color``, ``ocean.color``, ``lakes.color``         Face color for land, ocean, and lake patches.
``rivers.color``, ``rivers.linewidth``                   Line color and linewidth for river lines.
``borders.color``, ``borders.linewidth``                 Line color and linewidth for country border lines.
``innerborders.color``, ``innerborders.linewidth``       Line color and linewidth for internal border lines.
=======================================================  =============================================================================

There are two new additions to pre-existing ``image`` category, described
in the table below.

=================  =========================================================================================================================================================================================================================================================
Key                Description
=================  =========================================================================================================================================================================================================================================================
``image.levels``   Default number of levels for ``pcolormesh`` and ``contourf`` plots.
``image.edgefix``  Whether to fix the the `white-lines-between-filled-contours <https://stackoverflow.com/q/8263769/4970632>`__ and `white-lines-between-pcolor-rectangles <https://stackoverflow.com/q/27092991/4970632>`__ issues. This slows down figure rendering a bit.
=================  =========================================================================================================================================================================================================================================================

The ``colorbar`` category has settings that control the default colorbar
layout. It is analogous to the builtin ``legend`` category, but configures
both *inset* and *panel* colorbars (see the `~proplot.axes.BaseAxes`
`~proplot.axes.BaseAxes.colorbar` and `~proplot.axes.PanelAxes`
`~proplot.axes.PanelAxes.colorbar` methods for details).

========================  =========================================================================================================================
Key                       Description
========================  =========================================================================================================================
``colorbar.grid``         Boolean, indicates whether to draw "gridlines" between each level of the colorbar.
``colorbar.frameon``      Boolean, indicates whether to draw a frame behind inset colorbars.
``colorbar.framealpha``   The opacity of colorbar frame.
``colorbar.loc``          Default inset colorbar location, one of "upper right", "upper left", "lower left", or "lower right", or "center" options.
``colorbar.length``       Default length of inset colorbars.
``colorbar.width``        Default width of inset colorbars.
``colorbar.extend``       Length of rectangular or triangular "extensions" for panel colorbars.
``colorbar.extendinset``  Length of rectangular or triangular "extensions" for inset colorbars.
``colorbar.pad``          Default padding between figure edge of rectangular or triangular "extensions" for inset colorbars.
========================  =========================================================================================================================

The ``subplot`` category has settings that control the default figure
layout. As with all sizing arguments, if specified as a number, the units
are inches. If string, the units are interpreted by `~proplot.utils.units`.

=======================  ==================================================================
Key                      Description
=======================  ==================================================================
``subplots.axwidth``      Default width of each axes.
``subplots.panelwidth``   Width of side panels.
``subplots.cbarwidth``    Width of "colorbar" panels.
``subplots.legwidth``     Width of "legend" panels.
``subplots.borderpad``    Padding around edges for tight subplots.
``subplots.subplotpad``   Padding between main subplots.
``subplots.panelpad``     Padding between axes panels and their parents.
``subplots.titlespace``   Vertical space for titles.
``subplots.ylabspace``    Horizontal space between subplots alotted for *y*-labels.
``subplots.xlabspace``    Vertical space between subplots alotted for *x*-labels.
``subplots.nolabspace``   Space between subplots alotted for tick marks.
``subplots.innerspace``   Purely empty space between subplots.
``subplots.panelspace``   Purely empty space between main axes and side panels.
=======================  ==================================================================
"""
# TODO: Add 'style' setting that overrides .proplotrc
import re
import os
import sys
import yaml
import cycler
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
import matplotlib.pyplot as plt
from . import utils
try:
    import IPython
    get_ipython = IPython.get_ipython
except ModuleNotFoundError:
    get_ipython = lambda: None
__all__ = ['rc', 'rc_configurator', 'nb_setup']

# Initialize
from matplotlib import rcParams as _rcParams
_rcGlobals = {}
_rcExtraParams = {}

# Get default font
# WARNING: Had issues with Helvetica Neue on Linux, weirdly some characters
# failed to render/printed nonsense, but Helvetica was fine.
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
    'backend':      [],
    'autosave':     [],
    'autoreload':   [],
    # Subplots
    'tight':        [],
    'abc':          [],
    # Style
    'fontname':     ['font.family'], # specify family directly, so we can easily switch between serif/sans-serif; requires text.usetex = False; see below
    'cmap':         ['image.cmap'],
    'lut':          ['image.lut'],
    'cycle':        [],
    'rgbcycle':     [],
    'alpha':        ['axes.alpha'], # this is a custom setting
    'color':        ['axes.edgecolor', 'map.edgecolor', 'axes.labelcolor', 'tick.labelcolor', 'hatch.color', 'xtick.color', 'ytick.color'], # change the 'color' of an axes
    'facecolor':    ['axes.facecolor', 'map.facecolor'],
    'small':        ['font.size', 'tick.labelsize', 'xtick.labelsize', 'ytick.labelsize', 'axes.labelsize', 'legend.fontsize', 'geogrid.labelsize'], # the 'small' fonts
    'large':        ['abc.fontsize', 'figure.titlesize', 'axes.titlesize', 'suptitle.fontsize', 'title.fontsize', 'rowlabel.fontsize', 'collabel.fontsize'], # the 'large' fonts
    'linewidth':    ['axes.linewidth', 'map.linewidth', 'hatch.linewidth', 'xtick.major.width', 'ytick.major.width'],
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
_rc_names = {*_rcParams.keys()}
_rc_names_global = {*_rcGlobals_children.keys()}
_rc_names_custom = {
    'axes.formatter.zerotrim', 'axes.formatter.timerotation',
    'axes.gridminor', 'axes.geogrid', 'axes.alpha',
    'image.levels', 'image.edgefix',
    'map.linewidth', 'map.facecolor', 'map.edgecolor',
    'land.color', 'ocean.color', 'lakes.color', 'coast.color', 'coast.linewidth',
    'borders.color', 'borders.linewidth', 'innerborders.color', 'innerborders.linewidth', 'rivers.color', 'rivers.linewidth',
    'abc.fontsize', 'abc.weight', 'abc.color', 'abc.loc', 'abc.format', 'abc.border', 'abc.linewidth',
    'title.loc', 'title.color', 'title.border', 'title.linewidth', 'title.weight', 'title.fontsize',
    'suptitle.fontsize', 'suptitle.weight', 'suptitle.color',
    'rowlabel.fontsize', 'rowlabel.weight', 'rowlabel.color',
    'collabel.fontsize', 'collabel.weight', 'collabel.color',
    'gridminor.alpha', 'gridminor.color', 'gridminor.linestyle', 'gridminor.linewidth',
    'geogrid.labels', 'geogrid.alpha', 'geogrid.color', 'geogrid.labelsize', 'geogrid.linewidth', 'geogrid.linestyle', 'geogrid.latmax', 'geogrid.lonstep', 'geogrid.latstep',
    'subplots.subplotpad', 'subplots.panelpad', 'subplots.borderpad', 'subplots.titlespace',
    'subplots.innerspace',
    'tick.labelweight', 'tick.labelcolor', 'tick.labelsize',
    'subplots.legwidth', 'subplots.cbarwidth', 'subplots.ylabspace', 'subplots.xlabspace', 'subplots.nolabspace',
    'subplots.axwidth', 'subplots.panelwidth', 'subplots.panelspace',
    'colorbar.grid', 'colorbar.frameon', 'colorbar.framealpha', 'colorbar.length', 'colorbar.width', 'colorbar.loc', 'colorbar.extend', 'colorbar.extendinset', 'colorbar.axespad',
    }
# Used by BaseAxes.format, allows user to pass rc settings as keyword args,
# way less verbose. For example, landcolor='b' vs. rc_kw={'land.color':'b'}.
_rc_names_nodots = { # useful for passing these as kwargs
    name.replace('.', ''):name for names in
    (_rc_names_custom, _rc_names, _rc_names_global)
    for name in names
    }
# Categories for returning dict of subcategory properties
_rc_categories = {
    *(re.sub('\.[^.]*$', '', name) for names in (_rc_names_custom, _rc_names) for name in names),
    *(re.sub('\..*$', '', name) for names in (_rc_names_custom, _rc_names) for name in names)
    }

#-------------------------------------------------------------------------------
# Contextual settings management
# Adapted from seaborn; see: https://github.com/mwaskom/seaborn/blob/master/seaborn/rcmod.py
#-------------------------------------------------------------------------------
def _to_points(value, key=None):
    """Converts certain keys to the units "points". If "key" is passed, tests
    that key against possible keys that accept physical units."""
    # See: https://matplotlib.org/users/customizing.html, all props matching
    # the below strings use the units 'points'!
    if isinstance(value,str) and (
        key is None or
        (isinstance(key,str) and any(s in key for s in ('width','size','pad')))
        ):
        value = utils.units(value, 'pt')
    return value

def _set_cycler(name):
    """Sets the default color cycler."""
    # Draw from dictionary
    try:
        colors = mcm.cmap_d[name].colors
    except Exception:
        cycles = sorted(name for name,cmap in mcm.cmap_d.items() if isinstance(cmap, mcolors.ListedColormap))
        raise ValueError(f'Invalid cycle name "{name}". Options are: {", ".join(cycles)}')
    # Apply color name definitions
    if _rcGlobals['rgbcycle'] and name.lower()=='colorblind':
        regcolors = colors + [(0.1, 0.1, 0.1)]
    elif mcolors.to_rgb('r') != (1.0,0.0,0.0): # reset
        regcolors = [(0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.75, 0.75, 0.0), (0.75, 0.75, 0.0), (0.0, 0.75, 0.75), (0.0, 0.0, 0.0)]
    else:
        regcolors = [] # no reset necessary
    for code,color in zip('brgmyck', regcolors):
        rgb = mcolors.to_rgb(color)
        mcolors.colorConverter.colors[code] = rgb
        mcolors.colorConverter.cache[code]  = rgb
    # Pass to cycle constructor
    _rcParams['patch.facecolor'] = colors[0]
    _rcParams['axes.prop_cycle'] = cycler.cycler('color', colors)

def _get_globals(key=None, value=None):
    """Returns dictionaries for updating "child" properties in
    `rcParams` and `rcExtraParams` with global property."""
    kw = {} # builtin properties that global setting applies to
    kw_custom = {} # custom properties that global setting applies to
    if key is not None and value is not None:
        items = [(key,value)]
    else:
        items = _rcGlobals.items()
    for key,value in items:
        # Tick length/major-minor tick length ratio
        if key in ('ticklen', 'ticklenratio'):
            if key=='ticklen':
                ticklen = _to_points(value)
                ratio = _rcGlobals['ticklenratio']
            else:
                ticklen = _rcGlobals['ticklen']
                ratio = value
            kw['xtick.minor.size'] = ticklen*ratio
            kw['ytick.minor.size'] = ticklen*ratio
        # Spine width/major-minor tick width ratio
        elif key in ('linewidth', 'tickratio'):
            if key=='linewidth':
                tickwidth = _to_points(value)
                ratio = _rcGlobals['tickratio']
            else:
                tickwidth = _rcGlobals['linewidth']
                ratio = value
            kw['xtick.minor.width'] = tickwidth*ratio
            kw['ytick.minor.width'] = tickwidth*ratio
        # Grid line
        # NOTE: Changing minor width does not affect major width
        elif key in ('grid.linewidth', 'gridratio'):
            if key=='grid.linewidth':
                gridwidth = _to_points(value)
                ratio = _rcGlobals['gridratio']
            else:
                gridwidth = _rcParams['grid.linewidth']
                ratio = value
            kw_custom['gridminor.linewidth'] = gridwidth*ratio
        # Now update linked settings
        if key not in ('gridratio','tickratio','ticklenratio'):
            points = (key.split('.')[0] not in ('colorbar','subplot'))
            for name in _rcGlobals_children[key]:
                if points:
                    value = _to_points(value, key=name)
                if name in _rcExtraParams:
                    kw_custom[name] = value
                else:
                    kw[name] = value
            if key=='linewidth' and value==0:
                ikw, ikw_custom = _get_globals('ticklen', 0)
                kw.update(ikw)
                kw_custom.update(ikw_custom)
    return kw, kw_custom

class rc_configurator(object):
    _public_api = ('get', 'fill', 'category', 'reset', 'context', 'update') # getattr and setattr will not look for these items on underlying dictionary
    # @_counter # about 0.05s
    def __init__(self):
        """Magical abstract class for managing matplotlib `rcParams
        <https://matplotlib.org/users/customizing.html>`__ settings, ProPlot
        :ref:`rcExtraParams` settings, and :ref:`rcGlobals` "global" settings.
        See the `~proplot.rctools` documentation for details."""
        # Set the default style. Note that after first figure made, backend
        # is 'sticky', never changes! See: https://stackoverflow.com/a/48322150/4970632
        # TODO: Flexible styles, and register proplot as a style!
        plt.style.use('default')

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
            # Add keys to dictionaries
            gkeys, ckeys = {*()}, {*()}
            for key,value in data.items():
                if key in _rc_names_global:
                    _rcGlobals[key] = value
                    if i==0:
                        gkeys.add(key)
                else:
                    points = (key.split('.')[0] not in ('colorbar','subplot'))
                    if points:
                        value = _to_points(value, key=key)
                    if key in _rc_names_custom:
                        _rcExtraParams[key] = value
                        if i==0:
                            ckeys.add(key)
                    else:
                        try:
                            _rcParams[key] = value
                        except KeyError:
                            raise RuntimeError(('Default', 'User')[i] + f' .proplotrc file has invalid key "{key}".')
            # Make sure we did not miss anything
            if i==0:
                if gkeys!=_rc_names_global:
                    raise RuntimeError(f'Default .proplotrc file has incomplete or invalid global keys {_rc_names_global - gkeys}.')
                if ckeys!=_rc_names_custom:
                    raise RuntimeError(f'Default .proplotrc file has incomplete or invalid custom keys {_rc_names_custom - ckeys}.')

        # Set default fontname and cycler
        _set_cycler(_rcGlobals['cycle'])
        if _rcGlobals.get('fontname', None) is None:
            _rcGlobals['fontname'] = _default_font

        # Apply *global settings* to children settings
        rc, rc_new = _get_globals()
        for key,value in rc.items():
            _rcParams[key] = value
        for key,value in rc_new.items():
            _rcExtraParams[key] = value

        # Caching stuff
        self._init = True
        self._getitem_mode = 0
        self._context = {}
        self._cache = {}
        self._cache_orig = {}
        self._cache_restore = {}

    def __getitem__(self, key):
        """Returns `rcParams <https://matplotlib.org/users/customizing.html>`__,
        :ref:`rcExtraParams`, and :ref:`rcGlobals` settings. If we are in a
        `~rc_configurator.context` block, may return ``None`` if the setting
        is not cached (i.e. if it was not changed by the user).
        """
        # Can get a whole bunch of different things
        # Get full dictionary e.g. for rc[None]
        if not key:
            return {**_rcParams, **_rcExtraParams}

        # Standardize
        # NOTE: If key is invalid, raise error down the line.
        if '.' not in key and key not in _rcGlobals:
            key = _rc_names_nodots.get(key, key)

        # Allow for special time-saving modes where we *ignore _rcParams*
        # or even *ignore _rcExtraParams*.
        mode = self._getitem_mode
        if mode==0:
            kws = (self._cache, _rcGlobals, _rcExtraParams, _rcParams)
        elif mode==1:
            kws = (self._cache, _rcGlobals, _rcExtraParams) # custom only!
        elif mode==2:
            kws = (self._cache,) # changed only!
        else:
            raise KeyError(f'Invalid caching mode {mode}.')

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
        :ref:`rcExtraParams`, and :ref:`rcGlobals` settings."""
        # Check whether we are in context block
        # NOTE: Do not add key to cache until we are sure it is a valid key
        cache = self._cache
        context = bool(self._context) # test if context dict is non-empty
        if context:
            restore = self._cache_restore

        # Standardize
        # NOTE: If key is invalid, raise error down the line.
        if '.' not in key and key not in _rcGlobals:
            key = _rc_names_nodots.get(key, key)

        # Set the default cycler
        # If rgbcycle is changed, change it, and re-apply the cycler
        if key=='rgbcycle':
            cache[key] = value
            _rcGlobals[key] = value
            key, value = 'cycle', _rcGlobals['cycle']
        if key=='cycle':
            cache[key] = value
            if context:
                restore[key] = _rcGlobals[key]
                restore['axes.prop_cycle'] = _rcParams['axes.prop_cycle']
                restore['patch.facecolor'] = _rcParams['patch.facecolor']
            _set_cycler(value)

        # Gridline toggling, complicated because of the clunky way this is
        # implemented in matplotlib. There should be a gridminor setting!
        elif key in ('grid', 'gridminor'):
            cache[key] = value
            ovalue = _rcParams['axes.grid']
            owhich = _rcParams['axes.grid.which']
            if context:
                restore[key] = _rcGlobals[key]
                restore['axes.grid'] = ovalue
                restore['axes.grid.which'] = owhich
            # Instruction is to turn off gridlines
            if not value:
                # Gridlines are already off, or they are on for the particular
                # ones that we want to turn off. Instruct to turn both off.
                if not ovalue or (key=='grid' and owhich=='major') or (key=='gridminor' and owhich=='minor'):
                    which = 'both' # disable both sides
                # Gridlines are currently on for major and minor ticks, so we instruct
                # to turn on gridlines for the one we *don't* want off
                elif owhich=='both': # and ovalue is True, as we already tested
                    value = True
                    which = 'major' if key=='gridminor' else 'minor' # if gridminor=False, enable major, and vice versa
                # Gridlines are on for the ones that we *didn't* instruct to turn
                # off, and off for the ones we do want to turn off. This just
                # re-asserts the ones that are already on.
                else:
                    value = True
                    which = owhich
            # Instruction is to turn on gridlines
            else:
                # Gridlines are already both on, or they are off only for the ones
                # that we want to turn on. Turn on gridlines for both.
                if owhich=='both' or (key=='grid' and owhich=='minor') or (key=='gridminor' and owhich=='major'):
                    which = 'both'
                # Gridlines are off for both, or off for the ones that we
                # don't want to turn on. We can just turn on these ones.
                else:
                    which = owhich
            cache.update({'axes.grid':value, 'axes.grid.which':which})
            _rcParams.update({'axes.grid':value, 'axes.grid.which':which})

        # Ordinary settings
        elif key in _rcGlobals:
            # Update global setting
            cache[key] = value
            if context:
                restore[key] = _rcGlobals[key]
            _rcGlobals[key] = value
            # Update children of setting
            rc, rc_new = _get_globals(key, value)
            cache.update(rc)
            cache.update(rc_new)
            if context:
                restore.update({key:_rcParams[key] for key in rc})
                restore.update({key:_rcExtraParams[key] for key in rc_new})
            _rcParams.update(rc)
            _rcExtraParams.update(rc_new)
        else:
            points = (key.split('.')[0] not in ('colorbar','subplot'))
            if points:
                value = _to_points(value, key=key)
            if key in _rc_names_custom:
                cache[key] = value
                if context:
                    restore[key] = _rcExtraParams[key]
                _rcExtraParams[key] = value
            elif key in _rc_names:
                cache[key] = value
                if context:
                    restore[key] = _rcParams[key]
                _rcParams[key] = value # rcParams dict has key validation
            else:
                raise KeyError(f'Invalid key "{key}".')
        self._init = False # setitem was successful, we are no longer in initial state

    def __getattribute__(self, attr):
        """Invokes `~rc_configurator.__getitem__`."""
        if attr[:1]=='_' or attr in self._public_api:
            return object.__getattribute__(self, attr)
        else:
            return self[attr]

    def __setattr__(self, attr, value):
        """Invokes `~rc_configurator.__setitem__`."""
        if attr[:1]=='_':
            object.__setattr__(self, attr, value)
        else:
            self[attr] = value

    def __enter__(self):
        """Apply temporary user settings."""
        self._cache_orig = rc._cache.copy()
        for key,value in self._context.items():
            self[key] = value # applies globally linked and individual settings

    def __exit__(self, _type, _value, _traceback):
        """Restore configurator cache to initial state."""
        self._context = {}
        self._getitem_mode = 0
        for key,value in self._cache_restore.items():
            self[key] = value
        self._cache = self._cache_orig
        self._cache_restore = {}
        self._cache_orig = {}

    def __delitem__(self, *args):
        """Pseudo-immutability."""
        raise NotImplementedError

    def __delattr__(self, *args):
        """Pseudo-immutability."""
        raise NotImplementedError

    def __repr__(self):
        """Nice string representation."""
        length = 1 + max(len(key) for key in _rcGlobals.keys())
        string = 'rc = {\n' + '\n'.join(f'  {key}: {" "*(length-len(key))}{value}'
                                    for key,value in _rcGlobals.items()) + '\n}'
        return string

    def __str__(self):
        """Short string representation."""
        string = ', '.join(f'{key}: {value}' for key,value in _rcGlobals.items())
        return string

    # Internally used, but public methods.
    def context(self, *args, mode=0, **kwargs):
        """
        Temporarily modifies settings in a ``with...as`` block,
        used by ProPlot internally but may also be useful for power users.

        This function was invented to prevent successive calls to
        `~proplot.axes.BaseAxes.format` from constantly looking up and re-applying
        unchanged settings.  Testing showed that these gratuitous
        `rcParams <https://matplotlib.org/users/customizing.html>`__
        lookups and artist updates tended to increase runtime by seconds,
        even for relatively simple plots.

        Parameters
        ----------
        *args
            Dictionaries of setting names and values.
        **kwargs
            Setting names and values passed as keyword arguments.

        Other parameters
        ----------------
        mode : {0,1,2}, optional
            The `~rc_configurator.__getitem__` mode.
            Dictates the behavior of the `rc` object within a ``with...as``
            block when settings are requested with e.g. ``rc['setting']``. If
            you are using `~rc_configurator.context` manually, the `mode` is
            automatically set to ``0`` -- other input is ignored. Internally,
            ProPlot uses all of the three available modes.

            0. All settings (`rcParams <https://matplotlib.org/users/customizing.html>`__,
               :ref:`rcExtraParams`, and :ref:`rcGlobals`) are returned, whether
               or not `~rc_configurator.context` has changed them.
            1. Unchanged `rcParams <https://matplotlib.org/users/customizing.html>`__
               return ``None``. :ref:`rcExtraParams` and :ref:`rcGlobals` are
               returned whether or not `~rc_configurator.context` has changed them.
               This is used in the initial `~proplot.axes.BaseAxes.__init__`
               call to `~proplot.axes.BaseAxes.format`. When a setting lookup
               returns ``None``, `~proplot.axes.BaseAxes.format` does not apply it.
            2. All unchanged settings return ``None``. This is used during user
               calls to `~proplot.axes.BaseAxes.format`.

        Example
        -------
        .. code-block:: python

            import proplot as plot
            with plot.rc.context(linewidth=2, ticklen=5):
                f, ax = plot.subplots()
                ax.plot(data)

        """
        # Apply mode
        if mode not in range(3):
            raise ValueError(f'Invalid _getitem_mode {mode}.')
        for arg in args:
            if not isinstance(arg, dict):
                raise ValueError('Non-dictionary argument.')
            kwargs.update(arg)
        self._context = kwargs # could be empty
        self._getitem_mode = mode
        return self

    def get(self, key, cache=False):
        """
        Returns a setting.

        Parameters
        ----------
        key : str
            The setting name.
        cache : bool, optional
            If ``False``, the `~rc_configurator.__getitem__` mode is temporarily
            set to ``0`` (see `~rc_configurator.context`).
        """
        if not cache:
            orig = self._getitem_mode
            self._getitem_mode = 0
        item = self[key]
        if not cache:
            self._getitem_mode = orig
        return item

    def fill(self, props, cache=True):
        """
        Returns a dictionary filled with `rc` settings, used internally to build
        dictionaries for updating `~matplotlib.artist.Artist` instances.

        Parameters
        ----------
        props : dict-like
            Dictionary whose values are names of `rc` settings. The values
            are replaced with the corresponding property only if
            `~rc_configurator.__getitem__` does not return ``None``. Otherwise,
            that key, value pair is omitted from the output dictionary.
        cache : bool, optional
            If ``False``, the `~rc_configurator.__getitem__` mode is temporarily
            set to ``0`` (see `~rc_configurator.context`). Otherwise, if an `rc`
            lookup returns ``None``, the setting is omitted from the output
            dictionary.
        """
        if not cache:
            orig = self._getitem_mode
            self._getitem_mode = 0
        props_out = {}
        for key,value in props.items():
            item = self[value]
            if item is not None:
                props_out[key] = item
        if not cache:
            self._getitem_mode = orig
        return props_out

    def category(self, cat, cache=True):
        """
        Returns a dictionary of settings belonging to the indicated category,
        i.e. settings beginning with the substring ``cat + '.'``.

        Parameters
        ----------
        cat : str, optional
            The `rc` settings category.
        cache : bool, optional
            If ``False``, the `~rc_configurator.__getitem__` mode is temporarily
            set to ``0`` (see `~rc_configurator.context`).
        """
        # Check
        if cat not in _rc_categories:
            raise ValueError(f'RC category {cat} does not exist. Valid categories are {", ".join(_rc_categories)}.')
        if not cache:
            mode = 0
        else:
            mode = self._getitem_mode

        # Allow for special time-saving modes where we *ignore _rcParams*
        # or even *ignore _rcExtraParams*.
        if mode==0:
            kws = (self._cache, _rcGlobals, _rcExtraParams, _rcParams)
        elif mode==1:
            kws = (self._cache, _rcGlobals, _rcExtraParams)
        elif mode==2:
            kws = (self._cache, _rcGlobals)
        else:
            raise KeyError(f'Invalid caching mode {mode}.')

        # Return params dictionary
        params = {}
        for kw in kws:
            for category,value in kw.items():
                if re.search(f'^{cat}\.', category):
                    subcategory = re.sub(f'^{cat}\.', '', category)
                    if subcategory and '.' not in subcategory:
                        params[subcategory] = value
        return params

    def update(self, *args, **kwargs):
        """
        Bulk updates settings, usage is similar to python `dict` objects.

        Parameters
        ----------
        *args : str, dict, or (str, dict)
            Positional arguments can be a dictionary of `rc` settings and/or
            a "category" string name. If a category name is passed, all settings
            in the dictionary (if it was passed) and all keyword arg names
            (if they were passed) are prepended with the string ``cat + '.'``.
            For example, ``rc.update('axes', labelsize=20, titlesize=20)``
            changes the ``axes.labelsize`` and ``axes.titlesize`` properties.
        **kwargs
            `rc` settings passed as keyword args.
        """
        # Parse args
        kw = {}
        prefix = ''
        if len(args)>2:
            raise ValueError('Accepts 1-2 positional arguments. Use plot.rc.update(kw) to update a bunch of names, or plot.rc.update(category, kw) to update subcategories belonging to single category e.g. axes. All kwargs will be added to the dict.')
        elif len(args)==2:
            prefix = args[0]
            kw = args[1]
        elif len(args)==1:
            if isinstance(args[0], str):
                prefix = args[0]
            else:
                kw = args[0]
        # Apply settings
        if prefix:
            prefix = prefix + '.'
        kw.update(kwargs)
        for key,value in kw.items():
            self[prefix + key] = value

    def reset(self):
        """Restores settings to the default."""
        if not self._init: # save resources if rc is unchanged!
            return self.__init__()

# Declare rc object
# WARNING: Must be instantiated after ipython notebook setup! The default
# backend may change some rc settings!
rc = rc_configurator()
"""Instance of `rc_configurator`. This is used to change global settings.
See the `~proplot.rctools` documentation for details."""

# Ipython notebook behavior
def nb_setup():
    """
    Sets up your iPython workspace, called on import if ``rc['nbsetup']`` is
    ``True`` (the default). For all iPython sessions, passes ``rc['autoreload']``
    to the useful `autoreload <https://ipython.readthedocs.io/en/stable/config/extensions/autoreload.html>`__
    extension. For iPython *notebook* sessions, results in higher-quality inline figures
    and passes ``rc['autosave']`` to the `autosave <https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-matplotlib>`__
    extension.

    See the `~proplot.rctools` documentation for details.
    """
    # Make sure we are in session
    ipython = get_ipython() # save session
    if ipython is None:
        return

    # Only do this if not already loaded -- otherwise will get *recursive*
    # reloading, even with unload_ext command!
    if _rcGlobals['autoreload']:
        if 'autoreload' not in ipython.magics_manager.magics['line']:
            ipython.magic("reload_ext autoreload") # reload instead of load, to avoid annoying message
        ipython.magic("autoreload " + str(_rcGlobals['autoreload'])) # turn on expensive autoreloading

    # Initialize with default 'inline' settings
    # Reset rc object afterwards
    rc._init = False
    try:
        # For notebooks
        ipython.magic("matplotlib inline") # change print_figure_kwargs to see edges
        rc.reset()
    except Exception:
        # For ipython sessions
        ipython.magic("matplotlib qt")
        rc.reset()
    else:
        # Autosaving, capture the annoying message + 2 line breaks
        if _rcGlobals['autosave']:
            with IPython.utils.io.capture_output():
                ipython.magic("autosave " + str(_rcGlobals['autosave'])) # autosave every minute
        # Retina probably more space efficient (high-res bitmap), but svg is
        # prettiest and is only one preserving vector graphics
        ipython.magic("config InlineBackend.figure_formats = ['svg']")
        # Control all settings with 'rc' object, *no* notebook-specific overrides
        ipython.magic("config InlineBackend.rc = {}")
        # So don't have memory issues/have to keep re-closing them
        ipython.magic("config InlineBackend.close_figures = True")
        # Disable matplotlib tight layout, use proplot instead
        ipython.magic("config InlineBackend.print_figure_kwargs = {'bbox_inches':None}")

# Setup notebook
if _rcGlobals['nbsetup']:
    nb_setup()

