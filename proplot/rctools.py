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

#########
rcGlobals
#########

================  ====================================================================================================================================================================================================================
Key               Description
================  ====================================================================================================================================================================================================================
``nbsetup``       Whether to run `setup` on import. Can only be changed from the ``~/.proplotrc`` file.
``format``        The inline backend figure format, one of ``retina``, ``png``, ``jpeg``, ``pdf``, or ``svg``. Can only be changed from the ``~/.proplotrc`` file.
``autosave``      If non-empty and ``nbsetup`` is ``True``, passed to `%autosave <https://www.webucator.com/blog/2016/03/change-default-autosave-interval-in-ipython-notebook/>`__. Can only be changed from the ``~/.proplotrc`` file.
``autoreload``    If non-empty and ``nbsetup`` is ``True``, passed to `%autoreload <https://ipython.readthedocs.io/en/stable/config/extensions/autoreload.html#magic-autoreload>`__. Can only be changed from the ``~/.proplotrc`` file.
``abc``           Boolean, indicates whether to draw a-b-c labels by default.
``tight``         Boolean, indicates whether to auto-adjust figure bounds and subplot spacings.
``share``         The axis sharing level, one of ``0``, ``1``, ``2``, or ``3``.
``align``         Whether to align axis labels during draw; see `aligning labels <https://matplotlib.org/3.1.1/gallery/subplots_axes_and_figures/align_labels_demo.html>`__.
``span``          Boolean, toggles spanning axis labels.
``fontname``      Name of font used for all text in the figure. The default is ``Helvetica`` for Linux and ``Helvetica Neue`` for Windows/Mac. See `~proplot.fonttools` for details.
``cmap``          The default colormap.
``lut``           The number of colors to put in the colormap lookup table.
``cycle``         The default color cycle name, used e.g. for lines.
``rgbcycle``      Whether to register cycles names as ``'r'``, ``'b'``, ``'g'``, etc., like in `seaborn <https://seaborn.pydata.org/tutorial/color_palettes.html>`__.
``color``         The color of axis spines, tick marks, tick labels, and labels.
``alpha``         The opacity of the background axes patch.
``facecolor``     The color of the background axes patch.
``small``         Font size for legend text, tick labels, axis labels, and text generated with `~matplotlib.axes.Axes.text`.
``large``         Font size for titles, "super" titles, and a-b-c subplot labels.
``linewidth``     Thickness of axes spines and major tick lines.
``margin``        The margin of space between axes edges and objects plotted inside the axes, if ``xlim`` and ``ylim`` are unset.
``ticklen``       Length of major ticks in points.
``tickdir``       Major and minor tick direction. Must be one of ``out``, ``in``, or ``inout``.
``tickpad``       Padding between ticks and tick labels in points.
``grid``          Boolean, toggles major grid lines on and off.
``gridminor``     Boolean, toggles minor grid lines on and off.
``tickratio``     Ratio of minor tickline width to major tickline width.
``gridratio``     Ratio of minor gridline width to major gridline width.
``ticklenratio``  Ratio of minor tickline length to major tickline length.
``reso``          Resolution of geographic features, one of ``'lo'``, ``'med'``, or ``'hi'``
``geogrid``       Boolean, toggles meridian and parallel gridlines on and off.
``land``          Boolean, toggles land patches on and off.
``ocean``         Boolean, toggles ocean patches on and off.
``lakes``         Boolean, toggles lake patches on and off.
``coast``         Boolean, toggles coastline lines on and off.
``borders``       Boolean, toggles country border lines on and off.
``innerborders``  Boolean, toggles internal border lines on and off, e.g. for states and provinces.
``rivers``        Boolean, toggles river lines on and off.
================  ====================================================================================================================================================================================================================

#############
rcExtraParams
#############
For a-b-c labels and axes title settings, use the ``abc`` and ``title``
categories. For axis tick label settings, use the ``tick`` category.
For figure title, row label, and column label settings, use the new
``suptitle``, ``leftlabel``, ``toplabel``, ``rightlabel``, and ``bottomlabel``
categories.

=======================================================================  ================================================================================================================================================================================================
Key(s)                                                                   Description
=======================================================================  ================================================================================================================================================================================================
``abc.format``                                                           a-b-c label format (for options, see `~proplot.axes.Axes.format_partial`).
``title.pad``                                                            Alias for ``axes.titlepad``, the title offset in arbitrary units
``abc.loc``, ``title.loc``                                               a-b-c label or title position (for options, see `~proplot.axes.Axes.format_partial`).
``abc.border``, ``title.border``                                         Boolean, indicates whether to draw a white border around a-b-c labels or titles located inside an axes.
``abc.linewidth``, ``title.linewidth``                                   Width of the white border around a-b-c labels or titles.
``abc.color``, ``abc.fontsize``, ``abc.weight``                          Font color, size, and weight for a-b-c labels.
``title.color``, ``title.fontsize``, ``title.weight``                    Font color, size, and weight for subplot titles.
``leftlabel.color``, ``leftlabel.fontsize``, ``leftlabel.weight``        Font color, size, and weight for row labels on the left-hand side.
``rightlabel.color``, ``rightlabel.fontsize``, ``rightlabel.weight``     Font color, size, and weight for row labels on the right-hand side.
``toplabel.color``, ``toplabel.fontsize``, ``toplabel.weight``           Font color, size, and weight for column labels on the top of the figure.
``bottomlabel.color``, ``bottomlabel.fontsize``, ``bottomlabel.weight``  Font color, size, and weight for column labels on the bottom of the figure.
``suptitle.color``, ``suptitle.fontsize``, ``suptitle.weight``           Font color, size, and weight for the figure title.
``tick.labelcolor``, ``tick.labelsize``, ``tick.labelweight``            Font color, size, and weight for axis tick labels. These mirror the ``axes.labelcolor``, ``axes.labelsize``, and ``axes.labelweight`` `~matplotlib.rcParams` settings used for axes labels.
``axes.formatter.zerotrim``                                              Boolean, indicates whether trailing decimal zeros are trimmed on tick labels.
``axes.formatter.timerotation``                                          Float, indicates the default *x* axis tick label rotation for datetime tick labels.
=======================================================================  ================================================================================================================================================================================================

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

===================================================================  =============================================================================
Key                                                                  Description
===================================================================  =============================================================================
``geoaxes.facecolor``, ``geoaxes.edgecolor``, ``geoaxes.linewidth``  Face color, edge color, and edge width for the map outline patch.
``land.color``, ``ocean.color``, ``lakes.color``                     Face color for land, ocean, and lake patches.
``rivers.color``, ``rivers.linewidth``                               Line color and linewidth for river lines.
``borders.color``, ``borders.linewidth``                             Line color and linewidth for country border lines.
``innerborders.color``, ``innerborders.linewidth``                   Line color and linewidth for internal border lines.
===================================================================  =============================================================================

There are two new additions to pre-existing ``image`` category, described
in the table below.

=================  =========================================================================================================================================================================================================================================================
Key                Description
=================  =========================================================================================================================================================================================================================================================
``image.levels``   Default number of levels for ``pcolormesh`` and ``contourf`` plots.
``image.edgefix``  Whether to fix the the `white-lines-between-filled-contours <https://stackoverflow.com/q/8263769/4970632>`__ and `white-lines-between-pcolor-rectangles <https://stackoverflow.com/q/27092991/4970632>`__ issues. This slows down figure rendering a bit.
=================  =========================================================================================================================================================================================================================================================

The ``colorbar`` category controls default settings for *inset* and *panel*
colorbars (see `~proplot.axes.Axes.colorbar`). For the last 5 settings,
if float, units are inches. If string, units are interpreted by
`~proplot.utils.units`.

========================  =================================================================================
Key                       Description
========================  =================================================================================
``colorbar.loc``          Inset colorbar location, options are listed in `~proplot.axes.Axes.colorbar`.
``colorbar.grid``         Boolean, indicates whether to draw borders between each level of the colorbar.
``colorbar.frameon``      Boolean, indicates whether to draw a frame behind inset colorbars.
``colorbar.framealpha``   Opacity for inset colorbar frames.
``colorbar.length``       Length of outer colorbars.
``colorbar.lengthinset``  Length of inset colorbars.
``colorbar.width``        Width of outer colorbars.
``colorbar.widthinset``   Width of inset colorbars.
``colorbar.axespad``      Padding between axes edge and inset colorbars.
``colorbar.extend``       Length of rectangular or triangular "extensions" for panel colorbars.
``colorbar.extendinset``  Length of rectangular or triangular "extensions" for inset colorbars.
========================  =================================================================================

The ``subplots`` category controls default layout settings for the
`~proplot.subplots.subplots` command. If float, units are inches.
If string, the units are interpreted by `~proplot.utils.units`.

==========================  ==================================================================
Key                         Description
==========================  ==================================================================
``subplots.axwidth``        Default width of each axes.
``subplots.panelwidth``     Width of side panels.
``subplots.pad``            Padding around figure edge.
``subplots.axpad``          Padding between adjacent subplots.
``subplots.panelpad``       Padding between subplots and panels, and between stacked panels.
``subplots.titlespace``     Vertical space for titles.
``subplots.ylabspace``      Horizontal space between subplots alotted for *y*-labels.
``subplots.xlabspace``      Vertical space between subplots alotted for *x*-labels.
``subplots.innerspace``     Space between subplots alotted for tick marks.
``subplots.panelspace``     Purely empty space between main axes and side panels.
==========================  ==================================================================
"""
# TODO: Add 'style' setting that overrides .proplotrc
# Adapted from seaborn; see: https://github.com/mwaskom/seaborn/blob/master/seaborn/rcmod.py
from . import utils
from .utils import _counter, _timer, DEBUG
import re
import os
import sys
import yaml
import cycler
import warnings
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
if DEBUG:
    import time
    t = time.clock()
import matplotlib.pyplot as plt
if DEBUG:
    print(f'pyplot: {time.clock() - t}')
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
USER_RC = os.path.join(os.path.expanduser("~"), '.proplotrc')
DEFAULT_RC = os.path.join(os.path.dirname(__file__), '.proplotrc')
DEFAULT_FONT = 'Helvetica' if sys.platform == 'linux' else 'Helvetica Neue'
if not os.path.exists(DEFAULT_RC):
    raise ValueError('Default configuration file does not exist.')

# "Global" settings and the lower-level settings they change
# NOTE: This whole section, declaring dictionaries and sets, takes 1ms
RCGLOBALS_CHILDREN = {
    'abc':          [],
    'span':         [],
    'share':        [],
    'align':        [],
    'tight':        [],
    'fontname':     ['font.family'],
    'cmap':         ['image.cmap'],
    'lut':          ['image.lut'],
    'cycle':        [],
    'rgbcycle':     [],
    'alpha':        ['axes.alpha'], # this is a custom setting
    'facecolor':    ['axes.facecolor', 'geoaxes.facecolor'],
    'color':        ['axes.edgecolor', 'geoaxes.edgecolor', 'axes.labelcolor', 'tick.labelcolor', 'hatch.color', 'xtick.color', 'ytick.color'], # change the 'color' of an axes
    'small':        ['font.size', 'tick.labelsize', 'xtick.labelsize', 'ytick.labelsize', 'axes.labelsize', 'legend.fontsize', 'geogrid.labelsize'], # the 'small' fonts
    'large':        ['abc.fontsize', 'figure.titlesize', 'axes.titlesize', 'suptitle.fontsize', 'title.fontsize', 'leftlabel.fontsize', 'toplabel.fontsize', 'rightlabel.fontsize', 'bottomlabel.fontsize'], # the 'large' fonts
    'linewidth':    ['axes.linewidth', 'geoaxes.linewidth', 'hatch.linewidth', 'xtick.major.width', 'ytick.major.width'],
    'margin':       ['axes.xmargin', 'axes.ymargin'],
    'grid':         ['axes.grid'],
    'gridminor':    ['axes.gridminor'],
    'geogrid':      ['axes.geogrid'],
    'ticklen' :     ['xtick.major.size', 'ytick.major.size'],
    'tickdir':      ['xtick.direction',  'ytick.direction'],
    'tickpad':      ['xtick.major.pad', 'xtick.minor.pad', 'ytick.major.pad', 'ytick.minor.pad'],
    'tickratio':    [],
    'ticklenratio': [],
    'gridratio':    [],
    'reso':         [],
    'land':         [],
    'ocean':        [],
    'lakes':        [],
    'coast':        [],
    'borders':      [],
    'innerborders': [],
    'rivers':       [],
    'nbsetup':      [],
    'format':       [],
    'autosave':     [],
    'autoreload':   [],
    }

# Names of the new settings
RC_NAMES = {*_rcParams.keys()}
RC_NAMES_GLOBAL = {*RCGLOBALS_CHILDREN.keys()}
RC_NAMES_CUSTOM = {
    'axes.formatter.zerotrim', 'axes.formatter.timerotation',
    'axes.gridminor', 'axes.geogrid', 'axes.alpha',
    'image.levels', 'image.edgefix',
    'geoaxes.linewidth', 'geoaxes.facecolor', 'geoaxes.edgecolor',
    'land.color', 'ocean.color', 'lakes.color', 'coast.color', 'coast.linewidth',
    'borders.color', 'borders.linewidth', 'innerborders.color', 'innerborders.linewidth', 'rivers.color', 'rivers.linewidth',
    'abc.fontsize', 'abc.weight', 'abc.color', 'abc.loc', 'abc.format', 'abc.border', 'abc.linewidth',
    'title.loc', 'title.pad', 'title.color', 'title.border', 'title.linewidth', 'title.weight', 'title.fontsize',
    'suptitle.fontsize', 'suptitle.weight', 'suptitle.color',
    'leftlabel.fontsize', 'leftlabel.weight', 'leftlabel.color',
    'rightlabel.fontsize', 'rightlabel.weight', 'rightlabel.color',
    'toplabel.fontsize', 'toplabel.weight', 'toplabel.color',
    'bottomlabel.fontsize', 'bottomlabel.weight', 'bottomlabel.color',
    'gridminor.alpha', 'gridminor.color', 'gridminor.linestyle', 'gridminor.linewidth',
    'geogrid.labels', 'geogrid.alpha', 'geogrid.color', 'geogrid.labelsize', 'geogrid.linewidth', 'geogrid.linestyle', 'geogrid.latmax', 'geogrid.lonstep', 'geogrid.latstep',
    'tick.labelweight', 'tick.labelcolor', 'tick.labelsize',
    'subplots.pad', 'subplots.axpad', 'subplots.panelpad',
    'subplots.ylabspace', 'subplots.xlabspace', 'subplots.innerspace', 'subplots.titlespace',
    'subplots.axwidth', 'subplots.panelwidth', 'subplots.panelspace',
    'colorbar.grid', 'colorbar.frameon', 'colorbar.framealpha',
    'colorbar.loc', 'colorbar.length', 'colorbar.width', 'colorbar.lengthinset', 'colorbar.widthinset',
    'colorbar.extend', 'colorbar.extendinset', 'colorbar.axespad',
    }
# Used by Axes.format, allows user to pass rc settings as keyword args,
# way less verbose. For example, landcolor='b' vs. rc_kw={'land.color':'b'}.
RC_NAMES_NODOTS = { # useful for passing these as kwargs
    name.replace('.', ''):name for names in
    (RC_NAMES_CUSTOM, RC_NAMES, RC_NAMES_GLOBAL)
    for name in names
    }
# Categories for returning dict of subcategory properties
RC_CATEGORIES = {
    *(re.sub('\.[^.]*$', '', name) for names in (RC_NAMES_CUSTOM, RC_NAMES) for name in names),
    *(re.sub('\..*$', '', name) for names in (RC_NAMES_CUSTOM, RC_NAMES) for name in names)
    }

# Unit conversion
# See: https://matplotlib.org/users/customizing.html, all props matching
# the below strings use the units 'points', and my special categories are inches!
UNITS_REGEX = re.compile('^.*(width|space|size|pad|len|small|large)$') # len is for ticklen, space is for subplots props
def _convert_units(key, value):
    """Converts certain keys to the units "points". If "key" is passed, tests
    that key against possible keys that accept physical units."""
    # WARNING: Must keep colorbar and subplots units alive, so when user
    # requests em units, values change with respect to font size. The points
    # thing is a conveniene feature so not as important for them.
    if isinstance(value,str) and UNITS_REGEX.match(key) \
        and key.split('.')[0] not in ('colorbar','subplots'):
            value = utils.units(value, 'pt')
    return value

def _set_cycler(name):
    """Sets the default color cycler."""
    # Draw from dictionary
    try:
        colors = mcm.cmap_d[name].colors
    except Exception:
        cycles = sorted(name for name,cmap in mcm.cmap_d.items() if isinstance(cmap, mcolors.ListedColormap))
        raise ValueError(f'Invalid cycle name {name!r}. Options are: {", ".join(cycles)}')
    # Apply color name definitions
    if _rcGlobals['rgbcycle'] and name.lower() == 'colorblind':
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
            if key == 'ticklen':
                ticklen = _convert_units(key, value)
                ratio = _rcGlobals['ticklenratio']
            else:
                ticklen = _rcGlobals['ticklen']
                ratio = value
            kw['xtick.minor.size'] = ticklen*ratio
            kw['ytick.minor.size'] = ticklen*ratio
        # Spine width/major-minor tick width ratio
        elif key in ('linewidth', 'tickratio'):
            if key == 'linewidth':
                tickwidth = _convert_units(key, value)
                ratio = _rcGlobals['tickratio']
            else:
                tickwidth = _rcGlobals['linewidth']
                ratio = value
            kw['xtick.minor.width'] = tickwidth*ratio
            kw['ytick.minor.width'] = tickwidth*ratio
        # Grid line
        elif key in ('grid.linewidth', 'gridratio'):
            if key == 'grid.linewidth':
                gridwidth = _convert_units(key, value)
                ratio = _rcGlobals['gridratio']
            else:
                gridwidth = _rcParams['grid.linewidth']
                ratio = value
            kw_custom['gridminor.linewidth'] = gridwidth*ratio
        # Now update linked settings
        if key not in ('gridratio','tickratio','ticklenratio'):
            for name in RCGLOBALS_CHILDREN[key]:
                val = _convert_units(key, value)
                if name in _rcExtraParams:
                    kw_custom[name] = val
                else:
                    kw[name] = val
            if key == 'linewidth' and val == 0:
                ikw, ikw_custom = _get_globals('ticklen', 0)
                kw.update(ikw)
                kw_custom.update(ikw_custom)
    return kw, kw_custom

class rc_configurator(object):
    _public_api = (
        'get', 'fill', 'category', 'reset', 'context', 'update'
        ) # getattr and setattr will not look for these items
    def __str__(self):
        return type(_rcParams).__str__(_rcGlobals) # just show globals
    def __repr__(self):
        return type(_rcParams).__repr__(_rcGlobals)

    @_counter # about 0.05s
    def __init__(self):
        """Magical abstract class for managing matplotlib `rcParams
        <https://matplotlib.org/users/customizing.html>`__ settings, ProPlot
        :ref:`rcExtraParams` settings, and :ref:`rcGlobals` "global" settings.
        When initialized, this loads defaults settings plus any user overrides
        in the ``~/.proplotrc`` file. See the `~proplot.rctools` documentation
        for details."""
        # Set the default style. Note that after first figure made, backend
        # is 'sticky', never changes! See: https://stackoverflow.com/a/48322150/4970632
        plt.style.use('default')

        # Load the defaults from file
        for i,file in enumerate((DEFAULT_RC, USER_RC)):
            # Load
            if not os.path.exists(file):
                continue
            with open(file) as f:
                try:
                    data = yaml.safe_load(f)
                except yaml.YAMLError as err:
                    print('Error: Invalid .proplotrc file.')
                    raise err
            # Special duplicate keys
            val = data.get('title.pad', None)
            if val is not None:
                data['axes.titlepad'] = val
            # Add keys to dictionaries
            gkeys, ckeys = {*()}, {*()}
            for key,value in data.items():
                if key in RC_NAMES_GLOBAL:
                    _rcGlobals[key] = value
                    if i == 0:
                        gkeys.add(key)
                elif key in RC_NAMES_CUSTOM:
                    value = _convert_units(key, value)
                    _rcExtraParams[key] = value
                    if i == 0:
                        ckeys.add(key)
                elif key in RC_NAMES:
                    value = _convert_units(key, value)
                    _rcParams[key] = value
                else:
                    raise RuntimeError(('Default', 'User')[i] + f' .proplotrc file has invalid key {key!r}.')
            # Make sure we did not miss anything
            if i == 0:
                if gkeys != RC_NAMES_GLOBAL:
                    raise RuntimeError(f'Default .proplotrc file has incomplete or invalid global keys {RC_NAMES_GLOBAL - gkeys}.')
                if ckeys != RC_NAMES_CUSTOM:
                    raise RuntimeError(f'Default .proplotrc file has incomplete or invalid custom keys {RC_NAMES_CUSTOM - ckeys}.')

        # Set default fontname and cycler
        _set_cycler(_rcGlobals['cycle'])
        if _rcGlobals.get('fontname', None) is None:
            _rcGlobals['fontname'] = DEFAULT_FONT

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
        is not cached (i.e. if it was not changed by the user)."""
        # Can get a whole bunch of different things
        # Get full dictionary e.g. for rc[None]
        if not key:
            return {**_rcParams, **_rcExtraParams}

        # Standardize
        # NOTE: If key is invalid, raise error down the line.
        if '.' not in key and key not in _rcGlobals:
            key = RC_NAMES_NODOTS.get(key, key)

        # Allow for special time-saving modes where we *ignore _rcParams*
        # or even *ignore _rcExtraParams*.
        mode = self._getitem_mode
        if mode == 0:
            kws = (self._cache, _rcGlobals, _rcExtraParams, _rcParams)
        elif mode == 1:
            kws = (self._cache, _rcGlobals, _rcExtraParams) # custom only!
        elif mode == 2:
            kws = (self._cache,) # changed only!
        else:
            raise KeyError(f'Invalid caching mode {mode!r}.')

        # Get individual property. Will successively index a few different dicts
        # Try to return the value
        for kw in kws:
            try:
                return kw[key]
            except KeyError:
                continue
        # If we were in one of the exclusive modes, return None
        if mode == 0:
            raise KeyError(f'Invalid prop name {key!r}.')
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
            key = RC_NAMES_NODOTS.get(key, key)

        # Special keys
        if key == 'title.pad':
            key = 'axes.titlepad'
        if key == 'rgbcycle': # if must re-apply cycler afterward
            cache[key] = value
            _rcGlobals[key] = value
            key, value = 'cycle', _rcGlobals['cycle']

        # Set the default cycler
        if key == 'cycle':
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
                if not ovalue or (key == 'grid' and owhich == 'major') or (key == 'gridminor' and owhich == 'minor'):
                    which = 'both' # disable both sides
                # Gridlines are currently on for major and minor ticks, so we instruct
                # to turn on gridlines for the one we *don't* want off
                elif owhich == 'both': # and ovalue is True, as we already tested
                    value = True
                    which = 'major' if key == 'gridminor' else 'minor' # if gridminor=False, enable major, and vice versa
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
                if owhich == 'both' or (key == 'grid' and owhich == 'minor') or (key == 'gridminor' and owhich == 'major'):
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
        # Update normal settings
        elif key in RC_NAMES_CUSTOM:
            value = _convert_units(key, value)
            cache[key] = value
            if context:
                restore[key] = _rcExtraParams[key]
            _rcExtraParams[key] = value
        elif key in RC_NAMES:
            value = _convert_units(key, value)
            cache[key] = value
            if context:
                restore[key] = _rcParams[key]
            _rcParams[key] = value # rcParams dict has key validation
        else:
            raise KeyError(f'Invalid key {key!r}.')
        self._init = False # setitem was successful, we are no longer in initial state

    # Attributes same as items
    def __getattribute__(self, attr):
        """Invokes `~rc_configurator.__getitem__`."""
        if attr[:1] == '_' or attr in self._public_api:
            return object.__getattribute__(self, attr)
        else:
            return self[attr]

    def __setattr__(self, attr, value):
        """Invokes `~rc_configurator.__setitem__`."""
        if attr[:1] == '_':
            object.__setattr__(self, attr, value)
        else:
            self[attr] = value

    # Immutability
    def __delitem__(self, *args):
        """Pseudo-immutability."""
        raise RuntimeError('rc settings cannot be deleted.')

    def __delattr__(self, *args):
        """Pseudo-immutability."""
        raise RuntimeError('rc settings cannot be deleted.')

    # Context tools
    def __enter__(self):
        """Apply settings from configurator cache."""
        self._cache_orig = rc._cache.copy()
        for key,value in self._context.items():
            self[key] = value # applies globally linked and individual settings

    def __exit__(self, *args):
        """Restore configurator cache to initial state."""
        self._context = {}
        self._getitem_mode = 0
        for key,value in self._cache_restore.items():
            self[key] = value
        self._cache = self._cache_orig
        self._cache_restore = {}
        self._cache_orig = {}

    def context(self, *args, mode=0, **kwargs):
        """
        Temporarily modifies settings in a ``with...as`` block,
        used by ProPlot internally but may also be useful for power users.

        This function was invented to prevent successive calls to
        `~proplot.axes.Axes.format` from constantly looking up and
        re-applying unchanged settings. Testing showed that these gratuitous
        `rcParams <https://matplotlib.org/users/customizing.html>`__
        lookups and artist updates increased runtime by seconds, even for
        relatively simple plots.

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
               This is used in the `~proplot.axes.Axes.__init__` call to
               `~proplot.axes.Axes.format`. When a setting lookup returns
               ``None``, `~proplot.axes.Axes.format` does not apply it.
            2. All unchanged settings return ``None``. This is used during user
               calls to `~proplot.axes.Axes.format`.

        Example
        -------

        >>> import proplot as plot
        >>> with plot.rc.context(linewidth=2, ticklen=5):
        ...     f, ax = plot.subplots()
        ...     ax.plot(data)

        """
        if mode not in range(3):
            raise ValueError(f'Invalid _getitem_mode {mode}.')
        for arg in args:
            if not isinstance(arg, dict):
                raise ValueError('Non-dictionary argument.')
            kwargs.update(arg)
        self._context = kwargs # could be empty
        self._getitem_mode = mode
        return self

    # Other tools
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
        if cat not in RC_CATEGORIES:
            raise ValueError(f'RC category {cat} does not exist. Valid categories are {", ".join(RC_CATEGORIES)}.')
        if not cache:
            mode = 0
        else:
            mode = self._getitem_mode

        # Allow for special time-saving modes where we *ignore _rcParams*
        # or even *ignore _rcExtraParams*.
        if mode == 0:
            kws = (self._cache, _rcGlobals, _rcExtraParams, _rcParams)
        elif mode == 1:
            kws = (self._cache, _rcGlobals, _rcExtraParams)
        elif mode == 2:
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
        if len(args) > 2:
            raise ValueError('Accepts 1-2 positional arguments. Use plot.rc.update(kw) to update a bunch of names, or plot.rc.update(category, kw) to update subcategories belonging to single category e.g. axes. All kwargs will be added to the dict.')
        elif len(args) == 2:
            prefix = args[0]
            kw = args[1]
        elif len(args) == 1:
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
        """Restores settings to the initial state -- ProPlot defaults, plus
        any user overrides in the ``~/.proplotrc`` file."""
        if not self._init: # save resources if rc is unchanged!
            return self.__init__()

# Declare rc object
# WARNING: Must be instantiated after ipython notebook setup! The default
# backend may change some rc settings!
rc = rc_configurator()
"""Instance of `rc_configurator`. This is used to change global settings.
See the `~proplot.rctools` documentation for details."""

# Ipython notebook behavior
@_timer
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
        ipython.magic("matplotlib inline")
        rc.reset()
    except KeyError:
        # For terminals
        # TODO: Fix auto tight layout with osx backend -- currently has
        # standard issue where content adjusted but canvas size not adjusted
        # until second draw command, and other issue where window size does
        # not sync with figure size
        ipython.magic("matplotlib qt")
        rc.reset()
    else:
        # Choose svg vector or retina hi-res bitmap backends
        autosave = _rcGlobals['autosave']
        if autosave: # capture annoying message + line breaks
            with IPython.utils.io.capture_output():
                ipython.magic("autosave " + str(autosave))
        ipython.magic("config InlineBackend.figure_formats = ['" + _rcGlobals['format'] + "']")
        ipython.magic("config InlineBackend.rc = {}") # no notebook-specific overrides
        ipython.magic("config InlineBackend.close_figures = True") # memory issues
        ipython.magic("config InlineBackend.print_figure_kwargs = {'bbox_inches':None}") # use ProPlot tight layout

# Setup notebook and issue warning
# TODO: Add to list of incompatible backends?
if _rcGlobals['nbsetup']:
    nb_setup()
if _rcParams['backend'][:2] == 'nb' or _rcParams['backend'] in ('MacOSX',):
    warnings.warn(f'Due to automatic figure resizing, using ProPlot with the {_rcParams["backend"]!r} backend may result in unexpected behavior. Try using %matplotlib inline or %matplotlib qt, or just import ProPlot before specifying the backend and ProPlot will automatically load it.')

