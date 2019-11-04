#!/usr/bin/env python3
"""
A special object named `~proplot.rctools.rc`, belonging to the
`~proplot.rctools.rc_configurator` class, is created on import.
This is your one-stop shop for changing global settings belonging to any of
the following three categories.

1. Builtin matplotlib `rcParams <https://matplotlib.org/users/customizing.html>`__
   settings. These have the format ``x.y`` or ``x.y.z``.
2. ProPlot :ref:`rcParamsLong` settings. These also have the format ``x.y``
   (see below).
3. ProPlot :ref:`rcParamsShort` settings. These have no dots (see below).

You can change settings with the `~proplot.rctools.rc` object as follows.

* ``plot.rc.name = value``
* ``plot.rc['name'] = value``
* ``plot.rc.update(name1=value1, name2=value2)``
* ``plot.rc.update({'name1':value1, 'name2':value2})``

To temporarily change settings on a particular axes, use either of the
following.

* ``ax.format(name=value)``
* ``ax.format(rc_kw={'name':value})``

In all of these examples, if the setting name ``name`` contains
any dots, you can simply omit the dots. For example, to change the
:rcraw:`title.loc` property, use ``plot.rc.titleloc = value``,
``plot.rc.update(titleloc=value)``, or ``ax.format(titleloc=value)``.

#############
rcParamsShort
#############

These are **simple, short** names used to change multiple matplotlib and
ProPlot settings at once, as shorthands for settings with longer names, or
for special options. For example, :rcraw:`ticklen` changes the tick length for
the *x* and *y* axes in one go.

================  ====================================================================================================================================================================================================================================
Key               Description
================  ====================================================================================================================================================================================================================================
``nbsetup``       Whether to run `nb_setup` on import. Can only be changed from the ``~/.proplotrc`` file.
``format``        The inline backend figure format, one of ``retina``, ``png``, ``jpeg``, ``pdf``, or ``svg``. Can only be changed from the ``~/.proplotrc`` file.
``autosave``      If not empty or ``0`` and :rcraw:`nbsetup` is ``True``, passed to `%autosave <https://www.webucator.com/blog/2016/03/change-default-autosave-interval-in-ipython-notebook/>`__. Can only be changed from the ``~/.proplotrc`` file.
``autoreload``    If not empty or ``0`` and :rcraw:`nbsetup` is ``True``, passed to `%autoreload <https://ipython.readthedocs.io/en/stable/config/extensions/autoreload.html#magic-autoreload>`__. Can only be changed from the ``~/.proplotrc`` file.
``abc``           Boolean, indicates whether to draw a-b-c labels by default.
``tight``         Boolean, indicates whether to auto-adjust figure bounds and subplot spacings.
``share``         The axis sharing level, one of ``0``, ``1``, ``2``, or ``3``. See `~proplot.subplots.subplots` for details.
``align``         Whether to align axis labels during draw. See `aligning labels <https://matplotlib.org/3.1.1/gallery/subplots_axes_and_figures/align_labels_demo.html>`__.
``span``          Boolean, toggles spanning axis labels. See `~proplot.subplots.subplots` for details.
``fontname``      Name of font used for all text in the figure. The default is Helvetica Neue. See `~proplot.fonttools` for details.
``cmap``          The default colormap.
``lut``           The number of colors to put in the colormap lookup table.
``cycle``         The default color cycle name, used e.g. for lines.
``rgbcycle``      If ``True``, and ``colorblind`` is the current cycle, this registers the ``colorblind`` colors as ``'r'``, ``'b'``, ``'g'``, etc., like in `seaborn <https://seaborn.pydata.org/tutorial/color_palettes.html>`__.
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
================  ====================================================================================================================================================================================================================================

############
rcParamsLong
############
The ``subplots`` category controls the default layout for figures
and axes. The ``abc``, ``title``, and ``tick`` categories control
a-b-c label, title, and axis tick label settings. The
``suptitle``, ``leftlabel``, ``toplabel``, ``rightlabel``, and ``bottomlabel``
categories control figure title and edge label settings.

There are two new additions to the ``image`` category, and the new
``colorbar`` category controls *inset* and *outer*
`~proplot.axes.Axes.colorbar` properties.

The new ``gridminor`` category controls minor gridline settings,
and the new ``geogrid`` category controls meridian and parallel line settings
for `~proplot.axes.ProjectionAxes`. For both ``gridminor`` and ``geogrid``, if
a property is empty, the corresponding property from ``grid`` is used.

Finally, the ``geoaxes``, ``land``, ``ocean``, ``rivers``, ``lakes``,
``borders``, and ``innerborders`` categories control various
`~proplot.axes.ProjectionAxes` settings. These are used when the boolean
toggles for the corresponding :ref:`rcParamsShort` settings are turned on.

===================================================================  =========================================================================================================================================================================================================================================================
Key(s)                                                               Description
===================================================================  =========================================================================================================================================================================================================================================================
``abc.style``                                                        a-b-c label style. For options, see `~proplot.axes.Axes.format`.
``abc.loc``                                                          a-b-c label position. For options, see `~proplot.axes.Axes.format`.
``abc.border``                                                       Boolean, indicates whether to draw a white border around a-b-c labels inside an axes.
``abc.linewidth``                                                    Width of the white border around a-b-c labels.
``abc.color``, ``abc.size``, ``abc.weight``                          Font color, size, and weight for a-b-c labels.
``axes.formatter.zerotrim``                                          Boolean, indicates whether trailing decimal zeros are trimmed on tick labels.
``axes.formatter.timerotation``                                      Float, indicates the default *x* axis tick label rotation for datetime tick labels.
``borders.color``, ``borders.linewidth``                             Line color and linewidth for country border lines.
``bottomlabel.color``, ``bottomlabel.size``, ``bottomlabel.weight``  Font color, size, and weight for column labels on the bottom of the figure.
``colorbar.loc``                                                     Inset colorbar location, options are listed in `~proplot.axes.Axes.colorbar`.
``colorbar.grid``                                                    Boolean, indicates whether to draw borders between each level of the colorbar.
``colorbar.frameon``                                                 Boolean, indicates whether to draw a frame behind inset colorbars.
``colorbar.framealpha``                                              Opacity for inset colorbar frames.
``colorbar.length``                                                  Length of outer colorbars.
``colorbar.insetlength``                                             Length of inset colorbars. Units are interpreted by `~proplot.utils.units`.
``colorbar.width``                                                   Width of outer colorbars. Units are interpreted by `~proplot.utils.units`.
``colorbar.insetwidth``                                              Width of inset colorbars. Units are interpreted by `~proplot.utils.units`.
``colorbar.insetpad``                                                Padding between axes edge and inset colorbars. Units are interpreted by `~proplot.utils.units`.
``colorbar.extend``                                                  Length of rectangular or triangular "extensions" for panel colorbars. Units are interpreted by `~proplot.utils.units`.
``colorbar.insetextend``                                             Length of rectangular or triangular "extensions" for inset colorbars. Units are interpreted by `~proplot.utils.units`.
``geoaxes.facecolor``, ``geoaxes.edgecolor``, ``geoaxes.linewidth``  Face color, edge color, and edge width for the map outline patch.
``geogrid.labels``                                                   Boolean, indicates whether to label the parallels and meridians.
``geogrid.labelsize``                                                Font size for latitide and longitude labels. Inherits from ``small``.
``geogrid.latmax``                                                   Absolute latitude in degrees, poleward of which meridian gridlines are cut off.
``geogrid.lonstep``, ``geogrid.latstep``                             Interval for meridian and parallel gridlines, in degrees.
``gridminor.linewidth``, ``geogrid.linewidth``                       The line width.
``gridminor.linestyle``, ``geogrid.linestyle``                       The line style.
``gridminor.alpha``, ``geogrid.alpha``                               The line transparency.
``gridminor.color``, ``geogrid.color``                               The line color.
``image.levels``                                                     Default number of levels for ``pcolormesh`` and ``contourf`` plots.
``image.edgefix``                                                    Whether to fix the the `white-lines-between-filled-contours <https://stackoverflow.com/q/8263769/4970632>`__ and `white-lines-between-pcolor-rectangles <https://stackoverflow.com/q/27092991/4970632>`__ issues. This slows down figure rendering a bit.
``innerborders.color``, ``innerborders.linewidth``                   Line color and linewidth for internal border lines.
``land.color``, ``ocean.color``, ``lakes.color``                     Face color for land, ocean, and lake patches.
``leftlabel.color``, ``leftlabel.size``, ``leftlabel.weight``        Font color, size, and weight for row labels on the left-hand side.
``rightlabel.color``, ``rightlabel.size``, ``rightlabel.weight``     Font color, size, and weight for row labels on the right-hand side.
``rivers.color``, ``rivers.linewidth``                               Line color and linewidth for river lines.
``subplots.axwidth``                                                 Default width of each axes. Units are interpreted by `~proplot.utils.units`.
``subplots.panelwidth``                                              Width of side panels. Units are interpreted by `~proplot.utils.units`.
``subplots.pad``                                                     Padding around figure edge. Units are interpreted by `~proplot.utils.units`.
``subplots.axpad``                                                   Padding between adjacent subplots. Units are interpreted by `~proplot.utils.units`.
``subplots.panelpad``                                                Padding between subplots and panels, and between stacked panels. Units are interpreted by `~proplot.utils.units`.
``subplots.titlespace``                                              Vertical space for titles. Units are interpreted by `~proplot.utils.units`.
``subplots.ylabspace``                                               Horizontal space between subplots alotted for *y*-labels. Units are interpreted by `~proplot.utils.units`.
``subplots.xlabspace``                                               Vertical space between subplots alotted for *x*-labels. Units are interpreted by `~proplot.utils.units`.
``subplots.innerspace``                                              Space between subplots alotted for tick marks. Units are interpreted by `~proplot.utils.units`.
``subplots.panelspace``                                              Purely empty space between main axes and side panels. Units are interpreted by `~proplot.utils.units`.
``suptitle.color``, ``suptitle.size``, ``suptitle.weight``           Font color, size, and weight for the figure title.
``tick.labelcolor``, ``tick.labelsize``, ``tick.labelweight``        Font color, size, and weight for axis tick labels. These mirror the ``axes.labelcolor``, ``axes.labelsize``, and ``axes.labelweight`` `~matplotlib.rcParams` settings used for axes labels.
``title.loc``                                                        Title position. For options, see `~proplot.axes.Axes.format`.
``title.border``                                                     Boolean, indicates whether to draw a white border around titles inside an axes.
``title.linewidth``                                                  Width of the white border around titles.
``title.pad``                                                        Alias for ``axes.titlepad``, the title offset in arbitrary units
``title.color``, ``title.size``, ``title.weight``                    Font color, size, and weight for subplot titles.
``toplabel.color``, ``toplabel.size``, ``toplabel.weight``           Font color, size, and weight for column labels on the top of the figure.
===================================================================  =========================================================================================================================================================================================================================================================

##############
proplotrc file
##############

To modify the global settings, edit your
``~/.proplotrc`` file. To modify settings for a particular project,
create a ``.proplotrc`` file in the same directory as your ipython
notebook, or in an arbitrary parent directory.

As an example, the default ``.proplotrc`` file
is shown below. The syntax is roughly the same as that used for
``matplotlibrc`` files, although ``.proplotrc`` strictly adheres to
`YAML <https://en.wikipedia.org/wiki/YAML>`__.

.. include:: ../proplot/.proplotrc
   :literal:

"""
# TODO: Add 'style' setting that overrides .proplotrc
# Adapted from seaborn; see: https://github.com/mwaskom/seaborn/blob/master/seaborn/rcmod.py
from . import utils
from .utils import _counter, _timer, _benchmark
import re
import os
import yaml
import cycler
import warnings
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
with _benchmark('pyplot'):
    import matplotlib.pyplot as plt
try:
    import IPython
    get_ipython = IPython.get_ipython
except ModuleNotFoundError:
    get_ipython = lambda: None
__all__ = ['rc', 'rc_configurator', 'nb_setup']

# Initialize
from matplotlib import rcParams, RcParams
class RcParamsShort(dict):
    """Class for holding short-name settings. So far just `dict`."""
    pass
class RcParamsLong(dict):
    """Class for holding custom settings. So far just `dict`."""
    pass
defaultParamsShort = RcParamsShort({
    'nbsetup':      True,
    'format':       'retina',
    'autosave':     30,
    'autoreload':   2,
    'abc':          False,
    'share':        3,
    'align':        False,
    'span':         True,
    'tight':        True,
    'fontname':     'Helvetica Neue',
    'cmap':         'fire',
    'lut':          256,
    'cycle':        'colorblind',
    'rgbcycle':     False,
    'color':        'k',
    'alpha':        1,
    'facecolor':    'w',
    'small':        8,
    'large':        9,
    'linewidth':    0.6,
    'margin':       0.0,
    'grid':         True,
    'gridminor':    False,
    'ticklen':      4.0,
    'tickdir':      'out',
    'tickpad':      2.0,
    'tickratio':    0.8,
    'ticklenratio': 0.5,
    'gridratio':    0.5,
    'reso':         'lo',
    'geogrid':      True,
    'land':         False,
    'ocean':        False,
    'coast':        False,
    'rivers':       False,
    'lakes':        False,
    'borders':      False,
    'innerborders': False,
    })
defaultParamsLong = RcParamsLong({
    'abc.loc':                     'l', # left side above the axes
    'title.loc':                   'c', # centered above the axes
    'title.pad':                   3.0, # copy
    'abc.style':                   'a',
    'abc.size':                    None, # = large
    'abc.color':                   'k',
    'abc.weight':                  'bold',
    'abc.border':                  True,
    'abc.linewidth':               1.5,
    'tick.labelsize':              None, # = small
    'tick.labelcolor':             None, # = color
    'tick.labelweight':            'normal',
    'title.size':                  None, # = large
    'title.color':                 'k',
    'title.weight':                'normal',
    'title.border':                True,
    'title.linewidth':             1.5,
    'suptitle.size':               None, # = large
    'suptitle.color':              'k',
    'suptitle.weight':             'bold',
    'leftlabel.size':              None, # = large
    'leftlabel.weight':            'bold',
    'leftlabel.color':             'k',
    'toplabel.size':               None, # = large
    'toplabel.weight':             'bold',
    'toplabel.color':              'k',
    'rightlabel.size':             None, # = large
    'rightlabel.weight':           'bold',
    'rightlabel.color':            'k',
    'bottomlabel.size':            None, # = large
    'bottomlabel.weight':          'bold',
    'bottomlabel.color':           'k',
    'image.edgefix':               True,
    'image.levels':                11,
    'axes.alpha':                  None, # if empty, depends on 'savefig.transparent' setting
    'axes.formatter.zerotrim':     True,
    'axes.formatter.timerotation': 90,
    'axes.gridminor':              True,
    'axes.geogrid':                True,
    'gridminor.alpha':             None, # = grid.alpha
    'gridminor.color':             None, # = grid.color
    'gridminor.linestyle':         None, # = grid.linewidth
    'gridminor.linewidth':         None, # = grid.linewidth x gridratio
    'geogrid.labels':              False,
    'geogrid.labelsize':           None, # = small
    'geogrid.latmax':              90,
    'geogrid.lonstep':             30,
    'geogrid.latstep':             20,
    'geogrid.alpha':               0.5,
    'geogrid.color':               'k',
    'geogrid.linewidth':           1.0,
    'geogrid.linestyle':           ':                                                         ',
    'geoaxes.linewidth':           None, # = linewidth
    'geoaxes.facecolor':           None, # = facecolor
    'geoaxes.edgecolor':           None, # = color
    'land.color':                  'k',
    'ocean.color':                 'w',
    'lakes.color':                 'w',
    'coast.color':                 'k',
    'coast.linewidth':             0.6,
    'borders.color':               'k',
    'borders.linewidth':           0.6,
    'innerborders.color':          'k',
    'innerborders.linewidth':      0.6,
    'rivers.color':                'k',
    'rivers.linewidth':            0.6,
    'colorbar.loc':                'right',
    'colorbar.grid':               False,
    'colorbar.frameon':            True,
    'colorbar.framealpha':         0.8,
    'colorbar.insetpad':           '0.5em',
    'colorbar.extend':             '1.3em',
    'colorbar.insetextend':        '1em',
    'colorbar.length':             1,
    'colorbar.insetlength':        '8em',
    'colorbar.width':              '1.5em',
    'colorbar.insetwidth':         '1.2em',
    'subplots.axwidth':            '18em',
    'subplots.panelwidth':         '4em',
    'subplots.pad':                '0.5em',
    'subplots.axpad':              '1em',
    'subplots.panelpad':           '0.5em',
    'subplots.panelspace':         '1em',
    'subplots.innerspace':         '1.5em',
    'subplots.ylabspace':          '5.5em',
    'subplots.xlabspace':          '4em',
    'subplots.titlespace':         '2em',
    })
defaultParams = RcParams({
    'figure.dpi':              90,
    'figure.facecolor':        '#f2f2f2',
    'figure.autolayout':       False,
    'figure.titleweight':      'bold',
    'figure.max_open_warning': 0,
    'savefig.directory':       '',
    'savefig.dpi':             300,
    'savefig.facecolor':       'white',
    'savefig.transparent':     True,
    'savefig.format':          'pdf',
    'savefig.bbox':            'standard',
    'savefig.pad_inches':      0.0,
    'axes.titleweight':        'normal',
    'axes.xmargin':            0.0,
    'axes.ymargin':            0.0,
    'axes.grid':               True,
    'axes.labelpad':           3.0,
    'axes.titlepad':           3.0,
    'grid.color':              'k',
    'grid.alpha':              0.1,
    'grid.linewidth':          0.6,
    'grid.linestyle':          '-',
    'hatch.color':             'k',
    'hatch.linewidth':         0.6,
    'lines.linewidth':         1.3,
    'lines.markersize':        3.0,
    'legend.frameon':          True,
    'legend.framealpha':       0.8,
    'legend.fancybox':         False,
    'legend.labelspacing':     0.5,
    'legend.handletextpad':    0.5,
    'legend.handlelength':     1.5,
    'legend.columnspacing':    1.0,
    'legend.borderpad':        0.5,
    'legend.borderaxespad':    0,
    'xtick.minor.visible':     True,
    'ytick.minor.visible':     True,
    'mathtext.bf':             'sans:bold',
    'mathtext.it':             'sans:it',
    'mathtext.default':        'regular',
    })
rcParamsShort = RcParamsShort({})
rcParamsLong = RcParamsLong({})

# Initialize user file
_rc_file = os.path.join(os.path.expanduser('~'), '.proplotrc')
if not os.path.isfile(_rc_file):
    def _tabulate(rcdict):
        string = ''
        maxlen = max(map(len, rcdict))
        for key,value in rcdict.items():
            value = '' if value is None else repr(value)
            space = ' ' * (maxlen - len(key) + 1) * int(bool(value))
            string += f'#  {key}:{space}{value}\n'
    with open(_rc_file, 'x') as f:
        f.write(f"""
#------------------------------------------------------
# Use this file to customize your settings
# For descriptions of each key name see:
# https://proplot.readthedocs.io/en/latest/rctools.html
#------------------------------------------------------
# ProPlot short name settings
{_tabulate(defaultParamsShort)}
#
# ProPlot long name settings
{_tabulate(defaultParamsLong)}
#
# Matplotlib settings
{_tabulate(defaultParams)}
""")

# "Global" settings and the lower-level settings they change
# NOTE: This whole section, declaring dictionaries and sets, takes 1ms
RC_CHILDREN = {
    'fontname':     ('font.family',),
    'cmap':         ('image.cmap',),
    'lut':          ('image.lut',),
    'alpha':        ('axes.alpha',), # this is a custom setting
    'facecolor':    ('axes.facecolor', 'geoaxes.facecolor'),
    'color':        ('axes.edgecolor', 'geoaxes.edgecolor', 'axes.labelcolor', 'tick.labelcolor', 'hatch.color', 'xtick.color', 'ytick.color'), # change the 'color' of an axes
    'small':        ('font.size', 'tick.labelsize', 'xtick.labelsize', 'ytick.labelsize', 'axes.labelsize', 'legend.fontsize', 'geogrid.labelsize'), # the 'small' fonts
    'large':        ('abc.size', 'figure.titlesize', 'axes.titlesize', 'suptitle.size', 'title.size', 'leftlabel.size', 'toplabel.size', 'rightlabel.size', 'bottomlabel.size'), # the 'large' fonts
    'linewidth':    ('axes.linewidth', 'geoaxes.linewidth', 'hatch.linewidth', 'xtick.major.width', 'ytick.major.width'),
    'margin':       ('axes.xmargin', 'axes.ymargin'),
    'grid':         ('axes.grid',),
    'gridminor':    ('axes.gridminor',),
    'geogrid':      ('axes.geogrid',),
    'ticklen':      ('xtick.major.size', 'ytick.major.size'),
    'tickdir':      ('xtick.direction',  'ytick.direction'),
    'tickpad':      ('xtick.major.pad', 'xtick.minor.pad', 'ytick.major.pad', 'ytick.minor.pad'),
    }

# Used by Axes.format, allows user to pass rc settings as keyword args,
# way less verbose. For example, landcolor='b' vs. rc_kw={'land.color':'b'}.
RC_NODOTS = { # useful for passing these as kwargs
    name.replace('.', ''): name for names
    in (rcParams, rcParamsLong) for name in names
    }
# Categories for returning dict of subcategory properties
RC_CATEGORIES = {
    *(re.sub('\.[^.]*$', '', name) for names in (rcParams, rcParamsLong) for name in names),
    *(re.sub('\..*$', '', name) for names in (rcParams, rcParamsLong) for name in names)
    }

# Helper funcs
def _to_points(key, value):
    """Converts certain keys to the units "points". If "key" is passed, tests
    that key against possible keys that accept physical units."""
    # See: https://matplotlib.org/users/customizing.html, all props matching
    # the below strings use the units 'points', and custom categories are in
    if (isinstance(value,str) and key.split('.')[0] not in ('colorbar','subplots')
        and re.match('^.*(width|space|size|pad|len|small|large)$', key)):
        value = utils.units(value, 'pt')
    return value

def _get_config_paths():
    """Returns configuration file paths."""
    # Local configuration
    idir = os.getcwd()
    paths = []
    while idir: # not empty string
        ipath = os.path.join(idir, '.proplotrc')
        if os.path.exists(ipath):
            paths.append(ipath)
        ndir, _ = os.path.split(idir)
        if ndir == idir:
            break
        idir = ndir
    paths = paths[::-1] # sort from decreasing to increasing importantce
    # Home configuration
    ipath = os.path.join(os.path.expanduser('~'), '.proplotrc')
    if os.path.exists(ipath) and ipath not in paths:
        paths.insert(0, ipath)
    return paths

def _get_synced_params(key=None, value=None):
    """Returns dictionaries for updating "child" properties in
    `rcParams` and `rcParamsLong` with global property."""
    kw = {} # builtin properties that global setting applies to
    kw_custom = {} # custom properties that global setting applies to
    if key is not None and value is not None:
        items = [(key,value)]
    else:
        items = rcParamsShort.items()
    for key,value in items:
        # Cycler
        if key in ('cycle', 'rgbcycle'):
            if key == 'rgbcycle':
                cycle, rgbcycle = rcParamsShort['cycle'], value
            else:
                cycle, rgbcycle = value, rcParamsShort['rgbcycle']
            try:
                colors = mcm.cmap_d[cycle].colors
            except (KeyError, AttributeError):
                cycles = sorted(name for name,cmap in mcm.cmap_d.items() if isinstance(cmap, mcolors.ListedColormap))
                raise ValueError(f'Invalid cycle name {cycle!r}. Options are: {", ".join(map(repr, cycles))}')
            if rgbcycle and cycle.lower() == 'colorblind':
                regcolors = colors + [(0.1, 0.1, 0.1)]
            elif mcolors.to_rgb('r') != (1.0,0.0,0.0): # reset
                regcolors = [(0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.75, 0.75, 0.0), (0.75, 0.75, 0.0), (0.0, 0.75, 0.75), (0.0, 0.0, 0.0)]
            else:
                regcolors = [] # no reset necessary
            for code,color in zip('brgmyck', regcolors):
                rgb = mcolors.to_rgb(color)
                mcolors.colorConverter.colors[code] = rgb
                mcolors.colorConverter.cache[code]  = rgb
            kw['patch.facecolor'] = colors[0]
            kw['axes.prop_cycle'] = cycler.cycler('color', colors)

        # Zero linewidth almost always means zero tick length
        elif key == 'linewidth' and _to_points(key, value) == 0:
            ikw, ikw_custom = _get_synced_params('ticklen', 0)
            kw.update(ikw)
            kw_custom.update(ikw_custom)

        # Tick length/major-minor tick length ratio
        elif key in ('ticklen', 'ticklenratio'):
            if key == 'ticklen':
                ticklen = _to_points(key, value)
                ratio = rcParamsShort['ticklenratio']
            else:
                ticklen = rcParamsShort['ticklen']
                ratio = value
            kw['xtick.minor.size'] = ticklen*ratio
            kw['ytick.minor.size'] = ticklen*ratio

        # Spine width/major-minor tick width ratio
        elif key in ('linewidth', 'tickratio'):
            if key == 'linewidth':
                tickwidth = _to_points(key, value)
                ratio = rcParamsShort['tickratio']
            else:
                tickwidth = rcParamsShort['linewidth']
                ratio = value
            kw['xtick.minor.width'] = tickwidth*ratio
            kw['ytick.minor.width'] = tickwidth*ratio

        # Gridline width
        elif key in ('grid.linewidth', 'gridratio'):
            if key == 'grid.linewidth':
                gridwidth = _to_points(key, value)
                ratio = rcParamsShort['gridratio']
            else:
                gridwidth = rcParams['grid.linewidth']
                ratio = value
            kw_custom['gridminor.linewidth'] = gridwidth*ratio

        # Gridline toggling, complicated because of the clunky way this is
        # implemented in matplotlib. There should be a gridminor setting!
        elif key in ('grid', 'gridminor'):
            ovalue = rcParams['axes.grid']
            owhich = rcParams['axes.grid.which']
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
            kw['axes.grid'] = value
            kw['axes.grid.which'] = which

        # Now update linked settings
        value = _to_points(key, value)
        for name in RC_CHILDREN.get(key, ()):
            if name in rcParamsLong:
                kw_custom[name] = value
            else:
                kw[name] = value
    return kw, kw_custom

#-----------------------------------------------------------------------------#
# Class
#-----------------------------------------------------------------------------#
class _mode_mod(object):
    """Helper class that temporarily modifies the getitem mode."""
    def __init__(self, mode):
        self._mode = mode
    def __enter__(self):
        if self._mode is not None:
            self._mode_prev = rc._mode
            object.__setattr__(rc, '_mode', self._mode)
    def __exit__(self, *args):
        if self._mode is not None:
            object.__setattr__(rc, '_mode', self._mode_prev)

class rc_configurator(object):
    """
    Magical abstract class for managing matplotlib `rcParams
    <https://matplotlib.org/users/customizing.html>`__ and additional
    ProPlot :ref:`rcParamsLong` and :ref:`rcParamsShort` settings. When
    initialized, this loads defaults settings plus any user overrides in the
    ``~/.proplotrc`` file. See the `~proplot.rctools` documentation for
    details.
    """
    def __str__(self):
        return type(rcParams).__str__(rcParamsShort) # just show globals
    def __repr__(self):
        return type(rcParams).__repr__(rcParamsShort)
    def __contains__(self, key):
        return key in rcParamsShort or key in rcParamsLong or key in rcParams

    @_counter # about 0.05s
    def __init__(self, local=True):
        """
        Parameters
        ----------
        local : bool, optional
            Whether to load overrides from local and user ``.proplotrc``
            file(s). Default is ``True``.
        """
        # Attributes and style
        object.__setattr__(self, '_mode', 0)
        object.__setattr__(self, '_context', [])
        plt.style.use('default')

        # Update from defaults
        rcParamsLong.clear()
        rcParamsLong.update(defaultParamsLong)
        rcParamsShort.clear()
        rcParamsShort.update(defaultParamsShort)
        for key,value in rcParamsShort.items():
            rc, rc_long = _get_synced_params(key, value)
            rcParams.update(rc)
            rcParamsLong.update(rc_long)

        # Update from files
        if not local:
            return
        for i,file in enumerate(_get_config_paths()):
            if not os.path.exists(file):
                continue
            with open(file) as f:
                try:
                    data = yaml.safe_load(f)
                except yaml.YAMLError as err:
                    print('{file!r} has invalid YAML syntax.')
                    raise err
            for key,value in (data or {}).items():
                try:
                    self[key] = value
                except KeyError:
                    raise RuntimeError(f'{file!r} has invalid key {key!r}.')

    def __enter__(self):
        """Applies settings from the most recent context object."""
        # Get context information
        _, mode, kwargs, cache, restore = self._context[-1] # missing arg is previous mode
        object.__setattr__(self, '_mode', mode)
        def _set_item(rcdict, key, value):
            restore[key] = rcdict[key]
            cache[key] = rcdict[key] = value

        # Apply settings
        for key,value in kwargs.items():
            if key in rcParamsShort:
                rc, rc_long = _get_synced_params(key, value)
                _set_item(rcParamsShort, key, value)
                for ikey, ivalue in rc_long.items():
                    _set_item(rcParamsLong, ikey, ivalue)
                for ikey, ivalue in rc.items():
                    _set_item(rcParams, ikey, ivalue)
            elif key in rcParamsLong:
                _set_item(rcParamsLong, ikey, ivalue)
            elif key in rcParams:
                _set_item(rcParams, ikey, ivalue)

    def __exit__(self, *args):
        """Restores configurator cache to initial state."""
        mode, _, _, _, restore = self._context[-1]
        for key,value in restore.items():
            self[key] = value
        del self._context[-1]
        object.__setattr__(self, '_mode', mode)

    def __delitem__(self, *args):
        """Pseudo-immutability."""
        raise RuntimeError('rc settings cannot be deleted.')

    def __delattr__(self, *args):
        """Pseudo-immutability."""
        raise RuntimeError('rc settings cannot be deleted.')

    def __getattr__(self, attr):
        """Invokes `~rc_configurator.__getitem__`."""
        return self[attr]

    def __getitem__(self, key):
        """Returns `rcParams <https://matplotlib.org/users/customizing.html>`__,
        :ref:`rcParamsLong`, and :ref:`rcParamsShort` settings. If we are in a
        `~rc_configurator.context` block, may return ``None`` if the setting
        is not cached (i.e. if it was not changed by the user)."""
        # Can get a whole bunch of different things
        # Get full dictionary e.g. for rc[None]
        if not key:
            return {**rcParams, **rcParamsLong}

        # Standardize
        # NOTE: If key is invalid, raise error down the line.
        if '.' not in key and key not in rcParamsShort:
            key = RC_NODOTS.get(key, key)

        # Allow for special time-saving modes where we *ignore rcParams*
        # or even *ignore rcParamsLong*.
        mode = self._mode
        if mode == 0:
            kws = (rcParamsShort, rcParamsLong, rcParams)
        elif mode == 1:
            kws = (rcParamsShort, rcParamsLong) # custom only!
        elif mode == 2:
            kws = ()
        else:
            raise KeyError(f'Invalid caching mode {mode!r}.')
        if self._context:
            kws = (self._context[-1][-2], *kws)

        # Get individual property. Will successively index a few different dicts
        # Try to return the value
        for kw in kws:
            try:
                return kw[key]
            except KeyError:
                continue
        # If we were in one of the exclusive modes, return None
        if mode == 0:
            raise KeyError(f'Invalid property name {key!r}.')
        else:
            return None

    def __setattr__(self, attr, value):
        """Invokes `~rc_configurator.__setitem__`."""
        self[attr] = value

    def __setitem__(self, key, value):
        """Sets `rcParams <https://matplotlib.org/users/customizing.html>`__,
        :ref:`rcParamsLong`, and :ref:`rcParamsShort` settings."""
        if '.' not in key and key not in rcParamsShort:
            key = RC_NODOTS.get(key, key)
        if key == 'title.pad':
            key = 'axes.titlepad'
        if key in rcParamsShort:
            rc, rc_long = _get_synced_params(key, value)
            rcParamsShort[key] = value
            for ikey, ivalue in rc.items():
                rcParams[ikey] = ivalue
            for ikey, ivalue in rc_long.items():
                rcParamsLong[ikey] = ivalue
        elif key in rcParamsLong:
            rcParamsLong[key] = _to_points(key, value)
        elif key in rcParams:
            rcParams[key] = _to_points(key, value)
        else:
            raise KeyError(f'Invalid key {key!r}.')

    def _setattr(self, attr, value):
        """Helper function that sets attribute."""
        object.__setattr__(self, attr, value)

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
        # Check input mode
        if cat not in RC_CATEGORIES:
            raise ValueError(f'rc category {cat!r} does not exist. Valid categories are {", ".join(map(repr, RC_CATEGORIES))}.')
        mode = 0 if not cache else self._mode
        if mode == 0:
            kws = (rcParamsLong, rcParams)
        elif mode == 1:
            kws = (rcParamsLong,)
        elif mode == 2:
            kws = ()
        else:
            raise KeyError(f'Invalid caching mode {mode}.')
        if self._context:
            kws = (self._context[-1][-2], *kws)

        # Return params dictionary
        output = {}
        for kw in kws:
            for category,value in kw.items():
                if re.search(f'^{cat}\.', category):
                    subcategory = re.sub(f'^{cat}\.', '', category)
                    if subcategory and '.' not in subcategory:
                        output[subcategory] = value
        return output

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
            Dictionaries of `rc` names and values.
        **kwargs
            `rc` names and values passed as keyword arguments. If the
            name has dots, simply omit them.

        Other parameters
        ----------------
        mode : {0,1,2}, optional
            The `~rc_configurator.__getitem__` mode. Dictates the behavior of
            the `rc` object within a ``with...as`` block when settings are
            requested.

            0. All settings (`rcParams <https://matplotlib.org/users/customizing.html>`__,
               :ref:`rcParamsLong`, and :ref:`rcParamsShort`) are returned,
               whether or not `~rc_configurator.context` has changed them.
            1. Unchanged `rcParams <https://matplotlib.org/users/customizing.html>`__
               return ``None``. :ref:`rcParamsLong` and :ref:`rcParamsShort`
               are returned whether or not `~rc_configurator.context` has
               changed them.  This is used in the `~proplot.axes.Axes.__init__`
               call to `~proplot.axes.Axes.format`. When a lookup returns
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
            raise ValueError(f'Invalid mode {mode!r}.')
        for arg in args:
            if not isinstance(arg, dict):
                raise ValueError('Non-dictionary argument {arg!r}.')
            kwargs.update(arg)
        self._context.append((self._mode, mode, kwargs, {}, {}))

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
        mode = 0 if not cache else None
        with _mode_mod(mode):
            return self[key]

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
        mode = 0 if not cache else None
        output = {}
        with _mode_mod(mode):
            for key,value in props.items():
                item = self[value]
                if item is not None:
                    output[key] = item
        return output

    def update(self, *args, **kwargs):
        """
        Bulk updates settings.

        Parameters
        ----------
        *args : str, dict, or (str, dict)
            The first argument can optionally be a "category" string name,
            in which case all other setting names passed to this function are
            prepended with the string ``cat + '.'``. For example,
            ``rc.update('axes', labelsize=20, titlesize=20)`` changes the
            :rcraw:`axes.labelsize` and :rcraw:`axes.titlesize` properties.

            The first or second argument can also be a dictionary of `rc`
            names and values.
        **kwargs
            `rc` names and values passed as keyword arguments. If the
            name has dots, simply omit them.
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

    def reset(self, **kwargs):
        """
        Resets the configurator to its initial state.

        Parameters
        ----------
        **kwargs
            Passed to `rc_configurator`.
        """
        self.__init__(**kwargs)

# Declare rc object
# WARNING: Must be instantiated after ipython notebook setup! The default
# backend may change some rc settings!
rc = rc_configurator(True)
"""Instance of `rc_configurator`. This is used to change global settings.
See the `~proplot.rctools` documentation for details."""

# Ipython notebook behavior
@_timer
def nb_setup():
    """
    Sets up your iPython workspace, called on import if :rcraw:`nbsetup` is
    ``True``. For all iPython sessions, passes :rcraw:`autoreload` to the useful
    `autoreload <https://ipython.readthedocs.io/en/stable/config/extensions/autoreload.html>`__
    extension. For iPython *notebook* sessions, results in higher-quality inline figures
    and passes :rcraw:`autosave` to the `autosave <https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-matplotlib>`__
    extension.

    See the `~proplot.rctools` documentation for details.
    """
    # Make sure we are in session
    ipython = get_ipython()
    if ipython is None:
        return

    # Only do this if not already loaded -- otherwise will get *recursive*
    # reloading, even with unload_ext command!
    if rcParamsShort['autoreload']:
        if 'autoreload' not in ipython.magics_manager.magics['line']:
            ipython.magic("reload_ext autoreload") # reload instead of load, to avoid annoying message
        ipython.magic("autoreload " + str(rcParamsShort['autoreload'])) # turn on expensive autoreloading

    # Initialize with default 'inline' settings
    # Reset rc object afterwards
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
        autosave = rcParamsShort['autosave']
        if autosave: # capture annoying message + line breaks
            with IPython.utils.io.capture_output():
                ipython.magic("autosave " + str(autosave))
        ipython.magic("config InlineBackend.figure_formats = ['" + rcParamsShort['format'] + "']")
        ipython.magic("config InlineBackend.rc = {}") # no notebook-specific overrides
        ipython.magic("config InlineBackend.close_figures = True") # memory issues
        ipython.magic("config InlineBackend.print_figure_kwargs = {'bbox_inches':None}") # use ProPlot tight layout

# Setup notebook and issue warning
# TODO: Add to list of incompatible backends?
if rcParamsShort['nbsetup']:
    nb_setup()
if rcParams['backend'][:2] == 'nb' or rcParams['backend'] in ('MacOSX',):
    warnings.warn(f'Due to automatic figure resizing, using ProPlot with the {rcParams["backend"]!r} backend may result in unexpected behavior. Try using %matplotlib inline or %matplotlib qt, or just import ProPlot before specifying the backend and ProPlot will automatically load it.')

