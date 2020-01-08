#!/usr/bin/env python3
"""
Utilities for configuring matplotlib and ProPlot global settings.
See :ref:`Configuring proplot` for details.
"""
# NOTE: Make sure to add to docs/configuration.rst when updating or adding
# new settings! Much of this script was adapted from seaborn; see:
# https://github.com/mwaskom/seaborn/blob/master/seaborn/rcmod.py
import re
import os
import numpy as np
import cycler
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
from numbers import Number
from matplotlib import style, rcParams
try:
    import IPython
    from IPython import get_ipython
except ModuleNotFoundError:
    def get_ipython():
        return
from .utils import _warn_proplot, _counter, _benchmark, units

# Disable mathtext "missing glyph" warnings
import matplotlib.mathtext  # noqa
import logging
logger = logging.getLogger('matplotlib.mathtext')
logger.setLevel(logging.ERROR)  # suppress warnings!

__all__ = [
    'rc', 'rc_configurator', 'ipython_autosave',
    'ipython_autoreload', 'ipython_matplotlib',
]

# Dictionaries used to track custom proplot settings
rcParamsShort = {}
rcParamsLong = {}

# Dictionaries containing default settings
defaultParamsShort = {
    'abc': False,
    'align': False,
    'alpha': 1,
    'autoreload': 2,
    'autosave': 30,
    'borders': False,
    'cmap': 'fire',
    'coast': False,
    'color': 'k',
    'cycle': 'colorblind',
    'facecolor': 'w',
    'fontname': 'sans-serif',
    'inlinefmt': 'retina',
    'geogrid': True,
    'grid': True,
    'gridminor': False,
    'gridratio': 0.5,
    'innerborders': False,
    'lakes': False,
    'land': False,
    'large': 10,
    'linewidth': 0.6,
    'lut': 256,
    'margin': 0.0,
    'matplotlib': 'auto',
    'ocean': False,
    'reso': 'lo',
    'rgbcycle': False,
    'rivers': False,
    'share': 3,
    'small': 9,
    'span': True,
    'tickdir': 'out',
    'ticklen': 4.0,
    'ticklenratio': 0.5,
    'tickminor': True,
    'tickpad': 2.0,
    'tickratio': 0.8,
    'tight': True,
}
defaultParamsLong = {
    'abc.border': True,
    'abc.color': 'k',
    'abc.linewidth': 1.5,
    'abc.loc': 'l',  # left side above the axes
    'abc.size': None,  # = large
    'abc.style': 'a',
    'abc.weight': 'bold',
    'axes.facealpha': None,  # if empty, depends on 'savefig.transparent'
    'axes.formatter.timerotation': 90,
    'axes.formatter.zerotrim': True,
    'axes.geogrid': True,
    'axes.gridminor': True,
    'borders.color': 'k',
    'borders.linewidth': 0.6,
    'bottomlabel.color': 'k',
    'bottomlabel.size': None,  # = large
    'bottomlabel.weight': 'bold',
    'coast.color': 'k',
    'coast.linewidth': 0.6,
    'colorbar.extend': '1.3em',
    'colorbar.framealpha': 0.8,
    'colorbar.frameon': True,
    'colorbar.grid': False,
    'colorbar.insetextend': '1em',
    'colorbar.insetlength': '8em',
    'colorbar.insetpad': '0.5em',
    'colorbar.insetwidth': '1.2em',
    'colorbar.length': 1,
    'colorbar.loc': 'right',
    'colorbar.width': '1.5em',
    'geoaxes.edgecolor': None,  # = color
    'geoaxes.facealpha': None,  # = alpha
    'geoaxes.facecolor': None,  # = facecolor
    'geoaxes.linewidth': None,  # = linewidth
    'geogrid.alpha': 0.5,
    'geogrid.color': 'k',
    'geogrid.labels': False,
    'geogrid.labelsize': None,  # = small
    'geogrid.latmax': 90,
    'geogrid.latstep': 20,
    'geogrid.linestyle': ':',
    'geogrid.linewidth': 1.0,
    'geogrid.lonstep': 30,
    'gridminor.alpha': None,  # = grid.alpha
    'gridminor.color': None,  # = grid.color
    'gridminor.linestyle': None,  # = grid.linewidth
    'gridminor.linewidth': None,  # = grid.linewidth x gridratio
    'image.edgefix': True,
    'image.levels': 11,
    'innerborders.color': 'k',
    'innerborders.linewidth': 0.6,
    'lakes.color': 'w',
    'land.color': 'k',
    'leftlabel.color': 'k',
    'leftlabel.size': None,  # = large
    'leftlabel.weight': 'bold',
    'ocean.color': 'w',
    'rightlabel.color': 'k',
    'rightlabel.size': None,  # = large
    'rightlabel.weight': 'bold',
    'rivers.color': 'k',
    'rivers.linewidth': 0.6,
    'subplots.axpad': '1em',
    'subplots.axwidth': '18em',
    'subplots.pad': '0.5em',
    'subplots.panelpad': '0.5em',
    'subplots.panelwidth': '4em',
    'suptitle.color': 'k',
    'suptitle.size': None,  # = large
    'suptitle.weight': 'bold',
    'tick.labelcolor': None,  # = color
    'tick.labelsize': None,  # = small
    'tick.labelweight': 'normal',
    'title.border': True,
    'title.color': 'k',
    'title.linewidth': 1.5,
    'title.loc': 'c',  # centered above the axes
    'title.pad': 3.0,  # copy
    'title.size': None,  # = large
    'title.weight': 'normal',
    'toplabel.color': 'k',
    'toplabel.size': None,  # = large
    'toplabel.weight': 'bold',
}
defaultParams = {
    'axes.grid': True,
    'axes.labelpad': 3.0,
    'axes.titlepad': 3.0,
    'axes.titleweight': 'normal',
    'axes.xmargin': 0.0,
    'axes.ymargin': 0.0,
    'figure.autolayout': False,
    'figure.facecolor': '#f2f2f2',
    'figure.max_open_warning': 0,
    'figure.titleweight': 'bold',
    'font.serif': (
        'New Century Schoolbook',
        'Century Schoolbook L',
        'Utopia',
        'ITC Bookman',
        'Bookman',
        'Nimbus Roman No9 L',
        'Times New Roman',
        'Times',
        'Palatino',
        'Charter',
        'Computer Modern Roman',
        'DejaVu Serif',
        'Bitstream Vera Serif',
        'serif',
    ),
    'font.sans-serif': (
        'Helvetica',
        'Arial',
        'Lucida Grande',
        'Verdana',
        'Geneva',
        'Lucid',
        'Avant Garde',
        'TeX Gyre Heros',
        'DejaVu Sans',
        'Bitstream Vera Sans',
        'Computer Modern Sans Serif',
        'sans-serif',
    ),
    'font.monospace': (
        'Andale Mono',
        'Nimbus Mono L',
        'Courier New',
        'Courier',
        'Fixed',
        'Terminal',
        'Computer Modern Typewriter',
        'DejaVu Sans Mono',
        'Bitstream Vera Sans Mono',
        'monospace',
    ),
    'grid.alpha': 0.1,
    'grid.color': 'k',
    'grid.linestyle': '-',
    'grid.linewidth': 0.6,
    'hatch.color': 'k',
    'hatch.linewidth': 0.6,
    'legend.borderaxespad': 0,
    'legend.borderpad': 0.5,
    'legend.columnspacing': 1.0,
    'legend.fancybox': False,
    'legend.framealpha': 0.8,
    'legend.frameon': True,
    'legend.handlelength': 1.5,
    'legend.handletextpad': 0.5,
    'legend.labelspacing': 0.5,
    'lines.linewidth': 1.3,
    'lines.markersize': 3.0,
    'mathtext.fontset': 'custom',
    'mathtext.default': 'regular',
    'savefig.bbox': 'standard',
    'savefig.directory': '',
    'savefig.dpi': 300,
    'savefig.facecolor': 'white',
    'savefig.format': 'pdf',
    'savefig.pad_inches': 0.0,
    'savefig.transparent': True,
    'text.usetex': False,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
}

# "Global" settings and the lower-level settings they change
_rc_children = {
    'cmap': (
        'image.cmap',
    ),
    'lut': (
        'image.lut',
    ),
    'alpha': (  # this is a custom setting
        'axes.facealpha', 'geoaxes.facealpha',
    ),
    'facecolor': (
        'axes.facecolor', 'geoaxes.facecolor'
    ),
    'fontname': (
        'font.family',
    ),
    'color': (  # change the 'color' of an axes
        'axes.edgecolor', 'geoaxes.edgecolor', 'axes.labelcolor',
        'tick.labelcolor', 'hatch.color', 'xtick.color', 'ytick.color'
    ),
    'small': (  # the 'small' fonts
        'font.size', 'tick.labelsize', 'xtick.labelsize', 'ytick.labelsize',
        'axes.labelsize', 'legend.fontsize', 'geogrid.labelsize'
    ),
    'large': (  # the 'large' fonts
        'abc.size', 'figure.titlesize',
        'axes.titlesize', 'suptitle.size', 'title.size',
        'leftlabel.size', 'toplabel.size',
        'rightlabel.size', 'bottomlabel.size'
    ),
    'linewidth': (
        'axes.linewidth', 'geoaxes.linewidth', 'hatch.linewidth',
        'xtick.major.width', 'ytick.major.width'
    ),
    'margin': (
        'axes.xmargin', 'axes.ymargin'
    ),
    'grid': (
        'axes.grid',
    ),
    'gridminor': (
        'axes.gridminor',
    ),
    'geogrid': (
        'axes.geogrid',
    ),
    'ticklen': (
        'xtick.major.size', 'ytick.major.size'
    ),
    'tickdir': (
        'xtick.direction', 'ytick.direction'
    ),
    'labelpad': (
        'axes.labelpad',
    ),
    'titlepad': (
        'axes.titlepad',
    ),
    'tickpad': (
        'xtick.major.pad', 'xtick.minor.pad',
        'ytick.major.pad', 'ytick.minor.pad'
    ),
    'grid.color': (
        'gridminor.color',
    ),
    'grid.linewidth': (
        'gridminor.linewidth',
    ),
    'grid.linestyle': (
        'gridminor.linestyle',
    ),
    'grid.alpha': (
        'gridminor.alpha',
    ),
}

# Mapping of settings without "dots" to their full names. This lets us pass
# all settings as kwargs, e.g. ax.format(landcolor='b') instead of the much
# more verbose ax.format(rc_kw={'land.color':'b'}).
# WARNING: rcParamsShort has to be in here because Axes.format() only checks
# _rc_nodots to filter out the rc kwargs!
_rc_nodots = {
    name.replace('.', ''): name
    for names in (defaultParamsShort, defaultParamsLong, rcParams)
    for name in names.keys()
}

# Category names, used for returning dicts of subcategory properties
_rc_categories = {
    *(
        re.sub(r'\.[^.]*$', '', name)
        for names in (defaultParamsLong, rcParams)
        for name in names.keys()
    ),
    *(
        re.sub(r'\..*$', '', name)
        for names in (defaultParamsLong, rcParams)
        for name in names.keys()
    )
}


def _get_config_paths():
    """Return a list of configuration file paths in reverse order of
    precedence."""
    # Local configuration
    idir = os.getcwd()
    paths = []
    while idir:  # not empty string
        ipath = os.path.join(idir, '.proplotrc')
        if os.path.exists(ipath):
            paths.append(ipath)
        ndir = os.path.dirname(idir)
        if ndir == idir:  # root
            break
        idir = ndir
    paths = paths[::-1]  # sort from decreasing to increasing importantce
    # Home configuration
    ipath = os.path.join(os.path.expanduser('~'), '.proplotrc')
    if os.path.exists(ipath) and ipath not in paths:
        paths.insert(0, ipath)
    return paths


def _get_synced_params(key, value):
    """Return dictionaries for updating the `rcParamsShort`, `rcParamsLong`,
    and `rcParams` properties associated with this key."""
    kw = {}  # builtin properties that global setting applies to
    kw_long = {}  # custom properties that global setting applies to
    kw_short = {}  # short name properties

    # Skip full name keys
    key = _sanitize_key(key)
    if '.' in key:
        pass

    # Cycler
    elif key in ('cycle', 'rgbcycle'):
        if key == 'rgbcycle':
            cycle, rgbcycle = rcParamsShort['cycle'], value
        else:
            cycle, rgbcycle = value, rcParamsShort['rgbcycle']
        try:
            colors = mcm.cmap_d[cycle].colors
        except (KeyError, AttributeError):
            cycles = sorted(
                name for name, cmap in mcm.cmap_d.items()
                if isinstance(cmap, mcolors.ListedColormap)
            )
            raise ValueError(
                f'Invalid cycle name {cycle!r}. Options are: '
                ', '.join(map(repr, cycles)) + '.'
            )
        if rgbcycle and cycle.lower() == 'colorblind':
            regcolors = colors + [(0.1, 0.1, 0.1)]
        elif mcolors.to_rgb('r') != (1.0, 0.0, 0.0):  # reset
            regcolors = [
                (0.0, 0.0, 1.0),
                (1.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                (0.75, 0.75, 0.0),
                (0.75, 0.75, 0.0),
                (0.0, 0.75, 0.75),
                (0.0, 0.0, 0.0)
            ]
        else:
            regcolors = []  # no reset necessary
        for code, color in zip('brgmyck', regcolors):
            rgb = mcolors.to_rgb(color)
            mcolors.colorConverter.colors[code] = rgb
            mcolors.colorConverter.cache[code] = rgb
        kw['patch.facecolor'] = colors[0]
        kw['axes.prop_cycle'] = cycler.cycler('color', colors)

    # Zero linewidth almost always means zero tick length
    elif key == 'linewidth' and _to_points(key, value) == 0:
        _, ikw_long, ikw = _get_synced_params('ticklen', 0)
        kw.update(ikw)
        kw_long.update(ikw_long)

    # Tick length/major-minor tick length ratio
    elif key in ('ticklen', 'ticklenratio'):
        if key == 'ticklen':
            ticklen = _to_points(key, value)
            ratio = rcParamsShort['ticklenratio']
        else:
            ticklen = rcParamsShort['ticklen']
            ratio = value
        kw['xtick.minor.size'] = ticklen * ratio
        kw['ytick.minor.size'] = ticklen * ratio

    # Spine width/major-minor tick width ratio
    elif key in ('linewidth', 'tickratio'):
        if key == 'linewidth':
            tickwidth = _to_points(key, value)
            ratio = rcParamsShort['tickratio']
        else:
            tickwidth = rcParamsShort['linewidth']
            ratio = value
        kw['xtick.minor.width'] = tickwidth * ratio
        kw['ytick.minor.width'] = tickwidth * ratio

    # Gridline width
    elif key in ('grid.linewidth', 'gridratio'):
        if key == 'grid.linewidth':
            gridwidth = _to_points(key, value)
            ratio = rcParamsShort['gridratio']
        else:
            gridwidth = rcParams['grid.linewidth']
            ratio = value
        kw_long['gridminor.linewidth'] = gridwidth * ratio

    # Gridline toggling, complicated because of the clunky way this is
    # implemented in matplotlib. There should be a gridminor setting!
    elif key in ('grid', 'gridminor'):
        ovalue = rcParams['axes.grid']
        owhich = rcParams['axes.grid.which']
        # Instruction is to turn off gridlines
        if not value:
            # Gridlines are already off, or they are on for the particular
            # ones that we want to turn off. Instruct to turn both off.
            if not ovalue or (key == 'grid' and owhich == 'major') or (
                    key == 'gridminor' and owhich == 'minor'):
                which = 'both'  # disable both sides
            # Gridlines are currently on for major and minor ticks, so we
            # instruct to turn on gridlines for the one we *don't* want off
            elif owhich == 'both':  # and ovalue is True, as we already tested
                # if gridminor=False, enable major, and vice versa
                value = True
                which = 'major' if key == 'gridminor' else 'minor'
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
            if owhich == 'both' or (key == 'grid' and owhich == 'minor') or (
                    key == 'gridminor' and owhich == 'major'):
                which = 'both'
            # Gridlines are off for both, or off for the ones that we
            # don't want to turn on. We can just turn on these ones.
            else:
                which = owhich
        kw['axes.grid'] = value
        kw['axes.grid.which'] = which

    # Now update linked settings
    value = _to_points(key, value)
    if key in rcParamsShort:
        kw_short[key] = value
    elif key in rcParamsLong:
        kw_long[key] = value
    elif key in rcParams:
        kw[key] = value
    else:
        raise KeyError(f'Invalid key {key!r}.')
    for name in _rc_children.get(key, ()):
        if name in rcParamsLong:
            kw_long[name] = value
        else:
            kw[name] = value
    return kw_short, kw_long, kw


def _sanitize_key(key):
    """Ensure string and convert keys with omitted dots."""
    if not isinstance(key, str):
        raise KeyError(f'Invalid key {key!r}. Must be string.')
    if '.' not in key and key not in rcParamsShort:  # speedup
        key = _rc_nodots.get(key, key)
    return key.lower()


def _to_points(key, value):
    """Convert certain rc keys to the units "points"."""
    # TODO: Incorporate into more sophisticated validation system
    # See: https://matplotlib.org/users/customizing.html, all props matching
    # the below strings use the units 'points', except custom categories!
    if (
        isinstance(value, str)
        and key.split('.')[0] not in ('colorbar', 'subplots')
        and re.match('^.*(width|space|size|pad|len|small|large)$', key)
    ):
        value = units(value, 'pt')
    return value


def _update_from_file(file):
    """
    Apply updates from a file. This is largely copied from matplotlib.

    Parameters
    ----------
    file : str
        The path.
    """
    cnt = 0
    file = os.path.expanduser(file)
    added = set()
    with open(file, 'r') as fd:
        for line in fd:
            # Read file
            cnt += 1
            stripped = line.split('#', 1)[0].strip()
            if not stripped:
                continue
            pair = stripped.split(':', 1)
            if len(pair) != 2:
                _warn_proplot(
                    f'Illegal line #{cnt} in file {file!r}:\n{line!r}"'
                )
                continue
            key, value = pair
            key = key.strip()
            value = value.strip()
            if key in added:
                _warn_proplot(
                    f'Duplicate key {key!r} on line #{cnt} in file {file!r}.'
                )
            added.add(key)

            # *Very primitive* type conversion system. Just check proplot
            # settings (they are all simple/scalar) and leave rcParams alone.
            # TODO: Add built-in validation by making special RcParamsLong
            # and RcParamsShort classes just like matplotlib RcParams
            if key in rcParamsShort or key in rcParamsLong:
                if not value:
                    value = None  # older proplot versions supported this
                elif value in ('True', 'False', 'None'):
                    value = eval(value)  # rare case where eval is o.k.
                else:
                    try:
                        # int-float distinction does not matter in python3
                        value = float(value)
                    except ValueError:
                        pass

            # Add to dictionaries
            try:
                rc_short, rc_long, rc = _get_synced_params(key, value)
            except KeyError:
                _warn_proplot(
                    f'Invalid key {key!r} on line #{cnt} in file {file!r}.'
                )
            else:
                rcParamsShort.update(rc_short)
                rcParamsLong.update(rc_long)
                rcParams.update(rc)


def _write_defaults(filename, comment=True, overwrite=False):
    """
    Save a file to the specified path containing the default `rc` settings.

    Parameters
    ----------
    filename : str
        The path.
    comment : bool, optional
        Whether to "comment out" each setting.
    overwrite : bool, optional
        Whether to overwrite existing files.
    """
    def _tabulate(rcdict):
        string = ''
        prefix = '# ' if comment else ''
        maxlen = max(map(len, rcdict))
        NoneType = type(None)
        for key, value in rcdict.items():
            if isinstance(value, cycler.Cycler):  # special case!
                value = repr(value)
            elif isinstance(value, (str, Number, NoneType)):
                value = str(value)
            elif isinstance(value, (list, tuple)) and all(
                isinstance(val, (str, Number)) for val in value
            ):
                value = ', '.join(str(val) for val in value)
            else:
                raise ValueError(
                    f'Failed to write rc setting {key} = {value!r}. '
                    'Must be string, number, or list or tuple thereof, '
                    'or None or a cycler.'
                )
            space = ' ' * (maxlen - len(key) + 1)
            string += f'{prefix}{key}:{space}{value}\n'
        return string.strip()

    # Fill empty defaultParamsLong values with rcDefaultParamsShort
    # They are only allowed to be None in the *default dictionaries* because
    # they are immediately overwritten. However if users try to set them as
    # None in a .proplotrc file, may trigger error down the line.
    rc_parents = {
        child: parent
        for parent, children in _rc_children.items()
        for child in children
    }
    defaultParamsLong_filled = defaultParamsLong.copy()
    for key, value in defaultParamsLong.items():
        if value is None:
            try:
                defaultParamsLong_filled[key] = rc_parents[key]
            except KeyError:
                raise RuntimeError(
                    f'rcParamsLong param {key!r} has default value of None '
                    'but has no rcParmsShort parent!'
                )

            value = defaultParamsShort

    with open(filename, 'w') as f:
        f.write(f"""
#---------------------------------------------------------------------
# Use this file to change the default proplot and matplotlib settings
# The syntax is mostly the same as for matplotlibrc files
# For descriptions of each setting see:
# https://proplot.readthedocs.io/en/latest/rctools.html
# https://matplotlib.org/3.1.1/tutorials/introductory/customizing.html
#---------------------------------------------------------------------
# ProPlot short name settings
{_tabulate(defaultParamsShort)}

# ProPlot long name settings
{_tabulate(defaultParamsLong_filled)}

# Matplotlib settings
{_tabulate(defaultParams)}
""".strip())


class rc_configurator(object):
    """
    Magical abstract class for managing matplotlib
    `rcParams <https://matplotlib.org/users/customizing.html>`__
    and additional ProPlot :ref:`rcParamsLong` and :ref:`rcParamsShort`
    settings. When initialized, this loads defaults settings plus any user
    overrides in the ``~/.proplotrc`` file. See the `~proplot.rctools`
    documentation for details.
    """
    def __contains__(self, key):
        return key in rcParamsShort or key in rcParamsLong or key in rcParams

    def __iter__(self):
        for key in sorted((*rcParamsShort, *rcParamsLong, *rcParams)):
            yield key

    def __repr__(self):
        rcdict = type('rc', (dict,), {})(rcParamsShort)
        string = type(rcParams).__repr__(rcdict)
        indent = ' ' * 4  # indent is rc({
        return string.strip(
            '})') + f'\n{indent}... (rcParams) ...\n{indent}}})'

    def __str__(self):  # encapsulate params in temporary class
        rcdict = type('rc', (dict,), {})(rcParamsShort)
        string = type(rcParams).__str__(rcdict)
        return string + '\n... (rcParams) ...'

    @_counter  # about 0.05s
    def __init__(self, local=True):
        """
        Parameters
        ----------
        local : bool, optional
            Whether to load overrides from local and user ``.proplotrc``
            file(s). Default is ``True``.
        """
        # Attributes and style
        object.__setattr__(self, '_context', [])
        with _benchmark('  use'):
            style.use('default')

        # Update from defaults
        rcParams.update(defaultParams)
        rcParamsLong.clear()
        rcParamsLong.update(defaultParamsLong)
        rcParamsShort.clear()
        rcParamsShort.update(defaultParamsShort)
        for rcdict in (rcParamsShort, rcParamsLong):
            for key, value in rcdict.items():
                _, rc_long, rc = _get_synced_params(key, value)
                rcParamsLong.update(rc_long)
                rcParams.update(rc)

        # Update from files
        if not local:
            return
        for i, file in enumerate(_get_config_paths()):
            if not os.path.exists(file):
                continue
            _update_from_file(file)

    def __enter__(self):
        """Apply settings from the most recent context block."""
        if not self._context:
            raise RuntimeError(
                f'rc object must be initialized with rc.context().'
            )
        *_, kwargs, cache, restore = self._context[-1]

        def _update(rcdict, newdict):
            for key, value in newdict.items():
                restore[key] = rcdict[key]
                rcdict[key] = cache[key] = value
        for key, value in kwargs.items():
            rc_short, rc_long, rc = _get_synced_params(key, value)
            _update(rcParamsShort, rc_short)
            _update(rcParamsLong, rc_long)
            _update(rcParams, rc)

    def __exit__(self, *args):
        """Restore settings from the most recent context block."""
        if not self._context:
            raise RuntimeError(
                f'rc object must be initialized with rc.context().'
            )
        *_, restore = self._context[-1]
        for key, value in restore.items():
            rc_short, rc_long, rc = _get_synced_params(key, value)
            rcParamsShort.update(rc_short)
            rcParamsLong.update(rc_long)
            rcParams.update(rc)
        del self._context[-1]

    def __delitem__(self, *args):
        """Raise an error. This enforces pseudo-immutability."""
        raise RuntimeError('rc settings cannot be deleted.')

    def __delattr__(self, *args):
        """Raise an error. This enforces pseudo-immutability."""
        raise RuntimeError('rc settings cannot be deleted.')

    def __getattr__(self, attr):
        """Pass the attribute to `~rc_configurator.__getitem__` and return
        the result."""
        if attr[:1] == '_':
            return super().__getattr__(attr)
        else:
            return self[attr]

    def __getitem__(self, key):
        """Return an `rcParams \
<https://matplotlib.org/users/customizing.html>`__,
        :ref:`rcParamsLong`, or :ref:`rcParamsShort` setting."""
        key = _sanitize_key(key)
        for kw in (rcParamsShort, rcParamsLong, rcParams):
            try:
                return kw[key]
            except KeyError:
                continue
        raise KeyError(f'Invalid setting name {key!r}.')

    def __setattr__(self, attr, value):
        """Pass the attribute and value to `~rc_configurator.__setitem__`."""
        self[attr] = value

    def __setitem__(self, key, value):
        """Modify an `rcParams \
<https://matplotlibcorg/users/customizing.html>`__,
        :ref:`rcParamsLong`, and :ref:`rcParamsShort` setting(s)."""
        if key == 'matplotlib':
            return ipython_matplotlib(value)
        elif key == 'autosave':
            return ipython_autosave(value)
        elif key == 'autoreload':
            return ipython_autoreload(value)
        rc_short, rc_long, rc = _get_synced_params(key, value)
        rcParamsShort.update(rc_short)
        rcParamsLong.update(rc_long)
        rcParams.update(rc)

    def _get_item(self, key, mode=None):
        """As with `~rc_configurator.__getitem__` but the search is limited
        based on the context mode and ``None`` is returned if the key is not
        found in the dictionaries."""
        if mode is None:
            mode = min((context[0] for context in self._context), default=0)
        caches = (context[2] for context in self._context)
        if mode == 0:
            rcdicts = (*caches, rcParamsShort, rcParamsLong, rcParams)
        elif mode == 1:
            rcdicts = (*caches, rcParamsShort, rcParamsLong)  # custom only!
        elif mode == 2:
            rcdicts = (*caches,)
        else:
            raise KeyError(f'Invalid caching mode {mode!r}.')
        for rcdict in rcdicts:
            if not rcdict:
                continue
            try:
                return rcdict[key]
            except KeyError:
                continue
        if mode == 0:
            raise KeyError(f'Invalid setting name {key!r}.')
        else:
            return

    def category(self, cat, *, trimcat=True, context=False):
        """
        Return a dictionary of settings beginning with the substring
        ``cat + '.'``.

        Parameters
        ----------
        cat : str, optional
            The `rc` setting category.
        trimcat : bool, optional
            Whether to trim ``cat`` from the key names in the output
            dictionary. Default is ``True``.
        context : bool, optional
            If ``True``, then each category setting that is not found in the
            context mode dictionaries is omitted from the output dictionary.
            See `~rc_configurator.context`.
        """
        if cat not in _rc_categories:
            raise ValueError(
                f'Invalid rc category {cat!r}. Valid categories are '
                ', '.join(map(repr, _rc_categories)) + '.'
            )
        kw = {}
        mode = 0 if not context else None
        for rcdict in (rcParamsLong, rcParams):
            for key in rcdict:
                if not re.match(fr'\A{cat}\.[^.]+\Z', key):
                    continue
                value = self._get_item(key, mode)
                if value is None:
                    continue
                if trimcat:
                    key = re.sub(fr'\A{cat}\.', '', key)
                kw[key] = value
        return kw

    def context(self, *args, mode=0, **kwargs):
        """
        Temporarily modify the rc settings in a "with as" block.

        This is used by ProPlot internally but may also be useful for power
        users. It was invented to prevent successive calls to
        `~proplot.axes.Axes.format` from constantly looking up and
        re-applying unchanged settings. Testing showed that these gratuitous
        `rcParams <https://matplotlib.org/users/customizing.html>`__
        lookups and artist updates increased runtime by seconds, even for
        relatively simple plots. It also resulted in overwriting previous
        rc changes with the default values upon subsequent calls to
        `~proplot.axes.Axes.format`.

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
            The context mode. Dictates the behavior of `~rc_configurator.get`,
            `~rc_configurator.fill`, and `~rc_configurator.category` within a
            "with as" block when called with ``context=True``. The options are
            as follows.

            0. All settings (`rcParams \
<https://matplotlib.org/users/customizing.html>`__,
               :ref:`rcParamsLong`, and :ref:`rcParamsShort`) are returned,
               whether or not `~rc_configurator.context` has changed them.
            1. Unchanged `rcParams \
<https://matplotlib.org/users/customizing.html>`__
               return ``None``. :ref:`rcParamsLong` and :ref:`rcParamsShort`
               are returned whether or not `~rc_configurator.context` has
               changed them.  This is used in the `~proplot.axes.Axes.__init__`
               call to `~proplot.axes.Axes.format`. When a lookup returns
               ``None``, `~proplot.axes.Axes.format` does not apply it.
            2. All unchanged settings return ``None``. This is used during user
               calls to `~proplot.axes.Axes.format`.

        Example
        -------
        The below applies settings to axes in a specific figure using
        `~rc_configurator.context`.

        >>> import proplot as plot
        >>> with plot.rc.context(linewidth=2, ticklen=5):
        ...     f, ax = plot.subplots()
        ...     ax.plot(data)

        By contrast, the below applies settings to a specific axes using
        `~proplot.axes.Axes.format`.

        >>> import proplot as plot
        >>> f, ax = plot.subplots()
        >>> ax.format(linewidth=2, ticklen=5)

        """
        if mode not in range(3):
            raise ValueError(f'Invalid mode {mode!r}.')
        for arg in args:
            if not isinstance(arg, dict):
                raise ValueError('Non-dictionary argument {arg!r}.')
            kwargs.update(arg)
        self._context.append((mode, kwargs, {}, {}))
        return self

    def dict(self):
        """
        Return a raw dictionary of all settings.
        """
        output = {}
        for key in sorted((*rcParamsShort, *rcParamsLong, *rcParams)):
            output[key] = self[key]
        return output

    def get(self, key, *, context=False):
        """
        Return a single setting.

        Parameters
        ----------
        key : str
            The setting name.
        context : bool, optional
            If ``True``, then ``None`` is returned if the setting is not found
            in the context mode dictionaries. See `~rc_configurator.context`.
        """
        mode = 0 if not context else None
        return self._get_item(key, mode)

    def fill(self, props, *, context=False):
        """
        Return a dictionary filled with settings whose names match the
        string values in the input dictionary.

        Parameters
        ----------
        props : dict-like
            Dictionary whose values are `rc` setting names.
        context : bool, optional
            If ``True``, then each setting that is not found in the
            context mode dictionaries is omitted from the output dictionary.
            See `~rc_configurator.context`.
        """
        kw = {}
        mode = 0 if not context else None
        for key, value in props.items():
            item = self._get_item(value, mode)
            if item is not None:
                kw[key] = item
        return kw

    def items(self):
        """
        Return an iterator that loops over all setting names and values.
        Same as `dict.items`.
        """
        for key in self:
            yield key, self[key]

    def keys(self):
        """
        Return an iterator that loops over all setting names.
        Same as `dict.items`.
        """
        for key in self:
            yield key

    def update(self, *args, **kwargs):
        """
        Update several settings at once with a dictionary and/or
        keyword arguments.

        Parameters
        ----------
        *args : str, dict, or (str, dict), optional
            A dictionary containing `rc` keys and values. You can also
            pass a "category" name as the first argument, in which case all
            settings are prepended with ``'category.'``. For example,
            ``rc.update('axes', labelsize=20, titlesize=20)`` changes the
            :rcraw:`axes.labelsize` and :rcraw:`axes.titlesize` properties.
        **kwargs, optional
            `rc` keys and values passed as keyword arguments. If the
            name has dots, simply omit them.
        """
        # Parse args
        kw = {}
        prefix = ''
        if len(args) > 2:
            raise ValueError(
                f'rc.update() accepts 1-2 arguments, got {len(args)}. Usage '
                'is rc.update(kw), rc.update(category, kw), '
                'rc.update(**kwargs), or rc.update(category, **kwargs).'
            )
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
        for key, value in kw.items():
            self[prefix + key] = value

    def reset(self, **kwargs):
        """
        Reset the configurator to its initial state.

        Parameters
        ----------
        **kwargs
            Passed to `rc_configurator`.
        """
        self.__init__(**kwargs)

    def values(self):
        """
        Return an iterator that loops over all setting values.
        Same as `dict.values`.
        """
        for key in self:
            yield self[key]


def ipython_matplotlib(backend=None, fmt=None):
    """
    Set up the `matplotlib backend \
<https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-matplotlib>`__
    for ipython sessions and apply the following ``%config InlineBackend``
    magic commands.

    .. code-block:: ipython

        %config InlineBackend.figure_formats = fmt
        %config InlineBackend.rc = {}  # never override my rc settings!
        %config InlineBackend.close_figures = True  # memory issues
        %config InlineBackend.print_figure_kwargs = {'bbox_inches': None}  \
# we use our own tight layout algorithm

    This must be called *before drawing any figures*! For some ipython
    sessions (e.g. terminals) the backend can only be changed by adding
    ``matplotlib: backend`` to your ``.proplotrc`` file.
    See :ref:`Configuring proplot` for details.

    Parameters
    ----------
    backend : str, optional
        The backend name. The default is ``'auto'``, which applies
        ``%matplotlib inline`` for notebooks and ``%matplotlib qt`` for
        all other sessions.

        Note that when using the qt backend on macOS, you may want to prevent
        "tabbed" figure windows by navigating to Settings...Dock and changing
        "Prefer tabs when opening documents" to "Manually" (see \
`Issue #13164 <https://github.com/matplotlib/matplotlib/issues/13164>`__).
    fmt : str or list of str, optional
        The inline backend file format(s). Valid formats include ``'jpg'``,
        ``'png'``, ``'svg'``, ``'pdf'``, and ``'retina'``. This is ignored
        for non-inline backends.
    """  # noqa
    # Bail out
    ipython = get_ipython()
    backend = backend or rcParamsShort['matplotlib']
    if ipython is None or backend is None:
        return

    # Default behavior dependent on type of ipython session
    # See: https://stackoverflow.com/a/22424821/4970632
    ibackend = backend
    if backend == 'auto':
        if 'IPKernelApp' in getattr(get_ipython(), 'config', ''):
            ibackend = 'inline'
        else:
            ibackend = 'qt'
    try:
        ipython.magic('matplotlib ' + ibackend)
        if 'rc' in globals():  # should always be True, but just in case
            rc.reset()
    except KeyError:
        if backend != 'auto':
            _warn_proplot(f'{"%matplotlib " + backend!r} failed.')

    # Configure inline backend no matter what type of session this is
    # Should be silently ignored for terminal ipython sessions
    fmt = fmt or rcParamsShort['inlinefmt']
    if isinstance(fmt, str):
        fmt = [fmt]
    elif np.iterable(fmt):
        fmt = list(fmt)
    else:
        raise ValueError(
            f'Invalid inline backend format {fmt!r}. '
            'Must be string or list thereof.'
        )
    ipython.magic(f'config InlineBackend.figure_formats = {fmt!r}')
    ipython.magic('config InlineBackend.rc = {}')  # no notebook overrides
    ipython.magic('config InlineBackend.close_figures = True')  # memory issues
    ipython.magic(  # use ProPlot tight layout instead
        'config InlineBackend.print_figure_kwargs = {"bbox_inches":None}'
    )


def ipython_autoreload(autoreload=None):
    """
    Set up the `ipython autoreload utility \
<https://ipython.readthedocs.io/en/stable/config/extensions/autoreload.html>`__
    by running the following ipython magic.

    .. code-block:: ipython

        %autoreload autoreload

    This is called on import by default. Add ``autoreload:`` to your
    ``.proplotrc`` to disable. See :ref:`Configuring proplot` for details.

    Parameters
    ----------
    autoreload : float, optional
        The autoreload level. Default is :rc:`autoreload`.
    """  # noqa
    ipython = get_ipython()
    autoreload = autoreload or rcParamsShort['autoreload']
    if ipython is None or autoreload is None:
        return
    if 'autoreload' not in ipython.magics_manager.magics['line']:
        with IPython.utils.io.capture_output():  # capture annoying message
            ipython.magic('load_ext autoreload')
    ipython.magic('autoreload ' + str(autoreload))


def ipython_autosave(autosave=None):
    """
    Set up the `ipython autosave utility \
<https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-matplotlib>`__
    by running the following ipython magic.

    .. code-block:: ipython

        %autosave autosave

    This is called on import by default. Add ``autosave:`` to your
    ``.proplotrc`` to disable. See :ref:`Configuring proplot` for details.

    Parameters
    ----------
    autosave : float, optional
        The autosave interval in seconds. Default is :rc:`autosave`.
    """  # noqa
    ipython = get_ipython()
    autosave = autosave or rcParamsShort['autosave']
    if ipython is None or autosave is None:
        return
    with IPython.utils.io.capture_output():  # capture annoying message
        try:
            ipython.magic('autosave ' + str(autosave))
        except IPython.core.error.UsageError:
            pass


# Write defaults
_user_rc_file = os.path.join(os.path.expanduser('~'), '.proplotrc')
if not os.path.exists(_user_rc_file):
    _write_defaults(_user_rc_file)

#: Instance of `rc_configurator`. This is used to change global settings.
#: See :ref:`Configuring proplot` for details.
rc = rc_configurator()

# Manually call setup functions after rc has been instantiated
# We cannot call these inside rc.__init__ because ipython_matplotlib may
# need to reset the configurator to overwrite backend-imposed settings!
ipython_matplotlib()
ipython_autoreload()
ipython_autosave()
