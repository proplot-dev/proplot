#!/usr/bin/env python3
"""
Utilities for configuring matplotlib and ProPlot global settings.
See :ref:`Configuring proplot` for details.
"""
# TODO: Add 'style' option that overrides .proplotrc
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
from matplotlib import style, rcParams
try:
    import IPython
    from IPython import get_ipython
except ModuleNotFoundError:
    get_ipython = lambda: None
__all__ = [
    'rc', 'rc_configurator', 'autosave_setup',
    'autoreload_setup', 'backend_setup'
    ]

# Initialize
defaultParamsShort = {
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
    'tickminor':    True,
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
    }
defaultParamsLong = {
    'title.loc':                   'c', # centered above the axes
    'title.pad':                   3.0, # copy
    'abc.loc':                     'l', # left side above the axes
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
    'axes.facealpha':              None, # if empty, depends on 'savefig.transparent' setting
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
    'geoaxes.facealpha':           None, # = alpha
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
    }
defaultParams = {
    'axes.titleweight':        'normal',
    'axes.xmargin':            0.0,
    'axes.ymargin':            0.0,
    'axes.grid':               True,
    'axes.labelpad':           3.0,
    'axes.titlepad':           3.0,
    'figure.dpi':              90,
    'figure.facecolor':        '#f2f2f2',
    'figure.autolayout':       False,
    'figure.titleweight':      'bold',
    'figure.max_open_warning': 0,
    'grid.color':              'k',
    'grid.alpha':              0.1,
    'grid.linewidth':          0.6,
    'grid.linestyle':          '-',
    'hatch.color':             'k',
    'hatch.linewidth':         0.6,
    'legend.frameon':          True,
    'legend.framealpha':       0.8,
    'legend.fancybox':         False,
    'legend.labelspacing':     0.5,
    'legend.handletextpad':    0.5,
    'legend.handlelength':     1.5,
    'legend.columnspacing':    1.0,
    'legend.borderpad':        0.5,
    'legend.borderaxespad':    0,
    'lines.linewidth':         1.3,
    'lines.markersize':        3.0,
    'mathtext.bf':             'sans:bold',
    'mathtext.it':             'sans:it',
    'mathtext.default':        'regular',
    'savefig.directory':       '',
    'savefig.dpi':             300,
    'savefig.facecolor':       'white',
    'savefig.transparent':     True,
    'savefig.format':          'pdf',
    'savefig.bbox':            'standard',
    'savefig.pad_inches':      0.0,
    'xtick.minor.visible':     True,
    'ytick.minor.visible':     True,
    }
rcParamsShort = {}
rcParamsLong = {}

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
        return string.strip()
    with open(_rc_file, 'x') as f:
        f.write(f"""
#------------------------------------------------------
# Use this file to customize settings
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
""".strip())

# "Global" settings and the lower-level settings they change
# NOTE: This whole section, declaring dictionaries and sets, takes 1ms
RC_CHILDREN = {
    'fontname':  ('font.family',),
    'cmap':      ('image.cmap',),
    'lut':       ('image.lut',),
    'alpha':     ('axes.facealpha', 'geoaxes.facealpha'), # this is a custom setting
    'facecolor': ('axes.facecolor', 'geoaxes.facecolor'),
    'color':     ('axes.edgecolor', 'geoaxes.edgecolor', 'axes.labelcolor', 'tick.labelcolor', 'hatch.color', 'xtick.color', 'ytick.color'), # change the 'color' of an axes
    'small':     ('font.size', 'tick.labelsize', 'xtick.labelsize', 'ytick.labelsize', 'axes.labelsize', 'legend.fontsize', 'geogrid.labelsize'), # the 'small' fonts
    'large':     ('abc.size', 'figure.titlesize', 'axes.titlesize', 'suptitle.size', 'title.size', 'leftlabel.size', 'toplabel.size', 'rightlabel.size', 'bottomlabel.size'), # the 'large' fonts
    'linewidth': ('axes.linewidth', 'geoaxes.linewidth', 'hatch.linewidth', 'xtick.major.width', 'ytick.major.width'),
    'margin':    ('axes.xmargin', 'axes.ymargin'),
    'grid':      ('axes.grid',),
    'gridminor': ('axes.gridminor',),
    'geogrid':   ('axes.geogrid',),
    'ticklen':   ('xtick.major.size', 'ytick.major.size'),
    'tickdir':   ('xtick.direction', 'ytick.direction'),
    'tickminor': ('xtick.minor.visible', 'ytick.minor.visible'),
    'tickpad':   ('xtick.major.pad', 'xtick.minor.pad', 'ytick.major.pad', 'ytick.minor.pad'),
    'title.pad': ('axes.titlepad',),
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
    """Convert certain rc keys to the units "points"."""
    # See: https://matplotlib.org/users/customizing.html, all props matching
    # the below strings use the units 'points', and custom categories are in
    if (isinstance(value,str) and key.split('.')[0] not in ('colorbar','subplots')
        and re.match('^.*(width|space|size|pad|len|small|large)$', key)):
        value = utils.units(value, 'pt')
    return value

def _get_config_paths():
    """Return a list of configuration file paths."""
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

def _get_synced_params(key, value):
    """Return dictionaries for updating the `rcParamsShort`, `rcParamsLong`,
    and `rcParams` properties associted with this key."""
    kw = {} # builtin properties that global setting applies to
    kw_long = {} # custom properties that global setting applies to
    kw_short = {} # short name properties
    if '.' not in key and key not in rcParamsShort:
        key = RC_NODOTS.get(key, key)
    # Skip full name keys
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
        kw_long['gridminor.linewidth'] = gridwidth*ratio

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
    if key in rcParamsShort:
        kw_short[key] = value
    elif key in rcParamsLong:
        kw_long[key] = value
    elif key in rcParams:
        kw[key] = value
    else:
        raise KeyError(f'Invalid key {key!r}.')
    for name in RC_CHILDREN.get(key, ()):
        if name in rcParamsLong:
            kw_long[name] = value
        else:
            kw[name] = value
    return kw_short, kw_long, kw

#-----------------------------------------------------------------------------#
# Main class
#-----------------------------------------------------------------------------#
def _sanitize_key(key):
    """Convert the key to a palatable value."""
    if not isinstance(key, str):
        raise KeyError(f'Invalid key {key!r}. Must be string.')
    if '.' not in key and key not in rcParamsShort:
        key = RC_NODOTS.get(key, key)
    return key.lower()

class rc_configurator(object):
    """
    Magical abstract class for managing matplotlib `rcParams
    <https://matplotlib.org/users/customizing.html>`__ and additional
    ProPlot :ref:`rcParamsLong` and :ref:`rcParamsShort` settings. When
    initialized, this loads defaults settings plus any user overrides in the
    ``~/.proplotrc`` file. See the `~proplot.rctools` documentation for
    details.
    """
    def __contains__(self, key):
        return key in rcParamsShort or key in rcParamsLong or key in rcParams
    def __iter__(self):
        for key in sorted((*rcParamsShort, *rcParamsLong, *rcParams)):
            yield key
    def __repr__(self):
        rcdict = type('rc', (dict,), {})(rcParamsShort)
        string = type(rcParams).__repr__(rcdict)
        indent = ' ' * 4 # indent is rc({
        return string.strip('})') + f'\n{indent}... (rcParams) ...\n{indent}}})'
    def __str__(self): # encapsulate params in temporary class whose name is used by rcParams.__str__
        rcdict = type('rc', (dict,), {})(rcParamsShort)
        string = type(rcParams).__str__(rcdict)
        return string + '\n... (rcParams) ...'

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
        object.__setattr__(self, '_context', [])
        with _benchmark('  use'):
            style.use('default')

        # Update from defaults
        rcParamsLong.clear()
        rcParamsLong.update(defaultParamsLong)
        rcParamsShort.clear()
        rcParamsShort.update(defaultParamsShort)
        for rcdict in (rcParamsShort, rcParamsLong):
            for key,value in rcdict.items():
                _, rc_long, rc = _get_synced_params(key, value)
                rcParamsLong.update(rc_long)
                rcParams.update(rc)

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
        """Apply settings from the most recent context block."""
        *_, kwargs, cache, restore = self._context[-1] # missing arg is previous mode
        def _update(rcdict, newdict):
            for key,value in newdict.items():
                restore[key] = rcdict[key]
                rcdict[key] = cache[key] = value
        for key,value in kwargs.items():
            rc_short, rc_long, rc = _get_synced_params(key, value)
            _update(rcParamsShort, rc_short)
            _update(rcParamsLong, rc_long)
            _update(rcParams, rc)

    def __exit__(self, *args):
        """Restore settings from the most recent context block."""
        *_, restore = self._context[-1]
        for key,value in restore.items():
            self[key] = value
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
        return self[attr]

    def __getitem__(self, key):
        """Return the relevant `rcParams <https://matplotlib.org/users/customizing.html>`__,
        :ref:`rcParamsLong`, and :ref:`rcParamsShort` setting."""
        key = _sanitize_key(key)
        for kw in (rcParamsShort, rcParamsLong, rcParams):
            try:
                return kw[key]
            except KeyError:
                continue
        raise KeyError(f'Invalid property name {key!r}.')

    def __setattr__(self, attr, value):
        """Pass the attribute and value to `~rc_configurator.__setitem__`."""
        self[attr] = value

    def __setitem__(self, key, value):
        """Modify the relevant `rcParams <https://matplotlib.org/users/customizing.html>`__,
        :ref:`rcParamsLong`, and :ref:`rcParamsShort` setting(s)."""
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
            rcdicts = (*caches, rcParamsShort, rcParamsLong) # custom only!
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
            raise KeyError(f'Invalid property name {key!r}.')
        else:
            return None

    def category(self, cat, *, context=False):
        """
        Return a dictionary of settings beginning with the substring
        ``cat + '.'``.

        Parameters
        ----------
        cat : str, optional
            The `rc` settings category.
        context : bool, optional
            If ``True``, then each category setting that is not found in the
            context mode dictionaries is omitted from the output dictionary.
            See `~rc_configurator.context`.
        """
        if cat not in RC_CATEGORIES:
            raise ValueError(f'Invalid rc category {cat!r}. Valid categories are {", ".join(map(repr, RC_CATEGORIES))}.')
        kw = {}
        mode = 0 if not context else None
        for rcdict in (rcParamsLong, rcParams):
            for key in rcdict:
                if not re.search(f'^{cat}[.][^.]+$', key):
                    continue
                value = self._get_item(key, mode)
                if value is None:
                    continue
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
            Dictionary whose values are names of settings. The values
            are replaced with the corresponding property only if
            `~rc_configurator.__getitem__` does not return ``None``. Otherwise,
            that key, value pair is omitted from the output dictionary.
        context : bool, optional
            If ``True``, then each setting that is not found in the
            context mode dictionaries is omitted from the output dictionary.
            See `~rc_configurator.context`.
        """
        kw = {}
        mode = 0 if not context else None
        for key,value in props.items():
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
        Update multiple settings at once.

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

@_timer
def backend_setup(backend=None, fmt=None):
    """
    Set up the matplotlib backend for your iPython workspace.

    Parameters
    ----------
    backend : str, optional
        The backend name. Leave this empty or use ``'auto'`` to revert to the
        ProPlot defaults.
    fmt : str, optional
        The inline backend file format. Valid formats include ``'jpg'``,
        ``'png'``, ``'svg'``, ``'pdf'``, and ``'retina'``. This is ignored
        for non-inline backends.
    """
    # Initialize with default 'inline' settings
    # Reset rc object afterwards
    # TODO: Change nbsetup --> autobackend in add-subplot branch
    fmt = fmt or rcParamsShort['format']
    ipython = get_ipython()
    backend = backend or ('auto' if
        rcParamsShort.get('autobackend', rcParamsShort.get('nbsetup', True))
        else None) or rcParams['backend']
    if ipython is None or backend is None:
        return
    if backend[:2] == 'nb' or backend in ('MacOSX',):
        warnings.warn(f'Using ProPlot with the {backend!r} backend may result in unexpected behavior due to automatic figure resizing. Try using %matplotlib inline or %matplotlib qt, or just import proplot before specifying the backend and one of these will be automatically loaded.')
        backend = 'auto'

    # For notebooks
    try:
        ibackend = ('inline' if backend == 'auto' else backend)
        ipython.magic('matplotlib ' + ibackend)
        rc.reset()
    # For terminals
    # KeyError is evidently a subclass of the UnknownBackup exception
    except KeyError as err:
        if backend != 'auto':
            raise err
        ipython.magic('matplotlib qt') # use any available Qt backend
        rc.reset()

    # Configure inline backend no matter what type of session this is
    # Should be silently ignored for terminal ipython sessions
    ipython.magic("config InlineBackend.figure_formats = ['" + fmt + "']")
    ipython.magic("config InlineBackend.rc = {}") # no notebook-specific overrides
    ipython.magic("config InlineBackend.close_figures = True") # memory issues
    ipython.magic("config InlineBackend.print_figure_kwargs = {'bbox_inches':None}") # use ProPlot tight layout

def autoreload_setup(autoreload=None):
    """
    Set up the
    `autoreload <https://ipython.readthedocs.io/en/stable/config/extensions/autoreload.html>`__
    utility for ipython sessions.

    Parameters
    ----------
    autoreload : float, optional
        The autoreload level. Default is :rc:`autoreload`.
    """
    ipython = get_ipython()
    autoreload = autoreload or rcParamsShort['autoreload']
    if ipython is None or autoreload is None:
        return
    if 'autoreload' not in ipython.magics_manager.magics['line']:
        with IPython.utils.io.capture_output(): # capture annoying message
            ipython.magic("load_ext autoreload")
    ipython.magic("autoreload " + str(autoreload))

def autosave_setup(autosave=None):
    """
    Set up the
    `autosave <https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-matplotlib>`__
    utility for ipython notebook sessions.

    Parameters
    ----------
    autosave : float, optional
        The autosave interval in seconds. Default is :rc:`autosave`.
    """
    ipython = get_ipython()
    autosave = autosave or rcParamsShort['autosave']
    if ipython is None or autosave is None:
        return
    with IPython.utils.io.capture_output(): # capture annoying message
        try:
            ipython.magic("autosave " + str(autosave))
        except IPython.core.error.UsageError:
            pass

# Call setup functions and declare rc object
# WARNING: Must be instantiated after ipython notebook setup! The default
# backend may change some rc settings!
rc = rc_configurator()
"""Instance of `rc_configurator`. This is used to change global settings.
See :ref:`Configuring proplot` for details."""
backend_setup()
autoreload_setup()
autosave_setup()

