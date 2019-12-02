#!/usr/bin/env python3
"""
Utilities for configuring matplotlib and ProPlot global settings.
See :ref:`Configuring proplot` for details.
"""
# TODO: Add 'style' setting that overrides .proplotrc
# Adapted from seaborn; see:
# https://github.com/mwaskom/seaborn/blob/master/seaborn/rcmod.py
from matplotlib import rcParams as rcParams
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
    def get_ipython():
        return None
__all__ = [
    'rc', 'rc_configurator', 'autosave_setup',
    'autoreload_setup', 'backend_setup'
]

# Initialize
rcParamsShort = {}
rcParamsCustom = {}

# "Global" settings and the lower-level settings they change
# NOTE: This whole section, declaring dictionaries and sets, takes 1ms
RC_CHILDREN = {
    'fontname': ('font.family',),
    'cmap': ('image.cmap',),
    'lut': ('image.lut',),
    'alpha': ('axes.alpha',),  # this is a custom setting
    'facecolor': ('axes.facecolor', 'geoaxes.facecolor'),
    # change the 'color' of an axes
    'color': (
        'axes.edgecolor', 'geoaxes.edgecolor', 'axes.labelcolor',
        'tick.labelcolor', 'hatch.color', 'xtick.color', 'ytick.color'),
    # the 'small' fonts
    'small': (
        'font.size', 'tick.labelsize', 'xtick.labelsize', 'ytick.labelsize',
        'axes.labelsize', 'legend.fontsize', 'geogrid.labelsize'),
    # the 'large' fonts
    'large': (
        'abc.size', 'figure.titlesize',
        'axes.titlesize', 'suptitle.size', 'title.size',
        'leftlabel.size', 'toplabel.size',
        'rightlabel.size', 'bottomlabel.size'),
    'linewidth': (
        'axes.linewidth', 'geoaxes.linewidth', 'hatch.linewidth',
        'xtick.major.width', 'ytick.major.width'),
    'margin': ('axes.xmargin', 'axes.ymargin'),
    'grid': ('axes.grid',),
    'gridminor': ('axes.gridminor',),
    'geogrid': ('axes.geogrid',),
    'ticklen': ('xtick.major.size', 'ytick.major.size'),
    'tickdir': ('xtick.direction', 'ytick.direction'),
    'tickpad': (
        'xtick.major.pad', 'xtick.minor.pad',
        'ytick.major.pad', 'ytick.minor.pad'),
}

# Names of the new settings
RC_PARAMNAMES = {*rcParams.keys()}
RC_SHORTNAMES = {
    'abc', 'span', 'share', 'align', 'tight', 'fontname', 'cmap', 'lut',
    'cycle', 'rgbcycle', 'alpha', 'facecolor', 'color', 'small', 'large',
    'linewidth', 'margin', 'grid', 'gridminor', 'geogrid',
    'ticklen', 'tickdir', 'tickpad', 'tickratio', 'ticklenratio', 'gridratio',
    'reso', 'land', 'ocean', 'lakes', 'coast',
    'borders', 'innerborders', 'rivers',
    'nbsetup', 'format', 'autosave', 'autoreload'
}
RC_CUSTOMNAMES = {
    'axes.formatter.zerotrim', 'axes.formatter.timerotation',
    'axes.gridminor', 'axes.geogrid', 'axes.alpha',
    'image.levels', 'image.edgefix',
    'geoaxes.linewidth', 'geoaxes.facecolor', 'geoaxes.edgecolor',
    'land.color', 'ocean.color', 'lakes.color',
    'coast.color', 'coast.linewidth',
    'borders.color', 'borders.linewidth',
    'innerborders.color', 'innerborders.linewidth',
    'rivers.color', 'rivers.linewidth',
    'abc.size', 'abc.weight', 'abc.color',
    'abc.loc', 'abc.style', 'abc.border', 'abc.linewidth',
    'title.loc', 'title.pad', 'title.color', 'title.border',
    'title.linewidth', 'title.weight', 'title.size',
    'suptitle.size', 'suptitle.weight', 'suptitle.color',
    'leftlabel.size', 'leftlabel.weight', 'leftlabel.color',
    'rightlabel.size', 'rightlabel.weight', 'rightlabel.color',
    'toplabel.size', 'toplabel.weight', 'toplabel.color',
    'bottomlabel.size', 'bottomlabel.weight', 'bottomlabel.color',
    'gridminor.alpha', 'gridminor.color',
    'gridminor.linestyle', 'gridminor.linewidth',
    'geogrid.labels', 'geogrid.alpha',
    'geogrid.color', 'geogrid.labelsize',
    'geogrid.linewidth', 'geogrid.linestyle',
    'geogrid.latmax', 'geogrid.lonstep', 'geogrid.latstep',
    'tick.labelweight', 'tick.labelcolor', 'tick.labelsize',
    'subplots.pad', 'subplots.axpad', 'subplots.panelpad',
    'subplots.ylabspace', 'subplots.xlabspace',
    'subplots.innerspace', 'subplots.titlespace',
    'subplots.axwidth', 'subplots.panelwidth', 'subplots.panelspace',
    'colorbar.grid', 'colorbar.frameon', 'colorbar.framealpha',
    'colorbar.loc', 'colorbar.length',
    'colorbar.width', 'colorbar.insetlength', 'colorbar.insetwidth',
    'colorbar.extend', 'colorbar.insetextend', 'colorbar.axespad',
}
# Used by Axes.format, allows user to pass rc settings as keyword args,
# way less verbose. For example, landcolor='b' vs. rc_kw={'land.color':'b'}.
RC_NODOTSNAMES = {  # useful for passing these as kwargs
    name.replace('.', ''): name for names in
    (RC_CUSTOMNAMES, RC_PARAMNAMES, RC_SHORTNAMES)
    for name in names
}
# Categories for returning dict of subcategory properties
RC_CATEGORIES = {
    *(re.sub(r'\.[^.]*$', '', name) for names in
        (RC_CUSTOMNAMES, RC_PARAMNAMES) for name in names),
    *(re.sub(r'\..*$', '', name) for names in
        (RC_CUSTOMNAMES, RC_PARAMNAMES) for name in names)
}


def _convert_units(key, value):
    """Converts certain keys to the units "points". If "key" is passed, tests
    that key against possible keys that accept physical units."""
    # See: https://matplotlib.org/users/customizing.html, all props matching
    # the strings use the units 'points', and special categories are inches!
    # WARNING: Must keep colorbar and subplots units alive, so when user
    # requests em units, values change with respect to font size. The points
    # thing is a conveniene feature so not as important for them.
    if (isinstance(value, str)
            and key.split('.')[0] not in ('colorbar', 'subplots')
            and re.match('^.*(width|space|size|pad|len|small|large)$', key)):
        value = utils.units(value, 'pt')
    return value


def _set_cycler(name):
    """Sets the default color cycler."""
    # Draw from dictionary
    try:
        colors = mcm.cmap_d[name].colors
    except (KeyError, AttributeError):
        cycles = sorted(name for name, cmap in mcm.cmap_d.items()
                        if isinstance(cmap, mcolors.ListedColormap))
        raise ValueError(
            f'Invalid cycle name {name!r}. Options are: {", ".join(cycles)}')
    # Apply color name definitions
    if rcParamsShort['rgbcycle'] and name.lower() == 'colorblind':
        regcolors = colors + [(0.1, 0.1, 0.1)]
    elif mcolors.to_rgb('r') != (1.0, 0.0, 0.0):  # reset
        regcolors = [
            (0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0),
            (0.75, 0.75, 0.0), (0.75, 0.75, 0.0), (0.0, 0.75, 0.75),
            (0.0, 0.0, 0.0)]
    else:
        regcolors = []  # no reset necessary
    for code, color in zip('brgmyck', regcolors):
        rgb = mcolors.to_rgb(color)
        mcolors.colorConverter.colors[code] = rgb
        mcolors.colorConverter.cache[code] = rgb
    # Pass to cycle constructor
    rcParams['patch.facecolor'] = colors[0]
    rcParams['axes.prop_cycle'] = cycler.cycler('color', colors)


def _get_config_paths():
    """Returns configuration file paths."""
    # Get paths
    idir = os.getcwd()
    paths = []
    while idir:  # not empty string
        ipath = os.path.join(idir, '.proplotrc')
        if os.path.exists(ipath):
            paths.append(ipath)
        ndir, _ = os.path.split(idir)
        if ndir == idir:
            break
        idir = ndir
    paths = paths[::-1]  # sort from decreasing to increasing importantce
    # Home configuration
    ipath = os.path.join(os.path.expanduser('~'), '.proplotrc')
    if os.path.exists(ipath) and ipath not in paths:
        paths.insert(0, ipath)
    # Global configuration
    ipath = os.path.join(os.path.dirname(__file__), '.proplotrc')
    if ipath in paths:
        paths.remove(ipath)
    paths.insert(0, ipath)
    return paths


def _get_synced_params(key=None, value=None):
    """Returns dictionaries for updating "child" properties in
    `rcParams` and `rcParamsCustom` with global property."""
    kw = {}  # builtin properties that global setting applies to
    kw_custom = {}  # custom properties that global setting applies to
    if key is not None and value is not None:
        items = [(key, value)]
    else:
        items = rcParamsShort.items()
    for key, value in items:
        # Tick length/major-minor tick length ratio
        if key in ('ticklen', 'ticklenratio'):
            if key == 'ticklen':
                ticklen = _convert_units(key, value)
                ratio = rcParamsShort['ticklenratio']
            else:
                ticklen = rcParamsShort['ticklen']
                ratio = value
            kw['xtick.minor.size'] = ticklen * ratio
            kw['ytick.minor.size'] = ticklen * ratio
        # Spine width/major-minor tick width ratio
        elif key in ('linewidth', 'tickratio'):
            if key == 'linewidth':
                tickwidth = _convert_units(key, value)
                ratio = rcParamsShort['tickratio']
            else:
                tickwidth = rcParamsShort['linewidth']
                ratio = value
            kw['xtick.minor.width'] = tickwidth * ratio
            kw['ytick.minor.width'] = tickwidth * ratio
        # Grid line
        elif key in ('grid.linewidth', 'gridratio'):
            if key == 'grid.linewidth':
                gridwidth = _convert_units(key, value)
                ratio = rcParamsShort['gridratio']
            else:
                gridwidth = rcParams['grid.linewidth']
                ratio = value
            kw_custom['gridminor.linewidth'] = gridwidth * ratio
        # Now update linked settings
        val = None
        for name in RC_CHILDREN.get(key, ()):
            val = _convert_units(key, value)
            if name in rcParamsCustom:
                kw_custom[name] = val
            else:
                kw[name] = val
        if key == 'linewidth' and val == 0:
            ikw, ikw_custom = _get_synced_params('ticklen', 0)
            kw.update(ikw)
            kw_custom.update(ikw_custom)
    return kw, kw_custom


class rc_configurator(object):
    _public_api = (
        'get', 'fill', 'category', 'reset', 'context', 'update'
    )  # getattr and setattr will not look for these items

    def __str__(self):
        return type(rcParams).__str__(rcParamsShort)  # just show globals

    def __repr__(self):
        return type(rcParams).__repr__(rcParamsShort)

    def __contains__(self, key):
        return (key in RC_SHORTNAMES or key in RC_CUSTOMNAMES or key in
                RC_PARAMNAMES or key in RC_NODOTSNAMES)  # biggest lists last

    @_counter  # about 0.05s
    def __init__(self):
        """Magical abstract class for managing matplotlib `rcParams \
<https://matplotlib.org/users/customizing.html>`__
        settings, ProPlot :ref:`rcParamsCustom` settings, and
        :ref:`rcParamsShort` "global" settings. This starts with the
        default settings plus user ``.proplotrc`` overrides.
        See :ref:`Configuring proplot` for details."""
        # Set the default style. Note that after first figure made, backend
        # is 'sticky', never changes! See:
        # https://stackoverflow.com/a/48322150/4970632
        plt.style.use('default')

        # Load the defaults from file
        for i, file in enumerate(_get_config_paths()):
            # Load
            if not os.path.exists(file):
                continue
            with open(file) as f:
                try:
                    data = yaml.safe_load(f)
                except yaml.YAMLError as err:
                    print('{file!r} has invalid YAML syntax.')
                    raise err
            # Special duplicate keys
            if data is None:
                continue
            # Add keys to dictionaries
            gkeys, ckeys = {*()}, {*()}
            for key, value in data.items():
                if key in RC_SHORTNAMES:
                    rcParamsShort[key] = value
                    if i == 0:
                        gkeys.add(key)
                elif key in RC_CUSTOMNAMES:
                    value = _convert_units(key, value)
                    rcParamsCustom[key] = value
                    if i == 0:
                        ckeys.add(key)
                elif key in RC_PARAMNAMES:
                    value = _convert_units(key, value)
                    rcParams[key] = value
                else:
                    raise RuntimeError(f'{file!r} has invalid key {key!r}.')
            # Make sure we did not miss anything
            if i == 0:
                if gkeys != RC_SHORTNAMES:
                    raise RuntimeError(
                        f'{file!r} has incomplete or invalid global keys '
                        f'{RC_SHORTNAMES - gkeys}.')
                if ckeys != RC_CUSTOMNAMES:
                    raise RuntimeError(
                        f'{file!r} has incomplete or invalid custom keys '
                        f'{RC_CUSTOMNAMES - ckeys}.')

        # Apply *global settings* to children settings
        rcParams['axes.titlepad'] = rcParamsCustom['title.pad']
        _set_cycler(rcParamsShort['cycle'])
        rc, rc_new = _get_synced_params()
        for key, value in rc.items():
            rcParams[key] = value
        for key, value in rc_new.items():
            rcParamsCustom[key] = value

        # Caching stuff
        self._init = True
        self._getitem_mode = 0
        self._context = {}
        self._cache = {}
        self._cache_orig = {}
        self._cache_restore = {}

    def __getitem__(self, key):
        """Get `rcParams <https://matplotlib.org/users/customizing.html>`__,
        :ref:`rcParamsCustom`, and :ref:`rcParamsShort` settings. If we are in
        a `~rc_configurator.context` block, may return ``None`` if the setting
        is not cached (i.e. if it was not changed by the user)."""
        # Can get a whole bunch of different things
        # Get full dictionary e.g. for rc[None]
        if not key:
            return {**rcParams, **rcParamsCustom}

        # Standardize
        # NOTE: If key is invalid, raise error down the line.
        if '.' not in key and key not in rcParamsShort:
            key = RC_NODOTSNAMES.get(key, key)

        # Allow for special time-saving modes where we *ignore rcParams*
        # or even *ignore rcParamsCustom*.
        mode = self._getitem_mode
        if mode == 0:
            kws = (self._cache, rcParamsShort, rcParamsCustom, rcParams)
        elif mode == 1:
            kws = (self._cache, rcParamsShort, rcParamsCustom)  # custom only!
        elif mode == 2:
            kws = (self._cache,)  # changed only!
        else:
            raise KeyError(f'Invalid caching mode {mode!r}.')

        # Get individual property. Will successively index a few different
        # dicts. Try to return the value
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
        """Set `rcParams <https://matplotlib.org/users/customizing.html>`__,
        :ref:`rcParamsCustom`, and :ref:`rcParamsShort` settings."""
        # Check whether we are in context block
        # NOTE: Do not add key to cache until we are sure it is a valid key
        cache = self._cache
        context = bool(self._context)  # test if context dict is non-empty
        if context:
            restore = self._cache_restore

        # Standardize
        # NOTE: If key is invalid, raise error down the line.
        if '.' not in key and key not in rcParamsShort:
            key = RC_NODOTSNAMES.get(key, key)

        # Special keys
        if key == 'title.pad':
            key = 'axes.titlepad'
        if key == 'rgbcycle':  # if must re-apply cycler afterward
            cache[key] = value
            rcParamsShort[key] = value
            key, value = 'cycle', rcParamsShort['cycle']

        # Set the default cycler
        if key == 'cycle':
            cache[key] = value
            if context:
                restore[key] = rcParamsShort[key]
                restore['axes.prop_cycle'] = rcParams['axes.prop_cycle']
                restore['patch.facecolor'] = rcParams['patch.facecolor']
            _set_cycler(value)

        # Gridline toggling, complicated because of the clunky way this is
        # implemented in matplotlib. There should be a gridminor setting!
        elif key in ('grid', 'gridminor'):
            cache[key] = value
            ovalue = rcParams['axes.grid']
            owhich = rcParams['axes.grid.which']
            if context:
                restore[key] = rcParamsShort[key]
                restore['axes.grid'] = ovalue
                restore['axes.grid.which'] = owhich
            # Instruction is to turn off gridlines
            if not value:
                # Gridlines are already off, or they are on for the particular
                # ones that we want to turn off. Instruct to turn both off.
                if not ovalue or (key == 'grid' and owhich == 'major') or (
                        key == 'gridminor' and owhich == 'minor'):
                    which = 'both'  # disable both sides
                # Gridlines are currently on for major and minor ticks, so we
                # instruct to turn on gridlines for the one we *don't* want off
                elif owhich == 'both':  # and ovalue is True
                    value = True
                    # if gridminor=False, enable major, and vice versa
                    which = 'major' if key == 'gridminor' else 'minor'
                # Gridlines are on for the ones that we *didn't* instruct to
                # turn off, and off for the ones we do want to turn off. This
                # just re-asserts the ones that are already on.
                else:
                    value = True
                    which = owhich
            # Instruction is to turn on gridlines
            else:
                # Gridlines are already both on, or they are off only for the
                # ones that we want to turn on. Turn on gridlines for both.
                if owhich == 'both' or (key == 'grid' and owhich == 'minor') \
                        or (key == 'gridminor' and owhich == 'major'):
                    which = 'both'
                # Gridlines are off for both, or off for the ones that we
                # don't want to turn on. We can just turn on these ones.
                else:
                    which = owhich
            cache.update({'axes.grid': value, 'axes.grid.which': which})
            rcParams.update({'axes.grid': value, 'axes.grid.which': which})

        # Ordinary settings
        elif key in rcParamsShort:
            # Update global setting
            cache[key] = value
            if context:
                restore[key] = rcParamsShort[key]
            rcParamsShort[key] = value
            # Update children of setting
            rc, rc_new = _get_synced_params(key, value)
            cache.update(rc)
            cache.update(rc_new)
            if context:
                restore.update({key: rcParams[key] for key in rc})
                restore.update({key: rcParamsCustom[key] for key in rc_new})
            rcParams.update(rc)
            rcParamsCustom.update(rc_new)
        # Update normal settings
        elif key in RC_CUSTOMNAMES:
            value = _convert_units(key, value)
            cache[key] = value
            if context:
                restore[key] = rcParamsCustom[key]
            rcParamsCustom[key] = value
        elif key in RC_PARAMNAMES:
            value = _convert_units(key, value)
            cache[key] = value
            if context:
                restore[key] = rcParams[key]
            rcParams[key] = value  # rcParams dict has key validation
        else:
            raise KeyError(f'Invalid key {key!r}.')
        self._init = False  # we are no longer in initial state

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
        self._cache_restore = {}  # shouldn't be necessary but just in case
        self._cache = {}
        for key, value in self._context.items():
            self[key] = value  # applies linked and individual settings

    def __exit__(self, *args):
        """Restore configurator cache to initial state."""
        self._context = {}
        self._getitem_mode = 0
        for key, value in self._cache_restore.items():
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
            block when settings are requested with e.g. :rcraw:`setting`. If
            you are using `~rc_configurator.context` manually, the `mode` is
            automatically set to ``0`` -- other input is ignored. Internally,
            ProPlot uses all of the three available modes.

            0. All settings (`rcParams \
<https://matplotlib.org/users/customizing.html>`__,
               :ref:`rcParamsCustom`, and :ref:`rcParamsShort`) are returned,
               whether or not `~rc_configurator.context` has changed them.
            1. Unchanged `rcParams \
<https://matplotlib.org/users/customizing.html>`__
               return ``None``. :ref:`rcParamsCustom` and :ref:`rcParamsShort`
               are returned whether or not `~rc_configurator.context` has
               changed them.  This is used in the ``__init__`` call to
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
        self._context = kwargs  # could be empty
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
            If ``False``, the `~rc_configurator.__getitem__` mode is
            temporarily set to ``0`` (see `~rc_configurator.context`).
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
        Returns a dictionary filled with `rc` settings, used internally to
        build dictionaries for updating `~matplotlib.artist.Artist` instances.

        Parameters
        ----------
        props : dict-like
            Dictionary whose values are names of `rc` settings. The values
            are replaced with the corresponding property only if
            `~rc_configurator.__getitem__` does not return ``None``. Otherwise,
            that key, value pair is omitted from the output dictionary.
        cache : bool, optional
            If ``False``, the `~rc_configurator.__getitem__` mode is
            temporarily set to ``0`` (see `~rc_configurator.context`).
            Otherwise, if an `rc` lookup returns ``None``, the setting is
            omitted from the output dictionary.
        """
        if not cache:
            orig = self._getitem_mode
            self._getitem_mode = 0
        props_out = {}
        for key, value in props.items():
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
            If ``False``, the `~rc_configurator.__getitem__` mode is
            temporarily set to ``0`` (see `~rc_configurator.context`).
        """
        # Check
        if cat not in RC_CATEGORIES:
            raise ValueError(
                f'RC category {cat!r} does not exist. Valid categories are '
                f', '.join(map(repr, RC_CATEGORIES)) + '.')
        if not cache:
            mode = 0
        else:
            mode = self._getitem_mode

        # Allow for special time-saving modes where we *ignore rcParams*
        # or even *ignore rcParamsCustom*.
        if mode == 0:
            kws = (self._cache, rcParamsShort, rcParamsCustom, rcParams)
        elif mode == 1:
            kws = (self._cache, rcParamsShort, rcParamsCustom)
        elif mode == 2:
            kws = (self._cache, rcParamsShort)
        else:
            raise KeyError(f'Invalid caching mode {mode}.')

        # Return params dictionary
        params = {}
        for kw in kws:
            for category, value in kw.items():
                if re.search(rf'^{cat}\.', category):
                    subcategory = re.sub(rf'^{cat}\.', '', category)
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
            a "category" string name. If a category name is passed, all
            settings in the dictionary (if it was passed) and all keyword arg
            names (if they were passed) are prepended with the string
            ``cat + '.'``. For example,
            ``rc.update('axes', labelsize=20, titlesize=20)``
            changes the ``axes.labelsize`` and ``axes.titlesize`` properties.
        **kwargs
            `rc` settings passed as keyword args.
        """
        # Parse args
        kw = {}
        prefix = ''
        if len(args) > 2:
            raise ValueError(
                'Accepts 1-2 positional arguments. Use plot.rc.update(kw) '
                'to update a bunch of names, or plot.rc.update(category, kw) '
                'to update subcategories belonging to single category '
                'e.g. axes. Keyword args are added to the kw dict.')
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

    def reset(self):
        """Restores settings to the initial state -- ProPlot defaults, plus
        any user overrides in the ``~/.proplotrc`` file."""
        if not self._init:  # save resources if rc is unchanged!
            return self.__init__()


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
    backend = backend or (
        'auto' if rcParamsShort.get(
            'autobackend', rcParamsShort.get('nbsetup', True)
        ) else None) or rcParams['backend']
    if ipython is None or backend is None:
        return
    if backend[:2] == 'nb' or backend in ('MacOSX',):
        warnings.warn(
            f'Using ProPlot with the {backend!r} backend may result in '
            'unexpected behavior due to automatic figure resizing. '
            'Try using %matplotlib inline or %matplotlib qt, or just '
            'import proplot before specifying the backend and one of these '
            'will be automatically loaded.')
        backend = 'auto'

    # For notebooks
    rc._init = False
    ibackend = ('inline' if backend == 'auto' else backend)
    try:
        ipython.magic('matplotlib ' + ibackend)
        rc.reset()
    # For terminals (UnknownBackend is subclass of KeyError)
    except KeyError:
        if backend != 'auto':
            warnings.warn(f'Failed to import {backend!r} backend.')
        try:
            ipython.magic('matplotlib qt')  # use any available Qt backend
            rc.reset()
        except KeyError:  # should be impossible; pyplot needs Qt!
            warnings.warn('Failed to import \'qt\' backend.')

    # Configure inline backend no matter what type of session this is
    # Should be silently ignored for terminal ipython sessions
    ipython.magic("config InlineBackend.figure_formats = ['" + fmt + "']")
    ipython.magic('config InlineBackend.rc = {}')  # no notebook overrides
    ipython.magic('config InlineBackend.close_figures = True')  # memory issues
    # Use ProPlot tight layout
    ipython.magic(
        "config InlineBackend.print_figure_kwargs = {'bbox_inches':None}")


def autoreload_setup(autoreload=None):
    """
    Set up the `autoreload \
<https://ipython.readthedocs.io/en/stable/config/extensions/autoreload.html>`__
    utility for ipython sessions.

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


def autosave_setup(autosave=None):
    """
    Set up the `autosave \
<https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-matplotlib>`__
    utility for ipython notebook sessions.

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


#: Instance of `rc_configurator`. This is used to change global settings.
#: See :ref:`Configuring proplot` for details.
rc = rc_configurator()

# Call setup functions
backend_setup()
autoreload_setup()
autosave_setup()
