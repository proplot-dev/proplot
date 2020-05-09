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
import matplotlib as mpl
import matplotlib.font_manager as mfonts
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
import cycler
from collections import namedtuple
from . import colors as pcolors
from .utils import units, to_xyz
from .internals import ic  # noqa: F401
from .internals import defaults, docstring, timers, warnings, _not_none
try:
    from IPython import get_ipython
except ImportError:
    def get_ipython():
        return

__all__ = [
    'rc', 'rc_configurator',
    'register_cmaps', 'register_cycles', 'register_colors', 'register_fonts',
    'inline_backend_fmt',  # deprecated
]

# Disable mathtext "missing glyph" warnings
import matplotlib.mathtext  # noqa
import logging
logger = logging.getLogger('matplotlib.mathtext')
logger.setLevel(logging.ERROR)  # suppress warnings!

# Dictionaries used to track custom proplot settings
_Context = namedtuple('Context', ('mode', 'kwargs', 'rc_new', 'rc_old'))
rc_params = rcParams = mpl.rcParams  # PEP8 4 lyfe
rc_quick = rcParamsShort = defaults._rc_quick_default.copy()
rc_added = rcParamsLong = defaults._rc_added_default.copy()
defaultParams = defaults._rc_matplotlib_default
defaultParamsLong = defaults._rc_added_default
defaultParamsShort = default._rc_quick_defautl

# Misc constants
# TODO: Use explicit validators for specific settings like matplotlib.
REGEX_POINTS = re.compile(
    r'\A(?!colorbar|subplots|pdf|ps).*(width|space|size|pad|len)\Z'
)

ALWAYS_ADD = (
    *(  # common fancy names or natural names
        'charcoal', 'tomato', 'burgundy', 'maroon', 'burgundy', 'lavendar',
        'taupe', 'ocre', 'sand', 'stone', 'earth', 'sand brown', 'sienna',
        'terracotta', 'moss', 'crimson', 'mauve', 'rose', 'teal', 'forest',
        'grass', 'sage', 'pine', 'vermillion', 'russet', 'cerise', 'avocado',
        'wine', 'brick', 'umber', 'mahogany', 'puce', 'grape', 'blurple',
        'cranberry', 'sand', 'aqua', 'jade', 'coral', 'olive', 'magenta',
        'turquoise', 'sea blue', 'royal blue', 'slate blue', 'slate grey',
        'baby blue', 'salmon', 'beige', 'peach', 'mustard', 'lime', 'indigo',
        'cornflower', 'marine', 'cloudy blue', 'tangerine', 'scarlet', 'navy',
        'cool grey', 'warm grey', 'chocolate', 'raspberry', 'denim',
        'gunmetal', 'midnight', 'chartreuse', 'ivory', 'khaki', 'plum',
        'silver', 'tan', 'wheat', 'buff', 'bisque', 'cerulean',
    ),
    *(  # common combos
        'red orange', 'yellow orange', 'yellow green',
        'blue green', 'blue violet', 'red violet',
    ),
    *(  # common names
        prefix + color
        for color in (
            'red', 'orange', 'yellow', 'green', 'blue', 'indigo', 'violet',
            'brown', 'grey'
        )
        for prefix in (
            '', 'light ', 'dark ', 'medium ', 'pale ',
        )
    )
)
ALWAYS_REMOVE = (  # filter these out, let's try to be professional here...
    'shit', 'poop', 'poo', 'pee', 'piss', 'puke', 'vomit', 'snot',
    'booger', 'bile', 'diarrhea',
)
TRANSLATE_COLORS = (  # prevent registering similar-sounding names
    ('/', ' '),
    ("'s", ''),
    ('forrest', 'forest'),  # typo?
    ('reddish', 'red'),  # remove 'ish'
    ('purplish', 'purple'),
    ('bluish', 'blue'),
    ('ish ', ' '),
    ('grey', 'gray'),
    ('pinky', 'pink'),
    ('greeny', 'green'),
    ('bluey', 'blue'),
    ('purply', 'purple'),
    ('purpley', 'purple'),
    ('yellowy', 'yellow'),
    ('robin egg', 'robins egg'),
    ('egg blue', 'egg'),
    ('bluegray', 'blue gray'),
    ('grayblue', 'gray blue'),
    ('lightblue', 'light blue'),
    ('yellowgreen', 'yellow green'),
    ('yelloworange', 'yellow orange'),
)

OPEN_COLORS = {}  # populated during register_colors
XKCD_COLORS = {}  # populated during register_colors
BASE_COLORS = {
    **mcolors.BASE_COLORS,  # shorthand names like 'r', 'g', etc.
    'blue': (0, 0, 1),
    'green': (0, 0.5, 0),
    'red': (1, 0, 0),
    'cyan': (0, 0.75, 0.75),
    'magenta': (0.75, 0, 0.75),
    'yellow': (0.75, 0.75, 0),
    'black': (0, 0, 0),
    'white': (1, 1, 1),
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


def _get_data_paths(subfolder, user=True, default=True, reverse=False):
    """
    Return data folder paths in reverse order of precedence.
    """
    # When loading colormaps, cycles, and colors, files in the latter
    # directories overwrite files in the former directories. When loading
    # fonts, the resulting paths need to be *reversed*.
    paths = []
    if user:
        paths.append(os.path.join(os.path.dirname(__file__), subfolder))
    if default:
        paths.append(os.path.join(os.path.expanduser('~'), '.proplot', subfolder))
    if reverse:
        paths = paths[::-1]
    return paths


def _iter_data_paths(subfolder, **kwargs):
    """
    Iterate over all files in the data paths. Also yield an index indicating
    whether these are default ProPlot files or user files.
    """
    for i, path in enumerate(_get_data_paths(subfolder, **kwargs)):
        for dirname, dirnames, filenames in os.walk(path):
            for filename in filenames:
                if filename[0] == '.':  # UNIX-style hidden files
                    continue
                yield i, dirname, filename


def _write_defaults(filename, comment=True):
    """
    Save a file to the specified path containing the default `rc` settings.

    Parameters
    ----------
    filename : str
        The path.
    comment : bool, optional
        Whether to "comment out" each setting.
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
                parent = rc_parents[key]
            except KeyError:
                raise RuntimeError(
                    f'rcParamsLong param {key!r} has default value of None '
                    'but has no rcParmsShort parent!'
                )
            if parent in defaultParamsShort:
                value = defaultParamsShort[parent]
            elif parent in defaultParams:
                value = defaultParams[parent]
            else:
                value = rcParams[parent]
            defaultParamsLong_filled[key] = value

    with open(filename, 'w') as f:
        f.write(f"""
#---------------------------------------------------------------------
# Use this file to change the default proplot and matplotlib settings
# The syntax is mostly the same as for matplotlibrc files
# For descriptions of each setting see:
# https://proplot.readthedocs.io/en/latest/configuration.html
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
    settings. This loads the default ProPlot settings and the
    user ``.proplotrc`` overrides. See :ref:`Configuring proplot` for details.
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
        # Remove context objects
        object.__setattr__(self, '_context', [])

        # Set default style
        # NOTE: Previously, style.use would trigger first pyplot import because
        # rcParams.__getitem__['backend'] imports pyplot.switch_backend() so it
        # can determine the default backend.
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
        """
        Apply settings from the most recent context block.
        """
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

    def __exit__(self, *args):  # noqa: U100
        """
        Restore settings from the most recent context block.
        """
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

    def __delitem__(self, item):  # noqa: 100
        """
        Raise an error. This enforces pseudo-immutability.
        """
        raise RuntimeError('rc settings cannot be deleted.')

    def __delattr__(self, item):  # noqa: 100
        """
        Raise an error. This enforces pseudo-immutability.
        """
        raise RuntimeError('rc settings cannot be deleted.')

    def __getattr__(self, attr):
        """
        Pass the attribute to `~rc_configurator.__getitem__` and return
        the result.
        """
        if attr[:1] == '_':
            return super().__getattr__(attr)
        else:
            return self[attr]

    def __getitem__(self, key):
        """
        Return an `rcParams <https://matplotlib.org/users/customizing.html>`__,
        :ref:`rcParamsLong`, or :ref:`rcParamsShort` setting.
        """
        key = _sanitize_key(key)
        for kw in (rcParamsShort, rcParamsLong, rcParams):
            try:
                return kw[key]
            except KeyError:
                continue
        raise KeyError(f'Invalid setting name {key!r}.')

    def __setattr__(self, attr, value):
        """
        Pass the attribute and value to `~rc_configurator.__setitem__`.
        """
        self[attr] = value

    def __setitem__(self, key, value):
        """
        Modify an `rcParams <https://matplotlibcorg/users/customizing.html>`__,
        :ref:`rcParamsLong`, and :ref:`rcParamsShort` setting(s).
        """
        rc_short, rc_long, rc = _get_synced_params(key, value)
        rcParamsShort.update(rc_short)
        rcParamsLong.update(rc_long)
        rcParams.update(rc)

    def _get_item(self, key, mode=None):
        """
        As with `~rc_configurator.__getitem__` but the search is limited
        based on the context mode and ``None`` is returned if the key is not
        found in the dictionaries.
        """
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


def config_inline_backend(fmt=None):
    """
    Set up the `ipython inline backend \
<https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-matplotlib>`__
    format and ensure that inline figures always look the same as saved
    figures. This runs the following ipython magic commands:

    .. code-block:: ipython

        %config InlineBackend.figure_formats = rc['inlinefmt']
        %config InlineBackend.rc = {}  # never override rc settings
        %config InlineBackend.close_figures = True  \
# cells start with no active figures
        %config InlineBackend.print_figure_kwargs = {'bbox_inches': None}  \
# never override rc settings

    When the inline backend is inactive or unavailable, this has no effect.
    This function is called when you modify the :rcraw:`inlinefmt` property.

    Parameters
    ----------
    fmt : str or list of str, optional
        The inline backend file format(s). Default is :rc:`inlinefmt`.
        Valid formats include ``'jpg'``, ``'png'``, ``'svg'``, ``'pdf'``,
        and ``'retina'``.
    """
    # Note if inline backend is unavailable this will fail silently
    ipython = get_ipython()
    if ipython is None:
        return
    fmt = _not_none(fmt, rc_quick['inlinefmt'])
    if isinstance(fmt, str):
        fmt = [fmt]
    elif np.iterable(fmt):
        fmt = list(fmt)
    else:
        raise ValueError(
            f'Invalid inline backend format {fmt!r}. Must be string or list thereof.'
        )
    ipython.magic('config InlineBackend.figure_formats = ' + repr(fmt))
    ipython.magic('config InlineBackend.rc = {}')
    ipython.magic('config InlineBackend.close_figures = True')
    ipython.magic("config InlineBackend.print_figure_kwargs = {'bbox_inches': None}")
    ipython.magic(  # use ProPlot tight layout instead
        'config InlineBackend.print_figure_kwargs = {"bbox_inches": None}'
    )


@docstring.add_snippets
def register_cmaps(user=True, default=False):
    """
    Register colormaps packaged with ProPlot or saved to the
    ``~/.proplot/cmaps`` folder. This is called on import.
    Colormaps are registered according to their filenames -- for example,
    ``name.xyz`` will be registered as ``'name'``.

    %(register.ext_table)s

    To visualize the registered colormaps, use `~proplot.show.show_cmaps`.

    Parameters
    ----------
    %(register_cmaps.params)s
    """
    for i, dirname, filename in _iter_data_paths('cmaps', user=user, default=default):
        path = os.path.join(dirname, filename)
        cmap = pcolors.LinearSegmentedColormap.from_file(path, warn_on_failure=True)
        if not cmap:
            continue
        if i == 0 and cmap.name.lower() in (
            'phase', 'graycycle', 'romao', 'broco', 'corko', 'viko',
        ):
            cmap.set_cyclic(True)
        pcolors._cmapdict[cmap.name] = cmap


@docstring.add_snippets
def register_cycles(user=True, default=False):
    """
    Register color cycles packaged with ProPlot or saved to the
    ``~/.proplot/cycles`` folder. This is called on import. Color cycles
    are registered according to their filenames -- for example, ``name.hex``
    will be registered as ``'name'``.

    %(register.ext_table)s

    To visualize the registered color cycles, use `~proplot.show.show_cycles`.

    Parameters
    ----------
    %(register_cycles.params)s
    """
    for _, dirname, filename in _iter_data_paths('cycles', user=user, default=default):
        path = os.path.join(dirname, filename)
        cmap = pcolors.ListedColormap.from_file(path, warn_on_failure=True)
        if not cmap:
            continue
        pcolors._cmapdict[cmap.name] = cmap


@docstring.add_snippets
def register_colors(user=True, default=False, space='hcl', margin=0.10):
    """
    Register the `open-color <https://yeun.github.io/open-color/>`_ colors,
    XKCD `color survey <https://xkcd.com/color/rgb/>`_ colors, and colors
    saved to the ``~/.proplot/colors`` folder. This is called on import.
    The color survey colors are filtered to a subset that is "perceptually
    distinct" in the HCL colorspace. The user color names are loaded from
    ``.txt`` files saved in ``~/.proplot/colors``. Each file should contain
    one line per color in the format ``name : hex``. Whitespace is ignored.

    To visualize the registered colors, use `~proplot.show.show_colors`.

    Parameters
    ----------
    %(register_colors.params)s
    space : {'hcl', 'hsl', 'hpl'}, optional
        The colorspace used to detect "perceptually distinct" colors.
    margin : float, optional
        The margin by which a color's normalized hue, saturation, and
        luminance channel values must differ from the normalized channel
        values of the other colors to be deemed "perceptually distinct."
    """
    # Reset native colors dictionary
    mcolors.colorConverter.colors.clear()  # clean out!
    mcolors.colorConverter.cache.clear()  # clean out!

    # Add in base colors and CSS4 colors so user has no surprises
    for name, dict_ in (('base', BASE_COLORS), ('css', mcolors.CSS4_COLORS)):
        mcolors.colorConverter.colors.update(dict_)

    # Load colors from file and get their HCL values
    # NOTE: Colors that come *later* overwrite colors that come earlier.
    hex = re.compile(rf'\A{pcolors.HEX_PATTERN}\Z')  # match each string
    for i, dirname, filename in _iter_data_paths('colors', user=user, default=default):
        path = os.path.join(dirname, filename)
        cat, ext = os.path.splitext(filename)
        if ext != '.txt':
            raise ValueError(
                f'Unknown color data file extension ({path!r}). '
                'All files in this folder should have extension .txt.'
            )

        # Read data
        loaded = {}
        with open(path, 'r') as fh:
            for cnt, line in enumerate(fh):
                # Load colors from file
                stripped = line.strip()
                if not stripped or stripped[0] == '#':
                    continue
                pair = tuple(
                    item.strip().lower() for item in line.split(':')
                )
                if len(pair) != 2 or not hex.match(pair[1]):
                    warnings._warn_proplot(
                        f'Illegal line #{cnt + 1} in file {path!r}:\n'
                        f'{line!r}\n'
                        f'Lines must be formatted as "name: hexcolor".'
                    )
                    continue
                # Never overwrite "base" colors with xkcd colors.
                # Only overwrite with user colors.
                name, color = pair
                if i == 0 and name in BASE_COLORS:
                    continue
                loaded[name] = color

        # Add every user color and every opencolor color and ensure XKCD
        # colors are "perceptually distinct".
        if i == 1:
            mcolors.colorConverter.colors.update(loaded)
        elif cat == 'opencolor':
            mcolors.colorConverter.colors.update(loaded)
            OPEN_COLORS.update(loaded)
        elif cat == 'xkcd':
            # Always add these colors, but make sure not to add other
            # colors too close to them.
            hcls = []
            filtered = []
            for name in ALWAYS_ADD:
                color = loaded.pop(name, None)
                if color is None:
                    continue
                if 'grey' in name:
                    name = name.replace('grey', 'gray')
                hcls.append(to_xyz(color, space=space))
                filtered.append((name, color))
                mcolors.colorConverter.colors[name] = color
                XKCD_COLORS[name] = color

            # Get locations of "perceptually distinct" colors
            # WARNING: Unique axis argument requires numpy version >=1.13
            for name, color in loaded.items():
                for string, replace in TRANSLATE_COLORS:
                    if string in name:
                        name = name.replace(string, replace)
                if any(string in name for string in ALWAYS_REMOVE):
                    continue  # remove "unpofessional" names
                hcls.append(to_xyz(color, space=space))
                filtered.append((name, color))  # category name pair
            hcls = np.asarray(hcls)
            if not hcls.size:
                continue
            hcls = hcls / np.array([360, 100, 100])
            hcls = np.round(hcls / margin).astype(np.int64)
            _, idxs = np.unique(hcls, return_index=True, axis=0)

            # Register "distinct" colors
            for idx in idxs:
                name, color = filtered[idx]
                mcolors.colorConverter.colors[name] = color
                XKCD_COLORS[name] = color
        else:
            raise ValueError(f'Unknown proplot color database {path!r}.')


def register_fonts():
    """
    Add fonts packaged with ProPlot or saved to the ``~/.proplot/fonts``
    folder, if they are not already added. Detects ``.ttf`` and ``.otf`` files
    -- see `this link \
<https://gree2.github.io/python/2015/04/27/python-change-matplotlib-font-on-mac>`__
    for a guide on converting various other font file types to ``.ttf`` and
    ``.otf`` for use with matplotlib.

    To visualize the registered fonts, use `~proplot.show.show_fonts`.
    """
    # Find proplot fonts
    # WARNING: If you include a font file with an unrecognized style,
    # matplotlib may use that font instead of the 'normal' one! Valid styles:
    # 'ultralight', 'light', 'normal', 'regular', 'book', 'medium', 'roman',
    # 'semibold', 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black'
    # https://matplotlib.org/api/font_manager_api.html
    # For macOS the only fonts with 'Thin' in one of the .ttf file names
    # are Helvetica Neue and .SF NS Display Condensed. Never try to use these!
    paths_proplot = _get_data_paths('fonts', reverse=True)
    fnames_proplot = set(mfonts.findSystemFonts(paths_proplot))

    # Detect user-input ttc fonts and issue warning
    fnames_proplot_ttc = {
        file for file in fnames_proplot if os.path.splitext(file)[1] == '.ttc'
    }
    if fnames_proplot_ttc:
        warnings._warn_proplot(
            'Ignoring the following .ttc fonts because they cannot be '
            'saved into PDF or EPS files (see matplotlib issue #3135): '
            + ', '.join(map(repr, sorted(fnames_proplot_ttc)))
            + '. Please consider expanding them into separate .ttf files.'
        )

    # Rebuild font cache only if necessary! Can be >50% of total import time!
    fnames_all = {font.fname for font in mfonts.fontManager.ttflist}
    fnames_proplot -= fnames_proplot_ttc
    if not fnames_all >= fnames_proplot:
        warnings._warn_proplot('Rebuilding font cache.')
        if hasattr(mfonts.fontManager, 'addfont'):
            # New API lets us add font files manually
            for fname in fnames_proplot:
                mfonts.fontManager.addfont(fname)
            mfonts.json_dump(mfonts.fontManager, mfonts._fmcache)
        else:
            # Old API requires us to modify TTFPATH
            # NOTE: Previously we tried to modify TTFPATH before importing
            # font manager with hope that it would load proplot fonts on
            # initialization. But 99% of the time font manager just imports
            # the FontManager from cache, so this doesn't work.
            paths = ':'.join(paths_proplot)
            if 'TTFPATH' not in os.environ:
                os.environ['TTFPATH'] = paths
            elif paths not in os.environ['TTFPATH']:
                os.environ['TTFPATH'] += ':' + paths
            mfonts._rebuild()

    # Remove ttc files *after* rebuild
    mfonts.fontManager.ttflist = [
        font for font in mfonts.fontManager.ttflist
        if os.path.splitext(font.fname)[1] != '.ttc'
    ]


# Initialize .proplotrc file
_user_rc_file = os.path.join(os.path.expanduser('~'), '.proplotrc')
if not os.path.exists(_user_rc_file):
    _write_defaults(_user_rc_file)

# Initialize customization folders
_rc_folder = os.path.join(os.path.expanduser('~'), '.proplot')
if not os.path.isdir(_rc_folder):
    os.mkdir(_rc_folder)
for _rc_sub in ('cmaps', 'cycles', 'colors', 'fonts'):
    _rc_sub = os.path.join(_rc_folder, _rc_sub)
    if not os.path.isdir(_rc_sub):
        os.mkdir(_rc_sub)

# Convert colormaps that *should* be LinearSegmented from Listed
for _name in ('viridis', 'plasma', 'inferno', 'magma', 'cividis', 'twilight'):
    _cmap = pcolors._cmapdict.get(_name, None)
    if _cmap and isinstance(_cmap, pcolors.ListedColormap):
        del pcolors._cmapdict[_name]
        pcolors._cmapdict[_name] = pcolors.LinearSegmentedColormap.from_list(
            _name, _cmap.colors, cyclic=(_name == 'twilight')
        )

# Register objects and configure settings
with timers._benchmark('cmaps'):
    register_cmaps(default=True)

with timers._benchmark('cycles'):
    register_cycles(default=True)

with timers._benchmark('colors'):
    register_colors(default=True)

with timers._benchmark('fonts'):
    register_fonts()

with timers._benchmark('rc'):
    _ = rc_configurator()

#: Instance of `rc_configurator`. This is used to change global settings.
#: See :ref:`Configuring proplot` for details.
rc = rc_configurator()
