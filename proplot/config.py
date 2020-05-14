#!/usr/bin/env python3
"""
Tools for setting up ProPlot and configuring global settings.
See the :ref:`configuration guide <ug_config>` for details.
"""
# NOTE: The matplotlib analogue to this file is actually in __init__.py
# but it makes more sense to have all the setup actions in a separate file
# so the namespace of the top-level module is unpolluted.
# NOTE: Why also load colormaps and cycles in this file and not colors.py?
# Because I think it makes sense to have all the code that "runs" (i.e. not
# just definitions) in the same place, and I was having issues with circular
# dependencies and where import order of __init__.py was affecting behavior.
import re
import os
import numpy as np
import matplotlib as mpl
import matplotlib.font_manager as mfonts
import matplotlib.colors as mcolors
import matplotlib.style.core as mstyle
import matplotlib.cbook as cbook
import numbers
import cycler
from collections import namedtuple
from . import colors as pcolors
from .utils import units, to_xyz
from .internals import ic  # noqa: F401
from .internals import rcsetup, docstring, timers, warnings, _not_none
try:
    from IPython import get_ipython
except ImportError:
    def get_ipython():
        return

__all__ = [
    'rc', 'rc_configurator',
    'register_cmaps', 'register_cycles', 'register_colors', 'register_fonts',
    'config_inline_backend', 'use_style',
    'inline_backend_fmt',  # deprecated
]

# Disable mathtext "missing glyph" warnings
import matplotlib.mathtext  # noqa
import logging
logger = logging.getLogger('matplotlib.mathtext')
logger.setLevel(logging.ERROR)  # suppress warnings!

# Dictionaries used to track custom proplot settings
_Context = namedtuple('Context', ('mode', 'kwargs', 'rc_new', 'rc_old'))
rc_params = mpl.rcParams  # PEP8 4 lyfe
rc_quick = rcsetup._rc_quick_default.copy()
rc_added = rcsetup._rc_added_default.copy()

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

_config_docstring = """
user : bool, optional
    Whether to reload user {name}. Default is ``True``.
default : bool, optional
    Whether to reload default proplot {name}. Default is ``False``.
"""
docstring.snippets['register_cmaps.params'] = _config_docstring.format(name='colormaps')
docstring.snippets['register_cycles.params'] = _config_docstring.format(name='cycles')
docstring.snippets['register_colors.params'] = _config_docstring.format(name='colors')
docstring.snippets['rc.params'] = """
local : bool, optional
    Whether to reload ``.proplotrc`` settings in this directory and parent
    directories. Default is ``True``.
user : bool, optional
    Whether to reload ``~/.proplotrc`` user settings. Default is ``True``.
default : bool, optional
    Whether to reload default proplot settings. Default is ``True``.
"""

docstring.snippets['register.ext_table'] = """
Valid file extensions are as follows:

==================  =====================================================================================================================================================================================================================
Extension           Description
==================  =====================================================================================================================================================================================================================
``.hex``            List of HEX strings in any format (comma-separated, separate lines, with double quotes... anything goes).
``.xml``            XML files with ``<Point .../>`` tags specifying ``x``, ``r``, ``g``, ``b``, and (optionally) ``o`` parameters, where ``x`` is the coordinate and the rest are the red, blue, green, and opacity channel values.
``.rgb``, ``.txt``  3-4 column table of red, blue, green, and (optionally) opacity channel values, delimited by commas or spaces. If values larger than 1 are detected, they are assumed to be on the 0-255 scale and are divided by 255.
==================  =====================================================================================================================================================================================================================
"""  # noqa: E501


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
            elif isinstance(value, (str, numbers.Number, NoneType)):
                value = str(value)
            elif isinstance(value, (list, tuple)) and all(
                isinstance(val, (str, numbers.Number)) for val in value
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

    # Fill empty _rc_added_default values with _rc_quick_default
    # They are only allowed to be None in the *default dictionaries* because
    # they are immediately overwritten. However if users try to set them as
    # None in a .proplotrc file, may trigger error down the line.
    # TODO: Never set
    rc_parents = {
        child: parent
        for parent, children in rcsetup._rc_children.items()
        for child in children
    }
    rc_added_filled = rcsetup._rc_added_default.copy()
    for key, value in rcsetup._rc_added_default.items():
        if value is None:
            try:
                parent = rc_parents[key]
            except KeyError:
                raise RuntimeError(
                    f'rc_added param {key!r} has default value of None '
                    'but has no rcParmsShort parent!'
                )
            if parent in rcsetup._rc_quick_default:
                value = rcsetup._rc_quick_default[parent]
            elif parent in rcsetup._rc_params_default:
                value = rcsetup._rc_params_default[parent]
            else:
                value = rc_params[parent]
            rc_added_filled[key] = value

    with open(filename, 'w') as fh:
        fh.write(f"""
#---------------------------------------------------------------------
# Use this file to change the default proplot and matplotlib settings
# The syntax is mostly the same as for matplotlibrc files
# For descriptions of each setting see:
# https://proplot.readthedocs.io/en/latest/configuration.html
# https://matplotlib.org/3.1.1/tutorials/introductory/customizing.html
#---------------------------------------------------------------------
# ProPlot quick settings
{_tabulate(rcsetup._rc_quick_default)}

# ProPlot added settings
{_tabulate(rc_added_filled)}

# Matplotlib settings
{_tabulate(rcsetup._rc_params_default)}
""".strip())


class rc_configurator(object):
    """
    Magical abstract class for managing matplotlib's
    `builtin settings <rc_params>`_, ProPlot's
    :ref:`added settings <rc_added>`, and :ref:`quick settings <rc_quick>`.
    When ProPlot is imported, this class is instantiated as the `rc` object
    and the ProPlot default settings and ``.proplotrc`` user overrides
    are applied. To modify these settings, use the `rc` object.
    See the :ref:`configuration guide <ug_config>` for details.
    """
    def __repr__(self):  # encapsulate params in temporary class
        rcdict = type('rc', (dict,), {})(rc_quick)
        string = type(rc_params).__repr__(rcdict)
        return string.strip()[:-2] + ',\n    ... <rcParams> ...\n    })'

    def __str__(self):
        rcdict = type('rc', (dict,), {})(rc_quick)
        string = type(rc_params).__str__(rcdict)
        return string + '\n... <rcParams> ...'

    def __iter__(self):  # lets us build dict
        """
        Iterate over keys and values of matplotlib and proplot settings.
        """
        for key in sorted((*rc_quick, *rc_added, *rc_params)):
            yield key, self[key]

    def __contains__(self, key):
        """
        Test whether key exists as matplotlib or proplot setting.
        """
        return key in rc_quick or key in rc_added or key in rc_params

    @docstring.add_snippets
    def __init__(self, local=True, user=True, default=True):
        """
        Parameters
        ----------
        %(rc.params)s
        """
        self._context = []
        self.reset(local=local, user=user, default=default)

    def __enter__(self):
        """
        Apply settings from the most recent context block.
        """
        if not self._context:
            raise RuntimeError(
                'rc object must be initialized for context block '
                'using rc.context().'
            )
        context = self._context[-1]
        kwargs = context.kwargs
        rc_new = context.rc_new  # used for context-based _get_item
        rc_old = context.rc_old  # used to re-apply settings without copying whole dict
        for key, value in kwargs.items():
            kw_quick, kw_added, kw_params = self._get_param_dicts(key, value)
            for rc_dict, kw_new in zip(
                (rc_quick, rc_added, rc_params),
                (kw_quick, kw_added, kw_params),
            ):
                for key, value in kw_new.items():
                    rc_old[key] = rc_dict[key]
                    rc_new[key] = rc_dict[key] = value

    def __exit__(self, *args):  # noqa: U100
        """
        Restore settings from the most recent context block.
        """
        if not self._context:
            raise RuntimeError(
                'rc object must be initialized for context block '
                'using rc.context().'
            )
        context = self._context[-1]
        for key, value in context.rc_old.items():
            kw_quick, kw_added, kw_params = self._get_param_dicts(key, value)
            rc_quick.update(kw_quick)
            rc_added.update(kw_added)
            rc_params.update(kw_params)
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
            return super().__getattribute__(attr)
        else:
            return self[attr]

    def __getitem__(self, key):
        """
        Return a `builtin matplotlib setting <rc_params>`_,
        a ProPlot :ref:`added setting <rc_added>`,
        or a :ref:`quick setting <rc_quick>`.
        """
        key = self._sanitize_key(key)
        for kw in (rc_quick, rc_added, rc_params):
            try:
                return kw[key]
            except KeyError:
                continue
        raise KeyError(f'Invalid setting name {key!r}.')

    def __setattr__(self, attr, value):
        """
        Pass the attribute and value to `~rc_configurator.__setitem__`.
        """
        if attr[:1] == '_':
            super().__setattr__(attr, value)
        else:
            self.__setitem__(attr, value)

    def __setitem__(self, key, value):
        """
        Modify a `builtin matplotlib setting <rc_params>`_,
        a ProPlot :ref:`added setting <rc_added>`,
        or a :ref:`quick setting <rc_quick>`.
        """
        kw_quick, kw_added, kw_params = self._get_param_dicts(key, value)
        rc_quick.update(kw_quick)
        rc_added.update(kw_added)
        rc_params.update(kw_params)

    def _get_item(self, key, mode=None):
        """
        As with `~rc_configurator.__getitem__` but the search is limited
        based on the context mode and ``None`` is returned if the key is not
        found in the dictionaries.
        """
        if mode is None:
            mode = min((context.mode for context in self._context), default=0)
        cache = tuple(context.rc_new for context in self._context)
        if mode == 0:
            rcdicts = (*cache, rc_quick, rc_added, rc_params)
        elif mode == 1:
            rcdicts = (*cache, rc_quick, rc_added)  # custom only!
        elif mode == 2:
            rcdicts = (*cache,)
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

    def _get_param_dicts(self, key, value):
        """
        Return dictionaries for updating the `rc_quick`, `rc_added`,
        and `rc_params` properties associated with this key.
        """
        kw_params = {}  # builtin properties that global setting applies to
        kw_added = {}  # custom properties that global setting applies to
        kw_quick = {}  # short name properties
        key = self._sanitize_key(key)
        children = rcsetup._rc_children.get(key, ())

        # Permit arbitary units for builtin matplotlib params
        # See: https://matplotlib.org/users/customizing.html, props matching
        # the below strings use the units 'points'.
        # TODO: Incorporate into more sophisticated validation system
        if any(REGEX_POINTS.match(_) for _ in (children or (key,))):
            try:
                value = self._scale_font(value)
            except KeyError:
                value = units(value, 'pt')

        # Deprecations
        if key in rcsetup._rc_removed:
            version = rcsetup._rc_removed[key]
            warnings._warn_proplot(
                f'rc setting {key!r} was removed in version {version}.'
            )
            return {}, {}, {}
        if key in rcsetup._rc_renamed:
            key_new, version = rcsetup._rc_renamed[key]
            warnings._warn_proplot(
                f'rc setting {key!r} was renamed to {key_new} in version {version}.'
            )
            key = key_new

        # Special key: configure inline backend
        if key == 'inlinefmt':
            config_inline_backend(value)

        # Special key: apply stylesheet
        elif key == 'style':
            if value is not None:
                kw_params, kw_added = _get_style_dicts(value, infer=True)

        # Cycler
        elif key == 'cycle':
            colors = _get_cycle_colors(value)
            kw_params['patch.facecolor'] = colors[0]
            kw_params['axes.prop_cycle'] = cycler.cycler('color', colors)

        # Zero linewidth almost always means zero tick length
        # TODO: Document this feature
        elif key == 'linewidth' and value == 0:
            _, ikw_added, ikw_params = self._get_param_dicts('ticklen', 0)
            kw_added.update(ikw_added)
            kw_params.update(ikw_params)

        # Tick length/major-minor tick length ratio
        elif key in ('tick.len', 'tick.lenratio'):
            if key == 'tick.len':
                ticklen = value
                ratio = rc_added['tick.lenratio']
            else:
                ticklen = rc_added['tick.len']
                ratio = value
            kw_params['xtick.minor.size'] = ticklen * ratio
            kw_params['ytick.minor.size'] = ticklen * ratio

        # Spine width/major-minor tick width ratio
        elif key in ('linewidth', 'tick.ratio'):
            if key == 'linewidth':
                tickwidth = value
                ratio = rc_added['tick.ratio']
            else:
                tickwidth = rc_quick['linewidth']
                ratio = value
            kw_params['xtick.minor.width'] = tickwidth * ratio
            kw_params['ytick.minor.width'] = tickwidth * ratio

        # Gridline width
        elif key in ('grid.linewidth', 'grid.ratio'):
            if key == 'grid.linewidth':
                gridwidth = value
                ratio = rc_added['grid.ratio']
            else:
                gridwidth = rc_params['grid.linewidth']
                ratio = value
            kw_added['gridminor.linewidth'] = gridwidth * ratio

        # Gridline toggling, complicated because of the clunky way this is
        # implemented in matplotlib. There should be a gridminor setting!
        elif key in ('grid', 'gridminor'):
            ovalue = rc_params['axes.grid']
            owhich = rc_params['axes.grid.which']

            # Instruction is to turn off gridlines
            if not value:
                # Gridlines are already off, or they are on for the particular
                # ones that we want to turn off. Instruct to turn both off.
                if (
                    not ovalue
                    or (key == 'grid' and owhich == 'major')
                    or (key == 'gridminor' and owhich == 'minor')
                ):
                    which = 'both'  # disable both sides
                # Gridlines are currently on for major and minor ticks, so we
                # instruct to turn on gridlines for the one we *don't* want off
                elif owhich == 'both':  # and ovalue is True, as already tested
                    # if gridminor=False, enable major, and vice versa
                    value = True
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
                if (
                    owhich == 'both'
                    or (key == 'grid' and owhich == 'minor')
                    or (key == 'gridminor' and owhich == 'major')
                ):
                    which = 'both'
                # Gridlines are off for both, or off for the ones that we
                # don't want to turn on. We can just turn on these ones.
                else:
                    which = owhich

            # Finally apply settings
            kw_params['axes.grid'] = value
            kw_params['axes.grid.which'] = which

        # Update setting in dictionary, detect invalid keys
        if key in rc_quick:
            kw_quick[key] = value
        elif key in rc_added:
            kw_added[key] = value
        elif key in rc_params:
            kw_params[key] = value
        else:
            raise KeyError(f'Invalid rc key {key!r}.')

        # Update linked settings
        for key in children:
            if key in rc_added:
                kw_added[key] = value
            else:
                kw_params[key] = value
        return kw_quick, kw_added, kw_params

    @staticmethod
    def _get_local_paths():
        """
        Return locations of local proplotrc files in this directory
        and in parent directories.
        """
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
        return paths[::-1]  # sort from decreasing to increasing importantce

    @staticmethod
    def _get_user_path():
        """
        Return location of user proplotrc file.
        """
        return os.path.join(os.path.expanduser('~'), '.proplotrc')

    @staticmethod
    def _sanitize_key(key):
        """
        Ensure string and convert keys with omitted dots.
        """
        if not isinstance(key, str):
            raise KeyError(f'Invalid key {key!r}. Must be string.')
        if '.' not in key and key not in rc_quick:  # speedup
            key = rcsetup._rc_nodots.get(key, key)
        return key.lower()

    @staticmethod
    def _scale_font(size):
        """
        Translate font size to numeric.
        """
        if isinstance(size, str):
            size = rc_params['font.size'] * mfonts.font_scalings[size]
        return size

    def _update_from_file(self, path):
        """
        Apply updates from a file. This is largely copied from matplotlib.

        Parameters
        ----------
        path : str
            The path.
        """
        path = os.path.expanduser(path)
        added = set()
        with open(path, 'r') as fh:
            for cnt, line in enumerate(fh):
                # Parse line and ignore comments
                stripped = line.split('#', 1)[0].strip()
                if not stripped:
                    continue
                pair = stripped.split(':', 1)
                if len(pair) != 2:
                    warnings._warn_proplot(
                        f'Illegal line #{cnt + 1} in file {path!r}:\n{line!r}"'
                    )
                    continue

                # Get key value pair
                key, val = pair
                key = key.strip()
                val = val.strip()
                if key in added:
                    warnings._warn_proplot(
                        f'Duplicate key {key!r} on line #{cnt + 1} in file {path!r}.'
                    )
                added.add(key)

                # *Very primitive* type conversion system. Just check proplot
                # settings (they are all simple/scalar) and leave rc_params alone.
                # TODO: Add built-in validation just like matplotlib RcParams
                if key in rc_quick or key in rc_added:
                    if not val:
                        val = None  # older proplot versions supported this
                    elif val in ('True', 'False', 'None'):
                        val = eval(val)  # rare case where eval is o.k.
                    else:
                        try:
                            val = float(val) if '.' in val else int(val)
                        except ValueError:
                            pass  # retain string

                # Add to dictionaries
                try:
                    kw_quick, kw_added, kw_params = self._get_param_dicts(key, val)
                except KeyError:
                    warnings._warn_proplot(
                        f'Invalid key {key!r} on line #{cnt} in file {path!r}.'
                    )
                else:
                    rc_quick.update(kw_quick)
                    rc_added.update(kw_added)
                    rc_params.update(kw_params)

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
        if cat not in rcsetup._rc_categories:
            raise ValueError(
                f'Invalid rc category {cat!r}. Valid categories are '
                ', '.join(map(repr, rcsetup._rc_categories)) + '.'
            )
        kw = {}
        mode = 0 if not context else None
        for rcdict in (rc_added, rc_params):
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
            "with as" block when called with ``context=True``.

            The options are as follows:

            0. Matplotlib's `builtin settings <rc_params>`_,
               ProPlot's :ref:`added settings <rc_added>`,
               and the :ref:`quick settings <rc_quick>` are all returned,
               whether or not `~rc_configurator.context` has changed them.
            1. *Unchanged* `matplotlib settings <rc_params>`_
               return ``None``. All of ProPlot's
               :ref:`added settings <rc_added>`
               and the :ref:`quick settings <rc_quick>` are returned whether
               or not `~rc_configurator.context` has changed them.  This is
               used in the `~proplot.axes.Axes.__init__` call to
               `~proplot.axes.Axes.format`. When a lookup returns
               ``None``, `~proplot.axes.Axes.format` does not apply it.
            2. All unchanged settings return ``None``. This is used during
               user calls to `~proplot.axes.Axes.format`.

        Example
        -------
        The below applies settings to axes in a specific figure using
        `~rc_configurator.context`.

        >>> import proplot as plot
        >>> with plot.rc.context(linewidth=2, ticklen=5):
        >>>     fig, ax = plot.subplots()
        >>>     ax.plot(data)

        The below applies settings to a specific axes using `~proplot.axes.Axes.format`,
        which uses `~rc_configurator.context` internally.

        >>> import proplot as plot
        >>> fig, ax = plot.subplots()
        >>> ax.format(linewidth=2, ticklen=5)
        """
        if mode not in range(3):
            raise ValueError(f'Invalid mode {mode!r}.')
        for arg in args:
            if not isinstance(arg, dict):
                raise ValueError('Non-dictionary argument {arg!r}.')
            kwargs.update(arg)
        tup = _Context(mode=mode, kwargs=kwargs, rc_new={}, rc_old={})
        self._context.append(tup)
        return self

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
        Same as `dict.keys`.
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
            self.__setitem__(prefix + key, value)

    @docstring.add_snippets
    def reset(self, local=True, user=True, default=True):
        """
        Reset the configurator to its initial state.

        Parameters
        ----------
        %(rc.params)s
        """
        # Always remove context objects
        self._context.clear()

        # Update from default settings
        # NOTE: see _remove_blacklisted_style_params bugfix
        # TODO: Make proplot settings compatible with matplotlib stylesheets.
        # Tweak rc_added so they are *consistent* with properties in rcParams,
        # load from rcParamsOrig rather than rcParamsDefault, implement proplot
        # "style", and do not use rc_quick to apply any default settings -- this
        # should be for user convenience only and shold not be used internally.
        if default:
            rc_params.update(_get_style_dicts('original', infer=False))
            rc_params.update(rcsetup._rc_params_default)  # proplot changes
            rc_added.update(rcsetup._rc_added_default)  # proplot custom params
            rc_quick.update(rcsetup._rc_quick_default)  # proplot quick params
            for dict_ in (rc_quick, rc_added):
                for key, value in dict_.items():
                    _, kw_added, kw_params = self._get_param_dicts(key, value)
                    rc_added.update(kw_added)
                    rc_params.update(kw_params)

        # Update from user home
        user_path = None
        if user:
            user_path = self._get_user_path()
            if os.path.isfile(user_path):
                self._update_from_file(user_path)

        # Update from local paths
        if local:
            local_paths = self._get_local_paths()
            for path in local_paths:
                if path == user_path:  # local files always have precedence
                    continue
                self._update_from_file(path)

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


def _get_cycle_colors(cycle):
    """
    Update the color cycle.
    """
    try:
        colors = pcolors._cmap_database[cycle].colors
    except (KeyError, AttributeError):
        cycles = sorted(
            name for name, cmap in pcolors._cmap_database.items()
            if isinstance(cmap, pcolors.ListedColormap)
        )
        raise ValueError(
            f'Invalid cycle name {cycle!r}. Options are: '
            + ', '.join(map(repr, cycles)) + '.'
        )
    return colors


def _get_default_dict():
    """
    Get the default rc parameters dictionary with deprecated parameters filtered.
    """
    # NOTE: Use RcParams update to filter and translate deprecated settings
    # before actually applying them to rcParams down pipeline. This way we can
    # suppress warnings for deprecated default params but still issue warnings
    # when user-supplied stylesheets have deprecated params.
    # WARNING: Some deprecated rc params remain in dictionary as None so we
    # filter them out. Beware if hidden attribute changes.
    rcdict = _get_filtered_dict(mpl.rcParamsDefault, warn=False)
    with cbook._suppress_matplotlib_deprecation_warning():
        rcdict = dict(mpl.RcParams(rcdict))
    for attr in ('_deprecated_remain_as_none', '_deprecated_set'):
        if hasattr(mpl, attr):  # _deprecated_set is in matplotlib before v3
            for deprecated in getattr(mpl, attr):
                rcdict.pop(deprecated, None)
    return rcdict


def _get_filtered_dict(rcdict, warn=True):
    """
    Filter out blacklisted style parameters.
    """
    # NOTE: This implements bugfix: https://github.com/matplotlib/matplotlib/pull/17252
    # This fix is *critical* for proplot because we always run style.use()
    # when the configurator is made. Without fix backend is reset every time
    # you import proplot in jupyter notebooks. So apply retroactively.
    rcdict_filtered = {}
    for key in rcdict:
        if key in mstyle.STYLE_BLACKLIST:
            if warn:
                warnings._warn_proplot(
                    f'Dictionary includes a parameter, {key!r}, that is not related '
                    'to style. Ignoring.'
                )
        else:
            rcdict_filtered[key] = rcdict[key]
    return rcdict_filtered


def _get_style_dicts(style, infer=False):
    """
    Return a dictionary of settings belonging to the requested style(s). If `infer`
    is ``True``, two dictionaries are returned, where the second contains custom
    ProPlot settings "inferred" from the matplotlib settings.
    """
    # NOTE: This is adapted from matplotlib source for the following changes:
    # 1. Add 'original' option. Like rcParamsOrig except we also *reload*
    #    from user matplotlibrc file.
    # 2. When the style is changed we reset to the *default* state ignoring
    #    matplotlibrc. Matplotlib applies styles on top of current state
    #    (including matplotlibrc changes and runtime rcParams changes) but
    #    IMO the word 'style' implies a *rigid* static format.
    # 3. Add a separate function that returns lists of style dictionaries so
    #    that we can modify the active style in a context block. ProPlot context
    #    is more conservative than matplotlib's rc_context because it gets
    #    called a lot (e.g. every time you make an axes and every format() call).
    #    Instead of copying the entire rcParams dict we just track the keys
    #    that were changed.
    style_aliases = {
        '538': 'fivethirtyeight',
        'mpl20': 'default',
        'mpl15': 'classic',
        'original': mpl.matplotlib_fname(),
    }

    # Always apply the default style *first* so styles are rigid
    kw_params = _get_default_dict()
    if style == 'default' or style is mpl.rcParamsDefault:
        return kw_params

    # Apply "pseudo" default properties. Pretend some proplot settings are part of
    # the matplotlib specification so they propagate to other styles.
    kw_params['font.family'] = 'sans-serif'
    kw_params['font.sans-serif'] = rcsetup._rc_params_default['font.sans-serif']

    # Apply user input style(s) one by one
    # NOTE: Always use proplot fonts if style does not explicitly set them.
    if isinstance(style, str) or isinstance(style, dict):
        styles = [style]
    else:
        styles = style
    for style in styles:
        if isinstance(style, dict):
            kw = style
        elif isinstance(style, str):
            style = style_aliases.get(style, style)
            if style in mstyle.library:
                kw = mstyle.library[style]
            else:
                try:
                    kw = mpl.rc_params_from_file(style, use_default_template=False)
                except IOError:
                    raise IOError(
                        f'Style {style!r} not found in the style library and input is '
                        'not a valid URL or path. Available styles are: '
                        + ', '.join(map(repr, mstyle.available)) + '.'
                    )
        else:
            raise ValueError(f'Invalid style {style!r}. Must be string or dictionary.')
        kw = _get_filtered_dict(kw, warn=True)
        kw_params.update(kw)

    # Infer proplot params from stylesheet params
    if infer:
        kw_added = _infer_added_params(kw_params)
        return kw_params, kw_added
    else:
        return kw_params


def _infer_added_params(kw_params):
    """
    Infer values for proplot's "added" parameters from stylesheets.
    """
    kw_added = {}
    mpl_to_proplot = {
        'font.size': ('tick.labelsize',),
        'axes.titlesize': (
            'abc.size', 'suptitle.size', 'title.size',
            'leftlabel.size', 'rightlabel.size',
            'toplabel.size', 'bottomlabel.size',
        ),
        'axes.facecolor': ('geoaxes.facecolor',),
        'text.color': (
            'abc.color', 'suptitle.color', 'tick.labelcolor', 'title.color',
            'leftlabel.color', 'rightlabel.color',
            'toplabel.color', 'bottomlabel.color',
        ),
    }
    for key, params in mpl_to_proplot.items():
        if key in kw_params:
            value = kw_params[key]
            for param in params:
                kw_added[param] = value
    return kw_added


def use_style(style):
    """
    Apply the `matplotlib style(s) \
<https://matplotlib.org/tutorials/introductory/customizing.html>`__
    with `matplotlib.style.use`. This function is
    called when you modify the :rcraw:`style` property.

    Parameters
    ----------
    style : str, dict, or list thereof
        The matplotlib style name(s) or stylesheet filename(s), or dictionary(s)
        of settings. Use ``'default'`` to apply matplotlib default settings and
        ``'original'`` to include settings from your user ``matplotlibrc`` file.
    """
    # NOTE: This function is not really necessary but makes proplot's
    # stylesheet-supporting features obvious. Plus changing the style does
    # so much *more* than changing rc params or quick settings, so it is
    # nice to have dedicated function instead of just another rc_param name.
    kw_params, kw_added = _get_style_dicts(style, infer=True)
    rc_params.update(kw_params)
    rc_added.update(kw_added)


@docstring.add_snippets
def register_cmaps(user=True, default=False):
    """
    Register colormaps packaged with ProPlot or saved to the
    ``~/.proplot/cmaps`` folder. This is called on import.
    Colormaps are registered according to their filenames -- for example,
    ``name.xyz`` will be registered as ``'name'``.

    %(register.ext_table)s

    To visualize the registered colormaps, use `~proplot.demos.show_cmaps`.

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
        pcolors._cmap_database[cmap.name] = cmap


@docstring.add_snippets
def register_cycles(user=True, default=False):
    """
    Register color cycles packaged with ProPlot or saved to the
    ``~/.proplot/cycles`` folder. This is called on import. Color cycles
    are registered according to their filenames -- for example, ``name.hex``
    will be registered as ``'name'``.

    %(register.ext_table)s

    To visualize the registered color cycles, use `~proplot.demos.show_cycles`.

    Parameters
    ----------
    %(register_cycles.params)s
    """
    for _, dirname, filename in _iter_data_paths('cycles', user=user, default=default):
        path = os.path.join(dirname, filename)
        cmap = pcolors.ListedColormap.from_file(path, warn_on_failure=True)
        if not cmap:
            continue
        pcolors._cmap_database[cmap.name] = cmap


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

    To visualize the registered colors, use `~proplot.demos.show_colors`.

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

    To visualize the registered fonts, use `~proplot.demos.show_fonts`.
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

    # Remove ttc files and 'Thin' fonts *after* rebuild
    # NOTE: 'Thin' filter is ugly kludge but without this matplotlib picks up on
    # Roboto thin ttf files installed on the RTD server when compiling docs.
    mfonts.fontManager.ttflist = [
        font for font in mfonts.fontManager.ttflist
        if os.path.splitext(font.fname)[1] != '.ttc'
        or 'Thin' in os.path.basename(font.fname)
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
    _cmap = pcolors._cmap_database.get(_name, None)
    if _cmap and isinstance(_cmap, pcolors.ListedColormap):
        del pcolors._cmap_database[_name]
        pcolors._cmap_database[_name] = pcolors.LinearSegmentedColormap.from_list(
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
#: See the :ref:`configuration guide <ug_config>` for details.
rc = _

# Modify N of existing colormaps because ProPlot settings may have changed
# image.lut. We have to register colormaps and cycles first so that the 'cycle'
# property accepts named cycles registered by ProPlot. No performance hit here.
lut = rc['image.lut']
for cmap in pcolors._cmap_database.values():
    if isinstance(cmap, mcolors.LinearSegmentedColormap):
        cmap.N = lut

# Deprecated
inline_backend_fmt = warnings._rename_obj('inline_backend_fmt', config_inline_backend)
