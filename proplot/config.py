#!/usr/bin/env python3
"""
Tools for setting up ProPlot and configuring global settings.
See the :ref:`configuration guide <ug_config>` for details.
"""
# NOTE: The matplotlib analogue to this file is actually __init__.py
# but it makes more sense to have all the setup actions in a separate file
# so the namespace of the top-level module is unpolluted.
# NOTE: Why also load colormaps and cycles in this file and not colors.py?
# Because I think it makes sense to have all the code that "runs" (i.e. not
# just definitions) in the same place, and I was having issues with circular
# dependencies and where import order of __init__.py was affecting behavior.
import logging
import os
import re
from collections import namedtuple

import cycler
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.font_manager as mfonts
import matplotlib.mathtext  # noqa
import matplotlib.rcsetup as msetup
import matplotlib.style.core as mstyle
import numpy as np

from . import colors as pcolors
from .internals import ic  # noqa: F401
from .internals import _not_none, docstring, rcsetup, timers, warnings
from .utils import to_xyz, units

try:
    from matplotlib.cbook import _suppress_matplotlib_deprecation_warning  # mpl<3.4.0
except ImportError:
    from matplotlib._api import suppress_matplotlib_deprecation_warning as _suppress_matplotlib_deprecation_warning  # noqa: E501

try:
    from IPython import get_ipython
except ImportError:
    def get_ipython():
        return

__all__ = [
    'rc',
    'RcConfigurator',
    'register_cmaps',
    'register_cycles',
    'register_colors',
    'register_fonts',
    'config_inline_backend',
    'use_style',
    'inline_backend_fmt',  # deprecated
]

logger = logging.getLogger('matplotlib.mathtext')
logger.setLevel(logging.ERROR)  # suppress warnings!

# Dictionaries used to track settings
rc_proplot = rcsetup._rc_proplot_default.copy()
rc_matplotlib = mpl.rcParams  # PEP8 4 lyfe
_RcContext = namedtuple('RcContext', ('mode', 'kwargs', 'rc_new', 'rc_old'))
RcParams = mpl.RcParams  # the special class

# Misc constants
# TODO: Use explicit validators for specific settings like matplotlib.
REGEX_STRING = re.compile('\\A(\'.*\'|".*")\\Z')
REGEX_POINTS = re.compile(
    r'\A(?!font.mono|colorbar|subplots|pdf|ps).*(width|space|size|pad|len)\Z'
)
FONT_KEYS = (
    'abc.size', 'axes.labelsize', 'axes.titlesize', 'figure.titlesize',
    'xtick.labelsize', 'ytick.labelsize', 'tick.labelsize', 'grid.labelsize',
    'bottomlabel.size', 'leftlabel.size', 'rightlabel.size', 'toplabel.size',
    'suptitle.size', 'title.size', 'text.labelsize', 'text.titlesize',
    'legend.fontsize', 'legend.title_fontsize'
)
COLORS_OPEN = {}  # populated during register_colors
COLORS_XKCD = {}  # populated during register_colors
COLORS_BASE = {
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
COLORS_REMOVE = (  # filter these out, let's try to be professional here...
    'shit', 'poop', 'poo', 'pee', 'piss', 'puke', 'vomit', 'snot',
    'booger', 'bile', 'diarrhea',
)
COLORS_TRANSLATE = (  # prevent registering similar-sounding names
    ('/', ' '),
    ("'s", ''),
    ('forrest', 'forest'),  # survey typo?
    ('reddish', 'red'),  # remove 'ish'
    ('purplish', 'purple'),
    ('bluish', 'blue'),
    ('ish ', ' '),
    ('grey', 'gray'),  # 'Murica
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
COLORS_ADD = (
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
        for prefix in ('', 'light ', 'dark ', 'medium ', 'pale ')
    )
)


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


class RcConfigurator(object):
    """
    Dictionary-like class for managing matplotlib's `builtin settings <rc_matplotlib>`_
    and ProPlot's :ref:`added settings <rc_proplot>`. When ProPlot is imported,
    this class is instantiated as the `rc` object and the ProPlot default settings
    and ``.proplotrc`` user overrides are applied. See the
    :ref:`configuration guide <ug_config>` for details.
    """
    def __repr__(self):
        rcdict = type('rc', (dict,), {})({  # encapsulate params in temporary class
            key: value for key, value in rc_proplot.items()
            if '.' not in key  # show short names
        })
        string = type(rc_matplotlib).__repr__(rcdict)
        return string.strip()[:-2] + ',\n    ... <rcParams> ...\n    })'

    def __str__(self):
        rcdict = type('rc', (dict,), {})({
            key: value for key, value in rc_proplot.items()
            if '.' not in key  # show short names
        })
        string = type(rc_matplotlib).__str__(rcdict)
        return string + '\n... <rcParams> ...'

    def __iter__(self):  # lets us build dict
        """
        Iterate over keys and values of matplotlib and proplot settings.
        """
        for key in sorted((*rc_proplot, *rc_matplotlib)):
            yield key, self[key]

    def __contains__(self, key):
        """
        Test whether key exists as matplotlib or proplot setting.
        """
        return key in rc_proplot or key in rc_matplotlib

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
            kw_proplot, kw_matplotlib = self._get_synced_params(key, value)
            for rc_dict, kw_new in zip(
                (rc_proplot, rc_matplotlib),
                (kw_proplot, kw_matplotlib),
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
            kw_proplot, kw_matplotlib = self._get_synced_params(key, value)
            rc_proplot.update(kw_proplot)
            rc_matplotlib.update(kw_matplotlib)
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
        Pass the attribute to `~RcConfigurator.__getitem__` and return
        the result.
        """
        if attr[:1] == '_':
            return super().__getattribute__(attr)
        else:
            return self[attr]

    def __getitem__(self, key):
        """
        Return a `builtin matplotlib setting <rc_matplotlib>`_
        or a ProPlot :ref:`added setting <rc_proplot>`.
        """
        key = self._sanitize_key(key)
        if key is None:  # means key was *removed*, warnings was issued
            return None
        for kw in (rc_proplot, rc_matplotlib):
            try:
                return kw[key]
            except KeyError:
                continue
        raise KeyError(f'Invalid setting name {key!r}.')

    def __setattr__(self, attr, value):
        """
        Pass the attribute and value to `~RcConfigurator.__setitem__`.
        """
        if attr[:1] == '_':
            super().__setattr__(attr, value)
        else:
            self.__setitem__(attr, value)

    def __setitem__(self, key, value):
        """
        Modify a `builtin matplotlib setting <rc_matplotlib>`_ or
        a ProPlot :ref:`added setting <rc_proplot>`.
        """
        kw_proplot, kw_matplotlib = self._get_synced_params(key, value)
        rc_proplot.update(kw_proplot)
        rc_matplotlib.update(kw_matplotlib)

    def _get_item(self, key, mode=None):
        """
        As with `~RcConfigurator.__getitem__` but the search is limited
        based on the context mode and ``None`` is returned if the key is not
        found in the dictionaries.
        """
        if mode is None:
            mode = self._get_context_mode()
        cache = tuple(context.rc_new for context in self._context)
        if mode == 0:
            rcdicts = (*cache, rc_proplot, rc_matplotlib)
        elif mode == 1:
            rcdicts = (*cache, rc_proplot)  # custom only!
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

    def _get_context_mode(self):
        """
        Return highest (least permissive) context mode.
        """
        return max((context.mode for context in self._context), default=0)

    def _get_synced_params(self, key, value):
        """
        Return dictionaries for updating the `rc_proplot`
        and `rc_matplotlib` properties associated with this key.
        """
        key = self._sanitize_key(key)
        if key is None:  # means setting was removed
            return {}, {}, {}
        keys = (key,) + rcsetup._rc_children.get(key, ())  # settings to change
        value = self._sanitize_value(value)
        kw_proplot = {}  # custom properties that global setting applies to
        kw_matplotlib = {}  # builtin properties that global setting applies to

        # Permit arbitary units for builtin matplotlib params
        # Props matching the below strings use the units 'points'.
        # See: https://matplotlib.org/users/customizing.html
        # TODO: Incorporate into more sophisticated validation system
        if REGEX_POINTS.match(key):
            if key in FONT_KEYS and value in mfonts.font_scalings:
                pass
            elif key.startswith('legend') and not key.endswith('fontsize'):
                value = units(value, 'em')  # scaled by font size
            else:
                value = units(value, 'pt')  # untis points fontsize='10px'

        # Special key: configure inline backend
        if key == 'inlinefmt':
            config_inline_backend(value)

        # Special key: apply stylesheet
        elif key == 'style':
            if value is not None:
                kw_matplotlib, kw_proplot = _get_style_dicts(value, infer=True)

        # Cycler
        elif key == 'cycle':
            colors = _get_cycle_colors(value)
            kw_matplotlib['patch.facecolor'] = 'C0'
            kw_matplotlib['axes.prop_cycle'] = cycler.cycler('color', colors)

        # Turning bounding box on should turn border off and vice versa
        elif key in ('abc.bbox', 'title.bbox', 'abc.border', 'title.border'):
            if value:
                name, this = key.split('.')
                other = 'border' if this == 'bbox' else 'border'
                kw_proplot[name + '.' + other] = False

        # Fontsize
        # NOTE: Re-application of e.g. size='small' uses the updated 'font.size'
        elif key == 'font.size':
            kw_proplot.update({
                key: value for key, value in rc_proplot.items()
                if key in FONT_KEYS and value in mfonts.font_scalings
            })
            kw_matplotlib.update({
                key: value for key, value in rc_matplotlib.items()
                if key in FONT_KEYS and value in mfonts.font_scalings
            })

        # Zero linewidth almost always means zero tick length
        # TODO: Document this feature
        elif key == 'linewidth' and value == 0:
            ikw_proplot, ikw_matplotlib = self._get_synced_params('ticklen', 0)
            kw_proplot.update(ikw_proplot)
            kw_matplotlib.update(ikw_matplotlib)

        # Tick length/major-minor tick length ratio
        elif key in ('tick.len', 'tick.lenratio'):
            if key == 'tick.len':
                ticklen = value
                ratio = rc_proplot['tick.lenratio']
            else:
                ticklen = rc_proplot['tick.len']
                ratio = value
            kw_matplotlib['xtick.minor.size'] = ticklen * ratio
            kw_matplotlib['ytick.minor.size'] = ticklen * ratio

        # Spine width/major-minor tick width ratio
        elif key in ('linewidth', 'tick.ratio'):
            if key == 'linewidth':
                tickwidth = value
                ratio = rc_proplot['tick.ratio']
            else:
                tickwidth = rc_proplot['linewidth']
                ratio = value
            kw_matplotlib['xtick.minor.width'] = tickwidth * ratio
            kw_matplotlib['ytick.minor.width'] = tickwidth * ratio

        # Gridline width
        elif key in ('grid.linewidth', 'grid.ratio'):
            if key == 'grid.linewidth':
                gridwidth = value
                ratio = rc_proplot['grid.ratio']
            else:
                gridwidth = rc_matplotlib['grid.linewidth']
                ratio = value
            kw_proplot['gridminor.linewidth'] = gridwidth * ratio

        # Gridline toggling, complicated because of the clunky way this is
        # implemented in matplotlib. There should be a gridminor setting!
        elif key in ('grid', 'gridminor'):
            b = value
            ovalue = rc_matplotlib['axes.grid']
            owhich = rc_matplotlib['axes.grid.which']

            # Instruction is to turn off gridlines
            if not value:
                # Gridlines are already off, or they are on for the particular
                # ones that we want to turn off. Instruct to turn both off.
                if (
                    not ovalue
                    or key == 'grid' and owhich == 'major'
                    or key == 'gridminor' and owhich == 'minor'
                ):
                    which = 'both'  # disable both sides
                # Gridlines are currently on for major and minor ticks, so we
                # instruct to turn on gridlines for the one we *don't* want off
                elif owhich == 'both':  # and ovalue is True, as already tested
                    # if gridminor=False, enable major, and vice versa
                    b = True
                    which = 'major' if key == 'gridminor' else 'minor'
                # Gridlines are on for the ones that we *didn't* instruct to
                # turn off, and off for the ones we do want to turn off. This
                # just re-asserts the ones that are already on.
                else:
                    b = True
                    which = owhich

            # Instruction is to turn on gridlines
            else:
                # Gridlines are already both on, or they are off only for the
                # ones that we want to turn on. Turn on gridlines for both.
                if (
                    owhich == 'both'
                    or key == 'grid' and owhich == 'minor'
                    or key == 'gridminor' and owhich == 'major'
                ):
                    which = 'both'
                # Gridlines are off for both, or off for the ones that we
                # don't want to turn on. We can just turn on these ones.
                else:
                    which = owhich

            # Finally apply settings
            kw_matplotlib['axes.grid'] = b
            kw_matplotlib['axes.grid.which'] = which

        # Update original setting and linked settings
        for key in keys:
            if key in rc_proplot:
                kw_proplot[key] = value
            elif key in rc_matplotlib:
                kw_matplotlib[key] = value
            else:
                raise KeyError(f'Invalid rc key {key!r}.')
        return kw_proplot, kw_matplotlib

    @staticmethod
    def _get_user_path():
        """
        Return location of user proplotrc file.
        """
        return os.path.join(os.path.expanduser('~'), '.proplotrc')

    @staticmethod
    def _get_local_paths():
        """
        Return locations of local proplotrc files in this directory
        and in parent directories.
        """
        cdir = os.getcwd()
        paths = []
        # Loop until we reach root
        while cdir:
            # Look for hidden and unhidden proplotrc files
            path = os.path.join(cdir, 'proplotrc')
            if os.path.exists(path):
                paths.append(path)
            path = os.path.join(cdir, '.proplotrc')
            if os.path.exists(path):
                paths.append(path)
            # Move on to next parent directory
            ndir = os.path.dirname(cdir)
            if ndir == cdir:  # root
                break
            cdir = ndir
        return paths[::-1]  # sort from decreasing to increasing importantce

    @staticmethod
    def _sanitize_key(key):
        """
        Ensure string and convert keys with omitted dots.
        """
        if not isinstance(key, str):
            raise KeyError(f'Invalid key {key!r}. Must be string.')

        # Translate from nodots to 'full' version
        if '.' not in key:
            key = rcsetup._rc_nodots.get(key, key)

        # Handle deprecations
        if key in rcsetup._rc_removed:
            alternative, version = rcsetup._rc_removed[key]
            message = f'rc setting {key!r} was removed in version {version}.'
            if alternative:  # provide an alternative
                message = f'{message} {alternative}'
            warnings._warn_proplot(warnings)
            key = None
        if key in rcsetup._rc_renamed:
            key_new, version = rcsetup._rc_renamed[key]
            warnings._warn_proplot(
                f'rc setting {key!r} was renamed to {key_new} in version {version}.'
            )
            key = key_new

        return key.lower()

    @staticmethod
    def _sanitize_value(value):
        """
        Convert numpy ndarray to list.
        """
        if isinstance(value, np.ndarray):
            if value.size <= 1:
                value = value.item()
            else:
                value = value.tolist()
        return value

    @staticmethod
    def _scale_font(size):
        """
        Translate font size to numeric.
        """
        # NOTE: Critical this remains KeyError so except clause
        # in _get_synced_params works.
        if not isinstance(size, str):
            return size
        try:
            scale = mfonts.font_scalings[size]
        except KeyError:
            raise KeyError(
                f'Invalid font scaling {size!r}. Options are: '
                + ', '.join(
                    f'{key!r} ({value})'
                    for key, value in mfonts.font_scalings.items()
                ) + '.'
            )
        else:
            return rc_matplotlib['font.size'] * scale

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
            See `~RcConfigurator.context`.

        See also
        --------
        RcConfigurator.get
        RcConfigurator.fill
        """
        if cat not in rcsetup._rc_categories:
            raise ValueError(
                f'Invalid rc category {cat!r}. Valid categories are '
                + ', '.join(map(repr, rcsetup._rc_categories)) + '.'
            )
        kw = {}
        mode = 0 if not context else None
        for rcdict in (rc_proplot, rc_matplotlib):
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

    def context(self, *args, mode=0, file=None, **kwargs):
        """
        Temporarily modify the rc settings in a "with as" block.

        Parameters
        ----------
        *args
            Dictionaries of `rc` names and values.
        file : str, optional
            Filename from which settings should be loaded.
        **kwargs
            `rc` names and values passed as keyword arguments. If the
            name has dots, simply omit them.

        Other parameters
        ----------------
        mode : {0, 1, 2}, optional
            The context mode. Dictates the behavior of `~RcConfigurator.get`,
            `~RcConfigurator.fill`, and `~RcConfigurator.category` within a
            "with as" block when called with ``context=True``.

            The options are as follows:

            0. Matplotlib's `builtin settings <rc_matplotlib>`_ and ProPlot's
               :ref:`added settings <rc_proplot>` are all returned,
               whether or not `~RcConfigurator.context` has changed them.
            1. *Unchanged* `matplotlib settings <rc_matplotlib>`_ return ``None``.
               All of ProPlot's :ref:`added settings <rc_proplot>` are returned
               whether or not `~RcConfigurator.context` has changed them.
               This is used in the `~proplot.axes.Axes.__init__` call to
               `~proplot.axes.Axes.format`. When a lookup returns ``None``,
               `~proplot.axes.Axes.format` does not apply it.
            2. All unchanged settings return ``None``. This is used during
               user calls to `~proplot.axes.Axes.format`.

        Note
        ----
        This is used by ProPlot internally but may also be useful for power users.
        It was invented to prevent successive calls to `~proplot.axes.Axes.format`
        from constantly looking up and re-applying unchanged settings. These
        gratuitous lookups increased runtime significantly, and resulted in successive
        calls to `~proplot.axes.Axes.format` overwriting the previous calls.

        Example
        -------
        The below applies settings to axes in a specific figure using
        `~RcConfigurator.context`.

        >>> import proplot as pplt
        >>> with pplt.rc.context(linewidth=2, ticklen=5):
        >>>     fig, ax = pplt.subplots()
        >>>     ax.plot(data)

        The below applies settings to a specific axes using `~proplot.axes.Axes.format`,
        which uses `~RcConfigurator.context` internally.

        >>> import proplot as pplt
        >>> fig, ax = pplt.subplots()
        >>> ax.format(linewidth=2, ticklen=5)
        """
        # Add input dictionaries
        for arg in args:
            if not isinstance(arg, dict):
                raise ValueError('Non-dictionary argument {arg!r}.')
            kwargs.update(arg)

        # Add settings from file
        # TODO: Incoporate with matplotlib 'stylesheets'
        if file is not None:
            kw_proplot, kw_matplotlib = self._load_file(file)
            kwargs.update(kw_proplot)
            kwargs.update(kw_matplotlib)

        # Activate context object
        if mode not in range(3):
            raise ValueError(f'Invalid mode {mode!r}.')
        context = _RcContext(mode=mode, kwargs=kwargs, rc_new={}, rc_old={})
        self._context.append(context)
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
            in the context mode dictionaries. See `~RcConfigurator.context`.

        See also
        --------
        RcConfigurator.category
        RcConfigurator.fill
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
            See `~RcConfigurator.context`.

        See also
        --------
        RcConfigurator.category
        RcConfigurator.get
        """
        kw = {}
        mode = 0 if not context else None
        for key, value in props.items():
            item = self._get_item(value, mode)
            if item is not None:
                kw[key] = item
        return kw

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

        See also
        --------
        RcConfigurator.fill
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
        if default:
            rc_matplotlib.update(_get_style_dicts('original', filter=False))
            rc_matplotlib.update(rcsetup._rc_matplotlib_default)
            rc_proplot.update(rcsetup._rc_proplot_default)
            for key, value in rc_proplot.items():
                kw_proplot, kw_matplotlib = self._get_synced_params(key, value)
                rc_matplotlib.update(kw_matplotlib)
                rc_proplot.update(kw_proplot)

        # Update from user home
        user_path = None
        if user:
            user_path = self._get_user_path()
            if os.path.isfile(user_path):
                self.load_file(user_path)

        # Update from local paths
        if local:
            local_paths = self._get_local_paths()
            for path in local_paths:
                if path == user_path:  # local files always have precedence
                    continue
                self.load_file(path)

    def _load_file(self, path):
        """
        Return dictionaries of proplot and matplotlib settings loaded from the file.
        """
        added = set()
        path = os.path.expanduser(path)
        kw_proplot = {}
        kw_matplotlib = {}
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

                # *Very primitive* type conversion system for proplot settings.
                # Matplotlib does this automatically.
                if REGEX_STRING.match(val):  # also do this for matplotlib settings
                    val = val[1:-1]  # remove quotes from string
                if key in rc_proplot:
                    if not val or val == 'None':
                        val = None  # older proplot versions supported this
                    elif val == 'True':
                        val = True
                    elif val == 'False':
                        val = False
                    else:
                        try:
                            val = float(val) if '.' in val else int(val)
                        except ValueError:
                            pass

                # Add to dictionaries
                try:
                    ikw_proplot, ikw_matplotlib = self._get_synced_params(key, val)
                    kw_proplot.update(ikw_proplot)
                    kw_matplotlib.update(ikw_matplotlib)
                except KeyError:
                    warnings._warn_proplot(
                        f'Invalid key {key!r} on line #{cnt} in file {path!r}.'
                    )

        return kw_proplot, kw_matplotlib

    def load_file(self, path):
        """
        Load settings from the specified file.

        Parameters
        ----------
        path : str
            The file path.

        See also
        --------
        RcConfigurator.save
        """
        kw_proplot, kw_matplotlib = self._load_file(path)
        rc_proplot.update(kw_proplot)
        rc_matplotlib.update(kw_matplotlib)

    @staticmethod
    def _save_rst(path):
        """
        Used internally to create table for online docs.
        """
        string = rcsetup._gen_rst_table()
        with open(path, 'w') as fh:
            fh.write(string)

    @staticmethod
    def _save_proplotrc(path, comment=False):
        """
        Used internally to create initial proplotrc file and file for online docs.
        """
        self = object()  # self is unused when 'user' is False
        RcConfigurator.save(self, path, user=False, backup=False, comment=comment)

    def save(self, path=None, user=True, comment=None, backup=True, description=False):
        """
        Save the current settings to a ``.proplotrc`` file. This writes
        the default values commented out plus the values that *differ*
        from the defaults at the top of the file.

        Parameters
        ----------
        path : str, optional
            The path name. The default file name is ``.proplotrc`` and the default
            directory is the home directory. Use ``path=''`` to save to the current
            directory.
        user : bool, optional
            If ``True`` (the default), the settings you changed since importing
            proplot are shown uncommented at the very top of the file.
        backup : bool, optional
            If the file already exists and this is set to ``True``, it is moved
            to a backup file with the suffix ``.bak``.
        comment : bool, optional
            Whether to comment out the default settings. Default is the
            value of `user`.
        description : bool, optional
            Whether to include descriptions of each setting as comments.
            Default is ``False``.

        See also
        --------
        RcConfigurator.load_file
        """
        if path is None:
            path = '~'
        path = os.path.abspath(os.path.expanduser(path))
        if os.path.isdir(path):
            path = os.path.join(path, '.proplotrc')
        if os.path.isfile(path) and backup:
            os.rename(path, path + '.bak')
            warnings._warn_proplot(
                f'Existing proplotrc file {path!r} was moved to {path + ".bak"!r}.'
            )

        # Generate user-specific table, ignoring non-style related
        # settings that may be changed from defaults like 'backend'
        rc_user = ()
        if user:
            # Changed settings
            rcdict = {
                key: value for key, value in self
                if value != rcsetup._get_default_param(key)
            }

            # Special handling for certain settings
            # TODO: For now not sure how to detect if prop cycle changed since
            # we cannot load it from _cmap_database in rcsetup.
            rcdict.pop('interactive', None)  # changed by backend
            rcdict.pop('axes.prop_cycle', None)

            # Filter and get table
            rcdict = _get_filtered_dict(rcdict, warn=False)
            rc_user_table = rcsetup._gen_yaml_table(rcdict, comment=False)
            rc_user = ('# Settings changed by user', rc_user_table, '')  # + blank line

        # Generate tables and write
        comment = _not_none(comment, user)
        rc_proplot_table = rcsetup._gen_yaml_table(
            rcsetup._rc_proplot, comment=comment, description=description,
        )
        rc_matplotlib_table = rcsetup._gen_yaml_table(
            rcsetup._rc_matplotlib_default, comment=comment
        )
        with open(path, 'w') as fh:
            fh.write('\n'.join((
                '#--------------------------------------------------------------------',
                '# Use this file to change the default proplot and matplotlib settings',
                '# The syntax is identical to matplotlibrc syntax. For details see:',
                '# https://proplot.readthedocs.io/en/latest/configuration.html',
                '# https://matplotlib.org/stable/tutorials/introductory/customizing.html',  # noqa: E501
                '#--------------------------------------------------------------------',
                *rc_user,  # includes blank line
                '# ProPlot settings',
                rc_proplot_table,
                '\n# Matplotlib settings',
                rc_matplotlib_table,
            )))

    def items(self):
        """
        Return an iterator that loops over all setting names and values.
        Same as `dict.items`.

        See also
        --------
        RcConfigurator.keys
        RcConfigurator.values
        RcConfigurator.items
        """
        for key in self:
            yield key, self[key]

    def keys(self):
        """
        Return an iterator that loops over all setting names.
        Same as `dict.keys`.

        See also
        --------
        RcConfigurator.values
        RcConfigurator.items
        """
        for key in self:
            yield key

    def values(self):
        """
        Return an iterator that loops over all setting values.
        Same as `dict.values`.

        See also
        --------
        RcConfigurator.keys
        RcConfigurator.items
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

    See also
    --------
    RcConfigurator
    """
    # Note if inline backend is unavailable this will fail silently
    ipython = get_ipython()
    if ipython is None:
        return
    fmt = _not_none(fmt, rc_proplot['inlinefmt'])
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
    with _suppress_matplotlib_deprecation_warning():
        rcdict = dict(RcParams(rcdict))
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


def _get_style_dicts(style, infer=False, filter=True):
    """
    Return a dictionary of settings belonging to the requested style(s). If `infer`
    is ``True``, two dictionaries are returned, where the second contains custom
    ProPlot settings "inferred" from the matplotlib settings. If `filter` is ``True``,
    invalid style parameters like `backend` are filtered out.
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
    kw_matplotlib = _get_default_dict()
    if style == 'default' or style is mpl.rcParamsDefault:
        if infer:
            kw_proplot = _infer_added_params(kw_matplotlib)
            return kw_matplotlib, kw_proplot
        else:
            return kw_matplotlib

    # Apply limited deviations from the matplotlib style that we want to propagate to
    # other styles. Want users selecting stylesheets to have few surprises, so
    # currently just enforce the new aesthetically pleasing fonts.
    kw_matplotlib['font.family'] = 'sans-serif'
    for fmly in ('serif', 'sans-serif', 'monospace', 'cursive', 'fantasy'):
        kw_matplotlib['font.' + fmly] = rcsetup._rc_matplotlib_default['font.' + fmly]

    # Apply user input style(s) one by one
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
        if filter:
            kw = _get_filtered_dict(kw, warn=True)
        kw_matplotlib.update(kw)

    # Infer proplot params from stylesheet params
    if infer:
        kw_proplot = _infer_added_params(kw_matplotlib)
        return kw_matplotlib, kw_proplot
    else:
        return kw_matplotlib


def _infer_added_params(kw_params):
    """
    Infer values for proplot's "added" parameters from stylesheets.
    """
    kw_proplot = {}
    mpl_to_proplot = {
        'font.size': ('tick.labelsize',),
        'axes.titlesize': (
            'abc.size', 'suptitle.size', 'title.size',
            'leftlabel.size', 'rightlabel.size',
            'toplabel.size', 'bottomlabel.size',
        ),
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
                kw_proplot[param] = value
    return kw_proplot


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

    See also
    --------
    RcConfigurator
    matplotlib.style.use
    """
    # NOTE: This function is not really necessary but makes proplot's
    # stylesheet-supporting features obvious. Plus changing the style does
    # so much *more* than changing rc params or quick settings, so it is
    # nice to have dedicated function instead of just another rc_param name.
    kw_matplotlib, kw_proplot = _get_style_dicts(style, infer=True)
    rc_matplotlib.update(kw_matplotlib)
    rc_proplot.update(kw_proplot)


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

    See also
    --------
    register_cycles
    register_colors
    register_fonts
    proplot.demos.show_cmaps
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

    See also
    --------
    register_cmaps
    register_colors
    register_fonts
    proplot.demos.show_cycles
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

    See also
    --------
    register_cmaps
    register_cycles
    register_fonts
    proplot.demos.show_colors
    """
    # Reset native colors dictionary
    mcolors.colorConverter.colors.clear()  # clean out!
    mcolors.colorConverter.cache.clear()  # clean out!

    # Add in base colors and CSS4 colors so user has no surprises
    for name, dict_ in (('base', COLORS_BASE), ('css', mcolors.CSS4_COLORS)):
        mcolors.colorConverter.colors.update(dict_)

    # Load colors from file and get their HCL values
    # NOTE: Colors that come *later* overwrite colors that come earlier.
    hex = re.compile(rf'\A{pcolors.REGEX_HEX}\Z')  # match each string
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
                if i == 0 and name in COLORS_BASE:
                    continue
                loaded[name] = color

        # Add every user color and every opencolor color and ensure XKCD
        # colors are "perceptually distinct".
        if i == 1:
            mcolors.colorConverter.colors.update(loaded)
        elif cat == 'opencolor':
            mcolors.colorConverter.colors.update(loaded)
            COLORS_OPEN.update(loaded)
        elif cat == 'xkcd':
            # Always add these colors, but make sure not to add other
            # colors too close to them.
            hcls = []
            filtered = []
            for name in COLORS_ADD:
                color = loaded.pop(name, None)
                if color is None:
                    continue
                if 'grey' in name:
                    name = name.replace('grey', 'gray')
                hcls.append(to_xyz(color, space=space))
                filtered.append((name, color))
                mcolors.colorConverter.colors[name] = color
                COLORS_XKCD[name] = color

            # Get locations of "perceptually distinct" colors
            # WARNING: Unique axis argument requires numpy version >=1.13
            for name, color in loaded.items():
                for string, replace in COLORS_TRANSLATE:
                    if string in name:
                        name = name.replace(string, replace)
                if any(string in name for string in COLORS_REMOVE):
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
                COLORS_XKCD[name] = color
        else:
            raise ValueError(f'Unknown proplot color database {path!r}.')


def _patch_validators():
    """
    Fix the fontsize validators to allow for new font scalings.
    """
    # First define valdiators
    # NOTE: In the future will subclass RcParams directly and control the
    # validators ourselves.
    def _validate_fontsize(s):
        fontsizes = list(mfonts.font_scalings)
        if isinstance(s, str):
            s = s.lower()
        if s in fontsizes:
            return s
        try:
            return float(s)
        except ValueError:
            raise ValueError(
                f'{s!r} is not a valid font size. Valid sizes are: '
                + ', '.join(map(repr, fontsizes))
            )

    def _validate_fontsize_None(s):
        if s is None or s == 'None':
            return None
        else:
            return _validate_fontsize(s)

    _validate_fontsizelist = None
    if hasattr(msetup, '_listify_validator'):
        _validate_fontsizelist = msetup._listify_validator(_validate_fontsize)

    # Apply new functions
    validate = RcParams.validate
    for key in list(validate):  # modify in-place
        validator = validate[key]
        if validator is msetup.validate_fontsize:
            validate[key] = _validate_fontsize
        elif validator is getattr(msetup, 'validate_fontsize_None', None):
            validate[key] = _validate_fontsize_None
        elif validator is getattr(msetup, 'validate_fontsizelist', None):
            if _validate_fontsizelist is not None:
                validate[key] = _validate_fontsizelist


def register_fonts():
    """
    Add fonts packaged with ProPlot or saved to the ``~/.proplot/fonts``
    folder, if they are not already added. Detects ``.ttf`` and ``.otf`` files
    -- see `this link \
<https://gree2.github.io/python/2015/04/27/python-change-matplotlib-font-on-mac>`__
    for a guide on converting various other font file types to ``.ttf`` and
    ``.otf`` for use with matplotlib.

    To visualize the registered fonts, use `~proplot.demos.show_fonts`.

    See also
    --------
    register_cmaps
    register_cycles
    register_colors
    proplot.demos.show_fonts
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
            # New API lets us add font files manually and deprecates TTFPATH. However,
            # to cache fonts added this way, we must call json_dump explicitly.
            # NOTE: Previously, cache filename was specified as _fmcache variable, but
            # recently became inaccessible. Must reproduce mpl code instead! Annoying.
            # NOTE: Older mpl versions used fontList.json as the cache, but these
            # versions also did not have 'addfont', so makes no difference.
            for fname in fnames_proplot:
                mfonts.fontManager.addfont(fname)
            cache = os.path.join(
                mpl.get_cachedir(),
                f'fontlist-v{mfonts.FontManager.__version__}.json'
            )
            mfonts.json_dump(mfonts.fontManager, cache)
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
    RcConfigurator._save_proplotrc(_user_rc_file, comment=True)

# Initialize customization folders
_rc_folder = os.path.join(os.path.expanduser('~'), '.proplot')
if not os.path.isdir(_rc_folder):
    os.mkdir(_rc_folder)
for _rc_sub in ('cmaps', 'cycles', 'colors', 'fonts'):
    _rc_sub = os.path.join(_rc_folder, _rc_sub)
    if not os.path.isdir(_rc_sub):
        os.mkdir(_rc_sub)

# Add custom font scalings to font_manager and monkey patch rcParams validator
# NOTE: This is because we prefer large sizes
if hasattr(mfonts, 'font_scalings'):
    mfonts.font_scalings['med-small'] = 0.9
    mfonts.font_scalings['med-large'] = 1.1
    _patch_validators()

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
    _ = RcConfigurator()

#: Instance of `RcConfigurator`. This is used to change global settings.
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
inline_backend_fmt = warnings._rename_objs(
    '0.6', inline_backend_fmt=config_inline_backend
)
