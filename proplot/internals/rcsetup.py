#!/usr/bin/env python3
"""
Utilities for global configuration.
"""
import functools
import re
from collections.abc import MutableMapping
from numbers import Integral, Real

import matplotlib.rcsetup as msetup
import numpy as np
from cycler import Cycler
from matplotlib import RcParams
from matplotlib import rcParamsDefault as _rc_matplotlib_native
from matplotlib.colors import Colormap
from matplotlib.font_manager import font_scalings
from matplotlib.fontconfig_pattern import parse_fontconfig_pattern

from . import ic  # noqa: F401
from . import warnings
from .versions import _version_mpl

# Regex for "probable" unregistered named colors. Try to retain warning message for
# colors that were most likely a failed literal string evaluation during startup.
REGEX_NAMED_COLOR = re.compile(r'\A[a-zA-Z0-9:_ -]*\Z')

# Configurable validation settings
# NOTE: These are set to True inside __init__.py
# NOTE: We really cannot delay creation of 'rc' until after registration because
# colormap creation depends on rc['cmap.lut'] and rc['cmap.listedthresh'].
# And anyway to revoke that dependence would require other uglier kludges.
VALIDATE_REGISTERED_CMAPS = False
VALIDATE_REGISTERED_COLORS = False

# Initial synced properties
# NOTE: Important that LINEWIDTH is less than matplotlib default of 0.8.
# In general want axes lines to look about as thick as text.
# NOTE: Important that default values are equivalent to the *validated* values
# used in the RcParams dictionaries. Otherwise _user_settings() detects changes.
# NOTE: We *could* just leave some settings empty and leave it up to Configurator
# to sync them when proplot is imported... but also sync them here so that we can
# simply compare any Configurator state to these dictionaries and use save() to
# save only the settings changed by the user.
BLACK = 'black'
CYCLE = 'colorblind'
CMAPCYC = 'twilight'
CMAPDIV = 'BuRd'
CMAPSEQ = 'Fire'
CMAPCAT = 'colorblind10'
DIVERGING = 'div'
FRAMEALPHA = 0.8  # legend and colorbar
FONTNAME = 'sans-serif'
FONTSIZE = 9.0
GRIDALPHA = 0.1
GRIDBELOW = 'line'
GRIDPAD = 3.0
GRIDRATIO = 0.5  # differentiated from major by half size reduction
GRIDSTYLE = '-'
LABELPAD = 4.0  # default is 4.0, previously was 3.0
LARGESIZE = 'med-large'
LINEWIDTH = 0.6
MARGIN = 0.05
MATHTEXT = False
SMALLSIZE = 'medium'
TICKDIR = 'out'
TICKLEN = 4.0
TICKLENRATIO = 0.5  # very noticeable length reduction
TICKMINOR = True
TICKPAD = 2.0
TICKWIDTHRATIO = 0.8  # very slight width reduction
TITLEPAD = 5.0  # default is 6.0, previously was 3.0
WHITE = 'white'
ZLINES = 2  # default zorder for lines
ZPATCHES = 1

# Preset legend locations and aliases
LEGEND_LOCS = {
    'fill': 'fill',
    'inset': 'best',
    'i': 'best',
    0: 'best',  # offsetbox codes
    1: 'upper right',
    2: 'upper left',
    3: 'lower left',
    4: 'lower right',
    5: 'center right',  # weird, but see _get_anchored_bbox in matplotlib/offsetbox.py
    6: 'center left',
    7: 'center right',
    8: 'lower center',
    9: 'upper center',
    10: 'center',
    'l': 'left',  # custom shorthands
    'r': 'right',
    'b': 'bottom',
    't': 'top',
    'c': 'center',
    'ur': 'upper right',
    'ul': 'upper left',
    'll': 'lower left',
    'lr': 'lower right',
    'cl': 'center left',
    'cr': 'center right',
    'uc': 'upper center',
    'lc': 'lower center',
    'ne': 'upper right',  # BboxBase.anchored shorthands
    'nw': 'upper left',
    'sw': 'lower left',
    'se': 'lower right',
    'w': 'center left',
    'e': 'center right',
    's': 'lower center',
    'n': 'upper center',
}
for _loc in set(LEGEND_LOCS.values()):
    LEGEND_LOCS[_loc] = _loc  # identity assignments
TEXT_LOCS = {
    key: val for key, val in LEGEND_LOCS.items() if val in (
        'left', 'center', 'right',
        'upper left', 'upper center', 'upper right',
        'lower left', 'lower center', 'lower right',
    )
}
COLORBAR_LOCS = {
    key: val for key, val in LEGEND_LOCS.items() if val in (
        'fill', 'best',
        'left', 'right', 'top', 'bottom',
        'upper left', 'upper right', 'lower left', 'lower right',
    )
}
PANEL_LOCS = {
    key: val for key, val in LEGEND_LOCS.items() if val in (
        'left', 'right', 'top', 'bottom'
    )
}
ALIGN_LOCS = {
    key: val for key, val in LEGEND_LOCS.items() if isinstance(key, str) and val in (
        'left', 'right', 'top', 'bottom', 'center',
    )
}

# Matplotlib setting categories
EM_KEYS = (  # em-width units
    'legend.borderpad',
    'legend.labelspacing',
    'legend.handlelength',
    'legend.handleheight',
    'legend.handletextpad',
    'legend.borderaxespad',
    'legend.columnspacing',
)
PT_KEYS = (
    'font.size',  # special case
    'xtick.major.size',
    'xtick.minor.size',
    'ytick.major.size',
    'ytick.minor.size',
    'xtick.major.pad',
    'xtick.minor.pad',
    'ytick.major.pad',
    'ytick.minor.pad',
    'xtick.major.width',
    'xtick.minor.width',
    'ytick.major.width',
    'ytick.minor.width',
    'axes.labelpad',
    'axes.titlepad',
    'axes.linewidth',
    'grid.linewidth',
    'patch.linewidth',
    'hatch.linewidth',
    'lines.linewidth',
    'contour.linewidth',
)
FONT_KEYS = set()  # dynamically add to this below


def _get_default_param(key):
    """
    Get the default parameter from one of three places. This is used for
    the :rc: role when compiling docs and when saving proplotrc files.
    """
    sentinel = object()
    for dict_ in (
        _rc_proplot_default,
        _rc_matplotlib_default,  # imposed defaults
        _rc_matplotlib_native,  # native defaults
    ):
        value = dict_.get(key, sentinel)
        if value is not sentinel:
            return value
    raise KeyError(f'Invalid key {key!r}.')


def _validate_abc(value):
    """
    Validate a-b-c setting.
    """
    try:
        if np.iterable(value):
            return all(map(_validate_bool, value))
        else:
            return _validate_bool(value)
    except ValueError:
        pass
    if isinstance(value, str):
        if 'a' in value.lower():
            return value
    else:
        if all(isinstance(_, str) for _ in value):
            return tuple(value)
    raise ValueError(
        "A-b-c setting must be string containing 'a' or 'A' or sequence of strings."
    )


def _validate_belongs(*options):
    """
    Return a validator ensuring the item belongs in the list.
    """
    def _validate_belongs(value):  # noqa: E306
        for opt in options:
            if isinstance(value, str) and isinstance(opt, str):
                if value.lower() == opt.lower():  # noqa: E501
                    return opt
            elif value is True or value is False or value is None:
                if value is opt:
                    return opt
            elif value == opt:
                return opt
        raise ValueError(
            f'Invalid value {value!r}. Options are: '
            + ', '.join(map(repr, options))
            + '.'
        )
    return _validate_belongs


def _validate_cmap(subtype):
    """
    Validate the colormap or cycle. Possibly skip name registration check
    and assign the colormap name rather than a colormap instance.
    """
    def _validate_cmap(value):
        name = value
        if isinstance(value, str):
            if VALIDATE_REGISTERED_CMAPS:
                from ..colors import _get_cmap_subtype
                _get_cmap_subtype(name, subtype)  # may trigger useful error message
            return name
        elif isinstance(value, Colormap):
            name = getattr(value, 'name', None)
            if isinstance(name, str):
                from ..colors import _cmap_database  # avoid circular imports
                _cmap_database[name] = value
                return name
        raise ValueError(f'Invalid colormap or color cycle name {name!r}.')
    return _validate_cmap


def _validate_color(value, alternative=None):
    """
    Validate the color. Possibly skip name registration check.
    """
    if alternative and isinstance(value, str) and value.lower() == alternative:
        return value
    try:
        return msetup.validate_color(value)
    except ValueError:
        if (
            VALIDATE_REGISTERED_COLORS
            or not isinstance(value, str)
            or not REGEX_NAMED_COLOR.match(value)
        ):
            raise ValueError(f'{value!r} is not a valid color arg.') from None
        return value
    except Exception as error:
        raise error


def _validate_fontprops(s):
    """
    Parse font property with support for ``'regular'`` placeholder.
    """
    b = s.startswith('regular')
    if b:
        s = s.replace('regular', 'sans', 1)
    parse_fontconfig_pattern(s)
    if b:
        s = s.replace('sans', 'regular', 1)
    return s


def _validate_fontsize(value):
    """
    Validate font size with new scalings and permitting other units.
    """
    if value is None and None in font_scalings:  # has it always been this way?
        return
    if isinstance(value, str):
        value = value.lower()
        if value in font_scalings:
            return value
    try:
        return _validate_pt(value)  # note None is also a valid font size!
    except ValueError:
        pass
    raise ValueError(
        f'Invalid font size {value!r}. Can be points or one of the '
        'preset scalings: ' + ', '.join(map(repr, font_scalings)) + '.'
    )


def _validate_labels(labels, lon=True):
    """
    Convert labels argument to length-4 boolean array.
    """
    if labels is None:
        return [None] * 4
    which = 'lon' if lon else 'lat'
    if isinstance(labels, str):
        labels = (labels,)
    array = np.atleast_1d(labels).tolist()
    if all(isinstance(_, str) for _ in array):
        bool_ = [False] * 4
        opts = ('left', 'right', 'bottom', 'top')
        for string in array:
            if string in opts:
                string = string[0]
            elif set(string) - set('lrbt'):
                raise ValueError(
                    f'Invalid {which}label string {string!r}. Must be one of '
                    + ', '.join(map(repr, opts))
                    + " or a string of single-letter characters like 'lr'."
                )
            for char in string:
                bool_['lrbt'.index(char)] = True
        array = bool_
    if len(array) == 1:
        array.append(False)  # default is to label bottom or left
    if len(array) == 2:
        if lon:
            array = [False, False, *array]
        else:
            array = [*array, False, False]
    if len(array) != 4 or any(isinstance(_, str) for _ in array):
        raise ValueError(f'Invalid {which}label spec: {labels}.')
    return array


def _validate_or_none(validator):
    """
    Allow none otherwise pass to the input validator.
    """
    @functools.wraps(validator)
    def _validate_or_none(value):
        if value is None:
            return
        if isinstance(value, str) and value.lower() == 'none':
            return
        return validator(value)
    _validate_or_none.__name__ = validator.__name__ + '_or_none'
    return _validate_or_none


def _validate_rotation(value):
    """
    Valid rotation arguments.
    """
    if isinstance(value, str) and value.lower() in ('horizontal', 'vertical'):
        return value
    return _validate_float(value)


def _validate_units(dest):
    """
    Validate the input using the units function.
    """
    def _validate_units(value):
        if isinstance(value, str):
            from ..utils import units  # avoid circular imports
            value = units(value, dest)  # validation happens here
        return _validate_float(value)
    return _validate_units


def _rst_table():
    """
    Return the setting names and descriptions in an RST-style table.
    """
    # Initial stuff
    colspace = 2  # spaces between each column
    descrips = tuple(descrip for (_, _, descrip) in _rc_proplot_table.values())
    keylen = len(max((*_rc_proplot_table, 'Key'), key=len)) + 4  # literal backticks
    vallen = len(max((*descrips, 'Description'), key=len))
    divider = '=' * keylen + ' ' * colspace + '=' * vallen + '\n'
    header = 'Key' + ' ' * (keylen - 3 + colspace) + 'Description\n'

    # Build table
    string = divider + header + divider
    for key, (_, _, descrip) in _rc_proplot_table.items():
        spaces = ' ' * (keylen - (len(key) + 4) + colspace)
        string += f'``{key}``{spaces}{descrip}\n'

    string = string + divider
    return '.. rst-class:: proplot-rctable\n\n' + string.strip()


def _to_string(value):
    """
    Translate setting to a string suitable for saving.
    """
    # NOTE: Never safe hex strings with leading '#'. In both matplotlibrc
    # and proplotrc this will be read as comment character.
    if value is None or isinstance(value, (str, bool, Integral)):
        value = str(value)
        if value[:1] == '#':  # i.e. a HEX string
            value = value[1:]
    elif isinstance(value, Real):
        value = str(round(value, 6))  # truncate decimals
    elif isinstance(value, Cycler):
        value = repr(value)  # special case!
    elif isinstance(value, (list, tuple, np.ndarray)):
        value = ', '.join(map(_to_string, value))  # sexy recursion
    else:
        value = None
    return value


def _yaml_table(rcdict, comment=True, description=False):
    """
    Return the settings as a nicely tabulated YAML-style table.
    """
    prefix = '# ' if comment else ''
    data = []
    for key, args in rcdict.items():
        # Optionally append description
        includes_descrip = isinstance(args, tuple) and len(args) == 3
        if not description:
            descrip = ''
            value = args[0] if includes_descrip else args
        elif includes_descrip:
            value, validator, descrip = args
            descrip = '# ' + descrip  # skip the validator
        else:
            raise ValueError(f'Unexpected input {key}={args!r}.')

        # Translate object to string
        value = _to_string(value)
        if value is not None:
            data.append((key, value, descrip))
        else:
            warnings._warn_proplot(
                f'Failed to write rc setting {key} = {value!r}. Must be None, bool, '
                'string, int, float, a list or tuple thereof, or a property cycler.'
            )

    # Generate string
    string = ''
    keylen = len(max(rcdict, key=len))
    vallen = len(max((tup[1] for tup in data), key=len))
    for key, value, descrip in data:
        space1 = ' ' * (keylen - len(key) + 1)
        space2 = ' ' * (vallen - len(value) + 2) if descrip else ''
        string += f'{prefix}{key}:{space1}{value}{space2}{descrip}\n'

    return string.strip()


class _RcParams(MutableMapping, dict):
    """
    A simple dictionary with locked inputs and validated assignments.
    """
    # NOTE: By omitting __delitem__ in MutableMapping we effectively
    # disable mutability. Also disables deleting items with pop().
    def __init__(self, source, validate):
        self._validate = validate
        for key, value in source.items():
            self.__setitem__(key, value)  # trigger validation

    def __repr__(self):
        return RcParams.__repr__(self)

    def __str__(self):
        return RcParams.__repr__(self)

    def __len__(self):
        return dict.__len__(self)

    def __iter__(self):
        # NOTE: Proplot doesn't add deprecated args to dictionary so
        # we don't have to suppress warning messages here.
        yield from sorted(dict.__iter__(self))

    def __getitem__(self, key):
        key, _ = self._check_key(key)
        return dict.__getitem__(self, key)

    def __setitem__(self, key, value):
        key, value = self._check_key(key, value)
        if key not in self._validate:
            raise KeyError(f'Invalid rc key {key!r}.')
        try:
            value = self._validate[key](value)
        except (ValueError, TypeError) as error:
            raise ValueError(f'Key {key}: {error}') from None
        if key is not None:
            dict.__setitem__(self, key, value)

    @staticmethod
    def _check_key(key, value=None):
        # NOTE: If we assigned from the Configurator then the deprecated key will
        # still propagate to the same 'children' as the new key.
        # NOTE: This also translates values for special cases of renamed keys.
        # Currently the special cases are 'basemap' and 'cartopy.autoextent'.
        if key in _rc_renamed:
            key_new, version = _rc_renamed[key]
            warnings._warn_proplot(
                f'The rc setting {key!r} was deprecated in version {version} and may be '  # noqa: E501
                f'removed in {warnings._next_release()}. Please use {key_new!r} instead.'  # noqa: E501
            )
            if key == 'basemap':  # special case
                value = ('cartopy', 'basemap')[int(bool(value))]
            if key == 'cartopy.autoextent':
                value = ('globe', 'auto')[int(bool(value))]
            key = key_new
        if key in _rc_removed:
            info, version = _rc_removed[key]
            raise KeyError(
                f'The rc setting {key!r} was removed in version {version}.'
                + (info and ' ' + info)
            )
        return key, value

    def copy(self):
        source = {key: dict.__getitem__(self, key) for key in self}
        return _RcParams(source, self._validate)


# Borrow validators from matplotlib and construct some new ones
# WARNING: Instead of validate_fontweight matplotlib used validate_string
# until version 3.1.2. So use that as backup here.
# WARNING: We create custom 'or none' validators since their
# availability seems less consistent across matplotlib versions.
_validate_pt = _validate_units('pt')
_validate_em = _validate_units('em')
_validate_in = _validate_units('in')
_validate_bool = msetup.validate_bool
_validate_int = msetup.validate_int
_validate_float = msetup.validate_float
_validate_string = msetup.validate_string
_validate_fontname = msetup.validate_stringlist  # same as 'font.family'
_validate_fontweight = getattr(msetup, 'validate_fontweight', _validate_string)

# Special style validators
# See: https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.FancyBboxPatch.html
_validate_boxstyle = _validate_belongs(
    'square', 'circle', 'round', 'round4', 'sawtooth', 'roundtooth',
)
if hasattr(msetup, '_validate_linestyle'):  # fancy validation including dashes
    _validate_linestyle = msetup._validate_linestyle
else:  # no dashes allowed then but no big deal
    _validate_linestyle = _validate_belongs(
        '-', ':', '--', '-.', 'solid', 'dashed', 'dashdot', 'dotted', 'none', ' ', '',
    )

# Patch existing matplotlib validators.
# NOTE: validate_fontsizelist is unused in recent matplotlib versions and
# validate_colorlist is only used with prop cycle eval (which we don't care about)
font_scalings['med'] = 1.0  # consistent shorthand
font_scalings['med-small'] = 0.9  # add scaling
font_scalings['med-large'] = 1.1  # add scaling
if not hasattr(RcParams, 'validate'):  # not mission critical so skip
    warnings._warn_proplot('Failed to update matplotlib rcParams validators.')
else:
    _validate = RcParams.validate
    _validate['image.cmap'] = _validate_cmap('continuous')
    _validate['legend.loc'] = _validate_belongs(*LEGEND_LOCS)
    for _key, _validator in _validate.items():
        if _validator is getattr(msetup, 'validate_fontsize', None):  # should exist
            FONT_KEYS.add(_key)
            _validate[_key] = _validate_fontsize
        if _validator is getattr(msetup, 'validate_fontsize_None', None):
            FONT_KEYS.add(_key)
            _validate[_key] = _validate_or_none(_validate_fontsize)
        if _validator is getattr(msetup, 'validate_font_properties', None):
            _validate[_key] = _validate_fontprops
        if _validator is getattr(msetup, 'validate_color', None):  # should exist
            _validate[_key] = _validate_color
        if _validator is getattr(msetup, 'validate_color_or_auto', None):
            _validate[_key] = functools.partial(_validate_color, alternative='auto')
        if _validator is getattr(msetup, 'validate_color_or_inherit', None):
            _validate[_key] = functools.partial(_validate_color, alternative='inherit')
    for _keys, _validator_replace in ((EM_KEYS, _validate_em), (PT_KEYS, _validate_pt)):
        for _key in _keys:
            _validator = _validate.get(_key, None)
            if _validator is None:
                continue
            if _validator is msetup.validate_float:
                _validate[_key] = _validator_replace
            if _validator is getattr(msetup, 'validate_float_or_None'):
                _validate[_key] = _validate_or_none(_validator_replace)


# Proplot overrides of matplotlib default style
# WARNING: Critical to include every parameter here that can be changed by a
# "meta" setting so that _get_default_param returns the value imposed by *proplot*
# and so that "changed" settings detected by Configurator.save are correct.
_rc_matplotlib_default = {
    'axes.axisbelow': GRIDBELOW,
    'axes.formatter.use_mathtext': MATHTEXT,
    'axes.grid': True,  # enable lightweight transparent grid by default
    'axes.grid.which': 'major',
    'axes.edgecolor': BLACK,
    'axes.labelcolor': BLACK,
    'axes.labelpad': LABELPAD,  # more compact
    'axes.labelsize': SMALLSIZE,
    'axes.labelweight': 'normal',
    'axes.linewidth': LINEWIDTH,
    'axes.titlepad': TITLEPAD,  # more compact
    'axes.titlesize': LARGESIZE,
    'axes.titleweight': 'normal',
    'axes.xmargin': MARGIN,
    'axes.ymargin': MARGIN,
    'errorbar.capsize': 3.0,
    'figure.autolayout': False,
    'figure.figsize': (4.0, 4.0),  # for interactife backends
    'figure.dpi': 100,
    'figure.facecolor': '#f4f4f4',  # similar to MATLAB interface
    'figure.titlesize': LARGESIZE,
    'figure.titleweight': 'bold',  # differentiate from axes titles
    'font.serif': [
        'TeX Gyre Schola',  # Century lookalike
        'TeX Gyre Bonum',  # Bookman lookalike
        'TeX Gyre Termes',  # Times New Roman lookalike
        'TeX Gyre Pagella',  # Palatino lookalike
        'DejaVu Serif',
        'Bitstream Vera Serif',
        'Computer Modern Roman',
        'Bookman',
        'Century Schoolbook L',
        'Charter',
        'ITC Bookman',
        'New Century Schoolbook',
        'Nimbus Roman No9 L',
        'Noto Serif',
        'Palatino',
        'Source Serif Pro',
        'Times New Roman',
        'Times',
        'Utopia',
        'serif',
    ],
    'font.sans-serif': [
        'TeX Gyre Heros',  # Helvetica lookalike
        'DejaVu Sans',
        'Bitstream Vera Sans',
        'Computer Modern Sans Serif',
        'Arial',
        'Avenir',
        'Fira Math',
        'Fira Sans',
        'Frutiger',
        'Geneva',
        'Gill Sans',
        'Helvetica',
        'Lucid',
        'Lucida Grande',
        'Myriad Pro',
        'Noto Sans',
        'Roboto',
        'Source Sans Pro',
        'Tahoma',
        'Trebuchet MS',
        'Ubuntu',
        'Univers',
        'Verdana',
        'sans-serif',
    ],
    'font.cursive': [
        'TeX Gyre Chorus',  # Chancery lookalike
        'Apple Chancery',
        'Felipa',
        'Sand',
        'Script MT',
        'Textile',
        'Zapf Chancery',
        'cursive',
    ],
    'font.fantasy': [
        'TeX Gyre Adventor',  # Avant Garde lookalike
        'Avant Garde',
        'Charcoal',
        'Chicago',
        'Comic Sans MS',
        'Futura',
        'Humor Sans',
        'Impact',
        'Optima',
        'Western',
        'xkcd',
        'fantasy',
    ],
    'font.monospace': [
        'TeX Gyre Cursor',  # Courier lookalike
        'DejaVu Sans Mono',
        'Bitstream Vera Sans Mono',
        'Computer Modern Typewriter',
        'Andale Mono',
        'Courier New',
        'Courier',
        'Fixed',
        'Nimbus Mono L',
        'Terminal',
        'monospace',
    ],
    'font.family': FONTNAME,
    'font.size': FONTSIZE,
    'grid.alpha': GRIDALPHA,  # lightweight unobtrusive gridlines
    'grid.color': BLACK,  # lightweight unobtrusive gridlines
    'grid.linestyle': GRIDSTYLE,
    'grid.linewidth': LINEWIDTH,
    'hatch.color': BLACK,
    'hatch.linewidth': LINEWIDTH,
    'image.cmap': CMAPSEQ,
    'lines.linestyle': '-',
    'lines.linewidth': 1.5,
    'lines.markersize': 6.0,
    'legend.borderaxespad': 0,  # i.e. flush against edge
    'legend.borderpad': 0.5,  # a bit more roomy
    'legend.columnspacing': 1.5,  # a bit more compact (see handletextpad)
    'legend.edgecolor': BLACK,
    'legend.facecolor': WHITE,
    'legend.fancybox': False,  # i.e. BboxStyle 'square' not 'round'
    'legend.fontsize': SMALLSIZE,
    'legend.framealpha': FRAMEALPHA,
    'legend.handleheight': 1.0,  # default is 0.7
    'legend.handlelength': 2.0,  # default is 2.0
    'legend.handletextpad': 0.5,  # a bit more compact (see columnspacing)
    'mathtext.default': 'it',
    'mathtext.fontset': 'custom',
    'mathtext.bf': 'regular:bold',  # custom settings implemented above
    'mathtext.cal': 'cursive',
    'mathtext.it': 'regular:italic',
    'mathtext.rm': 'regular',
    'mathtext.sf': 'regular',
    'mathtext.tt': 'monospace',
    'patch.linewidth': LINEWIDTH,
    'savefig.bbox': None,  # do not use 'tight'
    'savefig.directory': '',  # use the working directory
    'savefig.dpi': 1000,  # use academic journal recommendation
    'savefig.facecolor': WHITE,  # use white instead of 'auto'
    'savefig.format': 'pdf',  # use vector graphics
    'savefig.transparent': False,
    'xtick.color': BLACK,
    'xtick.direction': TICKDIR,
    'xtick.labelsize': SMALLSIZE,
    'xtick.major.pad': TICKPAD,
    'xtick.major.size': TICKLEN,
    'xtick.major.width': LINEWIDTH,
    'xtick.minor.pad': TICKPAD,
    'xtick.minor.size': TICKLEN * TICKLENRATIO,
    'xtick.minor.width': LINEWIDTH * TICKWIDTHRATIO,
    'xtick.minor.visible': TICKMINOR,
    'ytick.color': BLACK,
    'ytick.direction': TICKDIR,
    'ytick.labelsize': SMALLSIZE,
    'ytick.major.pad': TICKPAD,
    'ytick.major.size': TICKLEN,
    'ytick.major.width': LINEWIDTH,
    'ytick.minor.pad': TICKPAD,
    'ytick.minor.size': TICKLEN * TICKLENRATIO,
    'ytick.minor.width': LINEWIDTH * TICKWIDTHRATIO,
    'ytick.minor.visible': TICKMINOR,
}
if 'mathtext.fallback' in _rc_matplotlib_native:
    _rc_matplotlib_default['mathtext.fallback'] = 'stixsans'

# Proplot pseudo-setting defaults, validators, and descriptions
# NOTE: Cannot have different a-b-c and title paddings because they are both controlled
# by matplotlib's _title_offset_trans transform and want to keep them aligned anyway.
_addendum_rotation = " Must be 'vertical', 'horizontal', or a float indicating degrees."
_addendum_em = ' Interpreted by `~proplot.utils.units`. Numeric units are em-widths.'
_addendum_in = ' Interpreted by `~proplot.utils.units`. Numeric units are inches.'
_addendum_pt = ' Interpreted by `~proplot.utils.units`. Numeric units are points.'
_addendum_font = (
    ' Must be a :ref:`relative font size <font_table>` or unit string '
    'interpreted by `~proplot.utils.units`. Numeric units are points.'
)
_rc_proplot_table = {
    # Stylesheet
    'style': (
        None,
        _validate_or_none(_validate_string),
        'The default matplotlib `stylesheet '
        '<https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html>`__ '  # noqa: E501
        'name. If ``None``, a custom proplot style is used. '
        "If ``'default'``, the default matplotlib style is used."
    ),

    # A-b-c labels
    'abc': (
        False,
        _validate_abc,
        'If ``False`` then a-b-c labels are disabled. If ``True`` the default label '
        'style ``a`` is used. If string this indicates the style and must contain the '
        "character ``a`` or ``A``, for example ``'a.'`` or ``'(A)'``."
    ),
    'abc.border': (
        True,
        _validate_bool,
        'Whether to draw a white border around a-b-c labels '
        'when :rcraw:`abc.loc` is inside the axes.'
    ),
    'abc.borderwidth': (
        1.5,
        _validate_pt,
        'Width of the white border around a-b-c labels.'
    ),
    'abc.box': (
        False,
        _validate_bool,
        'Whether to draw semi-transparent bounding boxes around a-b-c labels '
        'when :rcraw:`abc.loc` is inside the axes.'
    ),
    'abc.boxcolor': (
        WHITE,
        _validate_color,
        'a-b-c label bounding box color.'
    ),
    'abc.boxalpha': (
        0.5,
        _validate_float,
        'a-b-c label bounding box opacity.'
    ),
    'abc.boxstyle': (
        'square',
        _validate_boxstyle,
        'a-b-c label bounding box style.'
    ),
    'abc.boxpad': (
        None,
        _validate_or_none(_validate_pt),
        'Padding for the a-b-c label bounding box. By default this is scaled '
        'to make the box flush against the subplot edge.' + _addendum_pt
    ),
    'abc.color': (
        BLACK,
        _validate_color,
        'a-b-c label color.'
    ),
    'abc.loc': (
        'left',  # left side above the axes
        _validate_belongs(*TEXT_LOCS),
        'a-b-c label position. '
        'For options see the :ref:`location table <title_table>`.'
    ),
    'abc.size': (
        LARGESIZE,
        _validate_fontsize,
        'a-b-c label font size.' + _addendum_font
    ),
    'abc.titlepad': (
        LABELPAD,
        _validate_pt,
        'Padding separating the title and a-b-c label when in the same location.'
        + _addendum_pt
    ),
    'abc.weight': (
        'bold',
        _validate_fontweight,
        'a-b-c label font weight.'
    ),

    # Autoformatting
    'autoformat': (
        True,
        _validate_bool,
        'Whether to automatically apply labels from `pandas.Series`, '
        '`pandas.DataFrame`, and `xarray.DataArray` objects passed to '
        'plotting functions. See also :rcraw:`unitformat`.'
    ),

    # Axes additions
    'axes.alpha': (
        None,
        _validate_or_none(_validate_float),
        'Opacity of the background axes patch.'
    ),
    'axes.inbounds': (
        True,
        _validate_bool,
        'Whether to exclude out-of-bounds data when determining the default *y* (*x*) '
        'axis limits and the *x* (*y*) axis limits have been locked.'
    ),
    'axes.margin': (
        MARGIN,
        _validate_float,
        'The fractional *x* and *y* axis margins when limits are unset.'
    ),

    # Country borders
    'borders': (
        False,
        _validate_bool,
        'Toggles country border lines on and off.'
    ),
    'borders.alpha': (
        None,
        _validate_or_none(_validate_float),
        'Opacity for country border lines.',
    ),
    'borders.color': (
        BLACK,
        _validate_color,
        'Line color for country border lines.'
    ),
    'borders.linewidth': (
        LINEWIDTH,
        _validate_pt,
        'Line width for country border lines.'
    ),
    'borders.zorder': (
        ZLINES,
        _validate_float,
        'Z-order for country border lines.'
    ),

    # Bottom subplot labels
    'bottomlabel.color': (
        BLACK,
        _validate_color,
        'Font color for column labels on the bottom of the figure.'
    ),
    'bottomlabel.pad': (
        TITLEPAD,
        _validate_pt,
        'Padding between axes content and column labels on the bottom of the figure.'
        + _addendum_pt
    ),
    'bottomlabel.rotation': (
        'horizontal',
        _validate_rotation,
        'Rotation for column labels at the bottom of the figure.' + _addendum_rotation
    ),
    'bottomlabel.size': (
        LARGESIZE,
        _validate_fontsize,
        'Font size for column labels on the bottom of the figure.' + _addendum_font
    ),
    'bottomlabel.weight': (
        'bold',
        _validate_fontweight,
        'Font weight for column labels on the bottom of the figure.'
    ),

    # Coastlines
    'coast': (
        False,
        _validate_bool,
        'Toggles coastline lines on and off.'
    ),
    'coast.alpha': (
        None,
        _validate_or_none(_validate_float),
        'Opacity for coast lines',
    ),
    'coast.color': (
        BLACK,
        _validate_color,
        'Line color for coast lines.'
    ),
    'coast.linewidth': (
        LINEWIDTH,
        _validate_pt,
        'Line width for coast lines.'
    ),
    'coast.zorder': (
        ZLINES,
        _validate_float,
        'Z-order for coast lines.'
    ),

    # Colorbars
    'colorbar.edgecolor': (
        BLACK,
        _validate_color,
        'Color for the inset colorbar frame edge.'
    ),
    'colorbar.extend': (
        1.3,
        _validate_em,
        'Length of rectangular or triangular "extensions" for panel colorbars.'
        + _addendum_em
    ),
    'colorbar.fancybox': (
        False,
        _validate_bool,
        'Whether to use a "fancy" round bounding box for inset colorbar frames.'
    ),
    'colorbar.framealpha': (
        FRAMEALPHA,
        _validate_float,
        'Opacity for inset colorbar frames.'
    ),
    'colorbar.facecolor': (
        WHITE,
        _validate_color,
        'Color for the inset colorbar frame.'
    ),
    'colorbar.frameon': (
        True,
        _validate_bool,
        'Whether to draw a frame behind inset colorbars.'
    ),
    'colorbar.grid': (
        False,
        _validate_bool,
        'Whether to draw borders between each level of the colorbar.'
    ),
    'colorbar.insetextend': (
        0.9,
        _validate_em,
        'Length of rectangular or triangular "extensions" for inset colorbars.'
        + _addendum_em
    ),
    'colorbar.insetlength': (
        8,
        _validate_em,
        'Length of inset colorbars.' + _addendum_em
    ),
    'colorbar.insetpad': (
        0.7,
        _validate_em,
        'Padding between axes edge and inset colorbars.' + _addendum_em
    ),
    'colorbar.insetwidth': (
        1.2,
        _validate_em,
        'Width of inset colorbars.' + _addendum_em
    ),
    'colorbar.length': (
        1,
        _validate_em,
        'Length of outer colorbars.'
    ),
    'colorbar.loc': (
        'right',
        _validate_belongs(*COLORBAR_LOCS),
        'Inset colorbar location. '
        'For options see the :ref:`location table <colorbar_table>`.'
    ),
    'colorbar.width': (
        0.2,
        _validate_in,
        'Width of outer colorbars.' + _addendum_in
    ),
    'colorbar.rasterized': (
        False,
        _validate_bool,
        'Whether to use rasterization for colorbar solids.'
    ),
    'colorbar.shadow': (
        False,
        _validate_bool,
        'Whether to add a shadow underneath inset colorbar frames.'
    ),

    # Color cycle additions
    'cycle': (
        CYCLE,
        _validate_cmap('discrete'),
        'Name of the color cycle assigned to :rcraw:`axes.prop_cycle`.'
    ),

    # Colormap additions
    'cmap': (
        CMAPSEQ,
        _validate_cmap('continuous'),
        'Alias for :rcraw:`cmap.sequential` and :rcraw:`image.cmap`.'
    ),
    'cmap.autodiverging': (
        True,
        _validate_bool,
        'Whether to automatically apply a diverging colormap and '
        'normalizer based on the data.'
    ),
    'cmap.qualitative': (
        CMAPCAT,
        _validate_cmap('discrete'),
        'Default colormap for qualitative datasets.'
    ),
    'cmap.cyclic': (
        CMAPCYC,
        _validate_cmap('continuous'),
        'Default colormap for cyclic datasets.'
    ),
    'cmap.discrete': (
        None,
        _validate_or_none(_validate_bool),
        'If ``True``, `~proplot.colors.DiscreteNorm` is used for every colormap plot. '
        'If ``False``, it is never used. If ``None``, it is used for all plot types '
        'except `imshow`, `matshow`, `spy`, `hexbin`, and `hist2d`.'
    ),
    'cmap.diverging': (
        CMAPDIV,
        _validate_cmap('continuous'),
        'Default colormap for diverging datasets.'
    ),
    'cmap.inbounds': (
        True,
        _validate_bool,
        'If ``True`` and the *x* and *y* axis limits are fixed, only in-bounds data '
        'is considered when determining the default colormap `vmin` and `vmax`.'
    ),
    'cmap.levels': (
        11,
        _validate_int,
        'Default number of `~proplot.colors.DiscreteNorm` levels for plotting '
        'commands that use colormaps.'
    ),
    'cmap.listedthresh': (
        64,
        _validate_int,
        'Native `~matplotlib.colors.ListedColormap`\\ s with more colors than '
        'this are converted to `~proplot.colors.ContinuousColormap` rather than '
        '`~proplot.colors.DiscreteColormap`. This helps translate continuous '
        'colormaps from external projects.'
    ),
    'cmap.lut': (
        256,
        _validate_int,
        'Number of colors in the colormap lookup table. '
        'Alias for :rcraw:`image.lut`.'
    ),
    'cmap.robust': (
        False,
        _validate_bool,
        'If ``True``, the default colormap `vmin` and `vmax` are chosen using the '
        '2nd to 98th percentiles rather than the minimum and maximum.'
    ),
    'cmap.sequential': (
        CMAPSEQ,
        _validate_cmap('continuous'),
        'Default colormap for sequential datasets. Alias for :rcraw:`image.cmap`.'
    ),

    # Special setting
    'edgefix': (
        True,
        _validate_bool,
        'Whether to fix issues with "white lines" appearing between patches '
        'in saved vector graphics and with vector graphic backends. Applies '
        'to colorbar levels and bar, area, pcolor, and contour plots.'
    ),

    # Font settings
    'font.name': (
        FONTNAME,
        _validate_fontname,
        'Alias for :rcraw:`font.family`.'
    ),
    'font.small': (
        SMALLSIZE,
        _validate_fontsize,
        'Alias for :rcraw:`font.smallsize`.'
    ),
    'font.smallsize': (
        SMALLSIZE,
        _validate_fontsize,
        'Meta setting that changes the label-like sizes ``axes.labelsize``, '
        '``legend.fontsize``, ``tick.labelsize``, and ``grid.labelsize``. Default is '
        "``'medium'`` (equivalent to :rcraw:`font.size`)." + _addendum_font
    ),
    'font.large': (
        LARGESIZE,
        _validate_fontsize,
        'Alias for :rcraw:`font.largesize`.'
    ),
    'font.largesize': (
        LARGESIZE,
        _validate_fontsize,
        'Meta setting that changes the title-like sizes ``abc.size``, ``title.size``, '
        '``suptitle.size``, ``leftlabel.size``, ``rightlabel.size``, etc. Default is '
        "``'med-large'`` (i.e. 1.1 times :rcraw:`font.size`)." + _addendum_font
    ),

    # Formatter settings
    'formatter.timerotation': (
        'vertical',
        _validate_rotation,
        'Rotation for *x* axis datetime tick labels.' + _addendum_rotation
    ),
    'formatter.zerotrim': (
        True,
        _validate_bool,
        'Whether to trim trailing decimal zeros on tick labels.'
    ),
    'formatter.limits': (
        [-5, 6],  # must be list or else validated
        _validate['axes.formatter.limits'],
        'Alias for :rcraw:`axes.formatter.limits`.'
    ),
    'formatter.min_exponent': (
        0,
        _validate['axes.formatter.min_exponent'],
        'Alias for :rcraw:`axes.formatter.min_exponent`.'
    ),
    'formatter.offset_threshold': (
        4,
        _validate['axes.formatter.offset_threshold'],
        'Alias for :rcraw:`axes.formatter.offset_threshold`.'
    ),
    'formatter.use_locale': (
        False,
        _validate_bool,
        'Alias for :rcraw:`axes.formatter.use_locale`.'
    ),
    'formatter.use_mathtext': (
        MATHTEXT,
        _validate_bool,
        'Alias for :rcraw:`axes.formatter.use_mathtext`.'
    ),
    'formatter.use_offset': (
        True,
        _validate_bool,
        'Alias for :rcraw:`axes.formatter.useOffset`.'
    ),

    # Geographic axes settings
    'geo.backend': (
        'cartopy',
        _validate_belongs('cartopy', 'basemap'),
        'The backend used for `~proplot.axes.GeoAxes`. Must be '
        "either 'cartopy' or 'basemap'."
    ),
    'geo.extent': (
        'globe',
        _validate_belongs('globe', 'auto'),
        "If ``'globe'``, the extent of cartopy `~proplot.axes.GeoAxes` is always "
        "global. If ``'auto'``, the extent is automatically adjusted based on "
        "plotted content. Default is ``'globe'``."
    ),
    'geo.round': (
        True,
        _validate_bool,
        "If ``True`` (the default), polar `~proplot.axes.GeoAxes` like ``'npstere'`` "
        "and ``'spstere'`` are bounded with circles rather than squares."
    ),


    # Gridlines
    # NOTE: Here 'grid' and 'gridminor' or *not* aliases for native 'axes.grid' and
    # invented 'axes.gridminor' because native 'axes.grid' controls both major *and*
    # minor gridlines. Must handle it independently from these settings.
    'grid': (
        True,
        _validate_bool,
        'Toggle major gridlines on and off.'
    ),
    'grid.below': (
        GRIDBELOW,  # like axes.axisbelow
        _validate_belongs(False, 'line', True),
        'Alias for :rcraw:`axes.axisbelow`. If ``True``, draw gridlines below '
        "everything. If ``True``, draw them above everything. If ``'line'``, "
        'draw them above patches but below lines and markers.'
    ),
    'grid.checkoverlap': (
        True,
        _validate_bool,
        'Whether to have cartopy automatically check for and remove overlapping '
        '`~proplot.axes.GeoAxes` gridline labels.'
    ),
    'grid.dmslabels': (
        True,
        _validate_bool,
        'Whether to use degrees-minutes-seconds rather than decimals for '
        'cartopy `~proplot.axes.GeoAxes` gridlines.'
    ),
    'grid.geolabels': (
        True,
        _validate_bool,
        "Whether to include the ``'geo'`` spine in cartopy >= 0.20 when otherwise "
        'toggling left, right, bottom, or top `~proplot.axes.GeoAxes` gridline labels.'
    ),
    'grid.inlinelabels': (
        False,
        _validate_bool,
        'Whether to add inline labels for cartopy `~proplot.axes.GeoAxes` gridlines.'
    ),
    'grid.labels': (
        False,
        _validate_bool,
        'Whether to add outer labels for `~proplot.axes.GeoAxes` gridlines.'
    ),
    'grid.labelcolor': (
        BLACK,
        _validate_color,
        'Font color for `~proplot.axes.GeoAxes` gridline labels.'
    ),
    'grid.labelpad': (
        GRIDPAD,
        _validate_pt,
        'Padding between the map boundary and cartopy `~proplot.axes.GeoAxes` '
        'gridline labels.' + _addendum_pt
    ),
    'grid.labelsize': (
        SMALLSIZE,
        _validate_fontsize,
        'Font size for `~proplot.axes.GeoAxes` gridline labels.' + _addendum_font
    ),
    'grid.labelweight': (
        'normal',
        _validate_fontweight,
        'Font weight for `~proplot.axes.GeoAxes` gridline labels.'
    ),
    'grid.nsteps': (
        250,
        _validate_int,
        'Number of points used to draw cartopy `~proplot.axes.GeoAxes` gridlines.'
    ),
    'grid.pad': (
        GRIDPAD,
        _validate_pt,
        'Alias for :rcraw:`grid.labelpad`.'
    ),
    'grid.rotatelabels': (
        False,  # False limits projections where labels are available
        _validate_bool,
        'Whether to rotate cartopy `~proplot.axes.GeoAxes` gridline labels.'
    ),
    'grid.style': (
        '-',
        _validate_linestyle,
        'Major gridline style. Alias for :rcraw:`grid.linestyle`.'
    ),
    'grid.width': (
        LINEWIDTH,
        _validate_pt,
        'Major gridline width. Alias for :rcraw:`grid.linewidth`.'
    ),
    'grid.widthratio': (
        GRIDRATIO,
        _validate_float,
        'Ratio of minor gridline width to major gridline width.'
    ),

    # Minor gridlines
    'gridminor': (
        False,
        _validate_bool,
        'Toggle minor gridlines on and off.'
    ),
    'gridminor.alpha': (
        GRIDALPHA,
        _validate_float,
        'Minor gridline opacity.'
    ),
    'gridminor.color': (
        BLACK,
        _validate_color,
        'Minor gridline color.'
    ),
    'gridminor.linestyle': (
        GRIDSTYLE,
        _validate_linestyle,
        'Minor gridline style.'
    ),
    'gridminor.linewidth': (
        GRIDRATIO * LINEWIDTH,
        _validate_pt,
        'Minor gridline width.'
    ),
    'gridminor.style': (
        GRIDSTYLE,
        _validate_linestyle,
        'Minor gridline style. Alias for :rcraw:`gridminor.linestyle`.'
    ),
    'gridminor.width': (
        GRIDRATIO * LINEWIDTH,
        _validate_pt,
        'Minor gridline width. Alias for :rcraw:`gridminor.linewidth`.'
    ),

    # Backend stuff
    'inlineformat': (
        'retina',
        _validate_belongs('svg', 'pdf', 'retina', 'png', 'jpeg'),
        'The inline backend figure format. Valid formats include '
        "``'svg'``, ``'pdf'``, ``'retina'``, ``'png'``, and ``jpeg``."
    ),

    # Inner borders
    'innerborders': (
        False,
        _validate_bool,
        'Toggles internal political border lines (e.g. states and provinces) '
        'on and off.'
    ),
    'innerborders.alpha': (
        None,
        _validate_or_none(_validate_float),
        'Opacity for internal political border lines',
    ),
    'innerborders.color': (
        BLACK,
        _validate_color,
        'Line color for internal political border lines.'
    ),
    'innerborders.linewidth': (
        LINEWIDTH,
        _validate_pt,
        'Line width for internal political border lines.'
    ),
    'innerborders.zorder': (
        ZLINES,
        _validate_float,
        'Z-order for internal political border lines.'
    ),

    # Axis label settings
    'label.color': (
        BLACK,
        _validate_color,
        'Alias for :rcraw:`axes.labelcolor`.'
    ),
    'label.pad': (
        LABELPAD,
        _validate_pt,
        'Alias for :rcraw:`axes.labelpad`.'
        + _addendum_pt
    ),
    'label.size': (
        SMALLSIZE,
        _validate_fontsize,
        'Alias for :rcraw:`axes.labelsize`.' + _addendum_font
    ),
    'label.weight': (
        'normal',
        _validate_fontweight,
        'Alias for :rcraw:`axes.labelweight`.'
    ),

    # Lake patches
    'lakes': (
        False,
        _validate_bool,
        'Toggles lake patches on and off.'
    ),
    'lakes.alpha': (
        None,
        _validate_or_none(_validate_float),
        'Opacity for lake patches',
    ),
    'lakes.color': (
        WHITE,
        _validate_color,
        'Face color for lake patches.'
    ),
    'lakes.zorder': (
        ZPATCHES,
        _validate_float,
        'Z-order for lake patches.'
    ),

    # Land patches
    'land': (
        False,
        _validate_bool,
        'Toggles land patches on and off.'
    ),
    'land.alpha': (
        None,
        _validate_or_none(_validate_float),
        'Opacity for land patches',
    ),
    'land.color': (
        BLACK,
        _validate_color,
        'Face color for land patches.'
    ),
    'land.zorder': (
        ZPATCHES,
        _validate_float,
        'Z-order for land patches.'
    ),

    # Left subplot labels
    'leftlabel.color': (
        BLACK,
        _validate_color,
        'Font color for row labels on the left-hand side.'
    ),
    'leftlabel.pad': (
        TITLEPAD,
        _validate_pt,
        'Padding between axes content and row labels on the left-hand side.'
        + _addendum_pt
    ),
    'leftlabel.rotation': (
        'vertical',
        _validate_rotation,
        'Rotation for row labels on the left-hand side.' + _addendum_rotation
    ),
    'leftlabel.size': (
        LARGESIZE,
        _validate_fontsize,
        'Font size for row labels on the left-hand side.' + _addendum_font
    ),
    'leftlabel.weight': (
        'bold',
        _validate_fontweight,
        'Font weight for row labels on the left-hand side.'
    ),

    # Meta settings
    'margin': (
        MARGIN,
        _validate_float,
        'The fractional *x* and *y* axis data margins when limits are unset. '
        'Alias for :rcraw:`axes.margin`.'
    ),
    'meta.edgecolor': (
        BLACK,
        _validate_color,
        'Color of axis spines, tick marks, tick labels, and labels.'
    ),
    'meta.color': (
        BLACK,
        _validate_color,
        'Color of axis spines, tick marks, tick labels, and labels. '
        'Alias for :rcraw:`meta.edgecolor`.'
    ),
    'meta.linewidth': (
        LINEWIDTH,
        _validate_pt,
        'Thickness of axis spines and major tick lines.'
    ),
    'meta.width': (
        LINEWIDTH,
        _validate_pt,
        'Thickness of axis spines and major tick lines. '
        'Alias for :rcraw:`meta.linewidth`.'
    ),

    # For negative positive patches
    'negcolor': (
        'blue7',
        _validate_color,
        'Color for negative bars and shaded areas when using ``negpos=True``. '
        'See also :rcraw:`poscolor`.'
    ),
    'poscolor': (
        'red7',
        _validate_color,
        'Color for positive bars and shaded areas when using ``negpos=True``. '
        'See also :rcraw:`negcolor`.'
    ),

    # Ocean patches
    'ocean': (
        False,
        _validate_bool,
        'Toggles ocean patches on and off.'
    ),
    'ocean.alpha': (
        None,
        _validate_or_none(_validate_float),
        'Opacity for ocean patches',
    ),
    'ocean.color': (
        WHITE,
        _validate_color,
        'Face color for ocean patches.'
    ),
    'ocean.zorder': (
        ZPATCHES,
        _validate_float,
        'Z-order for ocean patches.'
    ),

    # Geographic resolution
    'reso': (
        'lo',
        _validate_belongs('lo', 'med', 'hi', 'x-hi', 'xx-hi'),
        'Resolution for `~proplot.axes.GeoAxes` geographic features. '
        "Must be one of ``'lo'``, ``'med'``, ``'hi'``, ``'x-hi'``, or ``'xx-hi'``."
    ),

    # Right subplot labels
    'rightlabel.color': (
        BLACK,
        _validate_color,
        'Font color for row labels on the right-hand side.'
    ),
    'rightlabel.pad': (
        TITLEPAD,
        _validate_pt,
        'Padding between axes content and row labels on the right-hand side.'
        + _addendum_pt
    ),
    'rightlabel.rotation': (
        'vertical',
        _validate_rotation,
        'Rotation for row labels on the right-hand side.' + _addendum_rotation
    ),
    'rightlabel.size': (
        LARGESIZE,
        _validate_fontsize,
        'Font size for row labels on the right-hand side.' + _addendum_font
    ),
    'rightlabel.weight': (
        'bold',
        _validate_fontweight,
        'Font weight for row labels on the right-hand side.'
    ),

    # River lines
    'rivers': (
        False,
        _validate_bool,
        'Toggles river lines on and off.'
    ),
    'rivers.alpha': (
        None,
        _validate_or_none(_validate_float),
        'Opacity for river lines.',
    ),
    'rivers.color': (
        BLACK,
        _validate_color,
        'Line color for river lines.'
    ),
    'rivers.linewidth': (
        LINEWIDTH,
        _validate_pt,
        'Line width for river lines.'
    ),
    'rivers.zorder': (
        ZLINES,
        _validate_float,
        'Z-order for river lines.'
    ),

    # Subplots settings
    'subplots.align': (
        False,
        _validate_bool,
        'Whether to align axis labels during draw. See `aligning labels '
        '<https://matplotlib.org/stable/gallery/subplots_axes_and_figures/align_labels_demo.html>`__.'  # noqa: E501
    ),
    'subplots.equalspace': (
        False,
        _validate_bool,
        'Whether to make the tight layout algorithm assign the same space for '
        'every row and the same space for every column.'
    ),
    'subplots.groupspace': (
        True,
        _validate_bool,
        'Whether to make the tight layout algorithm consider space between only '
        'adjacent subplot "groups" rather than every subplot in the row or column.'
    ),
    'subplots.innerpad': (
        1,
        _validate_em,
        'Padding between adjacent subplots.' + _addendum_em
    ),
    'subplots.outerpad': (
        0.5,
        _validate_em,
        'Padding around figure edge.' + _addendum_em
    ),
    'subplots.panelpad': (
        0.5,
        _validate_em,
        'Padding between subplots and panels, and between stacked panels.'
        + _addendum_em
    ),
    'subplots.panelwidth': (
        0.5,
        _validate_in,
        'Width of side panels.' + _addendum_in
    ),
    'subplots.refwidth': (
        2.5,
        _validate_in,
        'Default width of the reference subplot.' + _addendum_in
    ),
    'subplots.share': (
        True,
        _validate_belongs(0, 1, 2, 3, 4, False, 'labels', 'limits', True, 'all'),
        'The axis sharing level, one of ``0``, ``1``, ``2``, or ``3``, or the '
        "more intuitive aliases ``False``, ``'labels'``, ``'limits'``, or ``True``. "
        'See `~proplot.figure.Figure` for details.'
    ),
    'subplots.span': (
        True,
        _validate_bool,
        'Toggles spanning axis labels. See `~proplot.ui.subplots` for details.'
    ),
    'subplots.tight': (
        True,
        _validate_bool,
        'Whether to auto-adjust the subplot spaces and figure margins.'
    ),

    # Super title settings
    'suptitle.color': (
        BLACK,
        _validate_color,
        'Figure title color.'
    ),
    'suptitle.pad': (
        TITLEPAD,
        _validate_pt,
        'Padding between axes content and the figure super title.' + _addendum_pt
    ),
    'suptitle.size': (
        LARGESIZE,
        _validate_fontsize,
        'Figure title font size.' + _addendum_font
    ),
    'suptitle.weight': (
        'bold',
        _validate_fontweight,
        'Figure title font weight.'
    ),

    # Tick settings
    'tick.color': (
        BLACK,
        _validate_color,
        'Major and minor tick color.'
    ),
    'tick.dir': (
        TICKDIR,
        _validate_belongs('in', 'out', 'inout'),
        'Major and minor tick direction. Must be one of '
        "``'out'``, ``'in'``, or ``'inout'``."
    ),
    'tick.labelcolor': (
        BLACK,
        _validate_color,
        'Axis tick label color.'
    ),
    'tick.labelpad': (
        TICKPAD,
        _validate_pt,
        'Padding between ticks and tick labels.' + _addendum_pt
    ),
    'tick.labelsize': (
        SMALLSIZE,
        _validate_fontsize,
        'Axis tick label font size.' + _addendum_font
    ),
    'tick.labelweight': (
        'normal',
        _validate_fontweight,
        'Axis tick label font weight.'
    ),
    'tick.len': (
        TICKLEN,
        _validate_pt,
        'Length of major ticks in points.'
    ),
    'tick.lenratio': (
        TICKLENRATIO,
        _validate_float,
        'Ratio of minor tickline length to major tickline length.'
    ),
    'tick.linewidth': (
        LINEWIDTH,
        _validate_pt,
        'Major tickline width.'
    ),
    'tick.minor': (
        TICKMINOR,
        _validate_bool,
        'Toggles minor ticks on and off.',
    ),
    'tick.pad': (
        TICKPAD,
        _validate_pt,
        'Alias for :rcraw:`tick.labelpad`.'
    ),
    'tick.width': (
        LINEWIDTH,
        _validate_pt,
        'Major tickline width. Alias for :rcraw:`tick.linewidth`.'
    ),
    'tick.widthratio': (
        TICKWIDTHRATIO,
        _validate_float,
        'Ratio of minor tickline width to major tickline width.'
    ),

    # Title settings
    'title.above': (
        True,
        _validate_belongs(False, True, 'panels'),
        'Whether to move outer titles and a-b-c labels above panels, colorbars, or '
        "legends that are above the axes. If the string 'panels' then text is only "
        'redirected above axes panels. Otherwise should be boolean.'
    ),
    'title.border': (
        True,
        _validate_bool,
        'Whether to draw a white border around titles '
        'when :rcraw:`title.loc` is inside the axes.'
    ),
    'title.borderwidth': (
        1.5,
        _validate_pt,
        'Width of the border around titles.'
    ),
    'title.box': (
        False,
        _validate_bool,
        'Whether to draw semi-transparent bounding boxes around titles '
        'when :rcraw:`title.loc` is inside the axes.'
    ),
    'title.boxcolor': (
        WHITE,
        _validate_color,
        'Axes title bounding box color.'
    ),
    'title.boxalpha': (
        0.5,
        _validate_float,
        'Axes title bounding box opacity.'
    ),
    'title.boxstyle': (
        'square',
        _validate_boxstyle,
        'Axes title bounding box style.'
    ),
    'title.boxpad': (
        None,
        _validate_or_none(_validate_pt),
        'Padding for the title bounding box. By default this is scaled '
        'to make the box flush against the axes edge.' + _addendum_pt
    ),
    'title.color': (
        BLACK,
        _validate_color,
        'Axes title color. Alias for :rcraw:`axes.titlecolor`.'
    ),
    'title.loc': (
        'center',
        _validate_belongs(*TEXT_LOCS),
        'Title position. For options see the :ref:`location table <title_table>`.'
    ),
    'title.pad': (
        TITLEPAD,
        _validate_pt,
        'Padding between the axes edge and the inner and outer titles and '
        'a-b-c labels. Alias for :rcraw:`axes.titlepad`.' + _addendum_pt
    ),
    'title.size': (
        LARGESIZE,
        _validate_fontsize,
        'Axes title font size. Alias for :rcraw:`axes.titlesize`.' + _addendum_font
    ),
    'title.weight': (
        'normal',
        _validate_fontweight,
        'Axes title font weight. Alias for :rcraw:`axes.titleweight`.'
    ),

    # Top subplot label settings
    'toplabel.color': (
        BLACK,
        _validate_color,
        'Font color for column labels on the top of the figure.'
    ),
    'toplabel.pad': (
        TITLEPAD,
        _validate_pt,
        'Padding between axes content and column labels on the top of the figure.'
        + _addendum_pt
    ),
    'toplabel.rotation': (
        'horizontal',
        _validate_rotation,
        'Rotation for column labels at the top of the figure.' + _addendum_rotation
    ),
    'toplabel.size': (
        LARGESIZE,
        _validate_fontsize,
        'Font size for column labels on the top of the figure.' + _addendum_font
    ),
    'toplabel.weight': (
        'bold',
        _validate_fontweight,
        'Font weight for column labels on the top of the figure.'
    ),

    # Unit formatting
    'unitformat': (
        'L',
        _validate_string,
        'The format string used to format `pint.Quantity` default unit labels '
        'using ``format(units, unitformat)``. See also :rcraw:`autoformat`.'
    ),
}

# Child settings. Changing the parent changes all the children, but
# changing any of the children does not change the parent.
_rc_children = {
    'font.smallsize': (  # the 'small' fonts
        'tick.labelsize', 'xtick.labelsize', 'ytick.labelsize',
        'axes.labelsize', 'legend.fontsize', 'grid.labelsize'
    ),
    'font.largesize': (  # the 'large' fonts
        'abc.size', 'figure.titlesize', 'suptitle.size', 'axes.titlesize', 'title.size',
        'leftlabel.size', 'toplabel.size', 'rightlabel.size', 'bottomlabel.size'
    ),
    'meta.color': (  # change the 'color' of an axes
        'axes.edgecolor', 'axes.labelcolor', 'legend.edgecolor', 'colorbar.edgecolor',
        'tick.labelcolor', 'hatch.color', 'xtick.color', 'ytick.color'
    ),
    'meta.width': (  # change the tick and frame line width
        'axes.linewidth', 'tick.width', 'tick.linewidth', 'xtick.major.width',
        'ytick.major.width', 'grid.width', 'grid.linewidth',
    ),
    'axes.margin': ('axes.xmargin', 'axes.ymargin'),
    'grid.color': ('gridminor.color', 'grid.labelcolor'),
    'grid.alpha': ('gridminor.alpha',),
    'grid.linewidth': ('gridminor.linewidth',),
    'grid.linestyle': ('gridminor.linestyle',),
    'tick.color': ('xtick.color', 'ytick.color'),
    'tick.dir': ('xtick.direction', 'ytick.direction'),
    'tick.len': ('xtick.major.size', 'ytick.major.size'),
    'tick.minor': ('xtick.minor.visible', 'ytick.minor.visible'),
    'tick.pad': ('xtick.major.pad', 'xtick.minor.pad', 'ytick.major.pad', 'ytick.minor.pad'),  # noqa: E501
    'tick.width': ('xtick.major.width', 'ytick.major.width'),
    'tick.labelsize': ('xtick.labelsize', 'ytick.labelsize'),
}

# Recently added settings. Update these only if the version is recent enough
# NOTE: We don't make 'title.color' a child of 'axes.titlecolor' because
# the latter can take on the value 'auto' and can't handle that right now.
if _version_mpl >= '3.2':
    _rc_matplotlib_default['axes.titlecolor'] = BLACK
    _rc_children['title.color'] = ('axes.titlecolor',)
if _version_mpl >= '3.4':
    _rc_matplotlib_default['xtick.labelcolor'] = BLACK
    _rc_matplotlib_default['ytick.labelcolor'] = BLACK
    _rc_children['tick.labelcolor'] = ('xtick.labelcolor', 'ytick.labelcolor')
    _rc_children['grid.labelcolor'] = ('xtick.labelcolor', 'ytick.labelcolor')
    _rc_children['meta.color'] += ('xtick.labelcolor', 'ytick.labelcolor')

# Setting synonyms. Changing one setting changes the other. Also account for existing
# children. Most of these are aliased due to proplot settings overlapping with
# existing matplotlib settings.
_rc_synonyms = (
    ('cmap', 'image.cmap', 'cmap.sequential'),
    ('cmap.lut', 'image.lut'),
    ('font.name', 'font.family'),
    ('font.small', 'font.smallsize'),
    ('font.large', 'font.largesize'),
    ('formatter.limits', 'axes.formatter.limits'),
    ('formatter.use_locale', 'axes.formatter.use_locale'),
    ('formatter.use_mathtext', 'axes.formatter.use_mathtext'),
    ('formatter.min_exponent', 'axes.formatter.min_exponent'),
    ('formatter.use_offset', 'axes.formatter.useoffset'),
    ('formatter.offset_threshold', 'axes.formatter.offset_threshold'),
    ('grid.below', 'axes.axisbelow'),
    ('grid.labelpad', 'grid.pad'),
    ('grid.linewidth', 'grid.width'),
    ('grid.linestyle', 'grid.style'),
    ('gridminor.linewidth', 'gridminor.width'),
    ('gridminor.linestyle', 'gridminor.style'),
    ('label.color', 'axes.labelcolor'),
    ('label.pad', 'axes.labelpad'),
    ('label.size', 'axes.labelsize'),
    ('label.weight', 'axes.labelweight'),
    ('margin', 'axes.margin'),
    ('meta.width', 'meta.linewidth'),
    ('meta.color', 'meta.edgecolor'),
    ('tick.labelpad', 'tick.pad'),
    ('tick.labelsize', 'grid.labelsize'),
    ('tick.labelcolor', 'grid.labelcolor'),
    ('tick.labelweight', 'grid.labelweight'),
    ('tick.linewidth', 'tick.width'),
    ('title.pad', 'axes.titlepad'),
    ('title.size', 'axes.titlesize'),
    ('title.weight', 'axes.titleweight'),
)
for _keys in _rc_synonyms:
    for _key in _keys:
        _set = {_ for k in _keys for _ in {k, *_rc_children.get(k, ())}} - {_key}
        _rc_children[_key] = tuple(sorted(_set))

# Previously removed settings.
# NOTE: Initial idea was to defer deprecation warnings in Configurator to the
# subsequent RcParams indexing. However this turned out be complicated, because
# would have to detect deprecated keys in _get_item_dicts blocks, and need to
# validate values before e.g. applying 'tick.lenratio'. So the renamed parameters
# do not have to be added as _rc_children, since Configurator translates before
# retrieving the list of children in _get_item_dicts.
_rc_removed = {  # {key: (alternative, version)} dictionary
    'rgbcycle': ('', '0.6.0'),  # no alternative, we no longer offer this feature
    'geogrid.latmax': ('Please use ax.format(latmax=N) instead.', '0.6.0'),
    'geogrid.latstep': ('Please use ax.format(latlines=N) instead.', '0.6.0'),
    'geogrid.lonstep': ('Please use ax.format(lonlines=N) instead.', '0.6.0'),
    'gridminor.latstep': ('Please use ax.format(latminorlines=N) instead.', '0.6.0'),
    'gridminor.lonstep': ('Please use ax.format(lonminorlines=N) instead.', '0.6.0'),
}
_rc_renamed = {  # {old_key: (new_key, version)} dictionary
    'abc.format': ('abc', '0.5.0'),
    'align': ('subplots.align', '0.6.0'),
    'axes.facealpha': ('axes.alpha', '0.6.0'),
    'geoaxes.edgecolor': ('axes.edgecolor', '0.6.0'),
    'geoaxes.facealpha': ('axes.alpha', '0.6.0'),
    'geoaxes.facecolor': ('axes.facecolor', '0.6.0'),
    'geoaxes.linewidth': ('axes.linewidth', '0.6.0'),
    'geogrid.alpha': ('grid.alpha', '0.6.0'),
    'geogrid.color': ('grid.color', '0.6.0'),
    'geogrid.labels': ('grid.labels', '0.6.0'),
    'geogrid.labelpad': ('grid.pad', '0.6.0'),
    'geogrid.labelsize': ('grid.labelsize', '0.6.0'),
    'geogrid.linestyle': ('grid.linestyle', '0.6.0'),
    'geogrid.linewidth': ('grid.linewidth', '0.6.0'),
    'share': ('subplots.share', '0.6.0'),
    'small': ('font.smallsize', '0.6.0'),
    'large': ('font.largesize', '0.6.0'),
    'span': ('subplots.span', '0.6.0'),
    'tight': ('subplots.tight', '0.6.0'),
    'axes.formatter.timerotation': ('formatter.timerotation', '0.6.0'),
    'axes.formatter.zerotrim': ('formatter.zerotrim', '0.6.0'),
    'abovetop': ('title.above', '0.7.0'),
    'subplots.pad': ('subplots.outerpad', '0.7.0'),
    'subplots.axpad': ('subplots.innerpad', '0.7.0'),
    'subplots.axwidth': ('subplots.refwidth', '0.7.0'),
    'text.labelsize': ('font.smallsize', '0.8.0'),
    'text.titlesize': ('font.largesize', '0.8.0'),
    'alpha': ('axes.alpha', '0.8.0'),
    'facecolor': ('axes.facecolor', '0.8.0'),
    'edgecolor': ('meta.color', '0.8.0'),
    'color': ('meta.color', '0.8.0'),
    'linewidth': ('meta.width', '0.8.0'),
    'lut': ('cmap.lut', '0.8.0'),
    'image.levels': ('cmap.levels', '0.8.0'),
    'image.inbounds': ('cmap.inbounds', '0.8.0'),
    'image.discrete': ('cmap.discrete', '0.8.0'),
    'image.edgefix': ('edgefix', '0.8.0'),
    'tick.ratio': ('tick.widthratio', '0.8.0'),
    'grid.ratio': ('grid.widthratio', '0.8.0'),
    'abc.style': ('abc', '0.8.0'),
    'grid.loninline': ('grid.inlinelabels', '0.8.0'),
    'grid.latinline': ('grid.inlinelabels', '0.8.0'),
    'cmap.edgefix': ('edgefix', '0.9.0'),
    'basemap': ('geo.backend', '0.10.0'),
    'inlinefmt': ('inlineformat', '0.10.0'),
    'cartopy.circular': ('geo.round', '0.10.0'),
    'cartopy.autoextent': ('geo.extent', '0.10.0'),
    'colorbar.rasterize': ('colorbar.rasterized', '0.10.0'),
    'title.bbox': ('title.box', '0.11.0'),
    'title.bboxcolor': ('title.boxcolor', '0.11.0'),
    'title.bboxalpha': ('title.boxalpha', '0.11.0'),
    'title.bboxstyle': ('title.boxstyle', '0.11.0'),
    'title.bboxpad': ('title.boxpad', '0.11.0'),
    'abc.bbox': ('abc.box', '0.11.0'),
    'abc.bboxcolor': ('abc.boxcolor', '0.11.0'),
    'abc.bboxalpha': ('abc.boxalpha', '0.11.0'),
    'abc.bboxstyle': ('abc.boxstyle', '0.11.0'),
    'abc.bboxpad': ('abc.boxpad', '0.11.0'),
}

# Validate the default settings dictionaries using a custom proplot _RcParams
# and the original matplotlib RcParams. Also surreptitiously add proplot
# font settings to the font keys list (beoolean below always evalutes to True)
# font keys list whlie initializing.
_rc_proplot_default = {
    key: value for key, (value, _, _) in _rc_proplot_table.items()
}
_rc_proplot_validate = {
    key: validator for key, (_, validator, _) in _rc_proplot_table.items()
    if not (validator is _validate_fontsize and FONT_KEYS.add(key))
}
_rc_proplot_default = _RcParams(_rc_proplot_default, _rc_proplot_validate)
_rc_matplotlib_default = RcParams(_rc_matplotlib_default)

# Important joint matplotlib proplot constants
# NOTE: The 'nodots' dictionary should include removed and renamed settings
_rc_categories = {
    '.'.join(name.split('.')[:i + 1])
    for dict_ in (_rc_proplot_default, _rc_matplotlib_native)
    for name in dict_
    for i in range(len(name.split('.')) - 1)
}
_rc_nodots = {
    name.replace('.', ''): name
    for dict_ in (_rc_proplot_default, _rc_matplotlib_native, _rc_renamed, _rc_removed)
    for name in dict_.keys()
}
