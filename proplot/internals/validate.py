#!/usr/bin/env python3
"""
Configure validators for global settings.
"""
import functools
import re
from numbers import Integral, Real

import matplotlib.rcsetup as msetup
import numpy as np
from matplotlib import MatplotlibDeprecationWarning, RcParams
from matplotlib.colors import Colormap
from matplotlib.font_manager import font_scalings
from matplotlib.fontconfig_pattern import parse_fontconfig_pattern

from . import ic  # noqa: F401
from . import styles, warnings

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

# Preset legend locations and aliases
# TODO: Add additional inset colorbar locations.
LOCS_LEGEND = {
    'fill': 'fill',
    'inset': 'best',
    'i': 'best',
    0: 'best',
    1: 'upper right',
    2: 'upper left',
    3: 'lower left',
    4: 'lower right',
    5: 'center left',
    6: 'center right',
    7: 'lower center',
    8: 'upper center',
    9: 'center',
    'l': 'left',
    'r': 'right',
    'b': 'bottom',
    't': 'top',
    'c': 'center',
    'ur': 'upper right',
    'ul': 'upper left',
    'll': 'lower left',
    'lr': 'lower right',
    'cr': 'center right',
    'cl': 'center left',
    'uc': 'upper center',
    'lc': 'lower center',
}
LOCS_ALIGN = {
    key: val for key, val in LOCS_LEGEND.items()
    if isinstance(key, str)
    and val in ('left', 'right', 'top', 'bottom', 'center')
}
LOCS_PANEL = {
    key: val for key, val in LOCS_LEGEND.items()
    if isinstance(key, str)
    and val in ('left', 'right', 'top', 'bottom')
}
LOCS_COLORBAR = {
    key: val for key, val in LOCS_LEGEND.items()
    if val in (
        'fill', 'best', 'left', 'right', 'top', 'bottom',
        'upper left', 'upper right', 'lower left', 'lower right',
    )
}
LOCS_TEXT = {
    key: val for key, val in LOCS_LEGEND.items()
    if isinstance(key, str) and val in (
        'left', 'center', 'right',
        'upper left', 'upper center', 'upper right',
        'lower left', 'lower center', 'lower right',
    )
}


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


def _validate_belongs(*options, ignorecase=True):
    """
    Return a validator ensuring the item belongs in the list.
    """
    def _validate_belongs(value):  # noqa: E306
        for opt in options:
            if isinstance(value, str) and isinstance(opt, str):
                if ignorecase:
                    if value.lower() == opt.lower():
                        return opt
                else:
                    if value == opt:
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
        return value  # e.g. 'auto'
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
        f'Invalid font size {value!r}. Can be points or one of the preset scalings: '
        + ', '.join(map(repr, font_scalings))
        + '.'
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


def _validate_loc(loc, mode, **kwargs):
    """
    Validate the location specification.
    """
    # Create specific options dictionary
    if mode == 'align':
        loc_dict = LOCS_ALIGN
    elif mode == 'panel':
        loc_dict = LOCS_PANEL
    elif mode == 'legend':
        loc_dict = LOCS_LEGEND
    elif mode == 'colorbar':
        loc_dict = LOCS_COLORBAR
    elif mode == 'text':
        loc_dict = LOCS_TEXT
    else:
        raise ValueError(f'Invalid mode {mode!r}.')
    loc_dict = loc_dict.copy()
    loc_dict.update(kwargs)
    loc_dict.update({val: val for val in loc_dict.values()})

    # Translate the location
    # TODO: Implement 'best' colorbar location. Currently rely on a kludge.
    sentinel = object()
    default = kwargs.pop('default', sentinel)
    if default is not sentinel and (loc is None or loc is True):
        loc = default
    elif isinstance(loc, (str, Integral)):
        try:
            loc = loc_dict[loc]
        except KeyError:
            raise KeyError(
                f'Invalid {mode} location {loc!r}. Options are: '
                + ', '.join(map(repr, loc_dict))
                + '.'
            )
        if mode == 'colorbar' and loc == 'best':
            loc = 'lower right'
    elif (
        mode == 'legend'
        and np.iterable(loc)
        and len(loc) == 2
        and all(isinstance(l, Real) for l in loc)
    ):
        loc = tuple(loc)
    else:
        raise KeyError(f'Invalid {mode} location {loc!r}.')

    return loc


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
    Validate rotation arguments.
    """
    if isinstance(value, str) and value.lower() in ('horizontal', 'vertical'):
        return value
    return _validate_float(value)


def _validate_style(style):
    """
    Validate the style or stylesheet.
    """
    # NOTE: Important to ignore deprecation and blacklisted param warnings here. Only
    # show warnings once when we 'sync' the style in Configurator._get_item_dicts()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', MatplotlibDeprecationWarning)
        _, kw = styles._parse_style_spec(
            style, warn_blacklisted=False, allow_dictionary=False
        )
    return kw.pop('style')  # pull out from the returned 'rc_proplot' dictionary


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
    _validate['legend.loc'] = functools.partial(_validate_loc, mode='legend')
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
