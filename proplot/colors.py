#!/usr/bin/env python3
"""
Various colormap classes and colormap normalization classes.
"""
# NOTE: To avoid name conflicts between registered colormaps and colors, print
# set(pplt.colors._cmap_database) & set(pplt.colors._color_database) whenever
# you add new colormaps. v0.8 result is {'gray', 'marine', 'ocean', 'pink'} due
# to the MATLAB and GNUPlot colormaps. Want to minimize conflicts.
# NOTE: We feel that LinearSegmentedColormap should always be used for smooth color
# transitions while ListedColormap should always be used for qualitative color sets.
# Other sources use ListedColormap for dense "perceptually uniform" colormaps possibly
# seeking optimization. However testing reveals that initialization of even very
# dense 256-level colormaps is only 1.25ms vs. 0.25ms for a ListedColormap with the
# same data (+1ms). Also ListedColormap was designed for qualitative transitions
# because specifying N different from len(colors) will cyclically loop around the
# colors or truncate colors. So we translate the relevant ListedColormaps to
# LinearSegmentedColormaps for consistency. See :rc:`cmap.listedthresh`
import functools
import json
import os
import re
from collections.abc import MutableMapping, MutableSequence
from numbers import Integral, Number, Real
from xml.etree import ElementTree

import matplotlib.cm as mcm
import matplotlib.colors as mcolors
import numpy as np
import numpy.ma as ma

from .config import rc
from .externals import hsluv
from .internals import ic  # noqa: F401
from .internals import (
    _kwargs_to_args,
    _not_none,
    _pop_props,
    docstring,
    inputs,
    warnings,
)

__all__ = [
    'Color',
    'DiscreteColormap',
    'ContinuousColormap',
    'PerceptualColormap',
    'DiscreteNorm',
    'DivergingNorm',
    'SegmentedNorm',
    'ColorDatabase',
    'ColormapDatabase',
    'ListedColormap',  # deprecated
    'LinearSegmentedColormap',  # deprecated
    'PerceptuallyUniformColormap',  # deprecated
    'LinearSegmentedNorm',  # deprecated
]

# Default colormap properties
DEFAULT_NAME = '_no_name'
DEFAULT_SPACE = 'hsl'

# Color regexes
# NOTE: We do not compile hex regex because config.py needs this surrounded by \A\Z
_regex_hex = r'#(?:[0-9a-fA-F]{3,4}){2}'  # 6-8 digit hex
REGEX_HEX_MULTI = re.compile(_regex_hex)
REGEX_HEX_SINGLE = re.compile(rf'\A{_regex_hex}\Z')
REGEX_ADJUST = re.compile(r'\A(light|dark|medium|pale|charcoal)?\s*(gr[ea]y[0-9]?)?\Z')

# Colormap constants
CMAPS_CYCLIC = tuple(  # cyclic colormaps loaded from rgb files
    key.lower() for key in (
        'MonoCycle',
        'twilight',
        'Phase',
        'romaO',
        'brocO',
        'corkO',
        'vikO',
        'bamO',
    )
)
CMAPS_DIVERGING = {  # mirrored dictionary mapping for reversed names
    key.lower(): value.lower()
    for key1, key2 in (
        ('BR', 'RB'),
        ('NegPos', 'PosNeg'),
        ('CoolWarm', 'WarmCool'),
        ('ColdHot', 'HotCold'),
        ('DryWet', 'WetDry'),
        ('PiYG', 'GYPi'),
        ('PRGn', 'GnRP'),
        ('BrBG', 'GBBr'),
        ('PuOr', 'OrPu'),
        ('RdGy', 'GyRd'),
        ('RdBu', 'BuRd'),
        ('RdYlBu', 'BuYlRd'),
        ('RdYlGn', 'GnYlRd'),
    )
    for key, value in ((key1, key2), (key2, key1))
}
for _cmap_diverging in (  # remaining diverging cmaps (see PlotAxes._parse_cmap)
    'Div',
    'Vlag',
    'Spectral',
    'Balance',
    'Delta',
    'Curl',
    'roma',
    'broc',
    'cork',
    'vik',
    'bam',
    'lisbon',
    'tofino',
    'berlin',
    'vanimo',
):
    CMAPS_DIVERGING[_cmap_diverging.lower()] = _cmap_diverging.lower()
CMAPS_REMOVED = {
    'Blue0': '0.6.0',
    'Cool': '0.6.0',
    'Warm': '0.6.0',
    'Hot': '0.6.0',
    'Floral': '0.6.0',
    'Contrast': '0.6.0',
    'Sharp': '0.6.0',
    'Viz': '0.6.0',
}
CMAPS_RENAMED = {
    'GrayCycle': ('MonoCycle', '0.6.0'),
    'Blue1': ('Blues1', '0.7.0'),
    'Blue2': ('Blues2', '0.7.0'),
    'Blue3': ('Blues3', '0.7.0'),
    'Blue4': ('Blues4', '0.7.0'),
    'Blue5': ('Blues5', '0.7.0'),
    'Blue6': ('Blues6', '0.7.0'),
    'Blue7': ('Blues7', '0.7.0'),
    'Blue8': ('Blues8', '0.7.0'),
    'Blue9': ('Blues9', '0.7.0'),
    'Green1': ('Greens1', '0.7.0'),
    'Green2': ('Greens2', '0.7.0'),
    'Green3': ('Greens3', '0.7.0'),
    'Green4': ('Greens4', '0.7.0'),
    'Green5': ('Greens5', '0.7.0'),
    'Green6': ('Greens6', '0.7.0'),
    'Green7': ('Greens7', '0.7.0'),
    'Green8': ('Greens8', '0.7.0'),
    'Orange1': ('Yellows1', '0.7.0'),
    'Orange2': ('Yellows2', '0.7.0'),
    'Orange3': ('Yellows3', '0.7.0'),
    'Orange4': ('Oranges2', '0.7.0'),
    'Orange5': ('Oranges1', '0.7.0'),
    'Orange6': ('Oranges3', '0.7.0'),
    'Orange7': ('Oranges4', '0.7.0'),
    'Orange8': ('Yellows4', '0.7.0'),
    'Brown1': ('Browns1', '0.7.0'),
    'Brown2': ('Browns2', '0.7.0'),
    'Brown3': ('Browns3', '0.7.0'),
    'Brown4': ('Browns4', '0.7.0'),
    'Brown5': ('Browns5', '0.7.0'),
    'Brown6': ('Browns6', '0.7.0'),
    'Brown7': ('Browns7', '0.7.0'),
    'Brown8': ('Browns8', '0.7.0'),
    'Brown9': ('Browns9', '0.7.0'),
    'RedPurple1': ('Reds1', '0.7.0'),
    'RedPurple2': ('Reds2', '0.7.0'),
    'RedPurple3': ('Reds3', '0.7.0'),
    'RedPurple4': ('Reds4', '0.7.0'),
    'RedPurple5': ('Reds5', '0.7.0'),
    'RedPurple6': ('Purples1', '0.7.0'),
    'RedPurple7': ('Purples2', '0.7.0'),
    'RedPurple8': ('Purples3', '0.7.0'),
}

# Color constants
COLORS_OPEN = {}  # populated during register_colors
COLORS_XKCD = {}  # populated during register_colors
COLORS_KEEP = (
    *(  # always load these XKCD colors regardless of settings
        'charcoal', 'tomato', 'burgundy', 'maroon', 'burgundy', 'lavendar',
        'taupe', 'sand', 'stone', 'earth', 'sand brown', 'sienna',
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
    *(  # common combinations
        'red orange', 'yellow orange', 'yellow green',
        'blue green', 'blue violet', 'red violet',
        'bright red',  # backwards compatibility
    ),
    *(  # common names
        prefix + color
        for color in (
            'red', 'orange', 'yellow', 'green', 'blue', 'indigo', 'violet',
            'brown', 'grey', 'gray',
        )
        for prefix in ('', 'light ', 'dark ', 'medium ', 'pale ')
    )
)
COLORS_REMOVE = (
    # filter these out, let's try to be professional here...
    'shit',
    'poop',
    'poo',
    'pee',
    'piss',
    'puke',
    'vomit',
    'snot',
    'booger',
    'bile',
    'diarrhea',
    'icky',
    'sickly',
)
COLORS_REPLACE = (
    # prevent registering similar-sounding names
    # these can all be combined
    ('/', ' '),  # convert [color1]/[color2] to compound (e.g. grey/blue to grey blue)
    ("'s", 's'),  # robin's egg
    ('egg blue', 'egg'),  # robin's egg blue
    ('grey', 'gray'),  # 'Murica
    ('ochre', 'ocher'),  # ...
    ('forrest', 'forest'),  # ...
    ('ocre', 'ocher'),  # correct spelling
    ('kelley', 'kelly'),  # ...
    ('reddish', 'red'),  # remove [color]ish where it modifies the spelling of color
    ('purplish', 'purple'),  # ...
    ('pinkish', 'pink'),
    ('yellowish', 'yellow'),
    ('bluish', 'blue'),
    ('greyish', 'grey'),
    ('ish', ''),  # these are all [color]ish ('ish' substring appears nowhere else)
    ('bluey', 'blue'),  # remove [color]y trailing y
    ('greeny', 'green'),  # ...
    ('reddy', 'red'),
    ('pinky', 'pink'),
    ('purply', 'purple'),
    ('purpley', 'purple'),
    ('yellowy', 'yellow'),
    ('orangey', 'orange'),
    ('browny', 'brown'),
    ('minty', 'mint'),  # now remove [object]y trailing y
    ('grassy', 'grass'),  # ...
    ('mossy', 'moss'),
    ('dusky', 'dusk'),
    ('rusty', 'rust'),
    ('muddy', 'mud'),
    ('sandy', 'sand'),
    ('leafy', 'leaf'),
    ('dusty', 'dust'),
    ('dirty', 'dirt'),
    ('peachy', 'peach'),
    ('stormy', 'storm'),
    ('cloudy', 'cloud'),
    ('grayblue', 'gray blue'),  # separate merge compounds
    ('bluegray', 'gray blue'),  # ...
    ('lightblue', 'light blue'),
    ('yellowgreen', 'yellow green'),
    ('yelloworange', 'yellow orange'),
)

# Color snippets
_docstring_rgba = """
color : color-spec
    A color specification. Sanitized with `to_rgba`.
"""
_docstring_hex = """
color : str
    An 8-digit HEX string indicating red, green, blue, and alpha channel values.
"""
_docstring_convert = """
color : color-spec
    The color. Can be a 3-tuple or 4-tuple of channel values, a hex
    string, a registered color name, a cycle color like ``'C0'``, or
    a 2-tuple colormap coordinate specification like ``('magma', 0.5)``
    (see `~proplot.colors.ColorDatabase` for details).

    If `space` is ``'rgb'``, this is a tuple of RGB values, and any
    channels are larger than ``2``, the channels are assumed to be
    on the ``0`` to ``255`` scale and are divided by ``255``.
space : {'rgb', 'hsv', 'hcl', 'hpl', 'hsl'}, optional
    The colorspace for the input channel values. Ignored unless `color`
    is a tuple of numbers.
cycle : str, default: :rcraw:`cycle`
    The registered color cycle name used to interpret colors that
    look like ``'C0'``, ``'C1'``, etc.
clip : bool, default: True
    Whether to clip channel values into the valid ``0`` to ``1`` range.
    Setting this to ``False`` can result in invalid colors.
"""
_docstring_seealso = """
Color.set_hue
Color.set_saturation
Color.set_luminance
Color.set_alpha
Color.shift_hue
Color.scale_luminance
Color.scale_saturation
"""
docstring._snippet_manager['colors.color'] = _docstring_rgba
docstring._snippet_manager['colors.hex'] = _docstring_hex
docstring._snippet_manager['colors.convert'] = _docstring_convert
docstring._snippet_manager['colors.seealso'] = _docstring_seealso

# Colormap snippets
_N_docstring = """
N : int, default: :rc:`image.lut`
    Number of points in the colormap lookup table.
"""
_alpha_docstring = """
alpha : float, optional
    The opacity for the entire colormap. This overrides
    the input opacities.
"""
_cyclic_docstring = """
cyclic : bool, optional
    Whether the colormap is cyclic. If ``True``, this changes how the leftmost
    and rightmost color levels are selected, and `extend` can only be
    ``'neither'`` (a warning will be issued otherwise).
"""
_gamma_docstring = """
gamma : float, optional
    Set `gamma1` and `gamma2` to this identical value.
gamma1 : float, optional
    If greater than 1, make low saturation colors more prominent. If
    less than 1, make high saturation colors more prominent. Similar to
    the `HCLWizard <http://hclwizard.org:64230/hclwizard/>`_ option.
gamma2 : float, optional
    If greater than 1, make high luminance colors more prominent. If
    less than 1, make low luminance colors more prominent. Similar to
    the `HCLWizard <http://hclwizard.org:64230/hclwizard/>`_ option.
"""
_space_docstring = """
space : {'hsl', 'hpl', 'hcl', 'hsv'}, optional
    The hue, saturation, luminance-style colorspace to use for interpreting
    the channels. See `this page <http://www.hsluv.org/comparison/>`__ for
    a full description.
"""
_name_docstring = """
name : str, default: '_no_name'
    The colormap name. This can also be passed as the first
    positional string argument.
"""
_ratios_docstring = """
ratios : sequence of float, optional
    Relative extents of each color transition. Must have length
    ``len(colors) - 1``. Larger numbers indicate a slower
    transition, smaller numbers indicate a faster transition.
"""
docstring._snippet_manager['colors.N'] = _N_docstring
docstring._snippet_manager['colors.alpha'] = _alpha_docstring
docstring._snippet_manager['colors.cyclic'] = _cyclic_docstring
docstring._snippet_manager['colors.gamma'] = _gamma_docstring
docstring._snippet_manager['colors.space'] = _space_docstring
docstring._snippet_manager['colors.ratios'] = _ratios_docstring
docstring._snippet_manager['colors.name'] = _name_docstring

# List classmethod snippets
_from_list_docstring = """
colors : sequence of color-spec or tuple
    If a sequence of RGB[A] tuples or color strings, the colormap
    transitions evenly from ``colors[0]`` at the left-hand side
    to ``colors[-1]`` at the right-hand side.

    If a sequence of (float, color-spec) tuples, the float values are the
    coordinate of each transition and must range from 0 to 1. This
    can be used to divide  the colormap range unevenly.
%(colors.name)s
%(colors.ratios)s
    For example, ``('red', 'blue', 'green')`` with ``ratios=(2, 1)``
    creates a colormap with the transition from red to blue taking
    *twice as long* as the transition from blue to green.
"""
docstring._snippet_manager['colors.from_list'] = _from_list_docstring


def _clip_colors(colors, clip=True, gray=0.2, warn_invalid=False):
    """
    Clip impossible colors rendered in an HSL-to-RGB colorspace
    conversion. Used by `PerceptualColormap`.

    Parameters
    ----------
    colors : sequence of 3-tuple
        The RGB colors.
    clip : bool, optional
        If `clip` is ``True`` (the default), RGB channel values >1
        are clipped to 1. Otherwise, the color is masked out as gray.
    gray : float, optional
        The identical RGB channel values (gray color) to be used
        if `clip` is ``True``.
    """
    colors = np.asarray(colors)
    under = colors < 0
    over = colors > 1
    if clip:
        colors[under], colors[over] = 0, 1
    else:
        colors[under | over] = gray
    if warn_invalid:
        msg = 'Clipped' if clip else 'Invalid'
        for i, name in enumerate('rgb'):
            if np.any(under[:, i]) or np.any(over[:, i]):
                warnings._warn_proplot(f'{msg} {name!r} channel.')
    return colors


def _color_channel(color, channel, space='hcl'):
    """
    Get the hue, saturation, or luminance channel value from the input color. The
    color name `color` can optionally be a string with the format ``'color+x'``
    or ``'color-x'``, where `x` is the offset from the channel value.

    Parameters
    ----------
    color : color-spec
        The color. Sanitized with `to_rgba`.
    channel : optional
        The HCL channel to be retrieved.
    space : optional
        The colorspace for the corresponding channel value.

    Returns
    -------
    value : float
        The channel value.
    """
    # Interpret channel
    if callable(color) or isinstance(color, Number):
        return color
    if channel == 'hue':
        channel = 0
    elif channel in ('chroma', 'saturation'):
        channel = 1
    elif channel == 'luminance':
        channel = 2
    else:
        raise ValueError(f'Unknown channel {channel!r}.')
    # Interpret string or RGB tuple
    offset = 0
    if isinstance(color, str):
        m = re.search('([-+][0-9.]+)$', color)
        if m:
            offset = float(m.group(0))
            color = color[:m.start()]
    return offset + to_xyz(color, space)[channel]


def _make_segment_data(values, coords=None, ratios=None):
    """
    Return a segmentdata array or callable given the input colors
    and coordinates.

    Parameters
    ----------
    values : sequence of float
        The channel values.
    coords : sequence of float, optional
        The segment coordinates.
    ratios : sequence of float, optional
        The relative length of each segment transition.
    """
    # Allow callables
    if callable(values):
        return values
    values = np.atleast_1d(values)
    if len(values) == 1:
        value = values[0]
        return [(0, value, value), (1, value, value)]

    # Get coordinates
    if not np.iterable(values):
        raise TypeError('Colors must be iterable, got {values!r}.')
    if coords is not None:
        coords = np.atleast_1d(coords)
        if ratios is not None:
            warnings._warn_proplot(
                f'Segment coordinates were provided, ignoring '
                f'ratios={ratios!r}.'
            )
        if len(coords) != len(values) or coords[0] != 0 or coords[-1] != 1:
            raise ValueError(
                f'Coordinates must range from 0 to 1, got {coords!r}.'
            )
    elif ratios is not None:
        coords = np.atleast_1d(ratios)
        if len(coords) != len(values) - 1:
            raise ValueError(
                f'Need {len(values) - 1} ratios for {len(values)} colors, '
                f'but got {len(coords)} ratios.'
            )
        coords = np.concatenate(([0], np.cumsum(coords)))
        coords = coords / np.max(coords)  # normalize to 0-1
    else:
        coords = np.linspace(0, 1, len(values))

    # Build segmentdata array
    array = []
    for c, value in zip(coords, values):
        array.append((c, value, value))
    return array


def _make_lookup_table(N, data, gamma=1.0, inverse=False):
    r"""
    Generate lookup tables of HSL values given specified gradations. Similar to
    `~matplotlib.colors.makeMappingArray` but permits *circular* hue gradations,
    disables clipping of out-of-bounds values, and uses fancier "gamma" scaling.

    Parameters
    ----------
    N : int
        Number of points in the colormap lookup table.
    data : array-like
        Sequence of `(x, y_0, y_1)` tuples specifying channel jumps
        (from `y_0` to `y_1`) and `x` coordinate of those jumps
        (ranges between 0 and 1). See `~matplotlib.colors.LinearSegmentedColormap`.
    gamma : float or sequence of float, optional
        To obtain channel values between coordinates `x_i` and `x_{i+1}`
        in rows `i` and `i+1` of `data` we use the formula:

        .. math::

            y = y_{1,i} + w_i^{\gamma_i}*(y_{0,i+1} - y_{1,i})

        where `\gamma_i` corresponds to `gamma` and the weight `w_i` ranges from
        0 to 1 between rows ``i`` and ``i+1``. If `gamma` is float, it applies
        to every transition. Otherwise, its length must equal ``data.shape[0]-1``.

        This is similar to the `matplotlib.colors.makeMappingArray` `gamma` except
        it controls the weighting for transitions *between* each segment data
        coordinate rather than the coordinates themselves. This makes more sense
        for `PerceptualColormap`\ s because they usually contain just a
        handful of transitions representing chained segments.
    inverse : bool, optional
        If ``True``, `w_i^{\gamma_i}` is replaced with `1 - (1 - w_i)^{\gamma_i}` --
        that is, when `gamma` is greater than 1, this weights colors toward *higher*
        channel values instead of lower channel values.

        This is implemented in case we want to apply *equal* "gamma scaling"
        to different HSL channels in different directions. Usually, this
        is done to weight low data values with higher luminance *and* lower
        saturation, thereby emphasizing "extreme" data values.
    """
    # Allow for *callable* instead of linearly interpolating between segments
    gammas = np.atleast_1d(gamma)
    if np.any(gammas < 0.01) or np.any(gammas > 10):
        raise ValueError('Gamma can only be in range [0.01,10].')
    if callable(data):
        if len(gammas) > 1:
            raise ValueError('Only one gamma allowed for functional segmentdata.')
        x = np.linspace(0, 1, N)**gamma
        lut = np.array(data(x), dtype=float)
        return lut

    # Get array
    data = np.array(data)
    shape = data.shape
    if len(shape) != 2 or shape[1] != 3:
        raise ValueError('Mapping data must have shape N x 3.')
    if len(gammas) != 1 and len(gammas) != shape[0] - 1:
        raise ValueError(f'Expected {shape[0] - 1} gammas for {shape[0]} coords. Got {len(gamma)}.')  # noqa: E501
    if len(gammas) == 1:
        gammas = np.repeat(gammas, shape[:1])

    # Get indices
    x = data[:, 0]
    y0 = data[:, 1]
    y1 = data[:, 2]
    if x[0] != 0.0 or x[-1] != 1.0:
        raise ValueError('Data mapping points must start with x=0 and end with x=1.')
    if np.any(np.diff(x) < 0):
        raise ValueError('Data mapping points must have x in increasing order.')
    x = x * (N - 1)

    # Get distances from the segmentdata entry to the *left* for each requested
    # level, excluding ends at (0, 1), which must exactly match segmentdata ends.
    # NOTE: numpy.searchsorted returns where xq[i] must be inserted so it is
    # larger than x[ind[i]-1] but smaller than x[ind[i]].
    xq = (N - 1) * np.linspace(0, 1, N)
    ind = np.searchsorted(x, xq)[1:-1]
    offsets = (xq[1:-1] - x[ind - 1]) / (x[ind] - x[ind - 1])

    # Scale distances in each segment by input gamma
    # The ui are starting-points, the ci are counts from that point over which
    # segment applies (i.e. where to apply the gamma), the relevant 'segment'
    # is to the *left* of index returned by searchsorted
    _, uind, cind = np.unique(ind, return_index=True, return_counts=True)
    for ui, ci in zip(uind, cind):  # length should be N-1
        gamma = gammas[ind[ui] - 1]  # the relevant segment is *left* of this number
        if gamma == 1:
            continue
        if ci == 0:  # no lookup table coordinates fall inside this segment
            reverse = False
        else:  # reverse if we are transitioning to *lower* channel value
            reverse = (y0[ind[ui]] - y1[ind[ui] - 1]) < 0
        if inverse:  # reverse if we are transitioning to *higher* channel value
            reverse = not reverse
        if reverse:
            offsets[ui:ui + ci] = 1 - (1 - offsets[ui:ui + ci]) ** gamma
        else:
            offsets[ui:ui + ci] **= gamma

    # Perform successive linear interpolations rolled up into one equation
    lut = np.zeros((N,), float)
    lut[1:-1] = y1[ind - 1] + offsets * (y0[ind] - y1[ind - 1])
    lut[0] = y1[0]
    lut[-1] = y0[-1]
    return lut


def _load_colors(path, warn_on_failure=True):
    """
    Read colors from the input file.

    Parameters
    ----------
    warn_on_failure : bool, optional
        If ``True``, issue a warning when loading fails instead of raising an error.
    """
    # Warn or raise error (matches Colormap._from_file behavior)
    if not os.path.isfile(path):
        message = f'Failed to load color data file {path!r}. File not found.'
        if warn_on_failure:
            warnings._warn_proplot(message)
        else:
            raise FileNotFoundError(message)

    # Iterate through lines
    loaded = {}
    with open(path, 'r') as fh:
        for count, line in enumerate(fh):
            stripped = line.strip()
            if not stripped or stripped[0] == '#':
                continue
            pair = tuple(item.strip().lower() for item in line.split(':'))
            if len(pair) != 2 or not REGEX_HEX_SINGLE.match(pair[1]):
                warnings._warn_proplot(
                    f'Illegal line #{count + 1} in color file {path!r}:\n'
                    f'{line!r}\n'
                    f'Lines must be formatted as "name: hexcolor".'
                )
                continue
            loaded[pair[0]] = pair[1]

    return loaded


def _standardize_colors(input, space, margin):
    """
    Standardize the input colors.

    Parameters
    ----------
    input : dict
        The colors.
    space : optional
        The colorspace used to filter colors.
    margin : optional
        The proportional margin required for unique colors (e.g. 0.1
        is 36 hue units, 10 saturation units, 10 luminance units).
    """
    output = {}
    colors = []
    channels = []

    # Always add these colors and ignore other colors that are too close
    # We do this for colors with nice names or that proplot devs really like
    for name in COLORS_KEEP:
        color = input.pop(name, None)
        if color is None:
            continue
        if 'grey' in name:
            name = name.replace('grey', 'gray')
        colors.append((name, color))
        channels.append(to_xyz(color, space=space))
        output[name] = color  # required in case "kept" colors are close to each other

    # Translate remaining colors and remove bad names
    # WARNING: Unique axis argument requires numpy version >=1.13
    for name, color in input.items():
        for sub, rep in COLORS_REPLACE:
            if sub in name:
                name = name.replace(sub, rep)
        if any(sub in name for sub in COLORS_REMOVE):
            continue  # remove "unpofessional" names
        if name in output:
            continue  # prioritize names that come first
        colors.append((name, color))  # category name pair
        channels.append(to_xyz(color, space=space))

    # Get locations of "perceptually distinct" colors
    channels = np.asarray(channels)
    if not channels.size:
        return output
    channels = channels / np.array([360, 100, 100])
    channels = np.round(channels / margin).astype(np.int64)
    _, idxs = np.unique(channels, return_index=True, axis=0)

    # Return only "distinct" colors
    for idx in idxs:
        name, color = colors[idx]
        output[name] = color
    return output


class Color(MutableSequence, list):
    """
    A basic color class for storing HEX strings.
    """
    def _transform_color(self, func, space):
        """
        Standardize input for color transformation functions.
        """
        *color, opacity = self.to_rgba()
        color = self.to_xyz(color, space=space)
        color = func(list(color))  # apply transform
        return self.to_hex((*color, opacity), space=space)

    @docstring._snippet_manager
    def shift_hue(self, shift=0, space='hcl'):
        """
        Shift the hue channel of a color.

        Parameters
        ----------
        %(colors.color)s
        shift : float, optional
            The HCL hue channel is offset by this value.
        %(utils.space)s

        Returns
        -------
        %(colors.hex)s

        See also
        --------
        %(colors.seealso)s
        """
        def func(channels):
            channels[0] += shift
            channels[0] %= 360
            return channels
        return self._transform_color(func, space)

    @docstring._snippet_manager
    def scale_saturation(self, scale=1, space='hcl'):
        """
        Scale the saturation channel of a color.

        Parameters
        ----------
        %(colors.color)s
        scale : float, optional
            The HCL saturation channel is multiplied by this value.
        %(utils.space)s

        Returns
        -------
        %(colors.hex)s

        See also
        --------
        %(colors.seealso)s
        """
        def func(channels):
            channels[1] *= scale
            return channels
        return self._transform_color(func, space)

    @docstring._snippet_manager
    def scale_luminance(self, scale=1, space='hcl'):
        """
        Scale the luminance channel of a color.

        Parameters
        ----------
        %(colors.color)s
        scale : float, optional
            The luminance channel is multiplied by this value.
        %(utils.space)s

        Returns
        -------
        %(colors.hex)s

        See also
        --------
        %(colors.seealso)s
        """
        def func(channels):
            channels[2] *= scale
            return channels
        return self._transform_color(func, space)

    @docstring._snippet_manager
    def set_hue(self, hue, space='hcl'):
        """
        Return a color with a different hue and the same luminance and saturation
        as the input color.

        Parameters
        ----------
        %(colors.color)s
        hue : float, optional
            The new hue. Should lie between ``0`` and ``360`` degrees.
        %(utils.space)s

        Returns
        -------
        %(colors.hex)s

        See also
        --------
        %(colors.seealso)s
        """
        def func(channels):
            channels[0] = hue
            return channels
        return self._transform_color(func, space)

    @docstring._snippet_manager
    def set_saturation(self, saturation, space='hcl'):
        """
        Return a color with a different saturation and the same hue and luminance
        as the input color.

        Parameters
        ----------
        %(colors.color)s
        saturation : float, optional
            The new saturation. Should lie between ``0`` and ``360`` degrees.
        %(utils.space)s

        Returns
        -------
        %(colors.hex)s

        See also
        --------
        %(colors.seealso)s
        """
        def func(channels):
            channels[1] = saturation
            return channels
        return self._transform_color(func, space)

    @docstring._snippet_manager
    def set_luminance(self, luminance, space='hcl'):
        """
        Return a color with a different luminance and the same hue and saturation
        as the input color.

        Parameters
        ----------
        %(colors.color)s
        luminance : float, optional
            The new luminance. Should lie between ``0`` and ``100``.
        %(utils.space)s

        Returns
        -------
        %(colors.hex)s

        See also
        --------
        %(colors.seealso)s
        """
        def func(channels):
            channels[2] = luminance
            return channels
        return self._transform_color(func, space)

    @docstring._snippet_manager
    def set_alpha(self, alpha):
        """
        Return a color with the opacity channel set to the specified value.

        Parameters
        ----------
        %(colors.color)s
        alpha : float, optional
            The new opacity. Should be between ``0`` and ``1``.

        Returns
        -------
        %(colors.hex)s

        See also
        --------
        %(colors.seealso)s
        """
        def func(channels):
            channels[3] = alpha
            return channels
        return self._transform_color(func, alpha)

    def _translate_cycle_color(self, cycle=None):
        """
        Parse the input cycle color.
        """
        if isinstance(cycle, str):
            from .colors import _cmap_database
            try:
                cycle = _cmap_database[cycle].colors
            except (KeyError, AttributeError):
                cycles = sorted(
                    name
                    for name, cmap in _cmap_database.items()
                    if isinstance(cmap, mcolors.ListedColormap)
                )
                raise ValueError(
                    f'Invalid color cycle {cycle!r}. Options are: '
                    + ', '.join(map(repr, cycles))
                    + '.'
                )
        elif cycle is None:
            cycle = rc['axes.prop_cycle'].by_key()
            if 'color' not in cycle:
                cycle = ['k']
            else:
                cycle = cycle['color']
        else:
            raise ValueError(f'Invalid cycle {cycle!r}.')
        return cycle[int(self[-1]) % len(cycle)]

    @docstring._snippet_manager
    def to_hex(self, color, space='rgb', cycle=None, keep_alpha=True):
        """
        Translate the color from an arbitrary colorspace to a HEX string.
        This is a generalization of `matplotlib.colors.convert_hex`.

        Parameters
        ----------
        %(colors.convert)s
        keep_alpha : bool, default: True
            Whether to keep the opacity channel. If ``True`` an 8-digit HEX
            is returned. Otherwise a 6-digit HEX is returned.

        Returns
        -------
        %(colors.hex)s

        See also
        --------
        to_rgb
        to_rgba
        to_xyz
        to_xyza
        """
        rgba = self.to_rgba(color, space=space, cycle=cycle)
        return mcolors.convert_hex(rgba, keep_alpha=keep_alpha)

    @docstring._snippet_manager
    def to_rgb(self, color, space='rgb', cycle=None):
        """
        Translate the color from an arbitrary colorspace to an RGB tuple. This is
        a generalization of `matplotlib.colors.convert_rgb` and the inverse of `to_xyz`.

        Parameters
        ----------
        %(colors.convert)s

        Returns
        -------
        color : 3-tuple
            An RGB tuple.

        See also
        --------
        to_hex
        to_rgba
        to_xyz
        to_xyza
        """
        return self.to_rgba(color, space=space, cycle=cycle)[:3]

    @docstring._snippet_manager
    def to_rgba(self, color, space='rgb', cycle=None, clip=True):
        """
        Translate the color from an arbitrary colorspace to an RGBA tuple. This is a
        generalization of `matplotlib.colors.convert_rgba` and the inverse of `to_xyz`.

        Parameters
        ----------
        %(colors.convert)s

        Returns
        -------
        color : 4-tuple
            An RGBA tuple.

        See also
        --------
        to_hex
        to_rgb
        to_xyz
        to_xyza
        """
        # Translate color cycle strings
        if re.match(r'\AC[0-9]\Z', color):
            color = self._translate_cycle_color(cycle=cycle)

        # Translate RGB strings and (colormap, index) tuples
        # NOTE: Cannot use is_color_like because might have HSL channel values
        opacity = 1
        if (
            isinstance(color, str)
            or np.iterable(color) and len(color) == 2
        ):
            color = mcolors.convert_rgba(color)  # also enforced validity
        if (
            not np.iterable(color)
            or len(color) not in (3, 4)
            or not all(isinstance(c, Real) for c in color)
        ):
            raise ValueError(f'Invalid color-spec {color!r}.')
        if len(color) == 4:
            *color, opacity = color

        # Translate arbitrary colorspaces
        if space == 'rgb':
            if any(c > 2 for c in color):
                color = tuple(c / 255 for c in color)  # scale to within 0-1
            else:
                pass
        elif space == 'hsv':
            color = hsluv.hsl_to_rgb(*color)
        elif space == 'hcl':
            color = hsluv.hcl_to_rgb(*color)
        elif space == 'hsl':
            color = hsluv.hsluv_to_rgb(*color)
        elif space == 'hpl':
            color = hsluv.hpluv_to_rgb(*color)
        else:
            raise ValueError(f'Invalid colorspace {space!r}.')

        # Clip values. This should only be disabled when testing
        # translation functions.
        if clip:
            color = np.clip(color, 0, 1)  # clip to valid range

        # Return RGB or RGBA
        return (*color, opacity)

    @docstring._snippet_manager
    def to_xyz(self, color, space='hcl'):
        """
        Translate color in *any* format to a tuple of channel values in *any*
        colorspace. This is the inverse of `to_rgb`.

        Parameters
        ----------
        %(colors.color)s
        space : {'hcl', 'hpl', 'hsl', 'hsv', 'rgb'}, optional
            The colorspace for the output channel values.

        Returns
        -------
        color : 3-tuple
            Tuple of channel values for the colorspace `space`.

        See also
        --------
        to_hex
        to_rgb
        to_rgba
        to_xyza
        """
        return self.to_xyza(color, space)[:3]

    @docstring._snippet_manager
    def to_xyza(self, color, space='hcl'):
        """
        Translate color in *any* format to a tuple of channel values in *any*
        colorspace. This is the inverse of `to_rgba`.

        Parameters
        ----------
        %(colors.color)s
        space : {'hcl', 'hpl', 'hsl', 'hsv', 'rgb'}, optional
            The colorspace for the output channel values.

        Returns
        -------
        color : 3-tuple
            Tuple of channel values for the colorspace `space`.

        See also
        --------
        to_hex
        to_rgb
        to_rgba
        to_xyz
        """
        # Run tuple conversions
        # NOTE: Don't pass color tuple, because we may want to permit
        # out-of-bounds RGB values to invert conversion
        *color, opacity = self.to_rgba(color)
        if space == 'rgb':
            pass
        elif space == 'hsv':
            color = hsluv.rgb_to_hsl(*color)  # rgb_to_hsv would also work
        elif space == 'hcl':
            color = hsluv.rgb_to_hcl(*color)
        elif space == 'hsl':
            color = hsluv.rgb_to_hsluv(*color)
        elif space == 'hpl':
            color = hsluv.rgb_to_hpluv(*color)
        else:
            raise ValueError(f'Invalid colorspace {space}.')
        return (*color, opacity)


class _Colormap(object):
    """
    Mixin class used to add some helper methods.
    """
    def _get_data(self, ext, alpha=True):
        """
        Return a string containing the colormap colors for saving.

        Parameters
        ----------
        ext : {'hex', 'txt', 'rgb'}
            The filename extension.
        alpha : bool, optional
            Whether to include an opacity column.
        """
        # Get lookup table colors and filter out bad ones
        if not self._isinit:
            self._init()
        colors = self._lut[:-3, :]

        # Get data string
        if ext == 'hex':
            data = ', '.join(mcolors.convert_hex(color) for color in colors)
        elif ext in ('txt', 'rgb'):
            rgb = mcolors.convert_rgba if alpha else mcolors.convert_rgb
            data = [rgb(color) for color in colors]
            data = '\n'.join(' '.join(f'{num:0.6f}' for num in line) for line in data)
        else:
            raise ValueError(
                f'Invalid extension {ext!r}. Options are: '
                "'hex', 'txt', 'rgb', 'rgba'."
            )
        return data

    def _make_name(self, suffix=None):
        """
        Generate a default colormap name. Do not append more than one
        leading underscore or more than one identical suffix.
        """
        name = self.name
        name = name or ''
        if name[:1] != '_':
            name = '_' + name
        suffix = suffix or 'copy'
        suffix = '_' + suffix
        if name[-len(suffix):] != suffix:
            name = name + suffix
        return name

    def _parse_path(self, path, ext=None, subfolder=None):
        """
        Parse the user input path.

        Parameters
        ----------
        path : path-like, optional
            The file path.
        ext : str
            The default extension.
        subfolder : str, optional
            The subfolder.
        """
        # Get the folder
        folder = rc.user_folder(subfolder=subfolder)
        if path is not None:
            path = os.path.expanduser(path or '.')  # interpret empty string as '.'
            if os.path.isdir(path):
                folder, path = path, None
        # Get the filename
        if path is None:
            path = os.path.join(folder, self.name)
        if not os.path.splitext(path)[1]:
            path = path + '.' + ext  # default file extension
        return path

    @staticmethod
    def _pop_args(*args, names=None, **kwargs):
        """
        Pop the name as a first positional argument or keyword argument.
        Supports matplotlib-style ``Colormap(name, data, N)`` input
        algongside more intuitive ``Colormap(data, name, N)`` input.
        """
        names = names or ()
        if isinstance(names, str):
            names = (names,)
        names = ('name', *names)
        args, kwargs = _kwargs_to_args(names, *args, **kwargs)
        if args[0] is not None and args[1] is None:
            args[:2] = (None, args[0])
        if args[0] is None:
            args[0] = DEFAULT_NAME
        return (*args, kwargs)

    @classmethod
    def _from_file(cls, path, warn_on_failure=False):
        """
        Read generalized colormap and color cycle files.
        """
        path = os.path.expanduser(path)
        name, ext = os.path.splitext(os.path.basename(path))
        listed = issubclass(cls, mcolors.ListedColormap)
        reversed = name[-2:] == '_r'

        # Warn if loading failed during `register_cmaps` or `register_cycles`
        # but raise error if user tries to load a file.
        def _warn_or_raise(descrip, error=RuntimeError):
            prefix = f'Failed to load colormap or color cycle file {path!r}.'
            if warn_on_failure:
                warnings._warn_proplot(prefix + ' ' + descrip)
            else:
                raise error(prefix + ' ' + descrip)
        if not os.path.isfile(path):
            return _warn_or_raise('File not found.', FileNotFoundError)

        # Directly read segmentdata json file
        # NOTE: This is special case! Immediately return name and cmap
        ext = ext[1:]
        if ext == 'json':
            if listed:
                return _warn_or_raise('Cannot load cycles from JSON files.')
            try:
                with open(path, 'r') as fh:
                    data = json.load(fh)
            except json.JSONDecodeError:
                return _warn_or_raise('JSON decoding error.', json.JSONDecodeError)
            kw = {}
            for key in ('cyclic', 'gamma', 'gamma1', 'gamma2', 'space'):
                if key in data:
                    kw[key] = data.pop(key, None)
            if 'red' in data:
                cmap = ContinuousColormap(name, data)
            else:
                cmap = PerceptualColormap(name, data, **kw)
            if reversed:
                cmap = cmap.reversed(name[:-2])
            return cmap

        # Read .rgb and .rgba files
        if ext in ('txt', 'rgb'):
            # Load file
            # NOTE: This appears to be biggest import time bottleneck! Increases
            # time from 0.05s to 0.2s, with numpy loadtxt or with this regex thing.
            delim = re.compile(r'[,\s]+')
            data = [
                delim.split(line.strip())
                for line in open(path)
                if line.strip() and line.strip()[0] != '#'
            ]
            try:
                data = [[float(num) for num in line] for line in data]
            except ValueError:
                return _warn_or_raise(
                    'Expected a table of comma or space-separated floats.'
                )
            # Build x-coordinates and standardize shape
            data = np.array(data)
            if data.shape[1] not in (3, 4):
                return _warn_or_raise(
                    f'Expected 3 or 4 columns of floats. Got {data.shape[1]} columns.'
                )
            if ext[0] != 'x':  # i.e. no x-coordinates specified explicitly
                x = np.linspace(0, 1, data.shape[0])
            else:
                x, data = data[:, 0], data[:, 1:]

        # Load XML files created with scivizcolor
        # Adapted from script found here:
        # https://sciviscolor.org/matlab-matplotlib-pv44/
        elif ext == 'xml':
            try:
                doc = ElementTree.parse(path)
            except ElementTree.ParseError:
                return _warn_or_raise('XML parsing error.', ElementTree.ParseError)
            x, data = [], []
            for s in doc.getroot().findall('.//Point'):
                # Verify keys
                if any(key not in s.attrib for key in 'xrgb'):
                    return _warn_or_raise(
                        'Missing an x, r, g, or b key inside one or more <Point> tags.'
                    )
                # Get data
                color = []
                for key in 'rgbao':  # o for opacity
                    if key not in s.attrib:
                        continue
                    color.append(float(s.attrib[key]))
                x.append(float(s.attrib['x']))
                data.append(color)
            # Convert to array
            if not all(
                len(data[0]) == len(color) and len(color) in (3, 4) for color in data
            ):
                return _warn_or_raise(
                    'Unexpected channel number or mixed channels across <Point> tags.'
                )

        # Read hex strings
        elif ext == 'hex':
            # Read arbitrary format
            string = open(path).read()  # into single string
            data = REGEX_HEX_MULTI.findall(string)
            if len(data) < 2:
                return _warn_or_raise(
                    'Failed to find 6-digit or 8-digit HEX strings.'
                )
            # Convert to array
            x = np.linspace(0, 1, len(data))
            data = [Color(color) for color in data]

        # Invalid extension
        else:
            return _warn_or_raise(
                'Unknown colormap file extension {ext!r}. Options are: '
                + ', '.join(map(repr, ('json', 'txt', 'rgb', 'hex')))
                + '.'
            )

        # Standardize and reverse if necessary to cmap
        # TODO: Document the fact that filenames ending in _r return a reversed
        # version of the colormap stored in that file.
        x = np.array(x)
        x = (x - x.min()) / (x.max() - x.min())  # ensure they span 0-1
        data = np.array(data)
        if np.any(data > 2):  # from 0-255 to 0-1
            data = data / 255
        if reversed:
            name = name[:-2]
            data = data[::-1, :]
            x = 1 - x[::-1]
        if listed:
            return DiscreteColormap(data, name)
        else:
            data = [(x, color) for x, color in zip(x, data)]
            return ContinuousColormap.from_list(name, data)


class ContinuousColormap(mcolors.LinearSegmentedColormap, _Colormap):
    r"""
    Replacement for `~matplotlib.colors.LinearSegmentedColormap`.
    """
    def __str__(self):
        return type(self).__name__ + f'(name={self.name!r})'

    def __repr__(self):
        string = f" 'name': {self.name!r},\n"
        if hasattr(self, '_space'):
            string += f" 'space': {self._space!r},\n"
        if hasattr(self, '_cyclic'):
            string += f" 'cyclic': {self._cyclic!r},\n"
        for key, data in self._segmentdata.items():
            if callable(data):
                string += f' {key!r}: <function>,\n'
            else:
                stop = data[-1][1]
                start = data[0][2]
                string += f' {key!r}: [{start:.2f}, ..., {stop:.2f}],\n'
        return type(self).__name__ + '({\n' + string + '})'

    @docstring._snippet_manager
    def __init__(self, *args, gamma=1, alpha=None, cyclic=False, **kwargs):
        """
        Parameters
        ----------
        segmentdata : dict-like
            Dictionary containing the keys ``'red'``, ``'green'``, ``'blue'``, and
            (optionally) ``'alpha'``. The shorthands ``'r'``, ``'g'``, ``'b'``,
            and ``'a'`` are also acceptable. The key values can be callable
            functions that return channel values given a colormap index, or
            3-column arrays indicating the coordinates and channel transitions. See
            `matplotlib.colors.LinearSegmentedColormap` for a detailed explanation.
        %(colors.name)s
        %(colors.N)s
        gamma : float, optional
            Gamma scaling used for the *x* coordinates.
        %(colors.alpha)s
        %(colors.cyclic)s

        Other parameters
        ----------------
        **kwargs
            Passed to `matplotlib.colors.LinearSegmentedColormap`.

        See also
        --------
        DiscreteColormap
        matplotlib.colors.LinearSegmentedColormap
        proplot.constructor.Colormap
        """
        # NOTE: Additional keyword args should raise matplotlib error
        name, segmentdata, N, kwargs = self._pop_args(
            *args, names=('segmentdata', 'N'), **kwargs
        )
        if not isinstance(segmentdata, dict):
            raise ValueError(f'Invalid segmentdata {segmentdata}. Must be a dict.')
        N = _not_none(N, rc['image.lut'])
        data = _pop_props(segmentdata, 'rgba', 'hsla')
        if segmentdata:
            raise ValueError(f'Invalid segmentdata keys {tuple(segmentdata)}.')
        super().__init__(name, data, N=N, gamma=gamma, **kwargs)
        self._cyclic = cyclic
        if alpha is not None:
            self.set_alpha(alpha)

    def append(self, *args, ratios=None, name=None, N=None, **kwargs):
        """
        Return the concatenation of this colormap with the
        input colormaps.

        Parameters
        ----------
        *args
            Instances of `ContinuousColormap`.
        ratios : sequence of float, optional
            Relative extent of each component colormap in the
            merged colormap. Length must equal ``len(args) + 1``.
            For example, ``cmap1.append(cmap2, ratios=(2, 1))`` generates
            a colormap with the left two-thrids containing colors from
            ``cmap1`` and the right one-third containing colors from ``cmap2``.
        name : str, optional
            The colormap name. Default is to merge each name with underscores and
            prepend a leading underscore, for example ``_name1_name2``.
        N : int, optional
            The number of points in the colormap lookup table. Default is
            to sum the length of each lookup table.

        Other parameters
        ----------------
        **kwargs
            Passed to `ContinuousColormap.copy`
            or `PerceptualColormap.copy`.

        Returns
        -------
        ContinuousColormap
            The colormap.

        See also
        --------
        DiscreteColormap.append
        """
        # Parse input args
        if not args:
            return self
        if not all(isinstance(cmap, mcolors.LinearSegmentedColormap) for cmap in args):
            raise TypeError(f'Arguments {args!r} must be LinearSegmentedColormaps.')

        # PerceptualColormap --> ContinuousColormap conversions
        cmaps = [self, *args]
        spaces = {getattr(cmap, '_space', None) for cmap in cmaps}
        to_continuous = len(spaces) > 1  # mixed colorspaces *or* mixed types
        if to_continuous:
            for i, cmap in enumerate(cmaps):
                if isinstance(cmap, PerceptualColormap):
                    cmaps[i] = cmap.to_continuous()

        # Combine the segmentdata, and use the y1/y2 slots at merge points so
        # we never interpolate between end colors of different colormaps
        segmentdata = {}
        if name is None:
            name = '_' + '_'.join(cmap.name for cmap in cmaps)
        if not np.iterable(ratios):
            ratios = [1] * len(cmaps)
        ratios = np.asarray(ratios) / np.sum(ratios)
        x0 = np.append(0, np.cumsum(ratios))  # coordinates for edges
        xw = x0[1:] - x0[:-1]  # widths between edges
        for key in cmaps[0]._segmentdata.keys():  # not self._segmentdata
            # Callable segments
            # WARNING: If just reference a global 'funcs' list from inside the
            # 'data' function it can get overwritten in this loop. Must
            # embed 'funcs' into the definition using a keyword argument.
            datas = [cmap._segmentdata[key] for cmap in cmaps]
            if all(map(callable, datas)):  # expand range from x-to-w to 0-1
                def xyy(ix, funcs=datas):  # noqa: E306
                    ix = np.atleast_1d(ix)
                    kx = np.empty(ix.shape)
                    for j, jx in enumerate(ix.flat):
                        idx = max(np.searchsorted(x0, jx) - 1, 0)
                        kx.flat[j] = funcs[idx]((jx - x0[idx]) / xw[idx])
                    return kx

            # Concatenate segment arrays and make the transition at the
            # seam instant so we *never interpolate* between end colors
            # of different maps.
            elif not any(map(callable, datas)):
                datas = []
                for x, w, cmap in zip(x0[:-1], xw, cmaps):
                    xyy = np.array(cmap._segmentdata[key])
                    xyy[:, 0] = x + w * xyy[:, 0]
                    datas.append(xyy)
                for i in range(len(datas) - 1):
                    datas[i][-1, 2] = datas[i + 1][0, 2]
                    datas[i + 1] = datas[i + 1][1:, :]
                xyy = np.concatenate(datas, axis=0)
                xyy[:, 0] = xyy[:, 0] / xyy[:, 0].max(axis=0)  # fix fp errors

            else:
                raise TypeError(
                    'Cannot merge colormaps with mixed callable '
                    'and non-callable segment data.'
                )
            segmentdata[key] = xyy

            # Handle gamma values
            ikey = None
            if key == 'saturation':
                ikey = 'gamma1'
            elif key == 'luminance':
                ikey = 'gamma2'
            if not ikey or ikey in kwargs:
                continue
            gamma = []
            callable_ = all(map(callable, datas))
            for cmap in cmaps:
                igamma = getattr(cmap, '_' + ikey)
                if not np.iterable(igamma):
                    if callable_:
                        igamma = (igamma,)
                    else:
                        igamma = (igamma,) * (len(cmap._segmentdata[key]) - 1)
                gamma.extend(igamma)
            if callable_:
                if any(igamma != gamma[0] for igamma in gamma[1:]):
                    warnings._warn_proplot(
                        'Cannot use multiple segment gammas when concatenating '
                        f'callable segments. Using the first gamma of {gamma[0]}.'
                    )
                gamma = gamma[0]
            kwargs[ikey] = gamma

        # Return copy or merge mixed types
        if to_continuous and isinstance(self, PerceptualColormap):
            return ContinuousColormap(name, segmentdata, N, **kwargs)
        else:
            return self.copy(name, segmentdata, N, **kwargs)

    def cut(self, cut=None, name=None, left=None, right=None, **kwargs):
        """
        Return a version of the colormap with the center "cut out".
        This is great for making the transition from "negative" to "positive"
        in a diverging colormap more distinct.

        Parameters
        ----------
        cut : float, optional
            The proportion to cut from the center of the colormap. For example,
            ``cut=0.1`` cuts the central 10%, or ``cut=-0.1`` fills the ctranl 10%
            of the colormap with the current central color (usually white).
        name : str, default: '_name_copy'
            The new colormap name.
        left, right : float, default: 0, 1
            The colormap indices for the "leftmost" and "rightmost"
            colors. See `~ContinuousColormap.truncate` for details.
        right : float, optional
            The colormap index for the new "rightmost" color. Must fall between

        Other parameters
        ----------------
        **kwargs
            Passed to `ContinuousColormap.copy` or `PerceptualColormap.copy`.

        Returns
        -------
        ContinuousColormap
            The colormap.

        See also
        --------
        ContinuousColormap.truncate
        DiscreteColormap.truncate
        """
        # Parse input args
        left = max(_not_none(left, 0), 0)
        right = min(_not_none(right, 1), 1)
        cut = _not_none(cut, 0)
        offset = 0.5 * cut
        if offset < 0:  # add extra 'white' later on
            offset = 0
        elif offset == 0:
            return self.truncate(left, right)

        # Decompose cut into two truncations followed by concatenation
        if 0.5 - offset < left or 0.5 + offset > right:
            raise ValueError(f'Invalid cut={cut} for left={left} and right={right}.')
        if name is None:
            name = self._make_name()
        cmap_left = self.truncate(left, 0.5 - offset)
        cmap_right = self.truncate(0.5 + offset, right)

        # Permit adding extra 'white' to colormap center
        # NOTE: Rely on channel abbreviations to simplify code here
        args = []
        if cut < 0:
            ratio = 0.5 - 0.5 * abs(cut)  # ratio for flanks on either side
            space = getattr(self, '_space', None) or 'rgb'
            xyza = to_xyza(self(0.5), space=space)
            segmentdata = {
                key: _make_segment_data(x) for key, x in zip(space + 'a', xyza)
            }
            args.append(type(self)(DEFAULT_NAME, segmentdata, self.N))
            kwargs.setdefault('ratios', (ratio, abs(cut), ratio))
        args.append(cmap_right)

        return cmap_left.append(*args, name=name, **kwargs)

    def reversed(self, name=None, **kwargs):
        """
        Return a reversed copy of the colormap.

        Parameters
        ----------
        name : str, default: '_name_r'
            The new colormap name.

        Other parameters
        ----------------
        **kwargs
            Passed to `ContinuousColormap.copy`
            or `PerceptualColormap.copy`.

        See also
        --------
        matplotlib.colors.LinearSegmentedColormap.reversed
        """
        # Reverse segments
        segmentdata = {
            key: (
                (lambda x, func=data: func(x))
                if callable(data) else
                [(1.0 - x, y1, y0) for x, y0, y1 in reversed(data)]
            )
            for key, data in self._segmentdata.items()
        }

        # Reverse gammas
        if name is None:
            name = self._make_name(suffix='r')
        for key in ('gamma1', 'gamma2'):
            if key in kwargs:
                continue
            gamma = getattr(self, '_' + key, None)
            if gamma is not None and np.iterable(gamma):
                kwargs[key] = gamma[::-1]

        cmap = self.copy(name, segmentdata, **kwargs)
        cmap._rgba_under, cmap._rgba_over = cmap._rgba_over, cmap._rgba_under
        return cmap

    @docstring._snippet_manager
    def save(self, path=None, alpha=True):
        """
        Save the colormap data to a file.

        Parameters
        ----------
        path : path-like, optional
            The output filename. If not provided, the colormap is saved in the
            ``cmaps`` subfolder in `~proplot.config.Configurator.user_folder`
            under the filename ``name.json`` (where ``name`` is the colormap
            name). Valid extensions are shown in the below table.

        %(rc.cmap_exts)s

        alpha : bool, optional
            Whether to include an opacity column for ``.rgb``
            and ``.txt`` files.

        See also
        --------
        DiscreteColormap.save
        """
        # NOTE: We sanitize segmentdata before saving to json. Convert numpy float to
        # builtin float, np.array to list of lists, and callable to list of lists.
        # We tried encoding func.__code__ with base64 and marshal instead, but when
        # cmap.append() embeds functions as keyword arguments, this seems to make it
        # *impossible* to load back up the function with FunctionType (error message:
        # arg 5 (closure) must be tuple). Instead use this brute force workaround.
        filename = self._parse_path(path, ext='json', subfolder='cmaps')
        _, ext = os.path.splitext(filename)
        if ext[1:] != 'json':
            # Save lookup table colors
            data = self._get_data(ext[1:], alpha=alpha)
            with open(filename, 'w') as fh:
                fh.write(data)
        else:
            # Save segment data itself
            data = {}
            for key, value in self._segmentdata.items():
                if callable(value):
                    x = np.linspace(0, 1, rc['image.lut'])  # just save the transitions
                    y = np.array([value(_) for _ in x]).squeeze()
                    value = np.vstack((x, y, y)).T
                data[key] = np.asarray(value).astype(float).tolist()
            keys = ()
            if isinstance(self, PerceptualColormap):
                keys = ('cyclic', 'gamma1', 'gamma2', 'space')
            elif isinstance(self, ContinuousColormap):
                keys = ('cyclic', 'gamma')
            for key in keys:  # add all attrs to dictionary
                data[key] = getattr(self, '_' + key)
            with open(filename, 'w') as fh:
                json.dump(data, fh, indent=4)
        print(f'Saved colormap to {filename!r}.')

    def set_alpha(self, alpha, coords=None, ratios=None):
        """
        Set the opacity for the entire colormap or set up an opacity gradation.

        Parameters
        ----------
        alpha : float or sequence of float
            If float, this is the opacity for the entire colormap. If sequence of
            float, the colormap traverses these opacity values.
        coords : sequence of float, optional
            Colormap coordinates for the opacity values. The first and last
            coordinates must be ``0`` and ``1``. If `alpha` is not scalar, the
            default coordinates are ``np.linspace(0, 1, len(alpha))``.
        ratios : sequence of float, optional
            Relative extent of each opacity transition segment. Length should
            equal ``len(alpha) + 1``. For example
            ``cmap.set_alpha((1, 1, 0), ratios=(2, 1))`` creates a transtion from
            100 percent to 0 percent opacity in the right *third* of the colormap.

        See also
        --------
        DiscreteColormap.set_alpha
        """
        alpha = _make_segment_data(alpha, coords=coords, ratios=ratios)
        self._segmentdata['alpha'] = alpha
        self._isinit = False

    def set_cyclic(self, b):
        """
        Set whether this colormap is "cyclic". See `ContinuousColormap` for details.
        """
        self._cyclic = bool(b)
        self._isinit = False

    def shifted(self, shift=180, name=None, **kwargs):
        """
        Return a cyclicaly shifted version of the colormap. If the colormap
        cyclic property is set to ``False`` a warning will be raised.

        Parameters
        ----------
        shift : float, default: 180
            The number of degrees to shift, out of 360 degrees.
        name : str, default: '_name_s'
            The new colormap name.

        Other parameters
        ----------------
        **kwargs
            Passed to `ContinuousColormap.copy` or `PerceptualColormap.copy`.

        See also
        --------
        DiscreteColormap.shifted
        """
        shift = shift or 0
        shift %= 360
        shift /= 360
        if shift == 0:
            return self
        if name is None:
            name = self._make_name(suffix='s')
        if not self._cyclic:
            warnings._warn_proplot(
                f'Shifting non-cyclic colormap {self.name!r}. To suppress this '
                'warning use cmap.set_cyclic(True) or Colormap(..., cyclic=True).'
            )
            self._cyclic = True
        ratios = (1 - shift, shift)
        cmap_left = self.truncate(shift, 1)
        cmap_right = self.truncate(0, shift)
        return cmap_left.append(cmap_right, ratios=ratios, name=name, **kwargs)

    def truncate(self, left=None, right=None, name=None, **kwargs):
        """
        Return a truncated version of the colormap.

        Parameters
        ----------
        left : float, default: 0
            The colormap index for the new "leftmost" color. Must fall between ``0``
            and ``1``. For example, ``left=0.1`` cuts the leftmost 10%% of the colors.
        right : float, default: 1
            The colormap index for the new "rightmost" color. Must fall between ``0``
            and ``1``. For example, ``right=0.9`` cuts the leftmost 10%% of the colors.
        name : str, default: '_name_copy'
            The new colormap name.

        Other parameters
        ----------------
        **kwargs
            Passed to `ContinuousColormap.copy`
            or `PerceptualColormap.copy`.

        See also
        --------
        DiscreteColormap.truncate
        """
        # Bail out
        left = max(_not_none(left, 0), 0)
        right = min(_not_none(right, 1), 1)
        if left == 0 and right == 1:
            return self
        if name is None:
            name = self._make_name()

        # Resample the segmentdata arrays
        segmentdata = {}
        for key, data in self._segmentdata.items():
            # Callable array
            # WARNING: If just reference a global 'xyy' callable from inside
            # the lambda function it gets overwritten in the loop! Must embed
            # the old callable in the new one as a default keyword arg.
            if callable(data):
                def xyy(x, func=data):
                    return func(left + x * (right - left))

            # Slice
            # l is the first point where x > 0 or x > left, should be >= 1
            # r is the last point where r < 1 or r < right
            else:
                xyy = np.asarray(data)
                x = xyy[:, 0]
                l = np.searchsorted(x, left)  # first x value > left  # noqa
                r = np.searchsorted(x, right) - 1  # last x value < right
                xc = xyy[l:r + 1, :].copy()
                xl = xyy[l - 1, 1:] + (left - x[l - 1]) * (
                    (xyy[l, 1:] - xyy[l - 1, 1:]) / (x[l] - x[l - 1])
                )
                xr = xyy[r, 1:] + (right - x[r]) * (
                    (xyy[r + 1, 1:] - xyy[r, 1:]) / (x[r + 1] - x[r])
                )
                xyy = np.vstack(((left, *xl), xc, (right, *xr)))
                xyy[:, 0] = (xyy[:, 0] - left) / (right - left)

            # Retain the corresponding gamma *segments*
            segmentdata[key] = xyy
            if key == 'saturation':
                ikey = 'gamma1'
            elif key == 'luminance':
                ikey = 'gamma2'
            else:
                continue
            if ikey in kwargs:
                continue
            gamma = getattr(self, '_' + ikey)
            if np.iterable(gamma):
                if callable(xyy):
                    if any(igamma != gamma[0] for igamma in gamma[1:]):
                        warnings._warn_proplot(
                            'Cannot use multiple segment gammas when '
                            'truncating colormap. Using the first gamma '
                            f'of {gamma[0]}.'
                        )
                    gamma = gamma[0]
                else:
                    igamma = gamma[l - 1:r + 1]
                    if len(igamma) == 0:  # TODO: issue warning?
                        gamma = gamma[0]
                    else:
                        gamma = igamma
            kwargs[ikey] = gamma

        return self.copy(name, segmentdata, **kwargs)

    def copy(
        self, name=None, segmentdata=None, N=None, *,
        alpha=None, gamma=None, cyclic=None
    ):
        """
        Return a new colormap with relevant properties copied from this one
        if they were not provided as keyword arguments.

        Parameters
        ----------
        name : str, default: '_name_copy'
            The new colormap name.
        segmentdata, N, alpha, gamma, cyclic : optional
            See `ContinuousColormap`. If not provided, these are copied
            from the current colormap.

        See also
        --------
        DiscreteColormap.copy
        PerceptualColormap.copy
        """
        if name is None:
            name = self._make_name()
        if segmentdata is None:
            segmentdata = self._segmentdata.copy()
        if gamma is None:
            gamma = self._gamma
        if cyclic is None:
            cyclic = self._cyclic
        if N is None:
            N = self.N
        cmap = ContinuousColormap(
            name, segmentdata, N,
            alpha=alpha, gamma=gamma, cyclic=cyclic
        )
        cmap._rgba_bad = self._rgba_bad
        cmap._rgba_under = self._rgba_under
        cmap._rgba_over = self._rgba_over
        return cmap

    def to_discrete(self, samples=10, name=None, **kwargs):
        """
        Convert the `ContinuousColormap` to a `DiscreteColormap` by drawing
        samples from the colormap.

        Parameters
        ----------
        samples : int or sequence of float, optional
            If integer, draw samples at the colormap coordinates
            ``np.linspace(0, 1, samples)``. If sequence of float,
            draw samples at the specified points.
        name : str, default: '_name_copy'
            The new colormap name.

        Other parameters
        ----------------
        **kwargs
            Passed to `DiscreteColormap`.

        See also
        --------
        PerceptualColormap.to_continuous
        """
        if isinstance(samples, Integral):
            samples = np.linspace(0, 1, samples)
        elif not np.iterable(samples):
            raise TypeError('Samples must be integer or iterable.')
        samples = np.asarray(samples)
        colors = self(samples)
        if name is None:
            name = self._make_name()
        return DiscreteColormap(colors, name=name, **kwargs)

    @classmethod
    @docstring._snippet_manager
    def from_file(cls, path, *, warn_on_failure=False):
        """
        Load colormap from a file.

        Parameters
        ----------
        path : path-like
            The file path. Valid file extensions are shown in the below table.

        %(rc.cmap_exts)s

        warn_on_failure : bool, optional
            If ``True``, issue a warning when loading fails instead of
            raising an error.

        See also
        --------
        DiscreteColormap.from_file
        """
        return cls._from_file(path, warn_on_failure=warn_on_failure)

    @classmethod
    @docstring._snippet_manager
    def from_list(cls, *args, **kwargs):
        """
        Make a `ContinuousColormap` from a sequence of colors.

        Parameters
        ----------
        %(colors.from_list)s

        Other parameters
        ----------------
        **kwargs
            Passed to `ContinuousColormap`.

        Returns
        -------
        ContinuousColormap
            The colormap.

        See also
        --------
        matplotlib.colors.LinearSegmentedColormap.from_list
        PerceptualColormap.from_list
        """
        # Get coordinates
        name, colors, ratios, kwargs = cls._pop_args(
            *args, names=('colors', 'ratios'), **kwargs
        )
        coords = None
        if not np.iterable(colors):
            raise TypeError('Colors must be iterable.')
        if (
            np.iterable(colors[0])
            and len(colors[0]) == 2
            and not isinstance(colors[0], str)
        ):
            coords, colors = zip(*colors)
        colors = [to_rgba(color) for color in colors]

        # Build segmentdata
        keys = ('red', 'green', 'blue', 'alpha')
        cdict = {}
        for key, values in zip(keys, zip(*colors)):
            cdict[key] = _make_segment_data(values, coords, ratios)
        return cls(name, cdict, **kwargs)

    # Deprecated
    to_listed = warnings._rename_objs(
        '0.8.0',
        to_listed=to_discrete
    )
    concatenate, punched, truncated, updated = warnings._rename_objs(
        '0.6.0',
        concatenate=append,
        punched=cut,
        truncated=truncate,
        updated=copy,
    )


class DiscreteColormap(mcolors.ListedColormap, _Colormap):
    r"""
    Replacement for `~matplotlib.colors.ListedColormap`.
    """
    def __str__(self):
        return f'DiscreteColormap(name={self.name!r})'

    def __repr__(self):
        colors = [c if isinstance(c, str) else to_hex(c) for c in self.colors]
        string = 'DiscreteColormap({\n'
        string += f" 'name': {self.name!r},\n"
        string += f" 'colors': {colors!r},\n"
        string += '})'
        return string

    def __init__(self, colors, name=None, N=None, alpha=None, **kwargs):
        """
        Parameters
        ----------
        colors : sequence of color-spec, optional
            The colormap colors.
        name : str, default: '_no_name'
            The colormap name.
        N : int, default: ``len(colors)``
            The number of levels. The color list is truncated or wrapped
            to match this length.
        alpha : float, optional
            The opacity for the colormap colors. This overrides the
            input color opacities.

        Other parameters
        ----------------
        **kwargs
            Passed to `~matplotlib.colors.ListedColormap`.

        See also
        --------
        ContinuousColormap
        matplotlib.colors.ListedColormap
        proplot.constructor.Colormap
        """
        # NOTE: This also improves 'monochrome' detection to test all items
        # in the list. Otherwise ContourSet does not apply negative_linestyle
        # to monochromatic colormaps generated by passing a 'colors' keyword.
        # Also note that under the hood, just like proplot, ContourSet builds
        # identical monochromatic ListedColormaps when it receives scalar colors.
        N = _not_none(N, len(colors))
        name = _not_none(name, DEFAULT_NAME)
        super().__init__(colors, name=name, N=N, **kwargs)
        if alpha is not None:
            self.set_alpha(alpha)
        for i, color in enumerate(self.colors):
            if isinstance(color, np.ndarray):
                self.colors[i] = color.tolist()
        if self.colors and all(self.colors[0] == color for color in self.colors):
            self.monochrome = True  # for contour negative dash style

    def append(self, *args, name=None, N=None, **kwargs):
        """
        Append arbitrary colormaps onto this colormap.

        Parameters
        ----------
        *args
            Instances of `DiscreteColormap`.
        name : str, optional
            The new colormap name. Default is to merge each name with underscores and
            prepend a leading underscore, for example ``_name1_name2``.
        N : int, optional
            The number of points in the colormap lookup table. Default is
            the number of colors in the concatenated lists.

        Other parameters
        ----------------
        **kwargs
            Passed to `~DiscreteColormap.copy`.

        See also
        --------
        ContinuousColormap.append
        """
        if not args:
            return self
        if not all(isinstance(cmap, mcolors.ListedColormap) for cmap in args):
            raise TypeError(f'Arguments {args!r} must be DiscreteColormap.')
        cmaps = (self, *args)
        if name is None:
            name = '_' + '_'.join(cmap.name for cmap in cmaps)
        colors = [color for cmap in cmaps for color in cmap.colors]
        N = _not_none(N, len(colors))
        return self.copy(colors, name, N, **kwargs)

    @docstring._snippet_manager
    def save(self, path=None, alpha=True):
        """
        Save the colormap data to a file.

        Parameters
        ----------
        path : path-like, optional
            The output filename. If not provided, the colormap is saved in the
            ``cycles`` subfolder in `~proplot.config.Configurator.user_folder`
            under the filename ``name.hex`` (where ``name`` is the color cycle
            name). Valid extensions are described in the below table.

        %(rc.cycle_exts)s

        alpha : bool, optional
            Whether to include an opacity column for ``.rgb``
            and ``.txt`` files.

        See also
        --------
        ContinuousColormap.save
        """
        filename = self._parse_path(path, ext='hex', subfolder='cycles')
        _, ext = os.path.splitext(filename)
        data = self._get_data(ext[1:], alpha=alpha)
        with open(filename, 'w') as fh:
            fh.write(data)
        print(f'Saved colormap to {filename!r}.')

    def set_alpha(self, alpha):
        """
        Set the opacity for the entire colormap.

        Parameters
        ----------
        alpha : float
            The opacity.

        See also
        --------
        ContinuousColormap.set_alpha
        """
        self.colors = [set_alpha(color, alpha) for color in self.colors]
        self._init()

    def reversed(self, name=None, **kwargs):
        """
        Return a reversed version of the colormap.

        Parameters
        ----------
        name : str, default: '_name_r'
            The new colormap name.

        Other parameters
        ----------------
        **kwargs
            Passed to `DiscreteColormap.copy`

        See also
        --------
        matplotlib.colors.ListedColormap.reversed
        """
        if name is None:
            name = self._make_name(suffix='r')
        colors = self.colors[::-1]
        cmap = self.copy(colors, name, **kwargs)
        cmap._rgba_under, cmap._rgba_over = cmap._rgba_over, cmap._rgba_under
        return cmap

    def shifted(self, shift=1, name=None):
        """
        Return a cyclically shifted version of the colormap.

        Parameters
        ----------
        shift : float, default: 1
            The number of list indices to shift.
        name : str, eefault: '_name_s'
            The new colormap name.

        See also
        --------
        ContinuousColormap.shifted
        """
        if not shift:
            return self
        if name is None:
            name = self._make_name(suffix='s')
        shift = shift % len(self.colors)
        colors = list(self.colors)
        colors = colors[shift:] + colors[:shift]
        return self.copy(colors, name, len(colors))

    def truncate(self, left=None, right=None, name=None):
        """
        Return a truncated version of the colormap.

        Parameters
        ----------
        left : float, default: None
            The colormap index for the new "leftmost" color. Must fall between ``0``
            and ``self.N``. For example, ``left=2`` drops the first two colors.
        right : float, default: None
            The colormap index for the new "rightmost" color. Must fall between ``0``
            and ``self.N``. For example, ``right=4`` keeps the first four colors.
        name : str, default: '_name_copy'
            The new colormap name.

        See also
        --------
        ContinuousColormap.truncate
        """
        if left is None and right is None:
            return self
        if name is None:
            name = self._make_name()
        colors = self.colors[left:right]
        return self.copy(colors, name, len(colors))

    def copy(self, colors=None, name=None, N=None, *, alpha=None):
        """
        Return a new colormap with relevant properties copied from this one
        if they were not provided as keyword arguments.

        Parameters
        ----------
        name : str, default: '_name_copy'
            The new colormap name.
        colors, N, alpha : optional
            See `DiscreteColormap`. If not provided,
            these are copied from the current colormap.

        See also
        --------
        ContinuousColormap.copy
        PerceptualColormap.copy
        """
        if name is None:
            name = self._make_name()
        if colors is None:
            colors = list(self.colors)  # copy
        if N is None:
            N = self.N
        cmap = DiscreteColormap(colors, name, N=N, alpha=alpha)
        cmap._rgba_bad = self._rgba_bad
        cmap._rgba_under = self._rgba_under
        cmap._rgba_over = self._rgba_over
        return cmap

    @classmethod
    @docstring._snippet_manager
    def from_file(cls, path, *, warn_on_failure=False):
        """
        Load color cycle from a file.

        Parameters
        ----------
        path : path-like
            The file path. Valid file extensions are shown in the below table.

        %(rc.cycle_exts)s

        warn_on_failure : bool, optional
            If ``True``, issue a warning when loading fails instead of
            raising an error.

        See also
        --------
        ContinuousColormap.from_file
        """
        return cls._from_file(path, warn_on_failure=warn_on_failure)

    # Rename methods
    concatenate, truncated, updated = warnings._rename_objs(
        '0.6.0',
        concatenate=append,
        truncated=truncate,
        updated=copy,
    )


class PerceptualColormap(ContinuousColormap):
    """
    A `ContinuousColormap` with linear transitions across hue, saturation,
    and luminance rather than red, blue, and green.
    """
    @docstring._snippet_manager
    def __init__(
        self, *args, space=None, clip=True, gamma=None, gamma1=None, gamma2=None,
        **kwargs
    ):
        """
        Parameters
        ----------
        segmentdata : dict-like
            Dictionary containing the keys ``'hue'``, ``'saturation'``,
            ``'luminance'``, and (optionally) ``'alpha'``. The key ``'chroma'`` is
            treated as a synonym for ``'saturation'``. The shorthands ``'h'``,
            ``'s'``, ``'l'``, ``'a'``, and ``'c'`` are also acceptable. The key
            values can be callable functions that return channel values given a
            colormap index, or 3-column arrays indicating the coordinates and
            channel transitions. See `~matplotlib.colors.LinearSegmentedColormap`
            for a more detailed explanation.
        %(colors.name)s
        %(colors.N)s
        %(colors.space)s
        clip : bool, optional
            Whether to "clip" impossible colors (i.e. truncate HCL colors with
            RGB channels with values greater than 1) or mask them out as gray.
        %(colors.gamma)s
        %(colors.alpha)s
        %(colors.cyclic)s

        Other parameters
        ----------------
        **kwargs
            Passed to `matploitlib.colors.LinearSegmentedColormap`.

        Example
        -------
        The below example generates a `PerceptualColormap` from a
        `segmentdata` dictionary that uses color names for the hue data,
        instead of channel values between ``0`` and ``360``.

        >>> import proplot as pplt
        >>> data = {
        >>>     'h': [[0, 'red', 'red'], [1, 'blue', 'blue']],
        >>>     's': [[0, 100, 100], [1, 100, 100]],
        >>>     'l': [[0, 100, 100], [1, 20, 20]],
        >>> }
        >>> cmap = pplt.PerceptualColormap(data)

        See also
        --------
        ContinuousColormap
        proplot.constructor.Colormap
        """
        # Checks
        name, segmentdata, N, kwargs = self._pop_args(
            *args, names=('segmentdata', 'N'), **kwargs
        )
        data = _pop_props(segmentdata, 'hsla')
        if segmentdata:
            raise ValueError(f'Invalid segmentdata keys {tuple(segmentdata)}.')
        space = _not_none(space, DEFAULT_SPACE).lower()
        if space not in ('rgb', 'hsv', 'hpl', 'hsl', 'hcl'):
            raise ValueError(f'Unknown colorspace {space!r}.')
        # Convert color strings to channel values
        for key, array in data.items():
            if callable(array):  # permit callable
                continue
            for i, xyy in enumerate(array):
                xyy = list(xyy)  # make copy!
                for j, y in enumerate(xyy[1:]):  # modify the y values
                    xyy[j + 1] = _color_channel(y, key, space)
                data[key][i] = xyy
        # Initialize
        super().__init__(name, data, gamma=1.0, N=N, **kwargs)
        self._gamma1 = _not_none(gamma1, gamma, 1.0)
        self._gamma2 = _not_none(gamma2, gamma, 1.0)
        self._space = space
        self._clip = clip

    def _init(self):
        """
        As with `~matplotlib.colors.LinearSegmentedColormap`, but convert
        each value in the lookup table from ``self._space`` to RGB.
        """
        # First generate the lookup table
        channels = ('hue', 'saturation', 'luminance')
        inverses = (False, False, True)  # weight low chroma, high luminance
        gammas = (1.0, self._gamma1, self._gamma2)
        self._lut_hsl = np.ones((self.N + 3, 4), float)  # fill
        for i, (channel, gamma, inverse) in enumerate(zip(channels, gammas, inverses)):
            self._lut_hsl[:-3, i] = _make_lookup_table(
                self.N, self._segmentdata[channel], gamma, inverse
            )
        if 'alpha' in self._segmentdata:
            self._lut_hsl[:-3, 3] = _make_lookup_table(
                self.N, self._segmentdata['alpha']
            )
        self._lut_hsl[:-3, 0] %= 360

        # Make hues circular, set extremes i.e. copy HSL values
        self._lut = self._lut_hsl.copy()
        self._set_extremes()  # generally just used end values in segmentdata
        self._isinit = True

        # Now convert values to RGB and clip colors
        for i in range(self.N + 3):
            self._lut[i, :3] = to_rgb(self._lut[i, :3], self._space)
        self._lut[:, :3] = _clip_colors(self._lut[:, :3], self._clip)

    @docstring._snippet_manager
    def set_gamma(self, gamma=None, gamma1=None, gamma2=None):
        """
        Set the gamma value(s) for the luminance and saturation transitions.

        Parameters
        ----------
        %(colors.gamma)s
        """
        gamma1 = _not_none(gamma1, gamma)
        gamma2 = _not_none(gamma2, gamma)
        if gamma1 is not None:
            self._gamma1 = gamma1
        if gamma2 is not None:
            self._gamma2 = gamma2
        self._init()

    def copy(
        self, name=None, segmentdata=None, N=None, *,
        alpha=None, gamma=None, cyclic=None,
        clip=None, gamma1=None, gamma2=None, space=None
    ):
        """
        Return a new colormap with relevant properties copied from this one
        if they were not provided as keyword arguments.

        Parameters
        ----------
        name : str, default: '_name_copy'
            The new colormap name.
        segmentdata, N, alpha, clip, cyclic, gamma, gamma1, gamma2, space : optional
            See `PerceptualColormap`. If not provided,
            these are copied from the current colormap.

        See also
        --------
        DiscreteColormap.copy
        ContinuousColormap.copy
        """
        if name is None:
            name = self._make_name()
        if segmentdata is None:
            segmentdata = self._segmentdata.copy()
        if space is None:
            space = self._space
        if clip is None:
            clip = self._clip
        if gamma is not None:
            gamma1 = gamma2 = gamma
        if gamma1 is None:
            gamma1 = self._gamma1
        if gamma2 is None:
            gamma2 = self._gamma2
        if cyclic is None:
            cyclic = self._cyclic
        if N is None:
            N = self.N
        cmap = PerceptualColormap(
            name, segmentdata, N,
            alpha=alpha, clip=clip, cyclic=cyclic,
            gamma1=gamma1, gamma2=gamma2, space=space
        )
        cmap._rgba_bad = self._rgba_bad
        cmap._rgba_under = self._rgba_under
        cmap._rgba_over = self._rgba_over
        return cmap

    def to_continuous(self, name=None, **kwargs):
        """
        Convert the `PerceptualColormap` to a standard `ContinuousColormap`.
        This is used to merge such colormaps.

        Parameters
        ----------
        name : str, default: '_name_copy'
            The new colormap name.

        Other parameters
        ----------------
        **kwargs
            Passed to `ContinuousColormap`.

        See also
        --------
        ContinuousColormap.to_discrete
        """
        if not self._isinit:
            self._init()
        if name is None:
            name = self._make_name()
        return ContinuousColormap.from_list(name, self._lut[:-3, :], **kwargs)

    @classmethod
    @docstring._snippet_manager
    @warnings._rename_kwargs('0.7.0', fade='saturation', shade='luminance')
    def from_color(cls, *args, **kwargs):
        """
        Return a simple monochromatic "sequential" colormap that blends from white
        or near-white to the input color.

        Parameters
        ----------
        color : color-spec
            RGB tuple, hex string, or named color string.
        %(colors.name)s
        %(colors.space)s
        l, s, a, c
            Shorthands for `luminance`, `saturation`, `alpha`, and `chroma`.
        luminance : float or color-spec, default: 100
            If float, this is the luminance channel strength on the left-hand
            side of the colormap. If RGB[A] tuple, hex string, or named color
            string, the luminance is inferred from the color.
        saturation, alpha : float or color-spec, optional
            As with `luminance`, except the default `saturation` and the default
            `alpha` are the channel values taken from `color`.
        chroma
            Alias for `saturation`.

        Other parameters
        ----------------
        **kwargs
            Passed to `PerceptualColormap.from_hsl`.

        Returns
        -------
        PerceptualColormap
            The colormap.

        See also
        --------
        PerceptualColormap.from_hsl
        PerceptualColormap.from_list
        """
        name, color, space, kwargs = cls._pop_args(
            *args, names=('color', 'space'), **kwargs
        )
        space = _not_none(space, DEFAULT_SPACE).lower()
        props = _pop_props(kwargs, 'hsla')
        if props.get('hue', None) is not None:
            raise TypeError("from_color() got an unexpected keyword argument 'hue'")
        hue, saturation, luminance, alpha = to_xyza(color, space)
        alpha_fade = props.pop('alpha', 1)
        luminance_fade = props.pop('luminance', 100)
        saturation_fade = props.pop('saturation', saturation)
        return cls.from_hsl(
            name, hue=hue, space=space,
            alpha=(alpha_fade, alpha),
            saturation=(saturation_fade, saturation),
            luminance=(luminance_fade, luminance),
            **kwargs
        )

    @classmethod
    @docstring._snippet_manager
    def from_hsl(cls, *args, **kwargs):
        """
        Make a `~PerceptualColormap` by specifying the hue,
        saturation, and luminance transitions individually.

        Parameters
        ----------
        %(colors.space)s
        %(colors.name)s
        %(colors.ratios)s
            For example, ``luminance=(100, 50, 0)`` with ``ratios=(2, 1)`` results
            in a colormap with the transition from luminance ``100`` to ``50`` taking
            *twice as long* as the transition from luminance ``50`` to ``0``.
        h, s, l, a, c
            Shorthands for `hue`, `saturation`, `luminance`, `alpha`, and `chroma`.
        hue : float or color-spec or sequence, default: 0
            Hue channel value or sequence of values. The shorthand keyword `h` is also
            acceptable. Values can be any of the following.

            1. Numbers, within the range 0 to 360 for hue and 0 to 100 for
               saturation and luminance.
            2. Color string names or hex strings, in which case the channel
               value for that color is looked up.
        saturation : float or color-spec or sequence, default: 50
            As with `hue`, but for the saturation channel.
        luminance : float or color-spec or sequence, default: ``(100, 20)``
            As with `hue`, but for the luminance channel.
        alpha : float or color-spec or sequence, default: 1
            As with `hue`, but for the alpha (opacity) channel.
        chroma
            Alias for `saturation`.

        Other parameters
        ----------------
        **kwargs
            Passed to `PerceptualColormap`.

        Returns
        -------
        PerceptualColormap
            The colormap.

        See also
        --------
        PerceptualColormap.from_color
        PerceptualColormap.from_list
        """
        name, space, ratios, kwargs = cls._pop_args(
            *args, names=('space', 'ratios'), **kwargs
        )
        cdict = {}
        props = _pop_props(kwargs, 'hsla')
        for key, default in (
            ('hue', 0),
            ('saturation', 100),
            ('luminance', (100, 20)),
            ('alpha', 1),
        ):
            value = props.pop(key, default)
            cdict[key] = _make_segment_data(value, ratios=ratios)
        return cls(name, cdict, space=space, **kwargs)

    @classmethod
    @docstring._snippet_manager
    def from_list(cls, *args, adjust_grays=True, **kwargs):
        """
        Make a `PerceptualColormap` from a sequence of colors.

        Parameters
        ----------
        %(colors.from_list)s
        adjust_grays : bool, optional
            Whether to adjust the hues of grayscale colors (including ``'white'``,
            ``'black'``, and the ``'grayN'`` open-color colors) to the hues of the
            preceding and subsequent colors in the sequence. This facilitates the
            construction of diverging colormaps with monochromatic segments using
            e.g. ``PerceptualColormap.from_list(['blue', 'white', 'red'])``.

        Other parameters
        ----------------
        **kwargs
            Passed to `PerceptualColormap`.

        Returns
        -------
        PerceptualColormap
            The colormap.

        See also
        --------
        matplotlib.colors.LinearSegmentedColormap.from_list
        ContinuousColormap.from_list
        PerceptualColormap.from_color
        PerceptualColormap.from_hsl
        """
        # Get coordinates
        coords = None
        space = kwargs.get('space', DEFAULT_SPACE).lower()
        name, colors, ratios, kwargs = cls._pop_args(
            *args, names=('colors', 'ratios'), **kwargs
        )
        if not np.iterable(colors):
            raise ValueError(f'Colors must be iterable, got colors={colors!r}')
        if (
            np.iterable(colors[0]) and len(colors[0]) == 2
            and not isinstance(colors[0], str)
        ):
            coords, colors = zip(*colors)

        # Build segmentdata
        keys = ('hue', 'saturation', 'luminance', 'alpha')
        hslas = [to_xyza(color, space) for color in colors]
        cdict = {}
        for key, values in zip(keys, zip(*hslas)):
            cdict[key] = _make_segment_data(values, coords, ratios)

        # Adjust grays
        if adjust_grays:
            hues = cdict['hue']  # segment data
            for i, color in enumerate(colors):
                rgb = to_rgb(color)
                if isinstance(color, str) and REGEX_ADJUST.match(color):
                    pass
                elif not np.allclose(np.array(rgb), rgb[0]):
                    continue
                hues[i] = list(hues[i])  # enforce mutability
                if i > 0:
                    hues[i][1] = hues[i - 1][2]
                if i < len(hues) - 1:
                    hues[i][2] = hues[i + 1][1]

        return cls(name, cdict, **kwargs)

    # Deprecated
    to_linear_segmented = warnings._rename_objs(
        '0.8.0',
        to_linear_segmented=to_continuous
    )


def _interpolate_scalar(x, x0, x1, y0, y1):
    """
    Interpolate between two points.
    """
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def _interpolate_extrapolate_vector(xq, x, y):
    """
    Interpolate between two vectors. Similar to `numpy.interp` except this
    does not truncate out-of-bounds values (i.e. this is reversible).
    """
    # Follow example of _make_lookup_table for efficient, vectorized
    # linear interpolation across multiple segments.
    # * Normal test puts values at a[i] if a[i-1] < v <= a[i]; for
    #   left-most data, satisfy a[0] <= v <= a[1]
    # * searchsorted gives where xq[i] must be inserted so it is larger
    #   than x[ind[i]-1] but smaller than x[ind[i]]
    # yq = ma.masked_array(np.interp(xq, x, y), mask=ma.getmask(xq))
    x = np.asarray(x)
    y = np.asarray(y)
    xq = np.atleast_1d(xq)
    idx = np.searchsorted(x, xq)
    idx[idx == 0] = 1  # get normed value <0
    idx[idx == len(x)] = len(x) - 1  # get normed value >0
    distance = (xq - x[idx - 1]) / (x[idx] - x[idx - 1])
    yq = distance * (y[idx] - y[idx - 1]) + y[idx - 1]
    yq = ma.masked_array(yq, mask=ma.getmask(xq))
    return yq


def _sanitize_levels(levels, minsize=2):
    """
    Ensure the levels are monotonic. If they are descending, reverse them.
    """
    # NOTE: Matplotlib does not support datetime colormap levels as of 3.5
    levels = inputs._to_numpy_array(levels)
    if levels.ndim != 1 or levels.size < minsize:
        raise ValueError(f'Levels {levels} must be a 1D array with size >= {minsize}.')
    if isinstance(levels, ma.core.MaskedArray):
        levels = levels.filled(np.nan)
    if not inputs._is_numeric(levels) or not np.all(np.isfinite(levels)):
        raise ValueError(f'Levels {levels} does not support non-numeric cmap levels.')
    diffs = np.sign(np.diff(levels))
    if np.all(diffs == 1):
        descending = False
    elif np.all(diffs == -1):
        descending = True
        levels = levels[::-1]
    else:
        raise ValueError(f'Levels {levels} must be monotonic.')
    return levels, descending


class DiscreteNorm(mcolors.BoundaryNorm):
    """
    Meta-normalizer that discretizes the possible color values returned by
    arbitrary continuous normalizers given a sequence of level boundaries.
    """
    # See this post: https://stackoverflow.com/a/48614231/4970632
    # WARNING: Must be child of BoundaryNorm. Many methods in ColorBarBase
    # test for class membership, crucially including _process_values(), which
    # if it doesn't detect BoundaryNorm will try to use DiscreteNorm.inverse().
    @warnings._rename_kwargs(
        '0.7.0', extend='unique', descending='DiscreteNorm(descending_levels)'
    )
    def __init__(
        self, levels,
        norm=None, unique=None, step=None, clip=False, ticks=None, labels=None
    ):
        """
        Parameters
        ----------
        levels : sequence of float
            The level boundaries. Must be monotonically increasing or decreasing.
            If the latter then `~DiscreteNorm.descending` is set to ``True`` and the
            colorbar axis drawn with this normalizer will be reversed.
        norm : `~matplotlib.colors.Normalize`, optional
            The normalizer used to transform `levels` and data values passed to
            `~DiscreteNorm.__call__` before discretization. The ``vmin`` and ``vmax``
            of the normalizer are set to the minimum and maximum values in `levels`.
        unique : {'neither', 'both', 'min', 'max'}, optional
            Which out-of-bounds regions should be assigned unique colormap colors.
            Possible values are equivalent to the `extend` values. Internally, proplot
            sets this depending on the user-input `extend`, whether the colormap is
            cyclic, and whether `~matplotlib.colors.Colormap.set_under`
            or `~matplotlib.colors.Colormap.set_over` were called for the colormap.
        step : float, optional
            The intensity of the transition to out-of-bounds colors as a fraction
            of the adjacent step between in-bounds colors. Internally, proplot sets
            this to ``0.5`` for cyclic colormaps and ``1`` for all other colormaps.
            This only has an effect on lower colors when `unique` is ``'min'`` or
            ``'both'``, and on upper colors when `unique` is ``'max'`` or ``'both'``.
        clip : bool, optional
            Whether to clip values falling outside of the level bins. This only
            has an effect on lower colors when `unique` is ``'min'`` or ``'both'``,
            and on upper colors when `unique` is ``'max'`` or ``'both'``.

        Other parameters
        ----------------
        ticks : array-like, default: `levels`
            Default tick values to use for colorbars drawn with this normalizer. This
            is set to the level centers when `values` is passed to a plotting command.
        labels : array-like, optional
            Default tick labels to use for colorbars drawn with this normalizer. This
            is set to values when drawing on-the-fly colorbars.

        Note
        ----
        This normalizer makes sure that levels always span the full range of
        colors in the colormap, whether `extend` is set to ``'min'``, ``'max'``,
        ``'neither'``, or ``'both'``. In matplotlib, when `extend` is not ``'both'``,
        the most intense colors are cut off (reserved for "out of bounds" data),
        even though they are not being used.

        See also
        --------
        proplot.constructor.Norm
        proplot.colors.SegmentedNorm
        proplot.ticker.DiscreteLocator
        """
        # Parse input arguments
        # NOTE: This must be a subclass BoundaryNorm, so ColorbarBase will
        # detect it... even though we completely override it.
        if step is None:
            step = 1.0
        if unique is None:
            unique = 'neither'
        if not norm:
            norm = mcolors.Normalize()
        elif isinstance(norm, mcolors.BoundaryNorm):
            raise ValueError('Normalizer cannot be instance of BoundaryNorm.')
        elif not isinstance(norm, mcolors.Normalize):
            raise ValueError('Normalizer must be instance of Normalize.')
        uniques = ('min', 'max', 'both', 'neither')
        if unique not in uniques:
            raise ValueError(
                f'Unknown unique setting {unique!r}. Options are: '
                + ', '.join(map(repr, uniques))
                + '.'
            )

        # Process level boundaries and centers
        # NOTE: Currently there are no normalizers that reverse direction
        # of levels. Tried that with SegmentedNorm but colorbar ticks fail.
        # Instead user-reversed levels will always get passed here just as
        # they are passed to SegmentedNorm inside plot.py
        levels, descending = _sanitize_levels(levels)
        vcenter = getattr(norm, 'vcenter', None)
        vmin = norm.vmin = np.min(levels)
        vmax = norm.vmax = np.max(levels)
        bins, _ = _sanitize_levels(norm(levels))
        mids = np.zeros((levels.size + 1,))
        mids[1:-1] = 0.5 * (levels[1:] + levels[:-1])
        mids[0], mids[-1] = mids[1], mids[-2]

        # Adjust color coordinate for each bin
        # For same out-of-bounds colors, looks like [0 - eps, 0, ..., 1, 1 + eps]
        # For unique out-of-bounds colors, looks like [0 - eps, X, ..., 1 - X, 1 + eps]
        # NOTE: Critical that we scale the bin centers in "physical space" and *then*
        # translate to color coordinates so that nonlinearities in the normalization
        # stay intact. If we scaled the bin centers in *normalized space* to have
        # minimum 0 maximum 1, would mess up color distribution. However this is still
        # not perfect... get asymmetric color intensity either side of central point.
        # So we add special handling for diverging norms below to improve symmetry.
        if unique in ('min', 'both'):
            scale = levels[0] - levels[1] if len(levels) == 2 else mids[1] - mids[2]
            mids[0] += step * scale
        if unique in ('max', 'both'):
            scale = levels[-1] - levels[-2] if len(levels) == 2 else mids[-2] - mids[-3]
            mids[-1] += step * scale
        mmin = np.min(mids)
        mmax = np.max(mids)
        if np.isclose(mmin, mmax):
            mmin = mmin - (mmin or 1) * 1e-10
            mmax = mmax + (mmax or 1) * 1e-10
        if vcenter is None:  # not diverging norm or centered segmented norm
            mids = _interpolate_scalar(mids, mmin, mmax, vmin, vmax)
        else:
            mask1, mask2 = mids < vcenter, mids >= vcenter
            mids[mask1] = _interpolate_scalar(mids[mask1], mmin, vcenter, vmin, vcenter)
            mids[mask2] = _interpolate_scalar(mids[mask2], vcenter, mmax, vcenter, vmax)

        # Instance attributes
        # NOTE: If clip is True, we clip values to the centers of the end bins
        # rather than vmin/vmax to prevent out-of-bounds colors from getting an
        # in-bounds bin color due to landing on a bin edge.
        # NOTE: With unique='min' the minimimum in-bounds and out-of-bounds
        # colors are the same so clip=True will have no effect. Same goes
        # for unique='max' with maximum colors.
        eps = 1e-10
        dest = norm(mids)
        dest[0] -= eps  # dest guaranteed to be numpy.float64
        dest[-1] += eps
        self._ticks = _not_none(ticks, levels)
        self._labels = labels
        self._descending = descending
        self._bmin = np.min(mids)
        self._bmax = np.max(mids)
        self._bins = bins
        self._dest = dest
        self._norm = norm
        self.N = levels.size
        self.boundaries = levels
        mcolors.Normalize.__init__(self, vmin=vmin, vmax=vmax, clip=clip)

        # Add special clipping
        # WARNING: For some reason must clip manually for LogNorm, or end
        # up with unpredictable fill value, weird "out-of-bounds" colors
        self._norm_clip = None
        if isinstance(norm, mcolors.LogNorm):
            self._norm_clip = (1e-249, None)

    def __call__(self, value, clip=None):
        """
        Normalize data values to 0-1.

        Parameters
        ----------
        value : numeric
            The data to be normalized.
        clip : bool, default: ``self.clip``
            Whether to clip values falling outside of the level bins.
        """
        # Follow example of SegmentedNorm, but perform no interpolation,
        # just use searchsorted to bin the data.
        norm_clip = self._norm_clip
        if norm_clip:  # special extra clipping due to normalizer
            value = np.clip(value, *norm_clip)
        if clip is None:  # builtin clipping
            clip = self.clip
        if clip:  # note that np.clip can handle masked arrays
            value = np.clip(value, self._bmin, self._bmax)
        xq, is_scalar = self.process_value(value)
        xq = self._norm(xq)
        yq = self._dest[np.searchsorted(self._bins, xq)]
        yq = ma.array(yq, mask=ma.getmask(xq))
        if is_scalar:
            yq = np.atleast_1d(yq)[0]
        if self.descending:
            yq = 1 - yq
        return yq

    def inverse(self, value):  # noqa: U100
        """
        Raise an error.

        Raises
        ------
        ValueError
            Inversion after discretization is impossible.
        """
        raise ValueError('DiscreteNorm is not invertible.')

    @property
    def descending(self):
        """
        Boolean indicating whether the levels are descending.
        """
        return self._descending


class SegmentedNorm(mcolors.Normalize):
    """
    Normalizer that scales data linearly with respect to the
    interpolated index in an arbitrary monotonic level sequence.
    """
    def __init__(self, levels, vcenter=None, vmin=None, vmax=None, clip=None, fair=True):  # noqa: E501
        """
        Parameters
        ----------
        levels : sequence of float
            The level boundaries. Must be monotonically increasing or decreasing.
        vcenter : float, default: None
            The central colormap value. Default is to omit this.
        vmin : float, optional
            Ignored but included for consistency. Set to ``min(levels)``.
        vmax : float, optional
            Ignored but included for consistency. Set to ``max(levels)``.
        clip : bool, optional
            Whether to clip values falling outside of `vmin` and `vmax`.
        fair : bool, optional
            Whether to use fair scaling. See `DivergingNorm`.

        See also
        --------
        proplot.constructor.Norm
        proplot.colors.DiscreteNorm

        Note
        ----
        The algorithm this normalizer uses to select normalized values
        in-between level list indices is adapted from the algorithm
        `~matplotlib.colors.LinearSegmentedColormap` uses to select channel
        values in-between segment data points (hence the name `SegmentedNorm`).

        Example
        -------
        In the below example, unevenly spaced levels are passed to
        `~matplotlib.axes.Axes.contourf`, resulting in the automatic
        application of `SegmentedNorm`.

        >>> import proplot as pplt
        >>> import numpy as np
        >>> levels = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
        >>> data = 10 ** (3 * np.random.rand(10, 10))
        >>> fig, ax = pplt.subplots()
        >>> ax.contourf(data, levels=levels)
        """
        # WARNING: Tried using descending levels by adding 1 - yq to __call__() and
        # inverse() but then tick labels fail. Instead just silently reverse here and
        # the corresponding DiscreteLocator should enforce the descending axis.
        levels, _ = _sanitize_levels(levels)
        dest = np.linspace(0, 1, len(levels))
        vmin = np.min(levels)
        vmax = np.max(levels)
        if vcenter is not None:
            center = _interpolate_extrapolate_vector(vcenter, levels, dest)
            idxs, = np.where(np.isclose(vcenter, levels))
            if fair:
                delta = center - 0.5
                delta = max(-(dest[0] - delta), dest[-1] - delta - 1)
                dest = (dest - center) / (1 + 2 * delta) + 0.5
            elif idxs.size and idxs[0] > 0 and idxs[0] < len(levels) - 1:
                dest1 = np.linspace(0, 0.5, idxs[0] + 1)
                dest2 = np.linspace(0.5, 1, len(levels) - idxs[0])
                dest = np.append(dest1, dest2[1:])
            else:
                raise ValueError(f'Center {vcenter} not in level list {levels}.')
        super().__init__(vmin=vmin, vmax=vmax, clip=clip)
        self.vcenter = vcenter  # used for DiscreteNorm
        self._x = self.boundaries = levels  # 'boundaries' are used in PlotAxes
        self._y = dest

    def __call__(self, value, clip=None):
        """
        Normalize the data values to 0-1. Inverse of `~SegmentedNorm.inverse`.

        Parameters
        ----------
        value : numeric
            The data to be normalized.
        clip : bool, default: ``self.clip``
            Whether to clip values falling outside of the minimum and maximum levels.
        """
        if clip is None:  # builtin clipping
            clip = self.clip
        if clip:  # numpy.clip can handle masked arrays
            value = np.clip(value, self.vmin, self.vmax)
        xq, is_scalar = self.process_value(value)
        yq = _interpolate_extrapolate_vector(xq, self._x, self._y)
        if is_scalar:
            yq = np.atleast_1d(yq)[0]
        return yq

    def inverse(self, value):
        """
        Inverse of `~SegmentedNorm.__call__`.

        Parameters
        ----------
        value : numeric
            The data to be un-normalized.
        """
        yq, is_scalar = self.process_value(value)
        xq = _interpolate_extrapolate_vector(yq, self._y, self._x)
        if is_scalar:
            xq = np.atleast_1d(xq)[0]
        return xq


class DivergingNorm(mcolors.Normalize):
    """
    Normalizer that ensures some central data value lies at the central
    colormap color. The default central value is ``0``.
    """
    def __str__(self):
        return type(self).__name__ + f'(center={self.vcenter!r})'

    def __init__(self, vcenter=0, vmin=None, vmax=None, clip=None, fair=True):
        """
        Parameters
        ----------
        vcenter : float, default: 0
            The central data value.
        vmin : float, optional
            The minimum data value.
        vmax : float, optional
            The maximum data value.
        clip : bool, optional
            Whether to clip values falling outside of `vmin` and `vmax`.
        fair : bool, optional
            If ``True`` (default), the speeds of the color gradations on either side
            of the center point are equal, but colormap colors may be omitted. If
            ``False``, all colormap colors are included, but the color gradations on
            one side may be faster than the other side. ``False`` should be used with
            great care, as it may result in a misleading interpretation of your data.

        See also
        --------
        proplot.constructor.Norm
        """
        # NOTE: This post is an excellent summary of matplotlib's DivergingNorm history:
        # https://github.com/matplotlib/matplotlib/issues/15336#issuecomment-535291287
        # NOTE: This is a stale PR that plans to implement the same features.
        # https://github.com/matplotlib/matplotlib/pull/15333#issuecomment-537545430
        # Since proplot is starting without matplotlib's baggage we can just implement
        # a diverging norm like they would prefer if they didn't have to worry about
        # confusing users: single class, default "fair" scaling that can be turned off.
        super().__init__(vmin, vmax, clip)
        self.vcenter = vcenter
        self.fair = fair

    def __call__(self, value, clip=None):
        """
        Normalize the data values to 0-1.

        Parameters
        ----------
        value : numeric
            The data to be normalized.
        clip : bool, default: ``self.clip``
            Whether to clip values falling outside of `vmin` and `vmax`.
        """
        xq, is_scalar = self.process_value(value)
        self.autoscale_None(xq)  # sets self.vmin, self.vmax if None
        if clip is None:  # builtin clipping
            clip = self.clip
        if clip:  # note that np.clip can handle masked arrays
            value = np.clip(value, self.vmin, self.vmax)
        if self.vmin > self.vmax:
            raise ValueError('vmin must be less than or equal to vmax.')
        elif self.vmin == self.vmax:
            x = [self.vmin, self.vmax]
            y = [0.0, 0.0]
        elif self.vcenter >= self.vmax:
            x = [self.vmin, self.vcenter]
            y = [0.0, 0.5]
        elif self.vcenter <= self.vmin:
            x = [self.vcenter, self.vmax]
            y = [0.5, 1.0]
        elif self.fair:
            offset = max(abs(self.vcenter - self.vmin), abs(self.vmax - self.vcenter))
            x = [self.vcenter - offset, self.vcenter + offset]
            y = [0.0, 1.0]
        else:
            x = [self.vmin, self.vcenter, self.vmax]
            y = [0.0, 0.5, 1.0]
        yq = _interpolate_extrapolate_vector(xq, x, y)
        if is_scalar:
            yq = np.atleast_1d(yq)[0]
        return yq

    def autoscale_None(self, z):
        """
        Get vmin and vmax, and then clip at vcenter.
        """
        super().autoscale_None(z)
        if self.vmin > self.vcenter:
            self.vmin = self.vcenter
        if self.vmax < self.vcenter:
            self.vmax = self.vcenter


def _init_color_database():
    """
    Initialize the subclassed database.
    """
    database = mcolors._colors_full_map
    if not isinstance(database, ColorDatabase):
        database = mcolors._colors_full_map = ColorDatabase(database)
        if hasattr(mcolors, 'colorConverter'):  # suspect deprecation is coming soon
            mcolors.colorConverter.cache = database.cache
            mcolors.colorConverter.colors = database
    return database


def _init_cmap_database():
    """
    Initialize the subclassed database.
    """
    # WARNING: Skip over the matplotlib native duplicate entries
    # with suffixes '_r' and '_shifted'.
    attr = '_cmap_registry' if hasattr(mcm, '_cmap_registry') else 'cmap_d'
    database = getattr(mcm, attr)
    if mcm.get_cmap is not _get_cmap:
        mcm.get_cmap = _get_cmap
    if mcm.register_cmap is not _register_cmap:
        mcm.register_cmap = _register_cmap
    if not isinstance(database, ColormapDatabase):
        database = {
            key: value for key, value in database.items()
            if key[-2:] != '_r' and key[-8:] != '_shifted'
        }
        database = ColormapDatabase(database)
        setattr(mcm, attr, database)
    return database


_mpl_register_cmap = mcm.register_cmap
@functools.wraps(_mpl_register_cmap)  # noqa: E302
def _register_cmap(*args, **kwargs):
    """
    Monkey patch for `~matplotlib.cm.register_cmap`. Ignores warning
    message when re-registering existing colormaps. This is unnecessary
    and triggers 100 warnings when importing seaborn.
    """
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        return _mpl_register_cmap(*args, **kwargs)


@functools.wraps(mcm.get_cmap)
def _get_cmap(name=None, lut=None):
    """
    Monkey patch for `~matplotlib.cm.get_cmap`. Permits case-insensitive
    search of monkey-patched colormap database. This was broken in v3.2.0
    because matplotlib now uses _check_in_list with cmap dictionary keys.
    """
    if name is None:
        name = rc['image.cmap']
    if isinstance(name, mcolors.Colormap):
        return name
    cmap = _cmap_database[name]
    if lut is not None:
        cmap = cmap._resample(lut)
    return cmap


def _get_cmap_subtype(name, subtype):
    """
    Get a colormap belonging to a particular class. If none are found then raise
    a useful error message that omits colormaps from other classes.
    """
    # NOTE: Right now this is just used in rc validation but could be used elsewhere
    if subtype == 'discrete':
        cls = DiscreteColormap
    elif subtype == 'continuous':
        cls = ContinuousColormap
    elif subtype == 'perceptual':
        cls = PerceptualColormap
    else:
        raise RuntimeError(f'Invalid subtype {subtype!r}.')
    cmap = _cmap_database.get(name, None)
    if not isinstance(cmap, cls):
        names = sorted(k for k, v in _cmap_database.items() if isinstance(v, cls))
        raise ValueError(
            f'Invalid {subtype} colormap name {name!r}. Options are: '
            + ', '.join(map(repr, names))
            + '.'
        )
    return cmap


def _translate_cmap(cmap, lut=None, cyclic=None, listedthresh=None):
    """
    Translate the input argument to a proplot colormap subclass. Auto-detect
    cyclic colormaps based on names and re-apply default lookup table size.
    """
    # Parse args
    # WARNING: Apply default 'cyclic' property to native matplotlib colormaps
    # based on known names. Maybe slightly dangerous but cleanest approach
    lut = _not_none(lut, rc['image.lut'])
    cyclic = _not_none(cyclic, cmap.name and cmap.name.lower() in CMAPS_CYCLIC)
    listedthresh = _not_none(listedthresh, rc['cmap.listedthresh'])

    # Translate the colormap
    # WARNING: Here we ignore 'N' in order to respect proplotrc lut sizes
    # when initializing proplot.
    bad = cmap._rgba_bad
    under = cmap._rgba_under
    over = cmap._rgba_over
    name = cmap.name
    if isinstance(cmap, (DiscreteColormap, ContinuousColormap)):
        pass
    elif isinstance(cmap, mcolors.LinearSegmentedColormap):
        data = dict(cmap._segmentdata)
        cmap = ContinuousColormap(name, data, N=lut, gamma=cmap._gamma, cyclic=cyclic)
    elif isinstance(cmap, mcolors.ListedColormap):
        colors = list(cmap.colors)
        if len(colors) > listedthresh:  # see notes at top of file
            cmap = ContinuousColormap.from_list(name, colors, N=lut, cyclic=cyclic)
        else:
            cmap = DiscreteColormap(colors, name)
    elif isinstance(cmap, mcolors.Colormap):  # base class
        pass
    else:
        raise ValueError(
            f'Invalid colormap type {type(cmap).__name__!r}. '
            'Must be instance of matplotlib.colors.Colormap.'
        )

    # Apply hidden settings
    cmap._rgba_bad = bad
    cmap._rgba_under = under
    cmap._rgba_over = over

    return cmap


class _ColorCache(dict):
    """
    Replacement for the native color cache.
    """
    def __getitem__(self, key):
        """
        Get the standard color, colormap color, or color cycle color.
        """
        # NOTE: Matplotlib 'color' args are passed to to_rgba, which tries to read
        # directly from cache and if that fails, sanitizes input, which raises
        # error on receiving (colormap, idx) tuple. So we have to override cache.
        return self._get_rgba(*key)

    def _get_rgba(self, arg, alpha):
        """
        Try to get the color from the registered colormap or color cycle.
        """
        key = (arg, alpha)
        if isinstance(arg, str) or not np.iterable(arg) or len(arg) != 2:
            return dict.__getitem__(self, key)
        if not isinstance(arg[0], str) or not isinstance(arg[1], Number):
            return dict.__getitem__(self, key)
        # Try to get the colormap
        try:
            cmap = _cmap_database[arg[0]]
        except (KeyError, TypeError):
            return dict.__getitem__(self, key)
        # Read the colormap value
        if isinstance(cmap, DiscreteColormap):
            if not 0 <= arg[1] < len(cmap.colors):
                raise ValueError(
                    f'Color cycle sample for {arg[0]!r} cycle must be '
                    f'between 0 and {len(cmap.colors) - 1}, got {arg[1]}.'
                )
            rgba = cmap.colors[arg[1]]  # draw from list of colors
        else:
            if not 0 <= arg[1] <= 1:
                raise ValueError(
                    f'Colormap sample for {arg[0]!r} colormap must be '
                    f'between 0 and 1, got {arg[1]}.'
                )
            rgba = cmap(arg[1])  # get color selection
        # Return the colormap value
        rgba = to_rgba(rgba)
        a = _not_none(alpha, rgba[3])
        return (*rgba[:3], a)


class ColorDatabase(MutableMapping, dict):
    """
    Dictionary subclass used to replace the builtin matplotlib color database.
    See `~ColorDatabase.__getitem__` for details.
    """
    _colors_replace = (
        ('grey', 'gray'),  # British --> American synonyms
        ('ochre', 'ocher'),  # ...
        ('kelley', 'kelly'),  # backwards compatibility to correct spelling
    )

    def __iter__(self):
        yield from dict.__iter__(self)

    def __len__(self):
        return dict.__len__(self)

    def __delitem__(self, key):
        key = self._parse_key(key)
        dict.__delitem__(self, key)
        self.cache.clear()

    def __init__(self, mapping=None):
        """
        Parameters
        ----------
        mapping : dict-like, optional
            The colors.
        """
        # NOTE: Tested with and without standardization and speedup is marginal
        self._cache = _ColorCache()
        mapping = mapping or {}
        for key, value in mapping.items():
            self.__setitem__(key, value)

    def __getitem__(self, key):
        """
        Get a color. Translates ``grey`` into ``gray`` and supports retrieving
        colors "on-the-fly" from registered colormaps and color cycles.

        * For a colormap, use e.g. ``color=('Blues', 0.8)``.
          The number is the colormap index, and must be between 0 and 1.
        * For a color cycle, use e.g. ``color=('colorblind', 2)``.
          The number is the color list index.

        This works everywhere that colors are used in matplotlib, for
        example as `color`, `edgecolor', or `facecolor` keyword arguments
        passed to `~proplot.axes.PlotAxes` commands.
        """
        key = self._parse_key(key)
        return dict.__getitem__(self, key)

    def __setitem__(self, key, value):
        """
        Add a color. Translates ``grey`` into ``gray`` and clears the
        cache. The color must be a string.
        """
        # Always standardize assignments.
        key = self._parse_key(key)
        dict.__setitem__(self, key, value)
        self.cache.clear()

    def _parse_key(self, key):
        """
        Parse the color key. Currently this just translates grays.
        """
        if not isinstance(key, str):
            raise ValueError(f'Invalid color name {key!r}. Must be string.')
        if isinstance(key, str) and len(key) > 1:  # ignore base colors
            key = key.lower()
            for sub, rep in self._colors_replace:
                key = key.replace(sub, rep)
        return key

    @property
    def cache(self):
        # Matplotlib uses 'cache' but treat '_cache' as synonym
        # to guard against private API changes.
        return self._cache


class ColormapDatabase(MutableMapping, dict):
    """
    Dictionary subclass used to replace the matplotlib
    colormap registry. See `~ColormapDatabase.__getitem__` and
    `~ColormapDatabase.__setitem__` for details.
    """
    _regex_grays = re.compile(r'\A(grays)(_r|_s)*\Z', flags=re.IGNORECASE)
    _regex_suffix = re.compile(r'(_r|_s)*\Z', flags=re.IGNORECASE)

    def __iter__(self):
        yield from dict.__iter__(self)

    def __len__(self):
        return dict.__len__(self)

    def __delitem__(self, key):
        key = self._parse_key(key, mirror=True)
        dict.__delitem__(self, key)

    def __init__(self, kwargs):
        """
        Parameters
        ----------
        kwargs : dict-like
            The source dictionary.
        """
        for key, value in kwargs.items():
            self.__setitem__(key, value)

    def __getitem__(self, key):
        """
        Retrieve the colormap associated with the sanitized key name. The
        key name is case insensitive.

        * If the key ends in ``'_r'``, the result of ``cmap.reversed()`` is
          returned for the colormap registered under the preceding name.
        * If the key ends in ``'_s'``, the result of ``cmap.shifted(180)`` is
          returned for the colormap registered under the preceding name.
        * Reversed diverging colormaps can be requested with their "reversed"
          name -- for example, ``'BuRd'`` is equivalent to ``'RdBu_r'``.
        """
        return self._get_item(key)

    def __setitem__(self, key, value):
        """
        Store the colormap under its lowercase name. If the object is a
        `matplotlib.colors.ListedColormap` and ``cmap.N`` is smaller than
        :rc:`cmap.listedthresh`, it is converted to a `proplot.colors.DiscreteColormap`.
        Otherwise, it is converted to a `proplot.colors.ContinuousColormap`.
        """
        self._set_item(key, value)

    def _translate_deprecated(self, key):
        """
        Check if a colormap has been deprecated.
        """
        # WARNING: Must search only for case-sensitive *capitalized* names or we would
        # helpfully "redirect" user to SciVisColor cmap when they are trying to
        # generate open-color monochromatic cmaps and would disallow some color names
        if isinstance(key, str):
            test = self._regex_suffix.sub('', key)
        else:
            test = None
        if not self._has_item(test) and test in CMAPS_REMOVED:
            version = CMAPS_REMOVED[test]
            raise ValueError(
                f'The colormap name {key!r} was removed in version {version}.'
            )
        if not self._has_item(test) and test in CMAPS_RENAMED:
            test_new, version = CMAPS_RENAMED[test]
            warnings._warn_proplot(
                f'The colormap name {test!r} was deprecated in version {version} '
                f'and may be removed in {warnings._next_release()}. Please use '
                f'the colormap name {test_new!r} instead.'
            )
            key = re.sub(test, test_new, key, flags=re.IGNORECASE)
        return key

    def _translate_key(self, key, mirror=True):
        """
        Return the sanitized colormap name. Used for lookups and assignments.
        """
        # Sanitize key
        if not isinstance(key, str):
            raise KeyError(f'Invalid key {key!r}. Key must be a string.')
        key = key.lower()
        key = self._regex_grays.sub(r'greys\2', key)
        # Mirror diverging
        reverse = key[-2:] == '_r'
        if reverse:
            key = key[:-2]
        if mirror and not self._has_item(key):  # avoid recursion here
            key_mirror = CMAPS_DIVERGING.get(key, None)
            if key_mirror and self._has_item(key_mirror):
                reverse = not reverse
                key = key_mirror
        if reverse:
            key = key + '_r'
        return key

    def _has_item(self, key):
        """
        Redirect to unsanitized `dict.__contains__`.
        """
        return dict.__contains__(self, key)

    def _get_item(self, key):
        """
        Get the colormap with flexible input keys.
        """
        # Sanitize key
        key = self._translate_deprecated(key)
        key = self._translate_key(key, mirror=True)
        shift = key[-2:] == '_s' and not self._has_item(key)
        if shift:
            key = key[:-2]
        reverse = key[-2:] == '_r' and not self._has_item(key)
        if reverse:
            key = key[:-2]
        # Retrieve colormap
        try:
            value = dict.__getitem__(self, key)  # may raise keyerror
        except KeyError:
            raise KeyError(
                f'Invalid colormap or color cycle name {key!r}. Options are: '
                + ', '.join(map(repr, self))
                + '.'
            )
        # Modify colormap
        if reverse:
            value = value.reversed()
        if shift:
            value = value.shifted(180)
        return value

    def _set_item(self, key, value):
        """
        Add the colormap after validating and converting.
        """
        if not isinstance(key, str):
            raise KeyError(f'Invalid key {key!r}. Must be string.')
        if not isinstance(value, mcolors.Colormap):
            raise ValueError('Object is not a colormap.')
        key = self._translate_key(key, mirror=False)
        value = _translate_cmap(value)
        dict.__setitem__(self, key, value)


# Initialize databases
_cmap_database = _init_cmap_database()
_color_database = _init_color_database()

# Deprecated
(
    ListedColormap,
    LinearSegmentedColormap,
    PerceptuallyUniformColormap,
    LinearSegmentedNorm,
) = warnings._rename_objs(  # noqa: E501
    '0.8.0',
    ListedColormap=DiscreteColormap,
    LinearSegmentedColormap=ContinuousColormap,
    PerceptuallyUniformColormap=PerceptualColormap,
    LinearSegmentedNorm=SegmentedNorm,
)
