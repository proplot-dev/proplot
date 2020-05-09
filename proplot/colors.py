#!/usr/bin/env python3
"""
New colormap classes and colormap normalization classes.
"""
import os
import re
import json
import glob
import cycler
from xml.etree import ElementTree
from .utils import to_rgb, to_xyz
from numbers import Number, Integral
from matplotlib import rcParams
import numpy as np
from .internals import warnings
import numpy.ma as ma
import matplotlib.colors as mcolors
import matplotlib.cm as mcm
from .utils import _warn_proplot, _notNone, _timer
try:  # use this for debugging instead of print()!
    from icecream import ic
except ImportError:  # graceful fallback if IceCream isn't installed
    ic = lambda *a: None if not a else (a[0] if len(a) == 1 else a)  # noqa

__all__ = [
    'BinNorm', 'CmapDict', 'ColorDict',
    'LinearSegmentedNorm',
    'LinearSegmentedColormap',
    'ListedColormap',
    'MidpointNorm', 'PerceptuallyUniformColormap',
    'make_mapping_array',
    'Colormap', 'Colors', 'Cycle', 'Norm',
]

# Colormap stuff
CYCLES_TABLE = {
    'Matplotlib defaults': (
        'default', 'classic',
    ),
    'Matplotlib stylesheets': (
        'colorblind', 'colorblind10', 'ggplot', 'bmh', 'solarized', '538',
    ),
    'ColorBrewer2.0 qualitative': (
        'Accent', 'Dark2',
        'Paired', 'Pastel1', 'Pastel2',
        'Set1', 'Set2', 'Set3',
    ),
    'Other qualitative': (
        'FlatUI', 'Qual1', 'Qual2',
    ),
}
CMAPS_TABLE = {
    # Assorted origin, but these belong together
    'Grayscale': (
        'Grays', 'Mono', 'MonoCycle',
    ),
    # Builtin
    'Matplotlib sequential': (
        'viridis', 'plasma', 'inferno', 'magma', 'cividis',
    ),
    'Matplotlib cyclic': (
        'twilight',
    ),
    # Seaborn
    'Seaborn sequential': (
        'Rocket', 'Mako',
    ),
    'Seaborn diverging': (
        'IceFire', 'Vlag',
    ),
    # PerceptuallyUniformColormap
    'ProPlot sequential': (
        'Fire',
        'Stellar',
        'Boreal',
        'Marine',
        'Dusk',
        'Glacial',
        'Sunrise',
        'Sunset',
    ),
    'ProPlot diverging': (
        'Div', 'NegPos', 'DryWet',
    ),
    # Nice diverging maps
    'Other diverging': (
        'ColdHot', 'CoolWarm', 'BR',
    ),
    # cmOcean
    'cmOcean sequential': (
        'Oxy', 'Thermal', 'Dense', 'Ice', 'Haline',
        'Deep', 'Algae', 'Tempo', 'Speed', 'Turbid', 'Solar', 'Matter',
        'Amp',
    ),
    'cmOcean diverging': (
        'Balance', 'Delta', 'Curl',
    ),
    'cmOcean cyclic': (
        'Phase',
    ),
    # Fabio Crameri
    'Scientific colour maps sequential': (
        'batlow', 'oleron',
        'devon', 'davos', 'oslo', 'lapaz', 'acton',
        'lajolla', 'bilbao', 'tokyo', 'turku', 'bamako', 'nuuk',
        'hawaii', 'buda', 'imola',
    ),
    'Scientific colour maps diverging': (
        'roma', 'broc', 'cork', 'vik', 'berlin', 'lisbon', 'tofino',
    ),
    'Scientific colour maps cyclic': (
        'romaO', 'brocO', 'corkO', 'vikO',
    ),
    # ColorBrewer
    'ColorBrewer2.0 sequential': (
        'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
        'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
        'PuBu', 'PuBuGn', 'BuGn', 'GnBu', 'YlGnBu', 'YlGn'
    ),
    'ColorBrewer2.0 diverging': (
        'Spectral', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGY',
        'RdBu', 'RdYlBu', 'RdYlGn',
    ),
    # SciVisColor
    'SciVisColor blues': (
        'Blue1', 'Blue2', 'Blue3', 'Blue4', 'Blue5',
        'Blue6', 'Blue7', 'Blue8', 'Blue9', 'Blue10', 'Blue11',
    ),
    'SciVisColor greens': (
        'Green1', 'Green2', 'Green3', 'Green4', 'Green5',
        'Green6', 'Green7', 'Green8',
    ),
    'SciVisColor oranges': (
        'Orange1', 'Orange2', 'Orange3', 'Orange4', 'Orange5',
        'Orange6', 'Orange7', 'Orange8',
    ),
    'SciVisColor browns': (
        'Brown1', 'Brown2', 'Brown3', 'Brown4', 'Brown5',
        'Brown6', 'Brown7', 'Brown8', 'Brown9',
    ),
    'SciVisColor reds and purples': (
        'RedPurple1', 'RedPurple2', 'RedPurple3', 'RedPurple4',
        'RedPurple5', 'RedPurple6', 'RedPurple7', 'RedPurple8',
    ),
    # Builtin maps that will be deleted; categories are taken from comments in
    # matplotlib source code. Some of these are really bad, some are segmented
    # maps when the should be color cycles, and some are just uninspiring.
    'MATLAB': (
        'bone', 'cool', 'copper', 'autumn', 'flag', 'prism',
        'jet', 'hsv', 'hot', 'spring', 'summer', 'winter', 'pink', 'gray',
    ),
    'GNUplot': (
        'gnuplot', 'gnuplot2', 'ocean', 'afmhot', 'rainbow',
    ),
    'GIST': (
        'gist_earth', 'gist_gray', 'gist_heat', 'gist_ncar',
        'gist_rainbow', 'gist_stern', 'gist_yarg',
    ),
    'Other': (
        'binary', 'bwr', 'brg',  # appear to be custom matplotlib
        'cubehelix', 'Wistia', 'CMRmap',  # individually released
        'seismic', 'terrain', 'nipy_spectral',  # origin ambiguous
        'tab10', 'tab20', 'tab20b', 'tab20c',  # merged colormap cycles
    )
}
CMAPS_DIVERGING = tuple(
    (key1.lower(), key2.lower()) for key1, key2 in (
        ('PiYG', 'GYPi'),
        ('PRGn', 'GnRP'),
        ('BrBG', 'GBBr'),
        ('PuOr', 'OrPu'),
        ('RdGy', 'GyRd'),
        ('RdBu', 'BuRd'),
        ('RdYlBu', 'BuYlRd'),
        ('RdYlGn', 'GnYlRd'),
        ('BR', 'RB'),
        ('CoolWarm', 'WarmCool'),
        ('ColdHot', 'HotCold'),
        ('NegPos', 'PosNeg'),
        ('DryWet', 'WetDry')
    ))

# Named color filter props
COLORS_HEXPATTERN = r'#(?:[0-9a-fA-F]{3,4}){2}'  # 6-8 digit hex
COLORS_HCLTHRESH = 0.10  # bigger number equals fewer colors
COLORS_SPACE = 'hcl'  # color "distincness" is defined with this space
COLORS_ADD = (
    'charcoal', 'sky blue', 'eggshell', 'sea blue', 'coral', 'aqua',
    'tomato red', 'brick red', 'crimson',
    'red orange', 'yellow orange', 'yellow green', 'blue green',
    'blue violet', 'red violet',
)  # common names that should always be included
COLORS_BASE = {
    'blue': (0, 0, 1),
    'green': (0, 0.5, 0),
    'red': (1, 0, 0),
    'cyan': (0, 0.75, 0.75),
    'magenta': (0.75, 0, 0.75),
    'yellow': (0.75, 0.75, 0),
    'black': (0, 0, 0),
    'white': (1, 1, 1),
    **mcolors.BASE_COLORS,  # shorthand names like 'r', 'g', etc.
}
COLORS_RMREGEX = re.compile(
    '(' + '|'.join((
        'shit', 'poop', 'poo', 'pee', 'piss', 'puke', 'vomit', 'snot',
        'booger', 'bile', 'diarrhea',
    )) + ')'
)  # filter these out, let's try to be professional here...
COLORS_SUBREGEXES = tuple(
    (re.compile(regex), sub)
    for regex, sub in (
        ('/', ' '),
        ('\'s', ''),
        ('forrest', 'forest'),  # typo?
        (r'\s?majesty', ''),  # purple mountains majesty is too long
        ('reddish', 'red'),  # remove 'ish'
        ('purplish', 'purple'),
        ('bluish', 'blue'),
        (r'ish\b', ''),
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
)  # prevent registering similar-sounding names


def _get_channel(color, channel, space='hcl'):
    """
    Get the hue, saturation, or luminance channel value from the input color.
    The color name `color` can optionally be a string with the format
    ``'color+x'`` or ``'color-x'``, where `x` specifies the offset from the
    channel value.

    Parameters
    ----------
    color : color-spec
        The color. Sanitized with `to_rgb`.
    channel : {'hue', 'chroma', 'saturation', 'luminance'}
        The HCL channel to be retrieved.
    space : {'hcl', 'hpl', 'hsl', 'hsv', 'rgb'}, optional
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
        match = re.search('([-+][0-9.]+)$', color)
        if match:
            offset = float(match.group(0))
            color = color[:match.start()]
    return offset + to_xyz(color, space)[channel]


def _clip_colors(colors, clip=True, gray=0.2, warn=False):
    """
    Clip impossible colors rendered in an HSL-to-RGB colorspace conversion.
    Used by `PerceptuallyUniformColormap`. If `mask` is ``True``, impossible
    colors are masked out.

    Parameters
    ----------
    colors : list of length-3 tuples
        The RGB colors.
    clip : bool, optional
        If `clip` is ``True`` (the default), RGB channel values >1 are clipped
        to 1. Otherwise, the color is masked out as gray.
    gray : float, optional
        The identical RGB channel values (gray color) to be used if `mask`
        is ``True``.
    warn : bool, optional
        Whether to issue warning when colors are clipped.
    """
    colors = np.array(colors)
    over = colors > 1
    under = colors < 0
    if clip:
        colors[under] = 0
        colors[over] = 1
    else:
        colors[under | over] = gray
    if warn:
        msg = 'Clipped' if clip else 'Invalid'
        for i, name in enumerate('rgb'):
            if under[:, i].any():
                warnings._warn_proplot(f'{msg} {name!r} channel ( < 0).')
            if over[:, i].any():
                warnings._warn_proplot(f'{msg} {name!r} channel ( > 1).')
    return colors


def _make_segmentdata_array(values, coords=None, ratios=None):
    """
    Return a segmentdata array or callable given the input colors
    and coordinates.

    Parameters
    ----------
    values : list of float
        The channel values.
    coords : list of float, optional
        The segment coordinates.
    ratios : list of float, optional
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
            _warn_proplot(
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
                f'Need {len(values)-1} ratios for {len(values)} colors, '
                f'but got {len(ratios)} ratios.'
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


def make_mapping_array(N, data, gamma=1.0, inverse=False):
    r"""
    Similar to `~matplotlib.colors.makeMappingArray` but permits
    *circular* hue gradations along 0-360, disables clipping of
    out-of-bounds channel values, and uses fancier "gamma" scaling.

    Parameters
    ----------
    N : int
        Number of points in the colormap lookup table.
    data : 2D array-like
        List of :math:`(x, y_0, y_1)` tuples specifying the channel jump (from
        :math:`y_0` to :math:`y_1`) and the :math:`x` coordinate of that
        transition (ranges between 0 and 1).
        See `~matplotlib.colors.LinearSegmentedColormap` for details.
    gamma : float or list of float, optional
        To obtain channel values between coordinates :math:`x_i` and
        :math:`x_{i+1}` in rows :math:`i` and :math:`i+1` of `data`,
        we use the formula:

        .. math::

            y = y_{1,i} + w_i^{\gamma_i}*(y_{0,i+1} - y_{1,i})

        where :math:`\gamma_i` corresponds to `gamma` and the weight
        :math:`w_i` ranges from 0 to 1 between rows ``i`` and ``i+1``.
        If `gamma` is float, it applies to every transition. Otherwise,
        its length must equal ``data.shape[0]-1``.

        This is like the `gamma` used with matplotlib's
        `~matplotlib.colors.makeMappingArray`, except it controls the
        weighting for transitions *between* each segment data coordinate rather
        than the coordinates themselves. This makes more sense for
        `PerceptuallyUniformColormap`\\ s because they usually consist of just
        one linear transition for *sequential* colormaps and two linear
        transitions for *diverging* colormaps -- and in the latter case, it
        is often desirable to modify both "halves" of the colormap in the
        same way.
    inverse : bool, optional
        If ``True``, :math:`w_i^{\gamma_i}` is replaced with
        :math:`1 - (1 - w_i)^{\gamma_i}` -- that is, when `gamma` is greater
        than 1, this weights colors toward *higher* channel values instead
        of lower channel values.

        This is implemented in case we want to apply *equal* "gamma scaling"
        to different HSL channels in different directions. Usually, this
        is done to weight low data values with higher luminance *and* lower
        saturation, thereby emphasizing "extreme" data values with stronger
        colors.
    """
    # Allow for *callable* instead of linearly interpolating between segments
    gammas = np.atleast_1d(gamma)
    if (gammas < 0.01).any() or (gammas > 10).any():
        raise ValueError('Gamma can only be in range [0.01,10].')
    if callable(data):
        if len(gammas) > 1:
            raise ValueError(
                'Only one gamma allowed for functional segmentdata.')
        x = np.linspace(0, 1, N)**gamma
        lut = np.array(data(x), dtype=float)
        return lut

    # Get array
    data = np.array(data)
    shape = data.shape
    if len(shape) != 2 or shape[1] != 3:
        raise ValueError('Data must be nx3 format.')
    if len(gammas) != 1 and len(gammas) != shape[0] - 1:
        raise ValueError(
            f'Need {shape[0]-1} gammas for {shape[0]}-level mapping array, '
            f'but got {len(gamma)}.'
        )
    if len(gammas) == 1:
        gammas = np.repeat(gammas, shape[:1])

    # Get indices
    x = data[:, 0]
    y0 = data[:, 1]
    y1 = data[:, 2]
    if x[0] != 0.0 or x[-1] != 1.0:
        raise ValueError(
            'Data mapping points must start with x=0 and end with x=1.'
        )
    if (np.diff(x) < 0).any():
        raise ValueError(
            'Data mapping points must have x in increasing order.'
        )
    x = x * (N - 1)

    # Get distances from the segmentdata entry to the *left* for each requested
    # level, excluding ends at (0,1), which must exactly match segmentdata ends
    xq = (N - 1) * np.linspace(0, 1, N)
    # where xq[i] must be inserted so it is larger than x[ind[i]-1] but
    # smaller than x[ind[i]]
    ind = np.searchsorted(x, xq)[1:-1]
    distance = (xq[1:-1] - x[ind - 1]) / (x[ind] - x[ind - 1])

    # Scale distances in each segment by input gamma
    # The ui are starting-points, the ci are counts from that point
    # over which segment applies (i.e. where to apply the gamma), the relevant
    # 'segment' is to the *left* of index returned by searchsorted
    _, uind, cind = np.unique(ind, return_index=True, return_counts=True)
    for ui, ci in zip(uind, cind):  # length should be N-1
        # the relevant segment is to *left* of this number
        gamma = gammas[ind[ui] - 1]
        if gamma == 1:
            continue
        ireverse = False
        if ci > 1:  # i.e. more than 1 color in this 'segment'
            # by default want to weight toward a *lower* channel value
            ireverse = ((y0[ind[ui]] - y1[ind[ui] - 1]) < 0)
        if inverse:
            ireverse = not ireverse
        if ireverse:
            distance[ui:ui + ci] = 1 - (1 - distance[ui:ui + ci])**gamma
        else:
            distance[ui:ui + ci] **= gamma

    # Perform successive linear interpolations all rolled up into one equation
    lut = np.zeros((N,), float)
    lut[1:-1] = distance * (y0[ind] - y1[ind - 1]) + y1[ind - 1]
    lut[0] = y1[0]
    lut[-1] = y0[-1]
    return lut


_from_file_docstring = """
Valid file extensions are as follows:

==================  =====================================================================================================================================================================================================================
Extension           Description
==================  =====================================================================================================================================================================================================================
``.hex``            List of HEX strings in any format (comma-separated, separate lines, with double quotes... anything goes).
``.xml``            XML files with ``<Point .../>`` tags specifying ``x``, ``r``, ``g``, ``b``, and (optionally) ``o`` parameters, where ``x`` is the coordinate and the rest are the red, blue, green, and opacity channel values.
``.rgb``, ``.txt``  3-4 column table of red, blue, green, and (optionally) opacity channel values, delimited by commas or spaces. If values larger than 1 are detected, they are assumed to be on the 0-255 scale and are divided by 255.
==================  =====================================================================================================================================================================================================================

Parameters
----------
path : str
    The file path.
warn_on_failure : bool, optional
    If ``True``, issue a warning when loading fails rather than
    raising an error.
"""  # noqa


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
            data = ', '.join(mcolors.to_hex(color) for color in colors)
        elif ext in ('txt', 'rgb'):
            rgb = mcolors.to_rgba if alpha else mcolors.to_rgb
            data = [rgb(color) for color in colors]
            data = '\n'.join(
                ' '.join(f'{num:0.6f}' for num in line) for line in data
            )
        else:
            raise ValueError(
                f'Invalid extension {ext!r}. Options are: '
                "'hex', 'txt', 'rgb', 'rgba'."
            )
        return data

    def _parse_path(self, path, dirname='.', ext=''):
        """
        Parse the user input path.

        Parameters
        ----------
        dirname : str, optional
            The default directory.
        ext : str, optional
            The default extension.
        """
        path = os.path.expanduser(path or '')
        dirname = os.path.expanduser(dirname or '')
        if not path or os.path.isdir(path):
            path = os.path.join(path or dirname, self.name)  # default name
        dirname, basename = os.path.split(path)  # default to current directory
        path = os.path.join(dirname or '.', basename)
        if not os.path.splitext(path)[-1]:
            path = path + '.' + ext  # default file extension
        return path


class LinearSegmentedColormap(mcolors.LinearSegmentedColormap, _Colormap):
    r"""
    New base class for all `~matplotlib.colors.LinearSegmentedColormap`\ s.
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
                string += (f' {key!r}: [{data[0][2]:.3f}, '
                           f'..., {data[-1][1]:.3f}],\n')
        return type(self).__name__ + '({\n' + string + '})'

    def __init__(self, *args, cyclic=False, alpha=None, **kwargs):
        """
        Parameters
        ----------
        cyclic : bool, optional
            Whether the colormap is cyclic. If ``True``, this changes how the
            leftmost and rightmost color levels are selected, and `extend` can
            only be ``'neither'`` (a warning will be issued otherwise).
        alpha : float, optional
            The opacity for the entire colormap. Overrides the input
            segment data.

        Other parameters
        ----------------
        *args, **kwargs
            Passed to `~matplotlib.colors.LinearSegmentedColormap`.
        """
        super().__init__(*args, **kwargs)
        self._cyclic = cyclic
        if alpha is not None:
            self.set_alpha(alpha)

    def concatenate(self, *args, ratios=1, name=None, N=None, **kwargs):
        """
        Return the concatenation of this colormap with the
        input colormaps.

        Parameters
        ----------
        *args
            Instances of `LinearSegmentedColormap`.
        ratios : list of float, optional
            Relative extent of each component colormap in the merged colormap.
            Length must equal ``len(args) + 1``.

            For example, ``cmap1.concatenate(cmap2, ratios=[2,1])`` generates
            a colormap with the left two-thrids containing colors from
            ``cmap1`` and the right one-third containing colors from ``cmap2``.
        name : str, optional
            The colormap name. Default is
            ``'_'.join(cmap.name for cmap in args)``.
        N : int, optional
            The number of points in the colormap lookup table.
            Default is :rc:`image.lut` times ``len(args)``.

        Other parameters
        ----------------
        **kwargs
            Passed to `LinearSegmentedColormap.updated`
            or `PerceptuallyUniformColormap.updated`.

        Returns
        -------
        `LinearSegmentedColormap`
            The colormap.
        """
        # Parse input args
        if not args:
            raise ValueError(
                f'Got zero positional args, you must provide at least one.'
            )
        if not all(isinstance(cmap, type(self)) for cmap in args):
            raise ValueError(
                f'Colormaps {cmap.name + ": " + repr(cmap) for cmap in args} '
                f'must all belong to the same class.'
            )
        cmaps = (self, *args)
        spaces = {cmap.name: getattr(cmap, '_space', None) for cmap in cmaps}
        if len({*spaces.values(), }) > 1:
            raise ValueError(
                'Cannot merge PerceptuallyUniformColormaps that use '
                'different colorspaces: '
                + ', '.join(map(repr, spaces)) + '.'
            )
        N = N or len(cmaps) * rcParams['image.lut']
        if name is None:
            name = '_'.join(cmap.name for cmap in cmaps)

        # Combine the segmentdata, and use the y1/y2 slots at merge points so
        # we never interpolate between end colors of different colormaps
        segmentdata = {}
        ratios = ratios or 1
        if isinstance(ratios, Number):
            ratios = [1] * len(cmaps)
        ratios = np.array(ratios) / np.sum(ratios)
        x0 = np.concatenate([[0], np.cumsum(ratios)])  # coordinates for edges
        xw = x0[1:] - x0[:-1]  # widths between edges
        for key in self._segmentdata.keys():
            # Callable segments
            # WARNING: If just reference a global 'funcs' list from inside the
            # 'data' function it can get overwritten in this loop. Must
            # embed 'funcs' into the definition using a keyword argument.
            callable_ = [callable(cmap._segmentdata[key]) for cmap in cmaps]
            if all(callable_):  # expand range from x-to-w to 0-1
                funcs = [cmap._segmentdata[key] for cmap in cmaps]

                def xyy(ix, funcs=funcs):
                    ix = np.atleast_1d(ix)
                    kx = np.empty(ix.shape)
                    for j, jx in enumerate(ix.flat):
                        idx = max(np.searchsorted(x0, jx) - 1, 0)
                        kx.flat[j] = funcs[idx]((jx - x0[idx]) / xw[idx])
                    return kx
            # Concatenate segment arrays and make the transition at the
            # seam instant so we *never interpolate* between end colors
            # of different maps.
            elif not any(callable_):
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
                raise ValueError(
                    'Mixed callable and non-callable colormap values.'
                )
            segmentdata[key] = xyy
            # Handle gamma values
            if key == 'saturation':
                ikey = 'gamma1'
            elif key == 'luminance':
                ikey = 'gamma2'
            else:
                continue
            if ikey in kwargs:
                continue
            gamma = []
            for cmap in cmaps:
                igamma = getattr(cmap, '_' + ikey)
                if not np.iterable(igamma):
                    if all(callable_):
                        igamma = [igamma]
                    else:
                        igamma = (len(cmap._segmentdata[key]) - 1) * [igamma]
                gamma.extend(igamma)
            if all(callable_):
                if any(igamma != gamma[0] for igamma in gamma[1:]):
                    _warn_proplot(
                        'Cannot use multiple segment gammas when '
                        'concatenating callable segments. Using the first '
                        f'gamma of {gamma[0]}.'
                    )
                gamma = gamma[0]
            kwargs[ikey] = gamma

        # Return copy
        return self.updated(name=name, segmentdata=segmentdata, **kwargs)

    def punched(self, cut=None, name=None, **kwargs):
        """
        Return a version of the colormap with the center "punched out".
        This is great for making the transition from "negative" to "positive"
        in a diverging colormap more distinct.

        Parameters
        ----------
        cut : float, optional
            The proportion to cut from the center of the colormap.
            For example, ``center=0.1`` cuts the central 10%.
        name : str, optional
            The name of the new colormap. Default is
            ``self.name + '_punched'``.
        **kwargs
            Passed to `LinearSegmentedColormap.updated`
            or `PerceptuallyUniformColormap.updated`.

        Returns
        -------
        `LinearSegmentedColormap`
            The colormap.
        """
        cut = _notNone(cut, 0)
        if cut == 0:
            return self
        if name is None:
            name = self.name + '_punched'

        # Decompose cut into two truncations followed by concatenation
        left_center = 0.5 - cut / 2
        right_center = 0.5 + cut / 2
        cmap_left = self.truncated(0, left_center)
        cmap_right = self.truncated(right_center, 1)
        return cmap_left.concatenate(cmap_right, name=name, **kwargs)

    def reversed(self, name=None, **kwargs):
        """
        Return a reversed copy of the colormap, as in
        `~matplotlib.colors.LinearSegmentedColormap`.

        Parameters
        ----------
        name : str, optional
            The new colormap name. Default is ``self.name + '_r'``.
        **kwargs
            Passed to `LinearSegmentedColormap.updated`
            or `PerceptuallyUniformColormap.updated`.
        """
        if name is None:
            name = self.name + '_r'

        def factory(dat):
            def func_r(x):
                return dat(1.0 - x)
            return func_r
        segmentdata = {key:
                       factory(data) if callable(data) else
                       [(1.0 - x, y1, y0) for x, y0, y1 in reversed(data)]
                       for key, data in self._segmentdata.items()}
        for key in ('gamma1', 'gamma2'):
            if key in kwargs:
                continue
            gamma = getattr(self, '_' + key, None)
            if gamma is not None and np.iterable(gamma):
                kwargs[key] = gamma[::-1]
        return self.updated(name, segmentdata, **kwargs)

    def save(self, path=None, alpha=True):
        """
        Save the colormap data to a file.

        Parameters
        ----------
        path : str, optional
            The output filename. If not provided, the colormap
            is saved under ``~/.proplot/cmaps/name.json`` where ``name``
            is the colormap name. Valid extensions are described in
            the below table.

            =====================  ==========================================================
            Extension              Description
            =====================  ==========================================================
            ``.json`` (default)    JSON database of the channel segment data.
            ``.hex``               Comma-delimited list of HEX strings.
            ``.rgb``, ``.txt``     3-4 column table of channel values.
            =====================  ==========================================================

        alpha : bool, optional
            Whether to include an opacity column for ``.rgb``
            and ``.txt`` files.
        """  # noqa
        dirname = os.path.join('~', '.proplot', 'cmaps')
        filename = self._parse_path(path, dirname, 'json')

        # Save channel segment data in json file
        _, ext = os.path.splitext(filename)
        if ext[1:] == 'json':
            # Sanitize segmentdata values
            # Convert np.float to builtin float, np.array to list of lists,
            # and callable to list of lists. We tried encoding func.__code__
            # with base64 and marshal instead, but when cmap.concatenate()
            # embeds functions as keyword arguments, this seems to make it
            # *impossible* to load back up the function with FunctionType
            # (error message: arg 5 (closure) must be tuple). Instead use
            # this brute force workaround.
            data = {}
            for key, value in self._segmentdata.items():
                if callable(value):
                    x = np.linspace(0, 1, 256)  # just save the transitions
                    y = np.array([value(_) for _ in x]).squeeze()
                    value = np.vstack((x, y, y)).T
                data[key] = np.asarray(value).astype(float).tolist()
            # Add critical attributes to the dictionary
            keys = ()
            if isinstance(self, PerceptuallyUniformColormap):
                keys = ('cyclic', 'gamma1', 'gamma2', 'space')
            elif isinstance(self, LinearSegmentedColormap):
                keys = ('cyclic', 'gamma')
            for key in keys:
                data[key] = getattr(self, '_' + key)
            with open(filename, 'w') as file:
                json.dump(data, file, indent=4)

        # Save lookup table colors
        else:
            data = self._get_data(ext[1:], alpha=alpha)
            with open(filename, 'w') as f:
                f.write(data)
        print(f'Saved colormap to {filename!r}.')

    def set_alpha(self, alpha):
        """
        Set the opacity for the entire colormap.

        Parameters
        ----------
        alpha : float
            The opacity.
        """
        self._segmentdata['alpha'] = [(0, alpha, alpha), (1, alpha, alpha)]
        self._isinit = False

    def set_cyclic(self, b):
        """
        Set whether this colormap is "cyclic". See `LinearSegmentedColormap`
        for details.
        """
        self._cyclic = bool(b)
        self._isinit = False

    def shifted(self, shift=None, name=None, **kwargs):
        """
        Return a cyclicaly shifted version of the colormap. If the colormap
        cyclic property is set to ``False`` a warning will be raised.

        Parameters
        ----------
        shift : float, optional
            The number of degrees to shift, out of 360 degrees. If ``None``,
            the original colormap is returned.
        name : str, optional
            The name of the new colormap. Default is
            ``self.name + '_s'``.

        Other parameters
        ----------------
        **kwargs
            Passed to `LinearSegmentedColormap.copy`
            or `PerceptuallyUniformColormap.copy`.
        """
        shift = ((shift or 0) / 360) % 1
        if shift == 0:
            return self
        if name is None:
            name = self.name + '_s'
        if not self._cyclic:
            _warn_proplot(
                f'Shifting non-cyclic colormap {self.name!r}. '
                f'Use cmap.set_cyclic(True) or Colormap(..., cyclic=True) to '
                'suppress this warning.'
            )
            self._cyclic = True

        # Decompose shift into two truncations followed by concatenation
        cmap_left = self.truncated(shift, 1)
        cmap_right = self.truncated(0, shift)
        return cmap_left.concatenate(
            cmap_right, ratios=(1 - shift, shift), name=name, **kwargs
        )

    def truncated(self, left=None, right=None, name=None, **kwargs):
        """
        Return a truncated version of the colormap.

        Parameters
        ----------
        left : float, optional
            The colormap index for the new "leftmost" color. Must fall between
            ``0`` and ``1``. For example,
            ``left=0.1`` cuts the leftmost 10%% of the colors.
        right : float, optional
            The colormap index for the new "rightmost" color. Must fall between
            ``0`` and ``1``. For example,
            ``right=0.9`` cuts the leftmost 10%% of the colors.
        name : str, optional
            The name of the new colormap. Default is
            ``self.name + '_truncated'``.

        Other parameters
        ----------------
        **kwargs
            Passed to `LinearSegmentedColormap.updated`
            or `PerceptuallyUniformColormap.updated`.
        """
        # Bail out
        left = max(_notNone(left, 0), 0)
        right = min(_notNone(right, 1), 1)
        if left == 0 and right == 1:
            return self
        if name is None:
            name = self.name + '_truncated'

        # Resample the segmentdata arrays
        segmentdata = {}
        for key, xyy in self._segmentdata.items():
            # Callable array
            # WARNING: If just reference a global 'xyy' callable from inside
            # the lambda function it gets overwritten in the loop! Must embed
            # the old callable in the new one as a default keyword arg.
            if callable(xyy):
                def xyy(x, func=xyy):
                    return func(left + x * (right - left))
            # Slice
            # l is the first point where x > 0 or x > left, should be >= 1
            # r is the last point where r < 1 or r < right
            else:
                xyy = np.asarray(xyy)
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
            segmentdata[key] = xyy
            # Retain the corresponding gamma *segments*
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
                        _warn_proplot(
                            'Cannot use multiple segment gammas when '
                            'truncating colormap. Using the first gamma '
                            f'of {gamma[0]}.'
                        )
                    gamma = gamma[0]
                else:
                    gamma = gamma[l - 1:r + 1]
            kwargs[ikey] = gamma
        return self.updated(name, segmentdata, **kwargs)

    def updated(
        self, name=None, segmentdata=None, N=None, *,
        alpha=None, gamma=None, cyclic=None
    ):
        """
        Return a new colormap, with relevant properties copied from this one
        if they were not provided as keyword arguments.

        Parameters
        ----------
        name : str
            The colormap name. Default is ``self.name + '_updated'``.
        segmentdata, N, alpha, gamma, cyclic : optional
            See `LinearSegmentedColormap`. If not provided,
            these are copied from the current colormap.
        """
        if name is None:
            name = self.name + '_updated'
        if segmentdata is None:
            segmentdata = self._segmentdata
        if gamma is None:
            gamma = self._gamma
        if cyclic is None:
            cyclic = self._cyclic
        if N is None:
            N = self.N
        cmap = LinearSegmentedColormap(
            name, segmentdata, N,
            alpha=alpha, gamma=gamma, cyclic=cyclic
        )
        cmap._rgba_bad = self._rgba_bad
        cmap._rgba_under = self._rgba_under
        cmap._rgba_over = self._rgba_over
        return cmap

    @staticmethod
    def from_file(path, warn_on_failure=False):
        """Load colormap from a file."""
        return _from_file(path, listed=False, warn_on_failure=warn_on_failure)

    @staticmethod
    def from_list(name, colors, ratios=None, **kwargs):
        """
        Make a `LinearSegmentedColormap` from a list of colors.

        Parameters
        ----------
        name : str
            The colormap name.
        colors : list of color-spec or (float, color-spec) tuples, optional
            If list of RGB[A] tuples or color strings, the colormap transitions
            evenly from ``colors[0]`` at the left-hand side to
            ``colors[-1]`` at the right-hand side.

            If list of (float, color-spec) tuples, the float values are the
            coordinate of each transition and must range from 0 to 1. This
            can be used to divide  the colormap range unevenly.
        ratios : list of float, optional
            Relative extents of each color transition. Must have length
            ``len(colors) - 1``. Larger numbers indicate a slower
            transition, smaller numbers indicate a faster transition.

        Other parameters
        ----------------
        **kwargs
            Passed to `LinearSegmentedColormap`.

        Returns
        -------
        `LinearSegmentedColormap`
            The colormap.
        """
        # Get coordinates
        coords = None
        if not np.iterable(colors):
            raise TypeError('Colors must be iterable.')
        if (
            np.iterable(colors[0])
            and len(colors[0]) == 2
            and not isinstance(colors[0], str)
        ):
            coords, colors = zip(*colors)
        colors = [to_rgb(color, alpha=True) for color in colors]

        # Build segmentdata
        keys = ('red', 'green', 'blue', 'alpha')
        cdict = {}
        for key, values in zip(keys, zip(*colors)):
            cdict[key] = _make_segmentdata_array(values, coords, ratios)
        return LinearSegmentedColormap(name, cdict, **kwargs)

    # Fix docstrings
    # NOTE: Docstrings cannot be programatically altered in place e.g.
    # with f-strings. Can only be modified a posteriori.
    from_file.__func__.__doc__ += _from_file_docstring


class ListedColormap(mcolors.ListedColormap, _Colormap):
    r"""
    New base class for all `~matplotlib.colors.ListedColormap`\ s.
    """
    def __str__(self):
        return f'ListedColormap(name={self.name!r})'

    def __repr__(self):
        return (
            'ListedColormap({\n'
            f" 'name': {self.name!r},\n"
            f" 'colors': {[mcolors.to_hex(color) for color in self.colors]},\n"
            '})')

    def __init__(self, *args, alpha=None, **kwargs):
        """
        Parameters
        ----------
        alpha : float, optional
            The opacity for the entire colormap. Overrides the input
            colors.

        Other parameters
        ----------------
        *args, **kwargs
            Passed to `~matplotlib.colors.ListedColormap`.
        """
        super().__init__(*args, **kwargs)
        if alpha is not None:
            self.set_alpha(alpha)

    def concatenate(self, *args, name=None, N=None, **kwargs):
        """
        Append arbitrary colormaps onto this colormap.

        Parameters
        ----------
        *args
            Instances of `ListedColormap`.
        name : str, optional
            The colormap name. Default is
            ``'_'.join(cmap.name for cmap in args)``.
        N : int, optional
            The number of colors in the colormap lookup table. Default is
            the number of colors in the concatenated lists.

        Other parameters
        ----------------
        **kwargs
            Passed to `~ListedColormap.updated`.
        """
        if not args:
            raise ValueError(
                f'Got zero positional args, you must provide at least one.'
            )
        if not all(isinstance(cmap, type(self)) for cmap in args):
            raise ValueError(
                f'Input arguments {args} must all be ListedColormap.'
            )
        cmaps = (self, *args)
        if name is None:
            name = '_'.join(cmap.name for cmap in cmaps)
        colors = [color for cmap in cmaps for color in cmap.colors]
        return self.updated(colors, name, N or len(colors), **kwargs)

    def save(self, path=None, alpha=True):
        """
        Save the colormap data to a file.

        Parameters
        ----------
        path : str, optional
            The output filename. If not provided, the colormap
            is saved under ``~/.proplot/cycles/name.hex`` where ``name``
            is the colormap name. Valid extensions are described in
            the below table.

            =====================  ==========================================================
            Extension              Description
            =====================  ==========================================================
            ``.hex`` (default)     Comma-delimited list of HEX strings.
            ``.rgb``, ``.txt``     3-4 column table of channel values.
            =====================  ==========================================================

        alpha : bool, optional
            Whether to include an opacity column for ``.rgb``
            and ``.txt`` files.
        """  # noqa
        dirname = os.path.join('~', '.proplot', 'cycles')
        filename = self._parse_path(path, dirname, 'hex')

        # Save lookup table colors
        _, ext = os.path.splitext(filename)
        data = self._get_data(ext[1:], alpha=alpha)
        with open(filename, 'w') as f:
            f.write(data)
        print(f'Saved colormap to {filename!r}.')

    def set_alpha(self, alpha):
        """
        Set the opacity for the entire colormap.

        Parameters
        ----------
        alpha : float
            The opacity.
        """
        colors = [list(mcolors.to_rgba(color)) for color in self.colors]
        for color in colors:
            color[3] = alpha
        self.colors = colors
        self._init()

    def shifted(self, shift=None, name=None):
        """
        Return a cyclically shifted version of the colormap.

        Parameters
        ----------
        shift : float, optional
            The number of places to shift, between ``-self.N`` and ``self.N``.
            If ``None``, the original colormap is returned.
        name : str, optional
            The new colormap name. Default is ``self.name + '_s'``.
        """
        if not shift:
            return self
        if name is None:
            name = self.name + '_s'
        shift = shift % len(self.colors)
        colors = [*self.colors]  # ensure list
        colors = colors[shift:] + colors[:shift]
        return self.updated(colors, name, len(colors))

    def truncated(self, left=None, right=None, name=None):
        """
        Return a truncated version of the colormap.

        Parameters
        ----------
        left : float, optional
            The colormap index for the new "leftmost" color. Must fall between
            ``0`` and ``self.N``. For example,
            ``left=2`` deletes the two first colors.
        right : float, optional
            The colormap index for the new "rightmost" color. Must fall between
            ``0`` and ``self.N``. For example,
            ``right=4`` deletes colors after the fourth color.
        name : str, optional
            The new colormap name. Default is ``self.name + '_truncated'``.
        """
        if left is None and right is None:
            return self
        if name is None:
            name = self.name + '_truncated'
        colors = self.colors[left:right]
        return self.updated(colors, name, len(colors))

    def updated(self, colors=None, name=None, N=None, *, alpha=None):
        """
        Return a new colormap with relevant properties copied from this one
        if they were not provided as keyword arguments.

        Parameters
        ----------
        name : str
            The colormap name. Default is ``self.name + '_updated'``.
        colors, N, alpha : optional
            See `ListedColormap`. If not provided,
            these are copied from the current colormap.
        """
        if name is None:
            name = self.name + '_updated'
        if colors is None:
            colors = self.colors
        if N is None:
            N = self.N
        cmap = ListedColormap(colors, name, N, alpha=alpha)
        cmap._rgba_bad = self._rgba_bad
        cmap._rgba_under = self._rgba_under
        cmap._rgba_over = self._rgba_over
        return cmap

    @staticmethod
    def from_file(path, warn_on_failure=False):
        """
        Load color cycle from a file.
        """
        return _from_file(path, listed=True, warn_on_failure=warn_on_failure)

    # Fix docstrings
    # NOTE: Docstrings cannot be programatically altered in place e.g.
    # with f-strings. Can only be modified a posteriori.
    from_file.__func__.__doc__ += _from_file_docstring


class PerceptuallyUniformColormap(LinearSegmentedColormap, _Colormap):
    """
    Similar to `~matplotlib.colors.LinearSegmentedColormap`, but instead
    of varying the RGB channels, we vary hue, saturation, and luminance in
    either the HCL colorspace or the HSL or HPL scalings of HCL.
    """
    def __init__(
        self, name, segmentdata, N=None, space=None, clip=True,
        gamma=None, gamma1=None, gamma2=None,
        **kwargs
    ):
        """
        Parameters
        ----------
        name : str
            The colormap name.
        segmentdata : dict-like
            Mapping containing the keys ``'hue'``, ``'saturation'``, and
            ``'luminance'``. The key values can be callable functions that
            return channel values given a colormap index, or lists containing
            any of the following channel specifiers:

            1. Numbers, within the range 0-360 for hue and 0-100 for
               saturation and luminance.
            2. Color string names or hex tags, in which case the channel
               value for that color is looked up.

            See `~matplotlib.colors.LinearSegmentedColormap` for a more
            detailed explanation.
        N : int, optional
            Number of points in the colormap lookup table.
            Default is :rc:`image.lut`.
        space : {'hsl', 'hpl', 'hcl'}, optional
            The hue, saturation, luminance-style colorspace to use for
            interpreting the channels. See
            `this page <http://www.hsluv.org/comparison/>`__ for a description.
        clip : bool, optional
            Whether to "clip" impossible colors, i.e. truncate HCL colors
            with RGB channels with values >1, or mask them out as gray.
        cyclic : bool, optional
            Whether the colormap is cyclic. If ``True``, this changes how the
            leftmost and rightmost color levels are selected, and `extend` can
            only be ``'neither'`` (a warning will be issued otherwise).
        gamma : float, optional
            Sets `gamma1` and `gamma2` to this identical value.
        gamma1 : float, optional
            If >1, makes low saturation colors more prominent. If <1,
            makes high saturation colors more prominent. Similar to the
            `HCLWizard <http://hclwizard.org:64230/hclwizard/>`_ option.
            See `make_mapping_array` for details.
        gamma2 : float, optional
            If >1, makes high luminance colors more prominent. If <1,
            makes low luminance colors more prominent. Similar to the
            `HCLWizard <http://hclwizard.org:64230/hclwizard/>`_ option.
            See `make_mapping_array` for details.

        Other parameters
        ----------------
        **kwargs
            Passed to `LinearSegmentedColormap`.

        Example
        -------
        The following generates a `PerceptuallyUniformColormap` from a
        `segmentdata` dictionary that uses color names for the hue data,
        instead of channel values between ``0`` and ``360``.

        >>> import proplot as plot
        >>> data = {
        ...     'hue': [[0, 'red', 'red'], [1, 'blue', 'blue']],
        ...     'saturation': [[0, 100, 100], [1, 100, 100]],
        ...     'luminance': [[0, 100, 100], [1, 20, 20]],
        ... }
        >>> cmap = plot.PerceptuallyUniformColormap(data)

        """
        # Checks
        space = _notNone(space, 'hsl').lower()
        if space not in ('rgb', 'hsv', 'hpl', 'hsl', 'hcl'):
            raise ValueError(f'Unknown colorspace {space!r}.')
        keys = {*segmentdata.keys()}
        target = {'hue', 'saturation', 'luminance', 'alpha'}
        if not keys <= target:
            raise ValueError(
                f'Invalid segmentdata dictionary with keys {keys!r}.'
            )
        # Convert color strings to channel values
        for key, array in segmentdata.items():
            if callable(array):  # permit callable
                continue
            for i, xyy in enumerate(array):
                xyy = list(xyy)  # make copy!
                for j, y in enumerate(xyy[1:]):  # modify the y values
                    xyy[j + 1] = _get_channel(y, key, space)
                segmentdata[key][i] = xyy
        # Initialize
        N = N or rcParams['image.lut']
        super().__init__(name, segmentdata, N, gamma=1.0, **kwargs)
        # Custom properties
        self._gamma1 = _notNone(gamma1, gamma, 1.0)
        self._gamma2 = _notNone(gamma2, gamma, 1.0)
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
        for i, (channel, gamma, inverse) in enumerate(
                zip(channels, gammas, inverses)):
            self._lut_hsl[:-3, i] = make_mapping_array(
                self.N, self._segmentdata[channel], gamma, inverse)
        if 'alpha' in self._segmentdata:
            self._lut_hsl[:-3, 3] = make_mapping_array(
                self.N, self._segmentdata['alpha'])
        self._lut_hsl[:-3, 0] %= 360

        # Make hues circular, set extremes i.e. copy HSL values
        self._lut = self._lut_hsl.copy()
        self._set_extremes()  # generally just used end values in segmentdata
        self._isinit = True

        # Now convert values to RGB and clip colors
        for i in range(self.N + 3):
            self._lut[i, :3] = to_rgb(self._lut[i, :3], self._space)
        self._lut[:, :3] = _clip_colors(self._lut[:, :3], self._clip)

    @docstring.add_snippets
    def set_gamma(self, gamma=None, gamma1=None, gamma2=None):
        """
        Modify the gamma value(s) and refresh the lookup table.

        Parameters
        ----------
        %(cmap.gamma)s
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
        name : str
            The colormap name. Default is ``self.name + '_copy'``.
        segmentdata, N, alpha, clip, cyclic, gamma, gamma1, gamma2, space : \
optional
            See `PerceptuallyUniformColormap`. If not provided,
            these are copied from the current colormap.
        """
        if name is None:
            name = self.name + '_copy'
        if segmentdata is None:
            segmentdata = self._segmentdata
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
        cmap = PerceptuallyUniformColormap(
            name, segmentdata, N,
            alpha=alpha, clip=clip, cyclic=cyclic,
            gamma1=gamma1, gamma2=gamma2, space=space
        )
        cmap._rgba_bad = self._rgba_bad
        cmap._rgba_under = self._rgba_under
        cmap._rgba_over = self._rgba_over
        return cmap

    def to_linear_segmented(self, **kwargs):
        """
        Convert the `PerceptuallyUniformColormap` to a standard
        `LinearSegmentedColormap`. This is used to merge such colormaps.

        Parameters
        ----------
        **kwargs
            Passed to `LinearSegmentedColormap`.
        """
        if not self._isinit:
            self._init()
        return LinearSegmentedColormap.from_list(
            self.name, self._lut, **kwargs
        )

    @classmethod
    def from_color(cls, name, color, fade=None, space='hsl', **kwargs):
        """
        Return a monochromatic "sequential" colormap that blends from white
        or near-white to the input color.

        Parameters
        ----------
        name : str, optional
            The colormap name.
        color : color-spec
            RGB tuple, hex string, or named color string.
        fade : float or color-spec, optional
            If float, this is the luminance channel strength on the left-hand
            side of the colormap (default is ``100``), and the saturation
            channel is held constant throughout the colormap.

            If RGB tuple, hex string, or named color string, the luminance and
            saturation (but *not* the hue) from this color are used for the
            left-hand side of the colormap.
        space : {'hsl', 'hpl', 'hcl'}, optional
            The colorspace in which the luminance is varied.

        Other parameters
        ----------------
        **kwargs
            Passed to `PerceptuallyUniformColormap.from_hsl`.

        Returns
        -------
        `PerceptuallyUniformColormap`
            The colormap.
        """
        hue, saturation, luminance, alpha = to_xyz(color, space, alpha=True)
        if fade is None:
            fade = 100
        if isinstance(fade, Number):
            saturation_fade, luminance_fade = saturation, fade
        else:
            _, saturation_fade, luminance_fade = to_xyz(fade, space)
        return PerceptuallyUniformColormap.from_hsl(
            name, hue=hue, alpha=alpha, space=space,
            saturation=(saturation_fade, saturation),
            luminance=(luminance_fade, luminance),
            **kwargs
        )

    @staticmethod
    def from_hsl(
        name, hue=0, saturation=100, luminance=(100, 20), alpha=None,
        ratios=None, **kwargs
    ):
        """
        Make a `~PerceptuallyUniformColormap` by specifying the hue,
        saturation, and luminance transitions individually.

        Parameters
        ----------
        name : str, optional
            The colormap name.
        hue : float, str, or list thereof, optional
            Hue channel value or list of values. Values can be
            any of the following.

            1. Numbers, within the range 0-360 for hue and 0-100 for
               saturation and luminance.
            2. Color string names or hex strings, in which case the channel
               value for that color is looked up.

            If scalar, the hue does not change across the colormap.
        saturation, luminance, alpha : float, str, or list thereof, optional
            As with `hue`, but for the saturation, luminance, and alpha
            (opacity) channels, respectively.
        ratios : list of float, optional
            Relative extents of each color transition. Must have length
            ``len(colors) - 1``. Larger numbers indicate a slower
            transition, smaller numbers indicate a faster transition.

            For example, ``luminance=[100,50,0]`` with ``ratios=[2,1]``
            results in a colormap with the transition from luminance ``100``
            to ``50`` taking *twice as long* as the transition from luminance
            ``50`` to ``0``.

        Other parameters
        ----------------
        **kwargs
            Passed to `PerceptuallyUniformColormap`.

        Returns
        -------
        `PerceptuallyUniformColormap`
            The colormap.
        """
        cdict = {}
        alpha = _notNone(alpha, 1.0)
        for key, channel in zip(
            ('hue', 'saturation', 'luminance', 'alpha'),
            (hue, saturation, luminance, alpha)
        ):
            cdict[key] = _make_segmentdata_array(channel, ratios=ratios)
        return PerceptuallyUniformColormap(name, cdict, **kwargs)

    @staticmethod
    def from_list(name, colors, ratios=None, **kwargs):
        """
        Make a `PerceptuallyUniformColormap` from a list of colors.

        Parameters
        ----------
        name : str
            The colormap name.
        colors : list of color-spec or (float, color-spec) tuples, optional
            If list of RGB[A] tuples or color strings, the colormap transitions
            evenly from ``colors[0]`` at the left-hand side to
            ``colors[-1]`` at the right-hand side.

            If list of (float, color-spec) tuples, the float values are the
            coordinate of each transition and must range from 0 to 1. This
            can be used to divide  the colormap range unevenly.
        ratios : list of float, optional
            Relative extents of each color transition. Must have length
            ``len(colors) - 1``. Larger numbers indicate a slower
            transition, smaller numbers indicate a faster transition.

            For example, ``red=[1,0.5,0]`` with ``ratios=[2,1]``
            results in a colormap with the transition from red ``1``
            to ``0.5`` taking *twice as long* as the transition from red
            ``0.5`` to ``0``.

        Other parameters
        ----------------
        **kwargs
            Passed to `PerceptuallyUniformColormap`.

        Returns
        -------
        `PerceptuallyUniformColormap`
            The colormap.
        """
        # Get coordinates
        coords = None
        space = kwargs.get('space', 'hsl')  # use the builtin default
        if not np.iterable(colors):
            raise ValueError(f'Colors must be iterable, got colors={colors!r}')
        if (np.iterable(colors[0]) and len(colors[0]) == 2
                and not isinstance(colors[0], str)):
            coords, colors = zip(*colors)
        colors = [to_xyz(color, space, alpha=True) for color in colors]

        # Build segmentdata
        keys = ('hue', 'saturation', 'luminance', 'alpha')
        cdict = {}
        for key, values in zip(keys, zip(*colors)):
            cdict[key] = _make_segmentdata_array(values, coords, ratios)
        return cls(name, cdict, **kwargs)


def Colors(*args, **kwargs):
    """
    Pass all arguments to `Cycle` and return the list of colors from
    the resulting `~cycler.Cycler` object.
    """
    cycle = Cycle(*args, **kwargs)
    return [dict_['color'] for dict_ in cycle]


def Colormap(
    *args, name=None, listmode='perceptual', fade=None, cycle=None,
    shift=None, cut=None, left=None, right=None, reverse=False,
    save=False, save_kw=None, **kwargs
):
    """
    Generate or retrieve colormaps and optionally merge and manipulate
    them in a variety of ways. Used to interpret the `cmap` and `cmap_kw`
    arguments when passed to any plotting method wrapped by
    `~proplot.wrappers.cmap_changer`.

    Parameters
    ----------
    *args : colormap-spec
        Positional arguments that individually generate colormaps. If more than
        one argument is passed, the resulting colormaps are merged. Arguments
        are interpreted as follows.

        * If `~matplotlib.colors.Colormap` or a registered colormap name, the
          colormap is simply returned.
        * If a filename string with valid extension, the colormap data will
          be loaded. See `register_cmaps` and `register_cycles`.
        * If RGB tuple or color string, a `PerceptuallyUniformColormap` is
          generated with `~PerceptuallyUniformColormap.from_color`. If the
          string ends in ``'_r'``, the monochromatic map will be *reversed*,
          i.e. will go from dark to light instead of light to dark.
        * If list of RGB tuples or color strings, a
          `PerceptuallyUniformColormap` is generated with
          `~PerceptuallyUniformColormap.from_list`.
        * If dictionary containing the keys ``'hue'``, ``'saturation'``, and
          ``'luminance'``, a `PerceptuallyUniformColormap` is generated with
          `~PerceptuallyUniformColormap.from_hsl`.

    name : str, optional
        Name under which the final colormap is registered. It can then be
        reused by passing ``cmap='name'`` to plotting functions like
        `~matplotlib.axes.Axes.contourf`.
    fade : float, optional
        The maximum luminosity used when generating colormaps with
        `PerceptuallyUniformColormap.from_color`. Default is ``100`` when
        calling `Colormap` directly, and ``90`` when `Colormap` is called by
        `Cycle` (this prevents having pure white in the color cycle).

        For example, ``plot.Colormap('blue', fade=80)`` generates a blue
        colormap that fades to a pale blue with 80% luminance.
    cycle : str or list of color-spec, optional
        The registered cycle name or a list of colors used to interpret cycle
        color strings like ``'C0'`` and ``'C2'`` when generating colormaps
        with `PerceptuallyUniformColormap.from_color`. Default is colors
        from the currently active property cycler.
        For example, ``plot.Colormap('C0', 'C1', 'C2', cycle='538')``
        generates a colormap using colors from the ``'538'`` color cycle.
    listmode : {'perceptual', 'linear', 'listed'}, optional
        Controls how colormaps are generated when you input list(s) of colors.
        If ``'perceptual'``, a `PerceptuallyUniformColormap` is generated with
        `PerceptuallyUniformColormap.from_list`. If ``'linear'``,
        a `~matplotlib.colors.LinearSegmentedColormap` is generated with
        `~matplotlib.colors.LinearSegmentedColormap.from_list`. If
        ``'listed'``, the `~matplotlib.colors.ListedColormap` is generated.

        Default is ``'perceptual'`` when calling `Colormap` directly, and
        ``'listed'`` when `Colormap` is called by `Cycle`.
    cut : float, optional
        Passed to `LinearSegmentedColormap.punched`.
        This applies to the final *merged* colormap.
    left, right : float, optional
        Passed to `LinearSegmentedColormap.truncated` or
        `ListedColormap.truncated`. These apply to *each colormap*
        individually.
    reverse : bool, optional
        Passed to `LinearSegmentedColormap.reversed` or
        `ListedColormap.reversed`. This applies to *each colormap*
        individually.
    shift : float, optional
        Passed to `LinearSegmentedColormap.shifted` or
        `ListedColormap.shifted`. This applies to the final *merged* colormap.
    save : bool, optional
        Whether to call the colormap save method, i.e.
        `LinearSegmentedColormap.save` or
        `ListedColormap.save`.
    save_kw : dict-like, optional
        Ignored if `save` is ``False``. Passed to the colormap save method,
        i.e. `LinearSegmentedColormap.save` or
        `ListedColormap.save`.

    Other parameters
    ----------------
    **kwargs
        Passed to `LinearSegmentedColormap.concatenate` or
        `ListedColormap.concatenate`. Each of these functions accepts
        arbitrary colormap settings.

    Returns
    -------
    `~matplotlib.colors.Colormap`
        A `~matplotlib.colors.LinearSegmentedColormap` or
        `~matplotlib.colors.ListedColormap` instance.
    """
    # Parse input args
    # TODO: Play with using "qualitative" colormaps in realistic examples,
    # how to make colormaps cyclic.
    if not args:
        raise ValueError(
            f'Colormap() requires at least one positional argument.'
        )
    if listmode not in ('listed', 'linear', 'perceptual'):
        raise ValueError(
            f'Invalid listmode={listmode!r}. Options are: '
            "'listed', 'linear', 'perceptual'."
        )
    tmp = '_no_name'
    cmaps = []
    for i, cmap in enumerate(args):
        # Load registered colormaps and maps on file
        # TODO: Document how 'listmode' also affects loaded files
        if isinstance(cmap, str):
            if '.' in cmap:
                if listmode == 'listed':
                    cmap = ListedColormap.from_file(cmap)
                else:
                    cmap = LinearSegmentedColormap.from_file(cmap)
            else:
                try:
                    cmap = mcm.cmap_d[cmap]
                except KeyError:
                    pass
        # Convert matplotlib colormaps to subclasses
        if isinstance(cmap, (ListedColormap, LinearSegmentedColormap)):
            pass
        elif isinstance(cmap, mcolors.LinearSegmentedColormap):
            cmap = LinearSegmentedColormap(
                cmap.name, cmap._segmentdata, cmap.N, cmap._gamma)
        elif isinstance(cmap, mcolors.ListedColormap):
            cmap = ListedColormap(
                cmap.colors, cmap.name, cmap.N)
        # Dictionary of hue/sat/luminance values or 2-tuples representing
        # linear transition
        elif isinstance(cmap, dict):
            cmap = PerceptuallyUniformColormap.from_hsl(tmp, **cmap)
        # List of color tuples or color strings, i.e. iterable of iterables
        elif (not isinstance(cmap, str) and np.iterable(cmap)
              and all(np.iterable(color) for color in cmap)):
            try:
                colors = [to_rgb(color, cycle=cycle, alpha=True)
                          for color in cmap]
            except (ValueError, TypeError):
                raise ValueError(f'Invalid color(s) in list {cmap!r}.')
            if listmode == 'listed':
                cmap = ListedColormap(colors, tmp)
            elif listmode == 'linear':
                cmap = LinearSegmentedColormap.from_list(tmp, colors)
            else:
                cmap = PerceptuallyUniformColormap.from_list(tmp, colors)
        # Monochrome colormap from input color
        else:
            ireverse = (isinstance(cmap, str) and cmap[-2:] == '_r')
            if ireverse:
                cmap = cmap[:-2]
            try:
                color = to_rgb(cmap, cycle=cycle, alpha=True)
            except (ValueError, TypeError):
                msg = f'Invalid cmap, cycle, or color {cmap!r}.'
                if isinstance(cmap, str):
                    msg += (
                        f'\nValid cmap and cycle names: '
                        + ', '.join(map(repr, sorted(mcm.cmap_d))) + '.'
                        f'\nValid color names: '
                        + ', '.join(map(repr, sorted(
                            mcolors.colorConverter.colors))
                        ) + '.'
                    )
                raise ValueError(msg)
            cmap = PerceptuallyUniformColormap.from_color(tmp, color, fade)
            if ireverse:
                cmap = cmap.reversed()

        # Cut the edges and/or reverse the map
        if left is not None or right is not None:
            cmap = cmap.truncated(left, right)
        if reverse:
            cmap = cmap.reversed()
        cmaps.append(cmap)

    # Merge the result of this arbitrary user input
    if len(cmaps) > 1:  # more than one map?
        cmap = cmaps[0].concatenate(*cmaps[1:], **kwargs)
    elif kwargs:  # modify any props?
        cmap = cmaps[0].updated(**kwargs)

    # Cut the center and roate the colormap
    if cut is not None:
        cmap = cmap.punched(cut)
    if shift is not None:
        cmap = cmap.shifted(shift)

    # Initialize
    if not cmap._isinit:
        cmap._init()

    # Register and save the colormap
    if name is None:
        name = cmap.name  # may have been modified by e.g. .shifted()
    else:
        cmap.name = name
    mcm.cmap_d[name] = cmap
    if save:
        save_kw = save_kw or {}
        cmap.save(**save_kw)
    return cmap


def Cycle(
    *args, N=None, name=None,
    marker=None, alpha=None, dashes=None, linestyle=None, linewidth=None,
    markersize=None, markeredgewidth=None,
    markeredgecolor=None, markerfacecolor=None,
    save=False, save_kw=None, **kwargs
):
    """
    Generate and merge `~cycler.Cycler` instances in a variety of ways.
    Used to interpret the `cycle` and `cycle_kw` arguments when passed to
    any plotting method wrapped by
    `~proplot.wrappers.cycle_changer`.

    If you just want a list of colors instead of a `~cycler.Cycler` instance,
    use the `colors` function. If you want a `~cycler.Cycler` instance that
    imposes black as the default color and cycles through properties like
    ``linestyle`` instead, call this function without any positional arguments.

    Parameters
    ----------
    *args : colormap-spec or cycle-spec, optional
        Positional arguments control the *colors* in the `~cycler.Cycler`
        object. If more than one argument is passed, the resulting cycles are
        merged. Arguments are interpreted as follows.

        * If `~cycler.Cycler`, nothing more
          is done.
        * If list of RGB tuples or color strings, these
          colors are used.
        * If `~matplotlib.colors.ListedColormap`, colors from the ``colors``
          attribute are used.
        * If string color cycle name, that `~matplotlib.colors.ListedColormap`
          is looked up and its ``colors`` attribute is used. See `cycles`.
        * Otherwise, the argument is passed to `Colormap`, and colors
          from the resulting `~matplotlib.colors.LinearSegmentedColormap`
          are used. See the `N` argument.

        If the last positional argument is numeric, it is used for the `N`
        keyword argument.
    N : float or list of float, optional
        For `~matplotlib.colors.ListedColormap`\ s, this is the number of
        colors to select. For example, ``Cycle('538', 4)`` returns the first 4
        colors of the ``'538'`` color cycle.

        For `~matplotlib.colors.LinearSegmentedColormap`\ s, this is either
        a *list of sample coordinates* used to draw colors from the map, or an
        *integer number of colors* to draw. If the latter, the sample
        coordinates are ``np.linspace(0, 1, samples)``. For example,
        ``Cycle('Reds', 5)`` divides the ``'Reds'`` colormap into five evenly
        spaced colors.
    name : str, optional
        Name of the resulting `~matplotlib.colors.ListedColormap` used to
        register the color cycle. Default name is ``'no_name'``.
    marker, alpha, dashes, linestyle, linewidth, markersize, markeredgewidth, markeredgecolor, markerfacecolor : list of specs, optional
        Lists of `~matplotlib.lines.Line2D` properties that can be added to
        the `~cycler.Cycler` instance. If the lists have unequal length, they
        will be filled to match the length of the longest list.  See
        `~matplotlib.axes.Axes.set_prop_cycle` for more info on cyclers.
        Also see the `line style reference \
<https://matplotlib.org/gallery/lines_bars_and_markers/line_styles_reference.html>`__,
        the `marker reference \
<https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/marker_reference.html>`__,
        and the `custom dashes reference \
<https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/line_demo_dash_control.html>`__.
    save : bool, optional
        Whether to save the `ListedColormap` associated with this cycle.
        See `ListedColormap.save`.
    save_kw : dict-like, optional
        Ignored if `save` is ``False``. Passed to `ListedColormap.save`
        for the `ListedColormap` associated with this cycle.
    **kwargs
        Passed to `Colormap` when the input is not already a `~cycler.Cycler`
        instance.

    Returns
    -------
    `~cycler.Cycler`
        A cycler instance that can be passed to
        `~matplotlib.axes.Axes.set_prop_cycle`.
    """  # noqa
    # Add properties
    props = {}
    nprops = 0
    for key, value in (
        ('marker', marker),
        ('alpha', alpha),
        ('dashes', dashes),
        ('linestyle', linestyle),
        ('linewidth', linewidth),
        ('markersize', markersize),
        ('markeredgewidth', markeredgewidth),
        ('markeredgecolor', markeredgecolor),
        ('markerfacecolor', markerfacecolor),
    ):
        if value is not None:
            if isinstance(value, str) or not np.iterable(value):
                raise ValueError(
                    f'Invalid {key!r} property {value!r}. '
                    f'Must be list or tuple of properties.'
                )
            nprops = max(nprops, len(value))
            props[key] = [*value]  # ensure mutable list
    # If args is non-empty, means we want color cycle; otherwise is black
    if not args:
        props['color'] = ['k']  # ensures property cycler is non empty
        if kwargs:
            _warn_proplot(f'Ignoring Cycle() keyword arg(s) {kwargs}.')
    # Merge cycler objects
    elif all(isinstance(arg, cycler.Cycler) for arg in args):
        if kwargs:
            _warn_proplot(f'Ignoring Cycle() keyword arg(s) {kwargs}.')
        if len(args) == 1:
            return args[0]
        else:
            props = {}
            for arg in args:
                for key, value in arg.by_key():
                    if key not in props:
                        props[key] = []
                    props[key].extend([*value])
            return cycler.cycler(**props)
    # Build and register a ListedColormap
    else:
        # Collect samples
        if args and isinstance(args[-1], Number):
            # means we want to sample existing colormaps or cycles
            args, N = args[:-1], args[-1]
        kwargs.setdefault('fade', 90)
        kwargs.setdefault('listmode', 'listed')
        cmap = Colormap(*args, **kwargs)  # the cmap object itself
        if isinstance(cmap, ListedColormap):
            colors = cmap.colors[:N]  # if N is None, does nothing
        else:
            N = _notNone(N, 10)
            if isinstance(N, Integral):
                x = np.linspace(0, 1, N)  # from edge to edge
            elif np.iterable(N) and all(
                    isinstance(item, Number) for item in N):
                x = np.array(N)
            else:
                raise ValueError(f'Invalid samples {N!r}.')
            N = len(x)
            colors = cmap(x)

        # Register and save the samples as a ListedColormap
        name = name or '_no_name'
        cmap = ListedColormap(colors, name=name, N=N)
        mcm.cmap_d[name] = cmap
        if save:
            save_kw = save_kw or {}
            cmap.save(**save_kw)

        # Add to property dict
        nprops = max(nprops, len(colors))
        props['color'] = [
            tuple(color) if not isinstance(color, str) else color
            for color in cmap.colors
        ]  # save the tupled version!

    # Build cycler, make sure lengths are the same
    for key, value in props.items():
        if len(value) < nprops:
            value[:] = [
                value[i % len(value)] for i in range(nprops)
            ]  # make loop double back
    cycle = cycler.cycler(**props)
    cycle.name = name
    return cycle


def Norm(norm, *args, **kwargs):
    """
    Return an arbitrary `~matplotlib.colors.Normalize` instance.
    Used to interpret the `norm` and `norm_kw` arguments when passed to any
    plotting method wrapped by `~proplot.wrappers.cmap_changer`.

    Parameters
    ----------
    norm : str or `~matplotlib.colors.Normalize`
        The normalizer specification. If a `~matplotlib.colors.Normalize`
        instance already, the input argument is simply returned. Otherwise,
        `norm` should be a string corresponding to one of the "registered"
        colormap normalizers (see below table).

        If `norm` is a list or tuple and the first element is a "registered"
        normalizer name, subsequent elements are passed to the normalizer class
        as positional arguments.

        ===============================  ===============================
        Key(s)                           Class
        ===============================  ===============================
        ``'midpoint'``, ``'zero'``       `MidpointNorm`
        ``'segmented'``, ``'segments'``  `LinearSegmentedNorm`
        ``'null'``, ``'none'``           `~matplotlib.colors.NoNorm`
        ``'linear'``                     `~matplotlib.colors.Normalize`
        ``'log'``                        `~matplotlib.colors.LogNorm`
        ``'power'``                      `~matplotlib.colors.PowerNorm`
        ``'symlog'``                     `~matplotlib.colors.SymLogNorm`
        ===============================  ===============================

    Other parameters
    ----------------
    *args, **kwargs
        Passed to the `~matplotlib.colors.Normalize` initializer.
        See `this tutorial \
<https://matplotlib.org/tutorials/colors/colormapnorms.html>`__
        for more info.

    Returns
    -------
    `~matplotlib.colors.Normalize`
        A `~matplotlib.colors.Normalize` instance.
    """
    if isinstance(norm, mcolors.Normalize):
        return norm

    # Pull out extra args
    if np.iterable(norm) and not isinstance(norm, str):
        norm, args = norm[0], (*norm[1:], *args)
    if not isinstance(norm, str):
        raise ValueError(f'Invalid norm name {norm!r}. Must be string.')

    # Get class
    if norm not in normalizers:
        raise ValueError(
            f'Unknown normalizer {norm!r}. Options are: '
            + ', '.join(map(repr, normalizers.keys())) + '.'
        )
    if norm == 'symlog' and not args and 'linthresh' not in kwargs:
        kwargs['linthresh'] = 1  # special case, needs argument
    return normalizers[norm](*args, **kwargs)


class DiscreteNorm(mcolors.BoundaryNorm):
    """
    Meta-normalizer that discretizes the possible color values returned by
    arbitrary continuous normalizers given a list of level boundaries. This
    is applied to all colormap plots in ProPlot.
    """
    # See this post: https://stackoverflow.com/a/48614231/4970632
    # WARNING: Must be child of BoundaryNorm. Many methods in ColorBarBase
    # test for class membership, crucially including _process_values(), which
    # if it doesn't detect BoundaryNorm will try to use DiscreteNorm.inverse().
    def __init__(
        self, levels, norm=None, step=1.0, extend=None,
        clip=False, descending=False,
    ):
        """
        Parameters
        ----------
        levels : list of float
            The level boundaries.
        norm : `~matplotlib.colors.Normalize`, optional
            The normalizer used to transform `levels` and all data passed
            to `~DiscreteNorm.__call__` before discretization. The ``vmin``
            and ``vmax`` of the normalizer are set to the minimum and
            maximum values in `levels`.
        step : float, optional
            The intensity of the transition to out-of-bounds colors as a
            fraction of the adjacent step between in-bounds colors.
            Default is ``1``.
        extend : {'neither', 'both', 'min', 'max'}, optional
            Which out-of-bounds regions should be assigned unique colormap
            colors. The normalizer needs this information so it can ensure
            the colorbar always spans the full range of colormap colors.
        clip : bool, optional
            Whether to clip values falling outside of the level bins. This
            only has an effect on lower colors when extend is
            ``'min'`` or ``'both'``, and on upper colors when extend is
            ``'max'`` or ``'both'``.
        descending : bool, optional
            Whether the levels are meant to be descending. This will cause
            the colorbar axis to be reversed when it is drawn with a
            `~matplotlib.cm.ScalarMappable` that uses this normalizer.

        Note
        ----
        If you are using a diverging colormap with ``extend='max'`` or
        ``extend='min'``, the center will get messed up. But that is very
        strange usage anyway... so please just don't do that :)
        """
        # Parse input
        # NOTE: This must be a subclass BoundaryNorm, so ColorbarBase will
        # detect it... even though we completely override it.
        if not norm:
            norm = mcolors.Normalize()
        elif isinstance(norm, mcolors.BoundaryNorm):
            raise ValueError(f'Normalizer cannot be instance of BoundaryNorm.')
        elif not isinstance(norm, mcolors.Normalize):
            raise ValueError('Normalizer must be instance of Normalize.')
        extend = extend or 'neither'
        extends = ('both', 'min', 'max', 'neither')
        if extend not in extends:
            raise ValueError(
                f'Unknown extend option {extend!r}. Options are: '
                + ', '.join(map(repr, extends)) + '.'
            )

        # Ensure monotonically increasing levels
        levels, _ = _check_levels(levels, allow_descending=False)
        bins, _ = _check_levels(norm(levels), allow_descending=False)
        self.N = levels.size
        self.clip = clip
        self.boundaries = levels
        self.vmin = norm.vmin = vmin = np.min(levels)
        self.vmax = norm.vmax = vmax = np.max(levels)
        vcenter = getattr(norm, 'vcenter', None)

        # Get color coordinates corresponding to each bin, plus extra
        # 2 color coordinates for out-of-bounds color bins.
        # For *same* out-of-bounds colors, looks like [0, 0, ..., 1, 1]
        # For *unique* out-of-bounds colors, looks like [0, X, ..., 1 - X, 1]
        # NOTE: Critical that we scale the bin centers in "physical space"
        # and *then* translate to color coordinates so that nonlinearities in
        # the normalization stay intact. If we scaled the bin centers in
        # *normalized space* to have minimum 0 maximum 1, would mess up
        # color distribution. However this is still not perfect... get
        # asymmetric color intensity either side of central point. So add
        # special handling for diverging norms below to improve symmetry.
        mids = np.zeros((levels.size + 1,))
        mids[1:-1] = 0.5 * (levels[1:] + levels[:-1])
        mids[0], mids[-1] = mids[1], mids[-2]
        if extend in ('min', 'both'):
            mids[0] += step * (mids[1] - mids[2])
        if extend in ('max', 'both'):
            mids[-1] += step * (mids[-2] - mids[-3])
        if vcenter is None:
            mids = _interpolate_basic(
                mids, np.min(mids), np.max(mids), vmin, vmax
            )
        else:
            mids = mids.copy()
            mids[mids < vcenter] = _interpolate_basic(
                mids[mids < vcenter], np.min(mids), vcenter, vmin, vcenter,
            )
            mids[mids >= vcenter] = _interpolate_basic(
                mids[mids >= vcenter], vcenter, np.max(mids), vcenter, vmax,
            )
        dest = norm(mids)

        # Attributes
        # NOTE: If clip is True, we clip values to the centers of the end
        # bins rather than vmin/vmax to prevent out-of-bounds colors from
        # getting an in-bounds bin color due to landing on a bin edge.
        # NOTE: With extend='min' the minimimum in-bounds and out-of-bounds
        # colors are the same so clip=True will have no effect. Same goes
        # for extend='max' with maximum colors.
        # WARNING: For some reason must clip manually for LogNorm, or
        # end up with unpredictable fill value, weird "out-of-bounds" colors
        self._bmin = np.min(mids)
        self._bmax = np.max(mids)
        self._bins = bins
        self._dest = dest
        self._norm = norm
        self._norm_clip = None
        self._descending = descending
        if isinstance(norm, mcolors.LogNorm):
            self._norm_clip = (5e-249, None)

    def __call__(self, value, clip=None):
        """
        Normalize data values to 0-1.

        Parameters
        ----------
        value : numeric
            The data to be normalized.
        clip : bool, optional
            Whether to clip values falling outside of the level bins.
            Default is ``self.clip``.
        """
        # Follow example of LinearSegmentedNorm, but perform no interpolation,
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
        return yq

    def inverse(self, value):  # noqa: U100
        """
        Raise an error. Inversion after discretization is impossible.
        """
        raise ValueError('DiscreteNorm is not invertible.')


class LinearSegmentedNorm(mcolors.Normalize):
    """
    This is the default normalizer paired with `BinNorm` whenever `levels`
    are non-linearly spaced. The normalized value is linear with respect to
    its average index in the `levels` vector, allowing uniform color
    transitions across arbitrarily spaced monotonically increasing values.
    Can be explicitly used by passing ``norm='segmented'`` to any command
    accepting ``cmap``.
    """
    def __init__(self, levels, vmin=None, vmax=None, clip=False):
        """
        Parameters
        ----------
        levels : list of float
            The discrete data levels.
        vmin, vmax : None
            Ignored. `vmin` and `vmax` are set to the minimum and
            maximum of `levels`.
        clip : bool, optional
            Whether to clip values falling outside of the minimum and
            maximum levels.
        """
        levels = np.atleast_1d(levels)
        diffs = np.sign(np.diff(levels))
        y = np.linspace(0, 1, len(levels))
        self._descending = False
        if levels.ndim != 1:
            raise ValueError('Levels must be 1-dimensional.')
        elif levels.size < 2:
            raise ValueError('Need at least two levels.')
        elif all(diffs == -1):
            self._descending = True
            levels = levels[::-1]
            y = y[::-1]
        elif not all(diffs == 1):
            raise ValueError(
                f'Levels {levels!r} must be monotonically increasing.'
            )
        vmin, vmax = levels.min(), levels.max()
        super().__init__(vmin, vmax, clip=clip)  # second level superclass
        self._x = levels
        self._y = y

    def __call__(self, value, clip=None):
        """
        Normalize the data values to 0-1. Inverse
        of `~LinearSegmentedNorm.inverse`.

        Parameters
        ----------
        value : numeric
            The data to be normalized.
        clip : bool, optional
            Whether to clip values falling outside of the minimum and
            maximum levels. Default is ``self.clip``.
        """
        # Follow example of make_mapping_array for efficient, vectorized
        # linear interpolation across multiple segments.
        # * Normal test puts values at a[i] if a[i-1] < v <= a[i]; for
        #   left-most data, satisfy a[0] <= v <= a[1]
        # * searchsorted gives where xq[i] must be inserted so it is larger
        #   than x[ind[i]-1] but smaller than x[ind[i]]
        if clip is None:  # builtin clipping
            clip = self.clip
        if clip:  # note that np.clip can handle masked arrays
            value = np.clip(value, self.vmin, self.vmax)
        x = self._x  # from arbitrarily spaced monotonic levels
        y = self._y  # to linear range 0-1
        xq = np.atleast_1d(value)
        idx = np.searchsorted(x, xq)
        idx[idx == 0] = 1
        idx[idx == len(x)] = len(x) - 1
        distance = (xq - x[idx - 1]) / (x[idx] - x[idx - 1])
        yq = distance * (y[idx] - y[idx - 1]) + y[idx - 1]
        if self._descending:
            yq = 1 - yq
        mask = ma.getmaskarray(xq)
        return ma.array(yq, mask=mask)

    def inverse(self, value):
        """
        Inverse operation of `~LinearSegmentedNorm.__call__`.

        Parameters
        ----------
        value : numeric
            The data to be un-normalized.
        """
        x = self._x
        y = self._y
        yq = np.atleast_1d(value)
        idx = np.searchsorted(y, yq)
        idx[idx == 0] = 1
        idx[idx == len(y)] = len(y) - 1
        distance = (yq - y[idx - 1]) / (y[idx] - y[idx - 1])
        xq = distance * (x[idx] - x[idx - 1]) + x[idx - 1]
        mask = ma.getmaskarray(yq)
        return ma.array(xq, mask=mask)


class DivergingNorm(mcolors.Normalize):
    """
    Normalizer that ensures some central data value lies at the central
    colormap color.  The default central value is ``0``. Can be used by
    passing ``norm='diverging'`` to any command accepting ``cmap``.
    """
    def __init__(
        self, vcenter=0, vmin=None, vmax=None, fair=True, clip=None
    ):
        """
        Parameters
        ----------
        vcenter : float, optional
            The data value corresponding to the central position of the
            colormap. The default is ``0``.
        vmin, vmax : float, optional
            The minimum and maximum data values.
        fair : bool, optional
            If ``True`` (default), the speeds of the color gradations on
            either side of the center point are equal, but colormap colors may
            be omitted. If ``False``, all colormap colors are included, but
            the color gradations on one side may be faster than the other side.
            ``False`` should be used with great care, as it may result in
            a misleading interpretation of your data.
        clip : bool, optional
            Whether to clip values falling outside of `vmin` and `vmax`.
        """
        # NOTE: This post is an excellent summary of matplotlib's DivergingNorm history:
        # https://github.com/matplotlib/matplotlib/issues/15336#issuecomment-535291287
        # NOTE: This is a stale PR that plans to implement the same features.
        # https://github.com/matplotlib/matplotlib/pull/15333#issuecomment-537545430
        # Since proplot is starting without matplotlib's baggage we can just implement
        # DivergingNorm like they would prefer if they didn't have to worry about
        # confusing users: single class, default "fair" scaling that can be turned off.
        super().__init__(vmin, vmax, clip)
        self.vmin = vmin
        self.vmax = vmax
        self.vcenter = vcenter
        self.fair = fair

    def __call__(self, value, clip=None):
        """
        Normalize data values to 0-1.

        Parameters
        ----------
        value : numeric
            The data to be normalized.
        clip : bool, optional
            Whether to clip values falling outside of `vmin` and `vmax`.
            Default is ``self.clip``.
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
        elif not self.fair:
            x = [self.vmin, self.vcenter, self.vmax]
            y = [0, 0.5, 1.0]
        else:
            offset = max(
                np.abs(self.vcenter - self.vmin),
                np.abs(self.vmax - self.vcenter),
            )
            x = [self.vcenter - offset, self.vcenter + offset]
            y = [0, 1.0]
        yq = _interpolate_extrapolate(xq, x, y)
        if is_scalar:
            yq = np.atleast_1d(yq)[0]
        return yq

    def autoscale_None(self, z):
        """
        Get vmin and vmax, and then clip at vcenter
        """
        super().autoscale_None(z)
        if self.vmin > self.vcenter:
            self.vmin = self.vcenter
        if self.vmax < self.vcenter:
            self.vmax = self.vcenter


class ColorDatabase(dict):
    """
    Dictionary subclass used to replace the builtin matplotlib color
    database. This allows users to draw colors from named colormaps and color
    cycles for any plotting command that accepts a `color` keyword arg.
    See `~ColorDatabase.cache` for details.
    """
    def __init__(self, mapping):
        """
        Parameters
        ----------
        mapping : dict-like
            The colors.
        """
        super().__init__(mapping)
        self._cache = _ColorCache({})

    def __setitem__(self, key, value):
        """
        Add a color to the database and clear the cache.
        """
        if not isinstance(key, str):
            raise ValueError(f'Invalid color name {key!r}. Must be string.')
        super().__setitem__(key, value)
        self.cache.clear()

    def __delitem__(self, key):
        """
        Delete a color from the database and clear the cache.
        """
        super().__delitem__(key)
        self.cache.clear()

    @property
    def cache(self):
        """
        A special dictionary subclass capable of retrieving colors
        "on-the-fly" from registered colormaps and color cycles.

        * For a smooth colormap, usage is e.g. ``color=('Blues', 0.8)``. The
          number is the colormap index, and must be between 0 and 1.
        * For a color cycle, usage is e.g. ``color=('colorblind', 2)``. The
          number is the list index.

        These examples work with any matplotlib command that accepts a `color`
        keyword arg.
        """
        return self._cache


class _ColorCache(dict):
    def __getitem__(self, key):
        # Matplotlib 'color' args are passed to to_rgba, which tries to read
        # directly from cache and if that fails, sanitizes input, which
        # raises error on receiving (colormap, idx) tuple. So we *have* to
        # override cache instead of color dict itself.
        rgb, alpha = key
        if (
            not isinstance(rgb, str) and np.iterable(rgb) and len(rgb) == 2
            and isinstance(rgb[1], Number) and isinstance(rgb[0], str)
        ):
            try:
                cmap = _cmapdict[rgb[0]]
            except (TypeError, KeyError):
                pass
            else:
                if isinstance(cmap, ListedColormap):
                    if not 0 <= rgb[1] < len(cmap.colors):
                        raise ValueError(
                            f'Color cycle sample for {rgb[0]!r} cycle must be '
                            f'between 0 and {len(cmap.colors)-1}, '
                            f'got {rgb[1]}.'
                        )
                    rgb = cmap.colors[rgb[1]]  # draw from list of colors
                else:
                    if not 0 <= rgb[1] <= 1:
                        raise ValueError(
                            f'Colormap sample for {rgb[0]!r} colormap must be '
                            f'between 0 and 1, got {rgb[1]}.'
                        )
                    rgb = cmap(rgb[1])  # get color selection
                rgba = mcolors.to_rgba(rgb, alpha)
                return rgba


def _get_cmap(name=None, lut=None):
    """
    Monkey patch for matplotlib `~matplotlib.get_cmap`. Permits case-insensitive
    search of monkey-patched colormap database (which was broken in v3.2.0).
    """
    if name is None:
        name = rcParams['image.cmap']
    if isinstance(name, mcolors.Colormap):
        return name
    try:
        cmap = _cmapdict[name]
    except KeyError:
        raise KeyError(
            f'Invalid colormap name {name!r}. Valid names are: '
            + ', '.join(map(repr, _cmapdict)) + '.'
        )
    if lut is not None:
        cmap = cmap._resample(lut)
    return cmap


class ColormapDatabase(dict):
    """
    Dictionary subclass used to replace the `matplotlib.cm.cmap_d`
    colormap dictionary. See `~ColormapDatabase.__getitem__` and
    `~ColormapDatabase.__setitem__` for details.
    """
    def __init__(self, kwargs):
        """
        Parameters
        ----------
        kwargs : dict-like
            The source dictionary.
        """
        for key, value in kwargs.items():
            self.__setitem__(key, value)

    def __delitem__(self, key):
        """
        Delete the item from the database.
        """
        key = self._sanitize_key(key, mirror=True)
        super().__delitem__(key)

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
        key = self._sanitize_key(key, mirror=True)
        shift = key[-2:] == '_s'
        if shift:
            key = key[:-2]
        reverse = key[-2:] == '_r'
        if reverse:
            key = key[:-2]
        value = super().__getitem__(key)  # may raise keyerror
        if shift:
            if hasattr(value, 'shifted'):
                value = value.shifted(180)
            else:
                raise KeyError(
                    f'Item of type {type(value).__name__!r} '
                    'does not have shifted() method.'
                )
        if reverse:
            if hasattr(value, 'reversed'):
                value = value.reversed()
            else:
                raise KeyError(
                    f'Item of type {type(value).__name__!r} '
                    'does not have reversed() method.'
                )
        return value

    def __setitem__(self, key, item):
        """
        Store the colormap under its lowercase name. If the colormap is
        a matplotlib `~matplotlib.colors.ListedColormap` or
        `~matplotlib.colors.LinearSegmentedColormap`, it is converted to the
        ProPlot `ListedColormap` or `LinearSegmentedColormap` subclass.
        """
        if not isinstance(key, str):
            raise KeyError(f'Invalid key {key!r}. Must be string.')
        if isinstance(item, (ListedColormap, LinearSegmentedColormap)):
            pass
        elif isinstance(item, mcolors.LinearSegmentedColormap):
            item = LinearSegmentedColormap(
                item.name, item._segmentdata, item.N, item._gamma
            )
        elif isinstance(item, mcolors.ListedColormap):
            item = ListedColormap(
                item.colors, item.name, item.N
            )
        else:
            raise ValueError(
                f'Invalid colormap {item}. Must be instance of '
                'matplotlib.colors.ListedColormap or '
                'matplotlib.colors.LinearSegmentedColormap.'
            )
        key = self._sanitize_key(key, mirror=False)
        super().__setitem__(key, item)

    def __contains__(self, item):
        """
        Test for membership using the sanitized colormap name.
        """
        try:  # by default __contains__ ignores __getitem__ overrides
            self.__getitem__(item)
            return True
        except KeyError:
            return False

    def _sanitize_key(self, key, mirror=True):
        """
        Return the sanitized colormap name. This is used for lookups *and*
        assignments.
        """
        if not isinstance(key, str):
            raise KeyError(f'Invalid key {key!r}. Key must be a string.')
        key = key.lower()
        key = re.sub(r'\A(grays)(?:_r|_s)?\Z', 'greys', key)
        reverse = key[-2:] == '_r'
        if reverse:
            key = key[:-2]
        if mirror and not super().__contains__(key):  # search for mirrored key
            key_mirror = key
            for pair in CMAPS_DIVERGING:
                try:
                    idx = pair.index(key)
                    key_mirror = pair[1 - idx]
                except (ValueError, KeyError):
                    continue
            if super().__contains__(key_mirror):
                reverse = not reverse
                key = key_mirror
        if reverse:
            key = key + '_r'
        return key


# Replace color database with custom database
if not isinstance(mcolors._colors_full_map, ColorDatabase):
    _map = ColorDatabase(mcolors._colors_full_map)
    mcolors._colors_full_map = _map
    mcolors.colorConverter.cache = _map.cache
    mcolors.colorConverter.colors = _map

# Replace colormap database with custom database
if mcm.get_cmap is not _get_cmap:
    mcm.get_cmap = _get_cmap
if not isinstance(_cmapdict, ColormapDatabase):
    _cmapdict = {
        key: value for key, value in _cmapdict.items()
        if key[-2:] != '_r' and key[-8:] != '_shifted'
    }
    _cmapdict = ColormapDatabase(_cmapdict)
    setattr(mcm, _cmapdict_attr, _cmapdict)

#: Dictionary of possible normalizers. See `Norm` for a table.
normalizers = {
    'none': mcolors.NoNorm,
    'null': mcolors.NoNorm,
    'zero': MidpointNorm,
    'midpoint': MidpointNorm,
    'segments': LinearSegmentedNorm,
    'segmented': LinearSegmentedNorm,
    'log': mcolors.LogNorm,
    'linear': mcolors.Normalize,
    'power': mcolors.PowerNorm,
    'symlog': mcolors.SymLogNorm,
}

# Deprecations
CmapDict = warnings._rename_obj('CmapDict', ColormapDatabase)
ColorDict = warnings._rename_obj('ColorDict', ColorDatabase)
BinNorm = warnings._rename_obj('BinNorm', DiscreteNorm)
