#!/usr/bin/env python3
"""
New colormap classes and colormap normalization classes.
"""
# NOTE: Avoid colormap/color name conflicts by checking
# set(pplt.colors._cmap_database) & set(pplt.colors.mcolors._colors_full_map)
# whenever new default colormaps are added. Currently result is
# {'gray', 'marine', 'ocean', 'pink'} which correspond to MATLAB and GNUplot maps.
import json
import os
import re
from numbers import Integral, Number
from xml.etree import ElementTree

import matplotlib.cm as mcm
import matplotlib.colors as mcolors
import numpy as np
import numpy.ma as ma
from matplotlib import rcParams

from .internals import ic  # noqa: F401
from .internals import _not_none, _pop_props, docstring, warnings
from .utils import to_rgb, to_rgba, to_xyz, to_xyza

__all__ = [
    'ListedColormap',
    'LinearSegmentedColormap',
    'PerceptuallyUniformColormap',
    'DiscreteNorm',
    'DivergingNorm',
    'LinearSegmentedNorm',
    'ColorDatabase',
    'ColormapDatabase',
]

# Constants
# NOTE: Do not compile hex regex because config.py needs this with \A\Z surrounding
REGEX_HEX = r'#(?:[0-9a-fA-F]{3,4}){2}'  # 6-8 digit hex
CMAPS_DIVERGING = tuple(
    (key1.lower(), key2.lower())
    for key1, key2 in (
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
    )
)

# Deprecations
_cmaps_removed = {
    'Blue0': '0.6',
    'Cool': '0.6',
    'Warm': '0.6',
    'Hot': '0.6',
    'Floral': '0.6',
    'Contrast': '0.6',
    'Sharp': '0.6',
    'Viz': '0.6',
}
_cmaps_renamed = {
    'Blue1': ('Blues1', '0.7'),
    'Blue2': ('Blues2', '0.7'),
    'Blue3': ('Blues3', '0.7'),
    'Blue4': ('Blues4', '0.7'),
    'Blue5': ('Blues5', '0.7'),
    'Blue6': ('Blues6', '0.7'),
    'Blue7': ('Blues7', '0.7'),
    'Blue8': ('Blues8', '0.7'),
    'Blue9': ('Blues9', '0.7'),
    'Green1': ('Greens1', '0.7'),
    'Green2': ('Greens2', '0.7'),
    'Green3': ('Greens3', '0.7'),
    'Green4': ('Greens4', '0.7'),
    'Green5': ('Greens5', '0.7'),
    'Green6': ('Greens6', '0.7'),
    'Green7': ('Greens7', '0.7'),
    'Green8': ('Greens8', '0.7'),
    'Orange1': ('Yellows1', '0.7'),
    'Orange2': ('Yellows2', '0.7'),
    'Orange3': ('Yellows3', '0.7'),
    'Orange4': ('Oranges2', '0.7'),
    'Orange5': ('Oranges1', '0.7'),
    'Orange6': ('Oranges3', '0.7'),
    'Orange7': ('Oranges4', '0.7'),
    'Orange8': ('Yellows4', '0.7'),
    'Brown1': ('Browns1', '0.7'),
    'Brown2': ('Browns2', '0.7'),
    'Brown3': ('Browns3', '0.7'),
    'Brown4': ('Browns4', '0.7'),
    'Brown5': ('Browns5', '0.7'),
    'Brown6': ('Browns6', '0.7'),
    'Brown7': ('Browns7', '0.7'),
    'Brown8': ('Browns8', '0.7'),
    'Brown9': ('Browns9', '0.7'),
    'RedPurple1': ('Reds1', '0.7'),
    'RedPurple2': ('Reds2', '0.7'),
    'RedPurple3': ('Reds3', '0.7'),
    'RedPurple4': ('Reds4', '0.7'),
    'RedPurple5': ('Reds5', '0.7'),
    'RedPurple6': ('Purples1', '0.7'),
    'RedPurple7': ('Purples2', '0.7'),
    'RedPurple8': ('Purples3', '0.7'),
}

docstring.snippets['cmap.init'] = """
alpha : float, optional
    The opacity for the entire colormap. This overrides the input
    segment data.
cyclic : bool, optional
    Whether the colormap is cyclic. If ``True``, this changes how the leftmost
    and rightmost color levels are selected, and `extend` can only be
    ``'neither'`` (a warning will be issued otherwise).
"""

docstring.snippets['cmap.gamma'] = """
gamma : float, optional
    Sets `gamma1` and `gamma2` to this identical value.
gamma1 : float, optional
    If >1, makes low saturation colors more prominent. If <1,
    makes high saturation colors more prominent. Similar to the
    `HCLWizard <http://hclwizard.org:64230/hclwizard/>`_ option.
gamma2 : float, optional
    If >1, makes high luminance colors more prominent. If <1,
    makes low luminance colors more prominent. Similar to the
    `HCLWizard <http://hclwizard.org:64230/hclwizard/>`_ option.
"""


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


def _get_channel(color, channel, space='hcl'):
    """
    Get the hue, saturation, or luminance channel value from the input color.
    The color name `color` can optionally be a string with the format
    ``'color+x'`` or ``'color-x'``, where `x` specifies the offset from the
    channel value.

    Parameters
    ----------
    color : color-spec
        The color. Sanitized with `to_rgba`.
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


def _make_segment_data(values, coords=None, ratios=None):
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


def _make_lookup_table(N, data, gamma=1.0, inverse=False):
    r"""
    Used to generate lookup tables of HSL values given gradations specified
    by `PerceptuallyUniformColormap`. Similar to `~matplotlib.colors.makeMappingArray`
    but permits *circular* hue gradations along 0-360, disables clipping of
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
        `PerceptuallyUniformColormap`\ s because they usually consist of just
        one linear transition for *sequential* colormaps and two linear
        transitions for *diverging* colormaps.
    inverse : bool, optional
        If ``True``, :math:`w_i^{\gamma_i}` is replaced with
        :math:`1 - (1 - w_i)^{\gamma_i}` -- that is, when `gamma` is greater
        than 1, this weights colors toward *higher* channel values instead
        of lower channel values.

        This is implemented in case we want to apply *equal* "gamma scaling"
        to different HSL channels in different directions. Usually, this
        is done to weight low data values with higher luminance *and* lower
        saturation, thereby emphasizing "extreme" data values.
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

    @classmethod
    def _from_file(cls, filename, warn_on_failure=False):
        """
        Read generalized colormap and color cycle files.
        """
        filename = os.path.expanduser(filename)
        name, ext = os.path.splitext(os.path.basename(filename))
        listed = issubclass(cls, mcolors.ListedColormap)
        reversed = name[-2:] == '_r'

        # Warn if loading failed during `register_cmaps` or `register_cycles`
        # but raise error if user tries to load a file.
        def _warn_or_raise(msg, error=RuntimeError):
            if warn_on_failure:
                warnings._warn_proplot(msg)
            else:
                raise error(msg)
        if not os.path.exists(filename):
            return _warn_or_raise(f'File {filename!r} not found.', FileNotFoundError)

        # Directly read segmentdata json file
        # NOTE: This is special case! Immediately return name and cmap
        ext = ext[1:]
        if ext == 'json':
            if listed:
                raise TypeError(
                    f'Cannot load listed colormaps from json files ({filename!r}).'
                )
            try:
                with open(filename, 'r') as fh:
                    data = json.load(fh)
            except json.JSONDecodeError:
                return _warn_or_raise(
                    f'Failed to load {filename!r}.', json.JSONDecodeError
                )
            kw = {}
            for key in ('cyclic', 'gamma', 'gamma1', 'gamma2', 'space'):
                if key in data:
                    kw[key] = data.pop(key, None)
            if 'red' in data:
                cmap = LinearSegmentedColormap(name, data)
            else:
                cmap = PerceptuallyUniformColormap(name, data, **kw)
            if reversed:
                cmap = cmap.reversed(name[:-2])
            return cmap

        # Read .rgb and .rgba files
        if ext in ('txt', 'rgb'):
            # Load
            # NOTE: This appears to be biggest import time bottleneck! Increases
            # time from 0.05s to 0.2s, with numpy loadtxt or with this regex thing.
            delim = re.compile(r'[,\s]+')
            data = [
                delim.split(line.strip())
                for line in open(filename)
                if line.strip() and line.strip()[0] != '#'
            ]
            try:
                data = [[float(num) for num in line] for line in data]
            except ValueError:
                return _warn_or_raise(
                    f'Failed to load {filename!r}. Expected a table of comma '
                    'or space-separated values.'
                )
            # Build x-coordinates and standardize shape
            data = np.array(data)
            if data.shape[1] not in (3, 4):
                return _warn_or_raise(
                    f'Failed to load {filename!r}. Got {data.shape[1]} columns, '
                    f'but expected 3 or 4.'
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
                doc = ElementTree.parse(filename)
            except ElementTree.ParseError:
                return _warn_or_raise(
                    f'Failed to load {filename!r}. Parsing error.',
                    ElementTree.ParseError
                )
            x, data = [], []
            for s in doc.getroot().findall('.//Point'):
                # Verify keys
                if any(key not in s.attrib for key in 'xrgb'):
                    return _warn_or_raise(
                        f'Failed to load {filename!r}. Missing an x, r, g, or b '
                        'specification inside one or more <Point> tags.'
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
                len(data[0]) == len(color) and len(color) in (3, 4)
                for color in data
            ):
                return _warn_or_raise(
                    f'Failed to load {filename!r}. Unexpected number of channels '
                    'or mixed channels across <Point> tags.'
                )

        # Read hex strings
        elif ext == 'hex':
            # Read arbitrary format
            string = open(filename).read()  # into single string
            data = re.findall(REGEX_HEX, string)
            if len(data) < 2:
                return _warn_or_raise(
                    f'Failed to load {filename!r}. Hex strings not found.'
                )
            # Convert to array
            x = np.linspace(0, 1, len(data))
            data = [to_rgb(color) for color in data]

        # Invalid extension
        else:
            return _warn_or_raise(
                f'Colormap or cycle file {filename!r} has unknown extension.'
            )

        # Standardize and reverse if necessary to cmap
        # TODO: Document the fact that filenames ending in _r return a reversed
        # version of the colormap stored in that file.
        x, data = np.array(x), np.array(data)
        x = (x - x.min()) / (x.max() - x.min())  # ensure they span 0-1
        if np.any(data > 2):  # from 0-255 to 0-1
            data = data / 255
        if reversed:
            name = name[:-2]
            data = data[::-1, :]
            x = 1 - x[::-1]
        if listed:
            return ListedColormap(data, name)
        else:
            data = [(x, color) for x, color in zip(x, data)]
            return LinearSegmentedColormap.from_list(name, data)


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
                string += (
                    f' {key!r}: [{data[0][2]:.3f}, ..., {data[-1][1]:.3f}],\n'
                )
        return type(self).__name__ + '({\n' + string + '})'

    @docstring.add_snippets
    def __init__(
        self, name, segmentdata, N=None, gamma=1,
        cyclic=False, alpha=None,
    ):
        """
        Parameters
        ----------
        name : str
            The colormap name.
        segmentdata : dict-like
            Dictionary containing the keys ``'red'``, ``'blue'``, ``'green'``,
            and (optionally) ``'alpha'``. The shorthands ``'r'``, ``'g'``, ``'b'``,
            and ``'a'`` are also acceptable. The key values can be callable functions
            that return channel values given a colormap index, or 3-column arrays
            indicating the coordinates and channel transitions. See
            `matplotlib.colors.LinearSegmentedColormap` for a detailed explanation.
        N : int, optional
            Number of points in the colormap lookup table. Default is :rc:`image.lut`.
        gamma : float, optional
            Gamma scaling used for the *x* coordinates.
        %(cmap.init)s

        See also
        --------
        ListedColormap
        matplotlib.colors.LinearSegmentedColormap
        proplot.constructor.Colormap
        """
        N = _not_none(N, rcParams['image.lut'])
        data = _pop_props(segmentdata, 'rgba', 'hsla')
        if segmentdata:
            raise ValueError(f'Invalid segmentdata keys {tuple(segmentdata)}.')
        super().__init__(name, data, N=N, gamma=gamma)
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
            Instances of `LinearSegmentedColormap`.
        ratios : list of float, optional
            Relative extent of each component colormap in the merged colormap.
            Length must equal ``len(args) + 1``.

            For example, ``cmap1.append(cmap2, ratios=(2, 1))`` generates
            a colormap with the left two-thrids containing colors from
            ``cmap1`` and the right one-third containing colors from ``cmap2``.
        name : str, optional
            The name of the new colormap. Default is
            ``'_'.join(cmap.name for cmap in args)``.
        N : int, optional
            The number of points in the colormap lookup table.
            Default is :rc:`image.lut` times ``len(args)``.

        Other parameters
        ----------------
        **kwargs
            Passed to `LinearSegmentedColormap.copy`
            or `PerceptuallyUniformColormap.copy`.

        Returns
        -------
        `LinearSegmentedColormap`
            The colormap.

        See also
        --------
        ListedColormap.append
        """
        # Parse input args
        if not args:
            return self
        if not all(isinstance(cmap, mcolors.LinearSegmentedColormap) for cmap in args):
            raise TypeError(f'Arguments {args!r} must be LinearSegmentedColormaps.')

        # PerceptuallyUniformColormap --> LinearSegmentedColormap conversions
        cmaps = [self, *args]
        spaces = {getattr(cmap, '_space', None) for cmap in cmaps}
        to_linear_segmented = len(spaces) > 1  # mixed colorspaces *or* mixed types
        if to_linear_segmented:
            for i, cmap in enumerate(cmaps):
                if isinstance(cmap, PerceptuallyUniformColormap):
                    cmaps[i] = cmap.to_linear_segmented()

        # Combine the segmentdata, and use the y1/y2 slots at merge points so
        # we never interpolate between end colors of different colormaps
        segmentdata = {}
        if name is None:
            name = '_'.join(cmap.name for cmap in cmaps)
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
            callable_ = [callable(cmap._segmentdata[key]) for cmap in cmaps]
            if all(callable_):  # expand range from x-to-w to 0-1
                funcs = [cmap._segmentdata[key] for cmap in cmaps]
                def xyy(ix, funcs=funcs):  # noqa: E306
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
                raise TypeError('Mixed callable and non-callable colormap values.')
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
            for cmap in cmaps:
                igamma = getattr(cmap, '_' + ikey)
                if not np.iterable(igamma):
                    if all(callable_):
                        igamma = (igamma,)
                    else:
                        igamma = (igamma,) * (len(cmap._segmentdata[key]) - 1)
                gamma.extend(igamma)
            if all(callable_):
                if any(igamma != gamma[0] for igamma in gamma[1:]):
                    warnings._warn_proplot(
                        'Cannot use multiple segment gammas when concatenating '
                        f'callable segments. Using the first gamma of {gamma[0]}.'
                    )
                gamma = gamma[0]
            kwargs[ikey] = gamma

        # Return copy or merge mixed types
        if to_linear_segmented and isinstance(self, PerceptuallyUniformColormap):
            return LinearSegmentedColormap(name, segmentdata, N, **kwargs)
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
        name : str, optional
            The name of the new colormap. Default is ``self.name + '_copy'``.
        left, right : float, optional
            The colormap indices for the "leftmost" and "rightmost" colors.
            Defaults are ``0`` and ``1``. See
            `~LinearSegmentedColormap.truncate` for details.
        right : float, optional
            The colormap index for the new "rightmost" color. Must fall between

        Other parameters
        ----------------
        **kwargs
            Passed to `LinearSegmentedColormap.copy`
            or `PerceptuallyUniformColormap.copy`.

        Returns
        -------
        `LinearSegmentedColormap`
            The colormap.

        See also
        --------
        LinearSegmentedColormap.truncate
        ListedColormap.truncate
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
            raise ValueError(
                f'Invalid combination cut={cut} for left={left} and right={right}.'
            )
        if name is None:
            name = self.name + '_copy'
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
            args.append(type(self)('_no_name', segmentdata, self.N))
            kwargs.setdefault('ratios', (ratio, abs(cut), ratio))
        args.append(cmap_right)

        return cmap_left.append(*args, name=name, **kwargs)

    def reversed(self, name=None, **kwargs):
        """
        Return a reversed copy of the colormap.

        Parameters
        ----------
        name : str, optional
            The name of the new colormap. Default is ``self.name + '_r'``.

        Other parameters
        ----------------
        **kwargs
            Passed to `LinearSegmentedColormap.copy`
            or `PerceptuallyUniformColormap.copy`.

        See also
        --------
        matplotlib.colors.LinearSegmentedColormap.reversed
        """
        segmentdata = {
            key: (lambda x: data(1.0 - x)) if callable(data)
            else [(1.0 - x, y1, y0) for x, y0, y1 in reversed(data)]
            for key, data in self._segmentdata.items()
        }
        for key in ('gamma1', 'gamma2'):
            if key in kwargs:
                continue
            gamma = getattr(self, '_' + key, None)
            if gamma is not None and np.iterable(gamma):
                kwargs[key] = gamma[::-1]
        if name is None:
            name = self.name + '_r'
        return self.copy(name, segmentdata, **kwargs)

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

            ===================  ==========================================
            Extension            Description
            ===================  ==========================================
            ``.json`` (default)  JSON database of the channel segment data.
            ``.hex``             Comma-delimited list of HEX strings.
            ``.rgb``, ``.txt``   3-4 column table of channel values.
            ===================  ==========================================

        alpha : bool, optional
            Whether to include an opacity column for ``.rgb``
            and ``.txt`` files.

        See also
        --------
        ListedColormap.save
        """
        dirname = os.path.join('~', '.proplot', 'cmaps')
        filename = self._parse_path(path, dirname, 'json')

        # Save channel segment data in json file
        _, ext = os.path.splitext(filename)
        if ext[1:] == 'json':
            # Sanitize segmentdata values. Convert np.float to builtin float,
            # np.array to list of lists, and callable to list of lists. We
            # tried encoding func.__code__ with base64 and marshal instead, but
            # when cmap.concatenate() embeds functions as keyword arguments, this
            # seems to make it *impossible* to load back up the function with
            # FunctionType (error message: arg 5 (closure) must be tuple). Instead
            # use this brute force workaround.
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
            with open(filename, 'w') as fh:
                json.dump(data, fh, indent=4)

        # Save lookup table colors
        else:
            data = self._get_data(ext[1:], alpha=alpha)
            with open(filename, 'w') as fh:
                fh.write(data)
        print(f'Saved colormap to {filename!r}.')

    def set_alpha(self, alpha, coords=None, ratios=None):
        """
        Set the opacity for the entire colormap or set up an opacity gradation.

        Parameters
        ----------
        alpha : float or list of float
            If float, this is the opacity for the entire colormap. If list of
            float, the colormap traverses these opacity values.
        coords : list of float, optional
            Colormap coordinates for the opacity values. The first and last
            coordinates must be ``0`` and ``1``. If `alpha` is not scalar, the
            default coordinates are ``np.linspace(0, 1, len(alpha))``.
        ratios : list of float, optional
            Relative extent of each opacity transition segment. Length should
            equal ``len(alpha) + 1``. For example
            ``cmap.set_alpha((1, 1, 0), ratios=(2, 1))`` creates a transtion from
            100 percent to 0 percent opacity in the right *third* of the colormap.

        See also
        --------
        ListedColormap.set_alpha
        """
        alpha = _make_segment_data(alpha, coords=coords, ratios=ratios)
        self._segmentdata['alpha'] = alpha
        self._isinit = False

    def set_cyclic(self, b):
        """
        Set whether this colormap is "cyclic". See `LinearSegmentedColormap`
        for details.
        """
        self._cyclic = bool(b)
        self._isinit = False

    def shifted(self, shift=180, name=None, **kwargs):
        """
        Return a cyclicaly shifted version of the colormap. If the colormap
        cyclic property is set to ``False`` a warning will be raised.

        Parameters
        ----------
        shift : float, optional
            The number of degrees to shift, out of 360 degrees.
            The default is ``180``.
        name : str, optional
            The name of the new colormap. Default is ``self.name + '_s'``.

        Other parameters
        ----------------
        **kwargs
            Passed to `LinearSegmentedColormap.copy`
            or `PerceptuallyUniformColormap.copy`.

        See also
        --------
        ListedColormap.shifted
        """
        shift = ((shift or 0) / 360) % 1
        if shift == 0:
            return self
        if name is None:
            name = self.name + '_s'
        if not self._cyclic:
            warnings._warn_proplot(
                f'Shifting non-cyclic colormap {self.name!r}. '
                f'Use cmap.set_cyclic(True) or Colormap(..., cyclic=True) to '
                'suppress this warning.'
            )
            self._cyclic = True

        # Decompose shift into two truncations followed by concatenation
        cmap_left = self.truncate(shift, 1)
        cmap_right = self.truncate(0, shift)
        return cmap_left.append(
            cmap_right, ratios=(1 - shift, shift), name=name, **kwargs
        )

    def truncate(self, left=None, right=None, name=None, **kwargs):
        """
        Return a truncated version of the colormap.

        Parameters
        ----------
        left : float, optional
            The colormap index for the new "leftmost" color. Must fall between ``0``
            and ``1``. For example, ``left=0.1`` cuts the leftmost 10%% of the colors.
        right : float, optional
            The colormap index for the new "rightmost" color. Must fall between ``0``
            and ``1``. For example, ``right=0.9`` cuts the leftmost 10%% of the colors.
        name : str, optional
            The name of the new colormap. Default is ``self.name + '_copy'``.

        Other parameters
        ----------------
        **kwargs
            Passed to `LinearSegmentedColormap.copy`
            or `PerceptuallyUniformColormap.copy`.

        See also
        --------
        ListedColormap.truncate
        """
        # Bail out
        left = max(_not_none(left, 0), 0)
        right = min(_not_none(right, 1), 1)
        if left == 0 and right == 1:
            return self
        if name is None:
            name = self.name + '_copy'

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
        name : str
            The name of the new colormap. Default is ``self.name + '_copy'``.
        segmentdata, N, alpha, gamma, cyclic : optional
            See `LinearSegmentedColormap`. If not provided, these are copied from
            the current colormap.

        See also
        --------
        ListedColormap.copy
        PerceptuallyUniformColormap.copy
        """
        if name is None:
            name = self.name + '_copy'
        if segmentdata is None:
            segmentdata = self._segmentdata.copy()
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

    def to_listed(self, samples=10, name=None, **kwargs):
        """
        Convert the `LinearSegmentedColormap` to a `ListedColormap` by drawing
        samples from the colormap.

        Parameters
        ----------
        samples : int or list of float, optional
            If integer, draw samples at the colormap coordinates
            ``np.linspace(0, 1, samples)``. If list of float, draw samples
            at the specified points.
        name : str, optional
            The name of the new colormap. Default is ``self.name + '_copy'``.

        Other parameters
        ----------------
        **kwargs
            Passed to `ListedColormap`.

        See also
        --------
        ListedColormap.to_linear_segmented
        """
        if isinstance(samples, Integral):
            samples = np.linspace(0, 1, samples)
        elif not np.iterable(samples):
            raise TypeError('Samples must be integer or iterable.')
        samples = np.asarray(samples)
        colors = self(samples)
        if name is None:
            name = self.name + '_copy'
        return ListedColormap(colors, name=name, **kwargs)

    @classmethod
    def from_file(cls, path, warn_on_failure=False):
        """
        Load colormap from a file.

        Parameters
        ----------
        path : str
            The file path. The file extension should be one of the following:

            ===================  ==========================================
            Extension            Description
            ===================  ==========================================
            ``.json``            JSON database of the channel segment data.
            ``.hex``             Comma-delimited list of HEX strings.
            ``.rgb``, ``.txt``   3-4 column table of channel values.
            ===================  ==========================================

        warn_on_failure : bool, optional
            If ``True``, issue a warning when loading fails instead of
            raising an error.

        See also
        --------
        ListedColormap.from_file
        """
        return cls._from_file(path, warn_on_failure=warn_on_failure)

    @classmethod
    def from_list(cls, name, colors, ratios=None, **kwargs):
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

        See also
        --------
        matplotlib.colors.LinearSegmentedColormap.from_list
        PerceptuallyUniformColormap.from_list
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
        colors = [to_rgba(color) for color in colors]

        # Build segmentdata
        keys = ('red', 'green', 'blue', 'alpha')
        cdict = {}
        for key, values in zip(keys, zip(*colors)):
            cdict[key] = _make_segment_data(values, coords, ratios)
        return cls(name, cdict, **kwargs)

    # Rename methods
    concatenate, punched, truncated, updated = warnings._rename_objs(
        '0.6',
        concatenate=append,
        punched=cut,
        truncated=truncate,
        updated=copy,
    )


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
            '})'
        )

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

        See also
        --------
        LinearSegmentedColormap
        matplotlib.colors.ListedColormap
        proplot.constructor.Colormap
        """
        super().__init__(*args, **kwargs)
        if alpha is not None:
            self.set_alpha(alpha)

    def append(self, *args, name=None, N=None, **kwargs):
        """
        Append arbitrary colormaps onto this colormap.

        Parameters
        ----------
        *args
            Instances of `ListedColormap`.
        name : str, optional
            The name of the new colormap. Default is
            ``'_'.join(cmap.name for cmap in args)``.
        N : int, optional
            The number of colors in the colormap lookup table. Default is
            the number of colors in the concatenated lists.

        Other parameters
        ----------------
        **kwargs
            Passed to `~ListedColormap.copy`.

        See also
        --------
        LinearSegmentedColormap.append
        """
        if not args:
            return self
        if not all(isinstance(cmap, mcolors.ListedColormap) for cmap in args):
            raise TypeError(f'Arguments {args!r} must be ListedColormap.')
        cmaps = (self, *args)
        if name is None:
            name = '_'.join(cmap.name for cmap in cmaps)
        colors = [color for cmap in cmaps for color in cmap.colors]
        return self.copy(colors, name, N or len(colors), **kwargs)

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

            ==================  ====================================
            Extension           Description
            ==================  ====================================
            ``.hex`` (default)  Comma-delimited list of HEX strings.
            ``.rgb``, ``.txt``  3-4 column table of channel values.
            ==================  ====================================

        alpha : bool, optional
            Whether to include an opacity column for ``.rgb``
            and ``.txt`` files.

        See also
        --------
        LinearSegmentedColormap.save
        """
        dirname = os.path.join('~', '.proplot', 'cycles')
        filename = self._parse_path(path, dirname, 'hex')

        # Save lookup table colors
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
        LinearSegmentedColormap.set_alpha
        """
        colors = [list(mcolors.to_rgba(color)) for color in self.colors]
        for color in colors:
            color[3] = alpha
        self.colors = colors
        self._init()

    def shifted(self, shift=1, name=None):
        """
        Return a cyclically shifted version of the colormap.

        Parameters
        ----------
        shift : float, optional
            The number of places to shift, between ``-self.N`` and ``self.N``.
            The default is ``1``.
        name : str, optional
            The name of the new colormap. Default is ``self.name + '_s'``.

        See also
        --------
        LinearSegmentedColormap.shifted
        """
        if not shift:
            return self
        if name is None:
            name = self.name + '_s'
        shift = shift % len(self.colors)
        colors = list(self.colors)
        colors = colors[shift:] + colors[:shift]
        return self.copy(colors, name, len(colors))

    def truncate(self, left=None, right=None, name=None):
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
            The name of the new colormap. Default is ``self.name + '_copy'``.

        See also
        --------
        LinearSegmentedColormap.truncate
        """
        if left is None and right is None:
            return self
        if name is None:
            name = self.name + '_copy'
        colors = self.colors[left:right]
        return self.copy(colors, name, len(colors))

    def copy(self, colors=None, name=None, N=None, *, alpha=None):
        """
        Return a new colormap with relevant properties copied from this one
        if they were not provided as keyword arguments.

        Parameters
        ----------
        name : str
            The name of the new colormap. Default is ``self.name + '_copy'``.
        colors, N, alpha : optional
            See `ListedColormap`. If not provided,
            these are copied from the current colormap.

        See also
        --------
        LinearSegmentedColormap.copy
        PerceptuallyUniformColormap.copy
        """
        if name is None:
            name = self.name + '_copy'
        if colors is None:
            colors = list(self.colors)  # copy
        if N is None:
            N = self.N
        cmap = ListedColormap(colors, name, N, alpha=alpha)
        cmap._rgba_bad = self._rgba_bad
        cmap._rgba_under = self._rgba_under
        cmap._rgba_over = self._rgba_over
        return cmap

    @classmethod
    def from_file(cls, path, warn_on_failure=False):
        """
        Load color cycle from a file.

        Parameters
        ----------
        path : str
            The file path. The file extension should be one of the following:

            ==================  ==========================================
            Extension           Description
            ==================  ==========================================
            ``.hex``            Comma-delimited list of HEX strings.
            ``.rgb``, ``.txt``  3-4 column table of channel values.
            ==================  ==========================================

        warn_on_failure : bool, optional
            If ``True``, issue a warning when loading fails instead of
            raising an error.

        See also
        --------
        LinearSegmentedColormap.from_file
        """
        return cls._from_file(path, warn_on_failure=warn_on_failure)

    # Rename methods
    concatenate, truncated, updated = warnings._rename_objs(
        '0.6',
        concatenate=append,
        truncated=truncate,
        updated=copy,
    )


class PerceptuallyUniformColormap(LinearSegmentedColormap, _Colormap):
    """
    Similar to `~matplotlib.colors.LinearSegmentedColormap`, but instead
    of varying the RGB channels, we vary hue, saturation, and luminance in
    either the HCL colorspace or the HSL or HPL scalings of HCL.
    """
    @docstring.add_snippets
    def __init__(
        self, name, segmentdata, N=None,
        space=None, clip=True, gamma=None, gamma1=None, gamma2=None, **kwargs
    ):
        """
        Parameters
        ----------
        name : str
            The colormap name.
        segmentdata : dict-like
            Dictionary containing the keys ``'hue'``, ``'saturation'``,
            ``'luminance'``, and (optionally) ``'alpha'``. The key ``'chroma'`` is
            treated as a synonym for ``'saturation'``. The shorthands ``'h'``,
            ``'s'``, ``'l'``, ``'a'``, and ``'c'`` are also acceptable. The key
            values can be callable functions that return channel values given a
            colormap index, or 3-column arrays indicating the coordinates and
            channel transitions. See `~matplotlib.colors.LinearSegmentedColormap`
            for a more detailed explanation.
        space : {'hsl', 'hpl', 'hcl'}, optional
            The hue, saturation, luminance-style colorspace to use for
            interpreting the channels. See
            `this page <http://www.hsluv.org/comparison/>`__ for a description.
        clip : bool, optional
            Whether to "clip" impossible colors, i.e. truncate HCL colors
            with RGB channels with values >1, or mask them out as gray.
        %(cmap.init)s
        %(cmap.gamma)s

        Other parameters
        ----------------
        **kwargs
            Passed to `LinearSegmentedColormap`.

        Example
        -------
        The below example generates a `PerceptuallyUniformColormap` from a
        `segmentdata` dictionary that uses color names for the hue data,
        instead of channel values between ``0`` and ``360``.

        >>> import proplot as pplt
        >>> data = {
        >>>     'h': [[0, 'red', 'red'], [1, 'blue', 'blue']],
        >>>     's': [[0, 100, 100], [1, 100, 100]],
        >>>     'l': [[0, 100, 100], [1, 20, 20]],
        >>> }
        >>> cmap = pplt.PerceptuallyUniformColormap(data)

        See also
        --------
        LinearSegmentedColormap
        proplot.constructor.Colormap
        """
        # Checks
        data = _pop_props(segmentdata, 'hsla')
        if segmentdata:
            raise ValueError(f'Invalid segmentdata keys {tuple(segmentdata)}.')
        space = _not_none(space, 'hsl').lower()
        if space not in ('rgb', 'hsv', 'hpl', 'hsl', 'hcl'):
            raise ValueError(f'Unknown colorspace {space!r}.')
        # Convert color strings to channel values
        for key, array in data.items():
            if callable(array):  # permit callable
                continue
            for i, xyy in enumerate(array):
                xyy = list(xyy)  # make copy!
                for j, y in enumerate(xyy[1:]):  # modify the y values
                    xyy[j + 1] = _get_channel(y, key, space)
                data[key][i] = xyy
        # Initialize
        super().__init__(name, data, gamma=1.0, N=N, **kwargs)
        # Custom properties
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
            The name of the new colormap. Default is ``self.name + '_copy'``.
        segmentdata, N, alpha, clip, cyclic, gamma, gamma1, gamma2, space : optional
            See `PerceptuallyUniformColormap`. If not provided,
            these are copied from the current colormap.

        See also
        --------
        ListedColormap.copy
        LinearSegmentedColormap.copy
        """
        if name is None:
            name = self.name + '_copy'
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
        cmap = PerceptuallyUniformColormap(
            name, segmentdata, N,
            alpha=alpha, clip=clip, cyclic=cyclic,
            gamma1=gamma1, gamma2=gamma2, space=space
        )
        cmap._rgba_bad = self._rgba_bad
        cmap._rgba_under = self._rgba_under
        cmap._rgba_over = self._rgba_over
        return cmap

    def to_linear_segmented(self, name=None, **kwargs):
        """
        Convert the `PerceptuallyUniformColormap` to a standard
        `LinearSegmentedColormap`. This is used to merge such colormaps.

        Parameters
        ----------
        name : str, optional
            The name of the new colormap. Default is ``self.name + '_copy'``.

        Other parameters
        ----------------
        **kwargs
            Passed to `LinearSegmentedColormap`.

        See also
        --------
        ListedColormap.to_listed
        """
        if not self._isinit:
            self._init()
        if name is None:
            name = self.name + '_copy'
        return LinearSegmentedColormap.from_list(name, self._lut[:-3, :], **kwargs)

    @classmethod
    @warnings._rename_kwargs('0.7', shade='luminance')
    def from_color(cls, name, color, *, space='hsl', **kwargs):
        """
        Return a simple monochromatic "sequential" colormap that blends from white
        or near-white to the input color.

        Parameters
        ----------
        name : str, optional
            The colormap name.
        color : color-spec
            RGB tuple, hex string, or named color string.
        l, s, a, c
            Shorthands for `luminance`, `saturation`, `alpha`, and `chroma`.
        luminance : float or channel-spec, optional
            If float, this is the luminance channel strength on the left-hand
            side of the colormap (default is ``100``). If RGB[A] tuple, hex string,
            or named color string, the luminance is inferred from the color.
        saturation, alpha : float or channel-spec, optional
            As with `luminance`, except the default `saturation` and the default
            `alpha` are the channel values taken from `color`.
        chroma
            Alias for `saturation`.
        space : {'hsl', 'hpl', 'hcl'}, optional
            The colorspace in which luminance and/or saturation are varied.

        Other parameters
        ----------------
        **kwargs
            Passed to `PerceptuallyUniformColormap.from_hsl`.

        Returns
        -------
        `PerceptuallyUniformColormap`
            The colormap.

        See also
        --------
        PerceptuallyUniformColormap.from_hsl
        PerceptuallyUniformColormap.from_list
        """
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
    def from_hsl(cls, name, *, ratios=None, **kwargs):
        """
        Make a `~PerceptuallyUniformColormap` by specifying the hue,
        saturation, and luminance transitions individually.

        Parameters
        ----------
        name : str
            The colormap name.
        h, s, l, a, c
            Shorthands for `hue`, `saturation`, `luminance`, `alpha`, and `chroma`.
        hue : float, color-spec, or list thereof, optional
            Hue channel value or list of values. The shorthand keyword `h`
            is also acceptable. Values can be any of the following.

            1. Numbers, within the range 0 to 360 for hue and 0 to 100 for
               saturation and luminance.
            2. Color string names or hex strings, in which case the channel
               value for that color is looked up.

            If scalar, the hue does not change across the colormap.
            Default is ``0`` (i.e., red).
        saturation, luminance, alpha : float, color-spec, or list thereof, optional
            As with `hue`, but for the saturation, luminance, and alpha (opacity)
            channels, respectively. The default `saturation` is ``50``, luminance is
            ``(100, 20)``, and alpha is ``1`` (i.e., no transparency).
        chroma
            Alias for `saturation`.
        ratios : list of float, optional
            Relative extents of each color transition. Must have length
            ``len(colors) - 1``. Larger numbers indicate a slower
            transition, smaller numbers indicate a faster transition.

            For example, ``luminance=(100, 50, 0)`` with ``ratios=(2, 1)``
            results in a colormap with the transition from luminance ``100``
            to ``50`` taking *twice as long* as the transition from luminance
            ``50`` to ``0``.
        space : {'hsl', 'hpl', 'hcl'}, optional
            The colorspace in which hue, luminance, and/or saturation are varied.

        Other parameters
        ----------------
        **kwargs
            Passed to `PerceptuallyUniformColormap`.

        Returns
        -------
        `PerceptuallyUniformColormap`
            The colormap.

        See also
        --------
        PerceptuallyUniformColormap.from_color
        PerceptuallyUniformColormap.from_list
        """
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
        return cls(name, cdict, **kwargs)

    @classmethod
    def from_list(cls, name, colors, ratios=None, **kwargs):
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

            For example, ``red=(1, 0.5, 0)`` with ``ratios=(2, 1)``
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

        See also
        --------
        matplotlib.colors.LinearSegmentedColormap.from_list
        LinearSegmentedColormap.from_list
        PerceptuallyUniformColormap.from_color
        PerceptuallyUniformColormap.from_hsl
        """
        # Get coordinates
        coords = None
        space = kwargs.get('space', 'hsl')  # use the builtin default
        if not np.iterable(colors):
            raise ValueError(f'Colors must be iterable, got colors={colors!r}')
        if (
            np.iterable(colors[0]) and len(colors[0]) == 2
            and not isinstance(colors[0], str)
        ):
            coords, colors = zip(*colors)
        colors = [to_xyza(color, space) for color in colors]

        # Build segmentdata
        keys = ('hue', 'saturation', 'luminance', 'alpha')
        cdict = {}
        for key, values in zip(keys, zip(*colors)):
            cdict[key] = _make_segment_data(values, coords, ratios)
        return cls(name, cdict, **kwargs)


def _check_levels(levels, allow_descending=True):
    """
    Ensure the levels are monotonic. If they are descending, either
    reverse them or raise an error.
    """
    levels = np.atleast_1d(levels)
    if levels.ndim != 1 or levels.size < 2:
        raise ValueError(f'Levels {levels} must be 1d vector with size >= 2.')
    if isinstance(levels, ma.core.MaskedArray):
        levels = levels.filled(np.nan)
    if not np.all(np.isfinite(levels)):
        raise ValueError(f'Levels {levels} contain invalid values.')
    diffs = np.sign(np.diff(levels))
    if all(diffs == 1):
        descending = False
    elif all(diffs == -1) and allow_descending:
        levels = levels[::-1]
        descending = True
    elif allow_descending:
        raise ValueError(f'Levels {levels} must be monotonic.')
    else:
        raise ValueError(f'Levels {levels} must be monotonically increasing.')
    return levels, descending


def _interpolate_basic(x, x0, x1, y0, y1):
    """
    Basic interpolation between pairs of fixed points.
    """
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def _interpolate_extrapolate(xq, x, y):
    """
    Efficient vectorized linear interpolation. Similar to `numpy.interp`
    except this does not truncate out-of-bounds values (i.e. is reversible).
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
    @warnings._rename_kwargs('0.7', extend='unique')
    def __init__(
        self, levels, norm=None, cmap=None,
        unique=None, step=None, clip=False, descending=False,
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
        cmap : `matplotlib.colors.Colormap`, optional
            The colormap associated with this normalizer. This is used to
            apply default `unique` and `step` settings depending on whether the
            colormap is cyclic and whether distinct "extreme" colors have been
            designated with `~matplotlib.colors.Colormap.set_under` and/or
            `~matplotlib.colors.Colormap.set_over`.
        unique : {'neither', 'both', 'min', 'max'}, optional
            Which out-of-bounds regions should be assigned unique colormap
            colors. The normalizer needs this information so it can ensure
            the colorbar always spans the full range of colormap colors.
        step : float, optional
            The intensity of the transition to out-of-bounds colors as a
            fraction of the adjacent step between in-bounds colors.
            Default is ``1``.
        clip : bool, optional
            Whether to clip values falling outside of the level bins. This
            only has an effect on lower colors when unique is
            ``'min'`` or ``'both'``, and on upper colors when unique is
            ``'max'`` or ``'both'``.
        descending : bool, optional
            Whether the levels are meant to be descending. This will cause
            the colorbar axis to be reversed when it is drawn with a
            `~matplotlib.cm.ScalarMappable` that uses this normalizer.

        Note
        ----
        This normalizer also makes sure that levels always span the full range
        of colors in the colormap, whether `extend` is set to ``'min'``, ``'max'``,
        ``'neither'``, or ``'both'``. By default, when `extend` is not ``'both'``,
        matplotlib cuts off the most intense colors (reserved for "out of bounds"
        data), even though they are not being used. Note that this means
        using a diverging colormap with ``extend='max'`` or ``extend='min'``
        will shift the central color. But that is very strange usage anyway...
        so please just don't do that :)

        See also
        --------
        proplot.axes.apply_cmap
        proplot.constructor.Norm
        """
        # Special properties specific to colormap type
        if cmap is not None:
            over = cmap._rgba_over
            under = cmap._rgba_under
            cyclic = getattr(cmap, '_cyclic', None)
            if cyclic:
                # *Scale bins* as if extend='both' to omit end colors
                step = 0.5
                unique = 'both'

            else:
                # *Scale bins* as if unique='neither' because there may be *discrete
                # change* between minimum color and out-of-bounds colors.
                if over is not None and under is not None:
                    unique = 'both'
                elif over is not None:
                    # Turn off unique bin for over-bounds colors
                    if unique == 'both':
                        unique = 'min'
                    elif unique == 'max':
                        unique = 'neither'
                elif under is not None:
                    # Turn off unique bin for under-bounds colors
                    if unique == 'both':
                        unique = 'min'
                    elif unique == 'max':
                        unique = 'neither'

        # Validate input arguments
        # NOTE: This must be a subclass BoundaryNorm, so ColorbarBase will
        # detect it... even though we completely override it.
        if not norm:
            norm = mcolors.Normalize()
        elif isinstance(norm, mcolors.BoundaryNorm):
            raise ValueError('Normalizer cannot be instance of BoundaryNorm.')
        elif not isinstance(norm, mcolors.Normalize):
            raise ValueError('Normalizer must be instance of Normalize.')
        if unique is None:
            unique = 'neither'
        uniques = ('both', 'min', 'max', 'neither')
        if unique not in uniques:
            raise ValueError(
                f'Unknown unique option {unique!r}. Options are: '
                + ', '.join(map(repr, uniques)) + '.'
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
        # NOTE: Ensure out-of-bounds bin values correspond to out-of-bounds
        # coordinate in case user used set_under or set_over to apply discrete
        # change. See lukelbd/proplot#190
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
        if step is None:
            step = 1.0
        if unique in ('min', 'both'):
            mids[0] += step * (mids[1] - mids[2])
        if unique in ('max', 'both'):
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
        eps = 1e-10  # mids and dest are numpy.float64
        dest = norm(mids)
        dest[0] -= eps
        dest[-1] += eps

        # Attributes
        # NOTE: If clip is True, we clip values to the centers of the end
        # bins rather than vmin/vmax to prevent out-of-bounds colors from
        # getting an in-bounds bin color due to landing on a bin edge.
        # NOTE: With unique='min' the minimimum in-bounds and out-of-bounds
        # colors are the same so clip=True will have no effect. Same goes
        # for unique='max' with maximum colors.
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
    Normalizer that scales data linearly with respect to average *position*
    in an arbitrary monotonically increasing level lists. This is the same
    algorithm used by `~matplotlib.colors.LinearSegmentedColormap` to
    select colors in-between indices in the segment data tables.
    This is the default normalizer paired with `DiscreteNorm` whenever `levels`
    are non-linearly spaced. Can be explicitly used by passing ``norm='segmented'``
    to any command accepting ``cmap``.
    """
    def __init__(self, levels, vmin=None, vmax=None, clip=False):
        """
        Parameters
        ----------
        levels : list of float
            The level boundaries.
        vmin, vmax : None
            Ignored. `vmin` and `vmax` are set to the minimum and
            maximum of `levels`.
        clip : bool, optional
            Whether to clip values falling outside of the minimum and
            maximum levels.

        Example
        -------
        In the below example, unevenly spaced levels are passed to
        `~matplotlib.axes.Axes.contourf`, resulting in the automatic
        application of `LinearSegmentedNorm`.

        >>> import proplot as pplt
        >>> import numpy as np
        >>> levels = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
        >>> data = 10 ** (3 * np.random.rand(10, 10))
        >>> fig, ax = pplt.subplots()
        >>> ax.contourf(data, levels=levels)
        """
        levels = np.asarray(levels)
        levels, _ = _check_levels(levels, allow_descending=False)
        dest = np.linspace(0, 1, len(levels))
        vmin, vmax = np.min(levels), np.max(levels)
        super().__init__(vmin=vmin, vmax=vmax, clip=clip)
        self._x = levels
        self._y = dest

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
        if clip is None:  # builtin clipping
            clip = self.clip
        if clip:  # numpy.clip can handle masked arrays
            value = np.clip(value, self.vmin, self.vmax)
        xq, is_scalar = self.process_value(value)
        yq = _interpolate_extrapolate(xq, self._x, self._y)
        if is_scalar:
            yq = np.atleast_1d(yq)[0]
        return yq

    def inverse(self, value):
        """
        Inverse operation of `~LinearSegmentedNorm.__call__`.

        Parameters
        ----------
        value : numeric
            The data to be un-normalized.
        """
        yq, is_scalar = self.process_value(value)
        xq = _interpolate_extrapolate(yq, self._y, self._x)
        if is_scalar:
            xq = np.atleast_1d(xq)[0]
        return xq


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

        See also
        --------
        proplot.axes.apply_cmap
        proplot.constructor.Norm
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

        * For a smooth colormap, usage is e.g. ``color=('Blues', 0.8)``.
          The number is the colormap index, and must be between 0 and 1.
        * For a color cycle, usage is e.g. ``color=('colorblind', 2)``.
          The number is the list index.

        These examples work with any matplotlib command that accepts a `color`
        keyword arg.
        """
        return self._cache


class _ColorCache(dict):
    """
    Replacement for the native color cache.
    """
    # Matplotlib 'color' args are passed to to_rgba, which tries to read
    # directly from cache and if that fails, sanitizes input, which raises
    # error on receiving (colormap, idx) tuple. So we *have* to override cache.
    def __getitem__(self, key):
        # Split into RGB tuple and opacity
        # Matplotlib caches colors this way internally
        rgb, alpha = key
        if (
            not isinstance(rgb, str) and np.iterable(rgb) and len(rgb) == 2
            and isinstance(rgb[0], str) and isinstance(rgb[1], Number)
        ):
            # Try to get the colormap
            try:
                cmap = _cmap_database[rgb[0]]
            except (KeyError, TypeError):
                pass
            # Read the colormap value
            else:
                if isinstance(cmap, ListedColormap):
                    if not 0 <= rgb[1] < len(cmap.colors):
                        raise ValueError(
                            f'Color cycle sample for {rgb[0]!r} cycle must be '
                            f'between 0 and {len(cmap.colors) - 1}, got {rgb[1]}.'
                        )
                    rgb = cmap.colors[rgb[1]]  # draw from list of colors
                else:
                    if not 0 <= rgb[1] <= 1:
                        raise ValueError(
                            f'Colormap sample for {rgb[0]!r} colormap must be '
                            f'between 0 and 1, got {rgb[1]}.'
                        )
                    rgb = cmap(rgb[1])  # get color selection
                return mcolors.to_rgba(rgb, alpha)

        # Proceed as usual
        return super().__getitem__((rgb, alpha))


def _get_cmap(name=None, lut=None):
    """
    Return the registered colormap instance.

    Parameters
    ----------
    name : `matplotlib.colors.Colormap` or str or None, optional
        If a `~matplotlib.colors.Colormap` instance, it will be returned. Otherwise,
        the name of the registered colormap will be looked up and resampled by `lut`.
        If ``None``, the default colormap :rc:`image.cmap` is returned.
    lut : int or None, optional
        If `name` is not already a `~matplotlib.colors.Colormap` instance
        and `lut` is not None, the colormap will be resampled to have `lut`
        entries in the lookup table.
    """
    # Monkey patch for matplotlib `~matplotlib.get_cmap`. Permits case-insensitive
    # search of monkey-patched colormap database (which was broken in v3.2.0).
    if name is None:
        name = rcParams['image.cmap']
    if isinstance(name, mcolors.Colormap):
        return name
    cmap = _cmap_database[name]
    if lut is not None:
        cmap = cmap._resample(lut)
    return cmap


def _to_proplot_colormap(cmap):
    """
    Translate the input argument to a ProPlot colormap subclass.
    """
    cmap_orig = cmap
    if isinstance(cmap, (ListedColormap, LinearSegmentedColormap)):
        pass
    elif isinstance(cmap, mcolors.LinearSegmentedColormap):
        cmap = LinearSegmentedColormap(
            cmap.name, cmap._segmentdata, cmap.N, cmap._gamma
        )
    elif isinstance(cmap, mcolors.ListedColormap):
        cmap = ListedColormap(
            cmap.colors, cmap.name, cmap.N
        )
    elif isinstance(cmap, mcolors.Colormap):  # base class
        pass
    else:
        raise ValueError(
            f'Invalid colormap type {type(cmap).__name__!r}. '
            'Must be instance of matplotlib.colors.Colormap.'
        )
    cmap._rgba_bad = cmap_orig._rgba_bad
    cmap._rgba_under = cmap_orig._rgba_under
    cmap._rgba_over = cmap_orig._rgba_over
    return cmap


class ColormapDatabase(dict):
    """
    Dictionary subclass used to replace the matplotlib
    colormap registry. See `~ColormapDatabase.__getitem__` and
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
        # Deprecate previous names with support for '_r' and '_s' suffixes
        # NOTE: Must search only for case sensitive *capitalized* names or we would
        # helpfully "redirect" user to correct cmap when they are trying to generate
        # a monochromatic cmap in Colormap and would disallow some color names.
        test = re.sub(r'(_r(_s)?|_s)?\Z', '', key, flags=re.IGNORECASE)
        if not super().__contains__(test):
            if test in _cmaps_removed:
                version = _cmaps_removed[test]
                raise ValueError(
                    f'ProPlot colormap {key!r} was removed in version {version}.'
                )
            if test in _cmaps_renamed:
                test_new, version = _cmaps_renamed[test]
                warnings._warn_proplot(
                    f'Colormap {test!r} was renamed in version {version} and will be '
                    f'removed in the next major release. Please use {test_new!r} instead.'  # noqa: E501
                )
                key = re.sub(test, test_new, key, flags=re.IGNORECASE)

        # Sanitize key and handle suffixes
        key = self._sanitize_key(key, mirror=True)
        shift = key[-2:] == '_s'
        if shift:
            key = key[:-2]
        reverse = key[-2:] == '_r'
        if reverse:
            key = key[:-2]

        # Get item and issue nice error message
        try:
            value = super().__getitem__(key)  # may raise keyerror
        except KeyError:
            raise KeyError(
                f'Invalid colormap or cycle name {key!r}. Options are: '
                + ', '.join(map(repr, self)) + '.'
            )

        # Auto-reverse and auto-shift
        if reverse:
            if hasattr(value, 'reversed'):
                value = value.reversed()
            else:
                raise KeyError(
                    f'Item of type {type(value).__name__!r} '
                    'does not have reversed() method.'
                )
        if shift:
            if hasattr(value, 'shifted'):
                value = value.shifted(180)
            else:
                raise KeyError(
                    f'Item of type {type(value).__name__!r} '
                    'does not have shifted() method.'
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
        key = self._sanitize_key(key, mirror=False)
        item = _to_proplot_colormap(item)
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
        key = re.sub(r'\A(grays)(_r(_s)?|_s)?\Z', r'greys\2', key)
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
# WARNING: Skip over the matplotlib native duplicate entries with
# suffixes '_r' and '_shifted'.
_cmap_database_attr = '_cmap_registry' if hasattr(mcm, '_cmap_registry') else 'cmap_d'
_cmap_database = getattr(mcm, _cmap_database_attr)
if mcm.get_cmap is not _get_cmap:
    mcm.get_cmap = _get_cmap
if not isinstance(_cmap_database, ColormapDatabase):
    _cmap_database = {
        key: value for key, value in _cmap_database.items()
        if key[-2:] != '_r' and key[-8:] != '_shifted'
    }
    _cmap_database = ColormapDatabase(_cmap_database)
    setattr(mcm, _cmap_database_attr, _cmap_database)
