#!/usr/bin/env python3
"""
Various tools that may be useful while making plots.
"""
import re
from numbers import Integral, Number

import matplotlib.colors as mcolors
import matplotlib.font_manager as mfonts
import numpy as np
from matplotlib import rcParams

from .externals import hsluv
from .internals import ic  # noqa: F401
from .internals import docstring, warnings

__all__ = [
    'arange',
    'edges',
    'edges2d',
    'units',
    'set_hue',
    'set_saturation',
    'set_luminance',
    'set_alpha',
    'shift_hue',
    'scale_saturation',
    'scale_luminance',
    'to_hex',
    'to_rgb',
    'to_xyz',
    'to_rgba',
    'to_xyza',
    'shade',  # deprecated
    'saturate',  # deprecated
]

NUMBER = re.compile(r'\A([-+]?[0-9._]+(?:[eE][-+]?[0-9_]+)?)(.*)\Z')
UNIT_DICT = {
    'in': 1.0,
    'ft': 12.0,
    'yd': 36.0,
    'm': 39.37,
    'dm': 3.937,
    'cm': 0.3937,
    'mm': 0.03937,
    'pt': 1 / 72.0,
    'pc': 1 / 6.0,
    'ly': 3.725e+17,
}

# Shared parameters
docstring.snippets['param.rgba'] = """
color : color-spec
    The color. Sanitized with `to_rgba`.
"""
docstring.snippets['param.to_rgb'] = """
color : str, 3-tuple, or 4-tuple
    The color specification. Can be a tuple of channel values, a hex string,
    a registered color name, a cycle color like ``'C0'``, or a colormap color
    (see `~proplot.colors.ColorDatabase`).

    If `space` is ``'rgb'``, this is a tuple of RGB values, and if any
    channels are larger than ``2``, the channels are assumed to be on
    the ``0`` to ``255`` scale and are divided by ``255``.
space : {'rgb', 'hsv', 'hcl', 'hpl', 'hsl'}, optional
    The colorspace for the input channel values. Ignored unless `color` is
    a tuple of numbers.
cycle : str or list, optional
    The registered color cycle name used to interpret colors that
    look like ``'C0'``, ``'C1'``, etc. Default is :rc:`cycle`.
"""
docstring.snippets['param.space'] = """
space : {'hcl', 'hpl', 'hsl', 'hsv'}, optional
    The hue-saturation-luminance-like colorspace used to transform the color.
    Default is the perceptually uniform colorspace ``'hcl'``.
"""

# Shared return values
docstring.snippets['return.hex'] = """
color : str
    A HEX string.
"""


def arange(min_, *args):
    """
    Identical to `numpy.arange` but with inclusive endpoints. For
    example, ``pplt.arange(2, 4)`` returns ``np.array([2, 3, 4])`` instead
    of ``np.array([2, 3])``. This command is useful for generating lists of
    tick locations or colorbar level boundaries.

    See also
    --------
    proplot.axes.CartesianAxes.format
    proplot.constructor.Locator
    proplot.axes.apply_cmap
    """
    # Optional arguments just like np.arange
    if len(args) == 0:
        max_ = min_
        min_ = 0
        step = 1
    elif len(args) == 1:
        max_ = args[0]
        step = 1
    elif len(args) == 2:
        max_ = args[0]
        step = args[1]
    else:
        raise ValueError('Function takes from one to three arguments.')
    # All input is integer
    if all(isinstance(val, Integral) for val in (min_, max_, step)):
        min_, max_, step = np.int64(min_), np.int64(max_), np.int64(step)
        max_ += np.sign(step) * 1
    # Input is float or mixed, cast to float64
    # Don't use np.nextafter with np.finfo(np.dtype(np.float64)).max, because
    # round-off errors from continually adding step to min mess this up
    else:
        min_, max_, step = np.float64(min_), np.float64(max_), np.float64(step)
        max_ += np.sign(step) * (step / 2)
    return np.arange(min_, max_, step)


def edges(z, axis=-1):
    """
    Calculate the approximate "edge" values along an axis given "center"
    values. This is used internally to calculate graticule edges when
    you supply centers to `~matplotlib.axes.Axes.pcolor` or
    `~matplotlib.axes.Axes.pcolormesh`. It is also used to calculate colormap
    level boundaries when you supply centers to plotting methods wrapped by
    `~proplot.axes.apply_cmap`.

    Parameters
    ----------
    z : array-like
        Array of any shape or size.
    axis : int, optional
        The axis along which "edges" are calculated. The size of this axis
        will be increased by one.

    Returns
    -------
    `~numpy.ndarray`
        Array of "edge" coordinates.

    See also
    --------
    edges2d
    proplot.axes.standardize_2d
    proplot.axes.apply_cmap
    """
    z = np.asarray(z)
    z = np.swapaxes(z, axis, -1)
    *nextra, nx = z.shape
    zb = np.zeros((*nextra, nx + 1))

    # Inner edges
    zb[..., 1:-1] = 0.5 * (z[..., :-1] + z[..., 1:])

    # Left, right edges
    zb[..., 0] = 1.5 * z[..., 0] - 0.5 * z[..., 1]
    zb[..., -1] = 1.5 * z[..., -1] - 0.5 * z[..., -2]

    return np.swapaxes(zb, axis, -1)


def edges2d(z):
    """
    Calculate the approximate "edge" values given a 2d grid of "center"
    values. The size of both axes are increased by one. This is used
    internally to calculate graticule edges when you supply centers to
    `~matplotlib.axes.Axes.pcolor` or `~matplotlib.axes.Axes.pcolormesh`.

    Parameters
    ----------
    z : array-like
        A 2d array.

    Returns
    -------
    `~numpy.ndarray`
        Array of "edge" coordinates.

    See also
    --------
    edges
    proplot.axes.standardize_2d
    """
    z = np.asarray(z)
    if z.ndim != 2:
        raise ValueError(f'Input must be a 2d array, but got {z.ndim}d.')
    ny, nx = z.shape
    zb = np.zeros((ny + 1, nx + 1))

    # Inner edges
    zb[1:-1, 1:-1] = 0.25 * (
        z[1:, 1:] + z[:-1, 1:] + z[1:, :-1] + z[:-1, :-1]
    )

    # Left, right, top, bottom edges
    zb[0, :] += edges(1.5 * z[0, :] - 0.5 * z[1, :])
    zb[-1, :] += edges(1.5 * z[-1, :] - 0.5 * z[-2, :])
    zb[:, 0] += edges(1.5 * z[:, 0] - 0.5 * z[:, 1])
    zb[:, -1] += edges(1.5 * z[:, -1] - 0.5 * z[:, -2])
    zb[[0, 0, -1, -1], [0, -1, -1, 0]] *= 0.5  # corner correction

    return zb


def _transform_color(func, color, space):
    """
    Standardize input for color transformation functions.
    """
    *color, opacity = to_rgba(color)
    channels = list(to_xyz(color, space=space))
    channels = func(channels)  # apply transform
    color = to_rgb(channels, space=space)
    color = tuple(np.clip(color, 0, 1))  # clip to valid range
    return mcolors.to_hex((*color, opacity))


@docstring.add_snippets
def shift_hue(color, shift=0, space='hcl'):
    """
    Shift the hue channel of a color.

    Parameters
    ----------
    %(param.rgba)s
    shift : float, optoinal
        The HCL hue channel is offset by this value.
    %(param.space)s

    Returns
    -------
    %(return.hex)s

    See also
    --------
    set_hue
    set_saturation
    set_luminance
    set_alpha
    scale_saturation
    scale_luminance
    """
    def func(channels):
        channels[0] += shift
        channels[0] %= 360
        return channels

    return _transform_color(func, color, space)


@docstring.add_snippets
def scale_saturation(color, scale=1, space='hcl'):
    """
    Scale the saturation channel of a color.

    Parameters
    ----------
    %(param.rgba)s
    scale : float, optoinal
        The HCL saturation channel is multiplied by this value.
    %(param.space)s

    Returns
    -------
    %(return.hex)s

    See also
    --------
    set_hue
    set_saturation
    set_luminance
    set_alpha
    shift_hue
    scale_luminance
    """
    def func(channels):
        channels[1] *= scale
        return channels

    return _transform_color(func, color, space)


@docstring.add_snippets
def scale_luminance(color, scale=1, space='hcl'):
    """
    Scale the luminance channel of a color.

    Parameters
    ----------
    %(param.rgba)s
    scale : float, optoinal
        The luminance channel is multiplied by this value.
    %(param.space)s

    Returns
    -------
    %(return.hex)s

    See also
    --------
    set_hue
    set_saturation
    set_luminance
    set_alpha
    shift_hue
    scale_saturation
    """
    def func(channels):
        channels[2] *= scale
        return channels

    return _transform_color(func, color, space)


@docstring.add_snippets
def set_hue(color, hue, space='hcl'):
    """
    Return a color with a different hue and the same luminance and saturation
    as the input color.

    Parameters
    ----------
    %(param.rgba)s
    hue : float, optional
        The new hue. Should lie between ``0`` and ``360`` degrees.
    %(param.space)s

    Returns
    -------
    %(return.hex)s

    See also
    --------
    set_saturation
    set_luminance
    set_alpha
    shift_hue
    scale_saturation
    scale_luminance
    """
    def func(channels):
        channels[0] = hue
        return channels

    return _transform_color(func, color, space)


@docstring.add_snippets
def set_saturation(color, saturation, space='hcl'):
    """
    Return a color with a different saturation and the same hue and luminance
    as the input color.

    Parameters
    ----------
    %(param.rgba)s
    saturation : float, optional
        The new saturation. Should lie between ``0`` and ``360`` degrees.
    %(param.space)s

    Returns
    -------
    %(return.hex)s

    See also
    --------
    set_hue
    set_luminance
    set_alpha
    shift_hue
    scale_saturation
    scale_luminance
    """
    def func(channels):
        channels[1] = saturation
        return channels

    return _transform_color(func, color, space)


@docstring.add_snippets
def set_luminance(color, luminance, space='hcl'):
    """
    Return a color with a different luminance and the same hue and saturation
    as the input color.

    Parameters
    ----------
    %(param.rgba)s
    luminance : float, optional
        The new luminance. Should lie between ``0`` and ``100``.
    %(param.space)s

    Returns
    -------
    %(return.hex)s

    See also
    --------
    set_hue
    set_saturation
    set_alpha
    shift_hue
    scale_saturation
    scale_luminance
    """
    def func(channels):
        channels[2] = luminance
        return channels

    return _transform_color(func, color, space)


@docstring.add_snippets
def set_alpha(color, alpha):
    """
    Return a color with the opacity channel set to the specified value.

    Parameters
    ----------
    %(param.rgba)s
    alpha : float, optional
        The new opacity. Should be between ``0`` and ``1``.

    Returns
    -------
    %(return.hex)s

    See also
    --------
    set_hue
    set_saturation
    set_luminance
    shift_hue
    scale_saturation
    scale_luminance
    """
    color = list(to_rgba(color))
    color[3] = alpha
    return to_hex(color)


def to_hex(color, space='rgb', cycle=None, keep_alpha=True):
    """
    Translate the color in *any* format and from *any* colorspace
    to a HEX string. This is a generalization of `matplotlib.colors.to_hex`.

    Parameters
    ----------
    %(param.to_rgb)s
    keep_alpha : bool, optional
        Whether to keep the opacity channel. If ``True`` an 8-digit HEX
        is returned. Otherwise a 6-digit HEX is returned. Default is ``True``.

    Returns
    -------
    %(return.hex)s

    See also
    --------
    to_rgb
    to_rgba
    to_xyz
    to_xyza
    """
    rgba = to_rgba(color, space=space, cycle=cycle)
    return mcolors.to_hex(rgba, keep_alpha=keep_alpha)


@docstring.add_snippets
def to_rgb(color, space='rgb', cycle=None):
    """
    Translate the color in *any* format and from *any* colorspace to an RGB
    tuple. This is a generalization of `matplotlib.colors.to_rgb` and the
    inverse of `to_xyz`.

    Parameters
    ----------
    %(param.to_rgb)s

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
    return to_rgba(color, space=space, cycle=cycle)[:3]


@docstring.add_snippets
def to_rgba(color, space='rgb', cycle=None):
    """
    Translate the color in *any* format and from *any* colorspace to an RGB
    tuple. This is a generalization of `matplotlib.colors.to_rgba` and the
    inverse of `to_xyza`.

    Parameters
    ----------
    %(param.to_rgb)s

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
    # Convert color cycle strings
    if isinstance(color, str) and re.match(r'\AC[0-9]\Z', color):
        if isinstance(cycle, str):
            from .colors import _cmap_database
            try:
                cycle = _cmap_database[cycle].colors
            except (KeyError, AttributeError):
                cycles = sorted(
                    name for name, cmap in _cmap_database.items()
                    if isinstance(cmap, mcolors.ListedColormap)
                )
                raise ValueError(
                    f'Invalid cycle {cycle!r}. Options are: '
                    + ', '.join(map(repr, cycles)) + '.'
                )
        elif cycle is None:
            cycle = rcParams['axes.prop_cycle'].by_key()
            if 'color' not in cycle:
                cycle = ['k']
            else:
                cycle = cycle['color']
        else:
            raise ValueError(f'Invalid cycle {cycle!r}.')
        color = cycle[int(color[-1]) % len(cycle)]

    # Translate RGB strings and (colormap, index) tuples
    opacity = 1
    if isinstance(color, str) or np.iterable(color) and len(color) == 2:
        try:
            *color, opacity = mcolors.to_rgba(color)  # ensure is valid color
        except (ValueError, TypeError):
            raise ValueError(f'Invalid RGB argument {color!r}.')

    # Pull out alpha channel
    if len(color) == 4:
        *color, opacity = color
    elif len(color) != 3:
        raise ValueError(f'Invalid RGB argument {color!r}.')

    # Translate arbitrary colorspaces
    if space == 'rgb':
        try:
            if any(c > 2 for c in color):
                color = [c / 255 for c in color]  # scale to within 0-1
            color = tuple(color)
        except (ValueError, TypeError):
            raise ValueError(f'Invalid RGB argument {color!r}.')
    elif space == 'hsv':
        color = hsluv.hsl_to_rgb(*color)
    elif space == 'hcl':
        color = hsluv.hcl_to_rgb(*color)
    elif space == 'hsl':
        color = hsluv.hsluv_to_rgb(*color)
    elif space == 'hpl':
        color = hsluv.hpluv_to_rgb(*color)
    else:
        raise ValueError('Invalid color {color!r} for colorspace {space!r}.')

    # Return RGB or RGBA
    return (*color, opacity)


@docstring.add_snippets
def to_xyz(color, space='hcl'):
    """
    Translate color in *any* format to a tuple of channel values in *any*
    colorspace. This is the inverse of `to_rgb`.

    Parameters
    ----------
    %(param.rgba)s
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
    return to_xyza(color, space)[:3]


@docstring.add_snippets
def to_xyza(color, space='hcl'):
    """
    Translate color in *any* format to a tuple of channel values in *any*
    colorspace. This is the inverse of `to_rgba`.

    Parameters
    ----------
    %(param.rgba)s
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
    *color, opacity = to_rgba(color)
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


@warnings._rename_kwargs('0.6', units='dest')
def units(value, dest='in', axes=None, figure=None, width=True):
    """
    Convert values and lists of values between arbitrary physical units. This
    is used internally all over ProPlot, permitting flexible units for various
    keyword arguments.

    Parameters
    ----------
    value : float or str or list thereof
        A size specifier or *list thereof*. If numeric, nothing is done.
        If string, it is converted to the units `dest`. The string should look
        like ``'123.456unit'``, where the number is the magnitude and
        ``'unit'`` is one of the following.

        =========  =====================================================
        Key        Description
        =========  =====================================================
        ``'m'``    Meters
        ``'dm'``   Decimeters
        ``'cm'``   Centimeters
        ``'mm'``   Millimeters
        ``'yd'``   Yards
        ``'ft'``   Feet
        ``'in'``   Inches
        ``'pt'``   `Points <pt_>`_ (1/72 inches)
        ``'pc'``   `Pica <pc_>`_ (1/6 inches)
        ``'px'``   Pixels on screen, using dpi of :rcraw:`figure.dpi`
        ``'pp'``   Pixels once printed, using dpi of :rcraw:`savefig.dpi`
        ``'em'``   `Em square <em_>`_ for :rcraw:`font.size`
        ``'en'``   `En square <en_>`_ for :rcraw:`font.size`
        ``'Em'``   `Em square <em_>`_ for :rcraw:`axes.titlesize`
        ``'En'``   `En square <en_>`_ for :rcraw:`axes.titlesize`
        ``'ax'``   Axes-relative units (not always available)
        ``'fig'``  Figure-relative units (not always available)
        ``'ly'``   Light years ;)
        =========  =====================================================

        .. _pt: https://en.wikipedia.org/wiki/Point_(typography)
        .. _pc: https://en.wikipedia.org/wiki/Pica_(typography)
        .. _em: https://en.wikipedia.org/wiki/Em_(typography)
        .. _en: https://en.wikipedia.org/wiki/En_(typography)

    dest : str, optional
        The destination units. Default is inches, i.e. ``'in'``.
    axes : `~matplotlib.axes.Axes`, optional
        The axes to use for scaling units that look like ``'0.1ax'``.
    figure : `~matplotlib.figure.Figure`, optional
        The figure to use for scaling units that look like ``'0.1fig'``. If
        ``None`` we try to get the figure from ``axes.figure``.
    width : bool, optional
        Whether to use the width or height for the axes and figure relative
        coordinates.
    """
    # Font unit scales
    # NOTE: Delay font_manager import, because want to avoid rebuilding font
    # cache, which means import must come after TTFPATH added to environ
    # by register_fonts()!
    fontsize_small = rcParams['font.size']  # must be absolute
    fontsize_large = rcParams['axes.titlesize']
    if isinstance(fontsize_large, str):
        scale = mfonts.font_scalings.get(fontsize_large, 1)
        fontsize_large = fontsize_small * scale

    # Scales for converting physical units to inches
    unit_dict = UNIT_DICT.copy()
    unit_dict.update({
        'em': fontsize_small / 72.0,
        'en': 0.5 * fontsize_small / 72.0,
        'Em': fontsize_large / 72.0,
        'En': 0.5 * fontsize_large / 72.0,
    })

    # Scales for converting display units to inches
    # WARNING: In ipython shell these take the value 'figure'
    if not isinstance(rcParams['figure.dpi'], str):
        # once generated by backend
        unit_dict['px'] = 1 / rcParams['figure.dpi']
    if not isinstance(rcParams['savefig.dpi'], str):
        # once 'printed' i.e. saved
        unit_dict['pp'] = 1 / rcParams['savefig.dpi']

    # Scales relative to axes and figure objects
    if axes is not None and hasattr(axes, 'get_size_inches'):  # proplot axes
        unit_dict['ax'] = axes.get_size_inches()[1 - int(width)]
    if figure is None:
        figure = getattr(axes, 'figure', None)
    if figure is not None and hasattr(
            figure, 'get_size_inches'):  # proplot axes
        unit_dict['fig'] = figure.get_size_inches()[1 - int(width)]

    # Scale for converting inches to arbitrary other unit
    try:
        scale = unit_dict[dest]
    except KeyError:
        raise ValueError(
            f'Invalid destination units {dest!r}. Valid units are '
            + ', '.join(map(repr, unit_dict.keys())) + '.'
        )

    # Convert units for each value in list
    result = []
    singleton = not np.iterable(value) or isinstance(value, str)
    for val in ((value,) if singleton else value):
        if val is None or isinstance(val, Number):
            result.append(val)
            continue
        elif not isinstance(val, str):
            raise ValueError(
                f'Size spec must be string or number or list thereof. '
                f'Got {value!r}.'
            )
        regex = NUMBER.match(val)
        if not regex:
            raise ValueError(
                f'Invalid size spec {val!r}. Valid units are '
                + ', '.join(map(repr, unit_dict.keys())) + '.'
            )
        number, units = regex.groups()  # second group is exponential
        try:
            result.append(
                float(number) * (unit_dict[units] / scale if units else 1)
            )
        except (KeyError, ValueError):
            raise ValueError(
                f'Invalid size spec {val!r}. Valid units are '
                + ', '.join(map(repr, unit_dict.keys())) + '.'
            )
    if singleton:
        result = result[0]
    return result


# Deprecations
shade, saturate = warnings._rename_objs(
    '0.6',
    shade=scale_luminance,
    saturate=scale_saturation,
)
