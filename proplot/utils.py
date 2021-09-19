#!/usr/bin/env python3
"""
Various tools that may be useful while making plots.
"""
# WARNING: Cannot import 'rc' anywhere in this file or we get circular import
# issues. The rc param validators need functions in this file.
import functools
import re
from numbers import Integral, Real

import matplotlib.colors as mcolors
import matplotlib.font_manager as mfonts
import numpy as np
from matplotlib import rcParams as rc_matplotlib

from .externals import hsluv
from .internals import ic  # noqa: F401
from .internals import _not_none, docstring, warnings

__all__ = [
    'arange',
    'edges',
    'edges2d',
    'get_colors',
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
    'units',
    'shade',  # deprecated
    'saturate',  # deprecated
]

UNIT_REGEX = re.compile(
    r'\A([-+]?[0-9._]+(?:[eE][-+]?[0-9_]+)?)(.*)\Z'  # float with trailing units
)
UNIT_DICT = {
    'in': 1.0,
    'ft': 12.0,
    'yd': 36.0,
    'm': 39.37,
    'dm': 3.937,
    'cm': 0.3937,
    'mm': 0.03937,
    'pc': 1 / 6.0,
    'pt': 1 / 72.0,
    'ly': 3.725e17,
}


# Unit docstrings
# NOTE: Try to fit this into a single line. Cannot break up with newline as that will
# mess up docstring indentation since this is placed in indented param lines.
_units_docstring = 'If float, units are {units}. If string, interpreted by `~proplot.utils.units`.'  # noqa: E501
docstring._snippet_manager['units.pt'] = _units_docstring.format(units='points')
docstring._snippet_manager['units.in'] = _units_docstring.format(units='inches')
docstring._snippet_manager['units.em'] = _units_docstring.format(units='em-widths')


# Color docstrings
_docstring_rgba = """
color : color-spec
    The color. Sanitized with `to_rgba`.
"""
_docstring_to_rgb = """
color : color-spec
    The color. Can be a 3-tuple or 4-tuple of channel values, a hex
    string, a registered color name, a cycle color like ``'C0'``, or
    a 2-tuple colormap coordinate specification like ``('magma', 0.5)``
    (see `~proplot.colors.ColorDatabase` for details).

    If `space` is ``'rgb'``, this is a tuple of RGB values, and any
    channels are larger than ``2``, the channels are assumed to be
    on the ``0`` to ``255`` scale and are divided by ``255``.
space : {'rgb', 'hsv', 'hcl', 'hpl', 'hsl'}, optional
    The colorspace for the input channel values. Ignored unless `color` is
    a tuple of numbers.
cycle : str, optional
    The registered color cycle name used to interpret colors that
    look like ``'C0'``, ``'C1'``, etc. Default is :rc:`cycle`.
clip : bool, optional
    Whether to clip channel values into the valid ``0`` to ``1`` range.
    Default is ``True``.
"""
_docstring_space = """
space : {'hcl', 'hpl', 'hsl', 'hsv'}, optional
    The hue-saturation-luminance-like colorspace used to transform the color.
    Default is the perceptually uniform colorspace ``'hcl'``.
"""
_docstring_hex = """
color : str
    An 8-digit HEX string indicating the
    red, green, blue, and alpha channel values.
"""
docstring._snippet_manager['utils.color'] = _docstring_rgba
docstring._snippet_manager['utils.hex'] = _docstring_hex
docstring._snippet_manager['utils.space'] = _docstring_space
docstring._snippet_manager['utils.to'] = _docstring_to_rgb


def _keep_units(func):
    """
    Very simple decorator to strip and re-apply the same units.
    """
    # NOTE: Native UnitRegistry.wraps() is not sufficient since it enforces
    # unit types rather than arbitrary units. This wrapper is similar.
    @functools.wraps(func)
    def _with_stripped_units(data, *args, **kwargs):
        units = 1
        if hasattr(data, 'units') and hasattr(data, 'magnitude'):
            data, units = data.magnitude, data.units
        result = func(data, *args, **kwargs)
        return result * units
    return _with_stripped_units


def arange(min_, *args):
    """
    Identical to `numpy.arange` but with inclusive endpoints. For example,
    ``pplt.arange(2, 4)`` returns the numpy array ``[2, 3, 4]`` instead of
    ``[2, 3]``. This is useful for generating lists of tick locations or
    colormap levels, e.g. ``ax.format(xlocator=pplt.arange(0, 10))``
    or ``ax.pcolor(levels=pplt.arange(0, 10))``.

    Parameters
    ----------
    *args : float
        If three arguments are passed, these are the minimum, maximum, and step
        size. If fewer than three arguments are passed, the step size is ``1``.
        If one argument is passed, this is the maximum, and the minimum is ``0``.

    Returns
    -------
    numpy.ndarray
        Array of points.

    See also
    --------
    numpy.arange
    proplot.constructor.Locator
    proplot.axes.CartesianAxes.format
    proplot.axes.PolarAxes.format
    proplot.axes.GeoAxes.format
    proplot.axes.Axes.colorbar
    proplot.axes.PlotAxes
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
        max_ += np.sign(step)
    # Input is float or mixed, cast to float64
    # Don't use np.nextafter with np.finfo(np.dtype(np.float64)).max, because
    # round-off errors from continually adding step to min mess this up
    else:
        min_, max_, step = np.float64(min_), np.float64(max_), np.float64(step)
        max_ += 0.5 * step
    return np.arange(min_, max_, step)


@_keep_units
def edges(z, axis=-1):
    """
    Calculate the approximate "edge" values along an axis given "center" values. The
    size of the axis is increased by one. This is used internally to calculate graticule
    edges when you supply centers to pseudocolor commands.

    Parameters
    ----------
    z : array-like
        An array of any shape.
    axis : int, optional
        The axis along which "edges" are calculated. The size of this
        axis will be increased by one.

    Returns
    -------
    numpy.ndarray
        Array of "edge" coordinates.

    See also
    --------
    edges2d
    proplot.axes.PlotAxes.pcolor
    proplot.axes.PlotAxes.pcolormesh
    proplot.axes.PlotAxes.pcolorfast
    """
    z = np.asarray(z)
    z = np.swapaxes(z, axis, -1)
    *nextra, nx = z.shape
    zb = np.zeros((*nextra, nx + 1))

    # Inner edges
    zb[..., 1:-1] = 0.5 * (z[..., :-1] + z[..., 1:])

    # Outer edges
    zb[..., 0] = 1.5 * z[..., 0] - 0.5 * z[..., 1]
    zb[..., -1] = 1.5 * z[..., -1] - 0.5 * z[..., -2]

    return np.swapaxes(zb, axis, -1)


@_keep_units
def edges2d(z):
    """
    Calculate the approximate "edge" values given a 2D grid of "center"
    values. The size of both axes is increased by one. This is used
    internally to calculate graticule edges when you supply centers
    to pseudocolor plot commands.

    Parameters
    ----------
    z : array-like
        A 2D array.

    Returns
    -------
    numpy.ndarray
        Array of "edge" coordinates.

    See also
    --------
    edges
    proplot.axes.PlotAxes.pcolor
    proplot.axes.PlotAxes.pcolormesh
    proplot.axes.PlotAxes.pcolorfast
    """
    z = np.asarray(z)
    if z.ndim != 2:
        raise ValueError(f'Input must be a 2D array, but got {z.ndim}D.')
    ny, nx = z.shape
    zb = np.zeros((ny + 1, nx + 1))

    # Inner edges
    zb[1:-1, 1:-1] = 0.25 * (z[1:, 1:] + z[:-1, 1:] + z[1:, :-1] + z[:-1, :-1])

    # Outer edges
    zb[0, :] += edges(1.5 * z[0, :] - 0.5 * z[1, :])
    zb[-1, :] += edges(1.5 * z[-1, :] - 0.5 * z[-2, :])
    zb[:, 0] += edges(1.5 * z[:, 0] - 0.5 * z[:, 1])
    zb[:, -1] += edges(1.5 * z[:, -1] - 0.5 * z[:, -2])
    zb[[0, 0, -1, -1], [0, -1, -1, 0]] *= 0.5  # corner correction

    return zb


def get_colors(*args, **kwargs):
    """
    Get the colors associated with a registered or
    on-the-fly color cycle.

    Parameters
    ----------
    *args, **kwargs
        Passed to `~proplot.constructor.Cycle`.

    Returns
    -------
    colors : list of str
        A list of HEX strings.

    See also
    --------
    cycler.Cycler
    proplot.constructor.Cycle
    proplot.constructor.Colormap
    """
    from .constructor import Cycle  # delayed to avoid cyclic imports
    cycle = Cycle(*args, **kwargs)
    colors = [to_hex(dict_['color']) for dict_ in cycle]
    return colors


def _transform_color(func, color, space):
    """
    Standardize input for color transformation functions.
    """
    *color, opacity = to_rgba(color)
    color = to_xyz(color, space=space)
    color = func(list(color))  # apply transform
    return to_hex((*color, opacity), space=space)


@docstring._snippet_manager
def shift_hue(color, shift=0, space='hcl'):
    """
    Shift the hue channel of a color.

    Parameters
    ----------
    %(utils.color)s
    shift : float, optoinal
        The HCL hue channel is offset by this value.
    %(utils.space)s

    Returns
    -------
    %(utils.hex)s

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


@docstring._snippet_manager
def scale_saturation(color, scale=1, space='hcl'):
    """
    Scale the saturation channel of a color.

    Parameters
    ----------
    %(utils.color)s
    scale : float, optoinal
        The HCL saturation channel is multiplied by this value.
    %(utils.space)s

    Returns
    -------
    %(utils.hex)s

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


@docstring._snippet_manager
def scale_luminance(color, scale=1, space='hcl'):
    """
    Scale the luminance channel of a color.

    Parameters
    ----------
    %(utils.color)s
    scale : float, optoinal
        The luminance channel is multiplied by this value.
    %(utils.space)s

    Returns
    -------
    %(utils.hex)s

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


@docstring._snippet_manager
def set_hue(color, hue, space='hcl'):
    """
    Return a color with a different hue and the same luminance and saturation
    as the input color.

    Parameters
    ----------
    %(utils.color)s
    hue : float, optional
        The new hue. Should lie between ``0`` and ``360`` degrees.
    %(utils.space)s

    Returns
    -------
    %(utils.hex)s

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


@docstring._snippet_manager
def set_saturation(color, saturation, space='hcl'):
    """
    Return a color with a different saturation and the same hue and luminance
    as the input color.

    Parameters
    ----------
    %(utils.color)s
    saturation : float, optional
        The new saturation. Should lie between ``0`` and ``360`` degrees.
    %(utils.space)s

    Returns
    -------
    %(utils.hex)s

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


@docstring._snippet_manager
def set_luminance(color, luminance, space='hcl'):
    """
    Return a color with a different luminance and the same hue and saturation
    as the input color.

    Parameters
    ----------
    %(utils.color)s
    luminance : float, optional
        The new luminance. Should lie between ``0`` and ``100``.
    %(utils.space)s

    Returns
    -------
    %(utils.hex)s

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


@docstring._snippet_manager
def set_alpha(color, alpha):
    """
    Return a color with the opacity channel set to the specified value.

    Parameters
    ----------
    %(utils.color)s
    alpha : float, optional
        The new opacity. Should be between ``0`` and ``1``.

    Returns
    -------
    %(utils.hex)s

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


def _translate_cycle_color(color, cycle=None):
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
        cycle = rc_matplotlib['axes.prop_cycle'].by_key()
        if 'color' not in cycle:
            cycle = ['k']
        else:
            cycle = cycle['color']
    else:
        raise ValueError(f'Invalid cycle {cycle!r}.')

    return cycle[int(color[-1]) % len(cycle)]


@docstring._snippet_manager
def to_hex(color, space='rgb', cycle=None, keep_alpha=True):
    """
    Translate the color from an arbitrary colorspace to a HEX string.
    This is a generalization of `matplotlib.colors.to_hex`.

    Parameters
    ----------
    %(utils.to)s
    keep_alpha : bool, optional
        Whether to keep the opacity channel. If ``True`` an 8-digit HEX
        is returned. Otherwise a 6-digit HEX is returned. Default is ``True``.

    Returns
    -------
    %(utils.hex)s

    See also
    --------
    to_rgb
    to_rgba
    to_xyz
    to_xyza
    """
    rgba = to_rgba(color, space=space, cycle=cycle)
    return mcolors.to_hex(rgba, keep_alpha=keep_alpha)


@docstring._snippet_manager
def to_rgb(color, space='rgb', cycle=None):
    """
    Translate the color from an arbitrary colorspace to an RGB tuple. This is
    a generalization of `matplotlib.colors.to_rgb` and the inverse of `to_xyz`.

    Parameters
    ----------
    %(utils.to)s

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


@docstring._snippet_manager
def to_rgba(color, space='rgb', cycle=None, clip=True):
    """
    Translate the color from an arbitrary colorspace to an RGBA tuple. This is
    a generalization of `matplotlib.colors.to_rgba` and the inverse of `to_xyz`.

    Parameters
    ----------
    %(utils.to)s

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
    if isinstance(color, str) and re.match(r'\AC[0-9]\Z', color):
        color = _translate_cycle_color(color, cycle=cycle)

    # Translate RGB strings and (colormap, index) tuples
    # NOTE: Cannot use is_color_like because might have HSL channel values
    opacity = 1
    if (
        isinstance(color, str)
        or np.iterable(color) and len(color) == 2
    ):
        color = mcolors.to_rgba(color)  # also enforced validity
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
def to_xyz(color, space='hcl'):
    """
    Translate color in *any* format to a tuple of channel values in *any*
    colorspace. This is the inverse of `to_rgb`.

    Parameters
    ----------
    %(utils.color)s
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


@docstring._snippet_manager
def to_xyza(color, space='hcl'):
    """
    Translate color in *any* format to a tuple of channel values in *any*
    colorspace. This is the inverse of `to_rgba`.

    Parameters
    ----------
    %(utils.color)s
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


def _fontsize_to_pt(size):
    """
    Translate font preset size or unit string to points.
    """
    scalings = mfonts.font_scalings
    if not isinstance(size, str):
        return size
    if size in mfonts.font_scalings:
        return rc_matplotlib['font.size'] * scalings[size]
    try:
        return units(size, 'pt')
    except ValueError:
        raise KeyError(
            f'Invalid font size {size!r}. Can be points or one of the preset scalings: '
            + ', '.join(f'{key!r} ({value})' for key, value in scalings.items())
            + '.'
        )


@warnings._rename_kwargs('0.6', units='dest')
def units(
    value, numeric=None, dest=None, *, fontsize=None, figure=None, axes=None, width=None
):
    """
    Convert values between arbitrary physical units. This is used internally all
    over ProPlot, permitting flexible units for various keyword arguments.

    Parameters
    ----------
    value : float or str or sequence
        A size specifier or sequence of size specifiers. If numeric, units are
        converted from `numeric` to `dest`. If string, units are converted to
        `dest` according to the string specifier. The string should look like
        ``'123.456unit'``, where the number is the magnitude and ``'unit'``
        matches a key in the below table.

        .. _units_table:

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
        ``'pc'``   `Pica <pc_>`_ (1/6 inches)
        ``'pt'``   `Points <pt_>`_ (1/72 inches)
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

    numeric : str, optional
        The units associated with numeric input. Default is inches.
    dest : str, optional
        The destination units. Default is the same as `numeric`.
    fontsize : str or float, optional
        The font size in points used for scaling. Default is
        :rcraw:`font.size` for ``em`` and ``en`` units and
        :rcraw:`axes.titlesize` for ``Em`` and ``En`` units.
    axes : `~matplotlib.axes.Axes`, optional
        The axes to use for scaling units that look like ``'0.1ax'``.
    figure : `~matplotlib.figure.Figure`, optional
        The figure to use for scaling units that look like ``'0.1fig'``.
        If ``None`` we try to get the figure from ``axes.figure``.
    width : bool, optional
        Whether to use the width or height for the axes and figure
        relative coordinates.
    """
    # Scales for converting physical units to inches
    fontsize_small = _not_none(fontsize, rc_matplotlib['font.size'])  # always absolute
    fontsize_small = _fontsize_to_pt(fontsize_small)
    fontsize_large = _not_none(fontsize, rc_matplotlib['axes.titlesize'])
    fontsize_large = _fontsize_to_pt(fontsize_large)
    unit_dict = UNIT_DICT.copy()
    unit_dict.update(
        {
            'em': fontsize_small / 72.0,
            'en': 0.5 * fontsize_small / 72.0,
            'Em': fontsize_large / 72.0,
            'En': 0.5 * fontsize_large / 72.0,
        }
    )

    # Scales for converting display units to inches
    # WARNING: In ipython shell these take the value 'figure'
    if not isinstance(rc_matplotlib['figure.dpi'], str):
        unit_dict['px'] = 1 / rc_matplotlib['figure.dpi']  # once generated by backend
    if not isinstance(rc_matplotlib['savefig.dpi'], str):
        unit_dict['pp'] = 1 / rc_matplotlib['savefig.dpi']  # once 'printed' i.e. saved

    # Scales relative to axes and figure objects
    if axes is not None and hasattr(axes, '_get_size_inches'):  # proplot axes
        unit_dict['ax'] = axes._get_size_inches()[1 - int(width)]
    if figure is None:
        figure = getattr(axes, 'figure', None)
    if figure is not None and hasattr(figure, 'get_size_inches'):
        unit_dict['fig'] = figure.get_size_inches()[1 - int(width)]

    # Scale for converting inches to arbitrary other unit
    if numeric is None and dest is None:
        numeric = dest = 'in'
    elif numeric is None:
        numeric = dest
    elif dest is None:
        dest = numeric
    options = 'Valid units are ' + ', '.join(map(repr, unit_dict)) + '.'
    try:
        nscale = unit_dict[numeric]
    except KeyError:
        raise ValueError(f'Invalid numeric units {numeric!r}. ' + options)
    try:
        dscale = unit_dict[dest]
    except KeyError:
        raise ValueError(f'Invalid destination units {dest!r}. ' + options)

    # Convert units for each value in list
    result = []
    singleton = not np.iterable(value) or isinstance(value, str)
    for val in (value,) if singleton else value:
        # Silently pass None
        if val is None:
            result.append(val)
            continue
        # Get unit string
        if isinstance(val, Real):
            number, units = val, None
        elif isinstance(val, str):
            regex = UNIT_REGEX.match(val)
            if regex:
                number, units = regex.groups()  # second group is exponential
            else:
                raise ValueError(f'Invalid unit size spec {val!r}.')
        else:
            raise ValueError(f'Invalid unit size spec {val!r}.')
        # Convert with units
        if not units:
            result.append(float(number) * nscale / dscale)
        elif units in unit_dict:
            result.append(float(number) * unit_dict[units] / dscale)
        else:
            raise ValueError(f'Invalid input units {units!r}. ' + options)
    return result[0] if singleton else result


# Deprecations
shade, saturate = warnings._rename_objs(
    '0.6',
    shade=scale_luminance,
    saturate=scale_saturation,
)
