#!/usr/bin/env python3
"""
Various tools that may be useful while making plots.
"""
# WARNING: Cannot import 'rc' anywhere in this file or we get circular import
# issues. The rc param validators need functions in this file.
import functools
import re
from numbers import Integral, Real

import matplotlib.font_manager as mfonts
import numpy as np
from matplotlib import rcParams as rc_matplotlib

from .internals import ic  # noqa: F401
from .internals import _not_none, warnings

__all__ = [
    'arange',
    'edges',
    'edges2d',
    'units',
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
    Calculate the approximate "edge" values along an axis given "center" values.
    The size of the axis is increased by one. This is used internally to calculate
    coordinate edges when you supply coordinate centers to pseudocolor commands.

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
    *dims, n = z.shape
    zb = np.zeros((*dims, n + 1))

    # Inner edges
    zb[..., 1:-1] = 0.5 * (z[..., :-1] + z[..., 1:])

    # Outer edges
    zb[..., 0] = 1.5 * z[..., 0] - 0.5 * z[..., 1]
    zb[..., -1] = 1.5 * z[..., -1] - 0.5 * z[..., -2]

    return np.swapaxes(zb, axis, -1)


@_keep_units
def edges2d(z):
    """
    Calculate the approximate "edge" values given a 2D grid of "center" values.
    The size of both axes is increased by one. This is used internally to calculate
    coordinate edges when you supply coordinate to pseudocolor commands.

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
    on-the-fly color cycle or colormap.

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
    proplot.constructor.Cycle
    proplot.constructor.Colormap
    """
    from .constructor import Cycle  # delayed to avoid cyclic imports
    cycle = Cycle(*args, **kwargs)
    colors = [to_hex(dict_['color']) for dict_ in cycle]
    return colors


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


@warnings._rename_kwargs('0.6.0', units='dest')
def units(
    value, numeric=None, dest=None, *, fontsize=None, figure=None, axes=None, width=None
):
    """
    Convert values between arbitrary physical units. This is used internally all
    over proplot, permitting flexible units for various keyword arguments.

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

    numeric : str, default: 'in'
        The units associated with numeric input.
    dest : str, default: `numeric`
        The destination units.
    fontsize : str or float, default: :rc:`font.size` or :rc:`axes.titlesize`
        The font size in points used for scaling. Default is
        :rcraw:`font.size` for ``em`` and ``en`` units and
        :rcraw:`axes.titlesize` for ``Em`` and ``En`` units.
    axes : `~matplotlib.axes.Axes`, optional
        The axes to use for scaling units that look like ``'0.1ax'``.
    figure : `~matplotlib.figure.Figure`, optional
        The figure to use for scaling units that look like ``'0.1fig'``.
        If not provided we try to get the figure from ``axes.figure``.
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


# Deprecation function
def _deprecate_color_func(arg):
    """
    Deprecate the color transformation function.
    """
    if isinstance(arg, str):
        attr = dest = arg
    else:
        attr, dest = arg

    def _wrapper(color, *args, **kwargs):
        warnings._warn_proplot(
            f'pplt.{attr}() was deprecated in version 0.10.0 and will be removed '
            'in a future release. Please use pplt.Color(color).{dest}() instead. '
            'You may also be able to accomplish this transformation by passing keyword '
            'arguments to Color() (see the pplt.Color documentation for details).'
        )
        from .colors import Color
        color = Color(color)
        return getattr(color, dest)(*args, **kwargs)

    _wrapper.__name__ = attr
    __all__.append(attr)
    globals()[attr] = _wrapper


# Add the deprecated functions
for _color_func in (
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
    ('shade', 'scale_luminance'),
    ('saturate', 'scale_saturation'),
):
    _deprecate_color_func(_color_func)
