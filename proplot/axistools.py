#!/usr/bin/env python3
"""
Various axis `~matplotlib.ticker.Formatter` and `~matplotlib.scale.ScaleBase`
classes. Includes constructor functions so that these classes can be selected
with with a shorthand syntax.
"""
import re
from .utils import _notNone
from .rctools import rc
from numbers import Number
from fractions import Fraction
import numpy as np
import numpy.ma as ma
import warnings
import matplotlib.dates as mdates
import matplotlib.projections.polar as mpolar
import matplotlib.ticker as mticker
import matplotlib.scale as mscale
import matplotlib.transforms as mtransforms
__all__ = [
    'formatters', 'locators', 'scales',
    'Formatter', 'Locator', 'Scale',
    'AutoFormatter', 'CutoffScale', 'ExpScale',
    'FracFormatter', 'FuncScale',
    'InverseScale',
    'LinearScale',
    'LogitScale',
    'LogScale',
    'MercatorLatitudeScale', 'PowerScale', 'SimpleFormatter',
    'SineLatitudeScale',
    'SymmetricalLogScale',
]

# Scale preset names and positional args
SCALE_PRESETS = {
    'quadratic': ('power', 2,),
    'cubic': ('power', 3,),
    'quartic': ('power', 4,),
    'height': ('exp', np.e, -1 / 7, 1013.25, True),
    'pressure': ('exp', np.e, -1 / 7, 1013.25, False),
    'db': ('exp', 10, 1, 0.1, True),
    'idb': ('exp', 10, 1, 0.1, False),
    'np': ('exp', np.e, 1, 1, True),
    'inp': ('exp', np.e, 1, 1, False),
}


def Locator(locator, *args, **kwargs):
    """
    Returns a `~matplotlib.ticker.Locator` instance, used to interpret the
    `xlocator`, `xlocator_kw`, `ylocator`, `ylocator_kw`, `xminorlocator`,
    `xminorlocator_kw`, `yminorlocator`, and `yminorlocator_kw` arguments when
    passed to `~proplot.axes.XYAxes.format`, and the `locator`, `locator_kw`
    `minorlocator`, and `minorlocator_kw` arguments when passed to colorbar
    methods wrapped by `~proplot.wrappers.colorbar_wrapper`.

    Parameters
    ----------
    locator : `~matplotlib.ticker.Locator`, str, float, or list of float
        If `~matplotlib.ticker.Locator`, the object is returned.

        If number, specifies the *multiple* used to define tick separation.
        Returns a `~matplotlib.ticker.MultipleLocator` instance.

        If list of numbers, these points are ticked. Returns a
        `~matplotlib.ticker.FixedLocator` instance.

        If string, a dictionary lookup is performed (see below table).

        ======================  ============================================  =========================================================================================
        Key                     Class                                         Description
        ======================  ============================================  =========================================================================================
        ``'null'``, ``'none'``  `~matplotlib.ticker.NullLocator`              No ticks
        ``'auto'``              `~matplotlib.ticker.AutoLocator`              Major ticks at sensible locations
        ``'minor'``             `~matplotlib.ticker.AutoMinorLocator`         Minor ticks at sensible locations
        ``'date'``              `~matplotlib.dates.AutoDateLocator`           Default tick locations for datetime axes
        ``'log'``               `~matplotlib.ticker.LogLocator` preset        For log-scale axes, ticks on each power of the base
        ``'logminor'``          `~matplotlib.ticker.LogLocator` preset        For log-scale axes, ticks on the 1st through 9th multiples of each power of the base
        ``'maxn'``              `~matplotlib.ticker.MaxNLocator`              No more than ``N`` ticks at sensible locations
        ``'linear'``            `~matplotlib.ticker.LinearLocator`            Exactly ``N`` ticks encompassing the axis limits, spaced as ``numpy.linspace(lo, hi, N)``
        ``'multiple'``          `~matplotlib.ticker.MultipleLocator`          Ticks every ``N`` step away from zero
        ``'fixed'``             `~matplotlib.ticker.FixedLocator`             Ticks at these exact locations
        ``'index'``             `~matplotlib.ticker.IndexLocator`             Ticks on the non-negative integers
        ``'symlog'``            `~matplotlib.ticker.SymmetricalLogLocator`    Ticks for symmetrical log-scale axes
        ``'logit'``             `~matplotlib.ticker.LogitLocator`             Ticks for logit-scale axes
        ``'theta'``             `~matplotlib.projections.polar.ThetaLocator`  Like the base locator but default locations are every `numpy.pi`/8 radians
        ``'year'``              `~matplotlib.dates.YearLocator`               Ticks every ``N`` years
        ``'month'``             `~matplotlib.dates.MonthLocator`              Ticks every ``N`` months
        ``'weekday'``           `~matplotlib.dates.WeekdayLocator`            Ticks every ``N`` weekdays
        ``'day'``               `~matplotlib.dates.DayLocator`                Ticks every ``N`` days
        ``'hour'``              `~matplotlib.dates.HourLocator`               Ticks every ``N`` hours
        ``'minute'``            `~matplotlib.dates.MinuteLocator`             Ticks every ``N`` minutes
        ``'second'``            `~matplotlib.dates.SecondLocator`             Ticks every ``N`` seconds
        ``'microsecond'``       `~matplotlib.dates.MicrosecondLocator`        Ticks every ``N`` microseconds
        ======================  ============================================  =========================================================================================

    *args, **kwargs
        Passed to the `~matplotlib.ticker.Locator` class.

    Returns
    -------
    `~matplotlib.ticker.Locator`
        A `~matplotlib.ticker.Locator` instance.
    """  # noqa
    if isinstance(locator, mticker.Locator):
        return locator
    # Pull out extra args
    if np.iterable(locator) and not isinstance(locator, str) and not all(
            isinstance(num, Number) for num in locator):
        locator, args = locator[0], (*locator[1:], *args)
    # Get the locator
    if isinstance(locator, str):  # dictionary lookup
        # Shorthands and defaults
        if locator == 'logminor':
            locator = 'log'
            kwargs.setdefault('subs', np.arange(10))
        elif locator == 'index':
            args = args or (1,)
            if len(args) == 1:
                args = (*args, 0)
        # Lookup
        if locator not in locators:
            raise ValueError(
                f'Unknown locator {locator!r}. Options are '
                + ', '.join(map(repr, locators.keys())) + '.')
        locator = locators[locator](*args, **kwargs)
    elif isinstance(locator, Number):  # scalar variable
        locator = mticker.MultipleLocator(locator, *args, **kwargs)
    elif np.iterable(locator):
        locator = mticker.FixedLocator(
            np.sort(locator), *args, **kwargs)  # not necessary
    else:
        raise ValueError(f'Invalid locator {locator!r}.')
    return locator


def Formatter(formatter, *args, date=False, **kwargs):
    """
    Returns a `~matplotlib.ticker.Formatter` instance, used to interpret the
    `xformatter`, `xformatter_kw`, `yformatter`, and `yformatter_kw` arguments
    when passed to `~proplot.axes.XYAxes.format`, and the `formatter`
    and `formatter_kw` arguments when passed to colorbar methods wrapped by
    `~proplot.wrappers.colorbar_wrapper`.

    Parameters
    ----------
    formatter : `~matplotlib.ticker.Formatter`, str, list of str, or function
        If `~matplotlib.ticker.Formatter`, the object is returned.

        If list of strings, ticks are labeled with these strings. Returns a
        `~matplotlib.ticker.FixedFormatter` instance.

        If function, labels will be generated using this function. Returns a
        `~matplotlib.ticker.FuncFormatter` instance.

        If string, there are 4 possibilities:

        1. If string contains ``'%'`` and `date` is ``False``, ticks will be
           formatted using the C-notation ``string % number`` method. See
           `this page \
<https://docs.python.org/3/library/stdtypes.html#printf-style-string-formatting>`__
           for a review.
        2. If string contains ``'%'`` and `date` is ``True``, datetime
           ``string % number`` formatting is used. See
           `this page \
<https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes>`__
           for a review.
        3. If string contains ``{x}`` or ``{x:...}``, ticks will be
           formatted by calling ``string.format(x=number)``.
        4. In all other cases, a dictionary lookup is performed
           (see below table).

        ======================  ==============================================  ===================================================================================================================================
        Key                     Class                                           Description
        ======================  ==============================================  ===================================================================================================================================
        ``'null'``, ``'none'``  `~matplotlib.ticker.NullFormatter`              No tick labels
        ``'auto'``              `AutoFormatter`                                 New default tick labels for axes
        ``'simple'``            `SimpleFormatter`                               New default tick labels for e.g. contour labels
        ``'frac'``              `FracFormatter`                                 Rational fractions
        ``'date'``              `~matplotlib.dates.AutoDateFormatter`           Default tick labels for datetime axes
        ``'datestr'``           `~matplotlib.dates.DateFormatter`               Date formatting with C-style ``string % format`` notation
        ``'concise'``           `~matplotlib.dates.ConciseDateFormatter`        More concise date labels introduced in `matplotlib 3.1 <https://matplotlib.org/3.1.0/users/whats_new.html#concisedateformatter>`__
        ``'scalar'``            `~matplotlib.ticker.ScalarFormatter`            Old default tick labels for axes
        ``'strmethod'``         `~matplotlib.ticker.StrMethodFormatter`         From the ``string.format`` method
        ``'formatstr'``         `~matplotlib.ticker.FormatStrFormatter`         From C-style ``string % format`` notation
        ``'log'``, ``'sci'``    `~matplotlib.ticker.LogFormatterSciNotation`    For log-scale axes with scientific notation
        ``'math'``              `~matplotlib.ticker.LogFormatterMathtext`       For log-scale axes with math text
        ``'logit'``             `~matplotlib.ticker.LogitFormatter`             For logistic-scale axes
        ``'eng'``               `~matplotlib.ticker.EngFormatter`               Engineering notation
        ``'percent'``           `~matplotlib.ticker.PercentFormatter`           Trailing percent sign
        ``'fixed'``             `~matplotlib.ticker.FixedFormatter`             List of strings
        ``'index'``             `~matplotlib.ticker.IndexFormatter`             List of strings corresponding to non-negative integer positions along the axis
        ``'theta'``             `~matplotlib.projections.polar.ThetaFormatter`  Formats radians as degrees, with a degree symbol
        ``'pi'``                `FracFormatter` preset                          Fractions of :math:`\\pi`
        ``'e'``                 `FracFormatter` preset                          Fractions of *e*
        ``'deg'``               `SimpleFormatter` preset                        Trailing degree symbol
        ``'deglon'``            `SimpleFormatter` preset                        Trailing degree symbol and cardinal "WE" indicator
        ``'deglat'``            `SimpleFormatter` preset                        Trailing degree symbol and cardinal "SN" indicator
        ``'lon'``               `SimpleFormatter` preset                        Cardinal "WE" indicator
        ``'lat'``               `SimpleFormatter` preset                        Cardinal "SN" indicator
        ======================  ==============================================  ===================================================================================================================================

    date : bool, optional
        Toggles the behavior when `formatter` contains a ``'%'`` sign (see
        above).
    *args, **kwargs
        Passed to the `~matplotlib.ticker.Formatter` class.

    Returns
    -------
    `~matplotlib.ticker.Formatter`
        A `~matplotlib.ticker.Formatter` instance.
    """  # noqa
    if isinstance(formatter, mticker.Formatter):  # formatter object
        return formatter
    # Pull out extra args
    if np.iterable(formatter) and not isinstance(formatter, str) and not all(
            isinstance(item, str) for item in formatter):
        formatter, args = formatter[0], [*formatter[1:], *args]
    # Get the formatter
    if isinstance(formatter, str):  # assumption is list of strings
        # Format strings
        if re.search(r'{x?(:.+)?}', formatter):
            formatter = mticker.StrMethodFormatter(
                formatter, *args, **kwargs)  # new-style .format() form
        elif '%' in formatter:
            if date:
                formatter = mdates.DateFormatter(
                    formatter, *args, **kwargs)  # %-style, dates
            else:
                formatter = mticker.FormatStrFormatter(
                    formatter, *args, **kwargs)  # %-style, numbers
        else:
            # Fraction shorthands
            if formatter in ('pi', 'e'):
                if formatter == 'pi':
                    symbol, number = r'$\pi$', np.pi
                else:
                    symbol, number = '$e$', np.e
                kwargs.setdefault('symbol', symbol)
                kwargs.setdefault('number', number)
                formatter = 'frac'
            # Cartographic shorthands
            if formatter in ('deg', 'deglon', 'deglat', 'lon', 'lat'):
                negpos, suffix = None, None
                if 'deg' in formatter:
                    suffix = '\N{DEGREE SIGN}'
                if 'lat' in formatter:
                    negpos = 'SN'
                if 'lon' in formatter:
                    negpos = 'WE'
                kwargs.setdefault('suffix', suffix)
                kwargs.setdefault('negpos', negpos)
                formatter = 'simple'
            # Lookup
            if formatter not in formatters:
                raise ValueError(
                    f'Unknown formatter {formatter!r}. Options are '
                    + ', '.join(map(repr, formatters.keys())) + '.')
            formatter = formatters[formatter](*args, **kwargs)
    elif callable(formatter):
        formatter = mticker.FuncFormatter(formatter, *args, **kwargs)
    elif np.iterable(formatter):  # list of strings on the major ticks
        formatter = mticker.FixedFormatter(formatter)
    else:
        raise ValueError(f'Invalid formatter {formatter!r}.')
    return formatter


def Scale(scale, *args, **kwargs):
    """
    Returns a `~matplotlib.scale.ScaleBase` instance, used to interpret the
    `xscale`, `xscale_kw`, `yscale`, and `yscale_kw` arguments when passed to
    `~proplot.axes.XYAxes.format`.

    Parameters
    ----------
    scale : `~matplotlib.scale.ScaleBase`, str, (str, ...), or class
        If `~matplotlib.scale.ScaleBase`, the object is returned.

        If string, this is the registered scale name or scale "preset" (see
        below table). If an iterable is passed with the scale name as the
        first element, the subsequent items are passed to the scale class as
        positional arguments.

        =================  ===============================  =======================================================================
        Key                Class                            Description
        =================  ===============================  =======================================================================
        ``'linear'``       `~matplotlib.scale.LinearScale`  Linear
        ``'log'``          `LogScale`                       Logarithmic
        ``'symlog'``       `SymmetricalLogScale`            Logarithmic beyond finite space around zero
        ``'logit'``        `~matplotlib.scale.LogitScale`   Logistic
        ``'inverse'``      `InverseScale`                   Inverse
        ``'function'``     `FuncScale`                      Scale from arbitrary forward and backwards functions
        ``'sine'``         `SineLatitudeScale`              Sine function (in degrees)
        ``'mercator'``     `MercatorLatitudeScale`          Mercator latitude function (in degrees)
        ``'exp'``          `ExpScale`                       Arbitrary exponential function
        ``'power'``        `PowerScale`                     Arbitrary power function
        ``'cutoff'``       `CutoffScale`                    Arbitrary linear transformations
        ``'quadratic'``    `PowerScale` (preset)            Quadratic function
        ``'cubic'``        `PowerScale` (preset)            Cubic function
        ``'quartic'``      `PowerScale` (preset)            Cubic function
        ``'pressure'``     `ExpScale` (preset)              Height (in km) expressed linear in pressure
        ``'height'``       `ExpScale` (preset)              Pressure (in hPa) expressed linear in height
        ``'db'``           `ExpScale` (preset)              Ratio expressed as `decibels <https://en.wikipedia.org/wiki/Decibel>`__
        ``'np'``           `ExpScale` (preset)              Ratio expressed as `nepers <https://en.wikipedia.org/wiki/Neper>`__
        ``'idb'``          `ExpScale` (preset)              `Decibels <https://en.wikipedia.org/wiki/Decibel>`__ expressed as ratio
        ``'inp'``          `ExpScale` (preset)              `Nepers <https://en.wikipedia.org/wiki/Neper>`__ expressed as ratio
        =================  ===============================  =======================================================================

    *args, **kwargs
        Passed to the `~matplotlib.scale.ScaleBase` class.

    Returns
    -------
    `~matplotlib.scale.ScaleBase`
        The scale instance.
    """  # noqa
    if isinstance(scale, mscale.ScaleBase):
        return scale
    # Pull out extra args
    if np.iterable(scale) and not isinstance(scale, str):
        scale, args = scale[0], (*scale[1:], *args)
    if not isinstance(scale, str):
        raise ValueError(f'Invalid scale name {scale!r}. Must be string.')
    # Get scale preset
    if scale in SCALE_PRESETS:
        if args or kwargs:
            warnings.warn(
                f'Scale {scale!r} is a scale *preset*. Ignoring positional '
                'argument(s): {args} and keyword argument(s): {kwargs}. ')
        scale, *args = SCALE_PRESETS[scale]
    # Get scale
    scale = scale.lower()
    if scale in scales:
        scale = scales[scale]
    else:
        raise ValueError(
            f'Unknown scale or preset {scale!r}. Options are '
            + ', '.join(map(repr, list(scales) + list(SCALE_PRESETS))) + '.')
    axis = _dummy_axis()
    return scale(axis, *args, **kwargs)


class AutoFormatter(mticker.ScalarFormatter):
    """
    The new default formatter, a simple wrapper around
    `~matplotlib.ticker.ScalarFormatter`. Differs from
    `~matplotlib.ticker.ScalarFormatter` in the following ways:

    1. Trims trailing zeros if any exist.
    2. Allows user to specify *range* within which major tick marks
       are labelled.
    3. Allows user to add arbitrary prefix or suffix to every
       tick label string.

    """

    def __init__(self, *args,
                 zerotrim=None, precision=None, tickrange=None,
                 prefix=None, suffix=None, **kwargs):
        """
        Parameters
        ----------
        zerotrim : bool, optional
            Whether to trim trailing zeros.
            Default is :rc:`axes.formatter.zerotrim`.
        precision : float, optional
            The maximum number of digits after the decimal point.
        tickrange : (float, float), optional
            Range within which major tick marks are labelled.
        prefix, suffix : str, optional
            Optional prefix and suffix for all strings.
        *args, **kwargs
            Passed to `matplotlib.ticker.ScalarFormatter`.
        """
        tickrange = tickrange or (-np.inf, np.inf)
        super().__init__(*args, **kwargs)
        zerotrim = _notNone(zerotrim, rc.get('axes.formatter.zerotrim'))
        self._maxprecision = precision
        self._zerotrim = zerotrim
        self._tickrange = tickrange
        self._prefix = prefix or ''
        self._suffix = suffix or ''

    def __call__(self, x, pos=None):
        """
        Convert number to a string.

        Parameters
        ----------
        x : float
            The value.
        pos : float, optional
            The position.
        """
        # Tick range limitation
        eps = abs(x) / 1000
        tickrange = self._tickrange
        if (x + eps) < tickrange[0] or (x - eps) > tickrange[1]:
            return ''  # avoid some ticks
        # Normal formatting
        string = super().__call__(x, pos)
        if self._maxprecision is not None and '.' in string:
            head, tail = string.split('.')
            string = head + '.' + tail[:self._maxprecision]
        if self._zerotrim and '.' in string:
            string = string.rstrip('0').rstrip('.')
        if string == '-0' or string == '\N{MINUS SIGN}0':
            string = '0'
        # Prefix and suffix
        sign = ''
        string = string.replace('-', '\N{MINUS SIGN}')
        if string and string[0] == '\N{MINUS SIGN}':
            sign, string = string[0], string[1:]
        return sign + self._prefix + string + self._suffix


def SimpleFormatter(*args, precision=6,
                    prefix=None, suffix=None, negpos=None, zerotrim=True,
                    **kwargs):
    """
    Replicates features of `AutoFormatter`, but as a simpler
    `~matplotlib.ticker.FuncFormatter` instance. This is more suitable for
    arbitrary number formatting not necessarily associated with any
    `~matplotlib.axis.Axis` instance, e.g. labelling contours.

    Parameters
    ----------
    precision : int, optional
        Maximum number of digits after the decimal point.
    prefix, suffix : str, optional
        Optional prefix and suffix for all strings.
    negpos : str, optional
        Length-2 string that indicates suffix for "negative" and "positive"
        numbers, meant to replace the minus sign. This is useful for
        indicating cardinal geographic coordinates.
    zerotrim : bool, optional
        Whether to trim trailing zeros.
        Default is :rc:`axes.formatter.zerotrim`.
    """
    prefix = prefix or ''
    suffix = suffix or ''
    zerotrim = _notNone(zerotrim, rc['axes.formatter.zerotrim'])

    def f(x, pos):
        # Apply suffix if not on equator/prime meridian
        if not negpos:
            tail = ''
        elif x > 0:
            tail = negpos[1]
        else:
            x *= -1
            tail = negpos[0]
        # Finally use default formatter
        string = ('{:.%df}' % precision).format(x)
        if zerotrim and '.' in string:
            string = string.rstrip('0').rstrip('.')
        if string == '-0' or string == '\N{MINUS SIGN}0':
            string = '0'
        # Prefix and suffix
        sign = ''
        string = string.replace('-', '\N{MINUS SIGN}')
        if string and string[0] == '\N{MINUS SIGN}':
            sign, string = string[0], string[1:]
        return sign + prefix + string + suffix + tail
    return mticker.FuncFormatter(f)


def FracFormatter(symbol='', number=1):
    r"""
    Returns a `~matplotlib.ticker.FuncFormatter` that formats numbers as
    fractions or multiples of some value, e.g. a physical constant.

    This is powered by the python builtin `~fractions.Fraction` class.
    We account for floating point errors using the
    `~fractions.Fraction.limit_denominator` method.

    Parameters
    ----------
    symbol : str
        The symbol, e.g. ``r'$\pi$'``. Default is ``''``.
    number : float
        The value, e.g. `numpy.pi`. Default is ``1``.
    """
    def f(x, pos):  # must accept location argument
        frac = Fraction(x / number).limit_denominator()
        if x == 0:
            string = '0'
        elif frac.denominator == 1:  # denominator is one
            if frac.numerator == 1 and symbol:
                string = f'{symbol:s}'
            elif frac.numerator == -1 and symbol:
                string = f'-{symbol:s}'
            else:
                string = f'{frac.numerator:d}{symbol:s}'
        else:
            if frac.numerator == 1 and symbol:  # numerator is +/-1
                string = f'{symbol:s}/{frac.denominator:d}'
            elif frac.numerator == -1 and symbol:
                string = f'-{symbol:s}/{frac.denominator:d}'
            else:  # and again make sure we use unicode minus!
                string = f'{frac.numerator:d}{symbol:s}/{frac.denominator:d}'
        return string.replace('-', '\N{MINUS SIGN}')
    return mticker.FuncFormatter(f)


def _scale_factory(scale, axis, *args, **kwargs):
    """If `scale` is a `~matplotlib.scale.ScaleBase` instance, nothing is
    done. If it is a registered scale name, that scale is looked up and
    instantiated."""
    if isinstance(scale, mscale.ScaleBase):
        if args or kwargs:
            warnings.warn(f'Ignoring args {args} and keyword args {kwargs}.')
        return scale  # do nothing
    else:
        scale = scale.lower()
        if scale not in scales:
            raise ValueError(
                f'Unknown scale {scale!r}. Options are '
                + ', '.join(map(repr, scales.keys())) + '.')
        return scales[scale](axis, *args, **kwargs)


def _parse_logscale_args(kwargs, *keys):
    """Parses args for `LogScale` and `SymmetricalLogScale` that
    inexplicably require ``x`` and ``y`` suffixes by default."""
    for key in keys:
        value = _notNone(  # issues warning when multiple args passed!
            kwargs.pop(key, None),
            kwargs.pop(key + 'x', None),
            kwargs.pop(key + 'y', None),
            None, names=(key, key + 'x', key + 'y'),
        )
        if value is not None:  # _dummy_axis axis_name is 'x'
            kwargs[key + 'x'] = value
    return kwargs


class _dummy_axis(object):
    """Dummy axis used to initialize scales."""
    # See notes in source code for `~matplotlib.scale.ScaleBase`. All scales
    # accept 'axis' for backwards-compatibility reasons, but it is *virtually
    # unused* except to check for the `axis_name` attribute in log scales to
    # interpret input keyword args!
    # TODO: Submit matplotlib pull request! How has no one fixed this already!
    axis_name = 'x'


class _ScaleBase(object):
    """Mixin scale class that standardizes required methods."""

    def set_default_locators_and_formatters(self, axis, only_if_default=False):
        """
        Apply all locators and formatters defined as attributes on
        initialization, and define defaults for all scales.

        Parameters
        ----------
        axis : `~matplotlib.axis.Axis`
            The axis.
        only_if_default : bool, optional
            Whether to refrain from updating the locators and formatters if
            the axis is currently using non-default versions. Useful if we
            want to avoid overwriting user customization when the scale
            is changed.
        """
        # * We assert isDefault settings, because matplotlib does this in
        #   axis._set_scale but sometimes we need to bypass this method!
        # * Minor locator can be "non default" even when user has not changed
        #   it, due to "turning minor ticks" on and off, so always apply
        #   default if it is currently AutoMinorLocator.
        if getattr(self, '_smart_bounds', None):
            axis.set_smart_bounds(True)  # unnecessary?
        if not only_if_default or axis.isDefault_majloc:
            axis.set_major_locator(
                getattr(self, '_major_locator', None) or Locator('auto')
            )
            axis.isDefault_majloc = True
        if not only_if_default or axis.isDefault_majfmt:
            axis.set_major_formatter(
                getattr(self, '_major_formatter', None) or Formatter('auto')
            )
            axis.isDefault_majfmt = True
        if (not only_if_default or axis.isDefault_minloc or isinstance(
                axis.get_minor_locator(), mticker.AutoMinorLocator)):
            name = axis.axis_name if axis.axis_name in 'xy' else 'x'
            minor = 'minor' if rc.get(name + 'tick.minor.visible') else 'null'
            axis.set_minor_locator(
                getattr(self, '_minor_locator', None) or Locator(minor)
            )
            axis.isDefault_minloc = True
        if (not only_if_default or axis.isDefault_minfmt or isinstance(
                axis.get_minor_formatter(), mticker.NullFormatter)):
            axis.set_minor_formatter(
                getattr(self, '_minor_formatter', None) or Formatter('null')
            )
            axis.isDefault_minfmt = True

    def get_transform(self):
        """Returns the scale transform."""
        return getattr(self, '_transform', mtransforms.IdentityTransform())


class LinearScale(_ScaleBase, mscale.LinearScale):
    """
    As with `~matplotlib.scale.LinearScale`, but applies new default
    major formatter.
    """
    name = 'linear'
    """The registered scale name."""


class LogitScale(_ScaleBase, mscale.LogitScale):
    """
    As with `~matplotlib.scale.LogitScale`, but applies new default
    major formatter.
    """
    name = 'logit'
    """The registered scale name."""


class LogScale(_ScaleBase, mscale.LogScale):
    """
    As with `~matplotlib.scale.LogScale`, but applies new default major
    formatter and fixes the inexplicable choice to have separate "``x``" and
    "``y``" versions of each keyword argument.
    """
    name = 'log'
    """The registered scale name."""

    def __init__(self, axis, **kwargs):
        """
        Parameters
        ----------
        base : float, optional
            The base of the logarithm. Default is ``10``.
        nonpos : {'mask', 'clip'}, optional
            Non-positive values in *x* or *y* can be masked as
            invalid, or clipped to a very small positive number.
        subs : list of int, optional
            Default tick locations are on these multiples of each power
            of the base. For example, ``subs=(1,2,5)`` draws ticks on 1, 2, 5,
            10, 20, 50, etc.
        basex, basey, nonposx, nonposy, subsx, subsy
            Aliases for the above keywords. These used to be conditional
            on the *name* of the axis...... yikes.
        """
        kwargs = _parse_logscale_args(kwargs, 'base', 'nonpos', 'subs')
        super().__init__(axis, **kwargs)
        # self._major_formatter = Formatter('log')
        self._major_locator = Locator('log', base=self.base)
        self._minor_locator = Locator('log', base=self.base, subs=self.subs)


class SymmetricalLogScale(_ScaleBase, mscale.SymmetricalLogScale):
    """
    As with `~matplotlib.scale.SymmetricLogScale`, but applies new default
    major formatter and fixes the inexplicable choice to have separate "``x``"
    and "``y``" versions of each keyword argument.
    """
    name = 'symlog'
    """The registered scale name."""

    def __init__(self, axis, **kwargs):
        """
        Parameters
        ----------
        base : float, optional
            The base of the logarithm. Default is ``10``.
        linthresh : float, optional
            Defines the range ``(-linthresh, linthresh)``, within which the
            plot is linear.  This avoids having the plot go to infinity around
            zero. Defaults to 2.
        linscale : float, optional
            This allows the linear range ``(-linthresh, linthresh)`` to be
            stretched relative to the logarithmic range. Its value is the
            number of decades to use for each half of the linear range. For
            example, when `linscale` is ``1`` (the default), the space used
            for the positive and negative halves of the linear range will be
            equal to one decade in the logarithmic range.
        subs : sequence of int, optional
            Default minor tick locations are on these multiples of each power
            of the base. For example, ``subs=[1,2,5]`` draws ticks on 1, 2, 5,
            10, 20, 50, 100, 200, 500, etc.
        basex, basey, linthreshx, linthreshy, linscalex, linscaley, \
subsx, subsy
            Aliases for the above keywords. These used to be conditional
            on the *name* of the axis...... yikes.
        """
        kwargs = _parse_logscale_args(kwargs,
                                      'base', 'linthresh', 'linscale', 'subs')
        super().__init__(axis, **kwargs)
        # Note the symlog locator gets base and linthresh from the transform
        # self._major_formatter = Formatter('symlog'))
        self._major_locator = Locator('symlog', transform=self.get_transform())
        self._minor_locator = Locator('symlog', transform=self.get_transform(),
                                      subs=self.subs)


class FuncScale(_ScaleBase, mscale.ScaleBase):
    """
    Arbitrary scale with user-supplied forward and inverse functions and
    arbitrary additional transform applied thereafter. Input is a tuple
    of functions and, optionally, a `~matplotlib.transforms.Transform` or
    `~matplotlib.scale.ScaleBase` instance.
    """
    name = 'function'
    """The registered scale name."""

    def __init__(self, axis, functions, transform=None, scale=None,
                 major_locator=None, minor_locator=None,
                 major_formatter=None, minor_formatter=None,
                 ):
        """
        Parameters
        ----------
        axis : `~matplotlib.axis.Axis`
            The axis, required for compatibility reasons.
        functions : (function, function) or `~matplotlib.scale.ScaleBase`
            Length-2 tuple of forward and inverse functions, or another
            `~matplotlib.scale.ScaleBase` from which the functions are drawn.
        transform : `~matplotlib.transforms.Transform`, optional
            Additional transform applied after the forward function
            and before the inverse function.
        major_locator, minor_locator : `~matplotlib.ticker.Locator`, optional
            The default major and minor locator. By default these are the same
            as `~matplotlib.scale.LinearScale`.
        major_formatter, minor_formatter : `~matplotlib.ticker.Formatter`, \
optional
            The default major and minor formatter. By default these are the
            same as `~matplotlib.scale.LinearScale`.
        """
        if np.iterable(functions) and len(functions) == 2 and all(
                callable(ifunction) for ifunction in functions):
            forward, inverse = functions
        else:
            raise ValueError(
                f'scale needs length-2 list of forward and inverse '
                f'functions, not {functions!r}.')
        functransform = FuncTransform(forward, inverse)
        if transform is not None:
            if isinstance(transform, mtransforms.Transform):
                functransform = functransform + transform
            else:
                raise ValueError(
                    f'transform {transform!r} must be a Transform instance, '
                    f'not {type(transform)!r}.')
        self._transform = functransform
        if major_locator:
            self._major_locator = major_locator
        if minor_locator:
            self._minor_locator = minor_locator
        if major_formatter:
            self._major_formatter = major_formatter
        if minor_formatter:
            self._minor_formatter = minor_formatter


class FuncTransform(mtransforms.Transform):
    # Arbitrary forward and inverse transform
    # Mostly copied from matplotlib
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True

    def __init__(self, forward, inverse):
        super().__init__()
        if callable(forward) and callable(inverse):
            self._forward = forward
            self._inverse = inverse
        else:
            raise ValueError('arguments to FuncTransform must be functions')

    def inverted(self):
        return FuncTransform(self._inverse, self._forward)

    def transform_non_affine(self, values):
        return self._forward(values)


class PowerScale(_ScaleBase, mscale.ScaleBase):
    r"""
    Returns a "power scale" that performs the transformation

    .. math::

        x^{c}

    """
    name = 'power'
    """The registered scale name."""

    def __init__(self, axis,
                 power=1, inverse=False, *, minpos=1e-300, **kwargs):
        """
        Parameters
        ----------
        axis : `~matplotlib.axis.Axis`
            The axis, required for compatibility reasons.
        power : float, optional
            The power :math:`c` to which :math:`x` is raised.
        inverse : bool, optional
            If ``True``, the "forward" direction performs
            the inverse operation :math:`x^{1/c}`.
        minpos : float, optional
            The minimum permissible value, used to truncate negative values.
        """
        super().__init__(axis)
        if not inverse:
            self._transform = PowerTransform(power, minpos)
        else:
            self._transform = InvertedPowerTransform(power, minpos)

    def limit_range_for_scale(self, vmin, vmax, minpos):
        """Returns the range *vmin* and *vmax* limited to positive numbers."""
        return max(vmin, minpos), max(vmax, minpos)


class PowerTransform(mtransforms.Transform):
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True

    def __init__(self, power, minpos):
        super().__init__()
        self.minpos = minpos
        self._power = power

    def inverted(self):
        return InvertedPowerTransform(self._power, self.minpos)

    def transform(self, a):
        aa = np.array(a)
        aa[aa <= self.minpos] = self.minpos  # necessary
        return np.power(np.array(a), self._power)

    def transform_non_affine(self, a):
        return self.transform(a)


class InvertedPowerTransform(mtransforms.Transform):
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True

    def __init__(self, power, minpos):
        super().__init__()
        self.minpos = minpos
        self._power = power

    def inverted(self):
        return PowerTransform(self._power, self.minpos)

    def transform(self, a):
        aa = np.array(a)
        aa[aa <= self.minpos] = self.minpos  # necessary
        return np.power(np.array(a), 1 / self._power)

    def transform_non_affine(self, a):
        return self.transform(a)


class ExpScale(_ScaleBase, mscale.ScaleBase):
    r"""
    An "exponential scale". When `inverse` is ``False`` (the default), this
    performs the transformation

    .. math::

        Ca^{bx}

    where the constants :math:`a`, :math:`b`, and :math:`C` are set by the
    input (see below). When `inverse` is ``True``, this performs the inverse
    transformation

    .. math::

        (\log_a(x) - \log_a(C))/b

    which in appearence is equivalent to `LogScale` since it is just a linear
    transformation of the logarithm.
    """
    name = 'exp'
    """The registered scale name."""

    def __init__(
            self, axis,
            a=np.e, b=1, c=1, inverse=False, minpos=1e-300,
            **kwargs):
        """
        Parameters
        ----------
        axis : `~matplotlib.axis.Axis`
            The axis, required for compatibility reasons.
        a : float, optional
            The base of the exponential, i.e. the :math:`a` in :math:`Ca^{bx}`.
        b : float, optional
            The scale for the exponent, i.e. the :math:`b` in :math:`Ca^{bx}`.
        c : float, optional
            The coefficient of the exponential, i.e. the :math:`C`
            in :math:`Ca^{bx}`.
        minpos : float, optional
            The minimum permissible value, used to truncate negative values.
        inverse : bool, optional
            If ``True``, the "forward" direction performs the inverse
            operation.
        """
        super().__init__(axis)
        if not inverse:
            self._transform = ExpTransform(a, b, c, minpos)
        else:
            self._transform = InvertedExpTransform(a, b, c, minpos)

    def limit_range_for_scale(self, vmin, vmax, minpos):
        """Returns the range *vmin* and *vmax* limited to positive numbers."""
        return max(vmin, minpos), max(vmax, minpos)


class ExpTransform(mtransforms.Transform):
    # Arbitrary exponential function
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True

    def __init__(self, a, b, c, minpos):
        super().__init__()
        self.minpos = minpos
        self._a = a
        self._b = b
        self._c = c

    def inverted(self):
        return InvertedExpTransform(self._a, self._b, self._c, self.minpos)

    def transform(self, a):
        return self._c * np.power(self._a, self._b * np.array(a))

    def transform_non_affine(self, a):
        return self.transform(a)


class InvertedExpTransform(mtransforms.Transform):
    # Inverse exponential transform
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True

    def __init__(self, a, b, c, minpos):
        super().__init__()
        self.minpos = minpos
        self._a = a
        self._b = b
        self._c = c

    def inverted(self):
        return ExpTransform(self._a, self._b, self._c, self.minpos)

    def transform(self, a):
        aa = np.array(a)
        aa[aa <= self.minpos] = self.minpos  # necessary
        return np.log(aa / self._c) / (self._b * np.log(self._a))

    def transform_non_affine(self, a):
        return self.transform(a)


class CutoffScale(_ScaleBase, mscale.ScaleBase):
    """Axis scale with arbitrary cutoffs that "accelerate" parts of the
    axis, "decelerate" parts of the axes, or discretely jumps between
    numbers.

    If `upper` is not provided, you have the following two possibilities.

    1. If `scale` is greater than 1, the axis is "accelerated" to the right
       of `lower`.
    2. If `scale` is less than 1, the axis is "decelerated" to the right
       of `lower`.

    If `upper` is provided, you have the following three possibilities.

    1. If `scale` is `numpy.inf`, this puts a cliff between `lower` and
       `upper`. The axis discretely jumps from `lower` to `upper`.
    2. If `scale` is greater than 1, the axis is "accelerated" between `lower`
       and `upper`.
    3. If `scale` is less than 1, the axis is "decelerated" between `lower`
       and `upper`.
    """
    name = 'cutoff'
    """The registered scale name."""

    def __init__(self, axis, scale, lower, upper=None, **kwargs):
        """
        Parameters
        ----------
        axis : `~matplotlib.axis.Axis`
            The matplotlib axis. Required for compatibility reasons.
        scale : float
            Value satisfying ``0 < scale <= numpy.inf``. If `scale` is
            greater than ``1``, values to the right of `lower`, or
            between `lower` and `upper`, are "accelerated". Otherwise, values
            are "decelerated". Infinity represents a discrete jump.
        lower : float
            The first cutoff point.
        upper : float, optional
            The second cutoff point (optional, see above).

        Todo
        ----
        Add method for drawing diagonal "cutoff" strokes. See
        `this post <https://stackoverflow.com/a/5669301/4970632>`__
        for class-based and multi-axis solutions.
        """
        # Note the space between 1-9 in Paul's answer is because actual
        # cutoffs were 0.1 away (and tick locations are 0.2 apart).
        if scale < 0:
            raise ValueError('Scale must be a positive float.')
        if upper is None and scale == np.inf:
            raise ValueError(
                'For a discrete jump, need both lower and upper bounds. '
                'You just provided lower bounds.')
        super().__init__(axis)
        self._transform = CutoffTransform(scale, lower, upper)


class CutoffTransform(mtransforms.Transform):
    # Create transform object
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True

    def __init__(self, scale, lower, upper=None):
        super().__init__()
        self._scale = scale
        self._lower = lower
        self._upper = upper

    def inverted(self):
        return InvertedCutoffTransform(self._scale, self._lower, self._upper)

    def transform(self, a):
        a = np.array(a)  # very numpy array
        aa = a.copy()
        scale = self._scale
        lower = self._lower
        upper = self._upper
        if upper is None:  # just scale between 2 segments
            m = (a > lower)
            aa[m] = a[m] - (a[m] - lower) * (1 - 1 / scale)
        elif lower is None:
            m = (a < upper)
            aa[m] = a[m] - (upper - a[m]) * (1 - 1 / scale)
        else:
            m1 = (a > lower)
            m2 = (a > upper)
            m3 = (a > lower) & (a < upper)
            if scale == np.inf:
                aa[m1] = a[m1] - (upper - lower)
                aa[m3] = lower
            else:
                aa[m2] = a[m2] - (upper - lower) * (1 - 1 / scale)
                aa[m3] = a[m3] - (a[m3] - lower) * (1 - 1 / scale)
        return aa

    def transform_non_affine(self, a):
        return self.transform(a)


class InvertedCutoffTransform(mtransforms.Transform):
    # Inverse of cutoff transform
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True

    def __init__(self, scale, lower, upper=None):
        super().__init__()
        self._scale = scale
        self._lower = lower
        self._upper = upper

    def inverted(self):
        return CutoffTransform(self._scale, self._lower, self._upper)

    def transform(self, a):
        a = np.array(a)
        aa = a.copy()
        scale = self._scale
        lower = self._lower
        upper = self._upper
        if upper is None:
            m = (a > lower)
            aa[m] = a[m] + (a[m] - lower) * (1 - 1 / scale)
        elif lower is None:
            m = (a < upper)
            aa[m] = a[m] + (upper - a[m]) * (1 - 1 / scale)
        else:
            n = (upper - lower) * (1 - 1 / scale)
            m1 = (a > lower)
            m2 = (a > upper - n)
            m3 = (a > lower) & (a < (upper - n))
            if scale == np.inf:
                aa[m1] = a[m1] + (upper - lower)
            else:
                aa[m2] = a[m2] + n
                aa[m3] = a[m3] + (a[m3] - lower) * (1 - 1 / scale)
        return aa

    def transform_non_affine(self, a):
        return self.transform(a)


class MercatorLatitudeScale(_ScaleBase, mscale.ScaleBase):
    """
    Scales axis as with latitude in the `Mercator projection \
<http://en.wikipedia.org/wiki/Mercator_projection>`__.
    Adapted from `this example \
<https://matplotlib.org/examples/api/custom_scale_example.html>`__.
    The scale function is as follows.

    .. math::

        y = \\ln(\\tan(\\pi x/180) + \\sec(\\pi x/180))

    The inverse scale function is as follows.

    .. math::

        x = 180\\arctan(\\sinh(y))/\\pi

    """
    name = 'mercator'
    """The registered scale name."""

    def __init__(self, axis, *, thresh=85.0):
        """
        Parameters
        ----------
        axis : `~matplotlib.axis.Axis`
            The matplotlib axis. Required for compatibility reasons.
        thresh : float, optional
            Threshold between 0 and 90, used to constrain axis limits between
            ``-thresh`` and ``+thresh``.
        """
        super().__init__(axis)
        if thresh >= 90.0:
            raise ValueError('Threshold "thresh" must be <=90.')
        self._thresh = thresh
        self._transform = MercatorLatitudeTransform(thresh)
        self._major_formatter = Formatter('deg')
        self._smart_bounds = True

    def limit_range_for_scale(self, vmin, vmax, minpos):
        """Returns the range *vmin* and *vmax* limited to some range within
        +/-90 degrees (exclusive)."""
        return max(vmin, -self._thresh), min(vmax, self._thresh)


class MercatorLatitudeTransform(mtransforms.Transform):
    # Default attributes
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True

    def __init__(self, thresh):
        super().__init__()
        self._thresh = thresh

    def inverted(self):
        return InvertedMercatorLatitudeTransform(self._thresh)

    def transform_non_affine(self, a):
        # With safeguards
        # TODO: Can improve this?
        a = np.deg2rad(a)  # convert to radians
        m = ma.masked_where((a < -self._thresh) | (a > self._thresh), a)
        if m.mask.any():
            return ma.log(np.abs(ma.tan(m) + 1 / ma.cos(m)))
        else:
            return np.log(np.abs(np.tan(a) + 1 / np.cos(a)))


class InvertedMercatorLatitudeTransform(mtransforms.Transform):
    # As above, but for the inverse transform
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True

    def __init__(self, thresh):
        super().__init__()
        self._thresh = thresh

    def inverted(self):
        return MercatorLatitudeTransform(self._thresh)

    def transform_non_affine(self, a):
        # m = ma.masked_where((a < -self._thresh) | (a > self._thresh), a)
        # always assume in first/fourth quadrant, i.e. go from -pi/2 to pi/2
        return np.rad2deg(np.arctan2(1, np.sinh(a)))


class SineLatitudeScale(_ScaleBase, mscale.ScaleBase):
    r"""
    Scales axis to be linear in the *sine* of *x* in degrees.
    The scale function is as follows.

    .. math::

        y = \sin(\pi x/180)

    The inverse scale function is as follows.

    .. math::

        x = 180\arcsin(y)/\pi
    """
    name = 'sine'
    """The registered scale name."""

    def __init__(self, axis):
        """
        Parameters
        ----------
        axis : `~matplotlib.axis.Axis`
            The matplotlib axis. Required for compatibility reasons.
        """
        super().__init__(axis)
        self._transform = SineLatitudeTransform()
        self._major_formatter = Formatter('deg')
        self._smart_bounds = True

    def limit_range_for_scale(self, vmin, vmax, minpos):
        """Returns the range *vmin* and *vmax* limited to some range within
        +/-90 degrees (inclusive)."""
        return max(vmin, -90), min(vmax, 90)


class SineLatitudeTransform(mtransforms.Transform):
    # Default attributes
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True

    def __init__(self):
        # Initialize, declare attribute
        super().__init__()

    def inverted(self):
        return InvertedSineLatitudeTransform()

    def transform_non_affine(self, a):
        # With safeguards
        # TODO: Can improve this?
        with np.errstate(invalid='ignore'):  # NaNs will always be False
            m = (a >= -90) & (a <= 90)
        if not m.all():
            aa = ma.masked_where(~m, a)
            return ma.sin(np.deg2rad(aa))
        else:
            return np.sin(np.deg2rad(a))


class InvertedSineLatitudeTransform(mtransforms.Transform):
    # Inverse of SineLatitudeTransform
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True

    def __init__(self):
        super().__init__()

    def inverted(self):
        return SineLatitudeTransform()

    def transform_non_affine(self, a):
        # Clipping, instead of setting invalid
        # NOTE: Using ma.arcsin below caused super weird errors, dun do that
        aa = a.copy()
        return np.rad2deg(np.arcsin(aa))


class InverseScale(_ScaleBase, mscale.ScaleBase):
    r"""
    Scales axis to be linear in the *inverse* of *x*. The scale
    function and inverse scale function are as follows.

    .. math::

        y = x^{-1}

    """
    # Unlike log-scale, we can't just warp the space between
    # the axis limits -- have to actually change axis limits. Also this
    # scale will invert and swap the limits you provide. Weird!
    name = 'inverse'
    """The registered scale name."""

    def __init__(self, axis, **kwargs):
        """
        Parameters
        ----------
        axis : `~matplotlib.axis.Axis`
            The matplotlib axis. Required for compatibility reasons.
        """
        super().__init__(axis)
        self._transform = InverseTransform()
        self._major_locator = Locator('log', base=10, subs=(1, 2, 5))
        self._minor_locator = Locator('log', base=10, subs='auto')
        self._smart_bounds = True
        # self._minor_formatter = Fromatter('log')

    def limit_range_for_scale(self, vmin, vmax, minpos):
        """Returns the range *vmin* and *vmax* limited to positive numbers."""
        return max(vmin, minpos), max(vmax, minpos)


class InverseTransform(mtransforms.Transform):
    # Create transform object
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True

    def __init__(self):
        super().__init__()

    def inverted(self):
        return InverseTransform()

    def transform(self, a):
        a = np.array(a)
        # f = np.abs(a) <= self.minpos # attempt for negative-friendly
        # aa[f] = np.sign(a[f])*self.minpos
        with np.errstate(divide='ignore', invalid='ignore'):
            return 1.0 / a

    def transform_non_affine(self, a):
        return self.transform(a)


#: The registered scale names and their associated
#: `~matplotlib.scale.ScaleBase` classes. See `Scale` for a table.
scales = mscale._scale_mapping

#: Mapping of strings to `~matplotlib.ticker.Locator` classes. See
#: `Locator` for a table."""
locators = {
    'none': mticker.NullLocator,
    'null': mticker.NullLocator,
    'auto': mticker.AutoLocator,
    'log': mticker.LogLocator,
    'maxn': mticker.MaxNLocator,
    'linear': mticker.LinearLocator,
    'multiple': mticker.MultipleLocator,
    'fixed': mticker.FixedLocator,
    'index': mticker.IndexLocator,
    'symlog': mticker.SymmetricalLogLocator,
    'logit': mticker.LogitLocator,
    'minor': mticker.AutoMinorLocator,
    'date': mdates.AutoDateLocator,
    'microsecond': mdates.MicrosecondLocator,
    'second': mdates.SecondLocator,
    'minute': mdates.MinuteLocator,
    'hour': mdates.HourLocator,
    'day': mdates.DayLocator,
    'weekday': mdates.WeekdayLocator,
    'month': mdates.MonthLocator,
    'year': mdates.YearLocator,
}
if hasattr(mpolar, 'ThetaLocator'):
    locators['theta'] = mpolar.ThetaLocator

#: Mapping of strings to `~matplotlib.ticker.Formatter` classes. See
#: `Formatter` for a table.
formatters = {  # note default LogFormatter uses ugly e+00 notation
    'auto': AutoFormatter,
    'frac': FracFormatter,
    'simple': SimpleFormatter,
    'date': mdates.AutoDateFormatter,
    'datestr': mdates.DateFormatter,
    'scalar': mticker.ScalarFormatter,
    'none': mticker.NullFormatter,
    'null': mticker.NullFormatter,
    'strmethod': mticker.StrMethodFormatter,
    'formatstr': mticker.FormatStrFormatter,
    'log': mticker.LogFormatterSciNotation,
    'sci': mticker.LogFormatterSciNotation,
    'math': mticker.LogFormatterMathtext,
    'logit': mticker.LogitFormatter,
    'eng': mticker.EngFormatter,
    'percent': mticker.PercentFormatter,
    'index': mticker.IndexFormatter,
}
if hasattr(mdates, 'ConciseDateFormatter'):
    formatters['concise'] = mdates.ConciseDateFormatter
if hasattr(mpolar, 'ThetaFormatter'):
    formatters['theta'] = mpolar.ThetaFormatter

# Monkey patch. Force scale_factory to accept ScaleBase instances, so that
# set_xscale and set_yscale can accept scales returned by the Scale constructor
if mscale.scale_factory is not _scale_factory:
    mscale.scale_factory = _scale_factory

# Custom scales and overrides
mscale.register_scale(CutoffScale)
mscale.register_scale(ExpScale)
mscale.register_scale(LogScale)
mscale.register_scale(LinearScale)
mscale.register_scale(LogitScale)
mscale.register_scale(FuncScale)
mscale.register_scale(PowerScale)
mscale.register_scale(SymmetricalLogScale)
mscale.register_scale(InverseScale)
mscale.register_scale(SineLatitudeScale)
mscale.register_scale(MercatorLatitudeScale)
