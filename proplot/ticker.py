#!/usr/bin/env python3
"""
Various `~matplotlib.ticker.Locator` and `~matplotlib.ticker.Formatter`
classes.
"""
import numpy as np
import matplotlib.ticker as mticker
from fractions import Fraction
from .internals import ic  # noqa: F401
from .internals import _not_none

__all__ = [
    'AutoFormatter',
    'FracFormatter',
    'SigFigFormatter',
    'SimpleFormatter',
]



def Locator(locator, *args, **kwargs):
    """
    Return a `~matplotlib.ticker.Locator` instance. This function is used to
    interpret the `xlocator`, `xlocator_kw`, `ylocator`, `ylocator_kw`, `xminorlocator`,
    `xminorlocator_kw`, `yminorlocator`, and `yminorlocator_kw` arguments when
    passed to `~proplot.axes.XYAxes.format`, and the `locator`, `locator_kw`
    `minorlocator`, and `minorlocator_kw` arguments when passed to colorbar
    methods wrapped by `~proplot.wrappers.colorbar_wrapper`.

    Parameters
    ----------
    locator : `~matplotlib.ticker.Locator`, str, float, or list of float
        The axis locator specification, interpreted as follows:

        * If a `~matplotlib.ticker.Locator` instance already, the input
          argument is simply returned.
        * If a list of numbers, these points are ticked. Returns a
          `~matplotlib.ticker.FixedLocator`.
        * If number, this specifies the *step size* between tick locations.
          Returns a `~matplotlib.ticker.MultipleLocator`.

        Otherwise, `locator` should be a string corresponding to one
        of the "registered" locators (see below table). If `locator` is a
        list or tuple and the first element is a "registered" locator name,
        subsequent elements are passed to the locator class as positional
        arguments.

        ======================  =================================================  ================================================================================================
        Key                     Class                                              Description
        ======================  =================================================  ================================================================================================
        ``'null'``, ``'none'``  `~matplotlib.ticker.NullLocator`                   No ticks
        ``'auto'``              `~matplotlib.ticker.AutoLocator`                   Major ticks at sensible locations
        ``'minor'``             `~matplotlib.ticker.AutoMinorLocator`              Minor ticks at sensible locations
        ``'date'``              `~matplotlib.dates.AutoDateLocator`                Default tick locations for datetime axes
        ``'fixed'``             `~matplotlib.ticker.FixedLocator`                  Ticks at these exact locations
        ``'index'``             `~matplotlib.ticker.IndexLocator`                  Ticks on the non-negative integers
        ``'linear'``            `~matplotlib.ticker.LinearLocator`                 Exactly ``N`` ticks encompassing the axis limits, spaced as ``numpy.linspace(lo, hi, N)``
        ``'log'``               `~matplotlib.ticker.LogLocator`                    Ticks for log-scale axes
        ``'logminor'``          `~matplotlib.ticker.LogLocator` preset             Ticks for log-scale axes on the 1st through 9th multiples of each power of the base
        ``'logit'``             `~matplotlib.ticker.LogitLocator`                  Ticks for logit-scale axes
        ``'logitminor'``        `~matplotlib.ticker.LogitLocator` preset           Ticks for logit-scale axes with ``minor=True`` passed to `~matplotlib.ticker.LogitLocator`
        ``'maxn'``              `~matplotlib.ticker.MaxNLocator`                   No more than ``N`` ticks at sensible locations
        ``'multiple'``          `~matplotlib.ticker.MultipleLocator`               Ticks every ``N`` step away from zero
        ``'symlog'``            `~matplotlib.ticker.SymmetricalLogLocator`         Ticks for symmetrical log-scale axes
        ``'symlogminor'``       `~matplotlib.ticker.SymmetricalLogLocator` preset  Ticks for symmetrical log-scale axes on the 1st through 9th multiples of each power of the base
        ``'theta'``             `~matplotlib.projections.polar.ThetaLocator`       Like the base locator but default locations are every `numpy.pi`/8 radians
        ``'year'``              `~matplotlib.dates.YearLocator`                    Ticks every ``N`` years
        ``'month'``             `~matplotlib.dates.MonthLocator`                   Ticks every ``N`` months
        ``'weekday'``           `~matplotlib.dates.WeekdayLocator`                 Ticks every ``N`` weekdays
        ``'day'``               `~matplotlib.dates.DayLocator`                     Ticks every ``N`` days
        ``'hour'``              `~matplotlib.dates.HourLocator`                    Ticks every ``N`` hours
        ``'minute'``            `~matplotlib.dates.MinuteLocator`                  Ticks every ``N`` minutes
        ``'second'``            `~matplotlib.dates.SecondLocator`                  Ticks every ``N`` seconds
        ``'microsecond'``       `~matplotlib.dates.MicrosecondLocator`             Ticks every ``N`` microseconds
        ======================  =================================================  ================================================================================================

    Other parameters
    ----------------
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
        isinstance(num, Number) for num in locator
    ):
        locator, args = locator[0], (*locator[1:], *args)

    # Get the locator
    if isinstance(locator, str):  # dictionary lookup
        # Shorthands and defaults
        if locator in ('logminor', 'logitminor', 'symlogminor'):
            locator, _ = locator.split('minor')
            if locator == 'logit':
                kwargs.setdefault('minor', True)
            else:
                kwargs.setdefault('subs', np.arange(1, 10))
        elif locator == 'index':
            args = args or (1,)
            if len(args) == 1:
                args = (*args, 0)
        # Lookup
        if locator not in locators:
            raise ValueError(
                f'Unknown locator {locator!r}. Options are '
                + ', '.join(map(repr, locators.keys())) + '.'
            )
        locator = locators[locator](*args, **kwargs)
    elif isinstance(locator, Number):  # scalar variable
        locator = mticker.MultipleLocator(locator, *args, **kwargs)
    elif np.iterable(locator):
        locator = mticker.FixedLocator(np.sort(locator), *args, **kwargs)
    else:
        raise ValueError(f'Invalid locator {locator!r}.')
    return locator

from .internals import _not_none

def Formatter(formatter, *args, date=False, index=False, **kwargs):
    """
    Return a `~matplotlib.ticker.Formatter` instance. This function is used to
    interpret the `xformatter`, `xformatter_kw`, `yformatter`, and
    `yformatter_kw` arguments when passed to
    `~proplot.axes.XYAxes.format`, and the `formatter`
    and `formatter_kw` arguments when passed to colorbar methods wrapped by
    `~proplot.wrappers.colorbar_wrapper`.

    Parameters
    ----------
    formatter : `~matplotlib.ticker.Formatter`, str, list of str, or function
        The axis formatter specification, interpreted as follows:

        * If a `~matplotlib.ticker.Formatter` instance already, the input
          argument is simply returned.
        * If list of strings, the ticks are labeled with these strings. Returns
          a `~matplotlib.ticker.FixedFormatter` if `index` is ``False`` or an
          `~matplotlib.ticker.IndexFormatter` if `index` is ``True``.
        * If a function, the labels will be generated using this function.
          Returns a `~matplotlib.ticker.FuncFormatter`.
        * If a string containing ``{x}`` or ``{x:...}``, ticks will be
          formatted by calling ``string.format(x=number)``.
        * If a string containing ``'%'`` and `date` is ``False``, ticks will be
          formatted using the C-style ``string % number`` method. See
          `this page <https://docs.python.org/3/library/stdtypes.html#printf-style-string-formatting>`__
          for a review. Returns a `~matplotlib.ticker.FormatStrFormatter`.
        * If a string containing ``'%'`` and `date` is ``True``, *datetime*
          `string % number`` formatting is used. See
          `this page <https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes>`__
          for a review. Returns a `~matplotlib.ticker.DateFormatter`.

        Otherwise, `formatter` should be a string corresponding to one
        of the "registered" formatters (see below table). If `formatter` is
        a list or tuple and the first element is a "registered" formatter
        name, subsequent elements are passed to the formatter class as
        positional arguments.

        ======================  ==============================================  ===================================================================================================================================
        Key                     Class                                           Description
        ======================  ==============================================  ===================================================================================================================================
        ``'null'``, ``'none'``  `~matplotlib.ticker.NullFormatter`              No tick labels
        ``'auto'``              `AutoFormatter`                                 New default tick labels for axes
        ``'simple'``            `SimpleFormatter`                               New default tick labels for e.g. contour labels
        ``'frac'``              `FracFormatter`                                 Rational fractions
        ``'date'``              `~matplotlib.dates.AutoDateFormatter`           Default tick labels for datetime axes
        ``'concise'``           `~matplotlib.dates.ConciseDateFormatter`        More concise date labels introduced in `matplotlib 3.1 <https://matplotlib.org/3.1.0/users/whats_new.html#concisedateformatter>`__
        ``'datestr'``           `~matplotlib.dates.DateFormatter`               Date formatting with C-style ``string % format`` notation
        ``'eng'``               `~matplotlib.ticker.EngFormatter`               Engineering notation
        ``'fixed'``             `~matplotlib.ticker.FixedFormatter`             List of strings
        ``'formatstr'``         `~matplotlib.ticker.FormatStrFormatter`         From C-style ``string % format`` notation
        ``'func'``              `~matplotlib.ticker.FuncFormatter`              Use an arbitrary function
        ``'index'``             `~matplotlib.ticker.IndexFormatter`             List of strings corresponding to non-negative integer positions along the axis
        ``'log'``, ``'sci'``    `~matplotlib.ticker.LogFormatterSciNotation`    For log-scale axes with scientific notation
        ``'logit'``             `~matplotlib.ticker.LogitFormatter`             For logistic-scale axes
        ``'math'``              `~matplotlib.ticker.LogFormatterMathtext`       For log-scale axes with math text
        ``'percent'``           `~matplotlib.ticker.PercentFormatter`           Trailing percent sign
        ``'scalar'``            `~matplotlib.ticker.ScalarFormatter`            Old default tick labels for axes
        ``'strmethod'``         `~matplotlib.ticker.StrMethodFormatter`         From the ``string.format`` method
        ``'theta'``             `~matplotlib.projections.polar.ThetaFormatter`  Formats radians as degrees, with a degree symbol
        ``'e'``                 `FracFormatter` preset                          Fractions of *e*
        ``'pi'``                `FracFormatter` preset                          Fractions of :math:`\\pi`
        ``'deg'``               `AutoFormatter` preset                          Trailing degree symbol
        ``'deglat'``            `AutoFormatter` preset                          Trailing degree symbol and cardinal "SN" indicator
        ``'deglon'``            `AutoFormatter` preset                          Trailing degree symbol and cardinal "WE" indicator
        ``'lat'``               `AutoFormatter` preset                          Cardinal "SN" indicator
        ``'lon'``               `AutoFormatter` preset                          Cardinal "WE" indicator
        ======================  ==============================================  ===================================================================================================================================

    date : bool, optional
        Toggles the behavior when `formatter` contains a ``'%'`` sign (see
        above).
    index : bool, optional
        Controls the behavior when `formatter` is a list of strings (see
        above).

    Other parameters
    ----------------
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
        isinstance(item, str) for item in formatter
    ):
        formatter, args = formatter[0], (*formatter[1:], *args)

    # Get the formatter
    if isinstance(formatter, str):  # assumption is list of strings
        # Format strings
        if re.search(r'{x(:.+)?}', formatter):
            # string.format() formatting
            formatter = mticker.StrMethodFormatter(
                formatter, *args, **kwargs
            )
        elif '%' in formatter:
            # %-style formatting
            if date:
                formatter = mdates.DateFormatter(
                    formatter, *args, **kwargs
                )
            else:
                formatter = mticker.FormatStrFormatter(
                    formatter, *args, **kwargs
                )
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
                formatter = 'auto'
            # Lookup
            if formatter not in formatters:
                raise ValueError(
                    f'Unknown formatter {formatter!r}. Options are '
                    + ', '.join(map(repr, formatters.keys())) + '.'
                )
            formatter = formatters[formatter](*args, **kwargs)
    elif callable(formatter):
        # Function
        formatter = mticker.FuncFormatter(formatter, *args, **kwargs)
    elif np.iterable(formatter):
        # List of strings
        if index:
            formatter = mticker.IndexFormatter(formatter)
        else:
            formatter = mticker.FixedFormatter(formatter)
    else:
        raise ValueError(f'Invalid formatter {formatter!r}.')
    return formatter

from .internals import warnings

def Scale(scale, *args, **kwargs):
    """
    Return a `~matplotlib.scale.ScaleBase` instance. This function is used to
    interpret the `xscale`, `xscale_kw`, `yscale`, and `yscale_kw` arguments
    when passed to `~proplot.axes.CartesianAxes.format`.

    Parameters
    ----------
    scale : `~matplotlib.scale.ScaleBase`, str, or (str, ...)
        The axis scale specification. If a `~matplotlib.scale.ScaleBase`
        instance already, the input argument is simply returned. Otherwise,
        `scale` should be a string corresponding to one of the
        "registered" axis scales or axis scale presets (see below table).

        If `scale` is a list or tuple and the first element is a
        "registered" scale name, subsequent elements are passed to the
        scale class as positional arguments.

        =================  =======================  =======================================================================
        Key                Class                    Description
        =================  =======================  =======================================================================
        ``'linear'``       `LinearScale`            Linear
        ``'log'``          `LogScale`               Logarithmic
        ``'symlog'``       `SymmetricalLogScale`    Logarithmic beyond finite space around zero
        ``'logit'``        `LogitScale`             Logistic
        ``'inverse'``      `InverseScale`           Inverse
        ``'function'``     `FuncScale`              Scale from arbitrary forward and backwards functions
        ``'sine'``         `SineLatitudeScale`      Sine function (in degrees)
        ``'mercator'``     `MercatorLatitudeScale`  Mercator latitude function (in degrees)
        ``'exp'``          `ExpScale`               Arbitrary exponential function
        ``'power'``        `PowerScale`             Arbitrary power function
        ``'cutoff'``       `CutoffScale`            Arbitrary piecewise linear transformations
        ``'quadratic'``    `PowerScale` (preset)    Quadratic function
        ``'cubic'``        `PowerScale` (preset)    Cubic function
        ``'quartic'``      `PowerScale` (preset)    Cubic function
        ``'db'``           `ExpScale` (preset)      Ratio expressed as `decibels <https://en.wikipedia.org/wiki/Decibel>`__
        ``'np'``           `ExpScale` (preset)      Ratio expressed as `nepers <https://en.wikipedia.org/wiki/Neper>`__
        ``'idb'``          `ExpScale` (preset)      `Decibels <https://en.wikipedia.org/wiki/Decibel>`__ expressed as ratio
        ``'inp'``          `ExpScale` (preset)      `Nepers <https://en.wikipedia.org/wiki/Neper>`__ expressed as ratio
        ``'pressure'``     `ExpScale` (preset)      Height (in km) expressed linear in pressure
        ``'height'``       `ExpScale` (preset)      Pressure (in hPa) expressed linear in height
        =================  =======================  =======================================================================

    Other parameters
    ----------------
    *args, **kwargs
        Passed to the `~matplotlib.scale.ScaleBase` class.

    Returns
    -------
    `~matplotlib.scale.ScaleBase`
        The scale instance.
    """  # noqa
    # NOTE: Why not try to interpret FuncScale arguments, like when lists
    # of numbers are passed to Locator? Because FuncScale *itself* accepts
    # ScaleBase classes as arguments... but constructor functions cannot
    # do anything but return the class instance upon receiving one.
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
            _warn_proplot(
                f'Scale {scale!r} is a scale *preset*. Ignoring positional '
                'argument(s): {args} and keyword argument(s): {kwargs}. '
            )
        scale, *args = SCALE_PRESETS[scale]

    # Get scale
    scale = scale.lower()
    if scale in scales:
        scale = scales[scale]
    else:
        raise ValueError(
            f'Unknown scale or preset {scale!r}. Options are '
            + ', '.join(map(repr, list(scales) + list(SCALE_PRESETS))) + '.'
        )
    return scale(*args, **kwargs)


def _sanitize(string, zerotrim=False):
    """
    Sanitize tick label strings.
    """
    if zerotrim and '.' in string:
        string = string.rstrip('0').rstrip('.')
    string = string.replace('-', '\N{MINUS SIGN}')
    if string == '\N{MINUS SIGN}0':
        string = '0'
    return string


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
    def __init__(
        self, *args,
        zerotrim=None, tickrange=None,
        prefix=None, suffix=None, negpos=None, **kwargs
    ):
        """
        Parameters
        ----------
        zerotrim : bool, optional
            Whether to trim trailing zeros.
            Default is :rc:`axes.formatter.zerotrim`.
        tickrange : (float, float), optional
            Range within which major tick marks are labelled.
        prefix, suffix : str, optional
            Prefix and suffix for all strings.
        negpos : str, optional
            Length-2 string indicating the suffix for "negative" and "positive"
            numbers, meant to replace the minus sign. This is useful for
            indicating cardinal geographic coordinates.

        Other parameters
        ----------------
        *args, **kwargs
            Passed to `~matplotlib.ticker.ScalarFormatter`.

        Warning
        -------
        The matplotlib `~matplotlib.ticker.ScalarFormatter` determines the
        number of significant digits based on the axis limits, and therefore
        may *truncate* digits while formatting ticks on highly non-linear
        axis scales like `~proplot.scale.LogScale`. We try to correct
        this behavior with a patch.
        """
        tickrange = tickrange or (-np.inf, np.inf)
        super().__init__(*args, **kwargs)
        from .config import rc
        zerotrim = _not_none(zerotrim, rc['axes.formatter.zerotrim'])
        self._zerotrim = zerotrim
        self._tickrange = tickrange
        self._prefix = prefix or ''
        self._suffix = suffix or ''
        self._negpos = negpos or ''

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

        # Negative positive handling
        if not self._negpos or x == 0:
            tail = ''
        elif x > 0:
            tail = self._negpos[1]
        else:
            x *= -1
            tail = self._negpos[0]

        # Default string formatting
        string = super().__call__(x, pos)
        string = _sanitize_label(string, zerotrim=self._zerotrim)

        # Add just enough precision for small numbers. Default formatter is
        # only meant to be used for linear scales and cannot handle the wide
        # range of magnitudes in e.g. log scales. To correct this, we only
        # truncate if value is within one order of magnitude of the float
        # precision. Common issue is e.g. levels=plot.arange(-1, 1, 0.1).
        # This choice satisfies even 1000 additions of 0.1 to -100.
        # Example code:
        # def add(x, decimals=1, type_=np.float64):
        #     step = type_(10 ** -decimals)
        #     y = type_(x) + step
        #     if np.round(y, decimals) == 0:
        #         return y
        #     else:
        #         return add(y, decimals)
        # num = abs(add(-200, 1, float))
        # precision = abs(np.log10(num) // 1) - 1
        # ('{:.%df}' % precision).format(num)
        if string == '0' and x != 0:
            string = (
                '{:.%df}' % min(
                    int(abs(np.log10(abs(x)) // 1)),
                    np.finfo(type(x)).precision - 1
                )
            ).format(x)
            string = _sanitize_label(string, zerotrim=self._zerotrim)

        # Prefix and suffix
        sign = ''
        if string and string[0] == '\N{MINUS SIGN}':
            sign, string = string[0], string[1:]
        return sign + self._prefix + string + self._suffix + tail


def SimpleFormatter(precision=6, zerotrim=True):
    """
    Return a `~matplotlib.ticker.FuncFormatter` that replicates the
    `zerotrim` feature from `AutoFormatter`. This is more suitable for
    arbitrary number formatting not necessarily associated with any
    `~matplotlib.axis.Axis` instance, e.g. labeling contours.

    Parameters
    ----------
    precision : int, optional
        The maximum number of digits after the decimal point.
    zerotrim : bool, optional
        Whether to trim trailing zeros.
        Default is :rc:`axes.formatter.zerotrim`.
    """
    from .config import rc
    zerotrim = _not_none(zerotrim, rc['axes.formatter.zerotrim'])

    def func(x, pos):
        string = ('{:.%df}' % precision).format(x)
        if zerotrim and '.' in string:
            string = string.rstrip('0').rstrip('.')
        string = string.replace('-', '\N{MINUS SIGN}')
        if string == '\N{MINUS SIGN}0':
            string = '0'
        return string
    return mticker.FuncFormatter(func)


def FracFormatter(symbol='', number=1):
    r"""
    Return a `~matplotlib.ticker.FuncFormatter` that formats numbers as
    fractions or multiples of some arbitrary value.
    This is powered by the builtin `~fractions.Fraction` class
    and the `~fractions.Fraction.limit_denominator` method.

    Parameters
    ----------
    symbol : str
        The symbol, e.g. ``r'$\pi$'``. Default is ``''``.
    number : float
        The value, e.g. `numpy.pi`. Default is ``1``.
    """
    def func(x, pos):  # must accept location argument
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
    return mticker.FuncFormatter(func)

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
    'func': mticker.FuncFormatter,
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
#: The registered scale names and their associated
#: `~matplotlib.scale.ScaleBase` classes. See `Scale` for a table.
scales = mscale._scale_mapping
