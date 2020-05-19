#!/usr/bin/env python3
"""
Various `~matplotlib.ticker.Locator` and `~matplotlib.ticker.Formatter`
classes.
"""
import re
import numpy as np
import matplotlib.ticker as mticker
import locale
from fractions import Fraction
from .internals import ic  # noqa: F401
from .internals import docstring, _not_none

__all__ = [
    'AutoFormatter',
    'FracFormatter',
    'SigFigFormatter',
    'SimpleFormatter',
]

REGEX_ZERO = re.compile('\\A[-\N{MINUS SIGN}]?0(.0*)?\\Z')
REGEX_MINUS = re.compile('\\A[-\N{MINUS SIGN}]\\Z')
REGEX_MINUS_ZERO = re.compile('\\A[-\N{MINUS SIGN}]0(.0*)?\\Z')

docstring.snippets['formatter.params'] = """
zerotrim : bool, optional
    Whether to trim trailing zeros. Default is :rc:`axes.formatter.zerotrim`.
tickrange : (float, float), optional
    Range within which major tick marks are labelled. Default is ``(-np.inf, np.inf)``.
prefix, suffix : str, optional
    Prefix and suffix for all strings.
negpos : str, optional
    Length-2 string indicating the suffix for "negative" and "positive"
    numbers, meant to replace the minus sign. This is useful for indicating
    cardinal geographic coordinates.
"""


class _GeoLocator(mticker.MaxNLocator):
    """
    A locator for determining longitude and latitude gridlines.
    Adapted from cartopy.
    """
    # NOTE: Locator implementation is weird AF. __init__ just calls set_params with all
    # keyword args and fills in missing params with default_params class attribute.
    # Unknown params result in warning instead of error.
    default_params = mticker.MaxNLocator.default_params.copy()
    default_params.update(nbins=8, dms=False)

    def set_params(self, **kwargs):
        if 'dms' in kwargs:
            self._dms = kwargs.pop('dms')
        super().set_params(**kwargs)

    def _guess_steps(self, vmin, vmax):
        dv = abs(vmax - vmin)
        if dv > 180:
            dv -= 180
        if dv > 50:
            steps = np.array([1, 2, 3, 6, 10])
        elif not self._dms or dv > 3.0:
            steps = np.array([1, 1.5, 2, 2.5, 3, 5, 10])
        else:
            steps = np.array([1, 10 / 6.0, 15 / 6.0, 20 / 6.0, 30 / 6.0, 10])
        self.set_params(steps=np.array(steps))

    def _raw_ticks(self, vmin, vmax):
        self._guess_steps(vmin, vmax)
        return super()._raw_ticks(vmin, vmax)

    def bin_boundaries(self, vmin, vmax):  # matplotlib <2.2.0
        return self._raw_ticks(vmin, vmax)  # may call Latitude/LongitudeLocator copies


class _LongitudeLocator(mticker.MaxNLocator):
    """
    A locator for determining longitude gridlines. Includes considerations
    to prevent double gridlines on tick boundaries.

    Parameters
    ----------
    lonstep : float, optional
        The longitude step size. If ``None`` this is determined automatically.
    dms : bool, optional
        Allow the locator to stop on minutes and seconds. Default is ``False``.
    """
    # NOTE: Proplot ensures vmin, vmax are always the *actual* longitude range
    # accounting for central longitude position.
    def tick_values(self, vmin, vmax):
        ticks = super().tick_values(vmin, vmax)
        if np.isclose(ticks[0] + 360, ticks[-1]):
            eps = 1e-10
            if ticks[-1] % 360 > 0:
                # Make sure the label appears on *right*, not on
                # top of the leftmost label.
                ticks[-1] -= eps
            else:
                # Formatter formats label as 1e-10... so there is simply no way to
                # put label on right. Just shift this location off the map edge so
                # parallels still extend all the way to the edge, but label disappears.
                ticks[-1] += eps
        return ticks


class _LatitudeLocator(_GeoLocator):
    """
    A locator for latitudes that works even at very small scale.
    Copied  from cartopy so this can be used in basemap.

    Parameters
    ----------
    latstep : float, optional
        The latitude step size. If ``None`` this is determined automatically.
    latmax : float, optional
        Maximum absolute gridline latitude. Lines poleward of this latitude will be
        cut off. This can be set to e.g. ``80`` to create a "hole" around the poles
        for pole-centered projections.
    dms : bool, optional
        Allow the locator to stop on minutes and seconds. Default is ``False``.
    """
    default_params = _GeoLocator.default_params.copy()
    default_params.update(symmetric=True)  # always draw gridline on equator!

    def tick_values(self, vmin, vmax):
        latmax = self.axis.get_latmax()
        vmin = max(vmin, -latmax)
        vmax = min(vmax, latmax)
        return super().tick_values(vmin, vmax)

    def _guess_steps(self, vmin, vmax):
        latmax = self.axis.get_latmax()
        vmin = max(vmin, -latmax)
        vmax = min(vmax, latmax)
        super()._guess_steps(vmin, vmax)

    def _raw_ticks(self, vmin, vmax):
        latmax = self.axis.get_latmax()
        ticks = super()._raw_ticks(vmin, vmax)
        return [t for t in ticks if -latmax <= t <= latmax]


class AutoFormatter(mticker.ScalarFormatter):
    """
    The new default formatter. Differs from `~matplotlib.ticker.ScalarFormatter`
    in the following ways:

    1. Trims trailing decimal zeros by default.
    2. Permits specifying *range* within which major tick marks are labeled.
    3. Permits adding arbitrary prefix or suffix to every tick label string.
    4. Permits adding "negative" and "positive" indicator.
    """
    @docstring.add_snippets
    def __init__(
        self,
        zerotrim=None, tickrange=None,
        prefix=None, suffix=None, negpos=None,
        **kwargs
    ):
        """
        Parameters
        ----------
        %(formatter.params)s

        Other parameters
        ----------------
        **kwargs
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
        super().__init__(**kwargs)
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
        if self._outside_tick_range(x, self._tickrange):
            return ''

        # Negative positive handling
        x, tail = self._neg_pos_format(x, self._negpos)

        # Default string formatting
        string = super().__call__(x, pos)

        # Fix issue where non-zero string is formatted as zero
        string = self._fix_small_number(x, string)

        # Custom string formatting
        string = self._minus_format(string)
        if self._zerotrim:
            string = self._trim_trailing_zeros(string)

        # Prefix and suffix
        string = self._add_prefix_suffix(string, self._prefix, self._suffix)
        string = string + tail  # add negative-positive indicator
        return string

    @staticmethod
    def _add_prefix_suffix(string, prefix=None, suffix=None):
        """
        Add prefix and suffix to string.
        """
        sign = ''
        prefix = prefix or ''
        suffix = suffix or ''
        if string and REGEX_MINUS.match(string[0]):
            sign, string = string[0], string[1:]
        return sign + prefix + string + suffix

    def _fix_small_number(self, x, string, offset=2):
        """
        Fix formatting for non-zero number that gets formatted as zero. The `offset`
        controls the offset from the true floating point precision at which we want
        to limit maximum precision of the string.
        """
        # Add just enough precision for small numbers. Default formatter is
        # only meant to be used for linear scales and cannot handle the wide
        # range of magnitudes in e.g. log scales. To correct this, we only
        # truncate if value is within `offset` order of magnitude of the float
        # precision. Common issue is e.g. levels=plot.arange(-1, 1, 0.1).
        # This choice satisfies even 1000 additions of 0.1 to -100.
        match = REGEX_ZERO.match(string)
        decimal_point = self._get_decimal_point()

        if match and x != 0:
            # Get initial precision spit out by algorithm
            decimals, = match.groups()
            if decimals:
                precision_init = len(decimals.lstrip(decimal_point))
            else:
                precision_init = 0

            # Format with precision below floating point error
            precision_true = int(abs(np.log10(abs(x)) // 1))
            precision_max = np.finfo(type(x)).precision - offset
            precision = min(precision_true, precision_max)
            string = ('{:.%df}' % precision).format(x)

            # If number is zero after ignoring floating point error, generate
            # zero with precision matching original string.
            if REGEX_ZERO.match(string):
                string = ('{:.%df}' % precision_init).format(0)

            # Fix decimal point
            string = string.replace('.', decimal_point)

        return string

    def _get_decimal_point(self, use_locale=None):
        """
        Get decimal point symbol for current locale (e.g. in Europe will be comma).
        """
        from .config import rc
        use_locale = _not_none(
            use_locale, self.get_useLocale(), rc['axes.formatter.use_locale']
        )
        return locale.localeconv()['decimal_point'] if use_locale else '.'

    @staticmethod
    def _get_default_decimal_point(use_locale=None):
        """
        Get decimal point symbol for current locale. Called externally.
        """
        from .config import rc
        use_locale = _not_none(use_locale, rc['axes.formatter.use_locale'])
        return locale.localeconv()['decimal_point'] if use_locale else '.'

    @staticmethod
    def _minus_format(string):
        """
        Format the minus sign and avoid "negative zero," e.g. ``-0.000``.
        """
        from .config import rc
        if rc['axes.unicode_minus'] and not rc['text.usetex']:
            string = string.replace('-', '\N{MINUS SIGN}')
        if REGEX_MINUS_ZERO.match(string):
            string = string[1:]
        return string

    @staticmethod
    def _neg_pos_format(x, negpos):
        """
        Permit suffixes indicators for "negative" and "positive" numbers.
        """
        if not negpos or x == 0:
            tail = ''
        elif x > 0:
            tail = negpos[1]
        else:
            x *= -1
            tail = negpos[0]
        return x, tail

    @staticmethod
    def _outside_tick_range(x, tickrange):
        """
        Return whether point is outside tick range up to some precision.
        """
        eps = abs(x) / 1000
        return (x + eps) < tickrange[0] or (x - eps) > tickrange[1]

    def _trim_trailing_zeros(self, string):
        """
        Sanitize tick label strings.
        """
        decimal_point = self._get_decimal_point()
        if decimal_point in string:
            string = string.rstrip('0').rstrip(decimal_point)
        return string


def SigFigFormatter(sigfig=1, zerotrim=None):
    """
    Return a `~matplotlib.ticker.FuncFormatter` that rounds numbers
    to the specified number of *significant digits*.

    Parameters
    ----------
    sigfig : float, optional
        The number of significant digits.
    zerotrim : bool, optional
        Whether to trim trailing zeros.
    """
    from .config import rc
    zerotrim = _not_none(zerotrim, rc['axes.formatter.zerotrim'])

    def func(x, pos):
        # Limit digits to significant figures
        if x == 0:
            digits = 0
        else:
            digits = -int(np.log10(abs(x)) // 1)
        digits += sigfig - 1
        x = np.round(x, digits)
        string = ('{:.%df}' % max(0, digits)).format(x)
        string = string.replace('.', AutoFormatter._get_default_decimal_point())

        # Custom string formatting
        string = AutoFormatter._minus_format(string)
        if zerotrim:
            string = AutoFormatter._trim_trailing_zeros(string)
        return string

    return mticker.FuncFormatter(func)


@docstring.add_snippets
def SimpleFormatter(
    precision=6, zerotrim=None, tickrange=None,
    prefix=None, suffix=None, negpos=None,
):
    """
    Return a `~matplotlib.ticker.FuncFormatter` that replicates the
    `zerotrim` feature from `AutoFormatter`. This is more suitable for
    arbitrary number formatting not necessarily associated with any
    `~matplotlib.axis.Axis` instance, e.g. labeling contours.

    Parameters
    ----------
    precision : int, optional
        The maximum number of digits after the decimal point. Default is ``6``.
    %(formatter.params)s
    """
    from .config import rc
    zerotrim = _not_none(zerotrim, rc['axes.formatter.zerotrim'])
    tickrange = tickrange or (-np.inf, np.inf)
    prefix = prefix or ''
    suffix = suffix or ''
    negpos = negpos or ''

    def func(x, pos):
        # Tick range limitation
        if AutoFormatter._outside_tick_range(x, tickrange):
            return ''

        # Negative positive handling
        x, tail = AutoFormatter._neg_pos_format(x, negpos)

        # Default string formatting
        string = ('{:.%df}' % precision).format(x)
        string = string.replace('.', AutoFormatter._get_default_decimal_point())

        # Custom string formatting
        string = AutoFormatter._minus_format(string)
        if zerotrim:
            string = AutoFormatter._trim_trailing_zeros(string)

        # Prefix and suffix
        string = AutoFormatter._add_prefix_suffix(string, prefix, suffix)
        string = string + tail  # add negative-positive indicator
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
        string = AutoFormatter._minus_format(string)
        return string

    return mticker.FuncFormatter(func)
