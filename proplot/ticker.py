#!/usr/bin/env python3
"""
Various `~matplotlib.ticker.Locator` and `~matplotlib.ticker.Formatter` classes.
"""
import locale
import re
from fractions import Fraction

import matplotlib.ticker as mticker
import numpy as np

from .config import rc
from .internals import ic  # noqa: F401
from .internals import _not_none, context, docstring

try:
    import cartopy.crs as ccrs
    from cartopy.mpl.ticker import (
        _PlateCarreeFormatter, LatitudeFormatter, LongitudeFormatter
    )
except ModuleNotFoundError:
    ccrs = None
    _PlateCarreeFormatter = LatitudeFormatter = LongitudeFormatter = object

# NOTE: Keep IndexFormatter out of __all__ since we don't want it documented
# on website. Just represents a matplotlib replacement. However *do* keep
# it public so people can access it from module like all other classes.
__all__ = [
    'AutoFormatter',
    'FracFormatter',
    'SciFormatter',
    'SigFigFormatter',
    'SimpleFormatter',
    'DegreeFormatter',
    'LongitudeFormatter',
    'LatitudeFormatter',
    'DegreeLocator',
    'LongitudeLocator',
    'LatitudeLocator',
]

REGEX_ZERO = re.compile('\\A[-\N{MINUS SIGN}]?0(.0*)?\\Z')
REGEX_MINUS = re.compile('\\A[-\N{MINUS SIGN}]\\Z')
REGEX_MINUS_ZERO = re.compile('\\A[-\N{MINUS SIGN}]0(.0*)?\\Z')

_precision_docstring = """
precision : int, optional
    The maximum number of digits after the decimal point. Default is ``6``
    when `zerotrim` is ``True`` and ``2`` otherwise.
"""
_zerotrim_docstring = """
zerotrim : bool, optional
    Whether to trim trailing decimal zeros.
    Default is :rc:`formatter.zerotrim`.
"""
_auto_docstring = """
tickrange : 2-tuple of float, optional
    Range within which major tick marks are labelled. Default is
    ``(-np.inf, np.inf)``.
wraprange : 2-tuple of float, optional
    Range outside of which tick values are wrapped. For example,
    ``(-180, 180)`` will format a value of ``200`` as ``-160``.
prefix, suffix : str, optional
    Prefix and suffix for all tick strings. The suffix is added before
    the optional `negpos` suffix.
negpos : str, optional
    Length-2 string indicating the suffix for "negative" and "positive"
    numbers, meant to replace the minus sign.
"""
_formatter_call = """
Convert number to a string.

Parameters
----------
x : float
    The value.
pos : float, optional
    The position.
"""
docstring._snippet_manager['ticker.precision'] = _precision_docstring
docstring._snippet_manager['ticker.zerotrim'] = _zerotrim_docstring
docstring._snippet_manager['ticker.auto'] = _auto_docstring
docstring._snippet_manager['ticker.call'] = _formatter_call

_dms_docstring = """
Parameters
----------
dms : bool, optional
    Locate the ticks on clean degree-minute-second intervals and format the
    ticks with minutes and seconds instead of decimals. Default is ``False``.
"""
docstring._snippet_manager['ticker.dms'] = _dms_docstring


class DegreeLocator(mticker.MaxNLocator):
    """
    A locator for longitude and latitude gridlines. Adapted from cartopy.
    """
    # NOTE: Locator implementation is weird AF. __init__ just calls set_params with all
    # keyword args and fills in missing params with default_params class attribute.
    # Unknown params result in warning instead of error.
    default_params = mticker.MaxNLocator.default_params.copy()
    default_params.update(nbins=8, dms=False)

    @docstring._snippet_manager
    def __init__(self, *args, **kwargs):
        """
        %(ticker.dms)s
        """
        super().__init__(*args, **kwargs)

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


class LongitudeLocator(DegreeLocator):
    """
    A locator for longitude gridlines. Adapted from cartopy.
    """
    @docstring._snippet_manager
    def __init__(self, *args, **kwargs):
        """
        %(ticker.dms)s
        """
        super().__init__(*args, **kwargs)

    def tick_values(self, vmin, vmax):
        # NOTE: Proplot ensures vmin, vmax are always the *actual* longitude range
        # accounting for central longitude position.
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


class LatitudeLocator(DegreeLocator):
    """
    A locator for latitude gridlines. Adapted from cartopy.
    """
    @docstring._snippet_manager
    def __init__(self, *args, **kwargs):
        """
        %(ticker.dms)s
        """
        super().__init__(*args, **kwargs)

    def tick_values(self, vmin, vmax):
        vmin = max(vmin, -90)
        vmax = min(vmax, 90)
        return super().tick_values(vmin, vmax)

    def _guess_steps(self, vmin, vmax):
        vmin = max(vmin, -90)
        vmax = min(vmax, 90)
        super()._guess_steps(vmin, vmax)

    def _raw_ticks(self, vmin, vmax):
        ticks = super()._raw_ticks(vmin, vmax)
        return [t for t in ticks if -90 <= t <= 90]


class _CartopyFormatter(object):
    """
    Mixin class that fixes cartopy formatters.
    """
    # NOTE: Cartopy formatters pre 0.18 required axis, and *always* translated
    # input values from map projection coordinates to Plate Carrée coordinates.
    # After 0.18 you can avoid this behavior by not setting axis but really
    # dislike that inconsistency. Solution is temporarily change projection.
    def __init__(self, *args, **kwargs):
        import cartopy  # noqa: F401 (ensure available)
        super().__init__(*args, **kwargs)

    def __call__(self, value, pos=None):
        ctx = context._empty_context()
        if self.axis is not None:
            ctx = context._state_context(self.axis.axes, projection=ccrs.PlateCarree())
        with ctx:
            return super().__call__(value, pos)


class LongitudeFormatter(_CartopyFormatter, LongitudeFormatter):
    """
    A formatter for longitude gridlines.
    Adapted from `cartopy.mpl.ticker.LongitudeFormatter`.
    """
    @docstring._snippet_manager
    def __init__(self, *args, **kwargs):
        """
        %(ticker.dms)s
        """
        super().__init__(*args, **kwargs)


class LatitudeFormatter(_CartopyFormatter, LatitudeFormatter):
    """
    A formatter for latitude gridlines.
    Adapted from `cartopy.mpl.ticker.LatitudeFormatter`.
    """
    @docstring._snippet_manager
    def __init__(self, *args, **kwargs):
        """
        %(ticker.dms)s
        """
        super().__init__(*args, **kwargs)


class DegreeFormatter(_CartopyFormatter, _PlateCarreeFormatter):
    """
    A formatter for longitude and latitude gridlines. Adapted from cartopy.
    """
    @docstring._snippet_manager
    def __init__(self, *args, **kwargs):
        """
        %(ticker.dms)s
        """
        super().__init__(*args, **kwargs)

    def _apply_transform(self, value, *args, **kwargs):  # noqa: U100
        return value

    def _hemisphere(self, value, *args, **kwargs):  # noqa: U100
        return ''


class IndexFormatter(mticker.Formatter):
    """
    A duplicate of `~matplotlib.ticker.IndexFormatter`.
    """
    # NOTE: This was deprecated in matplotlib 3.3. For details check out
    # https://github.com/matplotlib/matplotlib/issues/16631 and bring some popcorn.
    def __init__(self, labels):
        self.labels = labels
        self.n = len(labels)

    def __call__(self, x, pos=None):  # noqa: U100
        i = int(round(x))
        if i < 0 or i >= self.n:
            return ''
        else:
            return self.labels[i]


def _default_precision_zerotrim(precision=None, zerotrim=None):
    """
    Return the default zerotrim and precision. Shared by several formatters.
    """
    zerotrim = _not_none(zerotrim, rc['formatter.zerotrim'])
    if precision is None:
        precision = 6 if zerotrim else 2
    return precision, zerotrim


class AutoFormatter(mticker.ScalarFormatter):
    """
    The new default tick label formatter.
    """
    @docstring._snippet_manager
    def __init__(
        self,
        zerotrim=None, tickrange=None, wraprange=None,
        prefix=None, suffix=None, negpos=None,
        **kwargs
    ):
        """
        Parameters
        ----------
        %(ticker.zerotrim)s
        %(ticker.auto)s

        Other parameters
        ----------------
        **kwargs
            Passed to `matplotlib.ticker.ScalarFormatter`.

        Note
        ----
        `matplotlib.ticker.ScalarFormatter` determines the number of
        significant digits based on the axis limits, and therefore may
        truncate digits while formatting ticks on highly non-linear axis
        scales like `~proplot.scale.LogScale`. `AutoFormatter` corrects
        this behavior, making it suitable for arbitrary axis scales. We
        therefore use `AutoFormatter` with every axis scale by default.
        """
        tickrange = tickrange or (-np.inf, np.inf)
        super().__init__(**kwargs)
        zerotrim = _not_none(zerotrim, rc['formatter.zerotrim'])
        self._zerotrim = zerotrim
        self._tickrange = tickrange
        self._wraprange = wraprange
        self._prefix = prefix or ''
        self._suffix = suffix or ''
        self._negpos = negpos or ''

    @docstring._snippet_manager
    def __call__(self, x, pos=None):
        """
        %(ticker.call)s
        """
        # Tick range limitation
        x = self._wrap_tick_range(x, self._wraprange)
        if self._outside_tick_range(x, self._tickrange):
            return ''

        # Negative positive handling
        x, tail = self._neg_pos_format(x, self._negpos, wraprange=self._wraprange)

        # Default string formatting
        string = super().__call__(x, pos)

        # Fix issue where non-zero string is formatted as zero
        string = self._fix_small_number(x, string)

        # Custom string formatting
        string = self._minus_format(string)
        if self._zerotrim:
            string = self._trim_trailing_zeros(string, self._get_decimal_point())

        # Prefix and suffix
        string = self._add_prefix_suffix(string, self._prefix, self._suffix)
        string = string + tail  # add negative-positive indicator
        return string

    def get_offset(self):
        """
        Get the offset but *always* use math text.
        """
        with context._state_context(self, _useMathText=True):
            return super().get_offset()

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

    def _fix_small_number(self, x, string, precision_offset=2):
        """
        Fix formatting for non-zero number that gets formatted as zero. The `offset`
        controls the offset from the true floating point precision at which we want
        to limit maximum precision of the string.
        """
        # Add just enough precision for small numbers. Default formatter is
        # only meant to be used for linear scales and cannot handle the wide
        # range of magnitudes in e.g. log scales. To correct this, we only
        # truncate if value is within `offset` order of magnitude of the float
        # precision. Common issue is e.g. levels=pplt.arange(-1, 1, 0.1).
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
            x -= getattr(self, 'offset', 0)  # guard against API change
            x /= 10 ** getattr(self, 'orderOfMagnitude', 0)  # guard against API change
            precision_true = max(0, self._decimal_place(x))
            precision_max = max(0, np.finfo(type(x)).precision - precision_offset)
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
        use_locale = _not_none(use_locale, self.get_useLocale())
        return self._get_default_decimal_point(use_locale)

    @staticmethod
    def _get_default_decimal_point(use_locale=None):
        """
        Get decimal point symbol for current locale. Called externally.
        """
        use_locale = _not_none(use_locale, rc['formatter.use_locale'])
        return locale.localeconv()['decimal_point'] if use_locale else '.'

    @staticmethod
    def _decimal_place(x):
        """
        Return the decimal place of the number (e.g., 100 is -2 and 0.01 is 2).
        """
        if x == 0:
            digits = 0
        else:
            digits = -int(np.log10(abs(x)) // 1)
        return digits

    @staticmethod
    def _minus_format(string):
        """
        Format the minus sign and avoid "negative zero," e.g. ``-0.000``.
        """
        if rc['axes.unicode_minus'] and not rc['text.usetex']:
            string = string.replace('-', '\N{MINUS SIGN}')
        if REGEX_MINUS_ZERO.match(string):
            string = string[1:]
        return string

    @staticmethod
    def _neg_pos_format(x, negpos, wraprange=None):
        """
        Permit suffixes indicators for "negative" and "positive" numbers.
        """
        # NOTE: If input is a symmetric wraprange, the value conceptually has
        # no "sign", so trim tail and format as absolute value.
        if not negpos or x == 0:
            tail = ''
        elif (
            wraprange is not None
            and np.isclose(-wraprange[0], wraprange[1])
            and np.any(np.isclose(x, wraprange))
        ):
            x = abs(x)
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

    @staticmethod
    def _trim_trailing_zeros(string, decimal_point='.'):
        """
        Sanitize tick label strings.
        """
        if decimal_point in string:
            string = string.rstrip('0').rstrip(decimal_point)
        return string

    @staticmethod
    def _wrap_tick_range(x, wraprange):
        """
        Wrap the tick range to within these values.
        """
        if wraprange is None:
            return x
        base = wraprange[0]
        modulus = wraprange[1] - wraprange[0]
        return (x - base) % modulus + base


class SciFormatter(mticker.Formatter):
    """
    Format numbers with scientific notation.
    """
    @docstring._snippet_manager
    def __init__(self, precision=None, zerotrim=None):
        """
        Parameters
        ----------
        %(ticker.precision)s
        %(ticker.zerotrim)s
        """
        precision, zerotrim = _default_precision_zerotrim(precision, zerotrim)
        self._precision = precision
        self._zerotrim = zerotrim

    @docstring._snippet_manager
    def __call__(self, x, pos=None):  # noqa: U100
        """
        %(ticker.call)s
        """
        # Get string
        decimal_point = AutoFormatter._get_default_decimal_point()
        string = ('{:.%de}' % self._precision).format(x)
        parts = string.split('e')

        # Trim trailing zeros
        significand = parts[0].rstrip(decimal_point)
        if self._zerotrim:
            significand = AutoFormatter._trim_trailing_zeros(significand, decimal_point)

        # Get sign and exponent
        sign = parts[1][0].replace('+', '')
        exponent = parts[1][1:].lstrip('0')
        if exponent:
            exponent = f'10^{{{sign}{exponent}}}'
        if significand and exponent:
            string = rf'{significand}{{\times}}{exponent}'
        else:
            string = rf'{significand}{exponent}'

        # Ensure unicode minus sign
        string = AutoFormatter._minus_format(string)

        # Return TeX string
        return f'${string}$'


class SigFigFormatter(mticker.Formatter):
    """
    Format numbers by retaining the specified number of significant digits.
    """
    @docstring._snippet_manager
    def __init__(self, sigfig=None, zerotrim=None, base=None):
        """
        Parameters
        ----------
        sigfig : float, optional
            The number of significant digits. Default is ``3``.
        %(ticker.zerotrim)s
        base : float, optional
            The base unit for rounding. Default is ``1``. For example
            ``SigFigFormatter(2, base=5)`` rounds to the nearest 5 with
            up to 2 significant digits (e.g., 87 --> 85, 8.7 --> 8.5).
        """
        self._sigfig = _not_none(sigfig, 3)
        self._zerotrim = _not_none(zerotrim, rc['formatter.zerotrim'])
        self._base = _not_none(base, 1)

    @docstring._snippet_manager
    def __call__(self, x, pos=None):  # noqa: U100
        """
        %(ticker.call)s
        """
        # Limit to significant figures
        digits = AutoFormatter._decimal_place(x) + self._sigfig - 1
        scale = self._base * 10 ** -digits
        x = scale * round(x / scale)

        # Create the string
        decimal_point = AutoFormatter._get_default_decimal_point()
        precision = max(0, digits) + max(0, AutoFormatter._decimal_place(self._base))
        string = ('{:.%df}' % precision).format(x)
        string = string.replace('.', decimal_point)

        # Custom string formatting
        string = AutoFormatter._minus_format(string)
        if self._zerotrim:
            string = AutoFormatter._trim_trailing_zeros(string, decimal_point)
        return string


class SimpleFormatter(mticker.Formatter):
    """
    A general purpose number formatter. This is similar to `AutoFormatter`
    but suitable for arbitrary formatting not necessarily associated with
    an `~matplotlib.axis.Axis` instance.
    """
    @docstring._snippet_manager
    def __init__(
        self, precision=None, zerotrim=None,
        tickrange=None, wraprange=None,
        prefix=None, suffix=None, negpos=None,
    ):
        """
        Parameters
        ----------
        %(ticker.precision)s
        %(ticker.zerotrim)s
        %(ticker.auto)s
        """
        precision, zerotrim = _default_precision_zerotrim(precision, zerotrim)
        self._precision = precision
        self._prefix = prefix or ''
        self._suffix = suffix or ''
        self._negpos = negpos or ''
        self._tickrange = tickrange or (-np.inf, np.inf)
        self._wraprange = wraprange
        self._zerotrim = zerotrim

    @docstring._snippet_manager
    def __call__(self, x, pos=None):  # noqa: U100
        """
        %(ticker.call)s
        """
        # Tick range limitation
        x = AutoFormatter._wrap_tick_range(x, self._wraprange)
        if AutoFormatter._outside_tick_range(x, self._tickrange):
            return ''

        # Negative positive handling
        x, tail = AutoFormatter._neg_pos_format(
            x, self._negpos, wraprange=self._wraprange
        )

        # Default string formatting
        decimal_point = AutoFormatter._get_default_decimal_point()
        string = ('{:.%df}' % self._precision).format(x)
        string = string.replace('.', decimal_point)

        # Custom string formatting
        string = AutoFormatter._minus_format(string)
        if self._zerotrim:
            string = AutoFormatter._trim_trailing_zeros(string, decimal_point)

        # Prefix and suffix
        string = AutoFormatter._add_prefix_suffix(string, self._prefix, self._suffix)
        string = string + tail  # add negative-positive indicator
        return string


class FracFormatter(mticker.Formatter):
    r"""
    Format numbers as fractions or multiples of some value.
    This is powered by the builtin `~fractions.Fraction` class
    and the `~fractions.Fraction.limit_denominator` method.
    """
    def __init__(self, symbol='', number=1):
        r"""
        Parameters
        ----------
        symbol : str
            The symbol, e.g. ``r'$\pi$'``. Default is ``''``.
        number : float
            The value, e.g. `numpy.pi`. Default is ``1``.
        """
        self._symbol = symbol
        self._number = number
        super().__init__()

    @docstring._snippet_manager
    def __call__(self, x, pos=None):  # noqa: U100
        """
        %(ticker.call)s
        """
        frac = Fraction(x / self._number).limit_denominator()
        symbol = self._symbol
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
