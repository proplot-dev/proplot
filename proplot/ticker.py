#!/usr/bin/env python3
"""
Various `~matplotlib.ticker.Locator` and `~matplotlib.ticker.Formatter` classes.
"""
import locale
import re
from fractions import Fraction

import matplotlib.axis as maxis
import matplotlib.ticker as mticker
import numpy as np

from .config import rc
from .internals import ic  # noqa: F401
from .internals import _not_none, context, docstring

try:
    import cartopy.crs as ccrs
    from cartopy.mpl.ticker import (
        LatitudeFormatter,
        LongitudeFormatter,
        _PlateCarreeFormatter,
    )
except ModuleNotFoundError:
    ccrs = None
    LatitudeFormatter = LongitudeFormatter = _PlateCarreeFormatter = object

__all__ = [
    "IndexLocator",
    "DiscreteLocator",
    "DegreeLocator",
    "LongitudeLocator",
    "LatitudeLocator",
    "AutoFormatter",
    "SimpleFormatter",
    "IndexFormatter",
    "SciFormatter",
    "SigFigFormatter",
    "FracFormatter",
    "DegreeFormatter",
    "LongitudeFormatter",
    "LatitudeFormatter",
]

REGEX_ZERO = re.compile("\\A[-\N{MINUS SIGN}]?0(.0*)?\\Z")
REGEX_MINUS = re.compile("\\A[-\N{MINUS SIGN}]\\Z")
REGEX_MINUS_ZERO = re.compile("\\A[-\N{MINUS SIGN}]0(.0*)?\\Z")

_precision_docstring = """
precision : int, default: {6, 2}
    The maximum number of digits after the decimal point. Default is ``6``
    when `zerotrim` is ``True`` and ``2`` otherwise.
"""
_zerotrim_docstring = """
zerotrim : bool, default: :rc:`format.zerotrim`
    Whether to trim trailing decimal zeros.
"""
_auto_docstring = """
tickrange : 2-tuple of float, optional
    Range within which major tick marks are labeled.
    All ticks are labeled by default.
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
docstring._snippet_manager["ticker.precision"] = _precision_docstring
docstring._snippet_manager["ticker.zerotrim"] = _zerotrim_docstring
docstring._snippet_manager["ticker.auto"] = _auto_docstring
docstring._snippet_manager["ticker.call"] = _formatter_call

_dms_docstring = """
Parameters
----------
dms : bool, default: False
    Locate the ticks on clean degree-minute-second intervals and format the
    ticks with minutes and seconds instead of decimals.
"""
docstring._snippet_manager["ticker.dms"] = _dms_docstring


def _default_precision_zerotrim(precision=None, zerotrim=None):
    """
    Return the default zerotrim and precision. Shared by several formatters.
    """
    zerotrim = _not_none(zerotrim, rc["formatter.zerotrim"])
    if precision is None:
        precision = 6 if zerotrim else 2
    return precision, zerotrim


class IndexLocator(mticker.Locator):
    """
    Format numbers by assigning fixed strings to non-negative indices. The ticks
    are restricted to the extent of plotted content when content is present.
    """

    def __init__(self, base=1, offset=0):
        self._base = base
        self._offset = offset

    def set_params(self, base=None, offset=None):
        if base is not None:
            self._base = base
        if offset is not None:
            self._offset = offset

    def __call__(self):
        # NOTE: We adapt matplotlib IndexLocator to support case where
        # the data interval is empty. Only restrict after data is plotted.
        dmin, dmax = self.axis.get_data_interval()
        vmin, vmax = self.axis.get_view_interval()
        min_ = max(dmin, vmin)
        max_ = min(dmax, vmax)
        return self.tick_values(min_, max_)

    def tick_values(self, vmin, vmax):
        base, offset = self._base, self._offset
        vmin = max(base * np.ceil(vmin / base), offset)
        vmax = max(base * np.floor(vmax / base), offset)
        locs = np.arange(vmin, vmax + 0.5 * base, base)
        return self.raise_if_exceeds(locs)


class DiscreteLocator(mticker.Locator):
    """
    A tick locator suitable for discretized colorbars. Adds ticks to some
    subset of the location list depending on the available space determined from
    `~matplotlib.axis.Axis.get_tick_space`. Zero will be used if it appears in the
    location list, and step sizes along the location list are restricted to "nice"
    intervals by default.
    """

    default_params = {
        "nbins": None,
        "minor": False,
        "steps": np.array([1, 2, 3, 4, 5, 6, 8, 10]),
        "min_n_ticks": 2,
    }

    @docstring._snippet_manager
    def __init__(self, locs, **kwargs):
        """
        Parameters
        ----------
        locs : array-like
            The tick location list.
        nbins : int, optional
            Maximum number of ticks to select. By default this is automatically
            determined based on the the axis length and tick label font size.
        minor : bool, default: False
            Whether this is for "minor" ticks. Setting to ``True`` will select more
            ticks with an index step that divides the index step used for "major" ticks.
        steps : array-like of int, default: ``[1 2 3 4 5 6 8]``
            Valid integer index steps when selecting from the tick list. Must fall
            between 1 and 9. Powers of 10 of these step sizes will also be permitted.
        min_n_ticks : int, default: 1
            The minimum number of ticks to select. See also `nbins`.
        """
        self.locs = np.array(locs)
        self._nbins = None  # otherwise unset
        self.set_params(**{**self.default_params, **kwargs})

    def __call__(self):
        """
        Return the locations of the ticks.
        """
        return self.tick_values(None, None)

    def set_params(self, steps=None, nbins=None, minor=None, min_n_ticks=None):
        """
        Set the parameters for this locator. See `DiscreteLocator` for details.
        """
        if steps is not None:
            steps = np.unique(np.array(steps, dtype=int))  # also sorts, makes 1D
            if np.any(steps < 1) or np.any(steps > 10):
                raise ValueError("Steps must fall between one and ten (inclusive).")
            if steps[0] != 1:
                steps = np.concatenate([[1], steps])
            if steps[-1] != 10:
                steps = np.concatenate([steps, [10]])
            self._steps = steps
        if nbins is not None:
            self._nbins = nbins
        if minor is not None:
            self._minor = bool(minor)  # needed to scale tick space
        if min_n_ticks is not None:
            self._min_n_ticks = int(min_n_ticks)  # compare to MaxNLocator

    def tick_values(self, vmin, vmax):  # noqa: U100
        """
        Return the locations of the ticks.
        """
        # NOTE: Critical that minor tick interval evenly divides major tick
        # interval. Otherwise get misaligned major and minor tick steps.
        # NOTE: This tries to select ticks that are integer steps away from zero (like
        # AutoLocator). The list minimum is used if this fails (like FixedLocator)
        # NOTE: This avoids awkward steps like '7' or '13' that produce awkward
        # jumps and have no integer divisors (and therefore eliminate minor ticks)
        # NOTE: We virtually always want to subsample the level list rather than
        # using continuous minor locators (e.g. LogLocator or SymLogLocator) because
        # _parse_autolev interpolates evenly in the norm-space (e.g. 1, 3.16, 10, 31.6
        # for a LogNorm) rather than in linear-space (e.g. 1, 5, 10, 15, 20).
        locs = self.locs
        axis = self.axis
        if axis is None:
            return locs
        nbins = self._nbins
        steps = self._steps
        if nbins is None:
            nbins = axis.get_tick_space()
            nbins = max((1, self._min_n_ticks - 1, nbins))
        step = max(1, int(np.ceil(locs.size / nbins)))
        fact = 10 ** max(0, -AutoFormatter._decimal_place(step))  # e.g. 2 for 100
        idx = min(len(steps) - 1, np.searchsorted(steps, step / fact))
        step = int(np.round(steps[idx] * fact))
        if self._minor:  # tick every half font size
            if isinstance(axis, maxis.XAxis):
                fact = 6  # unscale heuristic scaling of 3 em-widths
            elif isinstance(axis, maxis.YAxis):
                fact = 4  # unscale standard scaling of 2 em-widths
            else:
                fact = 2  # fall back to just one em-width
            for i in range(fact, 0, -1):
                if step % i == 0:
                    step = step // i
                    break
        diff = np.abs(np.diff(locs[: step + 1 : step]))
        (offset,) = np.where(np.isclose(locs % diff if diff.size else 0.0, 0.0))
        offset = offset[0] if offset.size else np.argmin(np.abs(locs))
        return locs[offset % step :: step]  # even multiples from zero or zero-close


class DegreeLocator(mticker.MaxNLocator):
    """
    Locate geographic gridlines with degree-minute-second support.
    Adapted from cartopy.
    """

    # NOTE: This is identical to cartopy except they only define LongitutdeLocator
    # for common methods whereas we use DegreeLocator. More intuitive this way in
    # case users need degree-minute-seconds for non-specific degree axis.
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
        if "dms" in kwargs:
            self._dms = kwargs.pop("dms")
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

    def bin_boundaries(self, vmin, vmax):  # matplotlib < 2.2.0
        return self._raw_ticks(vmin, vmax)  # may call Latitude/Longitude Locator copies


class LongitudeLocator(DegreeLocator):
    """
    Locate longitude gridlines with degree-minute-second support.
    Adapted from cartopy.
    """

    @docstring._snippet_manager
    def __init__(self, *args, **kwargs):
        """
        %(ticker.dms)s
        """
        super().__init__(*args, **kwargs)


class LatitudeLocator(DegreeLocator):
    """
    Locate latitude gridlines with degree-minute-second support.
    Adapted from cartopy.
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


class AutoFormatter(mticker.ScalarFormatter):
    """
    The default formatter used for proplot tick labels.
    Replaces `~matplotlib.ticker.ScalarFormatter`.
    """

    @docstring._snippet_manager
    def __init__(
        self,
        zerotrim=None,
        tickrange=None,
        wraprange=None,
        prefix=None,
        suffix=None,
        negpos=None,
        **kwargs,
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

        See also
        --------
        proplot.constructor.Formatter
        proplot.ticker.SimpleFormatter

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
        zerotrim = _not_none(zerotrim, rc["formatter.zerotrim"])
        self._zerotrim = zerotrim
        self._tickrange = tickrange
        self._wraprange = wraprange
        self._prefix = prefix or ""
        self._suffix = suffix or ""
        self._negpos = negpos or ""

    @docstring._snippet_manager
    def __call__(self, x, pos=None):
        """
        %(ticker.call)s
        """
        # Tick range limitation
        x = self._wrap_tick_range(x, self._wraprange)
        if self._outside_tick_range(x, self._tickrange):
            return ""

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
        sign = ""
        prefix = prefix or ""
        suffix = suffix or ""
        if string and REGEX_MINUS.match(string[0]):
            sign, string = string[0], string[1:]
        return sign + prefix + string + suffix

    def _fix_small_number(self, x, string, precision_offset=2):
        """
        Fix formatting for non-zero formatted as zero. The `offset` controls the offset
        from true floating point precision at which we want to limit string precision.
        """
        # Add just enough precision for small numbers. Default formatter is
        # only meant to be used for linear scales and cannot handle the wide
        # range of magnitudes in e.g. log scales. To correct this, we only
        # truncate if value is within `offset` order of magnitude of the float
        # precision. Common issue is e.g. levels=pplt.arange(-1, 1, 0.1).
        # This choice satisfies even 1000 additions of 0.1 to -100.
        m = REGEX_ZERO.match(string)
        decimal_point = self._get_decimal_point()

        if m and x != 0:
            # Get initial precision spit out by algorithm
            (decimals,) = m.groups()
            precision_init = len(decimals.lstrip(decimal_point)) if decimals else 0

            # Format with precision below floating point error
            x -= getattr(self, "offset", 0)  # guard against API change
            x /= 10 ** getattr(self, "orderOfMagnitude", 0)  # guard against API change
            precision_true = max(0, self._decimal_place(x))
            precision_max = max(0, np.finfo(type(x)).precision - precision_offset)
            precision = min(precision_true, precision_max)
            string = ("{:.%df}" % precision).format(x)

            # If zero ignoring floating point error then match original precision
            if REGEX_ZERO.match(string):
                string = ("{:.%df}" % precision_init).format(0)

            # Fix decimal point
            string = string.replace(".", decimal_point)

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
        use_locale = _not_none(use_locale, rc["formatter.use_locale"])
        return locale.localeconv()["decimal_point"] if use_locale else "."

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
        if rc["axes.unicode_minus"] and not rc["text.usetex"]:
            string = string.replace("-", "\N{MINUS SIGN}")
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
            tail = ""
        elif (
            wraprange is not None
            and np.isclose(-wraprange[0], wraprange[1])
            and np.any(np.isclose(x, wraprange))
        ):
            x = abs(x)
            tail = ""
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
    def _trim_trailing_zeros(string, decimal_point="."):
        """
        Sanitize tick label strings.
        """
        if decimal_point in string:
            string = string.rstrip("0").rstrip(decimal_point)
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


class SimpleFormatter(mticker.Formatter):
    """
    A general purpose number formatter. This is similar to `AutoFormatter`
    but suitable for arbitrary formatting not necessarily associated with
    an `~matplotlib.axis.Axis` instance.
    """

    @docstring._snippet_manager
    def __init__(
        self,
        precision=None,
        zerotrim=None,
        tickrange=None,
        wraprange=None,
        prefix=None,
        suffix=None,
        negpos=None,
    ):
        """
        Parameters
        ----------
        %(ticker.precision)s
        %(ticker.zerotrim)s
        %(ticker.auto)s

        See also
        --------
        proplot.constructor.Formatter
        proplot.ticker.AutoFormatter
        """
        precision, zerotrim = _default_precision_zerotrim(precision, zerotrim)
        self._precision = precision
        self._prefix = prefix or ""
        self._suffix = suffix or ""
        self._negpos = negpos or ""
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
            return ""

        # Negative positive handling
        x, tail = AutoFormatter._neg_pos_format(
            x, self._negpos, wraprange=self._wraprange
        )

        # Default string formatting
        decimal_point = AutoFormatter._get_default_decimal_point()
        string = ("{:.%df}" % self._precision).format(x)
        string = string.replace(".", decimal_point)

        # Custom string formatting
        string = AutoFormatter._minus_format(string)
        if self._zerotrim:
            string = AutoFormatter._trim_trailing_zeros(string, decimal_point)

        # Prefix and suffix
        string = AutoFormatter._add_prefix_suffix(string, self._prefix, self._suffix)
        string = string + tail  # add negative-positive indicator
        return string


class IndexFormatter(mticker.Formatter):
    """
    Format numbers by assigning fixed strings to non-negative indices. Generally
    paired with `IndexLocator` or `~matplotlib.ticker.FixedLocator`.
    """

    # NOTE: This was deprecated in matplotlib 3.3. For details check out
    # https://github.com/matplotlib/matplotlib/issues/16631 and bring some popcorn.
    def __init__(self, labels):
        self.labels = labels
        self.n = len(labels)

    def __call__(self, x, pos=None):  # noqa: U100
        i = int(round(x))
        if i < 0 or i >= self.n:
            return ""
        else:
            return self.labels[i]


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

        See also
        --------
        proplot.constructor.Formatter
        proplot.ticker.AutoFormatter
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
        string = ("{:.%de}" % self._precision).format(x)
        parts = string.split("e")

        # Trim trailing zeros
        significand = parts[0].rstrip(decimal_point)
        if self._zerotrim:
            significand = AutoFormatter._trim_trailing_zeros(significand, decimal_point)

        # Get sign and exponent
        sign = parts[1][0].replace("+", "")
        exponent = parts[1][1:].lstrip("0")
        if exponent:
            exponent = f"10^{{{sign}{exponent}}}"
        if significand and exponent:
            string = rf"{significand}{{\times}}{exponent}"
        else:
            string = rf"{significand}{exponent}"

        # Ensure unicode minus sign
        string = AutoFormatter._minus_format(string)

        # Return TeX string
        return f"${string}$"


class SigFigFormatter(mticker.Formatter):
    """
    Format numbers by retaining the specified number of significant digits.
    """

    @docstring._snippet_manager
    def __init__(self, sigfig=None, zerotrim=None, base=None):
        """
        Parameters
        ----------
        sigfig : float, default: 3
            The number of significant digits.
        %(ticker.zerotrim)s
        base : float, default: 1
            The base unit for rounding. For example ``SigFigFormatter(2, base=5)``
            rounds to the nearest 5 with up to 2 digits (e.g., 87 --> 85, 8.7 --> 8.5).

        See also
        --------
        proplot.constructor.Formatter
        proplot.ticker.AutoFormatter
        """
        self._sigfig = _not_none(sigfig, 3)
        self._zerotrim = _not_none(zerotrim, rc["formatter.zerotrim"])
        self._base = _not_none(base, 1)

    @docstring._snippet_manager
    def __call__(self, x, pos=None):  # noqa: U100
        """
        %(ticker.call)s
        """
        # Limit to significant figures
        digits = AutoFormatter._decimal_place(x) + self._sigfig - 1
        scale = self._base * 10**-digits
        x = scale * round(x / scale)

        # Create the string
        decimal_point = AutoFormatter._get_default_decimal_point()
        precision = max(0, digits) + max(0, AutoFormatter._decimal_place(self._base))
        string = ("{:.%df}" % precision).format(x)
        string = string.replace(".", decimal_point)

        # Custom string formatting
        string = AutoFormatter._minus_format(string)
        if self._zerotrim:
            string = AutoFormatter._trim_trailing_zeros(string, decimal_point)
        return string


class FracFormatter(mticker.Formatter):
    r"""
    Format numbers as integers or integer fractions. Optionally express the
    values relative to some constant like `numpy.pi`.
    """

    def __init__(self, symbol="", number=1):
        r"""
        Parameters
        ----------
        symbol : str, default: ''
            The constant symbol, e.g. ``r'$\pi$'``.
        number : float, default: 1
            The constant value, e.g. `numpy.pi`.

        Note
        ----
        The fractions shown by this formatter are resolved using the builtin
        `fractions.Fraction` class and `fractions.Fraction.limit_denominator`.

        See also
        --------
        proplot.constructor.Formatter
        proplot.ticker.AutoFormatter
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
            string = "0"
        elif frac.denominator == 1:  # denominator is one
            if frac.numerator == 1 and symbol:
                string = f"{symbol:s}"
            elif frac.numerator == -1 and symbol:
                string = f"-{symbol:s}"
            else:
                string = f"{frac.numerator:d}{symbol:s}"
        else:
            if frac.numerator == 1 and symbol:  # numerator is +/-1
                string = f"{symbol:s}/{frac.denominator:d}"
            elif frac.numerator == -1 and symbol:
                string = f"-{symbol:s}/{frac.denominator:d}"
            else:  # and again make sure we use unicode minus!
                string = f"{frac.numerator:d}{symbol:s}/{frac.denominator:d}"
        string = AutoFormatter._minus_format(string)
        return string


class _CartopyFormatter(object):
    """
    Mixin class for cartopy formatters.
    """

    # NOTE: Cartopy formatters pre 0.18 required axis, and *always* translated
    # input values from map projection coordinates to Plate Carrée coordinates.
    # After 0.18 you can avoid this behavior by not setting axis but really
    # dislike that inconsistency. Solution is temporarily assign PlateCarre().
    def __init__(self, *args, **kwargs):
        import cartopy  # noqa: F401 (ensure available)

        super().__init__(*args, **kwargs)

    def __call__(self, value, pos=None):
        ctx = context._empty_context()
        if self.axis is not None:
            ctx = context._state_context(self.axis.axes, projection=ccrs.PlateCarree())
        with ctx:
            return super().__call__(value, pos)


class DegreeFormatter(_CartopyFormatter, _PlateCarreeFormatter):
    """
    Formatter for longitude and latitude gridline labels.
    Adapted from cartopy.
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
        return ""


class LongitudeFormatter(_CartopyFormatter, LongitudeFormatter):
    """
    Format longitude gridline labels. Adapted from
    `cartopy.mpl.ticker.LongitudeFormatter`.
    """

    @docstring._snippet_manager
    def __init__(self, *args, **kwargs):
        """
        %(ticker.dms)s
        """
        super().__init__(*args, **kwargs)


class LatitudeFormatter(_CartopyFormatter, LatitudeFormatter):
    """
    Format latitude gridline labels. Adapted from
    `cartopy.mpl.ticker.LatitudeFormatter`.
    """

    @docstring._snippet_manager
    def __init__(self, *args, **kwargs):
        """
        %(ticker.dms)s
        """
        super().__init__(*args, **kwargs)
