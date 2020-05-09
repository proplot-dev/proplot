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


def _sanitize_label(string, zerotrim=False):
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


def SigFigFormatter(sigfig=1, zerotrim=False):
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
    def fmt(x, pos):
        if x == 0:
            digits = 0
        else:
            digits = -int(np.log10(abs(x)) // 1)
        digits += sigfig - 1
        x = np.round(x, digits)
        string = ('{:.%df}' % max(0, digits)).format(x)
        if zerotrim and '.' in string:
            string = string.rstrip('0').rstrip('.')
        string = string.replace('-', '\N{MINUS SIGN}')
        if string == '\N{MINUS SIGN}0':
            string = '0'
        return string
    return mticker.FuncFormatter(fmt)


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
