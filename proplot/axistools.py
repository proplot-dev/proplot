#!/usr/bin/env python3
"""
Defines various axis scales, locators, and formatters. Also "registers"
the locator and formatter names, so that they can be called selected with
the `~proplot.axes.CartesianAxes.format` method.
"""
################################################################################
# Developer notes section?
# Too damn long, comments are fine.
# .. raw:: html
#
#    <h1>Scales</h1>
#
#
# * These are complicated. See `~matplotlib.scale.ScaleBase`. Use existing ones
#   as inspiration -- e.g. `InverseScale` modeled after `~matplotlib.scale.LogScale`.
# * Way to think of these is that *every single value you see on an axes first
#   gets secretly converted through some equation*, e.g. logarithm, and plotted
#   linearly in that transformation space.
# * Methods:
#
#    - `get_transform`: Returns a `~matplotlib.transforms.Transform` instance.
#    - `set_default_locators_and_formatters`: Returns
#      the default locators and formatters.
#    - `limit_range_for_scale`: Can be used to raise errors or clip
#      stuff not within that range.
#
#      From the Mercator example: unlike the autoscaling provided by the
#      tick locators, this range limiting will always be adhered to, whether
#      the axis range is set manually, determined automatically or changed
#      through panning and zooming.
#
# * Notes on methods:
#
#     - When you use `set_xlim` or `set_ylim`, the `minpos` used is actually
#       the *data limits* `minpos` (i.e. the minimum coordinate for plotted
#       data). So don't try to e.g. clip data less than 0. That is job for
#       transform. If you use `minpos` in `limit_range_for_scale`, will get
#       wrong and weird results.
#     - Common to use `set_smart_bounds(True)` in the call to
#       `set_default_locators_and_formatters` call -- but this only draws ticks
#       where **data exists**. Often this may not be what we want. Check out
#       source code, see if we can develop own version smarter than this,
#       that still prevents these hanging ticks.
#
# * Note scales have to be *registered* unlike locators and formatters, which
#   can be passed to the setter methods directly.
#
# .. raw:: html
#
#    <h1>Transforms</h1>
#
#
# * These are complicted. See `the transforms module <https://matplotlib.org/_modules/matplotlib/transforms.html#Transform>`_.
# * Attributes:
#
#     - `input_dims`, `output_dims`, `is_separable`, and `has_inverse`. The
#       `dims` are because transforms can be N-D, but for *scales* they are
#       always 1. Note `is_separable` is true if the transform is separable
#       in the x/y dimensions.
#
# * Methods:
#
#     - `transform`: Transforms N-D coordinates, given M x N array of values. Can also
#       just declare `transform_affine` or `transform_non_affine`.
#     - `inverted`: If `has_inverse` is ``True``, performs the inverse transform.
#
# .. raw:: html
#
#    <h1>Locators</h1>
#
#
# * These are complicated. See `the ticker module <https://matplotlib.org/_modules/matplotlib/ticker.html#Locator>`_.
# * Special:
#
#     - `__init__` not defined on base class but *must* be defined for subclass.
#
# * Methods include:
#
#     - `tick_values`: Accepts vmin/vmax and returns values of located ticks
#     - `__call__`: Can return data limits, view limits, or
#       other stuff; not sure how this works or when it's invoked.
#     - `view_limits`: Changes the *view* limits from default `vmin`, `vmax`
#       to prevent singularities. Uses `~matplotlib.transforms.nonsingular`
#       method; for more info see the `matplotlib doc <https://matplotlib.org/_modules/matplotlib/transforms.html#nonsingular>`_.
#
# * Methods that usually can be left alone:
#
#     - `raise_if_exceeds`: Just tests if ticks exceed ``MAXTICKS`` number.
#     - `autoscale`: Calls the internal locator `view_limits` with
#       result of `~matplotlib.axis.Axis.get_view_interval`.
#     - `pan` and `zoom`: Interactive purposes.
#
# .. raw:: html
#
#    <h1>Formatters</h1>
#
# Some of these are easy to construct. Just pass `~matplotlib.ticker.FuncFormatter`
# a function that accepts two arguments the number and 'position', which maybe is used
# for offset or something (don't touch it, leave it default).
#
# The matplotlib default `~matplotlib.ticker.ScalarFormatter` is much more
# complex, but can be overridden in the typical way: adding stuff to the
# `__init__` and `__call__` methods.
#------------------------------------------------------------------------------#
import re
from .utils import _default
from .rctools import rc
from numbers import Number
from fractions import Fraction
import warnings
import numpy as np
import numpy.ma as ma
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import matplotlib.scale as mscale
import matplotlib.transforms as mtransforms
__all__ = [
    'Formatter', 'Locator', 'Scale',
    'AutoFormatter', 'SimpleFormatter', 'FracFormatter',
    'InverseScale', 'MercatorLatitudeScale', 'SineLatitudeScale',
    'CutoffScaleFactory', 'ExpScaleFactory', 'InvertedScaleFactory', 'PowerScaleFactory',
    'formatters', 'locators', 'scales',
    ]

#------------------------------------------------------------------------------#
# Helper functions for instantiating arbitrary Locator and Formatter classes
#------------------------------------------------------------------------------#
# When calling these functions, the format() method should automatically
# detect presence of date axis by testing if unit converter is on axis is
# DateConverter instance
# See: https://matplotlib.org/api/units_api.html
# And: https://matplotlib.org/api/dates_api.html
# Also see: https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/axis.py
# The axis_date() method just sets the converter to the date one
def Locator(locator, *args, **kwargs):
    """
    Returns a `~matplotlib.ticker.Locator` instance, used to interpret the
    `xlocator`, `xlocator_kw`, `ylocator`, `ylocator_kw`, `xminorlocator`,
    `xminorlocator_kw`, `yminorlocator`, and `yminorlocator_kw` arguments when
    passed to `~proplot.axes.CartesianAxes.format_partial`, and the `locator`, `locator_kw`
    `minorlocator`, and `minorlocator_kw` arguments when passed to colorbar
    methods wrapped by `~proplot.wrappers.colorbar_wrapper`.

    Parameters
    ----------
    locator : str, float, or list of float
        If number, specifies the *multiple* used to define tick separation.
        Returns a `~matplotlib.ticker.MultipleLocator` instance.

        If list of numbers, these points are ticked. Returns a
        `~matplotlib.ticker.FixedLocator` instance.

        If string, a dictionary lookup is performed (see below table).
    *args, **kwargs
        Passed to the `~matplotlib.ticker.Locator` class.


    For the `locator` dictionary lookup, options are as follows.

    ======================  ==========================================  =========================================================================================
    Key                     Class                                       Description
    ======================  ==========================================  =========================================================================================
    ``'null'``, ``'none'``  `~matplotlib.ticker.NullLocator`            No ticks
    ``'auto'``              `~matplotlib.ticker.AutoLocator`            Major ticks at sensible locations
    ``'minor'``             `~matplotlib.ticker.AutoMinorLocator`       Minor ticks at sensible locations
    ``'date'``              `~matplotlib.dates.AutoDateLocator`         Default tick locations for datetime axes
    ``'log'``               `~matplotlib.ticker.LogLocator` preset      For log-scale axes, ticks on each power of the base
    ``'logminor'``          `~matplotlib.ticker.LogLocator` preset      For log-scale axes, ticks on the 1st through 9th multiples of each power of the base
    ``'maxn'``              `~matplotlib.ticker.MaxNLocator`            No more than ``N`` ticks at sensible locations
    ``'linear'``            `~matplotlib.ticker.LinearLocator`          Exactly ``N`` ticks encompassing the axis limits, spaced as ``numpy.linspace(lo, hi, N)``
    ``'multiple'``          `~matplotlib.ticker.MultipleLocator`        Ticks every ``N`` step away from zero
    ``'fixed'``             `~matplotlib.ticker.FixedLocator`           Ticks at these exact locations
    ``'index'``             `~matplotlib.ticker.IndexLocator`           Ticks on the non-negative integers
    ``'symmetric'``         `~matplotlib.ticker.SymmetricalLogLocator`  Ticks for symmetrical log-scale axes
    ``'logit'``             `~matplotlib.ticker.LogitLocator`           Ticks for logit-scale axes
    ``'year'``              `~matplotlib.dates.YearLocator`             Ticks every ``N`` years
    ``'month'``             `~matplotlib.dates.MonthLocator`            Ticks every ``N`` months
    ``'weekday'``           `~matplotlib.dates.WeekdayLocator`          Ticks every ``N`` weekdays
    ``'day'``               `~matplotlib.dates.DayLocator`              Ticks every ``N`` days
    ``'hour'``              `~matplotlib.dates.HourLocator`             Ticks every ``N`` hours
    ``'minute'``            `~matplotlib.dates.MinuteLocator`           Ticks every ``N`` minutes
    ``'second'``            `~matplotlib.dates.SecondLocator`           Ticks every ``N`` seconds
    ``'microsecond'``       `~matplotlib.dates.MicrosecondLocator`      Ticks every ``N`` microseconds
    ======================  ==========================================  =========================================================================================

    Returns
    -------
    `~matplotlib.ticker.Locator`
        A `~matplotlib.ticker.Locator` instance.

    Note
    ----
    `~matplotlib.ticker.AutoLocator` has a useful ``nbins`` option. This
    limits the maximum number of ticks to ``nbins-1``.
    """
    if isinstance(locator, mticker.Locator):
        return locator
    # Pull out extra args
    if np.iterable(locator) and not isinstance(locator, str) and not all(isinstance(num, Number) for num in locator):
        locator, args = locator[0], (*locator[1:], *args)
    # Get the locator
    if isinstance(locator, str): # dictionary lookup
        # Shorthands and defaults
        if locator=='logminor':
            locator = 'log'
            kwargs.setdefault('subs', np.arange(10))
        elif locator=='index':
            args = args or (1,)
            if len(args)==1:
                args = (*args, 0)
        # Lookup
        if locator not in locators:
            raise ValueError(f'Unknown locator "{locator}". Options are {", ".join(locators.keys())}.')
        locator = locators[locator](*args, **kwargs)
    elif isinstance(locator, Number): # scalar variable
        locator = mticker.MultipleLocator(locator, *args, **kwargs)
    elif np.iterable(locator):
        locator = mticker.FixedLocator(np.sort(locator), *args, **kwargs) # not necessary
    else:
        raise ValueError(f'Invalid locator "{locator}".')
    return locator

def Formatter(formatter, *args, date=False, **kwargs):
    r"""
    Returns a `~matplotlib.ticker.Formatter` instance, used to interpret the
    `xformatter`, `xformatter_kw`, `yformatter`, and `yformatter_kw` arguments
    when passed to `~proplot.axes.CartesianAxes.format_partial`, and the `formatter`
    and `formatter_kw` arguments when passed to colorbar methods wrapped by
    `~proplot.wrappers.colorbar_wrapper`.

    Parameters
    ----------
    formatter : str, list of str, or function
        If list of strings, ticks are labeled with these strings. Returns a
        `~matplotlib.ticker.FixedFormatter` instance.

        If function, labels will be generated using this function. Returns a
        `~matplotlib.ticker.FuncFormatter` instance.

        If string, there are 4 possibilities:

        1. If string contains ``'%'`` and `date` is ``False``, ticks will be formatted
           using the C-notation ``string % number`` method. See `this page
           <https://docs.python.org/3.4/library/string.html#format-specification-mini-language>`__
           for a review.
        2. If string contains ``'%'`` and `date` is ``True``, datetime
           ``string % number`` formatting is used. See `this page
           <https://docs.python.org/3/library/datetime.html#strftime-strptime-behavior>`__
           for a review.
        3. If string contains ``{x}`` or ``{x:...}``, ticks will be
           formatted by calling ``string.format(x=number)``.
        4. In all other cases, a dictionary lookup is performed (see below table).

    date : bool, optional
        Toggles the behavior when `formatter` contains a ``'%'`` sign (see
        above).
    *args, **kwargs
        Passed to the `~matplotlib.ticker.Formatter` class.


    For the `formatter` dictionary lookup, options are as follows.

    =========================  ============================================  =================================================================
    Key                        Class                                         Description
    =========================  ============================================  =================================================================
    ``'null'``, ``'none'``     `~matplotlib.ticker.NullFormatter`            No tick labels
    ``'auto'``, ``'default'``  `AutoFormatter`                               New default tick labels for axes
    ``'simple'``               `SimpleFormatter`                             New default tick labels for e.g. contour labels
    ``'frac'``                 `FracFormatter`                               Rational fractions
    ``'date'``                 `~matplotlib.dates.AutoDateFormatter`         Default tick labels for datetime axes
    ``'scalar'``               `~matplotlib.ticker.ScalarFormatter`          Old default tick labels for axes
    ``'strmethod'``            `~matplotlib.ticker.StrMethodFormatter`       From the ``string.format`` method
    ``'formatstr'``            `~matplotlib.ticker.FormatStrFormatter`       From C-style ``string % format`` notation
    ``'log'``, ``'sci'``       `~matplotlib.ticker.LogFormatterSciNotation`  For log-scale axes with scientific notation
    ``'math'``                 `~matplotlib.ticker.LogFormatterMathtext`     For log-scale axes with math text
    ``'logit'``                `~matplotlib.ticker.LogitFormatter`           For logistic-scale axes
    ``'eng'``                  `~matplotlib.ticker.EngFormatter`             Engineering notation
    ``'percent'``              `~matplotlib.ticker.PercentFormatter`         Trailing percent sign
    ``'fixed'``                `~matplotlib.ticker.FixedFormatter`           List of strings
    ``'index'``                `~matplotlib.ticker.IndexFormatter`           List of strings corresponding to non-negative integer positions along the axis
    ``'pi'``                   `FracFormatter` preset                        Fractions of :math:`\pi`
    ``'e'``                    `FracFormatter` preset                        Fractions of *e*
    ``'deg'``                  `SimpleFormatter` preset                      Trailing degree symbol
    ``'deglon'``               `SimpleFormatter` preset                      Trailing degree symbol and cardinal "WE" indicator
    ``'deglat'``               `SimpleFormatter` preset                      Trailing degree symbol and cardinal "SN" indicator
    ``'lon'``                  `SimpleFormatter` preset                      Cardinal "WE" indicator
    ``'lat'``                  `SimpleFormatter` preset                      Cardinal "SN" indicator
    =========================  ============================================  =================================================================

    Returns
    -------
    `~matplotlib.ticker.Formatter`
        A `~matplotlib.ticker.Formatter` instance.
    """
    if isinstance(formatter, mticker.Formatter): # formatter object
        return formatter
    # Pull out extra args
    if np.iterable(formatter) and not isinstance(formatter, str) and not all(isinstance(item, str) for item in formatter):
        formatter, args = formatter[0], [*formatter[1:], *args]
    # Get the formatter
    if isinstance(formatter, str): # assumption is list of strings
        # Format strings
        if re.search(r'{x?(:.+)?}', formatter):
            formatter = mticker.StrMethodFormatter(formatter, *args, **kwargs) # new-style .format() form
        elif '%' in formatter:
            if date:
                formatter = mdates.DateFormatter(formatter, *args, **kwargs) # %-style, dates
            else:
                formatter = mticker.FormatStrFormatter(formatter, *args, **kwargs) # %-style, numbers
        else:
            # Fraction shorthands
            if formatter in ('pi', 'e'):
                if formatter=='pi':
                    kwargs.update({'symbol': r'$\pi$', 'number': np.pi})
                else:
                    kwargs.update({'symbol': '$e$', 'number': np.e})
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
                kwargs.update({'suffix':suffix, 'negpos':negpos})
                formatter = 'simple'
            # Lookup
            if formatter not in formatters:
                raise ValueError(f'Unknown formatter "{formatter}". Options are {", ".join(formatters.keys())}.')
            formatter = formatters[formatter](*args, **kwargs)
    elif callable(formatter):
        formatter = mticker.FuncFormatter(formatter, *args, **kwargs)
    elif np.iterable(formatter): # list of strings on the major ticks, wherever they may be
        formatter = mticker.FixedFormatter(formatter)
    else:
        raise ValueError(f'Invalid formatter "{formatter}".')
    return formatter

def Scale(scale, *args, **kwargs):
    """
    Returns the name for a registered or "on-the-fly" generated
    `~matplotlib.scale.ScaleBase` class, used to interpret the `xscale`,
    `xscale_kw`, `yscale`, and `yscale_kw` arguments when passed to
    `~proplot.axes.CartesianAxes.format_partial`.

    Parameters
    ----------
    scale : str or (str, ...)
        The string represents the registered scale name or corresponds to a
        scale factory function (see below table).

        If a tuple or list was passed, the items after the string are passed
        as positional arguments to the corresponding scale factory function.
        They are ignored if the string name does not correspond to a factory
        function.
    **kwargs
        Passed to the scale factory function. Ignored if the scale name does
        not correspond to a factory function.


    For the `scale` dictionary lookup, options are as follows.

    ===============  =======================================  ====================================================
    Key              Class or Factory                         Description
    ===============  =======================================  ====================================================
    ``'linear'``     `~matplotlib.scale.LinearScale`          Linear
    ``'log'``        `~matplotlib.scale.LogScale`             Logarithmic
    ``'symlog'``     `~matplotlib.scale.SymmetricalLogScale`  Logarithmic beyond finite space around zero
    ``'logit'``      `~matplotlib.scale.LogitScale`           Logistic
    ``'inverse'``    `InverseScale`                           Inverse
    ``'sine'``       `SineLatitudeScale`                      Sine function (in degrees)
    ``'mercator'``   `MercatorLatitudeScale`                  Mercator latitude function (in degrees)
    ``'quadratic'``  `PowerScaleFactory` preset               Quadratic function
    ``'cubic'``      `PowerScaleFactory` preset               Cubic function
    ``'height'``     `ExpScaleFactory` preset                 Pressure (in hPa) linear in height
    ``'pressure'``   `ExpScaleFactory` preset                 Height (in km) linear in pressure
    ``'power'``      `PowerScaleFactory`                      Arbitrary power function
    ``'exp'``        `ExpScaleFactory`                        Arbitrary exponential function
    ``'cutoff'``     `CutoffScaleFactory`                     Arbitrary linear transformations
    ===============  =======================================  ====================================================

    Returns
    -------
    `~matplotlib.scale.ScaleBase`
        A `~matplotlib.scale.ScaleBase` instance.
    """
    if isinstance(scale, mscale.ScaleBase):
        mscale.register_scale(scale) # ensure it is registered!
        return scale.name
    # Pull out extra args
    if np.iterable(scale) and not isinstance(scale, str):
        scale, args = scale[0], (*scale[1:], *args)
    # Lookup
    if scale in scales:
        if args:
            warnings.warn(f'Scale constructor ignored positional arguments {args}.')
        if kwargs:
            warnings.warn(f'Scale constructor ignored keyword arguments {kwargs}.')
        return scale # already registered
    # Build an on-the-fly scale
    # NOTE: The factories register the scales in one go
    if scale=='exp':
        return ExpScaleFactory(*args, **kwargs)
    elif scale=='power':
        return PowerScaleFactory(*args, **kwargs)
    elif scale=='cutoff':
        return CutoffScaleFactory(*args, **kwargs)
    else:
        raise ValueError(f'Unknown scale {scale}. Options are {", ".join(scales.keys())}.')

def InvertedScaleFactory(scale, name=None, **kwargs):
    """Returns name of newly registered *inverted* version of the
    `~matplotlib.scale.ScaleBase` corresponding to the scale name. The
    scale name `name` defaults to ``'{scale}_{inverted}'``."""
    scale = Scale(scale, **kwargs) # get the class
    name_ = name or f'{scale}_inverted' # name of inverted version
    class Inverted(scales[scale]):
        name = name_
        def get_transform(self):
            return super().get_transform().inverted() # that's all we need!
    mscale.register_scale(Inverted)
    return name_

#-------------------------------------------------------------------------------
# Formatting classes for mapping numbers (axis ticks) to formatted strings
# Create pseudo-class functions that actually return auto-generated formatting
# classes by passing function references to Funcformatter
#-------------------------------------------------------------------------------
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
    def __init__(self, *args, zerotrim=None, precision=None, tickrange=None,
                        prefix=None, suffix=None, **kwargs):
        """
        Parameters
        ----------
        zerotrim : bool, optional
            Whether to trim trailing zeros.
            Defaults to ``rc['axes.formatter.zerotrim']``.
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
        zerotrim = _default(zerotrim, rc.get('axes.formatter.zerotrim'))
        self._maxprecision = precision
        self._zerotrim = zerotrim
        self._tickrange = tickrange
        self._prefix = prefix or ''
        self._suffix = suffix or ''

    def __call__(self, x, pos=None):
        """
        Parameters
        ----------
        x : float
            The value.
        pos : float, optional
            The position.
        """
        # Tick range limitation
        eps = abs(x)/1000
        tickrange = self._tickrange
        if (x + eps) < tickrange[0] or (x - eps) > tickrange[1]:
            return '' # avoid some ticks
        # Normal formatting
        string = super().__call__(x, pos)
        if self._maxprecision is not None and '.' in string:
            head, tail = string.split('.')
            string = head + '.' + tail[:self._maxprecision]
        if self._zerotrim:
            string = re.sub(r'\.0+$', '', string)
            string = re.sub(r'^(.*\..*?)0+$', r'\1', string)
        string = re.sub(r'^[−-]0$', '0', string) # '-0' to '0'; necessary?
        # Prefix and suffix
        prefix = ''
        if string and string[0] in '−-': # unicode minus or hyphen
            prefix, string = string[0], string[1:]
        return prefix + self._prefix + string + self._suffix

def SimpleFormatter(*args, precision=6,
        prefix=None, suffix=None, negpos=None,
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
    """
    prefix = prefix or ''
    suffix = suffix or ''
    def f(x, pos):
        # Apply suffix if not on equator/prime meridian
        if not negpos:
            negpos_ = ''
        elif x>0:
            negpos_ = negpos[1]
        else:
            x *= -1
            negpos_ = negpos[0]
        # Finally use default formatter
        string = f'{{:.{precision}f}}'.format(x)
        string = re.sub(r'\.0+$', '', string)
        string = re.sub(r'^(.*\..*?)0+$', r'\1', string) # note the non-greedy secondary glob!
        if string=='-0':
            string = '0'
        string = string.replace('-', '\N{MINUS SIGN}')
        return prefix + string + suffix + negpos_
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
        The symbol, e.g. ``r'$\pi$'``. Defaults to ``''``.
    number : float
        The value, e.g. `numpy.pi`. Defaults to ``1``.
    """
    def f(x, pos): # must accept location argument
        frac = Fraction(x/number).limit_denominator()
        if x==0: # zero
            string = '0'
        elif frac.denominator==1: # denominator is one
            if frac.numerator==1 and symbol:
                string = f'{symbol:s}'
            elif frac.numerator==-1 and symbol:
                string = f'-{symbol:s}'
            else:
                string = f'{frac.numerator:d}{symbol:s}'
        elif frac.numerator==1 and symbol: # numerator is +/-1
            string = f'{symbol:s}/{frac.denominator:d}'
        elif frac.numerator==-1 and symbol:
            string = f'-{symbol:s}/{frac.denominator:d}'
        else: # and again make sure we use unicode minus!
            string = f'{frac.numerator:d}{symbol:s}/{frac.denominator:d}'
        return string.replace('-', '\N{MINUS SIGN}')
    # And create FuncFormatter class
    return mticker.FuncFormatter(f)

#------------------------------------------------------------------------------#
# Power axis scale
#------------------------------------------------------------------------------#
def PowerScaleFactory(power, inverse=False, name=None):
    r"""
    Returns a "power scale" that performs the transformation
    :math:`x^{c}`.

    Parameters
    ----------
    power : float
        The power :math:`c` to which :math:`x` is raised.
    inverse : bool, optional
        If ``True``, the "forward" direction performs
        the inverse operation :math:`x^{1/c}`.
    name : str, optional
        The registered scale name. Defaults to ``'polynomial_{power}'``.
    """
    name_ = _default(name, f'polynomial_{power:.1e}')
    class PowerScale(mscale.ScaleBase):
        # Declare name
        name = name_
        def __init__(self, axis, minpos=1e-300, **kwargs):
            super().__init__(axis)
            if not inverse:
                transform = _PowerTransform(power, minpos)
            else:
                transform = _InvertedPowerTransform(power, minpos)
            self._transform = transform
        def limit_range_for_scale(self, vmin, vmax, minpos):
            # return min(vmin, minpos), min(vmax, minpos)
            return vmin, vmax
        def set_default_locators_and_formatters(self, axis):
            axis.set_smart_bounds(True) # unnecessary?
            axis.set_major_formatter(Formatter('default'))
            axis.set_minor_formatter(Formatter('null'))
        def get_transform(self):
            return self._transform

    # Register and return
    mscale.register_scale(PowerScale)
    return name_

class _PowerTransform(mtransforms.Transform):
    """Power transform."""
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True
    def __init__(self, power, minpos):
        super().__init__()
        self.minpos = minpos
        self._power = power
    def transform(self, a):
        aa = np.array(a).copy()
        aa[aa<=self.minpos] = self.minpos # necessary
        return np.power(np.array(a), self._power)
    def transform_non_affine(self, a):
        return self.transform(a)
    def inverted(self):
        return _InvertedPowerTransform(self._power, self.minpos)

class _InvertedPowerTransform(mtransforms.Transform):
    """Inverse power transform."""
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True
    def __init__(self, power, minpos):
        super().__init__()
        self.minpos = minpos
        self._power = power
    def transform(self, a):
        aa = np.array(a).copy()
        aa[aa<=self.minpos] = self.minpos # necessary
        return np.power(np.array(a), 1/self._power)
    def transform_non_affine(self, a):
        return self.transform(a)
    def inverted(self):
        return _PowerTransform(self._power, self.minpos)

#------------------------------------------------------------------------------#
# Exp axis scale
# Why can't this just be a class that accepts args for changing the scale
# params? Because then would need kwargs for every set_xscale call.
#------------------------------------------------------------------------------#
def ExpScaleFactory(base, exp, scale=1, inverse=False, name=None):
    r"""
    Returns an "exponential scale" that performs the transformation
    :math:`Ca^{bx}`.

    This is used when adding a pressure coordinate axis
    for data that is plotted linearly w.r.t. height, or vice versa. Ignore
    this if you're not an atmospheric scientist.

    Parameters
    ----------
    base : float
        The base, i.e. the :math:`a` in :math:`Ca^{bx}`.
    exp : float
        The scale for the exonent, i.e. the :math:`b` in :math:`Ca^{bx}`.
    scale : float, optional
        The coefficient, i.e. the :math:`C` in :math:`Ca^{bx}`.
    inverse : bool, optional
        If ``True``, the "forward" direction performs
        the inverse operation :math:`(\log(x) - \log(C))/(b\log(a))`.
    name : str, optional
        The registered scale name. Defaults to ``'power_{base}_{exp}_{scale}'``.
    """
    name_ = _default(name, f'exp_{base:.1e}_{scale:.1e}_{exp:.1e}')
    class ExpScale(mscale.ScaleBase):
        # Declare name
        name = name_
        def __init__(self, axis, minpos=1e-300, **kwargs):
            super().__init__(axis)
            if not inverse:
                transform = _ExpTransform(base, exp, scale, minpos)
            else:
                transform = _InvertedExpTransform(base, exp, scale, minpos)
            self._transform = transform
        def limit_range_for_scale(self, vmin, vmax, minpos):
            # return min(vmin, minpos), min(vmax, minpos)
            return vmin, vmax
        def set_default_locators_and_formatters(self, axis):
            axis.set_smart_bounds(True) # unnecessary?
            axis.set_major_formatter(Formatter('default'))
            axis.set_minor_formatter(Formatter('null'))
        def get_transform(self):
            return self._transform

    # Register and return
    mscale.register_scale(ExpScale)
    return name_

class _ExpTransform(mtransforms.Transform):
    """Exponential coordinate transform."""
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True
    def __init__(self, base, exp, scale, minpos):
        super().__init__()
        self.minpos = minpos
        self._base = base
        self._exp = exp
        self._scale = scale
    def transform(self, a):
        return self._scale*np.power(self._base, self._exp*np.array(a))
    def transform_non_affine(self, a):
        return self.transform(a)
    def inverted(self):
        return _InvertedExpTransform(self._base, self._exp, self._scale, self.minpos)

class _InvertedExpTransform(mtransforms.Transform):
    """Inverse exponential coordinate transform."""
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True
    def __init__(self, base, exp, scale, minpos):
        super().__init__()
        self.minpos = minpos
        self._base = base
        self._exp = exp
        self._scale = scale
    def transform(self, a):
        aa = np.array(a).copy()
        aa[aa<=self.minpos] = self.minpos # necessary
        return (np.log(aa) - np.log(self._scale))/(np.log(self._base) * self._exp) # this one!
    def transform_non_affine(self, a):
        return self.transform(a)
    def inverted(self):
        return _ExpTransform(self._base, self._exp, self._scale, self.minpos)

#------------------------------------------------------------------------------#
# Cutoff axis
#------------------------------------------------------------------------------#
def CutoffScaleFactory(scale, lower, upper=None):
    """
    Constructs a scale with custom cutoffs.

    If `upper` is ``None``, you have the following two possibilities:

    1. If `scale` is >1, the axis is "faster" to the right
       of `lower`.
    2. If `scale` is <1, the axis is "slower" to the right
       of `lower`.

    If `upper` is not ``None``, you have the following three possibilities:

    1. If `scale` is `numpy.inf`, this puts a cliff between `lower` and
       `upper`. The axis cuts from `lower` to `upper` at a single point.
    2. If `scale` is >1, the axis is accelerated between `lower` and `upper`.
       So the axis is "slow" on the edges but "fast" in middle.
    3. If `scale` is <1, the axis is deccelerated between `lower` and `upper`.
       So the axis is "fast" on the edges but "slow" in middle.

    Todo
    ----
    Create method for drawing those diagonal "cutoff" strokes with whitespace
    between.  See `this post <https://stackoverflow.com/a/5669301/4970632>`_ for
    multi-axis solution and for this class-based solution. Note the space
    between 1-9 in Paul's answer is because actual cutoffs were 0.1 away
    (and tick locs are 0.2 apart).
    """
    if scale<0:
        raise ValueError('Scale must be a positive float.')
    if upper is None and scale==np.inf:
        raise ValueError('For infinite scale (i.e. discrete cutoff), need both lower and upper bounds.')
    name = f'cutoff_{scale:.1e}_{lower:.1e}'
    if upper is not None:
        name = f'{name}_{upper:.1e}'
    name_ = name # have to copy to different name

    class CutoffScale(mscale.ScaleBase):
        # Declare name
        name = name_
        def __init__(self, axis, **kwargs):
            super().__init__(axis)
            self._transform = _CutoffTransform()
        def get_transform(self):
            return self._transform
        def set_default_locators_and_formatters(self, axis):
            axis.set_major_formatter(Formatter('default'))
            axis.set_minor_formatter(Formatter('null'))
            axis.set_smart_bounds(True) # may prevent ticks from extending off sides

    class _CutoffTransform(mtransforms.Transform):
        # Create transform object
        input_dims = 1
        output_dims = 1
        has_inverse = True
        is_separable = True
        def __init__(self):
            super().__init__()
        def transform(self, a):
            a = np.array(a) # very numpy array
            aa = a.copy()
            if upper is None: # just scale between 2 segments
                m = (a > lower)
                aa[m] = a[m] - (a[m] - lower)*(1 - 1/scale)
            elif lower is None:
                m = (a < upper)
                aa[m] = a[m] - (upper - a[m])*(1 - 1/scale)
            else:
                m1 = (a > lower)
                m2 = (a > upper)
                m3 = (a > lower) & (a < upper)
                if scale==np.inf:
                    aa[m1] = a[m1] - (upper - lower)
                    aa[m3] = lower
                else:
                    aa[m2] = a[m2] - (upper - lower)*(1 - 1/scale)
                    aa[m3] = a[m3] - (a[m3] - lower)*(1 - 1/scale)
            return aa
        def transform_non_affine(self, a):
            return self.transform(a)
        def inverted(self):
            return _InvertedCutoffTransform()

    class _InvertedCutoffTransform(mtransforms.Transform):
        # Inverse of _CutoffTransform
        input_dims = 1
        output_dims = 1
        has_inverse = True
        is_separable = True
        def __init__(self):
            super().__init__()
        def transform(self, a):
            a = np.array(a)
            aa = a.copy()
            if upper is None:
                m = (a > lower)
                aa[m] = a[m] + (a[m] - lower)*(1 - 1/scale)
            elif lower is None:
                m = (a < upper)
                aa[m] = a[m] + (upper - a[m])*(1 - 1/scale)
            else:
                n = (upper-lower)*(1 - 1/scale)
                m1 = (a > lower)
                m2 = (a > upper - n)
                m3 = (a > lower) & (a < (upper - n))
                if scale==np.inf:
                    aa[m1] = a[m1] + (upper - lower)
                else:
                    aa[m2] = a[m2] + n
                    aa[m3] = a[m3] + (a[m3] - lower)*(1 - 1/scale)
            return aa
        def transform_non_affine(self, a):
            return self.transform(a)
        def inverted(self):
            return _CutoffTransform()

    # Register and return
    mscale.register_scale(CutoffScale)
    return name_

#------------------------------------------------------------------------------#
# Geographic axes
#------------------------------------------------------------------------------#
class MercatorLatitudeScale(mscale.ScaleBase):
    r"""
    Scales axis as with latitudes in the `Mercator projection
    <http://en.wikipedia.org/wiki/Mercator_projection>`__. Inspired by
    :cite:`barnes_rossby_2011`, and adapted from `this matplotlib example
    <https://matplotlib.org/examples/api/custom_scale_example.html>`__.

    The scale function is as follows.

    .. math::

        y = \ln(\tan(\pi x/180) + \sec(\pi x/180))

    The inverse scale function is as follows.

    .. math::

        x = 180\arctan(\sinh(y))/\pi

    Also uses a user-defined threshold :math:`\in (-90, 90)`, above and
    below which nothing will be plotted.

    .. bibliography:: ../refs.bib
    """
    name = 'mercator'
    """Registered scale name."""
    def __init__(self, axis, *, thresh=85.0, **kwargs):
        super().__init__(axis)
        if thresh >= 90.0:
            raise ValueError('Threshold "thresh" must be <=90.')
        self.thresh = thresh
    def get_transform(self):
        """See `~matplotlib.scale.ScaleBase`."""
        return _MercatorLatitudeTransform(self.thresh)
    def limit_range_for_scale(self, vmin, vmax, minpos):
        """See `~matplotlib.scale.ScaleBase`."""
        return max(vmin, -self.thresh), min(vmax, self.thresh)
    def set_default_locators_and_formatters(self, axis):
        """See `~matplotlib.scale.ScaleBase`."""
        axis.set_smart_bounds(True)
        axis.set_major_formatter(Formatter('deg'))
        axis.set_minor_formatter(Formatter('null'))

class _MercatorLatitudeTransform(mtransforms.Transform):
    # Default attributes
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True
    def __init__(self, thresh):
        super().__init__()
        self.thresh = thresh
    def transform_non_affine(self, a):
        # With safeguards
        # TODO: Can improve this?
        a = np.deg2rad(a) # convert to radians
        m = ma.masked_where((a < -self.thresh) | (a > self.thresh), a)
        if m.mask.any():
            return ma.log(np.abs(ma.tan(m) + 1/ma.cos(m)))
        else:
            return np.log(np.abs(np.tan(a) + 1/np.cos(a)))
    def inverted(self):
        return _InvertedMercatorLatitudeTransform(self.thresh)

class _InvertedMercatorLatitudeTransform(mtransforms.Transform):
    # As above, but for the inverse transform
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True
    def __init__(self, thresh):
        super().__init__()
        self.thresh = thresh
    def transform_non_affine(self, a):
        # m = ma.masked_where((a < -self.thresh) | (a > self.thresh), a)
        return np.rad2deg(np.arctan2(1, np.sinh(a))) # always assume in first/fourth quadrant, i.e. go from -pi/2 to pi/2
    def inverted(self):
        return _MercatorLatitudeTransform(self.thresh)

class SineLatitudeScale(mscale.ScaleBase):
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
    """Registered scale name."""
    def __init__(self, axis, **kwargs):
        super().__init__(axis)
    def get_transform(self):
        """See `~matplotlib.scale.ScaleBase`."""
        return _SineLatitudeTransform()
    def limit_range_for_scale(self, vmin, vmax, minpos):
        """See `~matplotlib.scale.ScaleBase`."""
        return vmin, vmax
    def set_default_locators_and_formatters(self, axis):
        """See `~matplotlib.scale.ScaleBase`."""
        axis.set_smart_bounds(True)
        axis.set_major_formatter(Formatter('deg'))
        axis.set_minor_formatter(Formatter('null'))

class _SineLatitudeTransform(mtransforms.Transform):
    # Default attributes
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True
    def __init__(self):
        # Initialize, declare attribute
        super().__init__()
    def transform_non_affine(self, a):
        # With safeguards
        # TODO: Can improve this?
        with np.errstate(invalid='ignore'): # NaNs will always be False
            m = (a >= -90) & (a <= 90)
        if not m.all():
            aa = ma.masked_where(~m, a)
            return ma.sin(np.deg2rad(aa))
        else:
            return np.sin(np.deg2rad(a))
    def inverted(self):
        return _InvertedSineLatitudeTransform()

class _InvertedSineLatitudeTransform(mtransforms.Transform):
    # Inverse of _SineLatitudeTransform
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True
    def __init__(self):
        super().__init__()
    def transform_non_affine(self, a):
        # Clipping, instead of setting invalid
        # NOTE: Using ma.arcsin below caused super weird errors, dun do that
        aa = a.copy()
        return np.rad2deg(np.arcsin(aa))
    def inverted(self):
        return _SineLatitudeTransform()

#------------------------------------------------------------------------------#
# Other transformations
#------------------------------------------------------------------------------#
class InverseScale(mscale.ScaleBase):
    r"""
    Scales axis to be linear in the *inverse* of *x*. The scale
    function and inverse scale function are as follows.

    .. math::

        y = x^{-1}

    """
    # Developer notes:
    # Unlike log-scale, we can't just warp the space between
    # the axis limits -- have to actually change axis limits. Also this
    # scale will invert and swap the limits you provide. Weird! But works great!
    # Declare name
    name = 'inverse'
    """Registered scale name."""
    def __init__(self, axis, minpos=1e-300, **kwargs):
        super().__init__(axis)
        self.minpos = minpos
    def get_transform(self):
        """See `~matplotlib.scale.ScaleBase`."""
        return _InverseTransform(self.minpos)
    def limit_range_for_scale(self, vmin, vmax, minpos):
        """See `~matplotlib.scale.ScaleBase`."""
        return min(vmin, minpos), min(vmax, minpos)
    def set_default_locators_and_formatters(self, axis):
        """See `~matplotlib.scale.ScaleBase`."""
        # TODO: Fix minor locator issue
        # NOTE: Log formatter can ignore certain major ticks! Why is that?
        axis.set_smart_bounds(True) # may prevent ticks from extending off sides
        axis.set_major_locator(mticker.LogLocator(base=10, subs=[1, 2, 5]))
        axis.set_minor_locator(mticker.LogLocator(base=10, subs='auto'))
        axis.set_major_formatter(Formatter('default')) # use 'log' instead?
        axis.set_minor_formatter(Formatter('null')) # use 'minorlog' instead?

class _InverseTransform(mtransforms.Transform):
    # Create transform object
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True
    def __init__(self, minpos):
        super().__init__()
        self.minpos = minpos
    def transform(self, a):
        a = np.array(a)
        aa = a.copy()
        # f = np.abs(a)<=self.minpos # attempt for negative-friendly
        # aa[f] = np.sign(a[f])*self.minpos
        aa[aa<=self.minpos] = self.minpos
        return 1.0/aa
    def transform_non_affine(self, a):
        return self.transform(a)
    def inverted(self):
        return _InverseTransform(self.minpos)

#------------------------------------------------------------------------------#
# Declare dictionaries
# Includes some custom classes, so has to go at end
#------------------------------------------------------------------------------#
scales = mscale._scale_mapping
"""The registered scale names and their associated
`~matplotlib.scale.ScaleBase` classes. See `Scale` for a table."""

locators = {
    'none':        mticker.NullLocator,
    'null':        mticker.NullLocator,
    'auto':        mticker.AutoLocator,
    'log':         mticker.LogLocator,
    'maxn':        mticker.MaxNLocator,
    'linear':      mticker.LinearLocator,
    'multiple':    mticker.MultipleLocator,
    'fixed':       mticker.FixedLocator,
    'index':       mticker.IndexLocator,
    'symmetric':   mticker.SymmetricalLogLocator,
    'logit':       mticker.LogitLocator,
    'minor':       mticker.AutoMinorLocator,
    'microsecond': mdates.MicrosecondLocator,
    'second':      mdates.SecondLocator,
    'minute':      mdates.MinuteLocator,
    'hour':        mdates.HourLocator,
    'day':         mdates.DayLocator,
    'weekday':     mdates.WeekdayLocator,
    'month':       mdates.MonthLocator,
    'year':        mdates.YearLocator,
    }
"""Mapping of strings to `~matplotlib.ticker.Locator` classes. See
`Locator` for a table."""

formatters = { # note default LogFormatter uses ugly e+00 notation
    'default':    AutoFormatter,
    'auto':       AutoFormatter,
    'frac':       FracFormatter,
    'simple':     SimpleFormatter,
    'scalar':     mticker.ScalarFormatter,
    'none':       mticker.NullFormatter,
    'null':       mticker.NullFormatter,
    'strmethod':  mticker.StrMethodFormatter,
    'formatstr':  mticker.FormatStrFormatter,
    'log':        mticker.LogFormatterSciNotation,
    'sci':        mticker.LogFormatterSciNotation,
    'math':       mticker.LogFormatterMathtext,
    'logit':      mticker.LogitFormatter,
    'eng':        mticker.EngFormatter,
    'percent':    mticker.PercentFormatter,
    'index':      mticker.IndexFormatter,
    }
"""Mapping of strings to `~matplotlib.ticker.Formatter` classes. See
`Formatter` for a table."""

# Register scale names, so user can set_xscale and set_yscale with strings.
# Misc
mscale.register_scale(InverseScale)
# Common powers
PowerScaleFactory(2, 'quadratic')
PowerScaleFactory(3, 'cubic')
# Geographic coordinates
mscale.register_scale(SineLatitudeScale)
mscale.register_scale(MercatorLatitudeScale)
# Height coordinates
# TODO: Some overlap maybe, since this sort-of duplicates a log scale?
ExpScaleFactory(np.e, -1/7, 1013.25, True, 'height')   # scale pressure so it matches a height axis
ExpScaleFactory(np.e, -1/7, 1013.25, False, 'pressure') # scale height so it matches a pressure axis

