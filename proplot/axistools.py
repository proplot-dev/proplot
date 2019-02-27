#!/usr/bin/env python3
"""
Defines various axis scales, locators, and formatters. Also "registers"
the locator and formatter names, so that they can be called selected with
the `~proplot.axes.XYAxes.format` method.

Below is the rough overview of matplotlib API.

.. raw:: html

   <h1>Review</h1>

We want to try to avoid using the Formatter to scale/transform values, and
passing the locator an array of scaled/transformed values.

Makes more sense
to instead define separate *axis transforms*, then can use locators and
formatters like normal, as they were intended to be used. This way, if e.g.
matching frequency-axis with wavelength-axis, just conver the *axis limits*
so they match, then you're good.

.. raw:: html

   <h2>Scales</h2>


* These are complicated. See `~matplotlib.scale.ScaleBase`. Use existing ones
  as inspiration -- e.g. `InverseScale` modeled after `~matplotlib.scale.LogScale`.
* Way to think of these is that *every single value you see on an axes first
  gets secretly converted through some equation*, e.g. logarithm, and plotted
  linearly in that transformation space.
* Methods are as follows:

   - `get_transform`: Returns a `~matplotlib.transforms.Transform` instance.
   - `set_default_locators_and_formatters`: Returns
     the default locators and formatters.
   - `limit_range_for_scale`: Can be used to raise errors or clip
     stuff not within that range.

     From the Mercator example: unlike the autoscaling provided by the
     tick locators, this range limiting will always be adhered to, whether
     the axis range is set manually, determined automatically or changed
     through panning and zooming.

* Important notes on methods:

    - When you use `set_xlim` or `set_ylim`, the `minpos` used is actually
      the *data limits* `minpos` (i.e. the minimum coordinate for plotted
      data). So don't try to e.g. clip data less than 0. That is job for
      transform. If you use `minpos` in `limit_range_for_scale`, will get
      wrong and weird results.
    - Common to use `set_smart_bounds(True)` in the call to
      `set_default_locators_and_formatters` call -- but this only draws ticks
      where **data exists**. Often this may not be what we want. Check out
      source code, see if we can develop own version smarter than this,
      that still prevents these hanging ticks.

* Note scales have to be *registered* unlike locators and formatters, which
  can be passed to the setter methods directly.

.. raw:: html

   <h2>Transforms</h2>


* These are complicted. See `the transforms module <https://matplotlib.org/_modules/matplotlib/transforms.html#Transform>`_.
* Attributes:
    - `input_dims`, `output_dims`, `is_separable`, and `has_inverse`. The
      `dims` are because transforms can be N-D, but for *scales* they are
      always 1. Note `is_separable` is true if the transform is separable
      in the x/y dimensions.

* Methods:
    - `transform`: Transforms N-D coordinates, given M x N array of values. Can also
      just declare `transform_affine` or `transform_non_affine`.
    - `inverted`: If `has_inverse` is ``True``, performs the inverse transform.

.. raw:: html

   <h2>Locators</h2>


* These are complicated. See `the ticker module <https://matplotlib.org/_modules/matplotlib/ticker.html#Locator>`_.
* Special:
    - `__init__` not defined on base class but *must* be defined for subclass.

* Methods include:
    - `tick_values`: Accepts vmin/vmax and returns values of located ticks
    - `__call__`: Can return data limits, view limits, or
      other stuff; not sure how this works or when it's invoked.
    - `view_limits`: Changes the *view* limits from default `vmin`, `vmax`
      to prevent singularities. Uses `~matplotlib.transforms.nonsingular`
      method; for more info see the `matplotlib doc <https://matplotlib.org/_modules/matplotlib/transforms.html#nonsingular>`_.

* Methods that usually can be left alone:
    - `raise_if_exceeds`: Just tests if ticks exceed ``MAXTICKS`` number.
    - `autoscale`: Calls the internal locator `view_limits` with
      result of ``axis.get_view_interval(...)``.
    - `pan` and `zoom`: Interactive purposes.

.. raw:: html

   <h2>Formatters</h2>


Easy to construct: just build with `~matplotlib.formatter.FuncFormatter`
a function that accepts the number and a 'position', which maybe is used
for offset or something (don't touch it, leave it default).
"""
#------------------------------------------------------------------------------#
# Imports
#------------------------------------------------------------------------------#
import re
from .utils import ic
from numbers import Number
from fractions import Fraction
from types import FunctionType
import numpy as np
import numpy.ma as ma
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import matplotlib.scale as mscale
import matplotlib.transforms as mtransforms

#------------------------------------------------------------------------------#
# Tick scales
#------------------------------------------------------------------------------#
def Scale(scale, **kwargs):
    """
    Returns a `~matplotlib.scale.ScaleBase` instance.

    Parameters
    ----------
    scale : str or (str, *args)
        If str, use the specified names.

        If tuple, ``args`` are passed to the class
        instantiator. Ignored uncless the str is ``'cutoff'``,
        ``'exp'``, ``'height'``, or ``'pressure'``.

        String name options are as follows:

        ============  =======================================  =========================================================
        Key           Class                                    Description
        ============  =======================================  =========================================================
        `'linear'`    `~matplotlib.scale.LinearScale`          Linear
        `'log'`       `~matplotlib.scale.LogScale`             Logarithmic
        `'symlog'`    `~matplotlib.scale.SymmetricalLogScale`  Logarithmic beyond space around zero
        `'logit'`     `~matplotlib.scale.LogitScale`           Logistic
        `'pressure'`  `ExpScale`                               Scale pressure coords to be linear in height
        `'height'`    `ExpScale`                               Scale height coords to be linear in pressure
        `'exp'`       `ExpScale`                               Scale with some exponential function
        `'sine'`      `SineLatitudeScale`                      Scale with sine function (in degrees)
        `'mercator'`  `MercatorLatitudeScale`                  Scale with Mercator latitude projection coords
        `'inverse'`   `InverseScale`                           Scale with the inverse
        ============  =======================================  =========================================================

    Other parameters
    ----------------
    **kwargs
        Passed to the ``Scale`` class on instantiation.
    """
    args = []
    if np.iterable(scale) and not isinstance(scale, str):
        scale, args = scale[0], scale[1:]
    if scale in scales and not args:
        pass # already registered
    elif scale=='cutoff':
        scale = CutoffScaleFactory(*args, **kwargs)
    elif scale in ('exp', 'height', 'pressure'): # note here args is non-zero
        if scale=='height':
            if len(args)!=1:
                raise ValueError('Only one non-keyword arg allowed.')
            args = [*args, True]
        if scale=='pressure':
            if len(args)!=1:
                raise ValueError('Only one non-keyword arg allowed.')
            args = [*args, False]
        scale = ExpScaleFactory(*args, **kwargs)
    else:
        raise ValueError(f'Unknown scale {scale}.')
    return scale

#------------------------------------------------------------------------------#
# Exp axis
#------------------------------------------------------------------------------#
def ExpScaleFactory(scale, to_exp=True, name='exp'):
    """
    Exponential scale, used e.g. when adding a pressure coordinate axis
    for data plotted linear w.r.t. height (this is atmospheric science stuff).
    """
    scale_num = scale
    scale_name = name # must make a copy
    class ExpScale(mscale.ScaleBase):
        name = scale_name # assigns as attribute
        """Registered scale name."""
        scale = scale_num
        forward = to_exp
        # Declare name
        def __init__(self, axis, thresh=1e-300, **kwargs):
            # Initialize
            mscale.ScaleBase.__init__(self)
            self.thresh = thresh

        def limit_range_for_scale(self, vmin, vmax, minpos):
            # Prevent conversion from inverting axis scale, which
            # happens when scale is less than zero
            # 99% of time this is want user will want I think
            # NOTE: Do not try to limit data to range above zero here, that
            # is transform's job; see header for info.
            if self.scale < 0:
                vmax, vmin = vmin, vmax
            return vmin, vmax

        def set_default_locators_and_formatters(self, axis):
            # Consider changing this
            # axis.set_smart_bounds(True) # may prevent ticks from extending off sides
            axis.set_major_formatter(Formatter('custom'))
            axis.set_minor_formatter(Formatter('null'))

        def get_transform(self):
            # Either sub into e(scale*z), the default, or invert
            # the exponential
            if self.forward:
                return _ExpTransform(self.scale, self.thresh)
            else:
                return _InvertedExpTransform(self.scale, self.thresh)

    # Register and return
    mscale.register_scale(ExpScale)
    return scale_name

class _ExpTransform(mtransforms.Transform):
    # Exponential coordinate transform
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True
    def __init__(self, scale, thresh):
        mtransforms.Transform.__init__(self)
        self.thresh = thresh
        self.scale = scale
    def transform(self, a):
        return np.exp(self.scale*np.array(a))
    def transform_non_affine(self, a):
        return self.transform(a)
    def inverted(self):
        return _InvertedExpTransform(self.scale, self.thresh)

class _InvertedExpTransform(mtransforms.Transform):
    # Inverse of _ExpTransform
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True
    def __init__(self, scale, thresh):
        mtransforms.Transform.__init__(self)
        self.thresh = thresh
        self.scale = scale
    def transform(self, a):
        a = np.array(a)
        aa = a.copy()
        aa[a<=self.thresh] = self.thresh
        return np.log(aa)/self.scale
    def transform_non_affine(self, a):
        return self.transform(a)
    def inverted(self):
        return _ExpTransform(self.scale, self.thresh)

#------------------------------------------------------------------------------#
# Cutoff axis
#------------------------------------------------------------------------------#
def CutoffScaleFactory(scale, lower, upper=None, name='cutoff'):
    """
    Constructer for scale with custom cutoffs. Three options here:

      1. Put a 'cliff' between two numbers (default).
      2. Accelerate the scale gradient between two numbers (scale>1).
      3. Deccelerate the scale gradient between two numbers (scale<1). So
         scale is fast on edges but slow in middle.

    Todo
    ----
    Alongside this, create method for drawing those cutoff diagonal marks
    with white space between.

    See `this post <https://stackoverflow.com/a/5669301/4970632>`_ for
    multi-axis solution and for this class-based solution. Note the space
    between 1-9 in Paul's answer is because actual cutoffs were 0.1 away
    (and tick locs are 0.2 apart).
    """
    scale_name = name # have to copy to different name
    if scale<0:
        raise ValueError('Scale must be a positive float.')
    if upper is None:
        if scale==np.inf:
            raise ValueError('For infinite scale (i.e. discrete cutoff), need both lower and upper bounds.')

    class CutoffScale(mscale.ScaleBase):
        # Declare name
        name = scale_name
        """Registered scale name."""
        def __init__(self, axis, **kwargs):
            mscale.ScaleBase.__init__(self)
            self.name = scale_name

        def get_transform(self):
            return _CutoffTransform()

        def set_default_locators_and_formatters(self, axis):
            axis.set_major_formatter(Formatter('custom'))
            axis.set_minor_formatter(Formatter('null'))
            axis.set_smart_bounds(True) # may prevent ticks from extending off sides

    # Register and return
    mscale.register_scale(CutoffScale)
    # print(f'Registered scale "{scale_name}".')
    return scale_name

class _CutoffTransform(mtransforms.Transform):
    # Create transform object
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True
    def __init__(self):
        mtransforms.Transform.__init__(self)
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
        mtransforms.Transform.__init__(self)
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

#------------------------------------------------------------------------------#
# Mercator axis
#------------------------------------------------------------------------------#
class MercatorLatitudeScale(mscale.ScaleBase):
    r"""
    Scales axis as with latitudes in the `Mercator projection
    <http://en.wikipedia.org/wiki/Mercator_projection>`_. Inspired by
    :cite:`barnes_rossby_2011`, and adapted from `this matplotlib example
    <https://matplotlib.org/examples/api/custom_scale_example.html>`_.

    The scale function is as follows:

    .. math::

        \ln(\tan(y) + \sec(y))

    The inverse scale function is as follows:

    .. math::

        \arctan(\sinh(y))

    Also uses a user-defined threshold :math:`\in (-90, 90)`, above and
    below which nothing will be plotted.

    .. bibliography:: ../refs.bib
    """
    name = 'mercator'
    """Registered scale name."""
    def __init__(self, axis, *, thresh=85.0, **kwargs):
        # Initialize
        mscale.ScaleBase.__init__(self)
        if thresh >= 90.0:
            raise ValueError('Threshold "thresh" must be <=90.')
        self.thresh = thresh

    def get_transform(self):
        # Return special transform object
        return _MercatorLatitudeTransform(self.thresh)

    def limit_range_for_scale(self, vmin, vmax, minpos):
        # *Hard* limit on axis boundaries
        return max(vmin, -self.thresh), min(vmax, self.thresh)

    def set_default_locators_and_formatters(self, axis):
        # Apply these
        axis.set_smart_bounds(True)
        axis.set_major_locator(Locator(20)) # every 20 degrees
        axis.set_major_formatter(Formatter('deg'))
        axis.set_minor_formatter(Formatter('null'))

class _MercatorLatitudeTransform(mtransforms.Transform):
    # Default attributes
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True
    def __init__(self, thresh):
        # Initialize, declare attribute
        mtransforms.Transform.__init__(self)
        self.thresh = thresh
    def transform_non_affine(self, a):
        # For M N-dimensional transform, transform MxN into result
        # So numbers stay the same, but data will then be linear in the
        # result of the math below.
        a = np.radians(a) # convert to radians
        m = ma.masked_where((a < -self.thresh) | (a > self.thresh), a)
        # m[m.mask] = np.nan
        # a[m.mask] = np.nan
        if m.mask.any():
            return ma.log(np.abs(ma.tan(m) + 1.0 / ma.cos(m)))
        else:
            return np.log(np.abs(np.tan(a) + 1.0 / np.cos(a)))
    def inverted(self):
        # Just call inverse transform class
        return _InvertedMercatorLatitudeTransform(self.thresh)

class _InvertedMercatorLatitudeTransform(mtransforms.Transform):
    # As above, but for the inverse transform
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True
    def __init__(self, thresh):
        mtransforms.Transform.__init__(self)
        self.thresh = thresh
    def transform_non_affine(self, a):
        # m = ma.masked_where((a < -self.thresh) | (a > self.thresh), a)
        return np.degrees(np.arctan2(1, np.sinh(a))) # always assume in first/fourth quadrant, i.e. go from -pi/2 to pi/2
    def inverted(self):
        return _MercatorLatitudeTransform(self.thresh)

#------------------------------------------------------------------------------#
# Sine axis
#------------------------------------------------------------------------------#
class SineLatitudeScale(mscale.ScaleBase):
    r"""
    The scale function is...

    .. math:

        \sin((y))

    The inverse scale function:

    .. math:

        (\arcsin(y))
    """
    name = 'sine'
    """Registered scale name."""
    def __init__(self, axis, **kwargs):
        # Initialize
        mscale.ScaleBase.__init__(self)

    def get_transform(self):
        # Return special transform object
        return _SineLatitudeTransform()

    def limit_range_for_scale(self, vmin, vmax, minpos):
        # *Hard* limit on axis boundaries
        return vmin, vmax
        # return max(vmin, -90), min(vmax, 90)

    def set_default_locators_and_formatters(self, axis):
        # Apply these
        axis.set_smart_bounds(True)
        axis.set_major_locator(Locator(20)) # every 20 degrees
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
        mtransforms.Transform.__init__(self)
    def transform_non_affine(self, a):
        # Transformation
        with np.errstate(invalid='ignore'): # NaNs will always be False
            m = (a >= -90) & (a <= 90)
        if not m.all():
            aa = ma.masked_where(~m, a)
            return ma.sin(np.deg2rad(aa))
        else:
            return np.sin(np.deg2rad(a))
    def inverted(self):
        # Just call inverse transform class
        return _InvertedSineLatitudeTransform()

class _InvertedSineLatitudeTransform(mtransforms.Transform):
    # Inverse of _SineLatitudeTransform
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True
    def __init__(self):
        mtransforms.Transform.__init__(self)
    def transform_non_affine(self, a):
        # Clipping, instead of setting invalid
        # NOTE: Using ma.arcsin below caused super weird errors, dun do that
        aa = a.copy()
        return np.rad2deg(np.arcsin(aa))
    def inverted(self):
        return _SineLatitudeTransform()

#------------------------------------------------------------------------------#
# Scale as the *inverse* of the axis coordinate
#------------------------------------------------------------------------------#
class InverseScale(mscale.ScaleBase):
    """
    Similar to `~matplotlib.scale.LogScale`, but this scales to be linear
    in the *inverse* of *x*. Very useful e.g. for plotting wavelengths
    on twin axis with wavenumbers.
    """
    # Developer notes:
    # Unlike log-scale, we can't just warp the space between
    # the axis limits -- have to actually change axis limits. Also this
    # scale will invert and swap the limits you provide. Weird! But works great!
    # Declare name
    name = 'inverse'
    """Registered scale name."""
    def __init__(self, axis, minpos=1e-2, **kwargs):
        # Initialize (note thresh is always needed)
        mscale.ScaleBase.__init__(self)
        self.minpos = minpos

    def get_transform(self):
        # Return transform class
        return _InverseTransform(self.minpos)

    def limit_range_for_scale(self, vmin, vmax, minpos):
        # *Hard* limit on axis boundaries
        if not np.isfinite(minpos):
            minpos = 1e-300
        return (minpos if vmin <= 0 else vmin,
                minpos if vmax <= 0 else vmax)

    def set_default_locators_and_formatters(self, axis):
        # TODO: fix minor locator issue
        # NOTE: log formatter can ignore certain major ticks! why is that?
        axis.set_smart_bounds(True) # may prevent ticks from extending off sides
        axis.set_major_locator(mticker.LogLocator(base=10, subs=[1, 2, 5]))
        axis.set_minor_locator(mticker.LogLocator(base=10, subs='auto'))
        axis.set_major_formatter(Formatter('custom'))
        axis.set_minor_formatter(Formatter('null'))
        # axis.set_major_formatter(mticker.LogFormatter())

class _InverseTransform(mtransforms.Transform):
    # Create transform object
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True
    def __init__(self, minpos):
        mtransforms.Transform.__init__(self)
        self.minpos = minpos
    def transform(self, a):
        a = np.array(a)
        aa = a.copy()
        aa[a<=0] = self.minpos
        # aa[a<=0] = np.nan # minpos
        return 1.0/aa
    def transform_non_affine(self, a):
        return self.transform(a)
    def inverted(self):
        return _InvertedInverseTransform(self.minpos)

class _InvertedInverseTransform(mtransforms.Transform):
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True
    def __init__(self, minpos):
        mtransforms.Transform.__init__(self)
        self.minpos = minpos
    def transform(self, a):
        a = np.array(a)
        aa = a.copy()
        aa[a<=0] = self.minpos
        # aa[a<=0] = np.nan # messes up automatic ylim setting
        return 1.0/aa
    def transform_non_affine(self, a):
        return self.transform(a)
    def inverted(self):
        return _InverseTransform(self.minpos)

# Register hard-coded scale names, so user can set_xscale and set_yscale with strings
mscale.register_scale(InverseScale)
mscale.register_scale(SineLatitudeScale)
mscale.register_scale(MercatorLatitudeScale)
ExpScaleFactory(-1.0/7, False, 'pressure') # scale pressure so it matches a height axis
ExpScaleFactory(-1.0/7, True,  'height') # scale height so it matches a pressure axis

#------------------------------------------------------------------------------#
# Helper functions for instantiating arbitrary Locator and Formatter classes
# When calling these functions, the format() method should automatically
# detect presence of date axis by testing if unit converter is on axis is
# DateConverter instance
# See: https://matplotlib.org/api/units_api.html
# And: https://matplotlib.org/api/dates_api.html
# Also see: https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/axis.py
# The axis_date() method just sets the converter to the date one
#------------------------------------------------------------------------------#
def Locator(locator, *args, minor=False, time=False, **kwargs):
    """
    Returns a `~matplotlib.ticker.Locator` instance.

    Parameters
    ----------
    locator : None, float, list of float, or str
        If None, returns the default `~matpltolib.ticker.AutoLocator`,
        unless `minor` or `time` are ``True`` (see below).

        If number, specifies the *multiple* used to define tick separation.
        Returns a `~matplotlib.ticker.MultipleLocator` instance.

        If list of numbers, these points are ticked. Returns a
        `~matplotlib.ticker.FixedLocator` instance.

        If str, a dictionary lookup is performed. Options are as follows:

        ===============  ==========================================
        Key              Class
        ===============  ==========================================
        `'none'`         `~matplotlib.ticker.NullLocator`
        `'null'`         `~matplotlib.ticker.NullLocator`
        `'log'`          `~matplotlib.ticker.LogLocator`
        `'maxn'`         `~matplotlib.ticker.MaxNLocator`
        `'linear'`       `~matplotlib.ticker.LinearLocator`
        `'log'`          `~matplotlib.ticker.LogLocator`
        `'multiple'`     `~matplotlib.ticker.MultipleLocator`
        `'fixed'`        `~matplotlib.ticker.FixedLocator`
        `'index'`        `~matplotlib.ticker.IndexLocator`
        `'symmetric'`    `~matplotlib.ticker.SymmetricalLogLocator`
        `'logit'`        `~matplotlib.ticker.LogitLocator`
        `'minor'`        `~matplotlib.ticker.AutoMinorLocator`
        `'microsecond'`  `~matplotlib.dates.MicrosecondLocator`
        `'second'`       `~matplotlib.dates.SecondLocator`
        `'minute'`       `~matplotlib.dates.MinuteLocator`
        `'hour'`         `~matplotlib.dates.HourLocator`
        `'day'`          `~matplotlib.dates.DayLocator`
        `'weekday'`      `~matplotlib.dates.WeekdayLocator`
        `'month'`        `~matplotlib.dates.MonthLocator`
        `'year'`         `~matplotlib.dates.YearLocator`
        ===============  ==========================================

    minor : bool, optional
        Ignored if `locator` is not ``None``. Otherwise, if ``True``, returns
        a `~matplotlib.ticker.AutoMinorLocator` instance.

    time : bool, optional
        Ignored if `locator` is not ``None``. Otherwise, if ``True``, returns
        a `~matplotlib.ticker.AutoDateLocator` instance.

    Other parameters
    ----------------
    **kwargs
        Passed to the ``Locator`` class on instantiation.

    Notes
    -----
    `~matpltolib.ticker.AutoLocator` has a useful ``nbins`` option; no more
    than ``nbins-1`` ticks are drawn.
    """
    # Do nothing, and return None if locator is None
    if isinstance(locator, mticker.Locator):
        return locator
    # Decipher user input
    if locator is None:
        if time:
            locator = mticker.AutoDateLocator(*args, **kwargs)
        elif minor:
            locator = mticker.AutoMinorLocator(*args, **kwargs)
        else:
            locator = mticker.AutoLocator(*args, **kwargs)
    elif type(locator) is str: # dictionary lookup
        if locator=='logminor':
            locator = 'log'
            kwargs.update({'subs':np.arange(0,10)})
        elif locator not in locators:
            raise ValueError(f'Unknown locator "{locator}". Options are {", ".join(locators.keys())}.')
        locator = locators[locator](*args, **kwargs)
    elif isinstance(locator, Number): # scalar variable
        locator = mticker.MultipleLocator(locator, *args, **kwargs)
    else:
        locator = mticker.FixedLocator(np.sort(locator), *args, **kwargs) # not necessary
    return locator

def Formatter(formatter, *args, time=False, tickrange=None, **kwargs):
    """
    Returns a `~matplotlib.ticker.Formatter` instance.

    Parameters
    ----------
    formatter : None, float, function, list of str, or str
        If None, returns `CustomFormatter`, unless `time` is ``True`` (see below).

        If float, returns `CustomFormatter` with precision `precision`.

        If function, function is used to output label string from numeric
        input. Returns a `~matplotlib.ticker.FuncFormatter` instance.

        If list of str, labels major ticks with these strings. Returns a
        `~matplotlib.ticker.FixedFormatter` instance.

        If str, there are 3 possibilities:

            1. If string contains ``{}``, ticks will be formatted by
               calling ``string.format(number)``.
            2. If string contains ``%``, ticks will be formatted
               using the C-notation ``string % number`` method.
            3. Otherwise, a dictionary lookup is performed.

        For the dictionary lookup, options are as follows:

        =============  ============================================
        Key            Class
        =============  ============================================
        `'none'`       `~matplotlib.ticker.NullFormatter`
        `'null'`       `~matplotlib.ticker.NullFormatter`
        `'strmethod'`  `~matplotlib.ticker.StrMethodFormatter`
        `'formatstr'`  `~matplotlib.ticker.FormatStrFormatter`
        `'scalar'`     `~matplotlib.ticker.ScalarFormatter`
        `'log'`        `~matplotlib.ticker.LogFormatterSciNotation`
        `'eng'`        `~matplotlib.ticker.LogFormatterMathtext`
        `'sci'`        `~matplotlib.ticker.LogFormatterSciNotation`
        `'logit'`      `~matplotlib.ticker.LogitFormatter`
        `'eng'`        `~matplotlib.ticker.EngFormatter`
        `'percent'`    `~matplotlib.ticker.PercentFormatter`
        `'index'`      `~matplotlib.ticker.IndexFormatter`
        `'default'`    `CustomFormatter`
        `'custom'`     `CustomFormatter`
        `'proplot'`    `CustomFormatter`
        `'$'`          `DollarFormatter`
        `'dollar'`     `DollarFormatter`
        `'euro'`       `EuroFormatter`
        `'pound'`      `PoundFormatter`
        `'pi'`         `PiFormatter`
        `'e'`          `eFormatter`
        `'deg'`        `CoordinateFormatter`
        `'lat'`        `LatFormatter`, without degree symbol
        `'lon'`        `LonFormatter`, without degree symbol
        `'deglat'`     `LatFormatter`, with degree symbol
        `'deglon'`     `LonFormatter`, with degree symbol
        =============  ============================================

    time : bool, optional
        Ignored if `formatter` is not ``None``. Otherwise, if ``True``, returns
        a `~matplotlib.ticker.AutoDateFormatter` instance.

    tickrange : (float, float), optional
        See `CustomFormatter`.

    Other parameters
    ----------------
    **kwargs
        Passed to the `~matplotlib.ticker.Formatter` class on instantiation.
    """
    # Already have a formatter object
    if isinstance(formatter, mticker.Formatter): # formatter object
        return formatter
    if np.iterable(formatter) and formatter[0]=='frac':
        args.append(formatter[1]) # the number
        formatter = formatter[0]
    # Interpret user input
    if formatter is None: # by default use my special super cool formatter, better than original
        if time:
            formatter = mdates.AutoDateFormatter(*args, **kwargs)
        else:
            formatter = CustomFormatter(*args, tickrange=tickrange, **kwargs)
    elif isinstance(formatter, Number): # interpret scalar number as *precision*
        formatter = CustomFormatter(formatter, *args, tickrange=tickrange, **kwargs)
    elif isinstance(formatter, FunctionType):
        formatter = mticker.FuncFormatter(formatter, *args, **kwargs)
    elif type(formatter) is str: # assumption is list of strings
        if '{}' in formatter:
            formatter = mticker.StrMethodFormatter(formatter, *args, **kwargs) # new-style .format() form
        elif '%' in formatter:
            if time:
                formatter = mdates.DateFormatter(formatter, *args, **kwargs) # %-style, dates
            else:
                formatter = mticker.FormatStrFormatter(formatter, *args, **kwargs) # %-style, numbers
        else:
            if formatter not in formatters:
                raise ValueError(f'Unknown formatter "{formatter}". Options are {", ".join(formatters.keys())}.')
            if formatter in ['deg','deglon','deglat','lon','lat']:
                kwargs.update({'deg':('deg' in formatter)})
            formatter = formatters[formatter](*args, **kwargs)
    else:
        formatter = mticker.FixedFormatter(formatter) # list of strings on the major ticks, wherever they may be
    return formatter

#-------------------------------------------------------------------------------
# Formatting classes for mapping numbers (axis ticks) to formatted strings
# Create pseudo-class functions that actually return auto-generated formatting
# classes by passing function references to Funcformatter
#-------------------------------------------------------------------------------
# First the default formatter
def CustomFormatter(precision=6, tickrange=[-np.inf, np.inf],
                    zerotrim=True):
    r"""
    The ProPlot default tick formatter. `CustomFormatter` differs from
    the default `~matplotlib.ticker.Formatter` in the following ways:

    1. Trims trailing zeros if any exist, and switches to exponential
       notation for big numbers instead of adding that little exponential
       tag at the top of the axes.
    2. Allows user to specify *range* within which major tick marks
       are labelled.

    Parameters
    ----------
    precision : int, optional
        For `precision` :math:`p`, switch to exponential notation if number
        is :math:`<10^{-p}` or :math:`>10^{p}`.
    tickrange : (float, float), optional
        Range within which major tick marks are labelled.
    zerotrim : bool, optional
        Whether to trim trailing zeros.
    """
    # Format definition
    if tickrange is None:
        tickrange = [-np.inf, np.inf]
    elif isinstance(tickrange, Number): # use e.g. -1 for no ticks
        tickrange = [-tickrange, tickrange]
    def f(value, location):
        # Exit if not in tickrange
        eps = abs(value)/1000
        if (value + eps) < tickrange[0] or (value - eps) > tickrange[1]:
            return '' # avoid some ticks
        # Return special string
        if value==0:
            string = '0'
        else:
            power = np.log10(abs(value))
            if power > precision or power < -precision + 1:
                string = '{:.1e}'.format(value)
                if zerotrim:
                    string = re.sub(r'\.0', '', string)
            elif power < 1:
                iprecision = 1 - int(power//1)
                string = f'{{:.{iprecision:d}f}}'.format(value) # f-string compiled, then format run
                if zerotrim:
                    string = re.sub(r'\.?0$', '', string)
            else:
                string = '{:.0f}'.format(value)
        # Use unicode minus instead of ASCII hyphen (which is default)
        minus = '−' # '${-}$', u'\u002d', r'\scalebox{0.75}[1.0]{$-$}'
        string = re.sub('-', minus, string) # pure unicode minus
        return string
    # And create object
    return mticker.FuncFormatter(f)

#------------------------------------------------------------------------------#
# Formatting with prefixes
#------------------------------------------------------------------------------#
def _PrefixSuffixFormatter(*args, prefix=None, suffix=None, **kwargs):
    """
    For adding arbitrary prefix/suffix to tick labels, inside of
    minus sign.
    """
    prefix = prefix or ''
    suffix = suffix or ''
    def f(value, location):
        string = CustomFormatter(*args, **kwargs)(value, location)
        if string[0] in '−-': # unicode minus or hyphen
            string = string[0] + prefix + string[1:] + suffix
        else:
            string = prefix + string + suffix
        return string
    return mticker.FuncFormatter(f)

def DollarFormatter(*args, **kwargs):
    """
    Dollar sign in front of values.
    """
    kwargs.update({'precision':2, 'zerotrim':False})
    return _PrefixSuffixFormatter(*args, prefix='$', **kwargs)

def PoundFormatter(*args, **kwargs):
    """
    Pound sign in front of values.
    """
    kwargs.update({'precision':2, 'zerotrim':False})
    return _PrefixSuffixFormatter(*args, prefix='£', **kwargs)

def EuroFormatter(*args, **kwargs):
    """
    Euro sign in front of values.
    """
    kwargs.update({'precision':2, 'zerotrim':False})
    return _PrefixSuffixFormatter(*args, prefix='€', **kwargs)

#------------------------------------------------------------------------------#
# Formatters for dealing with one axis in geographic coordinates
#------------------------------------------------------------------------------#
def CoordinateFormatter(*args, cardinal=None, deg=True, **kwargs):
    """
    For axes corresponding to geographic coordinates.

    Parameters
    ----------
    cardinal : None or length-2 str, optional
        If str, indicates the "negative" and "positive" coordinate. For
        example, ``cardinal='SN'``.
    deg : bool, optional
        Whether to add unicode degree symbol.
    """
    # Helper function for below formatters
    def f(value, location):
        # Optional degree symbol
        suffix = ''
        if deg:
            suffix = '\N{DEGREE SIGN}' # Unicode lookup by name
        # Apply suffix if not on equator/prime meridian
        if isinstance(cardinal,str):
            if value<0:
                value *= -1
                suffix += cardinal[0]
            elif value>0:
                suffix += cardinal[1]
        # Finally use default formatter
        func = CustomFormatter(*args, **kwargs)
        return func(value, location) + suffix
    # And create object
    return mticker.FuncFormatter(f)

def LatFormatter(*args, **kwargs):
    """
    Calls `CoordinateFormatter` with ``cardinal='SN'``.
    """
    return CoordinateFormatter(*args, cardinal='SN', **kwargs)

def LonFormatter(*args, **kwargs):
    """
    Calls `CoordinateFormatter` with ``cardinal='WE'``.
    """
    return CoordinateFormatter(*args, cardinal='WE', **kwargs)

#------------------------------------------------------------------------------#
# Formatters with fractions
#------------------------------------------------------------------------------#
def FracFormatter(symbol, number):
    r"""
    Format as fractions and/or multiples of some value, e.g. a
    physical constant.

    Parameters
    ----------
    symbol : str
        The symbol, e.g. ``r'$\pi$'``.
    number : float
        The value, e.g. `numpy.pi`.
    """
    def f(n, loc): # must accept location argument
        frac = Fraction(n/number).limit_denominator()
        if n==0: # zero
            string = '0'
        elif frac.denominator==1: # denominator is one
            if frac.numerator==1:
                string = f'${symbol}$'
            elif frac.numerator==-1:
                string = f'${{-}}{symbol:s}$'
            else:
                string = f'${frac.numerator:d}{symbol:s}$'
        elif frac.numerator==1: # numerator is +/-1
            string = f'${symbol:s}/{frac.denominator:d}$'
        elif frac.numerator==-1:
            string = f'${{-}}{symbol:s}/{frac.denominator:d}$'
        else: # and again make sure we use unicode minus!
            string = f'${frac.numerator:d}{symbol:s}/{frac.denominator:d}$'
        # string = re.sub('-', '−', string) # minus will be converted to unicode version since it's inside LaTeX math
        return string
    # And create FuncFormatter class
    return mticker.FuncFormatter(f)

def PiFormatter():
    r"""
    Calls `FracFormatter` with ``symbol=r'$\pi$'``, ``number=numpy.pi``.
    """
    return FracFormatter(r'\pi', np.pi)

def eFormatter():
    r"""
    Calls `FracFormatter` with ``symbol='e'``, ``number=numpy.exp(1)``.
    """
    return FracFormatter('e', np.exp(1))

# Declare dictionaries
# Includes some custom classes, so has to go at end
scales = ['linear',
          'log',
          'symlog',
          'logit', # builtin
          'pressure',
          'height',
          'exp',
          'sine',
          'mercator',
          'inverse'] # custom
"""
List of registered scales.
"""

locators = {
    'none':        mticker.NullLocator,
    'null':        mticker.NullLocator,
    'log':         mticker.LogLocator,
    'maxn':        mticker.MaxNLocator,
    'linear':      mticker.LinearLocator,
    'log':         mticker.LogLocator,
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
"""
Mapping of strings to `~matplotlib.ticker.Locator` classes.
"""

formatters = { # note default LogFormatter uses ugly e+00 notation
    'none':      mticker.NullFormatter,
    'null':      mticker.NullFormatter,
    'strmethod': mticker.StrMethodFormatter,
    'formatstr': mticker.FormatStrFormatter,
    'scalar':    mticker.ScalarFormatter,
    'log':       mticker.LogFormatterSciNotation,
    'eng':       mticker.LogFormatterMathtext,
    'sci':       mticker.LogFormatterSciNotation,
    'logit':     mticker.LogitFormatter,
    'eng':       mticker.EngFormatter,
    'percent':   mticker.PercentFormatter,
    'index':     mticker.IndexFormatter,
    'default':   CustomFormatter,
    'custom':    CustomFormatter,
    'proplot':   CustomFormatter,
    '$':         DollarFormatter,
    'dollar':    DollarFormatter,
    'euro':      EuroFormatter,
    'pound':     PoundFormatter,
    'pi':        PiFormatter,
    'e':         eFormatter,
    'deg':       CoordinateFormatter,
    'lat':       LatFormatter,
    'lon':       LonFormatter,
    'deglat':    LatFormatter,
    'deglon':    LonFormatter,
    }
"""
Mapping of strings to `~matplotlib.ticker.Formatter` classes.
"""
