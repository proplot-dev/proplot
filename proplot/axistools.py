#!/usr/bin/env python3
"""
Defines various axis scales, locators, and formatters. Also "registers"
the locator and formatter names, so that they can be called selected with
the `~proplot.axes.XYAxes.format` method.

This was pretty tricky. You can find my notes on overriding the
`~matplotlib.scale.ScaleBase`, `~matplotlib.transforms.Transform`,
`~matplotlib.ticker.Locator`, and
`~matplotlib.ticker.Formatter` classes below.

.. raw:: html

   <h1>Scales</h1>


* These are complicated. See `~matplotlib.scale.ScaleBase`. Use existing ones
  as inspiration -- e.g. `InverseScale` modeled after `~matplotlib.scale.LogScale`.
* Way to think of these is that *every single value you see on an axes first
  gets secretly converted through some equation*, e.g. logarithm, and plotted
  linearly in that transformation space.
* Methods:

   - `get_transform`: Returns a `~matplotlib.transforms.Transform` instance.
   - `set_default_locators_and_formatters`: Returns
     the default locators and formatters.
   - `limit_range_for_scale`: Can be used to raise errors or clip
     stuff not within that range.

     From the Mercator example: unlike the autoscaling provided by the
     tick locators, this range limiting will always be adhered to, whether
     the axis range is set manually, determined automatically or changed
     through panning and zooming.

* Notes on methods:

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

   <h1>Transforms</h1>


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

   <h1>Locators</h1>


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
      result of `~matplotlib.axis.Axis.get_view_interval`.
    - `pan` and `zoom`: Interactive purposes.

.. raw:: html

   <h1>Formatters</h1>


Some of these are easy to construct. Just pass `~matplotlib.ticker.FuncFormatter`
a function that accepts two arguments the number and 'position', which maybe is used
for offset or something (don't touch it, leave it default).

The matplotlib default `~matplotlib.ticker.ScalarFormatter` is much more
complex, but can be overridden in the typical way: adding stuff to the
`__init__` and `__call__` methods.
"""
#------------------------------------------------------------------------------#
# Imports
#------------------------------------------------------------------------------#
import re
from .utils import ic
from .rcmod import rc
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
def Scale(scale, *args, **kwargs):
    """
    Returns a `~matplotlib.scale.ScaleBase` instance.

    Parameters
    ----------
    scale : str or (str, *args)
        If tuple, ``args`` are passed to the class
        instantiator. Ignored unless the str is ``'cutoff'``,
        ``'exp'``, ``'height'``, or ``'pressure'``.

        If str, a dictionary lookup is performed. Options are as follows:

        ==============  =======================================  =========================================================
        Key             Class                                    Description
        ==============  =======================================  =========================================================
        ``'linear'``    `~matplotlib.scale.LinearScale`          Linear
        ``'log'``       `~matplotlib.scale.LogScale`             Logarithmic
        ``'symlog'``    `~matplotlib.scale.SymmetricalLogScale`  Logarithmic beyond space around zero
        ``'logit'``     `~matplotlib.scale.LogitScale`           Logistic
        ``'pressure'``  `ExpScale`                               Scale height coords to be linear in pressure.
        ``'height'``    `ExpScale`                               Scale pressure coords to be linear in height.
        ``'exp'``       `ExpScale`                               Scale with some exponential function
        ``'sine'``      `SineLatitudeScale`                      Scale with sine function (in degrees)
        ``'mercator'``  `MercatorLatitudeScale`                  Scale with Mercator latitude projection coords
        ``'inverse'``   `InverseScale`                           Scale with the inverse
        ==============  =======================================  =========================================================

    Other parameters
    ----------------
    **kwargs
        Passed to the `~matplotlib.scale.ScaleBase` class on instantiation.
    """
    # Existing scales
    args = []
    if isinstance(scale, mscale.ScaleBase):
        scale = scale.name
    if np.iterable(scale) and not isinstance(scale, str):
        scale, args = scale[0], [*scale[1:], *args]
    if scale in scales:
        return scale # already registered
    # Build an on-the-fly scale
    if scale=='cutoff':
        scale = CutoffScaleFactory(*args, **kwargs)
    elif scale in ('exp', 'height', 'pressure'): # note here args is non-zero
        scale = ExpScaleFactory(*args, to_exp=(scale!='pressure'), **kwargs)
    else:
        raise ValueError(f'Unknown scale {scale}. Options are {", ".join(scales.keys())}.')
    return scale

def InvertedScaleFactory(scale, **kwargs):
    """Returns name of newly registered *inverted* version of the
    `~matplotlib.scale.ScaleBase` corresponding to the scale name."""
    # Note I use private API here; no other way to get class directly
    scale = Scale(scale, **kwargs) # get the class
    name_ = f'{scale}_inverted' # name of inverted version
    class Inverted(scales[scale]):
        name = name_
        def get_transform(self):
            return super().get_transform().inverted() # that's all we need!
    mscale.register_scale(Inverted)
    return name_

#------------------------------------------------------------------------------#
# Exp axis
# Why can't this just be a class that accepts args for changing
# the scale params? Because then would need kwargs for every set_xscale
# call.
#------------------------------------------------------------------------------#
def ExpScaleFactory(expscale, mulscale, to_exp=True, name=None):
    r"""
    Exponential scale.

    This is used when adding a pressure coordinate axis
    for data that is plotted linearly w.r.t. height, or vice versa. Ignore
    this if you're not an atmospheric scientist.

    Parameters
    ----------
    expscale : float
        The scale for the exponent, i.e. the :math:`b` in :math:`Ae^{bx}`.
    mulscale : float
        The coefficient, i.e. the :math:`A` in :math:`Ae^{bx}`.
    to_exp : float
        Whether the "forward" direction performs exponentiation, or
        the inverse operation :math:`\log(x/A)/b`.
    name : None or str, optional
        The scale name. Defaults to ``'exp_{expscale}_{mulscale}'``.
    """
    if name is None:
        name = f'exp_{expscale:.1e}_{mulscale:.1e}'
    name_ = name # must make a copy
    class ExpScale(mscale.ScaleBase):
        # Declare name
        name = name_
        def __init__(self, axis, minpos=1e-300, **kwargs):
            # Initialize
            super().__init__()
            if to_exp:
                transform = _ExpTransform(expscale, mulscale, minpos)
            else:
                transform = _InvertedExpTransform(expscale, mulscale, minpos)
            self._transform = transform
        def limit_range_for_scale(self, vmin, vmax, minpos):
            # Do nothing
            return vmin, vmax
        def set_default_locators_and_formatters(self, axis):
            # Apply log spacing by default
            axis.set_smart_bounds(True) # unnecessary?
            axis.set_major_formatter(Formatter('default'))
            axis.set_minor_formatter(Formatter('null'))
        def get_transform(self):
            # Either sub into e(scale*z), the default, or invert
            # the exponential
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
    def __init__(self, expscale, mulscale, minpos):
        super().__init__()
        self.minpos = minpos
        self._expscale = expscale
        self._mulscale = mulscale
    def transform(self, a):
        return self._mulscale*np.power(np.e, self._expscale*np.array(a))
    def transform_non_affine(self, a):
        return self.transform(a)
    def inverted(self):
        return _InvertedExpTransform(self._expscale, self._mulscale, self.minpos)

class _InvertedExpTransform(mtransforms.Transform):
    """Inverse exponential coordinate transform."""
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True
    def __init__(self, expscale, mulscale, minpos):
        super().__init__()
        self.minpos = minpos
        self._expscale = expscale
        self._mulscale = mulscale
    def transform(self, a):
        aa = np.array(a).copy()
        aa[aa<=self.minpos] = self.minpos # necessary
        return (np.log(aa) - np.log(self._mulscale))/self._expscale # this one!
    def transform_non_affine(self, a):
        return self.transform(a)
    def inverted(self):
        return _ExpTransform(self._expscale, self._mulscale, self.minpos)

#------------------------------------------------------------------------------#
# Cutoff axis
#------------------------------------------------------------------------------#
def CutoffScaleFactory(scale, lower, upper=None):
    """
    Construct a scale with custom cutoffs.

    If `upper` is ``None``, you have the following two possibilities:

    1. If `scale` is >1, the axis is "faster" to the right
       of `lower`.
    2. If `scale` is <1, the axis is "slower" to the right
       of `lower`.

    If `upper` is not ``None``, you have the following three possibilities:

    1. If `scale` is >1, the axis is accelerated between `lower` and `upper`.
       So the axis is "slow" on the edges but "fast" in middle.
    2. If `scale` is >1, the axis is deccelerated between `lower` and `upper`.
       So the axis is "fast" on the edges but "slow" in middle.
    3. If `scale` is `numpy.inf`, this puts a cliff between `lower` and
       `upper`. The axis cuts from `lower` to `upper` at a single point.

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
    if upper is None:
        if scale==np.inf:
            raise ValueError('For infinite scale (i.e. discrete cutoff), need both lower and upper bounds.')
    name = f'cutoff_{scale:.1e}_{lower:.1e}'
    if upper is not None:
        name = f'{name}_{upper:.1e}'
    name_ = name # have to copy to different name

    class CutoffScale(mscale.ScaleBase):
        # Declare name
        name = name_
        def __init__(self, axis, **kwargs):
            super().__init__()
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
        super().__init__()
        if thresh >= 90.0:
            raise ValueError('Threshold "thresh" must be <=90.')
        self.thresh = thresh
    def get_transform(self):
        """See `~matplotlib.scale.ScaleBase`."""
        # Return special transform object
        return _MercatorLatitudeTransform(self.thresh)
    def limit_range_for_scale(self, vmin, vmax, minpos):
        """See `~matplotlib.scale.ScaleBase`."""
        # *Hard* limit on axis boundaries
        return max(vmin, -self.thresh), min(vmax, self.thresh)
    def set_default_locators_and_formatters(self, axis):
        """See `~matplotlib.scale.ScaleBase`."""
        # Apply these
        axis.set_smart_bounds(True)
        axis.set_major_locator(Locator(20)) # every 20 degrees
        axis.set_major_formatter(Formatter('deg'))
        axis.set_minor_formatter(Formatter('null'))

class _MercatorLatitudeTransform(mtransforms.Transform):
    """Mercator latitude coordinate transform."""
    # Default attributes
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True
    def __init__(self, thresh):
        # Initialize, declare attribute
        super().__init__()
        self.thresh = thresh
    def transform_non_affine(self, a):
        """See `~matplotlib.transforms.Transform`."""
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
        """See `~matplotlib.transforms.Transform`."""
        # Just call inverse transform class
        return _InvertedMercatorLatitudeTransform(self.thresh)

class _InvertedMercatorLatitudeTransform(mtransforms.Transform):
    """Inverse Mercator latitude coordinate transform."""
    # As above, but for the inverse transform
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True
    def __init__(self, thresh):
        super().__init__()
        self.thresh = thresh
    def transform_non_affine(self, a):
        """See `~matplotlib.transforms.Transform`."""
        # m = ma.masked_where((a < -self.thresh) | (a > self.thresh), a)
        return np.degrees(np.arctan2(1, np.sinh(a))) # always assume in first/fourth quadrant, i.e. go from -pi/2 to pi/2
    def inverted(self):
        """See `~matplotlib.transforms.Transform`."""
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
        super().__init__()
    def get_transform(self):
        """See `~matplotlib.scale.ScaleBase`."""
        return _SineLatitudeTransform()
    def limit_range_for_scale(self, vmin, vmax, minpos):
        """See `~matplotlib.scale.ScaleBase`."""
        # *Hard* limit on axis boundaries
        # return max(vmin, -90), min(vmax, 90)
        return vmin, vmax
    def set_default_locators_and_formatters(self, axis):
        """See `~matplotlib.scale.ScaleBase`."""
        axis.set_smart_bounds(True)
        axis.set_major_locator(Locator(20)) # every 20 degrees
        axis.set_major_formatter(Formatter('deg'))
        axis.set_minor_formatter(Formatter('null'))

class _SineLatitudeTransform(mtransforms.Transform):
    """Sine latitude transform."""
    # Default attributes
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True
    def __init__(self):
        # Initialize, declare attribute
        super().__init__()
    def transform_non_affine(self, a):
        """See `~matplotlib.transforms.Transform`."""
        with np.errstate(invalid='ignore'): # NaNs will always be False
            m = (a >= -90) & (a <= 90)
        if not m.all():
            aa = ma.masked_where(~m, a)
            return ma.sin(np.deg2rad(aa))
        else:
            return np.sin(np.deg2rad(a))
    def inverted(self):
        """See `~matplotlib.transforms.Transform`."""
        return _InvertedSineLatitudeTransform()

class _InvertedSineLatitudeTransform(mtransforms.Transform):
    """Inverse sine latitude transform."""
    # Inverse of _SineLatitudeTransform
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True
    def __init__(self):
        super().__init__()
    def transform_non_affine(self, a):
        """See `~matplotlib.transforms.Transform`."""
        # Clipping, instead of setting invalid
        # NOTE: Using ma.arcsin below caused super weird errors, dun do that
        aa = a.copy()
        return np.rad2deg(np.arcsin(aa))
    def inverted(self):
        """See `~matplotlib.transforms.Transform`."""
        return _SineLatitudeTransform()

#------------------------------------------------------------------------------#
# Scale as the *inverse* of the axis coordinate
#------------------------------------------------------------------------------#
class InverseScale(mscale.ScaleBase):
    """Similar to `~matplotlib.scale.LogScale`, but this scales to be linear
    in the *inverse* of *x*. Very useful e.g. for plotting wavelengths
    on twin axis with wavenumbers."""
    # Developer notes:
    # Unlike log-scale, we can't just warp the space between
    # the axis limits -- have to actually change axis limits. Also this
    # scale will invert and swap the limits you provide. Weird! But works great!
    # Declare name
    name = 'inverse'
    """Registered scale name."""
    def __init__(self, axis, minpos=1e-2, **kwargs):
        # Initialize
        super().__init__()
        self.minpos = minpos
    def get_transform(self):
        """See `~matplotlib.scale.ScaleBase`."""
        # Return transform class
        return _InverseTransform(self.minpos)
    def limit_range_for_scale(self, vmin, vmax, minpos):
        """See `~matplotlib.scale.ScaleBase`."""
        # *Hard* limit on axis boundaries
        if not np.isfinite(minpos):
            minpos = 1e-300
        return (minpos if vmin <= 0 else vmin,
                minpos if vmax <= 0 else vmax)
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
    """Inverse transform."""
    # Create transform object
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True
    def __init__(self, minpos):
        super().__init__()
        self.minpos = minpos
    def transform(self, a):
        """See `~matplotlib.scale.ScaleBase`."""
        a = np.array(a)
        aa = a.copy()
        # f = np.abs(a)<=self.minpos # attempt for negative-friendly
        # aa[f] = np.sign(a[f])*self.minpos
        aa[aa<=self.minpos] = self.minpos
        return 1.0/aa
    def transform_non_affine(self, a):
        """See `~matplotlib.scale.ScaleBase`."""
        return self.transform(a)
    def inverted(self):
        """See `~matplotlib.scale.ScaleBase`."""
        return _InverseTransform(self.minpos)

# Register hard-coded scale names, so user can set_xscale and set_yscale
# with strings
mscale.register_scale(InverseScale)
mscale.register_scale(SineLatitudeScale)
mscale.register_scale(MercatorLatitudeScale)
ExpScaleFactory(-1.0/7, 1013.25, False, 'height')   # scale pressure so it matches a height axis
ExpScaleFactory(-1.0/7, 1013.25, True,  'pressure') # scale height so it matches a pressure axis

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
        If ``None`` or ``'default'``, returns the default
        `~matplotlib.ticker.AutoLocator` unless `minor` or `time` are
        ``True`` (see below).

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
        Ignored if `locator` is not ``None`` or ``'default'``. Otherwise,
        if ``True``, returns an `~matplotlib.ticker.AutoMinorLocator` instance.
    time : bool, optional
        Ignored if `locator` is not ``None`` or ``'default'``. Otherwise,
        if ``True``, returns an `~matplotlib.dates.AutoDateLocator` instance.

    Other parameters
    ----------------
    *args, **kwargs
        Passed to the `~matplotlib.ticker.Locator` class on instantiation.

    Note
    ----
    `~matplotlib.ticker.AutoLocator` has a useful ``nbins`` option. This
    limits the maximum number of ticks to ``nbins-1``.
    """
    # Do nothing, and return None if locator is None
    if isinstance(locator, mticker.Locator):
        return locator
    # Decipher user input
    if np.iterable(locator) and not isinstance(locator, str) and not all(isinstance(num, Number) for num in locator):
        locator, args = locator[0], [*locator[1:], *args]
    if locator is None or (isinstance(locator, str) and locator=='default'):
        if time:
            locator = mdates.AutoDateLocator(*args, **kwargs)
        elif minor:
            locator = mticker.AutoMinorLocator(*args, **kwargs)
        else:
            locator = mticker.AutoLocator(*args, **kwargs)
    elif isinstance(locator, str): # dictionary lookup
        if locator=='logminor':
            locator = 'log'
            kwargs.update({'subs':np.arange(0,10)})
        elif locator not in locators:
            raise ValueError(f'Unknown locator "{locator}". Options are {", ".join(locators.keys())}.')
        if locator=='index':
            if not args:
                args = [1] # step
            if len(args)==1:
                args += [0]
        locator = locators[locator](*args, **kwargs)
    elif isinstance(locator, Number): # scalar variable
        locator = mticker.MultipleLocator(locator, *args, **kwargs)
    else:
        locator = mticker.FixedLocator(np.sort(locator), *args, **kwargs) # not necessary
    return locator

def Formatter(formatter, *args, time=False, tickrange=None, **kwargs):
    r"""
    Returns a `~matplotlib.ticker.Formatter` instance.

    Parameters
    ----------
    formatter : float, function, list of str, or str
        If ``None`` or ``'default'``, returns `ScalarFormatter` unless
        `time` is ``True`` (see below).

        If function, function is used to output label string from numeric
        input. Returns a `~matplotlib.ticker.FuncFormatter` instance.

        If list of str, labels major ticks with these strings. Returns a
        `~matplotlib.ticker.FixedFormatter` instance.

        If str, there are 3 possibilities:

            1. If string contains ``%``, ticks will be formatted
               using the C-notation ``string % number`` method. See `this link
               <https://docs.python.org/3.4/library/string.html#format-specification-mini-language>`__
               for a review.

               If the axis represents time (i.e. ``time=True``,
               passed automatically by `~proplot.axes.XYAxes.smart_update`),
               datetime %-formatting is used. See `this page
               <https://docs.python.org/3/library/datetime.html#strftime-strptime-behavior>`__
               for an overview.
            2. If string contains ``{x}`` or ``{x:...}``, ticks will be
               formatted by calling ``string.format(x=number)``.
            3. Otherwise, a dictionary lookup is performed.

        For the dictionary lookup, options are as follows:

        =============  ====================================================================
        Key            Class
        =============  ====================================================================
        `'default'`     `ScalarFormatter`
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
        `'pi'`         `FracFormatter`, with symbol :math:`\pi` and value `numpy.pi`
        `'e'`          `FracFormatter`, with symbol *e* and value `numpy.e`
        `'deg'`        `CoordFormatter`, with just a degree symbol
        `'lat'`        `CoordFormatter`, cardinal "SN" indicator without degree symbol
        `'lon'`        `CoordFormatter`, cardinal "WE" indicator without degree symbol
        `'deglat'`     `CoordFormatter`, cardinal "SN" indicator with degree symbol
        `'deglon'`     `CoordFormatter`, cardinal "WE" indicator with degree symbol
        =============  ====================================================================

    time : bool, optional
        If ``True``, returns an `~matplotlib.dates.AutoDateFormatter` when
        `formatter` is ``None``, or a `~matplotlib.dates.DateFormatter`
        instance when `formatter` is a string with a ``'%'`` character.
    tickrange : (float, float), optional
        See `ScalarFormatter`.

    Other parameters
    ----------------
    *args, **kwargs
        Passed to the `~matplotlib.ticker.Formatter` class on instantiation.
    """
    # Already have a formatter object
    if isinstance(formatter, mticker.Formatter): # formatter object
        return formatter
    if np.iterable(formatter) and not isinstance(formatter, str) and not all(isinstance(item, str) for item in formatter):
        formatter, args = formatter[0], [*formatter[1:], *args]
    # Interpret user input
    if formatter is None or (isinstance(formatter, str) and formatter=='default'): # by default use my special super cool formatter, better than original
        if not time:
            formatter = ScalarFormatter(*args, tickrange=tickrange, **kwargs)
        else:
            formatter = mdates.AutoDateFormatter(*args, **kwargs)
    elif isinstance(formatter, FunctionType):
        formatter = mticker.FuncFormatter(formatter, *args, **kwargs)
    elif isinstance(formatter, str): # assumption is list of strings
        if re.search(r'{x(:.+)?}', formatter):
            formatter = mticker.StrMethodFormatter(formatter, *args, **kwargs) # new-style .format() form
        elif '%' in formatter:
            if time:
                formatter = mdates.DateFormatter(formatter, *args, **kwargs) # %-style, dates
            else:
                formatter = mticker.FormatStrFormatter(formatter, *args, **kwargs) # %-style, numbers
        else:
            # Coordinate versions
            if formatter in ['deg', 'deglon', 'deglat', 'lon', 'lat']:
                kwargs.update({
                    'degrees':  ('deg' in formatter),
                    'cardinal': (None if formatter=='deg' else 'SN' if 'lat' in formatter else 'WE')
                    })
                formatter = 'deg'
            # Fractions
            if formatter=='pi':
                kwargs.update({'symbol': r'\pi', 'number': np.pi})
                formatter = 'frac'
            if formatter=='e':
                kwargs.update({'symbol': 'e', 'number': np.e})
                formatter = 'frac'
            if formatter not in formatters:
                raise ValueError(f'Unknown formatter "{formatter}". Options are {", ".join(formatters.keys())}.')
            formatter = formatters[formatter](*args, **kwargs)
    elif np.iterable(formatter):
        formatter = mticker.FixedFormatter(formatter) # list of strings on the major ticks, wherever they may be
    else:
        raise ValueError(f'Invalid formatter "{formatter}".')
    return formatter

#-------------------------------------------------------------------------------
# Formatting classes for mapping numbers (axis ticks) to formatted strings
# Create pseudo-class functions that actually return auto-generated formatting
# classes by passing function references to Funcformatter
#-------------------------------------------------------------------------------
# First the default formatter
class ScalarFormatter(mticker.ScalarFormatter):
    r"""
    The new default formatter, a simple wrapper around
    `~matplotlib.ticker.ScalarFormatter`. Differs from
    `~matplotlib.ticker.ScalarFormatter` in the following ways:

        1. Trims trailing zeros if any exist.
        2. Allows user to specify *range* within which major tick marks
           are labelled.
        3. Allows user to add arbitrary prefix or suffix to every
           tick label string.

    """
    def __init__(self, *args, zerotrim=None, tickrange=None,
                              prefix=None, suffix=None, **kwargs):
        """
        Parameters
        ----------
        tickrange : None or (float, float), optional
            Range within which major tick marks are labelled.
        zerotrim : bool, optional
            Whether to trim trailing zeros.
        prefix, suffix : None or str, optional
            Optional prefix and suffix for strings.
        *args, **kwargs
            Passed to `matplotlib.ticker.ScalarFormatter`.
        """
        tickrange = tickrange or (-np.inf, np.inf)
        super().__init__(*args, **kwargs)
        if zerotrim is None:
            zerotrim = rc.get('axes.formatter.zerotrim')
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
        pos : None or float, optional
            The position.
        """
        # Tick range limitation
        eps = abs(x)/1000
        tickrange = self._tickrange
        if (x + eps) < tickrange[0] or (x - eps) > tickrange[1]:
            return '' # avoid some ticks
        # Normal formatting
        string = super().__call__(x, pos)
        if x!=0 and string=='0': # happens sometimes with log locator and scalar formatter
            string = f'{x:.6f}' # TODO: more robust fix?
        if self._zerotrim:
            string = re.sub(r'\.0+$', '', string)
            string = re.sub(r'^(.*\..*?)0+$', r'\1', string)
        string = re.sub(r'^[−-]0$', '0', string) # '-0' to '0'; necessary?
        # Prefix and suffix
        prefix = ''
        if string and string[0] in '−-': # unicode minus or hyphen
            prefix, string = string[0], string[1:]
        return prefix + self._prefix + string + self._suffix

#------------------------------------------------------------------------------#
# Formatters for dealing with one axis in geographic coordinates
#------------------------------------------------------------------------------#
def CoordFormatter(*args, cardinal=None, degrees=True, **kwargs):
    """
    Returns a `~matplotlib.ticker.FuncFormatter` that formats numbers as
    geographic coordinates.

    Parameters
    ----------
    cardinal : None or length-2 str, optional
        If str, indicates the "negative" and "positive" coordinate to append
        to the end of the tick label string. For example, ``cardinal='SN'``.
    degrees : bool, optional
        Whether to add the unicode "degree sign" symbol.
    """
    def f(x, pos):
        # Optional degree symbol
        suffix = ''
        if degrees:
            suffix = '\N{DEGREE SIGN}' # Unicode lookup by name
        # Apply suffix if not on equator/prime meridian
        if cardinal:
            if x < 0:
                x *= -1
                suffix += cardinal[0]
            elif x > 0:
                suffix += cardinal[1]
        # Finally use default formatter
        string = '{:.6f}'.format(x).replace('-', '\N{MINUS SIGN}')
        string = re.sub(r'\.0+$', '', string)
        string = re.sub(r'^(.*\..*?)0+$', r'\1', string) # note the non-greedy secondary glob!
        return string + suffix
    return mticker.FuncFormatter(f)

#------------------------------------------------------------------------------#
# Formatters with fractions
#------------------------------------------------------------------------------#
def FracFormatter(symbol, number):
    r"""
    Returns a `~matplotlib.ticker.FuncFormatter` that formats numbers as
    fractions or multiples of some value, e.g. a physical constant.

    This is powerd by the python builtin `~fractions.Fraction` class.
    We account for floating point errors using the
    `~fractions.Fraction.limit_denominator` method.

    Parameters
    ----------
    symbol : str
        The symbol, e.g. ``r'$\pi$'``.
    number : float
        The value, e.g. `numpy.pi`.
    """
    def f(x, pos): # must accept location argument
        frac = Fraction(x/number).limit_denominator()
        if x==0: # zero
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

# Declare dictionaries
# Includes some custom classes, so has to go at end
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
"""Mapping of strings to `~matplotlib.ticker.Locator` classes. See
`Locator` for a table."""

formatters = { # note default LogFormatter uses ugly e+00 notation
    'default':   ScalarFormatter,
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
    'frac':      FracFormatter,
    'deg':       CoordFormatter,
    }
"""Mapping of strings to `~matplotlib.ticker.Formatter` classes. See
`Formatter` for a table."""
