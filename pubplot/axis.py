#!/usr/bin/env python3
"""
Define various axis scales, locators, and formatters. Also define normalizers
generally used for colormap scaling. Below is rough overview of API.
General Notes:
  * Want to try to avoid **using the Formatter to scale/transform values, and
    passing the locator an array of scaled/transformed values**. Makes more sense
    to instead define separate **axis transforms**, then can use locators and
    formatters like normal, as they were intended to be used. This way, if e.g.
    matching frequency-axis with wavelength-axis, just conver the **axis limits**
    so they match, then you're good.
Scales:
  * These are complicated. See: https://matplotlib.org/_modules/matplotlib/scale.html#ScaleBase
    Use existing ones as inspiration -- e.g. InverseScale modeled after LogScale.
  * Way to think of these is that *every single value you see on an axes first
    gets secretly converted through some equation*, e.g. logarithm, and plotted
    linearly in that transformation space.
  * Special:
      - __init__() not defined on base class but *must* be defined for subclass.
  * Methods include:
      - get_transform(), which should return an mtransforms.Transform instance
      - set_default_locators_and_formatters(), which should return
        locators and formatters
      - limit_range_for_scale(), which can be used to raise errors/clip
        stuff not within that range. From Mercator example: unlike the
        autoscaling provided by the tick locators, this range limiting will
        always be adhered to, whether the axis range is set manually,
        determined automatically or changed through panning and zooming.
  * Also, have to be 'registered' unlike locators and formatters, which
    can be passed to the 'set' methods. Or maybe not?
Transforms:
  * These are complicted. See: https://matplotlib.org/_modules/matplotlib/transforms.html#Transform
  * Attributes:
      - input_dims, output_dims, is_separable, and has_inverse; the dims are because
        transforms can be N-D, but for *scales* are always 1, 1. Note is_separable is
        true if transform is separable in x/y dimensions.
  * Methods:
      - transform(): transforms N-D coordinates, given M x N array of values. Can also
        just declare transform_affine or transform_non_affine.
      - inverted(): if has_inverse True, performs inverse transform.
Locators:
  * These are complicated. See: https://matplotlib.org/_modules/matplotlib/ticker.html#Locator
  * Special:
      - __init__() not defined on base class but *must* be defined for subclass.
  * Methods include:
      - tick_values(), which accepts vmin/vmax and returns values of located ticks
      - __call__(), which can return data limits, view limits, or
        other stuff; not sure how this works or when it's invoked.
      - view_limits(), which changes the *view* limits from default vmin, vmax
        to prevent singularities (uses mtransforms.nonsingular method; for
        more info on this see: https://matplotlib.org/_modules/matplotlib/transforms.html#nonsingular)
  * Methods that usually can be left alone:
      - raise_if_exceeds(), which just tests if ticks exceed MAXTICKS number
      - autoscale(), which calls the internal locator 'view_limits' with
        result of axis.get_view_interval()
      - pan() and zoom() for interactive purposes
Formatters:
  * Easy to construct: just build with FuncFormatter a function that accepts
    the number and a 'position', which maybe is used for offset or something
    but almost always don't touch it, leave it default.
Normalizers:
  * Generally these are used for colormaps, easy to construct: require
    only an __init__ method and a __call__ method.
  * The init method takes vmin, vmax, and clip, and can define custom
    attributes. The call method just returns a *masked array* to handle NaNs,
    and call transforms data from physical units to *normalized* units
    from 0-1, representing position in colormap.
"""
#------------------------------------------------------------------------------#
# Imports
#------------------------------------------------------------------------------#
import re
import numpy as np
import numpy.ma as ma
from fractions import Fraction
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import matplotlib.scale as mscale
import matplotlib.transforms as mtransforms

#------------------------------------------------------------------------------#
# Tick scales
#------------------------------------------------------------------------------#
def ScaleFactory(l, u, scale=np.inf):
    """
    Constructer for scale with custom cutoffs. Three options here:
      1. Put a 'cliff' between two numbers (default).
      2. Accelerate the scale gradient between two numbers (scale>1).
      3. Deccelerate the scale gradient between two numbers (scale<1). So
         scale is fast on edges but slow in middle.
    TODO:
      * Alongside this, create method for drawing those cutoff diagonal marks
        with white space between.
    See: https://stackoverflow.com/a/5669301/4970632 for multi-axis solution
    and for this class-based solution. Note the space between 1-9 in Paul's answer
    is because actual cutoffs were 0.1 away (and tick locs are 0.2 apart).
    """
    class CustomScale(mscale.ScaleBase):
        # Declare name
        name = 'cutoff'
        def __init__(self, axis, **kwargs):
            mscale.ScaleBase.__init__(self)

        def get_transform(self):
            return self.CustomTransform()

        def set_default_locators_and_formatters(self, axis):
            pass # for now default can be same as for ScaleBase

        class CustomTransform(mtransforms.Transform):
            # Create transform object
            input_dims = 1
            output_dims = 1
            is_separable = True
            lower = l
            upper = u
            def __init__(self):
                mtransforms.Transform.__init__(self)
            def transform(self, a):
                a = np.array(a) # very numpy array
                aa = a.copy()
                m1 = (a>self.lower)
                m2 = (a>self.upper)
                m3 = (a>self.lower) & (a<self.upper)
                if scale==np.inf:
                    aa[m1] = a[m1] - (self.upper - self.lower)
                    aa[m3] = self.lower
                else:
                    aa[m2] = a[m2] - (self.upper - self.lower)*(1 - 1/scale)
                    aa[m3] = a[m3] - (a[m3] - self.lower)*(1 - 1/scale)
                return aa
            def transform_non_affine(self, a):
                return self.transform(a)
            def inverted(self):
                return CustomScale.InvertedCustomTransform()

        class InvertedCustomTransform(mtransforms.Transform):
            input_dims = 1
            output_dims = 1
            is_separable = True
            lower = l
            upper = u
            def __init__(self):
                mtransforms.Transform.__init__(self)
            def transform(self, a):
                a = np.array(a)
                aa = a.copy()
                n = (self.upper-self.lower)*(1 - 1/scale)
                m1 = (a>self.lower)
                m2 = (a>self.upper - n)
                m3 = (a>self.lower) & (a<(self.upper - n))
                if scale==np.inf:
                    aa[m1] = a[m1] + (self.upper - self.lower)
                else:
                    aa[m2] = a[m2] + n
                    aa[m3] = a[m3] + (a[m3] - self.lower)*(1 - 1/scale)
                return aa
            def transform_non_affine(self, a):
                return self.transform(a)
            def inverted(self):
                return CustomScale.CustomTransform()

    # Register and return
    mscale.register_scale(CustomScale)
    print('Registered scale "cutoff".')
    return CustomScale

class MercatorLatitudeScale(mscale.ScaleBase):
    """
    The scale function:
        ln(tan(y) + sec(y))
    The inverse scale function:
        atan(sinh(y))
    Applies user-defined threshold below +/-90 degrees above and below which nothing
    will be plotted. See: http://en.wikipedia.org/wiki/Mercator_projection
    """
    name = 'mercator'
    def __init__(self, axis, *, thresh=np.deg2rad(85), **kwargs):
        # Initialize
        mscale.ScaleBase.__init__(self)
        if thresh >= np.pi/2:
            raise ValueError('Threshold "thresh" must be less than pi/2.')
        self.thresh = thresh

    def get_transform(self):
        # Return special transform object
        return self.MercatorLatitudeTransform(self.thresh)

    def limit_range_for_scale(self, vmin, vmax, minpos):
        # *Hard* limit on axis boundaries
        return max(vmin, -self.thresh), min(vmax, self.thresh)

    def set_default_locators_and_formatters(self, axis):
        # Apply these
        ticks = np.deg2rad(np.arange(-90, 90, 10))
        axis.set_major_locator(mticker.FixedLocator(ticks))
        axis.set_major_formatter(self.DegreeFormatter())
        axis.set_minor_formatter(mticker.NullFormatter())

    class DegreeFormatter(mticker.Formatter):
        # For the default formatter, easy
        def __call__(self, x, pos=None):
            return f'{np.rad2deg(x):.0f}\N{DEGREE SIGN}'

    class MercatorLatitudeTransform(mtransforms.Transform):
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
            masked = ma.masked_where((a < -self.thresh) | (a > self.thresh), a)
            if masked.mask.any():
                return ma.log(np.abs(ma.tan(masked) + 1.0 / ma.cos(masked)))
            else:
                return np.log(np.abs(np.tan(a) + 1.0 / np.cos(a)))
        def inverted(self):
            # Just call inverse transform class
            return MercatorLatitudeScale.InvertedMercatorLatitudeTransform(
                self.thresh)

    class InvertedMercatorLatitudeTransform(mtransforms.Transform):
        # As above, but for the inverse transform
        input_dims = 1
        output_dims = 1
        is_separable = True
        has_inverse = True
        def __init__(self, thresh):
            mtransforms.Transform.__init__(self)
            self.thresh = thresh
        def transform_non_affine(self, a):
            return np.arctan(np.sinh(a))
        def inverted(self):
            return MercatorLatitudeScale.MercatorLatitudeTransform(self.thresh)

class SineLatitudeScale(mscale.ScaleBase):
    """
    The scale function:
        sin(rad(y))
    The inverse scale function:
        deg(arcsin(y))
    """
    name = 'sine'
    def __init__(self, axis, **kwargs):
        # Initialize
        mscale.ScaleBase.__init__(self)

    def get_transform(self):
        # Return special transform object
        return self.SineLatitudeTransform()

    def limit_range_for_scale(self, vmin, vmax, minpos):
        # *Hard* limit on axis boundaries
        return max(vmin, -90), min(vmax, 90)

    def set_default_locators_and_formatters(self, axis):
        # Apply these
        ticks = np.arange(-80, 81, 20)
        axis.set_major_locator(mticker.FixedLocator(ticks))
        axis.set_major_formatter(self.DegreeFormatter())
        axis.set_minor_formatter(mticker.NullFormatter())

    class DegreeFormatter(mticker.Formatter):
        # For the default formatter, easy
        def __call__(self, x, pos=None):
            return f'{x:.0f}\N{DEGREE SIGN}'

    class SineLatitudeTransform(mtransforms.Transform):
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
            return SineLatitudeScale.InvertedSineLatitudeTransform()

    class InvertedSineLatitudeTransform(mtransforms.Transform):
        # As above, but for the inverse transform
        input_dims = 1
        output_dims = 1
        is_separable = True
        has_inverse = True
        def __init__(self):
            mtransforms.Transform.__init__(self)
        def transform_non_affine(self, a):
            # Clipping, instead of setting invalid, fixes weird issue where
            # ylim is inappropriately zoomed way in
            aa = a.copy()
            aa[a < -1] = -1
            aa[a >  1] = 1
            return np.rad2deg(ma.arcsin(aa))
        def inverted(self):
            return MercatorLatitudeScale.SineLatitudeTransform()

class InverseScale(mscale.ScaleBase):
    """
    Similar to LogScale, but this scales to be linear in *inverse* of x. Very
    useful e.g. to plot wavelengths on twin axis with wavenumbers.
    Important note: Unlike log-scale, we can't just warp the space between
    the axis limits -- have to actually change axis limits. This scale will
    invert and swap the limits you provide. Weird! But works great!
    """
    # Declare name
    name = 'inverse'
    def __init__(self, axis, minpos=1e-2, **kwargs):
        # Initialize (note thresh is always needed)
        mscale.ScaleBase.__init__(self)
        self.minpos = minpos

    def get_transform(self):
        # Return transform class
        return self.InverseTransform(self.minpos)

    def limit_range_for_scale(self, vmin, vmax, minpos):
        # *Hard* limit on axis boundaries
        # vmin = self.minpos if vmin<=0 else vmin
        # vmax = self.minpos if vmax<=0 else vmax
        # print(vmax, vmin)
        # return 1.0/vmax, 1.0/vmin # also *reverses* the axis direction, for some reason (not sure how, because we swap vmin/vmax here)
        # Above was dumb, see LogScale example; this should set the limits
        # in *scaled coordinates* (e.g. log uses -1000, meaning 10^(-1000))
        return vmin, vmax

    def set_default_locators_and_formatters(self, axis):
        # Consider changing this
        # TODO: fix minor locator issue
        axis.set_major_locator(mticker.LogLocator(base=10, subs=[1]))
        axis.set_major_formatter(mticker.LogFormatter())
        axis.set_minor_locator(mticker.LogLocator(base=10, subs='auto'))
        axis.set_minor_formatter(mticker.NullFormatter())

    class InverseTransform(mtransforms.Transform):
        # Create transform object
        input_dims = 1
        output_dims = 1
        is_separable = True
        def __init__(self, minpos):
            mtransforms.Transform.__init__(self)
            self.minpos = minpos
        def transform(self, a, eps=1e-2):
            a = np.array(a)
            aa = a.copy()
            aa[a<=0] = eps
            # aa[a<=0] = np.nan # minpos
            return 1.0/aa
        def transform_non_affine(self, a):
            return self.transform(a)
        def inverted(self):
            return InverseScale.InvertedInverseTransform(self.minpos)

    class InvertedInverseTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True
        def __init__(self, minpos):
            mtransforms.Transform.__init__(self)
            self.minpos = minpos
        def transform(self, a, eps=1e-2):
            a = np.array(a)
            aa = a.copy()
            aa[a<=0] = eps
            # aa[a<=0] = np.nan # messes up automatic ylim setting
            return 1.0/aa
        def transform_non_affine(self, a):
            return self.transform(a)
        def inverted(self):
            return InverseScale.InverseTransform(self.minpos)

# Register hard-coded scale names, so user can set_xscale and set_yscale with strings
mscale.register_scale(InverseScale)
mscale.register_scale(SineLatitudeScale)
mscale.register_scale(MercatorLatitudeScale)

#------------------------------------------------------------------------------#
# Tick locators
#------------------------------------------------------------------------------#
def Locator(locator, *args, axis='major', **kwargs):
    """
    Default locator. Can be instantiated with a bunch of
    different objects.
    Notes:
     * Default FixedLocator includes 'nbins' option to subsample
       the points passed so that no more than 'nbins' ticks are selected.
    """
    # Declare dictionary
    locators = {'log': mtransforms.LogLocator,
        'auto':      mtransforms.AutoLocator,
        'maxn':      mtransforms.MaxNLocator,
        'linear':    mtransforms.LinearLocator,
        'log':       mtransforms.LogLocator,
        'multiple':  mtransforms.MultipleLocator,
        'fixed':     mtransforms.FixedLocator,
        'index':     mtransforms.IndexLocator,
        'null':      mtransforms.NullLocator,
        'symmetric': mtransforms.SymmetricalLogLocator,
        'logit':     mtransforms.LogitLocator,
        'autominor': mtransforms.AutoMinorLocator,
        }
    # Return instance
    if locator is None or isinstance(locator, mticker.Locator):
        # Do nothing, and return None if locator is None
        pass
    elif type(locator) is str:
        # Lookup from dictionary
        locator = locators.get(locator,None)
        if locator is None:
            raise ValueError(f'Unknown locator "{locator}". Options are {", ".join(locators.keys())}.')
    elif not hasattr(locator,'__init__'):
        # Simple multiple of this value
        locator = mticker.MultipleLocator(locator)
    else:
        # Fixed tickmarks
        locator = mticker.FixedLocator(np.sort(locator)) # not necessary
        # locator = mticker.FixedLocator(np.array(locator))
    return locator

def AutoLocate(min_, max_, base=5):
    """
    Return auto-generated levels in the provided interval. Minimum and maximum
    values will be rounded to the nearest number divisible by the specified base.
    Notes:
     * This is *outdated*. Instead should just use MultipleLocator. Not
       necessary!
    """
    _round = lambda x: base*round(float(x)/base)
    return np.arange(_round(min_), _round(max_)+base/2, base)

def LatLocator():
    """
    Wrapper around MultipleLocator, chooses points that are somewhere.
    """
    return

def LonLocator():
    """
    Wrapper around MultipleLocator, chooses points that are somewhere.
    """
    return

class InverseLocator(mticker.MultipleLocator):
    """
    Locate ticks linear in 1/x.
    """
    def __init__(self, numticks=5):
        self.numticks = numticks
    def __call__(self):
        vmin, vmax = self.axis.get_view_interval()
        ticklocs = np.reciprocal(np.linspace(1/vmax, 1/vmin, self.numticks))
        return self.raise_if_exceeds(ticklocs)
    def tick_values(vmin, vmax):
        pass

#-------------------------------------------------------------------------------
# Formatting classes for mapping numbers (axis ticks) to formatted strings
# Create pseudo-class functions that actually return auto-generated formatting
# classes by passing function references to Funcformatter
#-------------------------------------------------------------------------------
def Formatter(precision=None, tickrange=None):
    """
    Format as a number, with N sigfigs, and trimming trailing zeros.
    Recall, must pass function in terms of n (number) and loc.
    For minus sign scaling, see: https://tex.stackexchange.com/a/79158/73149
    """
    # Format definition
    def f(value, location):
        # Exit if not in tickrange
        if tickrange is not None:
            eps = abs(value)/1000
            if (value+eps)<tickrange[0] or (value-eps)>tickrange[1]:
                return '' # avoid some ticks
        # Return special string
        # * Note CANNOT use "g" because "g" precision operator specifies count of
        #   significant digits, not places after decimal place.
        # * There is no format that specifies digits after decimal place AND trims trailing zeros.
        decimals = precision or 3
        string = f'{{{0}:.{decimals:d}f}}'.format(value) # f-string compiled, then format run
        if '.' in string: # g-style trimming
            string = string.rstrip('0').rstrip('.')
        if string=='-0': # special case
            string = '0'
        if value>0 and string=='0':
            pass
            # raise RuntimeError('Tried to round tick position label to zero. Add precision or use an exponential formatter.')
        string = re.sub('-', '−', string) # pure unicode minus
        # string = re.sub('-', '${-}$', string) # latex version
        # string = re.sub('-', u'\u002d', string) # unicode hyphen minus, looks same as hyphen
        # string = re.sub('-', r'\scalebox{0.75}[1.0]{$-$}', string)
        return string
    # And create object
    return mticker.FuncFormatter(f)

def LambdaFormatter(transform, *args, **kwargs):
    """
    Arbitrary formatter to run value through some operation, the
    transform lambda function, before passing to Formatter.
    """
    # Format defintion
    def f(value, location):
        value = transform(value) # take the inverse
        formatter = Formatter(*args, **kwargs)
        return formatter(value, location)
    return mticker.FuncFormatter(f)

def InverseFormatter(*args, **kwargs):
    """
    Just take the inverse, then run through Formatter.
    """
    # Format definition
    def f(value, location):
        value = 1.0/value # take the inverse
        formatter = Formatter(*args, **kwargs)
        return formatter(value, location)
    return mticker.FuncFormatter(f)

def LatFormatter(precision=None, cardinal=False, sine=False, degree=True):
    """
    Format latitude labels; can convert sine-lats back into lats (for
    areal weighting) and can apply N/S instead of postiive/negative.
    """
    def f(value, location):
        # Convert from sine to latitude number
        if sine:
            if abs(value)>1: raise ValueError("Sine latitudes must be in range [-1,1].")
            value = np.arcsin(value)*180/np.pi
        # Optional degree symbol
        suffix = ''
        if degree:
            suffix = '\N{DEGREE SIGN}' # Unicode lookup by name
        # Suffix to apply
        if cardinal:
            if value<0:
                value *= -1
                suffix += 'S'
            elif value>0:
                suffix += 'N'
        # string = formatstr % (value, string)
        # Return special string, as in Formatter method
        decimals = precision or 3
        string = f'{{{0}:.{decimals:d}f}}'.format(value) # f-string compiled, then call format
        if '.' in string: # g-style trimming
            string = string.rstrip('0').rstrip('.')
        if string=='-0': # special case
            string = '0'
        string = re.sub('-', '−', string) # pure unicode minus
        return string + suffix
    # And create object
    return mticker.FuncFormatter(f)

def LonFormatter(precision=None, cardinal=False, degree=True):
    """
    Format latitude labels; can convert sine-lats back into lats (for
    areal weighting) and can apply N/S instead of postiive/negative.
    """
    def f(value, location):
        # Optional degree symbol
        suffix = ''
        if degree:
            suffix = '\N{DEGREE SIGN}' # Unicode lookup by name
        # Suffix to apply
        if cardinal:
            if value<0:
                value *= -1
                suffix += 'W'
            elif value>0:
                suffix += 'E'
        # string = formatstr % (value, string)
        # Return special string, as in Formatter method
        decimals = precision or 3
        string = f'{{{0}:.{decimals:d}f}}'.format(value) # f-string compiled, then call format
        if '.' in string: # g-style trimming
            string = string.rstrip('0').rstrip('.')
        if string=='-0': # special case
            string = '0'
        string = re.sub('-', '−', string) # pure unicode minus
        return string + suffix
    # And create object
    return mticker.FuncFormatter(f)

def PiFormatter(number=np.pi, symbol=r'\pi'):
    """
    Format as fractions, multiples of some value.
    Note: Since everything is put inside LaTeX math mode the hyphens
    should be converted to unicode minus.
    """
    def f(n, loc): # must accept location argument
        frac = Fraction(n/number).limit_denominator()
        minus = '−' # pure unicode minus
        if n==0: # zero
            return '0'
        elif frac.denominator==1: # denominator is one
            if frac.numerator==1:
                return f'${symbol}$'
            elif frac.numerator==-1:
                return f'${{-}}{symbol:s}$'
            else:
                return f'${frac.numerator:d}{symbol:s}$'
        elif frac.numerator==1: # numerator is +/-1
            return f'${symbol:s}/{frac.denominator:d}$'
        elif frac.numerator==-1:
            return f'${{-}}{symbol:s}/{frac.denominator:d}$'
        else: # and again make sure we use unicode minus!
            return f'${frac.numerator:d}{symbol:s}/{frac.denominator:d}$'
    # And create FuncFormatter class
    return mticker.FuncFormatter(f)

#------------------------------------------------------------------------------
# Normalization classes for mapping data to colors (i.e. colormaps)
#------------------------------------------------------------------------------
class Norm(mcolors.Normalize):
    """
    Like the default BoundaryNorm, except instead of drawing color 'x' directly
    from RGB array, out of n colors for n+1 levels, this *interpolates* from an RGB
    array at point '(level_index+0.5)/num_levels' along the array.
    Example: Your levels edges are weirdly spaced [-1000, 100, 0, 100, 1000] or
    even [0, 10, 12, 20, 22], but center "colors" are always at colormap
    coordinates [.2, .4, .6, .8] no matter the spacing; levels just must be monotonic.
    """
    def __init__(self, levels, midpoint=None, clip=False, **kwargs):
        # Very simple
        try: iter(levels)
        except TypeError:
            raise ValueError("Must call Norm with your boundary vaues.")
        # mcolors.Normalize.__init__(self, min(levels), max(levels), clip, **kwargs)
        super().__init__(min(levels), max(levels), clip)
        self.midpoint = midpoint
        self.levels = np.array(levels)

    def __call__(self, value, clip=None):
        # TOTO: Add optional midpoint; this class will probably end up being one of
        # my most used if so; midpoint would just ensure <value> corresponds to 0.5 in cmap
        # Some checks (maybe not necessary)
        try: iter(value)
        except TypeError:
            value = np.array([value])
        if value.ndim>1:
            raise ValueError("Array is multi-dimensional... not sure what to do.")
        # Map data values in range (vmin,vmax) to color indices in colormap
        nvalues = np.empty(value.shape)
        for i,v in enumerate(value.flat):
            if np.isnan(v):
                continue
            locs = np.where(v>=self.levels)[0]
            if locs.size==0:
                nvalues[i] = 0
            elif locs.size==self.levels.size:
                nvalues[i] = 1
            else:
                # Define the boundary level integers
                # So if 11 levels and between ids 9 and 10, interpolant is between .9 and 1
                interpolee = self.levels[[locs[-1],locs[-1]+1]] # the boundary level values
                interpolant = np.array([locs[-1],locs[-1]+1])/(self.levels.size-1)
                nvalues[i] = np.interp(v, interpolee, interpolant)
                # nvalues[i] = max(0, nvalues[i]-.5/(self.levels.size-1))
                # print(self.vmin, self.vmax, min(self.levels), max(self.levels))
        return ma.masked_array(nvalues, np.isnan(value))

class StretchNorm(mcolors.Normalize):
    """
    Class that can 'stretch' and 'compress' either side of a colormap about
    some midpoint; proceeds exponentially (exp>0) or logarithmically (exp<0)
    down the linear colormap from the center point. Default midpoint is vmin, i.e.
    we just stretch to the right. For diverging colormaps, use midpoint 0.5.
    """
    def __init__(self, exp=0, extend='neither', midpoint=None, vmin=None, vmax=None, clip=None):
        # User will use -10 to 10 scale; converted to value used in equation
        if abs(exp) > 10: raise ValueError('Warping scale must be between -10 and 10.')
        super().__init__(vmin, vmax, clip)
        self.midpoint = midpoint
        self.exp = exp
        self.extend = extend
        # mcolors.Normalize.__init__(self, vmin, vmax, clip)

    # Function
    def warp(x, exp, exp_max=4):
        # Returns indices stretched so neutral/low values are sampled more heavily
        # Will artifically use exp to signify stretching away from neutral vals,
        # or compressing toward neutral vals
        if exp > 0:
            invert = True
        else:
            invert, exp = False, -exp
        exp = exp*(exp_max/10)
        # Apply function; approaches x=1 as a-->Inf, x=x as a-->0
        if invert: x = 1-x
        value =  (x-1+(np.exp(x)-x)**exp)/(np.e-1)**exp
        if invert:
            value = 1-value # flip on y-axis
        return value

    def __call__(self, value, clip=None):
        # Initial stuff
        if self.midpoint is None:
            midpoint = self.vmin
        else:
            midpoint = self.midpoint
        # Get middle point in 0-1 coords, and value
        midpoint_scaled = (midpoint - self.vmin)/(self.vmax - self.vmin)
        value_scaled    = (value - self.vmin)/(self.vmax - self.vmin)
        try: iter(value_scaled)
        except TypeError:
            value_scaled = np.arange(value_scaled)
        value_cmap = ma.empty(value_scaled.size)
        for i,v in enumerate(value_scaled):
            # Get values, accounting for midpoints
            if v < 0: v = 0
            if v > 1: v = 1
            if v >= midpoint_scaled:
                block_width = 1 - midpoint_scaled
                value_cmap[i] = (midpoint_scaled + 
                    block_width*warp((v - midpoint_scaled)/block_width, self.exp)
                    )
            else:
                block_width = midpoint_scaled
                value_cmap[i] = (midpoint_scaled - 
                        block_width*warp((midpoint_scaled - v)/block_width, self.exp)
                        )
        if self.extend=='both' or self.extend=='max':
            value_cmap[value_cmap>1] = 1
        if self.extend=='both' or self.extend=='min':
            value_cmap[value_cmap<0] = 0
        return value_cmap

class MidpointNorm(mcolors.Normalize):
    """
    Pass as norm=<instance>, when declaring new pcolor or contourf objects.
    Creates new normalization of existing registered cmap by changing midpoint
    away from (vmin+vmax)/2.
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        # Default None values so can declare by name in any order
        # mcolors.Normalize.__init__(self, vmin, vmax, clip)
        super().__init__(vmin, vmax, clip)
        self.midpoint = midpoint
        if any(x is None for x in [vmin,vmax,midpoint]):
            raise ValueError("Must declare vmin, vmax, and midpoint explicitly.")
        if self.midpoint<self.vmin or self.midpoint>self.vmax:
            raise ValueError('Midpoint {self.midpoint} is not between vmin {self.vmin} and vmax {self.vmax}.')

    def __call__(self, value, clip=None):
        # How to map data values in range (vmin,vmax) to color indices in colormap
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return ma.masked_array(np.interp(value, x, y), np.isnan(value))

