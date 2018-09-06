#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# Imports
#------------------------------------------------------------------------------#
import re
import numpy as np
from fractions import Fraction
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker

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
                interpolee = self.levels[[locs[-1],locs[-1]+1]] # the boundary level values
                interpolant = np.array([locs[-1],locs[-1]+1])/(self.levels.size-1) # the boundary level integers
                    # so if 11 levels and between ids 9 and 10, interpolant is between .9 and 1
                nvalues[i] = np.interp(v, interpolee, interpolant)
                # nvalues[i] = max(0, nvalues[i]-.5/(self.levels.size-1))
                # print(self.vmin, self.vmax, min(self.levels), max(self.levels))
        return np.ma.masked_array(nvalues, np.isnan(value))

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

    def __call__(self, value, clip=None):
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
            if invert: value = 1-value # flip on y-axis
            return value
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
        value_cmap = np.ma.empty(value_scaled.size)
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

    def __call__(self, value, clip=None):
        # How to map data values in range (vmin,vmax) to color indices in colormap
        # if self.midpoint<self.vmin or self.midpoint>self.vmax:
        #     raise ValueError("Midpoint {self.midpoint} is not between"\
        #             "vmin {self.vmin} and vmax {self.vmax}.")
        # print(x, np.interp(value, x, y))
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

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
    precision = 3 if precision is None else precision
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
        # print(string)
        string = f'{{{0}:.{precision:d}f}}'.format(value) # f-string compiled, then format run
        if '.' in string: # g-style trimming
            string = string.rstrip('0').rstrip('.')
        if string=='-0': # special case
            string = '0'
        string = re.sub('-', '−', string) # pure unicode minus
        # string = re.sub('-', '${-}$', string) # latex version
        # string = re.sub('-', u'\u002d', string) # unicode hyphen minus, looks same as hyphen
        # string = re.sub('-', r'\scalebox{0.75}[1.0]{$-$}', string)
        return string
    # And create object
    return mticker.FuncFormatter(f)

def LatFormatter(precision=None, sine=False, cardinal=True):
    """
    Format latitude labels; can convert sine-lats back into lats (for
    areal weighting) and can apply N/S instead of postiive/negative.
    """
    precision = 0 if precision is None else precision
    def f(value, location):
        # Convert from sine to latitude number
        if sine:
            if abs(value)>1: raise ValueError("Sine latitudes must be in range [-1,1].")
            value = np.arcsin(value)*180/np.pi
        # Suffix to apply
        if cardinal and value<0:
            value *= -1
            suffix = 'S'
        elif cardinal:
            suffix = 'N'
        else:
            suffix = ''
        # string = formatstr % (value, string)
        # Return special string, as in Formatter method
        string = f'{{{0}:.{precision:d}f}}'.format(value) # f-string compiled, then call format
        if '.' in string: # g-style trimming
            string = string.rstrip('0').rstrip('.')
        if string=='-0': # special case
            string = '0'
        return string+suffix
    # And create object
    return mticker.FuncFormatter(f)

def LonFormatter(precision=None, cardinal=True):
    """
    Format latitude labels; can convert sine-lats back into lats (for
    areal weighting) and can apply N/S instead of postiive/negative.
    """
    precision = 0 if precision is None else precision
    def f(value, location):
        # Suffix to apply
        if cardinal and value<0:
            value *= -1
            suffix = 'W'
        elif cardinal:
            suffix = 'E'
        else:
            suffix = ''
        # string = formatstr % (value, string)
        # Return special string, as in Formatter method
        string = f'{{{0}:.{precision:d}f}}'.format(value) # f-string compiled, then call format
        if '.' in string: # g-style trimming
            string = string.rstrip('0').rstrip('.')
        if string=='-0': # special case
            string = '0'
        return string+suffix
    # And create object
    return mticker.FuncFormatter(f)

def FracFormatter(fact=np.pi, symbol=r'\pi'):
    """
    Format as fractions, multiples of some value.
    """
    # Start with fraction definition
    def f(n, loc): # must accept location argument
        frac = Fraction(n/fact).limit_denominator()
        if n==0: # zero
            return '0'
        elif frac.denominator==1: # denominator is one
            if frac.numerator==1:
                return r'$%s$' % (symbol,)
            elif frac.numerator==-1:
                return r'${-}%s$' % (symbol,)
            else:
                return r'$%d%s$' % (frac.numerator, symbol)
        elif frac.numerator==1: # numerator is +/-1
            return r'$%s/%d$' % (symbol, frac.denominator)
        elif frac.numerator==-1:
            return r'${-}%s/%d$' % (symbol, frac.denominator)
        else: # otherwise
            return r'$%d%s/%d$' % (frac.numerator, symbol, frac.denominator)
    # And create FuncFormatter class
    return mticker.FuncFormatter(f)
