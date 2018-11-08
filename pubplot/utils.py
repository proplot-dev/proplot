#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# Just a few tools
# These aren't strictly plot-related functions, but will be useful for user
# in the context of making plots.
#------------------------------------------------------------------------------#
import numpy as np
from numbers import Number

#------------------------------------------------------------------------------#
# Helper class
#------------------------------------------------------------------------------#
class dot_dict(dict):
    """
    Simple class for accessing elements with dot notation.
    See: https://stackoverflow.com/a/23689767/4970632
    """
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

#------------------------------------------------------------------------------#
# Definitions
#------------------------------------------------------------------------------#
def isvector(item):
    """
    Just test if is iterable, but not a string (we almost never mean this).
    """
    # return hasattr(item, '__iter__') and not isinstance(item, str)
    return np.iterable(item) and not isinstance(item, str)

def isnumber(item):
    """
    Just test if number.
    """
    return isinstance(item, Number)

def isscalar(item):
    """
    Test if item is a number or a string, as opposed to list/tuple/array.
    See: https://stackoverflow.com/questions/4187185/how-can-i-check-if-my-python-object-is-a-number
    Note: Numpy numbers have __getitem__ attribute! So cannot test this.
    Why is this done? So they can be converted to ND singleton numpy arrays
    easily with number[None,None,...].
    """
    # return (hasattr(item,'__iter__') or hasattr(item,'__getitem__')) # integers have this!
    # return hasattr(item,'__iter__')
    return isinstance(item, Number) or isinstance(item, str)

def edges(values, axis=-1):
    """
    Get approximate edge values along arbitrary axis.
    """
    # First permute
    values = np.array(values)
    values = np.swapaxes(values, axis, -1)
    # Next operate
    flip = False
    if values[...,1]<values[...,0]:
        flip = True
        values = np.flip(values, axis=-1)
    values = np.concatenate((
        values[...,:1] - (values[...,1]-values[...,0])/2,
        (values[...,1:] + values[...,:-1])/2,
        values[...,-1:] + (values[...,-1]-values[...,-2])/2,
        ), axis=-1)
    if flip:
        values = np.flip(values, axis=-1)
    # Permute back and return
    values = np.swapaxes(values, axis, -1)
    return values

def arange(min_, *args):
    """
    Duplicate behavior of np.arange, except with inclusive endpoints; dtype is
    controlled very carefully, so should be 'most precise' among min/max/step args.
    Input:
        stop
        start, stop, [step]
        Just like np.arange!
    Output:
        The array sequence.
    """
    # Optional arguments just like np.arange
    if len(args)==0:
        max_ = min_
        min_ = 0 # this re-assignes the NAME "min_" to 0
        step = 1
    elif len(args)==1:
        max_ = args[0]
        step = 1
    elif len(args)==2:
        max_ = args[0]
        step = args[1]
    else:
        raise ValueError('Function takes from one to three arguments.')
    # All input is integer? Get new "max"
    if min_//1==min_ and max_//1==max_ and step//1==step:
        min_, max_, step = np.int64(min_), np.int64(max_), np.int64(step)
        max_ += 1
    # Input is float or mixed; cast all to float64, then get new "max"
    else:
        # Get the next FLOATING POINT, in direction of the second argument
        # Forget this; round-off errors from continually adding step to min mess this up
        # max_ = np.nextafter(max_, np.finfo(np.dtype(np.float64)).max)
        min_, max_, step = np.float64(min_), np.float64(max_), np.float64(step)
        max_ += step/2
    return np.arange(min_, max_, step)

