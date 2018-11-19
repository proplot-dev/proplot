#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# Just a few tools
# These aren't strictly plot-related functions, but will be useful for user
# in the context of making plots.
#------------------------------------------------------------------------------#
import time
import numpy as np
from numbers import Number
from functools import wraps
from inspect import cleandoc
try:
    from icecream import ic
except ImportError:  # graceful fallback if IceCream isn't installed.
    ic = lambda *a: None if not a else (a[0] if len(a) == 1 else a) # noqa

#------------------------------------------------------------------------------#
# Decorators
#------------------------------------------------------------------------------#
def docstring_fix(child):
    """
    Decorator function for appending documentation from overridden method
    onto the overriding method docstring.
    Adapted from: https://stackoverflow.com/a/8101598/4970632
    """
    for name,chfunc in vars(child).items(): # returns __dict__ object
        if not callable(chfunc): # better! see: https://stackoverflow.com/a/624939/4970632
        # if not isinstance(chfunc, FunctionType):
            continue
        for parent in getattr(child, '__bases__', ()):
            parfunc = getattr(parent, name, None)
            if not getattr(parfunc, '__doc__', None):
                continue
            if not getattr(chfunc, '__doc__', None):
                chfunc.__doc__ = '' # in case it's None
            cmessage = f'Full name: {parfunc.__qualname__}()'
            pmessage = f'Parent method (documentation below): {chfunc.__qualname__}()'
            chfunc.__doc__ = f'\n{cmessage}\n{cleandoc(chfunc.__doc__)}\n{pmessage}\n{cleandoc(parfunc.__doc__)}'
            break # only do this for the first parent class
    return child

def fancy_decorator(decorator):
    """
    Normally to make a decorator that accepts arguments, you have to create
    3 nested function definitions. This abstracts that away -- if you decorate
    your decorator-function declaration with this, the decorator will now accept arguments.
    See: https://stackoverflow.com/a/1594484/4970632
    """
    @wraps(decorator)
    def decorator_maker(*args, **kwargs):
        def decorator_wrapper(func):
            return decorator(func, *args, **kwargs)
        return decorator_wrapper
    return decorator_maker

def timer(func):
    """
    A decorator that prints the time a function takes to execute.
    See: https://stackoverflow.com/a/1594484/4970632
    """
    @wraps(func)
    def decorator(*args, **kwargs):
        t = time.clock()
        print(f'{func.__name__}()')
        res = func(*args, **kwargs)
        print(f'{func.__name__}() time: {time.clock()-t}s')
        return res
    return decorator

def logger(func):
    """
    A decorator that logs the activity of the script (it actually just prints it,
    but it could be logging!)
    See: https://stackoverflow.com/a/1594484/4970632
    """
    @wraps(func)
    def decorator(*args, **kwargs):
        res = func(*args, **kwargs)
        print(f'{func.__name__} called with: {args} {kwargs}')
        return res
    return decorator

def counter(func):
    """
    A decorator that counts and prints the number of times a function
    has been executed.
    See: https://stackoverflow.com/a/1594484/4970632
    """
    @wraps(func)
    def decorator(*args, **kwargs):
        # decorator.count += 1
        t = time.clock()
        res = func(*args, **kwargs)
        decorator.time += (time.clock() - t)
        decorator.count += 1
        print(f'{func.__name__} cumulative time: {decorator.time}s ({decorator.count} calls)')
        # print(f'{func.__name__} has been used: {decorator.count}x')
        return res
    decorator.time = 0
    decorator.count = 0 # initialize
    return decorator

#------------------------------------------------------------------------------#
# Helper stuff
#------------------------------------------------------------------------------#
_fill = (lambda x,y: x if x is not None else y)

class _dot_dict(dict):
    """
    Simple class for accessing elements with dot notation.
    See: https://stackoverflow.com/a/23689767/4970632
    """
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

def isnumber(item):
    """
    Just test if number.
    See: https://stackoverflow.com/questions/4187185/how-can-i-check-if-my-python-object-is-a-number
    Note: Numpy numbers have __getitem__ attribute! So cannot test this.
    Why is this done? So they can be converted to ND singleton numpy arrays
    easily with number[None,None,...].
    """
    return isinstance(item, Number)

def isvector(item):
    """
    Just test if is iterable, but not a string (we almost never mean this).
    """
    # return hasattr(item, '__iter__') and not isinstance(item, str)
    return np.iterable(item) and not isinstance(item, str)

#------------------------------------------------------------------------------#
# Accessible for user
#------------------------------------------------------------------------------#
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

