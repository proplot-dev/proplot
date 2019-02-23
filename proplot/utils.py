#!/usr/bin/env python3
"""
A couple simple tools. They aren't strictly plot-related functions, but
will be useful for user in the context of making plots.
"""
import time
import numpy as np
import functools
# from inspect import cleandoc
try:
    from icecream import ic
except ImportError:  # graceful fallback if IceCream isn't installed.
    ic = lambda *a: None if not a else (a[0] if len(a) == 1 else a) # noqa
_default = (lambda x,y: x if x is not None else y) # fill if not None

# Helper class
class _dot_dict(dict):
    """
    Simple class for accessing elements with dot notation.
    See: https://stackoverflow.com/a/23689767/4970632
    """
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

def units(value, error=True):
    """
    Flexible units! See `this link <http://iamvdo.me/en/blog/css-font-metrics-line-height-and-vertical-align#lets-talk-about-font-size-first>`_
    for info on the em square units.

    Parameters
    ----------
    value : float or str
        A size 'unit'. If numeric, assumed unit is inches.

        If string, we look for the format ``'123.456units'``, where the
        number is the value and `'units'` are one of the following:

        =====  ================================
        Name   Description
        =====  ================================
        em     Em-square 
        ex     Ex-square 
        lh     Line height, or 1.2 em-squares
        lem    Em-square for large-sized text
        lex    Ex-square for large-sized text
        llh    Line height, or 1.2 em-squares for large-sized text
        cm     Centimeters
        mm     Millimeters
        pt     Points, or 1/72 inches
        in     Inches
        =====  ================================

    error : bool, optional
        Raise error on failure?

    """
    if not isinstance(value, str):
        return value # assume int/float is in inches
    unit_dict = {
        'em': rc['small']/72.0,
        'ex': 0.5*rc['large']/72.0, # more or less; see URL
        'lh': 1.2*rc['small']/72.0, # line height units (default spacing is 1.2 em squares)
        'lem': rc['small']/72.0, # for large text
        'lex': 0.5*rc['large']/72.0,
        'llh': 1.2*rc['large']/72.0,
        'cm': 0.3937,
        'mm': 0.03937,
        'pt': 1/72.0,
        'in': 1.0, # already in inches
        }
    regex = re.match('^(.*)(' + '|'.join(unit_dict.keys()) + ')$', value)
    if not regex:
        if error:
            raise ValueError(f'Invalid size spec {value}.')
        else:
            return value
    num, unit = regex.groups()
    try:
        num = float(num)
    except ValueError:
        if error:
            raise ValueError(f'Invalid size spec {value}.')
        else:
            return value
    return num*unit_dict[unit] # e.g. cm / (in / cm)

#------------------------------------------------------------------------------#
# Decorators
# Below is very simple example that demonstrates how simple decorators work
# def decorator1(func):
#     def decorator():
#         print('decorator 1 called')
#         func()
#         print('decorator 1 finished')
#     return decorator
# def decorator2(func):
#     def decorator():
#         print('decorator 2 called')
#         func()
#         print('decorator 2 finished')
#     return decorator
# @decorator1
# @decorator2
# def hello():
#     print('hello world!')
# hello()
#------------------------------------------------------------------------------#
# Throw this one out, use kwarg interpolation from method.
# Why? Because matplotlib plot_directive sphinx extension will look for
# gallery images in the old documentation that do not exist for ProPlot.
# def _docstring_fix(child):
#     """
#     Decorator function for appending documentation from overridden method
#     onto the overriding method docstring.
#     Adapted from: https://stackoverflow.com/a/8101598/4970632
#     """
#     for name,chfunc in vars(child).items(): # returns __dict__ object
#         if not callable(chfunc): # better! see: https://stackoverflow.com/a/624939/4970632
#             continue
#         for parent in getattr(child, '__bases__', ()):
#             # Obtain documentation
#             parfunc = getattr(parent, name, None)
#             if not getattr(parfunc, '__doc__', None):
#                 continue
#             if not getattr(chfunc, '__doc__', None):
#                 chfunc.__doc__ = '' # in case it's None
#
#             # Ugly
#             # cmessage = f'Full name: {parfunc.__qualname__}()'
#             # pmessage = f'Parent method (documentation below): {chfunc.__qualname__}()'
#             # chfunc.__doc__ = f'\n{cmessage}\n{cleandoc(chfunc.__doc__)}\n{pmessage}\n{cleandoc(parfunc.__doc__)}'
#
#             # Simple
#             # chfunc.__doc__ = f'{cleandoc(chfunc.__doc__)}\n\n\n{cleandoc(parfunc.__doc__)}'
#
#             # Fails because numpydoc disallows custom subsections; see: https://developer.lsst.io/python/numpydoc.html#sections-are-restricted-to-the-numpydoc-section-set
#             # chfunc.__doc__ = f'ProPlot Override\n================\n{cleandoc(chfunc.__doc__)}\n\n\n' \
#             #                  f'Original Documentation\n======================\n{cleandoc(parfunc.__doc__)}'
#
#             # Differentiate with bold text
#             chfunc.__doc__ = f'**ProPlot Override**\n{cleandoc(chfunc.__doc__)}\n\n\n' \
#                              f'**Original Documentation**\n{cleandoc(parfunc.__doc__)}'
#             break # only do this for the first parent class
#     return child

def _timer(func):
    """
    A decorator that prints the time a function takes to execute.
    See: https://stackoverflow.com/a/1594484/4970632
    """
    @functools.wraps(func)
    def decorator(*args, **kwargs):
        t = time.clock()
        print(f'{func.__name__}()')
        res = func(*args, **kwargs)
        print(f'{func.__name__}() time: {time.clock()-t}s')
        return res
    return decorator

def _logger(func):
    """
    A decorator that logs the activity of the script (it actually just prints it,
    but it could be logging!)
    See: https://stackoverflow.com/a/1594484/4970632
    """
    @functools.wraps(func)
    def decorator(*args, **kwargs):
        res = func(*args, **kwargs)
        print(f'{func.__name__} called with: {args} {kwargs}')
        return res
    return decorator

def _counter(func):
    """
    A decorator that counts and prints the number of times a function
    has been executed.
    See: https://stackoverflow.com/a/1594484/4970632
    """
    @functools.wraps(func)
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
# Accessible for user
#------------------------------------------------------------------------------#
def arange(min_, *args):
    """
    Duplicate behavior of `numpy.arange`, except with **inclusive endpoints**.
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
    idxs = [[0] for _ in range(values.ndim-1)] # must be list because we use it twice
    if values[np.ix_(*idxs, [1])] < values[np.ix_(*idxs, [0])]:
        flip = True
        values = np.flip(values, axis=-1)
    values = np.concatenate((
        values[...,:1]  - (values[...,1]-values[...,0])/2,
        (values[...,1:] + values[...,:-1])/2,
        values[...,-1:] + (values[...,-1]-values[...,-2])/2,
        ), axis=-1)
    if flip:
        values = np.flip(values, axis=-1)
    # Permute back and return
    values = np.swapaxes(values, axis, -1)
    return values

