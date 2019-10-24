#!/usr/bin/env python3
"""
Simple tools used in various places across this package.
"""
import re
import os
import glob
import time
import numpy as np
import functools
import warnings
import matplotlib as mpl
from numbers import Number, Integral
rcParams = mpl.rcParams
try:
    from icecream import ic
except ImportError:  # graceful fallback if IceCream isn't installed
    ic = lambda *a: None if not a else (a[0] if len(a) == 1 else a) # noqa
__all__ = ['arange', 'edges', 'units', 'DEBUG']

# Debug decorators
DEBUG = False # debug mode, used for profiling and activating timer decorators
def _logger(func):
    """A decorator that logs the activity of the script (it actually just prints it,
    but it could be logging!). See: https://stackoverflow.com/a/1594484/4970632"""
    @functools.wraps(func)
    def decorator(*args, **kwargs):
        res = func(*args, **kwargs)
        if DEBUG:
            print(f'{func.__name__} called with: {args} {kwargs}')
        return res
    return decorator

def _timer(func):
    """A decorator that prints the time a function takes to execute.
    See: https://stackoverflow.com/a/1594484/4970632"""
    @functools.wraps(func)
    def decorator(*args, **kwargs):
        if DEBUG:
            t = time.clock()
        res = func(*args, **kwargs)
        if DEBUG:
            print(f'{func.__name__}() time: {time.clock()-t}s')
        return res
    return decorator

def _counter(func):
    """A decorator that counts and prints the cumulative time a function
    has benn running. See: https://stackoverflow.com/a/1594484/4970632"""
    @functools.wraps(func)
    def decorator(*args, **kwargs):
        if DEBUG:
            t = time.clock()
        res = func(*args, **kwargs)
        if DEBUG:
            decorator.time += (time.clock() - t)
            decorator.count += 1
            print(f'{func.__name__}() cumulative time: {decorator.time}s ({decorator.count} calls)')
        return res
    decorator.time = 0
    decorator.count = 0 # initialize
    return decorator

# Important private helper func
def _notNone(*args, names=None):
    """Returns the first non-``None`` value, used with keyword arg aliases and
    for setting default values. Ugly name but clear purpose. Pass the `names`
    keyword arg to issue warning if multiple args were passed. Must be list
    of non-empty strings."""
    if names is None:
        for arg in args:
            if arg is not None:
                return arg
        return arg # last one
    else:
        first = None
        kwargs = {}
        if len(names) != len(args) - 1:
            raise ValueError(f'Need {len(args)+1} names for {len(args)} args, but got {len(names)} names.')
        names = [*names, '']
        for name,arg in zip(names,args):
            if arg is not None:
                if first is None:
                    first = arg
                if name:
                    kwargs[name] = arg
        if len(kwargs)>1:
            warnings.warn(f'Got conflicting or duplicate keyword args, using the first one: {kwargs}')
        return first

# Accessible for user
def arange(min_, *args):
    """Identical to `numpy.arange`, but with inclusive endpoints. For
    example, ``plot.arange(2,4)`` returns ``np.array([2,3,4])`` instead
    of ``np.array([2,3])``."""
    # Optional arguments just like np.arange
    if len(args) == 0:
        max_ = min_
        min_ = 0
        step = 1
    elif len(args) == 1:
        max_ = args[0]
        step = 1
    elif len(args) == 2:
        max_ = args[0]
        step = args[1]
    else:
        raise ValueError('Function takes from one to three arguments.')
    # All input is integer
    if all(isinstance(val,Integral) for val in (min_,max_,step)):
        min_, max_, step = np.int64(min_), np.int64(max_), np.int64(step)
        max_ += 1
    # Input is float or mixed, cast to float64
    # Don't use np.nextafter with np.finfo(np.dtype(np.float64)).max, because
    # round-off errors from continually adding step to min mess this up
    else:
        min_, max_, step = np.float64(min_), np.float64(max_), np.float64(step)
        max_ += step/2
    return np.arange(min_, max_, step)

def edges(array, axis=-1):
    """
    Calculates approximate "edge" values given "center" values. This is used
    internally to calculate graitule edges when you supply centers to
    `~matplotlib.axes.Axes.pcolor` or `~matplotlib.axes.Axes.pcolormesh`, and
    in a few other places.

    Parameters
    ----------
    array : array-like
        Array of any shape or size. Generally, should be monotonically
        increasing or decreasing along `axis`.
    axis : int, optional
        The axis along which "edges" are calculated. The size of this axis
        will be augmented by one.

    Returns
    -------
    `~numpy.ndarray`
        Array of "edge" coordinates.
    """
    # First permute
    array = np.array(array)
    array = np.swapaxes(array, axis, -1)
    # Next operate
    flip = False
    idxs = [[0] for _ in range(array.ndim-1)] # must be list because we use it twice
    if array[np.ix_(*idxs, [1])] < array[np.ix_(*idxs, [0])]:
        flip = True
        array = np.flip(array, axis=-1)
    array = np.concatenate((
        array[...,:1]  - (array[...,1]-array[...,0])/2,
        (array[...,1:] + array[...,:-1])/2,
        array[...,-1:] + (array[...,-1]-array[...,-2])/2,
        ), axis=-1)
    if flip:
        array = np.flip(array, axis=-1)
    # Permute back and return
    array = np.swapaxes(array, axis, -1)
    return array

def _check_path(dir, *args):
    """Returns configuration directory and checks for unexpected extra files
    inside the directory. Issues helpful warning for new users."""
    dir = os.path.join(dir, '.proplot')
    paths = {os.path.basename(path) for path in glob.glob(os.path.join(dir, '*'))}
    paths_allowed = {'proplotrc', 'cmaps', 'colors', 'cycles', 'fonts'}
    if not paths <= paths_allowed:
        warnings.warn(f'Found extra files {", ".join(paths - paths_allowed)} in the ~/.proplot folder. Files must be placed in the .proplot/cmaps, .proplot/colors, .proplot/cycles, or .proplot/fonts subfolders.')
    return os.path.join(dir, *args)

def get_configpaths(sub=None, *, local=True):
    """
    Returns ProPlot configuration file and folder paths, or the paths of
    root configuration directories.

    Parameters
    ----------
    sub : {None, 'proplotrc', 'cmaps', 'cycles', 'fonts'}, optional
        The subfolder or file name. If empty, this function just returns
        the root configuration directories.
    local : bool, optional
        Whether to look for local configuration folders in current directory
        and parent directories.

    Returns
    -------
    paths : list
        List of paths.
    """
    # Checks
    if not sub:
        args = ()
    elif sub in ('proplotrc', 'cmaps', 'colors', 'cycles', 'fonts'):
        args = (sub,)
    else:
        raise ValueError(f'Invalid configuration location {sub!r}. Options are "proplotrc", "cmaps", "cycles", and "fonts".')
    if sub == 'proplotrc' and os.path.exists(os.path.join(os.path.expanduser('~'), '.proplotrc')):
        warnings.warn(f'Configuration file location is now "~/.proplot/proplotrc" to match matplotlib convention. Please move "~/.proplotrc" to "~/.proplot/proplotrc".')
    # Get paths
    paths = []
    if local:
        idir = os.getcwd()
        while idir: # not empty string
            ipath = _check_path(idir, *args)
            if os.path.exists(ipath):
                paths.append(ipath)
            ndir, _ = os.path.split(idir)
            if ndir == idir:
                break
            idir = ndir
        paths = paths[::-1] # sort from decreasing to increasing importantce
    # Home configuration
    ipath = _check_path(os.path.expanduser('~'), *args)
    if os.path.exists(ipath) and ipath not in paths:
        paths.insert(0, ipath)
    # Global configuration
    ipath = _check_path(os.path.dirname(__file__), *args)
    if not os.path.exists(ipath):
        raise ValueError(f'Default configuration location {ipath!r} does not exist.')
    elif ipath not in paths:
        paths.insert(0, ipath)
    return paths

def units(value, numeric='in'):
    """
    Flexible units -- this function is used internally all over ProPlot, so
    that you don't have to use "inches" or "points" for all sizing arguments.
    See `this link <http://iamvdo.me/en/blog/css-font-metrics-line-height-and-vertical-align#lets-talk-about-font-size-first>`_
    for info on the em square units.

    Parameters
    ----------
    value : float or str or list thereof
        A size "unit" or *list thereof*. If numeric, assumed unit is `numeric`.
        If string, we look for the format ``'123.456unit'``, where the
        number is the value and ``'unit'`` is one of the following.

        ======  ===================================================================
        Key     Description
        ======  ===================================================================
        ``m``   Meters
        ``cm``  Centimeters
        ``mm``  Millimeters
        ``ft``  Feet
        ``in``  Inches
        ``pt``  Points (1/72 inches)
        ``px``  Pixels on screen, uses dpi of ``rc['figure.dpi']``
        ``pp``  Pixels once printed, uses dpi of ``rc['savefig.dpi']``
        ``em``  Em-square for ``rc['font.size']``
        ``ex``  Ex-square for ``rc['font.size']``
        ``Em``  Em-square for ``rc['axes.titlesize']``
        ``Ex``  Ex-square for ``rc['axes.titlesize']``
        ======  ===================================================================

    numeric : str, optional
        The assumed unit for numeric arguments, and the output unit. Default
        is **inches**, i.e. ``'in'``.
    """
    # Loop through arbitrary list, or return None if input was None (this
    # is the exception).
    if not np.iterable(value) or isinstance(value, str):
        singleton = True
        values = (value,)
    else:
        singleton = False
        values = value

    # Font unit scales
    # NOTE: Delay font_manager import, because want to avoid rebuilding font
    # cache, which means import must come after TTFPATH added to environ,
    # i.e. inside styletools.register_fonts()!
    small = rcParams['font.size'] # must be absolute
    large = rcParams['axes.titlesize']
    if isinstance(large, str):
        import matplotlib.font_manager as mfonts
        scale = mfonts.font_scalings.get(large, 1) # error will be raised somewhere else if string name is invalid!
        large = small*scale

    # Dict of possible units
    unit_dict = {
        # Physical units
        'in': 1.0, # already in inches
        'm':  39.37,
        'ft': 12.0,
        'cm': 0.3937,
        'mm': 0.03937,
        'pt': 1/72.0,
        # Font units
        'em': small/72.0,
        'ex': 0.5*small/72.0, # more or less; see URL
        'Em': large/72.0, # for large text
        'Ex': 0.5*large/72.0,
        }
    # Display units
    # WARNING: In ipython shell these take the value 'figure'
    if not isinstance(rcParams['figure.dpi'], str):
        unit_dict['px'] = 1/rcParams['figure.dpi'] # on screen
    if not isinstance(rcParams['savefig.dpi'], str):
        unit_dict['pp'] = 1/rcParams['savefig.dpi'] # once 'printed', i.e. saved

    # Iterate
    try:
        scale = unit_dict[numeric]
    except KeyError:
        raise ValueError(f'Invalid numeric unit {numeric}. Valid units are {", ".join(unit_dict.keys())}.')
    result = []
    for value in values:
        if value is None or isinstance(value, Number):
            result.append(value)
            continue
        elif not isinstance(value, str):
            raise ValueError(f'Size spec must be string or number or list thereof, received {values}.')
        regex = re.match('^([-+]?[0-9.]*)(.*)$', value)
        num, unit = regex.groups()
        try:
            result.append(float(num)*unit_dict[unit]/scale)
        except (KeyError, ValueError):
            raise ValueError(f'Invalid size spec {value}. Valid units are {", ".join(unit_dict.keys())}.')
    if singleton:
        result = result[0]
    return result

