#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# Import everything into the top-level module namespace
# Have sepearate files for various categories, so we don't end up with a
# single enormous 12,000-line file
#------------------------------------------------------------------------------#
# Constants
name = 'ProPlot'
__version__ = '1.0'

# Monkey patch warnings format
# See: https://stackoverflow.com/a/2187390/4970632
# For internal warning call signature: https://docs.python.org/3/library/warnings.html#warnings.showwarning
# For default warning source code see: https://github.com/python/cpython/blob/master/Lib/warnings.py
# WARNING: Message must end with newline or will not be shown in ipython sessions
import warnings
def _warning_proplot(message, category, filename, lineno, line=None):
    if 'proplot' in filename:
        string = f'{filename}:{lineno}: ProPlotWarning: {message}'
    else:
        string = f'{filename}:{lineno}: {category.__name__}: {message}'
        if line is None:
            try:
                import linecache
                line = linecache.getline(filename, lineno)
                string = f'{string}\n{line}'
            except Exception:
                pass
    return string + '\n'
if warnings.formatwarning is not _warning_proplot:
    warnings.formatwarning = _warning_proplot

# Import stuff
# WARNING: Beware inter-dependencies! For example import styletools must come
# first because rctools may try to register a colormap that doesn't exist.
from .utils import * # misc stuff, debug mode
if _debug:
    import time
    t = time.clock()
    t0 = t
from .styletools import * # colors and fonts
if _debug:
    print(f'styletools: {time.clock() - t}')
    t = time.clock()
from .rctools import * # custom configuration implementation
if _debug:
    print(f'rctools: {time.clock() - t}')
    t = time.clock()
from .axistools import * # locators, normalizers, and formatters
if _debug:
    print(f'axistools: {time.clock() - t}')
    t = time.clock()
from .wrappers import * # wrappers
if _debug:
    print(f'wrappers: {time.clock() - t}')
    t = time.clock()
from .projs import * # projections and whatnot
if _debug:
    print(f'projs: {time.clock() - t}')
    t = time.clock()
from .axes import * # axes classes
if _debug:
    print(f'axes: {time.clock() - t}')
    t = time.clock()
from .subplots import * # subplots and figure class
if _debug:
    print(f'subplots: {time.clock() - t}')
    t = time.clock()
    print(f'total time: {time.clock() - t0}')

