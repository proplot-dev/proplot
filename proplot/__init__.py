#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# Import everything into the top-level module namespace
# Have sepearate files for various categories, so we don't end up with a
# single enormous 12,000-line file
#------------------------------------------------------------------------------#
# Constants
name = 'ProPlot'
__version__ = '1.0'

# Monkey patch warnings format for warnings issued by ProPlot, make sure to
# detect if this is just a matplotlib warning traced back to ProPlot code
# See: https://stackoverflow.com/a/2187390/4970632
# For internal warning call signature: https://docs.python.org/3/library/warnings.html#warnings.showwarning
# For default warning source code see: https://github.com/python/cpython/blob/master/Lib/warnings.py
import warnings
def _warning_proplot(message, category, filename, lineno, line=None):
    if line is None:
        try:
            import linecache
            line = linecache.getline(filename, lineno)
        except ModuleNotFoundError:
            pass
    if 'proplot' in filename and line is not None and 'warnings' in line:
        string = f'{filename}:{lineno}: ProPlotWarning: {message}'
    else:
        string = f'{filename}:{lineno}: {category.__name__}: {message}'
        if line is not None:
            string += ('\n' + line) # default behavior
    return (string + '\n') # must end in newline or not shown in IPython
if warnings.formatwarning is not _warning_proplot:
    warnings.formatwarning = _warning_proplot

# Import stuff
# WARNING: Import order is meaningful! Loads modules that are dependencies
# of other modules last, and loads styletools early so we can try to update
# TTFPATH before the fontManager is loaded by other matplotlib modules
from .utils import DEBUG
from .utils import *
if DEBUG:
    import time
    t = time.clock()
    t0 = t
from .styletools import * # colors and fonts
if DEBUG:
    print(f'styletools: {time.clock() - t}')
    t = time.clock()
from .rctools import * # custom configuration implementation
if DEBUG:
    print(f'rctools: {time.clock() - t}')
    t = time.clock()
from .axistools import * # locators, normalizers, and formatters
if DEBUG:
    print(f'axistools: {time.clock() - t}')
    t = time.clock()
from .wrappers import * # wrappers
if DEBUG:
    print(f'wrappers: {time.clock() - t}')
    t = time.clock()
from .projs import * # projections and whatnot
if DEBUG:
    print(f'projs: {time.clock() - t}')
    t = time.clock()
from .axes import * # axes classes
if DEBUG:
    print(f'axes: {time.clock() - t}')
    t = time.clock()
from .subplots import * # subplots and figure class
if DEBUG:
    print(f'subplots: {time.clock() - t}')
    t = time.clock()
    print(f'total time: {time.clock() - t0}')

