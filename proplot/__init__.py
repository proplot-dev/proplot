#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# Import everything in this folder into a giant module
# Files are segretated by function, so we don't end up with
# giant 5,000-line single file
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
# WARNING: Must import styletools and register names first, since rctools will
# try to look up e.g. 'sunset' in the colormap dictionary! Also must import
# subplots within demo functions.
from .utils import *      # misc stuff
from .styletools import * # color tools
from .rctools import *    # custom configuration implementation
from .axes import *       # everything, axes definitions
from .wrappers import *   # wrappers
from .axistools import *  # locators, normalizers, and formatters
from .subplots import *
from .gridspec import *   # gridspec objects
from .projs import *      # projections and whatnot

