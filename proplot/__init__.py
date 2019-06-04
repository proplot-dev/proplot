#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# Import everything in this folder into a giant module
# Files are segretated by function, so we don't end up with
# giant 5,000-line single file
#------------------------------------------------------------------------------#
# First set up notebook
name = 'ProPlot'
# Monkey patch warnings format
# See: https://stackoverflow.com/a/2187390/4970632
# For internal warning call signature: https://docs.python.org/3/library/warnings.html#warnings.showwarning
import warnings
_warning_default = warnings.formatwarning
def _warning_proplot(message, category, filename, lineno, line=None):
    if 'proplot' in filename or _warning_default is _warning_proplot:
        return f'{filename}:{lineno}: ProPlotWarning: {message}'
    else:
        return _warning_default(message, category, filename, lineno, line=line)
warnings.formatwarning = _warning_proplot
# Then import stuff
# WARNING: Must import colortools and register names first, since rcmod will
# try to look up e.g. 'sunset' in the colormap dictionary!
from .utils import *      # misc stuff
from .colortools import * # color tools
from .rcmod import *      # custom configuration implementation
from .axes import *       # everything, axes definitions
from .gridspec import *   # gridspec objects
from .wrappers import *   # wrappers
from .subplots import *
from .fonttools import *  # fonts
from .axistools import *  # locators, normalizers, and formatters
from .projs import *      # projections and whatnot
from .demos import *      # demonstrations
