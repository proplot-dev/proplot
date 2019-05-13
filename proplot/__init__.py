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
def _warning_no_line(message, category, filename, lineno, file=None, line=None):
    return f'{filename}:{lineno}: ProPlotWarning: {message}'
warnings.formatwarning = _warning_no_line
# Then import stuff
from .rcmod import *      # custom configuration implementation
from .utils import *      # misc stuff
from .axes import *       # everything, axes definitions
from .subplots import *
from .gridspec import *
from .colortools import * # color tools
from .fonttools import *  # fonts
from .axistools import *  # locators, normalizers, and formatters
from .projs import *      # projections and whatnot
from .demos import *      # demonstrations
from .notebook import *
