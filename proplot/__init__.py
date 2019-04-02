#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# Import everything in this folder into a giant module
# Files are segretated by function, so we don't end up with
# giant 5,000-line single file
#------------------------------------------------------------------------------#
# First set up notebook
name = 'ProPlot'
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
