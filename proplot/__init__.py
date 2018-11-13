#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# Import everything in this folder into a giant module
# Files are segretated by function, so we don't end up with
# giant 5,000-line single file
#------------------------------------------------------------------------------#
# First set up notebook
from .notebook import *
name = 'ProPlot'
# Then import stuff
from .base import *   # basic tools
from .subplots import *
from .gridspec import *
from .colortools import * # color tools
from .fonttools import * # fonts
from .proj import *   # cartopy projections and whatnot
from .axis import *   # locators, normalizers, and formatters
from .utils import *  # misc stuff
from .rcmod import *     # custom configuration implementation
