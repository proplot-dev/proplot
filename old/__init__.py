#------------------------------------------------------------------------------;
# Leave empty to keep basic module structure; usage should be as follows:
# from pyfuncs import tools
# from pyfuncs import plot
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Import everything into a single giant module; usage should be as follows:
# import pyfuncs as py
#------------------------------------------------------------------------------
# Imports
# from . import CONSTANTS as c # physical constants
# from .numerics import *
# from .ioutils import *
# from .sciplotlib import *

#------------------------------------------------------------------------------
# Import contents of __all__ from modules by name
#------------------------------------------------------------------------------
# Instead, import by name explicitly, and control __all__ in each file so that
# we don't re-import things like numpy, matplotlib as components of that module
# Also should NOT, say, import ncutils as nc (this *will* import extra stuff)
# __all__ = ['utils', 'numerics', 'constants', 'ncutils', 'geoutils', 'plotutils']
# from . import utils
# from . import plotutils
# from . import *
# __all__ = []
# from .utils import *
# del utils
# from .ncutils import *
# from .plotutils import *
# from .geoutils import *
# from .numerics import *
# from .constants import * # for this, will don't need __all__

# if 'utils' in globals():
#     del utils
# if 'ncutils' in globals():
#     del ncutils
# if 'plotutils' in globals():
#     del plotutils
# if 'geoutils' in globals():
#     del geoutils
# if 'numerics' in globals():
#     del numerics
# if 'constants' in globals():
#     del constants
# del utils, ncutils

#------------------------------------------------------------------------------
# Import contents of __all__ from *arbitrary* modules in directory
# Note that only difference between dir() and inspect.getmembers()
# is that getmembers returns list of tuple pairs of objects and their names,
# whereas dir() returns only a list of names
#------------------------------------------------------------------------------
# import pkgutil
# __all__ = []
# for loader, name, ispkg in pkgutil.walk_packages(__path__):
#     module = loader.find_module(name).load_module(name)
#     for member_name in module.__all__:
#         globals()[member_name] = getattr(module,member_name)
#         __all__.append(member_name)
# # Cleanup; prevent all this extra crap from being in dir(myfuncs)
# del loader, name, ispkg, pkgutil, module, member_name

# ## Option 3
# # From stack overflow; imports all CONTENTS of every module in directory
# # Adds them to globals; also import * will only import that stuff
# __all__ = []
# import pkgutil
# import inspect
# for loader, name, ispkg in pkgutil.walk_packages(__path__):
#     module = loader.find_module(name).load_module(name)
#     for member_name, member in inspect.getmembers(module):
#         if member_name.startswith('__') or not (
#                 inspect.isclass(member) or inspect.isfunction(member)
#                 ):
#             continue
#         globals()[member_name] = member
#         __all__.append(member_name)
# # Cleanup; prevent all this extra crap from being in dir(myfuncs)
# del pkgutil, inspect, loader, name, ispkg, member_name, member, module
