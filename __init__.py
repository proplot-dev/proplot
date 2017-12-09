# Note that only difference between dir() and inspect.getmembers()
# is that getmembers returns list of tuple pairs of objects and their names,
# whereas dir() returns only a list of names

#------------------------------------------------------------------------------;
# Most sane/common method; stop worrying about extra objects
# The import * will ignore any functions starting with _; use this property
#------------------------------------------------------------------------------
# Imports
# from . import CONSTANTS as c # physical constants
# from .numerics import *
# from .ioutils import *
# from .sciplotlib import *

#-------------------------------------------------------------------------------
# IPython session management
# TODO Now moved to notebook executable, maybe delete this
#-------------------------------------------------------------------------------
# def notebook(stem=None):
#     """
#     Setup up magic commands and other notebook properties.
#     Won't mess things up if we're not in a notebook.
#     Input...
#         stem: the directory stem from HOME directory.
#     """
#     import os
#     from IPython import get_ipython
#     if stem is not None:
#         stem = stem.replace('~',os.environ['HOME'])
#         os.chdir(stem)
#     if get_ipython() is not None:
#         # Basic configuration
#         # ipy.core.display.display(ipy.core.display.HTML(html_format))
#         get_ipython().magic('reload_ext autoreload') # reload instead of load, to avoid annoying message
#         get_ipython().magic('autoreload 2') # turn on expensive autoreloading
#         if getattr(get_ipython(), 'kernel', None) is None:
#             # Nothing else
#             print("Configured ipython session.")
#         else:
#             # Notebook configuration
#             # get_ipython().magic('config InlineBackend.figure_format=\'svg\'')
#             get_ipython().magic('autosave 120') # autosave every 120 seconds
#             get_ipython().magic('config InlineBackend.figure_format=\'retina\'')
#                 # retina probably more space efficient (high-res bitmap), but svg is prettiest
#                 # and is only one preserving vector graphics
#             get_ipython().magic('config InlineBackend.print_figure_kwargs=dict(bbox_inches=None)') #bbox_inches=\'tight\', pad_inches=0.1)')
#             get_ipython().magic('matplotlib inline') # change print_figure_kwargs to  see edges
#             print("Configured notebook.")
#
#             # CSS modifications; use here, or place in ~/.jupyter/custom/custom.css
#             # import IPython.core.display as display # use these to edit CSS
#             # ipy.core.display.display(ipy.core.display.HTML("<style>.container { width:100% !important; }</style>"))
#             # display.display(display.HTML("<style>"
#             #     + "#notebook { padding-top:0px !important; }" # hashtags refer to ids; dots refer to classes
#             #     + ".container { width:100% !important; } "
#             #     + ".end_space { min-height:0px !important; }"
#             #     + "</style>"))
#             # # icolor, ncolor = "#2b2b25"
#             # # icolor, ncolor = "#f2f2f2", "#f2f2f2"
#             # icolor, ncolor = "#f6f6f6", "#f6f6f6"
#             # display.display(display.HTML("<style>"
#             #     + ".edit_mode .cell.selected .CodeMirror-focused.cm-fat-cursor { background-color: " + ncolor + " !important; } " # normal mode
#             #     + ".edit_mode .cell.selected .CodeMirror-focused:not(.cm-fat-cursor) { background-color: " + icolor + " !important; }" # edit mode
#             #         # the above matches color of ITerm2 settings
#             #     + "</style>")) # <style> denotes CSS code block

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
# ...

#------------------------------------------------------------------------------
# Import contents of __all__ from *arbitrary* modules in directory
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
