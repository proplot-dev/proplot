#!/usr/bin/env python3
"""
Script to set up magic commands and other notebook properties.
Just call this in the first cell with "run pyfuncs/notebook". And voila!
Call with optional directory argument to make that the new working directory.

Unbelievably weird problem:
  * Warning on calling rcdefaults(): https://stackoverflow.com/q/48320804/4970632
    Apparently you can change the backend until the first plot is drawn, then
    it stays the same. So the rcdefault() command changes the backend to a
    non-inline version.

Notes on python figures:
  * Can use InlineBackend rc configuration to make inline figure properties
    different from figure.<subproperty> settings in rcParams.
  * Problem is, whenever rcParams are reset/pyfuncs module is reloaded, the previous
    InlineBackend properties disappear.
  * It is *also* necessary to maintain separate savefig options, including 'transparent'
    and 'facecolor' -- cannot just set these to use the figure properties.
    If transparent set to False, saved figure will have no transparency *even if* the
    default figure.facecolor has zero alpha. Will *only* be transparent if alpha explicitly
    changed by user command. Try playing with settings in plot.globals to see.
  * In conclusion: Best workflow is probably to set InlineBackend settings to empty, but
    control figure settings separate from savefig settings. Also need to test: if
    transparent=True but patches were set to have zero transparency manually, will they
    be made re-transparent when figure is saved?

Notes on jupyter configuration:
  * In .jupyter, the jupyter_nbconvert_config.json sets up locations of stuff; templates
      for nbextensions and formatting files for markdown/code cells.
  * In .jupyter, the jupyter_notebook_config.json installs the configurator extension
      for managing extra plugins.
  * In .jupyter, not sure yet how to successfully use jupyter_console_config.py and
      jupyter_notebook_config.py; couldn't get it to do what this function does on startup.
  * In .jupyter/custom, current_theme.txt lists the current jupyterthemes theme, custom.css
      contains CSS formatting for it, and fonts should contain font files -- note that there
      are not font files on my Mac, even though jupyterthemes works; sometimes may be empty
  * In .jupyter/nbconfig, tree.json loads the extra tab for the NBconfigurator, and 
      common.json gives option to hide incompatible plugs, and notebook.json contains all
      the new settings; just copy it over to current notebook to update
"""
# Imports
from IPython import get_ipython
from IPython.utils import io
import os
import sys
import socket
from .rcmod import rc

# @block_printing
def nbsetup(directory=None, backend='inline'):
    # Variables
    cd = os.getcwd()
    home = os.path.expanduser('~')
    hostname = socket.gethostname().split('.')[0]
    autosave = 30

    # Make sure we are in session
    ipython = get_ipython() # save session
    if ipython is None:
        print("Warning: IPython kernel not found.")
        return
        # exit()

    # Optional argument
    if directory:
        os.chdir(os.path.expanduser(directory)) # move to this directory
        print(f'Moved to directory {os.path.expanduser(directory)}.')

    # Reload modules, so can easily edit/run imported functions
    # with redirect_stdout(_null):
    if 'autoreload' not in ipython.magics_manager.magics['line']:
        # Only do this if not already loaded -- otherwise will get *recursive* 
        # reloading, even with unload_ext command!
        # with redirect_stdout(_null):
        ipython.magic("reload_ext autoreload") # reload instead of load, to avoid annoying message
        ipython.magic("autoreload 2") # turn on expensive autoreloading

    # Autosaving
    # with redirect_stdout(_null):
    # with redirect_stdout(_null):
    with io.capture_output() as captured:
        ipython.magic(f"autosave {autosave:d}") # autosave every minute
    # sys.stdout.write("\033[F") #back to previous line
    # sys.stdout.write("\033[K") #clear line

    # Initialize with default 'inline' settings
    ipython.magic("matplotlib " + backend) # change print_figure_kwargs to see edges

    # Below so don't have memory issues/have to keep re-closing them
    ipython.magic("config InlineBackend.close_figures = True")

    # Retina probably more space efficient (high-res bitmap), but svg is prettiest
    # and is only one preserving vector graphics
    ipython.magic("config InlineBackend.figure_formats = ['retina','svg']")

    # Control all rc settings directly with 'pubplot.rc' object, *no* notebook-
    # specific overrides.
    ipython.magic("config InlineBackend.rc = {}")

    # For some reason this is necessary, even with rc['savefig.bbox'] = 'standard'
    ipython.magic("config InlineBackend.print_figure_kwargs = {'bbox_inches':None}") #bbox_inches=\'tight\', pad_inches=0.1)')

    # Re-assert defaults (some get overwritten on inline initiatoin)
    rc.reset()

    # Message
    print("Configured ipython notebook.")

