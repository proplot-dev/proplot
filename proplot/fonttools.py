#!/usr/bin/env python3
"""
Script that simply lists the availble system fonts.
Add to this.
"""
import os
import re
import shutil
from glob import glob
from matplotlib import matplotlib_fname
from matplotlib import get_cachedir
import matplotlib.font_manager as mfonts
# from subprocess import Popen, PIPE
#------------------------------------------------------------------------------
# List the current font names, original version; works on Linux but not Mac,
# because can't find mac system fonts?
#------------------------------------------------------------------------------#
# List the system font names, smarter version
# See: https://github.com/olgabot/sciencemeetproductivity.tumblr.com/blob/master/posts/2012/11/how-to-set-helvetica-as-the-default-sans-serif-font-in.md
# Also see: https://olgabotvinnik.com/blog/2012-11-15-how-to-set-helvetica-as-the-default-sans-serif-font-in/
# Also see: https://stackoverflow.com/questions/18821795/how-can-i-get-list-of-font-familyor-name-of-font-in-matplotlib
_dir_data = re.sub('/matplotlibrc$', '', matplotlib_fname())
fonts_mpl_files = sorted(glob(f"{_dir_data}/fonts/ttf/*.[ot]tf"))
fonts_os_files  = sorted(mfonts.findSystemFonts(fontpaths=None, fontext='ttf')) # even with that fontext, will include otf! weird
fonts_os, fonts_mpl = set(), set()
for _file in fonts_os_files:
    try:
        fonts_os.add(mfonts.FontProperties(fname=_file).get_name())
    except Exception as err:
        pass # fails sometimes
for _file in fonts_mpl_files:
    try:
        fonts_mpl.add(mfonts.FontProperties(fname=_file).get_name())
    except Exception as err:
        pass # fails sometimes
fonts = {*fonts_os, *fonts_mpl}

# Missing fonts (add to this list whenever user requests one)
_missing_fonts = []

#------------------------------------------------------------------------------#
# Function that sets up any .ttf fonts contained in the <fonts> directory to
# be detected by matplotlib, regardless of the OS.
# See: https://olgabotvinnik.com/blog/2012-11-15-how-to-set-helvetica-as-the-default-sans-serif-font-in/
# Best fonts for dyslexia: http://dyslexiahelp.umich.edu/sites/default/files/good_fonts_for_dyslexia_study.pdf
#------------------------------------------------------------------------------#
# ***Notes on getting ttf files on Mac****
# /System/Library/Font *OR* /Library/Fonts
# * The .otf files work in addition to .ttf files; you can verify this by
#   looking at plot.fonts_files_os (results of findSystemFonts command) --
#   the list will include a bunch of .otf files.
# * Some system fonts are .dfont which are unreadable to matplotlib; use
#   fondu (download with Homebrew) to break down.
# * Sometimes fondu created .bdf files, not .ttf; use https://github.com/Tblue/mkttf
#   Requires FontForge and PoTrace
# * Install new fonts with: "brew cask install font-<name-of-font>" after using
#   "brew tap caskroom/fonts" to initialize; these will appear in ~/Library/Fonts.
# * For some conversion tools, run "pip install fonttools". Documentation
#   is here: https://github.com/fonttools/fonttools
#------------------------------------------------------------------------------#
# ***Notes on default files packaged in font directory.
# * Location will be something like: /lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf
# * 'STIX' fonts allow different LaTeX-like math modes e.g. blackboard bold
#   and caligraphy; see: https://matplotlib.org/gallery/text_labels_and_annotations/stix_fonts_demo.html
# * The 'cm'-prefix fonts seem to provide additional mathematical symbols
#   like integrals, and italized math-mode fonts.
# * We also have 'pdfcorefonts' in this directory, but I think since these
#   are afm matplotlib cannot use them? Don't know.
#------------------------------------------------------------------------------#
def install_fonts():
    """
    Install matplotlib fonts from ttf files located in the 'fonts' directory.
    May require restarting iPython session. Note font cache will be deleted
    in this process, which could cause delays.
    """
    # See: https://stackoverflow.com/a/2502883/4970632
    # Just print strings because in notebooks will get printed
    # to terminal; see: https://stackoverflow.com/q/38616077/4970632
    # _font_script = f'{os.path.dirname(__file__)}/fontscript.sh'
    # p = Popen(_font_script, stdout=PIPE, shell=False)
    # out, err = p.communicate()
    # print(out.decode('utf-8').strip()) # strip trailing newline

    # New method, just do this in python
    dir_source = f'{os.path.dirname(__file__)}/fonts' # should be in same place as scripts
    dir_dest = f'{_dir_data}/fonts/ttf'
    # print(f'Transfering .ttf and .otf files from {dir_source} to {dir_dest}.')
    for file in glob(f'{dir_source}/*.[ot]tf'):
        if not os.path.exists(f'{dir_dest}/{os.path.basename(file)}'):
            print(f'Adding font "{os.path.basename(file)}".')
            shutil.copy(file, dir_dest)

    # Delete cache
    dir_cache = get_cachedir()
    for file in glob(f'{dir_cache}/*.cache') + glob(f'{dir_cache}/font*'):
        if not os.path.isdir(file): # don't dump the tex.cache folder... because dunno why
            os.remove(file)
            print(f'Deleted font cache {file}.')

    # Rebuild
    mfonts._rebuild()
    print('Rebuilt font library.')

    # Message
    print('Fonts have been installed and font cache has been emptied. Please restart your iPython session.')
