#!/usr/bin/env python3
"""
Registers new fonts with `register_fonts`. Provides handy lists of available
font names. Makes Helvetica or Helvetica Neue the default font.
"""
# See: https://gree2.github.io/python/2015/04/27/python-change-matplotlib-font-on-mac
# Notes on getting ttf files on Mac
# * Location in /System/Library/Font, /Library/Fonts, or ~/Library/Fonts
# * To break down .dfont files, use fondu (homebrew download).
#   To break down .ttc files, use dfontsplitter (https://peter.upfold.org.uk/projects/dfontsplitter)
#   To break down .bdf files made by fondu, use mkttf (https://github.com/Tblue/mkttf; requires FontForge and PoTrace)
# * Install new fonts with "brew cask install font-<name-of-font>" after using
#   "brew tap caskroom/fonts" to initialize; appear in ~/Library/Fonts; see https://github.com/Homebrew/homebrew-cask-fonts
# * The .otf files work in addition to .ttf files. You can verify this by
#   looking at plot.fonts_files_os -- it includes .otf files.
# Notes on default files packaged in font directory:
# * Location will be something like:
#   /lib/python3.6/site-packages/matplotlib/mpl-data/fonts/ttf
# * 'STIX' fonts allow different LaTeX-like math modes e.g. blackboard bold
#   and caligraphy; see: https://matplotlib.org/gallery/text_labels_and_annotations/stix_fonts_demo.html
# * The 'cm'-prefix fonts seem to provide additional mathematical symbols
#   like integrals, and italized math-mode fonts.
# * We also have 'pdfcorefonts' in this directory, but I think since these
#   are afm matplotlib cannot use them? Don't know.
# Helvetica: https://olgabotvinnik.com/blog/2012-11-15-how-to-set-helvetica-as-the-default-sans-serif-font-in/
# Classic fonts: https://www.lifewire.com/classic-sans-serif-fonts-clean-appearance-1077406
# Good for downloading fonts: https://www.cufonfonts.com
import os
import re
import sys
import shutil
import glob
import matplotlib.font_manager as mfonts
from matplotlib import matplotlib_fname, get_cachedir
_data_fonts = os.path.join(os.path.dirname(__file__), 'fonts') # proplot fonts
_data_user = os.path.join(os.path.expanduser('~'), '.proplot')
_data_user_fonts = os.path.join(_data_user, 'fonts') # user fonts
_data_matplotlib = re.sub('/matplotlibrc$', '', matplotlib_fname())
_data_matplotlib_fonts = os.path.join(_data_matplotlib, 'fonts', 'ttf')
if not os.path.isdir(_data_user):
    os.mkdir(_data_user)
if not os.path.isdir(_data_user_fonts):
    os.mkdir(_data_user_fonts)

# Register fonts
# https://github.com/olgabot/sciencemeetproductivity.tumblr.com/blob/master/posts/2012/11/how-to-set-helvetica-as-the-default-sans-serif-font-in.md
# https://olgabotvinnik.com/blog/2012-11-15-how-to-set-helvetica-as-the-default-sans-serif-font-in/
# https://stackoverflow.com/questions/18821795/how-can-i-get-list-of-font-familyor-name-of-font-in-matplotlib
def register_fonts():
    """Adds fonts packaged with ProPlot or saved to the ``~/.proplot/fonts``
    folder. Also deletes the font cache, which may cause delays.
    Detects ``.ttf`` and ``.otf`` files -- see `this link
    <https://gree2.github.io/python/2015/04/27/python-change-matplotlib-font-on-mac>`__
    for a guide on converting various other font file types to ``.ttf`` and
    ``.otf`` for use with matplotlib."""
    # Populate file and font lists
    fonts[:] = []
    for (ifiles,ifonts) in ((fonts_os_files,fonts_os), (fonts_mpl_files,fonts_mpl)):
        ifonts[:] = []
        if ifiles is fonts_os_files:
            ifiles[:] = sorted(mfonts.findSystemFonts(fontpaths=None, fontext='ttf'))
        else:
            ifiles[:] = sorted(glob.glob(os.path.join(_data_matplotlib, 'fonts', 'ttf', '*.[ot]tf')))
        for file in ifiles:
            try:
                font = mfonts.FontProperties(fname=file).get_name()
            except Exception as err:
                pass
            else:
                ifonts.append(font)
        seen = {*()}
        ifonts[:] = [font for font in ifonts if font not in seen and not seen.add(font)]
        fonts.extend(ifonts)
    fonts[:] = sorted(fonts)

    # Transfer files to matplotlibrc
    fonts_new = []
    for filename in sorted(glob.glob(os.path.join(_data_fonts, '*.[ot]tf'))) + \
            sorted(glob.glob(os.path.join(_data_user_fonts, '*.[ot]tf'))):
        base = os.path.basename(filename)
        if os.path.exists(os.path.join(_data_matplotlib_fonts, base)):
            continue
        shutil.copy(filename, _data_matplotlib_fonts)
        fonts_new.append(base)

    # Delete and rebuild cache
    # TODO: Is rebuilding the cache sufficient?
    if fonts_new:
        print(f'Installed font(s): {", ".join(fonts_new)}.')
        caches = []
        dir_cache = get_cachedir()
        for file in glob.glob(os.path.join(dir_cache, '*.cache')) + glob.glob(os.path.join(dir_cache, '*.json')):
            if os.path.isdir(file): # important to dump even tex.cache, get weird errors with wrong font weight!
                shutil.rmtree(file)
            else:
                os.remove(file)
            caches.append(file)
        if caches:
            print(f'Deleted font cache(s): {", ".join(caches)}.')
        mfonts._rebuild()

# Font lists
fonts_mpl_files = []
"""Font filenames provided by matplotlib or ProPlot."""
fonts_os_files = []
"""Font filenames provided by your operating system."""
fonts_mpl = []
"""Registered font names provided by matplotlib or ProPlot."""
fonts_os = []
"""Registered font names provided by your operating system."""
fonts = []
"""All registered font names."""
_missing_fonts = {*()}
"""Missing fonts, filled whenever user requests a bad name."""

# Register fonts
register_fonts()

