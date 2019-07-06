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
# WARNING: Beware of unrecognized font variant 'Thin'! Matplotlib gets random
# unsorted list of font files when rebuilding the manager via os.walk(),
# which means if an unrecognized variant is added before the normal one,
# matplotlib uses that font! Found out by printing mfonts.fontManager.ttflist
# and noticing the 'Thin' files registered as 'normal' for all properties. See
# https://matplotlib.org/api/font_manager_api.html for valid weights, styles, etc.
# print(*[font for font in mfonts.fontManager.ttflist if 'HelveticaNeue' in font.fname], sep='\n')
# print(*[font.fname for font in mfonts.fontManager.ttflist if 'HelveticaNeue' in font.fname], sep='\n')
import os
import re
import sys
import shutil
import glob
import matplotlib.font_manager as mfonts
from matplotlib import get_cachedir, get_data_path
_data_fonts = os.path.join(os.path.dirname(__file__), 'fonts') # proplot fonts
_data_user = os.path.join(os.path.expanduser('~'), '.proplot')
_data_user_fonts = os.path.join(_data_user, 'fonts') # user fonts
_data_matplotlib_fonts = os.path.join(get_data_path(), 'fonts', 'ttf')
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
            ifiles[:] = sorted(glob.glob(os.path.join(_data_matplotlib_fonts, '*.[ot]tf')))
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

    # Transfer files to mpl-data
    fonts_new = []
    for filename in sorted(glob.glob(os.path.join(_data_fonts, '*.[ot]tf'))) + \
            sorted(glob.glob(os.path.join(_data_user_fonts, '*.[ot]tf'))):
        base = os.path.basename(filename)
        if os.path.exists(os.path.join(_data_matplotlib_fonts, base)):
            continue
        shutil.copy(filename, _data_matplotlib_fonts)
        fonts_new.append(base)

    # Rebuild cache
    if fonts_new:
        print(f'Installed font(s): {", ".join(fonts_new)}.')
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

#------------------------------------------------------------------------------#
# Fonts demo
#------------------------------------------------------------------------------#
def show_fonts(size=12):
    """Display nicely-formatted table of the fonts available in the matplotlib
    mpl-data folder."""
    from . import subplots
    fonts = ['DejaVu Sans', 'Arial', 'Avenir', 'Franklin Gothic Book', 'Frutiger', 'Futura',
            'Gotham', 'Helvetica', 'Helvetica Neue', 'Geneva', 'Gill Sans',
            'Lucida Grande', 'Noto Sans', 'Myriad Pro', 'Open Sans', 'Optima', 'Tahoma', 'Trebuchet MS', 'Univers', 'Verdana']
    math = r'(0) + {1} - [2] * <3> / 4,0 $\geq\gg$ 5.0 $\leq\ll$ ~6 $\times$ 7 $\equiv$ 8 $\approx$ 9 $\propto$'
    greek = r'$\alpha\beta$ $\Gamma\gamma$ $\Delta\delta$ $\epsilon\zeta\eta$ $\Theta\theta$ $\kappa\mu\nu$ $\Lambda\lambda$ $\Pi\pi$ $\xi\rho\tau\chi$ $\Sigma\sigma$ $\Phi\phi$ $\Psi\psi$ $\Omega\omega$ !?&#%'
    letters = 'the quick brown fox jumps over a lazy dog\nTHE QUICK BROWN FOX JUMPS OVER A LAZY DOG'
    # letters = 'Aa Bb Cc Dd Ee Ff Gg Hh Ii Jj Kk Ll Mm Nn Oo Pp Qq Rr Ss Tt Uu Vv Ww Xx Yy Zz'
    for weight in ('normal',):
        f, axs = subplots(ncols=1, nrows=len(fonts), flush=True, axwidth=4.5, axheight=5.5*size/72)
        axs.format(xloc='neither', yloc='neither', xlocator='null', ylocator='null', alpha=0)
        axs[0].format(title='Fonts demo', titlefontsize=size, titleloc='l', titleweight='bold')
        for i,ax in enumerate(axs):
            font = fonts[i]
            plot.rc.fontname = font
            ax.text(0, 0.5, f'{font}: {letters}\n{math}\n{greek}',
                    fontsize=size, weight=weight, ha='left', va='center')
    return f
