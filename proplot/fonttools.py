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
# Valid styles: https://matplotlib.org/api/font_manager_api.html for valid weights, styles, etc.
# Classic fonts: https://www.lifewire.com/classic-sans-serif-fonts-clean-appearance-1077406
# Good for downloading fonts: https://www.cufonfonts.com
# WARNING: Check out ttflist whenever adding new ttf files! For example, realized
# could dump all of the Gotham-Name.ttf files instead of GothamName files, and
# got Helvetica bug due to unrecognized 'thin' font style overwriting normal one.
# print(*[font for font in mfonts.fontManager.ttflist if 'HelveticaNeue' in font.fname], sep='\n')
# print(*[font.fname for font in mfonts.fontManager.ttflist if 'HelveticaNeue' in font.fname], sep='\n')
import os
import shutil
import glob
import matplotlib.font_manager as mfonts
from .utils import _check_data
from matplotlib import get_data_path
__all__ = ['clean_fonts', 'register_fonts', 'show_fonts', 'fonts', 'fonts_system', 'fonts_proplot']

# Data
_data_user = os.path.join(os.path.expanduser('~'), '.proplot')
_data_user_fonts = os.path.join(_data_user, 'fonts') # user fonts
_data_fonts = os.path.join(os.path.dirname(__file__), 'fonts') # proplot fonts
if not os.path.isdir(_data_user):
    os.mkdir(_data_user)
if not os.path.isdir(_data_user_fonts):
    os.mkdir(_data_user_fonts)

def register_fonts():
    """Adds fonts packaged with ProPlot or saved to the ``~/.proplot/fonts``
    folder. Also deletes the font cache, which may cause delays.
    Detects ``.ttf`` and ``.otf`` files -- see `this link
    <https://gree2.github.io/python/2015/04/27/python-change-matplotlib-font-on-mac>`__
    for a guide on converting various other font file types to ``.ttf`` and
    ``.otf`` for use with matplotlib."""
    # Add proplot path to TTFLIST and rebuild cache
    _check_data()
    paths = _data_fonts + ':' + _data_user_fonts
    if 'TTFPATH' not in os.environ:
        os.environ['TTFPATH'] = paths
    elif paths not in os.environ['TTFPATH']:
        os.environ['TTFPATH'] += (':' + paths)
    mfonts._rebuild()
    # Populate font lists
    fonts_system[:] = sorted({font.name for font in mfonts.fontManager.ttflist if not (_data_user_fonts in font.fname or _data_fonts in font.fname)})
    fonts_proplot[:] =  sorted({font.name for font in mfonts.fontManager.ttflist if (_data_user_fonts in font.fname or _data_fonts in font.fname)})
    fonts[:] = sorted((*fonts_system, *fonts_proplot))

def clean_fonts():
    """Remove fonts from ``mpl-data`` that were added by ProPlot."""
    # TODO: Delete this when enough time has passed, no longer copy stuff to mpl-data
    rm = []
    data_matplotlib_fonts = os.path.join(get_data_path(), 'fonts', 'ttf')
    fonts = {os.path.basename(font) for font in glob.glob(os.path.join(_data_fonts, '*'))}
    fonts_mpl = sorted(glob.glob(os.path.join(data_matplotlib_fonts, '*.[ot]tf')))
    for font_mpl in fonts_mpl:
        if os.path.basename(font_mpl) in fonts:
            rm.append(font_mpl)
            os.remove(font_mpl)
    print(f'Removed fonts {", ".join(os.path.basename(font) for font in rm)}.')
    mfonts._rebuild()

# Font lists
fonts_proplot = []
"""Names of fonts added by ProPlot."""
fonts_system = []
"""Names of fonts provided by matplotlib or your operating system."""
fonts = []
"""All registered font names."""

# Register fonts
register_fonts()

#------------------------------------------------------------------------------#
# Fonts demo
#------------------------------------------------------------------------------#
def show_fonts(fonts=None, size=12):
    """Displays table of the fonts installed by ProPlot or in the user-supplied
    `fonts` list. Use `size` to change the fontsize for fonts shown in the figure."""
    # letters = 'Aa Bb Cc Dd Ee Ff Gg Hh Ii Jj Kk Ll Mm Nn Oo Pp Qq Rr Ss Tt Uu Vv Ww Xx Yy Zz'
    from . import subplots
    fonts = ('DejaVu Sans', *fonts_proplot)
    math = r'(0) + {1} - [2] * <3> / 4,0 $\geq\gg$ 5.0 $\leq\ll$ ~6 $\times$ 7 $\equiv$ 8 $\approx$ 9 $\propto$'
    greek = r'$\alpha\beta$ $\Gamma\gamma$ $\Delta\delta$ $\epsilon\zeta\eta$ $\Theta\theta$ $\kappa\mu\nu$ $\Lambda\lambda$ $\Pi\pi$ $\xi\rho\tau\chi$ $\Sigma\sigma$ $\Phi\phi$ $\Psi\psi$ $\Omega\omega$ !?&#%'
    letters = 'the quick brown fox jumps over a lazy dog\nTHE QUICK BROWN FOX JUMPS OVER A LAZY DOG'
    for weight in ('normal',):
        f, axs = subplots(ncols=1, nrows=len(fonts), flush=True, axwidth=4.5, axheight=5.5*size/72)
        axs.format(xloc='neither', yloc='neither', xlocator='null', ylocator='null', alpha=0)
        axs[0].format(title='Fonts demo', titlefontsize=size, titleloc='l', titleweight='bold')
        for i,ax in enumerate(axs):
            font = fonts[i]
            ax.text(0, 0.5, f'{font}: {letters}\n{math}\n{greek}', fontfamily=font,
                    fontsize=size, weight=weight, ha='left', va='center')
    return f
