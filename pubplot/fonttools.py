#!/usr/bin/env python3
"""
Script that simply lists the availble system fonts.
Add to this.
"""
import os
from glob import glob
from matplotlib import matplotlib_fname
import matplotlib.font_manager as mfonts
#------------------------------------------------------------------------------
# List the current font names, original version; works on Linux but not Mac,
# because can't find mac system fonts?
#------------------------------------------------------------------------------#
# List the system font names, smarter version
# See: https://olgabotvinnik.com/blog/2012-11-15-how-to-set-helvetica-as-the-default-sans-serif-font-in/
# fonts = mfonts.get_fontconfig_fonts()
# fonts = [mfonts.FontProperties(fname=fname).get_name() for fname in flist]
fonts = [font.split('/')[-1].split('.')[0] for font in # system fonts
            mfonts.findSystemFonts(fontpaths=None, fontext='ttf')] + \
        [os.path.basename(font.rstrip('.ttf')) for font in # hack-installed fonts
            glob(f"{matplotlib_fname().rstrip('matplotlibrc')}/fonts/ttf/*.ttf")]
fonts = sorted(set(fonts)) # unique ones only
