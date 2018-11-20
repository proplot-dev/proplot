#!/usr/bin/env python3
"""
Script that simply lists the availble system fonts.
Add to this.
"""
import os
from glob import glob
from subprocess import Popen, PIPE
from matplotlib import matplotlib_fname
import matplotlib.font_manager as mfonts
#------------------------------------------------------------------------------
# List the current font names, original version; works on Linux but not Mac,
# because can't find mac system fonts?
#------------------------------------------------------------------------------#
# List the system font names, smarter version
# See: https://github.com/olgabot/sciencemeetproductivity.tumblr.com/blob/master/posts/2012/11/how-to-set-helvetica-as-the-default-sans-serif-font-in.md
# Also see: https://olgabotvinnik.com/blog/2012-11-15-how-to-set-helvetica-as-the-default-sans-serif-font-in/
# Also see: https://stackoverflow.com/questions/18821795/how-can-i-get-list-of-font-familyor-name-of-font-in-matplotlib
# fonts_mpl = sorted({os.path.basename(font.rstrip('.ttf')) for font in # hack-installed fonts
#     glob(f"{matplotlib_fname().rstrip('matplotlibrc')}/fonts/ttf/*.ttf")})
# fonts_os  = sorted({os.path.basename(font.split('.')[0]) for font in # system fonts
#     mfonts.findSystemFonts(fontpaths=None, fontext='ttf')})
# fonts_ttf = sorted({f.name for f in mfonts.fontManager.ttflist})
# fonts_afm = sorted({f.name for f in mfonts.fontManager.afmlist})
# fonts = mfonts.get_fontconfig_fonts() # same thing? but this lets us segregate them
fonts_mpl_files = sorted(glob(f"{matplotlib_fname().rstrip('matplotlibrc')}/fonts/ttf/*.ttf"))
fonts_os_files  = sorted(mfonts.findSystemFonts(fontpaths=None, fontext='ttf'))
fonts_os, fonts_mpl = set(), set()
for file in fonts_os_files:
    try:
        fonts_os.add(mfonts.FontProperties(fname=file).get_name())
    except Exception as err:
        pass # fails sometimes
for file in fonts_mpl_files:
    try:
        fonts_mpl.add(mfonts.FontProperties(fname=file).get_name())
    except Exception as err:
        pass # fails sometimes
fonts = {*fonts_os, *fonts_mpl}

# Missing fonts (add to this list whenever user requests one)
_missing_fonts = []

# Optionally update system fonts
_font_script = f'{os.path.dirname(__file__)}/fontscript.sh'
def install_fonts():
    """
    Install matplotlib fonts from ttf files located in the 'fonts' directory.
    May require restarting iPython session. Note font cache will be deleted
    in this process, which could cause delays.
    """
    # See: https://stackoverflow.com/a/17413045/4970632
    # os.system(_font_script)
    # p = Popen(_font_script, stdout=PIPE).wait()
    p = Popen(_font_script, shell=False).wait()
    # while p.poll() is None:
    #     print(p.stdout.readline()) # this blocks until it receives a newline.
    print('Fonts have ben installed and font cache has been emptied. Please restart your iPython session.')
