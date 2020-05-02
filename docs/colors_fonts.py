# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_colors_fonts:
#
# Colors and fonts
# ================
#
# ProPlot registers several new color names and font families, and includes
# tools for defining your own color names and adding your own font families.
# These features are described below.


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_colors:
#
# Included colors
# ---------------
#
# ProPlot adds new color names from the `XKCD color survey
# <https://blog.xkcd.com/2010/05/03/color-survey-results/>`__  and
# the `"Open color" <https://github.com/yeun/open-color>`__ Github project.
# This was inspired by `seaborn
# <https://seaborn.pydata.org/tutorial/color_palettes.html>`__. Use
# `~proplot.show.show_colors` to generate tables of these colors. Note that
# the matplotlib's native `X11 named colors
# <https://matplotlib.org/examples/color/named_colors.html>`__ are still
# registered, but we encourage using colors from the tables instead, and
# some of the X11 color names may be overwritten by the XKCD names. We prefer
# the XKCD color names because the selection is larger and the names are
# more likely to match your intuition for what a color "should" look like.
#
# To reduce the number of registered color names to a more manageable size,
# ProPlot filters the available XKCD colors so that they are *sufficiently
# distinct* in the :ref:`perceptually uniform colorspace <ug_perceptual>`.
# This makes it a bit easier to pick out colors from the table generated with
# `~proplot.show.show_colors`. Similar names were also cleaned up -- for
# example, ``'reddish'`` and ``'reddy'`` are changed to ``'red'``.

# %%
import proplot as plot
figs = plot.show_colors()


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_colors_cmaps:
#
# Colors from colormaps
# ---------------------
#
# If you want to draw an individual color from a colormap or a color cycle,
# use ``color=(cmap, coord)`` or ``color=(cycle, index)`` with any command
# that accepts the `color` keyword. The ``coord`` should be between ``0`` and
# ``1``, while the ``index`` is the index on the list of cycle colors. This
# feature is powered by the `~proplot.colors.ColorDatabase` class. This is
# useful if you spot a nice color in one of the available colormaps and want
# to use it for some arbitrary plot element.

# %%
import proplot as plot
import numpy as np
state = np.random.RandomState(51423)
fig, axs = plot.subplots(nrows=2, aspect=2, axwidth=3, share=0)

# Drawing from colormap
ax = axs[0]
cmap = 'Deep'
m = ax.pcolormesh([[0], [1]], cmap=cmap, N=1000)
idxs = plot.arange(0, 1, 0.2)
state.shuffle(idxs)
for idx in idxs:
    h = ax.plot(
        (np.random.rand(20) - 0.4).cumsum(), lw=5, color=(cmap, idx),
        label=f'idx {idx:.1f}', legend='r', legend_kw={'ncols': 1}
    )
ax.colorbar(m, loc='ul', locator=0.2, label='colormap')
ax.format(title='Drawing from the Solar colormap', grid=True)

# Drawing from color cycle
ax = axs[1]
idxs = np.arange(6)
state.shuffle(idxs)
for idx in idxs:
    h = ax.plot(
        (np.random.rand(20) - 0.4).cumsum(), lw=5, color=('qual1', idx),
        label=f'idx {idx:.0f}', legend='r', legend_kw={'ncols': 1}
    )
ax.format(title='Drawing from the ggplot color cycle')
axs.format(
    xlocator='null', abc=True, abcloc='ur', abcstyle='A.',
    suptitle='Getting individual colors from colormaps and cycles'
)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_colors_user:
#
# Using your own colors
# ---------------------
#
# You can register your own colors by adding ``.txt`` files to the
# ``~/.proplot/colors`` directory and calling
# `~proplot.config.register_colors`. This command is also called on import.
# Each file should contain lines that look like ``color: #xxyyzz`` where
# ``color`` is the registered color name and ``#xxyyzz`` is the HEX color
# value.


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_fonts:
#
# Included fonts
# --------------
#
# ProPlot adds several open source fonts, including the
# `TeX Gyre <http://www.gust.org.pl/projects/e-foundry/tex-gyre>`__ font
# series, and introduces a `~proplot.show.show_fonts` command to compare
# fonts. By default, this command displays the *sans-serif* fonts packaged
# with ProPlot and available on your system (see `~matplotlib.font_manager`).
# Generally speaking, sans-serif fonts are more appropriate for figures than
# serif fonts.
#
# ProPlot also changes the default font to the Helvetica-lookalike
# `TeX Gyre Heros <http://www.gust.org.pl/projects/e-foundry/tex-gyre/heros>`__.
# Matplotlib uses `DejaVu Sans <https://dejavu-fonts.github.io>`__ in part
# because it includes glyphs for a wider range of mathematical symbols (where
# you see the “¤” dummy symbol in the below table, that character is
# unavailable), but IMHO TeX Gyre Heros is much more aesthetically pleasing.
# If your plot has lots of symbols, you may want to switch to DejaVu Sans or
# `Fira Math <https://github.com/firamath/firamath>`__ (which is also
# packaged with ProPlot).

# %%
import proplot as plot
fig = plot.show_fonts()


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_fonts_user:
#
# Using your own fonts
# --------------------
#
# You can register your own fonts by adding files to the ``~/.proplot/fonts``
# directory and calling `~proplot.config.register_fonts`. This command is
# also called on import. To change the default font, use the
# `~proplot.config.rc` object or modify your ``~/.proplotrc``. See
# the :ref:`configuration section <ug_config>` for details.
#
# Sometimes the font you would like to use *is* installed, but the font file
# is not stored under the matplotlib-compatible ``.ttf``, ``.otf``, or
# ``.afm`` formats. For example, several macOS fonts are unavailable because
# they are stored as ``.dfont`` collections. Also, while matplotlib nominally
# supports ``.ttc`` collections, ProPlot manually removes them because
# figures with ``.ttc`` fonts
# `cannot be saved as PDFs <https://github.com/matplotlib/matplotlib/issues/3135>`__.
# You can get matplotlib to use these fonts by expanding the "collections"
# into individual ``.ttf`` files with the
# `DFontSplitter application <https://peter.upfold.org.uk/projects/dfontsplitter>`__,
# then saving the files in-place or in the ``~/.proplot/fonts`` folder.
#
# To find font files, check the paths listed in ``OSXFontDirectories``,
# ``X11FontDirectories``, ``MSUserFontDirectories``, and ``MSFontDirectories``
# under the `~matplotlib.font_manager` module. Note that if the font in question has
# a "thin" style, implied by file names with the word ``Thin``,
# `a matplotlib bug <https://github.com/matplotlib/matplotlib/issues/8788>`__
# may cause these styles to override the "normal" style.
