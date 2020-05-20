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
#
# Color names
# ===========
#
# ProPlot registers several new color names and includes tools for defining
# your own color names. These features are described below.


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_colors:
#
# Included colors
# ---------------
#
# ProPlot adds new color names from the `XKCD color survey
# <https://blog.xkcd.com/2010/05/03/color-survey-results/>`__  and
# the `Open Color <https://github.com/yeun/open-color>`__ UI design color
# palettes. You can use `~proplot.demos.show_colors` to generate a table of these
# colors. Note that the matplotlib's native `X11 named colors
# <https://matplotlib.org/examples/color/named_colors.html>`__ are still
# registered, but some of the X11 color names may be overwritten by the XKCD names,
# and we encourage choosing colors from the below tables instead. ProPlot
# registers XKCD colors because the selection is larger and the names are
# more likely to match your intuition for what a color "should" look like.
#
# To reduce the number of registered color names to a more manageable size,
# ProPlot filters the available XKCD colors so that they are *sufficiently
# distinct* in the :ref:`perceptually uniform colorspace <ug_perceptual>`.
# This makes it a bit easier to pick out colors from the table generated with
# `~proplot.demos.show_colors`. Similar names were also cleaned up -- for
# example, ``'reddish'`` and ``'reddy'`` are changed to ``'red'``.

# %%
import proplot as plot
fig, axs = plot.show_colors()


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
fig, axs = plot.subplots(nrows=2, axwidth=3.2, share=0)
axs.format(
    xformatter='null', yformatter='null', abc=True, abcloc='ul', abcstyle='A.',
    suptitle='Getting individual colors from colormaps and cycles'
)

# Drawing from colormap
ax = axs[0]
name = 'Deep'
cmap = plot.Colormap(name)
idxs = plot.arange(0, 1, 0.2)
state.shuffle(idxs)
for idx in idxs:
    data = (state.rand(20) - 0.4).cumsum()
    h = ax.plot(
        data, lw=5, color=(name, idx),
        label=f'idx {idx:.1f}', legend='r', legend_kw={'ncols': 1}
    )
ax.colorbar(cmap, loc='ur', label='colormap', length='12em')
ax.format(title='Drawing from the Solar colormap', grid=True)

# Drawing from color cycle
ax = axs[1]
idxs = np.arange(6)
state.shuffle(idxs)
for idx in idxs:
    data = (state.rand(20) - 0.4).cumsum()
    h = ax.plot(
        data, lw=5, color=('qual1', idx),
        label=f'idx {idx:.0f}', legend='r', legend_kw={'ncols': 1}
    )
ax.format(title='Drawing from the ggplot color cycle')


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
# value. Lines beginning with ``#`` are ignored as comments.
