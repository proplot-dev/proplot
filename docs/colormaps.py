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
# .. _ug_cmaps:
#
# Colormaps
# =========
#
# ProPlot defines **colormaps** as color palettes that sample some
# *continuous function* between two end colors. Colormaps are generally used
# to encode data values on a pseudo-third dimension. They are are implemented
# with the `~proplot.colors.LinearSegmentedColormap` and
# `~proplot.colors.PerceptuallyUniformColormap` classes, which are
# :ref:`subclassed from <ug_cmaps_new>`
# `matplotlib.colors.LinearSegmentedColormap`.
#
# ProPlot adds several features to help you use colormaps effectively in your
# figures. This section documents the new registered colormaps, explains how
# to make and modify colormaps, and shows how to apply them to your plots.


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_cmaps_included:
#
# Included colormaps
# ------------------
#
# On import, ProPlot registers a few sample
# :ref:`perceptually uniform colormaps <ug_perceptual>`, plus several
# colormaps from other online data viz projects. Use
# `~proplot.demos.show_cmaps` to generate a table of registered maps. The
# figure is broken down into the following sections:
#
# * "User" colormaps, i.e. colormaps saved to your ``~/.proplot/cmaps``
#   folder. You can use ProPlot to :ref:`make new colormaps <ug_cmaps_new>` and
#   add them to this folder.
# * Matplotlib and seaborn original colormaps.
# * ProPlot original :ref:`perceptually uniform colormaps <ug_perceptual>`.
# * `cmOcean <https://matplotlib.org/cmocean/>`__ colormaps designed for
#   oceanographic data but useful for everyone.
# * Fabio Crameri's
#   `"scientific colour maps" <http://www.fabiocrameri.ch/colourmaps.php>`__.
# * Cynthia Brewer's `ColorBrewer <http://colorbrewer2.org/>`__ colormaps, included
#   with matplotlib by default.
# * Colormaps from the `SciVisColor <https://sciviscolor.org/home/colormoves/>`__
#   online interactive tool. There are so many of these colormaps because they
#   are intended to be *merged* into more complex colormaps.
#
# ProPlot removes some default matplotlib colormaps with erratic color
# transitions.
#
# .. note::
#
#    Colormap and color cycle identification is more flexible in ProPlot. The names
#    are are case-insensitive (e.g. ``'Viridis'``, ``'viridis'``, and ``'ViRiDiS'``
#    are equivalent), diverging colormap names can be specified in their "reversed"
#    form (e.g. ``'BuRd'`` is equivalent to ``'RdBu_r'``), and appending ``'_r'``
#    or ``'_s'`` to *any* colormap name will return a
#    `~proplot.colors.LinearSegmentedColormap.reversed` or
#    `~proplot.colors.LinearSegmentedColormap.shifted` version of the colormap
#    or color cycle. See `~proplot.colors.ColormapDatabase` for more info.

# %%
import proplot as plot
fig, axs = plot.show_cmaps()


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_perceptual:
#
# Perceptually uniform colormaps
# ------------------------------
#
# ProPlot's custom colormaps are instances of the
# `~proplot.colors.PerceptuallyUniformColormap` class. These colormaps
# generate colors by interpolating between coordinates in any of the
# following three colorspaces:
#
# * **HCL** (a.k.a. `CIE LChuv <https://en.wikipedia.org/wiki/CIELUV>`__):
#   A purely perceptually uniform colorspace, where colors are broken down
#   into “hue” (color, range 0-360), “chroma” (saturation, range 0-100), and
#   “luminance” (brightness, range 0-100). This space is difficult to work
#   with due to *impossible colors* -- colors that, when translated back from
#   HCL to RGB, result in RGB channels greater than ``1``.
# * **HPL** (a.k.a. `HPLuv <http://www.hsluv.org/comparison>`__): Hue and
#   luminance are identical to HCL, but 100 saturation is set to the minimum
#   maximum saturation *across all hues* for a given luminance. HPL restricts
#   you to soft pastel colors, but is closer to HCL in terms of uniformity.
# * **HSL** (a.k.a. `HSLuv <http://www.hsluv.org/comparison>`__): Hue and
#   luminance are identical to HCL, but 100 saturation is set to the maximum
#   saturation *for a given hue and luminance*. HSL gives you access to the
#   entire RGB colorspace, but often results in sharp jumps in chroma.
#
# The colorspace used by each `~proplot.colors.PerceptuallyUniformColormap`
# is set with the `space` keyword arg. To plot arbitrary cross-sections of
# these colorspaces, use `~proplot.demos.show_colorspaces` (the black
# regions represent impossible colors). To see how colormaps vary with
# respect to each channel, use `~proplot.demos.show_channels`. Some examples
# are shown below.
#
# In theory, "uniform" colormaps should have *straight* lines in hue, chroma,
# and luminance (second figure, top row). In practice, this is
# difficult to accomplish due to impossible colors. Matplotlib and seaborn's
# ``'magma'`` and ``'Rocket'`` colormaps are fairly linear with respect to
# hue and luminance, but not chroma. ProPlot's ``'Fire'`` is linear in hue,
# luminance, and *HSL saturation* (bottom left), while ``'Dusk'`` is linear
# in hue, luminance, and *HPL saturation* (bottom right).

# %%
# Colorspace demo
import proplot as plot
fig, axs = plot.show_colorspaces(axwidth=1.6, luminance=50)
fig, axs = plot.show_colorspaces(axwidth=1.6, saturation=60)
fig, axs = plot.show_colorspaces(axwidth=1.6, hue=0)

# %%
# Compare colormaps
import proplot as plot
for cmaps in (('magma', 'rocket'), ('fire', 'dusk')):
    fig, axs = plot.show_channels(
        *cmaps, axwidth=1.5, minhue=-180, maxsat=400, rgb=False
    )


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_cmaps_new:
#
# Making new colormaps
# --------------------
#
# ProPlot doesn't just include new colormaps -- it provides tools for merging
# colormaps, modifying colormaps, making :ref:`perceptually uniform colormaps
# <ug_perceptual>` from scratch, and saving the results for future use. For
# your convenience, most of these features can be accessed via the
# `~proplot.constructor.Colormap` constructor function. Note that every
# plotting command that accepts a `cmap` keyword passes it through this
# function (see `~proplot.axes.cmap_changer`).
#
# To make `~proplot.colors.PerceptuallyUniformColormap`\ s from scratch, you
# have the following three options:
#
# * Pass a color name, hex string, or RGB tuple to
#   `~proplot.constructor.Colormap`. This builds a *monochromatic* (single
#   hue) colormap by calling the
#   `~proplot.colors.PerceptuallyUniformColormap.from_color` static method.
#   The colormap colors will vary from the specified color to pure white or
#   some shade *near* white (see the `fade` keyword arg).
# * Pass a *list of colors* to `~proplot.constructor.Colormap`. This calls
#   the `~proplot.colors.PerceptuallyUniformColormap.from_list` static method,
#   which linearly interpolates between each color in hue, saturation, and
#   luminance.
# * Pass a *dictionary* to `~proplot.constructor.Colormap`. This calls the
#   `~proplot.colors.PerceptuallyUniformColormap.from_hsl` static method,
#   which draws lines between channel values specified by the keyword arguments
#   `hue`, `saturation`, and `luminance`. The values can be numbers, color
#   strings, or lists thereof. Numbers indicate the channel value. For color
#   strings, the channel value is *inferred* from the specified color. You can
#   end any color string with ``'+N'`` or ``'-N'`` to *offset* the channel
#   value by the number ``N``.
#
# In the below example, we use all of these methods to make brand new
# `~proplot.colors.PerceptuallyUniformColormap`\ s in the ``'hsl'`` and
# ``'hpl'`` colorspaces.

# %%
import proplot as plot
import numpy as np
state = np.random.RandomState(51423)
data = state.rand(30, 30).cumsum(axis=1)

# Initialize figure
fig, axs = plot.subplots(
    [[0, 1, 1, 2, 2, 0], [3, 3, 4, 4, 5, 5]],
    ncols=2, axwidth=2, aspect=1
)
axs.format(
    xlabel='x axis', ylabel='y axis',
    suptitle='Building your own PerceptuallyUniformColormaps'
)

# Simple monochromatic colormap
axs[0].format(title='From single color')
m = axs[0].contourf(data, cmap='ocean blue', cmap_kw={'name': 'water'})
cmap1 = m.cmap

# Three merged monochromatic colormaps
# The trailing '_r' make the colormap go dark-to-light instead of light-to-dark
axs[1].format(title='From three colors')
cmap2 = plot.Colormap(
    'dark red_r', 'denim_r', 'warm gray_r',
    fade=90, name='tricolor'
)
axs[1].contourf(data, cmap=cmap2, levels=12)

# Colormaps from channel value dictionaries
axs[2:4].format(title='From channel values')
cmap3 = plot.Colormap(
    {
        'hue': ['red-90', 'red+90'],
        'saturation': [50, 70, 30],
        'luminance': [20, 100]
    },
    name='Matter',
    space='hcl',
)
axs[2].pcolormesh(data, cmap=cmap3)
cmap4 = plot.Colormap(
    {
        'hue': ['red', 'red-720'],
        'saturation': [80, 20],
        'luminance': [20, 100]
    },
    name='cubehelix',
    space='hpl',
)
axs[3].pcolormesh(data, cmap=cmap4)

# Colormap from lists
m = axs[4].pcolormesh(
    data,
    cmap=('maroon', 'ivory'),
    cmap_kw={'name': 'heat'}
)
cmap5 = m.cmap
axs[4].format(title='From list of colors')

# Display the channels
fig, axs = plot.show_channels(cmap1, cmap2, axwidth=1.5, rgb=False)
fig, axs = plot.show_channels(
    cmap3, cmap4, cmap5, minhue=-180, axwidth=1.5, rgb=False
)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_cmaps_merge:
#
# Merging colormaps
# -----------------
#
# To *merge* colormaps, simply pass multiple positional arguments to the
# `~proplot.constructor.Colormap` constructor. Each positional argument can
# be a colormap name, a colormap instance, or a
# :ref:`special argument <ug_cmaps_new>` that generates a new colormap
# on-the-fly. This lets you create new diverging colormaps and segmented
# `SciVisColor <https://sciviscolor.org/home/colormoves/>`__ style colormaps
# right inside ProPlot. Segmented colormaps are often desirable for complex
# datasets with complex statistical distributions.
#
# In the below example, we create a new divering colormap and reconstruct the
# colormap from `this SciVisColor example\
# <https://sciviscolor.org/wp-content/uploads/sites/14/2018/04/colormoves-icon-1.png>`__.
# We also *save* the results for future use by passing ``save=True`` to
# `~proplot.constructor.Colormap`.

# %%
import proplot as plot
import numpy as np
fig, axs = plot.subplots([[0, 1, 1, 0], [2, 2, 3, 3]], axwidth=2.4, span=False)
state = np.random.RandomState(51423)
data = state.rand(30, 30).cumsum(axis=1)

# Diverging colormap example
title1 = 'Custom diverging map'
cmap1 = plot.Colormap('Blue4_r', 'RedPurple3', name='Diverging', save=True)

# SciVisColor examples
title2 = 'Custom complex map'
cmap2 = plot.Colormap(
    'Green1_r', 'Orange5', 'Blue1_r', 'Blue6',
    name='Complex', save=True
)
title3 = 'SciVisColor example reproduction'
cmap3 = plot.Colormap(
    'Green1_r', 'Orange5', 'Blue1_r', 'Blue6',
    ratios=(1, 3, 5, 10), name='SciVisColor', save=True
)

# Plot examples
for ax, cmap, title in zip(axs, (cmap1, cmap2, cmap3), (title1, title2, title3)):
    func = (ax.pcolormesh if cmap is cmap1 else ax.contourf)
    m = func(data, cmap=cmap, levels=256)
    ax.colorbar(m, loc='b', locator='null', label=cmap.name)
    ax.format(title=title)
axs.format(
    xlabel='xlabel', ylabel='ylabel',
    suptitle='Merging existing colormaps'
)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_cmaps_mod:
#
# Modifying colormaps
# -------------------
#
# ProPlot allows you to create modified versions of *existing* colormaps
# using the `~proplot.constructor.Colormap` constructor and the new
# `~proplot.colors.LinearSegmentedColormap` and
# `~proplot.colors.ListedColormap` classes, which are used to replace the
# native matplotlib colormap classes. They can be modified in the following
# ways:
#
# * To remove colors from the left or right ends of a colormap, pass `left`
#   or `right` to `~proplot.constructor.Colormap`. This calls the
#   `~proplot.colors.LinearSegmentedColormap.truncate` method, and can be
#   useful when you want to use colormaps as :ref:`color cycles <ug_cycles>`
#   and need to remove the "white" part so that your lines stand out against
#   the background.
# * To remove central colors from a diverging colormap, pass `cut` to
#   `~proplot.constructor.Colormap`. This calls the
#   `~proplot.colors.LinearSegmentedColormap.cut` method, and can be used
#   to create a sharper cutoff between negative and positive values. This
#   should generally be used *without* a central level.
# * To rotate a cyclic colormap,  pass `shift` to
#   `~proplot.constructor.Colormap`. This calls the
#   `~proplot.colors.LinearSegmentedColormap.shifted` method. ProPlot ensures
#   the colors at the ends of "shifted" colormaps are *distinct* so that
#   levels never blur together.
# * To change the opacity of a colormap or add an opacity *gradation*, pass
#   `alpha` to `~proplot.constructor.Colormap`. This calls the
#   `~proplot.colors.LinearSegmentedColormap.set_alpha` method, and can be
#   useful when *layering* filled contour or mesh elements.
# * To change the "gamma" of a `~proplot.colors.PerceptuallyUniformColormap`,
#   pass `gamma` to `~proplot.constructor.Colormap`. This calls the
#   `~proplot.colors.PerceptuallyUniformColormap.set_gamma` method, and
#   controls how the luminance and saturation channels vary between colormap
#   segments. ``gamma > 1`` emphasizes high luminance, low saturation colors,
#   while ``gamma < 1`` emphasizes low luminance, high saturation colors.

# %%
import proplot as plot
import numpy as np
fig, axs = plot.subplots(
    [[1, 1, 2, 2, 3, 3], [0, 4, 4, 5, 5, 0], [0, 6, 6, 7, 7, 0]],
    axwidth=1.7, span=False
)
state = np.random.RandomState(51423)
data = state.rand(40, 40).cumsum(axis=0) - 12

# Cutting left and right
for ax, coord in zip(axs[:3], (None, 0.3, 0.7)):
    cmap = 'grays'
    if coord is None:
        title, cmap_kw = 'Original', {}
    elif coord < 0.5:
        title, cmap_kw = f'left={coord}', {'left': coord}
    else:
        title, cmap_kw = f'right={coord}', {'right': coord}
    ax.pcolormesh(
        data, cmap=cmap, cmap_kw=cmap_kw,
        colorbar='b', colorbar_kw={'locator': 'null'}
    )
    ax.format(xlabel='x axis', ylabel='y axis', title=title)

# Cutting central colors
levels = plot.arange(-10, 10, 2)
for i, (ax, cut) in enumerate(zip(axs[3:], (None, None, 0.1, 0.2))):
    if i == 0:
        title = 'With central level'
        levels = plot.edges(plot.arange(-10, 10, 2))
    else:
        title = 'Without central level'
        levels = plot.arange(-10, 10, 2)
    if cut is not None:
        title = f'cut = {cut}'
    m = ax.contourf(
        data, cmap='Div', cmap_kw={'cut': cut},
        extend='both', levels=levels,
    )
    ax.format(
        xlabel='x axis', ylabel='y axis', title=title,
        suptitle='Truncating sequential and diverging colormaps'
    )
    ax.colorbar(m, loc='b', locator='null')

# %%
import proplot as plot
import numpy as np
state = np.random.RandomState(51423)

# Rotating cyclic colormaps
fig, axs = plot.subplots(ncols=3, axwidth=1.7)
data = (state.rand(50, 50) - 0.48).cumsum(axis=1).cumsum(axis=0) - 50
for ax, shift in zip(axs, (0, 90, 180)):
    m = ax.contourf(data, cmap='twilight', cmap_kw={'shift': shift}, levels=12)
    ax.format(
        xlabel='x axis', ylabel='y axis', title=f'shift = {shift}',
        suptitle='Rotating cyclic colormaps'
    )
    ax.colorbar(m, loc='b', locator='null')

# Changing the colormap opacity
fig, axs = plot.subplots(ncols=3, axwidth=1.7)
data = state.rand(10, 10).cumsum(axis=1)
for ax, alpha in zip(axs, (1.0, 0.5, 0.0)):
    alpha = (alpha, 1 - 0.2 * (1 - alpha), 1.0)
    cmap = plot.Colormap('lajolla', alpha=alpha)
    m = ax.contourf(data, cmap=cmap, levels=10, extend='both', linewidth=0)
    ax.colorbar(m, loc='b', locator='none')
    ax.format(
        title=f'alpha = {alpha}', xlabel='x axis', ylabel='y axis',
        suptitle='Adding opacity gradations'
    )

# Changing the colormap gamma
fig, axs = plot.subplots(ncols=3, axwidth=1.7)
data = state.rand(10, 10).cumsum(axis=1)
for ax, gamma in zip(axs, (0.7, 1.0, 1.4)):
    cmap = plot.Colormap('boreal', gamma=gamma)
    m = ax.pcolormesh(data, cmap=cmap, levels=10, extend='both')
    ax.colorbar(m, loc='b', locator='none')
    ax.format(
        title=f'gamma = {gamma}', xlabel='x axis', ylabel='y axis',
        suptitle='Changing the PerceptuallyUniformColormap gamma'
    )


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_cmaps_dl:
#
# Downloading colormaps
# ---------------------
#
# There are plenty of online interactive tools for generating perceptually
# uniform colormaps, including
# `Chroma.js <https://gka.github.io/palettes/>`__,
# `HCLWizard <http://hclwizard.org:64230/hclwizard/>`__,
# `HCL picker <http://tristen.ca/hcl-picker/>`__,
# the `CCC-tool <https://ccctool.com>`__,
# and `SciVisColor <https://sciviscolor.org/home/colormaps/>`__.
#
# To add colormaps downloaded from any of these sources, save the colormap
# data to a file in your ``~/.proplot/cmaps`` folder and call
# `~proplot.config.register_cmaps` (or restart your python session), or use
# `~proplot.colors.LinearSegmentedColormap.from_file`. The file name is used
# as the registered colormap name. See
# `~proplot.colors.LinearSegmentedColormap.from_file` for a table of valid
# file extensions.
