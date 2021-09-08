# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _cmocean: https://matplotlib.org/cmocean/
#
# .. _fabio: http://www.fabiocrameri.ch/colourmaps.php
#
# .. _brewer: http://colorbrewer2.org/
#
# .. _sciviscolor: https://sciviscolor.org/home/colormoves/
#
# .. _matplotlib: https://matplotlib.org/stable/tutorials/colors/colormaps.html
#
# .. _seaborn: https://seaborn.pydata.org/tutorial/color_palettes.html
#
# .. _ug_cmaps:
#
# Colormaps
# =========
#
# ProPlot defines **continuous colormaps** as color palettes that sample some
# *continuous function* between two end colors. They are generally used
# to encode data values on a pseudo-third dimension. They are are implemented
# with the `~proplot.colors.ContinuousColormap` and
# `~proplot.colors.PerceptualColormap` classes, which are
# :ref:`subclassed from <ug_cmaps_new>`
# `matplotlib.colors.LinearSegmentedColormap`.
#
# ProPlot :ref:`adds several features <why_colormaps_cycles>` to help you use
# colormaps effectively in your figures. This section documents the new registered
# colormaps, explains how to make and modify colormaps, and shows how to apply them
# to your plots.


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
# * "User" colormaps created with `~proplot.constructor.Colormap`
#   or loaded from `~proplot.config.Configurator.user_folder`.
# * `Matplotlib <matplotlib_>`_ and `seaborn <seaborn_>`_ original colormaps.
# * ProPlot original :ref:`perceptually uniform colormaps <ug_perceptual>`.
# * The `cmOcean <cmocean_>`_ colormaps, designed for
#   oceanographic data but useful for everyone.
# * Fabio Crameri's `"scientific colour maps" <fabio_>`_.
# * Cynthia Brewer's `ColorBrewer <brewer_>`_ colormaps,
#   included with matplotlib by default.
# * Colormaps from the `SciVisColor <sciviscolor_>`_ project. There are so many
#   of these because they are intended to be merged into more complex colormaps.
#
# Matplotlib colormaps with erratic color transitions like ``'jet'`` are still
# registered, but they are hidden from this table, and their usage is discouraged.
#
# .. note::
#
#    Colormap and :ref:`color cycle <ug_cycles>` identification is more flexible in
#    ProPlot. The names are are case-insensitive (e.g., ``'Viridis'``, ``'viridis'``,
#    and ``'ViRiDiS'`` are equivalent), diverging colormap names can be specified in
#    their "reversed" form (e.g., ``'BuRd'`` is equivalent to ``'RdBu_r'``), and
#    appending ``'_r'`` or ``'_s'`` to *any* colormap name will return a
#    `~proplot.colors.ContinuousColormap.reversed` or
#    `~proplot.colors.ContinuousColormap.shifted` version of the colormap
#    or color cycle. See `~proplot.colors.ColormapDatabase` for more info.

# %%
import proplot as pplt
fig, axs = pplt.show_cmaps()


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_perceptual:
#
# Perceptually uniform colormaps
# ------------------------------
#
# ProPlot's custom colormaps are instances of the
# `~proplot.colors.PerceptualColormap` class. These colormaps
# generate colors by interpolating between coordinates in any
# of the following three hue-saturation-luminance colorspaces:
#
# * **HCL** (a.k.a. `CIE LChuv <https://en.wikipedia.org/wiki/CIELUV>`__):
#   A purely perceptually uniform colorspace, where colors are broken down
#   into “hue” (color, range 0-360), “chroma” (saturation, range 0-100), and
#   “luminance” (brightness, range 0-100). This colorspace is difficult to work
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
# The colorspace used by a `~proplot.colors.PerceptualColormap`
# is set with the `space` keyword arg. To plot arbitrary cross-sections of
# these colorspaces, use `~proplot.demos.show_colorspaces` (the black
# regions represent impossible colors). To see how colormaps vary with
# respect to each channel, use `~proplot.demos.show_channels`. Some examples
# are shown below.
#
# In theory, "uniform" colormaps should have *straight* lines in hue, chroma,
# and luminance (second figure, top row). In practice, this is
# difficult to accomplish due to impossible colors. Matplotlib's and seaborn's
# ``'magma'`` and ``'Rocket'`` colormaps are fairly linear with respect to
# hue and luminance, but not chroma. ProPlot's ``'Fire'`` is linear in hue,
# luminance, and *HSL saturation* (bottom left), while ``'Dusk'`` is linear
# in hue, luminance, and *HPL saturation* (bottom right).

# %%
# Colorspace demo
import proplot as pplt
fig, axs = pplt.show_colorspaces(refwidth=1.6, luminance=50)
fig, axs = pplt.show_colorspaces(refwidth=1.6, saturation=60)
fig, axs = pplt.show_colorspaces(refwidth=1.6, hue=0)

# %%
# Compare colormaps
import proplot as pplt
for cmaps in (('magma', 'rocket'), ('fire', 'dusk')):
    fig, axs = pplt.show_channels(
        *cmaps, refwidth=1.5, minhue=-180, maxsat=400, rgb=False
    )


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_cmaps_new:
#
# Making colormaps
# ----------------
#
# ProPlot doesn't just include new colormaps -- it provides tools
# for merging colormaps, modifying existing colormaps, making new
# :ref:`perceptually uniform colormaps <ug_perceptual>`, and saving colormaps
# for future use. Most of these features can be accessed via the
# `~proplot.constructor.Colormap` :ref:`constructor function <why_constructor>`.
# Note that every plotting command that accepts a `cmap` keyword passes
# it through this function (see the :ref:`2D plotting section <ug_apply_cmap>`).
#
# To make `~proplot.colors.PerceptualColormap`\ s from
# scratch, you have the following three options:
#
# * Pass a color name, HEX string, or RGB tuple to `~proplot.constructor.Colormap`.
#   This builds a monochromatic (single hue) colormap by calling
#   `~proplot.colors.PerceptualColormap.from_color`. The colormap colors will
#   progress from the specified color to a color with the same hue but changed
#   saturation or luminance. These can be set with the `saturation` and `luminance`
#   keywords (or their shorthands `s` and `l`). By default, the colormap will
#   progress to pure white.
# * Pass a list of color names, HEX strings, or RGB
#   tuples to `~proplot.constructor.Colormap`. This calls
#   `~proplot.colors.PerceptualColormap.from_list`, which linearly interpolates
#   between the hues, saturations, and luminances of the input colors. To facillitate
#   the construction of diverging colormaps, the hue channel values for nuetral
#   colors (i.e., white, black, and gray) are adjusted to the hues of the preceding
#   and subsequent colors in the list, with sharp hue cutoffs at the neutral colors.
#   This permits generating diverging colormaps with e.g. ``['blue', 'white', 'red']``.
# * Pass the keywords `hue`, `saturation`, or `luminance` (or their shorthands `h`,
#   `s`, and `l`) to `~proplot.constructor.Colormap` without any positional arguments
#   (or pass a dictionary containing these keys as a positional argument).
#   This calls `~proplot.colors.PerceptualColormap.from_hsl`, which
#   linearly interpolates between the specified channel values. Channel values can be
#   specified with numbers between ``0`` and ``100``, color strings, or lists thereof.
#   For color strings, the value is *inferred* from the specified color. You can
#   end any color string with ``'+N'`` or ``'-N'`` to *offset* the channel
#   value by the number ``N`` (e.g., ``hue='red+50'``).
#
# To change the :ref:`colorspace <ug_perceptual>` used to construct the colormap,
# use the `space` keyword. The default colorspace is ``'hsl'``. In the below example,
# we use all of these methods to make `~proplot.colors.PerceptualColormap`\ s
# in the ``'hsl'`` and ``'hpl'`` colorspaces.

# %%
# Sample data
import proplot as pplt
import numpy as np
state = np.random.RandomState(51423)
data = state.rand(30, 30).cumsum(axis=1)

# %%
# Colormap from a color
# The trailing '_r' makes the colormap go dark-to-light instead of light-to-dark
fig = pplt.figure(refwidth=2, span=False)
ax = fig.subplot(121)
ax.format(title='From single named color')
cmap1 = pplt.Colormap('prussian blue_r', l=100, name='Pacific', space='hpl')
m = ax.contourf(data, cmap=cmap1)
ax.colorbar(m, loc='b', ticks='none', label=cmap1.name)

# Colormap from lists
ax = fig.subplot(122)
ax.format(title='From list of colors')
cmap2 = pplt.Colormap(('maroon', 'light tan'), name='Heatwave')
m = ax.contourf(data, cmap=cmap2)
ax.colorbar(m, loc='b', ticks='none', label=cmap2.name)
fig.format(
    xticklabels='none',
    yticklabels='none',
    suptitle='Making PerceptualColormaps'
)

# Display the channels
fig, axs = pplt.show_channels(cmap1, cmap2, refwidth=1.5, rgb=False)

# %%
# Sequential colormap from channel values
cmap3 = pplt.Colormap(
    h=('red', 'red-720'), s=(80, 20), l=(20, 100), space='hpl', name='CubeHelix'
)
fig = pplt.figure(refwidth=2, span=False)
ax = fig.subplot(121)
ax.format(title='Sequential from channel values')
m = ax.contourf(data, cmap=cmap3)
ax.colorbar(m, loc='b', ticks='none', label=cmap3.name)

# Cyclic colormap from channel values
ax = fig.subplot(122)
ax.format(title='Cyclic from channel values')
cmap4 = pplt.Colormap(
    h=(0, 360), c=50, l=70, space='hcl', cyclic=True, name='Spectrum'
)
m = ax.contourf(data, cmap=cmap4)
ax.colorbar(m, loc='b', ticks='none', label=cmap4.name)
fig.format(
    xticklabels='none',
    yticklabels='none',
    suptitle='Making PerceptualColormaps'
)

# Display the channels
fig, axs = pplt.show_channels(cmap3, cmap4, refwidth=1.5, rgb=False)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_cmaps_merge:
#
# Merging colormaps
# -----------------
#
# To *merge* colormaps, you can pass multiple positional arguments to the
# `~proplot.constructor.Colormap` constructor function. This calls the
# `~proplot.colors.ContinuousColormap.append` method. Each positional
# argument can be a colormap name, a colormap instance, or a
# :ref:`special argument <ug_cmaps_new>` that generates a new colormap
# on-the-fly. This lets you create new diverging colormaps and segmented
# `SciVisColor <https://sciviscolor.org/home/colormoves/>`__ style colormaps
# right inside ProPlot. Segmented colormaps are often desirable for complex
# datasets with complex statistical distributions.
#
# In the below example, we create a new divering colormap and
# reconstruct the colormap from `this SciVisColor example
# <https://sciviscolor.org/media/filer_public/c7/27/c727e638-82eb-445b-bc96-e7b64c13efa2/colormoves.png>`__.
# We also *save* the results for future use by passing ``save=True`` to
# `~proplot.constructor.Colormap`.

# %%
import proplot as pplt
import numpy as np
state = np.random.RandomState(51423)
data = state.rand(30, 30).cumsum(axis=1)

# Generate figure
fig, axs = pplt.subplots([[0, 1, 1, 0], [2, 2, 3, 3]], refwidth=2.4, span=False)
axs.format(
    xlabel='xlabel', ylabel='ylabel',
    suptitle='Merging colormaps'
)

# Diverging colormap example
title1 = 'Diverging from sequential maps'
cmap1 = pplt.Colormap('Blues4_r', 'Reds3', name='Diverging', save=True)

# SciVisColor examples
title2 = 'SciVisColor example'
cmap2 = pplt.Colormap(
    'Greens1_r', 'Oranges1', 'Blues1_r', 'Blues6',
    ratios=(1, 3, 5, 10), name='SciVisColorUneven', save=True
)
title3 = 'SciVisColor with equal ratios'
cmap3 = pplt.Colormap(
    'Greens1_r', 'Oranges1', 'Blues1_r', 'Blues6',
    name='SciVisColorEven', save=True
)

# Plot examples
for ax, cmap, title in zip(axs, (cmap1, cmap2, cmap3), (title1, title2, title3)):
    m = ax.contourf(data, cmap=cmap, levels=500)
    ax.colorbar(m, loc='b', locator='null', label=cmap.name)
    ax.format(title=title)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_cmaps_mod:
#
# Modifying colormaps
# -------------------
#
# ProPlot lets you create modified versions of *existing* colormaps
# using the `~proplot.constructor.Colormap` constructor function and the
# new `~proplot.colors.ContinuousColormap` and
# `~proplot.colors.DiscreteColormap` classes, which replace the native
# matplotlib colormap classes. They can be modified in the following ways:
#
# * To remove colors from the left or right ends of a colormap, pass `left`
#   or `right` to `~proplot.constructor.Colormap`. This calls the
#   `~proplot.colors.ContinuousColormap.truncate` method, and can be
#   useful when you want to use colormaps as :ref:`color cycles <ug_cycles>`
#   and need to remove the light part so that your lines stand out
#   against the background.
# * To modify the central colors of a diverging colormap, pass `cut` to
#   `~proplot.constructor.Colormap`. This calls the
#   `~proplot.colors.ContinuousColormap.cut` method, and can be used
#   to create a sharper cutoff between negative and positive values or (when
#   `cut` is negative) to expand the "neutral" region of the colormap.
# * To rotate a cyclic colormap,  pass `shift` to
#   `~proplot.constructor.Colormap`. This calls the
#   `~proplot.colors.ContinuousColormap.shifted` method. ProPlot ensures
#   the colors at the ends of "shifted" colormaps are *distinct* so that
#   levels never blur together.
# * To change the opacity of a colormap or add an opacity *gradation*, pass
#   `alpha` to `~proplot.constructor.Colormap`. This calls the
#   `~proplot.colors.ContinuousColormap.set_alpha` method, and can be
#   useful when *layering* filled contour or mesh elements.
# * To change the "gamma" of a `~proplot.colors.PerceptualColormap`,
#   pass `gamma` to `~proplot.constructor.Colormap`. This calls the
#   `~proplot.colors.PerceptualColormap.set_gamma` method, and
#   controls how the luminance and saturation channels vary between colormap
#   segments. ``gamma > 1`` emphasizes high luminance, low saturation colors,
#   while ``gamma < 1`` emphasizes low luminance, high saturation colors. This
#   is similar to the effect of the `HCL wizard
#   <http://hclwizard.org:64230/hclwizard/>`__ "power" sliders.

# %%
import proplot as pplt
import numpy as np
state = np.random.RandomState(51423)
data = state.rand(40, 40).cumsum(axis=0)

# Generate figure
fig, axs = pplt.subplots([[0, 1, 1, 0], [2, 2, 3, 3]], refwidth=1.9, span=False)
axs.format(xlabel='y axis', ylabel='x axis', suptitle='Truncating sequential colormaps')

# Cutting left and right
cmap = 'Ice'
for ax, coord in zip(axs, (None, 0.3, 0.7)):
    if coord is None:
        title, cmap_kw = 'Original', {}
    elif coord < 0.5:
        title, cmap_kw = f'left={coord}', {'left': coord}
    else:
        title, cmap_kw = f'right={coord}', {'right': coord}
    ax.format(title=title)
    ax.contourf(
        data, cmap=cmap, cmap_kw=cmap_kw, colorbar='b', colorbar_kw={'locator': 'null'}
    )

# %%
import proplot as pplt
import numpy as np
state = np.random.RandomState(51423)
data = (state.rand(40, 40) - 0.5).cumsum(axis=0).cumsum(axis=1)

# Create figure
fig, axs = pplt.subplots(ncols=2, nrows=2, refwidth=1.7, span=False)
axs.format(
    xlabel='x axis', ylabel='y axis', xticklabels='none',
    suptitle='Modifying diverging colormaps',
)

# Cutting out central colors
titles = (
    'Negative-positive cutoff', 'Neutral-valued center',
    'Sharper cutoff', 'Expanded center'
)
for i, (ax, title, cut) in enumerate(zip(axs, titles, (None, None, 0.2, -0.1))):
    if i % 2 == 0:
        kw = {'levels': pplt.arange(-10, 10, 2)}  # negative-positive cutoff
    else:
        kw = {'values': pplt.arange(-10, 10, 2)}  # dedicated center
    if cut is not None:
        fmt = pplt.SimpleFormatter()  # a proper minus sign
        title = f'{title}\ncut = {fmt(cut)}'
    ax.format(title=title)
    m = ax.contourf(
        data, cmap='Div', cmap_kw={'cut': cut}, extend='both',
        colorbar='b', colorbar_kw={'locator': 'null'},
        **kw  # level edges or centers
    )

# %%
import proplot as pplt
import numpy as np
state = np.random.RandomState(51423)
data = (state.rand(50, 50) - 0.48).cumsum(axis=0).cumsum(axis=1) % 30

# Rotating cyclic colormaps
fig, axs = pplt.subplots(ncols=3, refwidth=1.7, span=False)
for ax, shift in zip(axs, (0, 90, 180)):
    m = ax.pcolormesh(data, cmap='romaO', cmap_kw={'shift': shift}, levels=12)
    ax.format(
        xlabel='x axis', ylabel='y axis', title=f'shift = {shift}',
        suptitle='Rotating cyclic colormaps'
    )
    ax.colorbar(m, loc='b', locator='null')

# %%
import proplot as pplt
import numpy as np
state = np.random.RandomState(51423)
data = state.rand(20, 20).cumsum(axis=1)

# Changing the colormap opacity
fig, axs = pplt.subplots(ncols=3, refwidth=1.7, span=False)
for ax, alpha in zip(axs, (1.0, 0.5, 0.0)):
    alpha = (alpha, 1.0)
    cmap = pplt.Colormap('batlow_r', alpha=alpha)
    m = ax.imshow(data, cmap=cmap, levels=10, extend='both')
    ax.colorbar(m, loc='b', locator='none')
    ax.format(
        title=f'alpha = {alpha}', xlabel='x axis', ylabel='y axis',
        suptitle='Adding opacity gradations'
    )

# %%
import proplot as pplt
import numpy as np
state = np.random.RandomState(51423)
data = state.rand(20, 20).cumsum(axis=1)

# Changing the colormap gamma
fig, axs = pplt.subplots(ncols=3, refwidth=1.7, span=False)
for ax, gamma in zip(axs, (0.7, 1.0, 1.4)):
    cmap = pplt.Colormap('boreal', gamma=gamma)
    m = ax.pcolormesh(data, cmap=cmap, levels=10, extend='both')
    ax.colorbar(m, loc='b', locator='none')
    ax.format(
        title=f'gamma = {gamma}', xlabel='x axis', ylabel='y axis',
        suptitle='Changing the PerceptualColormap gamma'
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
# To add colormaps downloaded from any of these sources, save the colormap data file
# to the ``cmaps`` subfolder inside `~proplot.config.Configurator.user_folder`
# and call `~proplot.config.register_cmaps` (or restart your python session). You
# can also use `~proplot.colors.ContinuousColormap.from_file` or manually pass
# continuous colormaps or file paths to `~proplot.config.register_cmaps`. See
# `~proplot.colors.ContinuousColormap.from_file` for a table of valid
# data file extensions.
