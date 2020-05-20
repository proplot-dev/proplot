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
# .. _ug_2dplots:
#
# Plotting 2D data
# ================
#
# ProPlot adds new features to various `~matplotlib.axes.Axes` plotting
# methods using a set of wrapper functions. When a plotting method like
# `~matplotlib.axes.Axes.contourf` is "wrapped" by one of these functions, it
# accepts the same parameters as the wrapper. These features are a strict
# *superset* of the matplotlib API.
# This section documents the features added by wrapper functions to 2D
# plotting commands like `~matplotlib.axes.Axes.contour`,
# `~matplotlib.axes.Axes.contourf`, `~matplotlib.axes.Axes.pcolor`, and
# `~matplotlib.axes.Axes.pcolormesh`.


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_cmap_changer:
#
# Colormaps and normalizers
# -------------------------
#
# It is often desirable to create ProPlot colormaps on-the-fly, without
# explicitly using the `~proplot.constructor.Colormap` constructor function.
# To enable this, the `~proplot.axes.cmap_changer` wrapper adds the
# `cmap` and `cmap_kw` arguments to every 2D plotting method. These
# arguments are passed to the `~proplot.constructor.Colormap` constructor
# function, and the resulting colormap is used for the input data. For
# example, to create and apply a monochromatic colormap, you can simply use
# ``cmap='color name'``.

# The `~proplot.axes.cmap_changer` wrapper also
# adds the `norm` and `norm_kw` arguments. They are passed to the
# `~proplot.constructor.Norm` constructor function, and the resulting
# normalizer is used for the input data. For more information on colormaps
# and normalizers, see the :ref:`colormaps section <ug_cmaps>` and `this
# matplotlib tutorial
# <https://matplotlib.org/tutorials/colors/colormapnorms.html>`__.

# %%
import proplot as plot
import numpy as np
N = 20
cmap = plot.Colormap(('orange0', 'blood'))
state = np.random.RandomState(51423)
data = 10 ** (0.25 * np.cumsum(state.rand(N, N), axis=0))
with plot.rc.context({'lines.linewidth': 3}):
    fig, axs = plot.subplots(ncols=2, span=False)
    axs.format(
        xlabel='xlabel', ylabel='ylabel',
        suptitle='On-the-fly colormaps and normalizers'
    )

    # On-the-fly colormaps and normalizers
    axs[0].pcolormesh(data, cmap=cmap, colorbar='b')
    axs[1].pcolormesh(data, norm='log', cmap=cmap, colorbar='b')
    axs[0].format(title='Linear normalizer')
    axs[1].format(title='Logarithmic normalizer')


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_discrete:
#
# Discrete colormap levels
# ------------------------
#
# The `~proplot.axes.cmap_changer` wrapper also applies the
# `~proplot.colors.DiscreteNorm` normalizer to every colormap plot.
# `~proplot.colors.DiscreteNorm` converts data values to colormap colors by (1)
# transforming data using an arbitrary *continuous* normalizer (e.g.
# `~matplotlib.colors.LogNorm`), then (2) mapping the normalized data to
# *discrete* colormap levels (just like `~matplotlib.colors.BoundaryNorm`).
#
# By applying `~proplot.colors.DiscreteNorm` to every plot, ProPlot permits
# distinct "levels" even for commands like `~matplotlib.axes.Axes.pcolor` and
# `~matplotlib.axes.Axes.pcolormesh`. Distinct levels can help the reader
# discern exact numeric values and tends to reveal qualitative structure in
# the figure. They are also critical for users that would *prefer* contours,
# but have complex 2D coordinate matrices that trip up the contouring
# algorithm.  `~proplot.colors.DiscreteNorm` also fixes the colormap
# end-colors by ensuring the following conditions are met (this may seem
# nitpicky, but it is crucial for plots with very few levels):
#
# #. All colormaps always span the *entire color range*, independent of the
#    `extend` setting.
# #. Cyclic colormaps always have *distinct color levels* on
#    either end of the colorbar.

# %%
import proplot as plot
import numpy as np

# Pcolor plot with and without distinct levels
fig, axs = plot.subplots(ncols=2, axwidth=2)
state = np.random.RandomState(51423)
data = (state.normal(0, 1, size=(33, 33))).cumsum(axis=0).cumsum(axis=1)
axs.format(suptitle='Pcolor plot with levels')
for ax, n, mode, side in zip(axs, (200, 10), ('Ambiguous', 'Discernible'), 'lr'):
    ax.pcolor(data, cmap='spectral_r', N=n, symmetric=True, colorbar=side)
    ax.format(title=f'{mode} level boundaries', yformatter='null')

# %%
import proplot as plot
import numpy as np
fig, axs = plot.subplots(
    [[0, 0, 1, 1, 0, 0], [2, 3, 3, 4, 4, 5]],
    wratios=(1.5, 0.5, 1, 1, 0.5, 1.5), axwidth=1.7, ref=1, right='2em'
)
axs.format(suptitle='DiscreteNorm color-range standardization')
levels = plot.arange(0, 360, 45)
state = np.random.RandomState(51423)
data = (20 * (state.rand(20, 20) - 0.4).cumsum(axis=0).cumsum(axis=1)) % 360

# Cyclic colorbar with distinct end colors
ax = axs[0]
ax.pcolormesh(
    data, levels=levels, cmap='phase', extend='neither',
    colorbar='b', colorbar_kw={'locator': 90}
)
ax.format(title='cyclic colormap\nwith distinct end colors')

# Colorbars with different extend values
for ax, extend in zip(axs[1:], ('min', 'max', 'neither', 'both')):
    ax.pcolormesh(
        data[:, :10], levels=levels, cmap='oxy',
        extend=extend, colorbar='b', colorbar_kw={'locator': 90}
    )
    ax.format(title=f'extend={extend!r}')


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_norm:
#
# Special normalizers
# -------------------
#
# The `~proplot.colors.LinearSegmentedNorm` colormap normalizer
# provides even color gradations with respect to *index* for an
# arbitrary monotonically increasing list of levels. This is automatically applied
# if you pass unevenly spaced `levels` to a plotting command, or it can be manually
# applied using e.g. ``norm='segmented'``.
#
# The `~proplot.colors.DivergingNorm` normalizer
# ensures the colormap midpoint lies on some *central* data value (usually 0),
# even if `vmin`, `vmax`, or `levels` are asymmetric with respect to the central
# value. This can be applied using e.g. ``norm='diverging'`` and be configured
# to scale colors "fairly" or "unfairly":
#
# * With fair scaling (the default), the gradations on either side of the midpoint
#   have equal intensity. If `vmin` and `vmax` are not symmetric about zero, the most
#   intense colormap colors on one side of the midpoint will be truncated.
# * With unfair scaling, the gradations on either side of the midpoint are warped
#   so that the full range of colormap colors is traversed. This configuration should
#   be used with care, as it may lead you to misinterpret your data!
#
# The below example demonstrates how these normalizers can be used for datasets
# with unusual statistical distributions.

# %%
import proplot as plot
import numpy as np

# Linear segmented norm
state = np.random.RandomState(51423)
data = 10**(2 * state.rand(20, 20).cumsum(axis=0) / 7)
fig, axs = plot.subplots(ncols=2, axwidth=2.4)
ticks = [5, 10, 20, 50, 100, 200, 500, 1000]
for i, (norm, title) in enumerate(zip(
    ('linear', 'segmented'),
    ('Linear normalizer', 'LinearSegmentedNorm')
)):
    m = axs[i].contourf(
        data, levels=ticks, extend='both',
        cmap='Mako', norm=norm,
        colorbar='b', colorbar_kw={'ticks': ticks},
    )
    axs[i].format(title=title)
axs.format(suptitle='Linear segmented normalizer demo')

# %%
import proplot as plot
import numpy as np

# Diverging norm
data1 = (state.rand(20, 20) - 0.43).cumsum(axis=0)
data2 = (state.rand(20, 20) - 0.57).cumsum(axis=0)
fig, axs = plot.subplots(nrows=2, ncols=2, axwidth=2.4, order='F')
cmap = plot.Colormap('DryWet', cut=0.1)
axs.format(suptitle='Diverging normalizer demo')
i = 0
for data, mode, fair in zip(
    (data1, data2),
    ('positive', 'negative'),
    ('fair', 'unfair')
):
    for fair in ('fair', 'unfair'):
        norm = plot.Norm('diverging', fair=(fair == 'fair'))
        ax = axs[i]
        m = ax.contourf(data, cmap=cmap, norm=norm)
        ax.colorbar(m, loc='b', locator=1)
        ax.format(title=f'Skewed {mode} data, {fair!r} scaling')
        i += 1


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_2dstd:
#
# Standardized arguments
# ----------------------
#
# The `~proplot.axes.standardize_2d` wrapper is used to standardize
# positional arguments across all 2D plotting methods.  Among other things,
# it guesses coordinate *edges* for `~matplotlib.axes.Axes.pcolor` and
# `~matplotlib.axes.Axes.pcolormesh` plots when you supply coordinate
# *centers*, and calculates coordinate *centers* for
# `~matplotlib.axes.Axes.contourf` and `~matplotlib.axes.Axes.contour` plots
# when you supply coordinate *edges*. Notice the locations of the rectangle
# edges in the ``pcolor`` plots shown below.

# %%
import proplot as plot
import numpy as np

# Figure and sample data
state = np.random.RandomState(51423)
x = y = np.array([-10, -5, 0, 5, 10])
xedges = plot.edges(x)
yedges = plot.edges(y)
data = state.rand(y.size, x.size)  # "center" coordinates
lim = (np.min(xedges), np.max(xedges))
with plot.rc.context({'image.cmap': 'Grays', 'image.levels': 21}):
    fig, axs = plot.subplots(ncols=2, nrows=2, share=False)
    axs.format(
        xlabel='xlabel', ylabel='ylabel',
        xlim=lim, ylim=lim, xlocator=5, ylocator=5,
        suptitle='Standardized input demonstration'
    )
    axs[0].format(title='Supplying coordinate centers')
    axs[1].format(title='Supplying coordinate edges')

    # Plot using both centers and edges as coordinates
    axs[0].pcolormesh(x, y, data)
    axs[1].pcolormesh(xedges, yedges, data)
    axs[2].contourf(x, y, data)
    axs[3].contourf(xedges, yedges, data)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_2dintegration:
#
# Pandas and xarray integration
# -----------------------------
#
# The `~proplot.axes.standardize_2d` wrapper also integrates 2D
# plotting methods with pandas `~pandas.DataFrame`\ s and xarray
# `~xarray.DataArray`\ s. When you pass a DataFrame or DataArray to any
# plotting command, the x-axis label, y-axis label, legend label, colorbar
# label, and/or title are configured from the metadata. This restores some of
# the convenience you get with the builtin `pandas
# <https://pandas.pydata.org>`__ and `xarray <https://pandas.pydata.org>`__
# plotting functions. This feature is *optional*; installation of pandas and
# xarray are not required.

# %%
import xarray as xr
import numpy as np
import pandas as pd

# DataArray
state = np.random.RandomState(51423)
linspace = np.linspace(0, np.pi, 20)
data = 50 * state.normal(1, 0.2, size=(20, 20)) * (
    np.sin(linspace * 2) ** 2
    * np.cos(linspace + np.pi / 2)[:, None] ** 2
)
lat = xr.DataArray(
    np.linspace(-90, 90, 20),
    dims=('lat',),
    attrs={'units': 'deg_north'}
)
plev = xr.DataArray(
    np.linspace(1000, 0, 20),
    dims=('plev',),
    attrs={'long_name': 'pressure', 'units': 'mb'}
)
da = xr.DataArray(
    data,
    name='u',
    dims=('plev', 'lat'),
    coords={'plev': plev, 'lat': lat},
    attrs={'long_name': 'zonal wind', 'units': 'm/s'}
)

# DataFrame
data = state.rand(12, 20)
df = pd.DataFrame(
    (data - 0.4).cumsum(axis=0).cumsum(axis=1),
    index=list('JFMAMJJASOND'),
)
df.name = 'temporal data'
df.index.name = 'month'
df.columns.name = 'variable (units)'

# %%
import proplot as plot
fig, axs = plot.subplots(nrows=2, axwidth=2.5, share=0)
axs.format(collabels=['Automatic subplot formatting'])

# Plot DataArray
cmap = plot.Colormap('RdPu', left=0.05)
axs[0].contourf(da, cmap=cmap, colorbar='l', linewidth=0.7, color='k')
axs[0].format(yreverse=True)

# Plot DataFrame
axs[1].contourf(df, cmap='YlOrRd', colorbar='r', linewidth=0.7, color='k')
axs[1].format(xtickminor=False)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_labels:
#
# Contour and gridbox labels
# --------------------------
#
# The `~proplot.axes.cmap_changer` wrapper also allows you to quickly add
# *labels* to `~proplot.axes.Axes.heatmap`, `~matplotlib.axes.Axes.pcolor`,
# `~matplotlib.axes.Axes.pcolormesh`, `~matplotlib.axes.Axes.contour`, and
# `~matplotlib.axes.Axes.contourf` plots by simply using ``labels=True``.
# The label text is colored black or white depending on the luminance of
# the underlying grid box or filled contour.
#
# `~proplot.axes.cmap_changer` draws contour labels with
# `~matplotlib.axes.Axes.clabel` and grid box labels with
# `~matplotlib.axes.Axes.text`. You can pass keyword arguments to these
# functions using the `labels_kw` dictionary keyword argument, and change the
# label precision with the `precision` keyword argument. See
# `~proplot.axes.cmap_changer` for details.

# %%
import proplot as plot
import pandas as pd
import numpy as np
fig, axs = plot.subplots(
    [[1, 1, 2, 2], [0, 3, 3, 0]],
    axwidth=2.2, share=1, span=False, hratios=(1, 0.9)
)
state = np.random.RandomState(51423)
data = state.rand(6, 6)
data = pd.DataFrame(data, index=pd.Index(['a', 'b', 'c', 'd', 'e', 'f']))
axs.format(xlabel='xlabel', ylabel='ylabel', suptitle='Labels demo')

# Heatmap with labeled boxes
ax = axs[0]
m = ax.heatmap(
    data, cmap='rocket', labels=True,
    precision=2, labels_kw={'weight': 'bold'}
)
ax.format(title='Heatmap plot with labels')

# Filled contours with labels
ax = axs[1]
m = ax.contourf(
    data.cumsum(axis=0), labels=True,
    cmap='rocket', labels_kw={'weight': 'bold'}
)
ax.format(title='Filled contour plot with labels')

# Line contours with labels
ax = axs[2]
ax.contour(
    data.cumsum(axis=1) - 2, color='gray8',
    labels=True, lw=2, labels_kw={'weight': 'bold'}
)
ax.format(title='Line contour plot with labels')


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_heatmap:
#
# Heatmap plots
# -------------
#
# The new `~proplot.axes.Axes.heatmap` command calls
# `~matplotlib.axes.Axes.pcolormesh` and configures the axes with settings
# that are suitable for heatmaps -- fixed aspect ratio, no gridlines, no minor ticks,
# and major ticks at the center of each box. Among other things, this is useful for
# displaying covariance and correlation matrices, as shown below.

# %%
import proplot as plot
import numpy as np
import pandas as pd

# Covariance data
state = np.random.RandomState(51423)
data = state.normal(size=(10, 10)).cumsum(axis=0)
data = (data - data.mean(axis=0)) / data.std(axis=0)
data = (data.T @ data) / data.shape[0]
data[np.tril_indices(data.shape[0], -1)] = np.nan  # fill half with empty boxes
data = pd.DataFrame(data, columns=list('abcdefghij'), index=list('abcdefghij'))

# Covariance matrix plot
fig, ax = plot.subplots(axwidth=4.5)
m = ax.heatmap(
    data, cmap='ColdHot', vmin=-1, vmax=1, N=100,
    lw=0.5, edgecolor='k', labels=True, labels_kw={'weight': 'bold'},
    clip_on=False,  # turn off clipping so box edges are not cut in half
)
ax.format(
    suptitle='Heatmap demo', title='Table of correlation coefficients', alpha=0,
    xloc='top', yloc='right', yreverse=True, ticklabelweight='bold', linewidth=0,
    ytickmajorpad=4,  # the ytick.major.pad rc setting; adds extra space
)
