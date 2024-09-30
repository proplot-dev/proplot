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
# .. _pandas: https://pandas.pydata.org
#
# .. _xarray: http://xarray.pydata.org/en/stable/
#
# .. _seaborn: https://seaborn.pydata.org
#
# .. _ug_2dplots:
#
# 2D plotting commands
# ====================
#
# Proplot adds :ref:`several new features <why_plotting>` to matplotlib's
# plotting commands using the intermediate `~proplot.axes.PlotAxes` class.
# For the most part, these additions represent a *superset* of matplotlib -- if
# you are not interested, you can use the plotting commands just like you would
# in matplotlib. This section documents the features added for 2D plotting commands
# like `~proplot.axes.PlotAxes.contour`, `~proplot.axes.PlotAxes.pcolor`,
# and `~proplot.axes.PlotAxes.imshow`.
#
# .. important::
#
#    By default, proplot automatically adjusts the width of
#    `~proplot.axes.PlotAxes.contourf` and `~proplot.axes.PlotAxes.pcolor` edges
#    to eliminate the appearance of `"white lines" in saved vector graphic files
#    <https://github.com/jklymak/contourfIssues>`__. However, this can significantly
#    slow down the drawing time for large datasets. To disable this feature,
#    pass ``edgefix=False`` to the relevant `~proplot.axes.PlotAxes` command,
#    or set :rcraw:`edgefix` to ``False`` to disable globally.

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_2dstd:
#
# Data arguments
# --------------
#
# The treatment of data arguments passed to the 2D `~proplot.axes.PlotAxes`
# commands is standardized. For each command, you can optionally omit the *x*
# and *y* coordinates, in which case they are inferred from the data
# (see :ref:`xarray and pandas integration <ug_2dintegration>`). If coordinates
# are string labels, they are converted to indices and tick labels using
# `~proplot.ticker.IndexLocator` and `~proplot.ticker.IndexFormatter`.
# If coordinates are descending and the axis limits are unset, the axis
# direction is automatically reversed. If coordinate *centers* are passed to commands
# like `~proplot.axes.PlotAxes.pcolor` and `~proplot.axes.PlotAxes.pcolormesh`, they
# are automatically converted to edges using `~proplot.utils.edges` or
# `~proplot.utils.edges2d`, and if coordinate *edges* are passed to commands like
# `~proplot.axes.PlotAxes.contour` and `~proplot.axes.PlotAxes.contourf`, they are
# automatically converted to centers (notice the locations of the rectangle edges
# in the ``pcolor`` plots below). All positional arguments can also be specified
# as keyword arguments (see the documentation for each plotting command).
#
# .. note::
#
#    By default, when choosing the colormap :ref:`normalization
#    range <ug_apply_cmap>`, proplot ignores data outside the *x* or *y* axis
#    limits if they were previously fixed by `~matplotlib.axes.Axes.set_xlim` or
#    `~matplotlib.axes.Axes.set_ylim` (or, equivalently, by passing `xlim` or
#    `ylim` to `proplot.axes.CartesianAxes.format`). This can be useful if you
#    wish to restrict the view along the *x* or *y* axis within a large dataset.
#    To disable this feature, pass ``inbounds=False`` to the plotting command or
#    set :rcraw:`cmap.inbounds` to ``False`` (see also the :rcraw:`axes.inbounds`
#    setting and the :ref:`user guide <ug_1dstd>`).

# %%
import proplot as pplt
import numpy as np

# Sample data
state = np.random.RandomState(51423)
x = y = np.array([-10, -5, 0, 5, 10])
xedges = pplt.edges(x)
yedges = pplt.edges(y)
data = state.rand(y.size, x.size)  # "center" coordinates
lim = (np.min(xedges), np.max(xedges))

with pplt.rc.context({"cmap": "Grays", "cmap.levels": 21}):
    # Figure
    fig = pplt.figure(refwidth=2.3, share=False)
    axs = fig.subplots(ncols=2, nrows=2)
    axs.format(
        xlabel="xlabel",
        ylabel="ylabel",
        xlim=lim,
        ylim=lim,
        xlocator=5,
        ylocator=5,
        suptitle="Standardized input demonstration",
        toplabels=("Coordinate centers", "Coordinate edges"),
    )

    # Plot using both centers and edges as coordinates
    axs[0].pcolormesh(x, y, data)
    axs[1].pcolormesh(xedges, yedges, data)
    axs[2].contourf(x, y, data)
    axs[3].contourf(xedges, yedges, data)

# %%
import proplot as pplt
import numpy as np

# Sample data
cmap = "turku_r"
state = np.random.RandomState(51423)
N = 80
x = y = np.arange(N + 1)
data = 10 + (state.normal(0, 3, size=(N, N))).cumsum(axis=0).cumsum(axis=1)
xlim = ylim = (0, 25)

# Plot the data
fig, axs = pplt.subplots(
    [[0, 1, 1, 0], [2, 2, 3, 3]],
    wratios=(1.3, 1, 1, 1.3),
    span=False,
    refwidth=2.2,
)
axs[0].fill_between(
    xlim,
    *ylim,
    zorder=3,
    edgecolor="red",
    facecolor=pplt.set_alpha("red", 0.2),
)
for i, ax in enumerate(axs):
    inbounds = i == 1
    title = f"Restricted lims inbounds={inbounds}"
    title += " (default)" if inbounds else ""
    ax.format(
        xlim=(None if i == 0 else xlim),
        ylim=(None if i == 0 else ylim),
        title=("Default axis limits" if i == 0 else title),
    )
    ax.pcolor(x, y, data, cmap=cmap, inbounds=inbounds)
fig.format(
    xlabel="xlabel",
    ylabel="ylabel",
    suptitle="Default vmin/vmax restricted to in-bounds data",
)

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_2dintegration:
#
# Pandas and xarray integration
# -----------------------------
#
# The 2D `~proplot.axes.PlotAxes` commands recognize `pandas`_
# and `xarray`_ data structures. If you omit *x* and *y* coordinates,
# the commands try to infer them from the `pandas.DataFrame` or
# `xarray.DataArray`. If you did not explicitly set the *x* or *y* axis label
# or :ref:`legend or colorbar <ug_guides_loc>` label(s), the commands
# try to retrieve them from the `pandas.DataFrame` or `xarray.DataArray`.
# The commands also recognize `pint.Quantity` structures and apply
# unit string labels with formatting specified by :rc:`unitformat`.
#
# These features restore some of the convenience you get with the builtin
# `pandas`_ and `xarray`_ plotting functions. They are also *optional* --
# installation of pandas and xarray are not required to use proplot. The
# automatic labels can be disabled by setting :rcraw:`autoformat` to ``False``
# or by passing ``autoformat=False`` to any plotting command.
#
# .. note::
#
#    For every plotting command, you can pass a `~xarray.Dataset`, `~pandas.DataFrame`,
#    or `dict` to the `data` keyword with strings as data arguments instead of arrays
#    -- just like matplotlib. For example, ``ax.plot('y', data=dataset)`` and
#    ``ax.plot(y='y', data=dataset)`` are translated to ``ax.plot(dataset['y'])``.
#    This is the preferred input style for most `seaborn`_ plotting commands.
#    Also, if you pass a `pint.Quantity` or `~xarray.DataArray`
#    containing a `pint.Quantity`, proplot will automatically call
#    `~pint.UnitRegistry.setup_matplotlib` so that the axes become unit-aware.

# %%
import xarray as xr
import numpy as np
import pandas as pd

# DataArray
state = np.random.RandomState(51423)
linspace = np.linspace(0, np.pi, 20)
data = (
    50
    * state.normal(1, 0.2, size=(20, 20))
    * (np.sin(linspace * 2) ** 2 * np.cos(linspace + np.pi / 2)[:, None] ** 2)
)
lat = xr.DataArray(
    np.linspace(-90, 90, 20), dims=("lat",), attrs={"units": "\N{DEGREE SIGN}N"}
)
plev = xr.DataArray(
    np.linspace(1000, 0, 20),
    dims=("plev",),
    attrs={"long_name": "pressure", "units": "hPa"},
)
da = xr.DataArray(
    data,
    name="u",
    dims=("plev", "lat"),
    coords={"plev": plev, "lat": lat},
    attrs={"long_name": "zonal wind", "units": "m/s"},
)

# DataFrame
data = state.rand(12, 20)
df = pd.DataFrame(
    (data - 0.4).cumsum(axis=0).cumsum(axis=1)[::1, ::-1],
    index=pd.date_range("2000-01", "2000-12", freq="MS"),
)
df.name = "temperature (\N{DEGREE SIGN}C)"
df.index.name = "date"
df.columns.name = "variable (units)"

# %%
import proplot as pplt

fig = pplt.figure(refwidth=2.5, share=False, suptitle="Automatic subplot formatting")

# Plot DataArray
cmap = pplt.Colormap("PuBu", left=0.05)
ax = fig.subplot(121, yreverse=True)
ax.contourf(da, cmap=cmap, colorbar="t", lw=0.7, ec="k")

# Plot DataFrame
ax = fig.subplot(122, yreverse=True)
ax.contourf(df, cmap="YlOrRd", colorbar="t", lw=0.7, ec="k")
ax.format(xtickminor=False, yformatter="%b", ytickminor=False)

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_apply_cmap:
#
# Changing the colormap
# ---------------------
#
# It is often useful to create custom colormaps on-the-fly,
# without explicitly calling the `~proplot.constructor.Colormap`
# :ref:`constructor function <why_constructor>`. You can do so using the `cmap`
# and `cmap_kw` keywords, available with most `~proplot.axes.PlotAxes` 2D plot
# commands. For example, to create and apply a monochromatic colormap, you can use
# ``cmap='color_name'`` (see the :ref:`colormaps section <ug_cmaps_new>` for more info).
# You can also create on-the-fly "qualitative" `~proplot.colors.DiscreteColormap`\ s
# by passing lists of colors to the keyword `c`, `color`, or `colors`.
#
# Proplot defines :ref:`global defaults <ug_rc>` for four different colormap
# types: :ref:`sequential <ug_cmaps_included>` (setting :rcraw:`cmap.sequential`),
# :ref:`diverging <ug_cmaps_included>` (setting :rcraw:`cmap.diverging`),
# :ref:`cyclic <ug_cmaps_included>` (setting :rcraw:`cmap.cyclic`),
# and :ref:`qualitative <ug_cycles_included>` (setting :rcraw:`cmap.qualitative`).
# To use the default colormap for a given type, pass ``sequential=True``,
# ``diverging=True``, ``cyclic=True``, or ``qualitative=True`` to any plotting
# command. If the colormap type is not explicitly specified, `sequential` is
# used with the default linear normalizer when data is strictly positive
# or negative, and `diverging` is used with the :ref:`diverging normalizer <ug_norm>`
# when the data limits or colormap levels cross zero (see :ref:`below <ug_autonorm>`).

# %%
import proplot as pplt
import numpy as np

# Sample data
N = 18
state = np.random.RandomState(51423)
data = np.cumsum(state.rand(N, N), axis=0)

# Custom defaults of each type
pplt.rc["cmap.sequential"] = "PuBuGn"
pplt.rc["cmap.diverging"] = "PiYG"
pplt.rc["cmap.cyclic"] = "bamO"
pplt.rc["cmap.qualitative"] = "flatui"

# Make plots. Note the default behavior is sequential=True or diverging=True
# depending on whether data contains negative values (see below).
fig = pplt.figure(refwidth=2.2, span=False, suptitle="Colormap types")
axs = fig.subplots(ncols=2, nrows=2)
axs.format(xformatter="none", yformatter="none")
axs[0].pcolor(data, sequential=True, colorbar="l", extend="max")
axs[1].pcolor(data - 5, diverging=True, colorbar="r", extend="both")
axs[2].pcolor(data % 8, cyclic=True, colorbar="l")
axs[3].pcolor(data, levels=pplt.arange(0, 12, 2), qualitative=True, colorbar="r")
types = ("sequential", "diverging", "cyclic", "qualitative")
for ax, typ in zip(axs, types):
    ax.format(title=typ.title() + " colormap")
pplt.rc.reset()

# %%
import proplot as pplt
import numpy as np

# Sample data
N = 20
state = np.random.RandomState(51423)
data = np.cumsum(state.rand(N, N), axis=1) - 6

# Continuous "diverging" colormap
fig = pplt.figure(refwidth=2.3, spanx=False)
ax = fig.subplot(121, title="Diverging colormap with 'cmap'", xlabel="xlabel")
ax.contourf(
    data,
    norm="div",
    cmap=("cobalt", "white", "violet red"),
    cmap_kw={"space": "hsl", "cut": 0.15},
    colorbar="b",
)

# Discrete "qualitative" colormap
ax = fig.subplot(122, title="Qualitative colormap with 'colors'")
ax.contourf(
    data,
    levels=pplt.arange(-6, 9, 3),
    colors=["red5", "blue5", "yellow5", "gray5", "violet5"],
    colorbar="b",
)

import proplot as pplt
import numpy as np

# Sample data
N = 20
state = np.random.RandomState(51423)
data = 11 ** (0.25 * np.cumsum(state.rand(N, N), axis=0))

# Create figure
gs = pplt.GridSpec(ncols=2)
fig = pplt.figure(refwidth=2.3, span=False, suptitle="Normalizer types")

# Different normalizers
ax = fig.subplot(gs[0], title="Default linear normalizer")
ax.pcolormesh(data, cmap="magma", colorbar="b")
ax = fig.subplot(gs[1], title="Logarithmic normalizer with norm='log'")
ax.pcolormesh(data, cmap="magma", norm="log", colorbar="b")


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_norm:
#
# Special normalizers
# -------------------
#
# Proplot includes two new :ref:`"continuous" normalizers <ug_norm>`. The
# `~proplot.colors.SegmentedNorm` normalizer provides even color gradations with respect
# to index for an arbitrary monotonically increasing or decreasing list of levels. This
# is automatically applied if you pass unevenly spaced `levels` to a plotting command,
# or it can be manually applied using e.g. ``norm='segmented'``. This can be useful for
# datasets with unusual statistical distributions or spanning many orders of magnitudes.
#
# The `~proplot.colors.DivergingNorm` normalizer ensures that colormap midpoints lie
# on some central data value (usually ``0``), even if `vmin`, `vmax`, or `levels`
# are asymmetric with respect to the central value. This is automatically applied
# if your data contains negative and positive values (see :ref:`below <ug_autonorm>`),
# or it can be manually applied using e.g. ``diverging=True`` or ``norm='diverging'``.
# It can also be configured to scale colors "fairly" or "unfairly":
#
# * With fair scaling (the default), gradations on either side of the midpoint
#   have equal intensity. If `vmin` and `vmax` are not symmetric about zero, the most
#   intense colormap colors on one side of the midpoint will be truncated.
# * With unfair scaling, gradations on either side of the midpoint are warped
#   so that the full range of colormap colors is always traversed. This configuration
#   should be used with care, as it may lead you to misinterpret your data.
#
# The below examples demonstrate how these normalizers
# affect the interpretation of different datasets.

# %%
import proplot as pplt
import numpy as np

# Sample data
state = np.random.RandomState(51423)
data = 11 ** (2 * state.rand(20, 20).cumsum(axis=0) / 7)

# Linear segmented norm
fig, axs = pplt.subplots(ncols=2, refwidth=2.4)
fig.format(suptitle="Segmented normalizer demo")
ticks = [5, 10, 20, 50, 100, 200, 500, 1000]
for ax, norm in zip(axs, ("linear", "segmented")):
    m = ax.contourf(
        data,
        levels=ticks,
        extend="both",
        cmap="Mako",
        norm=norm,
        colorbar="b",
        colorbar_kw={"ticks": ticks},
    )
    ax.format(title=norm.title() + " normalizer")
# %%
import proplot as pplt
import numpy as np

# Sample data
state = np.random.RandomState(51423)
data1 = (state.rand(20, 20) - 0.485).cumsum(axis=1).cumsum(axis=0)
data2 = (state.rand(20, 20) - 0.515).cumsum(axis=0).cumsum(axis=1)

# Figure
fig, axs = pplt.subplots(nrows=2, ncols=2, refwidth=2.2, order="F")
axs.format(suptitle="Diverging normalizer demo")
cmap = pplt.Colormap("DryWet", cut=0.1)

# Diverging norms
i = 0
for data, mode, fair in zip(
    (data1, data2),
    ("positive", "negative"),
    ("fair", "unfair"),
):
    for fair in ("fair", "unfair"):
        norm = pplt.Norm("diverging", fair=(fair == "fair"))
        ax = axs[i]
        m = ax.contourf(data, cmap=cmap, norm=norm)
        ax.colorbar(m, loc="b")
        ax.format(title=f"{mode.title()}-skewed + {fair} scaling")
        i += 1

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_discrete:
#
# Discrete levels
# ---------------
#
# By default, proplot uses `~proplot.colors.DiscreteNorm` to "discretize"
# the possible colormap colors for contour and pseudocolor `~proplot.axes.PlotAxes`
# commands (e.g., `~proplot.axes.PlotAxes.contourf`, `~proplot.axes.PlotAxes.pcolor`).
# This is analogous to `matplotlib.colors.BoundaryNorm`, except
# `~proplot.colors.DiscreteNorm` can be paired with arbitrary
# continuous normalizers specified by `norm` (see :ref:`above <ug_norm>`).
# Discrete color levels can help readers discern exact numeric values and
# tend to reveal qualitative structure in the data. `~proplot.colors.DiscreteNorm`
# also repairs the colormap end-colors by ensuring the following conditions are met:
#
# #. All colormaps always span the *entire color range*
#    regardless of the `extend` parameter.
# #. Cyclic colormaps always have *distinct color levels*
#    on either end of the colorbar.
#
# To explicitly toggle discrete levels on or off, change :rcraw:`cmap.discrete`
# or pass ``discrete=False`` or ``discrete=True`` to any plotting command
# that accepts a `cmap` argument. The level edges or centers used with
# `~proplot.colors.DiscreteNorm` can be explicitly specified using the `levels` or
# `values` keywords, respectively (`~proplot.utils.arange` and `~proplot.utils.edges`
# are useful for generating `levels` and `values` lists). You can also pass an integer
# to these keywords (or to the `N` keyword) to automatically generate approximately this
# many level edges or centers at "nice" intervals. The algorithm used to generate levels
# is similar to matplotlib's algorithm for generarting contour levels. The default
# number of levels is controlled by :rcraw:`cmap.levels`, and the level selection
# is constrainted by the keywords `vmin`, `vmax`, `locator`, and `locator_kw` -- for
# example, ``vmin=100`` ensures the minimum level is greater than or equal to ``100``,
# and ``locator=5`` ensures a level step size of 5 (see :ref:`this section
# <ug_locators>` for more on locators). You can also use the keywords `negative`,
# `positive`, or `symmetric` to ensure that your levels are strictly negative,
# positive, or symmetric about zero, or use the `nozero` keyword to remove
# the zero level (useful for single-color `~proplot.axes.PlotAxes.contour` plots).

# %%
import proplot as pplt
import numpy as np

# Sample data
state = np.random.RandomState(51423)
data = 10 + state.normal(0, 1, size=(33, 33)).cumsum(axis=0).cumsum(axis=1)

# Figure
fig, axs = pplt.subplots([[1, 1, 2, 2], [0, 3, 3, 0]], ref=3, refwidth=2.3)
axs.format(yformatter="none", suptitle="Discrete vs. smooth colormap levels")

# Pcolor
axs[0].pcolor(data, cmap="viridis", colorbar="l")
axs[0].set_title("Pcolor plot\ndiscrete=True (default)")
axs[1].pcolor(data, discrete=False, cmap="viridis", colorbar="r")
axs[1].set_title("Pcolor plot\ndiscrete=False")

# Imshow
m = axs[2].imshow(data, cmap="oslo", colorbar="b")
axs[2].format(title="Imshow plot\ndiscrete=False (default)", yformatter="auto")

# %%
import proplot as pplt
import numpy as np

# Sample data
state = np.random.RandomState(51423)
data = (20 * (state.rand(20, 20) - 0.4).cumsum(axis=0).cumsum(axis=1)) % 360
levels = pplt.arange(0, 360, 45)

# Figure
gs = pplt.GridSpec(nrows=2, ncols=4, hratios=(1.5, 1))
fig = pplt.figure(refwidth=2.4, right=2)
fig.format(suptitle="DiscreteNorm end-color standardization")

# Cyclic colorbar with distinct end colors
cmap = pplt.Colormap("twilight", shift=-90)
ax = fig.subplot(gs[0, 1:3], title='distinct "cyclic" end colors')
ax.pcolormesh(
    data,
    cmap=cmap,
    levels=levels,
    colorbar="b",
    colorbar_kw={"locator": 90},
)

# Colorbars with different extend values
for i, extend in enumerate(("min", "max", "neither", "both")):
    ax = fig.subplot(gs[1, i], title=f"extend={extend!r}")
    ax.pcolormesh(
        data[:, :10],
        levels=levels,
        cmap="oxy",
        extend=extend,
        colorbar="b",
        colorbar_kw={"locator": 180},
    )

# %% [raw] raw_mimetype="text/restructuredtext" tags=[]
# .. _ug_autonorm:
#
# Auto normalization
# ------------------
#
# By default, colormaps are normalized to span from roughly the minimum
# data value to the maximum data value. However in the presence of outliers,
# this is not desirable. Proplot adds the `robust` keyword to change this
# behavior, inspired by the `xarray keyword
# <http://xarray.pydata.org/en/stable/user-guide/plotting.html#robust>`__
# of the same name. Passing ``robust=True`` to a `~proplot.axes.PlotAxes`
# 2D plot command will limit the default colormap normalization between
# the 2nd and 98th data percentiles. This range can be customized by passing
# an integer to `robust` (e.g. ``robust=90`` limits the normalization range
# between the 5th and 95th percentiles) or by passing a 2-tuple to `robust`
# (e.g. ``robust=(0, 90)`` limits the normalization range between the
# data minimum and the 90th percentile). This can be turned on persistently
# by setting :rcraw:`cmap.robust` to ``True``.
#
# Additionally, `similar to xarray
# <http://xarray.pydata.org/en/stable/user-guide/plotting.html#colormaps>`__,
# proplot can automatically detect "diverging" datasets. By default,
# the 2D `~proplot.axes.PlotAxes` commands will apply the diverging colormap
# :rc:`cmap.diverging` (rather than :rc:`cmap.sequential`) and the diverging
# normalizer `~proplot.colors.DivergingNorm` (rather than `~matplotlib.colors.Normalize`
# -- see :ref:`above <ug_norm>`) if the following conditions are met:
#
# #. If discrete levels are enabled (see :ref:`above <ug_discrete>`) and the
#    level list includes at least 2 negative and 2 positive values.
# #. If discrete levels are disabled (see :ref:`above <ug_discrete>`) and the
#    normalization limits `vmin` and `vmax` are negative and positive.
# #. A colormap was not explicitly passed, or a colormap was passed but it
#    matches the name of a :ref:`known diverging colormap <ug_cmaps_included>`.
#
# The automatic detection of "diverging" datasets can be disabled by
# setting :rcraw:`cmap.autodiverging` to ``False``.

# %%
import proplot as pplt
import numpy as np

N = 20
state = np.random.RandomState(51423)
data = N * 2 + (state.rand(N, N) - 0.45).cumsum(axis=0).cumsum(axis=1) * 10
fig, axs = pplt.subplots(
    nrows=2, ncols=2, refwidth=2, suptitle="Auto normalization demo"
)

# Auto diverging
pplt.rc["cmap.sequential"] = "lapaz_r"
pplt.rc["cmap.diverging"] = "vik"
for i, ax in enumerate(axs[:2]):
    ax.pcolor(data - i * N * 6, colorbar="b")
    ax.format(title="Diverging " + ("on" if i else "off"))

# Auto range
pplt.rc["cmap.sequential"] = "lajolla"
data = data[::-1, :]
data[-1, 0] = 2e3
for i, ax in enumerate(axs[2:]):
    ax.pcolor(data, robust=bool(i), colorbar="b")
    ax.format(title="Robust " + ("on" if i else "off"))
pplt.rc.reset()

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_labels:
#
# Quick labels
# ------------
#
# You can now quickly add labels to `~proplot.axes.PlotAxes.contour`,
# `~proplot.axes.PlotAxes.contourf`, `~proplot.axes.PlotAxes.pcolor`,
# `~proplot.axes.PlotAxes.pcolormesh`, and `~proplot.axes.PlotAxes.heatmap`,
# plots by passing ``labels=True`` to the plotting command. The
# label text is colored black or white depending on the luminance of the underlying
# grid box or filled contour (see the section on :ref:`colorspaces <ug_perceptual>`).
# Contour labels are drawn with `~matplotlib.axes.Axes.clabel` and grid box
# labels are drawn with `~proplot.axes.Axes.text`. You can pass keyword arguments
# to these functions by passing a dictionary to `labels_kw`, and you can
# change the label precision using the `precision` keyword. See the plotting
# command documentation for details.

# %%
import proplot as pplt
import pandas as pd
import numpy as np

# Sample data
state = np.random.RandomState(51423)
data = state.rand(6, 6)
data = pd.DataFrame(data, index=pd.Index(["a", "b", "c", "d", "e", "f"]))

# Figure
fig, axs = pplt.subplots(
    [[1, 1, 2, 2], [0, 3, 3, 0]],
    refwidth=2.3,
    share="labels",
    span=False,
)
axs.format(xlabel="xlabel", ylabel="ylabel", suptitle="Labels demo")

# Heatmap with labeled boxes
ax = axs[0]
m = ax.heatmap(
    data,
    cmap="rocket",
    labels=True,
    precision=2,
    labels_kw={"weight": "bold"},
)
ax.format(title="Heatmap with labels")

# Filled contours with labels
ax = axs[1]
m = ax.contourf(
    data.cumsum(axis=0), cmap="rocket", labels=True, labels_kw={"weight": "bold"}
)
ax.format(title="Filled contours with labels")

# Line contours with labels and no zero level
data = 5 * (data - 0.45).cumsum(axis=0) - 2
ax = axs[2]
ax.contour(data, nozero=True, color="gray8", labels=True, labels_kw={"weight": "bold"})
ax.format(title="Line contours with labels")

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_heatmap:
#
# Heatmap plots
# -------------
#
# The `~proplot.axes.PlotAxes.heatmap` command can be used to draw "heatmaps" of
# 2-dimensional data. This is a convenience function equivalent to
# `~proplot.axes.PlotAxes.pcolormesh`, except the axes are configured with settings
# suitable for heatmaps: fixed aspect ratios (ensuring "square" grid boxes), no
# gridlines, no minor ticks, and major ticks at the center of each box. Among other
# things, this is useful for displaying covariance and correlation matrices, as shown
# below. `~proplot.axes.PlotAxes.heatmap` should generally only be used with
# `~proplot.axes.CartesianAxes`.

# %%
import proplot as pplt
import numpy as np
import pandas as pd

# Covariance data
state = np.random.RandomState(51423)
data = state.normal(size=(10, 10)).cumsum(axis=0)
data = (data - data.mean(axis=0)) / data.std(axis=0)
data = (data.T @ data) / data.shape[0]
data[np.tril_indices(data.shape[0], -1)] = np.nan  # fill half with empty boxes
data = pd.DataFrame(data, columns=list("abcdefghij"), index=list("abcdefghij"))

# Covariance matrix plot
fig, ax = pplt.subplots(refwidth=4.5)
m = ax.heatmap(
    data,
    cmap="ColdHot",
    vmin=-1,
    vmax=1,
    N=100,
    lw=0.5,
    ec="k",
    labels=True,
    precision=2,
    labels_kw={"weight": "bold"},
    clip_on=False,  # turn off clipping so box edges are not cut in half
)
ax.format(
    suptitle="Heatmap demo",
    title="Table of correlation coefficients",
    xloc="top",
    yloc="right",
    yreverse=True,
    ticklabelweight="bold",
    alpha=0,
    linewidth=0,
    tickpad=4,
)
