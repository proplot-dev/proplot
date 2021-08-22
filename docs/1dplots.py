# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _pandas: https://pandas.pydata.org
#
# .. _xarray: http://xarray.pydata.org/en/stable/
#
# .. _ug_1dplots:
#
# Plotting 1D data
# ================
#
# ProPlot adds :ref:`several new features <why_plotting>` to matplotlib's
# plotting commands using the intermediate `~proplot.axes.PlotAxes` subclass.
# For the most part, these additions represent a *superset* of matplotlib -- if
# you are not interested, you can use the plotting commands just like you always
# have. This section documents the features added for 1D plotting commands
# like `~proplot.axes.PlotAxes.plot`, `~proplot.axes.PlotAxes.scatter`,
# and `~proplot.axes.PlotAxes.bar`.


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_1dstd:
#
# Standardized arguments
# ----------------------
#
# Input arguments passed to 1D plotting commands are now uniformly
# standardized. For each command, you can optionally omit the dependent
# variable coordinates, in which case they are inferred from the data
# (see :ref:`xarray and pandas integration <ug_1dintegration>`), or pass
# 2D dependent or independent variable coordinates, in which case the
# plotting command is called for each column of the 2D array(s). If coordinates
# are string labels, they are converted to indices and tick labels using
# `~matplotlib.ticker.FixedLocator` and `~matplotlib.ticker.IndexFormatter`.
# All positional arguments can also be optionally specified as keyword
# arguments (see the individual command documentation).
#
# .. note::
#
#    By default, when just the *x* or *y* axis was explicitly fixed by
#    `~matplotlib.axes.Axes.set_xlim` or `~matplotlib.axes.Axes.set_ylim`
#    (or, equivalently, by passing `xlim` or `ylim` to
#    `proplot.axes.CartesianAxes.format`), ProPlot ignores the out of bounds
#    data when determining the other axis limits. This can be useful if
#    you wish to restrict the view within a large dataset. To disable
#    this feature, pass ``inbounds=False`` to the plotting command or set
#    :rcraw:`axes.inbounds` to ``False`` (see also the :rcraw:`cmap.inbounds`
#    setting and the :ref:`user guide <ug_apply_cmap>`).

# %%
import proplot as pplt
import numpy as np

N = 5
state = np.random.RandomState(51423)
with pplt.rc.context({'axes.prop_cycle': pplt.Cycle('Grays', N=N, left=0.3)}):
    # Sample data
    x = np.linspace(-5, 5, N)
    y = state.rand(N, 5)
    fig = pplt.figure(share=False)

    # Plot by passing both x and y coordinates
    ax = fig.subplot(121)
    ax.area(x, -1 * y / N, stack=True)
    ax.bar(x, y, linewidth=0, alpha=1, width=0.8)
    ax.plot(x, y + 1, linewidth=2)
    ax.scatter(x, y + 2, marker='s', markersize=5**2)
    ax.format(title='Manual x coordinates')

    # Plot by passing just y coordinates
    # Default x coordinates are inferred from DataFrame,
    # inferred from DataArray, or set to np.arange(0, y.shape[0])
    ax = fig.subplot(122)
    ax.area(-1 * y / N, stack=True)
    ax.bar(y, linewidth=0, alpha=1)
    ax.plot(y + 1, linewidth=2)
    ax.scatter(y + 2, marker='s', markersize=5**2)
    ax.format(title='Auto x coordinates')
    fig.format(xlabel='xlabel', ylabel='ylabel')
    fig.format(suptitle='Standardized input demonstration')

# %%
import proplot as pplt
import numpy as np

# Sample data
cycle = pplt.Cycle('davos', right=0.8)
state = np.random.RandomState(51423)
N, M = 400, 20
xmax = 20
x = np.linspace(0, 100, N)
y = 100 * (state.rand(N, M) - 0.42).cumsum(axis=0)

# Plot the data
fig = pplt.figure(refwidth=2.2, share=False)
axs = fig.subplots([[0, 1, 1, 0], [2, 2, 3, 3]], wratios=(2, 1, 1, 2))
axs[0].axvspan(
    0, xmax, zorder=3, edgecolor='red', facecolor=pplt.set_alpha('red', 0.2),
)
for i, ax in enumerate(axs):
    inbounds = i == 1
    title = f'Manual limits inbounds={inbounds}'
    title += ' (default)' if inbounds else ''
    ax.format(
        xmax=(None if i == 0 else xmax),
        title=('Auto x axis limits' if i == 0 else title),
    )
    ax.plot(x, y, cycle=cycle, inbounds=inbounds)
fig.format(
    xlabel='xlabel',
    ylabel='ylabel',
    suptitle='Auto y limits with in-bounds data'
)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_1dintegration:
#
# Pandas and xarray integration
# -----------------------------
#
# The `~proplot.axes.PlotAxes` plotting commands are seamlessly integrated
# with `pandas`_ and `xarray`_. If you omit dependent variable coordinates,
# the plotting command tries to infer them from the `pandas.DataFrame`
# or `xarray.DataArray`. If you did not explicitly set the *x* or *y* axis label
# or :ref:`legend or colorbar <ug_cbars_axes>` title, the plotting command tries to
# retrieve them from the `pandas.DataFrame` or `xarray.DataArray`. You can also pass
# a `~xarray.Dataset`, `~pandas.DataFrame`, or `dict` to any plotting command using the
# `data` keyword, then pass string keys as the data arguments rather than arrays (for
# example, ``ax.plot('y', data=dataset)`` is translated to ``ax.plot(dataset['y'])``).
# Finally, if you pass `pint.Quantity`\ s or `xarray.DataArray`\ s containing
# `pint.Quantity`\ s to a plotting command, ProPlot will automatically call
# `~pint.UnitRegistry.setup_matplotlib` and apply the unit string formatted as
# :rcraw:`unitformat` for the default content labels.
#
# These features restore some of the convenience you get with the builtin
# `pandas`_ and `xarray`_ plotting functions. They are also *optional* --
# installation of pandas and xarray are not required. All of these features
# can be disabled by setting :rcraw:`autoformat` to ``False`` or by passing
# ``autoformat=False`` to any plotting command.

# %%
import xarray as xr
import numpy as np
import pandas as pd

# DataArray
state = np.random.RandomState(51423)
data = (
    np.sin(np.linspace(0, 2 * np.pi, 20))[:, None]
    + state.rand(20, 8).cumsum(axis=1)
)
coords = {
    'x': xr.DataArray(
        np.linspace(0, 1, 20),
        dims=('x',),
        attrs={'long_name': 'distance', 'units': 'km'}
    ),
    'num': xr.DataArray(
        np.arange(0, 80, 10),
        dims=('num',),
        attrs={'long_name': 'parameter'}
    )
}
da = xr.DataArray(
    data, dims=('x', 'num'), coords=coords, name='energy', attrs={'units': 'kJ'}
)

# DataFrame
data = (
    (np.cos(np.linspace(0, 2 * np.pi, 20))**4)[:, None] + state.rand(20, 5) ** 2
)
ts = pd.date_range('1/1/2000', periods=20)
df = pd.DataFrame(data, index=ts, columns=['foo', 'bar', 'baz', 'zap', 'baf'])
df.name = 'data'
df.index.name = 'date'
df.columns.name = 'category'

# %%
import proplot as pplt
fig = pplt.figure(share=False)
fig.format(suptitle='Automatic subplot formatting')

# Plot DataArray
cycle = pplt.Cycle('dark blue', space='hpl', N=da.shape[1])
ax = fig.subplot(121)
ax.scatter(da, cycle=cycle, lw=3, colorbar='t', colorbar_kw={'locator': 20})

# Plot Dataframe
cycle = pplt.Cycle('dark green', space='hpl', N=df.shape[1])
ax = fig.subplot(122)
ax.plot(df, cycle=cycle, lw=3, legend='t', legend_kw={'frame': False})


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_apply_cycle:
#
# Property cycles
# ---------------
#
# It is often useful to create on-the-fly `property cycles
# <https://matplotlib.org/tutorials/intermediate/color_cycle.html#sphx-glr-tutorials-intermediate-color-cycle-py>`__
# and use different property cycles for different plot elements. You can create and
# apply property cycles on-the-fly using the `cycle` and `cycle_kw` keywords, available
# with most `~proplot.axes.PlotAxes` 1D plotting commands. `cycle` and `cycle_kw` are
# passed to the `~proplot.constructor.Cycle` :ref:`constructor function
# <why_constructor>`, and the resulting property cycle is used for the plot. You
# can specify `cycle` once with 2D input data (in which case each column is
# plotted in succession according to the property cycle) or call a plotting
# command multiple times with the same `cycle` argument each time (the property
# cycle is not reset). You can also disable property cycling with
# ``cycle=False``, ``cycle='none'``, or ``cycle=()`` and re-enable the default
# property cycle with ``cycle=True``. For more information on property cycling,
# see the :ref:`color cycles section <ug_cycles>` and `this matplotlib tutorial
# <https://matplotlib.org/tutorials/intermediate/color_cycle.html#sphx-glr-tutorials-intermediate-color-cycle-py>`__.

# %%
import proplot as pplt
import numpy as np

# Sample data
M, N = 50, 4
state = np.random.RandomState(51423)
data1 = (state.rand(M, N) - 0.5).cumsum(axis=0)
data2 = (state.rand(M, N) - 0.5).cumsum(axis=0) * 1.5
data1 += state.rand(M, N)
data2 += state.rand(M, N)

with pplt.rc.context({'lines.linewidth': 3}):
    # Use property cycle for columns of 2D input data
    fig = pplt.figure(share=False)
    ax = fig.subplot(121)
    ax.format(title='Grayscale cycle')
    ax.plot(
        data1 * data2,
        cycle='black',  # cycle from monochromatic colormap
        cycle_kw={'ls': ('-', '--', '-.', ':')}
    )

    # Use property cycle with successive plot() calls
    ax = fig.subplot(122)
    ax.format(title='Colorful cycle')
    for i in range(data1.shape[1]):
        ax.plot(data1[:, i], cycle='Reds', cycle_kw={'N': N, 'left': 0.3})
    for i in range(data1.shape[1]):
        ax.plot(data2[:, i], cycle='Blues', cycle_kw={'N': N, 'left': 0.3})
    fig.format(xlabel='xlabel', ylabel='ylabel', suptitle='Local property cycles demo')


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_lines:
#
# Line plots
# ----------
#
# Line plots can be drawn with `~proplot.axes.PlotAxes.plot` or
# `~proplot.axes.PlotAxes.plotx` (or their aliases, `~proplot.axes.PlotAxes.line`
# or `~proplot.axes.PlotAxes.linex`). For the ``x`` commands, positional
# arguments are interpreted as *x* coordinates or (*y*, *x*) pairs. This is
# analogous to `~proplot.axes.PlotAxes.barh` and `~proplot.axes.PlotAxes.areax`.
# Also, the default *x* bounds for lines drawn with `~proplot.axes.PlotAxes.plot`
# and *y* bounds for lines drawn with `~proplot.axes.PlotAxes.plotx` are now
# "sticky", i.e. there is no padding between the lines and axes edges by default.
#
# Step and stem plots can be drawn with `~proplot.axes.PlotAxes.step`,
# `~proplot.axes.PlotAxes.stepx`, `~proplot.axes.PlotAxes.stem`, and
# `~proplot.axes.PlotAxes.stemx`. Plots of parallel vertical and horizontal
# lines can be drawn with `~proplot.axes.PlotAxes.vlines` and
# `~proplot.axes.PlotAxes.hlines`. You can have different colors for
# "negative" and "positive" lines using ``negpos=True``, ``negcolor=color``, and
# ``poscolor=color`` (the default colors are :rc:`negcolor` and :rc:`poscolor`).

# %%
import proplot as pplt
import numpy as np
state = np.random.RandomState(51423)
gs = pplt.GridSpec(nrows=3, ncols=2)
fig = pplt.figure(refwidth=2.2, span=False, share='labels')

# Vertical vs. horizontal
data = (state.rand(10, 5) - 0.5).cumsum(axis=0)
ax = fig.subplot(gs[0])
ax.format(title='Dependent x-axis')
ax.line(data, lw=2.5, cycle='seaborn')
ax = fig.subplot(gs[1])
ax.format(title='Dependent y-axis')
ax.linex(data, lw=2.5, cycle='seaborn')

# Vertical lines
gray = 'gray7'
data = state.rand(20) - 0.5
ax = fig.subplot(gs[2])
ax.area(data, color=gray, alpha=0.2)
ax.vlines(data, negpos=True, lw=2)
ax.format(title='Vertical lines')

# Horizontal lines
ax = fig.subplot(gs[3])
ax.areax(data, color=gray, alpha=0.2)
ax.hlines(data, negpos=True, lw=2)
ax.format(title='Horizontal lines')

# Step
ax = fig.subplot(gs[4])
data = state.rand(20, 4).cumsum(axis=1).cumsum(axis=0)
cycle = ('gray6', 'blue7', 'red7', 'gray4')
ax.step(data, cycle=cycle, labels=list('ABCD'), legend='ul', legend_kw={'ncol': 2})
ax.format(title='Step plot')

# Stems
ax = fig.subplot(gs[5])
data = state.rand(20)
ax.stem(data)
ax.format(title='Stem plot')
fig.format(suptitle='Line plots demo', xlabel='xlabel', ylabel='ylabel')


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_scatter:
#
# Scatter plots
# -------------
#
# The `~proplot.axes.PlotAxes.scatter` command now permits omitting *x*
# coordinates and accepts 2D *y* coordinates, just like `~proplot.axes.PlotAxes.plot`.
# As with `~proplot.axes.PlotAxes.plotx`, the `~proplot.axes.PlotAxes.scatterx`
# command is just like `~proplot.axes.PlotAxes.scatter`, except positional
# arguments are interpreted as *x* coordinates and (*y*, *x*) pairs.
# `~proplot.axes.PlotAxes.scatter` also now accepts keywords
# that look like `~proplot.axes.PlotAxes.plot` keywords (e.g., `color` instead of
# `c` and `markersize` instead of `s`). This way, `~proplot.axes.PlotAxes.scatter`
# can be used simply to "plot markers, not lines" without changing the input
# arguments relative to `~proplot.axes.PlotAxes.plot`.
#
# The property cycler used by `~proplot.axes.PlotAxes.scatter` can be changed
# using the `cycle` keyword argument, and unlike matplotlib it can include
# properties like `marker` and `markersize`. The colormap `cmap` and normalizer
# `norm` used with the optional `c` color array are now passed through the
# `~proplot.constructor.Colormap` and `~proplot.constructor.Norm` constructor
# functions, and the the `s` marker size array can now be conveniently scaled using
# the keywords `smin` and `smax` (analogous to `vmin` and `vmax` used for colors).

# %%
import proplot as pplt
import numpy as np
import pandas as pd

# Sample data
state = np.random.RandomState(51423)
x = (state.rand(20) - 0).cumsum()
data = (state.rand(20, 4) - 0.5).cumsum(axis=0)
data = pd.DataFrame(data, columns=pd.Index(['a', 'b', 'c', 'd'], name='label'))

# Figure
gs = pplt.GridSpec(ncols=2, nrows=2)
fig = pplt.figure(refwidth=2.2, share='labels', span=False)

# Vertical vs. horizontal
ax = fig.subplot(gs[0])
ax.format(title='Dependent x-axis')
ax.scatter(data, cycle='538')
ax = fig.subplot(gs[1])
ax.format(title='Dependent y-axis')
ax.scatterx(data, cycle='538')

# Scatter plot with property cycler
ax = fig.subplot(gs[2])
ax.format(title='With property cycle')
obj = ax.scatter(
    x, data, legend='ul', legend_kw={'ncols': 2},
    cycle='Set2', cycle_kw={'m': ['x', 'o', 'x', 'o'], 'ms': [5, 10, 20, 30]}
)

# Scatter plot with colormap
ax = fig.subplot(gs[3])
ax.format(title='With colormap')
data = state.rand(2, 100)
obj = ax.scatter(
    *data,
    s=state.rand(100), smin=3, smax=60, marker='o',
    c=data.sum(axis=0), cmap='maroon',
    colorbar='lr', colorbar_kw={'label': 'label'},
)
fig.format(suptitle='Scatter plot demo', xlabel='xlabel', ylabel='ylabel')


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_bar:
#
# Bar plots and area plots
# ------------------------
#
# The `~proplot.axes.PlotAxes.bar` and `~proplot.axes.PlotAxes.barh` commands
# apply default *x* or *y* coordinates if you failed to provide them explicitly
# and can *group* or *stack* columns of data if you pass 2D arrays instead of
# 1D arrays -- just like `pandas`_. Similarly, `~proplot.axes.PlotAxes.fill_between`
# and `~proplot.axes.PlotAxes.fill_betweenx` also apply default *x* or *y* coordinates
# if you failed to provide them explicitly, and can *stack* or *overlay* columns of
# data by passing 2D arrays instead of 1D arrays. You can also use the shorthands
# `~proplot.axes.PlotAxes.area` and `~proplot.axes.PlotAxes.areax` instead of
# ``fill_between``.
#
# For both bar and area plots, you can have different colors for "negative" and
# "positive" regions using ``negpos=True``, ``negcolor=color``, and ``poscolor=color``
# (the default colors are :rc:`negcolor` and :rc:`poscolor`). Also, the default *x*
# bounds for shading drawn with `~proplot.axes.PlotAxes.area` and *y* bounds for
# shading drawn with `~proplot.axes.PlotAxes.areax` is now "sticky", i.e. there
# is no padding between the shading and axes edges by default.

# %%
import proplot as pplt
import numpy as np
import pandas as pd

# Sample data
state = np.random.RandomState(51423)
data = state.rand(5, 5).cumsum(axis=0).cumsum(axis=1)[:, ::-1]
data = pd.DataFrame(
    data, columns=pd.Index(np.arange(1, 6), name='column'),
    index=pd.Index(['a', 'b', 'c', 'd', 'e'], name='row idx')
)

# Figure
pplt.rc.abc = 'a.'
pplt.rc.titleloc = 'l'
gs = pplt.GridSpec(nrows=2, hratios=(3, 2))
fig = pplt.figure(refaspect=2, refwidth=4.8, share=False)

# Side-by-side bars
ax = fig.subplot(gs[0])
obj = ax.bar(
    data, cycle='Reds', edgecolor='red9', colorbar='ul', colorbar_kw={'frameon': False}
)
ax.format(xlocator=1, xminorlocator=0.5, ytickminor=False, title='Side-by-side')

# Stacked bars
ax = fig.subplot(gs[1])
obj = ax.barh(
    data.iloc[::-1, :], cycle='Blues', edgecolor='blue9', legend='ur', stack=True,
)
ax.format(title='Stacked')
fig.format(grid=False, suptitle='Bar plot demo')

# %%
import proplot as pplt
import numpy as np

# Sample data
state = np.random.RandomState(51423)
data = state.rand(5, 3).cumsum(axis=0)
cycle = ('gray3', 'gray5', 'gray7')

# Figure
fig = pplt.figure(refwidth=2.3, share=False)

# Overlaid area patches
ax = fig.subplot(121)
ax.area(
    np.arange(5), data, data + state.rand(5)[:, None], cycle=cycle, alpha=0.7,
    legend='uc', legend_kw={'center': True, 'ncols': 2, 'labels': ['z', 'y', 'qqqq']},
)
ax.format(title='Fill between columns')

# Stacked area patches
ax = fig.subplot(122)
ax.area(
    np.arange(5), data, stack=True, cycle=cycle, alpha=0.8,
    legend='ul', legend_kw={'center': True, 'ncols': 2, 'labels': ['z', 'y', 'qqqq']},
)
ax.format(title='Stack between columns')
fig.format(grid=False, xlabel='xlabel', ylabel='ylabel', suptitle='Area plot demo')

# %%
import proplot as pplt
import numpy as np

# Sample data
state = np.random.RandomState(51423)
data = 4 * (state.rand(40) - 0.5)

# Figure
fig, axs = pplt.subplots(nrows=2, refaspect=2, figwidth=5)
axs.format(
    xmargin=0, xlabel='xlabel', ylabel='ylabel', grid=True,
    suptitle='Positive and negative colors demo',
)
for ax in axs:
    ax.axhline(0, color='k', linewidth=1)  # zero line

# Bar plot
ax = axs[0]
ax.bar(data, width=1, negpos=True)
ax.format(title='Bar plot')

# Area plot
ax = axs[1]
ax.area(data, negpos=True, lw=0.5, edgecolor='k')
ax.format(title='Area plot')

# Reset title styles changed above
pplt.rc.reset()


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_errorbars:
#
# Shading and error bars
# ----------------------
#
# Error bars and error shading can be quickly added on-the-fly to
# `~proplot.axes.PlotAxes.line`, `~proplot.axes.PlotAxes.linex`
# (equivalently, `~proplot.axes.PlotAxes.plot`,
# `~proplot.axes.PlotAxes.plotx`), `~proplot.axes.PlotAxes.scatter`,
# `~proplot.axes.PlotAxes.scatterx`, `~proplot.axes.PlotAxes.bar`, and
# `~proplot.axes.PlotAxes.barh` plots using any of several keyword arguments.
#
# If you pass 2D arrays to these commands with ``mean=True``, ``means=True``,
# ``median=True``, or ``medians=True``, the means or medians of each column are
# drawn as lines, points, or bars, while *error bars* or *error shading*
# indicates the spread of the distribution in each column. Invalid data is
# ignored. You can also specify the error bounds *manually* with the `bardata`,
# `boxdata`, `shadedata`, and `fadedata` keywords. These commands can draw and
# style thin error bars (the ``bar`` keywords), thick "boxes" overlaid on top of
# these bars (the ``box`` keywords; think of them as miniature boxplots), a
# transparent primary shading region (the ``shade`` keywords), and a more
# transparent secondary shading region (the ``fade`` keywords). See the documentation
# on the plotting commands for details.


# %%
import numpy as np
import pandas as pd

# Sample data
# Each column represents a distribution
state = np.random.RandomState(51423)
data = state.rand(20, 8).cumsum(axis=0).cumsum(axis=1)[:, ::-1]
data = data + 20 * state.normal(size=(20, 8)) + 30
data = pd.DataFrame(data, columns=np.arange(0, 16, 2))
data.columns.name = 'column number'
data.name = 'variable'

# Calculate error data
# Passed to 'errdata' in the 3rd subplot example
means = data.mean(axis=0)
means.name = data.name  # copy name for formatting
fadedata = np.percentile(data, (5, 95), axis=0)  # light shading
shadedata = np.percentile(data, (25, 75), axis=0)  # dark shading

# %%
import proplot as pplt
import numpy as np

# Loop through "vertical" and "horizontal" versions
varray = [[1], [2], [3]]
harray = [[1, 1], [2, 3], [2, 3]]
for orientation, array in zip(('horizontal', 'vertical'), (harray, varray)):
    # Figure
    fig = pplt.figure(refwidth=4, refaspect=1.5, share=False)
    axs = fig.subplots(array, hratios=(2, 1, 1))
    axs.format(abc='A.', suptitle=f'Indicating {orientation} error bounds')

    # Medians and percentile ranges
    ax = axs[0]
    kw = dict(
        color='light red', legend=True,
        median=True, barpctile=90, boxpctile=True,
        # median=True, barpctile=(5, 95), boxpctile=(25, 75)  # equivalent
    )
    if orientation == 'horizontal':
        ax.barh(data, **kw)
    else:
        ax.bar(data, **kw)
    ax.format(title='Bar plot')

    # Means and standard deviation range
    ax = axs[1]
    kw = dict(
        color='denim', marker='x', markersize=8**2, linewidth=0.8,
        label='mean', shadelabel=True,
        mean=True, shadestd=1,
        # mean=True, shadestd=(-1, 1)  # equivalent
    )
    if orientation == 'horizontal':
        ax.scatterx(data, legend='b', legend_kw={'ncol': 1}, **kw)
    else:
        ax.scatter(data, legend='ll', **kw)
    ax.format(title='Marker plot')

    # User-defined error bars
    ax = axs[2]
    kw = dict(
        shadedata=shadedata, fadedata=fadedata,
        label='mean', shadelabel='50% CI', fadelabel='90% CI',
        color='ocean blue', barzorder=0, boxmarker=False,
    )
    if orientation == 'horizontal':
        ax.linex(means, legend='b', legend_kw={'ncol': 1}, **kw)
    else:
        ax.line(means, legend='ll', **kw)
    ax.format(title='Line plot')


# %% [raw] raw_mimetype="text/restructuredtext" tags=[]
# .. _ug_hist:
#
# Histogram plots
# ---------------
#
# Vertical and horizontal histograms can be drawn with
# `~proplot.axes.PlotAxes.hist` and `~proplot.axes.PlotAxes.histh`.
# As with the other plotting commands, multiple histograms can be
# drawn by passing 2D arrays instead of 1D arrays, and the color
# cycle used to color histograms can be changed on-the-fly using
# the `cycle` and `cycle_kw` keywords. Likewise, 2D histograms can
# be drawn with the `~proplot.axes.PlotAxes.hist2d`
# `~proplot.axes.PlotAxes.hexbin` commands, and their colormaps can
# be changed on-the-fly with the `cmap` and `cmap_kw` keywords (see
# the :ref:`2d plotting section <ug_apply_cmap>`). Marginal distributions
# for the 2D histograms can be added using :ref:`panel axes <ug_panels>`.
# In the future, ProPlot may include options for drawing "smooth"
# kernel density estimations with these commands.

# %%
import proplot as pplt
import numpy as np

# Sample data
M, N = 300, 3
state = np.random.RandomState(51423)
x = state.normal(size=(M, N)) + state.rand(M)[:, None] * np.arange(N) + 2 * np.arange(N)

# Sample overlayed histograms
fig, ax = pplt.subplots(refwidth=4, refaspect=(3, 2))
ax.format(suptitle='Overlaid histograms', xlabel='distribution', ylabel='count')
res = ax.hist(
    x, pplt.arange(-3, 8, 0.2), alpha=0.7,
    cycle=('indigo9', 'gray3', 'red9'), labels=list('abc'), legend='ul',
)

# %%
import proplot as pplt
import numpy as np

# Sample data
N = 500
state = np.random.RandomState(51423)
x = state.normal(size=(N,))
y = state.normal(size=(N,))
bins = pplt.arange(-3, 3, 0.25)

# Histogram with marginal distributions
fig, axs = pplt.subplots(ncols=2, refwidth=2.3)
axs.format(
    abc='A.', abcloc='l', titleabove=True,
    ylabel='y axis', suptitle='Histograms with marginal distributions'
)
colors = ('indigo9', 'red9')
titles = ('Group 1', 'Group 2')
for ax, which, color, title in zip(axs, 'lr', colors, titles):
    ax.hist2d(
        x, y, bins, vmin=0, vmax=10, levels=50,
        cmap=color, colorbar='b', colorbar_kw={'label': 'count'}
    )
    color = pplt.scale_luminance(color, 1.5)  # histogram colors
    px = ax.panel(which, space=0)
    px.hist(y, bins, lw=0, color=color, vert=False)  # or orientation='horizontal'
    px.format(grid=False, xlocator=[], xreverse=(which == 'l'))
    px = ax.panel('t', space=0)
    px.hist(x, bins, lw=0, color=color)
    px.format(grid=False, ylocator=[], title=title, titleloc='l')

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_boxplots:
#
# Box plots and violin plots
# --------------------------
#
# Vertical and horizontal box and violin plots can be drawn using
# `~proplot.axes.PlotAxes.boxplot`, `~proplot.axes.PlotAxes.violinplot`,
# `~proplot.axes.PlotAxes.boxploth`, and `~proplot.axes.PlotAxes.violinploth` (or
# their new shorthands, `~proplot.axes.PlotAxes.box`, `~proplot.axes.PlotAxes.violin`,
# `~proplot.axes.PlotAxes.boxh`, and `~proplot.axes.PlotAxes.violinh`). The
# ProPlot versions employ aesthetically pleasing defaults and permit flexible
# configuration using keywords like `color`, `barcolor`, and `fillcolor`.
# They also automatically apply axis labels based on the `~pandas.DataFrame`
# or `~xarray.DataArray` column labels.

# %%
import proplot as pplt
import numpy as np
import pandas as pd

# Sample data
N = 500
state = np.random.RandomState(51423)
data1 = state.normal(size=(N, 5)) + 2 * (state.rand(N, 5) - 0.5) * np.arange(5)
data1 = pd.DataFrame(data1, columns=pd.Index(list('abcde'), name='label'))
data2 = state.rand(100, 7)
data2 = pd.DataFrame(data2, columns=pd.Index(list('abcdefg'), name='label'))

# Figure
fig, axs = pplt.subplots([[1, 1, 2, 2], [0, 3, 3, 0]], span=False)
axs.format(
    abc='A.', titleloc='l', grid=False,
    suptitle='Boxes and violins demo'
)

# Box plots
ax = axs[0]
obj1 = ax.boxplot(data1, means=True, marker='x', meancolor='r', fillcolor='gray4')
ax.format(title='Box plots')

# Violin plots
ax = axs[1]
obj2 = ax.violinplot(data1, fillcolor='gray6', means=True, points=100)
ax.format(title='Violin plots')

# Boxes with different colors
ax = axs[2]
colors = pplt.get_colors('pastel2')  # list of colors from the cycle
ax.boxplot(data2, fillcolor=colors, orientation='horizontal')
ax.format(title='Multiple colors', ymargin=0.15)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_parametric:
#
# Parametric plots
# ----------------
#
# Parametric plots can be drawn using the new `~proplot.axes.PlotAxes.parametric`
# command. This creates `~matplotlib.collections.LineCollection`\ s that map
# individual line segments to individual colors, where each segment represents a
# "parametric" coordinate (e.g., time). The parametric coordinates are specified with
# a third positional argument or with the keywords `c`, `color`, `colors` or `values`.
# Representing parametric coordinates with colors instead of text labels can be
# cleaner. The below example makes a simple `~proplot.axes.PlotAxes.parametric`
# plot with a colorbar indicating the parametric coordinate.

# %%
import proplot as pplt
import numpy as np
import pandas as pd
gs = pplt.GridSpec(ncols=2, wratios=(2, 1))
fig = pplt.figure(figwidth='16cm', refaspect=(2, 1), share=False)
fig.format(suptitle='Parametric plots demo')
cmap = 'IceFire'

# Sample data
state = np.random.RandomState(51423)
N = 50
x = (state.rand(N) - 0.52).cumsum()
y = state.rand(N)
c = np.linspace(-N / 2, N / 2, N)  # color values
c = pd.Series(c, name='parametric coordinate')

# Parametric line with smooth gradations
ax = fig.subplot(gs[0])
m = ax.parametric(
    x, y, c, interp=10, capstyle='round', joinstyle='round',
    lw=7, cmap=cmap, colorbar='b', colorbar_kw={'locator': 5}
)
ax.format(xlabel='xlabel', ylabel='ylabel', title='Line with smooth gradations')

# Sample data
N = 12
radii = np.linspace(1, 0.2, N + 1)
angles = np.linspace(0, 4 * np.pi, N + 1)
x = radii * np.cos(1.4 * angles)
y = radii * np.sin(1.4 * angles)
c = np.linspace(-N / 2, N / 2, N + 1)

# Parametric line with stepped gradations
ax = fig.subplot(gs[1])
m = ax.parametric(x, y, c, cmap=cmap, lw=15)
ax.format(
    xlim=(-1, 1), ylim=(-1, 1), title='Step gradations',
    xlabel='cosine angle', ylabel='sine angle'
)
ax.colorbar(m, loc='b', maxn=10, label='parametric coordinate')
