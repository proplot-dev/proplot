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
# .. _ug_1dplots:
#
# Plotting 1-dimensional data
# ===========================
#
# ProPlot adds new features to various `~matplotlib.axes.Axes` plotting
# methods using a set of wrapper functions. When a plotting method like
# `~matplotlib.axes.Axes.plot` is "wrapped" by one of these functions, it
# accepts the same parameters as the wrapper. These features are a strict
# *superset* of the matplotlib API -- if you want, you can use the plotting
# methods exactly as you always have.
#
# This section documents the features added by wrapper functions to 1D
# plotting commands like `~matplotlib.axes.Axes.plot`,
# `~matplotlib.axes.Axes.scatter`, `~matplotlib.axes.Axes.bar`, and
# `~matplotlib.axes.Axes.barh`.


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_cycle_changer:
#
# Property cycles
# ---------------
#
# It is often desirable to use different `property cycles
# <https://matplotlib.org/tutorials/intermediate/color_cycle.html#sphx-glr-tutorials-intermediate-color-cycle-py>`__
# for different axes or different plot elements. To enable this, the
# `~proplot.wrappers.cycle_changer` adds the `cycle` and `cycle_kw` to the 1D
# plotting methods. These arguments are passed to the
# `~proplot.constructor.Cycle` constructor function, and the resulting property
# cycle is used to style the input data. ProPlot iterates through property
# cycle properties when (1) making multiple calls to a plotting command, or (2)
# plotting successive columns of 2-dimensional input data. For more information
# on property cycles, see the :ref:`color cycles section <ug_cycles>` and `this
# matplotlib tutorial
# <https://matplotlib.org/tutorials/intermediate/color_cycle.html#sphx-glr-tutorials-intermediate-color-cycle-py>`__.

# %%
import proplot as plot
import numpy as np
state = np.random.RandomState(51423)
data1 = state.rand(6, 4)
data2 = state.rand(6, 4) * 1.5
with plot.rc.context({'lines.linewidth': 3}):
    fig, axs = plot.subplots(ncols=2)
    axs.format(xlabel='xlabel', ylabel='ylabel', suptitle='Local property cycles demo')

    # Property cycles specific to datasets
    axs[0].plot(data1, cycle='Reds', cycle_kw={'left': 0.3})
    axs[0].plot(data2, cycle='Blues', cycle_kw={'left': 0.3})
    axs[1].plot(
        data1 * data2,
        cycle='black',
        cycle_kw={'linestyle': ('-', '--', '-.', ':')}
    )


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_1dstd:
#
# Standardized arguments
# ----------------------
#
# The `~proplot.wrappers.standardize_1d` wrapper is used to standardize
# positional arguments across all 1D plotting methods.
# `~proplot.wrappers.standardize_1d` allows you to optionally omit *x*
# coordinates, in which case they are inferred from the data. It also permits
# passing 2D *y* coordinate arrays to any plotting method, in which case the
# plotting method is called for each column of the array.

# %%
import proplot as plot
import numpy as np

# Figure and sample data
N = 5
state = np.random.RandomState(51423)
with plot.rc.context({'axes.prop_cycle': plot.Cycle('Grays', N=N, left=0.3)}):
    fig, axs = plot.subplots(ncols=2, share=False)
    x = np.linspace(-5, 5, N)
    y = state.rand(N, 5)
    axs.format(xlabel='xlabel', ylabel='ylabel', ymargin=0.05)
    axs.format(suptitle='Standardized arguments demonstration')

    # Plot by passing both x and y coordinates
    ax = axs[0]
    ax.area(x, -1 * y / N, stacked=True)
    ax.bar(x, y, linewidth=0, alpha=1, width=0.8 * (x[1] - x[0]))
    ax.plot(x, y + 1, linewidth=2)
    ax.scatter(x, y + 2, marker='s', markersize=5**2)
    ax.format(title='Manual x coordinates')

    # Plot by passing just y coordinates
    # Default x coordinates are inferred from DataFrame,
    # inferred from DataArray, or set to np.arange(0, y.shape[0])
    ax = axs[1]
    ax.area(-1 * y / N, stacked=True)
    ax.bar(y, linewidth=0, alpha=1)
    ax.plot(y + 1, linewidth=2)
    ax.scatter(y + 2, marker='s', markersize=5**2)
    ax.format(title='Auto x coordinates')


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_1dintegration:
#
# Pandas and xarray integration
# -----------------------------
#
# The `~proplot.wrappers.standardize_1d` wrapper integrates 1D plotting
# methods with pandas `~pandas.DataFrame`\ s and xarray `~xarray.DataArray`\ s.
# When you pass a DataFrame or DataArray to any plotting command, the x-axis
# label, y-axis label, legend label, colorbar label, and/or title are
# configured from the metadata. This restores some of the convenience you get
# with the builtin `pandas <https://pandas.pydata.org>`__ and `xarray
# <https://pandas.pydata.org>`__ plotting functions. This feature is
# *optional*; installation of pandas and xarray are not required.

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
da = xr.DataArray(data, dims=('x', 'cat'), coords={
    'x': xr.DataArray(
        np.linspace(0, 1, 20),
        dims=('x',),
        attrs={'long_name': 'distance', 'units': 'km'}
    ),
    'cat': xr.DataArray(
        np.arange(0, 80, 10),
        dims=('cat',),
        attrs={'long_name': 'parameter', 'units': 'K'}
    )
}, name='position series')

# DataFrame
data = (
    (np.cos(np.linspace(0, 2 * np.pi, 20))**4)[:, None] + state.rand(20, 5)**2
)
ts = pd.date_range('1/1/2000', periods=20)
df = pd.DataFrame(data, index=ts, columns=['foo', 'bar', 'baz', 'zap', 'baf'])
df.name = 'time series'
df.index.name = 'time (s)'
df.columns.name = 'columns'

# %%
import proplot as plot
fig, axs = plot.subplots(ncols=2, axwidth=2.2, share=0)
axs.format(suptitle='Automatic subplot formatting')

# Plot DataArray
cycle = plot.Cycle('dark blue', fade=90, space='hpl', N=da.shape[1])
axs[0].scatter(da, cycle=cycle, lw=3, colorbar='ul', colorbar_kw={'locator': 20})

# Plot Dataframe
cycle = plot.Cycle('dark green', fade=90, space='hpl', N=df.shape[1])
axs[1].plot(df, cycle=cycle, lw=3, legend='uc')


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_errorbars:
#
# Adding error bars
# -----------------
#
# The `~proplot.wrappers.add_errorbars` wrapper lets you draw error bars
# on-the-fly by passing certain keyword arguments to
# `~matplotlib.axes.Axes.plot`, `~matplotlib.axes.Axes.scatter`,
# `~matplotlib.axes.Axes.bar`, or `~matplotlib.axes.Axes.barh`.
#
# If you pass 2D arrays to these methods with ``means=True`` or
# ``medians=True``, the means or medians of each column are drawn as points,
# lines, or bars, and error bars are drawn to represent the spread in each
# column. `~proplot.wrappers.add_errorbars` lets you draw both thin error
# "bars" with optional whiskers, and thick error "boxes" overlayed on top of
# these bars (this can be used to represent different percentil ranges).
# Instead of using 2D arrays, you can also pass error bar coordinates
# *manually* with the `bardata` and `boxdata` keyword arguments. See
# `~proplot.wrappers.add_errorbars` for details.

# %%
import proplot as plot
import numpy as np
import pandas as pd
plot.rc['title.loc'] = 'uc'
plot.rc['axes.ymargin'] = plot.rc['axes.xmargin'] = 0.05
state = np.random.RandomState(51423)
data = (
    state.rand(20, 8).cumsum(axis=0).cumsum(axis=1)[:, ::-1]
    + 20 * state.normal(size=(20, 8)) + 30
)
fig, axs = plot.subplots(
    nrows=3, aspect=1.5, axwidth=4,
    share=0, hratios=(2, 1, 1)
)
axs.format(suptitle='Error bars with various plotting commands')
axs[1:].format(xlabel='column number', xticks=1, xgrid=False)

# Asking add_errorbars to calculate bars
ax = axs[0]
obj = ax.barh(data, color='red orange', means=True)
ax.format(title='Column statistics')
ax.format(ylabel='column number', title='Bar plot', ygrid=False)

# Showing a standard deviation range instead of percentile range
ax = axs[1]
ax.scatter(
    data, color='k', marker='_', markersize=50,
    medians=True, barstd=True, boxes=False, capsize=2,
    barcolor='gray6', barrange=(-1, 1), barzorder=0, barlw=1,
)
ax.format(title='Scatter plot')

# Supplying error bar data manually
ax = axs[2]
boxdata = np.percentile(data, (25, 75), axis=0)
bardata = np.percentile(data, (5, 95), axis=0)
ax.plot(
    data.mean(axis=0), boxes=True,
    edgecolor='k', color='gray9',
    boxdata=boxdata, bardata=bardata, barzorder=0,
    boxcolor='gray7', barcolor='gray7', boxmarker=False,
)
ax.format(title='Line plot')
plot.rc.reset()


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_bar:
#
# Bar plots and area plots
# ------------------------
#
# The `~matplotlib.axes.Axes.bar` and `~matplotlib.axes.Axes.barh` methods
# are wrapped by `~proplot.wrappers.bar_wrapper`,
# `~proplot.wrappers.cycle_changer`, and `~proplot.wrappers.standardize_1d`.
# You can now *group* or *stack* columns of data by passing 2D arrays to
# `~matplotlib.axes.Axes.bar` or `~matplotlib.axes.Axes.barh`, just like in
# `pandas`. Also, `~matplotlib.axes.Axes.bar` and `~matplotlib.axes.Axes.barh`
# now employ "default" *x* coordinates if you failed to provide them
# explicitly, just like `~matplotlib.axes.Axes.plot`.
#
# To make filled "area" plots, use the new `~proplot.axes.Axes.area` and
# `~proplot.axes.Axes.areax` methods. These are alises for
# `~matplotlib.axes.Axes.fill_between` and
# `~matplotlib.axes.Axes.fill_betweenx`, which are now wrapped by
# `~proplot.wrappers.fill_between_wrapper` and
# `~proplot.wrappers.fill_betweenx_wrapper`. You can now *stack* or *overlay*
# columns of data by passing 2D arrays to `~proplot.axes.Axes.area` and
# `~proplot.axes.Axes.areax`, just like in `pandas`. You can also now draw
# area plots that *change color* when the fill boundaries cross each other by
# passing ``negpos=True`` to `~matplotlib.axes.Axes.fill_between`. The most
# common use case for this is highlighting negative and positive areas with
# different colors, as shown below.

# %%
import proplot as plot
import numpy as np
import pandas as pd
plot.rc.titleloc = 'uc'
plot.rc.margin = 0.05
fig, axs = plot.subplots(nrows=2, aspect=2, axwidth=5, share=0, hratios=(3, 2))
state = np.random.RandomState(51423)
data = state.rand(5, 5).cumsum(axis=0).cumsum(axis=1)[:, ::-1]
data = pd.DataFrame(
    data, columns=pd.Index(np.arange(1, 6), name='column'),
    index=pd.Index(['a', 'b', 'c', 'd', 'e'], name='row idx')
)

# Side-by-side bars
ax = axs[0]
obj = ax.bar(
    data, cycle='Reds', colorbar='ul',
    edgecolor='red9', colorbar_kw={'frameon': False}
)
ax.format(
    xlocator=1, xminorlocator=0.5, ytickminor=False,
    title='Side-by-side', suptitle='Bar plot wrapper demo'
)

# Stacked bars
ax = axs[1]
obj = ax.barh(
    data.iloc[::-1, :], cycle='Blues',
    legend='ur', edgecolor='blue9', stacked=True
)
ax.format(title='Stacked')
axs.format(grid=False)
plot.rc.reset()

# %%
import proplot as plot
import numpy as np
plot.rc.margin = 0
fig, axs = plot.subplots(array=[[1, 2], [3, 3]], hratios=(1, 0.8), share=0)
axs.format(xlabel='xlabel', ylabel='ylabel', suptitle='Area plot demo')
state = np.random.RandomState(51423)
data = state.rand(5, 3).cumsum(axis=0)
cycle = ('gray3', 'gray5', 'gray7')

# Overlaid and stacked area patches
ax = axs[0]
ax.area(
    np.arange(5), data, data + state.rand(5)[:, None], cycle=cycle, alpha=0.5,
    legend='uc', legend_kw={'center': True, 'ncols': 2, 'labels': ['z', 'y', 'qqqq']},
)
ax.format(title='Fill between columns')
ax = axs[1]
ax.area(
    np.arange(5), data, stacked=True, cycle=cycle, alpha=0.8,
    legend='ul', legend_kw={'center': True, 'ncols': 2, 'labels': ['z', 'y', 'qqqq']},
)
ax.format(title='Stack between columns')

# Positive and negative color area patches
ax = axs[2]
data = 5 * (state.rand(20) - 0.5)
ax.area(data, negpos=True, negcolor='blue7', poscolor='red7')
ax.format(title='Positive and negative colors', xlabel='xlabel', ylabel='ylabel')
axs.format(grid=False)
plot.rc.reset()


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_boxplots:
#
# Box plots and violin plots
# --------------------------
#
# The `~matplotlib.axes.Axes.boxplot` and `~matplotlib.axes.Axes.violinplot`
# methods are now wrapped with `~proplot.wrappers.boxplot_wrapper`,
# `~proplot.wrappers.violinplot_wrapper`, `~proplot.wrappers.cycle_changer`,
# and `~proplot.wrappers.standardize_1d`. These wrappers add some useful
# options and apply aesthetically pleasing default settings. They also
# automatically apply axis labels based on the `~pandas.DataFrame` column
# labels or the input *x* coordinate labels.

# %%
import proplot as plot
import numpy as np
import pandas as pd
N = 500
state = np.random.RandomState(51423)
fig, axs = plot.subplots(ncols=2, axwidth=2.5)
data = state.normal(size=(N, 5)) + 2 * (state.rand(N, 5) - 0.5) * np.arange(5)
data = pd.DataFrame(
    data,
    columns=pd.Index(['a', 'b', 'c', 'd', 'e'], name='xlabel')
)
axs.format(
    ymargin=0.1, xmargin=0.1, grid=False,
    suptitle='Boxes and violins demo'
)

# Box plots
ax = axs[0]
obj1 = ax.boxplot(
    data, lw=0.7, marker='x', fillcolor='gray5',
    medianlw=1, mediancolor='k'
)
ax.format(title='Box plots', titleloc='uc')

# Violin plots
ax = axs[1]
obj2 = ax.violinplot(
    data, lw=0.7, fillcolor='gray7',
    points=500, bw_method=0.3, means=True
)
ax.format(title='Violin plots', titleloc='uc')


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_parametric:
#
# Parametric plots
# ----------------
#
# To make "parametric" plots, use the new `~proplot.axes.Axes.parametric`
# method. Parametric plots are
# `~matplotlib.collections.LineCollection`\ s that map individual line
# segments to individual colors, where each segment represents a "parametric"
# coordinate (e.g. time). The parametric coordinates are specified with the
# `values` keyword argument. See `~proplot.axes.Axes.parametric` for details.
# As shown below, it is also easy to build colorbars from the
# `~matplotlib.collections.LineCollection` returned by
# `~proplot.axes.Axes.parametric`.

# %%
import proplot as plot
import numpy as np
N = 50
cmap = 'IceFire'
values = np.linspace(-N / 2, N / 2, N)
fig, axs = plot.subplots(
    share=0, ncols=2, wratios=(2, 1),
    axwidth='7cm', aspect=(2, 1)
)
axs.format(suptitle='Parametric plots demo')

# Parametric line with smooth gradations
ax = axs[0]
state = np.random.RandomState(51423)
m = ax.parametric(
    (state.rand(N) - 0.5).cumsum(), state.rand(N),
    cmap=cmap, values=values, lw=7, extend='both'
)
ax.format(
    xlabel='xlabel', ylabel='ylabel',
    title='Line with smooth gradations'
)
ax.format(xlim=(-1, 5), ylim=(-0.2, 1.2))
ax.colorbar(m, loc='b', label='parametric coordinate', locator=5)

# Parametric line with stepped gradations
N = 12
ax = axs[1]
values = np.linspace(-N / 2, N / 2, N + 1)
radii = np.linspace(1, 0.2, N + 1)
angles = np.linspace(0, 4 * np.pi, N + 1)
x = radii * np.cos(1.4 * angles)
y = radii * np.sin(1.4 * angles)
m = ax.parametric(x, y, cmap=cmap, values=values, linewidth=15, interp=False)
ax.format(
    xlim=(-1, 1), ylim=(-1, 1), title='Step gradations',
    xlabel='cosine angle', ylabel='sine angle'
)
ax.colorbar(m, loc='b', maxn=10, label=f'parametric coordinate')


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_scatter:
#
# Scatter plots
# -------------
#
# The `~matplotlib.axes.Axes.scatter` method is now wrapped by
# `~proplot.wrappers.scatter_wrapper`, `~proplot.wrappers.cycle_changer`, and
# `~proplot.wrappers.standardize_1d`. This means that
# `~matplotlib.axes.Axes.scatter` now accepts 2D arrays, just like
# `~matplotlib.axes.Axes.plot`. Also, successive calls to
# `~matplotlib.axes.Axes.scatter` now use the property cycler properties
# (e.g.  `color`, `marker`, and `markersize`), and
# `~matplotlib.axes.Axes.scatter` now optionally accepts keywords that look
# like `~matplotlib.axes.Axes.plot` keywords (e.g. `color` instead of `c` and
# `markersize` instead of `s`).
#
# We are also considering supporting 2D array input and property cycle
# iteration for more obscure matplotlib plotting commands like
# `~matplotlib.axes.Axes.stem`, `~matplotlib.axes.Axes.step`,
# `~matplotlib.axes.Axes.vlines`, and `~matplotlib.axes.Axes.hlines`. Stay
# tuned.

# %%
import proplot as plot
import numpy as np
import pandas as pd
fig, axs = plot.subplots(ncols=2, share=1)
state = np.random.RandomState(51423)
x = (state.rand(20) - 0).cumsum()
data = (state.rand(20, 4) - 0.5).cumsum(axis=0)
data = pd.DataFrame(data, columns=pd.Index(['a', 'b', 'c', 'd'], name='label'))

# Scatter plot with property cycler
ax = axs[0]
ax.format(title='Extra prop cycle properties', suptitle='Scatter plot demo')
obj = ax.scatter(
    x, data, legend='ul', cycle='Set2', legend_kw={'ncols': 2},
    cycle_kw={'marker': ['x', 'o', 'x', 'o'], 'markersize': [5, 10, 20, 30]}
)

# Scatter plot with colormap
ax = axs[1]
ax.format(title='Scatter plot with cmap')
data = state.rand(2, 100)
obj = ax.scatter(
    *data, color=data.sum(axis=0), size=state.rand(100), smin=3, smax=30,
    marker='o', cmap='dark red', colorbar='lr', vmin=0, vmax=2,
    colorbar_kw={'label': 'label', 'locator': 0.5}
)
axs.format(xlabel='xlabel', ylabel='ylabel')
