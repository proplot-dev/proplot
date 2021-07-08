# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.3
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
# .. _ug_1dplots:
#
# Plotting 1D data
# ================
#
# ProPlot :ref:`adds new features <why_plotting>` to various `~matplotlib.axes.Axes`
# plotting methods using a set of "wrapper" functions. When a plotting method like
# `~matplotlib.axes.Axes.plot` is "wrapped" by one of these functions, it accepts
# the same parameters as the wrapper. These additions are a strict *superset* of
# matplotlib -- if you are not interested, you can use matplotlib's plotting methods
# just like you always have. This section documents the features added by wrapper
# functions to 1D plotting commands like `~matplotlib.axes.Axes.plot`,
# `~matplotlib.axes.Axes.scatter`, `~matplotlib.axes.Axes.bar`, and
# `~matplotlib.axes.Axes.barh`.


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_1dstd:
#
# Standardized arguments
# ----------------------
#
# The `~proplot.axes.standardize_1d` wrapper standardizes
# positional arguments across all 1D plotting methods.
# `~proplot.axes.standardize_1d` allows you to optionally omit *x*
# coordinates, in which case they are inferred from the data. It also permits
# passing 2D *y* coordinate arrays to any plotting method, in which case the
# plotting method is called for each column of the array.

# %%
import proplot as plot
import numpy as np

N = 5
state = np.random.RandomState(51423)
with plot.rc.context({'axes.prop_cycle': plot.Cycle('Grays', N=N, left=0.3)}):
    # Figure and sample data
    x = np.linspace(-5, 5, N)
    y = state.rand(N, 5)
    fig, axs = plot.subplots(ncols=2, share=False)
    axs.format(xlabel='xlabel', ylabel='ylabel')
    axs.format(suptitle='Standardized arguments demonstration')

    # Plot by passing both x and y coordinates
    ax = axs[0]
    ax.area(x, -1 * y / N, stack=True)
    ax.bar(x, y, linewidth=0, alpha=1, width=0.8)
    ax.plot(x, y + 1, linewidth=2)
    ax.scatter(x, y + 2, marker='s', markersize=5**2)
    ax.format(title='Manual x coordinates')

    # Plot by passing just y coordinates
    # Default x coordinates are inferred from DataFrame,
    # inferred from DataArray, or set to np.arange(0, y.shape[0])
    ax = axs[1]
    ax.area(-1 * y / N, stack=True)
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
# The `~proplot.axes.standardize_1d` wrapper integrates 1D plotting
# methods with pandas `~pandas.DataFrame`\ s and xarray `~xarray.DataArray`\ s.
# When you pass a DataFrame or DataArray to any plotting command, the x-axis
# label, y-axis label, legend label, colorbar label, and/or title are
# configured from the metadata. This restores some of the convenience you get
# with the builtin `pandas`_ and `xarray`_ plotting functions. This feature is
# *optional*. Installation of pandas and xarray are not required, and
# it can be disabled by setting :rcraw:`autoformat` to ``False``.

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
    (np.cos(np.linspace(0, 2 * np.pi, 20))**4)[:, None] + state.rand(20, 5)**2
)
ts = pd.date_range('1/1/2000', periods=20)
df = pd.DataFrame(data, index=ts, columns=['foo', 'bar', 'baz', 'zap', 'baf'])
df.name = 'data'
df.index.name = 'date'
df.columns.name = 'category'

# %%
import proplot as plot
fig, axs = plot.subplots(ncols=2, refwidth=2.2, share=0)
axs.format(suptitle='Automatic subplot formatting')

# Plot DataArray
cycle = plot.Cycle('dark blue', fade=90, space='hpl', N=da.shape[1])
axs[0].scatter(da, cycle=cycle, lw=3, colorbar='ul', colorbar_kw={'locator': 20})

# Plot Dataframe
cycle = plot.Cycle('dark green', fade=90, space='hpl', N=df.shape[1])
axs[1].plot(df, cycle=cycle, lw=3, legend='uc')


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_apply_cycle:
#
# Property cycles
# ---------------
#
# It is often useful to create on-the-fly `property cycles
# <https://matplotlib.org/tutorials/intermediate/color_cycle.html#sphx-glr-tutorials-intermediate-color-cycle-py>`__
# and use different property cycles for different plot elements. You can create and
# apply property cycles on-the-fly using the `cycle` and `cycle_kw` arguments, available
# with any plotting method wrapped by `~proplot.axes.apply_cycle`. `cycle` and
# `cycle_kw` are passed to the `~proplot.constructor.Cycle` :ref:`constructor function
# <why_constructor>`, and the resulting property cycle is used for the plot. You can
# specify `cycle` once with 2D input data (in which case each column is plotted in
# succession according to the property cycle) or call a plotting command multiple times
# with the same `cycle` argument each time (the property cycle is not reset). For more
# information on property cycles, see the :ref:`color cycles section <ug_cycles>` and
# `this matplotlib tutorial
# <https://matplotlib.org/tutorials/intermediate/color_cycle.html#sphx-glr-tutorials-intermediate-color-cycle-py>`__.

# %%
import proplot as plot
import numpy as np

# Sample data
N = 4
state = np.random.RandomState(51423)
data1 = state.rand(6, N)
data2 = state.rand(6, N) * 1.5

with plot.rc.context({'lines.linewidth': 3}):
    # Figure
    fig, axs = plot.subplots(ncols=2)
    axs.format(xlabel='xlabel', ylabel='ylabel', suptitle='Local property cycles demo')

    # Use property cycle for columns of 2D input data
    axs[0].plot(
        data1 * data2,
        cycle='black',
        cycle_kw={'linestyle': ('-', '--', '-.', ':')}
    )

    # Use property cycle with successive plot() calls
    for i in range(data1.shape[1]):
        axs[1].plot(data1[:, i], cycle='Reds', cycle_kw={'N': N, 'left': 0.3})
    for i in range(data1.shape[1]):
        axs[1].plot(data2[:, i], cycle='Blues', cycle_kw={'N': N, 'left': 0.3})


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_errorbars:
#
# Shading and error bars
# ----------------------
#
# The `~proplot.axes.indicate_error` wrapper lets you draw error bars
# and error shading on-the-fly by passing certain keyword arguments to
# `~matplotlib.axes.Axes.plot`, `~matplotlib.axes.Axes.scatter`,
# `~matplotlib.axes.Axes.bar`, or `~matplotlib.axes.Axes.barh`.
#
# If you pass 2D arrays to these methods with ``mean=True`` or ``median=True``,
# the means or medians of each column are drawn as points, lines, or bars, and
# *error bars* or *shading* is drawn to represent the spread of the distribution
# for each column. You can also specify the error bounds *manually* with the
# `bardata`, `boxdata`, `shadedata`, and `fadedata` keywords.
# `~proplot.axes.indicate_error` can draw thin error bars with optional whiskers,
# thick "boxes" overlayed on top of these bars (think of these as miniature boxplots),
# and up to 2 layers of shading. See `~proplot.axes.indicate_error` for details.


# %%
import proplot as plot
import numpy as np
import pandas as pd
plot.rc['title.loc'] = 'uc'

# Sample data
state = np.random.RandomState(51423)
data = state.rand(20, 8).cumsum(axis=0).cumsum(axis=1)[:, ::-1]
data = data + 20 * state.normal(size=(20, 8)) + 30
data = pd.DataFrame(data, columns=np.arange(0, 16, 2))
data.name = 'variable'

# Figure
fig, axs = plot.subplots(
    nrows=3, refaspect=1.5, refwidth=4,
    share=0, hratios=(2, 1, 1)
)
axs.format(suptitle='Indicating error bounds')
axs[1:].format(xlabel='column number', xticks=1, xgrid=False)

# Medians and percentile ranges
ax = axs[0]
obj = ax.barh(
    data, color='light red', legend=True,
    median=True, barpctile=90, boxpctile=True,
    # median=True, barpctile=(5, 95), boxpctile=(25, 75)  # equivalent
)
ax.format(title='Column statistics')
ax.format(ylabel='column number', title='Bar plot', ygrid=False)

# Means and standard deviation range
ax = axs[1]
ax.scatter(
    data, color='denim', marker='x', markersize=8**2, linewidth=0.8, legend='ll',
    label='mean', shadelabel=True,
    mean=True, shadestd=1,
    # mean=True, shadestd=(-1, 1)  # equivalent
)
ax.format(title='Marker plot')

# User-defined error bars
ax = axs[2]
means = data.mean(axis=0)
means.name = data.name
shadedata = np.percentile(data, (25, 75), axis=0)  # dark shading
fadedata = np.percentile(data, (5, 95), axis=0)  # light shading
ax.plot(
    means,
    shadedata=shadedata, fadedata=fadedata,
    label='mean', shadelabel='50% CI', fadelabel='90% CI',
    color='ocean blue', barzorder=0, boxmarker=False, legend='ll',
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
# are wrapped by `~proplot.axes.bar_extras`,
# `~proplot.axes.apply_cycle`, and `~proplot.axes.standardize_1d`.
# This means that `~matplotlib.axes.Axes.bar` and `~matplotlib.axes.Axes.barh` employ
# default *x* or *y* coordinates if you failed to provide them explicitly.
# You can now *group* or *stack* columns of data by passing 2D arrays to
# `~matplotlib.axes.Axes.bar` or `~matplotlib.axes.Axes.barh`, just like in
# `pandas`_. You can also use different colors for "negative" and "positive"
# bars by passing ``negpos=True`` (the default colors are :rc:`negcolor`
# and :rc:`poscolor`).
#
# The `~matplotlib.axes.Axes.fill_between` and `~matplotlib.axes.Axes.fill_betweenx`
# commands are wrapped by `~proplot.axes.fill_between_extras` and
# `~proplot.axes.fill_betweenx_extras`. They also now have the optional shorthands
# `~proplot.axes.Axes.area` and `~proplot.axes.Axes.areax`.
# You can now *stack* or *overlay* columns of data by passing 2D arrays to
# to these commands, just like in `pandas`_. You can also draw area plots that
# change color when the fill boundaries cross each other by passing ``negpos=True``
# (the default colors are :rc:`negcolor` and :rc:`poscolor`). The most common
# use case for this is highlighting *negative* and *positive* areas with different
# colors, as shown below.

# %%
import proplot as plot
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
plot.rc.titleloc = 'uc'
fig, axs = plot.subplots(nrows=2, refaspect=2, refwidth=4.8, share=0, hratios=(3, 2))

# Side-by-side bars
ax = axs[0]
obj = ax.bar(
    data, cycle='Reds', colorbar='ul',
    edgecolor='red9', colorbar_kw={'frameon': False}
)
ax.format(
    xlocator=1, xminorlocator=0.5, ytickminor=False,
    title='Side-by-side', suptitle='Bar plot demo'
)

# Stacked bars
ax = axs[1]
obj = ax.barh(
    data.iloc[::-1, :], cycle='Blues',
    legend='ur', edgecolor='blue9', stack=True
)
ax.format(title='Stacked')
axs.format(grid=False)
plot.rc.reset()

# %%
import proplot as plot
import numpy as np

# Sample data
state = np.random.RandomState(51423)
data = state.rand(5, 3).cumsum(axis=0)
cycle = ('gray3', 'gray5', 'gray7')

# Figure
fig, axs = plot.subplots(ncols=2, refwidth=2.3, share=0)
axs.format(grid=False, xlabel='xlabel', ylabel='ylabel', suptitle='Area plot demo')

# Overlaid area patches
ax = axs[0]
ax.area(
    np.arange(5), data, data + state.rand(5)[:, None], cycle=cycle, alpha=0.7,
    legend='uc', legend_kw={'center': True, 'ncols': 2, 'labels': ['z', 'y', 'qqqq']},
)
ax.format(title='Fill between columns')

# Stacked area patches
ax = axs[1]
ax.area(
    np.arange(5), data, stack=True, cycle=cycle, alpha=0.8,
    legend='ul', legend_kw={'center': True, 'ncols': 2, 'labels': ['z', 'y', 'qqqq']},
)
ax.format(title='Stack between columns')

# %%
import proplot as plot
import numpy as np

# Sample data
state = np.random.RandomState(51423)
data = 4 * (state.rand(40) - 0.5)

# Figure
fig, axs = plot.subplots(nrows=2, refaspect=2, figwidth=5)
axs.format(
    xmargin=0, xlabel='xlabel', ylabel='ylabel', grid=True,
    suptitle='Positive and negative colors demo',
)
axs.axhline(0, color='k', linewidth=1)  # zero line

# Bar plot
axs[0].bar(data, width=1, edgecolor='none', negpos=True)
axs[0].format(title='Bar plot')

# Area plot
axs[1].area(data, negpos=True)
axs[1].format(title='Area plot')


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_boxplots:
#
# Box plots and violin plots
# --------------------------
#
# The `~matplotlib.axes.Axes.boxplot` and `~matplotlib.axes.Axes.violinplot`
# commands are wrapped by `~proplot.axes.boxplot_extras`,
# `~proplot.axes.violinplot_extras`, `~proplot.axes.apply_cycle`,
# and `~proplot.axes.standardize_1d`. They also now have the optional shorthands
# `~proplot.axes.Axes.boxes` and `~proplot.axes.Axes.violins`. The wrappers
# apply aesthetically pleasing default settings and permit configuration using
# keyword arguments like ``color``, ``boxcolor``, and ``fillcolor``. They also
# automatically apply axis labels based on the `~pandas.DataFrame` or
# `~xarray.DataArray` column labels or the input *x* coordinate labels.

# %%
import proplot as plot
import numpy as np
import pandas as pd

# Sample data
N = 500
state = np.random.RandomState(51423)
data = state.normal(size=(N, 5)) + 2 * (state.rand(N, 5) - 0.5) * np.arange(5)
data = pd.DataFrame(
    data, columns=pd.Index(['a', 'b', 'c', 'd', 'e'], name='xlabel')
)

# Figure
fig, axs = plot.subplots([[1, 1, 2, 2], [0, 3, 3, 0]])
axs.format(grid=False, suptitle='Boxes and violins demo')

# Box plots
ax = axs[0]
obj1 = ax.boxplot(
    data, means=True, meancolor='red', marker='x', fillcolor='gray5',
)
ax.format(title='Box plots', titleloc='uc')

# Violin plots
ax = axs[1]
obj2 = ax.violinplot(
    data, fillcolor='gray7',
    points=500, bw_method=0.3, means=True,
)
ax.format(title='Violin plots', titleloc='uc')

# Boxes with different colors
ax = axs[2]
data = state.rand(100, 7)
colors = plot.Colors('pastel2')  # list of colors from the cycle
ax.boxplot(data, fillcolor=colors, orientation='horizontal')
ax.format(title='Multiple colors', titleloc='uc', ymargin=0.15)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_lines:
#
# Line plots
# ----------
#
# The `~matplotlib.axes.Axes.plot` command is wrapped by
# `~proplot.axes.apply_cycle` and `~proplot.axes.standardize_1d`.
# But in general, its behavior is the same -- ProPlot simply tries to expand
# the flexibility of this command to the rest of the 1D plotting commands.
# The new `~proplot.axes.Axes.plotx` command can be used just like
# `~matplotlib.axes.Axes.plotx`, except a single argument is interpreted
# as *x* coordinates (with default *y* coordinates inferred from the data),
# and multiple arguments are interpreted as *y* and *x* coordinates (in that order).
# This is analogous to `~matplotlib.axes.Axes.barh` and
# `~matplotlib.axes.Axes.fill_betweenx`.
#
# As with the other 1D plotting commands, `~matplotlib.axes.Axes.step`,
# `~matplotlib.axes.Axes.hlines`, `~matplotlib.axes.Axes.vlines`, and
# `~matplotlib.axes.Axes.stem` are wrapped by `~proplot.axes.standardize_1d`.
# `~proplot.axes.Axes.step` now use the property cycle, just like
# `~matplotlib.axes.Axes.plot`. `~matplotlib.axes.Axes.vlines` and
# `~matplotlib.axes.Axes.hlines` are also wrapped by `~proplot.axes.vlines_extras`
# and `~proplot.axes.hlines_extras`, which can use
# different colors for "negative" and "positive" lines using ``negpos=True``
# (the default colors are :rc:`negcolor` and :rc:`poscolor`).

# %%
import proplot as plot
import numpy as np
state = np.random.RandomState(51423)
fig, axs = plot.subplots(ncols=2, nrows=2, share=0)
axs.format(suptitle='Line plots demo', xlabel='xlabel', ylabel='ylabel')

# Step
ax = axs[0]
data = state.rand(20, 4).cumsum(axis=1).cumsum(axis=0)
cycle = ('blue7', 'gray5', 'red7', 'gray5')
ax.step(data, cycle=cycle, labels=list('ABCD'), legend='ul', legend_kw={'ncol': 2})
ax.format(title='Step plot')

# Stems
ax = axs[1]
data = state.rand(20)
ax.stem(data, linefmt='k-')
ax.format(title='Stem plot')

# Vertical lines
gray = 'gray7'
data = state.rand(20) - 0.5
ax = axs[2]
ax.area(data, color=gray, alpha=0.2)
ax.vlines(data, negpos=True, linewidth=2)
ax.format(title='Vertical lines')

# Horizontal lines
ax = axs[3]
ax.areax(data, color=gray, alpha=0.2)
ax.hlines(data, negpos=True, linewidth=2)
ax.format(title='Horizontal lines')

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_scatter:
#
# Scatter plots
# -------------
#
# The `~matplotlib.axes.Axes.scatter` command is wrapped by
# `~proplot.axes.scatter_extras`, `~proplot.axes.apply_cycle`, and
# `~proplot.axes.standardize_1d`. This means that
# `~matplotlib.axes.Axes.scatter` now accepts 2D *y* coordinates and permits
# omitting *x* coordinates, just like `~matplotlib.axes.Axes.plot`.
# `~matplotlib.axes.Axes.scatter` now also accepts keywords that look like
# `~matplotlib.axes.Axes.plot` keywords (e.g., `color` instead of `c` and
# `markersize` instead of `s`). This way, `~matplotlib.axes.Axes.scatter` can
# optionally be used simply to "plot markers, not lines" without changing the
# input arguments relative to `~matplotlib.axes.Axes.plot`.
#
# Just like `~matplotlib.axes.Axes.plot`, the property cycle is used
# with `~matplotlib.axes.Axes.scatter` plots by default. It can be changed
# using the `cycle` keyword argument, and it can include properties like `marker`
# and `markersize`. The colormap `cmap` and normalizer `norm` used with the
# optional `c` color array are now passed through the `~proplot.constructor.Colormap`
# and `~proplot.constructor.Norm` constructor functions, and the the `s` marker
# size array can now be conveniently scaled using the arguments `smin` and `smax`
# (analogous to `vmin` and `vmax` used for colors).

# %%
import proplot as plot
import numpy as np
import pandas as pd

# Sample data
state = np.random.RandomState(51423)
x = (state.rand(20) - 0).cumsum()
data = (state.rand(20, 4) - 0.5).cumsum(axis=0)
data = pd.DataFrame(data, columns=pd.Index(['a', 'b', 'c', 'd'], name='label'))

# Figure
fig, axs = plot.subplots(ncols=2, share=1)
axs.format(suptitle='Scatter plot demo')

# Scatter plot with property cycler
ax = axs[0]
ax.format(title='With property cycle')
obj = ax.scatter(
    x, data, legend='ul', cycle='Set2', legend_kw={'ncols': 2},
    cycle_kw={'marker': ['x', 'o', 'x', 'o'], 'markersize': [5, 10, 20, 30]}
)

# Scatter plot with colormap
ax = axs[1]
ax.format(title='With colormap')
data = state.rand(2, 100)
obj = ax.scatter(
    *data, color=data.sum(axis=0), size=state.rand(100), smin=3, smax=30,
    marker='o', cmap='dark red', cmap_kw={'fade': 90}, vmin=0, vmax=2,
    colorbar='lr', colorbar_kw={'label': 'label', 'locator': 0.5},
)
axs.format(xlabel='xlabel', ylabel='ylabel')

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_parametric:
#
# Parametric plots
# ----------------
#
# To make "parametric" plots, use the new `~proplot.axes.Axes.parametric`
# command. Parametric plots are `~matplotlib.collections.LineCollection`\ s that
# map individual line segments to individual colors, where each segment represents a
# "parametric" coordinate (e.g., time). The parametric coordinates are specified with
# the `values` keyword argument. See `~proplot.axes.Axes.parametric` for details. As
# shown below, it is also easy to build colorbars from the
# `~matplotlib.collections.LineCollection` returned by `~proplot.axes.Axes.parametric`.

# %%
import proplot as plot
import numpy as np
fig, axs = plot.subplots(
    share=0, ncols=2, wratios=(2, 1),
    figwidth='16cm', refaspect=(2, 1)
)
axs.format(suptitle='Parametric plots demo')
cmap = 'IceFire'

# Sample data
state = np.random.RandomState(51423)
N = 50
x = (state.rand(N) - 0.52).cumsum()
y = state.rand(N)
c = np.linspace(-N / 2, N / 2, N)  # color values

# Parametric line with smooth gradations
ax = axs[0]
m = ax.parametric(
    x, y, c, cmap=cmap, lw=7, interp=5, capstyle='round', joinstyle='round'
)
ax.format(xlabel='xlabel', ylabel='ylabel', title='Line with smooth gradations')
ax.colorbar(m, loc='b', label='parametric coordinate', locator=5)

# Sample data
N = 12
radii = np.linspace(1, 0.2, N + 1)
angles = np.linspace(0, 4 * np.pi, N + 1)
x = radii * np.cos(1.4 * angles)
y = radii * np.sin(1.4 * angles)
c = np.linspace(-N / 2, N / 2, N + 1)

# Parametric line with stepped gradations
ax = axs[1]
m = ax.parametric(x, y, c, cmap=cmap, lw=15)
ax.format(
    xlim=(-1, 1), ylim=(-1, 1), title='Step gradations',
    xlabel='cosine angle', ylabel='sine angle'
)
ax.colorbar(m, loc='b', maxn=10, label='parametric coordinate')
