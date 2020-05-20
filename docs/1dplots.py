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
# Plotting 1D data
# ================
#
# ProPlot adds new features to various `~matplotlib.axes.Axes` plotting
# methods using a set of wrapper functions. When a plotting method like
# `~matplotlib.axes.Axes.plot` is "wrapped" by one of these functions, it
# accepts the same parameters as the wrapper. These features are a strict
# *superset* of the matplotlib API.
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
# `~proplot.axes.cycle_changer` adds the `cycle` and `cycle_kw` to the 1D
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
# The `~proplot.axes.standardize_1d` wrapper is used to standardize
# positional arguments across all 1D plotting methods.
# `~proplot.axes.standardize_1d` allows you to optionally omit *x*
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
    axs.format(xlabel='xlabel', ylabel='ylabel')
    axs.format(suptitle='Standardized arguments demonstration')

    # Plot by passing both x and y coordinates
    ax = axs[0]
    ax.area(x, -1 * y / N, stacked=True)
    ax.bar(x, y, linewidth=0, alpha=1, width=0.8)
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
# The `~proplot.axes.standardize_1d` wrapper integrates 1D plotting
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
# Shading and error bars
# ----------------------
#
# The `~proplot.axes.indicate_error` wrapper lets you draw error bars
# and error shading on-the-fly by passing certain keyword arguments to
# `~matplotlib.axes.Axes.plot`, `~matplotlib.axes.Axes.scatter`,
# `~matplotlib.axes.Axes.bar`, or `~matplotlib.axes.Axes.barh`.
#
# If you pass 2D arrays to these methods with ``means=True`` or
# ``medians=True``, the means or medians of each column are drawn as points,
# lines, or bars, and *error bars* or *shading* is drawn to represent the spread
# of the distribution for each column. You can also specify the error bounds
# *manually* with the `bardata`, `boxdata`, `shadedata`, and `fadedata` keywords.
# `~proplot.axes.indicate_error` can draw thin error bars with optional whiskers,
# thick "boxes" overlayed on top of these bars (think of this as a miniature boxplot),
# and up to 2 regions of shading. See `~proplot.axes.indicate_error` for details.


# %%
import proplot as plot
import numpy as np
import pandas as pd
plot.rc['title.loc'] = 'uc'

# Generate sample data
state = np.random.RandomState(51423)
data = state.rand(20, 8).cumsum(axis=0).cumsum(axis=1)[:, ::-1]
data = data + 20 * state.normal(size=(20, 8)) + 30
data = pd.DataFrame(data, columns=np.arange(0, 16, 2))
data.name = 'variable'

# Generate figure
fig, axs = plot.subplots(
    nrows=3, aspect=1.5, axwidth=4,
    share=0, hratios=(2, 1, 1)
)
axs.format(suptitle='Indicating error bounds with various plotting commands')
axs[1:].format(xlabel='column number', xticks=1, xgrid=False)

# Automatically calculate medians and display default percentile range
ax = axs[0]
obj = ax.barh(
    data, color='light red', legend=True,
    medians=True, boxpctiles=True, barpctiles=(5, 95),
)
ax.format(title='Column statistics')
ax.format(ylabel='column number', title='Bar plot', ygrid=False)

# Automatically calculate means and display requested standard deviation range
ax = axs[1]
ax.scatter(
    data, color='denim', marker='x', markersize=8**2, linewidth=0.8,
    means=True, shadestds=(-1, 1), legend='ll',
)
ax.format(title='Scatter plot')

# Manually supply error bar data and legend labels
ax = axs[2]
means = data.mean(axis=0)
means.name = data.name
shadedata = np.percentile(data, (25, 75), axis=0)  # dark shading
fadedata = np.percentile(data, (5, 95), axis=0)  # light shading
ax.plot(
    means, shadedata=shadedata, fadedata=fadedata,
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
# are wrapped by `~proplot.axes.bar_wrapper`,
# `~proplot.axes.cycle_changer`, and `~proplot.axes.standardize_1d`.
# You can now *group* or *stack* columns of data by passing 2D arrays to
# `~matplotlib.axes.Axes.bar` or `~matplotlib.axes.Axes.barh`, just like in
# `pandas`, or use different colors for negative and positive bars by
# passing ``negpos=True``. Also, `~matplotlib.axes.Axes.bar` and
# `~matplotlib.axes.Axes.barh` now employ "default" *x* coordinates if you
# failed to provide them explicitly.
#
# To make filled "area" plots, use the new `~proplot.axes.Axes.area` and
# `~proplot.axes.Axes.areax` methods. These are alises for
# `~matplotlib.axes.Axes.fill_between` and
# `~matplotlib.axes.Axes.fill_betweenx`, which are wrapped by
# `~proplot.axes.fill_between_wrapper` and
# `~proplot.axes.fill_betweenx_wrapper`. You can now *stack* or *overlay*
# columns of data by passing 2D arrays to `~proplot.axes.Axes.area` and
# `~proplot.axes.Axes.areax`, just like in `pandas`. You can also now draw
# area plots that change color when the fill boundaries cross each other by
# passing ``negpos=True`` to `~matplotlib.axes.Axes.fill_between`. The most
# common use case for this is highlighting negative and positive areas with
# different colors, as shown below.

# %%
import proplot as plot
import numpy as np
import pandas as pd
plot.rc.titleloc = 'uc'
fig, axs = plot.subplots(nrows=2, aspect=2, axwidth=4.8, share=0, hratios=(3, 2))
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
    title='Side-by-side', suptitle='Bar plot demo'
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
fig, axs = plot.subplots(array=[[1, 2], [3, 3]], hratios=(1, 1.5), axwidth=2.3, share=0)
axs.format(grid=False, xlabel='xlabel', ylabel='ylabel', suptitle='Area plot demo')
state = np.random.RandomState(51423)
data = state.rand(5, 3).cumsum(axis=0)
cycle = ('gray3', 'gray5', 'gray7')

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
    np.arange(5), data, stacked=True, cycle=cycle, alpha=0.8,
    legend='ul', legend_kw={'center': True, 'ncols': 2, 'labels': ['z', 'y', 'qqqq']},
)
ax.format(title='Stack between columns')

# Positive and negative color bars and area patches
ax = axs[2]
data = 4 * (state.rand(20) - 0.5)
ax.bar(data, bottom=-2, width=1, edgecolor='none', negpos=True)
ax.area(data + 2, y2=2, negpos=True)
for offset in (-2, 2):
    ax.axhline(offset, color='k', linewidth=1, linestyle='--')
ax.format(
    xmargin=0, xlabel='xlabel', ylabel='ylabel', grid=True,
    title='Positive and negative colors demo', titleweight='bold',
)
plot.rc.reset()


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_boxplots:
#
# Box plots and violin plots
# --------------------------
#
# The `~matplotlib.axes.Axes.boxplot` and `~matplotlib.axes.Axes.violinplot`
# methods are now wrapped with `~proplot.axes.boxplot_wrapper`,
# `~proplot.axes.violinplot_wrapper`, `~proplot.axes.cycle_changer`,
# and `~proplot.axes.standardize_1d`. These wrappers add some useful
# options and apply aesthetically pleasing default settings. They also
# automatically apply axis labels based on the `~pandas.DataFrame` column
# labels or the input *x* coordinate labels.

# %%
import proplot as plot
import numpy as np
import pandas as pd

# Generate sample data
N = 500
state = np.random.RandomState(51423)
data = state.normal(size=(N, 5)) + 2 * (state.rand(N, 5) - 0.5) * np.arange(5)
data = pd.DataFrame(
    data,
    columns=pd.Index(['a', 'b', 'c', 'd', 'e'], name='xlabel')
)

# Generate figure
fig, axs = plot.subplots(ncols=2, axwidth=2.5)
axs.format(grid=False, suptitle='Boxes and violins demo')

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
fig, axs = plot.subplots(
    share=0, ncols=2, wratios=(2, 1),
    width='16cm', aspect=(2, 1)
)
axs.format(suptitle='Parametric plots demo')
cmap = 'IceFire'

# Parametric line with smooth gradations
ax = axs[0]
state = np.random.RandomState(51423)
N = 50
x = (state.rand(N) - 0.52).cumsum()
y = state.rand(N)
c = np.linspace(-N / 2, N / 2, N)  # color values
m = ax.parametric(
    x, y, c, cmap=cmap, lw=7, interp=5, capstyle='round', joinstyle='round'
)
ax.format(xlabel='xlabel', ylabel='ylabel', title='Line with smooth gradations')
ax.colorbar(m, loc='b', label='parametric coordinate', locator=5)

# Parametric line with stepped gradations
N = 12
ax = axs[1]
radii = np.linspace(1, 0.2, N + 1)
angles = np.linspace(0, 4 * np.pi, N + 1)
x = radii * np.cos(1.4 * angles)
y = radii * np.sin(1.4 * angles)
c = np.linspace(-N / 2, N / 2, N + 1)
m = ax.parametric(x, y, c, cmap=cmap, lw=15)
ax.format(
    xlim=(-1, 1), ylim=(-1, 1), title='Step gradations',
    xlabel='cosine angle', ylabel='sine angle'
)
ax.colorbar(m, loc='b', maxn=10, label='parametric coordinate')


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_scatter:
#
# Other plotting methods
# ----------------------
#
# The `~matplotlib.axes.Axes.scatter` method is now wrapped by
# `~proplot.axes.scatter_wrapper`, `~proplot.axes.cycle_changer`, and
# `~proplot.axes.standardize_1d`. This means that
# `~matplotlib.axes.Axes.scatter` now accepts 2D arrays, just like
# `~matplotlib.axes.Axes.plot`. Also, successive calls to
# `~matplotlib.axes.Axes.scatter` now use the property cycler properties
# (e.g.  `color`, `marker`, and `markersize`), and
# `~matplotlib.axes.Axes.scatter` now optionally accepts keywords that look
# like `~matplotlib.axes.Axes.plot` keywords (e.g. `color` instead of `c` and
# `markersize` instead of `s`).
#
# ProPlot also supports property cycling for `~proplot.axes.Axes.step` plots
# and wraps the `~matplotlib.axes.Axes.vlines` and `~matplotlib.axes.Axes.hlines`
# methods with `~proplot.axes.vlines_wrapper` and `~proplot.axes.hlines_wrapper`,
# which adds the ability to use different colors for "negative" and "positive" lines.

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
ax.format(suptitle='Scatter plot demo', title='Extra prop cycle properties')
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
