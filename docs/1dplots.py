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
# .. _ug_1dplots:
#
# 1D plotting
# ===========
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
# Input arguments passed to 1D plot commands are now uniformly
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
    title = f'Manual xlim inbounds={inbounds}'
    title += ' (default)' if inbounds else ''
    ax.format(
        xmax=(None if i == 0 else xmax),
        title=('Default xlim' if i == 0 else title),
    )
    ax.plot(x, y, cycle=cycle, inbounds=inbounds)
fig.format(
    xlabel='xlabel',
    ylabel='ylabel',
    suptitle='Default ylim restricted to in-bounds data'
)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_1dintegration:
#
# Pandas and xarray integration
# -----------------------------
#
# The `~proplot.axes.PlotAxes` plotting commands recognize `pandas`_ and
# `xarray`_ data structures. If you omit dependent variable coordinates,
# the plotting commands try to infer them from the `pandas.DataFrame` or
# `xarray.DataArray`. If you did not explicitly set the *x* or *y* axis label
# or :ref:`legend or colorbar <ug_cbars_axes>` label(s), the plotting commands
# try to retrieve them from the `pandas.DataFrame` or `xarray.DataArray`.
# The plotting commands also recognize `pint.Quantity` structures and apply
# unit string labels with formatting specified by :rc:`unitformat`.
#
# These features restore some of the convenience you get with the builtin
# `pandas`_ and `xarray`_ plotting functions. They are also *optional* --
# installation of pandas and xarray are not required to use ProPlot. All of these
# features can be disabled by setting :rcraw:`autoformat` to ``False`` or by
# passing ``autoformat=False`` to any plotting command.
#
# .. note::
#
#    For every plotting command, you can pass a `~xarray.Dataset`, `~pandas.DataFrame`,
#    or `dict` to the `data` keyword with strings as data arguments instead of arrays
#    -- just like matplotlib. For example, ``ax.plot('y', data=dataset)`` and
#    ``ax.plot(y='y', data=dataset)`` are translated to ``ax.plot(dataset['y'])``.
#    This is the preferred input style for most `seaborn`_ plotting commands.
#    Also, if you pass a `pint.Quantity` or `~xarray.DataArray`
#    containing a `pint.Quantity`, ProPlot will automatically call
#    `~pint.UnitRegistry.setup_matplotlib` so that the axes become unit-aware.

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
# <https://matplotlib.org/stable/tutorials/intermediate/color_cycle.html#sphx-glr-tutorials-intermediate-color-cycle-py>`__
# and use different property cycles for different plot elements. You can create and
# apply property cycles on-the-fly using the `cycle` and `cycle_kw` keywords, available
# with most `~proplot.axes.PlotAxes` 1D plot commands. `cycle` and `cycle_kw` are
# passed to the `~proplot.constructor.Cycle` :ref:`constructor function
# <why_constructor>`, and the resulting property cycle is used for the plot. You
# can specify `cycle` once with 2D input data (in which case each column is
# plotted in succession according to the property cycle) or call a plotting
# command multiple times with the same `cycle` argument (the property
# cycle is not reset). You can also disable property cycling with ``cycle=False``,
# ``cycle='none'``, or ``cycle=()`` and re-enable the default property cycle with
# ``cycle=True`` (note that as usual, you can also simply override the property cycle
# with relevant artist keywords like `color`). For more information on property cycles,
# see the :ref:`color cycles section <ug_cycles>` and `this matplotlib tutorial
# <https://matplotlib.org/stable/tutorials/intermediate/color_cycle.html#sphx-glr-tutorials-intermediate-color-cycle-py>`__.

# %%
import proplot as pplt
import numpy as np

# Sample data
M, N = 50, 5
state = np.random.RandomState(51423)
data1 = (state.rand(M, N) - 0.48).cumsum(axis=1).cumsum(axis=0)
data2 = (state.rand(M, N) - 0.48).cumsum(axis=1).cumsum(axis=0) * 1.5
data1 += state.rand(M, N)
data2 += state.rand(M, N)

with pplt.rc.context({'lines.linewidth': 3}):
    # Use property cycle for columns of 2D input data
    fig = pplt.figure(share=False)
    ax = fig.subplot(121)
    ax.format(title='Single plot call')
    ax.plot(
        2 * data1 + data2,
        cycle='black',  # cycle from monochromatic colormap
        cycle_kw={'ls': ('-', '--', '-.', ':')}
    )

    # Use property cycle with successive plot() calls
    ax = fig.subplot(122)
    ax.format(title='Multiple plot calls')
    for i in range(data1.shape[1]):
        ax.plot(data1[:, i], cycle='Reds', cycle_kw={'N': N, 'left': 0.3})
    for i in range(data1.shape[1]):
        ax.plot(data2[:, i], cycle='Blues', cycle_kw={'N': N, 'left': 0.3})
    fig.format(
        xlabel='xlabel', ylabel='ylabel', suptitle='On-the-fly property cycles'
    )


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_lines:
#
# Line plots
# ----------
#
# Line plots can be drawn with `~proplot.axes.PlotAxes.plot` or
# `~proplot.axes.PlotAxes.plotx` (or their aliases, `~proplot.axes.PlotAxes.line`
# or `~proplot.axes.PlotAxes.linex`). For the ``x`` commands, positional
# arguments are interpreted as *x* coordinates or (*y*, *x*) pairs. This is analogous
# to `~proplot.axes.PlotAxes.barh` and `~proplot.axes.PlotAxes.fill_betweenx`.
# Also, the default *x* bounds for lines drawn with `~proplot.axes.PlotAxes.plot`
# and *y* bounds for lines drawn with `~proplot.axes.PlotAxes.plotx` are now
# "sticky", i.e. there is no padding between the lines and axes edges by default.
#
# Step and stem plots can be drawn with `~proplot.axes.PlotAxes.step`,
# `~proplot.axes.PlotAxes.stepx`, `~proplot.axes.PlotAxes.stem`, and
# `~proplot.axes.PlotAxes.stemx`. Plots of parallel vertical and horizontal
# lines can be drawn with `~proplot.axes.PlotAxes.vlines` and
# `~proplot.axes.PlotAxes.hlines`. You can have different colors for "negative" and
# "positive" lines using ``negpos=True`` (see :ref:`below <ug_negpos>` for details).

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


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_bar:
#
# Bar plots and area plots
# ------------------------
#
# The `~proplot.axes.PlotAxes.bar` and `~proplot.axes.PlotAxes.barh` commands
# apply default *x* or *y* coordinates if you failed to provide them explicitly
# and can *group* or *stack* successive columns of data if you pass 2D arrays instead
# of 1D arrays -- just like `pandas`_. The widths of bars are expressed in step
# size-relative units by default, but matplotlib's behavior can be restored by
# passing ``absolute_width=True``. When bars are grouped, their widths are
# positions are adjusted according to the number of bars in the group.
#
# The `~proplot.axes.PlotAxes.fill_between` and `~proplot.axes.PlotAxes.fill_betweenx`
# have the new shorthands `~proplot.axes.PlotAxes.area`
# and `~proplot.axes.PlotAxes.areax`. Similar to `~proplot.axes.PlotAxes.bar` and
# `~proplot.axes.PlotAxes.barh`, they apply default *x* coordinates if you failed
# to provide them explicitly, and can *overlay* or *stack* successive columns of
# data if you pass 2D arrays instead of 1D arrays -- just like `pandas`_. Also, the
# default *x* bounds for shading drawn with `~proplot.axes.PlotAxes.area` and *y*
# bounds for shading drawn with `~proplot.axes.PlotAxes.areax` is now "sticky",
# i.e. there is no padding between the shading and axes edges by default.

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
pplt.rc.reset()

# %%
import proplot as pplt
import numpy as np

# Sample data
state = np.random.RandomState(51423)
data = state.rand(5, 3).cumsum(axis=0)
cycle = ('gray3', 'gray5', 'gray7')

# Figure
pplt.rc.abc = 'a.'
pplt.rc.titleloc = 'l'
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
pplt.rc.reset()

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_negpos:
#
# Negative and positive colors
# ----------------------------
#
# You can use different colors for "negative" and
# "positive" data by passing ``negpos=True`` to any of the
# `~proplot.axes.PlotAxes.fill_between`, `~proplot.axes.PlotAxes.fill_betweenx`
# (shorthands `~proplot.axes.PlotAxes.area`, `~proplot.axes.PlotAxes.areax`),
# `~proplot.axes.PlotAxes.vlines`, `~proplot.axes.PlotAxes.hlines`,
# `~proplot.axes.PlotAxes.bar`, or `~proplot.axes.PlotAxes.barh` plotting commands.
# The default negative and positive colors are controlled with :rcraw:`negcolor`
# and :rcraw:`poscolor` but the colors can be modified for particular plots by
# passing ``negcolor=color`` and ``poscolor=color`` to the plotting commands.

# %%
import proplot as pplt
import numpy as np

# Sample data
state = np.random.RandomState(51423)
data = 4 * (state.rand(40) - 0.5)

# Figure
pplt.rc.abc = 'a.'
pplt.rc.titleloc = 'l'
fig, axs = pplt.subplots(nrows=3, refaspect=2, figwidth=5)
axs.format(
    xmargin=0, xlabel='xlabel', ylabel='ylabel', grid=True,
    suptitle='Positive and negative colors demo',
)
for ax in axs:
    ax.axhline(0, color='k', linewidth=1)  # zero line

# Line plot
ax = axs[0]
ax.vlines(data, linewidth=3, negpos=True)
ax.format(title='Line plot')

# Bar plot
ax = axs[1]
ax.bar(data, width=1, negpos=True, edgecolor='k')
ax.format(title='Bar plot')

# Area plot
ax = axs[2]
ax.area(data, negpos=True, lw=0.5, edgecolor='k')
ax.format(title='Area plot')

# Reset title styles changed above
pplt.rc.reset()
