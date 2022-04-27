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
# .. _ug_basics:
#
# The basics
# ==========

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_intro:
#
# Creating figures
# ----------------
#
# Proplot works by `subclassing
# <https://docs.python.org/3/tutorial/classes.html#inheritance>`__
# three fundamental matplotlib classes: `proplot.figure.Figure` replaces
# `matplotlib.figure.Figure`, `proplot.axes.Axes` replaces `matplotlib.axes.Axes`,
# and `proplot.gridspec.GridSpec` replaces `matplotlib.gridspec.GridSpec`
# (see this `tutorial
# <https://matplotlib.org/stable/tutorials/intermediate/gridspec.html>`__
# for more on gridspecs).
#
# To make plots with these classes, you must start with the top-level commands
# `~proplot.ui.figure`, `~proplot.ui.subplot`, or `~proplot.ui.subplots`. These are
# modeled after the `~matplotlib.pyplot` commands of the same name. As in
# `~matplotlib.pyplot`, `~proplot.ui.subplot` creates a figure and a single
# subplot, `~proplot.ui.subplots` creates a figure and a grid of subplots, and
# `~proplot.ui.figure` creates an empty figure that can be subsequently filled
# with subplots. A minimal example with just one subplot is shown below.
#
# %% [raw] raw_mimetype="text/restructuredtext"
# .. note::
#
#    Proplot changes the default :rcraw:`figure.facecolor`
#    so that the figure backgrounds shown by the `matplotlib backend
#    <https://matplotlib.org/faq/usage_faq#what-is-a-backend>`__ are light gray
#    (the :rcraw:`savefig.facecolor` applied to saved figures is still white).
#    Proplot also controls the appearance of figures in Jupyter notebooks
#    using the new :rcraw:`inlineformat` setting, which is passed to
#    `~proplot.config.config_inline_backend` on import. This
#    imposes a higher-quality default `"inline" format
#    <https://ipython.readthedocs.io/en/stable/interactive/plotting.html>`__
#    and disables the backend-specific settings ``InlineBackend.rc`` and
#    ``InlineBackend.print_figure_kwargs``, ensuring that the figures you save
#    look like the figures displayed by the backend.
#
#    Proplot also changes the default :rcraw:`savefig.format`
#    from PNG to PDF for the following reasons:
#
#        #. Vector graphic formats are infinitely scalable.
#        #. Vector graphic formats are preferred by academic journals.
#        #. Nearly all academic journals accept figures in the PDF format alongside
#           the `EPS <https://en.wikipedia.org/wiki/Encapsulated_PostScript>`__ format.
#        #. The EPS format is outdated and does not support transparent graphic
#           elements.
#
#    In case you *do* need a raster format like PNG, proplot increases the
#    default :rcraw:`savefig.dpi` to 1000 dots per inch, which is
#    `recommended <https://www.pnas.org/page/authors/format>`__ by most journals
#    as the minimum resolution for figures containing lines and text. See the
#    :ref:`configuration section <ug_proplotrc>` for how to change these settings.
#

# %%
# Simple subplot
import numpy as np
import proplot as pplt
state = np.random.RandomState(51423)
data = 2 * (state.rand(100, 5) - 0.5).cumsum(axis=0)
fig, ax = pplt.subplots(suptitle='Single subplot', xlabel='x axis', ylabel='y axis')
# fig = pplt.figure(suptitle='Single subplot')  # equivalent to above
# ax = fig.subplot(xlabel='x axis', ylabel='y axis')
ax.plot(data, lw=2)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_subplot:
#
# Creating subplots
# -----------------
#
# Similar to matplotlib, subplots can be added to figures one-by-one
# or all at once. Each subplot will be an instance of
# `proplot.axes.Axes`. To add subplots all at once, use
# `proplot.figure.Figure.add_subplots` (or its shorthand,
# `proplot.figure.Figure.subplots`). Note that under the hood, the top-level
# proplot command `~proplot.ui.subplots` simply calls `~proplot.ui.figure`
# followed by `proplot.figure.Figure.add_subplots`.
#
# * With no arguments, `~proplot.figure.Figure.add_subplots` returns a subplot
#   generated from a 1-row, 1-column `~proplot.gridspec.GridSpec`.
# * With `ncols` or `nrows`, `~proplot.figure.Figure.add_subplots` returns a
#   simple grid of subplots from a `~proplot.gridspec.GridSpec` with
#   matching geometry in either row-major or column-major `order`.
# * With `array`, `~proplot.figure.Figure.add_subplots` returns an arbitrarily
#   complex grid of subplots from a `~proplot.gridspec.GridSpec` with matching
#   geometry. Here `array` is a 2D array representing a "picture" of the subplot
#   layout, where each unique integer indicates a `~matplotlib.gridspec.GridSpec`
#   slot occupied by the corresponding subplot and ``0`` indicates an empty space.
#   The returned subplots are contained in a `~proplot.gridspec.SubplotGrid`
#   (:ref:`see below <ug_subplotgrid>` for details).
#
# To add subplots one-by-one, use the `proplot.figure.Figure.add_subplot`
# command (or its shorthand `proplot.figure.Figure.subplot`).
#
# * With no arguments, `~proplot.figure.Figure.add_subplot` returns a subplot
#   generated from a 1-row, 1-column `~proplot.gridspec.GridSpec`.
# * With integer arguments, `~proplot.figure.Figure.add_subplot` returns
#   a subplot matching the corresponding `~proplot.gridspec.GridSpec` geometry,
#   as in matplotlib. Note that unlike matplotlib, the geometry must be compatible
#   with the geometry implied by previous `~proplot.figure.Figure.add_subplot` calls.
# * With a `~matplotlib.gridspec.SubplotSpec` generated by indexing a
#   `proplot.gridspec.GridSpec`, `~proplot.figure.Figure.add_subplot` returns a
#   subplot at the corresponding location. Note that unlike matplotlib, only
#   one `~proplot.figure.Figure.gridspec` can be used with each figure.
#
# As in matplotlib, to save figures, use `~matplotlib.figure.Figure.savefig` (or its
# shorthand `proplot.figure.Figure.save`). User paths in the filename are expanded
# with `os.path.expanduser`. In the following examples, we add subplots to figures
# with a variety of methods and then save the results to the home directory.
#
# .. warning::
#
#    Proplot employs :ref:`automatic axis sharing <ug_share>` by default. This lets
#    subplots in the same row or column share the same axis limits, scales, ticks,
#    and labels. This is often convenient, but may be annoying for some users. To
#    keep this feature turned off, simply :ref:`change the default settings <ug_rc>`
#    with e.g. ``pplt.rc.update('subplots', share=False, span=False)``. See the
#    :ref:`axis sharing section <ug_share>` for details.

# %%
# Simple subplot grid
import numpy as np
import proplot as pplt
state = np.random.RandomState(51423)
data = 2 * (state.rand(100, 5) - 0.5).cumsum(axis=0)
fig = pplt.figure()
ax = fig.subplot(121)
ax.plot(data, lw=2)
ax = fig.subplot(122)
fig.format(
    suptitle='Simple subplot grid', title='Title',
    xlabel='x axis', ylabel='y axis'
)
# fig.save('~/example1.png')  # save the figure
# fig.savefig('~/example1.png')  # alternative


# %%
# Complex grid
import numpy as np
import proplot as pplt
state = np.random.RandomState(51423)
data = 2 * (state.rand(100, 5) - 0.5).cumsum(axis=0)
array = [  # the "picture" (0 == nothing, 1 == subplot A, 2 == subplot B, etc.)
    [1, 1, 2, 2],
    [0, 3, 3, 0],
]
fig = pplt.figure(refwidth=1.8)
axs = fig.subplots(array)
axs.format(
    abc=True, abcloc='ul', suptitle='Complex subplot grid',
    xlabel='xlabel', ylabel='ylabel'
)
axs[2].plot(data, lw=2)
# fig.save('~/example2.png')  # save the figure
# fig.savefig('~/example2.png')  # alternative


# %%
# Really complex grid
import numpy as np
import proplot as pplt
state = np.random.RandomState(51423)
data = 2 * (state.rand(100, 5) - 0.5).cumsum(axis=0)
array = [  # the "picture" (1 == subplot A, 2 == subplot B, etc.)
    [1, 1, 2],
    [1, 1, 6],
    [3, 4, 4],
    [3, 5, 5],
]
fig, axs = pplt.subplots(array, figwidth=5, span=False)
axs.format(
    suptitle='Really complex subplot grid',
    xlabel='xlabel', ylabel='ylabel', abc=True
)
axs[0].plot(data, lw=2)
# fig.save('~/example3.png')  # save the figure
# fig.savefig('~/example3.png')  # alternative

# %%
# Using a GridSpec
import numpy as np
import proplot as pplt
state = np.random.RandomState(51423)
data = 2 * (state.rand(100, 5) - 0.5).cumsum(axis=0)
gs = pplt.GridSpec(nrows=2, ncols=2, pad=1)
fig = pplt.figure(span=False, refwidth=2)
ax = fig.subplot(gs[:, 0])
ax.plot(data, lw=2)
ax = fig.subplot(gs[0, 1])
ax = fig.subplot(gs[1, 1])
fig.format(
    suptitle='Subplot grid with a GridSpec',
    xlabel='xlabel', ylabel='ylabel', abc=True
)
# fig.save('~/example4.png')  # save the figure
# fig.savefig('~/example4.png')  # alternative

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_subplotgrid:
#
# Multiple subplots
# -----------------
#
# If you create subplots all-at-once with e.g. `~proplot.ui.subplots`,
# proplot returns a `~proplot.gridspec.SubplotGrid` of subplots. This list-like,
# array-like object provides some useful features and unifies the behavior of the
# three possible return types used by `matplotlib.pyplot.subplots`:
#
# * `~proplot.gridspec.SubplotGrid` behaves like a scalar when it is singleton.
#   In other words, if you make a single subplot with ``fig, axs = pplt.subplots()``,
#   then ``axs[0].method(...)`` is equivalent to ``axs.method(...)``.
# * `~proplot.gridspec.SubplotGrid` permits list-like 1D indexing, e.g. ``axs[1]``
#   to return the second subplot. The subplots in the grid are sorted by
#   `~proplot.axes.Axes.number` (see :ref:`this page <ug_abc>` for details
#   on changing the `~proplot.axes.Axes.number` order).
# * `~proplot.gridspec.SubplotGrid` permits array-like 2D indexing, e.g.
#   ``axs[1, 0]`` to return the subplot in the second row, first column, or
#   ``axs[:, 0]`` to return a `~proplot.gridspec.SubplotGrid` of every subplot
#   in the first column. The 2D indexing is powered by the underlying
#   `~proplot.gridspec.SubplotGrid.gridspec`.
#
# `~proplot.gridspec.SubplotGrid` includes methods for working
# simultaneously with different subplots. Currently, this includes
# the commands `~proplot.gridspec.SubplotGrid.format`,
# `~proplot.gridspec.SubplotGrid.panel_axes`,
# `~proplot.gridspec.SubplotGrid.inset_axes`,
# `~proplot.gridspec.SubplotGrid.altx`, and `~proplot.gridspec.SubplotGrid.alty`.
# In the below example, we use `proplot.gridspec.SubplotGrid.format` on the grid
# returned by `~proplot.ui.subplots` to format different subgroups of subplots
# (:ref:`see below <ug_format>` for more on the format command).
#
# .. note::
#
#    If you create subplots one-by-one with `~proplot.figure.Figure.subplot` or
#    `~proplot.figure.Figure.add_subplot`, a `~proplot.gridspec.SubplotGrid`
#    containing the numbered subplots is available via the
#    `proplot.figure.Figure.subplotgrid` property. As with subplots made
#    all-at-once, the subplots in the grid are sorted by `~proplot.axes.Axes.number`.

# %%
import proplot as pplt
import numpy as np
state = np.random.RandomState(51423)

# Selected subplots in a simple grid
fig, axs = pplt.subplots(ncols=4, nrows=4, refwidth=1.2, span=True)
axs.format(xlabel='xlabel', ylabel='ylabel', suptitle='Simple SubplotGrid')
axs.format(grid=False, xlim=(0, 50), ylim=(-4, 4))
axs[:, 0].format(facecolor='blush', edgecolor='gray7', linewidth=1)  # eauivalent
axs[:, 0].format(fc='blush', ec='gray7', lw=1)
axs[0, :].format(fc='sky blue', ec='gray7', lw=1)
axs[0].format(ec='black', fc='gray5', lw=1.4)
axs[1:, 1:].format(fc='gray1')
for ax in axs[1:, 1:]:
    ax.plot((state.rand(50, 5) - 0.5).cumsum(axis=0), cycle='Grays', lw=2)

# Selected subplots in a complex grid
fig = pplt.figure(refwidth=1, refnum=5, span=False)
axs = fig.subplots([[1, 1, 2], [3, 4, 2], [3, 4, 5]], hratios=[2.2, 1, 1])
axs.format(xlabel='xlabel', ylabel='ylabel', suptitle='Complex SubplotGrid')
axs[0].format(ec='black', fc='gray1', lw=1.4)
axs[1, 1:].format(fc='blush')
axs[1, :1].format(fc='sky blue')
axs[-1, -1].format(fc='gray4', grid=False)
axs[0].plot((state.rand(50, 10) - 0.5).cumsum(axis=0), cycle='Grays_r', lw=2)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_plots:
#
# Plotting stuff
# --------------
#
# Matplotlib includes `two different interfaces
# <https://matplotlib.org/stable/api/index.html>`__ for plotting stuff:
# a python-style object-oriented interface with axes-level commands
# like `matplotlib.axes.Axes.plot`, and a MATLAB-style `~matplotlib.pyplot` interface
# with global commands like `matplotlib.pyplot.plot` that track the "current" axes.
# Proplot builds upon the python-style interface using the `proplot.axes.PlotAxes`
# class. Since every axes used by proplot is a child of `~proplot.axes.PlotAxes`, we
# are able to add features directly to the axes-level commands rather than relying
# on a separate library of commands  (note that while some of these features may be
# accessible via `~matplotlib.pyplot` commands, this is not officially supported).
#
# For the most part, the features added by `~proplot.axes.PlotAxes` represent
# a *superset* of matplotlib. If you are not interested, you can use the plotting
# commands just like you would in matplotlib. Some of the core added features include
# more flexible treatment of :ref:`data arguments <ug_1dstd>`, recognition of
# :ref:`xarray and pandas <ug_1dintegration>` data structures, integration with
# proplot's :ref:`colormap <ug_apply_cmap>` and :ref:`color cycle <ug_apply_cycle>`
# tools, and on-the-fly :ref:`legend and colorbar generation <ug_guides_plot>`.
# In the below example, we create a 4-panel figure with the
# familiar "1D" plotting commands `~proplot.axes.PlotAxes.plot` and
# `~proplot.axes.PlotAxes.scatter`, along with the "2D" plotting commands
# `~proplot.axes.PlotAxes.pcolormesh` and `~proplot.axes.PlotAxes.contourf`.
# See the :ref:`1D plotting <ug_1dplots>` and :ref:`2D plotting <ug_2dplots>`
# sections for details on the features added by proplot.


# %%
import proplot as pplt
import numpy as np

# Sample data
N = 20
state = np.random.RandomState(51423)
data = N + (state.rand(N, N) - 0.55).cumsum(axis=0).cumsum(axis=1)

# Example plots
cycle = pplt.Cycle('greys', left=0.2, N=5)
fig, axs = pplt.subplots(ncols=2, nrows=2, figwidth=5, share=False)
axs[0].plot(data[:, :5], linewidth=2, linestyle='--', cycle=cycle)
axs[1].scatter(data[:, :5], marker='x', cycle=cycle)
axs[2].pcolormesh(data, cmap='greys')
m = axs[3].contourf(data, cmap='greys')
axs.format(
    abc='a.', titleloc='l', title='Title',
    xlabel='xlabel', ylabel='ylabel', suptitle='Quick plotting demo'
)
fig.colorbar(m, loc='b', label='label')


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_format:
#
# Formatting stuff
# ----------------
#
# Matplotlib includes `two different interfaces
# <https://matplotlib.org/stable/api/index.html>`__ for formatting stuff:
# a "python-style" object-oriented interface with instance-level commands
# like `matplotlib.axes.Axes.set_title`, and a "MATLAB-style" interface
# that tracks current axes and provides global commands like
# `matplotlib.pyplot.title`.
#
# Proplot provides the ``format`` command as an
# alternative "python-style" command for formatting a variety of plot elements.
# While matplotlib's one-liner commands still work, ``format`` only needs to be
# called once and tends to cut down on boilerplate code. You can call
# ``format`` manually or pass ``format`` parameters to axes-creation commands
# like `~proplot.figure.Figure.subplots`, `~proplot.figure.Figure.add_subplot`,
# `~proplot.axes.Axes.inset_axes`, `~proplot.axes.Axes.panel_axes`, and
# `~proplot.axes.CartesianAxes.altx` or `~proplot.axes.CartesianAxes.alty`. The
# keyword arguments accepted by ``format`` can be grouped as follows:
#
# * Figure settings. These are related to row labels, column labels, and
#   figure "super" titles -- for example, ``fig.format(suptitle='Super title')``.
#   See `proplot.figure.Figure.format` for details.
#
# * General axes settings. These are related to background patches,
#   a-b-c labels, and axes titles -- for example, ``ax.format(title='Title')``
#   See `proplot.axes.Axes.format` for details.
#
# * Cartesian axes settings (valid only for `~proplot.axes.CartesianAxes`).
#   These are related to *x* and *y* axis ticks, spines, bounds, and labels --
#   for example, ``ax.format(xlim=(0, 5))`` changes the x axis bounds.
#   See `proplot.axes.CartesianAxes.format` and
#   :ref:`this section <ug_cartesian>` for details.
#
# * Polar axes settings (valid only for `~proplot.axes.PolarAxes`).
#   These are related to azimuthal and radial grid lines, bounds, and labels --
#   for example, ``ax.format(rlim=(0, 10))`` changes the radial bounds.
#   See `proplot.axes.PolarAxes.format`
#   and :ref:`this section <ug_polar>` for details.
#
# * Geographic axes settings (valid only for `~proplot.axes.GeoAxes`).
#   These are related to map bounds, meridian and parallel lines and labels,
#   and geographic features -- for example, ``ax.format(latlim=(0, 90))``
#   changes the meridional bounds. See `proplot.axes.GeoAxes.format`
#   and :ref:`this section <ug_geoformat>` for details.
#
# * `~proplot.config.rc` settings. Any keyword matching the name
#   of an rc setting is locally applied to the figure and axes.
#   If the name has "dots", you can pass it as a keyword argument with
#   the "dots" omitted, or pass it to `rc_kw` in a dictionary. For example, the
#   default a-b-c label location is controlled by :rcraw:`abc.loc`. To change
#   this for an entire figure, you can use ``fig.format(abcloc='right')``
#   or ``fig.format(rc_kw={'abc.loc': 'right'})``.
#   See :ref:`this section <ug_config>` for more on rc settings.
#
# A ``format`` command is available on every figure and axes.
# `proplot.figure.Figure.format` accepts both figure and axes
# settings (applying them to each numbered subplot by default).
# Similarly, `proplot.axes.Axes.format` accepts both axes and figure
# settings. There is also a `proplot.gridspec.SubplotGrid.format`
# command that can be used to change settings for a subset of
# subplots -- for example, ``axs[:2].format(xtickminor=True)``
# turns on minor ticks for the first two subplots (see
# :ref:`this section <ug_subplotgrid>` for more on subplot grids).
# The below example shows the many keyword arguments accepted
# by ``format``, and demonstrates how ``format`` can be
# used to succinctly and efficiently customize plots.

# %%
import proplot as pplt
import numpy as np
fig, axs = pplt.subplots(ncols=2, nrows=2, refwidth=2, share=False)
state = np.random.RandomState(51423)
N = 60
x = np.linspace(1, 10, N)
y = (state.rand(N, 5) - 0.5).cumsum(axis=0)
axs[0].plot(x, y, linewidth=1.5)
axs.format(
    suptitle='Format command demo',
    abc='A.', abcloc='ul',
    title='Main', ltitle='Left', rtitle='Right',  # different titles
    ultitle='Title 1', urtitle='Title 2', lltitle='Title 3', lrtitle='Title 4',
    toplabels=('Column 1', 'Column 2'),
    leftlabels=('Row 1', 'Row 2'),
    xlabel='xaxis', ylabel='yaxis',
    xscale='log',
    xlim=(1, 10), xticks=1,
    ylim=(-3, 3), yticks=pplt.arange(-3, 3),
    yticklabels=('a', 'bb', 'c', 'dd', 'e', 'ff', 'g'),
    ytickloc='both', yticklabelloc='both',
    xtickdir='inout', xtickminor=False, ygridminor=True,
)

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_rc:
#
# Settings and styles
# -------------------
#
# A dictionary-like object named `~proplot.config.rc` is created when you import
# proplot. `~proplot.config.rc` is similar to the matplotlib `~matplotlib.rcParams`
# dictionary, but can be used to change both `matplotlib settings
# <https://matplotlib.org/stable/tutorials/introductory/customizing.html>`__ and
# :ref:`proplot settings <ug_rcproplot>`. The matplotlib-specific settings are
# stored in `~proplot.config.rc_matplotlib` (our name for `matplotlib.rcParams`) and
# the proplot-specific settings are stored in `~proplot.config.rc_proplot`.
# Proplot also includes a :rcraw:`style` setting that can be used to
# switch between `matplotlib stylesheets
# <https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html>`__.
# See the :ref:`configuration section <ug_config>` for details.
#
# To modify a setting for just one subplot or figure, you can pass it to
# `proplot.axes.Axes.format` or `proplot.figure.Figure.format`. To temporarily
# modify setting(s) for a block of code, use `~proplot.config.Configurator.context`.
# To modify setting(s) for the entire python session, just assign it to the
# `~proplot.config.rc` dictionary or use `~proplot.config.Configurator.update`.
# To reset everything to the default state, use `~proplot.config.Configurator.reset`.
# See the below example.


# %%
import proplot as pplt
import numpy as np

# Update global settings in several different ways
pplt.rc.metacolor = 'gray6'
pplt.rc.update({'fontname': 'Source Sans Pro', 'fontsize': 11})
pplt.rc['figure.facecolor'] = 'gray3'
pplt.rc.axesfacecolor = 'gray4'
# pplt.rc.save()  # save the current settings to ~/.proplotrc

# Apply settings to figure with context()
with pplt.rc.context({'suptitle.size': 13}, toplabelcolor='gray6', metawidth=1.5):
    fig = pplt.figure(figwidth=6, sharey='limits', span=False)
    axs = fig.subplots(ncols=2)

# Plot lines with a custom cycler
N, M = 100, 7
state = np.random.RandomState(51423)
values = np.arange(1, M + 1)
cycle = pplt.get_colors('grays', M - 1) + ['red']
for i, ax in enumerate(axs):
    data = np.cumsum(state.rand(N, M) - 0.5, axis=0)
    lines = ax.plot(data, linewidth=3, cycle=cycle)

# Apply settings to axes with format()
axs.format(
    grid=False, xlabel='xlabel', ylabel='ylabel',
    toplabels=('Column 1', 'Column 2'),
    suptitle='Rc settings demo',
    suptitlecolor='gray7',
    abc='[A]', abcloc='l',
    title='Title', titleloc='r', titlecolor='gray7'
)

# Reset persistent modifications from head of cell
pplt.rc.reset()


# %%
import proplot as pplt
import numpy as np
# pplt.rc.style = 'style'  # set the style everywhere

# Sample data
state = np.random.RandomState(51423)
data = state.rand(10, 5)

# Set up figure
fig, axs = pplt.subplots(ncols=2, nrows=2, span=False, share=False)
axs.format(suptitle='Stylesheets demo')
styles = ('ggplot', 'seaborn', '538', 'bmh')

# Apply different styles to different axes with format()
for ax, style in zip(axs, styles):
    ax.format(style=style, xlabel='xlabel', ylabel='ylabel', title=style)
    ax.plot(data, linewidth=3)
