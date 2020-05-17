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
# ProPlot works by subclassing the matplotlib `~matplotlib.figure.Figure` and
# `~matplotlib.axes.Axes` classes. You can generate grids of proplot
# `~proplot.axes.Axes` axes on a proplot `~proplot.figure.Figure` using the
# `~proplot.ui.subplots` command.
#
# .. code-block:: python
#
#   import proplot as plot
#   fig, axs = plot.subplots(...)
#
# Just like `matplotlib.pyplot.subplots`, you can use
# `~proplot.ui.subplots` without arguments to generate a single-axes figure
# or with `ncols` or `nrows` to set up simple grids of subplots.  You can
# *also* draw arbitrarily complex grids in ProPlot by passing 2D arrays of
# integers to `~proplot.ui.subplots`. Just think of this array as a
# "picture" of your figure, where each unique integer indicates a slot that
# is occupied by the corresponding subplot and ``0`` indicates an empty space.
#
# In the below examples, we create subplot grids with `~proplot.ui.subplots`
# and modify the axes using `~proplot.axes.Axes.format` and
# `~proplot.ui.SubplotsContainer`. See the :ref:`formatting guide <ug_format>`
# and :ref:`subplots container <ug_container>` sections for details.
#
# Please note that by default, ProPlot sets :rcraw:`figure.facecolor` to gray,
# :rcraw:`savefig.facecolor` to white, and :rcraw:`savefig.transparent` to ``True``.
# That is, the default display background is gray, the default background for
# saved figures is transparent, and the default background is white when you pass
# ``transparent=False`` to `~matplotlib.figure.Figure.savefig`.
# ProPlot also sets the default :rcraw:`savefig.format` to PDF, because
# (1) vector graphic formats are always more suitable for matplotlib figures than
# raster formats, (2) most academic journals these days accept PDF format figures
# alongside the older EPS format, (3) PDF figures are easy to embed in LaTeX documents,
# and (4) the EPS format does not support transparent graphic elements. If you *do*
# need raster graphics, ProPlot sets the default :rcraw:`savefig.dpi` to 1200 dots per
# inch, which is recommended by most journals as the minimum resolution for rasterized
# figures containing lines and text. See the :ref:`configuration section <ug_proplotrc>`
# for how to change these settings.

# %%
import proplot as plot
import numpy as np
state = np.random.RandomState(51423)
data = 2 * (state.rand(100, 5) - 0.5).cumsum(axis=0)

# Simple plot
fig, axs = plot.subplots(ncols=2)
axs[0].plot(data, lw=2)
axs[0].format(xticks=20, xtickminor=False)
axs.format(
    suptitle='Simple subplot grid', title='Title',
    xlabel='x axis', ylabel='y axis'
)

# Complex grid
array = [  # the "picture"; 1 == subplot A, 2 == subplot B, etc.
    [1, 1, 2, 2],
    [0, 3, 3, 0],
]
fig, axs = plot.subplots(array, axwidth=1.8)
axs.format(
    abc=True, abcloc='ul', suptitle='Complex subplot grid',
    xlabel='xlabel', ylabel='ylabel'
)
axs[2].plot(data, lw=2)

# Really complex grid
array = [  # the "picture"
    [1, 1, 2],
    [1, 1, 6],
    [3, 4, 4],
    [3, 5, 5],
]
fig, axs = plot.subplots(array, width=5, span=False)
axs.format(
    suptitle='Really complex subplot grid',
    xlabel='xlabel', ylabel='ylabel', abc=True
)
axs[0].plot(data, lw=2)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_plots:
#
# Plotting data
# -------------
#
# Matplotlib includes two APIs for creating plots: an object-oriented API,
# and a MATLAB-style `~matplotlib.pyplot` API (see matplotlib's
# `API guide <https://matplotlib.org/3.2.1/api/index.html>`__ for details).
# If you are already familiar with the object-oriented API, plotting in
# ProPlot will look exactly the same to you. This is because ProPlot's plotting
# features are a strict *superset* of matplotlib's features. Rather than creating
# a brand new interface, ProPlot simply builds upon the existing matplotlib constructs
# of the `~matplotlib.axes.Axes` and the `~matplotlib.figure.Figure`.
# This means a shallow learning curve for the average matplotlib user.
#
# In the below example, we create a 4-panel figure with the familiar matplotlib
# commands `~matplotlib.axes.Axes.plot`, `~matplotlib.axes.Axes.scatter`,
# `~matplotlib.axes.Axes.pcolormesh`, and `~matplotlib.axes.Axes.contourf`.
# See the :ref:`1d plotting <ug_1dplots>` and :ref:`2d plotting <ug_2dplots>`
# sections for details on the plotting features added by ProPlot.


# %%
import proplot as plot
import numpy as np

# Sample data
N = 20
state = np.random.RandomState(51423)
data = (state.rand(N, N) - 0.5).cumsum(axis=0).cumsum(axis=1)

# Example plots
cycle = plot.Cycle('greys', left=0.2, N=5)
fig, axs = plot.subplots(ncols=2, nrows=2, share=0, width=5)
axs[0].plot(data[:, :5], linewidth=2, linestyle='--', cycle=cycle)
axs[1].scatter(data[:, :5], marker='x', cycle=cycle)
axs[2].pcolormesh(data, cmap='greys')
axs[3].contourf(data, cmap='greys')
axs.format(abc=True, xlabel='xlabel', ylabel='ylabel', suptitle='Quick plotting demo')


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_format:
#
# Formatting plots
# ----------------
#
# Every `~matplotlib.axes.Axes` returned by `~proplot.ui.subplots` has a
# ``format`` method. This is your one-stop-shop for changing axes settings.
# Keyword arguments passed to ``format`` are interpreted as follows:
#
# #. Any keyword matching the name of an `~proplot.config.rc` setting
#    will be used to update the axes. If the name has "dots", you can omit them
#    (e.g. ``titleloc='left'`` to change the ``title.loc`` property). See the
#    :ref:`configuration section <ug_config>` for details.
# #. Valid keywords arguments are passed to the
#    `proplot.axes.CartesianAxes.format`, `proplot.axes.PolarAxes.format`, or
#    `proplot.axes.GeoAxes.format` methods. These change settings that are
#    specific to the axes type. For example:

#    * To change the *x* axis bounds on a `~proplot.axes.CartesianAxes`,
#      use e.g. ``xlim=(0, 5)``.
#    * To change the radial bounds on a `~proplot.axes.PolarAxes`, use e.g.
#      ``rlim=(0, 10)``.
#    * To change the meridional bounds on a `~proplot.axes.GeoAxes`,
#      use e.g. ``lonlim=(-90, 0)``.
#
#
# #. All remaining keyword arguments are passed to the base `proplot.axes.Axes.format`
#    method. `~proplot.axes.Axes` is the base class for all other axes classes.
#    This changes things that are the same for all axes types, like titles and
#    a-b-c subplot labels (for example, ``abcstyle='A.'``).
#
# ``format`` lets you use simple shorthands for changing all kinds of
# settings at once, instead of one-liner setter methods like
# ``ax.set_title()``, ``ax.set_xlabel()``, and ``ax.xaxis.tick_params()``. It
# is also integrated with the `~proplot.constructor.Locator`,
# `~proplot.constructor.Formatter`, and `~proplot.constructor.Scale`
# constructor functions (see the :ref:`x and y axis settings <ug_xy_axis>`
# section for details).
#
# The below example shows the many different keyword arguments accepted by
# ``format``, and demonstrates how ``format`` can be used to succinctly and
# efficiently customize your plots.

# %%
import proplot as plot
import numpy as np
fig, axs = plot.subplots(ncols=2, nrows=2, share=0, tight=True, axwidth=2)
state = np.random.RandomState(51423)
x = np.linspace(1, 10, 80)
y = (state.rand(80, 5) - 0.5).cumsum(axis=0)
axs[0].plot(x, y, linewidth=1.5)
axs.format(
    suptitle='Format command demo',
    abc=True, abcloc='ul', abcstyle='A.',
    title='Main', ltitle='Left', rtitle='Right',  # different titles
    urtitle='Title A', lltitle='Title B', lrtitle='Title C',  # extra titles
    collabels=['Column label 1', 'Column label 2'],
    rowlabels=['Row label 1', 'Row label 2'],
    xlabel='x-axis', ylabel='y-axis',
    xscale='log',
    xlim=(1, 10), xticks=1,
    ymargin=0.05, yticks=plot.arange(-2, 2),
    yticklabels=('a', 'bb', 'c', 'dd', 'e'),
    ytickloc='both', yticklabelloc='both',
    xtickdir='inout', xtickminor=False, ygridminor=True,
    linewidth=0.8, gridlinewidth=0.8, gridminorlinewidth=0.5,
)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_rc:
#
# Changing rc settings
# --------------------
#
# A special object named `~proplot.config.rc` is created whenever you import
# ProPlot. `~proplot.config.rc` is similar to the matplotlib
# `~matplotlib.rcParams` dictionary, but can be used to change (1)
# matplotlib's `builtin settings
# <https://matplotlib.org/tutorials/introductory/customizing.html>`_, (2)
# ProPlot's :ref:`added settings <rc_added>`, and (3) :ref:`quick settings
# <rc_quick>` that can be used to change lots of matplotlib and ProPlot
# settings at once. `~proplot.config.rc` also provides a ``style`` parameter
# that can be used to switch between `matplotlib stylesheets\
# <https://matplotlib.org/3.1.1/gallery/style_sheets/style_sheets_reference.html>`__.
# See the :ref:`configuration section <ug_config>` for details.
#
# To modify a setting for just one subplot, you can pass it to the
# `~proplot.axes.Axes` `~proplot.axes.Axes.format` method. To temporarily
# modify setting(s) for a block of code, use
# `~proplot.config.rc_configurator.context`. To modify setting(s) for the
# entire python session, just assign it to the `~proplot.config.rc` object or
# use `~proplot.config.rc_configurator.update`.  To reset everything to the
# default state, use `~proplot.config.rc_configurator.reset`. See the below
# example.

# %%
import proplot as plot
import numpy as np

# Update global settings in several different ways
plot.rc.cycle = 'colorblind'
plot.rc.color = 'gray6'
plot.rc.update({'fontname': 'Noto Sans'})
plot.rc['figure.facecolor'] = 'gray3'
plot.rc.axesfacecolor = 'gray4'
plot.rc.save()  # save the current settings to a ~/.proplotrc file

# Apply settings to figure with context()
with plot.rc.context({'suptitle.size': 11}, toplabelcolor='gray6', linewidth=1.5):
    fig, axs = plot.subplots(ncols=2, aspect=1, width=6, span=False, sharey=2)

# Plot lines
N, M = 100, 6
state = np.random.RandomState(51423)
values = np.arange(1, M + 1)
for i, ax in enumerate(axs):
    data = np.cumsum(state.rand(N, M) - 0.5, axis=0)
    lines = ax.plot(data, linewidth=3, cycle='Grays')

# Apply settings to axes with format()
axs.format(
    grid=False, xlabel='x label', ylabel='y label',
    collabels=['Column label 1', 'Column label 2'],
    suptitle='Rc settings demo',
    suptitlecolor='gray7',
    abc=True, abcloc='l', abcstyle='A)',
    title='Title', titleloc='r', titlecolor='gray7'
)
ay = axs[-1].twinx()
ay.format(ycolor='red', linewidth=1.5, ylabel='secondary axis')
ay.plot((state.rand(100) - 0.2).cumsum(), color='r', lw=3)

# Reset persistent modifications from head of cell
plot.rc.reset()

# %%
import proplot as plot
import numpy as np
# plot.rc.style = 'style'  # set the style everywhere

# Set up figure
styles = ('ggplot', 'seaborn', '538', 'bmh')
state = np.random.RandomState(51423)
data = state.rand(10, 5)
fig, axs = plot.subplots(ncols=2, nrows=2, span=False, share=False)

# Apply different styles to different axes with format()
axs.format(suptitle='Stylesheets demo')
for ax, style in zip(axs, styles):
    ax.format(style=style, xlabel='xlabel', ylabel='ylabel', title=style)
    ax.plot(data, linewidth=3)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_container:
#
# The subplots container
# ----------------------
#
# Instead of an `~numpy.ndarray` of axes, `~proplot.ui.subplots` returns a
# `~proplot.ui.SubplotsContainer` instance. This container behaves
# like a python *list*, but lets you call any arbitrary method on multiple
# axes at once. It supports both 2D indexing (e.g. ``axs[0, 1]``) and 1D
# indexing (e.g. ``axs[2]``), and is row-major by default. Further, slicing a
# subplot container (e.g. ``axs[:, 0]``) returns another subplot container.
#
# In the below example, the `~proplot.ui.SubplotsContainer` returned by
# `~proplot.ui.subplots` is used to call `~proplot.axes.Axes.format` on
# several axes at once.

# %%
import proplot as plot
import numpy as np
state = np.random.RandomState(51423)
fig, axs = plot.subplots(ncols=4, nrows=4, axwidth=1.2)
axs.format(
    xlabel='xlabel', ylabel='ylabel', suptitle='SubplotsContainer demo',
    grid=False, xlim=(0, 50), ylim=(-4, 4)
)

# Various ways to select subplots in the subplot grid
axs[:, 0].format(facecolor='blush', color='gray7', linewidth=1)
axs[0, :].format(facecolor='sky blue', color='gray7', linewidth=1)
axs[0].format(color='black', facecolor='gray5', linewidth=1.4)
axs[1:, 1:].format(facecolor='gray1')
for ax in axs[1:, 1:]:
    ax.plot((state.rand(50, 5) - 0.5).cumsum(axis=0), cycle='Grays', lw=2)
