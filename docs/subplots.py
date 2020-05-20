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
# .. _ug_subplots:
#
# Subplots
# ========
#
# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_autosize:
#
# Automatic figure size
# ---------------------
#
# By default, ProPlot automatically determines the suitable figure size given
# the geometry of your subplot grid and the physical dimensions of a "reference"
# subplot. ProPlot can also determine the suitable figure height given a fixed
# figure width, and figure width given a fixed figure height (which can be
# particularly useful when preparing publications).
#
# This algorithm is controlled by
# a variety of `~proplot.ui.subplots` keyword arguments:
#
# * The `ref` parameter sets the reference subplot number (default is ``1``,
#   i.e. the subplot in the upper left corner).
# * The `aspect` parameter sets the reference subplot aspect ratio (default
#   is ``1``). You can also use the built-in matplotlib
#   `~matplotlib.axes.Axes.set_aspect` method.
# * The `axwidth` and `axheight` parameters set the physical dimensions of
#   the *reference subplot* (default is ``axwidth=2``). If one is specified,
#   the other is calculated to satisfy `aspect`. If both are specified,
#   `aspect` is ignored. The physical dimensions of the *figure* are
#   determined automatically.
# * The `width` and `height` parameters set the physical dimensions of the
#   *figure*. If one is specified, the other is calculated to satisfy `aspect`
#   and the subplot spacing. If both are specified (or if the matplotlib
#   `figsize` parameter is specified), `aspect` is ignored.
# * The `journal` parameter constrains the physical dimensions of the figure
#   so it meets requirements for submission to an academic journal. For example,
#   figures created with ``journal='nat1'`` are sized as single-column
#   *Nature* figures. See :ref:`this table <journal_table>` for the list
#   of available journal specifications (feel free to add to this table by
#   submitting a PR).
#
# The below examples demonstrate the default behavior of the automatic figure
# sizing algorithm, and how it can be controlled with `~proplot.ui.subplots`
# keyword arguments.
#
# .. important::
#
#    The automatic figure size algorithm has the following notable properties:
#
#    * For very simple subplot grids (i.e. subplots created with the `ncols` and
#      `nrows` arguments), the arguments `aspect`, `axwidth`, and `axheight` apply
#      to every subplot in the figure -- not just the reference subplot.
#    * When the reference subplot `aspect ratio\
#      <https://matplotlib.org/2.0.2/examples/pylab_examples/equal_aspect_ratio.html>`__
#      has been fixed (e.g. with ``ax.set_aspect(1)``) or is set
#      to ``'equal'`` (as with :ref:`map projections <ug_geo>` and
#      `~matplotlib.axes.Axes.imshow` images), the fixed aspect ratio is used
#      and the `~proplot.ui.subplots` `aspect` parameter is ignored. This
#      is critical for getting the figure size right when working with grids
#      of images and grids of projections.
#    * When `~proplot.axes.Axes.colorbar`\ s and `~proplot.axes.Axes.panel`\ s
#      are present in the figure, their physical widths are *preserved* during
#      figure resizing. ProPlot specifies their widths in physical units to help
#      avoid colorbars that look "too skinny" or "too fat".

# %%
import proplot as plot
import numpy as np

# Auto sized grid of cartopy projections
fig, axs = plot.subplots(ncols=2, nrows=3, proj='robin')
axs.format(
    land=True, landcolor='k',
    suptitle='Auto figure sizing with grid of cartopy projections'
)

# Auto sized grid of images
state = np.random.RandomState(51423)
fig, axs = plot.subplots(ncols=3, nrows=2, axwidth=1.5)
colors = state.rand(13, 10, 3).cumsum(axis=2)
colors /= colors.max()
axs.imshow(colors)
axs.format(
    suptitle='Auto figure sizing with grid of images'
)

# %%
import proplot as plot

# Change the reference subplot width
suptitle = 'Effect of subplot width on figure size'
for axwidth in ('4cm', '6cm'):
    fig, axs = plot.subplots(ncols=2, axwidth=axwidth,)
    axs[0].format(
        suptitle=suptitle,
        title=f'axwidth = {axwidth}', titleweight='bold',
        titleloc='uc', titlecolor='red9',
    )

# Change the reference subplot aspect ratio
for aspect in (1, (3, 2)):
    fig, axs = plot.subplots(ncols=2, nrows=2, axwidth=1.6, aspect=aspect)
    axs[0].format(
        suptitle='Effect of subplot aspect ratio on figure size',
        title=f'aspect = {aspect}', titleweight='bold',
        titleloc='uc', titlecolor='red9',
    )

# %%
import proplot as plot

# Change the reference subplot in presence of unequal width/height ratios
for ref in (1, 2):
    fig, axs = plot.subplots(
        ref=ref, nrows=3, ncols=3, wratios=(3, 2, 2),
        axwidth=1.1,
    )
    axs[ref - 1].format(
        suptitle='Effect of reference subplot on figure size',
        title='reference', titleweight='bold',
        titleloc='uc', titlecolor='red9'
    )

# Change the reference subplot in a complex grid
for ref in (1, 2):
    fig, axs = plot.subplots(
        [[1, 2], [1, 3]],
        ref=ref, axwidth=1.8, span=False
    )
    axs[ref - 1].format(
        suptitle='Effect of reference subplot on figure size',
        title='reference', titleweight='bold',
        titleloc='uc', titlecolor='red9'
    )


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_tight:
#
# Automatic subplot spacing
# -------------------------
#
# In addition to automatic figure sizing, by default ProPlot applies a *tight layout*
# algorithm to every figure. This algorithm automatically adjusts the space between
# subplot rows and columns and the figure edge to accommodate labels.
# It can be disabled by passing ``tight=False`` to `~proplot.ui.subplots`.
# While matplotlib has `its own tight layout algorithm
# <https://matplotlib.org/3.1.1/tutorials/intermediate/tight_layout_guide.html>`__,
# ProPlot's algorithm may change the figure size to accommodate the correct spacing
# and permits variable spacing between subsequent subplot rows and columns (see the
# new `~proplot.gridspec.GridSpec` class for details).
#
# The tight layout algorithm can also be overridden. When you pass a
# spacing argument like `left`, `right`, `top`, `bottom`, `wspace`, or
# `hspace` to `~proplot.ui.subplots`, that value is always respected.
# For example:
#
# * ``left='2em'`` fixes the left margin, while the right, bottom, and
#   top margins are determined automatically.
# * ``wspace='1em'`` fixes the spaces between subplot columns, while the spaces
#   between subplot rows are determined automatically.
# * ``wspace=('3em', None)`` fixes the space between the first two columns of
#   a three-column plot, while the space between the second two columns is
#   determined automatically.
#
# The below examples demonstrate how the tight layout algorithm permits
# variable spacing between subplot rows and columns.

# %%
import proplot as plot

# Automatic spacing for all margins and between all columns and rows
fig, axs = plot.subplots(nrows=3, ncols=3, axwidth=1.1, share=0)

# Formatting that stress-tests the algorithm
axs[4].format(
    title='title\ntitle\ntitle',
    suptitle='Tight layout with variable row-column spacing'
)
axs[1].format(ylabel='ylabel\nylabel\nylabel')
axs[:4:2].format(xlabel='xlabel\nxlabel\nxlabel')
axs.format(
    rowlabels=['Row 1', 'Row 2', 'Row 3'],
    collabels=['Column 1', 'Column 2', 'Column 3']
)

# %%
import proplot as plot

# Manual spacing for certain margins and between certain columns and rows
fig, axs = plot.subplots(
    ncols=4, nrows=3, axwidth=1.1, span=False,
    bottom='5em', right='5em',  # margin spacing overrides
    wspace=(0, 0, None), hspace=(0, None),  # column and row spacing overrides
)

# Formatting that stress-tests the algorithm
axs.format(
    xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), xlocator=1, ylocator=1,
    suptitle='Tight layout with user overrides',
    rowlabels=['Row 1', 'Row 2', 'Row 3'],
    collabels=['Column 1', 'Column 2', 'Column 3', 'Column 4']
)
axs[0, :].format(xtickloc='top')
axs[2, :].format(xtickloc='both')
axs[:, 1].format(ytickloc='neither')
axs[:, 2].format(ytickloc='right')
axs[:, 3].format(ytickloc='both')
axs[-1, :].format(title='Title\nTitle\nTitle', xlabel='xlabel')
axs[:, 0].format(ylabel='ylabel\nylabel')


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_units:
#
# Arbitrary physical units
# ------------------------
#
# ProPlot supports arbitrary *physical units* for controlling the figure
# `width` and `height`, the reference subplot `axwidth` and `axheight`, the
# gridspec spacing values `left`, `right`, `bottom`, `top`, `wspace`, and
# `hspace`, and in a few other places, e.g. `~proplot.axes.Axes.panel` and
# `~proplot.axes.Axes.colorbar` widths. This feature is powered by the
# `~proplot.utils.units` function.
#
# If a sizing argument is numeric, the units are inches or points; if it is
# string, the units are converted to inches or points by
# `~proplot.utils.units`. A table of acceptable units is found in the
# `~proplot.utils.units` documentation. They include centimeters,
# millimeters, pixels,
# `em-heights <https://en.wikipedia.org/wiki/Em_(typography)>`__, and
# `points <https://en.wikipedia.org/wiki/Point_(typography)>`__.

# %%
import proplot as plot
import numpy as np
with plot.rc.context(fontsize='12px'):
    fig, axs = plot.subplots(
        ncols=3, width='15cm', height='2.5in',
        wspace=('10pt', '20pt'), right='10mm'
    )
    panel = axs[2].panel_axes('r', width='2em')
    panel.format(xlim=(0, 1))
axs.format(
    suptitle='Arguments with arbitrary units',
    xlabel='x axis', ylabel='y axis',
    xlim=(0, 1), ylim=(0, 1),
)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_abc:
#
# A-b-c subplot labels
# --------------------
#
# ProPlot can be used to add "a-b-c" labels to subplots. This is possible
# because `~proplot.ui.subplots` assigns unique `~proplot.axes.Axes.number`\ s
# to each subplot. If you :ref:`passed an array <ug_intro>` to
# `~proplot.ui.subplots`, the subplot numbers correspond to the numbers in
# the array. If you used the `ncols` and `nrows` keyword arguments, the number
# order is row-major by default but can be switched to column-major by passing
# ``order='C'`` to `~proplot.ui.subplots`. The number order also determines the
# subplot order in the `~proplot.ui.SubplotsContainer` returned by
# `~proplot.ui.subplots`.
#
# To turn on "a-b-c" labels, set :rcraw:`abc` to ``True`` or pass
# ``abc=True`` to `~proplot.axes.Axes.format` (see
# :ref:`the format command <ug_format>` for details). To change the label
# style, modify :rcraw:`abc.style` or pass e.g. ``abcstyle='A.'`` to
# `~proplot.axes.Axes.format`. You can also modify the "a-b-c" label
# location, weight, and size with the :rcraw:`abc.loc`, :rcraw:`abc.weight`,
# and :rcraw:`abc.size` settings.

# %%
import proplot as plot
fig, axs = plot.subplots(nrows=8, ncols=8, axwidth=0.7, space=0)
axs.format(
    abc=True, abcloc='ur', xlabel='x axis', ylabel='y axis',
    xticks=[], yticks=[], suptitle='Subplot labels demo'
)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_share:
#
# Axis sharing
# ------------
#
# Redundant labels are a common problem for figures with lots of subplots. To
# address this, `matplotlib.pyplot.subplots` includes `sharex` and `sharey`
# keywords that permit sharing axis limits, ticks, and tick labels between like
# rows and columns of subplots. However there is no convenient way to enable
# axis sharing between complex grids of subplots, nor is there a way to share axis
# labels between subplots in the same row or column.
#
# ProPlot expands upon this feature by introducing a new option for drawing
# labels that "span" subplots in the same row or column, controlled by the
# `spanx` and `spany` keywords, along with four axis-sharing "levels"
# rather than just one, controlled with the `sharex` and `sharey` keywords:
#
# * Level ``0`` disables axis sharing.
# * Level ``1`` shares duplicate *x* and *y* axis labels, but nothing else.
# * Level ``2`` is the same as ``1``, but the *x* and *y* axis limits, ticks,
#   and scales are also shared.
# * Level ``3`` is the same as ``2``, but the *x* and *y* tick labels are
#   also shared.
#
# These features are best illustrated by example (see below).
#
# .. note::
#
#    Since "shared" and "spanning" labels are determined automatically based on
#    position of each subplot in the `~proplot.gridspec.GridSpec`, these features
#    work for arbitrarily complex figures. This can be done without any ambiguity
#    because ProPlot ensures each figure has only one `~proplot.gridspec.GridSpec`.

# %%
import proplot as plot
import numpy as np
N = 50
M = 40
state = np.random.RandomState(51423)
colors = plot.Colors('grays_r', M, left=0.1, right=0.8)
datas = []
for scale in (1, 3, 7, 0.2):
    data = scale * (state.rand(N, M) - 0.5).cumsum(axis=0)[N // 2:, :]
    datas.append(data)

# Same plot with different sharing and spanning settings
for share in (0, 1, 2, 3):
    fig, axs = plot.subplots(
        ncols=4, aspect=1, axwidth=1.06,
        sharey=share, spanx=share // 2
    )
    for ax, data in zip(axs, datas):
        on = ['off', 'on'][share // 2]
        ax.plot(data, cycle=colors)
        ax.format(
            suptitle=f'Sharing level {share}, spanning labels {on}',
            grid=False, xlabel='spanning', ylabel='shared'
        )

# %%
import proplot as plot
import numpy as np
plot.rc.reset()
plot.rc.cycle = 'Set3'
state = np.random.RandomState(51423)
titles = ['With redundant labels', 'Without redundant labels']

# Same plot with and without default sharing settings
for mode in (0, 1):
    fig, axs = plot.subplots(
        nrows=4, ncols=4, share=3 * mode,
        span=1 * mode, axwidth=1
    )
    for ax in axs:
        ax.plot((state.rand(100, 20) - 0.4).cumsum(axis=0))
    axs.format(
        xlabel='xlabel', ylabel='ylabel', suptitle=titles[mode],
        abc=True, abcloc='ul',
        grid=False, xticks=25, yticks=5
    )
