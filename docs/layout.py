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
# .. _ug_layout:
#
# The layout
# ==========
#
# This section documents a variety of features related to ProPlot subplots,
# including automatic a-b-c subplot labels, axis sharing between subplots,
# automatic spacing between subplots, and a unique feature where the figure
# size is determined automatically from the subplot layout.
#
# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_abc:
#
# A-b-c labels
# ------------
#
# ProPlot can quickly add "a-b-c" labels to subplots. This is possible because
# ProPlot assigns a unique `~proplot.axes.Axes.number` to each subplot. The
# subplot number can be manually controlled by passing a `number` keyword to
# `~proplot.figure.Figure.add_subplot`. Otherwise, the subplot number is
# incremented by ``1`` each time you call `~proplot.figure.Figure.add_subplot`.
#
# If you draw all of your subplots at once with `~proplot.figure.Figure.add_subplots`,
# the subplot numbers depend on the input arguments. If you
# :ref:`passed an array <ug_intro>`, the subplot numbers correspond to the numbers
# in the array. But if you used the `ncols` and `nrows` keyword arguments, the
# number order is row-major by default and can be switched to column-major by
# passing ``order='F'``. The number order also determines the subplot order in
# the `~proplot.figure.SubplotGrid` returned by `~proplot.figure.Figure.add_subplots`.
#
# To turn on "a-b-c" labels, set :rcraw:`abc` to ``True`` or pass ``abc=True``
# to `~proplot.axes.Axes.format` (see :ref:`the format command <ug_format>`
# for details). To change the label style, set :rcraw:`abc` to e.g. ``'A.'`` or
# pass e.g. ``abc='A.'`` to `~proplot.axes.Axes.format`. You can also modify
# the "a-b-c" label location, weight, and size with the :rcraw:`abc.loc`,
# :rcraw:`abc.weight`, and :rcraw:`abc.size` settings. Also note that if the
# an "a-b-c" label and title are in the same position, they are automatically
# offset away from each other.
#
# .. note::
#
#    "Inner" a-b-c labels and titles are surrounded with a white border when
#    :rcraw:`abc.border` and :rcraw:`title.border` are ``True`` (the default).
#    White boxes can be used instead by setting :rcraw:`abc.bbox` and
#    :rcraw:`title.bbox` to ``True``. These options help labels stand out
#    against plotted content. Any text can be given "borders" or "boxes" by
#    passing ``border=True`` or ``bbox=True`` to `proplot.axes.Axes.text`.

# %%
import proplot as pplt
fig = pplt.figure(space=0, refwidth='10em')
axs = fig.subplots(nrows=3, ncols=3)
axs.format(
    abc='A.', abcloc='ul',
    xticks='null', yticks='null', facecolor='gray5',
    xlabel='x axis', ylabel='y axis',
    suptitle='A-b-c label offsetting, borders, and boxes',
)
axs[:3].format(abcloc='l', titleloc='l', title='Title')
axs[-3:].format(abcbbox=True)  # also disables abcborder
# axs[:-3].format(abcborder=True)  # this is already the default

# %%
import proplot as pplt
fig = pplt.figure(space=0, refwidth=0.7)
axs = fig.subplots(nrows=8, ncols=8)
axs.format(
    abc=True, abcloc='ur',
    xlabel='x axis', ylabel='y axis', xticks=[], yticks=[],
    suptitle='A-b-c label stress test'
)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_autosize:
#
# Automatic size
# --------------
#
# By default, ProPlot determines the suitable figure size given the
# geometry of the subplot grid and the size of a "reference" subplot.
# This "reference" subplot is specified with the `~proplot.figure.Figure`
# keyword `refnum` (default is ``1``, i.e. the first subplot added to the figure
# or the subplot in the upper-left corner when generated with `~proplot.ui.subplots`).
# ProPlot can also determine the suitable figure height given a fixed figure
# width, and the suitable figure width given a fixed figure height.
#
# The figure size is ultimately controlled by the following
# `~proplot.figure.Figure` keyword arguments:
#
# * `refwidth` and `refheight` set the physical dimensions of the reference subplot
#   (default is :rc:`subplots.refwidth`). If one is specified, the other is calculated
#   to satisfy the subplot aspect ratio `refaspect` (default is ``1``). If both are
#   specified, `refaspect` is ignored. When these keyword arguments are used, the
#   width and height of the figure are both determined automatically.
# * `figwidth` and `figheight` set the physical dimensions of the figure.
#   If one is specified, the other is calculated to satisfy `refaspect`
#   and the subplot spacing. If both are specified, or if the `figsize` parameter
#   is specified, the figure size is fixed and `refaspect` is ignored.
# * `journal` constrains the physical dimensions of the figure to meet requirements
#   for submission to an academic journal. For example, ``journal='nat1'``
#   results in a width suitable for single-column *Nature* figures. See
#   :ref:`this table <journal_table>` for the list of available journal
#   specifications (feel free to add to this table by submitting a pull request).
#
# The below examples show how these keyword arguments affect the figure size.
#
# .. important::
#
#    The automatic figure size algorithm has the following notable properties:
#
#    * For very simple subplot grids (e.g., subplots created with the `ncols` and
#      `nrows` arguments), the arguments `refaspect`, `refwidth`, and `refheight`
#      effectively apply to every subplot in the figure -- not just the
#      reference subplot.
#    * When the reference subplot `aspect ratio
#      <https://matplotlib.org/stable/examples/pylab_examples/equal_aspect_ratio.html>`__
#      has been fixed (e.g., using ``ax.set_aspect(1)``) or is set to ``'equal'`` (as
#      with :ref:`geographic projections <ug_geo>` and `~proplot.axes.PlotAxes.imshow`
#      images), the fixed aspect ratio is used and the `~proplot.ui.subplots`
#      `refaspect` parameter is ignored. This is critical for getting the figure
#      size right when working with grids of images and geographic projections.
#    * The physical widths of `~proplot.axes.Axes.colorbar`\ s and
#      `~proplot.axes.Axes.panel`\ s are always preserved during figure resizing.
#      ProPlot specifies their widths in physical units to help avoid colorbars
#      and panels that look "too skinny" or "too fat".

# %%
import proplot as pplt
import numpy as np

# Grid of images (note the square pixels)
state = np.random.RandomState(51423)
colors = np.tile(state.rand(8, 12, 1), (1, 1, 3))
fig, axs = pplt.subplots(ncols=3, nrows=2, refwidth=1.7)
fig.format(suptitle='Auto figure size for grid of images')
for ax in axs:
    ax.imshow(colors)

# Grid of cartopy projections
fig, axs = pplt.subplots(ncols=2, nrows=3, proj='robin')
axs.format(land=True, landcolor='k')
fig.format(suptitle='Auto figure size for grid of cartopy projections')


# %%
import proplot as pplt
pplt.rc.update(grid=False, titleloc='uc', titleweight='bold', titlecolor='red9')

# Change the reference subplot width
suptitle = 'Effect of subplot width on figure size'
for refwidth in ('3cm', '5cm'):
    fig, axs = pplt.subplots(ncols=2, refwidth=refwidth,)
    axs[0].format(title=f'refwidth = {refwidth}', suptitle=suptitle)

# Change the reference subplot aspect ratio
suptitle = 'Effect of subplot aspect ratio on figure size'
for refaspect in (1, 2):
    fig, axs = pplt.subplots(ncols=2, refwidth=1.6, refaspect=refaspect)
    axs[0].format(title=f'refaspect = {refaspect}', suptitle=suptitle)

# Change the reference subplot
suptitle = 'Effect of reference subplot on figure size'
for ref in (1, 2):  # with different width ratios
    fig, axs = pplt.subplots(ncols=3, wratios=(3, 2, 2), ref=ref, refwidth=1.1)
    axs[ref - 1].format(title='reference', suptitle=suptitle)
for ref in (1, 2):  # with complex subplot grid
    fig, axs = pplt.subplots([[1, 2], [1, 3]], refnum=ref, refwidth=1.8)
    axs[ref - 1].format(title='reference', suptitle=suptitle)

pplt.rc.reset()

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_tight:
#
# Automatic spaces
# ----------------
#
# By default, ProPlot automatically determines the suitable space between
# subplots using a tight layout algorithm. This algorithm automatically
# expands or contracts the space between subplots to accommodate labels.
# It can be disabled by passing ``tight=False`` to `~proplot.ui.subplots`
# or setting :rcraw:`subplots.tight` to ``False``. In contrast to
# `matplotlib's tight layout algorithm
# <https://matplotlib.org/stable/tutorials/intermediate/tight_layout_guide.html>`__,
# ProPlot's algorithm may change the figure size to accommodate the correct
# spacing and permits variable spacing between subsequent subplot rows and
# columns (see `proplot.gridspec.GridSpec` for details).
#
# The tight layout algorithm can also be completely or partly overridden. When
# you pass any of the spacing arguments `left`, `right`, `top`, `bottom`,
# `wspace`, or `hspace` to `~proplot.ui.figure`, `~proplot.ui.subplots`, or
# `~proplot.gridspec.GridSpec`, that value is always respected. For example:
#
# * ``left=2`` fixes the left margin at 2 em-widths, while the right,
#   bottom, and top margin widths are determined by the tight layout algorithm.
# * ``wspace=1`` fixes the spaces between subplot columns at 1 em-width, while the
#   spaces between subplot rows are determined by the tight layout algorithm.
# * ``wspace=(3, None)`` fixes the space between the first two columns of
#   a three-column plot at 3 em-widths, while the space between the second two
#   columns is determined by the tight layout algorithm.
#
# Alternatively, the padding used by the tight layout algorithm (rather than the
# absolute spaces between subplot edges) can be changed by passing `outerpad`,
# `innerpad`, or `panelpad` to `~proplot.ui.figure` or `~proplot.ui.subplots`.
# This padding can be set locally by passing an array of values to `wpad`
# and `hpad` (analogous to `wspace` and `hspace`), or by passing the `pad`
# keyword when creating :ref:`panel axes <ug_panels>` or :ref:`outer
# colorbars and legends <ug_cbars_axes>` (analogous to `space`). Finally,
# to constrain the tight layout algorithm to produce equal spacing between
# main subplot rows and columns, you can pass ``wequal=True``, ``hequal=True``
# or ``equal=True`` to `~proplot.ui.figure` or `~proplot.ui.subplots` (note that
# equal spacing is the default behavior when tight layout is disabled).

# All the spacing parameters described above can be specified with a
# :ref:`unit string <ug_units>` interpreted by `~proplot.utils.units`.
# The default unit assumed for numeric arguments is an "em-width" (i.e., a
# :rcraw:`font.size` width -- see the :ref:`units table <units_table>` for details).

# %%
import proplot as pplt

# Stress test of the tight layout algorithm
# Add large labels along the edge of one subplot
for equal, descrip in enumerate(('variable', 'equal')):
    fig, axs = pplt.subplots(
        nrows=3, ncols=3, refwidth=1.1, share=False, equal=bool(equal)
    )
    axs[1].format(
        xlabel='xlabel\nxlabel',
        ylabel='ylabel\nylabel\nylabel\nylabel'
    )
    axs.format(
        grid=False,
        toplabels=('Column 1', 'Column 2', 'Column 3'),
        leftlabels=('Row 1', 'Row 2', 'Row 3'),
        suptitle=f'Tight layout with {descrip} row-column spacing',
    )

# %%
import proplot as pplt

# Stress test of the tight layout algorithm
# This time override the algorithm between selected subplot rows/columns
fig, axs = pplt.subplots(
    ncols=4, nrows=3, refwidth=1.1, span=False,
    bottom='5em', right='5em',  # margin spacing overrides
    wspace=(0, 0, None), hspace=(0, None),  # column and row spacing overrides
)
axs.format(
    grid=False,
    xlocator=1, ylocator=1, tickdir='inout',
    xlim=(-1.5, 1.5), ylim=(-1.5, 1.5),
    suptitle='Tight layout with user overrides',
    toplabels=('Column 1', 'Column 2', 'Column 3', 'Column 4'),
    leftlabels=('Row 1', 'Row 2', 'Row 3'),
)
axs[0, :].format(xtickloc='top')
axs[2, :].format(xtickloc='both')
axs[:, 1].format(ytickloc='neither')
axs[:, 2].format(ytickloc='right')
axs[:, 3].format(ytickloc='both')
axs[-1, :].format(xlabel='xlabel', title='Title\nTitle\nTitle')
axs[:, 0].format(ylabel='ylabel')


# %% [raw] raw_mimetype="text/restructuredtext" tags=[]
# .. _ug_share:
#
# Axis sharing
# ------------
#
# Figures with lots of subplots often have :ref:`redundant labels <why_redundant>`.
# To help address this, `matplotlib.pyplot.subplots` includes the `sharex` and
# `sharey` keyword arguments that permit sharing axis limits, ticks, and tick labels
# between like rows and columns of subplots. ProPlot builds on this feature by...
#
# #. Automatically sharing axes between subplots and :ref:`panels <ug_panels>`
#    occupying the same rows or columns of the `~proplot.gridspec.GridSpec`. This
#    works for :ref:`aribtrarily complex subplot grids <ug_details>`. It also works
#    if subplots were generated one-by-one with `~proplot.figure.Figure.add_subplot`
#    rather than `~proplot.figure.Figure.subplots`. It is controlled by the `sharex`
#    and `sharey` keywords (default is :rc:`subplots.share`). You can use the
#    `share` keyword as a shorthand to set both `sharex` and `sharey`.
# #. Automatically sharing labels across subplots and :ref:`panels <ug_panels>`
#    with edges against the same row or column of the `~proplot.gridspec.GridSpec`.
#    This also works for complex grids and subplots generated one-by-one. It is
#    controlled by the `spanx` and `spany` keywords (default is :rc:`subplots.span`).
#    Use the `span` keyword as a shorthand to set both `spanx` and `spany`.
# #. Supporting five sharing "levels". These values can be passed to `sharex`,
#    `sharey`, or `share`, or assigned to :rcraw:`subplots.share`. The levels
#    are defined as follows:
#
#    * ``False`` or ``0``: Axis sharing is disabled.
#    * ``'labels'``, ``'labs'``, or ``1``: Axis labels are shared, but
#      nothing else. Labels will appear on the leftmost and bottommost subplots.
#    * ``'limits'``, ``'lims'``, or ``2``: Same as ``1``, but axis limits, axis
#      scales, and major and minor tick locations and formatting are also shared.
#    * ``True`` or ``3``: Same as ``2``, but axis tick labels are also shared.
#      Tick labels will appear on the leftmost and bottommost subplots.
#    * ``'all'`` or ``4``: Same as ``3``, but axis limits, axis scales, and
#      axis ticks are shared even between subplots not in the same row or column.
#
# The below examples demonstrate the effect of various axis and label sharing
# settings on the appearance of several subplot grids.

# %%
import proplot as pplt
import numpy as np
N = 50
M = 40
state = np.random.RandomState(51423)
cycle = pplt.Cycle('grays_r', M, left=0.1, right=0.8)
datas = []
for scale in (1, 3, 7, 0.2):
    data = scale * (state.rand(N, M) - 0.5).cumsum(axis=0)[N // 2:, :]
    datas.append(data)

# Same plot with different sharing and spanning settings
for i, share in enumerate((False, 'labels', 'limits', True)):
    fig = pplt.figure(refaspect=1, refwidth=1.06, sharey=share, spanx=i // 2)
    axs = fig.subplots(ncols=4)
    for ax, data in zip(axs, datas):
        on = ('off', 'on')[i // 2]
        ax.plot(data, cycle=cycle)
        ax.format(
            suptitle=f'Sharing mode {share!r} (level {i}) with spanning labels {on}',
            grid=False, xlabel='spanning axis', ylabel='shared axis'
        )

# %%
import proplot as pplt
import numpy as np
pplt.rc.reset()
pplt.rc.cycle = 'Set3'
state = np.random.RandomState(51423)

# Same plot with and without default sharing settings
titles = ('With redundant labels', 'Without redundant labels')
for b in (False, True):
    fig = pplt.figure(refwidth=1, share=b, span=b)
    axs = fig.subplots(nrows=4, ncols=4)
    for ax in axs:
        ax.plot((state.rand(100, 20) - 0.4).cumsum(axis=0))
    axs.format(
        abc=True, abcloc='ul', suptitle=titles[b],
        xlabel='xlabel', ylabel='ylabel',
        grid=False, xticks=25, yticks=5
    )


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_units:
#
# Physical units
# --------------
#
# ProPlot supports arbitrary physical units for controlling the figure
# `figwidth` and `figheight`, the reference subplot `refwidth` and `refheight`,
# the gridspec spacing and tight layout padding values `left`, `right`, `bottom`,
# `top`, `wspace`, `hspace`, `outerpad`, `innerpad`, `panelpad`, `wpad`, and `hpad`,
# the `~proplot.axes.Axes.panel_axes` and `~proplot.axes.Axes.colorbar` widths,
# and all applicable `~proplot.config.rc` settings (e.g., settings controlling
# legend spacing, label padding, and font size). This feature is powered by the
# `~proplot.utils.units` function.
#
# A table of acceptable physical units is found :ref:`here <units_table>`.
# They include centimeters, millimeters, pixels,
# `em-heights <https://en.wikipedia.org/wiki/Em_(typography)>`__,
# `en-heights <https://en.wikipedia.org/wiki/En_(typography)>`__,
# and `points <https://en.wikipedia.org/wiki/Point_(typography)>`__.
# The default physical unit (assumed when an argument is numeric) depends on the
# context. For legend and gridspec spaces, it is em-widths. For subplot and
# figure sizes, it is inches. For text padding and font sizes, it is points. See
# the relevant documentation in the :ref:`API reference <api>` for details.

# %%
import proplot as pplt
import numpy as np
with pplt.rc.context(fontsize='12px'):
    fig, axs = pplt.subplots(
        ncols=3, figwidth='15cm', figheight='3in',
        wspace=('10pt', '20pt'), right='10mm',
    )
    cmap = pplt.Colormap('Mono')
    cb = fig.colorbar(
        cmap, loc='b', extend='both', label='colorbar',
        width='2em', extendsize='3em', shrink=0.8,
    )
    pax = axs[2].panel_axes('r', width='5en')
axs.format(
    suptitle='Arguments with arbitrary units',
    xlabel='x axis', ylabel='y axis',
    xlim=(0, 1), ylim=(0, 1),
)
