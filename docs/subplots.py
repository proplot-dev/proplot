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
# .. _ug_layout:
#
# Subplots
# ========
#
# This section documents a variety of features related to proplot subplots,
# including a-b-c subplot labels, axis sharing between subplots, automatic
# "tight layout" spacing between subplots, and a unique feature where the figure
# width and/or height are automatically adjusted based on the subplot geometry.
#
# .. note::
#
#    Proplot only supports one `~proplot.gridspec.GridSpec` per figure
#    (see the section on :ref:`adding subplots <ug_subplot>`), and proplot
#    does not officially support the "nested" matplotlib structures
#    `~matplotlib.gridspec.GridSpecFromSubplotSpec` and `~matplotlib.figure.SubFigure`.
#    These restrictions have the advantage of 1) considerably simplifying the
#    :ref:`tight layout <ug_tight>` and :ref:`figure size <ug_autosize>`
#    algorithms and 2) reducing the ambiguity of :ref:`a-b-c label assignment <ug_abc>`
#    and :ref:`automatic axis sharing <ug_share>` between subplots. If you need the
#    features associated with "nested" matplotlib structures, some are reproducible
#    with proplot -- including :ref:`different spaces <ug_tight>` between distinct
#    subplot rows and columns and :ref:`different formatting <ug_subplotgrid>` for
#    distinct groups of subplots.
#
#
# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_abc:
#
# A-b-c labels
# ------------
#
# Proplot can quickly add "a-b-c" labels using the
# `~proplot.axes.Axes.number` assigned to each subplot.
# If you add subplots one-by-one with `~proplot.figure.Figure.add_subplot`, you can
# manually specify the number with the `number` keyword. By default, the subplot number
# is incremented by ``1`` each time you call `~proplot.figure.Figure.add_subplot`.
# If you draw all of your subplots at once with `~proplot.figure.Figure.add_subplots`,
# the numbers depend on the input arguments. If you
# :ref:`passed an array <ug_intro>`, the subplot numbers correspond to the numbers
# in the array. But if you used the `ncols` and `nrows` keyword arguments, the
# number order is row-major by default and can be switched to column-major by
# passing ``order='F'`` (note the number order also determines the list order in the
# `~proplot.gridspec.SubplotGrid` returned by `~proplot.figure.Figure.add_subplots`).
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

fig = pplt.figure(space=0, refwidth="10em")
axs = fig.subplots(nrows=3, ncols=3)
axs.format(
    abc="A.",
    abcloc="ul",
    xticks="null",
    yticks="null",
    facecolor="gray5",
    xlabel="x axis",
    ylabel="y axis",
    suptitle="A-b-c label offsetting, borders, and boxes",
)
axs[:3].format(abcloc="l", titleloc="l", title="Title")
axs[-3:].format(abcbbox=True)  # also disables abcborder
# axs[:-3].format(abcborder=True)  # this is already the default

# %%
import proplot as pplt

fig = pplt.figure(space=0, refwidth=0.7)
axs = fig.subplots(nrows=8, ncols=8)
axs.format(
    abc=True,
    abcloc="ur",
    xlabel="x axis",
    ylabel="y axis",
    xticks=[],
    yticks=[],
    suptitle="A-b-c label stress test",
)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_autosize:
#
# Figure width and height
# -----------------------
#
# Proplot automatically adjusts the figure width and height by default to
# respect the physical size of a "reference" subplot and the geometry of the
# `~proplot.figure.Figure.gridspec`. The "reference" subplot is the subplot whose
# `~proplot.axes.Axes.number` matches the `refnum` that was passed to
# `~proplot.figure.Figure` (the default `refnum` of ``1`` usually matches the subplot
# in the upper-left corner -- see :ref:`this section <ug_abc>` for more on subplot
# numbers). Alternatively, you can request a fixed figure width (height), and the
# algorithm will automatically adjusts the figure height (width) to respect
# the `~proplot.figure.Figure.gridspec` geometry.

# This algorithm is extremely powerful and generally produces more aesthetically
# pleasing subplot grids out-of-the-box, especially when they contain images or map
# projections (see below). It is constrained by the following `~proplot.figure.Figure`
# keyword arguments:
#
# * `refwidth` and `refheight` set the physical width and height of the reference
#   subplot (default is :rc:`subplots.refwidth`). If just the width (height) is
#   specified, then the height (width) is automatically adjusted to satisfy the
#   subplot spacing and the reference subplot aspect ratio `refaspect` (default
#   is ``1`` unless the data aspect ratio is fixed -- see below). If both the
#   width and height are specified, then `refaspect` is ignored.
# * `figwidth` and `figheight` set the physical width and height of the figure.
#   As in matplotlib, you can use `figsize` to set both at once. If just the width
#   (height) is specified, then the height (width) is automatically adjusted, just
#   like with `refwidth` and `refheight`. If both the width and height are specified
#   (e.g., using `figsize`), then `refaspect` is ignored and the figure size is fixed.
#   Note that `figwidth` and `figheight` always override `refwidth` and `refheight`.
# * `journal` sets the physical dimensions of the figure to meet requirements
#   for submission to an academic journal. For example, ``journal='nat1'`` sets
#   `figwidth` according to the `*Nature* standard for single-column figures
#   <https://www.nature.com/nature/for-authors/formatting-guide>`__ and
#   ``journal='aaas2'`` sets `figwidth` according to the
#   `*Science* standard for dual-column figures
#   <https://www.sciencemag.org/authors/instructions-preparing-initial-manuscript>`__.
#   See :ref:`this table <journal_table>` for the currently available journal
#   specifications (feel free to add to this list with a :ref:`pull request
#   <contrib_pr>`).
#
# The below examples demonstrate how different keyword arguments and
# subplot arrangements influence the figure size algorithm.
#
# .. important::
#
#    * If the `data aspect ratio
#      <https://matplotlib.org/stable/examples/pylab_examples/equal_aspect_ratio.html>`__
#      of the reference subplot is fixed (either due to calling
#      `~matplotlib.axes.Axes.set_aspect` or filling the subplot with a
#      :ref:`geographic projection <ug_geo>`, `~proplot.axes.PlotAxes.imshow`
#      plot, or `~proplot.axes.PlotAxes.heatmap` plot), then this is used as
#      the default value for the reference aspect ratio `refaspect`. This helps
#      minimize excess space between grids of subplots with fixed aspect ratios.
#    * For the simplest subplot grids (e.g., those created by passing integers to
#      `~proplot.figure.Figure.add_subplot` or passing `ncols` or `nrows` to
#      `~proplot.figure.Figure.add_subplots`) the keyword arguments `refaspect`,
#      `refwidth`, and `refheight` effectively apply to every subplot in the
#      figure -- not just the reference subplot.
#    * The physical widths of proplot `~proplot.axes.Axes.colorbar`\ s and
#      `~proplot.axes.Axes.panel`\ s are always independent of the figure size.
#      `~proplot.gridspec.GridSpec` specifies their widths in physical units to help
#      users avoid drawing colorbars and panels that look "too skinny" or "too fat"
#      depending on the number of subplots in the figure.

# %%
import proplot as pplt
import numpy as np

# Grid of images (note the square pixels)
state = np.random.RandomState(51423)
colors = np.tile(state.rand(8, 12, 1), (1, 1, 3))
fig, axs = pplt.subplots(ncols=3, nrows=2, refwidth=1.7)
fig.format(suptitle="Auto figure dimensions for grid of images")
for ax in axs:
    ax.imshow(colors)

# Grid of cartopy projections
fig, axs = pplt.subplots(ncols=2, nrows=3, proj="robin")
axs.format(land=True, landcolor="k")
fig.format(suptitle="Auto figure dimensions for grid of cartopy projections")


# %%
import proplot as pplt

pplt.rc.update(grid=False, titleloc="uc", titleweight="bold", titlecolor="red9")

# Change the reference subplot width
suptitle = "Effect of subplot width on figure size"
for refwidth in ("3cm", "5cm"):
    fig, axs = pplt.subplots(
        ncols=2,
        refwidth=refwidth,
    )
    axs[0].format(title=f"refwidth = {refwidth}", suptitle=suptitle)

# Change the reference subplot aspect ratio
suptitle = "Effect of subplot aspect ratio on figure size"
for refaspect in (1, 2):
    fig, axs = pplt.subplots(ncols=2, refwidth=1.6, refaspect=refaspect)
    axs[0].format(title=f"refaspect = {refaspect}", suptitle=suptitle)

# Change the reference subplot
suptitle = "Effect of reference subplot on figure size"
for ref in (1, 2):  # with different width ratios
    fig, axs = pplt.subplots(ncols=3, wratios=(3, 2, 2), ref=ref, refwidth=1.1)
    axs[ref - 1].format(title="reference", suptitle=suptitle)
for ref in (1, 2):  # with complex subplot grid
    fig, axs = pplt.subplots([[1, 2], [1, 3]], refnum=ref, refwidth=1.8)
    axs[ref - 1].format(title="reference", suptitle=suptitle)

pplt.rc.reset()

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_tight:
#
# Spacing and tight layout
# ------------------------
#
# Proplot automatically adjusts the spacing between subplots
# by default to accomadate labels using its own `"tight layout" algorithm
# <https://matplotlib.org/stable/tutorials/intermediate/tight_layout_guide.html>`__.
# In contrast to matplotlib's algorithm, proplot's algorithm can :ref:`change the
# figure size <ug_autosize>` and permits variable spacing between each subplot
# row and column (see `proplot.gridspec.GridSpec` for details).
# This algorithm can be disabled entirely by passing ``tight=False`` to
# `~proplot.figure.Figure` or by setting :rcraw:`subplots.tight` to ``False``, or
# it can be partly overridden by passing any of the spacing arguments `left`, `right`,
# `top`, `bottom`, `wspace`, or `hspace` to `~proplot.figure.Figure` or
# `~proplot.gridspec.GridSpec`. For example:
#
# * ``left=2`` fixes the left margin at 2 em-widths, while the right,
#   bottom, and top margin widths are determined by the tight layout algorithm.
# * ``wspace=1`` fixes the space between subplot columns at 1 em-width, while the
#   space between subplot rows is determined by the tight layout algorithm.
# * ``wspace=(3, None)`` fixes the space between the first two columns of
#   a three-column plot at 3 em-widths, while the space between the second two
#   columns is determined by the tight layout algorithm.
#
# The padding between the tight layout extents (rather than the absolute spaces
# between subplot edges) can also be changed by passing `outerpad`, `innerpad`,
# or `panelpad` to `~proplot.figure.Figure` or `~proplot.gridspec.GridSpec`.
# This padding can be set locally by passing an array of values to `wpad`
# and `hpad` (analogous to `wspace` and `hspace`), or by passing the `pad`
# keyword when creating :ref:`panel axes <ug_panels>` or :ref:`outer
# colorbars or legends <ug_guides_loc>` (analogous to `space`).
#
# All the subplot spacing arguments can be specified with a
# :ref:`unit string <ug_units>` interpreted by `~proplot.utils.units`.
# The default unit assumed for numeric arguments is an "em-width" (i.e., a
# :rcraw:`font.size` width -- see the :ref:`units table <units_table>` for details).
#
# .. note::

#    The core behavior of the tight layout algorithm can be modified with a few
#    keyword arguments and settings. Using ``wequal=True``, ``hequal=True``, or
#    ``equal=True`` (or setting :rcraw:`subplots.equalspace` to ``True``) constrains
#    the tight layout algorithm to produce equal spacing between main subplot columns
#    or rows (note that equal spacing is the default behavior when tight layout is
#    disabled). Similarly, using ``wgroup=False``, ``hgroup=False``, or ``group=False``
#    (or setting :rcraw:`subplots.groupspace` to ``False``) disables the default
#    behavior of only comparing subplot extent between adjacent subplot "groups"
#    and instead compares subplot extents across entire columns and rows
#    (note the spacing between the first and second row in the below example).

# %%
import proplot as pplt

# Stress test of the tight layout algorithm
# This time override the algorithm between selected subplot rows/columns
fig, axs = pplt.subplots(
    ncols=4,
    nrows=3,
    refwidth=1.1,
    span=False,
    bottom="5em",
    right="5em",  # margin spacing overrides
    wspace=(0, 0, None),
    hspace=(0, None),  # column and row spacing overrides
)
axs.format(
    grid=False,
    xlocator=1,
    ylocator=1,
    tickdir="inout",
    xlim=(-1.5, 1.5),
    ylim=(-1.5, 1.5),
    suptitle="Tight layout with user overrides",
    toplabels=("Column 1", "Column 2", "Column 3", "Column 4"),
    leftlabels=("Row 1", "Row 2", "Row 3"),
)
axs[0, :].format(xtickloc="top")
axs[2, :].format(xtickloc="both")
axs[:, 1].format(ytickloc="neither")
axs[:, 2].format(ytickloc="right")
axs[:, 3].format(ytickloc="both")
axs[-1, :].format(xlabel="xlabel", title="Title\nTitle\nTitle")
axs[:, 0].format(ylabel="ylabel")


# %%
import proplot as pplt

# Stress test of the tight layout algorithm
# Add large labels along the edge of one subplot
equals = [("unequal", False), ("unequal", False), ("equal", True)]
groups = [("grouped", True), ("ungrouped", False), ("grouped", True)]
for (name1, equal), (name2, group) in zip(equals, groups):
    suffix = " (default)" if group and not equal else ""
    suptitle = f'Tight layout with "{name1}" and "{name2}" row-column spacing{suffix}'
    fig, axs = pplt.subplots(
        nrows=3,
        ncols=3,
        refwidth=1.1,
        share=False,
        equal=equal,
        group=group,
    )
    axs[1].format(xlabel="xlabel\nxlabel", ylabel="ylabel\nylabel\nylabel\nylabel")
    axs[3:6:2].format(
        title="Title\nTitle",
        titlesize="med",
    )
    axs.format(
        grid=False,
        toplabels=("Column 1", "Column 2", "Column 3"),
        leftlabels=("Row 1", "Row 2", "Row 3"),
        suptitle=suptitle,
    )

# %% [raw] raw_mimetype="text/restructuredtext" tags=[]
# .. _ug_share:
#
# Axis label sharing
# ------------------
#
# Figures with lots of subplots often have :ref:`redundant labels <why_redundant>`.
# To help address this, the matplotlib command `matplotlib.pyplot.subplots` includes
# `sharex` and `sharey` keywords that permit sharing axis limits and ticks between
# like rows and columns of subplots. Proplot builds on this feature by:
#
# #. Automatically sharing axes between subplots and :ref:`panels <ug_panels>`
#    occupying the same rows or columns of the `~proplot.gridspec.GridSpec`. This
#    works for :ref:`aribtrarily complex subplot grids <ug_layout>`. It also works
#    for subplots generated one-by-one with `~proplot.figure.Figure.add_subplot`
#    rather than `~proplot.figure.Figure.subplots`. It is controlled by the `sharex`
#    and `sharey` `~proplot.figure.Figure` keywords (default is :rc:`subplots.share`).
#    Use the `share` keyword as a shorthand to set both `sharex` and `sharey`.
# #. Automatically sharing labels across subplots and :ref:`panels <ug_panels>`
#    with edges along the same row or column of the `~proplot.gridspec.GridSpec`.
#    This also works for complex subplot grids and subplots generated one-by-one.
#    It is controlled by the `spanx` and `spany` `~proplot.figure.Figure`
#    keywords (default is :rc:`subplots.span`). Use the `span` keyword
#    as a shorthand to set both `spanx` and `spany`. Note that unlike
#    `~matplotlib.figure.Figure.supxlabel` and `~matplotlib.figure.Figure.supylabel`,
#    these labels are aligned between gridspec edges rather than figure edges.
# #. Supporting five sharing "levels". These values can be passed to `sharex`,
#    `sharey`, or `share`, or assigned to :rcraw:`subplots.share`. The levels
#    are defined as follows:
#
#    * ``False`` or ``0``: Axis sharing is disabled.
#    * ``'labels'``, ``'labs'``, or ``1``: Axis labels are shared, but
#      nothing else. Labels will appear on the leftmost and bottommost subplots.
#    * ``'limits'``, ``'lims'``, or ``2``: Same as ``1``, but axis limits, axis
#      scales, and major and minor tick locations and formatting are also shared.
#    * ``True`` or ``3`` (default): Same as ``2``, but axis tick labels are also
#      shared. Tick labels will appear on the leftmost and bottommost subplots.
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
cycle = pplt.Cycle("grays_r", M, left=0.1, right=0.8)
datas = []
for scale in (1, 3, 7, 0.2):
    data = scale * (state.rand(N, M) - 0.5).cumsum(axis=0)[N // 2 :, :]
    datas.append(data)

# Plots with different sharing and spanning settings
# Note that span=True and share=True are the defaults
spans = (False, False, True, True)
shares = (False, "labels", "limits", True)
for i, (span, share) in enumerate(zip(spans, shares)):
    fig = pplt.figure(refaspect=1, refwidth=1.06, spanx=span, sharey=share)
    axs = fig.subplots(ncols=4)
    for ax, data in zip(axs, datas):
        on = ("off", "on")[int(span)]
        ax.plot(data, cycle=cycle)
        ax.format(
            grid=False,
            xlabel="spanning axis",
            ylabel="shared axis",
            suptitle=f"Sharing mode {share!r} (level {i}) with spanning labels {on}",
        )

# %%
import proplot as pplt
import numpy as np

state = np.random.RandomState(51423)

# Plots with minimum and maximum sharing settings
# Note that all x and y axis limits and ticks are identical
spans = (False, True)
shares = (False, "all")
titles = ("Minimum sharing", "Maximum sharing")
for span, share, title in zip(spans, shares, titles):
    fig = pplt.figure(refwidth=1, span=span, share=share)
    axs = fig.subplots(nrows=4, ncols=4)
    for ax in axs:
        data = (state.rand(100, 20) - 0.4).cumsum(axis=0)
        ax.plot(data, cycle="Set3")
    axs.format(
        abc=True,
        abcloc="ul",
        suptitle=title,
        xlabel="xlabel",
        ylabel="ylabel",
        grid=False,
        xticks=25,
        yticks=5,
    )


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_units:
#
# Physical units
# --------------
#
# Proplot supports arbitrary physical units for controlling the figure
# `figwidth` and `figheight`; the reference subplot `refwidth` and `refheight`;
# the gridspec spacing and tight layout padding keywords `left`, `right`, `bottom`,
# `top`, `wspace`, `hspace`, `outerpad`, `innerpad`, `panelpad`, `wpad`, and `hpad`;
# the `~proplot.axes.Axes.colorbar` and `~proplot.axes.Axes.panel` widths;
# various `~proplot.axes.Axes.legend` spacing and padding arguments; various
# `~proplot.axes.Axes.format` font size and padding arguments; the line width and
# marker size arguments passed to `~proplot.axes.PlotAxes` commands; and all
# applicable `~proplot.config.rc` settings, e.g. :rcraw:`subplots.refwidth`,
# :rcraw:`legend.columnspacing`, and :rcraw:`axes.labelpad`. This feature is
# powered by the physical units engine `~proplot.utils.units`.
#
# When one of these keyword arguments is numeric, a default physical unit is
# used. For subplot and figure sizes, the defult unit is inches. For gridspec and
# legend spaces, the default unit is `em-widths
# <https://en.wikipedia.org/wiki/Em_(typography)>`__.
# For font sizes, text padding, and
# line widths, the default unit is
# `points <https://en.wikipedia.org/wiki/Point_(typography)>`__.
# See the relevant documentation in the :ref:`API reference <api>` for details.
# A table of acceptable physical units is found :ref:`here <units_table>`
# -- they include centimeters, millimeters, pixels,
# `em-widths <https://en.wikipedia.org/wiki/Em_(typography)>`__,
# `en-heights <https://en.wikipedia.org/wiki/En_(typography)>`__,
# and `points <https://en.wikipedia.org/wiki/Point_(typography)>`__.

# %%
import proplot as pplt
import numpy as np

with pplt.rc.context(fontsize="12px"):  # depends on rc['figure.dpi']
    fig, axs = pplt.subplots(
        ncols=3,
        figwidth="15cm",
        figheight="3in",
        wspace=("10pt", "20pt"),
        right="10mm",
    )
    cb = fig.colorbar(
        "Mono",
        loc="b",
        extend="both",
        label="colorbar",
        width="2em",
        extendsize="3em",
        shrink=0.8,
    )
    pax = axs[2].panel_axes("r", width="5en")
axs.format(
    suptitle="Arguments with arbitrary units",
    xlabel="x axis",
    ylabel="y axis",
)
