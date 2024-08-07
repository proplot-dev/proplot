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
# .. _ug_guides:
#
# Colorbars and legends
# =====================
#
# Proplot includes some useful changes to the matplotlib API that make
# working with colorbars and legends :ref:`easier <why_colorbars_legends>`.
# Notable features include "inset" colorbars, "outer" legends,
# on-the-fly colorbars and legends, colorbars built from artists,
# and row-major and centered-row legends.


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_guides_loc:
#
# Outer and inset locations
# -------------------------
#
# Matplotlib supports drawing "inset" legends and "outer" colorbars using the `loc`
# and `location` keyword arguments. However, "outer" legends are only
# posssible using the somewhat opaque `bbox_to_anchor` keyword (see `here
# <https://matplotlib.org/stable/tutorials/intermediate/legend_guide.html#legend-location>`__)
# and "inset" colorbars are not possible without manually creating and positioning
# the associated axes. Proplot tries to improve this behavior:
#
# * `proplot.axes.Axes.legend` can draw both "inset" legends when you request an inset
#   location (e.g., ``loc='upper right'`` or the shorthand ``loc='ur'``) and "outer"
#   legends along a subplot edge when you request a :ref:`side location <legend_table>`
#   (e.g., ``loc='right'`` or the shorthand ``loc='r'``). If you draw multiple legends
#   or colorbars on one side, they are "stacked" on top of each other. Unlike using
#   `bbox_to_anchor`, the "outer" legend position is adjusted automatically when the
#   :ref:`tight layout algorithm <ug_tight>` is active.
# * Proplot adds the axes command `proplot.axes.Axes.colorbar`,
#   analogous to `proplot.axes.Axes.legend` and equivalent to
#   calling `proplot.figure.Figure.colorbar` with an `ax` keyword.
#   `~proplot.axes.Axes.colorbar` can draw both "outer" colorbars when you request
#   a side location (e.g., ``loc='right'`` or the shorthand ``loc='r'``) and "inset"
#   colorbars when you request an :ref:`inset location <colorbar_table>`
#   (e.g., ``loc='upper right'`` or the shorthand ``loc='ur'``). Inset
#   colorbars have optional background "frames" that can be configured
#   with various `~proplot.axes.Axes.colorbar` keywords.

# `~proplot.axes.Axes.colorbar` and `~proplot.axes.Axes.legend` also both accept
# `space` and `pad` keywords. `space` controls the absolute separation of the
# "outer" colorbar or legend from the parent subplot edge and `pad` controls the
# :ref:`tight layout <ug_tight>` padding relative to the subplot's tick and axis labels
# (or, for "inset" locations, the padding between the subplot edge and the inset frame).
# The below example shows a variety of arrangements of "outer" and "inset"
# colorbars and legends.
#
# .. important::
#
#    Unlike matplotlib, proplot adds "outer" colorbars and legends by allocating
#    new rows and columns in the `~proplot.gridspec.GridSpec` rather than
#    "stealing" space from the parent subplot (note that subsequently indexing
#    the `~proplot.gridspec.GridSpec` will ignore the slots allocated for
#    colorbars and legends). This approach means that "outer" colorbars and
#    legends :ref:`do not affect subplot aspect ratios <ug_autosize>`
#    and :ref:`do not affect subplot spacing <ug_tight>`, which lets
#    proplot avoid relying on complicated `"constrained layout" algorithms
#    <https://matplotlib.org/stable/tutorials/intermediate/constrainedlayout_guide.html>`__
#    and tends to improve the appearance of figures with even the most
#    complex arrangements of subplots, colorbars, and legends.

# %%
import proplot as pplt
import numpy as np

state = np.random.RandomState(51423)
fig = pplt.figure(share=False, refwidth=2.3)

# Colorbars
ax = fig.subplot(121, title="Axes colorbars")
data = state.rand(10, 10)
m = ax.heatmap(data, cmap="dusk")
ax.colorbar(m, loc="r")
ax.colorbar(m, loc="t")  # title is automatically adjusted
ax.colorbar(m, loc="ll", label="colorbar label")  # inset colorbar demonstration

# Legends
ax = fig.subplot(122, title="Axes legends", titlepad="0em")
data = (state.rand(10, 5) - 0.5).cumsum(axis=0)
hs = ax.plot(data, lw=3, cycle="ggplot", labels=list("abcde"))
ax.legend(loc="ll", label="legend label")  # automatically infer handles and labels
ax.legend(hs, loc="t", ncols=5, frame=False)  # automatically infer labels from handles
ax.legend(hs, list("jklmn"), loc="r", ncols=1, frame=False)  # manually override labels
fig.format(
    abc=True,
    xlabel="xlabel",
    ylabel="ylabel",
    suptitle="Colorbar and legend location demo",
)

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_guides_plot:
#
# On-the-fly colorbars and legends
# --------------------------------
#
# In proplot, you can add colorbars and legends on-the-fly by supplying keyword
# arguments to various `~proplot.axes.PlotAxes` commands. To plot data and
# draw a colorbar or legend in one go, pass a location (e.g., ``colorbar='r'``
# or ``legend='b'``) to the plotting command (e.g., `~proplot.axes.PlotAxes.plot`
# or `~proplot.axes.PlotAxes.contour`). To pass keyword arguments to the colorbar
# and legend commands, use the `legend_kw` and `colorbar_kw` arguments (e.g.,
# ``legend_kw={'ncol': 3}``). Note that `~proplot.axes.Axes.colorbar` can also
# build colorbars from lists of arbitrary matplotlib artists, for example the
# lines generated by `~proplot.axes.PlotAxes.plot` or `~proplot.axes.PlotAxes.line`
# (see :ref:`below <ug_colorbars>`).
#
# .. note::
#
#    Specifying the same `colorbar` location with multiple plotting calls will have
#    a different effect depending on the plotting command. For :ref:`1D commands
#    <ug_1dplots>`, this will add each item to a "queue" used to build colorbars
#    from a list of artists. For :ref:`2D commands <ug_2dplots>`, this will "stack"
#    colorbars in outer locations, or replace existing colorbars in inset locations.
#    By contrast, specifying the same `legend` location will always add items to
#    the same legend rather than creating "stacks".

# %%
import proplot as pplt

labels = list("xyzpq")
state = np.random.RandomState(51423)
fig = pplt.figure(share=0, refwidth=2.3, suptitle="On-the-fly colorbar and legend demo")

# Legends
data = (state.rand(30, 10) - 0.5).cumsum(axis=0)
ax = fig.subplot(121, title="On-the-fly legend")
ax.plot(  # add all at once
    data[:, :5],
    lw=2,
    cycle="Reds1",
    cycle_kw={"ls": ("-", "--"), "left": 0.1},
    labels=labels,
    legend="b",
    legend_kw={"title": "legend title"},
)
for i in range(5):
    ax.plot(  # add one-by-one
        data[:, 5 + i],
        label=labels[i],
        linewidth=2,
        cycle="Blues1",
        cycle_kw={"N": 5, "ls": ("-", "--"), "left": 0.1},
        colorbar="ul",
        colorbar_kw={"label": "colorbar from lines"},
    )

# Colorbars
ax = fig.subplot(122, title="On-the-fly colorbar")
data = state.rand(8, 8)
ax.contourf(
    data,
    cmap="Reds1",
    extend="both",
    colorbar="b",
    colorbar_kw={"length": 0.8, "label": "colorbar label"},
)
ax.contour(
    data,
    color="gray7",
    lw=1.5,
    label="contour",
    legend="ul",
    legend_kw={"label": "legend from contours"},
)

# %%
import proplot as pplt
import numpy as np

N = 10
state = np.random.RandomState(51423)
fig, axs = pplt.subplots(
    nrows=2,
    share=False,
    refwidth="55mm",
    panelpad="1em",
    suptitle="Stacked colorbars demo",
)

# Repeat for both axes
args1 = (0, 0.5, 1, 1, "grays", 0.5)
args2 = (0, 0, 0.5, 0.5, "reds", 1)
args3 = (0.5, 0, 1, 0.5, "blues", 2)
for j, ax in enumerate(axs):
    ax.format(xlabel="data", xlocator=np.linspace(0, 0.8, 5), title=f"Subplot #{j+1}")
    for i, (x0, y0, x1, y1, cmap, scale) in enumerate((args1, args2, args3)):
        if j == 1 and i == 0:
            continue
        data = state.rand(N, N) * scale
        x, y = np.linspace(x0, x1, N + 1), np.linspace(y0, y1, N + 1)
        m = ax.pcolormesh(x, y, data, cmap=cmap, levels=np.linspace(0, scale, 11))
        ax.colorbar(m, loc="l", label=f"dataset #{i + 1}")


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_guides_multi:
#
# Figure-wide colorbars and legends
# ---------------------------------
#
# In proplot, colorbars and legends can be added to the edge of figures using the
# figure methods `proplot.figure.Figure.colorbar` and `proplot.figure.Figure.legend`.
# These methods align colorbars and legends between the edges
# of the `~proplot.figure.Figure.gridspec` rather than the figure.
# As with :ref:`axes colorbars and legends <ug_guides_loc>`, if you
# draw multiple colorbars or legends on the same side, they are stacked on
# top of each other. To draw a colorbar or legend alongside particular row(s) or
# column(s) of the subplot grid, use the `row`, `rows`, `col`, or `cols` keyword
# arguments. You can pass an integer to draw the colorbar or legend beside a
# single row or column (e.g., ``fig.colorbar(m, row=1)``), or pass a tuple to
# draw the colorbar or legend along a range of rows or columns
# (e.g., ``fig.colorbar(m, rows=(1, 2))``). The space separation between the subplot
# grid edge and the colorbars or legends can be controlled with the `space` keyword,
# and the tight layout padding can be controlled with the `pad` keyword.

# %%
import proplot as pplt
import numpy as np

state = np.random.RandomState(51423)
fig, axs = pplt.subplots(ncols=3, nrows=3, refwidth=1.4)
for ax in axs:
    m = ax.pcolormesh(
        state.rand(20, 20), cmap="grays", levels=np.linspace(0, 1, 11), extend="both"
    )
fig.format(
    suptitle="Figure colorbars and legends demo",
    abc="a.",
    abcloc="l",
    xlabel="xlabel",
    ylabel="ylabel",
)
fig.colorbar(m, label="column 1", ticks=0.5, loc="b", col=1)
fig.colorbar(m, label="columns 2 and 3", ticks=0.2, loc="b", cols=(2, 3))
fig.colorbar(m, label="stacked colorbar", ticks=0.1, loc="b", minorticks=0.05)
fig.colorbar(m, label="colorbar with length <1", ticks=0.1, loc="r", length=0.7)

# %%
import proplot as pplt
import numpy as np

state = np.random.RandomState(51423)
fig, axs = pplt.subplots(
    ncols=2, nrows=2, order="F", refwidth=1.7, wspace=2.5, share=False
)

# Plot data
data = (state.rand(50, 50) - 0.1).cumsum(axis=0)
for ax in axs[:2]:
    m = ax.contourf(data, cmap="grays", extend="both")
hs = []
colors = pplt.get_colors("grays", 5)
for abc, color in zip("ABCDEF", colors):
    data = state.rand(10)
    for ax in axs[2:]:
        (h,) = ax.plot(data, color=color, lw=3, label=f"line {abc}")
    hs.append(h)

# Add colorbars and legends
fig.colorbar(m, length=0.8, label="colorbar label", loc="b", col=1, locator=5)
fig.colorbar(m, label="colorbar label", loc="l")
fig.legend(hs, ncols=2, center=True, frame=False, loc="b", col=2)
fig.legend(hs, ncols=1, label="legend label", frame=False, loc="r")
fig.format(abc="A", abcloc="ul", suptitle="Figure colorbars and legends demo")
for ax, title in zip(axs, ("2D {} #1", "2D {} #2", "Line {} #1", "Line {} #2")):
    ax.format(xlabel="xlabel", title=title.format("dataset"))


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_colorbars:
#
# Added colorbar features
# -----------------------
#
# The `proplot.axes.Axes.colorbar` and `proplot.figure.Figure.colorbar` commands are
# somehwat more flexible than their matplotlib counterparts. The following core
# features are unique to proplot:

# * Calling ``colorbar`` with a list of `~matplotlib.artist.Artist`\ s,
#   a `~matplotlib.colors.Colormap` name or object, or a list of colors
#   will build the required `~matplotlib.cm.ScalarMappable` on-the-fly. Lists
#   of `~matplotlib.artist.Artists`\ s are used when you use the `colorbar`
#   keyword with :ref:`1D commands <ug_1dplots>` like `~proplot.axes.PlotAxes.plot`.
# * The associated :ref:`colormap normalizer <ug_norm>` can be specified with the
#   `vmin`, `vmax`, `norm`, and `norm_kw` keywords. The `~proplot.colors.DiscreteNorm`
#   levels can be specified with `values`, or proplot will infer them from the
#   `~matplotlib.artist.Artist` labels (non-numeric labels will be applied to
#   the colorbar as tick labels). This can be useful for labeling discrete plot
#   elements that bear some numeric relationship to each other.
#
# Proplot also includes improvements for adding ticks and tick labels to colorbars.
# Similar to `proplot.axes.CartesianAxes.format`, you can flexibly specify
# major tick locations, minor tick locations, and major tick labels using the
# `locator`, `minorlocator`, `formatter`, `ticks`, `minorticks`, and `ticklabels`
# keywords. These arguments are passed through the `~proplot.constructor.Locator` and
# `~proplot.constructor.Formatter` :ref:`constructor functions <why_constructor>`.
# Unlike matplotlib, the default ticks for :ref:`discrete colormaps <ug_discrete>`
# are restricted based on the axis length using `~proplot.ticker.DiscreteLocator`.
# You can easily toggle minor ticks using ``tickminor=True``.
#
# Similar to :ref:`axes panels <ug_panels>`, the geometry of proplot colorbars is
# specified with :ref:`physical units <ug_units>` (this helps avoid the common issue
# where colorbars appear "too skinny" or "too fat" and preserves their appearance
# when the figure size changes). You can specify the colorbar width locally using the
# `width` keyword or globally using the :rcraw:`colorbar.width` setting (for outer
# colorbars) and the :rcraw:`colorbar.insetwidth` setting (for inset colorbars).
# Similarly, you can specify the colorbar length locally with the `length` keyword or
# globally using the :rcraw:`colorbar.insetlength` setting. The outer colorbar length
# is always relative to the subplot grid and always has a default of ``1``. You
# can also specify the size of the colorbar "extensions" in physical units rather
# than relative units using the `extendsize` keyword rather than matplotlib's
# `extendfrac`. The default `extendsize` values are :rcraw:`colorbar.extend` (for
# outer colorbars) and :rcraw:`colorbar.insetextend` (for inset colorbars).
# See `~proplot.axes.Axes.colorbar` for details.

# %%
import proplot as pplt
import numpy as np

fig = pplt.figure(share=False, refwidth=2)

# Colorbars from lines
ax = fig.subplot(121)
state = np.random.RandomState(51423)
data = 1 + (state.rand(12, 10) - 0.45).cumsum(axis=0)
cycle = pplt.Cycle("algae")
hs = ax.line(
    data,
    lw=4,
    cycle=cycle,
    colorbar="lr",
    colorbar_kw={"length": "8em", "label": "line colorbar"},
)
ax.colorbar(hs, loc="t", values=np.arange(0, 10), label="line colorbar", ticks=2)

# Colorbars from a mappable
ax = fig.subplot(122)
m = ax.contourf(data.T, extend="both", cmap="algae", levels=pplt.arange(0, 3, 0.5))
fig.colorbar(
    m, loc="r", length=1, label="interior ticks", tickloc="left"  # length is relative
)
ax.colorbar(
    m,
    loc="ul",
    length=6,  # length is em widths
    label="inset colorbar",
    tickminor=True,
    alpha=0.5,
)
fig.format(
    suptitle="Colorbar formatting demo",
    xlabel="xlabel",
    ylabel="ylabel",
    titleabove=False,
)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_legends:
#
# Added legend features
# ---------------------
#
# The `proplot.axes.Axes.legend` and `proplot.figure.Figure.legend` commands are
# somewhat more flexible than their matplotlib counterparts. The following core
# features are the same as matplotlib:

# * Calling ``legend`` without positional arguments will
#   automatically fill the legend with the labeled artist in the
#   the parent axes (when using `proplot.axes.Axes.legend`) or
#   or the parent figure (when using `proplot.figure.Figure.legend`).
# * Legend labels can be assigned early by calling plotting comamnds with
#   the `label` keyword (e.g., ``ax.plot(..., label='label')``) or on-the-fly by
#   passing two positional arguments to ``legend`` (where the first argument is the
#   "handle" list and the second is the "label" list).

# The following core features are unique to proplot:

# * Legend labels can be assigned for each column of a
#   :ref:`2D array passed to a 1D plotting command <ug_1dstd>`
#   using the `labels` keyword (e.g., ``labels=['label1', 'label2', ...]``).
# * Legend labels can be assigned to `~matplotlib.contour.ContourSet`\ s by passing
#   the `label` keyword to a contouring command (e.g., `~proplot.axes.PlotAxes.contour`
#   or `~proplot.axes.PlotAxes.contourf`).
# * A "handle" list can be passed to ``legend`` as the sole
#   positional argument and the labels will be automatically inferred
#   using `~matplotlib.artist.Artist.get_label`. Valid "handles" include
#   `~matplotlib.lines.Line2D`\ s returned by `~proplot.axes.PlotAxes.plot`,
#   `~matplotlib.container.BarContainer`\ s returned by `~proplot.axes.PlotAxes.bar`,
#   and `~matplotlib.collections.PolyCollection`\ s
#   returned by `~proplot.axes.PlotAxes.fill_between`.
# * A composite handle can be created by grouping the "handle"
#   list objects into tuples (see this `matplotlib guide
#   <https://matplotlib.org/stable/tutorials/intermediate/legend_guide.html#legend-handlers>`__
#   for more on tuple groups). The associated label will be automatically
#   inferred from the objects in the group. If multiple distinct
#   labels are found then the group is automatically expanded.
#
# `proplot.axes.Axes.legend` and `proplot.figure.Figure.legend` include a few other
# useful features. To draw legends with centered rows, pass ``center=True`` or
# a list of lists of "handles" to ``legend`` (this stacks several single-row,
# horizontally centered legends and adds an encompassing frame behind them).
# To switch between row-major and column-major order for legend entries,
# use the `order` keyword (the default ``order='C'`` is row-major,
# unlike matplotlib's column-major ``order='F'``). To alphabetize the legend
# entries, pass ``alphabetize=True`` to ``legend``. To modify the legend handles
# (e.g., `~proplot.axes.PlotAxes.plot` or `~proplot.axes.PlotAxes.scatter` handles)
# pass the relevant properties like `color`, `linewidth`, or `markersize` to ``legend``
# (or use the `handle_kw` keyword). See `proplot.axes.Axes.legend` for details.

# %%
import proplot as pplt
import numpy as np

pplt.rc.cycle = "538"
fig, axs = pplt.subplots(ncols=2, span=False, share="labels", refwidth=2.3)
labels = ["a", "bb", "ccc", "dddd", "eeeee"]
hs1, hs2 = [], []

# On-the-fly legends
state = np.random.RandomState(51423)
for i, label in enumerate(labels):
    data = (state.rand(20) - 0.45).cumsum(axis=0)
    h1 = axs[0].plot(
        data,
        lw=4,
        label=label,
        legend="ul",
        legend_kw={"order": "F", "title": "column major"},
    )
    hs1.extend(h1)
    h2 = axs[1].plot(
        data,
        lw=4,
        cycle="Set3",
        label=label,
        legend="r",
        legend_kw={"lw": 8, "ncols": 1, "frame": False, "title": "modified\n handles"},
    )
    hs2.extend(h2)

# Outer legends
ax = axs[0]
ax.legend(hs1, loc="b", ncols=3, title="row major", order="C", facecolor="gray2")
ax = axs[1]
ax.legend(hs2, loc="b", ncols=3, center=True, title="centered rows")
axs.format(xlabel="xlabel", ylabel="ylabel", suptitle="Legend formatting demo")
