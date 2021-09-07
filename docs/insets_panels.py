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
# .. _ug_insets_panels:
#
# Insets and panels
# =================

# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_panels:
#
# Panel axes
# ----------
#
# It is often useful to have narrow "panels" along the edge of a larger
# subplot for plotting secondary 1-dimensional datasets or summary statistics.
# In ProPlot, you can generate panels using the `~proplot.axes.Axes.panel_axes`
# command (or its shorthand, `~proplot.axes.Axes.panel`). The panel location
# is specified with a string, e.g. ``ax.panel('r')`` or ``ax.panel('right')``
# for a right-hand side panel, and the resulting panels are instances of
# `~proplot.axes.CartesianAxes`. By default, the panel shares its axis limits,
# axis labels, tick positions, and tick labels with the main subplot, but
# this can be disabled by passing ``share=False``. To generate "stacked" panels,
# call `~proplot.axes.Axes.panel_axes` more than once. To generate several
# panels at once, call `~proplot.gridspec.SubplotGrid.panel_axes` on
# the `~proplot.gridspec.SubplotGrid` returned by `~proplot.figure.Figure.subplots`.
# Note that panels :ref:`do not interfere with the tight layout algorithm <ug_tight>`
# and :ref:`do not affect the subplot aspect ratios <ug_autosize>`.
#
# In the first example below, the distances are automatically adjusted by the
# :ref:`tight layout algorithm <ug_tight>` according to the `pad` keyword
# (the default is :rcraw:`subplots.panelpad` -- this can be changed for an entire
# figure by passing `panelpad` to `~proplot.figure.Figure`). In the second example,
# the tight layout algorithm is overriden by manually setting the `space` to ``0``.
# Panel widths are specified in physical units, with the default controlled
# by :rcraw:`subplots.panelwidth`. This helps preserve the look of the
# figure if the figure size changes. Note that by default, panels are excluded
# when centering :ref:`spanning axis labels <ug_share>` and super titles.
# To include the panels, pass ``includepanels=True`` to `~proplot.figure.Figure`.

# %%
import proplot as pplt

# Demonstrate that complex arrangements preserve
# spacing, aspect ratios, and axis sharing
gs = pplt.GridSpec(nrows=2, ncols=2)
fig = pplt.figure(refwidth=1.5, share=False)
for ss, side in zip(gs, 'tlbr'):
    ax = fig.add_subplot(ss)
    px = ax.panel_axes(side, width='3em')
fig.format(
    xlim=(0, 1), ylim=(0, 1),
    xlabel='xlabel', ylabel='ylabel',
    xticks=0.2, yticks=0.2,
    title='Title', suptitle='Complex arrangement of panels',
    toplabels=('Column 1', 'Column 2'),
    abc=True, abcloc='ul', titleloc='uc', titleabove=False,
)

# %%
import proplot as pplt
import numpy as np
state = np.random.RandomState(51423)
data = (state.rand(20, 20) - 0.48).cumsum(axis=1).cumsum(axis=0)
data = 10 * (data - data.min()) / (data.max() - data.min())

# Stacked panels with outer colorbars
for cbarloc, ploc in ('rb', 'br'):
    # Create figure
    fig, axs = pplt.subplots(
        nrows=1, ncols=2, refwidth=1.8, panelpad=0.8,
        share=False, includepanels=True
    )
    axs.format(
        xlabel='xlabel', ylabel='ylabel', title='Title',
        suptitle='Using panels for summary statistics',
    )

    # Plot 2D dataset
    for ax in axs:
        ax.contourf(
            data, cmap='glacial', extend='both',
            colorbar=cbarloc, colorbar_kw={'label': 'colorbar'},
        )

    # Get summary statistics and settings
    axis = int(ploc == 'r')  # dimension along which stats are taken
    x1 = x2 = np.arange(20)
    y1 = data.mean(axis=axis)
    y2 = data.std(axis=axis)
    titleloc = 'upper center'
    if ploc == 'r':
        titleloc = 'center'
        x1, x2, y1, y2 = y1, y2, x1, x2

    # Panels for plotting the mean. Note SubplotGrid.panel() returns a SubplotGrid
    # of panel axes. We use this to call format() for all the panels at once.
    space = 0
    width = '4em'
    kwargs = {'titleloc': titleloc, 'xreverse': False, 'yreverse': False}
    pxs = axs.panel(ploc, space=space, width=width)
    pxs.format(title='Mean', **kwargs)
    for px in pxs:
        px.plot(x1, y1, color='gray7')

    # Panels for plotting the standard deviation
    pxs = axs.panel(ploc, space=space, width=width)
    pxs.format(title='Stdev', **kwargs)
    for px in pxs:
        px.plot(x2, y2, color='gray7', ls='--')


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_insets:
#
# Inset axes
# ----------
#
# `Inset axes
# <https://matplotlib.org/stable/gallery/subplots_axes_and_figures/zoom_inset_axes.html>`__
# can be generated with the `~proplot.axes.Axes.inset_axes` command (or its
# shorthand, `~proplot.axes.Axes.inset`). To generate several insets at once, call
# `~proplot.gridspec.SubplotGrid.inset_axes` on the `~proplot.gridspec.SubplotGrid`
# returned by `~proplot.figure.Figure.subplots`. By default, inset axes have the
# same projection as the parent axes, but you can also request a :ref:`different
# projection <ug_proj>` (e.g., ``ax.inset_axes(bounds, proj='polar')``). When the
# axes are both `~proplot.axes.Axes.CartesianAxes`, you can pass ``zoom=True``
# to `~proplot.axes.Axes.inset_axes` to quickly add a "zoom indication" box and
# lines (this uses `~matplotlib.axes.Axes.indicate_inset_zoom` internally). The box
# and line positions automatically follow the axis limits of the inset axes and parent
# axes. To modify the zoom line properties, you can pass a dictionary to `zoom_kw`.

# %%
import proplot as pplt
import numpy as np

# Sample data
N = 20
state = np.random.RandomState(51423)
x, y = np.arange(10), np.arange(10)
data = state.rand(10, 10).cumsum(axis=0)

# Plot data
fig, ax = pplt.subplots(refwidth=3)
m = ax.pcolormesh(data, cmap='Grays', levels=N)
ax.colorbar(m, loc='b', label='label')
ax.format(
    xlabel='xlabel', ylabel='ylabel',
    suptitle='"Zooming in" with an inset axes'
)

# Create inset axes representing a "zoom-in"
ix = ax.inset(
    [5, 5, 4, 4], transform='data', zoom=True,
    zoom_kw={'edgecolor': 'red3', 'lw': 2, 'ls': '--'}
)
ix.format(
    xlim=(2, 4), ylim=(2, 4), metacolor='red7',
    linewidth=1.5, ticklabelweight='bold'
)
ix.pcolormesh(data, cmap='Grays', levels=N, inbounds=False)
