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
# In matplotlib, there is no simple way to do this. In ProPlot, you can
# create panels by passing a location (e.g. ``loc='r'`` or ``loc='right'``)
# to the `~proplot.axes.Axes.panel` or `~proplot.axes.Axes.panel_axes` methods.
# The resulting axes are instances of `~proplot.axes.CartesianAxes`.
#
# To generate "stacked" panels, simply call `~proplot.axes.Axes.panel` more
# than once. To include panels when centering spanning axis labels and super
# titles, pass ``includepanels=True`` to `~proplot.figure.Figure`. Panels
# :ref:`do not interfere with the tight layout algorithm <ug_tight>` and
# :ref:`do not affect the subplot aspect ratios <ug_autosize>`.
#
# In the first example below, the panel distance from the main subplot is
# manually set to ``space=0``. In the second example, it is adjusted automatically
# by the tight layout algorithm.

# %%
import proplot as plot
import numpy as np
state = np.random.RandomState(51423)
data = (state.rand(20, 20) - 0.48).cumsum(axis=1).cumsum(axis=0)
data = 10 * (data - data.min()) / (data.max() - data.min())

# Stacked panels with outer colorbars
for cbarlocation, plocation in ('rb', 'br'):
    fig, axs = plot.subplots(
        axwidth=1.6, nrows=1, ncols=2,
        share=0, panelpad=0.1, includepanels=True
    )
    axs.contourf(
        data, cmap='glacial', extend='both',
        colorbar=cbarlocation, colorbar_kw={'label': 'colorbar'},
    )

    # Summary statistics and settings
    titleloc = 'upper center'
    axis = int(plocation == 'r')  # dimension along which stats are taken
    x1 = x2 = np.arange(20)
    y1 = data.mean(axis=axis)
    y2 = data.std(axis=axis)
    if plocation == 'r':
        titleloc = 'center'
        x1, x2, y1, y2 = y1, y2, x1, x2
    kwargs = {'titleloc': titleloc, 'xreverse': False, 'yreverse': False}
    space = 0
    width = '30pt'

    # Panels for plotting the mean
    panels = axs.panel(plocation, space=space, width=width)
    panels.plot(x1, y1, color='gray7')
    panels.format(title='Mean', **kwargs)

    # Panels for plotting the standard deviation
    panels = axs.panel(plocation, space=space, width=width)
    panels.plot(x2, y2, color='gray7', ls='--')
    panels.format(title='Stdev', **kwargs)

    # Apply formatting *after*
    axs.format(
        xlabel='xlabel', ylabel='ylabel', title='Title',
        suptitle='Using panels for summary statistics',
    )

# %%
import proplot as plot
fig, axs = plot.subplots(axwidth=1.5, nrows=2, ncols=2, share=0)

# Demonstration that complex arrangements of panels
# do not mess up tight layout algorithm
for ax, side in zip(axs, 'tlbr'):
    ax.panel(side, width='3em')
axs.format(
    xlim=(0, 1), ylim=(0, 1),
    xlabel='xlabel', ylabel='ylabel',
    yticks=plot.arange(0.2, 0.8, 0.2),
    xticks=plot.arange(0.2, 0.8, 0.2),
    title='Title', suptitle='Complex arrangement of panels',
    collabels=['Column 1', 'Column 2'],
    abc=True, abcloc='ul', titleloc='uc', abovetop=False,
)


# %% [raw] raw_mimetype="text/restructuredtext"
# .. _ug_insets:
#
# Inset axes
# ----------
#
# `Inset axes\
# <https://matplotlib.org/3.1.1/gallery/subplots_axes_and_figures/zoom_inset_axes.html>`__
# can be generated with the `~proplot.axes.Axes.inset` or
# `~proplot.axes.Axes.inset_axes` command. By defaut, the resulting axes
# use the same projection as the parent axes, but you can also specify
# a different projection (for example, ``ax.inset(bounds, proj='polar')``).
# Passing ``zoom=True`` to `~proplot.axes.Axes.inset` draws "zoom indication"
# lines with `~matplotlib.axes.Axes.indicate_inset_zoom` when the axes are both
# Cartesian, and ProPlot automatically updates the lines when the axis limits of
# the parent axes change. To modify the line properties, simply use the `zoom_kw`
# argument.

# %%
import proplot as plot
import numpy as np

# Generate sample data
N = 20
state = np.random.RandomState(51423)
x, y = np.arange(10), np.arange(10)
data = state.rand(10, 10)

# Plot sample data
fig, ax = plot.subplots(axwidth=3)
m = ax.pcolormesh(data, cmap='Grays', levels=N)
ax.colorbar(m, loc='b', label='label')
ax.format(
    xlabel='xlabel', ylabel='ylabel',
    suptitle='"Zooming in" with an inset axes'
)

# Create inset axes representing a "zoom-in"
iax = ax.inset(
    [5, 5, 4, 4], transform='data', zoom=True,
    zoom_kw={'color': 'red3', 'lw': 2, 'ls': '--'}
)
iax.format(
    xlim=(2, 4), ylim=(2, 4), color='red7',
    linewidth=1.5, ticklabelweight='bold'
)
iax.pcolormesh(data, cmap='Grays', levels=N)
