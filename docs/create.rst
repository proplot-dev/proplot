Figures and subplots
====================

The subplots command
--------------------
Figures and axes with special ProPlot features can be created with `~proplot.subplots.subplots` (or the alias `~proplot.figure.figure`). This function was inspired by `matplotlib.pyplot.subplots` and is packed with new features. It returns a figure instance and a special container "`~proplot.figure.axes_list`" for the axes.

To generate complex grids, pass a 2D array of numbers corresponding to unique subplots. Use zero to allot empty space. For example, the following command creates a grid with one tall plot on the left, and two smaller plots on the right:

.. code-block:: python

   import proplot as plot
   f, axs = plot.subplots([[1,2],[1,3]])

while this command creates a grid with one long plot on top and two smaller plots on the bottom, with a gap in the middle (the zero):

.. code-block:: python

   import proplot as plot
   f, axs = plot.subplots([[1,1,1],[2,0,4]])

The corresponding axes "numbers" are added as attributes, and can be used for a-b-c labeling with the `~proplot.BaseAxes.format` method.

For a single-axes figure (1 row, 1 column), simply call `~proplot.subplots` with no arguments. For simple grids, use the `nrows` and `ncols` keyword args; for example

.. code-block:: python

   import proplot as plot
   f, axs = plot.subplots(nrows=2)

The `~proplot.axes_list` container lets you invoke methods in every axes in ``axs`` simultaneously; for example:

.. code-block:: python

   import proplot as plot
   f, axs = plot.subplots()
   axs.format(xlim=(0,1)) # set to be identical everywhere
   axs[:3].format(color='red')  # make axis labels, spines, ticks red
   axs[3:].format(color='blue') # same, but make them blue

There are further arguments for drawing axes or figure "panels", setting particular axes
as map projections, sharing axis limits and tick/axis labels, scaling the
relative widths of rows/columns, setting axes aspect ratios and variable inter-axes
spacing, and setting average axes width/height instead of the full figure dimensions.

See `~proplot.subplots` for a full description of the arguments.

Shared axes
-----------
I've expanded the matplotlib `shared axis <https://matplotlib.org/examples/pylab_examples/shared_axis_demo.html>`_ setting and added a new "spanning axis label" feature. These
are best demonstrated by example; see :ref:`Examples` and `~proplot.subplots.subplots`.

Smarter subplot layout
----------------------

If you specify just one of `width`, `height`, `axwidth` (fixes
figure width), or `axheight` (fixes figure height), the unspecified
dimension will be **scaled** such that the top-left subplot has aspect
ratio `aspect`. **Spacing and panel widths are held fixed** during
this scaling.

This is accomplished with the new `~proplot.gridspec.FlexibleGridSpec` class, subclassed
from matplotlib’s `~matplotlib.gridspec.GridSpec` class. The “actual” `wspace` and
`hspace` passed to `~matplotlib.gridspec.GridSpec` are zero – the spaces you see in your
figure are empty subplot slots *masquerading* as spaces, whose widths
are controlled by `width_ratios` and `height_ratios`. Check out
`FlexibleGridSpec.__getitem__`.

Note this also means **inter-subplot spacing is now variable**. You can
specify `wspace` and `hspace` with:

1. A scalar constant, e.g.  ``wspace=0.2``.
2. A list of different spacings, e.g.  ``wspace=[0.1,0.5]`` to offset the 3rd column in a 3-column plot from the rest.

Smarter “tight” layout
----------------------

The subplot layout changes allowed me to create the
`~proplot.figure.FigureBase.smart_tight_layout` method. By default, this method is called
whenever the figure is drawn (i.e. when it is rendered by the
matplotlib backend or saved to file) – disable this behavior with
``plot.subplots(tight=False)``.

Previously, `~matplotlib.figure.Figure.tight_layout` could be used to fit the figure borders
over a box that perfectly encompasses all artists (i.e. text, subplots,
etc.). However, because `~matplotlib.gridspec.GridSpec` spaces are relative to the subplot
dimensions, changing the figure dimensions also changes the
inter-subplot spacings. Since your font size is specified in points
(i.e. a physical unit), *this can easily cause text to overlap with
other subplots where it didn’t before*.

The new `~proplot.figure.FigureBase.smart_tight_layout` method draws a tight bounding box that
**preserves inter-subplot spacing, panel widths, and subplot aspect
ratios**. It does so by letting either the height or width dimension of the figure vary;
by default, the height is allowed to vary. If you instead specify a fixed figure
size, the aspect ratios of subplots will vary -- but inter-subplot spacing and panel widths
will still be preserved.

