Figure creation
===============

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

The `~proplot.axes_list` container lets you invoke methods in every axes in `axs` simultaneously; for example:

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
This feature is demonstrated in the :ref:`examples`.

Smarter subplot layout
----------------------

If you specify just one of ``width``, ``height``, ``axwidth`` (fixes
figure width), or ``axheight`` (fixes figure height), the unspecified
dimension will be **scaled** such that the top-left subplot has aspect
ratio ``aspect``. **Spacing and panel widths are held fixed** during
this scaling.

This is accomplished with the new ``FlexibleGridSpec`` class, subclassed
from matplotlib’s ``GridSpec`` class. The “actual” ``wspace`` and
``hspace`` passed to ``GridSpec`` are zero – the spaces you see in your
figure are empty subplot slots *masquerading* as spaces, whose widths
are controlled by ``width_ratios`` and ``height_ratios``. Check out
``FlexibleGridSpec.__getitem__``.

Note this also means **inter-subplot spacing is now variable**. You can
specify ``wspace`` and ``hspace`` with: 1. A scalar constant, e.g.
``wspace=0.2``. 1. A list of different spacings, e.g.
``wspace=[0.1,0.5]`` to offset the 3rd column in a 3-column plot from
the rest.

Smarter “tight” layout
----------------------

The subplot layout changes allowed me to create the
``smart_tight_layout`` method. By default, this method is **called
whenever the figure is drawn** (i.e. when it is rendered by the
matplotlib backend or saved to file) – disable this behavior with
``plot.subplots(tight=False)``.

Previously, ``tight_layout`` could be used to fit the figure borders
over a box that perfectly encompasses all artists (i.e. text, subplots,
etc.). However, because ``GridSpec`` spaces are relative to the subplot
dimensions, changing the figure dimensions also changes the
inter-subplot spacings. Since your font size is specified in points
(i.e. a physical unit), *this can easily cause text to overlap with
other subplots where it didn’t before*.

The new ``smart_tight_layout`` method draws a tight bounding box that
**preserves inter-subplot spacing, panel widths, and subplot aspect
ratios**.

Academic journal standards
--------------------------

To create figures with dimensions that satisfy journal standards, use
the `journal` keyword argument.

Example:

.. code-block:: python

   f, axs = plot.subplots(ncols=3, nrows=2, journal='ams2') # medium-sized figure for AMS journal

The currently available specifiers are found in the `~proplot.gridspec.journal_size`
documentation.
