.. _cartopy: https://scitools.org.uk/cartopy/docs/latest/

.. _basemap: https://matplotlib.org/basemap/index.html

.. _seaborn: https://seaborn.pydata.org

.. _pandas: https://pandas.pydata.org

.. _xarray: http://xarray.pydata.org/en/stable/

.. _usage:

=============
Using ProPlot
=============

This page offers a condensed overview of ProPlot's features. It is populated
with links to the :ref:`API reference` and :ref:`User Guide <ug_basics>`.
For a more in-depth discussion, see :ref:`Why ProPlot?`

.. _usage_background:

Background
==========

ProPlot is an object-oriented matplotlib wrapper. The "wrapper" part means
that ProPlot's features are largely a *superset* of matplotlib.  You can use
plotting commands like `~matplotlib.axes.Axes.plot`, `~matplotlib.axes.Axes.scatter`,
`~matplotlib.axes.Axes.contour`, and `~matplotlib.axes.Axes.pcolor` like you always
have. The "object-oriented" part means that ProPlot's features are implemented with
*subclasses* of the `~matplotlib.figure.Figure` and `~matplotlib.axes.Axes` classes.

If you tend to use `~matplotlib.pyplot` and are not familiar with the figure and axes
classes, check out `this guide <https://matplotlib.org/stable/api/index.html>`__.
Directly working with matplotlib classes tends to be more clear and concise than
`~matplotlib.pyplot`, makes things easier when working with multiple figures and axes,
and is certainly more "`pythonic <https://www.python.org/dev/peps/pep-0020/>`__".
Therefore, although many ProPlot features may still work, we do not officially
support the `~matplotlib.pyplot` interface.

.. _usage_import:

Importing proplot
=================

Importing ProPlot immediately adds several
new :ref:`colormaps <ug_cmaps>`, :ref:`property cycles <ug_cycles>`,
:ref:`color names <ug_colors>`, and :ref:`fonts <ug_fonts>` to matplotlib.
If you are only interested in these features, you may want to
import ProPlot at the top of your script and do nothing else!
We recommend importing ProPlot as follows:

.. code-block:: python

   import proplot as pplt

This differentiates ProPlot from the usual ``plt`` abbreviation reserved for
the `~matplotlib.pyplot` module.

.. _usage_classes:

Figure and axes classes
=======================

Creating figures with ProPlot is very similar to
matplotlib. You can either create the figure and
all of its subplots at once:

.. code-block:: python

   fig, axs = pplt.subplots(...)

or create an empty figure
then fill it with subplots:

.. code-block:: python

   fig = pplt.figure(...)
   axs = fig.add_subplots(...)  # add several subplots
   ax = fig.add_subplot(...)  # add a single subplot
   # axs = fig.subplots(...)  # shorthand
   # ax = fig.subplot(...)  # shorthand

These commands are modeled after `matplotlib.pyplot.subplots` and
`matplotlib.pyplot.figure` and are :ref:`packed with new features <ug_layout>`.
One highlight is the `~proplot.figure.Figure.auto_layout` algorithm that
:ref:`automatically adjusts the space between subplots <ug_tight>` (similar to
matplotlib's `tight layout
<https://matplotlib.org/stable/tutorials/intermediate/tight_layout_guide.html>`__)
and :ref:`automatically adjusts the figure size <ug_autosize>` to preserve subplot
sizes and aspect ratios (particularly useful for grids of map projections
and images). All sizing arguments take :ref:`arbitrary units <ug_units>`,
including metric units like ``cm`` and ``mm``.

Instead of the native `matplotlib.figure.Figure` and `matplotlib.axes.Axes` classes,
ProPlot uses the `proplot.figure.Figure`, `proplot.axes.Axes`, and
`proplot.axes.PlotAxes` subclasses. ProPlot figures are saved with
`~proplot.figure.Figure.save` or `~matplotlib.figure.Figure.savefig`,
and ProPlot axes belong to one of the following three child classes:

* `proplot.axes.CartesianAxes`:
  For ordinary plots with *x* and *y* coordinates.
* `proplot.axes.GeoAxes`:
  For geographic plots with *longitude* and *latitude* coordinates.
* `proplot.axes.PolarAxes`:
  For polar plots with *azimuth* and *radius* coordinates.

Most of ProPlot's features are implemented using these subclasses.
They include several new figure and axes methods and added
functionality to existing figure and axes methods.

* The `proplot.axes.Axes.format` and `proplot.figure.Figure.format` commands fine-tunes
  various axes and figure settings.  Think of this as a dedicated
  `~matplotlib.artist.Artist.update` method for axes and figures. See
  :ref:`formatting subplots <ug_format>` for a broad overview, along with the
  individual sections on formatting :ref:`Cartesian plots <ug_cartesian>`,
  :ref:`geographic plots <ug_geoformat>`, and :ref:`polar plots <ug_polar>`.
* The `proplot.axes.Axes.colorbar` and `proplot.axes.Axes.legend` commands
  draw colorbars and legends inside of subplots or along the outside edges of
  subplots. The `proplot.figure.Figure.colorbar` and `proplot.figure.Figure.legend`
  commands draw colorbars or legends along the edges of figures (aligned by subplot
  boundaries). These commands considerably :ref:`simplify <ug_cbars_legends>`
  the process of drawing colorbars and legends.
* The `proplot.axes.PlotAxes` subclass (used for all ProPlot axes)
  adds many, many useful features to virtually every plotting command
  (including `~proplot.axes.PlotAxes.plot`, `~proplot.axes.PlotAxes.scatter`,
  `~proplot.axes.PlotAxes.bar`, `~proplot.axes.PlotAxes.area`,
  `~proplot.axes.PlotAxes.contour`, and `~proplot.axes.PlotAxes.pcolor`.
  See the :ref:`1D plotting <ug_1dplots>` and :ref:`2D plotting <ug_2dplots>`
  sections for details.

.. _usage_integration:

Integration features
====================

ProPlot includes *optional* integration features with four external
packages: the `pandas`_ and `xarray`_ packages, used for working with annotated
tables and arrays, and the `cartopy`_ and `basemap`_ geographic
plotting packages.

* If you pass a `pandas.Series`, `pandas.DataFrame`, or `xarray.DataArray`
  to any plotting command, the axis labels, tick labels, titles, colorbar
  labels, and legend labels are automatically applied from the metadata. If
  you did not supply the *x* and *y* coordinates, they are also inferred from
  the metadata. This works just like the native `xarray.DataArray.plot` and
  `pandas.DataFrame.plot` methods. A demonstration of this feature is given
  in the sections on :ref:`1D plotting <ug_1dintegration>` and
  :ref:`2D plotting <ug_2dintegration>`. This feature can be disabled by
  setting :rcraw:`autoformat` to ``False``.
* The `~proplot.axes.GeoAxes` class uses the `cartopy`_ or
  `basemap`_ packages to :ref:`plot geophysical data <ug_geoplot>`,
  :ref:`add geographic features <ug_geoformat>`, and
  :ref:`format projections <ug_geoformat>`. `~proplot.axes.GeoAxes` provides
  provides a simpler, cleaner interface than the original `cartopy`_ and `basemap`_
  interfaces. Figures can be filled with `~proplot.axes.GeoAxes` by passing the
  `proj` keyword to `~proplot.ui.subplots`.

Since these features are optional, ProPlot can be used without installing
any of these packages.

.. _usage_features:

Additional features
===================

Outside of the features provided by the `proplot.figure.Figure` and
`proplot.axes.Axes` subclasses, ProPlot includes several useful
classes and :ref:`constructor functions <why_constructor>`.

* The `~proplot.constructor.Colormap` and `~proplot.constructor.Cycle`
  constructor functions can be used to :ref:`slice <ug_cmaps_mod>`,
  and :ref:`merge <ug_cmaps_merge>` existing colormaps and color
  cycles. It can also :ref:`make new colormaps <ug_cmaps_new>`
  and :ref:`color cycles <ug_cycles_new>` from scratch.
* The `~proplot.colors.ContinuousColormap` and
  `~proplot.colors.DiscreteColormap` subclasses replace the default matplotlib
  colormap classes and add several methods. The new
  `~proplot.colors.PerceptualColormap` class is used to make
  colormaps with :ref:`perceptually uniform transitions <ug_perceptual>`.
* The `~proplot.demos.show_cmaps`, `~proplot.demos.show_cycles`,
  `~proplot.demos.show_colors`, `~proplot.demos.show_fonts`,
  `~proplot.demos.show_channels`, and `~proplot.demos.show_colorspaces`
  functions are used to visualize your :ref:`color scheme <ug_colors>`
  and :ref:`font options <ug_fonts>` and
  :ref:`inspect individual colormaps <ug_perceptual>`.
* The `~proplot.constructor.Norm` constructor function generates colormap
  normalizers from shorthand names. The new
  `~proplot.colors.SegmentedNorm` normalizer scales colors evenly
  w.r.t. index for arbitrarily spaced monotonic levels, and the new
  `~proplot.colors.DiscreteNorm` meta-normalizer is used to
  :ref:`break up colormap colors into discrete levels <ug_discrete>`.
* The `~proplot.constructor.Locator`, `~proplot.constructor.Formatter`, and
  `~proplot.constructor.Scale` constructor functions return corresponding class
  instances from flexible input types. These are used to interpret keyword
  arguments passed to `~proplot.axes.Axes.format`, and can be used to quickly
  and easily modify :ref:`x and y axis settings <ug_cartesian>`.
* The `~proplot.config.rc` object, an instance of
  `~proplot.config.Configurator`, is used for
  :ref:`modifying individual settings, changing settings in bulk, and
  temporarily changing settings in context blocks <ug_rc>`.
  It also introduces several :ref:`new setings <ug_config>`
  and sets up the inline plotting backend with `~proplot.config.inline_backend_fmt`
  so that your inline figures look the same as your saved figures.
