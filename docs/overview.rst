==============
Usage overview
==============

ProPlot is an *object-oriented* matplotlib wrapper, which means
most of its features derive from subclasses of the `~matplotlib.figure.Figure` and `~matplotlib.axes.Axes` classes. If you tend to use the `~matplotlib.pyplot` API and are not familiar with figure and axes "objects", you should first take a look at `this page <https://matplotlib.org/api/api_overview.html#the-pyplot-api>`__. Using the objects directly tends to be more clear and concise than `~matplotlib.pyplot` API, and makes life easier when working with multiple figures or axes.

..
   This page gives a condensed overview of these features, along with features
   outside of these classes.
..
   This page is meant as the starting point for new users. It is
   populated with links to the :ref:`API reference` and User Guide.
   For more in-depth descriptions, see :ref:`Why ProPlot?`.

Importing proplot
=================

We recommend importing ProPlot with

.. code-block:: python

   import proplot as plot

This differentiates ProPlot from the usual ``plt`` abbreviation used for the `~matplotlib.pyplot` module. Importing proplot immediately adds a bunch of new colormaps, property cyclers, color names, and fonts to matplotlib. See :ref:`Colormaps`, :ref:`Color cycles`, and :ref:`Colors and fonts` for details.

Figure and axes classes
=======================
ProPlot's features derive from `~proplot.subplots.subplots`, modeled after
the native `matplotlib.pyplot.subplots` command.
`~proplot.subplots.subplots` creates a `~proplot.subplots.Figure` subclass
populated with special `~proplot.axes.Axes` subclasses.
See :ref:`Creating figures`
and :ref:`Figure tight layout` for details.

Each `~proplot.axes.Axes` class also belongs to
one of the `~proplot.axes.XYAxes`, `~proplot.axes.PolarAxes`,
or `~proplot.axes.ProjAxes` parent classes, depending on the projection used. See
:ref:`Geographic and polar plots` for details.

The `~proplot.subplots.Figure` and `~proplot.axes.Axes` subclasses
include useful new methods and override several existing methods:

* The most important new method is `~proplot.axes.Axes.format`, whose behavior depends on whether the axes is an `~proplot.axes.XYAxes`, `~proplot.axes.PolarAxes`, or `~proplot.axes.ProjAxes`. This method fine-tunes various axes settings. See :ref:`Customizing figures` for details.
* The `~proplot.subplots.Figure` `~proplot.subplots.Figure.colorbar` and `~proplot.subplots.Figure.legend` and `~proplot.axes.Axes` `~proplot.axes.Axes.colorbar` and `~proplot.axes.Axes.legend` commands are used to add colorbars and legends *inside* of subplots, along the *outside edge* of subplots, and along the *edge of the figure*. See :ref:`Colorbars and legends` for details.
* There is a huge variety of new features for working with contour plots, pcolor plots, heatmaps, line plots, error bars, bar plots, area plots, and parametric plots. See :ref:`1d plotting` and :ref:`2d plotting` for details.

Integration with other packages
===============================
ProPlot includes integration with `xarray`, `pandas`, `cartopy`, and `~mpl_toolkits.basemap`:

* Axis labels, tick labels, titles, colorbar labels, and legend labels are automatically applied when you pass an `xarray.DataArray`, `pandas.DataFrame`, or `pandas.Series` object to any plotting command. This works just like the native `xarray.DataArray.plot` and `pandas.DataFrame.plot` methods. See :ref:`1d plotting` and :ref:`2d plotting` for details.
* The `~proplot.projs.Proj` function lets you make arbitrary grids of basemap `~mpl_toolkits.basemap.Basemap` and cartopy `~cartopy.crs.Projection` projections. It is used to interpret the `proj` keyword arg passed to `~proplot.subplots.subplots`. The resulting axes are instances of `~proplot.axes.ProjAxes` with `~proplot.axes.ProjAxes.format` methods that can be used to add geographic features and custom meridian and parallel gridlines. See :ref:`Geographic and polar plots` for details.

Other functions and classes
===========================
ProPlot includes a bunch of useful tools outside
of the `~proplot.subplots.Figure` and `~proplot.axes.Axes` subclasses:

* The `~proplot.styletools.Colormap` and `~proplot.styletools.Cycle` constructor functions. These can slice, merge, and modify colormaps and color cycles. See :ref:`Colormaps`, :ref:`Color cycles`, and :ref:`Colors and fonts` for details.
* The `~proplot.styletools.LinearSegmentedColormap` and  `~proplot.styletools.ListedColormap` subclasses, used to wrap the default matplotlib colormaps, and the new `~proplot.styletools.PerceptuallyUniformColormap` class, used for creating arbitrary colormaps with perceptually uniform transitions. See :ref:`Colormaps` for details.
* The `~proplot.styletools.Norm` constructor function, used to generated colormap normalizers from shorthand names; the `~proplot.styletools.LinearSegmentedNorm` normalizer, used to scale colors evenly w.r.t. index for arbitrarily spaced monotonic levels; and the `~proplot.styletools.BinNorm` meta-normalizer, used to discretized colormap colors. See :ref:`2d plotting` for details.
* The `~proplot.axistools.Locator`, `~proplot.axistools.Formatter`, and `~proplot.axistools.Scale` constructor functions, used to generate class instances from variable input types. These are used to interpret keyword arguments passed to `~proplot.axes.Axes.format` and `~proplot.subplots.Figure.colorbar`. See :ref:`X and Y axis settings` for details.
* The `~proplot.rctools.rc` object, an instance of `~proplot.rctools.rc_configurator`, for modifying global settings. You can also control settings with a ``~/.proplotrc`` file. See :ref:`Configuring proplot` for details.

