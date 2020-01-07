=============
Using ProPlot
=============

..
   This page gives a condensed overview of these features, along with features
   outside of these classes.
..
   This page is meant as the starting point for new users. It is
   populated with links to the :ref:`API reference` and User Guide.
   For more in-depth descriptions, see :ref:`Why ProPlot?`.

Background
==========

ProPlot is an object-oriented matplotlib wrapper. The "wrapper" part means that
ProPlot's features are largely a *superset* of matplotlib.
You can use your favorite plotting commands like
`~matplotlib.axes.Axes.plot`, `~matplotlib.axes.Axes.scatter`, `~matplotlib.axes.Axes.contour`, and `~matplotlib.axes.Axes.pcolor` like you always have.
The "object-oriented" part means that ProPlot's features are implemented with *subclasses* of the `~matplotlib.figure.Figure` and `~matplotlib.axes.Axes` classes.

If you tend to use `~matplotlib.pyplot` and are not familiar with figure and axes *classes*, check out `this guide from the matplotlib documentation <https://matplotlib.org/api/api_overview.html#the-pyplot-api>`__. Working with objects directly tends to be more clear and concise than `~matplotlib.pyplot`, makes things easier when working with multiple figures and axes, and is certainly more "`pythonic <https://www.python.org/dev/peps/pep-0020/>`__". Therefore, although some ProPlot features may still work, we do not officially support the `~matplotlib.pyplot` API.


Importing proplot
=================

We recommend importing ProPlot as follows:

.. code-block:: python

   import proplot as plot

This differentiates ProPlot from the usual ``plt`` abbreviation used for the `~matplotlib.pyplot` module.
Importing ProPlot immediately adds several new colormaps, property cyclers, color names, and fonts to matplotlib. See :ref:`Colormaps`, :ref:`Color cycles`, and :ref:`Colors and fonts` for details.

Importing ProPlot also configures your IPython environment by setting up the matplotlib backend and enabling the autoreload and autosave extensions (this feature can be disabled). See `~proplot.rctools.ipython_matplotlib`, `~proplot.rctools.ipython_autoreload`, and `~proplot.rctools.ipython_autosave` for details.

Figure and axes classes
=======================
Making figures in ProPlot always begins with a call to the
`~proplot.subplots.subplots` command:

.. code-block:: python

   f, axs = plot.subplots(...)

`~proplot.subplots.subplots` is modeled after
matplotlib's native `matplotlib.pyplot.subplots` command.
It creates an instance of ProPlot's
`~proplot.subplots.Figure` class
populated with instances of ProPlot's
`~proplot.axes.Axes` classes.
See :ref:`The basics`
and :ref:`Subplots features` for details.

Each `~proplot.axes.Axes` returned by `~proplot.subplots.subplots`
belongs to one of the following three child classes:

* `~proplot.axes.XYAxes`: For plotting simple data with *x* and *y* coordinates.
* `~proplot.axes.ProjAxes`: For geographic plots with *longitude* and *latitude* coordinates.
* `~proplot.axes.PolarAxes`: For polar plots with *radius* and *azimuth* coordinates.

See :ref:`X and Y axis settings` for details on working with `~proplot.axes.XYAxes` and
:ref:`Geographic and polar plots` for details on working with
`~proplot.axes.ProjAxes` and `~proplot.axes.PolarAxes`.

Figure and axes methods
=======================
The `~proplot.subplots.Figure` and `~proplot.axes.Axes` subclasses
include several *brand new* methods and add to the functionality of several *existing* methods.

* The new `~proplot.axes.Axes.format` method is used to fine-tune various axes settings.  Its behavior depends on whether the axes is an `~proplot.axes.XYAxes`, `~proplot.axes.PolarAxes`, or `~proplot.axes.ProjAxes`. Think of this as a dedicated `~matplotlib.artist.Artist.update` method for axes artists. See :ref:`Formatting subplots` and :ref:`Changing rc settings` for details.
* The `~proplot.subplots.Figure` `~proplot.subplots.Figure.colorbar` and `~proplot.subplots.Figure.legend` and `~proplot.axes.Axes` `~proplot.axes.Axes.colorbar` and `~proplot.axes.Axes.legend` commands are used to add colorbars and legends *inside* of subplots, along the *outside edge* of subplots, and along the *edge of the figure*. They considerably simplify the process of drawing colorbars and legends. See :ref:`Colorbars and legends` for details.
* ProPlot adds a huge variety of features for working with `~matplotlib.axes.Axes.contour` plots, `~matplotlib.axes.Axes.pcolor` plots, `~matplotlib.axes.Axes.plot` lines, `~proplot.axes.Axes.heatmap` plots, `~matplotlib.axes.Axes.errorbar` bars, `~matplotlib.axes.Axes.bar` plots, `~proplot.axes.Axes.area` plots, and `~proplot.axes.Axes.parametric` plots. See :ref:`1d plotting` and :ref:`2d plotting` for details.

Integration with other packages
===============================
ProPlot's features are integrated with the data containers
introduced by `xarray` and `pandas` and the
`cartopy` and `~mpl_toolkits.basemap` geographic
plotting toolkits.

* Axis labels, tick labels, titles, colorbar labels, and legend labels are automatically applied when you pass an `xarray.DataArray`, `pandas.DataFrame`, or `pandas.Series` object to any plotting command. This works just like the native `xarray.DataArray.plot` and `pandas.DataFrame.plot` methods. See :ref:`1d plotting` and :ref:`2d plotting` for details.
* The `~proplot.projs.Proj` function lets you make arbitrary grids of basemap `~mpl_toolkits.basemap.Basemap` and cartopy `~cartopy.crs.Projection` projections. It is used to interpret the `proj` keyword arg passed to `~proplot.subplots.subplots`. The resulting axes are instances of `~proplot.axes.ProjAxes` with `~proplot.axes.ProjAxes.format` methods that can be used to add geographic features and custom meridian and parallel gridlines. See :ref:`Geographic and polar plots` for details.

New functions and classes
=========================
ProPlot includes several useful *constructor functions*
and *subclasses* outside
of the `~proplot.subplots.Figure` and `~proplot.axes.Axes` subclasses.

* The `~proplot.styletools.Colormap` and `~proplot.styletools.Cycle` constructor functions can slice, merge, and modify colormaps and color cycles. See :ref:`Colormaps`, :ref:`Color cycles`, and :ref:`Colors and fonts` for details.
* The `~proplot.styletools.LinearSegmentedColormap` and  `~proplot.styletools.ListedColormap` subclasses replace the default matplotlib colormap classes and add several methods. The new `~proplot.styletools.PerceptuallyUniformColormap` class is used to make colormaps with perceptually uniform transitions. See :ref:`Colormaps` for details.
* The `~proplot.styletools.show_cmaps`, `~proplot.styletools.show_cycles`, `~proplot.styletools.show_colors`, `~proplot.styletools.show_fonts`, `~proplot.styletools.show_channels`, and `~~proplot.styletools.show_colorspaces` functions are used to visualize your color scheme and font options and inspect individual colormaps.
* The `~proplot.styletools.Norm` constructor function generates colormap normalizers from shorthand names. The new `~proplot.styletools.LinearSegmentedNorm` normalizer scales colors evenly w.r.t. index for arbitrarily spaced monotonic levels, and the new `~proplot.styletools.BinNorm` meta-normalizer is used to discretized colormap colors. See :ref:`2d plotting` for details.
* The `~proplot.axistools.Locator`, `~proplot.axistools.Formatter`, and `~proplot.axistools.Scale` constructor functions, used to generate class instances from variable input types. These are used to interpret keyword arguments passed to `~proplot.axes.Axes.format` and `~proplot.subplots.Figure.colorbar`. See :ref:`X and Y axis settings` for details.
* The `~proplot.rctools.rc` object, an instance of `~proplot.rctools.rc_configurator`, is used for modifying *individual* global settings, changing settings in *bulk*, and temporarily changing settings in *context blocks*. You can also control settings with a ``~/.proplotrc`` file. See :ref:`Configuring proplot` for details.
