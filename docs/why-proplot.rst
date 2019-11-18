============
Why ProPlot?
============

ProPlot is meant to improve upon the parts of matplotlib that
tend to be laborious or cumbersome
when trying to create
beautiful graphics. This page
enumerates these limitations and
describes how ProPlot addresses them.

Efficient modifications
=======================
Problem
-------
Power users often need to change lots of plot settings all at once. In matplotlib, this requires a bunch of one-liner setters and getters, like `~matplotlib.axes.Axes.set_title`. 

This workflow is cumbersome, verbose, and often confusing. It can be unclear whether settings can be changed from an `~matplotlib.axes.Axes` setter, an `~matplotlib.axes.Axes.xaxis` or `~matplotlib.axes.Axes.yaxis` setter, or a miscellaneous bulk function like `~matplotlib.axes.Axes.tick_params`.

Solution
--------
In ProPlot, the `~proplot.axes.Axes.format` command is introduced for changing **arbitrary** settings in bulk. This massively reduces the amount of code needed to create highly customized figures. For example, it is trivial to see that

.. code-block:: python

   import proplot as plot
   f, ax = plot.subplots()
   ax.format(linewidth=1, color='gray')
   ax.format(xticks=20, xtickminor=True, xlabel='x axis', ylabel='y axis')

is much more succinct than

.. code-block:: python

   import matplotlib.pyplot as plt
   import matplotlib.ticker as mticker
   from matplotlib import rcParams
   rcParams['axes.linewidth'] = 1
   rcParams['axes.color'] = 'gray'
   fig, ax = plt.subplots()
   ax.xaxis.set_major_locator(mticker.MultipleLocator(10))
   ax.tick_params(width=1, color='gray', labelcolor='gray')
   ax.tick_params(axis='x', which='minor', bottom=True)
   ax.set_xlabel('x axis', color='gray')
   ax.set_ylabel('y axis', color='gray')

Constructor functions
=====================
Problem
-------
Matplotlib and cartopy introduce a bunch of classes with verbose names like `~matplotlib.ticker.MultipleLocator`, `~matplotlib.ticker.FormatStrFormatter`, `~cartopy.crs.AlbersEqualArea`, and `~cartopy.crs.NorthPolarStereo`.
But matplotlib `axis scales <https://matplotlib.org/3.1.0/gallery/scales/scales.html>`__ and `colormaps <https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html>`__ are **registered** under string names, as are `map projections <https://matplotlib.org/basemap/users/mapsetup.html>`__ in basemap. So why not register other verbose class names too?

Solution
--------
In ProPlot, tick locators, tick formatters, axis scales, cartopy projections, colormaps, and property cyclers are all **registered**. This is done by creating **constructor functions** and passing various keyword argument through the functions, allowing users to use class names instead of instances.

This may seem "unpythonic" but it is absolutely invaluable when making plots. Plotting code has a half life of about 30 seconds, and typing out all of these extra class names and import statements gets to be a major drag.

The below table lists the constructor functions and the keyword arguments that use them.

==============================  =============================  ========================================================
Function                        Returns                        Interpreted by
==============================  =============================  ========================================================
`~proplot.axistools.Locator`    Axis locator                   ``locator=``, ``xlocator=``, ``ylocator=``
`~proplot.axistools.Formatter`  Axis formatter                 ``formatter=``, ``xformatter=``, ``yformatter=``
`~proplot.axistools.Scale`      Axis scale                     ``xscale=``, ``yscale=``
`~proplot.styletools.Colormap`  Colormap                       ``cmap=``
`~proplot.styletools.Cycle`     Property cycler                ``cycle=``
`~proplot.projs.Proj`           Cartopy or basemap projection  ``proj=``
==============================  =============================  ========================================================

Fluid figure dimensions
=======================
Problem
-------
In matplotlib, you have to specify the physical dimensions of the figure. However, the dimensions of the *individual subplots* are often more important:

#. The subplot aspect ratio is usually more relevant than the figure aspect ratio, e.g. for projections.
#. The subplot dimensions control the evident thickness of text and other content plotted inside the box.

Solution
--------
In ProPlot, you can specify the physical dimensions of *subplots* instead of the figure by passing ``axwidth`` or ``axheight`` to `~proplot.subplots.Figure`. The default behavior is ``axwidth=2``. Figure dimensions are then automatically calculated to accommodate the subplot geometry and the spacing adjustments.

Several matplotlib backends require figure dimensions to be fixed. When `~matplotlib.subplots.Figure.draw` changes the figure dimensions, this can "surprise" the backend and cause unexpected behavior. ProPlot fixes this issue for the static inline backend and the Qt popup backend. However, this issue is unfixable the "notebook" inline backend, the "macosx" popup backend, and possibly other untested backends.

The right layout every time
===========================
Problem
-------
In matplotlib, the tight layout algorithm is clunky, and cannot implement variable spacing between subplot rows and columns.

Solution
--------
In ProPlot, the tight layout algorithm is more accurate because:

#. The new `~proplot.subplots.GridSpec` class permits variable spacing between rows and columns.
#. The `~proplot.subplots.GridSpec` spacing parameters are specified in physical units instead of figure-relative units.

.. The `~matplotlib.gridspec.GridSpec` class is useful for creating figures with complex subplot geometry.
.. Users want to control axes positions with gridspecs.
.. * Matplotlib permits arbitrarily many `~matplotlib.gridspec.GridSpec`\ s per figure. This greatly complicates the tight layout algorithm for little evident gain.

Also, to simplify the tight layout algorithm, ProPlot permits only *one* `~matplotlib.gridspec.GridSpec` per figure.

* When a `~matplotlib.gridspec.SubplotSpec` is passed to `~proplot.subplots.Figure.add_subplot`, the figure is locked to the `~matplotlib.gridspec.GridSpec` associated with the `~matplotlib.gridspec.SubplotSpec`.
* When an integer or tuple is passed to `~proplot.subplots.Figure.add_subplot`, the geometry implied by subsequent calls must *divide* or *multiply* the initial geometry -- for example, two square subplots above a longer rectangle subplot can be drawn by passing ``221``, ``222``, and ``212`` to `~proplot.subplots.Figure.add_subplot`.

This introduces a marginal limitation (see discussion in :pr:`50`) but *considerably* simplifies the tight layout algorithm. 

Colorbars and legends made easy
===============================
Problem
-------
Drawing colorbars and legends is pretty clumsy in matplotlib -- especially when trying to draw them outside of the figure. They can be too narrow, too wide, and mess up your subplot aspect ratios.

Solution
--------
ProPlot introduces a brand new engine for drawing colorbars and legends along the outside of:

#. Individual subplots.
#. Arbitrarily many subplots in the same row or column.

And since `~proplot.subplots.GridSpec` permits variable spacing between subplot rows and columns, ProPlot can easily make room for these objects, where matplotlib would have to create giant gaps if e.g. there was a colorbar between columns 1 and 2 but nothing between columns 2 and 3.

The `~proplot.axes.Axes` `~proplot.axes.Axes.colorbar` and `~proplot.axes.Axes.legend` commands and `~proplot.subplots.Figure` `~proplot.subplots.Figure.colorbar` and `~proplot.subplots.Figure.legend` commands are also overridden to make this process simple. `~proplot.subplots.Figure` colorbars are always centered relative to the subplot grid.

It is often useful to have colorbars, legends, and text spanning along multiple subplots on the outside of a figure. This is really tricky to do with the matplotlib API, or can result in misalignment between the subplots and the object since the items are centered relative to the **figure**, not the **subplot grid**.

A handy axes container
======================
.. The `~matplotlib.pyplot.subplots` command is useful for generating a scaffolding of * axes all at once. This is generally faster than successive `~matplotlib.subplots.Figure.add_subplot` commands.

Problem
-------
In matplotlib, `~matplotlib.pyplot.subplots` returns a 2D `~numpy.ndarray`, a 1D `~numpy.ndarray`, or the axes itself. This variable output is cumbersome.

Solution
--------
In ProPlot, `~proplot.subplots.subplots` returns an `~proplot.subplots.axes_grid` of axes that unifies the behavior of these three possible return values:

* `~proplot.subplots.axes_grid` is a `list` subclass that behaves like a scalar when it contains just one element, supports row-major or column-major 1D indexing (e.g. ``axs[0]``), and permits 2D indexing (e.g. ``axs[1,0]``) no matter the geometry.
* Since `~proplot.subplots.subplots` can generate figures with arbitrarily complex subplot geometry, this 2D indexing is useful only when the arrangement happens to be a clean 2D matrix.
* Thanks to the `~proplot.subplots.axes_grid.__getattr__` override, `~proplot.subplots.axes_grid` allows you to call arbitrary methods on arbitrary axes all at once, e.g. ``axs.format(tickminor=False)``.

Arbitrary units
===============
.. * Configuring spaces and dimensions in matplotlib often requires physical units.
Problem
-------
Matplotlib uses "inches" for figure dimensions and figure-relative or axes-relative units almost everywhere else. The problem is:

* Inches are foreign to the world outside of the U.S.
* Figure-relative and axes-relative units encourage "tinkering" with meaningless numbers that change the subjective appearance when the physical dimensions change, since *text* and *lines* are specified in the physical units "points".

Solution
--------
ProPlot permits arbitrary physical units for almost all sizing arguments, e.g. ``left='0.5cm'``. This is done by passing various keyword arguments through the `~proplot.utils.units` engine.

* This prevents "tinkering" and encourages users to be aware of the physical dimensions describing their figure.
* You can also use font-relative units, e.g. ``left='1em'``. This is nice when you don't care about physical dimensions, but need something more intuitive than figure-relative units.

.. * You can still use axes-relative and figure-relative units for most arguments with e.g. ``left='0.1fig'`` or ``left='0.1ax'``.
Working with colormaps
======================
Problem
-------
In matplotlib, colormaps are implemented with the `~matplotlib.colors.ListedColormap` and `~matplotlib.colors.LinearSegmentedColormap` classes. They are very hard to modify and hard to create.

Colormap identification by string name is also suboptimal. The names are case-sensitive, and reversed versions of each colormap (names that end in ``'_r'``) are not guaranteed to exist.

Solution
--------
In ProPlot, it is easy to generate, combine, and modify colormaps using the `~proplot.styletools.Colormap` constructor function, and thanks to the new `~proplot.styletools.ListedColormap`, `~proplot.styletools.LinearSegmentedColormap`, and `~proplot.styletools.PerceptuallyUniformColormap`. This includes new tools for making colormaps that vary linearly in `perceptually uniform colorspaces <https://en.wikipedia.org/wiki/HCL_color_space>`__, so that they portray your data accurately (see :ref:`Perceptually uniform colormaps` for details).

The `~proplot.styletools.CmapDict` dictionary used to store colormaps also makes colormap identification a bit easier. All colormap names are case-insensitive, and reversed colormaps are automatically created when you request a name ending in ``'_r'``.

.. Also, "colormaps" and "color cycles" are now *fluid*, e.g. you can use a colormap as the color cycler for line plots. This is ProPlot's answer to seaborn's "palettes".

.. -- matplotlib's "colormaps" and "property cyclers" are sufficient.

Working with property cycles
============================
Problem
-------
Changing the property cycle is tricky in matplotlib. You have to work with the :rcraw:`axes. prop_cycle` setting and the `~cycler.Cycler` class directly.

Solution
--------
In ProPlot, you can create arbitrary property cycles with `~proplot.styletools.Cycle` and use them with arbitrary plotting commands with the `cycle` keyword argument. You can also create property cycles from arbitrary colormaps! See `~proplot.styletools.Cycle` for details.

.. Changing the property cycle is easy in ProPlot.

More accurate colorbars
=======================
Problem
-------
In matplotlib, when ``extend='min'``, ``extend='max'``, or ``extend='neither'`` is passed to `~matplotlib.figure.Figure.colorbar` , colormap colors reserved for "out-of-bounds" values are simply truncated. But most of the time, users want to utilize the full colormap range.

The problem is that matplotlib divides colormap colors up into "levels" by generating a low-resolution lookup table (see `~matplotolib.colors.LinearSegmentedColormap` for details).
This approach cannot be fine-tuned, creates an unnecessary copy of the colormap, and prevents you from using the resulting colormap for plots with different numbers of levels.

It is clear that the task of "discretizing" the possible colormap colors should be left to the **normalizer**, not the colormap itself. Matplotlib provides `~matplotlib.colors.BoundaryNorm` for this purpose, but it is seldom used and its features are limited.

Solution
--------
In ProPlot, all colormap visualizations are automatically discretized with the `~proplot.styletools.BinNorm` class. This reads the ``extend`` property passed to your plotting command and chooses colormap indices so that your colorbar levels *always* traverse the full range of colormap colors.

`~proplot.styletools.BinNorm` can also apply an arbitrary continuous normalizer, e.g. `~proplot.styletools.LogNorm`, before discretization. Think of it as a "meta-normalizer" -- other normalizers perform the continuous transformation step, while this performs the discretization step.

Working with fonts
==================
Problem
-------
In matplotlib, the default font is DejaVu Sans. In this developer's not-at-all humble opinion, DejaVu Sans is fugly AF. It is also really tricky to work with custom fonts in matplotlib.

.. This font is not very aesthetically pleasing.

Solution
--------
In ProPlot, the default font is Helvetica. Albeit somewhat overused, this is a tried and tested, aesthetically pleasing sans serif font.

ProPlot also makes it easier to work with custom fonts by making use of a completely undocumented feature: Matplotlib adds ``.ttf`` and ``.otf`` font files in folders indicated by the ``TTFPATH`` environmental variable. ProPlot simply adds ``~/.proplot/fonts`` to ``TTFPATH`` -- feel free to drop font files in that directory, and you're good to go.

Cartopy and basemap defaults
============================
Problem
-------
In basemap and cartopy, the default coordinate system is always **map projection coordinates**. For most of us, this choice is very frustrating. Geophysical data is almost always stored in longitude-latitude or "Plate Carrée" coordinates.

Solution
--------
ProPlot makes the following changes:

* ``latlon=True`` is the default for `~proplot.axes.BasemapAxes` plotting methods.
* ``transform=ccrs.PlateCarree()`` is the default for `~proplot.axes.CartopyAxes` plotting methods.


