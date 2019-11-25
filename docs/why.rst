============
Why ProPlot?
============

ProPlot's core mission
is to improve upon the parts of matplotlib that
tend to be cumbersome or repetitive
for power users.
This page
enumerates these limitations and
describes how ProPlot addresses them.
To start using these new features, see
see :ref:`Quick overview` and the User Guide.

No more boilerplate
===================

.. raw:: html

   <h3>Problem</h3>

Power users often need to change lots of plot settings all at once. In matplotlib, this requires a bunch of one-liner setters and getters, like `~matplotlib.axes.Axes.set_title`. 

This workflow is verbose and often confusing. It can be unclear whether settings can be changed from a `~matplotlib.figure.Figure` setter, an `~matplotlib.axes.Axes` setter, an `~matplotlib.axis.XAxis` or `~matplotlib.axis.YAxis` setter, or a miscellaneous bulk function like `~matplotlib.axes.Axes.tick_params`.

.. raw:: html

   <h3>Solution</h3>

ProPlot introduces the `~proplot.axes.Axes.format` method for changing arbitrary settings **in bulk**.
Think of this as an improved version of the `~matplotlib.artist.Artist` `~matplotlib.artist.Artist.update` method.
This massively reduces the amount of code needed to create highly customized figures. For example, it is trivial to see that

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
.. raw:: html

   <h3>Problem</h3>

Matplotlib and cartopy introduce a bunch of classes with verbose names like `~matplotlib.ticker.MultipleLocator`, `~matplotlib.ticker.FormatStrFormatter`, `~cartopy.crs.AlbersEqualArea`, and `~cartopy.crs.NorthPolarStereo`.
But matplotlib `axis scales <https://matplotlib.org/3.1.0/gallery/scales/scales.html>`__ and `colormaps <https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html>`__ are **registered** under string names, as are `map projections <https://matplotlib.org/basemap/users/mapsetup.html>`__ in basemap. So why not register other verbose class names too?

.. raw:: html

   <h3>Solution</h3>

In ProPlot, tick locators, tick formatters, axis scales, cartopy projections, colormaps, and property cyclers are all **registered**. This is done by creating **constructor functions** and passing various keyword argument through the functions, allowing users to use class names instead of instances.

This may seem "unpythonic" but it is absolutely invaluable when making plots. Plotting code has a half life of about 30 seconds, and typing out all of these extra class names and import statements gets to be a major drag.

The below table lists the constructor functions and the keyword arguments that use them.

==============================  =============================  ================================================================================================================================================================================================
Function                        Returns                        Interpreted by
==============================  =============================  ================================================================================================================================================================================================
`~proplot.axistools.Locator`    Axis locator                   ``locator=``, ``xlocator=``, ``ylocator=``, ``minorlocator=``, ``xminorlocator=``, ``yminorlocator=``, ``ticks=``, ``xticks=``, ``yticks=``, ``minorticks=``, ``xminorticks=``, ``yminorticks=``
`~proplot.axistools.Formatter`  Axis formatter                 ``formatter=``, ``xformatter=``, ``yformatter=``, ``ticklabels=``, ``xticklabels=``, ``yticklabels=``
`~proplot.axistools.Scale`      Axis scale                     ``xscale=``, ``yscale=``
`~proplot.styletools.Colormap`  Colormap                       ``cmap=``
`~proplot.styletools.Cycle`     Property cycler                ``cycle=``
`~proplot.styletools.Norm`      Colormap normalizer            ``norm=``
`~proplot.projs.Proj`           Cartopy or basemap projection  ``proj=``
==============================  =============================  ================================================================================================================================================================================================



Nice-looking figures out of the box
===================================

.. raw:: html

   <h3>Problem</h3>

Matplotlib plots tend to require lots of "tweaking" when you have more than one subplot in the figure. This is partly because you must specify the physical dimensions of the figure, while the dimensions of the *individual subplots* are more important:

#. The subplot aspect ratio is usually more relevant than the figure aspect ratio, e.g. for map projections.
#. The subplot width and height control the evident thickness of text and other content plotted inside the axes.

Matplotlib has a `tight layout <https://matplotlib.org/tutorials/intermediate/tight_layout_guide.html>`__ algorithm to keep you from having to "tweak" the spacing, but the algorithm cannot apply different amounts of spacing between different subplot row and column boundaries. This sometimes creates figures with unnecessary extra whitespace.

.. raw:: html

   <h3>Solution</h3>

In ProPlot, you can specify the physical dimensions of *subplots* instead of the figure by passing `axwidth` or `axheight` to `~proplot.subplots.Figure`. The default behavior is ``axwidth=2`` (inches). Figure dimensions are then automatically calculated to accommodate the subplot geometry and the spacing adjustments.

..
   Several matplotlib backends require figure dimensions to be fixed. When `~proplot.subplots.Figure.draw` changes the figure dimensions, this can "surprise" the backend and cause unexpected behavior. ProPlot fixes this issue for the static inline backend and the Qt popup backend. However, this issue is unfixable the "notebook" inline backend, the "macosx" popup backend, and possibly other untested backends.

The tight layout algorithm is also simpler and more accurate because:

#. The new `~proplot.subplots.FlexibleGridSpec` class permits variable spacing between rows and columns.
#. The `~proplot.subplots.FlexibleGridSpec` spacing parameters are specified in physical units instead of figure-relative units.
#. Figures are restricted to have only *one* `~proplot.subplots.FlexibleGridSpec` per figure. This is done by requiring users to draw all of their subplots at once with `~proplot.subplots.subplots`. This requirement *considerably* simplifies the algorithm (see :pr:`50` for details).

See :ref:`Figure tight layout` for details.

..
   The `~matplotlib.gridspec.FlexibleGridSpec` class is useful for creating figures with complex subplot geometry.
..
   Users want to control axes positions with gridspecs.
..
   * Matplotlib permits arbitrarily many `~matplotlib.gridspec.FlexibleGridSpec`\ s per figure. This greatly complicates the tight layout algorithm for little evident gain.
..
   ProPlot introduces a marginal limitation (see discussion in :pr:`50`) but *considerably* simplifies the tight layout algorithm.

Simpler colorbars and legends
=============================

.. raw:: html

   <h3>Problem</h3>

In matplotlib, it is hard to put colorbars and legends on the outside of subplots. It can end up messing up subplot aspect ratios, and colorbars tend to be too narrow or too wide.

It is *also* difficult to draw colorbars and legends that span the figure edge or serve as reference to more than one subplot.
Because this requires so much tinkering, most users just add identical colorbars
to every single subplot, which is incredibly repetitive!

..
   Drawing colorbars and legends is pretty clumsy in matplotlib -- especially when trying to draw them outside of the figure. They can be too narrow, too wide, and mess up your subplot aspect ratios.

.. raw:: html

   <h3>Solution</h3>

ProPlot introduces a brand new engine for drawing colorbars and legends along the outside of
individual subplots and along contiguous subplots on the edge of the figure:

* The `~proplot.axes.Axes` `~proplot.axes.Axes.legend` command and the `~proplot.subplots.Figure` `~proplot.subplots.Figure.colorbar` and `~proplot.subplots.Figure.legend` commands are overridden, adding various new features.
* There is a new `~proplot.axes.Axes` `~proplot.axes.Axes.colorbar` method for drawing *inset* colorbars or adding colorbars along the outer edge of axes.
* The `~proplot.subplots.Figure` `~proplot.subplots.Figure.colorbar` and `~proplot.subplots.Figure.legend` commands draw colorbars and legends that are centered relative to the *subplot grid*, not the axes. This is critical if your left-right or top-bottom border padding is asymmetric.
* You can put colorbars and legends along the edge of axes or along the edge of the whole figure by passing ``loc='l'``, ``loc='r'``, ``loc='b'``, or ``loc='t'`` to the colorbar and legend commands.
* Outer colorbars and legends don't mess up the subplot layout or subplot aspect ratios, since the new `~proplot.subplots.FlexibleGridSpec` class permits variable spacing between subplot rows and columns. This is critical e.g. if you have a colorbar between columns 1 and 2 but nothing between columns 2 and 3.
* The width of colorbars are now specified in physical units. This makes it easier to get the thickness just right, and makes thickness independent of figure size.

A useful axes container
=======================

..
   The `~matplotlib.pyplot.subplots` command is useful for generating a scaffolding of * axes all at once. This is generally faster than successive `~matplotlib.subplots.Figure.add_subplot` commands.

.. raw:: html

   <h3>Problem</h3>

In matplotlib, `~matplotlib.pyplot.subplots` returns a 2D `~numpy.ndarray`, a 1D `~numpy.ndarray`, or the axes itself. This variable output is cumbersome.

.. raw:: html

   <h3>Solution</h3>

In ProPlot, `~proplot.subplots.subplots` returns an `~proplot.subplots.axes_grid` of axes that unifies the behavior of these three possible return values:

* `~proplot.subplots.axes_grid` is a `list` subclass that behaves like a scalar when it contains just one element.
* `~proplot.subplots.axes_grid` supports row-major or column-major 1D indexing, e.g. ``axs[0]``. The order can be changed by passing ``order='F'`` or ``order='C'`` to `~proplot.subplots.subplots`.
* `~proplot.subplots.axes_grid` permits 2D indexing, e.g. ``axs[1,0]``. Since `~proplot.subplots.subplots` can generate figures with arbitrarily complex subplot geometry, this 2D indexing is useful only when the arrangement happens to be a clean 2D matrix.

Further, thanks to the `~proplot.subplots.axes_grid.__getattr__` override, `~proplot.subplots.axes_grid` allows you to call arbitrary methods on arbitrary axes all at once, e.g. ``axs.format(tickminor=False)``.

Xarray and pandas integration
=============================

.. raw:: html

   <h3>Problem</h3>

Matplotlib strips metadata from the array-like `xarray` `~xarray.DataArray` container and the `pandas` `~pandas.DataFrame` and `~pandas.Series` containers. To create plots that
are automatically labeled with metadata from these containers, you need to use
the dedicated `xarray plotting <http://xarray.pydata.org/en/stable/plotting.html>`__ and `pandas plotting <https://pandas.pydata.org/pandas-docs/stable/user_guide/visualization.html>`__ tools.

This approach is somewhat cumbersome -- plotting methods should be invoked on the axes, not on the data container! It also requires learning a slightly different syntax, and tends to encourage using the `~matplotlib.pyplot` API rather than the object-oriented API.

.. raw:: html

   <h3>Solution</h3>

ProPlot *reproduces* most of the `xarray.DataArray.plot` and `pandas.DataFrame.plot` features on the `~proplot.axes.Axes`
plotting methods themselves! Axis tick labels, axis labels, subplot titles, and colorbar and legend labels are automatically applied
when a `~xarray.DataArray`, `~pandas.DataFrame`, or `~pandas.Series` is passed through
a plotting method instead of a `~numpy.ndarray`.
This is accomplished by passing positional arguments through the
`~proplot.wrappers.standardize_1d` and `~proplot.wrappers.standardize_2d`
wrappers.

Various plotting improvements
=============================

.. raw:: html

   <h3>Problem</h3>

Certain plotting tasks are quite difficult to accomplish
with the default matplotlib API. The `seaborn`, `xarray`, and `pandas`
packages offer improvements, but it would be nice
to have this functionality build right into matplotlib.

.. raw:: html

   <h3>Solutions</h3>

The ProPlot `~proplot.axes.Axes` class
wraps various plotting methods to reproduce
certain features from `seaborn`, `xarray`, and `pandas`:

* The new `~proplot.axes.Axes.heatmap` command draws `~matplotlib.axes.Axes.pcolormesh` plots and puts ticks at the center of each box.
* The `~matplotlib.axes.Axes.bar` and `~matplotlib.axes.Axes.barh` commands now accept 2D arrays, and can *stack* or *group* successive columns of data, thanks to `~proplot.wrappers.bar_wrapper`.
* The new `~proplot.axes.Axes.area` and `~proplot.axes.Axes.areax` commands mimic the `~proplot.axes.Axes.fill_between` and `~proplot.axes.Axes.fill_betweenx` commands, but also support drawing *stacked* area plots for 2D arrays.

`~proplot.axes.Axes` also includes the following
new plotting features:

* `~matplotlib.axes.Axes.pcolor` and `~matplotlib.axes.Axes.pcolormesh` plots use auto-generated coordinate *edges* if you pass coordinate *centers*.
* `~proplot.axes.Axes.area` plots can be assigned different colors for negative and positive values. This will also be added to `~matplotlib.axes.Axes.bar` soon.
* `~matplotlib.axes.Axes.pcolor`, `~matplotlib.axes.Axes.pcolormesh`, `~proplot.axes.Axes.heatmap`,  `~matplotlib.axes.Axes.contour` and `~matplotlib.axes.Axes.contourf` plots can be assigned contour and box labels by simply passing ``labels=True`` to the plotting command.
* `~matplotlib.axes.Axes.plot`, `~matplotlib.axes.Axes.scatter`, and `~matplotlib.axes.Axes.bar` plots can be assigned error bars using a variety of `~proplot.wrappers.errorbar_wrapper` keyword args.
* `~proplot.axes.Axes.parametric` plots can be made with colormap colors marking the parametric coordinates rather than text annotations.
* `~matplotlib.axes.Axes.pcolor`, `~matplotlib.axes.Axes.pcolormesh`, `~proplot.axes.Axes.heatmap`,  `~matplotlib.axes.Axes.contour` and `~matplotlib.axes.Axes.contourf` plots on geographic axes can be inteprolated to global coverage by passing ``globe=True`` tot he plotting command.

See :ref:`1d plotting commands` and :ref:`2d plotting commands`
for details.

Cartopy and basemap integration
===============================

.. raw:: html

   <h3>Problem</h3>

There are two widely-used engines
for plotting geophysical data with matplotlib: `cartopy` and `~mpl_toolkits.basemap`.
Using cartopy tends to be quite verbose and involve lots of boilerplate code,
while basemap is outdated and requires you to use plotting commands on a separate `~mpl_toolkits.basemap.Basemap` object.

Also, `cartopy` and `~mpl_toolkits.basemap` plotting commands assume *map projection coordinates* unless specified otherwise. For most of us, this choice is very frustrating, since geophysical data are usually stored in longitude-latitude or "Plate Carrée" coordinates.

.. raw:: html

   <h3>Solution</h3>

ProPlot includes various `cartopy` and `~mpl_toolkits.basemap` features
using the `~proplot.axes.ProjAxes` class. The corresponding `~proplot.axes.ProjAxes.format` command lets you apply all kinds of geographic plot settings, like coastlines, continents, political boundaries, and meridian and parallel gridlines.
It also makes longitude-latitude coordinates the *default*:

* ``latlon=True`` is the default for `~proplot.axes.BasemapAxes` plotting methods.
* ``transform=ccrs.PlateCarree()`` is the default for `~proplot.axes.GeoAxes` plotting methods.

Note that the basemap developers plan to `halt active development after 2020 <https://matplotlib.org/basemap/users/intro.html#cartopy-new-management-and-eol-announcement>`__, since cartopy is integrated more closely with the matplotlib API and has more room for growth. For now, cartopy is `missing several features <https://matplotlib.org/basemap/api/basemap_api.html#module-mpl_toolkits.basemap>`__ offered by basemap -- namely, flexible meridian and parallel gridline labels, drawing physical map scales, and convenience features for adding background images like the "blue marble". But once these are added to cartopy, ProPlot support for basemap may be removed.


Colormaps and property cycles
=============================

.. raw:: html

   <h3>Problem</h3>

In matplotlib, colormaps are implemented with the `~matplotlib.colors.ListedColormap` and `~matplotlib.colors.LinearSegmentedColormap` classes. They are very hard to modify and hard to create. Colormap identification by string name is also suboptimal. The names are case-sensitive, and reversed versions of each colormap (i.e. names that end in ``'_r'``) are not guaranteed to exist.

.. raw:: html

   <h3>Solution</h3>

In ProPlot, it is easy to generate, combine, and modify colormaps using the `~proplot.styletools.Colormap` constructor function, and thanks to the new `~proplot.styletools.ListedColormap`, `~proplot.styletools.LinearSegmentedColormap`, and `~proplot.styletools.PerceptuallyUniformColormap`. This includes new tools for making colormaps that are `perceptually uniform <https://en.wikipedia.org/wiki/HCL_color_space>`__ (see :ref:`Perceptually uniform colormaps` for details).
The `~proplot.styletools.CmapDict` dictionary used to store colormaps also makes colormap identification a bit easier. All colormap names are case-insensitive, and reversed colormaps are automatically created when you request a name ending in ``'_r'``.

In ProPlot, you can also create arbitrary property cycles with `~proplot.styletools.Cycle` and use them with arbitrary plotting commands with the `cycle` keyword argument. You can also create property cycles from arbitrary colormaps! See `~proplot.styletools.Cycle` for details.

Improved colormap normalization
===============================
.. raw:: html

   <h3>Problem</h3>

In matplotlib, when ``extend='min'``, ``extend='max'``, or ``extend='neither'`` is passed to `~matplotlib.figure.Figure.colorbar` , colormap colors reserved for "out-of-bounds" values are truncated. The problem is that matplotlib discretizes colormaps by generating a low-resolution lookup table (see `~matplotlib.colors.LinearSegmentedColormap` for details).
This approach cannot be fine-tuned, creates an unnecessary copy of the colormap, and prevents you from using the resulting colormap for plots with different numbers of levels.

It is clear that the task discretizing colormap colors should be left to the **normalizer**, not the colormap itself. Matplotlib provides `~matplotlib.colors.BoundaryNorm` for this purpose, but it is seldom used and its features are limited.

.. raw:: html

   <h3>Solution</h3>

In ProPlot, all colormap visualizations are automatically discretized with the `~proplot.styletools.BinNorm` class. This reads the `extend` property passed to your plotting command and chooses colormap indices so that your colorbar levels *always* traverse the full range of colormap colors.

`~proplot.styletools.BinNorm` can also apply an arbitrary continuous normalizer, e.g. `~matplotlib.colors.LogNorm`, before discretization. Think of it as a "meta-normalizer" -- other normalizers perform the continuous transformation step, while this performs the discretization step.

Working with fonts
==================
.. raw:: html

   <h3>Problem</h3>

In matplotlib, the default font is DejaVu Sans. In this developer's humble opinion, DejaVu Sans is fugly AF. It is also really tricky to work with custom fonts in matplotlib.

..
   This font is not very aesthetically pleasing.

.. raw:: html

   <h3>Solution</h3>

In ProPlot, the default font is Helvetica. Albeit somewhat overused, this is a tried and tested, aesthetically pleasing sans serif font.

ProPlot also makes it easier to work with custom fonts by making use of a completely undocumented feature: ``$TTFPATH``. Matplotlib adds ``.ttf`` and ``.otf`` font files in folders listed in the ``$TTFPATH`` environment variable, so ProPlot simply populates that variable. Feel free to  drop your own font files into ``~/.proplot/fonts``, and you're good to go.


