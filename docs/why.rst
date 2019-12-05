============
Why ProPlot?
============

ProPlot's core mission
is to improve upon the parts of matplotlib that
tend to be cumbersome or repetitive
for power users.
This page enumerates the stickiest of these limitations
and describes how ProPlot addresses them.

..
   This page is not comprehensive --
   see the User Guide for a comprehensive overview
   with worked examples.

..
   To start using these new features, see
   see :ref:`Usage overview` and the User Guide.

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

...is much more succinct than

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

Matplotlib and cartopy introduce a bunch of classes with verbose names like `~matplotlib.ticker.MultipleLocator`, `~matplotlib.ticker.FormatStrFormatter`, and
`~cartopy.crs.LambertAzimuthalEqualArea`. Since plotting code has a half life of about 30 seconds, typing out all of these extra class names and import statements can be a *major* drag.

*Other* parts of the matplotlib API were designed with this in mind.
`Native axes projections <https://matplotlib.org/3.1.1/api/projections_api.html>`__,
`axis scale classes <https://matplotlib.org/3.1.0/gallery/scales/scales.html>`__,
`box style classes <https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.patches.FancyBboxPatch.html?highlight=boxstyle>`__, `arrow style classes <https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.patches.FancyArrowPatch.html?highlight=arrowstyle>`__,
`arc style classes <https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.patches.ConnectionStyle.html?highlight=class%20name%20attrs>`__, and
`backends <https://matplotlib.org/faq/usage_faq.html#what-is-a-backend>`__
are referenced with "registered" string names,
as are `basemap projection types <https://matplotlib.org/basemap/users/mapsetup.html>`__.
If these are already "registered", why not "register" everything else?

.. raw:: html

   <h3>Solution</h3>

In ProPlot, tick locators, tick formatters, axis scales, cartopy projections, colormaps, and property cyclers are all "registered". This is done by creating several **constructor functions** and passing various keyword argument *through* the constructor functions.
This may seem "unpythonic" but it is absolutely invaluable when making plots.

The constructor functions also accept *other* types of input: For
example, scalar numbers passed to `~proplot.axistools.Locator` returns
a `~matplotlib.ticker.MultipleLocator` instance, lists of strings passed
to `~proplot.axistools.Formatter` returns a `~matplotlib.ticker.FixedFormatter` instance, and `~proplot.styletools.Colormap` and `~proplot.styletools.Cycle` accept all kinds of input. And if an instance of the corresponding class is passed to any constructor function, it is simply returned. See :ref:`X and Y axis settings`, :ref:`Colormaps`, and :ref:`Color cycles` for details.

The below table lists the constructor functions and the keyword arguments that
use them.

==============================  ============================================================  =============================================================  =================================================================================================================================================================================================
Function                        Returns                                                       Used by                                                        Keyword argument(s)
==============================  ============================================================  =============================================================  =================================================================================================================================================================================================
`~proplot.axistools.Locator`    Axis `~matplotlib.ticker.Locator`                             `~proplot.axes.Axes.format` and `~proplot.axes.Axes.colorbar`  ``locator=``, ``xlocator=``, ``ylocator=``, ``minorlocator=``, ``xminorlocator=``, ``yminorlocator=``, ``ticks=``, ``xticks=``, ``yticks=``, ``minorticks=``, ``xminorticks=``, ``yminorticks=``
`~proplot.axistools.Formatter`  Axis `~matplotlib.ticker.Formatter`                           `~proplot.axes.Axes.format` and `~proplot.axes.Axes.colorbar`  ``formatter=``, ``xformatter=``, ``yformatter=``, ``ticklabels=``, ``xticklabels=``, ``yticklabels=``
`~proplot.axistools.Scale`      Axis `~matplotlib.scale.ScaleBase`                            `~proplot.axes.Axes.format`                                    ``xscale=``, ``yscale=``
`~proplot.styletools.Cycle`     Property `~cycler.Cycler`                                     1d plotting methods                                            ``cycle=``
`~proplot.styletools.Colormap`  `~matplotlib.colors.Colormap` instance                        2d plotting methods                                            ``cmap=``
`~proplot.styletools.Norm`      `~matplotlib.colors.Normalize` instance                       2d plotting methods                                            ``norm=``
`~proplot.projs.Proj`           `~cartopy.crs.Projection` or `~mpl_toolkits.basemap.Basemap`  `~proplot.subplots.subplots`                                   ``proj=``
==============================  ============================================================  =============================================================  =================================================================================================================================================================================================

Note that `~matplotlib.axes.Axes.set_xscale` and `~matplotlib.axes.Axes.set_yscale`
now accept instances of `~matplotlib.scale.ScaleBase` thanks to a monkey patch
applied by ProPlot.

Automatic dimensions and spacing
================================

.. raw:: html

   <h3>Problem</h3>

Matplotlib plots tend to require lots of "tweaking" when you have more than one subplot in the figure. This is partly because you must specify the physical dimensions of the figure, while the dimensions of the *individual subplots* are more important:

#. The subplot aspect ratio is usually more relevant than the figure aspect ratio, e.g. for map projections.
#. The subplot width and height control the evident thickness of text and other content plotted inside the axes.

Matplotlib has a `tight layout <https://matplotlib.org/tutorials/intermediate/tight_layout_guide.html>`__ algorithm to keep you from having to "tweak" the spacing, but the algorithm cannot apply different amounts of spacing between different subplot row and column boundaries. This is a silly limitation that often results in unnecessary whitespace, and can be a major problem when you want to put e.g. a legend on the outside of a subplot.

.. raw:: html

   <h3>Solution</h3>

In ProPlot, you can specify the physical dimensions of a *reference subplot* instead of the figure by passing `axwidth`, `axheight`, and/or `aspect` to `~proplot.subplots.Figure`. The default behavior is ``aspect=1`` and ``axwidth=2`` (inches). If the `aspect ratio mode <https://matplotlib.org/2.0.2/examples/pylab_examples/equal_aspect_ratio.html>`__ for the reference subplot is set to ``'equal'``, as with :ref:`Geographic and polar plots` and `~matplotlib.axes.Axes.imshow` plots, the existing aspect will be used instead.

Figure dimensions are calculated to accommodate subplot geometry and tight layout spacing as follows:

* If `axwidth` or `axheight` are used, the figure width and height are calculated automatically.
* If `width` is used, the figure height is calculated automatically.
* If `height` is used, the figure width is calculated automatically.
* If `width` *and* `height` or `figsize` is used, the figure dimensions are fixed.

..
   Several matplotlib backends require figure dimensions to be fixed. When `~proplot.subplots.Figure.draw` changes the figure dimensions, this can "surprise" the backend and cause unexpected behavior. ProPlot fixes this issue for the static inline backend and the Qt popup backend. However, this issue is unfixable the "notebook" inline backend, the "macosx" popup backend, and possibly other untested backends.

By default, ProPlot also uses a custom tight layout algorithm that automatically determines the `left`, `right`, `bottom`, `top`, `wspace`, and `hspace` `~matplotlib.gridspec.GridSpec` parameters. This algorithm is simpler and more precise because:

#. The new `~proplot.subplots.GridSpec` class permits variable spacing between rows and columns. It turns out this is *critical* for putting :ref:`Colorbars and legends` on the outside of subplots.
#. Figures are restricted to have only *one* `~proplot.subplots.GridSpec` per figure. This is done by requiring users to draw all of their subplots at once with `~proplot.subplots.subplots`, and it *considerably* simplifies the algorithm (see :pr:`50` for details).

See :ref:`Subplots features` for details.

..
   #. The `~proplot.subplots.GridSpec` spacing parameters are specified in physical units instead of figure-relative units.

..
   The `~matplotlib.gridspec.GridSpec` class is useful for creating figures with complex subplot geometry.
..
   Users want to control axes positions with gridspecs.
..
   * Matplotlib permits arbitrarily many `~matplotlib.gridspec.GridSpec`\ s per figure. This greatly complicates the tight layout algorithm for little evident gain.
..
   ProPlot introduces a marginal limitation (see discussion in :pr:`50`) but *considerably* simplifies the tight layout algorithm.

Outer colorbars and legends
===========================

.. raw:: html

   <h3>Problem</h3>

In matplotlib, it is difficult to draw `~matplotlib.figure.Figure.colorbar`\ s and
`~matplotlib.axes.Axes.legend`\ s on the outside of subplots. It is very easy to mess up the subplot aspect ratios and the colorbar widths. It is even *more* difficult to draw `~matplotlib.figure.Figure.colorbar`\ s and `~matplotlib.figure.Figure.legend`\ s that reference more than one subplot:

* Matplotlib has no capacity for drawing colorbar axes that span multiple plots -- you have to create the axes yourself. This requires so much tinkering that most users just add identical colorbars to every single subplot!
* Legends that span multiple plots tend to require *manual* positioning and tinkering with the `~matplotlib.gridspec.GridSpec` spacing, just like legends placed outside of individual subplots.

..
   The matplotlib example for `~matplotlib.figure.Figure` legends is `not pretty <https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/figlegend_demo.html>`__.

..
   Drawing colorbars and legends is pretty clumsy in matplotlib -- especially when trying to draw them outside of the figure. They can be too narrow, too wide, and mess up your subplot aspect ratios.

.. raw:: html

   <h3>Solution</h3>

ProPlot introduces a *brand new engine* for drawing colorbars and legends along the outside of
individual subplots and along contiguous subplots on the edge of the figure:

* Passing ``loc='l'``, ``loc='r'``, ``loc='b'``, or ``loc='t'`` to `~proplot.axes.Axes` `~proplot.axes.Axes.colorbar` or `~proplot.axes.Axes` `~proplot.axes.Axes.legend` draws the colorbar or legend along the outside of the axes.
* Passing ``loc='l'``, ``loc='r'``, ``loc='b'``, or ``loc='t'`` to `~proplot.subplots.Figure` `~proplot.subplots.Figure.colorbar` and `~proplot.subplots.Figure.legend` draws the colorbar or legend along the edge of the figure, centered relative to the *subplot grid* rather than figure coordinates.
* Outer colorbars and legends don't mess up the subplot layout or subplot aspect ratios, since `~proplot.subplots.GridSpec` permits variable spacing between subplot rows and columns. This is critical e.g. if you have a colorbar between columns 1 and 2 but nothing between columns 2 and 3.
* `~proplot.subplots.Figure` and `~proplot.axes.Axes` colorbar widths are specified in *physical* units rather than relative units. This makes colorbar thickness independent of figure size and easier to get just right.

The colorbar and legend commands also add several new features, like colorbars-from-lines and centered-row legends. And to make `~proplot.axes.Axes` `~proplot.axes.Axes.colorbar` consistent with `~proplot.axes.Axes` `~proplot.axes.Axes.legend`, you can also now draw *inset* colorbars. See :ref:`Colorbars and legends` for details.

The axes container class
========================

..
   The `~matplotlib.pyplot.subplots` command is useful for generating a scaffolding of * axes all at once. This is generally faster than successive `~matplotlib.subplots.Figure.add_subplot` commands.

.. raw:: html

   <h3>Problem</h3>

In matplotlib, `~matplotlib.pyplot.subplots` returns a 2D `~numpy.ndarray`, a 1D `~numpy.ndarray`, or the axes itself. This inconsistent behavior can be confusing.

.. raw:: html

   <h3>Solution</h3>

In ProPlot, `~proplot.subplots.subplots` returns an `~proplot.subplots.subplot_grid` container filled with axes objects. This container lets you call arbitrary methods on arbitrary subplots all at once, which can be useful when you want to style your subplots identically (e.g. ``axs.format(tickminor=False)``).
See :ref:`The basics` for details.

The `~proplot.subplots.subplot_grid` class also
unifies the behavior of the three possible `matplotlib.pyplot.subplots` return values:

* `~proplot.subplots.subplot_grid` permits 2d indexing, e.g. ``axs[1,0]``. Since `~proplot.subplots.subplots` can generate figures with arbitrarily complex subplot geometry, this 2d indexing is useful only when the arrangement happens to be a clean 2d matrix.
* Since `~proplot.subplots.subplot_grid` is a `list` subclass, it also supports 1d indexing, e.g. ``axs[1]``. The default order can be switched from row-major to column-major by passing ``order='F'`` to `~proplot.subplots.subplots`.
* `~proplot.subplots.subplot_grid` behaves like a scalar when it contains just one element. So if you just made a single axes with ``f, axs = plot.subplots()``, calling ``axs[0].command`` is equivalent to ``axs.command``.

..
   This goes with ProPlot's theme of preserving the object-oriented spirit, but making things easier for users.

New and improved plotting methods
=================================

.. raw:: html

   <h3>Problem</h3>

Certain plotting tasks are quite difficult to accomplish
with the default matplotlib API. The `seaborn`, `xarray`, and `pandas`
packages offer improvements, but it would be nice
to have this functionality build right into matplotlib.
There is also room for improvement that none of these packages
address.

..
   Matplotlib also has some finicky plotting issues
   that normally requires
..
   For example, when you pass coordinate *centers* to `~matplotlib.axes.Axes.pcolor` and `~matplotlib.axes.Axes.pcolormesh`, they are interpreted as *edges* and the last column and row of your data matrix is ignored. Also, to add labels to `~matplotlib.axes.Axes.contour` and `~matplotlib.axes.Axes.contourf`, you need to call a dedicated `~matplotlib.axes.Axes.clabel` method instead of just using a keyword argument.


.. raw:: html

   <h3>Solution</h3>


The ProPlot `~proplot.axes.Axes` class
overrides various plotting methods to make
your life easier.

..
  ProPlot also provides
  *constistent behavior* when
  switching between different commands, for
  example `~matplotlib.axes.Axes.plot` and `~matplotlib.axes.Axes.scatter`
  or `~matplotlib.axes.Axes.contourf` and `~matplotlib.axes.Axes.pcolormesh`.

..
   ProPlot also uses wrappers to *unify* the behavior of various
   plotting methods.

..
  All positional arguments for "1d" plotting methods are standardized by `~proplot.wrappers.standardize_1d`. All positional arguments for "2d" plotting methods are standardized by `~proplot.wrappers.standardize_2d`. See :ref:`1d plotting` and :ref:`2d plotting` for details.

* Just like `~proplot.axes.Axes.format`, certain keyword arguments for 1d and 2d plotting commands are passed through :ref:`Constructor functions`.

      * All 1d plotting methods accept a `cycle` keyword argument interpreted by `~proplot.styletools.Cycle`. See :ref:`Color cycles` for details.
      * All 2d plotting methods accept a `cmap` keyword argument interpreted by `~proplot.styletools.Colormap` and a `norm` keyword argument interpreted by `~proplot.styletools.Norm`. See :ref:`Colormaps` for details.

* All 2d plotting methods now accept a `labels` keyword argument. This is used to draw contour labels or grid box labels on heatmap plots. See :ref:`2d plotting` for details.
* ProPlot fixes the irritating `white-lines-between-filled-contours <https://stackoverflow.com/q/8263769/4970632>`__, `white-lines-between-pcolor-patches <https://stackoverflow.com/q/27092991/4970632>`__, and `white-lines-between-colorbar-patches <https://stackoverflow.com/q/15003353/4970632>`__ observed when saving vector graphics.
* Matplotlib requires coordinate *centers* for contour plots and *edges* for pcolor plots. If you pass *centers* to pcolor, matplotlib treats them as *edges* and silently trims one row/column of your data. ProPlot changes this behavior:

    * If edges are passed to `~matplotlib.axes.Axes.contour` or `~matplotlib.axes.Axes.contourf`, centers are *calculated* from the edges
    * If centers are passed to `~matplotlib.axes.Axes.pcolor` or `~matplotlib.axes.Axes.pcolormesh`, edges are *estimated* from the centers.


ProPlot also duplicates various `seaborn`, `xarray`, and `pandas` features
using *new* axes methods and by overriding *old* axes
methods.

=====================================  ====  =============================================================================================================================================================================================
Plotting method                        New?  Description
=====================================  ====  =============================================================================================================================================================================================
`~proplot.axes.Axes.area`              ✓     Alias for `~matplotlib.axes.Axes.fill_between`.
`~proplot.axes.Axes.areax`             ✓     Alias for `~matplotlib.axes.Axes.fill_betweenx`.
`~proplot.axes.Axes.parametric`        ✓     Draws *parametric* line plots, where the parametric coordinate is denoted with colormap colors.
`~matplotlib.axes.Axes.bar`            ✗     Now accepts 2D arrays, *stacks* or *groups* successive columns. Soon will be able to use different colors for positive/negative data.
`~matplotlib.axes.Axes.barh`           ✗     As with `~matplotlib.axes.Axes.bar`, but for horizontal bars.
`~matplotlib.axes.Axes.fill_between`   ✗     Now accepts 2D arrays, *stacks* or *overlays* successive columns. Also can use different colors for positive/negative data.
`~matplotlib.axes.Axes.fill_betweenx`  ✗     As with `~matplotlib.axes.Axes.fill_between`, but for horizontal fills.
`~proplot.axes.Axes.heatmap`           ✓     Invokes `~matplotlib.axes.Axes.pcolormesh` and puts ticks at the center of each box.
`~matplotlib.axes.Axes.contour`        ✗     Add contour labels by passing ``labels=True``, interpolate to global coverage for `~proplot.axes.ProjAxes` by passing ``globe=True``.
`~matplotlib.axes.Axes.contourf`       ✗     As with `~matplotlib.axes.Axes.contour`. Labels are colored black or white according to the filled contour's luminance.
`~matplotlib.axes.Axes.pcolor`         ✗     Add gridbox labels by passing ``labels=True``, interpolate to global coverage by passing ``globe=True``.  Labels are colored black or white according to the underlying box's luminance.
`~matplotlib.axes.Axes.pcolormesh`     ✗     As with `~matplotlib.axes.Axes.pcolor`.
=====================================  ====  =============================================================================================================================================================================================


..
   * `~proplot.axes.Axes.area` plots can be assigned different colors for negative and positive values. This will also be added to `~matplotlib.axes.Axes.bar` soon.

Xarray and pandas integration
=============================

.. raw:: html

   <h3>Problem</h3>

When you pass the array-like `xarray.DataArray`, `pandas.DataFrame`, and `pandas.Series` containers to matplotlib plotting commands, the metadata is ignored. To create plots that are automatically labeled with this metadata, you must use
the dedicated `xarray.DataArray.plot`, `pandas.DataFrame.plot`, and `pandas.Series.plot`
tools instead.

This approach is not ideal -- plotting methods should be invoked on the `~proplot.axes.Axes`, not on the data container! It also requires learning a slightly different syntax, and tends to encourage using the `~matplotlib.pyplot` API rather than the object-oriented API.

.. raw:: html

   <h3>Solution</h3>

ProPlot *reproduces* most of the `xarray.DataArray.plot`, `pandas.DataFrame.plot`, and `pandas.Series.plot` features on the `~proplot.axes.Axes` methods themselves!

When passing an `xarray.DataArray`, `pandas.DataFrame`, or `pandas.Series` through any
plotting method instead of a `numpy.ndarray`, the
axis tick labels, axis labels, subplot titles, and colorbar and legend labels are
automatically updated. This can be disabled by passing
``autoformat`` to the plotting method or to `~proplot.subplots.subplots`.
See :ref:`1d plotting` and :ref:`2d plotting` for details.

Cartopy and basemap integration
===============================

.. raw:: html

   <h3>Problem</h3>

There are two widely-used engines
for plotting geophysical data with matplotlib: `cartopy` and `~mpl_toolkits.basemap`.
Using cartopy tends to be quite verbose and involve lots of boilerplate code,
while basemap is outdated and requires you to use plotting commands on a separate `~mpl_toolkits.basemap.Basemap` object.

Also, `cartopy` and `~mpl_toolkits.basemap` plotting commands assume
*map projection coordinates* unless specified otherwise. For most of us, this
choice is very frustrating, since geophysical data are usually stored in
longitude-latitude or "Plate Carrée" coordinates.

.. raw:: html

   <h3>Solution</h3>

ProPlot integrates various `cartopy` and `~mpl_toolkits.basemap` features
into the `~proplot.axes.ProjAxes` `~proplot.axes.ProjAxes.format` method.
This lets you apply all kinds of geographic plot settings, like coastlines, continents, political boundaries, and meridian and parallel gridlines.
`~proplot.axes.ProjAxes` also
overrides various plotting methods:

* ``transform=ccrs.PlateCarree()`` is the new default for all `~proplot.axes.GeoAxes` plotting methods.
* ``latlon=True`` is the new default for all `~proplot.axes.BasemapAxes` plotting methods.
* ``globe=True`` can be passed to any 2D plotting command to enforce *global* coverage over the poles and across the longitude boundaries.

See :ref:`Geographic and polar plots` for details.
Note that active development on basemap will `halt after 2020 <https://matplotlib.org/basemap/users/intro.html#cartopy-new-management-and-eol-announcement>`__.
For now, cartopy is
`missing several features <https://matplotlib.org/basemap/api/basemap_api.html#module-mpl_toolkits.basemap>`__
offered by basemap -- namely, flexible meridian and parallel gridline labels,
drawing physical map scales, and convenience features for adding background images like
the "blue marble". But once these are added to cartopy, ProPlot may remove the `~mpl_toolkits.basemap` integration features.

..
  This is the right decision: Cartopy is integrated more closely with the matplotlib API
  and is more amenable to further development.

Colormaps and property cycles
=============================

.. raw:: html

   <h3>Problem</h3>

In matplotlib, colormaps are implemented with the `~matplotlib.colors.ListedColormap` and `~matplotlib.colors.LinearSegmentedColormap` classes.
They are hard to edit and hard to create from scratch.

..
   Colormap identification is also suboptimal, since the names are case-sensitive, and reversed versions of each colormap are not guaranteed to exist.

.. raw:: html

   <h3>Solution</h3>

In ProPlot, it is easy to manipulate colormaps and property cycles:

* The `~proplot.styletools.Colormap` constructor function can be used to slice and merge existing colormaps and/or generate brand new ones.
* The `~proplot.styletools.Cycle` constructor function can be used to make *color cycles* from *colormaps*! These cycles can be applied by passing the `cycle` keyword argument to plotting commands or changing the :rcraw:`cycle` setting. See :ref:`Color cycles` for details.
* The new `~proplot.styletools.ListedColormap` and `~proplot.styletools.LinearSegmentedColormap` classes include several new methods, e.g. `~proplot.styletools.LinearSegmentedColormap.save` and `~proplot.styletools.LinearSegmentedColormap.concatenate`, and have a much nicer REPL representation.
* The `~proplot.styletools.PerceptuallyUniformColormap` class is used to make :ref:`Perceptually uniform colormaps`. These have smooth, aesthetically pleasing color transitions represent your data *accurately*.

Importing ProPlot also makes colormap names *case-insensitive* and *reviersible* by appending ``'_r'`` tot he colormap name. This is powered by the `~proplot.styletools.CmapDict` dictionary, which replaces matplotlib's native colormap database.

Smarter colormap normalization
==============================
.. raw:: html

   <h3>Problem</h3>

In matplotlib, when ``extend='min'``, ``extend='max'``, or ``extend='neither'`` is passed to `~matplotlib.figure.Figure.colorbar` , the colormap colors reserved for "out-of-bounds" values are truncated.

The problem is that matplotlib discretizes colormaps by generating a low-resolution lookup table (see `~matplotlib.colors.LinearSegmentedColormap` for details).
This approach cannot be fine-tuned and creates an unnecessary copy of the colormap.

..
   and prevents you from using the resulting colormap for plots with different numbers of levels.

It is clear that the task discretizing colormap colors should be left to the **normalizer**, not the colormap itself. Matplotlib provides `~matplotlib.colors.BoundaryNorm` for this purpose, but it is seldom used and its features are limited.

.. raw:: html

   <h3>Solution</h3>

In ProPlot, all colormap visualizations are automatically discretized with the `~proplot.styletools.BinNorm` class. This reads the `extend` property passed to your plotting command and chooses colormap indices so that your colorbar levels *always* traverse the full range of colormap colors.

`~proplot.styletools.BinNorm` also applies arbitrary continuous normalizer requested by the user, e.g. `~matplotlib.colors.Normalize` or `~matplotlib.colors.LogNorm`, before discretization. Think of `~proplot.styletools.BinNorm` as a "meta-normalizer" -- other normalizers perform the continuous transformation step, while this performs the discretization step.

Bulk global settings
====================
.. raw:: html

   <h3>Problem</h3>

In matplotlib, there are several `~matplotlib.rcParams` that you often
want to set *all at once*, like the tick lengths and spine colors.
It is also often desirable to change these settings for *individual subplots*
or *individual blocks of code* rather than globally.

.. raw:: html

   <h3>Solution</h3>

In ProPlot, you can use the `~proplot.rctools.rc` object to
change lots of settings at once with convenient shorthands.
This is meant to replace matplotlib's `~matplotlib.rcParams`.
dictionary. Settings can be changed with ``plot.rc.key = value``, ``plot.rc[key] = value``,
``plot.rc.update(...)``, with the `~proplot.axes.Axes.format` method, or with the
`~proplot.rctools.rc_configurator.context` method.

For details, see :ref:`Configuring proplot`.
The most notable bulk settings are described below.

=============  =============================================  ===========================================================================================================================================================================
Key            Description                                    Children
=============  =============================================  ===========================================================================================================================================================================
``color``      The color for axes bounds, ticks, and labels.  ``axes.edgecolor``, ``geoaxes.edgecolor``, ``axes.labelcolor``, ``tick.labelcolor``, ``hatch.color``, ``xtick.color``, ``ytick.color``
``linewidth``  The width of axes bounds and ticks.            ``axes.linewidth``, ``geoaxes.linewidth``, ``hatch.linewidth``, ``xtick.major.width``, ``ytick.major.width``
``small``      Font size for "small" labels.                  ``font.size``, ``tick.labelsize``, ``xtick.labelsize``, ``ytick.labelsize``, ``axes.labelsize``, ``legend.fontsize``, ``geogrid.labelsize``
``large``      Font size for "large" labels.                  ``abc.size``, ``figure.titlesize``, ``axes.titlesize``, ``suptitle.size``, ``title.size``, ``leftlabel.size``, ``toplabel.size``, ``rightlabel.size``, ``bottomlabel.size``
``tickpad``    Padding between ticks and labels.              ``xtick.major.pad``, ``xtick.minor.pad``, ``ytick.major.pad``, ``ytick.minor.pad``
``tickdir``    Tick direction.                                ``xtick.direction``, ``ytick.direction``
``ticklen``    Tick length.                                   ``xtick.major.size``, ``ytick.major.size``, ``ytick.minor.size * tickratio``, ``xtick.minor.size * tickratio``
``tickratio``  Ratio between major and minor tick lengths.    ``xtick.major.size``, ``ytick.major.size``, ``ytick.minor.size * tickratio``, ``xtick.minor.size * tickratio``
``margin``     Margin width when limits not explicitly set.    ``axes.xmargin``, ``axes.ymargin``
=============  =============================================  ===========================================================================================================================================================================

Physical units engine
=====================
.. raw:: html

   <h3>Problem</h3>

Matplotlib requires users to use
*inches* for the figure size `figsize`. This must be confusing for users outside
of the U.S.

Matplotlib also uses figure-relative units for the margins
`left`, `right`, `bottom`, and `top`, and axes-relative units
for the column and row spacing `wspace` and `hspace`.
Relative units tend to require "tinkering" with meaningless numbers until you find the
right one... and then if your figure size changes, you have to adjust them again.

.. raw:: html

   <h3>Solution</h3>

ProPlot introduces the physical units engine `~proplot.utils.units`
for interpreting `figsize`, `width`, `height`, `axwidth`, `axheight`,
`left`, `right`, `top`, `bottom`, `wspace`, `hspace`, and arguments
in a few other places. Acceptable units include inches, centimeters,
millimeters, pixels, `points <https://en.wikipedia.org/wiki/Point_(typography)>`__,
`picas <https://en.wikipedia.org/wiki/Pica_(typography)>`__, `em-heights <https://en.wikipedia.org/wiki/Em_(typography)>`__, and `light years <https://en.wikipedia.org/wiki/Light-year>`__ (because why not?).
Em-heights are particularly useful, as labels already
present can be useful *rulers* for figuring out the amount
of space needed.

`~proplot.utils.units` is also used to convert settings
passed to `~proplot.rctools.rc` from arbitrary physical units
to *points* -- for example, :rcraw:`linewidth`, :rcraw:`ticklen`,
:rcraw:`axes.titlesize`, and :rcraw:`axes.titlepad`.
See :ref:`Configuring proplot` for details.

Working with fonts
==================
.. raw:: html

   <h3>Problem</h3>

In matplotlib, the default font is DejaVu Sans. In this developer's humble opinion, DejaVu Sans is fugly AF. It is also really tricky to work with custom fonts in matplotlib.

..
   This font is not very aesthetically pleasing.

.. raw:: html

   <h3>Solution</h3>

ProPlot comes packaged with several additional fonts. The new default font is Helvetica. Albeit somewhat overused, this is a tried and tested, aesthetically pleasing sans serif font.

ProPlot adds fonts to matplotlib by making use of a completely undocumented feature: the ``$TTFPATH`` environment variable (matplotlib adds ``.ttf`` and ``.otf`` font files from folders listed in ``$TTFPATH``). You can also use *your own* font files by dropping them in ``~/.proplot/fonts``.

..
   ...and much more!
   =================
   This page is not comprehensive -- it just
   illustrates how ProPlot addresses
   some of the stickiest matplotlib limitations
   that bug your average power user.
   See the User Guide for a more comprehensive overview.
