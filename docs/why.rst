============
Why ProPlot?
============

Matplotlib is an extremely powerful plotting package used
by academics, engineers, and data scientists far and wide. However, certain
plotting tasks can be cumbersome or repetitive for its heaviest users, i.e. users who...

* ...make very rich, complex figures.
* ...want to finely tune their figure annotations and aesthetics.
* ...create new figures nearly every day.

ProPlot's core mission is to provide a smoother plotting experience
for heavy matplotlib users.
We do this by expanding upon the object-oriented matplotlib API.
ProPlot makes changes that would be hard to justify or difficult
to incorporate into matplotlib itself, owing to design choices and backwards
compatibility considerations.
This page enumerates these changes and explains how they
address limitations of the matplotlib API.

..
   This page is not comprehensive --
   see the User Guide for a comprehensive overview
   with worked examples.

..
   To start using these new features, see
   see :ref:`Usage overview` and the User Guide.

Less typing, more plotting
==========================

.. rubric:: Problem

Matplotlib users often need to change lots of plot settings all at once. With the default API, this requires calling a series of one-liner setter methods.

This workflow is quite verbose -- it tends to require "boilerplate code" that gets copied and pasted a hundred times. It can also be confusing -- it is often unclear whether properties are applied from an `~matplotlib.axes.Axes` setter (e.g. `~matplotlib.axes.Axes.set_title`, `~matplotlib.axes.Axes.set_xlabel` and `~matplotlib.axes.Axes.set_xticks`), an `~matplotlib.axis.XAxis` or `~matplotlib.axis.YAxis` setter (e.g. `~matplotlib.axis.Axis.set_major_locator` and `~matplotlib.axis.Axis.set_major_formatter`), a `~matplotlib.spines.Spine` setter (e.g. `~matplotlib.spines.Spine.set_bounds`), a miscellaneous "bulk" setter (e.g. `~matplotlib.axes.Axes.tick_params`), or whether they require tinkering with several different objects. Also, one often needs to *loop through* lists of subplots to apply identical settings to each subplot.

..
   This is perhaps one reason why many users prefer the `~matplotlib.pyplot` API to the object-oriented API (see :ref:`Using ProPlot`).

.. rubric:: Solution

ProPlot introduces the `~proplot.axes.Axes.format` method for changing arbitrary settings *in bulk*. Think of this as an expanded and thoroughly documented version of the
`~matplotlib.artist.Artist` `~matplotlib.artist.Artist.update` method.
`~proplot.axes.Axes.format` can
also be used to update :ref:`Bulk global settings` and various :ref:`other rc settings <Changing rc settings>` for a particular subplot, and to concisely work with verbose classes using the :ref:`Class constructor functions`. Further, :ref:`The subplot container class` can be used to invoke `~proplot.axes.Axes.format` on several subplots at once.

Together, these features significantly reduce
the amount of code needed to create highly customized figures.
As an example, it is trivial to see that

.. code-block:: python

   import proplot as plot
   f, axs = plot.subplots(ncols=2)
   axs.format(linewidth=1, color='gray')
   axs.format(xticks=20, xtickminor=True, xlabel='x axis', ylabel='y axis')

...is much more succinct than

.. code-block:: python

   import matplotlib.pyplot as plt
   import matplotlib.ticker as mticker
   from matplotlib import rcParams
   rcParams['axes.linewidth'] = 1
   rcParams['axes.color'] = 'gray'
   fig, axs = plt.subplots(ncols=2)
   for ax in axs:
      ax.xaxis.set_major_locator(mticker.MultipleLocator(10))
      ax.tick_params(width=1, color='gray', labelcolor='gray')
      ax.tick_params(axis='x', which='minor', bottom=True)
      ax.set_xlabel('x axis', color='gray')
      ax.set_ylabel('y axis', color='gray')
   plt.style.use('default')  # restore


Class constructor functions
===========================
.. rubric:: Problem

Matplotlib and cartopy introduce a bunch of classes with verbose names like `~matplotlib.ticker.MultipleLocator`, `~matplotlib.ticker.FormatStrFormatter`, and
`~cartopy.crs.LambertAzimuthalEqualArea`. Since plotting code has a half life of about 30 seconds, typing out all of these extra class names and import statements can be a *major* drag.

Parts of the matplotlib API were actually designed with this in mind.
`Backend classes <https://matplotlib.org/faq/usage_faq.html#what-is-a-backend>`__,
`native axes projections <https://matplotlib.org/3.1.1/api/projections_api.html>`__,
`axis scales <https://matplotlib.org/3.1.0/gallery/scales/scales.html>`__,
`box styles <https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.patches.FancyBboxPatch.html?highlight=boxstyle>`__, `arrow styles <https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.patches.FancyArrowPatch.html?highlight=arrowstyle>`__, and
`arc styles <https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.patches.ConnectionStyle.html?highlight=class%20name%20attrs>`__
are referenced with "registered" string names,
as are `basemap projection types <https://matplotlib.org/basemap/users/mapsetup.html>`__.
So, why not "register" everything else?

.. rubric:: Solution

In ProPlot, tick locators, tick formatters, axis scales, cartopy projections, colormaps, and property cyclers are all "registered". This is done by creating several **constructor functions** and passing various keyword argument *through* the constructor functions.
This may seem "unpythonic" but it is absolutely invaluable when writing
plotting code.

Each constructor function accepts various *other* input types for your convenience. For
example, scalar numbers passed to `~proplot.axistools.Locator` returns
a `~matplotlib.ticker.MultipleLocator` instance, lists of strings passed
to `~proplot.axistools.Formatter` returns a `~matplotlib.ticker.FixedFormatter` instance, and `~proplot.styletools.Colormap` and `~proplot.styletools.Cycle` accept colormap names, individual colors, and lists of colors. When a *class instance* is passed to the relevant constructor function, it is simply returned. See :ref:`X and Y axis settings`, :ref:`Colormaps`, and :ref:`Color cycles` for details.

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

.. rubric:: Problem

Matplotlib plots tend to require lots of "tweaking" when you have more than one subplot in the figure. This is partly because you must specify the physical dimensions of the figure, while the dimensions of the *individual subplots* are more important:

#. The subplot aspect ratio is usually more relevant than the figure aspect ratio, e.g. for map projections.
#. The subplot width and height control the evident thickness of text and other content plotted inside the axes.

Matplotlib has a `tight layout <https://matplotlib.org/tutorials/intermediate/tight_layout_guide.html>`__ algorithm to keep you from having to "tweak" the spacing, but the algorithm cannot apply different amounts of spacing between different subplot row and column boundaries. This limitation often results in unnecessary whitespace, and can be a major problem when you want to put e.g. a legend on the outside of a subplot.

.. rubric:: Solution

In ProPlot, you can specify the physical dimensions of a *reference subplot* instead of the figure by passing `axwidth`, `axheight`, and/or `aspect` to `~proplot.subplots.Figure`. The default behavior is ``aspect=1`` and ``axwidth=2`` (inches). If the `aspect ratio mode <https://matplotlib.org/2.0.2/examples/pylab_examples/equal_aspect_ratio.html>`__ for the reference subplot is set to ``'equal'``, as with :ref:`Geographic and polar plots` and `~matplotlib.axes.Axes.imshow` plots, the *imposed* aspect ratio will be used instead.
Figure dimensions are constrained as follows:

* When `axwidth` or `axheight` are specified, the figure width and height are calculated automatically.
* When `width` is specified, the figure height is calculated automatically.
* When `height` is specified, the figure width is calculated automatically.
* When `width` *and* `height` or `figsize` is specified, the figure dimensions are fixed.

..
   Several matplotlib backends require figure dimensions to be fixed. When `~proplot.subplots.Figure.draw` changes the figure dimensions, this can "surprise" the backend and cause unexpected behavior. ProPlot fixes this issue for the static inline backend and the Qt popup backend. However, this issue is unfixable the "notebook" inline backend, the "macosx" popup backend, and possibly other untested backends.

ProPlot also uses a custom tight layout algorithm that automatically determines the `left`, `right`, `bottom`, `top`, `wspace`, and `hspace` `~matplotlib.gridspec.GridSpec` parameters. This algorithm is simpler because:

* The new `~proplot.subplots.GridSpec` class permits variable spacing between rows and columns. It turns out this is critical for putting :ref:`Colorbars and legends` on the outside of subplots.
* Figures are restricted to have only *one* `~proplot.subplots.GridSpec` per figure. This is done by requiring users to draw all of their subplots at once with `~proplot.subplots.subplots` (see :pr:`50`).

See :ref:`Automatic layout` for details.

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

Eliminating redundancies
========================

.. rubric:: Problem

For many of us, figures with just one subplot are a rarity. We tend to need multiple
subplots for comparing different datasets and illustrating complex concepts.
Unfortunately, it is easy to end up with *redundant* figure elements
when drawing multiple subplots; namely:

* Repeated axis tick labels.
* Repeated axis labels.
* Repeated colorbars.
* Repeated legends.

These sorts of redundancies are extremely common even in publications, where
they waste valuable page space. They arise because this is the path of least
resistance for the default API -- removing redundancies
tends to require a fair amount of extra work.

.. rubric:: Solution

ProPlot seeks to eliminate redundant elements
to help you make clear, concise figures.
We tackle this issue using
:ref:`Shared and spanning labels` and :ref:`Figure colorbars and legends`.

* By default, axis tick labels and axis labels are *shared* between subplots in the same row or column. This is controlled by the `sharex`, `sharey`, `spanx`, and `spany` `~proplot.subplots.subplots` keyword args.
* The new `~proplot.subplots.Figure` `~proplot.subplots.Figure.colorbar` and `~proplot.subplots.Figure.legend` methods make it easy to draw colorbars and legends intended to reference more than one subplot. For details, see the next section.

Outer colorbars and legends
===========================

.. rubric:: Problem

In matplotlib, it is difficult to draw `~matplotlib.figure.Figure.colorbar`\ s and
`~matplotlib.axes.Axes.legend`\ s intended to reference more than one subplot or
along the outside of subplots:

* To draw legends outside of subplots, you usually need to position the legend manually and adjust various `~matplotlib.gridspec.GridSpec` spacing properties to make *room* for the legend.
* To make colorbars that span multiple subplots, you have to supply `~matplotlib.figure.Figure.colorbar` with a `cax` you drew yourself. This requires so much tinkering that most users just add identical colorbars to every single subplot!

Furthermore, drawing colorbars with ``fig.colorbar(..., ax=ax)`` tends to mess up subplot aspect ratios since the space allocated for the colorbar is "stolen" from the parent axes.

..
   And since colorbar widths are specified in *axes relative* coordinates, they often look "too skinny" or "too fat" after the first draw.


..
   The matplotlib example for `~matplotlib.figure.Figure` legends is `not pretty <https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/figlegend_demo.html>`__.

..
   Drawing colorbars and legends is pretty clumsy in matplotlib -- especially when trying to draw them outside of the figure. They can be too narrow, too wide, and mess up your subplot aspect ratios.

.. rubric:: Solution

ProPlot introduces a brand new framework for drawing :ref:`Axes colorbars and legends`
(colorbars and legends inside or along the outside edge of a subplot)
and :ref:`Figure colorbars and legends`
(colorbars and legends sapnning contiguous subplots along the edge of the figure):

* Passing an "outer" location to `~proplot.axes.Axes` `~proplot.axes.Axes.colorbar` or `~proplot.axes.Axes` `~proplot.axes.Axes.legend` (e.g. ``loc='l'`` or ``loc='left'``) draws the colorbar or legend along the outside of the axes. Passing an "inner" location (e.g. ``loc='ur'`` or ``loc='upper right'``) draws an *inset* colorbar or legend. And yes, that's right, you can now draw inset colorbars!
* To draw a colorbar or legend along the edge of the figure, use `~proplot.subplots.Figure` `~proplot.subplots.Figure.colorbar` and `~proplot.subplots.Figure.legend`. The `col`, `row`, and `span` keyword args control which `~matplotlib.gridspec.GridSpec` rows and columns are spanned by the colorbar or legend.
* Since `~proplot.subplots.GridSpec` permits variable spacing between subplot rows and columns, "outer" colorbars and legends do not mess up subplot spacing or add extra whitespace. This is critical e.g. if you have a colorbar between columns 1 and 2 but nothing between columns 2 and 3.
* `~proplot.subplots.Figure` and `~proplot.axes.Axes` colorbar widths are specified in *physical* units rather than relative units. This makes colorbar thickness independent of subplot size and easier to get just right.

There are also several :ref:`New colorbar features` and :ref:`New legend features`.

The subplot container class
===========================

..
   The `~matplotlib.pyplot.subplots` command is useful for generating a scaffolding of * axes all at once. This is generally faster than successive `~matplotlib.subplots.Figure.add_subplot` commands.

.. rubric:: Problem

In matplotlib, `~matplotlib.pyplot.subplots` returns a 2d `~numpy.ndarray` for figures with more than one column and row, a 1d `~numpy.ndarray` for single-row or single-column figures, or just an `~matplotlib.axes.Axes` instance for single-subplot figures.

.. rubric:: Solution

In ProPlot, `~proplot.subplots.subplots` returns a `~proplot.subplots.subplot_grid`
filled with `~proplot.axes.Axes` instances.
This container lets you call arbitrary methods on arbitrary subplots all at once, which can be useful when you want to style your subplots identically (e.g. ``axs.format(tickminor=False)``).
The `~proplot.subplots.subplot_grid` class also
unifies the behavior of the three possible `matplotlib.pyplot.subplots` return values:

* `~proplot.subplots.subplot_grid` permits 2d indexing, e.g. ``axs[1,0]``. Since `~proplot.subplots.subplots` can generate figures with arbitrarily complex subplot geometry, this 2d indexing is useful only when the arrangement happens to be a clean 2d matrix.
* `~proplot.subplots.subplot_grid` also permits 1d indexing, e.g. ``axs[0]``, since it is a `list` subclass. The default order can be switched from row-major to column-major by passing ``order='F'`` to `~proplot.subplots.subplots`.
* When it is singleton, `~proplot.subplots.subplot_grid` behaves like a scalar. So when you make a single axes with ``f, axs = plot.subplots()``, ``axs[0].method(...)`` is equivalent to ``axs.method(...)``.

See :ref:`Subplot grids` for details.

..
   This goes with ProPlot's theme of preserving the object-oriented spirit, but making things easier for users.

New and improved plotting methods
=================================

.. rubric:: Problem

Certain plotting tasks are quite difficult to accomplish
with the default matplotlib API. The `seaborn`, `xarray`, and `pandas`
packages offer improvements, but it would be nice
to have this functionality build right into matplotlib.
There is also room for improvement of the native matplotlib plotting methods
that none of these packages address.

..
   Matplotlib also has some finicky plotting issues
   that normally requires
..
   For example, when you pass coordinate *centers* to `~matplotlib.axes.Axes.pcolor` and `~matplotlib.axes.Axes.pcolormesh`, they are interpreted as *edges* and the last column and row of your data matrix is ignored. Also, to add labels to `~matplotlib.axes.Axes.contour` and `~matplotlib.axes.Axes.contourf`, you need to call a dedicated `~matplotlib.axes.Axes.clabel` method instead of just using a keyword argument.


.. rubric:: Solution


ProPlot adds various
`seaborn`, `xarray`, and `pandas` features
to the `~proplot.axes.Axes` plotting methods
along with several *brand new* features designed to
make your life easier.

* The new `~proplot.axes.Axes.area` and `~proplot.axes.Axes.areax` methods call `~matplotlib.axes.Axes.fill_between` and `~matplotlib.axes.Axes.fill_betweenx`. These methods now accept 2D arrays and *stack* or *overlay* successive columns, and a `negpos` keyword argument that can be used to assign separate colors for negative and positive data.
* The new `~proplot.axes.Axes.parametric` method draws *parametric* line plots, where the parametric coordinate is denoted with a colorbar rather than text annotations. This is much cleaner and more aesthetically pleasing than the conventional approach.
* The new `~proplot.axes.Axes.heatmap` method invokes `~matplotlib.axes.Axes.pcolormesh` and draws ticks at the center of each box. This is more convenient for things like covariance matrices.
* The `~matplotlib.axes.Axes.bar` and `~matplotlib.axes.Axes.barh` methods accept 2D arrays and *stack* or *group* successive columns. Just like `~matplotlib.axes.Axes.fill_between` and `~matplotlib.axes.Axes.fill_betweenx`, you will be able to use different colors for positive/negative bars.
* All :ref:`1d plotting` can be used to draw :ref:`On-the-fly error bars` using the `means`, `medians`, `boxdata`, and `bardata` keyword arguments. You no longer have to work with `~matplotlib.axes.Axes.add_errobar` method directly.
* All :ref:`1d plotting` methods accept a `cycle` keyword argument interpreted by `~proplot.styletools.Cycle` and optional `legend` and `colorbar` keyword arguments for populating legends and colorbars at the specified location with the result of the plotting command. See :ref:`Color cycles` and :ref:`Colorbars and legends`.
* All :ref:`2d plotting` methods accept a `cmap` keyword argument interpreted by `~proplot.styletools.Colormap`, a `norm` keyword argument interpreted by `~proplot.styletools.Norm`, and an optional `colorbar` keyword argument for drawing on-the-fly colorbars with the resulting mappable. See :ref:`Colormaps` and :ref:`Colorbars and legends`.
* All :ref:`2d plotting` methods accept a `labels` keyword argument. This is used to draw contour labels or grid box labels on heatmap plots. Labels are colored black or white according to the luminance of the underlying filled contour or grid box color. See :ref:`2d plotting` for details.
* ProPlot fixes the irritating `white-lines-between-filled-contours <https://stackoverflow.com/q/8263769/4970632>`__, `white-lines-between-pcolor-patches <https://stackoverflow.com/q/27092991/4970632>`__, and `white-lines-between-colorbar-patches <https://stackoverflow.com/q/15003353/4970632>`__ vector graphic issues.
* Matplotlib requires coordinate *centers* for contour plots and *edges* for pcolor plots. If you pass *centers* to pcolor, matplotlib treats them as *edges* and silently trims one row/column of your data. Most people don't realize this! ProPlot changes this behavior: If edges are passed to `~matplotlib.axes.Axes.contour` or `~matplotlib.axes.Axes.contourf`, centers are *calculated* from the edges; if centers are passed to `~matplotlib.axes.Axes.pcolor` or `~matplotlib.axes.Axes.pcolormesh`, edges are *estimated* from the centers.

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

Xarray and pandas integration
=============================

.. rubric:: Problem

When you pass the array-like `xarray.DataArray`, `pandas.DataFrame`, and `pandas.Series` containers to matplotlib plotting commands, the metadata is ignored. To create plots that are automatically labeled with this metadata, you must use
the dedicated `xarray.DataArray.plot`, `pandas.DataFrame.plot`, and `pandas.Series.plot`
tools instead.

This approach is fine for quick plots, but not ideal for complex ones.
It requires learning a different syntax from matplotlib, and tends to encourage using the `~matplotlib.pyplot` API rather than the object-oriented API.
These tools also introduce features that would be useful additions to matplotlib
in their *own* right, without requiring special data containers and
an entirely separate API.

.. rubric:: Solution

ProPlot *reproduces* most of the `xarray.DataArray.plot`, `pandas.DataFrame.plot`, and `pandas.Series.plot` features on the `~proplot.axes.Axes` plotting methods themselves.
Passing an `~xarray.DataArray`, `~pandas.DataFrame`, or `~pandas.Series` through
any plotting method automatically updates the
axis tick labels, axis labels, subplot titles, and colorbar and legend labels
from the metadata.  This can be disabled by passing
``autoformat=False`` to the plotting method or to `~proplot.subplots.subplots`.

Also, as described in :ref:`New and improved plotting methods`, ProPlot implements certain
features like grouped bar plots, layered area plots, heatmap plots,
and on-the-fly colorbars and legends from the
`xarray` and `pandas` APIs directly on the `~proplot.axes.Axes` class.

Cartopy and basemap integration
===============================

.. rubric:: Problem

There are two widely-used engines
for plotting geophysical data with matplotlib: `cartopy` and `~mpl_toolkits.basemap`.
Using cartopy tends to be verbose and involve boilerplate code,
while using basemap requires you to use plotting commands on a
separate `~mpl_toolkits.basemap.Basemap` object rather than an axes object.

Also, `cartopy` and `~mpl_toolkits.basemap` plotting commands assume
*map projection coordinates* unless specified otherwise. For most of us, this
choice is very frustrating, since geophysical data are usually stored in
longitude-latitude or "Plate Carrée" coordinates.

.. rubric:: Solution

ProPlot integrates various `cartopy` and `~mpl_toolkits.basemap` features
into the `~proplot.axes.ProjAxes` `~proplot.axes.ProjAxes.format` method.
This lets you apply all kinds of geographic plot settings, like coastlines, continents, political boundaries, and meridian and parallel gridlines.
`~proplot.axes.ProjAxes` also
overrides various plotting methods:

* The new default for all `~proplot.axes.GeoAxes` plotting methods is ``transform=ccrs.PlateCarree()``.
* The new default for all `~proplot.axes.BasemapAxes` plotting methods is ``latlon=True``.
* *Global* coverage over the poles and across the matrix longitude boundaries can be enforced by passing ``globe=True`` to any 2d plotting command, e.g. `~matplotlib.axes.Axes.pcolormesh` and `~matplotlib.axes.Axes.contourf`.

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

.. rubric:: Problem

In matplotlib, colormaps are implemented with the `~matplotlib.colors.ListedColormap` and `~matplotlib.colors.LinearSegmentedColormap` classes.
They are hard to edit and hard to create from scratch.

..
   Colormap identification is also suboptimal, since the names are case-sensitive, and reversed versions of each colormap are not guaranteed to exist.

.. rubric:: Solution

In ProPlot, it is easy to manipulate colormaps and property cycles:

* The `~proplot.styletools.Colormap` constructor function can be used to slice and merge existing colormaps and/or generate brand new ones.
* The `~proplot.styletools.Cycle` constructor function can be used to make *color cycles* from *colormaps*! Color cycles can be applied to plots in a variety of ways; see :ref:`Color cycles` for details.
* The new `~proplot.styletools.ListedColormap` and `~proplot.styletools.LinearSegmentedColormap` classes include several convenient methods and have a much nicer REPL string representation.
* The `~proplot.styletools.PerceptuallyUniformColormap` class is used to make :ref:`Perceptually uniform colormaps`. These have smooth, aesthetically pleasing color transitions represent your data *accurately*.

Importing ProPlot also makes all colormap names *case-insensitive*, and colormaps can be *reversed* or *cyclically shifted* by 180 degrees simply by appending ``'_r'`` or ``'_shifted'`` to the colormap name. This is powered by the `~proplot.styletools.CmapDict` dictionary, which replaces matplotlib's native colormap database.

Smarter colormap normalization
==============================
.. rubric:: Problem

In matplotlib, when ``extend='min'``, ``extend='max'``, or ``extend='neither'`` is passed to `~matplotlib.figure.Figure.colorbar` , the colormap colors reserved for "out-of-bounds" values are truncated. This can be irritating for plots with very few colormap levels, which are often more desirable (see :ref:`Discrete colormap levels`).

The problem is that matplotlib "discretizes" colormaps by generating low-resolution lookup tables (see `~matplotlib.colors.LinearSegmentedColormap`).
While straightforward,
this approach has limitations and results in unnecessary
plot-specific copies of the colormap.
Ideally, the task of discretizing colormap colors should be left to the *normalizer*; matplotlib provides `~matplotlib.colors.BoundaryNorm` for this purpose, but it is seldom used and its features are limited.

.. rubric:: Solution

In ProPlot, all colormaps retain a high-resolution lookup table and the `~proplot.styletools.BinNorm` class is applied to every plot. `~proplot.styletools.BinNorm` restricts your plot to a *subset* of lookup table colors matching the number of requested levels. It chooses indices such that the colorbar levels *always* traverse the full range of colors, no matter the `extend` setting, and makes sure the end colors on *cyclic* colormaps are distinct.

Also, before discretization, `~proplot.styletools.BinNorm` passes values through the *continuous* normalizer requested by the user with the `norm` keyword argument (e.g. `~matplotlib.colors.LogNorm` or `~proplot.styletools.MidpointNorm`). You can thus think of `~proplot.styletools.BinNorm` as a "meta-normalizer": `~proplot.styletools.BinNorm` simply discretizes the result of any arbitrary continuous transformation.

Bulk global settings
====================
.. rubric:: Problem

In matplotlib, there are several `~matplotlib.rcParams` that you often
want to set *all at once*, like the tick lengths and spine colors.
It is also often desirable to change these settings for *individual subplots*
or *individual blocks of code* rather than globally.

.. rubric:: Solution

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
.. rubric:: Problem

Matplotlib requires users to use
inches for the figure size `figsize`. This may be confusing for users outside
of the U.S.

Matplotlib also uses figure-relative units for the margins
`left`, `right`, `bottom`, and `top`, and axes-relative units
for the column and row spacing `wspace` and `hspace`.
Relative units tend to require "tinkering" with numbers until you find the
right one. And since they are *relative*, if you decide to change your
figure size or add a subplot, they will have to be readjusted.

.. rubric:: Solution

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


The .proplot folder
===================
.. rubric:: Problem

In matplotlib, it can be difficult to design your
own colormaps and color cycles, and there is no builtin
way to *save* them for future use. It is also
difficult to get matplotlib to use custom ``.ttc``, ``.ttf``,
and ``.otf`` font files, which may be desirable when you are
working on Linux servers with limited font selections.


.. rubric:: Solution

ProPlot automatically adds colormaps, color cycles, and font files
saved in the ``.proplot/cmaps``,  ``.proplot/cycles``, and ``.proplot/fonts``
folders in your home directory.
You can save colormaps and color
cycles to these folders simply by passing ``save=True`` to
`~proplot.styletools.Colormap` and `~proplot.styletools.Cycle`.
To *manually* load from these folders, e.g. if you have added
files to these folders but you do not want to restart your
ipython session, simply call
`~proplot.styletools.register_cmaps`,
`~proplot.styletools.register_cycles`, and
`~proplot.styletools.register_fonts`.

..
   As mentioned above,
   ProPlot introduces the `~proplot.styletools.Colormap` and  `~proplot.styletools.Cycle`.
   functions for designing your own colormaps and color cycles.

..
   ...and much more!
   =================
   This page is not comprehensive -- it just
   illustrates how ProPlot addresses
   some of the stickiest matplotlib limitations
   that bug your average power user.
   See the User Guide for a more comprehensive overview.
