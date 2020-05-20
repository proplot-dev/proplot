============
Why ProPlot?
============

Matplotlib is an extremely powerful plotting package used
by academics, engineers, and data scientists far and wide. However,
the default matplotlib API can be cumbersome or repetitive for
users who...

* ...make very complex figures with multiple subplots.
* ...want to finely tune their figure annotations and aesthetics.
* ...need to make new figures nearly every day.

ProPlot's core mission is to provide a smoother plotting experience
for matplotlib's heaviest users. We accomplish this by *expanding upon*
the object-oriented matplotlib API. ProPlot makes changes that would be
hard to justify or difficult to incorporate into matplotlib itself, owing
to differing design choices and backwards compatibility considerations.

This page enumerates these changes and explains how they
address the limitations of the matplotlib API.

..
   This page is not comprehensive --
   see the User Guide for a comprehensive overview
   with worked examples.

..
   To start using these new features, see
   see :ref:`Usage overview` and the User Guide.

.. _why_less_typing:

Less typing, more plotting
==========================

.. rubric:: Problem

Matplotlib users often need to change lots of plot settings all at once. With
the default API, this requires calling a series of one-liner setter methods.

This workflow is quite verbose -- it tends to require "boilerplate code" that
gets copied and pasted a hundred times. It can also be confusing -- it is
often unclear whether properties are applied from an `~matplotlib.axes.Axes`
setter (e.g. `~matplotlib.axes.Axes.set_xlabel` and
`~matplotlib.axes.Axes.set_xticks`), an `~matplotlib.axis.XAxis` or
`~matplotlib.axis.YAxis` setter (e.g.
`~matplotlib.axis.Axis.set_major_locator` and
`~matplotlib.axis.Axis.set_major_formatter`), a `~matplotlib.spines.Spine`
setter (e.g. `~matplotlib.spines.Spine.set_bounds`), or a "bulk" property
setter (e.g. `~matplotlib.axes.Axes.tick_params`), or whether one must dig
into the figure architecture and apply settings to several different objects.
While this is in the spirit of object-oriented design, it seems like there
should be a more unified, straightforward way to change settings for
day-to-day matplotlib usage.

..
   This is perhaps one reason why many users prefer the `~matplotlib.pyplot`
   API to the object-oriented API (see :ref:`Using ProPlot`).

.. rubric:: Solution

ProPlot introduces the `~proplot.axes.Axes.format` method for changing
arbitrary settings all at once. Think of this as an expanded and thoroughly
documented version of the `matplotlib.artist.Artist.update` method.  It can
also be used to update so-called :ref:`quick settings <why_quick>` and
various :ref:`other rc settings <Changing rc settings>` for a particular
subplot, and to concisely work with verbose classes using the
:ref:`constructor functions <why_constructor>`. Further, the :ref:`subplots
container class <why_container>` can be used to invoke
`~proplot.axes.Axes.format` on several subplots at once.

Together, these features significantly reduce
the amount of code needed to create highly customized figures.
As an example, it is trivial to see that

.. code-block:: python

   import proplot as plot
   fig, axs = plot.subplots(ncols=2)
   axs.format(linewidth=1, color='gray')
   axs.format(xlim=(0, 100), xticks=10, xtickminor=True, xlabel='x axis', ylabel='y axis')

...is much more succinct than

.. code-block:: python

   import matplotlib.pyplot as plt
   import matplotlib.ticker as mticker
   import matplotlib as mpl
   with mpl.rc_context(rc={'axes.linewidth': 1, 'axes.color': 'gray'}):
       fig, axs = plt.subplots(ncols=2, sharey=True)
       axs[0].set_ylabel('y axis', color='gray')
       for ax in axs:
           ax.set_xlim(0, 100)
           ax.xaxis.set_major_locator(mticker.MultipleLocator(10))
           ax.tick_params(width=1, color='gray', labelcolor='gray')
           ax.tick_params(axis='x', which='minor', bottom=True)
           ax.set_xlabel('x axis', color='gray')


.. _why_constructor:

Class constructor functions
===========================

.. rubric:: Problem

Matplotlib and cartopy introduce a bunch of classes with verbose names like `~matplotlib.ticker.MultipleLocator`, `~matplotlib.ticker.FormatStrFormatter`, and `~cartopy.crs.LambertAzimuthalEqualArea`. Since plotting code has a half life of about 30 seconds, typing out these extra class names and import statements can be a *major* drag.

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

In ProPlot, tick locators, tick formatters, axis scales, cartopy projections, colormaps, and property cyclers are all "registered". ProPlot does this by introducing several *constructor functions* and passing various keyword argument through the constructor functions. This may seem "unpythonic" but it is absolutely invaluable when writing plotting code.

The constructor functions also accept other input types for your convenience.
For example, scalar numbers passed to `~proplot.constructor.Locator` returns
a `~matplotlib.ticker.MultipleLocator` instance, lists of strings passed to
`~proplot.constructor.Formatter` returns a
`~matplotlib.ticker.FixedFormatter` instance, and
`~proplot.constructor.Colormap` and `~proplot.constructor.Cycle` accept
colormap names, individual colors, and lists of colors. Passing the relevant
class instance to a constructor function simply returns the instance.

See the user guide sections on :ref:`x and y axis settings <ug_xy_axis>`,
:ref:`colormaps <ug_cmaps>`, and :ref:`color cycles <ug_cycles>` for
details. The below table lists the constructor functions and the keyword
arguments that use them.

================================  ============================================================  =============================================================  =================================================================================================================================================================================================
Function                          Return type                                                   Used by                                                        Keyword argument(s)
================================  ============================================================  =============================================================  =================================================================================================================================================================================================
`~proplot.constructor.Locator`    `~matplotlib.ticker.Locator`                                  `~proplot.axes.Axes.format` and `~proplot.axes.Axes.colorbar`  ``locator=``, ``xlocator=``, ``ylocator=``, ``minorlocator=``, ``xminorlocator=``, ``yminorlocator=``, ``ticks=``, ``xticks=``, ``yticks=``, ``minorticks=``, ``xminorticks=``, ``yminorticks=``
`~proplot.constructor.Formatter`  `~matplotlib.ticker.Formatter`                                `~proplot.axes.Axes.format` and `~proplot.axes.Axes.colorbar`  ``formatter=``, ``xformatter=``, ``yformatter=``, ``ticklabels=``, ``xticklabels=``, ``yticklabels=``
`~proplot.constructor.Scale`      `~matplotlib.scale.ScaleBase`                                 `~proplot.axes.Axes.format`                                    ``xscale=``, ``yscale=``
`~proplot.constructor.Cycle`      `~cycler.Cycler`                                              :ref:`1D plotting methods <ug_1dplots>`                        ``cycle=``
`~proplot.constructor.Colormap`   `~matplotlib.colors.Colormap`                                 :ref:`2D plotting methods <ug_2dplots>`                        ``cmap=``
`~proplot.constructor.Norm`       `~matplotlib.colors.Normalize`                                :ref:`2D plotting methods <ug_2dplots>`                        ``norm=``
`~proplot.constructor.Proj`       `~cartopy.crs.Projection` or `~mpl_toolkits.basemap.Basemap`  `~proplot.ui.subplots`                                         ``proj=``
================================  ============================================================  =============================================================  =================================================================================================================================================================================================

Note that `~matplotlib.axes.Axes.set_xscale` and `~matplotlib.axes.Axes.set_yscale` now accept instances of `~matplotlib.scale.ScaleBase` thanks to a monkey patch
applied by ProPlot.

.. _why_spacing:

Automatic dimensions and spacing
================================

.. rubric:: Problem

Matplotlib plots tend to require lots of "tweaking" when you have more than
one subplot in the figure. This is partly because you must specify the
physical dimensions of the figure, despite the fact that...

#. ...the subplot aspect ratio is generally more relevant than the figure
   aspect ratio. An aspect ratio of ``1`` is desirable for most plots, and
   the aspect ratio must be held fixed for
   :ref:`geographic and polar <ug_proj>` projections and most
   `~matplotlib.axes.Axes.imshow` plots.
#. ...the physical width and height of the subplot controls the *evident*
   thickness of text, lines, and other content plotted inside the subplot.
   By comparison, the effect of the figure size on the *evident* thickness
   depends on the number of subplots in the figure.

Also, while matplotlib has a `tight layout
<https://matplotlib.org/tutorials/intermediate/tight_layout_guide.html>`__
algorithm to keep you from having to "tweak" the *spacing*, the algorithm
cannot apply different amounts of spacing between different subplot row and
column boundaries. This limitation often results in unnecessary whitespace,
and can be a major problem when you want to put e.g. a legend on the outside
of a subplot.

.. rubric:: Solution

In ProPlot, you can specify the physical dimensions of a *reference subplot*
instead of the figure by passing `axwidth`, `axheight`, and/or `aspect` to
`~proplot.figure.Figure`. The default behavior is ``aspect=1`` and
``axwidth=2`` (inches). If the `aspect ratio mode
<https://matplotlib.org/2.0.2/examples/pylab_examples/equal_aspect_ratio.html>`__
for the reference subplot is set to ``'equal'``, as with
:ref:`geographic and polar <ug_proj>` plots and
`~matplotlib.axes.Axes.imshow` plots, the *imposed* aspect ratio will be used
instead.  Figure dimensions are constrained as follows:

* When `axwidth` or `axheight` are specified, the figure dimensions are
  determined automatically.
* When `width` is specified, the figure height is determined automatically.
* When `height` is specified, the figure width is determined automatically.
* When `width` *and* `height` or `figsize` is specified, the figure
  dimensions are fixed.

..
   Several matplotlib backends require figure dimensions to be fixed. When `~proplot.figure.Figure.draw` changes the figure dimensions, this can "surprise" the backend and cause unexpected behavior. ProPlot fixes this issue for the static inline backend and the Qt popup backend. However, this issue is unfixable the "notebook" inline backend, the "macosx" popup backend, and possibly other untested backends.

ProPlot also uses a custom tight layout algorithm that automatically
determines the `left`, `right`, `bottom`, `top`, `wspace`, and `hspace`
`~matplotlib.gridspec.GridSpec` parameters. This algorithm is simpler
and more accurate for the following reasons:

* The new `~proplot.gridspec.GridSpec` class permits variable spacing between
  rows and columns. This is critical for putting :ref:`colorbars
  and legends <ug_cbars_legends>` outside of subplots without "stealing
  space" from the parent subplot.
* Figures are restricted to have only *one* `~proplot.gridspec.GridSpec` per
  figure. This is done by requiring users to draw all of their subplots at
  once with `~proplot.ui.subplots` (although in a :pr:`future version <50>`,
  there will be a `~proplot.ui.figure` function that allows users to add
  subplots one-by-one while retaining the gridspec restriction).

See the :ref:`user guide <ug_subplots>` for details.

..
   #. The `~proplot.gridspec.GridSpec` spacing parameters are specified in
   physical units instead of figure-relative units.

..
   The `~matplotlib.gridspec.GridSpec` class is useful for creating figures
   with complex subplot geometry.

..
   Users want to control axes positions with gridspecs.

..
   * Matplotlib permits arbitrarily many `~matplotlib.gridspec.GridSpec`\ s
   per figure. This greatly complicates the tight layout algorithm for
   little evident gain.

..
   ProPlot introduces a marginal limitation (see discussion in :pr:`50`) but
   *considerably* simplifies the tight layout algorithm.


.. _why_redundancies:

Eliminating redundancies
========================

.. rubric:: Problem

For many of us, figures with just one subplot are a rarity. We tend to need
multiple subplots for comparing different datasets and illustrating complex
concepts.  Unfortunately, it is easy to end up with *redundant* figure
elements when drawing multiple subplots, namely...

* ...repeated axis tick labels.
* ...repeated axis labels.
* ...repeated colorbars.
* ...repeated legends.

These sorts of redundancies are very common even in publications, where
they waste valuable page space. They arise because this is often the path of
least resistance -- removing redundancies tends to require extra work.

.. rubric:: Solution

ProPlot seeks to eliminate redundant elements to help you make clear, concise
figures.  We tackle this issue using :ref:`shared and spanning axis labels
<ug_share>` and :ref:`figure-spanning colorbars and legends
<ug_cbars_figure>`.

* Axis tick labels and axis labels are *shared* between subplots in the
  same row or column by default. This is controlled by the `sharex`, `sharey`,
  `spanx`, and `spany` `~proplot.ui.subplots` keyword args.
* The new `~proplot.figure.Figure` `~proplot.figure.Figure.colorbar` and
  `~proplot.figure.Figure.legend` methods make it easy to draw colorbars and
  legends intended to reference more than one subplot. For details, see the
  next section.


.. _why_colorbars_legends:

Outer colorbars and legends
===========================

.. rubric:: Problem

In matplotlib, it can be difficult to draw `~matplotlib.figure.Figure.legend`\ s
along the outside of subplots. Generally, you need to position the legend
manually and adjust various `~matplotlib.gridspec.GridSpec` spacing
properties to make *room* for the legend. And while colorbars can be drawn
along the outside of subplots with ``fig.colorbar(..., ax=ax)``, this can
cause asymmetry in plots with more than one subplot, since the space allocated
for the colorbar is "stolen" from the parent axes.

It is also notoriously difficult to make colorbars that span multiple
subplots in matplotlib. You have to supply
`~matplotlib.figure.Figure.colorbar` with an axes that you drew yourself,
generally using an elaborate `~matplotlib.gridspec.GridSpec` specification.
This requires so much tinkering that most users just add identical colorbars
to every single subplot!


..
   And since colorbar widths are specified in *axes relative* coordinates,
   they often look "too skinny" or "too fat" after the first draw.

..
   The matplotlib example for `~matplotlib.figure.Figure` legends is `not
   pretty
   <https://matplotlib.org/3.1.1/gallery/text_labels_and_annotations/figlegend_demo.html>`__.

..
   Drawing colorbars and legends is pretty clumsy in matplotlib -- especially
   when trying to draw them outside of the figure. They can be too narrow,
   too wide, and mess up your subplot aspect ratios.

.. rubric:: Solution

ProPlot introduces a brand new framework for drawing colorbars and legends
referencing :ref:`individual subplots <ug_cbars_axes>` and
:ref:`multiple contiguous subplots <ug_cbars_figure>`. This framework makes
the process of drawing colorbars and legends much cleaner.

* To draw a colorbar or legend on the outside of a specific subplot, pass an
  "outer" location (e.g. ``loc='l'`` or ``loc='left'``)
  to `proplot.axes.Axes.colorbar` or `proplot.axes.Axes.legend`.
* To draw a colorbar or legend on the inside of a specific subplot, pass an
  "inner" location (e.g. ``loc='ur'`` or ``loc='upper right'``)
  to `proplot.axes.Axes.colorbar` or `proplot.axes.Axes.legend`.
* To draw a colorbar or legend along the edge of the figure, use
  `proplot.figure.Figure.colorbar` and `proplot.figure.Figure.legend`.
  The `col`, `row`, and `span` keyword args control which
  `~matplotlib.gridspec.GridSpec` rows and columns are spanned by the
  colorbar or legend.

Since `~proplot.gridspec.GridSpec` permits variable spacing between subplot
rows and columns, "outer" colorbars and legends do not mess up subplot
spacing or add extra whitespace. This is critical e.g. if you have a
colorbar between columns 1 and 2 but nothing between columns 2 and 3.
Also, `~proplot.figure.Figure` and `~proplot.axes.Axes` colorbar widths are
now specified in *physical* units rather than relative units, which makes
colorbar thickness independent of subplot size and easier to get just right.

There are also several new :ref:`colorbar <ug_cbars>` and
:ref:`legend <ug_legends>` features described in greater detail in the user
guide.


.. _why_plotting:

Enhanced plotting methods
=========================

.. rubric:: Problem

Certain common plotting tasks take a lot of work when using the default
matplotlib API.  The `seaborn`, `xarray`, and `pandas` packages offer
improvements, but it would be nice to have this functionality build right
into matplotlib.  There is also room for improvement of the native matplotlib
plotting methods that none of these packages address.

..
   Matplotlib also has some finicky plotting issues
   that normally requires
..
   For example, when you pass coordinate *centers* to
   `~matplotlib.axes.Axes.pcolor` and `~matplotlib.axes.Axes.pcolormesh`,
   they are interpreted as *edges* and the last column and row of your data
   matrix is ignored. Also, to add labels to `~matplotlib.axes.Axes.contour`
   and `~matplotlib.axes.Axes.contourf`, you need to call a dedicated
   `~matplotlib.axes.Axes.clabel` method instead of just using a keyword
   argument.


.. rubric:: Solution


ProPlot adds various `seaborn`, `xarray`, and `pandas` features to the
`~proplot.axes.Axes` plotting methods along with several *brand new* features
designed to make your life easier.

* The new `~proplot.axes.Axes.heatmap` method invokes
  `~matplotlib.axes.Axes.pcolormesh` and draws ticks at the center of each
  box. This is more convenient for things like covariance matrices.
* The new `~proplot.axes.Axes.parametric` method draws *parametric* line
  plots, where the parametric coordinate is denoted with a colorbar and
  colormap colors rather than text annotations.
* The `~matplotlib.axes.Axes.bar` and `~matplotlib.axes.Axes.barh` methods
  accept 2D arrays and can *stack* or *group* successive columns. Similarly,
  the new `~proplot.axes.Axes.area` and `~proplot.axes.Axes.areax` methods
  (aliases for `~matplotlib.axes.Axes.fill_between` and
  `~matplotlib.axes.Axes.fill_betweenx`) also accept 2D arrays
  and can *stack* or *overlay* successive columns.
* The `~matplotlib.axes.Axes.bar`, `~matplotlib.axes.Axes.barh`,
  `~matplotlib.axes.Axes.vlines`, `~matplotlib.axes.Axes.hlines`,
  `~proplot.axes.Axes.area`, and `~proplot.axes.Axes.areax` commands
  all accept a `negpos` keyword argument that can be used to assign
  "negative" and "positive" colors to different regions.
* You can now :ref:`add error bars or error shading <ug_errorbars>`
  to `~matplotlib.axes.Axes.bar`, `~matplotlib.axes.Axes.barh`,
  `~matplotlib.axes.Axes.plot`, and `~matplotlib.axes.Axes.scatter` plots
  by passing special keyword arguments to these functions. You no longer have
  to work with the `~matplotlib.axes.Axes.errorbar` method directly.
* All :ref:`1D plotting methods <ug_1dplots>` accept a
  :ref:`"cycle" keyword argument <ug_cycle_changer>`
  interpreted by `~proplot.constructor.Cycle` and optional
  :ref:`"colorbar" and "legend" keyword arguments <ug_cbars_axes>`
  for populating legends and colorbars at the specified location with the
  result of the plotting command.
* All :ref:`2D plotting methods <ug_2dplots>` methods accept
  :ref:`"cmap" and "norm" keyword arguments <ug_cmap_changer>`
  interpreted by `~proplot.constructor.Colormap` and
  `~proplot.constructor.Norm`, along with an optional
  :ref:`"colorbar" keyword argument <ug_cbars_axes>`
  for drawing on-the-fly colorbars. They also accept a
  :ref:`"labels" keyword argument <ug_labels>` used to
  draw contour labels or grid box labels on-the-fly, and labels are colored
  black or white according to the luminance of the underlying filled
  contour or grid box color.
* Matplotlib requires coordinate *centers* for contour plots and *edges* for
  pcolor plots. If you pass *centers* to pcolor, matplotlib treats them as
  *edges* and silently trims one row/column of your data. ProPlot
  :ref:`changes this behavior <ug_2dstd>` so that your data is no longer
  trimmed.
* ProPlot fixes an irritating issue with saved vector graphics where white
  lines appear between `filled contours
  <https://stackoverflow.com/q/8263769/4970632>`__, `pcolor patches
  <https://stackoverflow.com/q/27092991/4970632>`__, and `colorbar patches
  <https://stackoverflow.com/q/15003353/4970632>`__.

..
  ProPlot also provides *constistent behavior* when switching between
  different commands, for example `~matplotlib.axes.Axes.plot` and
  `~matplotlib.axes.Axes.scatter` or `~matplotlib.axes.Axes.contourf`
  and `~matplotlib.axes.Axes.pcolormesh`.

..
   ProPlot also uses wrappers to *unify* the behavior of various
   plotting methods.

..
  All positional arguments for 1D plotting methods are standardized by
  `~proplot.axes.standardize_1d`. All positional arguments for 2D
  plotting methods are standardized by `~proplot.axes.standardize_2d`.
  See :ref:`1D plotting methods <1d_plots>` and :ref:`2D plotting methods <2d_plots>`
  for details.


.. _why_xarray_pandas:

Xarray and pandas integration
=============================

.. rubric:: Problem

When you pass the array-like `xarray.DataArray`, `pandas.DataFrame`, and
`pandas.Series` containers to matplotlib plotting commands, the metadata is
ignored. To create plots that are automatically labeled with this metadata,
you must use the dedicated `xarray.DataArray.plot`, `pandas.DataFrame.plot`,
and `pandas.Series.plot` tools instead.

This approach is fine for quick plots, but not ideal for complex ones.  It
requires learning a different syntax from matplotlib, and tends to encourage
using the `~matplotlib.pyplot` API rather than the object-oriented API.
These tools also include features that would be useful additions to matplotlib
in their own right, without requiring special data containers and an entirely
separate API.

.. rubric:: Solution

ProPlot reproduces many of the `xarray.DataArray.plot`,
`pandas.DataFrame.plot`, and `pandas.Series.plot` features on the
`~proplot.axes.Axes` plotting methods themselves.  Passing a
`~xarray.DataArray`, `~pandas.DataFrame`, or `~pandas.Series` through any
plotting method automatically updates the axis tick labels, axis labels,
subplot titles, and colorbar and legend labels from the metadata.  This can
be disabled by passing ``autoformat=False`` to the plotting method or to
`~proplot.ui.subplots`.

Also, as described in the section on :ref:`plotting methods <why_plotting>`,
ProPlot implements certain features like grouped bar plots, layered area
plots, heatmap plots, and on-the-fly colorbars and legends from the `xarray`
and `pandas` APIs directly on the `~proplot.axes.Axes` class.


.. _why_cartopy_basemap:

Cartopy and basemap integration
===============================

.. rubric:: Problem

There are two widely-used engines for plotting geophysical data with
matplotlib: `cartopy` and `~mpl_toolkits.basemap`.  Using cartopy tends to be
verbose and involve boilerplate code, while using basemap requires you to use
plotting commands on a separate `~mpl_toolkits.basemap.Basemap` object rather
than an axes object.

Furthermore, when you use `cartopy` and `~mpl_toolkits.basemap` plotting
commands, the assumed coordinate system is *map projection* coordinates
rather than longitude-latitude coordinates. For many users, this choice is
confusing, since the vast majority of geophysical data are stored in
longitude-latitude or "Plate Carrée" coordinates.

.. rubric:: Solution

ProPlot integrates various `cartopy` and `~mpl_toolkits.basemap` features
into the `proplot.axes.GeoAxes.format` method.  This lets you apply all kinds
of geographic plot settings, like continents, coastlines, political
boundaries, and meridian and parallel gridlines.

`~proplot.axes.GeoAxes` also overrides various plotting methods as follows:

* The new default for all `~proplot.axes.CartopyAxes` plotting methods is
  ``transform=ccrs.PlateCarree()``.
* The new default for all `~proplot.axes.BasemapAxes` plotting methods is
  ``latlon=True``.
* *Global* coverage over the poles and across the matrix longitude boundaries
  can be enforced by passing ``globe=True`` to any 2D plotting command, e.g.
  `~matplotlib.axes.Axes.pcolormesh` and `~matplotlib.axes.Axes.contourf`.

See the :ref:`user guide <ug_proj>` for details.  Note that active
development on basemap will `halt after 2020
<https://matplotlib.org/basemap/users/intro.html#cartopy-new-management-and-eol-announcement>`__.
For now, cartopy is `missing several features
<https://matplotlib.org/basemap/api/basemap_api.html#module-mpl_toolkits.basemap>`__
offered by basemap -- namely, flexible meridian and parallel gridline labels,
drawing physical map scales, and convenience features for adding background
images like the "blue marble". But once these are added to cartopy, ProPlot
may deprecate the `~mpl_toolkits.basemap` integration features.

..
  This is the right decision: Cartopy is integrated more closely with the matplotlib API
  and is more amenable to further development.


.. _why_colormaps_cycles:

Colormaps and property cycles
=============================

.. rubric:: Problem

In matplotlib, colormaps are implemented with the
`~matplotlib.colors.ListedColormap` and
`~matplotlib.colors.LinearSegmentedColormap` classes. They are generally
difficult to edit or create from scratch. The `seaborn` package introduces
"color palettes" to make this easier, but it would be nice to have similar
features built right into the matplotlib API.

..
   Colormap identification is also suboptimal, since the names are case-sensitive, and reversed versions of each colormap are not guaranteed to exist.

.. rubric:: Solution

In ProPlot, it is easy to manipulate colormaps and property cycles:

* The `~proplot.constructor.Colormap` constructor function can be used to
  slice and merge existing colormaps and/or generate brand new ones.
* The `~proplot.constructor.Cycle` constructor function can be used to make
  property cycles from *colormaps*! Property cycles can be applied to plots
  in a variety of ways -- see the :ref:`user guide <ug_cycles>` for details.
* The new `~proplot.colors.ListedColormap` and
  `~proplot.colors.LinearSegmentedColormap` classes include several
  convenient methods and have a much nicer REPL string representation.
* The `~proplot.colors.PerceptuallyUniformColormap` class is used to make
  :ref:`perceptually uniform colormaps <ug_perceptual>`. These have smooth,
  aesthetically pleasing color transitions  that represent your data
  *accurately*.

Importing ProPlot also makes all colormap names *case-insensitive*, and
colormaps can be *reversed* or *cyclically shifted* by 180 degrees simply by
appending ``'_r'`` or ``'_s'`` to the colormap name. This is powered by the
`~proplot.colors.ColormapDatabase` dictionary, which replaces matplotlib's
native colormap database.

.. _why_normalization:

Cleaner colormap normalization
==============================

.. rubric:: Problem

In matplotlib, when ``extend='min'``, ``extend='max'``, or
``extend='neither'`` is passed to `~matplotlib.figure.Figure.colorbar` , the
colormap colors reserved for "out-of-bounds" values are truncated. This can
be irritating for plots with very few colormap levels, which are often more
desirable (see the :ref:`user guide <ug_discrete>` for an example).

The problem is that matplotlib "discretizes" colormaps by generating
low-resolution lookup tables (see
`~matplotlib.colors.LinearSegmentedColormap`).  While straightforward, this
approach has limitations and results in unnecessary plot-specific copies of
the colormap.  Ideally, the task of discretizing colormap colors should be
left to the *normalizer*. While matplotlib provides
`~matplotlib.colors.BoundaryNorm` for this purpose, it is seldom used and
its features are limited.

.. rubric:: Solution

In ProPlot, all colormaps retain a high-resolution lookup table and the
`~proplot.colors.DiscreteNorm` class is used to make distinct "levels".
`~proplot.colors.DiscreteNorm` first normalizes the levels with a
"continuous" normalizer like `~matplotlib.colors.Normalize` or
`~matplotlib.colors.LogNorm`, then restricts the plot to a *subset* of lookup
table colors at the indices of the normalized level centers. It also adjusts
the indices such that the colorbar colors always traverse the full range of
the colormap, regardless of the `extend` setting, and ensures that the end
colors on cyclic colormaps are distinct.

.. _why_container:

The subplot container class
===========================

..
   The `~matplotlib.pyplot.subplots` command is useful for generating a scaffolding of * axes all at once. This is generally faster than successive `~matplotlib.figure.Figure.add_subplot` commands.

.. rubric:: Problem

In matplotlib, `~matplotlib.pyplot.subplots` returns a 2D `~numpy.ndarray`
for figures with more than one column and row, a 1D `~numpy.ndarray` for
single-row or single-column figures, or just an `~matplotlib.axes.Axes`
instance for single-subplot figures.

.. rubric:: Solution

In ProPlot, `~proplot.ui.subplots` returns a `~proplot.ui.SubplotsContainer`
filled with `~proplot.axes.Axes` instances.  This container lets you call
arbitrary methods on arbitrary subplots all at once, which can be useful when
you want to style your subplots identically (e.g.
``axs.format(tickminor=False)``).  The `~proplot.ui.SubplotsContainer` class
also unifies the behavior of the three possible `matplotlib.pyplot.subplots`
return values:

* `~proplot.ui.SubplotsContainer` permits 2D indexing, e.g. ``axs[1,0]``.
  Since `~proplot.ui.subplots` can generate figures with arbitrarily complex
  subplot geometry, this 2D indexing is useful only when the arrangement
  happens to be a clean 2D matrix.
* `~proplot.ui.SubplotsContainer` also permits 1D indexing, e.g. ``axs[0]``.
  The default order can be switched from row-major to column-major by passing
  ``order='F'`` to `~proplot.ui.subplots`.
* When it is singleton, `~proplot.ui.SubplotsContainer` behaves like a
  scalar. So when you make a single axes with ``fig, axs = plot.subplots()``,
  ``axs[0].method(...)`` is equivalent to ``axs.method(...)``.

See the :ref:`user guide <ug_container>` for details.

..
   This goes with ProPlot's theme of preserving the object-oriented spirit,
   but making things easier for users.

.. _why_quick:

Quick global settings
=====================

.. rubric:: Problem

In matplotlib, there are several `~matplotlib.rcParams` that you often want
to set *all at once*, like the tick lengths and spine colors.  It is also
often desirable to change these settings for *individual subplots* rather
than globally.

.. rubric:: Solution

In ProPlot, you can use the `~proplot.config.rc` object to change lots of
settings at once with convenient shorthands.  This is meant to replace
matplotlib's `~matplotlib.rcParams`.  dictionary. Settings can be changed
with ``plot.rc.key = value``, ``plot.rc[key] = value``,
``plot.rc.update(...)``, with the `~proplot.axes.Axes.format` method, or with
the `~proplot.config.RcConfigurator.context` method.

The most notable bulk settings are described below.
See the :ref:`user guide <ug_config>` for details.

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


.. _why_units:

Physical units engine
=====================

.. rubric:: Problem

Matplotlib uses figure-relative units for the margins `left`, `right`,
`bottom`, and `top`, and axes-relative units for the column and row spacing
`wspace` and `hspace`.  Relative units tend to require "tinkering" with
numbers until you find the right one. And since they are *relative*, if you
decide to change your figure size or add a subplot, they will have to be
readjusted.

Matplotlib also requires users to set the figure size `figsize` in inches.
This may be confusing for users outside of the United States.


.. rubric:: Solution

ProPlot introduces the physical units engine `~proplot.utils.units` for
interpreting `figsize`, `width`, `height`, `axwidth`, `axheight`, `left`,
`right`, `top`, `bottom`, `wspace`, `hspace`, and arguments in a few other
places. Acceptable units include inches, centimeters, millimeters, pixels,
`points <https://en.wikipedia.org/wiki/Point_(typography)>`__, `picas
<https://en.wikipedia.org/wiki/Pica_(typography)>`__, `em-heights
<https://en.wikipedia.org/wiki/Em_(typography)>`__, and `light years
<https://en.wikipedia.org/wiki/Light-year>`__.  Em-heights
are particularly useful, as labels already present can be useful *rulers* for
figuring out the amount of space needed.

`~proplot.utils.units` is also used to convert settings passed to
`~proplot.config.rc` from arbitrary physical units to *points* -- for
example, :rcraw:`linewidth`, :rcraw:`ticklen`, :rcraw:`axes.titlesize`, and
:rcraw:`axes.titlepad`.  See the :ref:`user guide <ug_config>` for
details.

.. _why_dotproplot:

The .proplot folder
===================

.. rubric:: Problem

In matplotlib, it can be difficult to design your own colormaps and color
cycles, and there is no builtin way to *save* them for future use. It is also
difficult to get matplotlib to use custom ``.ttc``, ``.ttf``, and ``.otf``
font files, which may be desirable when you are working on Linux servers with
limited font selections.


.. rubric:: Solution

ProPlot automatically adds colormaps, color cycles, and font files saved in
the ``.proplot/cmaps``,  ``.proplot/cycles``, and ``.proplot/fonts`` folders
in your home directory.  You can save colormaps and color cycles to these
folders simply by passing ``save=True`` to `~proplot.constructor.Colormap`
and `~proplot.constructor.Cycle`.  To *manually* load from these folders,
e.g. if you have added files to these folders but you do not want to restart
your ipython session, simply call `~proplot.config.register_cmaps`,
`~proplot.config.register_cycles`, and `~proplot.config.register_fonts`.

..
   As mentioned above, ProPlot introduces the `~proplot.constructor.Colormap`
   and `~proplot.constructor.Cycle` functions for designing your own
   colormaps and color cycles.

..
   ...and much more!
   =================
   This page is not comprehensive -- it just illustrates how ProPlot
   addresses some of the stickiest matplotlib limitations that bug your
   average power user.  See the User Guide for a more comprehensive overview.
